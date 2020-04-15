/*
 * neo_localization_node.cpp
 *
 *  Created on: Apr 8, 2020
 *      Author: mad
 */

#include "../include/GridMap.h"
#include "../include/Solver.h"


#include <ros/ros.h>
#include <angles/angles.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/OccupancyGrid.h>
#include <sensor_msgs/LaserScan.h>
#include <geometry_msgs/TransformStamped.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>

#include <mutex>


template<typename T>
Matrix<T, 4, 4> rotate25_z(T rad) {
	return {std::cos(rad), -std::sin(rad), 0, 0,
			std::sin(rad), std::cos(rad), 0, 0,
			0, 0, 1, rad,
			0, 0, 0, 1};
}

template<typename T>
Matrix<T, 4, 4> rotate3_z(T rad) {
	return {std::cos(rad), -std::sin(rad), 0, 0,
			std::sin(rad), std::cos(rad), 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1};
}

template<typename T>
Matrix<T, 4, 4> translate25(T x, T y) {
	return {1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, 0,
			0, 0, 0, 1};
}

inline
Matrix<double, 4, 4> convert_transform_25(const tf::Transform& trans)
{
	Matrix<double, 4, 4> res;
	res(0, 0) = trans.getBasis()[0][0];
	res(1, 0) = trans.getBasis()[1][0];
	res(0, 1) = trans.getBasis()[0][1];
	res(1, 1) = trans.getBasis()[1][1];
	res(0, 3) = trans.getOrigin()[0];
	res(1, 3) = trans.getOrigin()[1];
	res(2, 3) = tf::getYaw(trans.getRotation());
	res(2, 2) = 1;
	res(3, 3) = 1;
	return res;
}

inline
Matrix<double, 4, 4> convert_transform_3(const tf::Transform& trans)
{
	Matrix<double, 4, 4> res;
	for(int j = 0; j < 3; ++j) {
		for(int i = 0; i < 3; ++i) {
			res(i, j) = trans.getBasis()[i][j];
		}
	}
	res(0, 3) = trans.getOrigin()[0];
	res(1, 3) = trans.getOrigin()[1];
	res(2, 3) = trans.getOrigin()[2];
	res(3, 3) = 1;
	return res;
}


class NeoLocalizationNode {
public:
	NeoLocalizationNode()
	{
		m_node_handle.param("broadcast_tf", m_broadcast_tf, true);
		m_node_handle.param<std::string>("base_frame", m_base_frame, "base_link");
		m_node_handle.param<std::string>("odom_frame", m_odom_frame, "odom");
		m_node_handle.param<std::string>("map_frame", m_map_frame, "map");
		m_node_handle.param("update_gain", m_update_gain, 0.1);
		m_node_handle.param("map_sub_sample", m_map_sub_sample, 1);
		m_node_handle.param("num_smooth", m_num_smooth, 20);
		m_node_handle.param("solver_iterations", m_solver_iterations, 10);
		m_node_handle.param("solver_gain", m_solver.gain, 0.1);
		m_node_handle.param("solver_damping", m_solver.damping, 100.);

		m_sub_scan_topic = m_node_handle.subscribe("/scan", 10, &NeoLocalizationNode::scan_callback, this);
		m_sub_map_topic = m_node_handle.subscribe("/map", 1, &NeoLocalizationNode::map_callback, this);
		m_sub_pose_estimate = m_node_handle.subscribe("/initialpose", 1, &NeoLocalizationNode::pose_callback, this);
	}

protected:
	void scan_callback(const sensor_msgs::LaserScan::ConstPtr& scan)
	{
		std::lock_guard<std::mutex> lock(m_node_mutex);

		if(!m_map_fine || !m_map_coarse) {
			return;
		}

		tf::StampedTransform sensor_to_base;
		try {
			m_tf.lookupTransform(m_base_frame, scan->header.frame_id, scan->header.stamp, sensor_to_base);
		} catch(...) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(scan->header.frame_id, m_base_frame) failed!");
			return;
		}

		tf::StampedTransform base_to_odom;
		try {
			m_tf.lookupTransform(m_odom_frame, m_base_frame, scan->header.stamp, base_to_odom);
		} catch(...) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(m_base_frame, m_odom_frame) failed!");
			return;
		}

		const Matrix<double, 4, 4> S = convert_transform_3(sensor_to_base);
		const Matrix<double, 4, 4> L = convert_transform_25(base_to_odom);
		const Matrix<double, 4, 4> T = translate25(m_offset_x, m_offset_y) * rotate25_z(m_offset_yaw);		// odom to map

		const Matrix<double, 3, 1> grid_pose = (m_grid_to_map.inverse() * T * L * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

		// set initial guess to odometry prediction
		m_solver.pose_x = grid_pose[0];
		m_solver.pose_y = grid_pose[1];
		m_solver.pose_yaw = grid_pose[2];

		std::vector<scan_point_t> points;

		for(size_t i = 0; i < scan->ranges.size(); ++i)
		{
			if(scan->ranges[i] <= 0) {
				continue;	// no measurement
			}

			// transform sensor points into base coordinate system
			const Matrix<double, 3, 1> scan_pos = (S * rotate3_z<double>(scan->angle_min + i * scan->angle_increment)
													* Matrix<double, 4, 1>{scan->ranges[i], 0, 0, 1}).project();
			scan_point_t point;
			point.x = scan_pos[0];
			point.y = scan_pos[1];
			point.w = 1;
			points.emplace_back(point);
		}

		// coarse localization
		for(int iter = 0; iter < m_solver_iterations; ++iter)
		{
			ROS_INFO_STREAM("Coarse iter " << iter << ": pose_x=" << m_solver.pose_x << ", pose_y=" << m_solver.pose_y << ", pose_yaw=" << m_solver.pose_yaw
					<< ", r_norm=" << m_solver.r_norm);

			m_solver.solve<float>(*m_map_coarse, points, 1);
		}

		// fine localization
		for(int iter = 0; iter < m_solver_iterations; ++iter)
		{
			ROS_INFO_STREAM("Fine iter " << iter << ": pose_x=" << m_solver.pose_x << ", pose_y=" << m_solver.pose_y << ", pose_yaw=" << m_solver.pose_yaw
					<< ", r_norm=" << m_solver.r_norm);

			m_solver.solve<float>(*m_map_fine, points, 1);
		}

		// get new pose from solver
		const Matrix<double, 4, 4> grid_pose_new = translate25(m_solver.pose_x, m_solver.pose_y) * rotate25_z(m_solver.pose_yaw);

		// compute new odom to map offset from new pose
		const Matrix<double, 3, 1> new_offset =
				(m_grid_to_map * grid_pose_new * L.inverse() * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

		// apply new offset with an exponential low pass filter
		m_offset_x = new_offset[0] * m_update_gain + m_offset_x * (1 - m_update_gain);
		m_offset_y = new_offset[1] * m_update_gain + m_offset_y * (1 - m_update_gain);
		m_offset_yaw += angles::shortest_angular_distance(m_offset_yaw, new_offset[2]) * m_update_gain;
		m_offset_time = scan->header.stamp;

		// publish new transform
		broadcast();
	}

	void pose_callback(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr& pose)
	{
		std::lock_guard<std::mutex> lock(m_node_mutex);

		if(pose->header.frame_id != m_map_frame) {
			ROS_WARN_STREAM("NeoLocalizationNode: Invalid pose estimate frame: " << pose->header.frame_id);
			return;
		}

		tf::Transform map_pose;
		tf::poseMsgToTF(pose->pose.pose, map_pose);

		tf::StampedTransform base_to_odom;
		try {
			m_tf.lookupTransform(m_odom_frame, m_base_frame, ros::Time(), base_to_odom);
		} catch(...) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(m_base_frame, m_odom_frame) failed!");
			return;
		}

		const Matrix<double, 4, 4> L = convert_transform_25(base_to_odom);

		// compute new odom to map offset
		const Matrix<double, 3, 1> new_offset =
				(convert_transform_25(map_pose) * L.inverse() * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

		m_offset_x = new_offset[0];
		m_offset_y = new_offset[1];
		m_offset_yaw = new_offset[2];
		m_offset_time = pose->header.stamp;

		broadcast();

		ROS_INFO_STREAM("NeoLocalizationNode: Got new pose estimate!");
	}

	void map_callback(const nav_msgs::OccupancyGrid::ConstPtr& ros_map)
	{
		ROS_INFO_STREAM("NeoLocalizationNode: Got map with dimensions " << ros_map->info.width << " x " << ros_map->info.height
				<< " and cell size " << ros_map->info.resolution);

		if(ros_map->info.width != ros_map->info.height)
		{
			ROS_WARN_STREAM("NeoLocalizationNode: Invalid map dimensions!");
			return;
		}

		auto map = std::make_shared<GridMap<float>>(ros_map->info.width, ros_map->info.resolution);

		// convert map to our format (occupancy between 0 and 1)
		for(int y = 0; y < map->size(); ++y) {
			for(int x = 0; x < map->size(); ++x) {
				const auto cell = ros_map->data[y * map->size() + x];
				if(cell >= 0) {
					(*map)(x, y) = fminf(cell / 100.f, 1.f);
				} else {
					(*map)(x, y) = 0;
				}
			}
		}

		auto fine = map;
		auto coarse = map->downscale()->downscale();		// coarse is 4x lower resolution

		for(int i = 0; i < m_num_smooth; ++i) {
			ROS_INFO_STREAM("Smooth iter " << i);
			fine->smooth_33_1();
			coarse->smooth_33_1();
		}

		{
			std::lock_guard<std::mutex> lock(m_node_mutex);
			{
				tf::Transform tmp;
				tf::poseMsgToTF(ros_map->info.origin, tmp);
				m_grid_to_map = convert_transform_25(tmp);
			}
			m_map_fine = fine;
			m_map_coarse = coarse;
		}
	}

	void broadcast()
	{
		if(m_broadcast_tf)
		{
			// compose and publish transform for tf package
			geometry_msgs::TransformStamped pose;
			// compose header
			pose.header.stamp = m_offset_time;
			pose.header.frame_id = m_map_frame;
			pose.child_frame_id = m_odom_frame;
			// compose data container
			pose.transform.translation.x = m_offset_x;
			pose.transform.translation.y = m_offset_y;
			pose.transform.translation.z = 0;
			tf::quaternionTFToMsg(tf::createQuaternionFromYaw(m_offset_yaw), pose.transform.rotation);

			// publish the transform
			m_tf_broadcaster.sendTransform(pose);
		}
	}

private:
	std::mutex m_node_mutex;

	ros::NodeHandle m_node_handle;

	ros::Publisher m_pub_smoothed_map;

	ros::Subscriber m_sub_map_topic;
	ros::Subscriber m_sub_scan_topic;
	ros::Subscriber m_sub_pose_estimate;

	tf::TransformListener m_tf;
	tf::TransformBroadcaster m_tf_broadcaster;

	bool m_broadcast_tf = false;
	std::string m_base_frame;
	std::string m_odom_frame;
	std::string m_map_frame;

	double m_update_gain = 0;
	int m_map_sub_sample = 0;
	int m_num_smooth = 0;
	int m_solver_iterations = 0;

	double m_offset_x = 0;			// current x offset between odom and map
	double m_offset_y = 0;			// current y offset between odom and map
	double m_offset_yaw = 0;		// current yaw offset between odom and map
	ros::Time m_offset_time;

	Matrix<double, 4, 4> m_grid_to_map;

	Solver m_solver;
	std::shared_ptr<GridMap<float>> m_map_fine;
	std::shared_ptr<GridMap<float>> m_map_coarse;

};


int main(int argc, char** argv)
{
	// initialize ROS
	ros::init(argc, argv, "neo_localization_node");

	try {
		NeoLocalizationNode node;

		ros::spin();
	}
	catch(const std::exception& ex) {
		ROS_ERROR_STREAM("NeoOmniDriveNode: " << ex.what());
		return -1;
	}

	return 0;
}

