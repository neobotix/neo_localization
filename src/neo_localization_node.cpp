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
#include <sensor_msgs/LaserScan.h>
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
Matrix<T, 4, 4> translate25(T x, T y) {
	return {1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, 0,
			0, 0, 0, 1};
}

inline
Matrix<double, 4, 4> convert_transform(const tf::Transform& trans)
{
	Matrix<double, 4, 4> res;
	res(0, 0) = trans.getBasis()[0];
	res(1, 0) = trans.getBasis()[1];
	res(0, 1) = trans.getBasis()[3];
	res(1, 1) = trans.getBasis()[4];
	res(0, 3) = trans.getOrigin()[0];
	res(1, 3) = trans.getOrigin()[1];
	res(2, 3) = tf::getYaw(trans.getRotation());
	res(2, 2) = 1;
	res(3, 3) = 1;
	return res;
}


class NeoLocalizationNode {
public:
	NeoLocalizationNode()
	{
		// TODO
	}

protected:
	void scan_callback(const sensor_msgs::LaserScan::ConstPtr& scan)
	{
		std::lock_guard<std::mutex> lock(m_node_mutex);

		if(!m_map) {
			return;
		}

		tf::StampedTransform sensor_to_base;
		try {
			m_tf.lookupTransform(scan->header.frame_id, m_base_frame, scan->header.stamp, sensor_to_base);
		} catch(...) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(scan->header.frame_id, m_base_frame) failed!");
			return;
		}

		tf::StampedTransform base_to_odom;
		try {
			m_tf.lookupTransform(m_base_frame, m_odom_frame, scan->header.stamp, base_to_odom);
		} catch(...) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(m_base_frame, m_odom_frame) failed!");
			return;
		}

		const Matrix<double, 4, 4> L = convert_transform(base_to_odom);
		const Matrix<double, 4, 4> T = translate25(m_offset_x, m_offset_y) * rotate25_z(m_offset_yaw);

		const Matrix<double, 3, 1> map_pose = (T * L * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

		m_solver.pos_x = map_pose[0];
		m_solver.pos_y = map_pose[1];
		m_solver.pos_yaw = map_pose[2];

		std::vector<scan_point_t> points(scan->ranges.size());

		for(size_t i = 0; i < points.size(); ++i)
		{
			// TODO: convert points to base link
		}

		for(int iter = 0; iter < m_num_iterations; ++iter)
		{
			m_solver.solve<float>(*m_map, points, 1);
		}

		const Matrix<double, 3, 1> new_offset =
				(L.inverse() * Matrix<double, 4, 1>{m_solver.pos_x, m_solver.pos_y, m_solver.pos_yaw, 1}).project();

		const double delta_yaw = angles::shortest_angular_distance(m_offset_yaw, new_offset[2]);

		m_offset_x = new_offset[0] * m_update_gain + m_offset_x * (1 - m_update_gain);
		m_offset_y = new_offset[1] * m_update_gain + m_offset_y * (1 - m_update_gain);
		m_offset_yaw += delta_yaw * m_update_gain;
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

	double m_update_gain = 0;
	int m_sub_sample = 1;
	int m_num_iterations = 0;

	double m_offset_x = 0;			// current x offset between odom and map
	double m_offset_y = 0;			// current y offset between odom and map
	double m_offset_yaw = 0;		// current yaw offset between odom and map

	Solver m_solver;
	std::shared_ptr<GridMap<float>> m_map;

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
	}

	return 0;
}

