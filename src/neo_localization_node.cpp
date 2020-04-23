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
#include <geometry_msgs/PoseArray.h>
#include <geometry_msgs/TransformStamped.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>

#include <mutex>
#include <thread>
#include <random>
#include <cmath>
#include <array>


/*
 * Creates a 2.5D (x, y, yaw) rotation matrix.
 */
template<typename T>
Matrix<T, 4, 4> rotate25_z(T rad) {
	return {std::cos(rad), -std::sin(rad), 0, 0,
			std::sin(rad), std::cos(rad), 0, 0,
			0, 0, 1, rad,
			0, 0, 0, 1};
}

/*
 * Creates a 3D rotation matrix for a yaw rotation around Z axis.
 */
template<typename T>
Matrix<T, 4, 4> rotate3_z(T rad) {
	return {std::cos(rad), -std::sin(rad), 0, 0,
			std::sin(rad), std::cos(rad), 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1};
}

/*
 * Creates a 2.5D (x, y, yaw) translation matrix.
 */
template<typename T>
Matrix<T, 4, 4> translate25(T x, T y) {
	return {1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, 0,
			0, 0, 0, 1};
}

/*
 * Converts ROS 3D Transform to a 2.5D matrix.
 */
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

/*
 * Converts ROS 3D Transform to a 3D matrix.
 */
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

/*
 * Computes 1D variance.
 */
template<typename T>
T compute_variance(const std::vector<double>& values, T& mean)
{
	if(values.size() < 2) {
		throw std::logic_error("values.size() < 2");
	}
	mean = 0;
	for(auto v : values) {
		mean += v;
	}
	mean /= T(values.size());

	double var = 0;
	for(auto v : values) {
		var += std::pow(v - mean, 2);
	}
	var /= T(values.size() - 1);
	return var;
}

/*
 * Computes ND covariance matrix.
 */
template<typename T, size_t N, size_t M>
Matrix<T, N, N> compute_covariance(const std::vector<Matrix<T, M, 1>>& points, Matrix<T, N, 1>& mean)
{
	if(M < N) {
		throw std::logic_error("M < N");
	}
	if(points.size() < 2) {
		throw std::logic_error("points.size() < 2");
	}
	mean = Matrix<T, N, 1>();
	for(auto point : points) {
		mean += point.template get<N, 1>();
	}
	mean /= T(points.size());

	Matrix<T, N, N> mat;
	for(auto point : points) {
		for(int j = 0; j < N; ++j) {
			for(int i = 0; i < N; ++i) {
				mat(i, j) += (point[i] - mean[i]) * (point[j] - mean[j]);
			}
		}
	}
	mat /= T(points.size() - 1);
	return mat;
}

/*
 * Computes 1D variance along a 2D axis (given by direction unit vector), around given mean position.
 */
template<typename T, size_t N>
T compute_variance_along_direction_2(	const std::vector<Matrix<T, N, 1>>& points,
										const Matrix<T, 2, 1>& mean,
										const Matrix<T, 2, 1>& direction)
{
	if(N < 2) {
		throw std::logic_error("N < 2");
	}
	if(points.size() < 2) {
		throw std::logic_error("points.size() < 2");
	}
	double var = 0;
	for(auto point : points)
	{
		const auto delta = Matrix<T, 2, 1>{point[0], point[1]} - mean;
		var += std::pow(direction.dot(delta), 2);
	}
	var /= T(points.size() - 1);
	return var;
}

// See: http://croninprojects.org/Vince/Geodesy/FindingEigenvectors.pdf
// See: https://math.stackexchange.com/questions/395698/fast-way-to-calculate-eigen-of-2x2-matrix-using-a-formula
template<typename T>
Matrix<T, 2, 1> compute_eigenvectors_2(	const Matrix<T, 2, 2>& mat,
										std::array<Matrix<T, 2, 1>, 2>& eigen_vectors)
{
	Matrix<T, 2, 1> eigen_values;
	const T tmp_0 = std::sqrt(std::pow(mat(0, 0) + mat(1, 1), T(2)) - T(4) * (mat(0, 0) * mat(1, 1) - mat(1, 0) * mat(0, 1)));
	eigen_values[0] = (mat(0, 0) + mat(1, 1) + tmp_0) / T(2);
	eigen_values[1] = (mat(0, 0) + mat(1, 1) - tmp_0) / T(2);
	if(eigen_values[1] > eigen_values[0]) {
		std::swap(eigen_values[0], eigen_values[1]);
	}

	eigen_vectors[0] = Matrix<T, 2, 1>{mat(0, 1), eigen_values[0] - mat(0, 0)};
	eigen_vectors[1] = Matrix<T, 2, 1>{eigen_values[1] - mat(1, 1), mat(1, 0)};
	eigen_vectors[0] /= eigen_vectors[0].norm();
	eigen_vectors[1] /= eigen_vectors[1].norm();
	return eigen_values;
}


/*
 * Coordinate systems:
 * - Sensor in [meters, rad], aka. "laserX"
 * - Base Link in [meters, rad], aka. "base_link"
 * - Odometry in [meters, rad], aka. "odom"
 * - Map in [meters, rad], aka. "map"
 * - World Grid in [meters, rad], aka. "world"
 * - Tile Grid in [meters, rad], aka. "grid"
 * - World Grid in [pixels]
 * - Tile Grid in [pixels]
 *
 */
class NeoLocalizationNode {
public:
	NeoLocalizationNode()
	{
		m_node_handle.param("broadcast_tf", m_broadcast_tf, true);
		m_node_handle.param<std::string>("base_frame", m_base_frame, "base_link");
		m_node_handle.param<std::string>("odom_frame", m_odom_frame, "odom");
		m_node_handle.param<std::string>("map_frame", m_map_frame, "map");
		m_node_handle.param("map_size", m_map_size, 1000);
		m_node_handle.param("map_downscale", m_map_downscale, 0);
		m_node_handle.param("map_update_rate", m_map_update_rate, 0.5);
		m_node_handle.param("num_smooth", m_num_smooth, 5);
		m_node_handle.param("solver_gain", m_solver.gain, 0.1);
		m_node_handle.param("solver_damping", m_solver.damping, 1000.);
		m_node_handle.param("solver_iterations", m_solver_iterations, 20);
		m_node_handle.param("sample_rate", m_sample_rate, 10);
		m_node_handle.param("min_points", m_min_points, 10);
		m_node_handle.param("update_gain", m_update_gain, 0.5);
		m_node_handle.param("confidence_gain", m_confidence_gain, 0.01);
		m_node_handle.param("confidence_loss", m_confidence_loss, 0.1);
		m_node_handle.param("max_confidence", m_max_confidence, 0.95);
		m_node_handle.param("sample_std_xy", m_sample_std_xy, 0.5);
		m_node_handle.param("sample_std_yaw", m_sample_std_yaw, 0.5);
		m_node_handle.param("constrain_ratio", m_constrain_ratio, 0.2);
		m_node_handle.param("constrain_threshold", m_constrain_threshold, 0.3);
		m_node_handle.param("disable_threshold", m_disable_threshold, 1.);

		m_sub_scan_topic = m_node_handle.subscribe("/scan", 10, &NeoLocalizationNode::scan_callback, this);
		m_sub_map_topic = m_node_handle.subscribe("/map", 1, &NeoLocalizationNode::map_callback, this);
		m_sub_pose_estimate = m_node_handle.subscribe("/initialpose", 1, &NeoLocalizationNode::pose_callback, this);

		m_pub_map_tile = m_node_handle.advertise<nav_msgs::OccupancyGrid>("/map_tile", 1);
		m_pub_loc_pose = m_node_handle.advertise<geometry_msgs::PoseWithCovarianceStamped>("/amcl_pose", 10);
		m_pub_loc_pose_2 = m_node_handle.advertise<geometry_msgs::PoseWithCovarianceStamped>("/map_pose", 10);
		m_pub_pose_array = m_node_handle.advertise<geometry_msgs::PoseArray>("/particlecloud", 10);

		m_update_thread = std::thread(&NeoLocalizationNode::update_loop, this);
	}

	~NeoLocalizationNode()
	{
		if(m_update_thread.joinable()) {
			m_update_thread.join();
		}
	}

protected:
	/*
	 * Computes localization update for a single laser scan.
	 */
	void scan_callback(const sensor_msgs::LaserScan::ConstPtr& scan)
	{
		std::lock_guard<std::mutex> lock(m_node_mutex);
		if(!m_map) {
			return;
		}

		tf::StampedTransform sensor_to_base;
		try {
			m_tf.waitForTransform(m_base_frame, scan->header.frame_id, scan->header.stamp, ros::Duration(0.2));
			m_tf.lookupTransform(m_base_frame, scan->header.frame_id, scan->header.stamp, sensor_to_base);
		} catch(const std::exception& ex) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(scan->header.frame_id, m_base_frame) failed: " << ex.what());
			return;
		}

		tf::StampedTransform base_to_odom;
		try {
			m_tf.waitForTransform(m_odom_frame, m_base_frame, scan->header.stamp, ros::Duration(0.2));
			m_tf.lookupTransform(m_odom_frame, m_base_frame, scan->header.stamp, base_to_odom);
		} catch(const std::exception& ex) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(m_base_frame, m_odom_frame) failed: " << ex.what());
			return;
		}

		const Matrix<double, 4, 4> S = convert_transform_3(sensor_to_base);
		const Matrix<double, 4, 4> L = convert_transform_25(base_to_odom);
		const Matrix<double, 4, 4> T = translate25(m_offset_x, m_offset_y) * rotate25_z(m_offset_yaw);		// odom to map

		const Matrix<double, 3, 1> odom_pose = (L * Matrix<double, 4, 1>{0, 0, 0, 1}).project();
		const double dist_moved = (odom_pose - m_last_odom_pose).get<2, 1>().norm();

		std::vector<scan_point_t> points;

		for(size_t i = 0; i < scan->ranges.size(); ++i)
		{
			if(scan->ranges[i] <= scan->range_min || scan->ranges[i] >= scan->range_max) {
				continue;	// no actual measurement
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

		// check for number of points
		if(points.size() < m_min_points)
		{
			ROS_WARN_STREAM("NeoLocalizationNode: Number of points too low: " << points.size());
			return;
		}

		auto pose_array = boost::make_shared<geometry_msgs::PoseArray>();
		pose_array->header.stamp = scan->header.stamp;
		pose_array->header.frame_id = m_map_frame;

		// calc predicted grid pose based on odometry
		const Matrix<double, 3, 1> grid_pose = (m_grid_to_map.inverse() * T * L * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

		// setup distributions
		const double rel_std_dev = fmax(1 - m_confidence, 0);
		std::normal_distribution<double> dist_x(grid_pose[0], m_sample_std_xy * rel_std_dev);
		std::normal_distribution<double> dist_y(grid_pose[1], m_sample_std_xy * rel_std_dev);
		std::normal_distribution<double> dist_yaw(grid_pose[2], m_sample_std_yaw * rel_std_dev);

		// solve odometry prediction first
		m_solver.pose_x = grid_pose[0];
		m_solver.pose_y = grid_pose[1];
		m_solver.pose_yaw = grid_pose[2];

		for(int iter = 0; iter < m_solver_iterations; ++iter) {
			m_solver.solve<float>(*m_map, points, 1);
		}

		double best_x = m_solver.pose_x;
		double best_y = m_solver.pose_y;
		double best_yaw = m_solver.pose_yaw;
		double best_score = m_solver.r_norm;

		std::vector<Matrix<double, 4, 1>> seeds(m_sample_rate);
		std::vector<Matrix<double, 4, 1>> samples(m_sample_rate);
		std::vector<double> sample_errors(m_sample_rate);

		for(int i = 0; i < m_sample_rate; ++i)
		{
			// generate new sample
			m_solver.pose_x = dist_x(m_generator);
			m_solver.pose_y = dist_y(m_generator);
			m_solver.pose_yaw = dist_yaw(m_generator);

			seeds[i] = Matrix<double, 4, 1>{m_solver.pose_x, m_solver.pose_y, m_solver.pose_yaw, 1};

			// solve sample
			for(int iter = 0; iter < m_solver_iterations; ++iter) {
				m_solver.solve<float>(*m_map, points, 1);
			}

			// save sample
			const auto sample = Matrix<double, 4, 1>{m_solver.pose_x, m_solver.pose_y, m_solver.pose_yaw, 1};
			samples[i] = sample;
			sample_errors[i] = m_solver.r_norm;

			// check if sample is better
			if(m_solver.r_norm < best_score) {
				best_x = m_solver.pose_x;
				best_y = m_solver.pose_y;
				best_yaw = m_solver.pose_yaw;
				best_score = m_solver.r_norm;
			}

			// add to visualization
			{
				const Matrix<double, 3, 1> map_pose = (m_grid_to_map * sample).project();
				tf::Pose pose;
				pose.setOrigin(tf::Vector3(map_pose[0], map_pose[1], 0));
				pose.setRotation(tf::createQuaternionFromYaw(map_pose[2]));
				geometry_msgs::Pose tmp;
				tf::poseTFToMsg(pose, tmp);
				pose_array->poses.push_back(tmp);
			}
		}

		// compute covariances
		double mean_error = 0;
		double mean_yaw = 0;
		Matrix<double, 3, 1> mean_xyw;
		Matrix<double, 3, 1> seed_mean_xyw;
		const double var_error = compute_variance(sample_errors, mean_error);
		const Matrix<double, 3, 3> var_xyw = compute_covariance(samples, mean_xyw);
		const Matrix<double, 3, 3> seed_var_xyw = compute_covariance(seeds, seed_mean_xyw);

		// compute "estimated" error characteristic
		std::array<Matrix<double, 2, 1>, 2> eigen_vectors;
		const Matrix<double, 2, 1> eigen_values = compute_eigenvectors_2(var_xyw.get<2, 2>(), eigen_vectors);
		const Matrix<double, 2, 1> sigma_uv {sqrt(fabs(eigen_values[0])), sqrt(fabs(eigen_values[1]))};

		// compute seed variance along "estimated" eigen vectors
		const double seed_sigma_u = sqrt(compute_variance_along_direction_2(seeds, seed_mean_xyw.get<2, 1>(), eigen_vectors[0]));
		const double seed_sigma_v = sqrt(compute_variance_along_direction_2(seeds, seed_mean_xyw.get<2, 1>(), eigen_vectors[1]));

		// decide if we have 2D, 1D or 0D localization
		const double factor_0 = fmax(sqrt(var_error) / 0.01, 0.1);
		const double factor_1 = sigma_uv[0] / seed_sigma_u;
		const double factor_2 = sigma_uv[1] / seed_sigma_v;

		int mode = -1;
		if(factor_2 / factor_0 > m_disable_threshold) {
			mode = 0;		// high sigma_v plus low var_error indicates no error gradient = nothing to localize with
		} else if(factor_1 > m_constrain_threshold && factor_2 < m_constrain_ratio * m_constrain_threshold) {
			mode = 1;		// sigma imbalance (above the threshold) means we can/should only localize in one direction
		} else {
			mode = 2;		// all good
		}

		if(mode > 0)
		{
			double new_grid_x = best_x;
			double new_grid_y = best_y;

			if(mode == 1)
			{
				// constrain update to the good direction (ie. in direction of the eigen vector with the smaller sigma)
				const auto delta = Matrix<double, 2, 1>{best_x, best_y} - Matrix<double, 2, 1>{grid_pose[0], grid_pose[1]};
				const auto dist = eigen_vectors[1].dot(delta);
				new_grid_x = grid_pose[0] + dist * eigen_vectors[1][0];
				new_grid_y = grid_pose[1] + dist * eigen_vectors[1][1];
			}

			// use best sample for update
			Matrix<double, 4, 4> grid_pose_new = translate25(new_grid_x, new_grid_y) * rotate25_z(best_yaw);

			// compute new odom to map offset from new grid pose
			const Matrix<double, 3, 1> new_offset =
					(m_grid_to_map * grid_pose_new * L.inverse() * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

			// apply new offset with an exponential low pass filter
			const double gain_factor = double(points.size()) / scan->ranges.size();
			const double gain = fmin(m_update_gain * gain_factor, 1);
			m_offset_x += (new_offset[0] - m_offset_x) * gain;
			m_offset_y += (new_offset[1] - m_offset_y) * gain;
			m_offset_yaw += angles::shortest_angular_distance(m_offset_yaw, new_offset[2]) * gain;

			if(scan->header.stamp > m_offset_time || (m_offset_time - scan->header.stamp).toSec() > 1) {
				m_offset_time = scan->header.stamp;
			}

			if(mode > 1) {
				// apply time based confidence gain
				m_confidence += (m_max_confidence - m_confidence) * m_confidence_gain * gain_factor;
			} else {
				// apply confidence loss based on distance moved
				m_confidence = fmax(m_confidence - dist_moved * m_confidence_loss, 0);
			}
		}

		// publish new transform
		broadcast();

		const Matrix<double, 3, 1> new_map_pose = (translate25(m_offset_x, m_offset_y) * rotate25_z(m_offset_yaw) *
													L * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

		// publish localization pose
		auto loc_pose = boost::make_shared<geometry_msgs::PoseWithCovarianceStamped>();
		loc_pose->header.stamp = m_offset_time;
		loc_pose->header.frame_id = m_map_frame;
		loc_pose->pose.pose.position.x = new_map_pose[0];
		loc_pose->pose.pose.position.y = new_map_pose[1];
		loc_pose->pose.pose.position.z = 0;
		tf::quaternionTFToMsg(tf::createQuaternionFromYaw(new_map_pose[2]), loc_pose->pose.pose.orientation);
		for(int j = 0; j < 3; ++j) {
			for(int i = 0; i < 3; ++i) {
				const int i_ = (i == 2 ? 5 : i);
				const int j_ = (j == 2 ? 5 : j);
				loc_pose->pose.covariance[j_ * 6 + i_] = var_xyw(i, j);
			}
		}
		m_pub_loc_pose.publish(loc_pose);
		m_pub_loc_pose_2.publish(loc_pose);

		// publish visualization
		m_pub_pose_array.publish(pose_array);

		// keep last odom pose
		m_last_odom_pose = odom_pose;

		ROS_INFO_STREAM("NeoLocalizationNode: r_norm=" << best_score << ", sigma_uv=[" << sigma_uv[0] << ", " << sigma_uv[1]
				<< "] m, confidence=" << m_confidence << ", mode=" << mode << "D");
	}

	/*
	 * Resets localization to given position.
	 */
	void pose_callback(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr& pose)
	{
		{
			std::lock_guard<std::mutex> lock(m_node_mutex);

			if(pose->header.frame_id != m_map_frame) {
				ROS_WARN_STREAM("NeoLocalizationNode: Invalid pose estimate frame: " << pose->header.frame_id);
				return;
			}

			tf::Transform map_pose;
			tf::poseMsgToTF(pose->pose.pose, map_pose);

			ROS_INFO_STREAM("NeoLocalizationNode: Got new map pose estimate: x=" << map_pose.getOrigin()[0]
							<< ", y=" <<  map_pose.getOrigin()[1] << ", yaw=" << tf::getYaw(map_pose.getRotation()));

			tf::StampedTransform base_to_odom;
			try {
				m_tf.lookupTransform(m_odom_frame, m_base_frame, ros::Time(), base_to_odom);
			} catch(const std::exception& ex) {
				ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(m_base_frame, m_odom_frame) failed: " << ex.what());
				return;
			}

			const Matrix<double, 4, 4> L = convert_transform_25(base_to_odom);

			// compute new odom to map offset
			const Matrix<double, 3, 1> new_offset =
					(convert_transform_25(map_pose) * L.inverse() * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

			// set new offset based on given position
			m_offset_x = new_offset[0];
			m_offset_y = new_offset[1];
			m_offset_yaw = new_offset[2];

			// reset confidence to zero
			m_confidence = 0;

			broadcast();
		}

		// get a new map tile immediately
		update_map();
	}

	/*
	 * Stores the given map.
	 */
	void map_callback(const nav_msgs::OccupancyGrid::ConstPtr& ros_map)
	{
		std::lock_guard<std::mutex> lock(m_node_mutex);

		ROS_INFO_STREAM("NeoLocalizationNode: Got new map with dimensions " << ros_map->info.width << " x " << ros_map->info.height
				<< " and cell size " << ros_map->info.resolution);

		{
			tf::Transform tmp;
			tf::poseMsgToTF(ros_map->info.origin, tmp);
			m_world_to_map = convert_transform_25(tmp);
		}
		m_world = ros_map;
		m_confidence = 0;
	}

	/*
	 * Extracts a new map tile around current position.
	 */
	void update_map()
	{
		Matrix<double, 4, 4> world_to_map;			// transformation from original grid map (integer coords) to "map frame"
		Matrix<double, 3, 1> world_pose;			// pose in the original (integer coords) grid map (not map tile)
		tf::StampedTransform base_to_odom;
		nav_msgs::OccupancyGrid::ConstPtr world;
		{
			std::lock_guard<std::mutex> lock(m_node_mutex);
			if(!m_world) {
				return;
			}

			try {
				m_tf.lookupTransform(m_odom_frame, m_base_frame, ros::Time(), base_to_odom);
			} catch(const std::exception& ex) {
				ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(m_base_frame, m_odom_frame) failed: " << ex.what());
				return;
			}

			const Matrix<double, 4, 4> L = convert_transform_25(base_to_odom);
			const Matrix<double, 4, 4> T = translate25(m_offset_x, m_offset_y) * rotate25_z(m_offset_yaw);		// odom to map
			world_pose = (m_world_to_map.inverse() * T * L * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

			world = m_world;
			world_to_map = m_world_to_map;
		}

		// compute tile origin in pixel coords
		const double world_scale = world->info.resolution;
		const int tile_x = int(world_pose[0] / world_scale) - m_map_size / 2;
		const int tile_y = int(world_pose[1] / world_scale) - m_map_size / 2;

		auto map = std::make_shared<GridMap<float>>(m_map_size, world_scale);

		// extract tile and convert to our format (occupancy between 0 and 1)
		for(int y = 0; y < map->size(); ++y) {
			for(int x = 0; x < map->size(); ++x) {
				const int x_ = std::min(std::max(tile_x + x, 0), int(world->info.width) - 1);
				const int y_ = std::min(std::max(tile_y + y, 0), int(world->info.height) - 1);
				const auto cell = world->data[y_ * world->info.width + x_];
				if(cell >= 0) {
					(*map)(x, y) = fminf(cell / 100.f, 1.f);
				} else {
					(*map)(x, y) = 0;
				}
			}
		}

		// optionally downscale map
		for(int i = 0; i < m_map_downscale; ++i) {
			map = map->downscale();
		}

		// smooth map
		for(int i = 0; i < m_num_smooth; ++i) {
			map->smooth_33_1();
		}

		// update map
		{
			std::lock_guard<std::mutex> lock(m_node_mutex);
			m_map = map;
			m_grid_to_map = world_to_map * translate25<double>(tile_x * world_scale, tile_y * world_scale);
		}

		const auto tile_origin = (m_grid_to_map * Matrix<double, 4, 1>{0, 0, 0, 1}).project();
		const auto tile_center = (m_grid_to_map * Matrix<double, 4, 1>{	map->scale() * m_map_size / 2,
																		map->scale() * m_map_size / 2, 0, 1}).project();

		// publish new map tile for visualization
		auto ros_grid = boost::make_shared<nav_msgs::OccupancyGrid>();
		ros_grid->info.resolution = map->scale();
		ros_grid->info.width = map->size();
		ros_grid->info.height = map->size();
		ros_grid->info.origin.position.x = tile_origin[0];
		ros_grid->info.origin.position.y = tile_origin[1];
		tf::quaternionTFToMsg(tf::createQuaternionFromYaw(tile_origin[2]), ros_grid->info.origin.orientation);
		ros_grid->data.resize(map->num_cells());
		for(int y = 0; y < map->size(); ++y) {
			for(int x = 0; x < map->size(); ++x) {
				ros_grid->data[y * map->size() + x] = (*map)(x, y) * 100.f;
			}
		}
		m_pub_map_tile.publish(ros_grid);

		ROS_INFO_STREAM("NeoLocalizationNode: Got new grid at offset (" << tile_x << ", " << tile_y << ") [iworld], "
				"center = (" << tile_center[0] << ", " << tile_center[1] << ") [map]");
	}

	/*
	 * Asynchronous map update loop, running in separate thread.
	 */
	void update_loop()
	{
		ros::Rate rate(m_map_update_rate);
		while(ros::ok()) {
			try {
				update_map();	// get a new map tile periodically
			}
			catch(const std::exception& ex) {
				ROS_WARN_STREAM("NeoLocalizationNode: update_map() failed: " << ex.what());
			}
			rate.sleep();
		}
	}

	/*
	 * Publishes "map" frame on tf.
	 */
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

	ros::Publisher m_pub_map_tile;
	ros::Publisher m_pub_loc_pose;
	ros::Publisher m_pub_loc_pose_2;
	ros::Publisher m_pub_pose_array;

	ros::Subscriber m_sub_map_topic;
	ros::Subscriber m_sub_scan_topic;
	ros::Subscriber m_sub_pose_estimate;

	tf::TransformListener m_tf;
	tf::TransformBroadcaster m_tf_broadcaster;

	bool m_broadcast_tf = false;
	std::string m_base_frame;
	std::string m_odom_frame;
	std::string m_map_frame;

	int m_map_size = 0;
	int m_map_downscale = 0;
	int m_num_smooth = 0;
	int m_solver_iterations = 0;
	int m_sample_rate = 0;
	int m_min_points = 0;
	double m_update_gain = 0;
	double m_confidence_gain = 0;
	double m_confidence_loss = 0;
	double m_max_confidence = 0;
	double m_sample_std_xy = 0;
	double m_sample_std_yaw = 0;
	double m_constrain_ratio = 0;
	double m_constrain_threshold = 0;
	double m_disable_threshold = 0;
	double m_map_update_rate = 0;

	double m_offset_x = 0;			// current x offset between odom and map
	double m_offset_y = 0;			// current y offset between odom and map
	double m_offset_yaw = 0;		// current yaw offset between odom and map
	double m_confidence = 0;		// current localization confidence
	ros::Time m_offset_time;

	Matrix<double, 3, 1> m_last_odom_pose;
	Matrix<double, 4, 4> m_grid_to_map;
	Matrix<double, 4, 4> m_world_to_map;
	std::shared_ptr<GridMap<float>> m_map;			// map tile
	nav_msgs::OccupancyGrid::ConstPtr m_world;		// whole map

	Solver m_solver;
	std::mt19937 m_generator;
	std::thread m_update_thread;

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

