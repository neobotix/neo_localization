/*
 * neo_localization_node.cpp
 *
 *  Created on: Apr 8, 2020
 *      Author: mad
 */

#include <neo_localization/Util.h>
#include <neo_localization/Convert.h>
#include <neo_localization/Solver.h>
#include <neo_localization/GridMap.h>

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
		m_node_handle.param("loc_update_rate", m_loc_update_rate, 10.);
		m_node_handle.param("num_smooth", m_num_smooth, 5);
		m_node_handle.param("min_score", m_min_score, 0.2);
		m_node_handle.param("solver_gain", m_solver.gain, 0.1);
		m_node_handle.param("solver_damping", m_solver.damping, 1000.);
		m_node_handle.param("solver_iterations", m_solver_iterations, 20);
		m_node_handle.param("sample_rate", m_sample_rate, 10);
		m_node_handle.param("min_points", m_min_points, 20);
		m_node_handle.param("update_gain", m_update_gain, 0.5);
		m_node_handle.param("confidence_gain", m_confidence_gain, 0.01);
		m_node_handle.param("odometry_std_xy", m_odometry_std_xy, 0.01);
		m_node_handle.param("odometry_std_yaw", m_odometry_std_yaw, 0.01);
		m_node_handle.param("min_sample_std_xy", m_min_sample_std_xy, 0.025);
		m_node_handle.param("min_sample_std_yaw", m_min_sample_std_yaw, 0.025);
		m_node_handle.param("max_sample_std_xy", m_max_sample_std_xy, 0.5);
		m_node_handle.param("max_sample_std_yaw", m_max_sample_std_yaw, 0.5);
		m_node_handle.param("constrain_threshold", m_constrain_threshold, 0.1);
		m_node_handle.param("constrain_threshold_yaw", m_constrain_threshold_yaw, 0.2);
		m_node_handle.param("transform_timeout", m_transform_timeout, 0.2);

		m_sub_scan_topic = m_node_handle.subscribe("/scan", 10, &NeoLocalizationNode::scan_callback, this);
		m_sub_map_topic = m_node_handle.subscribe("/map", 1, &NeoLocalizationNode::map_callback, this);
		m_sub_pose_estimate = m_node_handle.subscribe("/initialpose", 1, &NeoLocalizationNode::pose_callback, this);

		m_pub_map_tile = m_node_handle.advertise<nav_msgs::OccupancyGrid>("/map_tile", 1);
		m_pub_loc_pose = m_node_handle.advertise<geometry_msgs::PoseWithCovarianceStamped>("/amcl_pose", 10);
		m_pub_loc_pose_2 = m_node_handle.advertise<geometry_msgs::PoseWithCovarianceStamped>("/map_pose", 10);
		m_pub_pose_array = m_node_handle.advertise<geometry_msgs::PoseArray>("/particlecloud", 10);

		m_loc_update_timer = m_node_handle.createTimer(	ros::Rate(m_loc_update_rate),
														&NeoLocalizationNode::loc_update, this);

		m_map_update_thread = std::thread(&NeoLocalizationNode::update_loop, this);
	}

	~NeoLocalizationNode()
	{
		if(m_map_update_thread.joinable()) {
			m_map_update_thread.join();
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
		m_scan_buffer[scan->header.frame_id] = scan;
	}

	/*
	 * Convert/Transform a scan from ROS format to a specified base frame.
	 */
	std::vector<scan_point_t> convert_scan(const sensor_msgs::LaserScan::ConstPtr& scan, const Matrix<double, 4, 4>& odom_to_base)
	{
		std::vector<scan_point_t> points;

		tf::StampedTransform sensor_to_base;
		try {
			m_tf.lookupTransform(m_base_frame, scan->header.frame_id, ros::Time(0), sensor_to_base);
		} catch(const std::exception& ex) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(scan->header.frame_id, m_base_frame) failed: " << ex.what());
			return points;
		}

		tf::StampedTransform base_to_odom;
		try {
			m_tf.waitForTransform(m_odom_frame, m_base_frame, scan->header.stamp, ros::Duration(m_transform_timeout));
			m_tf.lookupTransform(m_odom_frame, m_base_frame, scan->header.stamp, base_to_odom);
		} catch(const std::exception& ex) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(m_base_frame, m_odom_frame) failed: " << ex.what());
			return points;
		}

		const Matrix<double, 4, 4> S = convert_transform_3(sensor_to_base);
		const Matrix<double, 4, 4> L = convert_transform_25(base_to_odom);

		// precompute transformation matrix from sensor to requested base
		const Matrix<double, 4, 4> T = odom_to_base * L * S;

		for(size_t i = 0; i < scan->ranges.size(); ++i)
		{
			if(scan->ranges[i] <= scan->range_min || scan->ranges[i] >= scan->range_max) {
				continue;	// no actual measurement
			}

			// transform sensor points into base coordinate system
			const Matrix<double, 3, 1> scan_pos = (T * rotate3_z<double>(scan->angle_min + i * scan->angle_increment)
													* Matrix<double, 4, 1>{scan->ranges[i], 0, 0, 1}).project();
			scan_point_t point;
			point.x = scan_pos[0];
			point.y = scan_pos[1];
			points.emplace_back(point);
		}
		return points;
	}

	void loc_update(const ros::TimerEvent& event)
	{
		std::lock_guard<std::mutex> lock(m_node_mutex);

		if(!m_map || m_scan_buffer.empty()) {
			return;
		}

		tf::StampedTransform base_to_odom;
		try {
			m_tf.lookupTransform(m_odom_frame, m_base_frame, ros::Time(0), base_to_odom);
		} catch(const std::exception& ex) {
			ROS_WARN_STREAM("NeoLocalizationNode: lookupTransform(m_base_frame, m_odom_frame) failed: " << ex.what());
			return;
		}

		const Matrix<double, 4, 4> L = convert_transform_25(base_to_odom);
		const Matrix<double, 4, 4> T = translate25(m_offset_x, m_offset_y) * rotate25_z(m_offset_yaw);		// odom to map

		const Matrix<double, 3, 1> odom_pose = (L * Matrix<double, 4, 1>{0, 0, 0, 1}).project();
		const double dist_moved = (odom_pose - m_last_odom_pose).get<2>().norm();
		const double rad_rotated = fabs(angles::normalize_angle(odom_pose[2] - m_last_odom_pose[2]));

		std::vector<scan_point_t> points;

		// convert all scans to current base frame
		for(const auto& scan : m_scan_buffer)
		{
			auto scan_points = convert_scan(scan.second, L.inverse());
			points.insert(points.end(), scan_points.begin(), scan_points.end());
		}

		// check for number of points
		if(points.size() < m_min_points)
		{
			ROS_WARN_STREAM("NeoLocalizationNode: Number of points too low: " << points.size());
			return;
		}

		auto pose_array = boost::make_shared<geometry_msgs::PoseArray>();
		pose_array->header.stamp = base_to_odom.stamp_;
		pose_array->header.frame_id = m_map_frame;

		// calc predicted grid pose based on odometry
		const Matrix<double, 3, 1> grid_pose = (m_grid_to_map.inverse() * T * L * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

		// setup distributions
		std::normal_distribution<double> dist_x(grid_pose[0], m_sample_std_xy);
		std::normal_distribution<double> dist_y(grid_pose[1], m_sample_std_xy);
		std::normal_distribution<double> dist_yaw(grid_pose[2], m_sample_std_yaw);

		// solve odometry prediction first
		m_solver.pose_x = grid_pose[0];
		m_solver.pose_y = grid_pose[1];
		m_solver.pose_yaw = grid_pose[2];

		for(int iter = 0; iter < m_solver_iterations; ++iter) {
			m_solver.solve<float>(*m_map, points);
		}

		double best_x = m_solver.pose_x;
		double best_y = m_solver.pose_y;
		double best_yaw = m_solver.pose_yaw;
		double best_score = m_solver.r_norm;

		std::vector<Matrix<double, 3, 1>> seeds(m_sample_rate);
		std::vector<Matrix<double, 3, 1>> samples(m_sample_rate);
		std::vector<double> sample_errors(m_sample_rate);

		for(int i = 0; i < m_sample_rate; ++i)
		{
			// generate new sample
			m_solver.pose_x = dist_x(m_generator);
			m_solver.pose_y = dist_y(m_generator);
			m_solver.pose_yaw = dist_yaw(m_generator);

			seeds[i] = Matrix<double, 3, 1>{m_solver.pose_x, m_solver.pose_y, m_solver.pose_yaw};

			// solve sample
			for(int iter = 0; iter < m_solver_iterations; ++iter) {
				m_solver.solve<float>(*m_map, points);
			}

			// save sample
			const auto sample = Matrix<double, 3, 1>{m_solver.pose_x, m_solver.pose_y, m_solver.pose_yaw};
			samples[i] = sample;
			sample_errors[i] = m_solver.r_norm;

			// check if sample is better
			if(m_solver.r_norm > best_score) {
				best_x = m_solver.pose_x;
				best_y = m_solver.pose_y;
				best_yaw = m_solver.pose_yaw;
				best_score = m_solver.r_norm;
			}

			// add to visualization
			{
				const Matrix<double, 3, 1> map_pose = (m_grid_to_map * sample.extend()).project();
				tf::Pose pose;
				pose.setOrigin(tf::Vector3(map_pose[0], map_pose[1], 0));
				pose.setRotation(tf::createQuaternionFromYaw(map_pose[2]));
				geometry_msgs::Pose tmp;
				tf::poseTFToMsg(pose, tmp);
				pose_array->poses.push_back(tmp);
			}
		}

		// compute covariances
		double mean_score = 0;
		Matrix<double, 3, 1> mean_xyw;
		Matrix<double, 3, 1> seed_mean_xyw;
		const double var_error = compute_variance(sample_errors, mean_score);
		const Matrix<double, 3, 3> var_xyw = compute_covariance(samples, mean_xyw);
		const Matrix<double, 3, 3> grad_var_xyw =
				compute_virtual_scan_covariance_xyw(m_map, points, Matrix<double, 3, 1>{best_x, best_y, best_yaw});

		// compute gradient characteristic
		std::array<Matrix<double, 2, 1>, 2> grad_eigen_vectors;
		const Matrix<double, 2, 1> grad_eigen_values = compute_eigenvectors_2(grad_var_xyw.get<2, 2>(), grad_eigen_vectors);
		const Matrix<double, 3, 1> grad_std_uvw{sqrt(grad_eigen_values[0]), sqrt(grad_eigen_values[1]), sqrt(grad_var_xyw(2, 2))};

		// decide if we have 3D, 2D, 1D or 0D localization
		int mode = 0;
		if(best_score > m_min_score) {
			if(grad_std_uvw[0] > m_constrain_threshold) {
				if(grad_std_uvw[1] > m_constrain_threshold) {
					mode = 3;	// 2D position + rotation
				} else if(grad_std_uvw[2] > m_constrain_threshold_yaw) {
					mode = 2;	// 1D position + rotation
				} else {
					mode = 1;	// 1D position only
				}
			}
		}

		if(mode > 0)
		{
			double new_grid_x = best_x;
			double new_grid_y = best_y;
			double new_grid_yaw = best_yaw;

			if(mode < 3)
			{
				// constrain update to the good direction (ie. in direction of the eigen vector with the smaller sigma)
				const auto delta = Matrix<double, 2, 1>{best_x, best_y} - Matrix<double, 2, 1>{grid_pose[0], grid_pose[1]};
				const auto dist = grad_eigen_vectors[0].dot(delta);
				new_grid_x = grid_pose[0] + dist * grad_eigen_vectors[0][0];
				new_grid_y = grid_pose[1] + dist * grad_eigen_vectors[0][1];
			}
			if(mode < 2) {
				new_grid_yaw = grid_pose[2];	// keep old orientation
			}

			// use best sample for update
			Matrix<double, 4, 4> grid_pose_new = translate25(new_grid_x, new_grid_y) * rotate25_z(new_grid_yaw);

			// compute new odom to map offset from new grid pose
			const Matrix<double, 3, 1> new_offset =
					(m_grid_to_map * grid_pose_new * L.inverse() * Matrix<double, 4, 1>{0, 0, 0, 1}).project();

			// apply new offset with an exponential low pass filter
			m_offset_x += (new_offset[0] - m_offset_x) * m_update_gain;
			m_offset_y += (new_offset[1] - m_offset_y) * m_update_gain;
			m_offset_yaw += angles::shortest_angular_distance(m_offset_yaw, new_offset[2]) * m_update_gain;
		}
		m_offset_time = base_to_odom.stamp_;

		// update particle spread depending on mode
		if(mode >= 3) {
			m_sample_std_xy *= (1 - m_confidence_gain);
		} else {
			m_sample_std_xy += dist_moved * m_odometry_std_xy;
		}
		if(mode >= 2) {
			m_sample_std_yaw *= (1 - m_confidence_gain);
		} else {
			m_sample_std_yaw += rad_rotated * m_odometry_std_yaw;
		}

		// limit particle spread
		m_sample_std_xy = fmin(fmax(m_sample_std_xy, m_min_sample_std_xy), m_max_sample_std_xy);
		m_sample_std_yaw = fmin(fmax(m_sample_std_yaw, m_min_sample_std_yaw), m_max_sample_std_yaw);

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

		if(update_counter++ % 10 == 0) {
			ROS_INFO_STREAM("NeoLocalizationNode: score=" << float(best_score) << ", grad_uvw=[" << float(grad_std_uvw[0]) << ", " << float(grad_std_uvw[1])
					<< ", " << float(grad_std_uvw[2]) << "], std_xy=" << float(m_sample_std_xy) << " m, std_yaw=" << float(m_sample_std_yaw)
					<< " rad, mode=" << mode << "D, " << m_scan_buffer.size() << " scans");
		}

		// clear scan buffer
		m_scan_buffer.clear();
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
							<< " m, y=" <<  map_pose.getOrigin()[1] << " m, yaw=" << tf::getYaw(map_pose.getRotation()) << " rad");

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

			// reset particle spread to maximum
			m_sample_std_xy = m_max_sample_std_xy;
			m_sample_std_yaw = m_max_sample_std_yaw;

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

		// reset particle spread to maximum
		m_sample_std_xy = m_max_sample_std_xy;
		m_sample_std_yaw = m_max_sample_std_yaw;
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

		auto map = std::make_shared<GridMap<float>>(m_map_size, m_map_size, world_scale);

		// extract tile and convert to our format (occupancy between 0 and 1)
		for(int y = 0; y < map->size_y(); ++y) {
			for(int x = 0; x < map->size_x(); ++x) {
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
		const auto tile_center = (m_grid_to_map * Matrix<double, 4, 1>{	map->scale() * map->size_x() / 2,
																		map->scale() * map->size_y() / 2, 0, 1}).project();

		// publish new map tile for visualization
		auto ros_grid = boost::make_shared<nav_msgs::OccupancyGrid>();
		ros_grid->info.resolution = map->scale();
		ros_grid->info.width = map->size_x();
		ros_grid->info.height = map->size_y();
		ros_grid->info.origin.position.x = tile_origin[0];
		ros_grid->info.origin.position.y = tile_origin[1];
		tf::quaternionTFToMsg(tf::createQuaternionFromYaw(tile_origin[2]), ros_grid->info.origin.orientation);
		ros_grid->data.resize(map->num_cells());
		for(int y = 0; y < map->size_y(); ++y) {
			for(int x = 0; x < map->size_x(); ++x) {
				ros_grid->data[y * map->size_x() + x] = (*map)(x, y) * 100.f;
			}
		}
		m_pub_map_tile.publish(ros_grid);

		ROS_DEBUG_STREAM("NeoLocalizationNode: Got new grid at offset (" << tile_x << ", " << tile_y << ") [iworld], "
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

	ros::Timer m_loc_update_timer;

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
	double m_min_score = 0;
	double m_odometry_std_xy = 0;			// odometry xy error in meter per meter driven
	double m_odometry_std_yaw = 0;			// odometry yaw error in rad per rad rotated
	double m_min_sample_std_xy = 0;
	double m_min_sample_std_yaw = 0;
	double m_max_sample_std_xy = 0;
	double m_max_sample_std_yaw = 0;
	double m_constrain_threshold = 0;
	double m_constrain_threshold_yaw = 0;
	double m_loc_update_rate = 0;
	double m_map_update_rate = 0;
	double m_transform_timeout = 0;

	ros::Time m_offset_time;
	double m_offset_x = 0;					// current x offset between odom and map
	double m_offset_y = 0;					// current y offset between odom and map
	double m_offset_yaw = 0;				// current yaw offset between odom and map
	double m_sample_std_xy = 0;				// current sample spread in xy
	double m_sample_std_yaw = 0;			// current sample spread in yaw

	Matrix<double, 3, 1> m_last_odom_pose;
	Matrix<double, 4, 4> m_grid_to_map;
	Matrix<double, 4, 4> m_world_to_map;
	std::shared_ptr<GridMap<float>> m_map;			// map tile
	nav_msgs::OccupancyGrid::ConstPtr m_world;		// whole map

	int64_t update_counter = 0;
	std::map<std::string, sensor_msgs::LaserScan::ConstPtr> m_scan_buffer;

	Solver m_solver;
	std::mt19937 m_generator;
	std::thread m_map_update_thread;

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
		ROS_ERROR_STREAM("NeoLocalizationNode: " << ex.what());
		return -1;
	}

	return 0;
}

