/*
 * Convert.h
 *
 *  Created on: Apr 30, 2020
 *      Author: mad
 */

#ifndef INCLUDE_NEO_LOCALIZATION_CONVERT_H_
#define INCLUDE_NEO_LOCALIZATION_CONVERT_H_

#include <neo_localization/GridMap.h>

#include <ros/ros.h>
#include <tf/transform_datatypes.h>
#include <nav_msgs/OccupancyGrid.h>


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
 * Converts a grid map to a ROS occupancy map.
 */
inline
nav_msgs::OccupancyGrid::Ptr convert_to_ros(	std::shared_ptr<GridMap<float>> map,
												Matrix<double, 3, 1> origin,
												ros::Time stamp = ros::Time())
{
	auto grid = boost::make_shared<nav_msgs::OccupancyGrid>();
	grid->header.stamp = stamp;
	grid->info.resolution = map->scale();
	grid->info.width = map->size_x();
	grid->info.height = map->size_y();
	grid->info.origin.position.x = origin[0];
	grid->info.origin.position.y = origin[1];
	tf::quaternionTFToMsg(tf::createQuaternionFromYaw(origin[2]), grid->info.origin.orientation);
	grid->data.resize(map->num_cells());
	for(int y = 0; y < map->size_y(); ++y) {
		for(int x = 0; x < map->size_x(); ++x) {
			grid->data[y * map->size_x() + x] = (*map)(x, y) * 100.f;
		}
	}
	return grid;
}

/*
 * Converts a grid map to a binary (-1, 0 or 100) ROS occupancy map.
 */
inline
nav_msgs::OccupancyGrid::Ptr convert_to_ros_binary(	std::shared_ptr<GridMap<float>> map,
													Matrix<double, 3, 1> origin,
													float threshold,
													ros::Time stamp = ros::Time())
{
	static const float coeff_55[5][5] = {
			{0.00118231, 0.01357, 0.0276652, 0.01357, 0.00118231},
			{0.01357, 0.06814, -0, 0.06814, 0.01357},
			{0.0276652, -0, -0.50349, -0, 0.0276652},
			{0.01357, 0.06814, -0, 0.06814, 0.01357},
			{0.00118231, 0.01357, 0.0276652, 0.01357, 0.00118231},
	};

	auto grid = boost::make_shared<nav_msgs::OccupancyGrid>();
	grid->header.stamp = stamp;
	grid->info.resolution = map->scale();
	grid->info.width = map->size_x();
	grid->info.height = map->size_y();
	grid->info.origin.position.x = origin[0];
	grid->info.origin.position.y = origin[1];
	tf::quaternionTFToMsg(tf::createQuaternionFromYaw(origin[2]), grid->info.origin.orientation);
	grid->data.resize(map->num_cells());
	for(int y = 0; y < map->size_y(); ++y) {
		for(int x = 0; x < map->size_x(); ++x) {
			// compute LoG filter
			float sum = 0;
			for(int j = -2; j <= 2; ++j) {
				const int y_ = std::min(std::max(y + j, 0), map->size_y() - 1);
				for(int i = -2; i <= 2; ++i) {
					const int x_ = std::min(std::max(x + i, 0), map->size_x() - 1);
					sum += coeff_55[j+2][i+2] * (*map)(x_, y_);
				}
			}
			if(-1 * sum > threshold) {
				grid->data[y * map->size_x() + x] = 100;
			} else if((*map)(x, y) < 0) {
				grid->data[y * map->size_x() + x] = -1;
			} else {
				grid->data[y * map->size_x() + x] = 0;
			}
		}
	}
	return grid;
}


#endif /* INCLUDE_NEO_LOCALIZATION_CONVERT_H_ */
