/*
 * Util.h
 *
 *  Created on: Apr 27, 2020
 *      Author: mad
 */

#ifndef INCLUDE_NEO_LOCALIZATION_UTIL_H_
#define INCLUDE_NEO_LOCALIZATION_UTIL_H_

#include <neo_common/Matrix.h>
#include <neo_localization/GridMap.h>

#include <ros/ros.h>
#include <tf/transform_datatypes.h>
#include <nav_msgs/OccupancyGrid.h>


template<typename T>
Matrix<T, 3, 3> rotate2_z(T rad) {
	return {std::cos(rad), -std::sin(rad), 0,
			std::sin(rad), std::cos(rad), 0,
			0, 0, 1};
}

template<typename T>
Matrix<T, 3, 3> translate2(T x, T y) {
	return {1, 0, x,
			0, 1, y,
			0, 0, 1};
}

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

inline
nav_msgs::OccupancyGrid::Ptr convert_to_ros(std::shared_ptr<GridMap<float>> map, Matrix<double, 3, 1> origin, ros::Time stamp = ros::Time())
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



#endif /* INCLUDE_NEO_LOCALIZATION_UTIL_H_ */
