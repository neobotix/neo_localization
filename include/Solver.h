/*
 * Solver.h
 *
 *  Created on: Apr 7, 2020
 *      Author: mad
 */

#ifndef INCLUDE_SOLVER_H_
#define INCLUDE_SOLVER_H_

#include "GridMap.h"

#include <neo_common/Matrix.h>

#include <vector>


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

struct scan_point_t
{
	float x = 0;		// [m]
	float y = 0;		// [m]
	float w = 1;		// weight [1]
};


class Solver {
public:
	double pose_x = 0;
	double pose_y = 0;
	double pose_yaw = 0;

	double gain = 0.1;					// how fast to converge (0 to 1)
	double damping = 1;					// numerical hessian damping
	double r_norm = 0;					// current error norm

	template<typename T>
	void solve(	const GridMap<T>& map,
				const std::vector<scan_point_t>& points,
				const float target)
	{
		Matrix<double, 3, 1> G;			// gradient vector
		Matrix<double, 3, 3> H;			// hessian matrix

		const float inv_scale = map.inv_scale();
		const Matrix<double, 3, 3> P = translate2(pose_x, pose_y) * rotate2_z(pose_yaw);

		r_norm = 0;

		for(const auto& point : points)
		{
			const auto q = (P * Matrix<double, 3, 1>{point.x, point.y, 1}).project();
			const float grid_x = q[0] * inv_scale - 0.5f;
			const float grid_y = q[1] * inv_scale - 0.5f;

			const double r_i = map.bilinear_lookup(grid_x, grid_y) - target;
			r_norm += r_i * r_i;

			float dx, dy;
			map.calc_gradient(grid_x, grid_y, dx, dy);

			const double J_x = dx * 1.f;
			const double J_y = dy * 1.f;
			const double J_yaw =  dx * (-sin(pose_yaw) * point.x - cos(pose_yaw) * point.y)
								+ dy * ( cos(pose_yaw) * point.x - sin(pose_yaw) * point.y);

			G[0] += J_x * r_i;
			G[1] += J_y * r_i;
			G[2] += J_yaw * r_i;

			H(0, 0) += J_x * J_x;
			H(1, 1) += J_y * J_y;
			H(2, 2) += J_yaw * J_yaw;

			H(0, 1) += J_x * J_y;
			H(1, 0) += J_x * J_y;

			H(0, 2) += J_x * J_yaw;
			H(2, 0) += J_x * J_yaw;

			H(1, 2) += J_y * J_yaw;
			H(2, 1) += J_y * J_yaw;
		}

		r_norm /= points.size();

		H(0, 0) += damping;
		H(1, 1) += damping;
		H(2, 2) += damping;

		const auto X = H.inverse() * G;
		pose_x -= gain * X[0];
		pose_y -= gain * X[1];
		pose_yaw -= gain * X[2];
	}

};



#endif /* INCLUDE_SOLVER_H_ */
