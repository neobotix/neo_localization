/*
 * Solver.h
 *
 *  Created on: Apr 7, 2020
 *      Author: mad
 */

#ifndef INCLUDE_SOLVER_H_
#define INCLUDE_SOLVER_H_

#include <neo_common/Matrix.h>
#include <neo_localization/Util.h>
#include <neo_localization/GridMap.h>

#include <vector>


struct scan_point_t
{
	double x = 0;		// [m]
	double y = 0;		// [m]
	double w = 1;		// weight [1]
};


class Solver {
public:
	double pose_x = 0;					// initial guess / output grid pose
	double pose_y = 0;					// initial guess / output grid pose
	double pose_yaw = 0;				// initial guess / output grid pose

	double gain = 0.1;					// how fast to converge (0 to 1)
	double damping = 1;					// numerical hessian damping
	double r_norm = 0;					// current error norm

	template<typename T>
	void solve(	const GridMap<T>& map,
				const std::vector<scan_point_t>& points)
	{
		Matrix<double, 3, 1> G;			// gradient vector
		Matrix<double, 3, 3> H;			// Hessian matrix

		// compute transformation matrix first
		const float inv_scale = map.inv_scale();
		const Matrix<double, 3, 3> P = translate2(pose_x, pose_y) * rotate2_z(pose_yaw);

		r_norm = 0;

		for(const auto& point : points)
		{
			// transform sensor point to grid coordinates
			const auto q = (P * Matrix<double, 3, 1>{point.x, point.y, 1}).project();
			const float grid_x = q[0] * inv_scale - 0.5f;
			const float grid_y = q[1] * inv_scale - 0.5f;

			// compute error based on grid
			const double r_i = map.bilinear_lookup(grid_x, grid_y);
			r_norm += r_i * r_i;

			// compute error gradient based on grid
			float dx, dy;
			map.calc_gradient(grid_x, grid_y, dx, dy);

			// compute Jacobian row for this point / equation
			const double J_x = dx * 1.f;
			const double J_y = dy * 1.f;
			const double J_yaw =  dx * (-sin(pose_yaw) * point.x - cos(pose_yaw) * point.y)
								+ dy * ( cos(pose_yaw) * point.x - sin(pose_yaw) * point.y);

			// direct gradient vector summation
			G[0] += J_x * r_i;
			G[1] += J_y * r_i;
			G[2] += J_yaw * r_i;

			// direct Hessian matrix summation
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

		// we want average r_norm
		r_norm = sqrt(r_norm / points.size());

		// add Hessian damping
		H(0, 0) += damping;
		H(1, 1) += damping;
		H(2, 2) += damping;

		// solve Gauss-Newton step
		const auto X = H.inverse() * G;

		// apply new solution with a gain (optimize max. r_norm)
		pose_x += gain * X[0];
		pose_y += gain * X[1];
		pose_yaw += gain * X[2];
	}

};



#endif /* INCLUDE_SOLVER_H_ */
