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
	float x = 0;		// [m]
	float y = 0;		// [m]
};

struct scan_point_ex_t : public scan_point_t
{
	float w = 1;		// weight [1]
	int layer = 0;		// layer index
};


class Solver {
public:
	double pose_x = 0;					// initial guess / output grid pose
	double pose_y = 0;					// initial guess / output grid pose
	double pose_yaw = 0;				// initial guess / output grid pose

	double gain = 0.1;					// how fast to converge (0 to 1)
	double damping = 1;					// numerical hessian damping
	double r_norm = 0;					// current error norm

	Matrix<double, 3, 1> G;				// gradient vector
	Matrix<double, 3, 3> H;				// Hessian matrix

	template<typename T>
	void solve(	const GridMap<T>& grid,
				const std::vector<scan_point_t>& points)
	{
		reset();

		// compute transformation matrix first
		const Matrix<double, 3, 3> P = transform2(pose_x, pose_y, pose_yaw);

		for(const auto& point : points)
		{
			// transform sensor point to grid coordinates
			const auto q = (P * Matrix<double, 3, 1>{point.x, point.y, 1}).project();
			const float grid_x = grid.world_to_grid(q[0]);
			const float grid_y = grid.world_to_grid(q[1]);

			// compute error based on grid
			const float r_i = grid.bilinear_lookup(grid_x, grid_y);
			r_norm += r_i * r_i;

			// compute error gradient based on grid
			float dx, dy;
			grid.calc_gradient(grid_x, grid_y, dx, dy);

			integrate(point.x, point.y, r_i, dx, dy);
		}

		// we want average r_norm
		r_norm = sqrt(r_norm / points.size());

		update();
	}

	template<typename T>
	void solve(	const MultiGridMap<T>& multi_grid,
				const std::vector<scan_point_ex_t>& points)
	{
		reset();

		// compute transformation matrix first
		const Matrix<double, 3, 3> P = transform2(pose_x, pose_y, pose_yaw);

		for(const auto& point : points)
		{
			auto& grid = multi_grid.layers[point.layer];

			// transform sensor point to grid coordinates
			const auto q = (P * Matrix<double, 3, 1>{point.x, point.y, 1}).project();
			const float grid_x = grid.world_to_grid(q[0]);
			const float grid_y = grid.world_to_grid(q[1]);

			// compute error based on grid
			const float r_i = grid.bilinear_lookup(grid_x, grid_y);
			r_norm += r_i * r_i;

			// compute error gradient based on grid
			float dx, dy;
			grid.calc_gradient(grid_x, grid_y, dx, dy);

			integrate(point.x, point.y, r_i, dx, dy);
		}

		// we want average r_norm
		r_norm = sqrt(r_norm / points.size());

		update();
	}

protected:
	void reset()
	{
		G = Matrix<double, 3, 1>();
		H = Matrix<double, 3, 3>();
		r_norm = 0;
	}

	void integrate(const float p_x, const float p_y, const float r_i, const float dx, const float dy)
	{
		const float J_x = dx * 1.f;
		const float J_y = dy * 1.f;
		const float J_yaw =   dx * (-sinf(pose_yaw) * p_x - cosf(pose_yaw) * p_y)
							+ dy * ( cosf(pose_yaw) * p_x - sinf(pose_yaw) * p_y);

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

	void update()
	{
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


/*
 * Computes a "virtual" covariance matrix based on second order gradients.
 * A higher covariance means a larger gradient, so the meaning of "covariance" is inverted here.
 * A higher gradient is better for localization accuracy.
 */
inline
Matrix<double, 3, 3> compute_virtual_scan_covariance_xyw(	std::shared_ptr<const GridMap<float>> grid,
															const std::vector<scan_point_t>& points,
															const Matrix<double, 3, 1>& pose)
{
	Matrix<double, 3, 3> var_xyw;
	const Matrix<double, 3, 3> P = transform2(pose);	// pre-compute transformation matrix

	for(const auto& point : points)
	{
		// transform sensor point to grid coordinates
		const auto q = (P * Matrix<double, 3, 1>{point.x, point.y, 1}).project();
		const float grid_x = grid->world_to_grid(q[0]);
		const float grid_y = grid->world_to_grid(q[1]);

		float ddx, ddy;
		grid->calc_gradient2(grid_x, grid_y, ddx, ddy);

		const float dir = pose[2] + atan2f(point.y, point.x);
		const float length = hypotf(point.x, point.y);
		const float ddyaw = (sinf(dir) * ddx + cosf(dir) * ddy) * length;

		var_xyw(0, 0) += ddx * ddx;
		var_xyw(1, 0) += ddx * ddy;
		var_xyw(0, 1) += ddy * ddx;
		var_xyw(1, 1) += ddy * ddy;
		var_xyw(2, 2) += ddyaw * ddyaw;
	}
	var_xyw *= 1. / points.size();
	return var_xyw;
}



#endif /* INCLUDE_SOLVER_H_ */
