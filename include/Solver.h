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


struct scan_point_t
{
	float x = 0;		// [m]
	float y = 0;		// [m]
	float w = 1;		// weight [1]
};


class Solver {
public:
	double pos_x = 0;
	double pos_y = 0;
	double pos_yaw = 0;

	double gain = 1;					// how fast to converge (0 to 1)
	double damping = 0;					// numerical hessian damping
	double r_norm = 0;					// current error norm

	template<typename T>
	void solve(	const GridMap<T>& map,
				const std::vector<scan_point_t>& points,
				const T target)
	{
		Matrix<double, 3, 1> G;			// gradient vector
		Matrix<double, 3, 3> H;			// hessian matrix

		for(const auto& point : points)
		{
			// TODO
		}

		H(0, 0) += damping;
		H(1, 1) += damping;
		H(2, 2) += damping;

		auto X = H.inverse() * G;
		pos_x -= gain * X[0];
		pos_y -= gain * X[1];
		pos_yaw -= gain * X[2];
	}

};



#endif /* INCLUDE_SOLVER_H_ */
