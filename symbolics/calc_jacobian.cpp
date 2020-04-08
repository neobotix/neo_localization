/*
 * calc_jacobian.cpp
 *
 *  Created on: Apr 8, 2020
 *      Author: mad
 */


#include <vector>
#include <array>
#include <string>
#include <iostream>

#include <neo_common/Matrix.h>

#include <ginac/ginac.h>

using namespace GiNaC;


Matrix<ex, 3, 3> rotate_z(ex radians) {
	return {cos(radians), -sin(radians), 0,
			sin(radians), cos(radians), 0,
			0, 0, 1};
}

Matrix<ex, 3, 3> translate(ex x, ex y) {
	return {1, 0, x,
			0, 1, y,
			0, 0, 1};
}


int main()
{
	symbol p_x("pos_x");
	symbol p_y("pos_y");
	symbol p_yaw("pos_yaw");

	symbol s_x("scan_x");
	symbol s_y("scan_y");

	Matrix<ex, 3, 3> T = translate(p_x, p_y) * rotate_z(p_yaw);

	Matrix<ex, 3, 1> P = T * Matrix<ex, 3, 1>{s_x, s_y, 1};

	std::cout << csrc << "dxdx = " << diff(P[0], p_x) << std::endl;
	std::cout << csrc << "dxdy = " << diff(P[0], p_y) << std::endl;
	std::cout << csrc << "dxdyaw = " << diff(P[0], p_yaw) << std::endl;
	std::cout << csrc << "dydx = " << diff(P[1], p_x) << std::endl;
	std::cout << csrc << "dydy = " << diff(P[1], p_y) << std::endl;
	std::cout << csrc << "dydyaw = " << diff(P[1], p_yaw) << std::endl;

}

