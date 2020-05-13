/*
 * Util.h
 *
 *  Created on: Apr 27, 2020
 *      Author: mad
 */

#ifndef INCLUDE_NEO_LOCALIZATION_UTIL_H_
#define INCLUDE_NEO_LOCALIZATION_UTIL_H_

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

/*
 * Creates a 2D (x, y, yaw) transformation matrix.
 */
template<typename T>
Matrix<T, 3, 3> transform2(T x, T y, T rad) {
	return translate2(x, y) * rotate2_z(rad);
}

/*
 * Creates a 2D (x, y, yaw) transformation matrix.
 */
template<typename T>
Matrix<T, 3, 3> transform2(const Matrix<T, 3, 1>& pose) {
	return translate2(pose[0], pose[1]) * rotate2_z(pose[2]);
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
 * Creates a 2.5D (x, y, yaw) transformation matrix.
 */
template<typename T>
Matrix<T, 4, 4> transform25(T x, T y, T rad) {
	return translate25(x, y) * rotate25_z(rad);
}

/*
 * Creates a 2.5D (x, y, yaw) transformation matrix.
 */
template<typename T>
Matrix<T, 4, 4> transform25(const Matrix<T, 3, 1>& pose) {
	return translate25(pose[0], pose[1]) * rotate25_z(pose[2]);
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
// See: http://people.math.harvard.edu/~knill/teaching/math21b2004/exhibits/2dmatrices/index.html
// See: http://math.colgate.edu/~wweckesser/math312Spring06/handouts/IMM_2x2linalg.pdf
// Returns eigenvalues in descending order (with matching eigenvector order)
template<typename T>
Matrix<T, 2, 1> compute_eigenvectors_2(	const Matrix<T, 2, 2>& mat,
										std::array<Matrix<T, 2, 1>, 2>& eigen_vectors)
{
	Matrix<T, 2, 1> eigen_values;
	const T tmp_0 = std::sqrt(std::pow(mat(0, 0) + mat(1, 1), T(2)) - T(4) * (mat(0, 0) * mat(1, 1) - mat(1, 0) * mat(0, 1)));
	eigen_values[0] = (mat(0, 0) + mat(1, 1) + tmp_0) / T(2);
	eigen_values[1] = (mat(0, 0) + mat(1, 1) - tmp_0) / T(2);

	if(std::abs(eigen_values[0] - eigen_values[1]) > 1e-6)
	{
		for(int i = 0; i < 2; ++i)
		{
			const Matrix<T, 2, 1> vector_a {-1 * mat(0, 1), mat(0, 0) - eigen_values[i]};
			const Matrix<T, 2, 1> vector_b {mat(1, 1) - eigen_values[i], -1 * mat(1, 0)};

			if(vector_a.norm() > vector_b.norm()) {
				eigen_vectors[i] = vector_a;
			} else{
				eigen_vectors[i] = vector_b;
			}
			eigen_vectors[i].normalize();
		}
		if(eigen_values[1] > eigen_values[0]) {
			std::swap(eigen_values[0], eigen_values[1]);
			std::swap(eigen_vectors[0], eigen_vectors[1]);
		}
	} else {
		eigen_vectors[0] = Matrix<T, 2, 1>{1, 0};
		eigen_vectors[1] = Matrix<T, 2, 1>{0, 1};
	}
	return eigen_values;
}


#endif /* INCLUDE_NEO_LOCALIZATION_UTIL_H_ */
