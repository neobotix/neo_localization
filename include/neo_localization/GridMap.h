/*
 * GridMap.h
 *
 *  Created on: Apr 7, 2020
 *      Author: mad
 */

#ifndef INCLUDE_GRIDMAP_H_
#define INCLUDE_GRIDMAP_H_

#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdexcept>
#include <memory>
#include <vector>


/*
 * Class for a rectangular grid map.
 */
template<typename T>
class GridMap {
public:
	/*
	 * @param size_x Size of map in pixels
	 * @param size_y Size of map in pixels
	 * @param scale Size of one pixel in meters
	 */
	GridMap(int size_x, int size_y, float scale)
		:	m_size_x(size_x),
			m_size_y(size_y),
			m_scale(scale),
			m_inv_scale(1 / scale)
	{
		m_map = new T[size_t(size_x) * size_y];
	}

	/*
	 * Deep copy constructor.
	 */
	GridMap(const GridMap& other)
		:	GridMap(other.m_size_x, other.m_size_y, other.m_scale)
	{
		*this = other;
	}

	~GridMap()
	{
		delete [] m_map;
		m_map = 0;
	}

	/*
	 * Deep assignment operator.
	 */
	GridMap& operator=(const GridMap& other)
	{
		if(m_size_x != other.m_size_x || m_size_y != other.m_size_y) {
			throw std::logic_error("grid size mismatch");
		}
		m_scale = other.m_scale;
		m_inv_scale = other.m_inv_scale;
		::memcpy(m_map, other.m_map, num_cells() * sizeof(T));
		return *this;
	}

	int size_x() const {
		return m_size_x;
	}

	int size_y() const {
		return m_size_y;
	}

	float scale() const {
		return m_scale;
	}

	float inv_scale() const {
		return m_inv_scale;
	}

	size_t num_cells() const {
		return size_t(m_size_x) * m_size_y;
	}

	/*
	 * Sets all pixels to given value.
	 */
	void clear(const T& value) const {
		for(size_t i = 0; i < num_cells(); ++i) {
			m_map[i] = value;
		}
	}

	/*
	 * Access a given cell by index.
	 */
	T& operator[](size_t index)
	{
		return m_map[index];
	}

	const T& operator[](size_t index) const
	{
		return m_map[index];
	}

	/*
	 * Access a given cell by coordinate.
	 */
	T& operator()(int x, int y)
	{
		return m_map[size_t(y) * m_size_x + x];
	}

	const T& operator()(int x, int y) const
	{
		return m_map[size_t(y) * m_size_x + x];
	}

	/*
	 * Converts a coordinate in meters to a grid coordinate in pixels.
	 */
	float world_to_grid(float meters) const {
		return meters * m_inv_scale - 0.5f;
	}

	/*
	 * Bilinear interpolation at given pixel position.
	 * A coordinate of (0, 0) gives the exact value of the first pixel.
	 */
	float bilinear_lookup(float x, float y) const
	{
		const float a = x - floorf(x);
		const float b = y - floorf(y);

		return bilinear_lookup_ex(x, y, a, b);
	}

	/*
	 * Same as bilinear_lookup() but with pre-computed offsets a and b.
	 */
	float bilinear_lookup_ex(int x, int y, float a, float b) const
	{
		const int x0 = std::min(std::max(x, 0), m_size_x - 1);
		const int x1 = std::min(x0 + 1, m_size_x - 1);
		const int y0 = std::min(std::max(y, 0), m_size_y - 1);
		const int y1 = std::min(y0 + 1, m_size_y - 1);

		return		(*this)(x0, y0) * ((1.f - a) * (1.f - b))
				+	(*this)(x1, y0) * (a * (1.f - b))
				+	(*this)(x0, y1) * ((1.f - a) * b)
				+	(*this)(x1, y1) * (a * b);
	}

	/*
	 * Bilinear summation at a given pixel position.
	 * A coordinate of (0, 0) will sum at the first pixel exclusively.
	 */
	void bilinear_summation(float x, float y, const T& value)
	{
		const float a = x - floorf(x);
		const float b = y - floorf(y);

		bilinear_summation_ex(x, y, a, b, value);
	}

	/*
	 * Same as bilinear_summation() but with pre-computed offsets a and b.
	 */
	void bilinear_summation_ex(int x, int y, float a, float b, const T& value)
	{
		const int x0 = std::min(std::max(x, 0), m_size_x - 1);
		const int x1 = std::min(x0 + 1, m_size_x - 1);
		const int y0 = std::min(std::max(y, 0), m_size_y - 1);
		const int y1 = std::min(y0 + 1, m_size_y - 1);

		(*this)(x0, y0) += value * ((1.f - a) * (1.f - b));
		(*this)(x1, y0) += value * (a * (1.f - b));
		(*this)(x0, y1) += value * ((1.f - a) * b);
		(*this)(x1, y1) += value * (a * b);
	}

	/*
	 * Computes gauss-filtered and bilinear-interpolated first-order x and y gradient
	 * at given pixel position.
	 */
	void calc_gradient(float x, float y, float& dx, float& dy) const
	{
		static const float coeff_33_dxy[3][3] = {
				{-0.139505, 0, 0.139505},
				{-0.220989, 0, 0.220989},
				{-0.139505, 0, 0.139505}
		};

		const int x0 = x;
		const int y0 = y;

		const float a = x - floorf(x);
		const float b = y - floorf(y);

		dx = 0;
		dy = 0;

		for(int j = -1; j <= 1; ++j) {
			for(int i = -1; i <= 1; ++i) {
				const float value = bilinear_lookup_ex(x0 + i, y0 + j, a, b);
				dx += coeff_33_dxy[j+1][i+1] * value;
				dy += coeff_33_dxy[i+1][j+1] * value;
			}
		}

		dx /= 2 * m_scale;
		dy /= 2 * m_scale;
	}

	/*
	 * Computes gauss-filtered and bilinear-interpolated second-order x and y gradient
	 * at given pixel position.
	 */
	void calc_gradient2(float x, float y, float& ddx, float& ddy) const
	{
		static const float coeff_33_ddxy[3][3] = {
				{0.069752, -0.139504, 0.069752},
				{0.110494, -0.220989, 0.110494},
				{0.069752, -0.139504, 0.069752}
		};

		const int x0 = x;
		const int y0 = y;

		const float a = x - floorf(x);
		const float b = y - floorf(y);

		ddx = 0;
		ddy = 0;

		for(int j = -1; j <= 1; ++j) {
			for(int i = -1; i <= 1; ++i) {
				const float value = bilinear_lookup_ex(x0 + i, y0 + j, a, b);
				ddx += coeff_33_ddxy[j+1][i+1] * value;
				ddy += coeff_33_ddxy[i+1][j+1] * value;
			}
		}

		ddx /= 2 * m_scale;
		ddy /= 2 * m_scale;
	}

	/*
	 * Applies one smoothing iteration using a 3x3 gaussian kernel with sigma 1.
	 */
	void smooth_33_1()
	{
		static const float coeff_33_1[3][3] = {
				{0.077847, 0.123317, 0.077847},
				{0.123317, 0.195346, 0.123317},
				{0.077847, 0.123317, 0.077847}
		};

		GridMap<T> tmp(m_size_x, m_size_y, m_scale);

		for(int y = 0; y < m_size_y; ++y) {
			for(int x = 0; x < m_size_x; ++x) {
				float sum = 0;
				for(int j = -1; j <= 1; ++j) {
					const int y_ = std::min(std::max(y + j, 0), m_size_y - 1);
					for(int i = -1; i <= 1; ++i) {
						const int x_ = std::min(std::max(x + i, 0), m_size_x - 1);
						sum += coeff_33_1[j+1][i+1] * (*this)(x_, y_);
					}
				}
				tmp(x, y) = sum;
			}
		}
		*this = tmp;
	}

	/*
	 * Applies one smoothing iteration using a 5x5 gaussian kernel with sigma 2.
	 */
	void smooth_55_2()
	{
		static const float coeff_55_2[5][5] = {
				{0.0232468, 0.033824, 0.0383276, 0.033824, 0.0232468},
				{0.033824, 0.0492136, 0.0557663, 0.0492136, 0.033824},
				{0.0383276, 0.0557663, 0.0631915, 0.0557663, 0.0383276},
				{0.033824, 0.0492136, 0.0557663, 0.0492136, 0.033824},
				{0.0232468, 0.033824, 0.0383276, 0.033824, 0.0232468}
		};

		GridMap<T> tmp(m_size_x, m_size_y, m_scale);

		for(int y = 0; y < m_size_y; ++y) {
			for(int x = 0; x < m_size_x; ++x) {
				float sum = 0;
				for(int j = -2; j <= 2; ++j) {
					const int y_ = std::min(std::max(y + j, 0), m_size_y - 1);
					for(int i = -2; i <= 2; ++i) {
						const int x_ = std::min(std::max(x + i, 0), m_size_x - 1);
						sum += coeff_55_2[j+2][i+2] * (*this)(x_, y_);
					}
				}
				tmp(x, y) = sum;
			}
		}
		*this = tmp;
	}

	/*
	 * Returns a 2x downscaled map using a 4x4 gaussian filter with sigma 1.
	 */
	std::shared_ptr<GridMap<T>> downscale()
	{
		static const float coeff_44_1[4][4] {
			{0.0180824, 0.049153, 0.049153, 0.0180824},
			{0.049153, 0.133612, 0.133612, 0.049153},
			{0.049153, 0.133612, 0.133612, 0.049153},
			{0.0180824, 0.049153, 0.049153, 0.0180824}
		};

		auto res = std::make_shared<GridMap<T>>(m_size_x / 2, m_size_y / 2, m_scale * 2);

		for(int y = 0; y < m_size_y / 2; ++y) {
			for(int x = 0; x < m_size_x / 2; ++x) {
				float sum = 0;
				for(int j = -1; j <= 2; ++j) {
					const int y_ = std::min(std::max(y * 2 + j, 0), m_size_y - 1);
					for(int i = -1; i <= 2; ++i) {
						const int x_ = std::min(std::max(x * 2 + i, 0), m_size_x - 1);
						sum += coeff_44_1[j+1][i+1] * (*this)(x_, y_);
					}
				}
				(*res)(x, y) = sum;
			}
		}
		return res;
	}

private:
	int m_size_x = 0;
	int m_size_y = 0;

	float m_scale = 0;
	float m_inv_scale = 0;

	T* m_map = 0;

};


template<typename T>
class MultiGridMap {
public:
	std::vector<std::shared_ptr<GridMap<T>>> layers;

	/*
	 * @param size_x_m Size of map in meters
	 * @param size_y_m Size of map in meters
	 * @param scale Size of one pixel in meters
	 */
	MultiGridMap(float size_x_m, float size_y_m, float scale, int num_layers)
	{
		layers.resize(num_layers);
		const int base_size_x = size_x_m / (scale * (1 << (num_layers - 1))) + 0.5f;
		const int base_size_y = size_y_m / (scale * (1 << (num_layers - 1))) + 0.5f;

		for(int i = 0; i < num_layers; ++i) {
			layers[i] = std::make_shared<GridMap<T>>(base_size_x * (1 << (num_layers - i - 1)),
													 base_size_y * (1 << (num_layers - i - 1)),
													 scale * (1 << i));
		}
	}


};



#endif /* INCLUDE_GRIDMAP_H_ */
