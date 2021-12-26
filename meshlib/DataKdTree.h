/////////////////////////////////////////////////////////////////////////////
// DKdTree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2016-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DATAKDTREE_H__INCLUDED
#define DATAKDTREE_H__INCLUDED

#include "common.h"
#include <iostream>
#include <functional>
#include <vector>

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"


template<typename T, typename KdTreeElement, size_t SamplesSize>
struct Kd3dRegularErrorEstimator
{
	static double estimateError(const DBox &box, const KdTreeElement * kd_te,
		const std::function<T(const DPoint3d & pt)> &f)
	{
		double error = 0.0;
		double fx = 1.0 / (double)SamplesSize;
		double dx = box.getDX() * fx;
		double dy = box.getDY() * fx;
		double dz = box.getDZ() * fx;

		DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
		for (; pt.x < box.x1; pt.x += dx) {
			for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
				for (pt.z = box.z0 + 0.5*dz; pt.z < box.z1; pt.z += dz) {
					error = std::max(diffKdValue(f(pt), kd_te->getValue(pt, box)),
						error);
				}
			}
		}

		return error;
	}
};

template<typename T, typename KdTreeElement, size_t SamplesSize>
struct Kd2dRegularErrorEstimator
{
	static double estimateError(const DRect &box, const KdTreeElement * kd_te,
		const std::function<T(const DPoint2d & pt)> &f)
	{
		double error = 0.0;
		double fx = 1.0 / (double)SamplesSize;
		double dx = box.getDX() * fx;
		double dy = box.getDY() * fx;

		DPoint2d pt(box.x0 + 0.5*dx, box.y0);
		for (; pt.x < box.x1; pt.x += dx) {
			for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
				error = std::max(diffKdValue(f(pt), kd_te->getValue(pt, box)),
					error);
			}
		}

		return error;
	}
};

/**
* An example implementation of the KdSplitter.
*
* Splits the box in the middle along the longest axis.
*/
template<typename T>
struct Kd3dHalfLongestSplitter
{
	static void whereToSplit(const DBox & box, const std::function<T(const DPoint3d & pt)> &,
		Axis & split_axis, double & split_value)
	{
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
	static void whereToSplit(const DBox & box, Axis & split_axis, double & split_value)
	{
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
};

template<typename T>
struct Kd2dHalfLongestSplitter
{
	static void whereToSplit(const DRect & box, const std::function<T(const DPoint2d & pt)> &,
		Axis & split_axis, double & split_value)
	{
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
	static void whereToSplit(const DRect & box, Axis & split_axis, double & split_value)
	{
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
};

/**
* An example implementation of the KdSplitter.
*
* Splits the box with golden-ratio, along the longest axis (in metric), closest to the given point.
*/
template<typename T>
struct Kd3dGoldenRatioLongestSplitter
{
	static void whereToSplit(const DBox & box, const std::function<T(const DPoint3d & pt)> &,
		Axis & split_axis, double & split_value)
	{
		assert(false); // just for interface rule, but maybe can be written to work...
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
	static void whereToSplit(const DBox & box, const DPoint3d& pt, const DVector3d & dlen,
		Axis & split_axis, double & split_value)
	{
		// longest axis, but with metric (metric lengths stored in dlen)
		if (dlen.x > dlen.y && dlen.x > dlen.z) split_axis = Axis::X;
		else if (dlen.y > dlen.z) split_axis = Axis::Y;
		else split_axis = Axis::Z;
		//split_axis = box.getLongestAxis();
		const double GR1 = 2.0 / (1 + sqrt(5.0));
		const double GR0 = 1.0 - GR1;
		double dp = pt[split_axis];
		double d0 = box.getD0(split_axis);
		double d10 = box.getDLen(split_axis);
		double t = (dp - d0) / d10;
		t = (abs(t - GR0) < abs(t - GR1)) ? GR0 : GR1;
		split_value = d0 + t * d10;
	}
};

template<typename T>
struct Kd2dGoldenRatioLongestSplitter
{
	static void whereToSplit(const DRect & box, const std::function<T(const DPoint2d & pt)> &,
		Axis & split_axis, double & split_value)
	{
		assert(false); // just for interface rule, but maybe can be written to work...
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
	static void whereToSplit(const DRect & box, const DPoint2d& pt, const DVector2d & dlen,
		Axis & split_axis, double & split_value)
	{
		// longest axis, but with metric (metric lengths stored in dlen)
		split_axis = (dlen.x > dlen.y) ? Axis::X : Axis::Y;
		//split_axis = box.getLongestAxis();
		const double GR1 = 2.0 / (1 + sqrt(5.0));
		const double GR0 = 1.0 - GR1;
		double dp = pt[split_axis];
		double d0 = box.getD0(split_axis);
		double d10 = box.getDLen(split_axis);
		double t = (dp - d0) / d10;
		t = (abs(t - GR0) < abs(t - GR1)) ? GR0 : GR1;
		split_value = d0 + t * d10;
	}
};
/**
* An example implementation of the KdSplitter.
*
* Splits the box with golden-ratio, along the longest axis (in metric), closest to the given point.
*/
template<typename T>
struct Kd3dMinLongestSplitter
{
	static void whereToSplit(const DBox & box, const std::function<T(const DPoint3d & pt)> &,
		Axis & split_axis, double & split_value)
	{
		assert(false); // just for interface ...
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
	static void whereToSplit(const DBox & box, const DPoint3d& pt, const DVector3d & dlen2,
		Axis & split_axis, double & split_value)
	{
		// longest axis, but with metric (metric lengths stored in dlen)
		if (dlen2.x > dlen2.y && dlen2.x > dlen2.z) split_axis = Axis::X;
		else if (dlen2.y > dlen2.z) split_axis = Axis::Y;
		else split_axis = Axis::Z;
		double len = sqrt(dlen2[split_axis]);
		//split_axis = box.getLongestAxis();
		double dp = pt[split_axis];
		double d0 = box.getD0(split_axis);
		double d10 = box.getDLen(split_axis);
		double t = (dp - d0) / d10;
		double t_len = t * len;
		if (t_len <= 2.0)
			split_value = d0 + t*d10;
		else if (len - t_len <= 2.0)
			split_value = d0 + (1.0 - t)*d10;
		else if( t < 0.5)
			split_value = d0 + (t*(t_len-1.0)/t_len)*d10; // len==1.0 (metric) before _t_
		else
			split_value = d0 + (1.0 - (1.0-t)*(t_len - 1.0) / t_len)*d10; // len==1.0 (metric) after _t_
	}
};

template<typename T>
struct Kd2dMinLongestSplitter
{
	static void whereToSplit(const DRect & box, const std::function<T(const DPoint2d & pt)> &,
		Axis & split_axis, double & split_value)
	{
		assert(false); // just for interface ...
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
	static void whereToSplit(const DRect & box, const DPoint2d& pt, const DVector2d & dlen2,
		Axis & split_axis, double & split_value)
	{
		// longest axis, but with metric (metric lengths stored in dlen)
		split_axis = (dlen2.x > dlen2.y) ? Axis::X : Axis::Y;
		double len = sqrt(dlen2[split_axis]);
		//split_axis = box.getLongestAxis();
		double dp = pt[split_axis];
		double d0 = box.getD0(split_axis);
		double d10 = box.getDLen(split_axis);
		double t = (dp - d0) / d10;
		double t_len = t * len;
		if (t_len <= 2.0)
			split_value = d0 + t*d10;
		else if (len - t_len <= 2.0)
			split_value = d0 + (1.0 - t)*d10;
		else if (t < 0.5)
			split_value = d0 + (t*(t_len - 1.0) / t_len)*d10; // len==1.0 (metric) before _t_
		else
			split_value = d0 + (1.0 - (1.0 - t)*(t_len - 1.0) / t_len)*d10; // len==1.0 (metric) after _t_
	}
};
/**
* An example implementation of the BoxSplitter (see \ref kd_tree).
*
* Splits the box according to gradient of the represented function.
*
* The current box is split into SamplesSize segments along the axes.
* Let's introduce split boxes numbering, \f$ S_{i,j,k} \f$
* where \f$ i,j,k \in [1, \, \text{SamplesSize}] = [1,\,N] \f$.
* Then we will consider gradient to be the difference of function value
* between two faces of the box \f$ S_{i,j,k} \f$, i.e.
* \f{equation}{
*   \nabla^x_{i,j,k} = f(x_m+\frac{1}{2}d_x, y_m, z_m) - f(x_m-\frac{1}{2}d_x, y_m, z_m)
* \f}
* where \f$ d_x \f$ is the size of the box along the x axis
*   and \f$ x_m, y_m, z_m \f$ are the coordinates of the box's center.
*
* For each dimension (x, y and z) the sums of absolute values of
* gradient are computed as
* \f{equation}{
*   \Sigma_x = \sum\limits_{i=1}^N \sum\limits_{j=1}^N \sum\limits_{k=1}^N \left| \nabla^x_{i,j,k} \right|
*   \text{ .}
* \f}
* The axis with the largest sum of gradients is then chosen to be split along.
*
* In order to determine the point of split, a value of T is obtained
* \f{equation}{
*   T = \min \{ t \mid
*         \sum\limits_{i=1}^t \sum\limits_{j=1}^N \sum\limits_{k=1}^N
*           \left| \nabla^x_{i,j,k} \right| > \Sigma_x \}
* \f}
* The box is split at the position \f$ x_m \f$ of the box \f$ S_{T,j,k} \f$
* where \f$ j, k \in [1,\,N] \f$.
*
*/
template<typename T, size_t SamplesSize>
struct KdMaxGradientSumSplitter
{
	static_assert(SamplesSize > 0, "number of samples must be nonzero");

	static void whereToSplit(const DBox &box, const std::function<T(const DPoint3d & pt)> &f,
		Axis & split_axis, double & split_value)
	{
		std::vector<T> cache;
		cache.reserve(SamplesSize*SamplesSize*SamplesSize);

		double fx = 1.0 / (double)SamplesSize;
		double dx = box.getDX() * fx;
		double dy = box.getDY() * fx;
		double dz = box.getDZ() * fx;

		DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
		for (; pt.x < box.x1; pt.x += dx) {
			for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
				for (pt.z = box.z0 + 0.5*dz; pt.z < box.z1; pt.z += dz) {
					cache.push_back(f(pt));
				}
			}
		}

		size_t step_i[] = { SamplesSize * SamplesSize,	SamplesSize, 1 };
		auto ind = [&step_i](size_t i, size_t j, size_t k) {
			return step_i[0] * i + step_i[1] * j + step_i[2] * k;
		};

		double sum[] = { 0, 0, 0 };
		for (size_t i = 0; i < SamplesSize; i++) {
			for (size_t j = 0; j < SamplesSize; j++) {
				for (size_t k = 0; k < SamplesSize; k++) {
					if (i > 0) sum[0] += diffKdValue(cache[ind(i, j, k)], cache[ind(i - 1, j, k)]);
					if (j > 0) sum[1] += diffKdValue(cache[ind(i, j, k)], cache[ind(i, j - 1, k)]);
					if (k > 0) sum[2] += diffKdValue(cache[ind(i, j, k)], cache[ind(i, j, k - 1)]);
				}
			}
		}

		double d0[] = { box.x0, box.y0, box.z0 };
		double step[] = { dx, dy, dz };
		if (sum[0] > sum[1]  && sum[0] > sum[2]) split_axis = Axis::X;
		else split_axis = (sum[1] > sum[2]) ? Axis::Y : Axis::Z;
		int dim = static_cast<int>(split_axis);

		auto ind_ijk = [&ind, dim](size_t i, size_t j, size_t k) {
			switch (dim) {
			case 0:	return ind(i, j, k);
			case 1:	return ind(k, i, j);
			case 2:	return ind(j, k, i);
			}
			assert(!"unreachable");
			return ind(i, j, k);
		};

		double sum_half = 0;
		split_value = d0[dim] + step[dim] * SamplesSize / 2;
		for (size_t i = 1; i < SamplesSize; i++) {
			for (size_t j = 0; j < SamplesSize; j++) {
				for (size_t k = 0; k < SamplesSize; k++) {
					sum_half += diffKdValue(cache[ind_ijk(i, j, k)], cache[ind_ijk(i - 1, j, k)]);
				}
			}
			if (2 * sum_half > sum[dim]) {
				split_value = d0[dim] + step[dim] * i;
				break;
			}
		}
	}
};

template<typename T, size_t SamplesSize>
struct KdMaxGradientSplitter
{
	static_assert(SamplesSize > 0, "number of samples must be nonzero");

	static void whereToSplit(const DBox &box, const std::function<T(const DPoint3d & pt)> &f,
		Axis & split_axis, double & split_value)
	{
		T cache[SamplesSize][SamplesSize][SamplesSize];

		double fx = 1.0 / (double)SamplesSize;
		double dx = box.getDX() * fx;
		double dy = box.getDY() * fx;
		double dz = box.getDZ() * fx;

		DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
		for (size_t i = 0; i < SamplesSize; i++, pt.x += dx) {
			pt.y = box.y0 + 0.5*dy;
			for (size_t j = 0; j < SamplesSize; j++, pt.y += dy) {
				pt.z = box.z0 + 0.5*dz;
				for (size_t k = 0; k < SamplesSize; k++, pt.z += dz) {
					cache[i][j][k] = f(pt);
				}
			}
		}

		double max_diff = 0.0;
		for (size_t i = 0; i < SamplesSize; i++) {
			for (size_t j = 0; j < SamplesSize; j++) {
				for (size_t k = 0; k < SamplesSize; k++) {
					if (i > 0) {
						double diff = diffKdValue(cache[i][j][k], cache[i - 1][j][k]);
						if (diff > max_diff) {
							max_diff = diff;
							split_axis = Axis::X;
							split_value = box.x0 + i * dx;
						}
					}
					if (j > 0) {
						double diff = diffKdValue(cache[i][j][k], cache[i][j - 1][k]);
						if (diff > max_diff) {
							max_diff = diff;
							split_axis = Axis::Y;
							split_value = box.y0 + j * dy;
						}
					}
					if (k > 0) {
						double diff = diffKdValue(cache[i][j][k], cache[i][j][k - 1]);
						if (diff > max_diff) {
							max_diff = diff;
							split_axis = Axis::Z;
							split_value = box.z0 + k * dz;
						}
					}
				}
			}
		}
	}
};

template<typename T, size_t SamplesSize>
struct Kd2dMaxGradientSplitter
{
	static_assert(SamplesSize > 0, "number of samples must be nonzero");

	static void whereToSplit(const DRect &box, const std::function<T(const DPoint2d & pt)> &f,
		Axis & split_axis, double & split_value)
	{
		T cache[SamplesSize][SamplesSize];

		double fx = 1.0 / (double)SamplesSize;
		double dx = box.getDX() * fx;
		double dy = box.getDY() * fx;

		DPoint2d pt(box.x0 + 0.5*dx, box.y0);
		for (size_t i = 0; i < SamplesSize; i++, pt.x += dx) {
			pt.y = box.y0 + 0.5*dy;
			for (size_t j = 0; j < SamplesSize; j++, pt.y += dy) {
				cache[i][j] = f(pt);
			}
		}

		double max_diff = 0.0;
		for (size_t i = 0; i < SamplesSize; i++) {
			for (size_t j = 0; j < SamplesSize; j++) {
				if (i > 0) {
					double diff = diffKdValue(cache[i][j], cache[i - 1][j]);
					if (diff > max_diff) {
						max_diff = diff;
						split_axis = Axis::X;
						split_value = box.x0 + i * dx;
					}
				}
				if (j > 0) {
					double diff = diffKdValue(cache[i][j], cache[i][j - 1]);
					if (diff > max_diff) {
						max_diff = diff;
						split_axis = Axis::Y;
						split_value = box.y0 + j * dy;
					}
				}
			}
		}
	}
};

template<typename T, size_t SamplesSize>
struct KdMaxDiffGradientSplitter
{
	static_assert(SamplesSize > 1, "number of samples must be greater than 1");

	static void whereToSplit(const DBox &box, const std::function<T(const DPoint3d & pt)> &f,
		Axis & split_axis, double & split_value)
	{
		T cache[SamplesSize][SamplesSize][SamplesSize];

		double fx = 1.0 / (double)SamplesSize;
		double dx = box.getDX() * fx;
		double dy = box.getDY() * fx;
		double dz = box.getDZ() * fx;

		DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
		for (size_t i = 0; i < SamplesSize; i++, pt.x += dx) {
			pt.y = box.y0 + 0.5*dy;
			for (size_t j = 0; j < SamplesSize; j++, pt.y += dy) {
				pt.z = box.z0 + 0.5*dz;
				for (size_t k = 0; k < SamplesSize; k++, pt.z += dz) {
					cache[i][j][k] = f(pt);
				}
			}
		}

		double diff, max_avg_diff = 0.0;
		for (size_t i = 0; i < SamplesSize; i++) {
			for (size_t j = 0; j < SamplesSize; j++) {
				for (size_t k = 0; k < SamplesSize; k++) {
					if (i > 0) {
						diff = diffKdValue(cache[i][j][k], cache[i - 1][j][k]);
						if (diff > max_avg_diff) {
							max_avg_diff = diff;
							split_axis = Axis::X;
							split_value = box.x0 + (i - 0.5) * dx;
						}
					}
					if (j > 0) {
						diff = diffKdValue(cache[i][j][k], cache[i][j - 1][k]);
						if (diff > max_avg_diff) {
							max_avg_diff = diff;
							split_axis = Axis::Y;
							split_value = box.y0 + (j - 0.5) * dy;
						}
					}
					if (k > 0) {
						diff = diffKdValue(cache[i][j][k], cache[i][j][k - 1]);
						if (diff > max_avg_diff) {
							max_avg_diff = diff;
							split_axis = Axis::Z;
							split_value = box.z0 + (k - 0.5) * dz;
						}
					}
				}
			}
		}
	}
};

/**
 * This class implements a kd-tree in a 3-dimensional space
 *  (values in leaves only) plus some basic operations.
 */
template<typename T,
		typename KdSplitter = Kd3dHalfLongestSplitter<T>,
		size_t SamplingSize = 4,
		size_t MaxLevel = 20>
class DataKdTreeL3d
{
public:
	DataKdTreeL3d(const DBox & bbox) : m_box(bbox) {}
	~DataKdTreeL3d() {	clear(); }
private:
	DataKdTreeL3d(const DataKdTreeL3d& tree) = delete;
	DataKdTreeL3d& operator=(const DataKdTreeL3d& tree) = delete;
public:
	int adaptToField(const std::function<T(const DPoint3d & pt)> & f, double max_error) {
		T v = f(m_box.getMiddlePoint());
		if (m_top != nullptr) clear();
		m_top = new KdElement(v);
		return 1 + m_top->adaptToField(f, m_box, max_error);
	}
	void clear() {
		if (m_top != nullptr) {
			m_top->clear();
			delete m_top;
			m_top = nullptr;
		}
	}
	const T& getValue(const DPoint3d& pt) const {
		KdElement* top = m_top;
		while (!top->is_leaf) {
			top = (pt[top->m_data.split.axis] <= top->m_data.split.value) ? 
				top->m_data.split.element_0 : top->m_data.split.element_1;
		}
		return top->m_data.value;
	}
	inline const DBox& getBBox() const { return m_box; }
public: // stat
	void forRegularGrid(int grid_size, const std::function<void(const DPoint3d& pt)>& fg) {
		double fx = 1.0 / (double)grid_size;
		double dx = m_box.getDX() * fx;
		double dy = m_box.getDY() * fx;
		double dz = m_box.getDZ() * fx;

		DPoint3d pt(m_box.x0 + 0.5*dx, m_box.y0, m_box.z0);
		for (; pt.x < m_box.x1; pt.x += dx) {
			for (pt.y = m_box.y0 + 0.5*dy; pt.y < m_box.y1; pt.y += dy) {
				for (pt.z = m_box.z0 + 0.5*dz; pt.z < m_box.z1; pt.z += dz) {
					fg(pt);
				}
			}
		}
	}
	int getElementCount() const { return m_top ? m_top->getElementCount() : 0; }
	int getValueCount() const { return m_top ? m_top->getValueCount() : 0; }
	int getMaxDepth() const { return m_top ? m_top->getMaxDepth() : 0; }
	int getMaxBalance() const { 
		int max_balance = 0;
		if(m_top) m_top->getMaxDepthAndBalance(max_balance);
		return max_balance;
	}
	int getTotalBytes() const {
		return (int)sizeof(*this) // main struct
			+ (int)sizeof(KdElement) * getElementCount(); // tree nodes
	}
private:
	class KdElement {
	public:
		KdElement(const T & v) : m_data(v) { }
		int adaptToField(const std::function<T(const DPoint3d & pt)> & f, 
			const DBox & box, double max_error, int level = 0) 
		{

			//static size_t counter = 0;
			//if (++counter % 100 == 0) {
			//	LOG4CPLUS_INFO(MeshLog::logger_console, "** counter", counter);
			//	LOG4CPLUS_INFO(MeshLog::logger_console, "** box-diam", box.getDiameter());
			//}

			assert(is_leaf);
			if (level >= MaxLevel) return 0;
			double curr_error = estimateError(box, f);
			if (curr_error <= max_error) return 0;
			is_leaf = false;
			KdSplitter::whereToSplit(box, f, m_data.split.axis, m_data.split.value);
			DBox box0, box1;
			box.split(m_data.split.axis, m_data.split.value, box0, box1);
			T v0 = f(box0.getMiddlePoint());
			T v1 = f(box1.getMiddlePoint());
			m_data.split.element_0 = new KdElement(v0);
			m_data.split.element_1 = new KdElement(v1);
			return 2 + 
				m_data.split.element_0->adaptToField(f, box0, max_error, level+1) +
				m_data.split.element_1->adaptToField(f, box1, max_error, level+1);
		}
		double estimateError(const DBox &box, const std::function<T(const DPoint3d & pt)> &f)
		{
			double error = 0.0;
			double fx = 1.0 / (double)SamplingSize;
			double dx = box.getDX() * fx;
			double dy = box.getDY() * fx;
			double dz = box.getDZ() * fx;

			DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
			for (; pt.x < box.x1; pt.x += dx) {
				for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
					for (pt.z = box.z0 + 0.5*dz; pt.z < box.z1; pt.z += dz) {
						error = std::max(diffKdValue(f(pt), m_data.value),
							error);
					}
				}
			}

			return error;
		}
		void clear() {
			if (!is_leaf) {
				m_data.split.element_0->clear();
				delete m_data.split.element_0;
				m_data.split.element_1->clear();
				delete m_data.split.element_1;
				is_leaf = true;
			}
		}
	public:
		int getElementCount() const {
			return is_leaf ? 1 :
				1 + m_data.split.element_0->getElementCount() + m_data.split.element_1->getElementCount();
		}
		int getValueCount() const {
			return is_leaf ? 1 :
				m_data.split.element_0->getValueCount() + m_data.split.element_1->getValueCount();
		}
		int getMaxDepth() const {
			return is_leaf ? 1 :
				1+std::max(m_data.split.element_0->getMaxDepth(), m_data.split.element_1->getMaxDepth());
		}
		int getMaxDepthAndBalance(int & max_balance) const {
			if (is_leaf) return 1;
			int left_depth  = m_data.split.element_0->getMaxDepthAndBalance(max_balance);
			int right_depth = m_data.split.element_1->getMaxDepthAndBalance(max_balance);
			int balance = std::abs(left_depth - right_depth);
			if (balance > max_balance) max_balance = balance;
			return 1+std::max(left_depth, right_depth);
		}
	public:
		bool is_leaf = true;
		union KdData {
			KdData(const T & v) : value(v) {}
			struct KdSplitElement{
				Axis		axis		= Axis::X;
				double		value		= 0.0;
				KdElement * element_0	= nullptr;
				KdElement * element_1	= nullptr;
			} split;
			T value;
		} m_data;
	};
private:
	DBox m_box;
	KdElement* m_top = nullptr;
};

/**
* This class implements a kd-tree in a 3-dimensional space
*  (values in vertices) plus some basic operations.
*/
template<typename T,
	typename KdSplitter = Kd3dHalfLongestSplitter<T>,
	size_t SamplingSize = 4,
	size_t MaxLevel = 20>
	class DataKdTreeV3d
{
public:
	DataKdTreeV3d(const DBox & bbox) : m_box(bbox), m_values(100) {}
	~DataKdTreeV3d() { clear(); }
private:
	DataKdTreeV3d(const DataKdTreeV3d& tree) = delete;
	DataKdTreeV3d& operator=(const DataKdTreeV3d& tree) = delete;
public:
	int adaptToField(const std::function<T(const DPoint3d & pt)> & f, double max_error) {
		if (m_top != nullptr) clear();
		typename KdElement::KdVertices kv;

		kv.vx_lsw = &m_values.get(m_values.add(f(DPoint3d(m_box.x0, m_box.y0, m_box.z0))));
		kv.vx_lse = &m_values.get(m_values.add(f(DPoint3d(m_box.x1, m_box.y0, m_box.z0))));
		kv.vx_lnw = &m_values.get(m_values.add(f(DPoint3d(m_box.x0, m_box.y1, m_box.z0))));
		kv.vx_lne = &m_values.get(m_values.add(f(DPoint3d(m_box.x1, m_box.y1, m_box.z0))));
		kv.vx_hsw = &m_values.get(m_values.add(f(DPoint3d(m_box.x0, m_box.y0, m_box.z1))));
		kv.vx_hse = &m_values.get(m_values.add(f(DPoint3d(m_box.x1, m_box.y0, m_box.z1))));
		kv.vx_hnw = &m_values.get(m_values.add(f(DPoint3d(m_box.x0, m_box.y1, m_box.z1))));
		kv.vx_hne = &m_values.get(m_values.add(f(DPoint3d(m_box.x1, m_box.y1, m_box.z1))));

		m_top = new KdElement(kv);
		return 1 + m_top->adaptToField(f, m_box, max_error, m_values);
	}
	void clear() {
		if (m_top != nullptr) {
			m_top->clear();
			delete m_top;
			m_top = nullptr;
			m_values.clear();
		}
	}
	T getValue(const DPoint3d& pt) const {
		KdElement* top = m_top;
		DBox box = m_box;
		while (!top->is_leaf) {
			bool split_lower = pt[top->m_data.split.axis] <= top->m_data.split.value;
			box.splitAndUpdate(top->m_data.split.axis, top->m_data.split.value, split_lower);
			assert(box.contains(pt));
			top = split_lower ?	top->m_data.split.element_0 : top->m_data.split.element_1;
		}

		return top->getValue(pt, box);
	}
	inline const DBox& getBBox() const { return m_box; }
public: // stat
	void forRegularGrid(int grid_size, const std::function<void(const DPoint3d& pt)>& fg) {
		double fx = 1.0 / (double)grid_size;
		double dx = m_box.getDX() * fx;
		double dy = m_box.getDY() * fx;
		double dz = m_box.getDZ() * fx;

		DPoint3d pt(m_box.x0 + 0.5*dx, m_box.y0, m_box.z0);
		for (; pt.x < m_box.x1; pt.x += dx) {
			for (pt.y = m_box.y0 + 0.5*dy; pt.y < m_box.y1; pt.y += dy) {
				for (pt.z = m_box.z0 + 0.5*dz; pt.z < m_box.z1; pt.z += dz) {
					fg(pt);
				}
			}
		}
	}
	int getElementCount() const { return m_top ? m_top->getElementCount() : 0; }
	int getValueCount() const { return (int)m_values.countInt(); }
	int getMaxDepth() const { return m_top ? m_top->getMaxDepth() : 0; }
	int getMaxBalance() const {
		int max_balance = 0;
		if (m_top) m_top->getMaxDepthAndBalance(max_balance);
		return max_balance;
	}
	int getTotalBytes() const {
		return (int)sizeof(*this) // main struct
			+ (int)sizeof(KdElement) * getElementCount() // tree nodes
			+ (int)sizeof(T) * getValueCount(); // values
	}
private:
	class KdElement {
	public:
		struct KdVertices {
			T* vx_lsw = nullptr;
			T* vx_lse = nullptr;
			T* vx_lnw = nullptr;
			T* vx_lne = nullptr;
			T* vx_hsw = nullptr;
			T* vx_hse = nullptr;
			T* vx_hnw = nullptr;
			T* vx_hne = nullptr;
		};
	public:
		KdElement(const KdVertices& kv) : m_data(kv) { }
		int adaptToField(const std::function<T(const DPoint3d & pt)> & f,
			const DBox & box, double max_error, DataVector<T>& kvalues, int level = 0)
		{

			//static size_t counter = 0;
			//if (++counter % 100 == 0) {
			//	LOG4CPLUS_INFO(MeshLog::logger_console, "** counter", counter);
			//	LOG4CPLUS_INFO(MeshLog::logger_console, "** box-diam", box.getDiameter());
			//}

			assert(is_leaf);
			if (level >= MaxLevel) return 0;
			double curr_error = estimateError(box, f);
			if (curr_error <= max_error) return 0;

			KdVertices kv0 = m_data.values;
			KdVertices kv1 = m_data.values;

			is_leaf = false;
			KdSplitter::whereToSplit(box, f, m_data.split.axis, m_data.split.value);
			DBox box0, box1;
			box.split(m_data.split.axis, m_data.split.value, box0, box1);
			// update values_0 and values_1
			switch (m_data.split.axis) {
			case Axis::X:
				{
					assert(box0.x1 == m_data.split.value);
					T* v0 = &kvalues.get(
						kvalues.add(f(DPoint3d(box0.x1, box.y0, box.z0))));
					T* v1 = &kvalues.get(
						kvalues.add(f(DPoint3d(box0.x1, box.y1, box.z0))));
					T* v2 = &kvalues.get(
						kvalues.add(f(DPoint3d(box0.x1, box.y0, box.z1))));
					T* v3 = &kvalues.get(
						kvalues.add(f(DPoint3d(box0.x1, box.y1, box.z1))));
					// east for kv0
					kv0.vx_lse = v0;
					kv0.vx_lne = v1;
					kv0.vx_hse = v2;
					kv0.vx_hne = v3;
					// west for kv1
					kv1.vx_lsw = v0;
					kv1.vx_lnw = v1;
					kv1.vx_hsw = v2;
					kv1.vx_hnw = v3;
					break;
				}
			case Axis::Y:
				{
					assert(box0.y1 == m_data.split.value);
					T* v0 = &kvalues.get(
						kvalues.add(f(DPoint3d(box.x0, box0.y1, box.z0))));
					T* v1 = &kvalues.get(
						kvalues.add(f(DPoint3d(box.x1, box0.y1, box.z0))));
					T* v2 = &kvalues.get(
						kvalues.add(f(DPoint3d(box.x0, box0.y1, box.z1))));
					T* v3 = &kvalues.get(
						kvalues.add(f(DPoint3d(box.x1, box0.y1, box.z1))));
					// north for kv0
					kv0.vx_lnw = v0;
					kv0.vx_lne = v1;
					kv0.vx_hnw = v2;
					kv0.vx_hne = v3;
					// south for kv1
					kv1.vx_lsw = v0;
					kv1.vx_lse = v1;
					kv1.vx_hsw = v2;
					kv1.vx_hse = v3;
					break;
				}
			case Axis::Z:
				{
					assert(box0.z1 == m_data.split.value);
					T* v0 = &kvalues.get(
						kvalues.add(f(DPoint3d(box.x0, box.y0, box0.z1))));
					T* v1 = &kvalues.get(
						kvalues.add(f(DPoint3d(box.x1, box.y0, box0.z1))));
					T* v2 = &kvalues.get(
						kvalues.add(f(DPoint3d(box.x0, box.y1, box0.z1))));
					T* v3 = &kvalues.get(
						kvalues.add(f(DPoint3d(box.x1, box.y1, box0.z1))));
					// high for kv0
					kv0.vx_hsw = v0;
					kv0.vx_hse = v1;
					kv0.vx_hnw = v2;
					kv0.vx_hne = v3;
					// low for kv1
					kv1.vx_lsw = v0;
					kv1.vx_lse = v1;
					kv1.vx_lnw = v2;
					kv1.vx_lne = v3;
					break;
				}
			}
			m_data.split.element_0 = new KdElement(kv0);
			m_data.split.element_1 = new KdElement(kv1);
			return 2 +
				m_data.split.element_0->adaptToField(f, box0, max_error, kvalues, level+1) +
				m_data.split.element_1->adaptToField(f, box1, max_error, kvalues, level+1);
		}
		double estimateError(const DBox &box, const std::function<T(const DPoint3d & pt)> &f)
		{
			double error, max_error = 0.0;
			double fx = 1.0 / (double)SamplingSize;
			double dx = box.getDX() * fx;
			double dy = box.getDY() * fx;
			double dz = box.getDZ() * fx;

			DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
			for (; pt.x < box.x1; pt.x += dx) {
				for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
					for (pt.z = box.z0 + 0.5*dz; pt.z < box.z1; pt.z += dz) {
						error = diffKdValue(f(pt), getValue(pt, box));
						if (error > max_error) max_error = error;
					}
				}
			}

			return max_error;
		}
		T getValue(const DPoint3d& pt, const DBox& box) const
		{
			assert(box.contains(pt));

			double tx = (pt.x - box.x0) / box.getDX();
			double ty = (pt.y - box.y0) / box.getDY();
			double tz = (pt.z - box.z0) / box.getDZ();

			double coeff[8] = {
				(1 - tx) * (1 - ty) * (1 - tz),	// VX_LSW
				tx	 * (1 - ty) * (1 - tz),	// VX_LSE
				(1 - tx) *    ty	* (1 - tz),	// VX_LNW
				tx	 *    ty	* (1 - tz),	// VX_LNE
				(1 - tx) * (1 - ty) *    tz ,	// VX_HSW
				tx	 * (1 - ty) *    tz ,	// VX_HSE
				(1 - tx) *    ty	*    tz ,	// VX_HNW
				tx	 *    ty	*    tz 	// VX_HNE
			};
			return
				*(m_data.values.vx_lsw) * coeff[0] +
				*(m_data.values.vx_lse) * coeff[1] +
				*(m_data.values.vx_lnw) * coeff[2] +
				*(m_data.values.vx_lne) * coeff[3] +
				*(m_data.values.vx_hsw) * coeff[4] +
				*(m_data.values.vx_hse) * coeff[5] +
				*(m_data.values.vx_hnw) * coeff[6] +
				*(m_data.values.vx_hne) * coeff[7];
		}
		void clear() {
			if (!is_leaf) {
				m_data.split.element_0->clear();
				delete m_data.split.element_0;
				m_data.split.element_1->clear();
				delete m_data.split.element_1;
				is_leaf = true;
			}
		}
	public:
		int getElementCount() const {
			return is_leaf ? 1 :
				1 + m_data.split.element_0->getElementCount() + m_data.split.element_1->getElementCount();
		}
		int getMaxDepth() const {
			return is_leaf ? 1 :
				1 + std::max(m_data.split.element_0->getMaxDepth(), m_data.split.element_1->getMaxDepth());
		}
		int getMaxDepthAndBalance(int & max_balance) const {
			if (is_leaf) return 1;
			int left_depth = m_data.split.element_0->getMaxDepthAndBalance(max_balance);
			int right_depth = m_data.split.element_1->getMaxDepthAndBalance(max_balance);
			int balance = std::abs(left_depth - right_depth);
			if (balance > max_balance) max_balance = balance;
			return 1 + std::max(left_depth, right_depth);
		}
	public:
		bool is_leaf = true;
		union KdData {
			KdData(const KdVertices & kv) : values(kv) {}
			struct KdSplitElement {
				Axis		axis = Axis::X;
				double		value = 0.0;
				KdElement * element_0 = nullptr;
				KdElement * element_1 = nullptr;
			} split;
			KdVertices values;
		} m_data;
	};
private:
	DBox m_box;
	KdElement* m_top = nullptr;
	DataVector<T> m_values;
};

/**
* This class implements a kd-tree in a 3-dimensional space
*  (values in leaves + interpolation from neighbours) plus some basic operations.
*/
template<typename T,
	typename KdSplitter = Kd3dHalfLongestSplitter<T>,
	size_t SamplingSize = 4,
	size_t MaxLevel = 20>
	class DataKdTreeLi3d
{
public:
	DataKdTreeLi3d(const DBox & bbox) : m_box(bbox) {}
	~DataKdTreeLi3d() { clear(); }
private:
	DataKdTreeLi3d(const DataKdTreeLi3d& tree) = delete;
	DataKdTreeLi3d& operator=(const DataKdTreeLi3d& tree) = delete;
public:
	int adaptToField(const std::function<T(const DPoint3d & pt)> & f, double max_error) {
		T v = f(m_box.getMiddlePoint());
		if (m_top != nullptr) clear();
		typename KdElement::KdNeighbors kdn;
		m_top = new KdElement(v, 0, m_box, kdn, false);
		return 1 + m_top->adaptToField(f, max_error);
	}
	void clear() {
		if (m_top != nullptr) {
			m_top->clear();
			delete m_top;
			m_top = nullptr;
		}
	}
	T getValue(const DPoint3d& pt) const {
		KdElement* top = m_top;

		while (!top->is_leaf) {
			bool split_lower = pt[top->m_data.split.axis] <= top->m_data.split.value;
			top = split_lower ? top->m_data.split.element_0 : top->m_data.split.element_1;
		}

		return top->getValue(pt);
	}
	inline const DBox& getBBox() const { return m_box; }
public: // stat
	void forRegularGrid(int grid_size, const std::function<void(const DPoint3d& pt)>& fg) {
		double fx = 1.0 / (double)grid_size;
		double dx = m_box.getDX() * fx;
		double dy = m_box.getDY() * fx;
		double dz = m_box.getDZ() * fx;

		DPoint3d pt(m_box.x0 + 0.5*dx, m_box.y0, m_box.z0);
		for (; pt.x < m_box.x1; pt.x += dx) {
			for (pt.y = m_box.y0 + 0.5*dy; pt.y < m_box.y1; pt.y += dy) {
				for (pt.z = m_box.z0 + 0.5*dz; pt.z < m_box.z1; pt.z += dz) {
					fg(pt);
				}
			}
		}
	}
	int getElementCount() const { return m_top ? m_top->getElementCount() : 0; }
	int getValueCount() const { return m_top ? m_top->getValueCount() : 0; }
	int getMaxDepth() const { return m_top ? m_top->getMaxDepth() : 0; }
	int getMaxBalance() const {
		int max_balance = 0;
		if (m_top) m_top->getMaxDepthAndBalance(max_balance);
		return max_balance;
	}
	int getTotalBytes() const {
		return (int)sizeof(*this) // main struct
			+ (int)sizeof(KdElement) * getElementCount(); // tree nodes
	}
private:
	class KdElement {
	public:
		struct KdNeighbors {
			KdElement* low = nullptr;
			KdElement* high = nullptr;
			KdElement* south = nullptr;
			KdElement* north = nullptr;
			KdElement* east = nullptr;
			KdElement* west = nullptr;
		};
	public:
		KdElement(const T & v, int lvl, const DBox& bx, const KdNeighbors& kn, bool ready = true) 
			: m_data(v, lvl, bx, kn, ready) { }
		int adaptToField(const std::function<T(const DPoint3d & pt)> & f,
			double max_error, bool recursive = true)
		{
			//static size_t counter = 0;
			//if (++counter % 100 == 0) {
			//	LOG4CPLUS_INFO(MeshLog::logger_console, "** counter", counter);
			//	LOG4CPLUS_INFO(MeshLog::logger_console, "** box-diam", box.getDiameter());
			//}

			int result = 0;

			if (is_leaf) {
				if (m_data.leaf.is_adapted) return result;
				if (m_data.leaf.level >= MaxLevel) return result;

				double curr_error = estimateError(f);
				if (curr_error <= max_error) {
					m_data.leaf.is_adapted = true;
					return result;
				}

				int split_level = m_data.leaf.level;
				DBox split_box = m_data.leaf.box;
				KdNeighbors split_kn = m_data.leaf.kn;

				is_leaf = false;
				KdSplitter::whereToSplit(split_box, f, m_data.split.axis, m_data.split.value);
				DBox box0, box1;
				split_box.split(m_data.split.axis, m_data.split.value, box0, box1);

				KdNeighbors kn0 = split_kn;
				KdNeighbors kn1 = split_kn;

				T v0 = f(box0.getMiddlePoint());
				T v1 = f(box1.getMiddlePoint());

				m_data.split.element_0 = new KdElement(v0, split_level + 1, box0, kn0, false);
				m_data.split.element_1 = new KdElement(v1, split_level + 1, box1, kn1, false);
				result += 2;

				switch (m_data.split.axis) {
				case Axis::X:
					kn0.east = m_data.split.element_1;
					kn1.west = m_data.split.element_0;
					break;
				case Axis::Y:
					kn0.north = m_data.split.element_1;
					kn1.south = m_data.split.element_0;
					break;
				case Axis::Z:
					kn0.high = m_data.split.element_1;
					kn1.low = m_data.split.element_0;
					break;
				}

				m_data.split.element_0->setKN(kn0); // updated
				m_data.split.element_1->setKN(kn1); // updated

				KdElement* kd_elements[6] = { 
					split_kn.east, split_kn.west, 
					split_kn.north, split_kn.south,
					split_kn.low, split_kn.high };

				DPoint3d pt = split_box.getMiddlePoint();
				for (KdElement* kde : kd_elements) {
					if (kde) {
						KdElement* nleaf = kde->getNearestLeaf(pt);
						assert(nleaf && nleaf->is_leaf);
						if (!nleaf->m_data.leaf.is_adapted && nleaf->m_data.leaf.level < split_level)
							result += nleaf->adaptToField(f, max_error, false);
					}
				}
			}

			if (recursive) {
				result += m_data.split.element_0->adaptToField(f, max_error);
				result += m_data.split.element_1->adaptToField(f, max_error);
			}

			return result;
		}
		double estimateError(const std::function<T(const DPoint3d & pt)> &f)
		{
			assert(is_leaf);
			double error, max_error = 0.0;
			double fx = 1.0 / (double)SamplingSize;
			double dx = m_data.leaf.box.getDX() * fx;
			double dy = m_data.leaf.box.getDY() * fx;
			double dz = m_data.leaf.box.getDZ() * fx;

			DPoint3d pt(m_data.leaf.box.x0 + 0.5*dx, m_data.leaf.box.y0, m_data.leaf.box.z0);
			for (; pt.x < m_data.leaf.box.x1; pt.x += dx) {
				for (pt.y = m_data.leaf.box.y0 + 0.5*dy; pt.y < m_data.leaf.box.y1; pt.y += dy) {
					for (pt.z = m_data.leaf.box.z0 + 0.5*dz; pt.z < m_data.leaf.box.z1; pt.z += dz) {
						error = diffKdValue(f(pt), getValue(pt));
						if (error > max_error) max_error = error;
					}
				}
			}

			return max_error;
		}
		KdElement* getNearestLeaf(const DPoint3d& pt) {
			KdElement* element = this;
			while (!element->is_leaf) {
				element = (pt[element->m_data.split.axis] <= element->m_data.split.value) ?
					element->m_data.split.element_0 : element->m_data.split.element_1;
			}
			return element;
		}
		T getValue(const DPoint3d& pt) const
		{
			assert(is_leaf);
			assert(m_data.leaf.box.contains(pt));

			double tx = (pt.x - m_data.leaf.box.x0) / m_data.leaf.box.getDX();
			double ty = (pt.y - m_data.leaf.box.y0) / m_data.leaf.box.getDY();
			double tz = (pt.z - m_data.leaf.box.z0) / m_data.leaf.box.getDZ();

			double kx = 2 * tx - 1.0; // [-1, 1]
			double ky = 2 * ty - 1.0; // [-1, 1]
			double kz = 2 * tz - 1.0; // [-1, 1]
			double kx2 = kx*kx;
			double ky2 = ky*ky;
			double kz2 = kz*kz;

			double N[28] = { 0.0 }; // shape functions for 27-node hexahedral

			N[1] = 0.125 * (1 - tx)*(1 - ty)*(1 - tz);
			N[2] = 0.125 * (1 + tx)*(1 - ty)*(1 - tz);
			N[3] = 0.125 * (1 + tx)*(1 + ty)*(1 - tz);
			N[4] = 0.125 * (1 - tx)*(1 + ty)*(1 - tz);
			N[5] = 0.125 * (1 - tx)*(1 - ty)*(1 + tz);
			N[6] = 0.125 * (1 + tx)*(1 - ty)*(1 + tz);
			N[7] = 0.125 * (1 + tx)*(1 + ty)*(1 + tz);
			N[8] = 0.125 * (1 - tx)*(1 + ty)*(1 + tz);

			N[9]  = 0.25 * (1 - kx2) * (1 - ky)  * (1 - kz);
			N[10] = 0.25 * (1 + kx)  * (1 - ky2) * (1 - kz);
			N[11] = 0.25 * (1 - kx2) * (1 + ky)  * (1 - kz);
			N[12] = 0.25 * (1 - kx)  * (1 - ky2) * (1 - kz);
			N[13] = 0.25 * (1 - kx2) * (1 - ky)  * (1 + kz);
			N[14] = 0.25 * (1 + kx)  * (1 - ky2) * (1 + kz);
			N[15] = 0.25 * (1 - kx2) * (1 + ky)  * (1 + kz);
			N[16] = 0.25 * (1 - kx)  * (1 - ky2) * (1 + kz);
			N[17] = 0.25 * (1 - kx)  * (1 - ky)  * (1 - kz2);
			N[18] = 0.25 * (1 + kx)  * (1 - ky)  * (1 - kz2);
			N[19] = 0.25 * (1 + kx)  * (1 + ky)  * (1 - kz2);
			N[20] = 0.25 * (1 - kx)  * (1 + ky)  * (1 - kz2);

			N[1] -= 0.5 * (N[12] + N[9]  + N[17]);
			N[2] -= 0.5 * (N[9]  + N[10] + N[18]);
			N[3] -= 0.5 * (N[10] + N[11] + N[19]);
			N[4] -= 0.5 * (N[11] + N[12] + N[20]);
			N[5] -= 0.5 * (N[16] + N[13] + N[17]);
			N[6] -= 0.5 * (N[13] + N[14] + N[18]);
			N[7] -= 0.5 * (N[14] + N[15] + N[19]);
			N[8] -= 0.5 * (N[15] + N[16] + N[20]);

			N[21] = 0.5 * (1 - kx2) * (1 - ky2) * (1 - kz);
			N[22] = 0.5 * (1 - kx2) * (1 - ky2) * (1 + kz);
			N[23] = 0.5 * (1 - kx2) * (1 - ky)  * (1 - kz2);
			N[24] = 0.5 * (1 + kx)  * (1 - ky2) * (1 - kz2);
			N[25] = 0.5 * (1 - kx2) * (1 + ky)  * (1 - kz2);
			N[26] = 0.5 * (1 - kx)  * (1 - ky2) * (1 - kz2);

			N[27] = (1 - kx2)*(1 - ky2)*(1 - kz2);

			N[1] += 0.25 * (N[21] + N[26] + N[23]);
			N[2] += 0.25 * (N[21] + N[23] + N[24]);
			N[3] += 0.25 * (N[21] + N[24] + N[25]);
			N[4] += 0.25 * (N[21] + N[25] + N[26]);
			N[5] += 0.25 * (N[22] + N[26] + N[23]);
			N[6] += 0.25 * (N[22] + N[23] + N[24]);
			N[7] += 0.25 * (N[22] + N[24] + N[25]);
			N[8] += 0.25 * (N[22] + N[25] + N[26]);

			N[9]  -= 0.5 * (N[21] + N[23]);
			N[10] -= 0.5 * (N[21] + N[24]);
			N[11] -= 0.5 * (N[21] + N[25]);
			N[12] -= 0.5 * (N[21] + N[26]);
			N[13] -= 0.5 * (N[22] + N[23]);
			N[14] -= 0.5 * (N[22] + N[24]);
			N[15] -= 0.5 * (N[22] + N[25]);
			N[16] -= 0.5 * (N[22] + N[26]);
			N[17] -= 0.5 * (N[26] + N[23]);
			N[18] -= 0.5 * (N[23] + N[24]);
			N[19] -= 0.5 * (N[24] + N[25]);
			N[20] -= 0.5 * (N[25] + N[26]);

			for (int i = 1; i <= 8; i++) N[i] -= 0.125 * N[27];
			for (int i = 9; i <= 20; i++) N[i] += 0.25 * N[27];
			for (int i = 21; i <= 26; i++) N[i] -= 0.5 * N[27];

			const T& vc = m_data.leaf.value;
			const T& vx0 = m_data.leaf.kn.west ? m_data.leaf.kn.west->getNearestLeaf(pt)->m_data.leaf.value : vc;
			const T& vx1 = m_data.leaf.kn.east ? m_data.leaf.kn.east->getNearestLeaf(pt)->m_data.leaf.value : vc;
			const T& vy0 = m_data.leaf.kn.south ? m_data.leaf.kn.south->getNearestLeaf(pt)->m_data.leaf.value : vc;
			const T& vy1 = m_data.leaf.kn.north ? m_data.leaf.kn.north->getNearestLeaf(pt)->m_data.leaf.value : vc;
			const T& vz0 = m_data.leaf.kn.low ? m_data.leaf.kn.low->getNearestLeaf(pt)->m_data.leaf.value : vc;
			const T& vz1 = m_data.leaf.kn.high ? m_data.leaf.kn.high->getNearestLeaf(pt)->m_data.leaf.value : vc;

			T P[28];

			//P[27] = vc;

			P[21] = (vz0 + vc) * 0.5;
			P[22] = (vz1 + vc) * 0.5;
			P[23] = (vy0 + vc) * 0.5;
			P[24] = (vx1 + vc) * 0.5;
			P[25] = (vy1 + vc) * 0.5;
			P[26] = (vx0 + vc) * 0.5;

			P[9]  = (P[21] + P[23]) * 0.5;
			P[10] = (P[21] + P[24]) * 0.5;
			P[11] = (P[21] + P[25]) * 0.5;
			P[12] = (P[21] + P[26]) * 0.5;
			P[13] = (P[22] + P[23]) * 0.5;
			P[14] = (P[22] + P[24]) * 0.5;
			P[15] = (P[22] + P[25]) * 0.5;
			P[16] = (P[22] + P[26]) * 0.5;
			P[17] = (P[26] + P[23]) * 0.5;
			P[18] = (P[23] + P[24]) * 0.5;
			P[19] = (P[24] + P[25]) * 0.5;
			P[20] = (P[25] + P[26]) * 0.5;

			P[1] = (P[21] + P[26] + P[23]) * (1.0 / 3.0);
			P[2] = (P[21] + P[23] + P[24]) * (1.0 / 3.0);
			P[3] = (P[21] + P[24] + P[25]) * (1.0 / 3.0);
			P[4] = (P[21] + P[26] + P[23]) * (1.0 / 3.0);
			P[5] = (P[22] + P[26] + P[23]) * (1.0 / 3.0);
			P[6] = (P[22] + P[23] + P[24]) * (1.0 / 3.0);
			P[7] = (P[22] + P[24] + P[25]) * (1.0 / 3.0);
			P[8] = (P[22] + P[25] + P[26]) * (1.0 / 3.0);


			T result = vc * N[27];
			for (int i = 1; i <= 26; i++)
				result += P[i] * N[i];
			return result;
			/*
			T result = m_data.leaf.value * 0.5;

			T ax = m_data.leaf.value;
			T ay = m_data.leaf.value;
			T az = m_data.leaf.value;

			if (tx < 0.5) {
				if (m_data.leaf.kn.west) {
					tx *= 2;
					ax *= tx;
					ax += m_data.leaf.kn.west->getNearestLeaf(pt)->m_data.leaf.value * (1.0 - tx);
				}
			}
			else {
				if (m_data.leaf.kn.east) {
					tx = (1.0 - tx) * 2;
					ax *= tx;
					ax += m_data.leaf.kn.east->getNearestLeaf(pt)->m_data.leaf.value * (1.0 - tx);
				}

			}

			if (ty < 0.5) {
				if (m_data.leaf.kn.south) {
					ty *= 2;
					ay *= ty;
					ay += m_data.leaf.kn.south->getNearestLeaf(pt)->m_data.leaf.value * (1.0 - ty);
				}
			}
			else {
				if (m_data.leaf.kn.north) {
					ty = (1.0 - ty) * 2;
					ay *= ty;
					ay += m_data.leaf.kn.north->getNearestLeaf(pt)->m_data.leaf.value * (1.0 - ty);
				}
			}

			if (tz < 0.5) {
				if (m_data.leaf.kn.low) {
					tz *= 2;
					az *= tz;
					az += m_data.leaf.kn.low->getNearestLeaf(pt)->m_data.leaf.value * (1.0 - tz);
				}
			}
			else {
				if (m_data.leaf.kn.high) {
					tz = (1.0 - tz) * 2;
					az *= tz;
					az += m_data.leaf.kn.high->getNearestLeaf(pt)->m_data.leaf.value * (1.0 - tz);
				}

			}

			result += (ax + ay + az) / 6;
			*/

			//T v_east = kn.east ? kn.east->getValueFromNearestLeaf(pt) : m_data.value;
			//T v_west = kn.west ? kn.west->getValueFromNearestLeaf(pt) : m_data.value;
			//T v_south = kn.south ? kn.south->getValueFromNearestLeaf(pt) : m_data.value;
			//T v_north = kn.north ? kn.north->getValueFromNearestLeaf(pt) : m_data.value;
			//T v_low  = kn.low  ? kn.low->getValueFromNearestLeaf(pt)  : m_data.value;
			//T v_high = kn.high ? kn.high->getValueFromNearestLeaf(pt) : m_data.value;

			return result;
		}
		void clear() {
			if (!is_leaf) {
				m_data.split.element_0->clear();
				delete m_data.split.element_0;
				m_data.split.element_1->clear();
				delete m_data.split.element_1;
				is_leaf = true;
			}
		}
		void setKN(const KdNeighbors& kn) {
			assert(is_leaf);
			m_data.leaf.kn = kn;
		}
	public:
		int getElementCount() const {
			return is_leaf ? 1 :
				1 + m_data.split.element_0->getElementCount() + m_data.split.element_1->getElementCount();
		}
		int getValueCount() const {
			return is_leaf ? 1 :
				m_data.split.element_0->getValueCount() + m_data.split.element_1->getValueCount();
		}
		int getMaxDepth() const {
			return is_leaf ? 1 :
				1 + std::max(m_data.split.element_0->getMaxDepth(), m_data.split.element_1->getMaxDepth());
		}
		int getMaxDepthAndBalance(int & max_balance) const {
			if (is_leaf) return 1;
			int left_depth = m_data.split.element_0->getMaxDepthAndBalance(max_balance);
			int right_depth = m_data.split.element_1->getMaxDepthAndBalance(max_balance);
			int balance = std::abs(left_depth - right_depth);
			if (balance > max_balance) max_balance = balance;
			return 1 + std::max(left_depth, right_depth);
		}
	public:
		bool is_leaf = true;
		union KdData {
			KdData(const T & v, int lvl, const DBox& box, const KdNeighbors& kn, bool adapted) 
				: leaf(v, lvl, box, kn, adapted) {}
			struct KdSplitElement {
				Axis		axis = Axis::X;
				double		value = 0.0;
				KdElement * element_0 = nullptr;
				KdElement * element_1 = nullptr;
			} split;
			struct KdLeaf {
				KdLeaf(const T& v, int lvl, const DBox& bx, const KdNeighbors& _kn, bool adapted = true) 
					: value(v), level(lvl), box(bx), kn(_kn), is_adapted(adapted) {}
				T value;
				int level = 0;
				DBox box;
				KdNeighbors kn;
				bool is_adapted = true;
			} leaf;
		} m_data;
	};
private:
	DBox m_box;
	KdElement* m_top = nullptr;
};

#endif // !defined(DATAKDTREE_H__INCLUDED)
