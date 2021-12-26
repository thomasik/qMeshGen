/////////////////////////////////////////////////////////////////////////////
// DataKdSplit3dContField.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DATAKD3DSPLITCF_H__INCLUDED
#define DATAKD3DSPLITCF_H__INCLUDED

#include "common.h"

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "DMetric3d.h"
#include "DataKdSplit3d.h"

class Kd3dSplitterContField : public Kd3dSplitter
{
public:
	Kd3dSplitterContField(int _csize = 0) : Kd3dSplitter(_csize) { }
public:
	virtual void whereToSplit(const DBox & box,
		const std::function<CDM3d(const DPoint3d & pt)> & f,
		Axis & split_axis, double & split_value) = 0;
};

/**
* An example implementation of the KdSplitter.
*
* Splits the box in the middle along the longest axis.
*/
class Kd3dCFHalfLongestSplitter : public Kd3dSplitterContField
{
public:
	void whereToSplit(const DBox & box, 
		const std::function<CDM3d(const DPoint3d & pt)> & /* f */,
		Axis & split_axis, double & split_value) override
	{
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
	string getLabel() const override { return "HLong"; }
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
class Kd3dCFMaxGradientSumSplitter : public Kd3dSplitterContField
{
private:
public:
	Kd3dCFMaxGradientSumSplitter(int _samples_size = 4) : Kd3dSplitterContField(_samples_size) {
		assert(_samples_size > 1);
	}
public:
	void whereToSplit(const DBox &box, 
		const std::function<CDM3d(const DPoint3d & pt)> & f,
		Axis & split_axis, double & split_value) override;
	string getLabel() const override;
};

class Kd3dCFMaxGradientSplitter : public Kd3dSplitterContField
{
public:
	Kd3dCFMaxGradientSplitter(int _samples_size = 4) : Kd3dSplitterContField(_samples_size) {
		assert(_samples_size > 1);
	}
public:
	void whereToSplit(const DBox &box, 
		const std::function<CDM3d(const DPoint3d & pt)> & f,
		Axis & split_axis, double & split_value) override;
	string getLabel() const override;
};

class Kd3dCFMaxDiffGradientSplitter : public Kd3dSplitterContField
{
public:
	Kd3dCFMaxDiffGradientSplitter(int _samples_size = 4) : Kd3dSplitterContField(_samples_size) {
		assert(_samples_size > 1);
	}
public:
	void whereToSplit(const DBox &box, 
		const std::function<CDM3d(const DPoint3d & pt)> & f,
		Axis & split_axis, double & split_value) override;
	string getLabel() const override;
};

#endif // !defined(DATAKD3DSPLITCF_H__INCLUDED)
