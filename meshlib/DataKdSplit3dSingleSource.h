/////////////////////////////////////////////////////////////////////////////
// DataKdSplit3dSingleSource.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DATAKD3DSPLITSS_H__INCLUDED
#define DATAKD3DSPLITSS_H__INCLUDED

#include "common.h"

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "DMetric3d.h"
#include "DataKdSplit3d.h"

class Kd3dSplitterSingleSource : public Kd3dSplitter
{
public:
	Kd3dSplitterSingleSource(int _csize = 0) : Kd3dSplitter(_csize) { }
public:
	virtual void whereToSplit(const DBox & box,
		const DPoint3d& pt, const DVector3d& dlen2,
		Axis & split_axis, double & split_value) = 0;
};

/**
* An example implementation of the KdSplitter.
* Splits the box in the middle along the longest axis.
*/
class Kd3dSSHalfLongestSplitter : public Kd3dSplitterSingleSource
{
public:
	void whereToSplit(const DBox & box, 
		const DPoint3d& /*pt*/, const DVector3d& /*dlen2*/,
		Axis & split_axis, double & split_value) override
	{
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
	string getLabel() const override { return "HLong"; }
};

class Kd3dSSGoldenRatioLongestSplitter : public Kd3dSplitterSingleSource
{
public:
	void whereToSplit(const DBox & box,
		const DPoint3d& pt, const DVector3d& dlen2,
		Axis & split_axis, double & split_value) override;
	string getLabel() const override { return "GoldenRatio"; }
};

class Kd3dSSMinLongestSplitter : public Kd3dSplitterSingleSource
{
public:
	void whereToSplit(const DBox & box,
		const DPoint3d& pt, const DVector3d& dlen2,
		Axis & split_axis, double & split_value) override;
	string getLabel() const override { return "MinLong"; }
};
#endif // !defined(DATAKD3DSPLITSS_H__INCLUDED)
