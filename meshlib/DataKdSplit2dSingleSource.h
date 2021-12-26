/////////////////////////////////////////////////////////////////////////////
// DataKdSplit2dSingleSource.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DATAKD2DSPLITSS_H__INCLUDED
#define DATAKD2DSPLITSS_H__INCLUDED

#include "common.h"

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "DMetric2d.h"
#include "DataKdSplit2d.h"

class Kd2dSplitterSingleSource : public Kd2dSplitter
{
public:
	Kd2dSplitterSingleSource(int _csize = 0) : Kd2dSplitter(_csize) { }
public:
	virtual void whereToSplit(const DRect & box,
		const DPoint2d& pt, const DVector2d& dlen2,
		Axis & split_axis, double & split_value) = 0;
};

/**
* An example implementation of the KdSplitter.
* Splits the box in the middle along the longest axis.
*/
class Kd2dSSHalfLongestSplitter : public Kd2dSplitterSingleSource
{
public:
	void whereToSplit(const DRect & box,
		const DPoint2d& /*pt*/, const DVector2d& /*dlen2*/,
		Axis & split_axis, double & split_value) override
	{
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
};

class Kd2dSSGoldenRatioLongestSplitter : public Kd2dSplitterSingleSource
{
public:
	void whereToSplit(const DRect & box,
		const DPoint2d& pt, const DVector2d& dlen2,
		Axis & split_axis, double & split_value) override;
};

class Kd2dSSMinLongestSplitter : public Kd2dSplitterSingleSource
{
public:
	void whereToSplit(const DRect & box,
		const DPoint2d& pt, const DVector2d& dlen2,
		Axis & split_axis, double & split_value) override;
};

#endif // !defined(DATAKD2DSPLITSS_H__INCLUDED)
