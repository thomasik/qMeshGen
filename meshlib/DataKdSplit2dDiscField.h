/////////////////////////////////////////////////////////////////////////////
// DataKdSplit2dDiscField.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DATAKD2DSPLITDF_H__INCLUDED
#define DATAKD2DSPLITDF_H__INCLUDED

#include "common.h"

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "DMetric2d.h"
#include "DataKdSplit2d.h"

class Kd2dSplitterDiscField : public Kd2dSplitter
{
public:
	Kd2dSplitterDiscField(int _csize = 0) : Kd2dSplitter(_csize) { }
public:
	virtual void whereToSplit(const DRect & box,
		/* ??? */
		Axis & split_axis, double & split_value) = 0;
};

/**
* An example implementation of the KdSplitter.
*
* Splits the box in the middle along the longest axis.
*/
class Kd2dDFHalfLongestSplitter : public Kd2dSplitterDiscField
{
public:
	void whereToSplit(const DRect & box, 
		/* ??? */
		Axis & split_axis, double & split_value)
	{
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
};

#endif // !defined(DATAKD2DSPLITDF_H__INCLUDED)
