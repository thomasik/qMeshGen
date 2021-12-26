/////////////////////////////////////////////////////////////////////////////
// DataKdSplit3dDiscField.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DATAKD3DSPLITDF_H__INCLUDED
#define DATAKD3DSPLITDF_H__INCLUDED

#include "common.h"

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "DMetric3d.h"
#include "DataKdSplit3d.h"

class Kd3dSplitterDiscField : public Kd3dSplitter
{
public:
	Kd3dSplitterDiscField(int _csize = 0) : Kd3dSplitter(_csize) { }
public:
	virtual void whereToSplit(const DBox & box,
		/* ??? */
		Axis & split_axis, double & split_value) = 0;
};

/**
* An example implementation of the KdSplitter.
*
* Splits the box in the middle along the longest axis.
*/
class Kd3dDFHalfLongestSplitter : public Kd3dSplitterDiscField
{
public:
	void whereToSplit(const DBox & box, 
		/* ??? */
		Axis & split_axis, double & split_value)
	{
		split_axis = box.getLongestAxis();
		split_value = box.getMiddleValue(split_axis);
	}
	string getLabel() const override { return "HLong"; }
};

#endif // !defined(DATAKD3DSPLITDF_H__INCLUDED)
