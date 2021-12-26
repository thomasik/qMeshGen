/////////////////////////////////////////////////////////////////////////////
// DPlane.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DPLANE_H__INCLUDED
#define DPLANE_H__INCLUDED

#include "DPoint.h"
#include "DVector.h"
#include "DOrientedBox.h"
#include "DataMatrix.h"

/**
 * This class implements a plane representation
 */
class DPlane
{
public:
	/// Standard constructor
	DPlane();
	DPlane(const DPoint3d& pt, const DVector3d& e0, const DVector3d& e1);
	DPlane(const DPoint3d& pt, const DVector3d& vn);
public:
	DOrientedBox getOrientedBox() const;
	DOrientedBox getOrientedBox( const DataVector<DPoint3d>& points ) const;
	DOrientedBox getOrientedBoxOpt( const DataVector<DPoint3d>& points ) const;
	void countNormal();
	bool switchOrientation();
	DPoint2d projectToPlane(const DPoint3d& pt) const;
	bool projectToPlane(const DPoint3d& pt, DPoint2d& param, double& z) const;
	DPoint3d projectToSpace(const DPoint2d& pt) const;
	/// Stores the parameters of this plane into the stream
	friend ostream& operator<<(ostream& os, const DPlane& plane);
	/// Restores the parameters of this plane from the stream
	friend istream& operator>>(istream& is, DPlane& plane);
public:
	/// coefficients
	DPoint3d p0;
	DVector3d e0, e1;
	DVector3d vn;
};

#endif // !defined(DPLANE_H__INCLUDED)
