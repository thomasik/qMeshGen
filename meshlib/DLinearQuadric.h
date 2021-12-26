/////////////////////////////////////////////////////////////////////////////
// DLinearQuadric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DLINEARQUADRIC_H__INCLUDED
#define DLINEARQUADRIC_H__INCLUDED

#include "DVectorN.h"
#include "DLine.h"
#include "DataMatrix.h"

/**
 * This class implements a quadric curve representation for line: y(x)
 */
class DLinearQuadric
{
public:
	static const int N = 3;
	//static const int SKETCH_LINES = 10;
public:
	/// Standard constructor
	DLinearQuadric();
	// v=(1,y,yy)
public:
	/// point on the curve
	DPoint2d getPoint(const double& t) const;
	/// count distance of point and quadric
	double distance(const DPoint2d& pt) const;
	/// Stores the parameters of this linear quadric into the stream
	friend ostream& operator<<(ostream& os, const DLinearQuadric& lq) {
		return os << lq.line << " " << lq.vq;
	}
	/// Restores the parameters of this linear quadric from the stream
	friend istream& operator>>(istream& is, DLinearQuadric& lq) {
		return is >> lq.line >> ws >> lq.vq;
	}
public:
	/// coefficients
	DVectorN<N> vq;
	DLine2d line;
};

/**
 * This class implements a quadric curve representation for line: y1(x), y2(x)
 */
class DLinearQuadric3d
{
public:
	static const int N = 3;
	//static const int SKETCH_LINES = 10;
public:
	/// Standard constructor
	DLinearQuadric3d();
	// v=(1,y,yy)
public:
	/// point on the curve
	DPoint3d getPoint(const double& t) const;
	double getParameter(const DPoint3d& pt) const { return line.paramForPoint(pt); }
	/// count distance of point and quadric
	double distance(const DPoint3d& pt) const;
	/// Stores the parameters of this linear quadric into the stream
	friend ostream& operator<<(ostream& os, const DLinearQuadric3d& lq) {
		return os << lq.line << " " << lq.e1 << " " << lq.vq1 << " " << lq.e2 << " " << lq.vq2;
	}
	/// Restores the parameters of this linear quadric from the stream
	friend istream& operator>>(istream& is, DLinearQuadric3d& lq) {
		return is >> lq.line >> ws >> lq.e1 >> ws >> lq.vq1 >> ws >> lq.e2 >> ws >> lq.vq2;
	}
public:
	/// coefficients
	DVectorN<N> vq1, vq2;
	DVector3d e1,e2;
	DLine3d line;
};

#endif // !defined(DLINEARQUADRIC_H__INCLUDED)
