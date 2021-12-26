/////////////////////////////////////////////////////////////////////////////
// DLinearQuadric.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DLinearQuadric.h"
#include "DVectorN.h"
#include "DPoint.h"

/// Standard constructor
DLinearQuadric::DLinearQuadric() : vq(0.0) { }

/// count distance of point and quadric
double DLinearQuadric::distance(const DPoint2d& pt) const
{
	double t = line.paramForPoint(pt);
	return pt.distance( getPoint(t) );
}

DPoint2d DLinearQuadric::getPoint(const double& t) const
{
	const double v[N] = { 1,  t, t*t };

	double dist = 0.0;
	for(int i = 0; i < N; i++) dist += v[i] * vq[i];
	return line.getPoint(t) + line.m_vn * dist;
}

/// Standard constructor
DLinearQuadric3d::DLinearQuadric3d() : vq1(0.0), vq2(0.0) { }

/// count distance of point and quadric
double DLinearQuadric3d::distance(const DPoint3d& pt) const
{
	double t = line.paramForPoint(pt);
	return pt.distance( getPoint(t) );
}

DPoint3d DLinearQuadric3d::getPoint(const double& t) const
{
	const double v[N] = { 1,  t, t*t };

	double d1 = 0.0, d2 = 0.0;
	for(int i = 0; i < N; i++){
		d1 += v[i] * vq1[i];
		d2 += v[i] * vq2[i];
	}
	return line.getPoint(t) + e1 * d1 + e2 * d2;
}
