/////////////////////////////////////////////////////////////////////////////
// DQuadric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DQUADRIC_H__INCLUDED
#define DQUADRIC_H__INCLUDED

#include "DVectorN.h"

class DPoint3d;
class DVector3d;

/**
 * This class implements a quadric surface representation
 */
class DQuadric
{
public:
	static const int N = 10;
public:
	/// Standard constructor
	DQuadric();
	DQuadric(const DVectorN<N> & v);
	// v=(1,x,y,z,xx,yy,zz,xy,xz,yz)
public:
	/// count distance of point and quadric
	double distance(const DPoint3d& pt) const;
	/// solve Q( pt + t * vt) == 0, with possible two solutions (or some other cases)
	bool solve(const DPoint3d& pt, const DVector3d& vt, double& t1, double &t2) const;
public:
	/// coefficients
	DVectorN<N> vq; 
};

#endif // !defined(DQUADRIC_H__INCLUDED)
