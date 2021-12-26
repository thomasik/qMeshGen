/////////////////////////////////////////////////////////////////////////////
// DTetrahedron.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DTETRAHEDRON_H__INCLUDED
#define DTETRAHEDRON_H__INCLUDED

#include "DPoint.h"
#include "DTriangle.h"

/**
 * This class implements a tetrahedron (four points) in a three-dimensional space
 *  plus some basic operations
 */
class DTetrahedron
{
public:
	/// Standard constructor
	DTetrahedron(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d) 
		: pt_a(a), pt_b(b), pt_c(c), pt_d(d) {}
public:
	/// Insphere 3d
	static double insphere(const DPoint3d& pt, const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d);
	/// Insphere 3d
	static double insphere(const DMPoint3d& pt, const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d);
	/// Volume of tetrahedron
	static double volume(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d){
		return (1.0/6.0) * DTriangle3d::orient3d(a, b, c, d);
	}
	/// Volume of tetrahedron
	static double volume(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d){
		return (1.0/6.0) * DMTriangle3d::orient3d(a, b, c, d);
	}

	/// Returns the quality coefficient of the tetrahedron
	static double aspectRatio(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d);
	/// Returns the quality coefficient of the tetrahedron
	static double aspectRatio(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d);

	/// Returns the length of the radius of the circumscribed sphere (for a tetrahedron given by four points)
	static double outerSphereRadius(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d);
	/// Returns the squared length of the radius of the circumscribed sphere (for a tetrahedron given by four points)
	static double outerSphereRadius2(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d);
	/// Returns the length of the radius of the circumscribed sphere (for a tetrahedron given by four points)
	static double outerSphereRadius(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d);
	/// Returns the squared length of the radius of the circumscribed sphere (for a tetrahedron given by four points)
	static double outerSphereRadius2(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d);
	/// Returns the middle of the circumscribed sphere (for a tetrahedron given by four points)
	static DPoint3d outerSphereCenter(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d);
	/// Returns the middle of the circumscribed sphere (for a tetrahedron given by four points)
	static DMPoint3d outerSphereCenter(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d);

public:
	/// points
	DPoint3d pt_a, pt_b, pt_c, pt_d; 
};

#endif // !defined(DTETRAHEDRON_H__INCLUDED)
