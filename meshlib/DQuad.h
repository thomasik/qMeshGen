/////////////////////////////////////////////////////////////////////////////
// DQuad.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DQUAD_H__INCLUDED
#define DQUAD_H__INCLUDED

#include "DPoint.h"

/**
 * This class implements a quad (four points) in a two-dimensional space
 *  plus some basic operations.
 */
class DQuad2d
{
public:
	/// Standard constructor
	DQuad2d(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c, const DPoint2d& d) 
		: pt_a(a), pt_b(b), pt_c(c), pt_d(d) {}
public:
	/// Returns the area of the triangle
//	double area() const { 
//		return 0.5 * ((pt_b.x - pt_a.x)*(pt_c.y - pt_a.y) - ((pt_b.y - pt_a.y)*(pt_c.x - pt_a.x))); }
	/// Returns the area of the triangle
//	static double area(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c) { 
//		return 0.5 * ((b.x - a.x)*(c.y - a.y) - ((b.y - a.y)*(c.x - a.x))); }
	/// Returns the area of the triangle
//	static double area(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c) { 
//		return 0.5 * ((b.x - a.x)*(c.y - a.y) - ((b.y - a.y)*(c.x - a.x))); }
	/// Returns the area of the triangle
//	static double det(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c) { 
//		return (b.x - a.x)*(c.y - a.y) - ((b.y - a.y)*(c.x - a.x)); }
	/// Returns the area of the triangle
//	static double det(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c) { 
//		return (b.x - a.x)*(c.y - a.y) - ((b.y - a.y)*(c.x - a.x)); }

	/// Returns the "alpha"-quality coefficient of the quad
	static double alphaQuality(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c, const DPoint2d& d);
	/// Returns the "alpha"-quality coefficient of the quad
	static double alphaQuality(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c, const DMPoint2d& d);

	/// Whether is valid (not inverted or self-crossing)
	static bool valid(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c, const DPoint2d& d);

public:
	/// points
	DPoint2d pt_a, pt_b, pt_c, pt_d; 
};

/**
 * This class implements a quad (four points) in a three-dimensional space
 *  plus some basic operations
 */
class DQuad3d
{
public:
	/// Standard constructor
	DQuad3d(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d) 
		: pt_a(a), pt_b(b), pt_c(c), pt_d(d) {}
public:
	/// Returns the area of the triangle
//	double area() const { 
//		return 0.5 * (pt_b-pt_a).crossProduct(pt_c-pt_a).length(); }
	/// Returns the area of the triangle
//	static double area(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c) { 
//		return 0.5 * (b-a).crossProduct(c-a).length(); }

	// Returns the "alpha"-quality coefficient of the quad
	static double alphaQuality(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d);
	// Returns the "shape"-quality coefficient of the quad
	static double shapeQuality(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d);
public:
	/// points
	DPoint3d pt_a, pt_b, pt_c, pt_d; 
};

#endif // !defined(DQUAD_H__INCLUDED)
