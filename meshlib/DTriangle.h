/////////////////////////////////////////////////////////////////////////////
// DTriangle.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DTRIANGLE_H__INCLUDED
#define DTRIANGLE_H__INCLUDED

#include "DPoint.h"

class ITriangle
{
public:
	ITriangle(int _a, int _b, int _c) : a(_a), b(_b), c(_c) { assert( a >= 0 && b >= 0 && c >= 0 ); }
	ITriangle() : a(-1), b(-1), c(-1) {}
public:
	int a, b, c;
};

/**
 * This class implements a triangle (three points) in a two-dimensional space
 *  plus some basic operations.
 */
class DTriangle2d
{
public:
	/// Standard constructor
	DTriangle2d(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c) : pt_a(a), pt_b(b), pt_c(c) {}
public:
	/// Returns the area of the triangle
	double area() const { 
		return 0.5 * ((pt_b.x - pt_a.x)*(pt_c.y - pt_a.y) - ((pt_b.y - pt_a.y)*(pt_c.x - pt_a.x))); }
	/// Returns the area of the triangle
	static double area(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c) { 
		return 0.5 * ((b.x - a.x)*(c.y - a.y) - ((b.y - a.y)*(c.x - a.x))); }
	/// Returns the area of the triangle
	static double area(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c) { 
		return 0.5 * ((b.x - a.x)*(c.y - a.y) - ((b.y - a.y)*(c.x - a.x))); }
	/// Returns the area of the triangle
	static double det(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c) { 
		return (b.x - a.x)*(c.y - a.y) - ((b.y - a.y)*(c.x - a.x)); }
	/// Returns the area of the triangle
	double det() const { 
		return (pt_b.x - pt_a.x)*(pt_c.y - pt_a.y) - ((pt_b.y - pt_a.y)*(pt_c.x - pt_a.x)); }
	/// Returns the area of the triangle
	static double det(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c) { 
		return (b.x - a.x)*(c.y - a.y) - ((b.y - a.y)*(c.x - a.x)); }
	/// Whether triangle is valid (not inverted)
	static bool valid(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c) { return det(a,b,c) > 0.0; }

	/// Check if point is within the circumscribed circle for a triangle given by three points
	double inSphereCheck(const DPoint2d& pt) const;
	/// Check if point is within the circumscribed circle for a triangle given by three points
	static double inSphereCheck(const DPoint2d& pt, const DPoint2d& a, const DPoint2d& b, const DPoint2d& c);
	/// Check if point is within the circumscribed circle for a triangle given by three points
	static double inSphereCheck(const DMPoint2d& pt, const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c);

	/// Whether contains the given point
	static bool containsPoint(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c, 
		const DPoint2d& pt, double eps = 0.0);

	/// Returns the length of the radius of the inscribed circle (for a triangle given by three points)
	static double innerCircleRadius(const DPoint2d &a, const DPoint2d &b, const DPoint2d &c);
	/// Returns the length of the radius of the inscribed circle (for a triangle given by three points)
	static double innerCircleRadius(const DMPoint2d &a, const DMPoint2d &b, const DMPoint2d &c);
	/// Returns the length of the radius of the circumscribed circle (for a triangle given by three points)
	static double outerCircleRadius(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c);
	/// Returns the length of the radius of the circumscribed circle (for a triangle given by three points)
	static double outerCircleRadius(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c);
	/// Returns the length of the radius of the circumscribed circle (for a triangle given by three points)
	static double outerCircleRadius2(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c);
	/// Returns the length of the radius of the circumscribed circle (for a triangle given by three points)
	static double outerCircleRadius2(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c);

	/// Returns middle of the circumscribed circle (for a triangle given by three points)
	static const DPoint2d outerCircleCenter(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c);
	/// Returns middle of the circumscribed circle (for a triangle given by three points)
	static const DMPoint2d outerCircleCenter(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c);

	/// Returns the "alpha"-quality coefficient of the triangle
	static double alphaQuality(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c);
	/// Returns the "alpha"-quality coefficient of the triangle
	static double alphaQuality(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c);

	/// Distance (squared) from triangle to point
	static double distance2ToPoint(const DPoint2d& pt, const DPoint2d& a, const DPoint2d& b, const DPoint2d& c );
	/// Distance (squared) from triangle to point
	double distance2ToPoint(const DPoint2d& pt) const;
	/// Distance (squared) from triangle to point
	static double distance2ToPoint(const DPoint2d& pt, const DPoint2d& a, const DVector2d& ab, const DVector2d& ac);

	/// Distance (squared) from triangle to point
	static double distance2ToPoint(const DMPoint2d& pt, const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c );
	/// Distance (squared) from triangle to point
	static double distance2ToPoint(const DMPoint2d& pt, const DMPoint2d& a, const DMVector2d& ab, const DMVector2d& ac);

public:
	/// points
	DPoint2d pt_a, pt_b, pt_c; 
};

/**
 * This class implements a segment (two points) in a three-dimensional space
 *  plus some basic operations
 */
class DTriangle3d
{
public:
	/// Standard constructor
	DTriangle3d() {}
	DTriangle3d(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c) : pt_a(a), pt_b(b), pt_c(c) {}
public:
	/// Return normal vector
	DVector3d normalVector() const;
	/// Returns the det of the triangle
	double det() const { 
		return (pt_b-pt_a).crossProduct(pt_c-pt_a).length(); }
	/// Returns the area of the triangle
	double area() const { 
		return 0.5 * (pt_b-pt_a).crossProduct(pt_c-pt_a).length(); }
	/// Returns the area of the triangle
	static double area(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c) { 
		return 0.5 * (b-a).crossProduct(c-a).length(); }

	/// Orient 3d
	double orient3d(const DPoint3d& d) const;
	/// Orient 3d
	static double orient3d(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d);

	/// Distance (squared) from triangle to point
	static double distance2ToPoint(const DPoint3d& pt, const DPoint3d& a, const DPoint3d& b, const DPoint3d& c );
	/// Distance (squared) from triangle to point
	double distance2ToPoint(const DPoint3d& pt) const;
	/// Distance (squared) from triangle to point
	static double distance2ToPoint(const DPoint3d& pt, const DPoint3d& a, const DVector3d& ab, const DVector3d& ac);

	/// Whether the point lies within the triangle after planar projection
	static bool containsProjectedPoint(const DPoint3d& pt, const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, double eps = 0.0 );
	/// Whether the point lies within the triangle after planar projection
	static bool containsProjectedPoint(const DPoint3d& pt, const DPoint3d& a, const DVector3d& ab, const DVector3d& ac, double eps = 0.0 );

	/// Distance (squared) from triangle to point - but only if projected onto triangle, the point is within
	static double distance2ToPointOrthogonal(const DPoint3d& pt, const DPoint3d& a, const DPoint3d& b, const DPoint3d& c);

	/// Check if segment |p1p2| is crossing triangle |p3p4p5|
	static bool crossSegment(const DPoint3d& sp1, const DPoint3d& sp2, 
					const DPoint3d& a, const DPoint3d& b, const DPoint3d& c);
	/// Check if segment |p1p2| is crossing triangle |p3p4p5|
	bool crossSegment(const DPoint3d& sp1, const DPoint3d& sp2) const;

	/// Returns the "alpha"-quality coefficient of the triangle
	static double alphaQuality(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c);

public:
	/// points
	DPoint3d pt_a, pt_b, pt_c; 
};

/**
 * This class implements a segment (two points) in a three-dimensional space
 *  plus some basic operations
 */
class DMTriangle3d
{
public:
	/// Standard constructor
	DMTriangle3d() {}
	DMTriangle3d(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c) : pt_a(a), pt_b(b), pt_c(c) {}
public:
	/// Return normal vector
	DMVector3d normalVector() const;
	/// Returns the det of the triangle
	double det() const { 
		return (pt_b-pt_a).crossProduct(pt_c-pt_a).length(); }
	/// Returns the area of the triangle
	double area() const { 
		return 0.5 * (pt_b-pt_a).crossProduct(pt_c-pt_a).length(); }
	/// Returns the area of the triangle
	static double area(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c) { 
		return 0.5 * (b-a).crossProduct(c-a).length(); }

	/// Orient 3d
	double orient3d(const DMPoint3d& d) const;
	/// Orient 3d
	static double orient3d(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d);

	/// Distance (squared) from triangle to point
	double distance2ToPoint(const DMPoint3d& pt) const;
	/// Distance (squared) from triangle to point
	static double distance2ToPoint(const DMPoint3d& pt, const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c);

	/// Check if segment |p1p2| is crossing triangle |p3p4p5|
	bool crossSegment(const DMPoint3d& sp1, const DMPoint3d& sp2) const;
	/// Check if segment |p1p2| is crossing triangle |p3p4p5|
	static bool crossSegment(const DMPoint3d& sp1, const DMPoint3d& sp2, 
					const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c);

	/// Returns the "alpha"-quality coefficient of the triangle
	static double alphaQuality(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c);
public:
	/// points
	DMPoint3d pt_a, pt_b, pt_c; 
};

#endif // !defined(DTRIANGLE_H__INCLUDED)
