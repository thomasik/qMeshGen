/////////////////////////////////////////////////////////////////////////////
// DSegment.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DSEGMENT_H__INCLUDED
#define DSEGMENT_H__INCLUDED

#include "DPoint.h"


class ISegment
{
public:
	ISegment(int _a, int _b) : a(_a), b(_b) { assert( a >= 0 && b >= 0 ); }
	ISegment() : a(-1), b(-1) {}
public:
	int a, b;
};

/**
 * This class implements a segment (two points) in a two-dimensional space
 *  plus some basic operations.
 */
class DSegment2d
{
public:
	/// Standard constructor
	DSegment2d(const DPoint2d& a, const DPoint2d& b) : pt_a(a), pt_b(b) {}
public:
	/// Returns the length of the segment
	double length() const { return pt_a.distance(pt_b); }
	/// Returns the squared length of the segment
	double length2() const { return pt_a.distance2(pt_b); }

	/// Returns distance (squared) from segment to point
	static double distanceToPoint2(const DPoint2d& p1, const DPoint2d& p2, const DPoint2d& pt);
	/// Returns distance from segment to point
	static double distanceToPoint(const DPoint2d& p1, const DPoint2d& p2, const DPoint2d& pt);

	/// Check, whether the two segments are crossing ...
	bool crossingSegment(const DSegment2d& other) const;
public:
	/// first point
	DPoint2d pt_a; 
	/// second point
	DPoint2d pt_b; 
};

/**
 * This class implements a segment (two points) in a three-dimensional space
 *  plus some basic operations
 */
class DSegment3d
{
public:
	/// Standard constructor
	DSegment3d(const DPoint3d& a, const DPoint3d& b) : pt_a(a), pt_b(b) {}
public:
	/// Returns the length of the segment
	double length() const { return pt_a.distance(pt_b); }
	/// Returns the squared length of the segment
	double length2() const { return pt_a.distance2(pt_b); }

	/// Returns squared distance between two segments
	double distance2ToSegment(const DPoint3d& spt_a, const DPoint3d& spt_b) const {
		return distance2ToSegment(pt_a, pt_b, spt_a, spt_b);
	}

	/// Returns squared distance between two segments
	static double distance2ToSegment(const DPoint3d& spt0_a, const DPoint3d& spt0_b, 
		const DPoint3d& spt1_a, const DPoint3d& spt1_b);

public:
	/// first point
	DPoint3d pt_a; 
	/// second point
	DPoint3d pt_b; 
};

class DMSegment3d
{
public:
	/// Standard constructor
	DMSegment3d(const DMPoint3d& a, const DMPoint3d& b) : pt_a(a), pt_b(b) {}
public:
	/// Returns the length of the segment
	double length() const { return pt_a.distance(pt_b); }
	/// Returns the squared length of the segment
	double length2() const { return pt_a.distance2(pt_b); }

	/// Returns squared distance between two segments
	double distance2ToSegment(const DMPoint3d& spt_a, const DMPoint3d& spt_b) const {
		return distance2ToSegment(pt_a, pt_b, spt_a, spt_b);
	}

	/// Returns squared distance between two segments
	static double distance2ToSegment(const DMPoint3d& spt0_a, const DMPoint3d& spt0_b, 
		const DMPoint3d& spt1_a, const DMPoint3d& spt1_b);

public:
	/// first point
	DMPoint3d pt_a; 
	/// second point
	DMPoint3d pt_b; 
};

#endif // !defined(DSEGMENT_H__INCLUDED)
