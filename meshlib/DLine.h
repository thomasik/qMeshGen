/////////////////////////////////////////////////////////////////////////////
// DLine.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DLINE_H__INCLUDED
#define DLINE_H__INCLUDED

#include "DPoint.h"
#include "DVector.h"

/**
 * This class implements a line (point + vector) in a two-dimensional space
 *  plus some basic operations.
 */
class DLine2d
{
public:
	/// Standard constructor
	DLine2d(const DPoint2d& line_p = DPoint2d::zero, const DVector2d& line_v = DVector2d::v_ox);
	/// Standard constructor
	DLine2d(const DPoint2d& p1, const DPoint2d& p2);
public:
	/// set line attributes
	void setPtVn(const DPoint2d& pt, const DVector2d& vn);
	/// set line attributes
	void setPtVt(const DPoint2d& pt, const DVector2d& vt);
	/// Returns distance (squared) from line to point
	double distanceToPoint2(const DPoint2d& pt) const;
	/// Returns distance from line to point
	double distanceToPoint(const DPoint2d& pt) const;

	/// Returns distance (squared) from line to point
	static double distanceToPoint2(const DPoint2d& line_p, const DVector2d& line_v, const DPoint2d& pt);
	/// Returns distance from line to point
	static double distanceToPoint(const DPoint2d& line_p, const DVector2d& line_v, const DPoint2d& pt);

	/// Sets this point as the crossing point of lines |p1p2| and |p3p4|
	static const DPoint2d crossPoint(const DPoint2d& p1, const DPoint2d& p2, const DPoint2d& p3, const DPoint2d& p4);
	/// Sets this point as the crossing point of lines |p1v1| and |p2v2|
	static const DPoint2d crossPoint(const DPoint2d& p1, const DVector2d& v1, const DPoint2d& p2, const DVector2d& v2);

	/// Returns point on the line for the given parameter
	DPoint2d getPoint(double t) const;
	/// Returns parameter for point on the line, nearest to pt
	double paramForPoint(const DPoint2d& pt) const;
	/// Calculates parameters for point on the line, and distance
	bool projectToLine(const DPoint2d& pt, double& t, double& z) const;

	/// Stores the parameters of this line into the stream
	friend ostream& operator<<(ostream& os, const DLine2d& line);
	/// Restores the parameters of this line from the stream
	friend istream& operator>>(istream& is, DLine2d& line);
public:
	/// point
	DPoint2d m_pt; 
	/// vector (normalized)
	DVector2d m_vt; 
	// normal vector (normalized)
	DVector2d m_vn;
};

/**
 * This class implements a segment (two points) in a three-dimensional space
 *  plus some basic operations
 */
class DLine3d
{
public:
	/// Standard constructor
	DLine3d(const DPoint3d& line_p = DPoint3d::zero, const DVector3d& line_v = DVector3d::v_ox );
	/// Standard constructor
	DLine3d(const DPoint3d& p1, const DPoint3d& p2);
public:
	/// set line attributes
	void setPtVt(const DPoint3d& pt, const DVector3d& vt);

	/// Returns distance (squared) from line to point
	double distanceToPoint2(const DPoint3d& pt) const;
	/// Returns distance from line to point
	double distanceToPoint(const DPoint3d& pt) const;

	/// Returns parameter for point on the line, nearest to pt
	double paramForPoint(const DPoint3d& pt) const;

	/// Returns point on the line for the given parameter
	DPoint3d getPoint(double t) const;
	/// Returns distance (squared) from line to point
	static double distanceToPoint2(const DPoint3d& line_p, const DVector3d& line_v, const DPoint3d& pt, bool line_v_normalized);
	/// Returns distance from line to point
	static double distanceToPoint(const DPoint3d& line_p, const DVector3d& line_v, const DPoint3d& pt, bool line_v_normalized);

	/// Stores the parameters of this line into the stream
	friend ostream& operator<<(ostream& os, const DLine3d& line);
	/// Restores the parameters of this line from the stream
	friend istream& operator>>(istream& is, DLine3d& line);

public:
	/// first point
	DPoint3d m_pt; 
	/// second point (normalized)
	DVector3d m_vt; 
};

#endif // !defined(DLINE_H__INCLUDED)
