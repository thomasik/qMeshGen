/////////////////////////////////////////////////////////////////////////////
// Curve3dSegment.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVE3DSEGMENT_H__INCLUDED)
#define CURVE3DSEGMENT_H__INCLUDED

#include "Curve3dParametric.h"
#include "DPoint.h"

class Metric3dContext;

/**
 * This class implements an constructional straight segment
 *  with required methods.
 */
class Curve3dSegment : public Curve3dParametric
{
public:
	/// Standard constructor
	Curve3dSegment(const DPoint3d& pt0, const DPoint3d& pt1 ) 
		: Curve3dParametric(0.0, 1.0 ), m_pt(pt0), m_vt(pt1-pt0) {}
	Curve3dSegment(const DPoint3d& pt0, const DVector3d& vt )
		: Curve3dParametric(0.0, 1.0 ), m_pt(pt0), m_vt(vt) {}
public:
	/// Returns basic description of this curve type
	virtual string getSimpleDescription() const { return "segment3d"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return CURVE_3D_SEGMENT; }
	/// Returns a point 2D for the given parameter t
	virtual DPoint3d getPoint(double t) const { return m_pt + m_vt * t; }
	/// Returns a parameter t for the given point 2D (t > ts), in range
	virtual double getParameterInRange(const DPoint3d& pt, double ts, double t_min, double t_max ) const;
	/// Returns a parameter t for the given point 2D (t > ts)
	virtual double getParameter(const DPoint3d& pt, double ts ) const;
	/// Returns derivatives for the given parameter (dx/dt, dy/dt)
	virtual DVector3d getDerivative(double) const { return m_vt; }
	/// Returns second derivatives for the given parameter (d2x/dt2, d2y/dt2)
	virtual DVector3d getSecondDerivative(double) const { return DVector3d(0.0, 0.0, 0.0); }
	/// Returns third derivatives for the given parameter (d3x/dt3, d3y/dt3)
	virtual DVector3d getThirdDerivative(double) const { return DVector3d(0.0, 0.0, 0.0); }
	/// Returns length of the segment of the curve (marked via t0 - t1, with metric if needed)
	virtual double getLength(Metric3dContext& mc, double t0, double t1, bool local_metric = true) const;
	/// Returns polyline (array of points and count) approximating a segment (t1-t2) of this curve
	virtual void getPolyLineInRange(double t0, double t1, DataVector<DPoint3d>& polyline) const;
	/// Returns polyline (array of points and count) approximating a segment (t1-t2) of this curve
	virtual void getPolyLineInRange(double t0, double t1, DataVector<double>& polyline) const;
	/// Returns rectangle containing this curve (if possible)
	virtual DBox getBoundingBox(double t0, double t1) const;
	/// Returns curvature of the 3D curve for the given parameter ( ~ 1/R)
	virtual double getCurvature(double t, double* cdt_len = nullptr) const;
	/// Checks wheterh this segment is valid (i.e. is properly initialized)
	virtual bool isValid() const { return m_vt.length2() > 0.0; }
protected:
	/// Beginning of the segment
	DPoint3d m_pt;	
	/// Direction vector
	DVector3d m_vt;
};

#endif // !defined(CURVE3DSEGMENT_H__INCLUDED)
