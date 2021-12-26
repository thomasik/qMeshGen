/////////////////////////////////////////////////////////////////////////////
// Curve2dSegment.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2005-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(CURVESEGMENT_H__INCLUDED)
#define CURVESEGMENT_H__INCLUDED

#include "Curve2dParametric.h"
#include "DPoint.h"
#include "Metric2dContext.h"

class SurfaceParametric;

/**
 * This class implements an constructional straight segment
 *  with required methods.
 */
class Curve2dSegment : public Curve2dParametric
{
public:
	/// Standard constructor
	Curve2dSegment(const DPoint2d& pt0, const DPoint2d& pt1) 
		: Curve2dParametric(0.0, 1.0), m_pt(pt0), m_vt(pt1-pt0) {}
	Curve2dSegment(const DPoint2d& pt0, const DVector2d& vt)
		: Curve2dParametric(0.0, 1.0), m_pt(pt0), m_vt(vt) {}
public:
	virtual string getSimpleDescription() const override { return "segment2d"; }
	/// for storing in XML file as a name of the curve-type
	virtual string typenameXML() const { return "segment"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return CURVE_SEGMENT; }
	/// Returns a point 2D for the given parameter t
	virtual DPoint2d getPoint(double t) const { return m_pt + m_vt * t; }
	/// Returns a parameter t for the given point 2D (t > ts), in range
	virtual double getParameterInRange( const DPoint2d& pt, double ts, double t_min, double t_max ) const;
	/// Returns a parameter t for the given point 2D (t > ts)
	virtual double getParameter( const DPoint2d& pt, double ts ) const;
	/// Returns derivatives for the given parameter (dx/dt, dy/dt)
	virtual DVector2d getDerivative(double) const { return m_vt; }
	/// Returns second derivatives for the given parameter (d2x/dt2, d2y/dt2)
	virtual DVector2d getSecondDerivative(double) const { return DVector2d(0.0, 0.0); }
	/// Returns third derivatives for the given parameter (d3x/dt3, d3y/dt3)
	virtual DVector2d getThirdDerivative(double) const { return DVector2d(0.0, 0.0); }
	/// Returns length of the segment of the curve (marked via t0 - t1, with metric if needed)
	virtual double getLength(Metric2dContext& mc, double t0, double t1, bool local_metric = true) const;
	/// Returns polyline (array of points and count) approximating a segment (t1-t2) of this curve
	virtual DPoint2d* getPolyLineInRange(int & ct, double t1, double t2) const;
	/// Returns polyline (array of points and count) approximating a segment (t1-t2) of this curve
	virtual void getPolyLineInRange(double t1, double t2, DataVector<double>& polyline) const;
	/// Returns rectangle containing this curve (if possible)
	virtual DRect getBoundingRect(double t0, double t1) const;
	/// Returns curvature of the 2D curve for the given parameter ( ~ 1/R)
	virtual double getPlanarCurvature(double) const { return 0.0; }
	/// Returns curvature of the 3D curve for the given parameter ( ~ 1/R)
	virtual double getNonPlanarCurvature(std::shared_ptr<const SurfaceParametric> surface, 
		double t, double* cdt_len = nullptr) const;
	/// Checks wheterh this segment is valid (i.e. is properly initialized)
	virtual bool isValid() const { return m_vt.length2() > 0.0; }
protected:
	/// Beginning of the segment
	DPoint2d m_pt;	
	/// Direction vector
	DVector2d m_vt;
};

#endif // !defined(CURVESEGMENT_H__INCLUDED)
