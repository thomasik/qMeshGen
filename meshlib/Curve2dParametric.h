/////////////////////////////////////////////////////////////////////////////
// Curve2dParametric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVEPARAMERIC_H__INCLUDED)
#define CURVEPARAMERIC_H__INCLUDED

#include <memory>

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "TagExtended.h"
#include "MeshData.h"

class SurfaceParametric;
class Metric2dContext;

/**
 * This class implements an abstract constructional curve
 *  and required methods.
 */
class Curve2dParametric : public TagExtended
{
public:
	/// Standard constructor
	Curve2dParametric(const double& t_min = -LARGE_NUMBER, const double& t_max = LARGE_NUMBER )
		: m_param_min( t_min ), m_param_max( t_max ), m_fixed( false ) {}
	/// Virtual destructor
	virtual ~Curve2dParametric() {}
public:
	virtual string getSimpleDescription() const { return "unknown"; }
	/// for storing in XML file as a name of the surface-type
	virtual string typenameXML() const { return "unknown"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the object-specific type 
	virtual ElementType getType() const { return CURVE_UNKNOWN; }
	/// Returns a point 2D for the given parameter t
	virtual DPoint2d getPoint(double t) const = 0;
	/// Returns a parameter t for the given point 2D (t > ts), in range
	virtual double getParameterInRange( const DPoint2d& pt, double ts, double t_min, double t_max ) const;
	/// Returns a parameter t for the given point 2D (t > ts)
	virtual double getParameter( const DPoint2d& pt, double ts ) const;
	/// Returns derivatives for the given parameter (dx/dt, dy/dt)
	virtual DVector2d getDerivative(double t) const = 0;
	/// Returns second derivatives for the given parameter (d2x/dt2, d2y/dt2)
	virtual DVector2d getSecondDerivative(double t) const = 0;
	/// Returns third derivatives for the given parameter (d3x/dt3, d3y/dt3)
	virtual DVector2d getThirdDerivative(double t) const = 0;
	/// Returns length of the segment of the curve (marked via t0 - t1, in parametric or metric space)
	virtual double getLength(Metric2dContext& mc, double t0, double t1, bool local_metric = true) const;
	/// Returns length of the segment of the curve (marked via t0 - t1, in real space)
	virtual double getLengthOnSurface(double t0, double t1, std::shared_ptr<const SurfaceParametric> surface) const;
	/// Returns length of the segment of the curve (marked via t0 - t1, with metric if needed)
	virtual double checkAndGetLength(Metric2dContext& mc, double t0, double& t1, double max_len = 1.0, bool local_metric = false) const;
	/// Returns polyline (array of points and count) approximating a segment (t_min-t_max) of this curve
	virtual void getPolyLine( DataVector<double>& polyline ) const { getPolyLineInRange( m_param_min, m_param_max, polyline); }
	/// Returns polyline (array of points and count) approximating a segment (t_min-t_max) of this curve
	virtual DPoint2d* getPolyLine(int & ct) const { return getPolyLineInRange( ct, m_param_min, m_param_max); }
	/// Returns polyline (array of points and count) approximating a segment (t1-t2) of this curve
	virtual DPoint2d* getPolyLineInRange(int & ct, double t1, double t2) const;
	/// Returns polyline (array of points and count) approximating a segment (t1-t2) of this curve
	virtual void getPolyLineInRange(double t0, double t1, DataVector<double>& polyline) const;
	/// Returns rectangle containing this curve (if possible)
	virtual DRect getBoundingRect(double t0, double t1) const;
	/// Returns curvature of the 2D curve for the given parameter ( ~ 1/R)
	virtual double getPlanarCurvature(double t) const;
	/// Returns curvature of the 3D curve for the given parameter ( ~ 1/R)
	virtual double getNonPlanarCurvature(std::shared_ptr<const SurfaceParametric> surface, double t, double* cdt_len = nullptr) const;
	/// min/max
	void setMinParam(const double & t_min) { m_param_min = t_min; }
	void setMaxParam(const double & t_max) { m_param_max = t_max; }
	void setMinMaxParam(const double & t_min, const double & t_max) { m_param_min = t_min; m_param_max = t_max; }
	const double& getMinParam() const { return m_param_min; }
	const double& getMaxParam() const { return m_param_max; }
	const double& addMinParam(const double & t_min) { if( t_min < m_param_min ) m_param_min = t_min; return m_param_min; }
	const double& addMaxParam(const double & t_max) { if( t_max > m_param_max ) m_param_max = t_max; return m_param_max; }
	void addMinMaxParam(const double & t) { if( t > m_param_max ) m_param_max = t; else if( t < m_param_min) m_param_min = t; }
	/// set fixed - not mutable
	void setFixed( bool fixed = true ) { m_fixed = fixed; }
	/// whether is fixed
	bool isFixed() const { return m_fixed; }
protected:
	double m_param_min;
	double m_param_max;
protected:
	/// Whether this surface is "mutable" or "final-and-permanent"
	bool m_fixed;
	/// Step used for the approximation of this curve by polyline
	static double m_probe_step;	
};

typedef std::shared_ptr<Curve2dParametric> Curve2dPtr;
typedef std::shared_ptr<const Curve2dParametric> Curve2dConstPtr;

#endif // !defined(CURVEPARAMERIC_H__INCLUDED)
