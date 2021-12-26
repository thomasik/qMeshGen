/////////////////////////////////////////////////////////////////////////////
// Curve3dParametric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVE3DPARAMETRIC_H__INCLUDED)
#define CURVE3DPARAMETRIC_H__INCLUDED

#include <memory>

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "MeshData.h"
#include "DataMatrix.h"
#include "SurfaceParametric.h"

class Metric3dContext;

/**
 * This class implements an abstract constructional curve
 *  and required methods.
 */
class Curve3dParametric
{
public:
	/// Standard constructor
	Curve3dParametric( const double & t_min = -LARGE_NUMBER, const double & t_max = LARGE_NUMBER )
		: m_param_min(t_min), m_param_max(t_max), m_fixed( false ) {}
	/// Virtual destructor
	virtual ~Curve3dParametric();
public:
	/// Fit 2d curve on parametric surface
	static std::shared_ptr<Curve3dParametric> fitCurveOnSurface(Metric3dContext& mc, SurfaceConstPtr surface,
		const DataVector<DPoint2d> & points, DataVector<double> & params, double tolerance);
public:
	/// for storing in XML file as a name of the surface-type
	virtual string typenameXML() const { return "unknown"; }
	/// Returns the object-specific type 
	virtual ElementType getType() const { return CURVE_3D_UNKNOWN; }
	/// Returns basic description of this curve type
	virtual string getSimpleDescription() const { return "unknown"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// domain check
	virtual bool withinDomain( const DPoint3d& pt, double & t ) const;
	/// Returns a point 2D for the given parameter t
	virtual DPoint3d getPoint(double t) const = 0;
	/// Returns a parameter t for the given point 3D (t > ts), in range
	virtual double getParameterInRange( const DPoint3d& pt, double ts, double t_min, double t_max ) const;
	/// Returns a parameter t for the given point 3D (t > ts)
	virtual double getParameter( const DPoint3d& pt, double ts ) const;
	/// Returns derivatives for the given parameter (dx/dt, dy/dt)
	virtual DVector3d getDerivative(double t) const = 0;
	/// Returns second derivatives for the given parameter (d2x/dt2, d2y/dt2)
	virtual DVector3d getSecondDerivative(double t) const = 0;
	/// Returns third derivatives for the given parameter (d3x/dt3, d3y/dt3)
	virtual DVector3d getThirdDerivative(double t) const = 0;
	/// Returns length of the segment of the curve (marked via t0 - t1, in parametric or metric space)
	virtual double getLength(Metric3dContext& mc, double t0, double t1, bool local_metric = true) const;
	/// Returns length of the segment of the curve (marked via t0 - t1, with metric if needed)
	virtual double checkAndGetLength(Metric3dContext& mc, double t0, double& t1, double max_len = 1.0, bool local_metric = false) const;
	/// Returns polyline (array of points and count) approximating a segment (t_min-t_max) of this curve
	virtual void getPolyLine( DataVector<DPoint3d>& polyline ) const { getPolyLineInRange( m_param_min, m_param_max, polyline); }
	/// Returns polyline (array of points and count) approximating a segment (t1-t2) of this curve
	virtual void getPolyLineInRange(double t0, double t1, DataVector<DPoint3d>& polyline) const;
	/// Returns polyline (array of points and count) approximating a segment (t_min-t_max) of this curve
	virtual void getPolyLine(DataVector<double>& polyline) const { getPolyLineInRange( m_param_min, m_param_max, polyline); }
	/// Returns polyline (array of points and count) approximating a segment (t1-t2) of this curve
	virtual void getPolyLineInRange(double t0, double t1, DataVector<double>& polyline) const;
	/// Returns rectangle containing this curve (if possible)
	virtual DBox getBoundingBoxInRange(double t0, double t1) const;
	/// Returns rectangle containing this curve (if possible)
	virtual DBox getBoundingBox() const { return getBoundingBoxInRange( m_param_min, m_param_max ); }
	/// Returns curvature of the 3D curve for the given parameter ( ~ 1/R)
	virtual double getCurvature(double t, double* cdt_len = nullptr) const;
	/// min/max
	void setMinParam(const double & t_min) { m_param_min = t_min; }
	void setMaxParam(const double & t_max) { m_param_max = t_max; }
	void setMinMaxParam(const double & t_min, const double & t_max) { m_param_min = t_min; m_param_max = t_max; }
	const double& getMinParam() const { return m_param_min; }
	const double& getMaxParam() const { return m_param_max; }
	double getMiddleParam() const { return 0.5 * ( m_param_min + m_param_max ); }
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
	/// Whether this surface is "mutable" or "final-and-permanent"
	bool m_fixed;
protected:
	/// Step used for the approximation of this curve by polyline
	static double m_probe_step;	
};

typedef std::shared_ptr<Curve3dParametric> Curve3dPtr;
typedef std::shared_ptr<const Curve3dParametric> Curve3dConstPtr;

class Curve3dParametricSet
{
public:
	Curve3dParametricSet(Curve3dConstPtr curve = nullptr);
	Curve3dParametricSet(std::shared_ptr<const Curve3dParametricSet> cset, Curve3dConstPtr curve = nullptr);
public:
	int addCurve(Curve3dConstPtr curve);
	int count() const;
	inline int countInt() const { return count(); }
	const Curve3dConstPtr& getCurve(int i) const { return m_curves[i]; }
//	Curve3dConstPtr* selectCurveForPoint( double t ) const;
private:
	DataVector< Curve3dConstPtr > m_curves;
};

typedef std::shared_ptr<Curve3dParametricSet> Curve3dSetPtr;
typedef std::shared_ptr<const Curve3dParametricSet> Curve3dSetConstPtr;

#endif // !defined(CURVE3DPARAMETRIC_H__INCLUDED)
