/////////////////////////////////////////////////////////////////////////////
// Curve3dSurfaceParametric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVE3DSURFACEPARAMETRIC_H__INCLUDED)
#define CURVE3DSURFACEPARAMETRIC_H__INCLUDED

#include "Curve2dParametric.h"
#include "Curve3dParametric.h"
#include "SurfaceParametric.h"

class Metric3dContext;

/**
 * This class implements an constructional straight segment
 *  with required methods.
 */
class Curve3dSurfaceParametric : public Curve3dParametric
{
public:
	/// Standard constructor
	Curve3dSurfaceParametric(SurfaceConstPtr surface, Curve2dConstPtr curve);
public:
	/// Returns basic description of this curve type
	virtual string getSimpleDescription() const override;
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return CURVE_3D_SURFACE_PARAMETRIC; }
	/// Returns a point 2D for the given parameter t
	virtual DPoint3d getPoint(double t) const;
	/// Returns a parameter t for the given point 2D (t > ts) in Range
	virtual double getParameterInRange( const DPoint3d& pt, double ts, double t_min, double t_max ) const;
	/// Returns derivatives for the given parameter (dx/dt, dy/dt)
	virtual DVector3d getDerivative(double) const;
	/// Returns second derivatives for the given parameter (d2x/dt2, d2y/dt2)
	virtual DVector3d getSecondDerivative(double) const;
	/// Returns third derivatives for the given parameter (d3x/dt3, d3y/dt3)
	virtual DVector3d getThirdDerivative(double) const;
	/// Checks wheterh this curve is valid (i.e. is properly initialized)
	virtual bool isValid() const;
protected:
	/// base surface
	SurfaceConstPtr m_surface;	
	/// planar curve
	Curve2dConstPtr m_curve;
};

#endif // !defined(CURVE3DSURFACEPARAMETRIC_H__INCLUDED)
