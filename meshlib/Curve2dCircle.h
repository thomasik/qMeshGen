/////////////////////////////////////////////////////////////////////////////
// Curve2dCircle.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVECIRCLE_H__INCLUDED)
#define CURVECIRCLE_H__INCLUDED

#include "DPoint.h"
#include "Curve2dParametric.h"
#include "DSphere.h"

/**
 * This class implements a 2D base circle
 *	with several required methods.
 */
class Curve2dCircle : public Curve2dParametric
{
public:
	/// Standard constructor
	Curve2dCircle(double x=0.0, double y=0.0, double r=0.0) : Curve2dParametric(0.0, 1.0), m_middle(x, y), m_radius(r) {}
	/// Standard constructor
	Curve2dCircle(const DPoint2d& pt, double r) : Curve2dParametric(0.0, 1.0), m_middle(pt), m_radius(r) {}
	/// Standard constructor
	Curve2dCircle(const DCircle& circle) : Curve2dParametric(0.0, 1.0), m_middle(circle.m_center), m_radius(circle.m_radius) {}
public:
	virtual string getSimpleDescription() const override { return "circle"; }
	/// for storing in XML file as a name of the surface-type
	virtual string typenameXML() const { return "circle"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the object-specific type 
	virtual ElementType getType() const { return CURVE_CIRCLE; }
	/// Returns the point of this curve for the given parameter
	virtual DPoint2d getPoint(double t) const;
	/// Returns a parameter t for the given point 2D (t > ts), in range
	virtual double getParameterInRange( const DPoint2d& pt, double ts, double t_min, double t_max ) const;
	/// Returns the parameter for this curve associated with the given point (numerical approximation)
	virtual double getParameter( const DPoint2d& pt, double ts ) const;
	/// Returns the parameter for this curve associated with the given point (numerical approximation)
	static double getCircleParameter( const DPoint2d& mid, double r, const DPoint2d& pt, double ts );
	/// Returns the derivatives (dx/dt, dy/dt) of the curve for the given parameter
	virtual DVector2d getDerivative(double t) const;
	/// Returns the second derivatives (d2x/dt2, d2y/dt2) of the curve for the given parameter
	virtual DVector2d getSecondDerivative(double t) const;
	/// Returns the third derivatives (d3x/dt3, d3y/dt3) of the curve for the given parameter
	virtual DVector2d getThirdDerivative(double t) const;
	/// Returns curvature of the curve for the given parameter ( ~ 1/R)
	virtual double getPlanarCurvature(double) const { return 1.0 / m_radius; }
	/// Checks wheterh this curve is valid (i.e. is properly initialized)
	virtual bool isValid() const { return m_radius > 0; }
	/// Stores the parameters of this curve into the stream
	friend ostream& operator<<(ostream& os, const Curve2dCircle* circle);
	/// Restores the parameters of this curve from the stream
	friend istream& operator>>(istream& is, Curve2dCircle* circle);
private:
	/// Middle point
	DPoint2d	m_middle;
	/// Radius
	double		m_radius;
};

#endif // !defined(CURVECIRCLE_H__INCLUDED)
