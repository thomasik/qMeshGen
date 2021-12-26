/////////////////////////////////////////////////////////////////////////////
// Curve2dLinearQuadric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVELINEARQUADRIC_H__INCLUDED)
#define CURVELINEARQUADRIC_H__INCLUDED

#include "DLinearQuadric.h"
#include "Curve2dParametric.h"

/**
 * This class implements a 2D linear quadric
 *	with several required methods.
 */
class Curve2dLinearQuadric : public Curve2dParametric
{
public:
	/// Standard constructor
	Curve2dLinearQuadric(const DLinearQuadric& lq) : m_lq(lq) { }
public:
	virtual string getSimpleDescription() const override { return "linear-quadric"; }
	/// for storing in XML file as a name of the surface-type
	virtual string typenameXML() const { return "linear-quadric"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the object-specific type 
	virtual ElementType getType() const { return CURVE_LQUADRIC; }
	/// Returns the point of this curve for the given parameter
	virtual DPoint2d getPoint(double t) const;
	/// Returns a parameter t for the given point 2D (t > ts), in range
	virtual double getParameterInRange( const DPoint2d& pt, double ts, double t_min, double t_max ) const;
	/// Returns the derivatives (dx/dt, dy/dt) of the curve for the given parameter
	virtual DVector2d getDerivative(double t) const;
	/// Returns the second derivatives (d2x/dt2, d2y/dt2) of the curve for the given parameter
	virtual DVector2d getSecondDerivative(double t) const;
	/// Returns the third derivatives (d3x/dt3, d3y/dt3) of the curve for the given parameter
	virtual DVector2d getThirdDerivative(double t) const;
	/// Checks wheterh this curve is valid (i.e. is properly initialized)
	virtual bool isValid() const { return true; }
	/// Stores the parameters of this curve into the stream
	friend ostream& operator<<(ostream& os, const Curve2dLinearQuadric* lq);
	/// Restores the parameters of this curve from the stream
	friend istream& operator>>(istream& is, Curve2dLinearQuadric* lq);
private:
	/// linear quadric data
	DLinearQuadric m_lq;
};

#endif // !defined(CURVELINEARQUADRIC_H__INCLUDED)
