/////////////////////////////////////////////////////////////////////////////
// Curve2dAnalytic.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVEANALYTIC_H__INCLUDED)
#define CURVEANALYTIC_H__INCLUDED

#include "DPoint.h"
#include "Curve2dParametric.h"
#include "DEquation.h"

/**
 * This class implements a parametrically-defined curve in 2D
 *	with several required methods.
 */
class Curve2dAnalytic : public Curve2dParametric
{
public:
	/// Standard constructor
	Curve2dAnalytic(string x_str = "", string y_str = "", DEquationConstTable *ctable = nullptr);
	/// Standard destructor
	virtual ~Curve2dAnalytic();
public:
	virtual string getSimpleDescription() const override { return "analytic"; }
	/// for storing in XML file as a name of the surface-type
	virtual string typenameXML() const { return "analytic"; }
	/// Returns the object-specific type
	virtual ElementType getType() const { return CURVE_ANALYTIC; }
	/// Returns the point of this curve for the given parameter
	virtual DPoint2d getPoint(double t) const;
	/// Returns the derivatives (dx/dt, dy/dt) of the curve for the given parameter
	virtual DVector2d getDerivative(double t) const;
	/// Returns the second derivatives (d2x/dt2, d2y/dt2) of the curve for the given parameter
	virtual DVector2d getSecondDerivative(double t) const;
	/// Returns the third derivatives (d3x/dt3, d3y/dt3) of the curve for the given parameter
	virtual DVector2d getThirdDerivative(double t) const;
	/// Checks wheterh this curve is valid (i.e. is properly initialized)
	virtual bool isValid() const { return m_valid; }
	/// Stores the parameters of this curve into the stream
	friend ostream& operator<<(ostream& os, const Curve2dAnalytic* figure);
	/// Restores the parameters of this curve from the stream
	friend istream& operator>>(istream& is, Curve2dAnalytic* figure);
private:
	/// Equation object for x(t) formula for this curve
	DEquation m_eq_x;
	/// Equation object for y(t) formula for this curve
	DEquation m_eq_y;
	/// Equation object for x'(t) formula for this curve
	DEquation* m_eq_xt;
	/// Equation object for y'(t) formula for this curve
	DEquation* m_eq_yt;
	/// Equation object for x''(t) formula for this curve
	DEquation* m_eq_xtt;
	/// Equation object for y''(t) formula for this curve
	DEquation* m_eq_ytt;
	/// Equation object for x'''(t) formula for this curve
	DEquation* m_eq_xttt;
	/// Equation object for y'''(t) formula for this curve
	DEquation* m_eq_yttt;
	/// Textual representaion of x(t) formula
	string m_str_x;
	/// Textual representaion of y(t) formula
	string m_str_y;
	/// Validity flag
	bool m_valid;
};

#endif // !defined(CURVEANALYTIC_H__INCLUDED)
