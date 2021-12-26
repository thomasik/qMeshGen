/////////////////////////////////////////////////////////////////////////////
// Curve3dLinearQuadric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVE3DLINEARQUADRIC_H__INCLUDED)
#define CURVE3DLINEARQUADRIC_H__INCLUDED

#include "DLinearQuadric.h"
#include "Curve3dParametric.h"

/**
 * This class implements a 3D linear quadric
 *	with several required methods.
 */
class Curve3dLinearQuadric : public Curve3dParametric
{
public:
	/// Standard constructor
	Curve3dLinearQuadric(const DLinearQuadric3d& lq) : m_lq(lq) { }
public:
	virtual string getSimpleDescription() const override { return "linear-quadric3d"; }
	/// for storing in XML file as a name of the surface-type
	virtual string typenameXML() const override { return "linear-quadric3d"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const override;
	/// Returns the object-specific type 
	virtual ElementType getType() const override { return CURVE_3D_LQUADRIC; }
	/// Returns the point of this curve for the given parameter
	virtual DPoint3d getPoint(double t) const override ;
	/// Returns a parameter t for the given point 2D (t > ts), in range
	virtual double getParameterInRange( const DPoint3d& pt, double ts, double t_min, double t_max ) const override;
	/// Returns the derivatives (dx/dt, dy/dt) of the curve for the given parameter
	virtual DVector3d getDerivative(double t) const override;
	/// Returns the second derivatives (d2x/dt2, d2y/dt2) of the curve for the given parameter
	virtual DVector3d getSecondDerivative(double t) const override;
	/// Returns the third derivatives (d3x/dt3, d3y/dt3) of the curve for the given parameter
	virtual DVector3d getThirdDerivative(double t) const override;
	/// Checks wheterh this curve is valid (i.e. is properly initialized)
	virtual bool isValid() const { return true; }
	/// Stores the parameters of this curve into the stream
	friend ostream& operator<<(ostream& os, const Curve3dLinearQuadric* lq);
	/// Restores the parameters of this curve from the stream
	friend istream& operator>>(istream& is, Curve3dLinearQuadric* lq);
private:
	/// linear quadric data
	DLinearQuadric3d m_lq;
};

#endif // !defined(CURVE3DLINEARQUADRIC_H__INCLUDED)
