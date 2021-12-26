/////////////////////////////////////////////////////////////////////////////
// SurfaceMonoAnalytic.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEMONOANALYTIC_H__INCLUDED)
#define SURFACEMONOANALYTIC_H__INCLUDED

#include "SurfaceParametric.h"
#include "DEquation.h"

/**
 * This class implements a surface defined by three equations x(t), y(t), z(t)
 *	with several required methods.
 */
class SurfaceMonoAnalytic : public SurfaceParametric  
{
public:
	enum FTYPE { FXY, FXZ, FYZ };
public:
	/// Standard constructor
	SurfaceMonoAnalytic(const char *formula = "", FTYPE ftype = FXY, DEquationConstTable *ctable = nullptr);
	/// Standard destructor
	virtual ~SurfaceMonoAnalytic() { clear(); }
public:
	virtual string getSimpleDescription() const override { return "mono-analytic"; }
	/// Returns coordinates defined by parameters x(s,t), y(s,t), z(s,t)
	virtual const DPoint3d getPoint(const DPoint2d& param) const;
	/// Returns the (specified) derivative vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const;
	/// Returns the parameters of the surface for the given point with starting point
	virtual const DPoint2d getParametersNear(const DPoint3d& point, const DPoint2d& near_point) const;
	/// Returns the parameters of the surface for the given point
	virtual const DPoint2d getParameters(const DPoint3d& point) const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return SURFACE_MONOANALYTIC; }
public:
	/// Stores the parameters of this surface into the stream
	friend ostream& operator<<(ostream& os, const SurfaceMonoAnalytic* surf);
	/// Restores the parameters of this surface from the stream
	friend istream& operator>>(istream& is, SurfaceMonoAnalytic* surf);
	/// Sets new formulas defining this parametrical surface
	void setData(const char *formula_x = "", FTYPE ftype = FXY, DEquationConstTable *ctable = nullptr);
protected:
	/// Clears parameters of this surface (makes it invalid)
	void clear();
protected:
	/// Type of formula
	FTYPE m_ftype;
	/// Formula-object for f(u,v)
	DEquation m_eq;
	/// Array of formula-objects for derivatives [xs,ys,zs, xt,yt,zt, xss,yss,zss, xst,yst,zst, xtt, ytt, ztt]
	DEquation* m_eq_deriv[15];
	/// String formula for f(u,v)
	string m_str;
};

#endif // !defined(SURFACEMONOANALYTIC_H__INCLUDED)
