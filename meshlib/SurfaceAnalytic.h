/////////////////////////////////////////////////////////////////////////////
// SurfaceAnalytic.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEANALYTIC_H__INCLUDED)
#define SURFACEANALYTIC_H__INCLUDED

#include "SurfaceParametric.h"
#include "DEquation.h"

/**
 * This class implements a surface defined by three equations x(t), y(t), z(t)
 *	with several required methods.
 */
class SurfaceAnalytic : public SurfaceParametric  
{
public:
	/// Standard constructor
	SurfaceAnalytic(const char *formula_x = "", const char *formula_y = "", 
		const char *formula_z = "", DEquationConstTable *ctable = nullptr);
	/// Standard destructor
	virtual ~SurfaceAnalytic() { clear(); }
public:
	virtual string getSimpleDescription() const override { return "analytic"; }
	/// Returns coordinates defined by parameters x(s,t), y(s,t), z(s,t)
	virtual const DPoint3d getPoint(const DPoint2d& param) const;
	/// Returns the (specified) derivative vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return SURFACE_ANALYTIC; }
public:
	/// Stores the parameters of this surface into the stream
	friend ostream& operator<<(ostream& os, const SurfaceAnalytic* surf);
	/// Restores the parameters of this surface from the stream
	friend istream& operator>>(istream& is, SurfaceAnalytic* surf);
	/// Sets new formulas defining this parametrical surface
	void setData(const char *formula_x = "", const char *formula_y = "", 
		const char *formula_z = "", DEquationConstTable *ctable = nullptr);
protected:
	/// Clears parameters of this surface (makes it invalid)
	void clear();
protected:
	/// Formula-object for x(s,t)
	DEquation m_eq_x;
	/// Formula-object for y(s,t)
	DEquation m_eq_y;
	/// Formula-object for z(s,t)
	DEquation m_eq_z;
	/// Array of formula-objects for derivatives [xs,ys,zs, xt,yt,zt, xss,yss,zss, xst,yst,zst, xtt, ytt, ztt]
	DEquation* m_eq_deriv[15];
	/// Array of formula-objects for extra derivatives [xsss,ysss,zsss, xsst,ysst,zsst, xstt,ystt,zstt, xttt,yttt,zttt]
	DEquation* m_eq_deriv_ext[12];
	/// String formula for x(s,t)
	string m_str_x;
	/// String formula for y(s,t)
	string m_str_y;
	/// String formula for z(s,t)
	string m_str_z;
};

#endif // !defined(SURFACEANALYTIC_H__INCLUDED)
