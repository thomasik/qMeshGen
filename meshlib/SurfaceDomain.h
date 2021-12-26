/////////////////////////////////////////////////////////////////////////////
// SurfaceDomain.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEDOMAIN_H__INCLUDED)
#define SURFACEDOMAIN_H__INCLUDED

class DPoint2d;
class MeshViewSet;
class SurfaceParametric;

#include <iostream>
#include "DOrientedBox.h"

class SurfaceDomain
{
public:
	SurfaceDomain() {} // standard, empty constructor -- infinite domain 
	virtual ~SurfaceDomain() {}
public:
	// [0-1] inside, less than 0 (usually set as -1) - outside
	virtual double getInsideQuality( const DPoint3d& /* point */, const DPoint2d& /* param */) const {	return AQ_VALID_MAX; }
	virtual bool isInsideOBBox(const DPoint3d& pt) const;
	void setOBBox( const DOrientedBox& box ) { m_box = box; }
	const DOrientedBox& getOBBox() const { return m_box; }
	virtual DPoint2d getMiddleParam() const { return DPoint2d::zero; }
public:
	static double w_inside_ratio;
	static double w_approx_error;
public:
	virtual void draw(MeshViewSet* set = nullptr, const SurfaceParametric * surface = nullptr, int id = -2) const;
	virtual std::ostream& storeXML(std::ostream& os, const std::string& /* prefix */ = "") const { return os; }
protected:
	DOrientedBox m_box;
};

#endif // SURFACEDOMAIN_H__INCLUDED
