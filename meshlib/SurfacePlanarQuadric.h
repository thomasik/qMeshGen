/////////////////////////////////////////////////////////////////////////////
// SurfacePlanarQuadric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEPLANARQUADRIC_H__INCLUDED)
#define SURFACEPLANARQUADRIC_H__INCLUDED

#include "DPlanarQuadric.h"
#include "SurfaceParametric.h"

/**
 * This class implements a planar quadric with several required methods.
 */
class SurfacePlanarQuadric : public SurfaceParametric
{
public:
	/// Standard constructor
	SurfacePlanarQuadric(const DPlanarQuadric& pquadric) : m_pquadric(pquadric) { m_valid = true; }
public:
	virtual string getSimpleDescription() const override { return "planar-quadric"; }
	/// return oriented bbox for this surface (invalid)
	virtual DOrientedBox getOrientedBox() const;
	/// return oriented bbox for this surface initialized with the given set of points
	virtual DOrientedBox getOrientedBox( const DataVector<DPoint3d>& points ) const;
	/// return oriented bbox for this surface initialized with the given set of points
	virtual DOrientedBox getOrientedBoxOpt( const DataVector<DPoint3d>& points ) const;
	/// invert the orientation of the surface (change diretion of normal vector) if possible
	virtual bool invertOrientation();
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the point of the planar quadric for the given parameters
	virtual const DPoint3d getPoint(const DPoint2d& param) const;
	/// Returns the parameters of the planar quadric for the given point
	virtual const DPoint2d getParameters(const DPoint3d& point) const;
	/// Returns the parameters of the planar quadric for the given point
	virtual const DPoint2d getParametersNear(const DPoint3d& point, const DPoint2d& /* near_point */) const
		{ return getParameters(point); }
	/// Returns the derivative (specified) vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return SURFACE_PLANAR_QUADRIC; }
public:
	/// Stores the parameters of this planar quadric into the stream
	friend ostream& operator<<(ostream& os, const SurfacePlanarQuadric* surf);
	/// Restores the parameters of this planar quadric from the stream
	friend istream& operator>>(istream& is, SurfacePlanarQuadric* surf);
public:
	const DPoint3d& baseP() const { return m_pquadric.plane.p0; }
	const DVector3d& baseEu() const { return m_pquadric.plane.e0; }
	const DVector3d& baseEv() const { return m_pquadric.plane.e1; }
protected:
	/// planar quadric data
	DPlanarQuadric m_pquadric;
};

#endif // !defined(SURFACEPLANARQUADRIC_H__INCLUDED)
