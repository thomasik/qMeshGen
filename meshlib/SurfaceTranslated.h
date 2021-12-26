/////////////////////////////////////////////////////////////////////////////
// SurfaceTranslated.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2007-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACETRANSLATED_H__INCLUDED)
#define SURFACETRANSLATED_H__INCLUDED

#include "DPoint.h"
#include "DMetric2d.h"
#include "SurfaceParametric.h"

/**
 * Some parametric surface translated by given vector
 */
class SurfaceTranslated : public SurfaceParametric
{
public:
	/// Standard constructor
	SurfaceTranslated(std::shared_ptr<const SurfaceParametric> surface, const DVector3d& dv)
			: m_surface(surface), m_dv(dv) { assert(surface); }
public:
	virtual string getSimpleDescription() const override { return string("translated") + m_surface->getSimpleDescription(); }
	/// Returns the point of the surface for the given parameters
	virtual const DPoint3d getPoint(const DPoint2d& param) const { 
		return m_surface->getPoint(param) + m_dv; }
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVector(const DPoint2d& param) const { 
		return m_surface->getNormalVector(param); }
	/// Returns the derivative (specified) vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const { 
		return m_surface->getDerivative(deriv, param); }
	/// Returns the object-specific type 
	virtual ElementType getType() const { return SURFACE_TRANSLATED; }
protected:
	/// Base surface
	std::shared_ptr<const SurfaceParametric> m_surface;
	/// Translation vector
	DVector3d  m_dv;
};

#endif // !defined(SURFACETRANSLATED_H__INCLUDED)
