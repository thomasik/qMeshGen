/////////////////////////////////////////////////////////////////////////////
// SurfaceCorrected.h 
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACECORRECTED_H__INCLUDED)
#define SURFACECORRECTED_H__INCLUDED

#include "DPoint.h"
#include "SurfaceParametric.h"

/**
 *  Surface patch with additional (vector-)corrections
 */
class SurfaceCorrected : public SurfaceParametric
{
public:
	/// Standard constructor
	SurfaceCorrected(std::shared_ptr<const SurfaceParametric> base_surface);
private:
	SurfaceCorrected(const SurfaceCorrected& ) = delete; // not available
	SurfaceCorrected& operator=(const SurfaceCorrected& ) = delete; // not available
public:
	virtual string getSimpleDescription() const override { return string("corrected ")+m_base_surface->getSimpleDescription(); }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the point of the surface for the given parameters
	virtual const DPoint3d getPoint(const DPoint2d& param) const;
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVector(const DPoint2d& param) const;
	/// Returns the derivative (specified) vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const; 
	/// Returns the parameters of the surface for the given point (numerical approximation)
	virtual const DPoint2d getParameters(const DPoint3d& point) const;
	/// Returns the parameters of the surface for the given point with starting point (numerical approximation)
	virtual const DPoint2d getParametersNear(const DPoint3d& point, const DPoint2d& near_point) const;
	/// Returns the object-specific type 
	virtual ElementType getType() const { return SURFACE_CORRECTED; }
public:
	/// Insert an additional correction vector
	void insertCorrectionVector(const DPoint2d& center, double radius, const DVector3d& translation);
	/// Returns the base surface
	std::shared_ptr<const SurfaceParametric> getBaseSurface() const { return m_base_surface; }
protected:
	// Correction vector
	class CorrectionVector{
	public:
		CorrectionVector() {}
		CorrectionVector(const DPoint2d& _center, double _radius, const DMetric2d& dm, const DVector3d& _translation);
	public:
		DPoint2d center;
		double radius;
		DMetric2d dm;
		DVector3d dv;
	};
protected:
	/// Base surface
	std::shared_ptr<const SurfaceParametric> m_base_surface;
	/// Correction-vectors
	DataVector<CorrectionVector> m_corrections;
};

#endif // !defined(SURFACECORRECTED_H__INCLUDED)
