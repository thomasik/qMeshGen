/////////////////////////////////////////////////////////////////////////////
// SurfaceSphere.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACESPHERE_H__INCLUDED)
#define SURFACESPHERE_H__INCLUDED

#include "SurfaceParametric.h"
#include "DSphere.h"

/**
 * This class implements a sphere with several required methods.
 */
class SurfaceSphere : public SurfaceParametric  
{
public:
	/// Standard constructor
	SurfaceSphere(const DPoint3d& pnt, double radius);
	SurfaceSphere(const DSphere& sphere);
public:
	virtual string getSimpleDescription() const override { return "sphere"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const override;
	/// Returns the principle curvature for this sphere and parameters [s,t]
	virtual const SurfaceCurvature getCurvature(const DPoint2d& /* pt */, double & g_ratio) const override;
	/// Returns the point of the plane for the given parameters
	virtual const DPoint3d getPoint(const DPoint2d& param) const override;
	/// Returns the parameters of the surface for the given point (numerical approximation) and signed distance (according to normal vecor)
	virtual bool getParametersAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z ) const override;
	/// Returns the parameters of the plane for the given point
	virtual const DPoint2d getParameters(const DPoint3d& point) const override;
	/// Returns the parameters of the plane for the given point
	virtual const DPoint2d getParametersNear(const DPoint3d& point, const DPoint2d& /* near_point */) const override
		{ return getParameters(point); }
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVector(const DPoint2d& param) const override;
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVectorForPoint3d(const DPoint3d& pt) const override;
	/// Returns the derivative (specified) vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const override;
	/// Returns the object-specific type
	virtual ElementType getType() const override { return SURFACE_SPHERE; }
	/// invert the orientation of the surface (change diretion of normal vector) if possible
	virtual bool invertOrientation();
	std::shared_ptr<SurfaceParametric> adjustedForFaces( const DataVector< MeshFace* > & mfaces,
		const DataVector< MeshPoint3d* > & mpoints, MeshFace* central_face = nullptr  ) const override;
	// fix center and range (for parameters) to fix periodic surfaces
	virtual bool fixCenterAndRangeForPeriodic(const DPoint3d& center_point ) override;
	virtual bool withinParamRange( const DPoint2d& param ) const override;
public:
	const DPoint3d& getCenter() const { return m_sphere.m_center; }
	const double getRadius() const { return m_sphere.getRadius(); }
protected:
	DPoint2d shiftedParams(const DPoint2d& param) const;
	DPoint2d unshiftedParams(const DPoint2d& param) const;
protected:
	DSphere m_sphere;
protected:
	bool m_non_periodic;
	double m_param_shift_x;
};

#endif // !defined(SURFACESPHERE_H__INCLUDED)
