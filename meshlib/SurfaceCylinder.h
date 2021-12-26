/////////////////////////////////////////////////////////////////////////////
// SurfaceCylinder.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACECYLINDER_H__INCLUDED)
#define SURFACECYLINDER_H__INCLUDED

#include "SurfaceParametric.h"
#include "DSphere.h"

/**
 * This class implements a cylinder with several required methods.
 */
class SurfaceCylinder : public SurfaceParametric  
{
public:
	/// Standard constructor
	SurfaceCylinder(const DPoint3d& pnt0, const DPoint3d& pnt1, double radius);
	SurfaceCylinder(const DPoint3d& pnt0, const DVector3d& vt,  double radius);
	SurfaceCylinder(const DCylinder & cylinder);	
	/// Standard destructor
	virtual ~SurfaceCylinder(){}
public:
	virtual string getSimpleDescription() const override { return "cylinder"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const override;
	/// Returns the principle curvature for this plane and parameters [s,t]
	virtual const SurfaceCurvature getCurvature(const DPoint2d& /* pt */, double & g_ratio) const override;
	/// Returns the point of the plane for the given parameters
	virtual const DPoint3d getPoint(const DPoint2d& param) const override;
	/// Returns the parameters of the plane for the given point
	virtual const DPoint2d getParameters(const DPoint3d& point) const override;
	/// Returns the parameters of the surface for the given point (numerical approximation) and signed distance (according to normal vecor)
	virtual bool getParametersAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z ) const override;
	/// Returns the parameters of the plane for the given point
	virtual const DPoint2d getParametersNear(const DPoint3d& point, const DPoint2d& /* near_point */) const override
		{ return getParameters(point); }
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVector(const DPoint2d& param) const override;
	/// Returns the derivative of normal vector for the given parameters
	virtual const DVector3d getNormalVectorDerivative(int deriv, const DPoint2d& param) const override;
		/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVectorForPoint3d(const DPoint3d& pt) const override;
	/// Returns the derivative (specified) vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const override;
	/// Returns the object-specific type
	virtual ElementType getType() const override { return SURFACE_CYLINDER; }
	/// invert the orientation of the surface (change diretion of normal vector) if possible
	virtual bool invertOrientation();
	std::shared_ptr<SurfaceParametric> adjustedForFaces( const DataVector< MeshFace* > & mfaces,
		const DataVector< MeshPoint3d* > & mpoints, MeshFace* central_face = nullptr  ) const override;
	// fix center and range (for parameters) to fix periodic surfaces
	virtual bool fixCenterAndRangeForPeriodic(const DPoint3d& center_point ) override;
	virtual bool withinParamRange( const DPoint2d& param ) const override;
public:
	const DPoint3d & getCenter() const { return m_cylinder.m_center; }
	const DVector3d & getAxisVt() const { return m_cylinder.m_axis_vt; }
	double getRadius() const { return m_cylinder.m_radius; }
protected:
	DPoint2d shiftedParams(const DPoint2d& param) const;
	DPoint2d unshiftedParams(const DPoint2d& param) const;
protected:
	/// Surface data
	DCylinder m_cylinder;
protected:
	bool m_non_periodic;
	double m_param_shift_x;
private:
	/// Parameterization ratio
	double m_dratio;
};

#endif // !defined(SURFACECYLINDER_H__INCLUDED)
