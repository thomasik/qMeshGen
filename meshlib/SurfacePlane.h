/////////////////////////////////////////////////////////////////////////////
// SurfacePlane.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEPLANE_H__INCLUDED)
#define SURFACEPLANE_H__INCLUDED

class DPlane;
class MeshFace;
class MeshPoint3d;

#include "SurfaceParametric.h"
#include "DPoint.h"
#include "DataVector.h"

/**
 * This class implements a planar surface with several required methods.
 */
class SurfacePlane : public SurfaceParametric  
{
public:
	/// Standard constructor - three points
	SurfacePlane(const DPoint3d& pnt0, const DPoint3d& pnt1, const DPoint3d& pnt2);
	/// Standard constructor - point and two base vectors
	SurfacePlane(const DPoint3d& pnt0, const DVector3d& v01, const DVector3d& v02);
	/// Standard constructor - point and normal vector
	SurfacePlane(const DPoint3d& pnt0, const DVector3d& vn);
	/// Standard consturctor - from DPlane
	SurfacePlane(const DPlane& plane);
	/// Standard constructor - two base vectors and zero-point
	SurfacePlane(const DVector3d& v01, const DVector3d& v02);
	/// Standard constructor (empty)
	SurfacePlane() : SurfaceParametric() {}
	/// Standard destructor
	virtual ~SurfacePlane(){}
public:
	virtual string getSimpleDescription() const override { return "plane"; }
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
	/// Returns the principle curvature for this plane and parameters [s,t] (infinitezimal, since it is a plane)
	virtual const SurfaceCurvature getCurvature(const DPoint2d& /* pt */, double & g_ratio) const override { 
		g_ratio = m_dratio; return SurfaceCurvature(0.0, 0.0, 0.0); }
	virtual const DVector3d getCurvatureDirection(const DPoint2d& /* pt */, double & g_ratio) const override { 
		g_ratio = m_dratio; return m_v01; }
	/// Returns the point of the plane for the given parameters
	virtual const DPoint3d getPoint(const DPoint2d& param) const override;
	/// Returns the parameters of the plane for the given point
	virtual const DPoint2d getParameters(const DPoint3d& point) const override;
	/// Returns the parameters of the surface for the given point (numerical approximation) and signed distance (according to normal vecor)
	virtual bool getParametersAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z ) const override;
	/// Returns the parameters of the plane for the given point
	virtual const DPoint2d getParametersNear(const DPoint3d& point, const DPoint2d& /* near_point */) const override
		{ return getParameters(point); }
	/// Returns the shape (curve) parameter of the surface for the given point with starting point (numerical approximation)
	virtual double getShapeParameters(const DPoint3d& point, 
		std::shared_ptr<const Curve2dParametric> shape, double near_t, double min_t, double max_t) const override;
	/// Returns the segment parameter of the surface for the given point with starting point (numerical approximation)
	virtual double getSegmentParameters(const DPoint3d& point, const DPoint2d& pt0, const DPoint2d& pt1, double near_t, double min_t, double max_t) const override;
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVector(const DPoint2d& /* param */) const override { return m_vn; }
	/// Returns the derivative of normal vector for the given parameters
	virtual const DVector3d getNormalVectorDerivative(int /* deriv*/ , const DPoint2d& /* param */) const override {
		return DVector3d::zero;
	}
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVectorForPoint3d(const DPoint3d& /* pt */) const override { return m_vn; }
	/// Returns the derivative (specified) vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const override;
	/// Returns the length of a segment on the plane
	virtual double segmentLength(const DPoint2d & param_0, const DPoint2d & param_1) const override;
	/// Returns the approximation of a segment on surfaces via polyline (array of points)
	virtual void getPolyLine(DataVector<DPoint3d> & polyline, const DPoint2d& a, const DPoint2d& b) const override;
	/// Returns the approximation of a segment on surfaces via polyline (array of points)
	virtual void getPolyLine(DataVector<DPoint3d> & polyline, 
		std::shared_ptr<const Curve2dParametric> shape, double t0, double t1) const override;
	/// Returns the object-specific type
	virtual ElementType getType() const override { return SURFACE_PLANE; }
	std::shared_ptr<SurfaceParametric> adjustedForFaces( const DataVector< MeshFace* > & mfaces,
		const DataVector< MeshPoint3d* > & mpoints, MeshFace* central_face = nullptr ) const override;
public:
	/// Stores the parameters of this plane into the stream
	friend ostream& operator<<(ostream& os, const SurfacePlane* surf);
	/// Restores the parameters of this plane from the stream
	friend istream& operator>>(istream& is, SurfacePlane* surf);
	/// Sets plane parameters (three non-linear points)
	//bool setData(const DPoint3d& pnt0, const DPoint3d& pnt1, const DPoint3d& pnt2);
	/// Returns parameters of point (segment crossing plane)
	//bool getSegmentCrossingPointParam(const DPoint3d& pnt0, const DPoint3d& pnt1, DPoint2d& pnt) const;
public:
	const DPoint3d& baseP() const { return m_p0; }
	const DVector3d& baseEu() const { return m_v01; }
	const DVector3d& baseEv() const { return m_v02; }
protected:
	/// A point lying on the plane
	DPoint3d m_p0;
	/// First vector defining the plane
	DVector3d m_v01;
	/// Second vector defining the plane
	DVector3d m_v02;
	/// Parameterization coefficient
	double m_dratio;
private:
	/// Normal vector
	DVector3d m_vn; 
};

#endif // !defined(SURFACEPLANE_H__INCLUDED)
