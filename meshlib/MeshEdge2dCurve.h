/////////////////////////////////////////////////////////////////////////////
// MeshEdge2dCurve.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef	MESHEDGECURVE_H_INCLUDED
#define MESHEDGECURVE_H_INCLUDED

#include "DPoint.h"
#include "MeshEdge2d.h"
#include "Curve2dParametric.h"

class SurfaceParametric;

/**
 * This class implements a curvilinear mesh edge in 2D 
 */
class MeshEdge2dCurve : public MeshEdge2d
{
public:
	/// Standard constructor
	MeshEdge2dCurve(MeshPoint2d *p1, MeshPoint2d *p2, 
		std::shared_ptr<const Curve2dParametric> new_shape, double t0 = 0.0, double t1 = 1.0,
		MeshElement* e1 = nullptr, MeshElement* e2 = nullptr);
public:
	/// Returns the parameter ksi for a given point (numerical approximation)
	virtual double getParameter(const DPoint2d& pt) const;
	/// Returns curvature of edge shape on plane
	virtual double getPlanarCurvature(double ksi) const;
	/// Returns curvature of edge shape on non-planar surface
	virtual double getNonPlanarCurvature(std::shared_ptr<const SurfaceParametric> surface, double ksi, double* cdt_len = nullptr) const;
	/// Update given ACS for curvature of edge contour
	virtual bool updateACSwithCurvature(std::shared_ptr<ControlSpace2dAdaptive> space, std::shared_ptr<const SurfaceParametric> surface, double d2_threshold) const;
	/// Create new edge with similar geometry (but possibly adjusted points topology and location)
	virtual MeshEdge2d* cloneGeometric(MeshPoint2d *p0, MeshPoint2d *p1, double ksi0 = 0.0, double ksi1 = 1.0, MeshElement* e1 = nullptr, MeshElement* e2 = nullptr) const;
	/// numerical approximation -> calculate parameter for point on surface+curve
	virtual DPoint2d surfaceParameters(std::shared_ptr<const SurfaceParametric> surface, const DPoint3d& point, double & ksi, bool same_direction) const;
	/// Returns the type of this object
	virtual ElementType getType() const { return EDGE_CURVE; }
	/// Includes this edge into bounding rectangle
	virtual void addToBoundingRect(DRect& rect) const;
	/// Returns the coordinates of point within this edge, described by parameter ksi [0,1]
	virtual DPoint2d getPoint(double ksi) const;
	/// Returns the length of a subsegment of this edge (current metric can be used)
	virtual double getLength(Metric2dContext& mc, double ksi0, double ksi1, bool local_metric = true) const;
	/// Returns the length of a subsegment of this edge (current metric can be used)
	virtual double getLengthMax(Metric2dContext& mc, double ksi0, double ksi1) const;
	/// Returns the length of a subsegment of this edge (current metric can be used)
	virtual double checkAndGetLength(Metric2dContext& mc, double ksi0, double& ksi1, double max_len = 1.0, bool local_metric = true) const;
	/// Returns the length of this edge (current metric can be used)
	virtual double getLength(Metric2dContext& mc, bool local_metric = true) const;
	/// Returns the approximation of this edge via polyline
	virtual void getPolyLine(DataVector<DPoint2d> & polyline) const;
	/// Returns the approximation of this edge via polyline (array of points)
	virtual void getPolyLine(DataVector<DPoint3d> & polyline, 
		std::shared_ptr<const SurfaceParametric> surface) const;
protected:
	/// Reverses the directionality of the edge
	virtual void switchSide();
protected:
	/// Reference to the class providing the extended information about shape
	std::shared_ptr<const Curve2dParametric>	m_shape;
	/// Parameteres (beginning and end) of the extended shape function
	double			m_t0, m_t1;
};

#endif // MESHEDGECURVE_H_INCLUDED
