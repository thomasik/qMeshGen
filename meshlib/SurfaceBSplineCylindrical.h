/////////////////////////////////////////////////////////////////////////////
// SurfaceBSplineCylindrical.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(SURFACEBSPLINECYLINDRICAL_H__INCLUDED)
#define SURFACEBSPLINECYLINDRICAL_H__INCLUDED

#include "SurfaceParametric.h"
#include "SurfaceCylinder.h"
#include "DRect.h"
#include "DataMatrix.h"

class SurfaceMulti;
class MeshViewSet;

/**
 * This class implements a BSpline surface with several required methods.
 */
class SurfaceBSplineCylindrical : public SurfaceParametric  
{
public:
	/// Standard constructor
	SurfaceBSplineCylindrical(const SurfaceCylinder& cylinder, const DRect& rect);
	/// Standard constructor
	SurfaceBSplineCylindrical(const SurfaceCylinder& cylinder, DataMatrix<DPoint3d> *points);
	/// Standard destructor
	virtual ~SurfaceBSplineCylindrical();
public:
	virtual string getSimpleDescription() const override { return "bspline-cylindrical"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the point of the surface for the given parameters
	virtual const DPoint3d getPoint(const DPoint2d& param) const;
	/// Returns the parameters of the surface for the given point with starting point
	virtual const DPoint2d getParametersNear(const DPoint3d& point, const DPoint2d& near_point) const;
	/// Returns the parameters of the surface for the given point
	virtual const DPoint2d getParameters(const DPoint3d& point) const;
	/// Returns the normal vector to surface for the given parameters
	virtual const DVector3d getNormalVector(const DPoint2d& param) const;
	/// Returns the derivative (specified) vector for the given parameters
	virtual const DVector3d getDerivative(int deriv, const DPoint2d& param) const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return SURFACE_BSPLINE; }
public:
	static const int FIT_PROBE_PTS = 20;
	static const int FIT_PROBE_MIN = 10;
public:
	/// additional fit
	double fitAdditionally(const DataVector<DPoint3d> & points);
	/// fit to point-cloud using direct fit, return max-dist (ortogonal to plane)
	double fitToPoints(const DataVector<DPoint3d> & points, int init_res = 2);
	/// get view set for visualization
	MeshViewSet* getViewSet(MeshViewSet* set = nullptr, bool with_knots = false) const;
private:
	/// fit to point-cloud - border adjust
	void fitToPointsAdjustBorder(int level);
private:
	/// Calculate functional knots for this B-spline
	bool calculateNodes();
	/// Calculate difference (max square distance) between the surface and set of points
	double calculateDiff2(const DataVector<DPoint3d>& points) const;
	/// Calculates the number of B-Spline row and column and sets ui/vi as a parameter inside this segment
	bool countSplineIndices(double& u, double& v, int& ui, int& vi) const;
protected:
	/// structure for point-ordering
	struct NodeDist {
		NodeDist(int _id = 0, double _dist2 = 0.0) : id(_id), dist2(_dist2) {}
		int id;
		double dist2;
	};
protected:
	/// A base plane (e0 and e1 -> rectangular boundary)
	SurfaceCylinder m_cylinder;
	/// Knots (control)
	DataMatrix<DPoint3d> * m_nodes;
	/// Knots (real)
	DataMatrix<DPoint3d> * m_points;
};

#endif // !defined(SURFACEBSPLINECYLINDRICAL_H__INCLUDED)
