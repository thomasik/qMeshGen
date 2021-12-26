/////////////////////////////////////////////////////////////////////////////
// DPlanarQuadric.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DPLANARQUADRIC_H__INCLUDED
#define DPLANARQUADRIC_H__INCLUDED

#include "DVectorN.h"
#include "DPlane.h"
#include "DOrientedBox.h"
#include "DataMatrix.h"
#include "SurfaceParametric.h"

class DPoint3d;
class DVector3d;
class MeshViewSet;

/**
 * This class implements a quadric surface representation for plane: z(x,y)
 */
class DPlanarQuadric
{
public:
	static const int N = 6;
	static const int SKETCH_LINES = 10;
public:
	/// Standard constructor
	DPlanarQuadric();
	// v=(1,x,y,xx,yy,xy)
public:
	DOrientedBox getOrientedBox() const;
	DOrientedBox getOrientedBox( const DataVector<DPoint3d>& points ) const;
	DOrientedBox getOrientedBoxOpt( const DataVector<DPoint3d>& points ) const;
	static bool invertSurfaceQuadricOrientation( DVectorN<N> & vq );
	/// invert the orientation of the surface (change diretion of normal vector)
	bool invertOrientation();
	/// create sketchy representation of this surface for the area with the given points
	MeshViewSet * createViewSetForPoints(MeshViewSet* set, const DataVector<DPoint3d> & points) const;
	/// point on the surface
	DPoint3d getPoint(const DPoint2d& param) const;
	/// point on the surface
	static DPoint3d getPoint(const DPoint2d& param, const DPlane& _plane, const DVectorN<N> & _vq);
	/// point on the surface
	DVector3d getDerivative(int deriv, const DPoint2d& param) const;
	/// point on the surface
	static DVector3d getDerivative(int deriv, const DPoint2d& param, const DPlane& _plane, const DVectorN<N> & _vq);
	/// count distance of point and quadric
	double distance(const DPoint3d& pt) const;
	/// solve Q( pt + t * vt) == 0
	bool solve(const DPoint3d& pt, const DVector3d& vt, double& t) const;
	/// Stores the parameters of this planar quadric into the stream
	friend ostream& operator<<(ostream& os, const DPlanarQuadric& pq) {
		return os << pq.plane << " " << pq.vq;
	}
	/// Restores the parameters of this planar quadric from the stream
	friend istream& operator>>(istream& is, DPlanarQuadric& pq) {
		return is >> pq.plane >> ws >> pq.vq;
	}
public:
	/// coefficients
	DVectorN<N> vq;
	DPlane plane;
};

/**
 * This class implements a quadric surface representation for surface z(x,y)
 */
class DQuadricOnSurface
{
public:
	static const int N = 6;
	static const int SKETCH_LINES = 10;
public:
	/// Standard constructor
	DQuadricOnSurface( SurfacePtr _base_surface = nullptr );
	// v=(1,x,y,xx,yy,xy)
public:
	DOrientedBox getOrientedBox() const;
	DOrientedBox getOrientedBox( const DataVector<DPoint3d>& points ) const;
	DOrientedBox getOrientedBoxOpt( const DataVector<DPoint3d>& points ) const;
	/// invert the orientation of the surface (change diretion of normal vector)
	bool invertOrientation();
	/// create sketchy representation of this surface for the area with the given points
	MeshViewSet * createViewSetForPoints(MeshViewSet* set, const DataVector<DPoint3d> & points, int id = -2) const;
	/// create sketchy representation of this surface for the area with the given points
	static MeshViewSet * createViewSetForPoints(MeshViewSet* set, const DataVector<DPoint3d> & points,
		SurfaceConstPtr surface, const DVectorN<N> & _vq, int id = -2);
	/// point on the surface
	DPoint3d getPoint(const DPoint2d& param) const;
	/// point on the surface
	static DPoint3d getPoint(const DPoint2d& param, SurfaceConstPtr surface, const DVectorN<N> & _vq);
	/// point on the surface
	DVector3d getDerivative(int deriv, const DPoint2d& param) const;
	/// point on the surface
	static DVector3d getDerivative(int deriv, const DPoint2d& param, SurfaceConstPtr surface, const DVectorN<N> & _vq);
	/// count distance of point and quadric
	double distance(const DPoint3d& pt) const;
	/// solve Q( pt + t * vt) == 0
	bool solve(const DPoint3d& pt, const DVector3d& vt, double& t) const;
	/// Stores the parameters of this planar quadric into the stream
	friend ostream& operator<<(ostream& os, const DQuadricOnSurface& qs) {
		qs.base_surface->storeXML(os);
		return os << " " << qs.vq;
	}
	/// Restores the parameters of this planar quadric from the stream
	friend istream& operator>>(istream& is, DQuadricOnSurface& qs) {
		assert(false);
		//qs.base_surface->loadXML(is);
		return is >> ws >> qs.vq;
	}
public:
	/// coefficients
	DVectorN<N> vq;
	SurfacePtr base_surface;
};

#endif // !defined(DPLANARQUADRIC_H__INCLUDED)
