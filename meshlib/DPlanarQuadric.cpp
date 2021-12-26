/////////////////////////////////////////////////////////////////////////////
// DPlanarQuadric.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DPlanarQuadric.h"
#include "DVectorN.h"
#include "DPoint.h"
#include "MeshViewSet.h"
#include "DEquation.h"

/// Standard constructor
DPlanarQuadric::DPlanarQuadric() : vq(0.0) { }

/// count distance of point and quadric
double DPlanarQuadric::distance(const DPoint3d& pt) const
{
	DPoint2d pt2d = plane.projectToPlane(pt);
	return pt.distance( getPoint(pt2d) );
}

DPoint3d DPlanarQuadric::getPoint(const DPoint2d& pt2d) const
{
	return getPoint( pt2d, plane, vq );
}

DPoint3d DPlanarQuadric::getPoint(const DPoint2d& param, const DPlane& _plane, const DVectorN<N> & _vq)
{
	double dist = _vq[0] +
		param.x * _vq[1] +
		param.y * _vq[2] +
		param.x * param.x * _vq[3] +
		param.y * param.y * _vq[4] +
		param.x * param.y * _vq[5];
	return _plane.projectToSpace(param) + _plane.vn * dist;
}

DVector3d DPlanarQuadric::getDerivative(int deriv, const DPoint2d& pt2d) const
{
	return getDerivative( deriv, pt2d, plane, vq );
}

DVector3d DPlanarQuadric::getDerivative(int deriv, const DPoint2d& param, const DPlane& _plane, const DVectorN<N> & _vq)
{
	switch(deriv){
	case DEquation::deriv_ds:
		return _plane.e0 + _plane.vn * 
			(_vq[1] + 2 * param.x * _vq[3] + param.y * _vq[5]);
	case DEquation::deriv_dt:
		return _plane.e1 + _plane.vn * 
			(_vq[2] + 2 * param.y * _vq[4] + param.x * _vq[5]);
	case DEquation::deriv_dss:
		return _plane.vn * (2 * _vq[3]);
	case DEquation::deriv_dst:
		return _plane.vn * _vq[5];
	case DEquation::deriv_dtt:
		return _plane.vn * (2 * _vq[4]);
	case DEquation::deriv_dsss:
	case DEquation::deriv_dsst:
	case DEquation::deriv_dstt:
	case DEquation::deriv_dttt:
		return DVector3d::zero;
	default:
		assert(false);
	}
	return DVector3d::zero;
}

/// solve Q( pt + t * vt) == 0
bool DPlanarQuadric::solve(const DPoint3d& pt, const DVector3d& vt, double& t) const
{
	assert(false);
/*
	double a =  vq[4] * vt.x * vt.x + vq[7] * vt.x * vt.y + vq[5] * vt.y * vt.y + 
				vq[8] * vt.x * vt.z + vq[9] * vt.y * vt.z + vq[6] * vt.z * vt.z;
	double b =  vq[1] * vt.x + 2 * vq[4] * pt.x * vt.x + vq[7] * vt.x * pt.y + 
				vq[2] * vt.y + vq[7] * pt.x * vt.y + 2 * vq[5] * pt.y * vt.y + 
				vq[8] * vt.x * pt.z + vq[9] * vt.y * pt.z + vq[3] * vt.z + 
				vq[8] * pt.x * vt.z + vq[9] * pt.y * vt.z + 2 * vq[6] * pt.z * vt.z;
	double c =  vq[0] + vq[1] * pt.x + vq[4] * pt.x * pt.x + vq[2] * pt.y + 
				vq[7] * pt.x * pt.y + vq[5] * pt.y * pt.y + vq[3] * pt.z + 
				vq[8] * pt.x * pt.z + vq[9] * pt.y * pt.z + vq[6] * pt.z * pt.z;

	const double EPS = SMALL_NUMBER * std::max(std::max(abs(a), abs(b)), abs(c));

	if(abs(a) < EPS){ // b x + c = 0
		if(abs(b) < EPS) return false;
		t1 = t2 = -c / b;
		return true;
	}else{ 
		// quadratic
		double delta = b*b - 4*a*c;
		if(delta < -EPS) return false; // no real solutions
		if(delta < EPS){ // delta close to 0
			t1 = t2 = -b / (2*a);
		}else{
			delta = sqrt(delta);
			t1 = (-b-delta) / (2*a);
			t2 = (-b+delta) / (2*a);
		}
	}
	return true;
*/
	return false;
}

/// create sketchy representation of this surface for the area with the given points
MeshViewSet * DPlanarQuadric::createViewSetForPoints(MeshViewSet* set, const DataVector<DPoint3d> & points) const
{
	size_t pct = points.countInt();
	DataVector<DPoint2d> points2d(pct);
	DRect box;
	for(size_t i = 0; i < pct; i++){
		set->addPoint(points[i]);
		points2d.add( plane.projectToPlane(points[i]));
		box.addPoint(points2d[i]);
	}

	DPoint2d ptx = box.getX0Y0();
	DPoint2d pty = box.getX0Y0();
	DVector2d dx(box.getDX() / (SKETCH_LINES-1), 0.0);
	DVector2d dy(0.0, box.getDY() / (SKETCH_LINES-1));


	for(int i = 0; i < SKETCH_LINES; i++) {
		DPoint2d ptxi = ptx;
		DPoint2d ptyi = pty;
		for(int j = 1; j < SKETCH_LINES; j++) {
			set->addEdge( getPoint(ptxi), getPoint(ptxi + dy));
			set->addEdge( getPoint(ptyi), getPoint(ptyi + dx));
			ptxi += dy;
			ptyi += dx;
		}
		ptx += dx;
		pty += dy;
	}

	return set;
}

bool DPlanarQuadric::invertSurfaceQuadricOrientation( DVectorN<N> & vq )
{
	for(int i = 0; i < 6; i++)
		vq[i] = -vq[i];
	vq[1] = -vq[1]; // x
	vq[5] = -vq[5]; // xy
	return true;
}

/// invert the orientation of the surface (change diretion of normal vector)
bool DPlanarQuadric::invertOrientation()
{
	return plane.switchOrientation() && invertSurfaceQuadricOrientation( vq );
}

DOrientedBox DPlanarQuadric::getOrientedBox() const
{
	return plane.getOrientedBox();
}

DOrientedBox DPlanarQuadric::getOrientedBox( const DataVector<DPoint3d>& points ) const
{
	return plane.getOrientedBox( points );
}

DOrientedBox DPlanarQuadric::getOrientedBoxOpt( const DataVector<DPoint3d>& points ) const
{
	return plane.getOrientedBoxOpt( points );
}

/// Standard constructor
DQuadricOnSurface::DQuadricOnSurface( SurfacePtr _base_surface) : vq(0.0), base_surface(_base_surface) { }

/// count distance of point and quadric
double DQuadricOnSurface::distance(const DPoint3d& pt) const
{
	assert( base_surface );
	DPoint2d pt2d = base_surface->getParameters(pt);
	return pt.distance( getPoint(pt2d) );
}

DPoint3d DQuadricOnSurface::getPoint(const DPoint2d& pt2d) const
{
	assert( base_surface );
	return getPoint( pt2d, base_surface, vq );
}

DPoint3d DQuadricOnSurface::getPoint(const DPoint2d& param, SurfaceConstPtr surface, const DVectorN<N> & _vq)
{
	assert( surface );
	double dist = _vq[0] +
		param.x * _vq[1] +
		param.y * _vq[2] +
		param.x * param.x * _vq[3] +
		param.y * param.y * _vq[4] +
		param.x * param.y * _vq[5];

	return surface->getPoint(param) + surface->getNormalVector(param) * dist;
}

DVector3d DQuadricOnSurface::getDerivative(int deriv, const DPoint2d& pt2d) const
{
	assert( base_surface );
	return getDerivative( deriv, pt2d, base_surface, vq );
}

DVector3d DQuadricOnSurface::getDerivative(int deriv, const DPoint2d& param, SurfaceConstPtr surface, const DVectorN<N> & _vq)
{
	double dist = _vq[0] + 
		param.x * _vq[1] + 
		param.y * _vq[2] +
		param.x * param.x * _vq[3] + 
		param.y * param.y * _vq[4] +
		param.x * param.y * _vq[5];

	switch(deriv){
	case DEquation::deriv_ds:
	{
		double dist_ds = (_vq[1] + 2 * param.x * _vq[3] + param.y * _vq[5]);
		return surface->getDerivative(deriv, param) +
			surface->getNormalVector(param) * dist_ds +
			surface->getNormalVectorDerivative(deriv, param) * dist;
	}
	case DEquation::deriv_dt:
	{
		double dist_dt = (_vq[2] + 2 * param.y * _vq[4] + param.x * _vq[5]);
		return surface->getDerivative(deriv, param) +
			surface->getNormalVector(param) * dist_dt +
			surface->getNormalVectorDerivative(deriv, param) * dist;
	}
	case DEquation::deriv_dss:
	{
		double dist_ds = (_vq[1] + 2 * param.x * _vq[3] + param.y * _vq[5]);
		double dist_dss = (2 * _vq[3]);
		return surface->getDerivative(deriv, param) +
			surface->getNormalVectorDerivative(DEquation::deriv_dss, param) * dist +
			surface->getNormalVectorDerivative(DEquation::deriv_ds, param) * (2 * dist_ds) +
			surface->getNormalVector(param) * dist_dss;
	}
	case DEquation::deriv_dst:
	{
		double dist_ds = (_vq[1] + 2 * param.x * _vq[3] + param.y * _vq[5]);
		double dist_dt = (_vq[2] + 2 * param.y * _vq[4] + param.x * _vq[5]);
		double dist_dst = _vq[5];
		return surface->getDerivative(deriv, param) +
			surface->getNormalVectorDerivative(DEquation::deriv_dst, param) * dist +
			surface->getNormalVectorDerivative(DEquation::deriv_ds, param) * dist_dt +
			surface->getNormalVectorDerivative(DEquation::deriv_dt, param) * dist_ds +
			surface->getNormalVector(param) * dist_dst;
	}
	case DEquation::deriv_dtt:
	{
		double dist_dt = (_vq[2] + 2 * param.y * _vq[4] + param.x * _vq[5]);
		double dist_dtt = (2 * _vq[4]);
		return surface->getDerivative(deriv, param) +
			surface->getNormalVectorDerivative(DEquation::deriv_dtt, param) * dist +
			surface->getNormalVectorDerivative(DEquation::deriv_dt, param) * (2 * dist_dt) +
			surface->getNormalVector(param) * dist_dtt;
	}
	case DEquation::deriv_dsss:
	case DEquation::deriv_dsst:
	case DEquation::deriv_dstt:
	case DEquation::deriv_dttt:
		return surface->getDerivative(deriv, param);
	default:
		assert(false);
	}
	return DVector3d::zero;
}

/// solve Q( pt + t * vt) == 0
bool DQuadricOnSurface::solve(const DPoint3d& pt, const DVector3d& vt, double& t) const
{
	assert(false);
/*
	double a =  vq[4] * vt.x * vt.x + vq[7] * vt.x * vt.y + vq[5] * vt.y * vt.y + 
				vq[8] * vt.x * vt.z + vq[9] * vt.y * vt.z + vq[6] * vt.z * vt.z;
	double b =  vq[1] * vt.x + 2 * vq[4] * pt.x * vt.x + vq[7] * vt.x * pt.y + 
				vq[2] * vt.y + vq[7] * pt.x * vt.y + 2 * vq[5] * pt.y * vt.y + 
				vq[8] * vt.x * pt.z + vq[9] * vt.y * pt.z + vq[3] * vt.z + 
				vq[8] * pt.x * vt.z + vq[9] * pt.y * vt.z + 2 * vq[6] * pt.z * vt.z;
	double c =  vq[0] + vq[1] * pt.x + vq[4] * pt.x * pt.x + vq[2] * pt.y + 
				vq[7] * pt.x * pt.y + vq[5] * pt.y * pt.y + vq[3] * pt.z + 
				vq[8] * pt.x * pt.z + vq[9] * pt.y * pt.z + vq[6] * pt.z * pt.z;

	const double EPS = SMALL_NUMBER * std::max(std::max(abs(a), abs(b)), abs(c));

	if(abs(a) < EPS){ // b x + c = 0
		if(abs(b) < EPS) return false;
		t1 = t2 = -c / b;
		return true;
	}else{ 
		// quadratic
		double delta = b*b - 4*a*c;
		if(delta < -EPS) return false; // no real solutions
		if(delta < EPS){ // delta close to 0
			t1 = t2 = -b / (2*a);
		}else{
			delta = sqrt(delta);
			t1 = (-b-delta) / (2*a);
			t2 = (-b+delta) / (2*a);
		}
	}
	return true;
*/
	return false;
}

/// create sketchy representation of this surface for the area with the given points
MeshViewSet * DQuadricOnSurface::createViewSetForPoints(MeshViewSet* set, const DataVector<DPoint3d> & points,
	SurfaceConstPtr surface, const DVectorN<N> & _vq, int id)
{
	if( set == nullptr) set = new MeshViewSet;

	size_t pct = points.countInt();
	DataVector<DPoint2d> points2d(pct);
	DRect box;
	for(size_t i = 0; i < pct; i++){
		set->addPoint(points[i]);
		points2d.add( surface->getParameters(points[i]));
		box.addPoint(points2d[i]);
	}

	DPoint2d ptx = box.getX0Y0();
	DPoint2d pty = box.getX0Y0();
	DVector2d dx(box.getDX() / (SKETCH_LINES-1), 0.0);
	DVector2d dy(0.0, box.getDY() / (SKETCH_LINES-1));


	for(int i = 0; i < SKETCH_LINES; i++) {
		DPoint2d ptxi = ptx;
		DPoint2d ptyi = pty;
		for(int j = 1; j < SKETCH_LINES; j++) {
			set->addEdge( getPoint(ptxi, surface, _vq), getPoint(ptxi + dy, surface, _vq), id);
			set->addEdge( getPoint(ptyi, surface, _vq), getPoint(ptyi + dx, surface, _vq), id);
			ptxi += dy;
			ptyi += dx;
		}
		ptx += dx;
		pty += dy;
	}

	return set;
}

/// create sketchy representation of this surface for the area with the given points
MeshViewSet * DQuadricOnSurface::createViewSetForPoints(MeshViewSet* set, 
						const DataVector<DPoint3d> & points, int id) const
{
	return createViewSetForPoints( set, points, base_surface, vq, id);
}

/// invert the orientation of the surface (change diretion of normal vector)
bool DQuadricOnSurface::invertOrientation()
{
	return base_surface->invertOrientation() && DPlanarQuadric::invertSurfaceQuadricOrientation( vq );
}

DOrientedBox DQuadricOnSurface::getOrientedBox() const
{
	return base_surface->getOrientedBox();
}

DOrientedBox DQuadricOnSurface::getOrientedBox( const DataVector<DPoint3d>& points ) const
{
	return base_surface->getOrientedBox( points );
}

DOrientedBox DQuadricOnSurface::getOrientedBoxOpt( const DataVector<DPoint3d>& points ) const
{
	return base_surface->getOrientedBoxOpt( points );
}
