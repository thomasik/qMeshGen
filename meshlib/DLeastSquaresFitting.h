/////////////////////////////////////////////////////////////////////////////
// DLeastSquaresFitting.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DLEASTSQUARESFITTING_H__INCLUDED
#define DLEASTSQUARESFITTING_H__INCLUDED

class DPoint2d;
class DPoint3d;
class DVector2d;
class DVector3d;
class DPlane;
class DQuadric;
class DPlanarQuadric;
class DQuadricOnSurface;
class DLinearQuadric;
class DLinearQuadric3d;
class DLine2d;
class DCircle;
class DSphere;
class DCylinder;
class SurfaceParametric;

#include <memory>
#include "DataVector.h"
#include "DataMatrix.h"

#define LS_FIT_ERROR -1.0

/**
 * This class implements a least square approximation.
 */
class DLeastSquaresFitting
{
public:
	/// Least Squares linear fitting of points using orthogonal regression, returns max distance
	static double fitLine(const DataVector<DPoint2d> & points, DLine2d& line, bool ret_max_dist = true);
	/// Least Squares fitting of orthogonal vector, returns max distance
	static double fitVectorOrthonormal(const DataVector<DVector3d> & vectors, DVector3d& vn);
	/// Least Squares hyperplanar fitting of points using orthogonal regression, returns max distance
	static double fitHyperplaneOrthogonal(const DataVector<DPoint3d> & points,
		DPlane& plane, bool oriented_contour = false, bool ret_max_dist = true );
	/// Least Squares fitting of points to circle, returns max distance
	static double fitCircleMinArea(const DataVector<DPoint2d> & points, 
				DCircle& circle, bool ret_max_dist = true);
	/// Least Squares fitting of points to circle, returns max distance
	static double fitCircle(const DataVector<DPoint2d> & points, DCircle& circle, bool ret_max_dist = true) { 
					return fitCircleBullock( points, circle, ret_max_dist); }
	/// Least Squares fitting of points to circle, returns max distance
	static double fitCircleBullock(const DataVector<DPoint2d> & points, 
				DCircle& circle, bool ret_max_dist = true);
	/// Least Squares fitting of points to circle, iterative, returns max distance
	static double fitCircleIter(const DataVector<DPoint2d> & points, 
				DCircle& circle, bool ret_max_dist = true);
	/// Least Squares fitting of points to sphere, returns max distance
	static double fitSphere(const DataVector<DPoint3d> & points, DSphere& sphere, bool ret_max_dist = true) { 
					return fitSphereIter( points, sphere, true, ret_max_dist); }
	/// Least Squares fitting of points to sphere, iterative, returns max distance
	static double fitSphereIter(const DataVector<DPoint3d> & points, 
				DSphere& sphere, bool calculate_initial = true, bool ret_max_dist = true);
	/// Least Squares fitting of points to cylinder, using projection on plane, returns max distance
	static double fitCylinder(const DataVector<DPoint3d> & points, const DPlane& plane,
				DCylinder& cylinder);
	static double fitCylinder(const DataVector<DPoint3d> & points, const DVector3d& caxis,
				DCylinder& cylinder);
	/// Least Squares quadric fitting of points, returns max distance
	static double fitQuadric(const DataVector<DPoint3d> & points, DQuadric& quadric);
	/// Least Squares planar-quadric fitting of points for the given base surface, returns max distance
	static double fitQuadricOnSurface(const DataVector<DPoint3d> & points, 
		std::shared_ptr<SurfaceParametric> surface, DQuadricOnSurface& squadric);
	/// Least Squares planar-quadric fitting of points for the given plane, returns max distance
	static double fitPlanarQuadric(const DataVector<DPoint3d> & points, const DPlane& plane, DPlanarQuadric& pquadric);
	/// Least Squares planar-quadric fitting of points, returns max distance
	static double fitPlanarQuadric(const DataVector<DPoint3d> & points, DPlanarQuadric& pquadric);
	/// Least Squares linear-quadric fitting of points, returns max distance
	static double fitLinearQuadric(const DataVector<DPoint2d> & points, DLinearQuadric& lquadric, bool ret_max_dist = true);
	/// Least Squares linear-quadric fitting of points for the given line, returns max distance
	static double fitLinearQuadric(const DataVector<DPoint2d> & points, const DLine2d& line, 
				DLinearQuadric& lquadric, bool ret_max_dist = true);
	/// Least Squares linear-quadric fitting of points, returns max distance
	static double fitLinearQuadric(const DataVector<DPoint3d> & points, DLinearQuadric3d& lquadric, bool ret_max_dist = true);
	/// Least Squares linear-quadric fitting of points for the given plane and line, returns max distance
	static double fitLinearQuadric(const DataVector<DPoint3d> & points, 
				const DPlane & plane, DLinearQuadric3d& lquadric, bool ret_max_dist = true);
	/// Least Squares fitting of grid, returns max distance
	static double fitRegularGrid(const DataVector<DPoint3d> & points, 
		DataMatrix<DPoint3d> & grid);
	/// Least Squares fitting of point to several planes (in form <vn, d>), returns max distance
	static double fitPointToPlanes(DPoint3d& pt, 
		const DataVector<DVector3d> & vns, const DataVector<double> & ds);
	/// test fitting
	static void testFit();
};

#endif // !defined(DLEASTSQUARESFITTING_H__INCLUDED)
