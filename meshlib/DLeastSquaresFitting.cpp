/////////////////////////////////////////////////////////////////////////////
// DLeastSquaresFitting.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DLeastSquaresFitting.h"
#include "DRect.h"
#include "DMatrix.h"
#include "DataMatrix.h"
#include "DVectorN.h"
#include "DMatrixN.h"
#include "DQuadric.h"
#include "DLine.h"
#include "DPlane.h"
#include "DPlanarQuadric.h"
#include "DLinearQuadric.h"
#include "SurfacePlane.h"
#include "SurfaceCylinder.h"
#include "DSphere.h"

#include "MeshViewSet.h"

/// Least Squares hyperplanar fitting of points using orthogonal regression, returns max distance
double DLeastSquaresFitting::fitHyperplaneOrthogonal(
	const DataVector<DPoint3d> & points,
	DPlane & plane, bool oriented_contour, bool ret_max_dist)
{
	size_t pct = points.countInt();
	// -- ptA -> average of the sample points
	DPoint3d ptA = DPoint3d::zero;
	// -- box -> bounding box
	DBox box;
	for(size_t i = 0; i < pct; i++){
		const DPoint3d& pt = points[i];
		ptA.add(pt);
		box.addPoint(pt);
	}
	double dd = box.getDiameter();
	if(dd == 0.0){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "fitHyperplaneOrthogonal - degenerate points");
		return LS_FIT_ERROR;
	}
	ptA /= (double)pct;

	// M(A) -- second form of the energy function
	DMatrix3d mA;
	for(size_t i = 0; i < pct; i++){
		const DVector3d dv = points[i] - ptA;
		mA.m[0][0] += dv.x * dv.x;
		mA.m[0][1] += dv.x * dv.y;
		mA.m[0][2] += dv.x * dv.z;
		mA.m[1][0] += dv.y * dv.x;
		mA.m[1][1] += dv.y * dv.y;
		mA.m[1][2] += dv.y * dv.z;
		mA.m[2][0] += dv.z * dv.x;
		mA.m[2][1] += dv.z * dv.y;
		mA.m[2][2] += dv.z * dv.z;
	}
	DMatrix3d mA_vectors;
	double mA_values[3];
	mA.eigensystem(mA_vectors, mA_values);

	if( (abs(mA_values[0]) < SMALL_NUMBER * dd) &&
		(abs(mA_values[1]) < SMALL_NUMBER * dd) &&
		(abs(mA_values[2]) < SMALL_NUMBER * dd))
	{
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "fitHyperplaneOrthogonal - degenerate solution");
		return LS_FIT_ERROR;
	}

	// -- the smallest eigenvalue of M(A)
	int min_i = ( abs(mA_values[0]) < abs(mA_values[1])) ? 0 : 1;
	min_i = ( abs(mA_values[min_i]) < abs(mA_values[2])) ? min_i : 2;

	// -- solution
	DVector3d n = mA_vectors.column(min_i).normalized();
	double D =  -(n.x * ptA.x + n.y * ptA.y + n.z * ptA.z);

	plane.p0 = ptA;

	// -- minimum of energy function (compared to the diameter of the point cloud)
	// LOG4CPLUS_INFO(MeshLog::logger_mesh, "min eigenvalue (energy function) => " << mA_values[min_i]);
	double max_dist = 0.0;
	if( ret_max_dist ){
		for(size_t i = 0; i < pct; i++){
			const DPoint3d& pt = points[i];
			double dist = abs(n.x * pt.x + n.y * pt.y + n.z * pt.z + D);
			if(dist > max_dist){
				max_dist = dist;
			}
		}
	}else{
		max_dist = mA_values[min_i];
	}

	if(oriented_contour){
		// check normal-vector orientation
		const DVector3d contour_n = DVector3d::crossProduct(points[0], points[1], points[pct-1]);
		if(contour_n.length2() > VERY_SMALL_NUMBER){
			if(n.scalarProduct(contour_n) > 0.0) n *= -1.0;
			// associate e0 with first edge of contour
			plane.e0 = n.crossProduct(points[1] - points[0]);
		}else{
			plane.e0 = n.crossProduct(box.getMaxPoint() - ptA);
		}
		plane.e1 = n.crossProduct(plane.e0).normalized() * (0.5 * dd);
		plane.e0 = n.crossProduct(plane.e1).normalized() * (0.5 * dd);
		plane.vn = -n;
		// check point-sequence orientation
		DataVector<DPoint2d> vpoints(pct);
		for(size_t i = 0; i < pct; i++)
			vpoints.add(plane.projectToPlane(points[i]));

		DPoint2d::clearZeroEdges(vpoints);
		if(!DPoint2d::properOrientation(vpoints)){
			// reverse the vectors
			plane.switchOrientation();
		}
	}else{
		plane.e0 = n.crossProduct(box.getMaxPoint() - ptA);
		plane.e1 = n.crossProduct(plane.e0).normalized() * (0.5 * dd);
		plane.e0 = n.crossProduct(plane.e1).normalized() * (0.5 * dd);
		plane.vn = -n;
	}
	return max_dist;
}

/// Least Squares quadric fitting of points, returns max distance
double DLeastSquaresFitting::fitQuadric(const DataVector<DPoint3d> & points, DQuadric& quadric)
{
	const int N = 10;
	// M(A) -- second form of the energy function
	DMatrixN<N> mA(0.0);
	size_t pct = points.countInt();
	for(size_t i = 0; i < pct; i++){
		const DPoint3d& pt = points[i];
		const double v[N] = { 1, pt.x, pt.y, pt.z,
			pt.x * pt.x, pt.y * pt.y, pt.z * pt.z,
			pt.x * pt.y, pt.x * pt.z, pt.y * pt.z };
		for(int j = 0; j < N; j++)
			for(int k = 0; k < N; k++)
				mA(j,k) += v[j] * v[k];
	}

	double min_energy = mA.smallestEigenvectorNormalized(quadric.vq);
	if(min_energy < 0.0) return LS_FIT_ERROR;

	// -- minimum of energy function (compared to the diameter of the point cloud)
	// LOG4CPLUS_INFO(MeshLog::logger_mesh, "min eigenvalue (energy function) => " << mA_values[min_i]);
	double max_dist = 0.0;
	for(size_t i = 0; i < pct; i++){
		double dist = quadric.distance(points[i]);
		if(dist > max_dist) max_dist = dist;
	}

	return max_dist;
}

/// Least Squares quadric fitting of points for the given plane, returns max distance
double DLeastSquaresFitting::fitPlanarQuadric(const DataVector<DPoint3d> & points, 
			const DPlane& plane, DPlanarQuadric& pquadric)
{
	pquadric.plane = plane;

	// M(A) -- second form of the energy function
	DMatrixN<DPlanarQuadric::N> mA(0.0);
	DVectorN<DPlanarQuadric::N> mb(0.0);
	size_t pct = points.countInt();
	DPoint2d pt;
	double z;
	for(size_t i = 0; i < pct; i++){
		if(!pquadric.plane.projectToPlane(points[i], pt, z)) 
			continue;
		const double v[DPlanarQuadric::N] = { 1, pt.x, pt.y, 
			pt.x * pt.x, pt.y * pt.y, pt.x * pt.y };
		for(int j = 0; j < DPlanarQuadric::N; j++){
			for(int k = 0; k < DPlanarQuadric::N; k++)
				mA(j,k) += v[j] * v[k];
			mb[j] += v[j] * z;
		}
	}

	if(! mA.solveLU(mb, pquadric.vq)) return LS_FIT_ERROR;

	// -- minimum of energy function (compared to the diameter of the point cloud)
	// LOG4CPLUS_INFO(MeshLog::logger_mesh, "min eigenvalue (energy function) => " << mA_values[min_i]);
	double dist, max_dist = 0.0;
	for(size_t i = 0; i < pct; i++){
		dist = pquadric.distance(points[i]);
		if(dist > max_dist) max_dist = dist;
	}

#ifdef _DEBUG
	if(  pct < 10 ){
		MeshViewSet* set = new MeshViewSet;
		pquadric.createViewSetForPoints(set, points);
		set->addInfo("pct", (int)pct);
		set->addInfo("max_dist", max_dist);
		SHOW_MESH( " LSF::planar-quadric with few points", set );
	}
#endif

	return max_dist;
}

//#define FQOS_SHOW

/// Least Squares quadric fitting of points for the given base surface, returns max distance
double DLeastSquaresFitting::fitQuadricOnSurface(
		const DataVector<DPoint3d> & points, 
		SurfacePtr surface, 
		DQuadricOnSurface& squadric)
{
	squadric.base_surface = surface;
	int pct = points.countInt();

#ifdef FQOS_SHOW
	MeshViewSet* set = new MeshViewSet;
	DataVector<DPoint2d> params(pct);
#endif

	// M(A) -- second form of the energy function
	DMatrixN<DPlanarQuadric::N> mA(0.0);
	DVectorN<DPlanarQuadric::N> mb(0.0);
	DPoint2d pt;
	double z;
	for(size_t i = 0; i < pct; i++){
		bool ok = squadric.base_surface->getParametersAndSignedDistance( points[i], pt, z );
		//pt = squadric.base_surface->getParameters(points[i]);
		//z = points[i].distance( squadric.base_surface->getPoint( pt ) );
		const double v[DPlanarQuadric::N] = { 1, pt.x, pt.y, 
			pt.x * pt.x, pt.y * pt.y, pt.x * pt.y };
		for(int j = 0; j < DPlanarQuadric::N; j++){
			for(int k = 0; k < DPlanarQuadric::N; k++)
				mA(j,k) += v[j] * v[k];
			mb[j] += v[j] * z;
		}

#ifdef FQOS_SHOW
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, i << "\t" << pt.x << "\t" << pt.y << "\t" << z 
			<< "\t" << points[i].x << "\t" << points[i].y << "\t" << points[i].z << endl;
		params.add( pt );
		set->addLabel( points[i], to_string( z ) + " ["+ to_string(pt.x) + "," + to_string(pt.y) + "]" );
#endif
	}

	if(! mA.solveLU(mb, squadric.vq)) return LS_FIT_ERROR;

#ifdef FQOS_SHOW
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "============================");
	for(int PI = 0; PI < pct; PI++){
		DPoint2d pt2d = params[PI];
		const double v[squadric.N] = { 1,  pt2d.x, pt2d.y, pt2d.x * pt2d.x, pt2d.y * pt2d.y, pt2d.x * pt2d.y };
		double dist = 0.0;
		for(int si = 0; si < squadric.N; si++) dist += v[si] * squadric.vq[si];		
		DPoint3d base_pt = squadric.base_surface->getPoint(pt2d);
		DPoint3d quad_pt = squadric.getPoint(pt2d);
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, PI << "\t" << pt2d.x << "\t" << pt2d.y << "\t" << dist 
			<< "\t" << base_pt.x << "\t" << base_pt.y << "\t" << base_pt.z
			<< "\t" << quad_pt.x << "\t" << quad_pt.y << "\t" << quad_pt.z
			<< "\t" << base_pt.distance(quad_pt) << endl;
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "============================");
	surface->createViewSetForPoints(set, params);
	DQuadricOnSurface::createViewSetForPoints(set, points, 
		surface, squadric.vq, 3 );
	//m_base_surface->createViewSetForPoints(set, params);						
	for(int i = 0; i < pct; i++){
		DVector3d normal = squadric.base_surface->getNormalVector( params[i] );
		DPoint3d spt = squadric.base_surface->getPoint( params[i] );
		//set->addEdge( spt, spt + normal * 0.04, 2 );
		set->addEdge( spt, squadric.getPoint( params[i] ), 1 );
	}
	SHOW_MESH("DLSF::fitQuadricOnSurface", set);
#endif

	// -- minimum of energy function (compared to the diameter of the point cloud)
	// LOG4CPLUS_INFO(MeshLog::logger_mesh, "min eigenvalue (energy function) => " << mA_values[min_i]);
	double dist, max_dist = 0.0;
	for(size_t i = 0; i < pct; i++){
		dist = squadric.distance(points[i]);
		if(dist > max_dist) max_dist = dist;
	}

#ifdef _DEBUG
	if(  pct < 10 ){
		MeshViewSet* set = new MeshViewSet;
		squadric.createViewSetForPoints(set, points);
		set->addInfo("pct", (int)pct);
		set->addInfo("max_dist", max_dist);
		SHOW_MESH( " LSF::quadric-on-surface with few points", set );
	}
#endif

	return max_dist;
}
/// Least Squares quadric fitting of points, returns max distance
double DLeastSquaresFitting::fitPlanarQuadric(const DataVector<DPoint3d> & points, DPlanarQuadric& pquadric)
{
	DPlane plane;
	double dist = fitHyperplaneOrthogonal(points, plane);
	if(dist < 0.0) return dist; // error fitting

	return fitPlanarQuadric(points, plane, pquadric);
}

/// Least Squares quadric fitting of points for the given line, returns max distance
double DLeastSquaresFitting::fitLinearQuadric( const DataVector<DPoint2d> & points, 
			const DLine2d& line, DLinearQuadric& lquadric, bool ret_max_dist )
{
	lquadric.line = line;

	// M(A) -- second form of the energy function
	DMatrixN<DLinearQuadric::N> mA(0.0);
	DVectorN<DLinearQuadric::N> mb(0.0);
	size_t pct = points.countInt();
	double t;
	double z;
	for(size_t i = 0; i < pct; i++){
		if(!lquadric.line.projectToLine(points[i], t, z)) 
			continue;
		const double v[DLinearQuadric::N] = { 1, t, t*t };
		for(int j = 0; j < DLinearQuadric::N; j++){
			for(int k = 0; k < DLinearQuadric::N; k++)
				mA(j,k) += v[j] * v[k];
			mb[j] += v[j] * z;
		}
	}

	if(! mA.solveLU(mb, lquadric.vq)) return LS_FIT_ERROR;

	if( ret_max_dist ) {
		double dist, max_dist = 0.0;
		for(size_t i = 0; i < pct; i++){
			dist = lquadric.distance(points[i]);
			if(dist > max_dist) max_dist = dist;
		}
		return max_dist;
	}else{
		// -- minimum of energy function (compared to the diameter of the point cloud)
		// mA_values[min_i];
		return 0;
	}
}

/// Least Squares quadric fitting of points, returns max distance
double DLeastSquaresFitting::fitLinearQuadric(const DataVector<DPoint2d> & points, 
			DLinearQuadric& lquadric, bool ret_max_dist )
{
	DLine2d line;
	double dist = fitLine( points, line, ret_max_dist );
	if(dist < 0.0) return dist; // error fitting

	return fitLinearQuadric(points, line, lquadric, ret_max_dist );
}

/// Least Squares quadric fitting of points for the given line, returns max distance
double DLeastSquaresFitting::fitLinearQuadric( const DataVector<DPoint3d> & points, 
			const DPlane& plane, DLinearQuadric3d& lquadric, bool ret_max_dist )
{
	size_t pct = points.countInt();
	DataVector<DPoint2d> spoints( pct );
	for(size_t i = 0; i < pct; i++)
		spoints.add( plane.projectToPlane( points[i] ) );

	DLine2d sline;
	double dist = fitLine( spoints, sline, false );
	if(dist < 0.0) return dist; // error fitting
	
	DPoint3d ppt1 = plane.projectToSpace( sline.m_pt );
	DPoint3d ppt2 = plane.projectToSpace( sline.m_pt + sline.m_vt );
	DPoint3d pptn = plane.projectToSpace( sline.m_pt + sline.m_vn );

	DPlane plane1( ppt1, (ppt2-ppt1).normalized(), (pptn-ppt1).normalized() );

	lquadric.line.setPtVt( plane.p0, plane1.e0 );
	lquadric.e1 = plane1.e1;

	spoints.clear();
	for(size_t i = 0; i < pct; i++)
		spoints.add( plane1.projectToPlane( points[i] ) );
	sline.m_pt = DPoint2d::zero;
	sline.m_vt = DVector2d::v_ox;
	sline.m_vn = DVector2d::v_oy;

	DLinearQuadric slquadric;
	dist = fitLinearQuadric(spoints, sline, slquadric, false );
	if(dist < 0.0) return dist; // error fitting
	lquadric.vq1 = slquadric.vq;

	DPlane plane2( ppt1, plane1.e0, plane1.vn );
	lquadric.e2 = plane2.e1;

	spoints.clear();
	for(size_t i = 0; i < pct; i++)
		spoints.add( plane2.projectToPlane( points[i] ) );

	dist = fitLinearQuadric(spoints, sline, slquadric, false );
	if(dist < 0.0) return dist; // error fitting
	lquadric.vq2 = slquadric.vq;

	if( ret_max_dist ) {
		double dist, max_dist = 0.0;
		for(size_t i = 0; i < pct; i++){
			dist = lquadric.distance(points[i]);
			if(dist > max_dist) max_dist = dist;
		}
		return max_dist;
	}else{
		return 0;
	}
}

/// Least Squares quadric fitting of points, returns max distance
double DLeastSquaresFitting::fitLinearQuadric(const DataVector<DPoint3d> & points, 
			DLinearQuadric3d& lquadric, bool ret_max_dist )
{
	DPlane plane;
	double dist = fitHyperplaneOrthogonal( points, plane, false, false );
	if( dist < 0.0) return dist; // error fitting

	return fitLinearQuadric( points, plane, lquadric, ret_max_dist );
}

/// Least Squares fitting of grid, returns max distance
double DLeastSquaresFitting::fitRegularGrid(
			const DataVector<DPoint3d> & points, 
			DataMatrix<DPoint3d> & /* grid */)
{
	size_t pct = points.countInt();
/*
	const int N = grid.rows();
	const int M = grid.countInt() / N;
	if(N < 2 || M < 2) return LS_FIT_ERROR; // impossible
	DMatrixNV mA(grid.countInt(), 0.0);

	DPoint3d pt_first = grid(0,0);
	DPoint3d pt_last  = grid(N-1, M-1);
	double dx = (pt_last.x - pt_first.x) / (N-1);
	double dy = (pt_last.y - pt_first.y) / (M-1);
	// -> form M(A)
	for(int i = 0; i < pct; i++){
		const DPoint3d& pt = points[i];
		// -> identify leaf containing this point and its vertices
		double kdx = (pt.x - pt_first.x) / dx;
		double kdy = (pt.y - pt_first.y) / dy;
		int kx = (int)kdx;
		int ky = (int)kdy;
		assert(kx >= 0 && kx < N);
		assert(ky >= 0 && ky < M);
		kdx -= kx;
		kdy -= ky;
		// mA(j,k) += ;
	}

	double min_energy = mA.smallestEigenvectorNormalized(quadric.vq);
	if(min_energy < 0.0) return LS_FIT_ERROR;

	// -> remove zero-rows/columns

	// -> solve

	// -> adjust data for "zero-row" variables
*/
	// -> count maximum distance
	double max_dist = 0.0;
	for(size_t i = 0; i < pct; i++){
		double dist = /* TODO */ 0.0;
		if(dist > max_dist) max_dist = dist;
	}

	return max_dist;
}

/// Least Squares fitting of orthogonal vector, returns max distance
double DLeastSquaresFitting::fitVectorOrthonormal(const DataVector<DVector3d> & vectors, DVector3d& vn)
{
	size_t vct = vectors.countInt();
	DataVector<DPoint3d> points( vct );
	for(size_t i = 0; i < vct; i++) {
		const DVector3d& v = vectors[i];
		points.add( DPoint3d(v.x, v.y, v.z) );
	}
	DPlane plane;
	double res = fitHyperplaneOrthogonal( points, plane, false, false );
	vn = plane.vn;
	return res;
}

/// Least Squares linear fitting of points using orthogonal regression, returns max distance
double DLeastSquaresFitting::fitLine(const DataVector<DPoint2d> & points, DLine2d& line, bool ret_max_dist) 
{
	assert( points.countInt() > 1 );
	// square root approximation
	DPoint2d pt_mid = DPoint2d::zero;
	size_t pct = points.countInt();
	for(size_t k = 0; k < pct; k++) pt_mid.add(points[k]);
	pt_mid /= (double)pct;

	// M(A) -- second form of the energy function
	DMatrix2d mA;
	for(size_t k = 0; k < pct; k++){
		const DVector2d dv = points[k] - pt_mid;
		mA.m[0][0] += dv.x * dv.x;
		mA.m[0][1] += dv.x * dv.y;
		mA.m[1][0] += dv.y * dv.x;
		mA.m[1][1] += dv.y * dv.y;
	}
	DMatrix2d mA_vectors;
	double mA_values[2];
	mA.eigensystem(mA_vectors, mA_values);

	// -- the smallest eigenvalue of M(A)
	int min_i = ( abs(mA_values[0]) < abs(mA_values[1])) ? 0 : 1;

	// -- solution
	line.setPtVn( pt_mid, mA_vectors.column(min_i).normalized() );

	// -- minimum of energy function (compared to the diameter of the point cloud)
	//LOG4CPLUS_INFO(MeshLog::logger_console, "min eigenvalue (energy function)", mA_values[min_i]);
	if( ret_max_dist ){
		double c = line.m_vn.x * line.m_pt.x + line.m_vn.y * line.m_pt.y;
		double max_dist = 0.0;
		for(size_t k = 0; k < pct; k++){
			const DPoint2d& pt = points[k];
			double dist = abs(line.m_vn.x * line.m_pt.x + line.m_vn.y * line.m_pt.y - c);
			if(dist > max_dist) max_dist = dist;
		}
		return max_dist;
	}else{
		return mA_values[min_i];
	}
}

/// Least Squares fitting of point to several planes (in form <vn, d>), returns max distance
double DLeastSquaresFitting::fitPointToPlanes(DPoint3d& pt, 
	const DataVector<DVector3d> & vns, const DataVector<double> & ds)
{
	DMatrix3d A;
	DVector3d b;
	for(size_t j = 0; j < vns.countInt(); j++){
		const DVector3d & vn = vns[j];
		A.m[0][0] += vn.x * vn.x;
		A.m[1][1] += vn.y * vn.y;
		A.m[2][2] += vn.z * vn.z;
		A.m[0][1] += vn.x * vn.y;
		A.m[0][2] += vn.x * vn.z;
		A.m[1][2] += vn.y * vn.z;
		double d = ds[j];
		b.x += d * vn.x;
		b.y += d * vn.y;
		b.z += d * vn.z;
	}
	A.m[1][0] = A.m[0][1];
	A.m[2][0] = A.m[0][2];
	A.m[2][1] = A.m[1][2];

	// solve
	const DPoint3d old_pt = pt;

	bool deg_x = abs(b.x) < SMALL_NUMBER && abs(A.m[0][0]) < SMALL_NUMBER && 
		abs(A.m[0][1]) < SMALL_NUMBER && abs(A.m[0][2]) < SMALL_NUMBER;
	bool deg_y = abs(b.y) < SMALL_NUMBER && abs(A.m[1][0]) < SMALL_NUMBER && 
		abs(A.m[1][1]) < SMALL_NUMBER && abs(A.m[1][2]) < SMALL_NUMBER;
	bool deg_z = abs(b.z) < SMALL_NUMBER && abs(A.m[2][0]) < SMALL_NUMBER && 
		abs(A.m[2][1]) < SMALL_NUMBER && abs(A.m[2][2]) < SMALL_NUMBER;

	if(deg_x && deg_y && deg_z){
		// NOP
	}else if(deg_x && deg_y){
		// count z
		pt = DPoint3d(old_pt.x, old_pt.y, b.z / A.m[2][2]);
	}else if(deg_x && deg_z){
		// count y
		pt = DPoint3d(old_pt.x, b.y / A.m[1][1], old_pt.z);
	}else if(deg_y && deg_z){
		// count x
		pt = DPoint3d(b.x / A.m[0][0], old_pt.y, old_pt.z);
	}else if(deg_x){
		// count y and z
		DMatrix2d A2(A.m[1][1], A.m[1][2], A.m[2][1], A.m[2][2]);
		DVector2d b2(b.y, b.z), res;
		if(A2.solve(b2, res)){
			pt = DPoint3d(old_pt.x, res.x, res.y);
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Error solving (reduced for x) optimization system for node");
		}
	}else if(deg_y){
		// count x and z
		DMatrix2d A2(A.m[0][0], A.m[0][2], A.m[2][0], A.m[2][2]);
		DVector2d b2(b.x, b.z), res;
		if(A2.solve(b2, res)){
			pt = DPoint3d(res.x, old_pt.y, res.y);
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Error solving (reduced for y) optimization system for node");
		}
	}else if(deg_z){
		// count x and y
		DMatrix2d A2(A.m[0][0], A.m[0][1], A.m[1][0], A.m[1][1]);
		DVector2d b2(b.x, b.y), res;
		if(A2.solve(b2, res)){
			pt = DPoint3d(res.x, res.y, old_pt.z);
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Error solving (reduced for z) optimization system for node");
		}
	}else{
		// else - count all
		DVector3d res;
		if(A.solve(b, res)){
			pt = DPoint3d::zero + res;
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Error solving optimization system for node");
		}
	}

	// stat
	double max_dist = 0.0;
	for(size_t i = 0; i < vns.countInt(); i++){
		double dist = vns[i].x * old_pt.x + vns[i].y * old_pt.y + vns[i].z * old_pt.z - ds[i];
		if(dist > max_dist) max_dist = dist;
	}
	return max_dist;
	//return pt.distance(old_pt);
}

/// Least Squares fitting of points to cylinder, using projection on plane, returns max distance
double DLeastSquaresFitting::fitCylinder(
			const DataVector<DPoint3d> & points,
			const DPlane& plane,
			DCylinder& cylinder)
{
	const DVector3d pl_e0 = plane.e0.normalized();
	const DVector3d pl_e1 = plane.e1.normalized();
	const DVector3d& pl_n = plane.vn;
	assert(abs(1.0 - pl_e0.length2()) < SMALL_NUMBER);	// should be normalized
	assert(abs(1.0 - pl_e1.length2()) < SMALL_NUMBER);
	assert(abs(1.0 - pl_n.length2()) < SMALL_NUMBER);

	// check for two orthogonal planes
	DPlane planeA(plane.p0, pl_e0, pl_n);
	DPlane planeB(plane.p0, pl_e1, pl_n);
	// project points
	size_t pct = points.countInt();
	DataVector<DPoint2d> pointsA(pct);
	DataVector<DPoint2d> pointsB(pct);
	for(size_t i = 0; i < pct; i++){
		pointsA.add(planeA.projectToPlane(points.get(i)));
		pointsB.add(planeB.projectToPlane(points.get(i)));
	}
	DCircle ca, cb;
	double max_dist_2d_A = fitCircle(pointsA, ca);
	double max_dist_2d_B = fitCircle(pointsB, cb);

	// check real distance for 3d points / cylinder
	DCylinder cylinderA(planeA.projectToSpace(ca.m_center), planeA.vn, ca.m_radius);
	DCylinder cylinderB(planeB.projectToSpace(cb.m_center), planeB.vn, cb.m_radius);
	double max_dist2_A = 0.0;
	double max_dist2_B = 0.0;
	for(size_t i = 0; i < pct; i++){
		const DPoint3d& pt = points.get(i);
		double dist2 = pt.distance2(cylinderA.getPoint(cylinderA.getParam(pt)));
		if(max_dist2_A < dist2) max_dist2_A = dist2;
		dist2 = pt.distance2(cylinderB.getPoint(cylinderB.getParam(pt)));
		if(max_dist2_B < dist2) max_dist2_B = dist2;
	}

	if(max_dist2_A < max_dist2_B){
		cylinder = cylinderA;
		return sqrt(max_dist2_A);
	}else{
		cylinder = cylinderB;
		return sqrt(max_dist2_B);
	}
}

/// Least Squares fitting of points to cylinder, using projection on plane, returns max distance
double DLeastSquaresFitting::fitCylinder(
			const DataVector<DPoint3d> & points,
			const DVector3d& caxis,
			DCylinder& cylinder)
{
	DPlane cplane(DPoint3d::zero, caxis.normalized());
	// project points
	size_t pct = points.countInt();
	DataVector<DPoint2d> cplane_points(pct);
	for(size_t i = 0; i < pct; i++)
		cplane_points.add(cplane.projectToPlane(points[i]));
	DCircle cyl_circle;
	double max_dist_2d_A = fitCircle(cplane_points, cyl_circle);

	// check real distance for 3d points / cylinder
	cylinder.set( cplane.projectToSpace(cyl_circle.m_center), cplane.vn, cyl_circle.m_radius);
	double max_dist = 0.0;
	for(size_t i = 0; i < pct; i++){
		const DPoint3d& pt = points[i];
		double dist = cylinder.distanceFromSurface(pt);
		if(max_dist < dist) max_dist = dist;
	}
	return max_dist;
}

/// Least Squares fitting of points to circle, iterative, returns max distance
double DLeastSquaresFitting::fitCircleBullock(const DataVector<DPoint2d> & points, 
			DCircle& circle, bool ret_max_dist)
{
	size_t pct = points.countInt();
	assert(pct > 2);
	double fpct = 1.0 / pct;

	// average point
	DPoint2d pt_ave = DPoint2d::zero;
	for(size_t i = 0; i < pct; i++)
		pt_ave.add(points[i]);
	pt_ave *= fpct;

	double Suu = 0.0;
	double Suv = 0.0;
	double Svv = 0.0;
	double Suuu = 0.0;
	double Suuv = 0.0;
	double Suvv = 0.0;
	double Svvv = 0.0;
	for(size_t i = 0; i < pct; i++){
		DVector2d dv = points[i] - pt_ave;
		double uu = dv.x * dv.x;
		double uv = dv.x * dv.y;
		double vv = dv.y * dv.y;
		Suu += uu;
		Suv += uv;
		Svv += vv;
		Suuu += uu * dv.x;
		Suuv += uu * dv.y;
		Suvv += vv * dv.x;
		Svvv += vv * dv.y;
	}

	DMatrix2d mA( Suu, Suv, Suv, Svv);
	DVector2d mb( 0.5 * ( Suuu + Suvv ), 0.5 * ( Suuv + Svvv ) ), mx;

	if( !mA.solve(mb, mx) ) return LS_FIT_ERROR;

	circle.m_center = pt_ave + mx;
	circle.m_radius = sqrt( mx.length2() + (Suu+Svv)*fpct );

	if( ret_max_dist ){
		double max_dist = 0;
		// calculate max
		for(size_t i = 0; i < pct; i++){
			double dist = abs(circle.m_center.distance(points[i]) - circle.m_radius);
			if(dist > max_dist) max_dist = dist;
		}

		return max_dist;
	}else{
		return 0;
	}
}

/// Least Squares fitting of points to sphere, iterative, returns max distance
double DLeastSquaresFitting::fitSphereIter(const DataVector<DPoint3d> & points, 
			DSphere& sphere, bool calculate_initial, bool ret_max_dist)
{
	size_t pct = points.countInt();
	assert(pct > 3);
	double fpct = 1.0 / pct;
	double fpct2 = fpct * fpct;

	// average point
	DPoint3d pt_ave = DPoint3d::zero;
	for(size_t i = 0; i < pct; i++)
		pt_ave.add(points[i]);
	pt_ave *= fpct;

	// iterate
	if( calculate_initial )
		sphere.m_center = pt_ave;

	for(size_t k = 0; k < 20; k++){
		DPoint3d p0 = sphere.m_center;
		double L = 0.0;
		double La = 0.0;
		double Lb = 0.0;
		double Lc = 0.0;
		for(size_t i = 0; i < pct; i++){
			const DPoint3d& pt = points[i];
			double dx = p0.x - pt.x;
			double dy = p0.y - pt.y;
			double dz = p0.z - pt.z;
			double Li = sqrt(dx*dx + dy*dy + dz*dz);
			L += Li;
			La += dx / Li;
			Lb += dy / Li;
			Lc += dz / Li;
		}
		sphere.m_center = pt_ave + DVector3d(La, Lb, Lc) * (L * fpct2);
		sphere.m_radius = L * fpct;

		if( p0.distance2(sphere.m_center) < VERY_SMALL_NUMBER ) break;
	}
	if( ret_max_dist ){
		double max_dist = 0;
		// calculate max
		for(size_t i = 0; i < pct; i++){
			double dist = abs(sphere.m_center.distance(points[i]) - sphere.m_radius);
			if(dist > max_dist) max_dist = dist;
		}

		return max_dist;
	}else{
		return 0;
	}
}

/// Least Squares fitting of points to circle, iterative, returns max distance
double DLeastSquaresFitting::fitCircleIter(const DataVector<DPoint2d> & points, 
			DCircle& circle, bool ret_max_dist)
{
	size_t pct = points.countInt();
	assert(pct > 2);
	double fpct = 1.0 / pct;
	double fpct2 = fpct * fpct;

	// average point
	DPoint2d pt_ave = DPoint2d::zero;
	for(size_t i = 0; i < pct; i++)
		pt_ave.add(points[i]);
	pt_ave *= fpct;

	// iterate
	circle.m_center = pt_ave;
	for(size_t k = 0; k < 20; k++){
		DPoint2d p0 = circle.m_center;
		double L = 0.0;
		double La = 0.0;
		double Lb = 0.0;
		for(size_t i = 0; i < pct; i++){
			const DPoint2d& pt = points[i];
			double dx = p0.x - pt.x;
			double dy = p0.y - pt.y;
			double Li = sqrt(dx*dx + dy*dy);
			L += Li;
			La += dx / Li;
			Lb += dy / Li;
		}
		circle.m_center.x = pt_ave.x + L * La * fpct2;
		circle.m_center.y = pt_ave.y + L * Lb * fpct2;
		circle.m_radius = L * fpct;

		if( p0.distance2(circle.m_center) < VERY_SMALL_NUMBER ) break;
	}
	if( ret_max_dist ){
		double max_dist = 0;
		// calculate max
		for(size_t i = 0; i < pct; i++){
			double dist = abs(circle.m_center.distance(points[i]) - circle.m_radius);
			if(dist > max_dist) max_dist = dist;
		}

		return max_dist;
	}else{
		return 0;
	}
}

/*
1. Function to minimize is (sum of the area differences):
Q = Sum { [ (xi - X0)^2 + (yi -Y0)^2 - R^2 ]^2 , i=1..N }

2. A necessary condition is the system of 3 equations with 3 unknowns
X0, Y0, R. Calculate the partial derivatives of Q, with respect to 
X0, Y0, R. (all d's are partial)

dQ / dX0 = 0
dQ / dY0 = 0
dQ / dR = 0

3. DeveloPIng we get the linear least-squares problem:

| x1 y1 1 | | a |    | -x1^2-y1^2 |
| x2 y2 1 | | b | =~ | -x2^2-y2^2 |
| x3 y3 1 | | c |    | -x3^2-y3^2 |
| x4 y4 1 |          | -x4^2-y4^2 |
| x5 y5 1 |          | -x5^2-y5^2 |
  .....               .........
(for example, for 5 points)

where  a = 2 X0; b = 2 Y0 ; c = X0^2 + Y0^2 - R^2. 
Take any good least-squares algorithm to solve it, yielding a,b,c.  
So the final circle solution will be given with
X0 = a/2; Y0 = b/2; R^2 = X0^2+Y0^2 - c. 
*/
/// Least Squares fitting of points to circle, iterative, returns max distance
double DLeastSquaresFitting::fitCircleMinArea(const DataVector<DPoint2d> & points, 
			DCircle& circle, bool ret_max_dist)
{
	size_t pct = points.countInt();
	assert(pct > 2);

	// M(A) -- second form of the energy function
	DMatrix3d mA;
	DVector3d mb,mx;
	for(size_t i = 0; i < pct; i++){
		const DPoint2d& pt = points[i];
		double d = pt.x*pt.x + pt.y*pt.y;
		const double v[3] = { pt.x, pt.y, 1 };
		for(int j = 0; j < 3; j++){
			for(int k = 0; k < 3; k++)
				mA(j,k) += v[j] * v[k];
			mb[j] -= v[j]*d;
		}
	}

	if( mA.solve(mb, mx ) ) {
		// -- solution
		circle.m_center.x = -0.5*mx.x;
		circle.m_center.y = -0.5*mx.y;
		circle.m_radius = sqrt( sqr(circle.m_center.x)+sqr(circle.m_center.y) - mx.z); 
	}else{
		return LS_FIT_ERROR;
	}

	// -- minimum of energy function (compared to the diameter of the point cloud)
	// LOG4CPLUS_INFO(MeshLog::logger_mesh, "min eigenvalue (energy function) => " << mA_values[min_i]);
	// -- the smallest eigenvalue of M(A)

	if( ret_max_dist ){
		double max_dist = 0;
		// calculate max
		for(size_t i = 0; i < pct; i++){
			double dist = abs(circle.m_center.distance(points[i]) - circle.m_radius);
			if(dist > max_dist) max_dist = dist;
		}
		return max_dist;
	}else{
		return 0;
	}
}

/// test fitting
void DLeastSquaresFitting::testFit()
{
	DCircle circle1, circle2, circle3;
	static const int N = 100;
	RandomGen rg(1);
	DataVector<DPoint2d> points(100);

	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"=============================================");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"res1\tres2\tres3\tres_min\tres1/min\tres2/min\tres3/min"
		<< "\tr\tr1\tr2\tr3\tdc1\tdc2\tdc3");

	for(int i = 0; i < N; i++ ) {
		double r = rg.doub( 10, 100 );
		DPoint2d center( rg.doub(-100, 100 ), rg.doub(-100, 100 ) );

		unsigned int np = rg.int32(10, 100);
		points.clear();
		for(unsigned int j = 0; j < np; j++) {
			double a = rg.doub( 0, 2*PI );
			points.add( center + DVector2d( r*sin(a), r*cos(a) ) );
		}

		double res1 = fitCircleIter(points, circle1, true);
		double res2 = fitCircleMinArea(points, circle2, true);
		double res3 = fitCircleBullock(points, circle3, true);

		double res_min = std::min( std::min( res1, res2), res3);

		LOG4CPLUS_INFO(MeshLog::logger_mesh, 
			res1 << "\t" << res2 << "\t" << res3 
			<< "\t" << res_min << "\t" << res1/res_min << "\t" << res2/res_min << "\t" << res3/res_min
			<< "\t"	<< r << "\t" << circle1.m_radius << "\t" << circle2.m_radius << "\t" << circle3.m_radius 
			<< "\t" << center.distance(circle1.m_center)
			<< "\t" << center.distance(circle2.m_center)
			<< "\t" << center.distance(circle3.m_center));

	}

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "=============================================");
}
