// SurfaceParametric.cpp: implementation of the SurfaceParametric class.
//
//////////////////////////////////////////////////////////////////////

#include "SurfaceParametric.h"

#include "common.h"
#include "DEquation.h"
#include "MeshData.h"
#include "Curve2dParametric.h"
#include "DMetric3d.h"
#include "Curve2dSegment.h"
#include "DMatrix.h"
#include "Metric3dContext.h"
#include "ControlSpace2dAdaptive.h"
#include "ControlSpace3dAdaptive.h"
#include "MeshViewSet.h"
#include "DLeastSquaresFitting.h"
#include "DPlane.h"
#include "DPlanarQuadric.h"
#include "SurfacePlane.h"
#include "SurfacePlanarQuadric.h"
#include "SurfaceBSplinePlanar.h"
#include "SurfaceDomain.h"
#include "SurfaceCylinder.h"
#include "SurfaceSphere.h"
#include "MeshGenerator3dSurface.h"
#include "MeshPoint3d.h"
#include "MeshFace.h"
#include "DataStatistics.h"

#define HEAP_SIZE 50

SurfaceParametric::~SurfaceParametric()
{
	if(m_domain) delete m_domain;
}

/// domain check
double SurfaceParametric::withinDomainQuality( const DPoint3d& pt, DPoint2d & param ) const
{
	if( m_domain != nullptr ) { // no domain -> infinite
		if( !m_domain->isInsideOBBox( pt ) ) 
			return AQ_INVALID;
	}
	param = getParametersNear( pt, param );

	if( !withinParamRange( param ) ) return AQ_INVALID;

	return (m_domain == nullptr) ? AQ_VALID_MAX : m_domain->getInsideQuality( pt, param );
}

/// return domain middle param
DPoint2d SurfaceParametric::getDomainMiddleParam() const
{
	if( m_domain == nullptr ) return DPoint2d::zero;
	else return m_domain->getMiddleParam();
}

/*
bool SurfaceParametric::withinDomain(Metric3dContext& mc, MeshPoint3d * point, DPoint2d * param, double * quality) const
{
	SurfaceParametric* psurf = point->getLocalSurface();
	const double& pquality = point->getLocalSurfaceQuality();
	if( psurf == this && pquality > 0.0 ) {
		if(param) *param = point->getLocalSurfaceParam();
		if(quality) *quality = pquality;
		return true;
	}

	if(m_domain){
		const DPoint3d& pt_3d = point->getCoordinates();
		if( ! m_domain->isInsideOBBox( pt_3d ) ) { 
			assert( psurf != this ); 
			return false; 
		}
		DPoint2d pt_2d = ( psurf == this ) ? point->getLocalSurfaceParam() : getParameters(pt_3d);
		if( param ) *param = pt_2d;
		bool within = m_domain->isInside(pt_3d, pt_2d, quality);
		if( within && quality ) { 
			if( SurfaceDomain::w_surf_dist >= 0.0 ) {
				mc.countMetricAtPoint( pt_3d, true );
				double dist_q = mc.transformRStoMS(pt_3d - getPoint( pt_2d ) ).length(); // / MeshGenerator3dSurface::param_local_shape_tolerance;
				if( dist_q > 1.0 ) {
					//static int local_counter = 0;
					//local_counter++;
					if(false){
						MeshViewSet* set = new MeshViewSet;
						m_domain->draw(set, this);
						set->addPoint( pt_3d, 1 );
						set->addLabel( getPoint( pt_2d ), "surf-pt" );
						set->addInfo("isInside", within? "true" : "false");
						set->addInfo("dist_q", dist_q);
						SHOW_MESH( "SP::withinDomain failed", set );
					}
					return false;
				}else if ( SurfaceDomain::w_surf_dist > 0.0 )
					*quality += SurfaceDomain::w_surf_dist * (1.0 - dist_q );
			}
			if( !point->isBorder() && SurfaceDomain::w_normal_sp >= 0.0 ) {
				double sp = point->getBaseNormal().scalarProduct( getNormalVector( pt_2d ));
				if( sp <= 0.0 ) {
					if(false){
						MeshViewSet* set = new MeshViewSet;
						m_domain->draw(set, this);
						set->addPoint( pt_3d, 1 );
						set->addInfo("isInside", within? "true" : "false");
						set->addInfo("normal_sp", sp);
						SHOW_MESH( "SP::withinDomain failed", set );
					}
					return false;
				}else if( SurfaceDomain::w_normal_sp > 0.0 )
					*quality += SurfaceDomain::w_normal_sp * sp;
			}
		}
		return within;
	}else{
		if(quality != nullptr) *quality = 1.0;
		if(param) *param = point->getLocalSurfaceParam( this );
		return true; 
	}
}
*/

/// domain check
/*
bool SurfaceParametric::withinDomain(Metric3dContext& mc, const DPoint3d& pt, const DVector3d& pt_normal, DPoint2d& param, double * quality) const 
{ 
	if(m_domain){
		if( ! m_domain->isInsideOBBox( pt ) ) return false;
		param = getParameters(pt);
		bool within = m_domain->isInside(pt, param, quality);
		if( within && quality ) { 
			if( SurfaceDomain::w_surf_dist >= 0.0 ) {
				mc.countMetricAtPoint( pt, true );
				double dist_q = mc.transformRStoMS(pt - getPoint( param ) ).length(); // / MeshGenerator3dSurface::param_local_shape_tolerance;
				if( dist_q > 1.0 ) {
					//static int local_counter = 0;
					//local_counter++;
					if(false){
						MeshViewSet* set = new MeshViewSet;
						m_domain->draw(set, this);
						set->addPoint( pt, 1 );
						set->addLabel( getPoint( param ), "surf-pt" );
						set->addInfo("isInside", within? "true" : "false");
						set->addInfo("dist_q", dist_q);
						SHOW_MESH( "SP::withinDomain failed", set );
					}
					return false;
				}else if ( SurfaceDomain::w_surf_dist > 0.0 )
					*quality += SurfaceDomain::w_surf_dist * (1.0 - dist_q );
			}
			if( SurfaceDomain::w_normal_sp >= 0.0 ) {
				double sp = pt_normal.scalarProduct( getNormalVector( param ));
				if( sp <= 0.0 ) {
					if(false){
						MeshViewSet* set = new MeshViewSet;
						m_domain->draw(set, this);
						set->addPoint( pt, 1 );
						set->addInfo("isInside", within? "true" : "false");
						set->addInfo("normal_sp", sp);
						SHOW_MESH( "SP::withinDomain failed", set );
					}
					return false;
				}else if( SurfaceDomain::w_normal_sp > 0.0 )
					*quality += SurfaceDomain::w_normal_sp * sp;
			}
		}
		return within;
	}else{
		param = getParameters(pt);
		return true; 
	}
}
*/
void SurfaceParametric::setDomain(SurfaceDomain* domain) 
{
	if(m_domain) delete m_domain;
	m_domain = domain;
}

/// store domain to XML file
ostream& SurfaceParametric::storeDomainXML(ostream& os, const string& prefix) const
{
	if(m_domain) m_domain->storeXML(os, prefix);
	return os;
}

/// draw/add viewset data for the domain
void SurfaceParametric::drawDomain(MeshViewSet* set, int id) const
{
	if(m_domain && set) m_domain->draw(set, this, id);
}

const DVector3d SurfaceParametric::getNormalVector(const DPoint2d& param) const 
{
	const DVector3d vs = getDerivative(DEquation::deriv_ds, param);
	const DVector3d vt = getDerivative(DEquation::deriv_dt, param);
	return vs.crossProduct(vt).normalized();
}

const DVector3d SurfaceParametric::getNormalVectorForPoint3d(const DPoint3d& pt) const 
{
	return getNormalVector( getParameters(pt) );
}

double SurfaceParametric::getShapeParameters(const DPoint3d& point, 
	Curve2dConstPtr shape, double near_t, double min_t, double max_t) const
{
	const double F_ERR2	= sqr(mesh_data.relative_small_number);
	const int PROBE = 10;
	const int STEPS = 20;
	double best_t = near_t;
	double best_dist = point.distance2(getPoint(shape->getPoint(best_t)));
	double dt = (max_t - min_t) / PROBE;
	double t = min_t;
	for(int i = 0; i < STEPS; i++){
		for(int j = 0; (j <= PROBE) && (t <= max_t); j++){
			double dist = point.distance2(getPoint(shape->getPoint(t)));
			if(dist < best_dist){
				if(dist < F_ERR2) return t;
				best_t = t;
				best_dist = dist;
			}
			t += dt;
		}
		dt /= (0.45*PROBE);
		t = best_t - (PROBE/2 * dt);
		if(t < min_t) t = min_t;
	}

//	DPoint2d param = getParametersNear(point, shape->getPoint(near_t));
//	double best_t = shape->getParameter(param, near_t, min_t, max_t);
	return best_t;
}

double SurfaceParametric::getSegmentParameters(const DPoint3d& point, const DPoint2d& pt0, const DPoint2d& pt1, double near_t, double min_t, double max_t) const
{
	DPoint2d param = getParametersNear(point, DPoint2d(pt0, pt1, near_t));
	return Curve2dSegment(pt0, pt1).getParameterInRange(param, near_t, min_t, max_t);
}

const DPoint2d SurfaceParametric::getParametersNear(const DPoint3d& point, const DPoint2d& near_point) const
{
	DPoint2d t = near_point;

	double	F_ERR2		= sqr(mesh_data.relative_small_number);	// Maksymalny b³¹d wyniku
	double	T_ERR2		= sqr(1e-5 * mesh_data.relative_small_number);	// Maksymalny b³¹d parametru
	int		MAX_STEPS	= 30;	// Maksymalna iloœæ kroków

	DPoint2d last_t = t;
	double last_diff = (getPoint(t) - point).length2();
	for(int k = 0; k < MAX_STEPS; k++){
		const DVector3d fd = getPoint(t) - point;
		double f_diff = fd.length2();
		if(f_diff < F_ERR2) return t;
		if(k > 0 && f_diff > last_diff) return last_t;
		last_diff = f_diff;
		last_t = t;
		DVector3d dfu = getDerivative(DEquation::deriv_ds, t);
		DVector3d dfv = getDerivative(DEquation::deriv_dt, t);
		DVector3d dfn = dfu.crossProduct(dfv);
		DMatrix3d A(dfu, dfv, dfn);
		DVector3d w;
		bool res = A.solve(fd, w);
		assert(res);
		if(!res) return t;
		t.x -= w.x;
		t.y -= w.y;
		if(w.x * w.x + w.y * w.y < T_ERR2) return t;
	}
	return t;
}

const SurfaceCurvature SurfaceParametric::getCurvature(const DPoint2d &pt, double & g_ratio) const
{
	const DVector3d fs  = getDerivative(DEquation::deriv_ds, pt);
	const DVector3d ft  = getDerivative(DEquation::deriv_dt, pt);
	const DVector3d fss = getDerivative(DEquation::deriv_dss, pt);
	const DVector3d fst = getDerivative(DEquation::deriv_dst, pt);
	const DVector3d ftt = getDerivative(DEquation::deriv_dtt, pt);
	const DVector3d fn  = fs.crossProduct(ft).normalized();

	double g11 = fs.scalarProduct(fs);
	double g12 = fs.scalarProduct(ft);
	double g22 = ft.scalarProduct(ft);

	g_ratio = std::min(g11, g22) / std::max(g11, g22);
		
	double b11 = fn.scalarProduct(fss);
	double b12 = fn.scalarProduct(fst);
	double b22 = fn.scalarProduct(ftt);

	double A = g11*g22 - sqr(g12);
	double B = 2*b12*g12 - b22*g11 - b11*g22;
	double C = b11*b22 - sqr(b12);

	double kx = std::max(std::max(abs(A), abs(B)), abs(C));
	if(kx == 0.0) return SurfaceCurvature(false);

	A /= kx;
	B /= kx;
	C /= kx;

	if(abs(A) > SMALL_NUMBER){
		double delta = B*B - 4.0*A*C;
		if(delta > SMALL_NUMBER){
			delta = sqrt(delta);
			double x1 = (-B - delta) / (2.0*A);
			double x2 = (-B + delta) / (2.0*A);

			// eigenvector
			double c11 = b11-x1*g11;
			double c12 = b12-x1*g12;
			double c22 = b22-x1*g22;
			if(abs(c12) > SMALL_NUMBER)
				return SurfaceCurvature(x1, x2, atan(-c11 / c12));
			else if(abs(c22) > SMALL_NUMBER)
				return SurfaceCurvature(x1, x2, atan(-c12 / c22));
			else
				return SurfaceCurvature(x2, x1, 0.0);
		}else if(delta > -SMALL_NUMBER){
			double x1 = -B / (2.0*A);
			return SurfaceCurvature(x1, x1, 0.0);
		}
	}else if(abs(B) > SMALL_NUMBER){
		double x1 = - C / B;
		return SurfaceCurvature(x1, x1, 0.0);
	}

	return SurfaceCurvature( false );
}

const DVector3d SurfaceParametric::getCurvatureDirection(const DPoint2d &pt, double & g_ratio) const
{
	const DVector3d fs  = getDerivative(DEquation::deriv_ds, pt);
	const DVector3d ft  = getDerivative(DEquation::deriv_dt, pt);
	const DVector3d fss = getDerivative(DEquation::deriv_dss, pt);
	const DVector3d fst = getDerivative(DEquation::deriv_dst, pt);
	const DVector3d ftt = getDerivative(DEquation::deriv_dtt, pt);
	const DVector3d fn  = fs.crossProduct(ft).normalized();

	double g11 = fs.scalarProduct(fs);
	double g12 = fs.scalarProduct(ft);
	double g22 = ft.scalarProduct(ft);

	g_ratio = std::min(g11, g22) / std::max(g11, g22);
		
	double b11 = fn.scalarProduct(fss);
	double b12 = fn.scalarProduct(fst);
	double b22 = fn.scalarProduct(ftt);

	double A = g11*g22 - sqr(g12);
	double B = 2*b12*g12 - b22*g11 - b11*g22;
	double C = b11*b22 - sqr(b12);

	double kx = std::max(std::max(abs(A), abs(B)), abs(C));
	if(kx == 0.0) return fs;

	A /= kx;
	B /= kx;
	C /= kx;

	if(abs(A) > SMALL_NUMBER){
		double delta = B*B - 4.0*A*C;
		if(delta > SMALL_NUMBER){
			delta = sqrt(delta);
			double x1 = (-B - delta) / (2.0*A);
			//double x2 = (-B + delta) / (2.0*A);

			// eigenvector
			double c11 = b11-x1*g11;
			double c12 = b12-x1*g12;
			double c22 = b22-x1*g22;
			if(abs(c12) > SMALL_NUMBER)
				return getPoint( DPoint2d(pt.x + c12, pt.y - c11) ) - getPoint( pt );
			else if(abs(c22) > SMALL_NUMBER)
				return getPoint( DPoint2d(pt.x + c22, pt.y - c12) ) - getPoint( pt );
		}
	}
	return fs;
}

ControlDataMatrix2d SurfaceParametric::countParameterizationMatrix(const DPoint2d& pt, double & g_ratio) const
{
	const DVector3d ds = getDerivative(DEquation::deriv_ds, pt);
	const DVector3d dt = getDerivative(DEquation::deriv_dt, pt);

	double g11 = ds.scalarProduct(ds);
	double g22 = dt.scalarProduct(dt);
	g_ratio = std::min(g11, g22) / std::max(g11, g22);

	return DMetric2d::matrixSquareRoot(
		ControlDataMatrix2d(std::max(g11, MIN_PARAM_GRATIO), ds.scalarProduct(dt), std::max(g22, MIN_PARAM_GRATIO)));
}

double SurfaceParametric::segmentLength(const DPoint2d & param_0, const DPoint2d & param_1) const 
{
	double len = 0.0;
	
	double factor = 1.0 - MEASURE_PRECISION;
	DPoint2d heap[HEAP_SIZE];
	DPoint2d p0 = param_0;
	DPoint2d p1 = param_1;
	DPoint3d pt0 = getPoint(p0);
	DPoint3d pt1 = getPoint(p1);
	heap[0] = p1;
	int top = 0;
	while(true){
		DPoint2d p_middle(p0, p1, 0.5);
		DPoint3d pt_middle = getPoint(p_middle);
		double distance = pt0.distance(pt1);
		double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
		if((distance >= factor * real_distance) || (top == HEAP_SIZE-1)){
			len += real_distance;
			if(top == 0) break;
			p0 = p1;
			pt0 = pt1;
			pt1 = getPoint(p1 = heap[--top]);
		}else{						
			heap[++top] = p1 = p_middle;
			pt1 = pt_middle;
		}
	}	
	return len;
}

ControlDataMatrix2d SurfaceParametric::projectTransformationTensor(
	const DPoint2d& pt, const ControlDataMatrix3d& cdm3d) const
{
	const DVector3d dpu = getDerivative(DEquation::deriv_ds, pt).normalized();
	const DVector3d dpv = getDerivative(DEquation::deriv_dt, pt).normalized();

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "dpu= " << dpu);
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "dpv= " << dpv);

	const ControlDataMatrix3d cdmm3d = cdm3d.inverse().transformationToTensor();

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "cdm3d    = " << cdm3d);
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "cdm3dinv = " << cdm3d.inverse());
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "toTensor = " << cdm3d.inverse().transformationToTensor());

	DVector3d vect_u = cdmm3d * dpu;
	DVector3d vect_v = cdmm3d * dpv;
	ControlDataMatrix2d cdmm(
		vect_u.scalarProduct(dpu),
		vect_u.scalarProduct(dpv),
		vect_v.scalarProduct(dpv));

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "ret= " << cdmm.tensorToTransformation().inverse());

	return cdmm.tensorToTransformation().inverse();
}

/// Returns the approximation of a segment on surfaces via polyline (array of points)
void SurfaceParametric::getPolyLine(DataVector<DPoint3d> & polyline, const DPoint2d& a, const DPoint2d& b) const
{
	const double FACTOR = 1.0 - MEASURE_PRECISION;
	const int MAIN_GRID = 10;

	polyline.add(getPoint(a));

	for(int k = 0; k < MAIN_GRID; k++){
		DPoint2d heap[HEAP_SIZE];
		DPoint2d p0(a, b,   (k)/(double)MAIN_GRID);
		DPoint2d p1(a, b, (k+1)/(double)MAIN_GRID);

		DPoint3d pt0 = getPoint(p0);
		DPoint3d pt1 = getPoint(p1);
		heap[0] = p1;
		int top = 0;
		while(true){
			DPoint2d p_middle(p0, p1, 0.5);
			DPoint3d pt_middle = getPoint(p_middle);
			double distance = pt0.distance(pt1);
			double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
			if((distance >= FACTOR * real_distance) || (top == HEAP_SIZE-1)){
				polyline.add(pt_middle);
				polyline.add(pt1);
				if(top == 0) break;
				p0 = p1;
				pt0 = pt1;
				pt1 = getPoint(p1 = heap[--top]);
			}else{						
				heap[++top] = p1 = p_middle;
				pt1 = pt_middle;
			}
		}
	}	
}

/// Returns the approximation of a segment on surfaces via polyline (array of points)
void SurfaceParametric::getPolyLine(DataVector<DPoint3d> & polyline, 
	Curve2dConstPtr shape, double mt0, double mt1) const
{
	const double FACTOR = 1.0 - MEASURE_PRECISION;
	const int MAIN_GRID = 10;

	polyline.add(getPoint(shape->getPoint(mt0)));

	for(int k = 0; k < MAIN_GRID; k++){
		double heap[HEAP_SIZE];
		double t0 = mt0 +   (k)*(mt1-mt0)/MAIN_GRID;
		double t1 = mt0 + (k+1)*(mt1-mt0)/MAIN_GRID;

		DPoint3d pt0 = getPoint(shape->getPoint(t0));
		DPoint3d pt1 = getPoint(shape->getPoint(t1));
		heap[0] = t1;
		int top = 0;
		while(true){
			double t_middle = (t0+t1)*0.5;
			DPoint3d pt_middle = getPoint(shape->getPoint(t_middle));
			double distance = pt0.distance(pt1);
			double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
			if((distance >= FACTOR * real_distance) || (top == HEAP_SIZE-1)){
				polyline.add(pt_middle);
				polyline.add(pt1);
				if(top == 0) break;
				t0 = t1;
				pt0 = pt1;
				pt1 = getPoint(shape->getPoint(t1 = heap[--top]));
			}else{						
				heap[++top] = t1 = t_middle;
				pt1 = pt_middle;
			}
		}
	}	
}

/// Store XML description to stream
ostream& SurfaceParametric::storeXML(ostream& os, const string& prefix) const
{
	return os << prefix << "<data>not available<data/>\n";
}

/// create sketchy representation of this surface for the area with the given points
MeshViewSet * SurfaceParametric::createViewSetForPoints(MeshViewSet* set, const DataVector<DPoint2d> & points) const
{
	static const int SKETCH_LINES = 20;

	size_t pct = points.countInt();
	DRect box;
	for(size_t i = 0; i < pct; i++){
		set->addPoint( getPoint(points[i]) );
		box.addPoint(points[i]);
	}

	DPoint2d ptx = box.getX0Y0();
	DPoint2d pty = box.getX0Y0();
	DVector2d dx(box.getDX() / (SKETCH_LINES-1), 0.0);
	DVector2d dy(0.0, box.getDY() / (SKETCH_LINES-1));


	for(int i = 0; i < SKETCH_LINES; i++) {
		DPoint2d ptxi = ptx;
		DPoint2d ptyi = pty;
		for(int j = 1; j < SKETCH_LINES; j++) {
			DPoint2d ptxi2 = ptxi + dy;
			if( withinParamRange(ptxi) && withinParamRange(ptxi2) )
				set->addEdge( getPoint(ptxi), getPoint(ptxi2), -1);
			DPoint2d ptyi2 = ptyi + dx;
			if( withinParamRange(ptyi) && withinParamRange(ptyi2) )
				set->addEdge( getPoint(ptyi), getPoint(ptyi2), -1);
			ptxi += dy;
			ptyi += dx;
		}
		ptx += dx;
		pty += dy;
	}

	return set;
}

bool SurfaceParametric::checkSurfaceFit( 
		Metric3dContext& mc, 
		SurfacePtr &surface, 
		const DataVector<DPoint3d> & points, 
		const DataVector<DVector3d> & normals, 
		DataVector<DPoint2d> & params, double tolerance, 
		DataVector<double> * approx_quality,
		bool require_within_tolerance)
{
	DataStatistics sp_stats;
	double metric_max_dist = 0.0;
	int pct = points.countInt();
	bool count_params = params.empty();
	if( count_params ) params.setExactly( pct );
	assert(params.countInt() == pct );

	if(approx_quality) approx_quality->setExactly( pct );

	for(int i = 0; i < pct; i++){
		mc.countMetricAtPoint(points[i]);
		if( count_params ) 
			params[i] = surface->getParameters( points[i] );
		const DPoint2d& pt_2d = params[i];
		double dist = mc.transformRStoMS( points[i] - surface->getPoint( pt_2d ) ).length(); 
		if(require_within_tolerance && (dist > tolerance) ) {
			surface.reset();
			return false;
		}else if (dist > metric_max_dist) metric_max_dist = dist;

		if(approx_quality)
			approx_quality->set(i, std::max(0.0, (1.0 - dist / tolerance) )); // possibly also including normal_sp_quality?

		const DVector3d& dn = normals[i];
		if( !dn.isZero() ) {
			sp_stats.add( dn.scalarProduct( surface->getNormalVector( pt_2d ) ) );
		}
	}

	if( true ) {
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < pct; i++ )
			set->addPoint( points[i] );
		if( approx_quality ) {
			for(int i = 0; i < pct; i++ ) {
				double q = approx_quality->get(i);
				if( q < 0.9 ) 
					set->addLabel( points[i], to_string( q ) );
			}
		}
		set = surface->createViewSetForPoints( set, params );
		set->addInfo("metric_max_dist", metric_max_dist );
		SHOW_MESH("checkSurfaceFit", set);
	}

	if( sp_stats.calculate() ) {
		double min_sp = 0.0;
		if( sp_stats.average() < 0.0 ) {
			surface->invertOrientation();
			min_sp = -sp_stats.maximum();
			// recalculate params
			assert( params.countInt() == points.countInt() );
			for(int i = 0; i < points.countInt(); i++)
				params[i] = surface->getParameters( points[i] );
		}else{
			min_sp = sp_stats.minimum();
		}
		LOG4CPLUS_INFO(MeshLog::logger_mesh,
			"Local surface for " << points.countInt() << "pts (" 
			<< surface->getSimpleDescription() << "), max_dist="
			<< metric_max_dist << ", min_sp=" << min_sp);
	}else{
		LOG4CPLUS_INFO(MeshLog::logger_mesh, 
			"Local surface for " << points.countInt() << "pts ("
			<< surface->getSimpleDescription() << "), max_dist=" 
			<< metric_max_dist << ", min_sp= xxx");
	}

	return true;
}

SurfacePtr SurfaceParametric::fitSurface( Metric3dContext& mc, 
	const DataVector<DPoint3d> & points, 
	const DataVector<DVector3d> & normals,
	DataVector<DPoint2d> & params, double tolerance, MeshFace* base_face, 
	double plane_fit_ratio, DataVector<double> * approx_quality)
{

	// ====> separate function for calculating error
	// ====?    + including normals evaluation....

	// a) plane
	DPlane plane;
	/* double plane_max_dist = */ DLeastSquaresFitting::fitHyperplaneOrthogonal(points, plane);
	SurfacePtr surface(new SurfacePlane(plane));

	params.clear();
	if( checkSurfaceFit( mc, surface, points, normals, params, plane_fit_ratio*tolerance, approx_quality ) )
		return surface;

	// b) planar quadric
	DPlanarQuadric pquadric;
	/* double pquadric_max_dist = */ 
	DLeastSquaresFitting::fitPlanarQuadric(points, plane, pquadric);
	surface = std::make_shared<SurfacePlanarQuadric>(pquadric);

	params.clear();
	if( checkSurfaceFit( mc, surface, points, normals, params, tolerance, approx_quality ) )
		return surface;

	if( base_face != nullptr ) {
		assert( ! base_face->getCurvatureDirection0().isZero() );
		DPlane fplane( base_face->getMiddlePoint(), base_face->getCurvatureDirection0(), 
			base_face->getCurvatureDirection1() );
		/* pquadric_max_dist = */ 
		DLeastSquaresFitting::fitPlanarQuadric(points, fplane, pquadric);
		surface = std::make_shared<SurfacePlanarQuadric>(pquadric);

		params.clear();
		if( checkSurfaceFit( mc, surface, points, normals, params, tolerance, approx_quality ) )
			return surface;
	}

/*
	// c) planar quadric, rotated to reflect the main curvature direction
	SurfacePlanarQuadric spq(pquadric);
	double gratio;
	DVector3d rot_e0 = spq.getCurvatureDirection(DPoint2d::zero, gratio);
	DPlane rot_plane( pquadric.getPoint(DPoint2d::zero),
		rot_e0, spq.getNormalVector(DPoint2d::zero).crossProduct(rot_e0));

	DPlanarQuadric rot_pquadric;
	double rot_pquadric_max_dist = DLeastSquaresFitting::fitPlanarQuadric(points, rot_plane, rot_pquadric);

	max_dist2 = 0.0;
	for(int i = 0; i < points.countInt(); i++){
		if(local_metric) mc.countMetricAtPoint(points[i]);
		double dist2 = mc.transformRStoMS( points[i] - rot_pquadric.getPoint( rot_plane.projectToPlane( points[i] ) ) ).length2(); 
		if(dist2 > max_dist2) max_dist2 = dist2;
	}
	metric_max_dist = sqrt(max_dist2);

	if(metric_max_dist <= tolerance){
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Local surface for " << points.countInt() << "pts (rot planar quadric), max_dist = " << metric_max_dist);
		return new SurfacePlanarQuadric(rot_pquadric);
	}

	// d) planar bsplines
	SurfacePlane splane(rot_plane);
	DRect brect;
	DataVector<DPoint2d> rot_params(points.countInt());
	for(int j = 0; j < points.countInt(); j++){
		rot_params.add(splane.getParameters(points[j])); 
		brect.addPoint(rot_params[j]);
	}
	brect.inflate(0.05);

	SurfaceBSplinePlanar* bsurf2 = new SurfaceBSplinePlanar(splane, brect);
	double bsurf_max_dist2 = bsurf2->fitToPoints(points, points.countInt(), 2);

	max_dist2 = 0.0;
	for(int i = 0; i < points.countInt(); i++){
		if(local_metric) mc.countMetricAtPoint(points[i]);
		double dist2 = mc.transformRStoMS( points[i] - bsurf2->getPoint( rot_params[i] ) ).length2(); 
		if(dist2 > max_dist2) max_dist2 = dist2;
	}
	metric_max_dist = sqrt(max_dist2);

	if(metric_max_dist <= tolerance){
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Local surface for " << points.countInt() << "pts (bspline-planar, 2), max_dist = " << metric_max_dist);
		return bsurf2;
	}else{
		delete bsurf2;
	}

	SurfaceBSplinePlanar* bsurf3 = new SurfaceBSplinePlanar(splane, brect);
	double bsurf_max_dist3 = bsurf3->fitToPoints(points, points.countInt(), 3);

	max_dist2 = 0.0;
	for(int i = 0; i < points.countInt(); i++){
		if(local_metric) mc.countMetricAtPoint(points[i]);
		double dist2 = mc.transformRStoMS( points[i] - bsurf3->getPoint( rot_params[i] ) ).length2(); 
		if(dist2 > max_dist2) max_dist2 = dist2;
	}
	metric_max_dist = sqrt(max_dist2);

	if(metric_max_dist <= tolerance){
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Local surface for " << points.countInt() << "pts (bspline-planar, 3), max_dist = " << metric_max_dist);
		return bsurf3;
	}else{
		delete bsurf3;
	}
*/
	//	DQuadric quadric;
	//	double quadric_max_dist = DLeastSquaresFitting::fitQuadric(local_points, quadric);

	// b) cylinder
	//DPoint3d cyl_pt;
	//DVector3d cyl_e0;
	//double cyl_radius;
	//double cyl_max_dist = DLeastSquaresFitting::fitCylinder(local_points,
	//						plane_pt, plane_e0, plane_e1,
	//						cyl_pt, cyl_e0, cyl_radius);
	//metric_max_dist = cyl_max_dist * metric_f;
	//if(cyl_max_dist <= tolerance){
	//	//surface = new SurfaceCylinder(cyl_pt, cyl_e0, cyl_radius));
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Local surface (cylinder), max_dist = " << metric_max_dist);
	//}

	//}
	// d) cylindrical bsplines
	//if(false){
	//	SurfaceCylinder cylinder(cyl_pt, cyl_e0, cyl_radius);
	//	DRect brect_cyl;
	//	for(int j = 0; j < local_points.countInt(); j++)
	//		brect_cyl.addPoint(cylinder.getParameters(local_points[j]));
	//	brect_cyl.inflate(0.05);

	//	SurfaceBSplineCylindrical* bcsurf2 = new SurfaceBSplineCylindrical(cylinder, brect_cyl);
	//	double bcsurf_max_dist2 = bcsurf2->fitToPoints(local_points, 2);
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Local surface (bsline-cyl, 2), max_dist = " << bcsurf_max_dist2*metric_f);

	//	SurfaceBSplineCylindrical* bcsurf3 = new SurfaceBSplineCylindrical(cylinder, brect_cyl);
	//	double bcsurf_max_dist3 = bcsurf3->fitToPoints(local_points, 3);
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Local surface (bsline-cyl, 3), max_dist = " << bcsurf_max_dist3*metric_f);

	//	delete bcsurf2;
	//	delete bcsurf3;
	//}

	// z) failed ...
	return nullptr;
}


SurfaceParametricSet::SurfaceParametricSet(SurfaceConstPtr surface)
{
	if(surface) m_surfaces.add(surface);
}

SurfaceParametricSet::SurfaceParametricSet(SurfaceSetConstPtr sset, SurfaceConstPtr surface)
{
	if(sset) // copy from other set
		for(int i = 0; i < sset->countInt(); i++)
			m_surfaces.add(sset->getSurface(i));
	if(surface) // and add an additional surface
		m_surfaces.add(surface);
}

int SurfaceParametricSet::addSurface(SurfaceConstPtr surface) {
	return (int)m_surfaces.add(surface);
}

bool SurfaceParametricSet::containsSurface(SurfaceConstPtr surface) const
{
	return m_surfaces.contains( surface );
}

/*
DPoint3d SurfaceParametricSet::fitToNearest(Metric3dContext & mc, const DPoint3d& pt, const DVector3d& vn, SurfaceParametric** nearest_surface) const {
	if(nearest_surface) *nearest_surface = nullptr;
	if(m_surfaces.empty()) return pt;
	DPoint3d best_pt;
	double min_dist2 = -1.0;
	mc.countMetricAtPoint( pt );
	DPoint2d sparam;
	for(int i = 0; i < m_surfaces.countInt(); i++){
		SurfaceParametric* surf = m_surfaces[i];
		if( surf->withinDomain( mc, pt,  vn, sparam ) ) { 
			const DPoint3d spt = surf->getPoint( sparam );
			double dist2 = mc.transformRStoMS(spt - pt).length2();
			if(min_dist2 < 0.0 || dist2 < min_dist2){
				min_dist2 = dist2;
				best_pt = spt;
				if(nearest_surface) *nearest_surface = surf;
			}
		}
	}
	return best_pt;
}

DPoint3d SurfaceParametricSet::fitToAll(Metric3dContext & mc, const DPoint3d& pt, const DVector3d& vn) const {
	if(m_surfaces.empty()) return pt;
	DPoint3d spt = pt;
	DPoint2d sparam;
	for(int i = 0; i < m_surfaces.countInt(); i++){
		if( m_surfaces[i]->withinDomain( mc, spt, vn, sparam ) )
			spt = m_surfaces[i]->getPoint( sparam );
	}
	return spt;
}
*/

int SurfaceParametricSet::count() const { 
	return (int)m_surfaces.countInt(); 
}

bool SurfaceParametricSet::countCurvatureMetric(const DPoint3d& pt, double model_diameter, ControlDataMatrix3d& cdm) const
{
	if(m_surfaces.empty()) return false;

	double max_len = model_diameter * std::min(
		ControlSpace2dAdaptive::param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
	double min_len = std::max(model_diameter * ControlSpace2dAdaptive::param_min_diameter_ratio, 
		ControlSpace2dAdaptive::param_min_length);

	double g_ratio;
	int counter = 0;
	for(int i = 0; i < m_surfaces.countInt(); i++){
		DPoint2d spt2d = m_surfaces[i]->getParameters(pt);
		//DPoint3d spt = m_surfaces[i]->getPoint( spt2d );
		auto curvature = m_surfaces[i]->getCurvature(spt2d, g_ratio);
		if(g_ratio > MIN_PARAM_GRATIO){ 
			// metric 2d
			ControlDataStretch2d data = ControlSpace2dAdaptive::adjustCurvatureData(
				curvature, ControlSpace2dAdaptive::param_curvature_ratio,
				model_diameter, min_len, max_len, 
				ControlSpace2dAdaptive::param_stretch_max_ratio);
			// Convert 2d -> 3d
			DMatrix2d e;
			double d[3];
			bool success = DMetric2d::stretchToMatrix(data).eigensystem(e, d);
			assert(success); if(!success) continue;
			const DVector3d pu = m_surfaces[i]->getDerivative(DEquation::deriv_du, spt2d);
			const DVector3d pv = m_surfaces[i]->getDerivative(DEquation::deriv_dv, spt2d);
			const DVector2d pt_u = e.column(0);
			const DVector2d pt_v = e.column(1);
			const DVector3d e0 = DVector3d(pt_u.x * pu.x + pt_u.y * pv.x,
					pt_u.x * pu.y + pt_u.y * pv.y, pt_u.x * pu.z + pt_u.y * pv.z).normalized();
			const DVector3d e1 = DVector3d(pt_v.x * pu.x + pt_v.y * pv.x,
					pt_v.x * pu.y + pt_v.y * pv.y, pt_v.x * pu.z + pt_v.y * pv.z);
			// after surface transformation -> e0 and e1 are not necessarily orthogonal
			const DVector3d e2 = e0.crossProduct(e1).normalized();
			const DVector3d e1x = e2.crossProduct(e0).normalized(); 
			d[2] = ControlSpace2dAdaptive::param_stretch_max_ratio * std::min(d[0], d[1]);

			if(++counter == 1){
				cdm = ControlDataMatrix3d(e0, e1x, e2, d);
			}else{
				cdm.setMinimum(ControlDataMatrix3d(e0, e1x, e2, d));
			}
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Improper parameterization " << g_ratio);
		}
	}

	return counter > 0;
}

SurfacePtr SurfaceParametric::calculateBestBaseSurface( 
		const DataVector<MeshFace*> & faces, double cmin_eps )
{
	int fct = faces.countInt();

	assert( fct > 0 );
	if( fct == 0 ) return nullptr;

	if( fct == 1 ) { 
		return calculateBestBaseSurface( faces[0], cmin_eps, faces[0] );
	} else {
		int hc_i = -1;
		double hc_max = -1.0;
		for(int i = 0; i < fct; i++ ) {
			double cmax = std::max( faces[i]->getCurvature0(), faces[i]->getCurvature1() );
			if( hc_i == -1 || cmax > hc_max ) {
				hc_max = cmax;
				hc_i = i;
			}
		}
		if( hc_i != -1 ) {
			SurfaceConstPtr surface = calculateBestBaseSurface( faces[hc_i], cmin_eps );
			DataHashTable< MeshPoint3d* > hpoints( 10*fct, nullptr );
			for(int i = 0; i < fct; i++) {
				MeshFace* face = faces[i];
				int fpct = face->getPointCount();
				for(int j = 0; j < fpct; j++)
					hpoints.insert( face->getPoint(j) );
			}
			DataVector< MeshPoint3d* > mpoints( hpoints.countInt() );
			hpoints.getValues( mpoints );
			auto adjusted_surface = surface->adjustedForFaces( faces, mpoints, faces[0] );
			adjusted_surface->fixCenterAndRangeForPeriodic( faces[0]->getMiddlePoint() );
			return adjusted_surface;
		}else
			return nullptr;
		// some averaging or choosing best, or least squares... ??? 
	}
	
	assert( false );
	return nullptr;
}

SurfacePtr SurfaceParametric::calculateBestBaseSurface( 
		MeshFace* face, double cmin_eps, MeshFace* central_face )
{
	// from curvature, or sth...
	double c0 = face->getCurvature0();
	double c1 = face->getCurvature1();

	bool c0smaller = abs(c0) < cmin_eps;
	bool c1smaller = abs(c1) < cmin_eps;

	DPoint3d fmp = face->getMiddlePoint();
	const DVector3d& v0 = face->getCurvatureDirection0();
	const DVector3d& v1 = face->getCurvatureDirection1();

	SurfacePtr  base_surface;

	if(  c0smaller && c1smaller ) { // plane
		if( v0.isZero() || v1.isZero() ) {
			DVector3d vn = face->getBaseNormal();
			assert( ! vn.isZero() );
			base_surface = std::make_shared<SurfacePlane>( fmp, vn );
		}else{
			base_surface = std::make_shared<SurfacePlane>( fmp, v0, v1 );
		}
	/////  -----> or if max/min ratio greater than 1e5 ??????
	}else{
		assert( !v0.isZero() && !v1.isZero() );
		const DVector3d& vn = face->getBaseNormal();

		if( c0smaller ) { // cylinder along direction of c0
			double r_signed = 1.0/c1;
			DPoint3d cyl_mid = fmp - vn * r_signed; // = ...
			base_surface = std::make_shared<SurfaceCylinder>( cyl_mid, v0.normalized(), abs(r_signed) );
		}else if( c1smaller ) { // cylinder along direction of c1
			double r_signed = 1.0/c0;
			DPoint3d cyl_mid = fmp - vn * r_signed; // = ...
			base_surface = std::make_shared<SurfaceCylinder>( cyl_mid, v1.normalized(), abs(r_signed) );
		}else { // sphere / (rotated and stretched)
			double r0_signed = 1.0/c0;
			double r1_signed = 1.0/c1;
			if( sgn( r0_signed ) != sgn( r1_signed ) ) { // saddle surface
				base_surface = std::make_shared<SurfacePlane>( fmp, v0, v1 );
			}else {
				//double r_signed =  (abs(r0_signed) > abs(r1_signed) ) ? r0_signed : r1_signed;
				double r_signed =  0.5 * ( r0_signed + r1_signed );
				DPoint3d sphere_mid = fmp - vn * r_signed;
				base_surface = std::make_shared<SurfaceSphere>( sphere_mid, abs(r_signed) );
			}
		}
	}


	if( central_face != nullptr ) {
		base_surface->setOrientationLikeFace( central_face );
		base_surface->fixCenterAndRangeForPeriodic( central_face->getMiddlePoint() );
	}

	return base_surface;
}

bool SurfaceParametric::setOrientationLikeFace( MeshFace* face )
{
	// check orientation
	DPoint3d fmp = face->getMiddlePoint();
	double sp = face->getBaseNormal().scalarProduct(
		getNormalVector( getParameters( fmp ) ) );

	bool inverted = false;
	if( sp < 0.0 ) {
		inverted = invertOrientation();
		assert(inverted);
		assert( (sp = face->getBaseNormal().scalarProduct( getNormalVector( getParameters( fmp ) ) ) ) > 0.0 );
	}

	return inverted;
}
