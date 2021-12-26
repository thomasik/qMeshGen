/////////////////////////////////////////////////////////////////////////////
// Curve3dParametric.cpp
// Abstract class describing curve for construction of 3D edges
//	[Curve3dParametric]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "Curve3dParametric.h"
#include "MeshData.h"
#include "DRect.h"
#include "DEquation.h"
#include "Metric3dContext.h"
#include "DLeastSquaresFitting.h"
#include "DLine.h"
#include "DSphere.h"
#include "DLinearQuadric.h"
#include "SurfaceParametric.h"
#include "Curve2dSegment.h"
#include "Curve2dCircle.h"
#include "Curve2dLinearQuadric.h"
#include "Curve3dSurfaceParametric.h"

double Curve3dParametric::m_probe_step = 1.0 / 40;

#define HEAP_DEPTH 20

/// domain check
bool Curve3dParametric::withinDomain( const DPoint3d& pt, double & t ) const
{ 
	t = getParameter( pt, t );

	return (t >= m_param_min) && (t <= m_param_max);
}

/// Virtual destructor
Curve3dParametric::~Curve3dParametric() { }

/////////////////////////////////////////////////////////////////////////////
// Returns array of points (in polyline structure)
void Curve3dParametric::getPolyLineInRange(double t0, double t1, DataVector<DPoint3d> & polyline) const
{
	double factor = 1.0 - MEASURE_PRECISION;
	double heap[HEAP_DEPTH];
	DPoint3d pt0 = getPoint(t0);
	polyline.add(pt0);
	DPoint3d pt1 = getPoint(t1);
	heap[0] = t1;
	int top = 0;
	while(true){
		double t_middle = (t0+t1)*0.5;
		DPoint3d pt_middle = getPoint(t_middle);
		double distance = pt0.distance(pt1);
		double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
		if((distance >= factor * real_distance) || (top == HEAP_DEPTH-1)){
			polyline.add(pt1);
			if(top == 0) return;
			t0 = t1;
			pt0 = pt1;
			pt1 = getPoint(t1 = heap[--top]);
		}else{						
			heap[++top] = t1 = t_middle;
			pt1 = pt_middle;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
// Returns array of point parameters (in polyline structure)
void Curve3dParametric::getPolyLineInRange(double t0, double t1, DataVector<double>& polyline) const
{
	double factor = 1.0 - MEASURE_PRECISION;
	double heap[HEAP_DEPTH];
	DPoint3d pt0 = getPoint(t0);
	polyline.add(t0);
	DPoint3d pt1 = getPoint(t1);
	heap[0] = t1;
	int top = 0;
	while(true){
		double t_middle = (t0+t1)*0.5;
		DPoint3d pt_middle = getPoint(t_middle);
		double distance = pt0.distance(pt1);
		double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
		if((distance >= factor * real_distance) || (top == HEAP_DEPTH-1)){
			polyline.add(t1);
			if(top == 0) break;
			t0 = t1;
			pt0 = pt1;
			pt1 = getPoint(t1 = heap[--top]);
		}else{						
			heap[++top] = t1 = t_middle;
			pt1 = pt_middle;
		}
	}
}

double Curve3dParametric::getLength(Metric3dContext& mc, double t0, double t1, bool local_metric) const
{
	double len = 0.0;
	
	double factor = 1.0 - MEASURE_PRECISION;
	double heap[20];
	if(local_metric) mc.countMetricAtPoint(getPoint((t0+t1)*0.5));
	DMPoint3d pt0 = mc.transformRStoMS(getPoint(t0));
	DMPoint3d pt1 = mc.transformRStoMS(getPoint(t1));
	heap[0] = t1;
	int top = 0;
	while(true){
		double t_middle = (t0+t1)*0.5;
		DMPoint3d pt_middle = mc.transformRStoMS(getPoint(t_middle));
		double distance = pt0.distance(pt1);
		double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
		if((distance >= factor * real_distance) || (top == 19)){
			len += real_distance;
			if(top == 0) break;
			t0 = t1;
			pt0 = pt1;
			pt1 = mc.transformRStoMS(getPoint(t1 = heap[--top]));
		}else{						
			heap[++top] = t1 = t_middle;
			pt1 = pt_middle;
		}
	}	
	return len;
}

double Curve3dParametric::checkAndGetLength(Metric3dContext& mc, double t0, double& end_t1, double max_len, bool local_metric) const
{
	double len = 0.0;
	double t1 = end_t1;
	
	double factor = 1.0 - MEASURE_PRECISION;
	double heap[20];
	if(local_metric) mc.countMetricAtPoint(getPoint((t0+t1)*0.5));
	DMPoint3d pt0 = mc.transformRStoMS(getPoint(t0));
	DMPoint3d pt1 = mc.transformRStoMS(getPoint(t1));
	heap[0] = t1;
	int top = 0;
	while(true){
		double t_middle = (t0+t1)*0.5;
		DMPoint3d pt_middle = mc.transformRStoMS(getPoint(t_middle));
		double distance = pt0.distance(pt1);
		double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
		if((distance >= factor * real_distance) || (top == 19)){
			if(len + real_distance > max_len){
				end_t1 = t0 + (max_len - len)/real_distance * (t1-t0);
				return max_len;
			}
			len += real_distance;
			if(top == 0) break;
			t0 = t1;
			pt0 = pt1;
			pt1 = mc.transformRStoMS(getPoint(t1 = heap[--top]));
		}else{						
			heap[++top] = t1 = t_middle;
			pt1 = pt_middle;
		}
	}	
	return len;
}

double Curve3dParametric::getCurvature(double t, double* cdt_len) const
{
	const DVector3d cdt   = getDerivative(t);
	const DVector3d cdtt  = getSecondDerivative(t);

	// curvature
	double ca = abs(cdt.crossProduct(cdtt).length());
	double cb = cdt.length2();
	assert(cb > 0.0);
	if(cdt_len) cb *= (*cdt_len = sqrt(cb));
	else cb *= sqrt(cb);
	return /* cr= */ ca / cb;
}

double Curve3dParametric::getParameterInRange(const DPoint3d& pt, double ts, double t_min, double t_max) const {

	const double	F_ERR		= METRIC_SMALL_NUMBER;	// Maksymalny b³¹d wyniku
	const double	F_ERR2		= F_ERR*F_ERR;	// Maksymalny b³¹d parametru
	const int		MAX_STEPS	= 10;	// Maksymalna iloœæ kroków

	DVector3d dft = getPoint(t_min) - pt;
	if(abs(dft.x) + abs(dft.y + abs(dft.z)) < F_ERR) return t_min;
	dft = getPoint(t_max) - pt;
	if(abs(dft.x) + abs(dft.y) + abs(dft.z) < F_ERR) return t_max;

	double dt = (t_max - t_min) / 100;
	bool found = false;
	// first - find bracketing segment
	double t0, t1;
	for(int i = 0; (i < MAX_STEPS) && !found; i++){
		t0 = t1 = ts;
		DVector3d dft0 = getPoint(t0) - pt;
		if(abs(dft0.x) + abs(dft0.y) + abs(dft0.z) < F_ERR) return ts;
		DVector3d dft1 = dft0;
		while(t0 > t_min || t1 < t_max){
			if(t0 > t_min){
				double next_t0 = std::max(t_min, t0-dt);
				DVector3d next_dft0 = getPoint(next_t0) - pt;
				if(abs(next_dft0.x) + abs(next_dft0.y) + abs(next_dft0.z) < F_ERR) return next_t0;
				if( (next_dft0.x * dft0.x < F_ERR2) && 
					(next_dft0.y * dft0.y < F_ERR2) &&
					(next_dft0.z * dft0.z < F_ERR2) )
				{	
					// found bracketing segment
					t1 = t0; t0 = next_t0;	found = true; break;
				}else{
					t0 = next_t0; dft0 = next_dft0;
				}
			}
			if(t1 < t_max){
				double next_t1 = std::min(t_max, t1+dt);
				DVector3d next_dft1 = getPoint(next_t1) - pt;
				if(abs(next_dft1.x) + abs(next_dft1.y) + abs(next_dft1.z) < F_ERR) return next_t1;
//				if((next_dft1.x * dft1.x < 0.0) && (next_dft1.y * dft1.y < 0.0)){	// found bracketing segment
				if( (next_dft1.x * dft1.x < F_ERR2) && 
					(next_dft1.y * dft1.y < F_ERR2) &&
					(next_dft1.z * dft1.z < F_ERR2) )
				{	
					// found bracketing segment
					t0 = t1; t1 = next_t1;	found = true; break;
				}else{
					t1 = next_t1; dft1 = next_dft1;
				}
			}
		}
		if(!found)
			dt *= 0.1;
	}
	assert(found);
	if(!found) return -1.0;

	// polish root
	for(int j = 0; j < MAX_STEPS; j++){
		// -> try bisection
		const DVector3d dft0 = getPoint(t0)-pt;
		for(int i=0; i < MAX_STEPS; i++){
			double t = 0.5 * (t0+t1);
			const DVector3d _dft = getPoint(t)-pt;
			if( (dft0.x * _dft.x < F_ERR2) && 
				(dft0.y * _dft.y < F_ERR2) &&
				(dft0.z * _dft.z < F_ERR2) )	// left
				t1 = t;
			else
				t0 = t;
		}
		if((t1-t0) < F_ERR) return 0.5*(t1+t0);
		// Newton method
		double t = 0.5 * (t0+t1);
		for(int i=0; i < MAX_STEPS; i++){
			const DPoint3d ft = getPoint(t);
			const DVector3d fd = ft - pt;
			if(abs(fd.x) + abs(fd.y) + abs(fd.z) < F_ERR) return t;
			const DVector3d _dft = getDerivative(t);
			double dtx = (abs(_dft.x) < VERY_SMALL_NUMBER) ? 0.0 : (fd.x / _dft.x);
			double dty = (abs(_dft.y) < VERY_SMALL_NUMBER) ? 0.0 : (fd.y / _dft.y);
			double dtz = (abs(_dft.z) < VERY_SMALL_NUMBER) ? 0.0 : (fd.z / _dft.z);
			if(abs(dtx) > abs(dty) && abs(dtx) > abs(dtz))
				t -= dtx;
			else if(abs(dty) > abs(dtz))
				t -= dty;
			else
				t -= dtz;
			if(abs(dtx) + abs(dty) + abs(dtz) < F_ERR) return t;
		}
	}

	LOG4CPLUS_WARN(MeshLog::logger_console, 
		"getParameter search failed - switching to linear");

	dt = (t_max - t_min) / 100;
	t0 = t_min;
	t1 = t_max;
	for(int i=0; i < MAX_STEPS; i++){
		double min_diff = SMALL_NUMBER;
		double best_t = ts;
		for(double t = t0+dt; t < t1; t += dt){
			const DVector3d _dft = getPoint(t) - pt;
			double diff = abs(_dft.x) + abs(_dft.y) + abs(_dft.z);
			if(diff < F_ERR) return t;
			else if(diff < min_diff){
				min_diff = diff;
				best_t = t;
			}
		}
		t0 = best_t - dt;
		t1 = best_t + dt;
		dt /= 10;
	}

	LOG4CPLUS_WARN(MeshLog::logger_console, 
		"getParameter linear search failed - give up...");
	assert(false);
	return ts;
}

double Curve3dParametric::getParameter( const DPoint3d& pt, double ts ) const
{
	return getParameterInRange( pt, ts, m_param_min, m_param_max );
}

DBox Curve3dParametric::getBoundingBoxInRange(double t0, double t1) const {
	DataVector<DPoint3d> polyline;
	getPolyLineInRange(t0, t1, polyline);
	DBox box; 
	for(size_t i = 0; i < polyline.countInt(); i++)
		box.addPoint(polyline[i]);
	return box;
}


Curve3dParametricSet::Curve3dParametricSet(Curve3dConstPtr curve)
{
	if(curve) m_curves.add(curve);
}

Curve3dParametricSet::Curve3dParametricSet(std::shared_ptr<const Curve3dParametricSet> cset, Curve3dConstPtr curve)
{
	if(cset) // copy from other set
		for(int i = 0; i < cset->countInt(); i++)
			m_curves.add(cset->getCurve(i));
	if(curve) // and add an additional curve
		m_curves.add(curve);
}

int Curve3dParametricSet::addCurve(Curve3dConstPtr curve) {
	return (int)m_curves.add(curve);
}

int Curve3dParametricSet::count() const { 
	return (int)m_curves.countInt(); 
}

//Curve3dConstPtr Curve3dParametricSet::selectCurveForPoint( double t ) const
//{
//	for(int i = 0; i < m_curves.countInt(); i++){
//		if( m_curves[i]->withinDomain(t) ) return m_curves[i];
//	}
//	return nullptr;
//}

/// Store XML description to stream
ostream& Curve3dParametric::storeXML(ostream& os, const string& prefix) const
{
	return os 
		<< prefix << "<param-min>" << m_param_min << "</param-min>\n"
		<< prefix << "<param-max>" << m_param_max << "</param-max>\n";
}

/// Fit 2d curve on parametric surface
Curve3dPtr Curve3dParametric::fitCurveOnSurface(
		Metric3dContext& mc, SurfaceConstPtr surface, 
		const DataVector<DPoint2d> & points, DataVector<double> & params, double tolerance)
{

	//DLeastSquaresFitting::testFit();

	int pct = points.countInt();
	// -- create candidate curves
	DataVector<Curve2dPtr> cand_curves;

	DLine2d line;
	if( pct > 1 ) {
		// a) line
		DLeastSquaresFitting::fitLine( points, line, false );
		auto curve = std::make_shared<Curve2dSegment>(line.m_pt, line.m_vt);
		cand_curves.add(curve);
		curve->setMinMaxParam( -LARGE_NUMBER, LARGE_NUMBER );
	}

	if( pct > 2 ) {
		// b) quadric-on-line ?
		DLinearQuadric lquadric;
		DLeastSquaresFitting::fitLinearQuadric( points, line, lquadric, false );
		auto curve = std::make_shared<Curve2dLinearQuadric>(lquadric);
		cand_curves.add(curve);
		curve->setMinMaxParam( -LARGE_NUMBER, LARGE_NUMBER );
	}

	if( pct > 3 ) {
		// c) circle
		DCircle circle;
		DLeastSquaresFitting::fitCircle( points, circle, false ); // best/default - Bullock
		auto curve = std::make_shared<Curve2dCircle>(circle);
		cand_curves.add(curve);
		curve->setMinMaxParam( -1.0, 2.0 ); // [0,1] +- 1
	}

	size_t n = cand_curves.countInt();
	assert( n > 0 );
	DataVector< DataVector<double> > cand_params(n, DataVector<double>());
	DataVector< double > cand_max_dist2(n, 0.0); 
	DataVector< double > cand_min_param(n, -LARGE_NUMBER );
	DataVector< double > cand_max_param(n,  LARGE_NUMBER );


	for(size_t k = 0; k < pct; k++){
		const DPoint3d pt_3d = surface->getPoint( points[k] );
		mc.countMetricAtPoint( pt_3d );

		for(size_t j = 0; j < n; j++ ){
			double last_t = cand_params[j].empty() ? 0.0 : cand_params[j].last();
			double t = cand_curves[j]->getParameter( points[k], last_t );
			if( k == 0 ) {
				cand_min_param[j] = t; 
				cand_max_param[j] = t;
			}else {
				if( t < cand_min_param[j] ) cand_min_param[j] = t;
				if( t > cand_max_param[j] ) cand_max_param[j] = t;
			}
			cand_params[j].add(t); // for circle -> enforce monotonous??!!
			double d = mc.transformRStoMS( surface->getPoint( cand_curves[j]->getPoint( t ) ) - pt_3d ).length2();
			if(d > cand_max_dist2[j]) cand_max_dist2[j] = d;
		}
	}

	int best_curve = 0;
	// -> select best ?
	double tol2 = tolerance*tolerance;
	if( cand_max_dist2[0] > 0.01*tol2 ) {
		for(int i = 1; i < n; i++)
			if( cand_max_dist2[i] < cand_max_dist2[best_curve] ) best_curve = i;

		if( cand_max_dist2[best_curve] > tol2 ) best_curve = -1; // none is good enough
	}// else leave 0 as simple line which is good enough

	if( best_curve >= 0 ) {
		params.clear();
		for(int k = 0; k < pct; k++){
			params.add( cand_params[best_curve][k] );
		}
	}

	if( best_curve >= 0 ){
		cand_curves[best_curve]->setMinMaxParam( cand_min_param[best_curve], cand_max_param[best_curve] );
		LOG4CPLUS_INFO(MeshLog::logger_mesh, 
			" * Local curve for " << pct << "pts (" 
			<< cand_curves[best_curve]->getSimpleDescription() 
			<< "), max_dist=" << sqrt(cand_max_dist2[best_curve]));
		return std::make_shared<Curve3dSurfaceParametric>( surface, cand_curves[best_curve] );
	}else 
		return nullptr;
}
