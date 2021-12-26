/////////////////////////////////////////////////////////////////////////////
// Curve2dParametric.cpp
// Abstract class for general parametrized curves 2D: u(t), v(t)
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve2dParametric.h"
#include "MeshData.h"
#include "ControlSpace2d.h"
#include "DRect.h"
#include "SurfaceParametric.h"
#include "DEquation.h"
#include "Metric2dContext.h"

double Curve2dParametric::m_probe_step = 1.0 / 40;

#define MAX_POINTS 300
#define HEAP_DEPTH 20

/////////////////////////////////////////////////////////////////////////////
// Zwraca tablicê punktów definiuj¹cych liniê ³aman¹ oPIsuj¹c¹ kszta³t tej 
//	krzywej
DPoint2d* Curve2dParametric::getPolyLineInRange(int &ct, double t0, double t1) const
{
	DPoint2d poly_array[MAX_POINTS];
	int threshold_factor = 1;
	while(true){
		double factor = 1.0 - threshold_factor * MEASURE_PRECISION;
		double heap[HEAP_DEPTH];
		DPoint2d pt0 = getPoint(t0);
		poly_array[0] = pt0;
		ct = 1;
		DPoint2d pt1 = getPoint(t1);
		heap[0] = t1;
		int top = 0;
		while(true){
			double t_middle = (t0+t1)*0.5;
			DPoint2d pt_middle = getPoint(t_middle);
			double distance = pt0.distance(pt1);
			double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
			if((distance >= factor * real_distance) || (top == HEAP_DEPTH-1)){
				poly_array[ct++] = pt1;
				if(ct >= MAX_POINTS) break;
				if(top == 0){
					DPoint2d* plot_points = new DPoint2d[ct];
					for(int k = 0; k < ct; k++)
						plot_points[k] = poly_array[k];
					return plot_points;
				}
				t0 = t1;
				pt0 = pt1;
				pt1 = getPoint(t1 = heap[--top]);
			}else{						
				heap[++top] = t1 = t_middle;
				pt1 = pt_middle;
			}
		}
		threshold_factor *= 2;
	}

	ct = 0;
	return nullptr;
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca tablicê punktów definiuj¹cych liniê ³aman¹ oPIsuj¹c¹ kszta³t tej 
//	krzywej
void Curve2dParametric::getPolyLineInRange(double t0, double t1, DataVector<double>& polyline) const
{
	double factor = 1.0 - MEASURE_PRECISION;
	double heap[HEAP_DEPTH];
	DPoint2d pt0 = getPoint(t0);
	polyline.add(t0);
	DPoint2d pt1 = getPoint(t1);
	heap[0] = t1;
	int top = 0;
	while(true){
		double t_middle = (t0+t1)*0.5;
		DPoint2d pt_middle = getPoint(t_middle);
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

double Curve2dParametric::getLength(Metric2dContext& mc, double t0, double t1, bool local_metric) const 
{
	double len = 0.0;
	
	double factor = 1.0 - MEASURE_PRECISION;
	double heap[20];
	if(local_metric) mc.countMetricAtPoint(getPoint((t0+t1)*0.5));
	DMPoint2d pt0 = mc.transformPStoMS(getPoint(t0));
	DMPoint2d pt1 = mc.transformPStoMS(getPoint(t1));
	heap[0] = t1;
	int top = 0;
	while(true){
		double t_middle = (t0+t1)*0.5;
		DMPoint2d pt_middle = mc.transformPStoMS(getPoint(t_middle));
		double distance = pt0.distance(pt1);
		double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
		if((distance >= factor * real_distance) || (top == 19)){
			len += real_distance;
			if(top == 0) break;
			t0 = t1;
			pt0 = pt1;
			pt1 = mc.transformPStoMS(getPoint(t1 = heap[--top]));
		}else{						
			heap[++top] = t1 = t_middle;
			pt1 = pt_middle;
		}
	}	
	return len;
}

double Curve2dParametric::getLengthOnSurface(double t0, double t1, SurfaceConstPtr surface) const 
{
	double len = 0.0;
	
	double factor = 1.0 - MEASURE_PRECISION;
	double heap[20];
	DPoint3d pt0 = surface->getPoint(getPoint(t0));
	DPoint3d pt1 = surface->getPoint(getPoint(t1));
	heap[0] = t1;
	int top = 0;
	while(true){
		double t_middle = (t0+t1)*0.5;
		DPoint3d pt_middle = surface->getPoint(getPoint(t_middle));
		double distance = pt0.distance(pt1);
		double real_distance = pt0.distance(pt_middle)+ pt_middle.distance(pt1);
		if((distance >= factor * real_distance) || (top == 19)){
			len += real_distance;
			if(top == 0) break;
			t0 = t1;
			pt0 = pt1;
			pt1 = surface->getPoint(getPoint(t1 = heap[--top]));
		}else{						
			heap[++top] = t1 = t_middle;
			pt1 = pt_middle;
		}
	}	
	return len;
}

double Curve2dParametric::checkAndGetLength(Metric2dContext& mc, double t0, double& end_t1, double max_len, bool local_metric) const 
{
	double len = 0.0;
	double t1 = end_t1;
	
	double factor = 1.0 - MEASURE_PRECISION;
	double heap[20];
	if(local_metric) mc.countMetricAtPoint(getPoint((t0+t1)*0.5));
	DMPoint2d pt0 = mc.transformPStoMS(getPoint(t0));
	DMPoint2d pt1 = mc.transformPStoMS(getPoint(t1));
	heap[0] = t1;
	int top = 0;
	while(true){
		double t_middle = (t0+t1)*0.5;
		DMPoint2d pt_middle = mc.transformPStoMS(getPoint(t_middle));
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
			pt1 = mc.transformPStoMS(getPoint(t1 = heap[--top]));
		}else{						
			heap[++top] = t1 = t_middle;
			pt1 = pt_middle;
		}
	}	
	return len;
}

double Curve2dParametric::getPlanarCurvature(double t) const
{
	const DVector2d dt = getDerivative(t);
	const DVector2d dtt = getSecondDerivative(t);
	double a = abs(dt.x*dtt.y - dt.y*dtt.x);
	double b = dt.length2();
	assert(b > 0.0);
	b *= sqrt(b);

	return a / b;
}

double Curve2dParametric::getNonPlanarCurvature(SurfaceConstPtr surface, double t, double* cdt_len) const
{
	const DPoint2d pt  = getPoint(t);
	const DVector3d fs  = surface->getDerivative(DEquation::deriv_ds, pt);
	const DVector3d ft  = surface->getDerivative(DEquation::deriv_dt, pt);

	// check
	double l1 = fs.length2();
	double l2 = ft.length2();
	double g_ratio = std::min(l1,l2) / std::max(l1,l2);
	if(g_ratio < MIN_PARAM_GRATIO) return 0.0;

	const DVector3d fss = surface->getDerivative(DEquation::deriv_dss, pt);
	const DVector3d fst = surface->getDerivative(DEquation::deriv_dst, pt);
	const DVector3d ftt = surface->getDerivative(DEquation::deriv_dtt, pt);
//	const DPoint3d fsss = surface->getDerivative(DEquation::deriv_dsss, pt);
//	const DPoint3d fsst = surface->getDerivative(DEquation::deriv_dsst, pt);
//	const DPoint3d fstt = surface->getDerivative(DEquation::deriv_dstt, pt);
//	const DPoint3d fttt = surface->getDerivative(DEquation::deriv_dttt, pt);
	const DVector2d dt   = getDerivative(t);
	const DVector2d dtt  = getSecondDerivative(t);
//	const DPoint2d dttt = getThirdDerivative(t);

	const DVector3d cdt = fs * dt.x + ft * dt.y;
	const DVector3d cdtt = fss * (dt.x*dt.x) + fst * (2*dt.x*dt.y) + ftt * (dt.y*dt.y)
		+ fs * dtt.x + ft * dtt.y;
//	const DPoint3d cdttt = 
//		((fsss*dt.x + fsst*(3*dt.y))*dt.x + (fss*dtt.x+fst*dtt.y)*3)*dt.x + fs*dttt.x + 
//		((fttt*dt.y + fstt*(3*dt.x))*dt.y + (ftt*dtt.y+fst*dtt.x)*3)*dt.y + ft*dttt.y;

	// curvature
	double ca = abs(cdt.crossProduct(cdtt).length());
	double cb = cdt.length2();
	assert(cb > 0.0);
	if(cdt_len) cb *= (*cdt_len = sqrt(cb));
	else cb *= sqrt(cb);
	return /* cr= */ ca / cb;

/*
	// torsion
	double tb = ca*ca;
	if(tb > mesh_data.relative_small_number){
		double ta = DPoint3d::det33(cdt.x, cdtt.x, cdttt.x, cdt.y, cdtt.y, cdttt.y, cdt.z, cdtt.z, cdttt.z);
		double tr = ta / tb;

		DPoint3d surf_cr = surface->getCurvature(pt);
		MESHLOG.precision(4);
		MESHLOG.width(10);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, t << '\t' << cr << '\t' << tr << '\t' << surf_cr);
	}else{
		DPoint3d surf_cr = surface->getCurvature(pt);
		MESHLOG.precision(4);
		MESHLOG.width(10);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, t << '\t' << cr << '\t' << '-' << '\t' << surf_cr);
	}
	return cr;
*/
}

double Curve2dParametric::getParameterInRange(const DPoint2d& pt, double ts, double t_min, double t_max) const {

	const double	F_ERR		= METRIC_SMALL_NUMBER;	// maximum error of result
	const double	F_ERR2		= F_ERR*F_ERR;	// maximum error of parameter
	const int		MAX_STEPS	= 10;	// maximum number of steps

	DVector2d dft = getPoint(t_min) - pt;
	if(abs(dft.x) + abs(dft.y) < F_ERR) return t_min;
	dft = getPoint(t_max) - pt;
	if(abs(dft.x) + abs(dft.y) < F_ERR) return t_max;

	double dt = (t_max - t_min) / 10;
	bool found = false;
	// first - find bracketing segment
	double t0, t1;
	for(int i = 0; (i < MAX_STEPS) && !found; i++){
		t0 = t1 = ts;
		DVector2d dft0 = getPoint(t0) - pt;
		if(abs(dft0.x) + abs(dft0.y) < F_ERR) return ts;
		DVector2d dft1 = dft0;
		while(t0 > t_min || t1 < t_max){
			if(t0 > t_min){
				double next_t0 = std::max(t_min, t0-dt);
				DVector2d next_dft0 = getPoint(next_t0) - pt;
				if(abs(dft0.x) + abs(dft0.y) < F_ERR) return next_t0;
//				if((next_dft0.x * dft0.x < 0.0) && (next_dft0.y * dft0.y < 0.0)){	// found bracketing segment
				if((next_dft0.x * dft0.x < F_ERR2) && 
					(next_dft0.y * dft0.y < F_ERR2)){	// found bracketing segment
					t1 = t0; t0 = next_t0;	found = true; break;
				}else{
					t0 = next_t0; dft0 = next_dft0;
				}
			}
			if(t1 < t_max){
				double next_t1 = std::min(t_max, t1+dt);
				DVector2d next_dft1 = getPoint(next_t1) - pt;
				if(abs(dft1.x) + abs(dft1.y) < F_ERR) return next_t1;
//				if((next_dft1.x * dft1.x < 0.0) && (next_dft1.y * dft1.y < 0.0)){	// found bracketing segment
				if((next_dft1.x * dft1.x < F_ERR2) && 
					(next_dft1.y * dft1.y < F_ERR2)){	// found bracketing segment
					t0 = t1; t1 = next_t1;	found = true; break;
				}else{
					t1 = next_t1; dft1 = next_dft1;
				}
			}
		}
		if(!found)
			dt *= 0.5;
	}
//	assert(found);
	if(found){
		// polish root
		for(int j = 0; j < MAX_STEPS; j++){
			// -> try bisection
			const DVector2d dft0 = getPoint(t0)-pt;
			for(int i=0; i < MAX_STEPS; i++){
				double t = 0.5 * (t0+t1);
				const DVector2d _dft = getPoint(t)-pt;
				if((dft0.x * _dft.x < F_ERR2) && 
					(dft0.y * _dft.y < F_ERR2))	// left
					t1 = t;
				else
					t0 = t;
			}
			if((t1-t0) < F_ERR) return 0.5*(t1+t0);
			// Newton method
			double t = 0.5 * (t0+t1);
			for(int i=0; i < MAX_STEPS; i++){
				const DPoint2d ft = getPoint(t);
				const DVector2d fd = ft - pt;
				if(abs(fd.x) + abs(fd.y) < F_ERR) return t;
				const DVector2d _dft = getDerivative(t);
				double dtx = (abs(_dft.x) < VERY_SMALL_NUMBER) ? 0.0 : (fd.x / _dft.x);
				double dty = (abs(_dft.y) < VERY_SMALL_NUMBER) ? 0.0 : (fd.y / _dft.y);
				t -= (abs(dtx) > abs(dty)) ? dtx : dty;
				if(abs(dtx) + abs(dty) < F_ERR) return t;
			}
		}
	}

	LOG4CPLUS_WARN(MeshLog::logger_console, 
		"getParameter search failed - switching to linear");

	dt = (t_max - t_min) / 100;
	t0 = t_min;
	t1 = t_max;
	double min_diff = VERY_SMALL_NUMBER;
	double best_t = ts;
	for(int i=0; i < MAX_STEPS; i++){
		for(double t = t0+dt; t < t1; t += dt){
			const DVector2d _dft = getPoint(t) - pt;
			double diff = abs(_dft.x) + abs(_dft.y);
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
		"getParameter linear search last resort, min_diff=" << pt.distance(getPoint(best_t)));
//	assert(false);
	return best_t;
}

double Curve2dParametric::getParameter( const DPoint2d& pt, double ts ) const {
	return getParameterInRange( pt, ts, m_param_min, m_param_max );
}

DRect Curve2dParametric::getBoundingRect(double t0, double t1) const {
	DRect rect; 
	int step_count = 10;
	double dt = (t1-t0)/10;
	rect.addPoint(getPoint(t0));
	rect.addPoint(getPoint(t1));
	for(int i = 1; i < step_count; i++)
		rect.addPoint(getPoint(t0+i*dt));
	return rect;
}

/// Store XML description to stream
ostream& Curve2dParametric::storeXML(ostream& os, const string& prefix) const
{
	return os 
		<< prefix << "<param-min>" << m_param_min << "</param-min>\n"
		<< prefix << "<param-max>" << m_param_max << "</param-max>\n";
}
