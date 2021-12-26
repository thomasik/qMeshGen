/////////////////////////////////////////////////////////////////////////////
// Curve2dSegment.cpp
//	[Curve2dParametric->Curve2dSegment]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2005-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve2dSegment.h"
#include "MeshData.h"
#include "ControlSpace2d.h"
#include "DRect.h"
#include "SurfaceParametric.h"
#include "DEquation.h"
#include "DMatrix.h"

/////////////////////////////////////////////////////////////////////////////
// Zwraca tablicê punktów definiuj¹cych liniê ³aman¹ oPIsuj¹c¹ kszta³t tej 
//	krzywej
DPoint2d* Curve2dSegment::getPolyLineInRange(int &ct, double t0, double t1) const
{
	DPoint2d* plot_points = new DPoint2d[ct=2];
	plot_points[0] = getPoint(t0);
	plot_points[1] = getPoint(t1);
	return plot_points;
}

void Curve2dSegment::getPolyLineInRange(double t1, double t2, DataVector<double>& polyline) const
{
	polyline.add(t1);
	polyline.add(t2);
}

double Curve2dSegment::getLength(Metric2dContext& mc, double t0, double t1, bool local_metric) const 
{
	if(local_metric) mc.countMetricAtPoint(getPoint((t0+t1)*0.5));

	return mc.transformPStoMS(getPoint(t0) - getPoint(t1)).length();
}

double Curve2dSegment::getNonPlanarCurvature(SurfaceConstPtr surface, double t, double* cdt_len) const
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
	const DVector2d dt  = m_vt;
//	const DPoint2d dtt = getSecondDerivative(t);	== 0.0
//	const DPoint2d dttt = getThirdDerivative(t);	== 0.0

	const DVector3d cdt = fs * dt.x + ft * dt.y;
	const DVector3d cdtt = fss * (dt.x*dt.x) + fst * (2*dt.x*dt.y) + ftt * (dt.y*dt.y);
//	const DPoint3d cdttt = 
//		(fsss*dt.x + fsst*(3*dt.y))*dt.x*dt.x +  
//		(fttt*dt.y + fstt*(3*dt.x))*dt.y*dt.y;

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
*/
}

double Curve2dSegment::getParameterInRange(const DPoint2d& pt, double ts, double t_min, double t_max) const 
{
	double t = getParameter( pt, ts );
	return std::max(t_min, std::min(t_max, t));
}

double Curve2dSegment::getParameter(const DPoint2d& pt, double /*ts*/ ) const 
{
	const DVector2d vn(m_vt.y, -m_vt.x);
	DMatrix2d m;
	m.setColumn(0, m_vt);
	m.setColumn(1, vn);
	DVector2d res;
	if(m.solve(pt - m_pt, res)){
		return res.x;
	}else{ // old method
		if(abs(m_vt.x) > abs(m_vt.y)){
			return (pt.x - m_pt.x) / m_vt.x;
		}else{
			return (pt.y - m_pt.y) / m_vt.y;
		}
	}
}

DRect Curve2dSegment::getBoundingRect(double t0, double t1) const {
	DRect rect; 
	rect.addPoint(getPoint(t0));
	rect.addPoint(getPoint(t1));
	return rect;
}

/// Store XML description to stream
ostream& Curve2dSegment::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<segment>\n";
	Curve2dParametric::storeXML(os, prefix + "\t");
	os << prefix << "\t<pt0> " << m_pt << " </pt0>\n";
	os << prefix << "\t<pt1> " << m_pt+m_vt << " </pt1>\n";
	return os << prefix << "</segment>\n";
}
