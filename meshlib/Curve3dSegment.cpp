/////////////////////////////////////////////////////////////////////////////
// Curve3dSegment.cpp
//	[Curve3dParametric->Curve3dSegment]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve3dSegment.h"
#include "MeshData.h"
#include "DRect.h"
#include "Metric3dContext.h"

/////////////////////////////////////////////////////////////////////////////
// Zwraca tablicê punktów definiuj¹cych liniê ³aman¹ oPIsuj¹c¹ kszta³t tej 
//	krzywej
void Curve3dSegment::getPolyLineInRange(double t0, double t1, DataVector<DPoint3d> & polyline) const
{
	polyline.add(getPoint(t0));
	polyline.add(getPoint(t1));
}

void Curve3dSegment::getPolyLineInRange(double t1, double t2, DataVector<double>& polyline) const
{
	polyline.add(t1);
	polyline.add(t2);
}

double Curve3dSegment::getLength(Metric3dContext& mc, double t0, double t1, bool local_metric) const 
{
	if(local_metric) mc.countMetricAtPoint(getPoint((t0+t1)*0.5));

	return mc.transformRStoMS(getPoint(t0) - getPoint(t1)).length();
}

double Curve3dSegment::getCurvature(double /* t */, double* cdt_len) const
{
	if(cdt_len) *cdt_len = m_vt.length();
	return 0.0;
}

double Curve3dSegment::getParameterInRange(const DPoint3d& pt, double ts, double t_min, double t_max) const 
{
	double t = getParameter( pt, ts );
	return std::max(t_min, std::min(t_max, t));
}

double Curve3dSegment::getParameter(const DPoint3d& pt, double /*ts*/ ) const 
{
	double fv[] = { abs(m_vt.x), abs(m_vt.y), abs(m_vt.z) };
	if(fv[0] > fv[1] && fv[0] > fv[2]){
		return (pt.x - m_pt.x) / m_vt.x;
	}else if (fv[1] > fv[2]){
		return (pt.y - m_pt.y) / m_vt.y;
	}else{
		return (pt.z - m_pt.z) / m_vt.z;
	}
}

DBox Curve3dSegment::getBoundingBox(double t0, double t1) const {
	DBox box; 
	box.addPoint(getPoint(t0));
	box.addPoint(getPoint(t1));
	return box;
}

/// Store XML description to stream
ostream& Curve3dSegment::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<segment3d>\n";
	Curve3dParametric::storeXML(os, prefix + "\t");
	os << prefix << "\t<pt0> " << m_pt << " </pt0>\n";
	os << prefix << "\t<pt1> " << m_pt+m_vt << " </pt1>\n";
	return os << prefix << "</segment3d>\n";
}
