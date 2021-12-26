/////////////////////////////////////////////////////////////////////////////
// Curve3dLinearQuadric.cpp
// Linear quadric 3D
//	[Curve3dParametric->Curve3dLinearQuadric]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve3dLinearQuadric.h"
#include "DPoint.h"
#include "DEquation.h"

/// t in [0,1]
DPoint3d Curve3dLinearQuadric::getPoint(double t) const {
	return m_lq.getPoint(t);
}
	
double Curve3dLinearQuadric::getParameterInRange( const DPoint3d& pt, double ts, double t_min, double t_max ) const {
	double t0 = m_lq.line.paramForPoint(pt);
	if( t0 < t_min) t0 = t_min;
	else if(t0 > t_max) t0 = t_max;
	//return Curve2dParametric::getParameterInRange( pt, t0, t_min, t_max );
	return t0;
}

/// t in [0,1]
DVector3d Curve3dLinearQuadric::getDerivative(double t) const {
	return m_lq.line.m_vt
		+ m_lq.e1 * ( m_lq.vq1[1] + 2 * t * m_lq.vq1[2] )
		+ m_lq.e2 * ( m_lq.vq2[1] + 2 * t * m_lq.vq2[2] );
}

/// t in [0,1]
DVector3d Curve3dLinearQuadric::getSecondDerivative(double t) const {
	return m_lq.e1 * ( 2 * m_lq.vq1[2] ) + m_lq.e2 * ( 2 * m_lq.vq2[2] );
}

/// t in [0,1]
DVector3d Curve3dLinearQuadric::getThirdDerivative(double t) const {
	return DVector3d::zero;
}

ostream& operator<<(ostream& os, const Curve3dLinearQuadric *lq){
	return os << lq->m_lq;
}

istream& operator>>(istream& is, Curve3dLinearQuadric *lq){
	return is >> lq->m_lq;
}

/// Store XML description to stream
ostream& Curve3dLinearQuadric::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<linearquadric3d>\n";
	Curve3dParametric::storeXML(os, prefix + "\t");
	os << prefix << "\t<pt0> " << m_lq.line.m_pt << " </pt0>\n";
	os << prefix << "\t<vt> "  << m_lq.line.m_vt << " </vt>\n";
	os << prefix << "\t<e1> "  << m_lq.e1 << " </e1>\n";
	os << prefix << "\t<e2> "  << m_lq.e2 << " </e2>\n";
	os << prefix << "\t<vq1> "  << m_lq.vq1 << " </vq1>\n";
	os << prefix << "\t<vq2> "  << m_lq.vq2 << " </vq2>\n";
	return os << prefix << "</linearquadric3d>\n";
}
