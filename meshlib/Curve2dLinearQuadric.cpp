/////////////////////////////////////////////////////////////////////////////
// Curve2dLinearQuadric.cpp
// Linear quadric 2D
//	[Curve2dParametric->Curve2dLinearQuadric]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve2dLinearQuadric.h"
#include "DPoint.h"
#include "DEquation.h"

/// t in [0,1]
DPoint2d Curve2dLinearQuadric::getPoint(double t) const {
	return m_lq.getPoint(t);
}
	
double Curve2dLinearQuadric::getParameterInRange( const DPoint2d& pt, double ts, double t_min, double t_max ) const {
	double t0 = m_lq.line.paramForPoint(pt);
	if( t0 < t_min) t0 = t_min;
	else if(t0 > t_max) t0 = t_max;
	//return Curve2dParametric::getParameterInRange( pt, t0, t_min, t_max );
	return t0;
}

/// t in [0,1]
DVector2d Curve2dLinearQuadric::getDerivative(double t) const {
	return m_lq.line.m_vt + m_lq.line.m_vn * ( m_lq.vq[1] + 2 * t * m_lq.vq[2] );
}

/// t in [0,1]
DVector2d Curve2dLinearQuadric::getSecondDerivative(double t) const {
	return m_lq.line.m_vn * ( 2 * m_lq.vq[2] );
}

/// t in [0,1]
DVector2d Curve2dLinearQuadric::getThirdDerivative(double t) const {
	return DVector2d::zero;
}

ostream& operator<<(ostream& os, const Curve2dLinearQuadric *lq){
	return os << lq->m_lq;
}

istream& operator>>(istream& is, Curve2dLinearQuadric *lq){
	return is >> lq->m_lq;
}

/// Store XML description to stream
ostream& Curve2dLinearQuadric::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<linearquadric>\n";
	Curve2dParametric::storeXML(os, prefix + "\t");
	os << prefix << "\t<pt0> " << m_lq.line.m_pt << " </pt0>\n";
	os << prefix << "\t<vt> "  << m_lq.line.m_vt << " </vt>\n";
	os << prefix << "\t<vq> "  << m_lq.vq << " </vq>\n";
	return os << prefix << "</linearquadric>\n";
}
