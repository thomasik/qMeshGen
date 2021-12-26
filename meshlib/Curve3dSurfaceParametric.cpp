/////////////////////////////////////////////////////////////////////////////
// Curve3dSurfaceParametric.cpp
//	[Curve3dParametric->Curve3dSurfaceParametric]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve3dSurfaceParametric.h"
#include "MeshData.h"
#include "DRect.h"
#include "Metric3dContext.h"
#include "Curve2dParametric.h"
#include "SurfaceParametric.h"
#include "DEquation.h"

/// Standard constructor
Curve3dSurfaceParametric::Curve3dSurfaceParametric(SurfaceConstPtr surface, Curve2dConstPtr curve )
	: Curve3dParametric( curve->getMinParam(), curve->getMaxParam() ),
		m_surface(surface), m_curve(curve)
{
}

DPoint3d Curve3dSurfaceParametric::getPoint(double t) const 
{
	return m_surface->getPoint( m_curve->getPoint(t) );
}

/// Returns derivatives for the given parameter (dx/dt, dy/dt)
DVector3d Curve3dSurfaceParametric::getDerivative(double t) const
{
	const DPoint2d pt  = m_curve->getPoint(t);
	const DVector3d fs  = m_surface->getDerivative(DEquation::deriv_ds, pt);
	const DVector3d ft  = m_surface->getDerivative(DEquation::deriv_dt, pt);

	// check
	//double l1 = fs.length2();
	//double l2 = ft.length2();
	//double g_ratio = std::min(l1,l2) / std::max(l1,l2);
	//assert(g_ratio > MIN_PARAM_GRATIO);

	const DVector2d dt = m_curve->getDerivative(t);
	return fs * dt.x + ft * dt.y;
}

/// Returns second derivatives for the given parameter (d2x/dt2, d2y/dt2)
DVector3d Curve3dSurfaceParametric::getSecondDerivative(double t) const
{
	const DPoint2d pt  = m_curve->getPoint(t);
	const DVector3d fs  = m_surface->getDerivative(DEquation::deriv_ds, pt);
	const DVector3d ft  = m_surface->getDerivative(DEquation::deriv_dt, pt);

	// check
	//double l1 = fs.length2();
	//double l2 = ft.length2();
	//double g_ratio = std::min(l1,l2) / std::max(l1,l2);
	//assert(g_ratio > MIN_PARAM_GRATIO);

	const DVector3d fss = m_surface->getDerivative(DEquation::deriv_dss, pt);
	const DVector3d fst = m_surface->getDerivative(DEquation::deriv_dst, pt);
	const DVector3d ftt = m_surface->getDerivative(DEquation::deriv_dtt, pt);
	const DVector2d dt   = m_curve->getDerivative(t);
	const DVector2d dtt  = m_curve->getSecondDerivative(t);

	return fss * (dt.x*dt.x) + fst * (2*dt.x*dt.y) + ftt * (dt.y*dt.y)
		+ fs * dtt.x + ft * dtt.y;
}

/// Returns third derivatives for the given parameter (d3x/dt3, d3y/dt3)
DVector3d Curve3dSurfaceParametric::getThirdDerivative(double t) const
{
	const DPoint2d pt  = m_curve->getPoint(t);
	const DVector3d fs  = m_surface->getDerivative(DEquation::deriv_ds, pt);
	const DVector3d ft  = m_surface->getDerivative(DEquation::deriv_dt, pt);

	// check
	//double l1 = fs.length2();
	//double l2 = ft.length2();
	//double g_ratio = std::min(l1,l2) / std::max(l1,l2);
	//assert(g_ratio > MIN_PARAM_GRATIO);

	const DVector3d fss = m_surface->getDerivative(DEquation::deriv_dss, pt);
	const DVector3d fst = m_surface->getDerivative(DEquation::deriv_dst, pt);
	const DVector3d ftt = m_surface->getDerivative(DEquation::deriv_dtt, pt);
	const DVector3d fsss = m_surface->getDerivative(DEquation::deriv_dsss, pt);
	const DVector3d fsst = m_surface->getDerivative(DEquation::deriv_dsst, pt);
	const DVector3d fstt = m_surface->getDerivative(DEquation::deriv_dstt, pt);
	const DVector3d fttt = m_surface->getDerivative(DEquation::deriv_dttt, pt);

	const DVector2d dt   = m_curve->getDerivative(t);
	const DVector2d dtt  = m_curve->getSecondDerivative(t);
	const DVector2d dttt  = m_curve->getThirdDerivative(t);

	return
		((fsss*dt.x + fsst*(3*dt.y))*dt.x + (fss*dtt.x+fst*dtt.y)*3)*dt.x + fs*dttt.x + 
		((fttt*dt.y + fstt*(3*dt.x))*dt.y + (ftt*dtt.y+fst*dtt.x)*3)*dt.y + ft*dttt.y;
}

/// Checks wheterh this curve is valid (i.e. is properly initialized)
bool Curve3dSurfaceParametric::isValid() const
{
	return m_surface && m_curve && m_surface->isValid();
}

double Curve3dSurfaceParametric::getParameterInRange(const DPoint3d& pt, double ts, double t_min, double t_max) const 
{
	return m_curve->getParameterInRange( m_surface->getParameters(pt), ts, t_min, t_max );
}

/// Store XML description to stream
ostream& Curve3dSurfaceParametric::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<curve3d-surface>\n";
	Curve3dParametric::storeXML(os, prefix + "\t");
	os << prefix << "\t<surface>" << endl;
	m_surface->storeXML(os, prefix+"\t\t");
	os << prefix << "\t</surface>" << endl;
	os << prefix << "\t<curve>" << endl;
	m_curve->storeXML(os, prefix+"\t\t");
	os << prefix << "\t</curve>" << endl;
	return os << prefix << "</curve3d-surface>\n";
}

string Curve3dSurfaceParametric::getSimpleDescription() const 
{ 
	return m_curve->getSimpleDescription() + " on " + m_surface->getSimpleDescription();
}
