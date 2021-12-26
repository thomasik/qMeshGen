/////////////////////////////////////////////////////////////////////////////
// Curve2dAnalytic.cpp
// Parametric curve 2D with { u(t), v(t) } given directly through analytic formulations
//	[Curve2dParametric->Curve2dAnalytic]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve2dAnalytic.h"
#include "DEquation.h"
#include "MeshLog.h"
#include "MeshData.h"

Curve2dAnalytic::Curve2dAnalytic(string str_x, string str_y, DEquationConstTable *ctable)
{
	m_valid  = m_eq_x.parse((m_str_x = str_x).c_str(), ctable);
	m_valid &= m_eq_y.parse((m_str_y = str_y).c_str(), ctable);
	if(m_valid){
		m_eq_xt = m_eq_x.getDerivative(DEquation::op_x);
		m_eq_yt = m_eq_y.getDerivative(DEquation::op_x);
		m_eq_xtt = m_eq_xt->getDerivative(DEquation::op_x);
		m_eq_ytt = m_eq_yt->getDerivative(DEquation::op_x);
//		m_eq_xttt = m_eq_xtt->getDerivative(DEquation::op_x);
//		m_eq_yttt = m_eq_ytt->getDerivative(DEquation::op_x);
		m_eq_xttt = m_eq_yttt = nullptr;	// third derivative is not really used now
	}else{
		m_eq_xt = m_eq_yt = nullptr;
		m_eq_xtt = m_eq_ytt = nullptr;
		m_eq_xttt = m_eq_yttt = nullptr;
	}
}

Curve2dAnalytic::~Curve2dAnalytic()
{
	if(m_eq_xt) delete m_eq_xt;
	if(m_eq_yt) delete m_eq_yt;
	if(m_eq_xtt) delete m_eq_xtt;
	if(m_eq_ytt) delete m_eq_ytt;
	if(m_eq_xttt) delete m_eq_xttt;
	if(m_eq_yttt) delete m_eq_yttt;
}

DPoint2d Curve2dAnalytic::getPoint(double t) const {
	return DPoint2d(m_eq_x.getValue(t), m_eq_y.getValue(t));
}

DVector2d Curve2dAnalytic::getDerivative(double t) const {
	assert(m_eq_xt && m_eq_yt);
	return DVector2d(m_eq_xt->getValue(t), m_eq_yt->getValue(t));
}

DVector2d Curve2dAnalytic::getSecondDerivative(double t) const {
	assert(m_eq_xtt && m_eq_ytt);
	return DVector2d(m_eq_xtt->getValue(t), m_eq_ytt->getValue(t));
}

DVector2d Curve2dAnalytic::getThirdDerivative(double t) const {
	assert(m_eq_xttt && m_eq_yttt);
	return DVector2d(m_eq_xttt->getValue(t), m_eq_yttt->getValue(t));
}

ostream& operator<<(ostream& os, const Curve2dAnalytic *figure){
	os << "X=" << figure->m_str_x << " Y=" << figure->m_str_y;
	return os;
}

istream& operator>>(istream& is, Curve2dAnalytic *figure){
	is >> ws >> figure->m_str_x >> ws >> figure->m_str_y;
	if(figure->m_eq_xt) { delete figure->m_eq_xt; figure->m_eq_xt = nullptr; }
	if(figure->m_eq_yt) { delete figure->m_eq_yt; figure->m_eq_yt = nullptr; }
	if(figure->m_eq_xtt) { delete figure->m_eq_xtt; figure->m_eq_xtt = nullptr; }
	if(figure->m_eq_ytt) { delete figure->m_eq_ytt; figure->m_eq_ytt = nullptr; }
	if(figure->m_eq_xttt) { delete figure->m_eq_xttt; figure->m_eq_xttt = nullptr; }
	if(figure->m_eq_yttt) { delete figure->m_eq_yttt; figure->m_eq_yttt = nullptr; }
	if(is){
		figure->m_str_x[0] = figure->m_str_x[1] = ' ';
		figure->m_valid  = figure->m_eq_x.parse(figure->m_str_x.c_str()+2);
		figure->m_str_y[0] = figure->m_str_y[1] = ' ';
		figure->m_valid &= figure->m_eq_y.parse(figure->m_str_y.c_str()+2);
		if(figure->m_valid){
			figure->m_eq_xt = figure->m_eq_x.getDerivative(DEquation::op_x);
			figure->m_eq_yt = figure->m_eq_y.getDerivative(DEquation::op_x);
			figure->m_eq_xtt = figure->m_eq_xt->getDerivative(DEquation::op_x);
			figure->m_eq_ytt = figure->m_eq_yt->getDerivative(DEquation::op_x);
//			figure->m_eq_xttt = figure->m_eq_xtt->getDerivative(DEquation::op_x);
//			figure->m_eq_yttt = figure->m_eq_ytt->getDerivative(DEquation::op_x);
		}else{
			figure->m_eq_xt = figure->m_eq_yt = nullptr;
			figure->m_eq_xtt = figure->m_eq_ytt = nullptr;
			figure->m_eq_xttt = figure->m_eq_yttt = nullptr;
		}
	}else{
		figure->m_eq_x.clear();
		figure->m_eq_y.clear();
		figure->m_valid = false;
	}
	return is;
}
