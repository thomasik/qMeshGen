// SurfaceAnalytic.cpp: implementation of the SurfaceAnalytic class.
//
//////////////////////////////////////////////////////////////////////

#include "SurfaceAnalytic.h"
#include "MeshLog.h"

SurfaceAnalytic::SurfaceAnalytic(const char* formula_x, const char *formula_y, 
			const char *formula_z, DEquationConstTable *ctable) : SurfaceParametric()
{
	for(int i = 0; i < 15; i++) m_eq_deriv[i] = nullptr;
	for(int i = 0; i < 12; i++) m_eq_deriv_ext[i] = nullptr;
	setData(formula_x, formula_y, formula_z, ctable);
}

const DPoint3d SurfaceAnalytic::getPoint(const DPoint2d& param) const
{
	return DPoint3d(
		m_eq_x.getValue(param.x, param.y), 
		m_eq_y.getValue(param.x, param.y), 
		m_eq_z.getValue(param.x, param.y));
}

const DVector3d SurfaceAnalytic::getDerivative(int deriv, const DPoint2d& param) const 
{
	assert(m_valid);
	switch(deriv){
	case DEquation::deriv_ds:
		return DVector3d(
			m_eq_deriv[0]->getValue(param.x, param.y), 
			m_eq_deriv[1]->getValue(param.x, param.y), 
			m_eq_deriv[2]->getValue(param.x, param.y));
	case DEquation::deriv_dt:
		return DVector3d(
			m_eq_deriv[3]->getValue(param.x, param.y), 
			m_eq_deriv[4]->getValue(param.x, param.y), 
			m_eq_deriv[5]->getValue(param.x, param.y));
	// second derivative
	case DEquation::deriv_dss:
		return DVector3d(
			m_eq_deriv[6]->getValue(param.x, param.y), 
			m_eq_deriv[7]->getValue(param.x, param.y), 
			m_eq_deriv[8]->getValue(param.x, param.y));
	case DEquation::deriv_dst:
		return DVector3d(
			m_eq_deriv[9]->getValue(param.x, param.y), 
			m_eq_deriv[10]->getValue(param.x, param.y), 
			m_eq_deriv[11]->getValue(param.x, param.y));
	case DEquation::deriv_dtt:
		return DVector3d(
			m_eq_deriv[12]->getValue(param.x, param.y), 
			m_eq_deriv[13]->getValue(param.x, param.y), 
			m_eq_deriv[14]->getValue(param.x, param.y));
	// third derivative
/*
	case DEquation::deriv_dsss:
		return DPoint3d(m_eq_deriv_ext[0]->getValue(s, t), 
			m_eq_deriv_ext[1]->getValue(s, t), m_eq_deriv_ext[2]->getValue(s, t));
	case DEquation::deriv_dsst:
		return DPoint3d(m_eq_deriv_ext[3]->getValue(s, t), 
			m_eq_deriv_ext[4]->getValue(s, t), m_eq_deriv_ext[5]->getValue(s, t));
	case DEquation::deriv_dstt:
		return DPoint3d(m_eq_deriv_ext[6]->getValue(s, t), 
			m_eq_deriv_ext[7]->getValue(s, t), m_eq_deriv_ext[8]->getValue(s, t));
	case DEquation::deriv_dttt:
		return DPoint3d(m_eq_deriv_ext[9]->getValue(s, t), 
			m_eq_deriv_ext[10]->getValue(s, t), m_eq_deriv_ext[11]->getValue(s, t));
*/
	default:
		assert(false);
	}
	return DVector3d(0.0, 0.0, 0.0);
}

void SurfaceAnalytic::setData(const char* formula_x, const char *formula_y, 
		const char *formula_z, DEquationConstTable *ctable)
{
	clear();

	m_str_x = formula_x;
	m_str_y = formula_y;
	m_str_z = formula_z;

	m_valid  = m_eq_x.parse(formula_x, ctable);
	m_valid &= m_eq_y.parse(formula_y, ctable);
	m_valid &= m_eq_z.parse(formula_z, ctable);

	if(m_valid){
		m_eq_deriv[0] = m_eq_x.getDerivative(DEquation::op_s);
		m_eq_deriv[1] = m_eq_y.getDerivative(DEquation::op_s);
		m_eq_deriv[2] = m_eq_z.getDerivative(DEquation::op_s);

		m_eq_deriv[3] = m_eq_x.getDerivative(DEquation::op_t);
		m_eq_deriv[4] = m_eq_y.getDerivative(DEquation::op_t);
		m_eq_deriv[5] = m_eq_z.getDerivative(DEquation::op_t);

		m_eq_deriv[6] = m_eq_deriv[0]->getDerivative(DEquation::op_s);	// dss
		m_eq_deriv[7] = m_eq_deriv[1]->getDerivative(DEquation::op_s);
		m_eq_deriv[8] = m_eq_deriv[2]->getDerivative(DEquation::op_s);

		m_eq_deriv[9]  = m_eq_deriv[0]->getDerivative(DEquation::op_t);	// dst
		m_eq_deriv[10] = m_eq_deriv[1]->getDerivative(DEquation::op_t);
		m_eq_deriv[11] = m_eq_deriv[2]->getDerivative(DEquation::op_t);

		m_eq_deriv[12] = m_eq_deriv[3]->getDerivative(DEquation::op_t);	// dtt
		m_eq_deriv[13] = m_eq_deriv[4]->getDerivative(DEquation::op_t);
		m_eq_deriv[14] = m_eq_deriv[5]->getDerivative(DEquation::op_t);

/*
		m_eq_deriv_ext[0] = m_eq_deriv[6]->getDerivative(DEquation::op_s); // dsss
		m_eq_deriv_ext[1] = m_eq_deriv[7]->getDerivative(DEquation::op_s);
		m_eq_deriv_ext[2] = m_eq_deriv[8]->getDerivative(DEquation::op_s);

		m_eq_deriv_ext[3] = m_eq_deriv[6]->getDerivative(DEquation::op_t); // dsst
		m_eq_deriv_ext[4] = m_eq_deriv[7]->getDerivative(DEquation::op_t);
		m_eq_deriv_ext[5] = m_eq_deriv[8]->getDerivative(DEquation::op_t);

		m_eq_deriv_ext[6] = m_eq_deriv[9]->getDerivative(DEquation::op_t); // dstt
		m_eq_deriv_ext[7] = m_eq_deriv[10]->getDerivative(DEquation::op_t);
		m_eq_deriv_ext[8] = m_eq_deriv[11]->getDerivative(DEquation::op_t);

		m_eq_deriv_ext[9]  = m_eq_deriv[12]->getDerivative(DEquation::op_t); // dttt
		m_eq_deriv_ext[10] = m_eq_deriv[13]->getDerivative(DEquation::op_t);
		m_eq_deriv_ext[11] = m_eq_deriv[14]->getDerivative(DEquation::op_t);
*/
	}
}

ostream& operator<<(ostream& os, const SurfaceAnalytic *surf){
	os << "X=" << surf->m_str_x << " Y=" << surf->m_str_y << " Z=" << surf->m_str_z;
	return os;
}

istream& operator>>(istream& is, SurfaceAnalytic *surf){
	string str_x, str_y, str_z;
	is >> ws >> str_x >> ws >> str_y >> ws >> str_z;
	if(is){
		str_x[0] = str_x[1] = ' ';
		str_y[0] = str_y[1] = ' ';
		str_z[0] = str_z[1] = ' ';
		surf->setData(str_x.c_str(), str_y.c_str(), str_z.c_str());
	}
	return is;
}

void SurfaceAnalytic::clear()
{
	for(int i = 0; i < 15; i++){
		if(m_eq_deriv[i]){
			delete m_eq_deriv[i];
			m_eq_deriv[i] = nullptr;
		}
	}
	for(int i = 0; i < 12; i++){
		if(m_eq_deriv_ext[i]){
			delete m_eq_deriv_ext[i];
			m_eq_deriv_ext[i] = nullptr;
		}
	}
	m_valid = false;
}
