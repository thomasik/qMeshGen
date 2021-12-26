// SurfaceMonoAnalytic.cpp: implementation of the SurfaceMonoAnalytic class.
//
//////////////////////////////////////////////////////////////////////

#include "SurfaceMonoAnalytic.h"
#include "MeshLog.h"

SurfaceMonoAnalytic::SurfaceMonoAnalytic(const char* formula, FTYPE ftype, DEquationConstTable *ctable)
	: SurfaceParametric()
{
	for(int i = 0; i < 15; i++) m_eq_deriv[i] = nullptr;
	setData(formula, ftype, ctable);
}

const DPoint3d SurfaceMonoAnalytic::getPoint(const DPoint2d& param) const
{
	switch(m_ftype){
	case FXY: return DPoint3d(param.x, param.y, m_eq.getValue(param.x, param.y));
	case FXZ: return DPoint3d(param.x, m_eq.getValue(param.x, param.y), param.y);
	case FYZ: return DPoint3d(m_eq.getValue(param.x, param.y), param.x, param.y);
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "SurfaceMonoAnalytic::getPoint, unknown ftype");
		return DPoint3d::zero;
	}
}

const DVector3d SurfaceMonoAnalytic::getDerivative(int deriv, const DPoint2d& param) const 
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
	default:
		assert(false);
	}
	return DVector3d(0.0, 0.0, 0.0);
}

void SurfaceMonoAnalytic::setData(const char* formula, FTYPE ftype, DEquationConstTable *ctable)
{
	clear();

	m_ftype = ftype;
	m_str = formula;

	m_valid  = m_eq.parse(formula, ctable);

	if(m_valid){
		switch(ftype){
		case FXY:
			m_eq_deriv[0] = new DEquation("1");
			m_eq_deriv[1] = new DEquation("0");
			m_eq_deriv[2] = m_eq.getDerivative(DEquation::op_s);

			m_eq_deriv[3] = new DEquation("0");
			m_eq_deriv[4] = new DEquation("1");
			m_eq_deriv[5] = m_eq.getDerivative(DEquation::op_t);
			break;
		case FXZ:
			m_eq_deriv[0] = new DEquation("1");
			m_eq_deriv[1] = m_eq.getDerivative(DEquation::op_s);
			m_eq_deriv[2] = new DEquation("0");

			m_eq_deriv[3] = new DEquation("0");
			m_eq_deriv[4] = m_eq.getDerivative(DEquation::op_t);
			m_eq_deriv[5] = new DEquation("1");
			break;
		case FYZ:
			m_eq_deriv[0] = m_eq.getDerivative(DEquation::op_s);
			m_eq_deriv[1] = new DEquation("1");
			m_eq_deriv[2] = new DEquation("0");

			m_eq_deriv[3] = m_eq.getDerivative(DEquation::op_t);
			m_eq_deriv[4] = new DEquation("0");
			m_eq_deriv[5] = new DEquation("1");
			break;
		default:
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "SurfaceMonoAnalytic::setData, unknown ftype");
			return;
		}

		m_eq_deriv[6] = m_eq_deriv[0]->getDerivative(DEquation::op_s);	// dss
		m_eq_deriv[7] = m_eq_deriv[1]->getDerivative(DEquation::op_s);
		m_eq_deriv[8] = m_eq_deriv[2]->getDerivative(DEquation::op_s);

		m_eq_deriv[9]  = m_eq_deriv[0]->getDerivative(DEquation::op_t);	// dst
		m_eq_deriv[10] = m_eq_deriv[1]->getDerivative(DEquation::op_t);
		m_eq_deriv[11] = m_eq_deriv[2]->getDerivative(DEquation::op_t);

		m_eq_deriv[12] = m_eq_deriv[3]->getDerivative(DEquation::op_t);	// dtt
		m_eq_deriv[13] = m_eq_deriv[4]->getDerivative(DEquation::op_t);
		m_eq_deriv[14] = m_eq_deriv[5]->getDerivative(DEquation::op_t);
	}
}

ostream& operator<<(ostream& os, const SurfaceMonoAnalytic *surf){
	os << "FTYPE=" << surf->m_ftype << " F=" << surf->m_str;
	return os;
}

istream& operator>>(istream& is, SurfaceMonoAnalytic *surf){
	string str;
	int ft;
	is >> ws >> ft >> ws >> str;
	if(is){
		str[0] = str[1] = ' ';
		surf->setData(str.c_str(), (SurfaceMonoAnalytic::FTYPE)ft);
	}
	return is;
}

void SurfaceMonoAnalytic::clear()
{
	for(int i = 0; i < 15; i++){
		if(m_eq_deriv[i]){
			delete m_eq_deriv[i];
			m_eq_deriv[i] = nullptr;
		}
	}
	m_valid = false;
}

/// Returns the parameters of the surface for the given point with starting point
const DPoint2d SurfaceMonoAnalytic::getParametersNear(const DPoint3d& point, const DPoint2d& /* near_point */) const
{
	switch(m_ftype){
	case FXY: return DPoint2d(point.x, point.y);
	case FXZ: return DPoint2d(point.x, point.z);
	case FYZ: return DPoint2d(point.y, point.z);
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "SurfaceMonoAnalytic::getParametersNear, unknown ftype");
		return DPoint2d::zero;
	}
}

const DPoint2d SurfaceMonoAnalytic::getParameters(const DPoint3d& point) const
{
	switch(m_ftype){
	case FXY: return DPoint2d(point.x, point.y);
	case FXZ: return DPoint2d(point.x, point.z);
	case FYZ: return DPoint2d(point.y, point.z);
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "SurfaceMonoAnalytic::getParameters, unknown ftype");
		return DPoint2d::zero;
	}
}

