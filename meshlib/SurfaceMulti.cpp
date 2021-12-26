// SurfaceMulti.cpp: implementation of the SurfaceMulti class.
// Tomasz Jurczyk, 2009-
// Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include "SurfaceMulti.h"
#include "SurfacePlane.h"

/// Standard constructor
SurfaceMulti::SurfaceMulti(std::shared_ptr<const SurfaceParametric> base_surface)
	: m_base_surface(base_surface)
{
	assert(m_base_surface);
}

const DPoint3d SurfaceMulti::getPoint(const DPoint2d& param) const
{
	DPoint3d result = DPoint3d::zero;
	double wtotal = 0.0;
	// check subsurfaces:
	for(size_t i = 0; i < m_surfaces.countInt(); i++){
		auto ss = m_surfaces[i];
		double w = ss->withinDomain(param);
		if(w > 0.0){
			DPoint2d local_param = ss->localParam(param);
			DPoint3d local_point =  ss->getSurface()->getPoint(local_param);
			result.add(local_point, w);
			wtotal += w;
		}
	}
	if(wtotal > 1.0){
		result /= wtotal;
	}else if(wtotal < 1.0){
		// adjust with the base surface
		result.add(m_base_surface->getPoint(param), 1.0 - wtotal);
	}

	return result;
}

const DVector3d SurfaceMulti::getNormalVector(const DPoint2d& param) const
{
	DVector3d result = DVector3d::zero;
	double wtotal = 0.0;
	// check subsurfaces:
	for(size_t i = 0; i < m_surfaces.countInt(); i++){
		auto ss = m_surfaces[i];
		double w = ss->withinDomain(param);
		if(w > 0.0){
			DPoint2d local_param = ss->localParam(param);
			DVector3d local_vector =  ss->getSurface()->getNormalVector(local_param);
			result += local_vector * w;
			wtotal += w;
		}
	}
	if(wtotal > 1.0){
		result /= wtotal;
	}else if(wtotal < 1.0){
		// adjust with the base surface
		result += m_base_surface->getNormalVector(param) * (1.0 - wtotal);
	}

	return result;
}

const DVector3d SurfaceMulti::getDerivative(int deriv, const DPoint2d& param) const
{
	DVector3d result = DVector3d::zero;
	double wtotal = 0.0;
	// check subsurfaces:
	for(size_t i = 0; i < m_surfaces.countInt(); i++){
		auto ss = m_surfaces[i];
		double w = ss->withinDomain(param);
		if(w > 0.0){
			DPoint2d local_param = ss->localParam(param);
			DVector3d local_vector =  ss->getSurface()->getDerivative(deriv, local_param);
			result += local_vector * w;
			wtotal += w;
		}
	}
	if(wtotal > 1.0){
		result /= wtotal;
	}else if(wtotal < 1.0){
		// adjust with the base surface
		result += m_base_surface->getDerivative(deriv, param) * (1.0 - wtotal);
	}

	return result;
}

/*
SurfaceParametric* SurfaceMulti::selectSubSurface(const DPoint2d& param, DPoint2d& local_param) const
{
	for(int i = 0; i < m_surfaces.countInt(); i++){
		SubSurface* ss = m_surfaces.getconst(i);
		if(ss->domain.getValue(param.x, param.y) > 0.0){
			local_param.x = ss->trans_param_u.getValue(param.x, param.y);
			local_param.y = ss->trans_param_v.getValue(param.x, param.y);
			return ss->surface;
		}
	}
	// base
	local_param = param;
	return m_base_surface;
}
*/

SurfaceMulti::SubSurface::~SubSurface()
{
	if(domain) delete domain;
	if(reparam) delete reparam;
}

double SurfaceMulti::DomainRect::testWithin(const DPoint2d& param) const
{
	double dx = abs(param.x - middle.x);
	double dy = abs(param.y - middle.y);
	if(dx < width.x && dy < width.y) return 1; // fully within
	if(dx > mwidth.x || dy > mwidth.y) return 0; // fully out

	double mx = (mwidth.x - dx) / (mwidth.x - width.x);
	double my = (mwidth.y - dy) / (mwidth.y - width.y);
	return std::min(mx, my);
}

SurfaceMulti::ReparameterizationPlanar::ReparameterizationPlanar(
	const SurfacePlane& source, const SurfacePlane& destination)
{
	DMatrix3d A;
	A.setColumn(0, destination.baseEu());
	A.setColumn(1, destination.baseEv());
	A.setColumn(2, source.getNormalVector(DPoint2d::zero));
	DVector3d res;
	A.solve(source.baseP() - destination.baseP(), res);
	transf_p = DPoint2d(res.x, res.y);
	A.solve(source.baseP() + source.baseEu() - destination.baseP(), res);
	transf_e0 = DPoint2d(res.x, res.y) - transf_p;
	A.solve(source.baseP() + source.baseEv() - destination.baseP(), res);
	transf_e1 = DPoint2d(res.x, res.y) - transf_p;
}

DPoint2d SurfaceMulti::ReparameterizationPlanar::transform(const DPoint2d& param) const
{
	return transf_p + transf_e0 * param.x + transf_e1 * param.y;
}

/// Insert an additional sub-surface
void SurfaceMulti::insertSubSurface(std::shared_ptr<const SurfaceParametric> surface, 
	Domain* domain, Reparameterization* reparam)
{
	m_surfaces.add(std::make_shared<SubSurface>(surface, domain, reparam));
}

/// Store XML description to stream
ostream& SurfaceMulti::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<multisurface>" << endl;
	os << prefix << "\t<basesurface>" << endl;
	m_base_surface->storeXML(os, prefix+"\t\t");
	os << prefix << "\t</basesurface>\n";
	if(m_surfaces.countInt() == 0) return os;

	os << prefix << "\t<subsurfaces>\n";
	for(size_t i = 0; i < m_surfaces.countInt(); i++){
		os << prefix << "\t\t<subsurface>\n";
		m_surfaces[i]->storeXML(os, prefix+"\t\t\t");
		os << prefix << "\t\t</subsurface>\n";
	}
	os << prefix << "\t</subsurfaces>" << endl;
	return os << prefix << "</multisurface>" << endl;
}

/// Store XML description to stream
ostream& SurfaceMulti::SubSurface::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<domain> to be done </domain>\n";
	os << prefix << "<reparameterization> to be done </reparameterization>\n";
	os << prefix << "<surface>\n";
	surface->storeXML(os, prefix+"\t");
	os << prefix << "</surface>\n";
	return os;
}
