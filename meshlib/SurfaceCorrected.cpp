// SurfaceCorrected.cpp: implementation of the SurfaceCorrected class.
// Tomasz Jurczyk, 2009-
// Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include "SurfaceCorrected.h"
#include "DEquation.h"

/// Standard constructor
SurfaceCorrected::SurfaceCorrected(std::shared_ptr<const SurfaceParametric> base_surface)
	: m_base_surface(base_surface)
{
	assert(m_base_surface);
}

const DPoint3d SurfaceCorrected::getPoint(const DPoint2d& param) const
{
	DPoint3d result = m_base_surface->getPoint(param);

	// check corrections
	if(!m_corrections.empty()){
		DVector3d cdv = DVector3d::zero;
		double wtotal = 0.0;
		double ksi_max = 0.0;
		for(size_t i = 0; i < m_corrections.countInt(); i++){
			const CorrectionVector& cv = m_corrections[i];
			double dist2 = cv.dm.transformPStoRS(param - cv.center).length2();
			if(dist2 > 1.0) continue;
			//double w = 1.0 - sqrt(dist2);
			//double ksi = w > 0.5 ? 1.0 : 2*w;
			double w = dist2 < 0.25 ? 1.0 : (2.0/(1+dist2)-1.0) * (5.0/3.0);
			//double w = 2.0 / (1 + dist2) - 1.0;
			double ksi = w;
			wtotal += w;
			cdv += cv.dv * w;
			if(ksi > ksi_max) ksi_max = ksi;
		}
		if(wtotal > 0.0)
			result += cdv * (ksi_max/wtotal);
	}

	return result;
}

const DVector3d SurfaceCorrected::getNormalVector(const DPoint2d& param) const
{
	return m_base_surface->getNormalVector(param);
}

const DVector3d SurfaceCorrected::getDerivative(int deriv, const DPoint2d& param) const
{
	return m_base_surface->getDerivative(deriv, param);
}

/// Store XML description to stream
ostream& SurfaceCorrected::storeXML(ostream& os, const string& prefix) const
{
	if(m_corrections.empty()){
		m_base_surface->storeXML(os, prefix);
		return os;
	}

	os << prefix << "<corrected-surface>" << endl;
	m_base_surface->storeXML(os, prefix+"\t");

	os << prefix << "\t<corrections>\n";
	for(size_t k = 0; k < m_corrections.countInt(); k++){
		os << prefix << "\t\t<correction>\n";
		os << prefix << "\t\t\t<center> " << m_corrections[k].center << " </center>\n";
		os << prefix << "\t\t\t<radius> " << DEquation::doubleToString(m_corrections[k].radius, DEquation::v_length) << " </radius>\n";
		os << prefix << "\t\t\t<translation> " << m_corrections[k].dv << " </translation>\n";
		os << prefix << "\t\t</correction>\n";
	}
	os << prefix << "\t</corrections>" << endl;
	return os << prefix << "</corrected-surface>" << endl;
}

/// Insert an additional correction vector
void SurfaceCorrected::insertCorrectionVector(const DPoint2d& center, double radius, const DVector3d& translation)
{
	double len = 1.0 / radius;
	m_corrections.add(
		CorrectionVector(center, radius, 
			DMetric2d(ControlDataMatrix2d(len, 0.0, len), this, center), translation));					
}

SurfaceCorrected::CorrectionVector::CorrectionVector(
	const DPoint2d& _center, double _radius, const DMetric2d& _dm, const DVector3d& _translation)
	: center(_center), radius(_radius), dm(_dm), dv(_translation) { }

/// Returns the parameters of the surface for the given point with starting point
const DPoint2d SurfaceCorrected::getParametersNear(const DPoint3d& point, const DPoint2d& near_point) const
{
	return m_base_surface->getParametersNear(point, near_point);
}

const DPoint2d SurfaceCorrected::getParameters(const DPoint3d& point) const
{
	return m_base_surface->getParameters(point);
}

