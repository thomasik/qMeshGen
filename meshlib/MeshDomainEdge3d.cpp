// MeshDomainEdge3d.cpp: implementation of the MeshDomainEdge3d class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshDomainEdge3d.h"
#include "ControlSpace3dAdaptive.h"
#include "DMetric3d.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MeshDomainEdge3d::MeshDomainEdge3d(MeshPoint3d* point1, MeshPoint3d* point2, MeshFace* face)
	: MeshEdge3d(point1, point2, face), m_discretization_valid(false)
{
}

MeshDomainEdge3d::~MeshDomainEdge3d() 
{
	for(size_t i = 0; i < m_discretization_points.countInt(); i++)
		m_discretization_points[i]->preDeleteAll();
	for(size_t i = 0; i < m_freepoints.countInt(); i++)
		m_freepoints[i]->preDeleteAll();
}

/// add free-point for this domain-edge
void MeshDomainEdge3d::addFreePoint(const std::shared_ptr<MeshPoint3d>& fpoint)
{
	m_freepoints.add(fpoint);
}

bool MeshDomainEdge3d::updateACS(CS3dPtr  cs) const
{
	if(!cs || !cs->isAdaptive()) return false;
	if(m_discretization_points.empty()) return false;

	auto acs = cs->getAsAdaptive();
	bool any_change = false;
	const MeshPoint3d* last_point = points[0];
	int dpct = m_discretization_points.countInt();
	for(int i = 0; i <= dpct; i++){
		const MeshPoint3d* point = (i < dpct) ? m_discretization_points[i].get() : points[1];
		const DPoint3d& pt = last_point->getCoordinates();
		const DVector3d dv = point->getCoordinates() - pt;
		double dv_len = dv.length();
		double len = 1.2*dv_len;
		DVector3d nv1, nv2;
		dv.orthonormalVectors(nv1, nv2);
		if(len < ControlSpace2dAdaptive::param_min_length) len = ControlSpace2dAdaptive::param_min_length;
		double nlen = len*ControlSpace2dAdaptive::param_stretch_max_ratio;
		double d[3] = {len, nlen, nlen};
		ControlDataMatrix3d cdm(dv / dv_len, nv1, nv2, d);
		any_change |= acs->setMinControl(ControlDataExtMatrix3dSegment(pt, dv, len, cdm));
		last_point = point;
	}
	return any_change;
}

