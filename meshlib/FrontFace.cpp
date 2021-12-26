/////////////////////////////////////////////////////////////////////////////
// FrontFace.cpp
// Class storing data for front faces used during frontal boundary-constrained meshing.
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "FrontFace.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshFace.h"
#include "MeshBlock.h"
#include "Metric3dContext.h"
#include "DTriangle.h"
#include "MeshViewSet.h"

FrontFace::FrontFace(MeshFace* face, bool oriented, int lev)
		: m_face(face), m_oriented(oriented), m_neighbors(3), m_angles(3), m_level(lev), m_min_angle_index(-1)
{
	if(m_face) m_face->setPtrTag(m_oriented ? TagExtended::TAG_FRONT_1 : TagExtended::TAG_FRONT_2, this);
}

FrontFace::~FrontFace()
{
	if(m_face) m_face->removeTag(m_oriented ? TagExtended::TAG_FRONT_1 : TagExtended::TAG_FRONT_2);
}

//////////////////////////////////////////////////////////////////////
/// Gathers incidency data for this face (neighbours and angles)
FrontFace* FrontFace::init(Metric3dContext& mc)
{
	assert(m_face);
	if(!m_face) return nullptr;

	mc.countMetricAtPoint(m_face->getMiddlePoint());

	// find neighbours
	int ect = m_face->getEdgeCount();
	m_neighbors.clear();
	m_angles.clear();
	for(int j = 0; j < ect; j++){
		MeshEdge3d* edge = m_face->getEdge(j);
		// find next face
		double best_angle = 3*PI;
		FrontFace* best_fface = nullptr;
		for(int k = 0; k < edge->getFaceCount(); k++){
			// -> find face
			MeshFace* other_mface = edge->getFaceAt(k);
			// ... different than this one
			if(other_mface == m_face) continue;
			// ... part of front
			FrontFace* fface = (FrontFace*)other_mface->getPtrTag(TagExtended::TAG_FRONT_1);
			if(fface){ // oriented == true
				double angle = calculateDihedralAngle(fface, mc);
				if(angle < best_angle){
					best_angle = angle;
					best_fface = fface;
				}
			}
			fface = (FrontFace*)other_mface->getPtrTag(TagExtended::TAG_FRONT_2);
			if(fface){ // oriented == false
				double angle = calculateDihedralAngle(fface, mc);
				if(angle < best_angle){
					best_angle = angle;
					best_fface = fface;
				}
			}
		}
		m_neighbors.add(best_fface);
		m_angles.add(best_angle);
	}
	
	m_min_angle_index = 0;
	for(int j = 1; j < ect; j++)
		if(m_angles[j] < m_angles[m_min_angle_index])
			m_min_angle_index = j;

	return this;
}

//////////////////////////////////////////////////////////////////////
// Compares two front edges, necessary for heap ordering
short FrontFace::compareTo(const FrontFace *item) const
{
	assert(m_min_angle_index > -1 && item->m_min_angle_index > -1);
	double key1 = m_angles[m_min_angle_index];
	double key2 = item->m_angles[item->m_min_angle_index];

	if(key1 > key2) return 1;
	else return (key1 < key2) ? -1 : 0;
}

void FrontFace::clear()
{
	assert(m_face);
	if(m_face) m_face->removeTag(m_oriented ? TagExtended::TAG_FRONT_1 : TagExtended::TAG_FRONT_2);
	m_face = nullptr;
}

void FrontFace::setFace(MeshFace *face, bool oriented)
{
	if(m_face) clear();

	m_face = face;
	m_oriented = oriented;

	if(m_face) face->setPtrTag(m_oriented ? TagExtended::TAG_FRONT_1 : TagExtended::TAG_FRONT_2, this);
}

MeshPoint3d* FrontFace::getPoint(int i) const
{
	assert(m_face);
	if(!m_face) return nullptr;

	return m_oriented ? m_face->getPoint(i) : m_face->getPoint(m_face->getEdgeCount() + 1 - i);
}

/// Returns other point from the (min-angle) incident front face
MeshPoint3d* FrontFace::getIncidentPoint(int i) const
{
	assert(i >= 0 && i < m_angles.countInt());
	FrontFace* fface = m_neighbors[i];
	MeshEdge3d* edge = m_face->getEdge(i);
	return fface->m_face->getOtherPoint(edge->getMeshPoint(0), edge->getMeshPoint(1));
}

/// Returns face of the i-th incident front-face
MeshFace* FrontFace::getIncidentFace(int i) const
{
	assert(i >= 0 && i < m_angles.countInt());
	return m_neighbors[i]->getFace();
}

/// Check if the given front face is among the neighbors
bool FrontFace::adjacentTo(FrontFace* fface) const
{
	return m_neighbors.contains(fface);
}

void FrontFace::gatherNeighbors(DataVector<FrontFace*> &ffaces)
{
	for(size_t i = 0; i < m_neighbors.countInt(); i++)
		ffaces.addIfNew(m_neighbors[i]);
}

/// Calculates dihedral angle between front faces
double FrontFace::calculateDihedralAngle(FrontFace* ff, Metric3dContext& mc) const
{
	MeshPoint3d* mpts[] = {
		m_face->getPoint(m_oriented ? 0 : 1),
		m_face->getPoint(m_oriented ? 1 : 0),
		m_face->getPoint(2)
	};

	const DMVector3d nv0 = DMVector3d::crossProduct(
		mpts[0]->getMetricCoordinates(mc),
		mpts[1]->getMetricCoordinates(mc),
		mpts[2]->getMetricCoordinates(mc));

	MeshPoint3d* other_mpts[] = {
		ff->getFace()->getPoint(ff->isOriented() ? 0 : 1),
		ff->getFace()->getPoint(ff->isOriented() ? 1 : 0),
		ff->getFace()->getPoint(2)
	};

	const DMVector3d nv1 = DMVector3d::crossProduct(
		other_mpts[0]->getMetricCoordinates(mc),
		other_mpts[1]->getMetricCoordinates(mc),
		other_mpts[2]->getMetricCoordinates(mc));

	double angle = PI - nv0.getAngle(nv1);

	MeshPoint3d* other_mpt = nullptr;
	for(int i = 0; i < 3; i++){
		if(other_mpts[i] != mpts[0] && other_mpts[i] != mpts[1] && other_mpts[i] != mpts[2]){
			other_mpt = other_mpts[i];
			break;
		}
	}
	assert(other_mpt);
	bool not_oriented = DMTriangle3d::orient3d(
		mpts[0]->getMetricCoordinates(mc),
		mpts[1]->getMetricCoordinates(mc),
		mpts[2]->getMetricCoordinates(mc),
		other_mpt->getMetricCoordinates(mc)) < 0.0;

	if(not_oriented) angle += PI;

	if(not_oriented){
		MeshViewSet* set = new MeshViewSet;
		// ... points
		for(int i = 0; i < 3; i++)
			set->addPoint(mpts[i]->getCoordinates(), 0, i);
		set->addPoint(other_mpt->getCoordinates(), 0, 4);
		// ... faces
		set->addFace(m_face, 0, 0.9, isOriented());
		set->addFace(ff->getFace(), 1, 0.9, ff->isOriented());
		// ... tetrahedra
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Dihedral angle=" << angle);
		SHOW_MESH("Dihedral angle", set);
	}

	return angle;

}

double FrontFace::getAngle(const FrontFace* ff) const { 
	for(int i = 0; i < m_neighbors.countInt(); i++)
		if(m_neighbors[i] == ff) return m_angles[i]; 
	return 0.0;
}
