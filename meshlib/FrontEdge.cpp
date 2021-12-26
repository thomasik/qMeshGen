/////////////////////////////////////////////////////////////////////////////
// FrontEdge.cpp
// Class storing data for front edges used durign triangle->quad conversion
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"
#include "FrontEdge.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "MeshElement.h"
#include "ControlSpace2d.h"

FrontEdge::FrontEdge(MeshEdge2d* edge, MeshPoint2d* point, int lev)
	: m_edge(edge), level(lev), type(2), metric_gradation(1.0), length_squared(1.0)
{
	fedges[EDGE_LEFT] = fedges[EDGE_RIGHT] = nullptr;
	angles[EDGE_LEFT] = angles[EDGE_RIGHT] = 0.0;
	assert(edge);
	edge->getMeshPoint(0)->incIntTag(TagExtended::TAG_QUAD);
	edge->getMeshPoint(1)->incIntTag(TagExtended::TAG_QUAD);
	assert(point);
	index_left = edge->getPointIndex(point);
	edge->setSideTag(point);
}

FrontEdge::~FrontEdge()
{
	if(m_edge){
		m_edge->clearSideTag(index_left);
		m_edge->getMeshPoint(0)->decIntTag(TagExtended::TAG_QUAD);
		assert(m_edge->getMeshPoint(0)->getIntTag(TagExtended::TAG_QUAD) >= 0);
		m_edge->getMeshPoint(1)->decIntTag(TagExtended::TAG_QUAD);		
		assert(m_edge->getMeshPoint(1)->getIntTag(TagExtended::TAG_QUAD) >= 0);
	}
}

//////////////////////////////////////////////////////////////////////
// Classifies type of the front edge (depending on angles between side edges)
FrontEdge* FrontEdge::classify(Metric2dContext& mc)
{
	type = 2;	// 00

	LOG_ASSERT(m_edge!=nullptr);

	// *** angles
	MeshPoint2d *pts[2] = {m_edge->getMeshPoint(index_left), m_edge->getMeshPoint(1-index_left) };
	const DMPoint2d dpts[2] = { 
		pts[0]->getMetricCoordinates(mc), 
		pts[1]->getMetricCoordinates(mc)};

	for(int i = 0; i < 2; i++){
		if(fedges[i]->m_edge == m_edge) angles[i] = 4; // (2*PI);
		else{
			LOG_ASSERT(fedges[i]->m_edge != nullptr);
			// Angle pt . pt_one_side . pt_other_side
			MeshPoint2d* other_point = fedges[i]->m_edge->getOtherPoint(pts[i]);
			angles[i] = (i==0) ? (dpts[1]-dpts[0]).getAngleCos(
				other_point->getMetricCoordinates(mc)-dpts[0])
					: (other_point->getMetricCoordinates(mc)-dpts[1]).getAngleCos(dpts[0]-dpts[1]);
		}
		if(angles[i] < ANGLE_THRESHOLD) type--;
	}
	if(type == 0) --type;
	if(angles[0] < ANGLE_2 || angles[1] < ANGLE_2) type -= 10;

	// *** length penalty (delay too large edges - counted in metric space)
	length_squared = pts[0]->getMetricCoordinates(mc).distance2(pts[1]->getMetricCoordinates(mc));
	if(length_squared > 2.0) type += (int)length_squared;

	//*** control gradation (favour areas with less raPId gradation)
	metric_gradation = mc.getMetricGradation(DPoint2d::average(pts[0]->getCoordinates(), pts[1]->getCoordinates())); // >= 1.0
	type += (int)(metric_gradation - 1.0);

	return this;
}

//////////////////////////////////////////////////////////////////////
// Compares two front edges, necessary for heap ordering
short FrontEdge::compareTo(const FrontEdge *item) const
{
	int key1 = level + type;
	int key2 = item->level + item->type;

	if(key1 > key2) return 1;
	else if(key1 < key2) return -1;
	else return 0;
}

void FrontEdge::clear()
{
	assert(m_edge);
	m_edge->clearSideTag(index_left);
	m_edge->getMeshPoint(0)->decIntTag(TagExtended::TAG_QUAD);
	assert(m_edge->getMeshPoint(0)->getIntTag(TagExtended::TAG_QUAD) >= 0);
	m_edge->getMeshPoint(1)->decIntTag(TagExtended::TAG_QUAD);		
	assert(m_edge->getMeshPoint(1)->getIntTag(TagExtended::TAG_QUAD) >= 0);
	m_edge = nullptr;
}

void FrontEdge::setEdge(MeshEdge2d *edge, MeshPoint2d *point)
{
	if(m_edge) clear();

	assert(edge);
	assert(point);
	m_edge = edge;
	index_left = edge->getPointIndex(point);
	edge->setSideTag(point);
	edge->getMeshPoint(0)->incIntTag(TagExtended::TAG_QUAD);
	edge->getMeshPoint(1)->incIntTag(TagExtended::TAG_QUAD);
}

void FrontEdge::setSideFrontEdge(FrontEdge* v, int side)
{ 
	assert(!m_edge || !v->getEdge() || v->getEdge()->incidentTo(m_edge->getMeshPoint((index_left + side)%2)));
	fedges[side] = v; 
}

ostream& operator<<(ostream& os, const FrontEdge& fe)
{
	os << "FE=[";
	if(!fe.m_edge) return os << "-]";
	else{
		assert(fe.m_edge->getMeshPoint(0)->nonZeroIntTag(TagExtended::TAG_QUAD));
		assert(fe.m_edge->getMeshPoint(1)->nonZeroIntTag(TagExtended::TAG_QUAD));
		return os << fe.m_edge->getMeshPoint(fe.index_left)->getIndex() << ',' 
			<< fe.m_edge->getMeshPoint(1-fe.index_left)->getIndex() << "], lev="
			<< fe.level << ", t=" << fe.type << ", mg=" << fe.metric_gradation 
			<<", ang=[" << fe.angles[0] << ',' << fe.angles[1] << ']';
	}
}
