/////////////////////////////////////////////////////////////////////////////
// FrontLines.cpp
// Klasa przechowuj¹ca zbiór krawêdzi frontu (podzielonych wed³ug typów)
//	dla algorytmu QMorph
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////

#include "FrontLines.h"
#include "FrontEdge.h"
#include "MeshEdge2d.h"
#include "MeshPoint2d.h"

FrontLines::FrontLines(int part_size)
{
	m_front_edges = new DataContainer<FrontEdge>(part_size);
}

FrontLines::~FrontLines()
{
	delete m_front_edges;
}

//////////////////////////////////////////////////////////////////////
// Funkcja dodaje now¹ krawêdŸ frontu
void FrontLines::addFrontEdge(FrontEdge *fedge)
{
//	fedge->classify();
	m_front_edges->addDataItem(fedge);
}

FrontEdge* FrontLines::getBestFrontEdge() const
{
	if(m_front_edges->countInt() > 0){
		return m_front_edges->getDataAt(0);
	}else{
		return nullptr;
	}
}

void FrontLines::classifyAllEdges(Metric2dContext& mc)
{
	m_front_edges->setHeapOrder(false);
	int count = m_front_edges->countInt();
	for(int i = 0; i < count; i++){
		m_front_edges->getDataAt(i)->classify(mc);
	}
	m_front_edges->setHeapOrder(true);
}

void FrontLines::removeFrontEdge(FrontEdge *fedge)
{
	delete m_front_edges->removeDataItem(fedge->getIndex());
}

FrontEdge* FrontLines::findFrontEdge(const MeshEdge2d *edge, int index) const
{
	int fct = m_front_edges->countInt();
	for(int i = 0; i < fct; i++){
		FrontEdge* fedge = m_front_edges->getDataAt(i);
		if(fedge->getEdge() == edge && fedge->getLeftIndex() == index)
			return fedge;
	}
	return nullptr;
}

bool FrontLines::isValid() const
{
	int fct = m_front_edges->countInt();
	for(int i = 0; i < fct; i++){
		const FrontEdge* fedge = m_front_edges->getDataAt(i);
		assert(fedge->getIndex() == i);
		if(fedge->getSideFrontEdge(0)->getSideFrontEdge(1) != fedge) return false;
		if(fedge->getSideFrontEdge(1)->getSideFrontEdge(0) != fedge) return false;
		MeshEdge2d* edge = fedge->getEdge();
		int index = fedge->getLeftIndex();
		MeshPoint2d* point = edge->getMeshPoint(index);
		if(!edge->isSideTagged(point)) return false;
		if(edge->getMeshElement(index) == nullptr) return false;
		if(edge->getMeshPoint(0)->zeroIntTag(TagExtended::TAG_QUAD)) return false;
		if(edge->getMeshPoint(1)->zeroIntTag(TagExtended::TAG_QUAD)) return false;
	}
	return true;
}

void FrontLines::postpone(FrontEdge* fe)
{
	fe->incLevel(5);
	fe->setType(3);
	updateFrontEdgePosition(fe);
}

