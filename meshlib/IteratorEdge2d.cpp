/////////////////////////////////////////////////////////////////////////////
// IteratorEdge2d.cpp
// Class describing iterator for browsing 2d edges in the mesh
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2005-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "IteratorEdge2d.h"
#include "MeshContainer2d.h"
#include "MeshEdge2d.h"

IteratorEdge2d::IteratorEdge2d(const MeshContainer2d* mesh) :
	current_edge(nullptr), current_point(nullptr), domain_mesh(mesh), 
   	PIndex(-1), ect(0), eindex(-1)
{
	assert(domain_mesh);
	pct = mesh->getPointsCount();
	nextEdge();
}

IteratorEdge2d& IteratorEdge2d::nextEdge()
{
	do{
		if(++eindex < ect){	// next edge for point
			assert(current_point);
			current_edge = current_point->getEdge(eindex);
		}else{
			while(++PIndex < pct){ // next point, first edge
				current_point = domain_mesh->getPointAt(PIndex);
				ect = current_point->getRank();
				if(ect > 0){
					current_edge = current_point->getEdge(eindex = 0);
					break;
				}
			}
			if(PIndex >= pct){ // end of edges
				current_point = nullptr;
				current_edge = nullptr;
				return *this;
			}
		}
	}while(current_edge->getPointIndex(current_point) != 0);

	return *this;
}
