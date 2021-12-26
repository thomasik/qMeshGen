/////////////////////////////////////////////////////////////////////////////
// IteratorEdge3d.cpp
// Class describing iterator for browsing 3d edges in the mesh
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2007-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "IteratorEdge3d.h"
#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshEdge3d.h"

IteratorEdge3d::IteratorEdge3d(const MeshContainer3d* mesh, 
							   TagExtended::TagType tag_type, int tag_value) :
	current_edge(nullptr), current_point(nullptr), 
	domain_mesh(mesh), domain_mesh_surface(nullptr),
   	PIndex(-1), ect(0), eindex(-1),
	tag_type(tag_type), tag_value(tag_value), tag_available(false)
{
	assert(domain_mesh);
	pct = mesh->getPointsCount();
	nextEdge();
}

IteratorEdge3d::IteratorEdge3d(const MeshContainer3dSurface* mesh, 
							   TagExtended::TagType tag_type, int tag_value) :
	current_edge(nullptr), current_point(nullptr), 
	domain_mesh(nullptr), domain_mesh_surface(mesh),
   	PIndex(-1), ect(0), eindex(-1),
	tag_type(tag_type), tag_value(tag_value), tag_available(false)
{
	assert(domain_mesh_surface);
	pct = mesh->getPointsCount();
	nextEdge();
}

IteratorEdge3d& IteratorEdge3d::nextEdge()
{
	if(tag_type == TagExtended::TAG_NONE){
		// next edge, don't mind tags
		do{
			if(++eindex < ect){	// next edge for point
				assert(current_point);
				current_edge = current_point->getEdge(eindex);
			}else{
				do{
					if(++PIndex < pct){ // next point, first edge
						current_point = domain_mesh ? domain_mesh->getPointAt(PIndex) : domain_mesh_surface->getPointAt(PIndex);
						ect = current_point->getRank();
						current_edge = (ect > 0) ? current_point->getEdge(eindex = 0) : nullptr;
					}else{ // end of edges
						current_point = nullptr;
						current_edge = nullptr;
						return *this;
					}
				}while(!current_edge);
			}
		}while(current_edge->getPointIndex(current_point) != 0);
	}else{
		// next edge, with at least one point marked with the given tag
		while(true){
			if(++eindex < ect){	// next edge for point
				assert(current_point);
				current_edge = current_point->getEdge(eindex);
			}else{
				do{
					if(++PIndex < pct){ // next point, first edge
						current_point = domain_mesh ? domain_mesh->getPointAt(PIndex) : domain_mesh_surface->getPointAt(PIndex);
						if(current_point->checkIntTag(tag_type, tag_value)){ // ok
							ect = current_point->getRank();
							current_edge = (ect > 0) ? current_point->getEdge(eindex = 0) : nullptr;
						}else{ // ... go to the next point
							current_edge = nullptr;
						}
					}else{ // end of edges
						current_point = nullptr;
						current_edge = nullptr;
						return *this;
					}
				}while(!current_edge);
			}
			if(current_edge->getOtherPoint(current_point)->checkIntTag(tag_type, tag_value)){
				if(current_edge->getPointIndex(current_point) == 0) break; // if both tagged, finish in only one case (out of two)
			}else break;
		}
	}

	return *this;
}

