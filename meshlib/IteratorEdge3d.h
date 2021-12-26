/////////////////////////////////////////////////////////////////////////////
// IteratorEdge3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2007-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(ITERATOREDGE3D_H__INCLUDED)
#define ITERATOREDGE3D_H__INCLUDED

#include "TagExtended.h"

class MeshContainer3d;
class MeshContainer3dSurface;
class MeshEdge3d;
class MeshPoint3d;

class IteratorEdge3d
{
public:
	/// Constructor
	IteratorEdge3d(const MeshContainer3d* mesh, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Constructor
	IteratorEdge3d(const MeshContainer3dSurface* mesh, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Mesh for current iterator
	MeshEdge3d* getEdge() { return current_edge; }
	/// Get next valid mesh2d
	IteratorEdge3d& nextEdge();
	/// for end-checking
	bool isValid() const { return PIndex < pct; }
private:
	MeshEdge3d* current_edge;
	MeshPoint3d* current_point;
	const MeshContainer3d* domain_mesh;
	const MeshContainer3dSurface* domain_mesh_surface;
	int pct;
	int PIndex;
	int ect;
	int eindex;
	TagExtended::TagType tag_type;
	int tag_value;
	bool tag_available;
};

#endif //!defined(ITERATOREDGE3D_H__INCLUDED)
