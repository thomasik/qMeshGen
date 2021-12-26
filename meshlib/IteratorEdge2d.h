/////////////////////////////////////////////////////////////////////////////
// IteratorEdge2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2005-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(ITERATOREDGE2D_H__INCLUDED)
#define ITERATOREDGE2D_H__INCLUDED

class MeshContainer2d;
class MeshEdge2d;
class MeshPoint2d;

class IteratorEdge2d
{
public:
	/// Constructor
	IteratorEdge2d(const MeshContainer2d* mesh);
	/// Mesh for current iterator
	MeshEdge2d* getEdge() { return current_edge; }
	/// Get next valid mesh2d
	IteratorEdge2d& nextEdge();
	/// for end-checking
	bool isValid() const { return PIndex < pct; }
private:
	MeshEdge2d* current_edge;
	MeshPoint2d* current_point;
	const MeshContainer2d* domain_mesh;
	int pct;
	int PIndex;
	int ect;
	int eindex;
};

#endif //!defined(ITERATOREDGE2D_H__INCLUDED)
