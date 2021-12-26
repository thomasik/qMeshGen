/////////////////////////////////////////////////////////////////////////////
// IteratorFace.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2007-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(ITERATORFACE_H__INCLUDED)
#define ITERATORFACE_H__INCLUDED

class MeshContainer3d;
class MeshFace;
class MeshBlock;

class IteratorFace
{
public:
	/// Constructor
	IteratorFace(const MeshContainer3d* mesh);
	/// Mesh for current iterator
	MeshFace* getFace() { return current_face; }
	/// Get next valid face
	IteratorFace& nextFace();
	/// for end-checking
	bool isValid() const { return findex < fct; }
private:
	MeshFace* current_face;
	MeshBlock* current_block;
	const MeshContainer3d* domain_mesh;
	int bct;
	int bindex;
	int fct;
	int findex;
};

#endif //!defined(ITERATORFACE_H__INCLUDED)
