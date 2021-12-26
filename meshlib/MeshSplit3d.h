/////////////////////////////////////////////////////////////////////////////
// MeshSplit3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2010-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHSPLIT3D_H__INCLUDED)
#define MESHSPLIT3D_H__INCLUDED

class MeshContainer2d;
class MeshContainer3d;
class MeshPoint3d;

#include "TagExtended.h"

/**
 * This class implements several methods for split/join 3d tetrahedral meshes
 */
class MeshSplit3d
{
public:
	/// Partition mesh (remove marked blocks and transfer them to other mesh)
	static MeshContainer3d* splitMeshByBlocks(MeshContainer3d* mesh, TagExtended::TagType tag_type, int tag_value);
	/// Merge two meshes into one (mesh_add will be erased!)
	static bool mergeMeshes(MeshContainer3d* mesh_base, MeshContainer3d* &mesh_add);
private:
	/// Structure for storing basic face information
	struct FaceInfo {
	public:
		FaceInfo(MeshPoint3d* _pt0, MeshPoint3d* _pt1, MeshPoint3d* _pt2, char _btype) 
			: pt0(_pt0), pt1(_pt1), pt2(_pt2), btype(_btype) {}
	public:
		MeshPoint3d* pt0;
		MeshPoint3d* pt1;
		MeshPoint3d* pt2;
		char btype;
	};
	/// Structure for storing basic edge information
	struct Edge3dInfo {
	public:
		Edge3dInfo(MeshPoint3d* _pt0, MeshPoint3d* _pt1) : pt0(_pt0), pt1(_pt1) {}
	public:
		MeshPoint3d* pt0;
		MeshPoint3d* pt1;
	};
};

#endif // !defined(MESHSPLIT3D_H__INCLUDED)
