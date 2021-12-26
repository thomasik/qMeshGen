/////////////////////////////////////////////////////////////////////////////
// MeshGenerator3dAdapt.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2011-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR3DADAPT_H__INCLUDED)
#define MESHGENERATOR3DADAPT_H__INCLUDED

#include "TagExtended.h"
#include "DataVector.h"
#include "DataList.h"
#include "DataContainer.h"
#include "DataMatrix.h"
#include "MeshEdge3d.h"

class MeshPoint3d;
class MeshEdge3d;
class MeshFace;
class Metric3dContext;
class MeshContainer3d;
class ControlSpace3dAdaptive;
class SurfaceParametric;
class ControlSpace3d;

/**
 * This class gathers several procedures responsible for adaptation of volume meshes
 * (mesh modification).
 */
class MeshGenerator3dAdapt  
{
protected:
	struct BorderFace {
		MeshPoint3d *p0, *p1, *p2;
		char border;
		TagExtended tags;
		BorderFace(MeshPoint3d *_p0 = 0, MeshPoint3d *_p1 = 0, MeshPoint3d *_p2 = 0, char _border = 0)
			: p0(_p0), p1(_p1), p2(_p2), border(_border) {}
	};
	struct BorderEdge {
		MeshPoint3d *p0, *p1;
		char border;
		BorderEdge(MeshPoint3d *_p0 = 0, MeshPoint3d *_p1 = 0, char _border = 0)
			: p0(_p0), p1(_p1), border(_border) {}
	};
public:
	/// Insert nodes at boundary, for too long edges, return number of inserted points
	static int addBoundaryNodesSplit(Metric3dContext& mc, MeshContainer3d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Collapse inner edge, return resulting point or nullptr - if impossible
	static MeshPoint3d* collapseInnerEdge(Metric3dContext& mc, MeshContainer3d* mesh,
		MeshEdge3d* edge, DataContainer<MeshEdge3d::ActiveEdge> * active_edges = nullptr);
	/// Collapse too short edges, return number of removed edges
	static int collapseInnerEdges(Metric3dContext& mc, MeshContainer3d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Collapse border edge, return resulting point or nullptr - if impossible
	static MeshPoint3d* collapseBoundaryEdge(Metric3dContext& mc, MeshContainer3d* mesh,
		MeshEdge3d* edge, DataContainer<MeshEdge3d::ActiveEdge> * active_edges = nullptr);
	/// Collapse too short edges, return number of removed edges
	static int collapseBoundaryEdges(Metric3dContext& mc, MeshContainer3d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Create new adaptive control space with sizing taken from mesh blocks
	static std::shared_ptr<ControlSpace3dAdaptive> createACSfromMeshBlocks(MeshContainer3d* mesh);
	/// Modify tetrahedral mesh according to the given CS3d, with local operators (edge collapsing and smoothing)
	static bool remeshWithLocalTransformations(Metric3dContext& mc, MeshContainer3d* mesh);
	/// Modify tetrahedral mesh according to the given CS3d, with local mesh cut & replace
	static bool remeshWithLocalCutRegions(Metric3dContext& mc, MeshContainer3d* mesh);
	/// returns total number of blocks in sub-regions with blocks deviating from the given control space
	static int findAdaptRegions(Metric3dContext& mc, MeshContainer3d* mesh,
		DataSimpleList< DataVector<MeshBlock*> > & regions,
		double minq, int min_region_size, int layers);
	/// Search for local reparameterization surfaces (for whole mesh, or only for points with the given tag)
	static int identifyLocalSurfaces(MeshContainer3d* mesh, double tolerance = 0.5, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Search for local reparameterization surface for the given point, returns number of ascribed faces
	static int approximateLocalSurface(MeshContainer3d* mesh, MeshFace* face, double tolerance = 0.2, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Fit surface to cloud of points
	static std::shared_ptr<SurfaceParametric> fitLocalSurface(const DataVector<DPoint3d> & points, double eps2);
};

#endif // !defined(MESHGENERATOR3DADAPT_H__INCLUDED)
