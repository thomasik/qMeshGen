/////////////////////////////////////////////////////////////////////////////
// MeshGenerator3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR3D_H__INCLUDED)
#define MESHGENERATOR3D_H__INCLUDED

#include <memory>

#include "DPoint.h"
#include "MeshData.h"
#include "DataVector.h"
#include "Metric3dContext.h"
#include "MeshContainer3d.h"
#include "OctTree.h"
#include "TagExtended.h"

class MeshPoint3d;
class MeshEdge3d;
class MeshFace;
class MeshBlock;
class MeshTetrahedron;
class MeshViewSet;
/**
 * This class gathers several procedures responsible for discretizing volumes 
 * into 3D tetrahedral meshes (triangulation and smoothing).
 */
class MeshGenerator3d  
{
protected:
	struct TFace{
		TFace() : face(nullptr), face_index(0), tetrahedron(nullptr) {}
		MeshFace *face;                 // common face
		int face_index;                 // orientation of this face for new tetrahedron
		MeshTetrahedron* tetrahedron;   // new tetrahedron
	};
public:
	/// Split the given edge (and all adjacent tetrahedra) in a given point (within the edge)
	static MeshPoint3d* splitEdgeSimple(MeshContainer3d* mesh, MeshEdge3d* edge, const DPoint3d& pt);
	/// Split the given face (and two adjacent tetrahedra) in a given point (within the face)
	static MeshPoint3d* splitFaceSimple(MeshContainer3d* mesh, MeshFace* face, const DPoint3d& pt);
	/// Create triangulation of a set of points
	static MeshContainer3d* triangulatePoints(
		std::shared_ptr<DataVector<std::shared_ptr<MeshPoint3d>>> vec);
	/// Creates special decomposition mesh
	static bool createDecompositionMesh(MeshContainer3d* mesh);
	/// Auto-triangulate modeled domain
	static bool autoTriangulate(MeshContainer3d* boundary, int sm_count = 2);
	/// Auto-triangulate modeled domain
	static bool autoTriangulateRemesh(MeshContainer3d* boundary, int sm_count = 2);
	/// Performs the second stage of triangulation by adding inner nodes in order to satisfy the control space requirements
	static int addInnerNodes(Metric3dContext& mc, MeshContainer3d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Transforms 2 tetrahedra into 2 tetrahedra by switching a common edge (on boundary), returns swapped edge if successful
	static MeshEdge3d* swap22(Metric3dContext& mc, MeshContainer3d* mesh, MeshEdge3d* edge, 
		MeshTetrahedron** created_tetrahedrons = nullptr);
	/// Transforms 4 tetrahedra into 4 tetrahedra by switching a common edge
	static bool swap44(Metric3dContext& mc, MeshContainer3d* mesh, MeshFace* face1, MeshEdge3d* edge, int phase);
	/// Checks whether the 3-2 transformation is possible for this configuration
	static bool swap32possible(Metric3dContext& mc, MeshEdge3d* edge);
	/// Transforms 3 tetrahedra into 2 tetrahedra by switching a common edge
	static bool swap32(Metric3dContext& mc, MeshContainer3d* mesh, MeshEdge3d* edge, 
		MeshTetrahedron** created_tetrahedrons = nullptr, 
		int quality_condition = MeshData::SWAP3_ALWAYS);
	/// Checks whether the 3-2 transformation is advantageous for this configuration and quality condition
	static bool swap32Advantageous(Metric3dContext& mc, MeshEdge3d* edge, 
		int quality_condition = MeshData::SWAP3_MIN_VOLUME);
	/// Transforms 2 tetrahedra into 3 tetrahedra by switching a common face
	static bool swap23(Metric3dContext& mc, MeshContainer3d* mesh, MeshFace* face, 
		MeshTetrahedron** created_tetrahedrons = nullptr, 
		int quality_condition = MeshData::SWAP3_ALWAYS);
	/// Checks whether the 2-3 transformation is possible for this configuration
	static bool swap23possible(Metric3dContext& mc, MeshFace* face);
	/// Checks whether the 2-3 transformation is advantageous for this configuration and quality condition
	static bool swap23Advantageous(Metric3dContext& mc, MeshFace* face, 
		int quality_condition = MeshData::SWAP3_MIN_VOLUME);
	static bool addPointToTriangulation(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* point, MeshTetrahedron* containing_tetrahedron = nullptr,
		MeshTetrahedron* improved_tetrahedron = nullptr, bool improving_only = false);
	/// Performs Delaunay3D retriangulation using the swap criterion
	static bool addPointToTriangulationBySwap(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* point, MeshTetrahedron* containing_tetrahedron = nullptr, 
		bool allow_swap = true);
	/// Performs Delaunay3D retriangulation using the empty cavity criterion
	static bool addPointToTriangulationBySphere(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* point, MeshTetrahedron* containing_tetrahedron = nullptr,
		MeshTetrahedron* improved_tetrahedron = nullptr, bool improving_only = false);
	/// Calls the triangulation (tetrahedralization) procedures for all domain-blocks
	static int createTetrahedra(MeshContainer3d* boundary);
	/// Prepares the requested boundary description for Delaunay 3D-tetrahedralization procedures
	static int prepareBoundaryMesh(MeshContainer3d* boundary);
	/// Prepares the requested boundary description for Delaunay 3D-tetrahedralization procedures
	static int prepareBoundarySurfaceMesh(MeshContainer3d* boundary);
	/// Restores normal metric configuration after special metric transformation
	static void iterativeTetrahedraSwap(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d *init_point, int max_layers, int max_swaps = 10000, bool local_metric = true,
		bool recalculate_quality = false);
public:
	/// Maximum number of retriangulation for automatic control space 3d adjustment
	static int param_max_auto_retriangulation_count;
	/// Method used for Delaunay retriangulation
	static int param_triangulation_type;
	/// Method of tetrahedra improvement (where to insert new node)
	static double param_quality_improvement;
	/// Quality threshold for tetrahedra improvement
	static double param_quality_threshold;
	/// Volume threshold for tetrahedra (swap, insertion)
	static double param_volume_threshold;
	/// Special metric ratio (for non-desired edges)
	static double param_special_metric_ratio;
	/// Method used for boundary triangulation
	static int param_boundary_triangulation_method;
};

#endif // !defined(MESHGENERATOR3D_H__INCLUDED)
