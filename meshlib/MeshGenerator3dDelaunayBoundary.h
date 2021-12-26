/////////////////////////////////////////////////////////////////////////////
// MeshGenerator3dDelaunayBoundary.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR3DDELAUNAYBOUNDARY_H__INCLUDED)
#define MESHGENERATOR3DDELAUNAYBOUNDARY_H__INCLUDED

#include <memory>

#include "DPoint.h"
#include "MeshData.h"
#include "DataVector.h"
#include "Metric3dContext.h"
#include "MeshContainer3d.h"
#include "OctTree.h"
#include "DataHashTable.h"

class DMTriangle3d;
class MeshPoint3d;
class MeshEdge3d;
class MeshFace;
class MeshBlock;
class MeshTetrahedron;
class MeshViewSet;
class MeshContainer3dSurface;

/**
 * This class gathers several procedures responsible for discretizing volumes 
 * into 3D tetrahedral meshes (triangulation and smoothing).
 */
class MeshGenerator3dDelaunayBoundary  
{
protected:
	enum PipeType { PIPE_START_POINT, PIPE_FACE_POINT, PIPE_BLOCK_POINT, 
		PIPE_FACE_EDGE, PIPE_BLOCK_EDGE, PIPE_BLOCK_FACE };
	enum ReductionType { REDUCE_FAILED, REDUCE_CLEAN, REDUCE_TRANSFORM, 
		REDUCE_POINT_CLEAN, REDUCE_POINT_TRANSFORM, REDUCE_CUT_SPLIT};
	struct PipeNode {
		/// Simple constructor
		PipeNode(PipeType etype, void* elem1, void* elem2, 
				 const DPoint3d& cross_pt, PipeNode* n = nullptr) : 
			element_type(etype), element1(elem1), element2(elem2), pt(cross_pt), next(n) {}
		PipeType element_type;		
		/// Pointer to the "containing" element (face or block)
		void* element1;
		/// Pointer to the "ending" element (face, edge or point)
		void* element2;
		/// 3D Point (e.g crossing point if the "ending" element is an edge or a face)
		DPoint3d pt;
		/// Pointer to the next Pipe-node
		PipeNode* next;
	};
	struct PointNode {
		PointNode(const DMPoint3d& _pt = DMPoint3d::zero, double v = 1.0) : pt(_pt), visibility(v) {}
		DMPoint3d pt;
		double visibility;
		bool operator<(const PointNode& rhs) const { return visibility > rhs.visibility; } // decreasing
	};
	struct FaceNode {
		FaceNode() {}
		FaceNode(const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2)
			: pt(pt0), v0(pt1-pt0), v1(pt2-pt0)
		{
			middle = pt + v0*0.5 + v1*0.5;
			double dist_0 = v0.length2();
			double dist_1 = v1.length2();
			double dist_2 = (pt2-pt1).length2();
			max_dist2 = 2 * std::max(std::max(dist_0, dist_1), dist_2);
			min_dist2 = 0.5 * std::min(std::min(dist_0, dist_1), dist_2);
		}
		DPoint3d pt;
		DVector3d v0, v1;
		DPoint3d middle;
		double max_dist2;
		double min_dist2;
	};
	struct CandidateEdge {
		MeshEdge3d* edge;
		double direction;
		bool operator>(const CandidateEdge& ce) const { return direction > ce.direction; }
	};
	struct CandidateFace {
		MeshFace* face;
		double direction;
		bool operator>(const CandidateFace& cf) const { return direction > cf.direction; }
	};
	struct MissingFace{
		MeshFace* face;
		MeshPoint3d *m_pt[3];
		ControlDataMatrix3d special_metric;
		bool valid_metric;
		int depth;
	public:
		MissingFace(MeshFace* f = nullptr, int d = 0) 
			: face(f), valid_metric(false), depth(d) 
		{ m_pt[0] = m_pt[1] = m_pt[2] = nullptr; }
		MissingFace(MeshFace* f, MeshPoint3d* mpt0, MeshPoint3d* mpt1, MeshPoint3d* mpt2, int d = 0) 
			: face(f), valid_metric(false), depth(d) 
		{ m_pt[0] = mpt0; m_pt[1] = mpt1; m_pt[2] = mpt2; }
		const ControlDataMatrix3d& getSpecialMetric();
		bool adjacentToPoints(MeshPoint3d* pt1, MeshPoint3d* pt2) const;
	};
	struct MissingEdge{
		MeshEdge3d* edge;
		MeshPoint3d *m_pt0, *m_pt1;
		ControlDataMatrix3d special_metric;
		bool valid_metric;
		int depth;
		DataVector<MeshPoint3d*> inserted_points;
	public:
		MissingEdge(MeshEdge3d* e = nullptr, int d = 0) : edge(e), m_pt0(nullptr), m_pt1(nullptr), valid_metric(false), depth(d) {}
		MissingEdge(MeshEdge3d* e, MeshPoint3d* mpt0, MeshPoint3d* mpt1, int d = 0) 
			: edge(e), m_pt0(mpt0), m_pt1(mpt1), valid_metric(false), depth(d) {}
		const ControlDataMatrix3d& getSpecialMetric();
	};
public:
	/// Try to retriangulate the cavity for recovery boundary edge
	static bool retriangulateCavity(Metric3dContext& mc, MeshContainer3d* mesh, 
			const DataVector<MeshBlock*> & crossing_blocks,
			const DataHashTable<MeshBlock*> & hblocks,
			const DataHashTable<MeshPoint3d*> & cavity_fpoints,
			const DataVector<MeshEdge3d*> & cavity_fedges,
			const DataVector<MeshGenerator3dDelaunayBoundary::MissingEdge> & local_edges,
			const DataVector<int> local_faces,
			const DataVector<MeshGenerator3dDelaunayBoundary::MissingFace> & missing_faces);
	/// Try recover boundary edge (and adjacent faces) by retrieving local cavity and remeshing
	static bool recoverBoundaryWithLocalCavityRemeshing(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* pt1, MeshPoint3d* pt2, const DataVector<MissingFace> & missing_faces);
	/// Try recover boundary face (single) by retrieving local cavity and remeshing
	static bool recoverBoundaryWithLocalCavityRemeshing(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* pt1, MeshPoint3d* pt2, MeshPoint3d* pt3);
	/// Match missing edge (and adjacent faces) by matching with border-entities (with transformations on both sides)
	static int matchBoundaryEdge(Metric3dContext& mc, const MissingEdge& me,
		DataHashTableKeyValue<MeshPoint3d*,MeshPoint3d*> & bmpoints);
	/// Generate boundary constrained tetrahedral mesh
	static MeshContainer3d* createBoundaryConstrainedMesh(Metric3dContext& mc, 
		MeshContainer3dSurface* surface_mesh, MeshDomainVolume* mdv);
	/// Tries to reduce the rank of the edge for swap44
	static bool reduceEdgeRankForSwap44(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshEdge3d* edge, MeshFace* first_face, MeshFace* second_face);
	/// Tries to recover the missing boundary face with these three points
	static int recoverBoundaryFace(Metric3dContext& mc, MeshContainer3d* mesh, 
		const DataVector<MissingFace> & missing_faces, 
		MeshPoint3d* pt1, MeshPoint3d* pt2, MeshPoint3d* pt3, int phase = 0);
	static MeshPoint3d* cutFirstPipeNode(Metric3dContext& mc, MeshContainer3d* mesh,
		MeshPoint3d* pt1, MeshPoint3d* pt2, PipeNode* pipe);
	/// Tries to recover the missing boundary edge between both points
	static int recoverBoundaryEdge(Metric3dContext& mc, MeshContainer3d* mesh, 
		MissingEdge& me, const DataVector<MissingFace> & missing_faces, int phase = 0);
	/// Constrains the Delaunay triangulation of boundary nodes (restores all boundary edges and faces, identifies valid tetrahedra and removes the obsolete ones)
	static bool constrainToBorder(Metric3dContext& mc, MeshContainer3d* mesh, 
		const MeshDomainVolume* mdv, MeshContainer3dSurface* surface_mesh);
	/// Insert auxiliary points for boundary meshing
	static int refineMeshForBoundaryTriangulation(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* point, const MeshContainer3d* boundary_mesh, int bindex);
	/// Insert auxiliary points for boundary meshing
	static int refineMeshForBoundaryTriangulationByEdges(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* point, const OctTreeMeshPoints& octree_points);
	/// Insert auxiliary points for boundary meshing
	static int refineMeshForBoundaryTriangulation(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* point, const OctTreeMeshPoints& octree_points);
	/// Performs the first stage of Delaunay tetrahedralization of boundary nodes (only)
	static MeshContainer3d* createInitialMesh(int points_count, const DBox& bbox, CS3dPtr cs);
	/// Tries to transform (improve locally) mesh by inserting new nodes within the worst block from Pipe
	static int insertPointInWorstBlock(Metric3dContext& mc, MeshContainer3d* mesh,
		MeshPoint3d* start_point, PipeNode* Pipe, const DataVector<FaceNode> &missing_faces);
	/// Tries to transform (improve locally) mesh by inserting new nodes in the vicinity of missing edge
	static int insertPointAtLongestPipeEdges(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* start_point, PipeNode* Pipe, const DataVector<MissingFace> & missing_faces);
	/// Tries to transform (improve locally) mesh by inserting new nodes in the vicinity of missing edge
	static bool insertPointAtLongestAdjacentEdge(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d* point, const DataVector<MissingFace> & missing_faces);
	/// Tries to fix (remove from boundary) additional boundary nodes
	static bool fixAdditionalBoundaryNodes(Metric3dContext& mc, MeshContainer3d* mesh,
		DataVector<MeshContainer3d::AdditionalBoundaryNode> & inserted_boundary_points);
private:
	/// Find good center of star-shaped cavity (if possible)
	static bool findStarCenter(const DataVector<DMTriangle3d> & orient_faces, DMPoint3d& center);
	/// Tries to recover the missing boundary edge using special metric approach
	static int recoverBoundaryEdgeWithSpecialMetricSwap(Metric3dContext& mc, 
		MeshContainer3d* mesh, MeshPoint3d* pt1, MeshPoint3d* pt2, int max_layers);
	/// Gathers segments and blocks crossing the missing face
	static void gatherCrossingEdgesForFace(Metric3dContext& mc, MeshEdge3d* edges[], 
		MeshPoint3d* points[], DataVector<MeshEdge3d*> & cross_edges);
	/// Tries to recover the missing boundary face using special metric approach
	static int recoverBoundaryFaceWithSpecialMetricSwap(Metric3dContext& mc, 
		MeshContainer3d* mesh, MeshPoint3d* pt1, MeshPoint3d* pt2, MeshPoint3d* pt3,
		int max_layers);
	/// Tries to recover boundary edge by introducing additional (boundary) nodes
	static bool insertBoundaryPointForEdgeRecovery(Metric3dContext& mc, MeshContainer3d* mesh, MissingEdge& me,
		DataVector<MissingFace> & missing_faces, DataVector<MissingEdge> & new_missing_edges, 
		DataVector<MeshContainer3d::AdditionalBoundaryNode> & inserted_boundary_points, int phase);
	/// Tries to recover boundary face by introducing additional (boundary) nodes
	static bool insertBoundaryPointForFaceRecovery(Metric3dContext& mc, 
		MeshContainer3d* mesh, MissingFace& mf, DataVector<MissingEdge> & missing_edges, 
		const DataVector<MissingFace> & missing_faces, DataVector<MissingFace> & new_missing_faces, 
		DataVector<MeshContainer3d::AdditionalBoundaryNode> & inserted_boundary_points, int phase);
	/// Tries to transform (improve locally) mesh by inserting new node (near face)
	static int insertPointForFaceRecovery(Metric3dContext& mc, MeshContainer3d* mesh,
		const DPoint3d& pt, const DataVector<MissingFace> &missing_faces, MeshTetrahedron* near_block);
	/// Creates ViewSet representing the given recovery-Pipe
	static MeshViewSet* createPipeViewSet(PipeNode* Pipe, MeshPoint3d* start_point, MeshPoint3d* last_point);
	/// Returns list of Pipe-nodes (starting from the given node, nullptr => from STARTING_POINT)
	static PipeNode* buildPipeWithInsertingNodes(Metric3dContext& mc, MeshContainer3d* mesh, 
		const DataVector<MissingFace> & missing_faces, 
		MeshPoint3d* &end_point, MeshPoint3d* &start_point, int phase = 0);
	/// Returns list of Pipe-nodes (starting from the given node, nullptr => from STARTING_POINT)
	static PipeNode* buildPipe(Metric3dContext& mc, MeshPoint3d* end_point, PipeNode* last_node = nullptr, 
		MeshPoint3d* start_point = nullptr);
	/// Tries to join two neighbouring Pipe-segments (through the local reorganization of the mesh)
	static bool reducePipeNodes(Metric3dContext& mc, MeshContainer3d* mesh, 
		PipeNode* last_node, PipeNode* node, int phase, bool& invalid_Pipe);
	/// Tries to transform vicinity of Pipe (through the local reorganization of the mesh)
	static bool transformPipeVicinity(Metric3dContext& mc, MeshContainer3d* mesh, PipeNode* Pipe,
		DataVector<int> &old_faces, DataVector<int> &old_edges);
	/// Tries to transform two neighbouring Pipe-segments (through the local reorganization of the mesh)
	static bool transformPipeNodes(Metric3dContext& mc, MeshContainer3d* mesh, PipeNode* node);
	/// Returns next fragment of the mesh crossing the edge to be restored (starting from the face)
	static PipeNode* nextPipeFromFace(Metric3dContext& mc, MeshBlock* last_block, 
		MeshFace* last_face, const DPoint3d& cross_point, MeshPoint3d* end_point);
	/// Returns next fragment of the mesh crossing the edge to be restored (starting from the edge)
	static PipeNode* nextPipeFromEdge(Metric3dContext& mc, MeshEdge3d* last_edge, 
		const DPoint3d& cross_point, MeshPoint3d* end_point);
	/// Returns next fragment of the mesh crossing the edge to be restored (starting from the vertex)
	static PipeNode* nextPipeFromPoint(Metric3dContext& mc, MeshPoint3d* last_point, 
		MeshPoint3d* end_point);
public:
	/// Whether to insert additional boundary nodes for recovery of boundary edges
	static int param_boundary_nodes_for_recovery;
	/// Dimension ratio for the initial mesh size
	static double param_initial_mesh_size_ratio;
	/// Volume threshold for block-face-edge checking during boundary recovery
	static double param_recovery_volume_threshold;
	/// Special metric ratio (for non-desired edges)
	static double param_special_metric_ratio;
	/// Refined algorithm for incremental constrained Delaunay - RefCD
	static double param_refined_ICD_ratio;
};

#endif // !defined(MESHGENERATOR3DDELAUNAYBOUNDARY_H__INCLUDED)
