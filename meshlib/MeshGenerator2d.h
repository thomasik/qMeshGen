/////////////////////////////////////////////////////////////////////////////
// MeshGenerator2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR2D_H__INCLUDED)
#define MESHGENERATOR2D_H__INCLUDED

#include "MeshData.h"
#include "DataStatistics.h"
#include "Metric2dContext.h"
#include "TagBorder.h"
#include "TagExtended.h"
#include "SurfaceParametric.h"

class MeshContainer3d;
class MeshContainer2d;
class MeshPoint2d;
class MeshEdge2d;
class Curve2dParametric;
class MeshTriangle2d;
class DRect;
class ControlSpace2d;
class ControlSpace2dAdaptive;

/**
 * This class gathers several procedures responsible for discretizing surfaces 
 * into 2D triangular meshes (triangulation and smoothing).
 */
class MeshGenerator2d  
{
public:
	/// Create triangulation of polygon
	static MeshContainer2d* triangulatePoly( const DataVector<DPoint2d> & poly, SurfaceConstPtr surface);
	/// Create triangulation of planar 3d polygon
	static MeshContainer2d* triangulatePoly( const DataVector<DPoint3d> & poly );
	/// Scatter to blocks/patches/contours, returns number of points remaining...
	static int scatterFreePoints(MeshContainer3d* boundary);
	/// Collapse too short edges
	static bool collapseEdges(Metric2dContext& mc, MeshContainer2d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Creates special decomposition mesh
	static bool createDecompositionMesh(MeshContainer2d* mesh);
	/// Creates new adaptive control space of selected type
	static std::shared_ptr<ControlSpace2dAdaptive> createNewControlSpace(std::shared_ptr<const SurfaceParametric> surface, const DRect& rect);
	/// Checks if control space is consistent with close boundary edges (and adjusts if necessary)
	static bool checkControlForCloseBoundaryEdges(Metric2dContext& mc, MeshContainer2d* mesh);
	/// Checks if control space is consistent with boundary edges (and adjusts if necessary)
	static bool checkControlAtBoundary(Metric2dContext& mc, MeshContainer2d* mesh);
	/// Count coordinates for new node insertion for longest edge of a triangle
	static bool countNewNodeLongestEdge(Metric2dContext& mc, const DMPoint2d dpoints[], 
		DPoint2d& dnew, MeshTriangle2d* triangle);
	/// Counts quality of mesh elements according to metric
	static bool statMetricQuality(Metric2dContext& mc, MeshContainer2d* mesh, 
		DataStatistics & stats, int criterion = MeshData::QUALITY_SPACE);
	/// Auto-triangulate modeled domain
	static bool autoTriangulate(MeshContainer3d* boundary,
		int sm_count = 2, bool start_clean = true, bool use_cs3d = false);
	/// Auto-triangulate modeled domain
	static bool triangulateWithCS3d(MeshContainer3d* boundary, int sm_count = 2);
	/// Auto-triangulate modeled domain
	static bool autoTriangulateRemesh(MeshContainer3d* boundary,
		int sm_count = 2, bool start_clean = true, bool use_cs3d = false);
	/// Returns statistical information about size of incidention tables of vertices in this mesh
	static bool getIncidenceInfo(MeshContainer2d* mesh, MeshData::StatData & stats);
	/// Tries to improve the quality of triangulation by swapPIng of common edge for neighbouring triangles
	static bool optimizeTriangleBySwapPIng(Metric2dContext& mc, 
		MeshTriangle2d* triangle);
	/// Removes the given point from mesh and re-validates the mesh in the vicinity
	static bool removeTriangulationPoint(Metric2dContext& mc, MeshContainer2d* mesh, 
		MeshPoint2d* point);
	/// Calls smoothing procedures for all faces in the given description of domain several times
	static bool smoothenFaces(MeshContainer3d* boundary, int steps = 1,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0, 
		int method = MeshData::SM_LAPLACE_MIXED | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE);
	/// Calls several smoothing-triangle procedures ("steps" times) for mesh
	static bool smoothen(Metric2dContext& mc, MeshContainer2d* mesh, int steps = 1,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0,
		int method = MeshData::SM_LAPLACE_MIXED | MeshData::SM_TOP_SWAP 
		| MeshData::SM_DEL_SWAP_COMPLETE);
	/// Call topologically oriented procedure of mesh smoothing
	static bool smoothenTopologicalSwap(MeshContainer2d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Smoothen mesh by edge swapPIng (using Delaunay criterion in metric space)
	static bool smoothenDelaunaySwapComplete(Metric2dContext& mc, MeshContainer2d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Smoothen mesh by edge swapPIng (using Delaunay criterion in metric space)
	static bool smoothenDelaunaySwap(Metric2dContext& mc, MeshContainer2d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Tries to move point to new coordinates (gradually, checking for inverted blocks)
	static bool tryMovingPoint(MeshPoint2d *point, const DPoint2d& new_pt);
	/// Tries to improve the mesh locally by moving the given point to the barycenter of incident points
	static bool movePointByLaplace(MeshPoint2d* point);
	/// Tries to improve the mesh locally by moving the given point to the barycenter of incident points
	static bool movePointByLaplaceForVariableMetric(Metric2dContext& mc, MeshPoint2d* point);
	/// Calls the metric based procedure for all points in the given mesh (with metric)
	static bool smoothenMetric(Metric2dContext& mc, MeshContainer2d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Calls the Laplacian-smoothing procedure for all points in the given mesh (with metric)
	static bool smoothenLaplace(Metric2dContext& mc, MeshContainer2d* mesh, bool variable_metric,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Calls the post-smoothing check in the given mesh (with metric)
	static bool smoothenPostCheck(Metric2dContext& mc, MeshContainer2d* mesh);
	/// Calls the Laplacian-smoothing procedure for all points in the given mesh (with metric)
	static bool smoothenLaplaceMixed(Metric2dContext& mc, MeshContainer2d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Performs the first stage of triangulation by adding boundary points into the triangulation
	static int addBoundaryNodes(Metric2dContext& mc, MeshContainer2d* mesh, MeshContainer2d* boundary_mesh);
	/// Performs the second stage of triangulation by adding inner nodes in order to satisfy the control space requirements
	static int addInnerNodes(Metric2dContext& mc, MeshContainer2d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Restores the missing edge (after the unconstrained Delaunay triangulation) by sequential edge swapPIng
	static bool makeEdgeByEdgesSwapPIng(Metric2dContext& mc, MeshPoint2d* pt1, MeshPoint2d* pt2, int dir, 
		char border_type = TagBorder::NONE, int area1 = -1, int area2 = -1);
	/// Automatically creates control space for the given domain-surface based on the curvature of surface and boundary
//	static CS2dPtr createUniformControlSpace(SurfaceConstPtr surface, const MeshContainer2d* boundary);
	/// Restores the Delaunay property in the vicinity of the inserted point by edge swapPIng
	static bool addPointToTriangulationByAngles(Metric2dContext& mc, MeshContainer2d *mesh, MeshPoint2d *point, 
		MeshTriangle2d *containing_triangle = nullptr, bool allow_swap = true,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Performs Delaunay retriangulation using the empty cavity criterion
	static bool addPointToTriangulationByCircle(Metric2dContext& mc, MeshContainer2d *mesh, MeshPoint2d *point, 
		MeshTriangle2d *containing_triangle = nullptr,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Constrains the Delaunay triangulation of boundary nodes (restores all boundary edges, identifies valid triangles and removes the obsolete ones)
	static bool constrainToBorder(Metric2dContext& mc, const MeshContainer2d* boundary, MeshContainer2d* mesh);
	/// Inserts the new point into triangular mesh, finds the containing triangle (if needed) and calls adequate reconfiguration procedure
	static bool addPointToTriangulation(Metric2dContext& mc, MeshContainer2d* mesh, MeshPoint2d* point, 
		MeshTriangle2d* containing_triangle = nullptr,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Calls the triangulation procedures for all domain-surfaces
	static int triangulateFaces(MeshContainer3d* boundary);
	/// Prepares the requested boundary description for Delaunay triangulation procedures
	static int prepareBoundaryMesh(MeshContainer3d* boundary);
	/// Performs the first stage of Delaunay triangulation of boundary nodes (only)
	static MeshContainer2d* createInitialMesh(const MeshContainer2d* boundary, const DRect* bounding_area = nullptr);
	/// Create initial triangulation
	static MeshContainer2d* createInitialMesh(int pct, std::shared_ptr<const SurfaceParametric> surface,
		const DRect& brect, std::shared_ptr<ControlSpace2d> space = nullptr);
public:
	struct ActiveEdge{
		ActiveEdge(MeshEdge2d* _edge = nullptr, double _len = 0.0) : edge(_edge), len(_len) {}
		MeshEdge2d* edge;
		double len;
		//
		short compareTo(ActiveEdge* ae) const {
			if(len == ae->len) return 0;
			else if(len > ae->len) return 1;
			else return -1;
		}
		void preDeleteAll() {}
		int getIndex() const { return index; }
		void setIndex(int i) { index = i; }
		int index;
	};
public:
	/// Maximum number of retriangulation for automatic control space adjustment
	static int param_max_auto_retriangulation_count;
	/// Whether to insert inner nodes during 2D Delaunay triangulation
	static int param_triangulate_with_inner_nodes;
	/// Method used for Delaunay retriangulation
	static int param_triangulation_type;
	/// Method of triangle improvement (where to insert new node)
	static int param_quality_improvement;
	/// Quality threshold for triangle improvement
	static double param_quality_threshold;
	/// Dimension ratio for the initial mesh size
	static double param_initial_mesh_size_ratio;
	/// Criterion used for swapPIng method of quality improvement
	static int param_swap_criterion;
	/// Maximum number of edge swap per single point retriangulation with inner angles criterion
	static int param_swap_maximum;
	/// Whether to show prediction
	static bool show_prediction;
	/// Extra decomposition mode for inner node insertion
	static int param_mesh_decomposition;
	/// Dense layer width for extra decomposition mode for inner node insertion
	static double param_mesh_decomposition_width;
	/// Metric gradation for extra decomposition mode for inner node insertion
	static double param_mesh_decomposition_gradation;
};

#endif // !defined(MESHGENERATOR2D_H__INCLUDED)
