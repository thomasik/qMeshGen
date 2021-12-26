/////////////////////////////////////////////////////////////////////////////
// MeshGenerator3dSurface.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2010-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR3DSURFACE_H__INCLUDED)
#define MESHGENERATOR3DSURFACE_H__INCLUDED

#include "TagExtended.h"
#include "MeshData.h"
#include "DataList.h"
#include "DataMatrix.h"
#include "MeshEdge3d.h"
#include "DataVector.h"
#include "DataPtrMatrix.h"
#include "DataContainer.h"

class MeshContainer3d;
class MeshContainer3dSurface;
class Metric3dContext;
class MeshPoint3d;
class DPoint3d;
class MeshFace;
class ControlSpace3d;
class ControlSpace3dAdaptive;
class MeshDomainVolume;
class SurfaceParametric;

/**
 * This class gathers procedures responsible for processing 3d (triangular) surface meshes
 */
class MeshGenerator3dSurface
{
public:
	static void remeshSurfaceMesh( MeshContainer3dSurface* mesh, double min_dist = -1.0 );
	/// get local surface valid for the given set of (local) mesh points
	//static Curve3dConstPtr getLocalCurveForPoints( Metric3dContext& mc, 
	//	MeshContainer3dSurface* mesh, const DataVector<MeshEdge3d*> & medges,
	//	const DataVector< MeshPoint3d* > &mpoints );
	/// get local surface valid for the given set of (local) mesh points
	static std::shared_ptr<const SurfaceParametric> getLocalSurfaceForPoints( Metric3dContext& mc,
		MeshContainer3dSurface* mesh, const DataVector<MeshFace*> & mfaces, 
		const DataVector< MeshPoint3d* > &mpoints );
	/// improve local surface approximation quality (returns number of reassigned surfaces)
	static int fixLocalSurfaceSeams(Metric3dContext& mc, MeshContainer3dSurface* mesh );
	/// improve local surface approximation quality (returns number of reassigned surfaces)
	static bool identifyLocalCurves(Metric3dContext& mc, MeshContainer3dSurface* mesh, double tolerance );
	/// improve local surface approximation quality (returns number of reassigned surfaces)
	//static int optimizeLocalSurfaces(Metric3dContext& mc, MeshContainer3dSurface* mesh );
	/// Check local surface information for the given point
	static bool checkPointForLocalSurfaces( Metric3dContext& mc, MeshContainer3dSurface* mesh, MeshPoint3d* point );
	/// Create additional local surfaces for points with low locals-surface-inside-ratio
	//static int complementLocalSurfaces(Metric3dContext& mc, MeshContainer3dSurface* mesh );
	/// try to split domains with open borders
	static int fixOpenBorders(MeshContainer3dSurface* mesh, 
		DataVector< std::shared_ptr<DataVector<MeshFace*>> > & sub_domains, int sub_tag_init);
	/// Search for local reparameterization surfaces (for whole mesh, or only faces with the given tag)
	static void identifyLocalSurfaces(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Split surface mesh into closed sub-domains (using boundary and tag info
	static int splitToSubdomains(MeshContainer3dSurface* mesh, 
		DataVector< std::shared_ptr<DataVector<MeshFace*>> > & sub_domains, 
		int sub_tag_init, TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Identify local surfaces using greedy approach
	static int identifyLocalSurfacesGreedy(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		std::shared_ptr<const DataVector<MeshFace*>> sub_faces );
	/// Identify local surfaces using higher-curvature first
	static int identifyLocalSurfacesHighCurvatureFirst(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		std::shared_ptr<const DataVector<MeshFace*>> sub_faces );
	/// Identify local surfaces using normal-clustering, then quadtree of quadrics
	static int identifyLocalSurfacesViaNormals(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		std::shared_ptr<DataVector<MeshFace*>> sub_faces );
	/// Check and fix faces using topological orientation info
	static bool checkAndFixInvertedFacesTopological(MeshContainer3dSurface* mesh);
	/// Calculate base normals for surface faces
	static bool calculateNormals(MeshContainer3dSurface* mesh);
	/// Calculate base normals for surface faces in the given subdomain
	static bool updateNormals(MeshContainer3dSurface* mesh, 
		std::shared_ptr<const DataVector<MeshFace*>> sub_faces);
	/// Try improving quality of the given face
	static void improveFace(Metric3dContext & mc, MeshContainer3dSurface* mesh, MeshFace* face);
	/// Rearrange boundary nodes for optimum length of border edges
	static void optimizeBoundaryEdges(Metric3dContext& mc, MeshContainer3dSurface* mesh);
	/// Transforms 2 face-triangles by switching a common edge (on boundary), returns swapped edge if successful
	static MeshEdge3d* swap22(Metric3dContext& mc, MeshContainer3dSurface* mesh, MeshEdge3d* edge);
	/// Split the given edge (normal or boundary) by inserting a new point in the middle
	static MeshPoint3d* splitEdgeWLT(Metric3dContext& mc, MeshContainer3dSurface* mesh, MeshEdge3d* edge,
		DataContainer<MeshEdge3d::ActiveEdge> * act_edges = nullptr, TagExtended::TagType act_tag = TagExtended::TAG_NONE);
	/// Split the given edge (normal or boundary) by inserting a new point in the middle
	static MeshPoint3d* splitEdgeSimple(MeshContainer3dSurface* mesh, MeshEdge3d* edge);
	/// Insert nodes at boundary, for too long edges, return number of inserted points
	static int splitEdges(Metric3dContext& mc, MeshContainer3dSurface* mesh);
	/// Create ACS using approximation of curvature from surface mesh faces
	static std::shared_ptr<ControlSpace3d> createACSFromApproxCurvatureDirect(MeshContainer3dSurface* surface_mesh);
	/// Create ACS using approximation of curvature from surface mesh faces
	static std::shared_ptr<ControlSpace3d> createACSFromApproxCurvatureAverage(MeshContainer3dSurface* surface_mesh);
	/// Check and fix all surface meshes for 
	static int validateSurfaceMesh(MeshContainer3dSurface* surface_mesh, double min_dist = 0.0);
	/// Create surface mesh from point/face data, with small fixes
	static MeshContainer3dSurface* createMeshFromFaces(DataVector<DPoint3d> & points, DataSimpleList< DataVector<int> > & faces);
	/// Update sizing data in the acs with curvature information aproximated from local parameterization surfaces (if available)
	static bool updateACSwithLocalCurvature(std::shared_ptr<ControlSpace3d> acs, MeshContainer3dSurface* mesh);
	/// Create new adaptive control space with sizing taken from mesh faces
	static std::shared_ptr<ControlSpace3dAdaptive> createACSfromMeshFaces(MeshContainer3dSurface* mesh);
	/// Creates surface mesh from volume mesh (copy boundary faces and entities)
	static MeshContainer3dSurface* copySurfaceMeshFromVolumeMesh(MeshContainer3d* mesh3d, MeshDomainVolume* mdv = nullptr);
	/// Check and fix faces with normals inverted with respect to base normal
	static int checkAndFixInvertedFacesForBaseNormal(MeshContainer3dSurface* mesh);
	/// Check and fix faces with normals inverted with respect to local surfaces
	static int checkAndFixInvertedFacesForLocalSurfaces(Metric3dContext& mc, MeshContainer3dSurface* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Smoothen mesh by edge swapPIng (using quality criterion in metric space)
	static bool smoothenQualitySwapComplete(Metric3dContext& mc, MeshContainer3dSurface* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Smoothen mesh by edge swapPIng (near border for badly aligned faces)
	static bool improveNearBorder(Metric3dContext& mc, MeshContainer3dSurface* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Smoothen mesh by edge swapPIng (using edge lengths in metric space)
	static bool smoothenEdgeLengthSwap(Metric3dContext& mc, MeshContainer3dSurface* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Smoothen mesh by edge swapPIng (using mean ratio quality in metric space)
	static bool smoothenQualitySwap(Metric3dContext& mc, MeshContainer3dSurface* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Calls the Laplacian-smoothing procedure for all points in the given mesh
	static bool smoothenLaplace(Metric3dContext& mc, MeshContainer3dSurface* mesh, bool variable_metric,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Calls the Laplacian-smoothing procedure for all points in the given mesh
	static bool smoothenLaplaceBoundary(Metric3dContext& mc, MeshContainer3dSurface* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0, int forbid_tag_value = 0);
	/// Calls the Laplacian-smoothing procedure for all points in the given mesh (with metric or without, depending on gradation)
	static bool smoothenLaplaceMixed(Metric3dContext& mc, MeshContainer3dSurface* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Call topologically oriented procedure of mesh smoothing
	static bool smoothenTopologicalSwap(Metric3dContext& mc, MeshContainer3dSurface* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Calls several smoothing-triangle procedures ("steps" times) for mesh
	static bool smoothen(Metric3dContext& mc, MeshContainer3dSurface* mesh, int steps = 1,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0,
		int method = MeshData::SM_LAPLACE_MIXED | MeshData::SM_TOP_SWAP 
		| MeshData::SM_DEL_SWAP_COMPLETE);
	/// Test Laplace methods
	static bool testLaplace(Metric3dContext& mc,  MeshContainer3dSurface* mesh, MeshPoint3d* point);
	/// Tries to improve the mesh locally by moving the given point to the mass-center of adjacent faces
	static bool movePointByMassCenter(MeshContainer3dSurface* mesh, Metric3dContext& mc, MeshPoint3d* point);
	/// Tries to improve the mesh locally by moving the given point to the barycenter of incident points
	static bool movePointByLaplace(MeshContainer3dSurface* mesh, Metric3dContext& mc, 
		MeshPoint3d* point, std::shared_ptr<const SurfaceParametric> surface = nullptr);
	/// Tries to improve the mesh locally by moving the given point to the barycenter of incident points
	static bool movePointByLaplaceForVariableMetric(MeshContainer3dSurface* mesh, Metric3dContext& mc, MeshPoint3d* point);
	/// Tries to improve the mesh locally by moving the given (boundary) point to the barycenter of incident points
	static bool moveBoundaryPointByLaplace( MeshContainer3dSurface* mesh, Metric3dContext& mc, MeshPoint3d* point );
	/// Tries to improve the mesh locally by moving the given (boundary) point to the barycenter of incident points
	static bool moveBoundaryPointByLaplaceForVariableMetric(MeshContainer3dSurface* mesh, Metric3dContext& mc, MeshPoint3d* point);
	/// Modify surface mesh according to the given CS3d, with local operators (edge collapsing and smoothing)
	static bool remeshSurfaceMeshWithLocalTransformations(Metric3dContext& mc, MeshContainer3dSurface* mesh);
	/// Modify crack surface mesh according to the given CS3d, with local operators (edge collapsing and smoothing)
	static bool remeshCrackSurfaceMeshWithLocalTransformations(Metric3dContext& mc, MeshContainer3dSurface* mesh);
	/// Collapse (all) vertices rank two and (maybe) three 
	static int collapseAllInnerVerticesRankLow(Metric3dContext& mc, MeshContainer3dSurface* mesh);
	/// Collapse too short edge, special case for poly-mesh and vertex rank two
	static bool collapseInnerVertexRankTwo(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		MeshPoint3d* removed_point, bool no_border_points = true, bool refill_edges = true, 
		DataContainer<MeshEdge3d::ActiveEdge> * act_edges = nullptr, TagExtended::TagType act_tag = TagExtended::TAG_NONE);
	/// Fix inverted triangle faces, return number of fixed faces
	static int fixInvertedTriangles(Metric3dContext& mc, MeshContainer3dSurface* mesh);
	/// Collapse too thin triangle faces, return number of fixed faces
	static int fixThinTriangles(Metric3dContext& mc, MeshContainer3dSurface* mesh);
	/// Collapse too short edges, return number of removed edges
	static int collapseInnerEdges(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		bool no_border_points = true, bool refill_edges = true,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Collapse too short edges, return number of removed edges
	static int collapseBoundaryEdges(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		bool no_corner_points = true, bool refill_edges = true,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value1 = 1, int tag_value2 = 1);
	/// Collapse too short edge
	static bool collapseBoundaryEdge(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		MeshEdge3d* col_edge, bool check_length = true, 
		bool no_corner_points = true, bool refill_edges = true, 
		DataContainer<MeshEdge3d::ActiveEdge> * act_edges = nullptr, TagExtended::TagType act_tag = TagExtended::TAG_NONE,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value1 = 1, int tag_value2 = 1);
	/// Removes the given point from mesh (for low-rank points)
	static bool removeTriangulationPointSimple(MeshContainer3dSurface* mesh, MeshPoint3d* point, 
		DataContainer<MeshEdge3d::ActiveEdge> * act_edges = nullptr, TagExtended::TagType act_tag = TagExtended::TAG_NONE);
private:
	static const double LOCAL_SURFACE_TOL_F[]; // tolerance coefficient
	static const int    LOCAL_SURFACE_LAY_F[]; // number of topological layers
	static const double LOCAL_SURFACE_MTR_F[]; // metric radius for geometrical layer
	static const int	LOCAL_SURFACE_TOL_F_CT;	// array lengths
	static const int	LOCAL_SURFACE_TOL_F_CT_WHOLE;	// for whole-surface approximation
public:
	static const double MIN_BEDGE_LEN;
	static const double MAX_BEDGE_LEN_JOINED;
	static const double MIN_EDGE_LEN;
	static const double MAX_BEDGE_LEN;
	static const double MAX_EDGE_LEN;
public:
	enum LocalSurfaceIdentification { LSI_GREEDY = 0, LSI_HIGH_CURVATURE_FIRST = 1, LSI_VIA_NORMALS = 2 };
public:
	/// Parameter for sharp edges identification (maximum scalar product value for normals)
	static double param_sharp_edge_threshold; 
	/// Parameter for approximation of local surfaces and curves
	static double param_local_shape_tolerance;
	/// Method of local surface identification
	static int param_local_surface_identification_method;
};

#endif // !defined(MESHGENERATOR3DSURFACE_H__INCLUDED)
