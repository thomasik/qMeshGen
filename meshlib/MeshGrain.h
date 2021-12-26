/////////////////////////////////////////////////////////////////////////////
// MeshGrain.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGRAIN_H__INCLUDED)
#define MESHGRAIN_H__INCLUDED

#include "MeshData.h"
#include "DataVector.h"
#include "DPoint.h"
#include "DRect.h"
#include "DQuadric.h"
#include "SurfaceParametric.h"
#include "Curve2dParametric.h"

class MeshViewSet;
class DMatrix3d;

/**
 * Grain mesh processing
 */
class MeshGrain
{
public:
	struct LineInfo{
		LineInfo(int n0 = -1, int n1 = -1, bool _valid=false) 
			: node0(n0), node1(n1), valid(_valid), trans_line(-1), adjacent_nonplanar(0) {}
		int node0, node1;
		DataVector<DPoint3d> inner_points;
		DataVector<int> incident_faces;
		bool valid;
		int trans_line;
		int adjacent_nonplanar;
	};
	struct FaceInfo{
	public:
		FaceInfo(int b0 = -1, int b1 = -1, bool _valid=false) 
			: block0(b0), block1(b1), valid(_valid), surf_id(-1) {}
	public:
		bool planar() const;
	public:
		DataVector<int> lines;
		DataVector<bool> linear;
		DataVector< std::shared_ptr<const Curve2dParametric> > curves;
		DataVector<int> nodes;
		DataVector<DPoint2d> local_node_params;
		int block0, block1;
		DataVector<DPoint3d> inner_points;
		bool valid;
		int surf_id;
		std::shared_ptr<SurfaceParametric> surface;
		DataVector<DPoint2d> correction_centers;
		DataVector<double> correction_radia;
		DataVector<DVector3d> correction_vectors;
	};
	struct BlockInfo{
		BlockInfo() {}
		DataVector<int> faces;
		DataVector<bool> inverted;
	};
public:
	// check for duplicate lines (due to node-collapsing - but not necessary!)
	static bool checkForDuplicateLines( 
		DataVector<LineInfo> & lines, 
		DataVector<DPoint3d> & nodes);
	// adjust surfaces to better fit face-vertices
	static void adjustSurfacesForNodes(
		DataVector<FaceInfo> & faces, 
		DataVector<DPoint3d> & nodes);
	/// Check vertices for planar faces
	static int recheckPlanarFaceVertices(
		DataVector<BlockInfo> & blocks,
		DataVector<FaceInfo> & faces, 
		DataVector<LineInfo> & lines, 
		DataVector<DPoint3d> & nodes);
	/// Check close nodes after optimization
	static int recheckCloseNodes(
		DataVector<BlockInfo> & blocks,
		DataVector<FaceInfo> & faces, 
		DataVector<LineInfo> & lines, 
		DataVector<DPoint3d> & nodes, 
		DataVector<int> & trans_nodes);
	/// Close line with too close vertices
	static bool closeLine(int i,
		DataVector<BlockInfo> & blocks,
		DataVector<FaceInfo> & faces, 
		DataVector<LineInfo> & lines, 
		DataVector<int> & trans_nodes);
	/// Check topology of model entities
	static MeshViewSet* getFaceViewSet(const FaceInfo& face,
		const DataVector<DPoint3d> & nodes, 
		const DataVector<LineInfo> & lines,
		MeshViewSet* set = nullptr);
	/// Transforms set of text files into .msh format
	static bool parseGrainFiles(const string& dir, const string& msh_file);
	/// Load list of nodes
	static int parseGrainNodeFile(const string& dir, 
		DataVector<DPoint3d> & nodes, 
		DataVector<int> & trans_nodes);
	/// Load list of nodes (short)
	static int parseGrainNodesShortFile(const string& dir, 
		DataVector<DPoint3d> & nodes, 
		DataVector<int> & trans_nodes);
	/// Load list of lines
	static int parseGrainLineFile(const string& dir, 
		DataVector<DPoint3d> & nodes, 
		const DataVector<int> & trans_nodes, 
		DataVector<LineInfo> & lines, 
		DataVector<FaceInfo> & faces);
	/// Load list of lines (short)
	static int parseGrainLinesShortFile(const string& dir, 
		DataVector<DPoint3d> & nodes, 
		const DataVector<int> & trans_nodes, 
		DataVector<LineInfo> & lines, 
		DataVector<FaceInfo> & faces);
	/// Load list of faces
	static int parseGrainFaceFile(const string& dir, 
		const DataVector<DPoint3d> & nodes, 
		const DataVector<int> & trans_nodes, 
		const DataVector<LineInfo> & lines, 
		DataVector<FaceInfo> & faces, 
		DataVector<BlockInfo> & blocks,
		const string& filename = "Planes.txt");
	/// Check topology of model entities
	static bool checkTopology(const DataVector<DPoint3d> & nodes, 
		DataVector<LineInfo> & lines, 
		DataVector<FaceInfo> & faces, 
		DataVector<BlockInfo> & blocks);
	/// Classify blocks (orientation of faces)
	static void classifyGrainBlocks(DataVector<BlockInfo> & blocks,
		const DataVector<FaceInfo> & faces, 
		const DataVector<DPoint3d> & nodes);
	/// Classify lines (straight, arcs, other - approximation)
	static void classifyGrainLines(DataVector<FaceInfo> & faces, 
		DataVector<LineInfo> & lines, 
		const DataVector<DPoint3d> & nodes);
	/// Classify surfaces (planar, quadratic, other - approximation)
	static void classifyGrainSurfaces(DataVector<FaceInfo> & faces, 
		DataVector<LineInfo> & lines, 
		const DataVector<DPoint3d> & nodes);
	/// Create surface meshes from points, for each face separately
	static void createSurfaceMeshes(DataVector<FaceInfo> & faces, 
		DataVector<LineInfo> & lines, 
		const DataVector<DPoint3d> & nodes);
	/// Optimize placement of nodes according to face-surfaces
	static void optimizeNodePlacement(DataVector<FaceInfo> & faces, 
		DataVector<DPoint3d> & nodes);
	/// Solve optimization problem for node placement
	static DPoint3d countOptimizedPoint(const DPoint3d& old_pt, const DMatrix3d& A, const DVector3d& b);
	/// Store parsed grain description into mesh XML file
	static bool storeGrainXml(const string& dir, const DataVector<DPoint3d> & nodes,
		const DataVector<int> & trans_nodes, const DataVector<LineInfo> & lines,
		const DataVector<FaceInfo> & faces, const DataVector<BlockInfo> & blocks);
	/// Store parsed grain description into mesh XML file (with planes)
	//static bool storeGrainXmlWithPlanes(const string& dir, const DataVector<DPoint3d> & nodes,
	//	const DataVector<int> & trans_nodes, const DataVector<LineInfo> & lines,
	//	const DataVector<FaceInfo> & faces, const DataVector<BlockInfo> & blocks);
public:
	// tolerance for geometric approximation (plane, quadric, etc.)
	static double param_grain_tolerance;
	// tolerance for how close can nodes be to be treated as different vertices
	static double param_grain_node_identity_tolerance;
	// what to show
	static int param_grain_view;
	// preset bounding box
	static DBox param_bounding_box;
};

#endif // !defined(MESHGRAIN_H__INCLUDED)
