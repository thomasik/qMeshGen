/////////////////////////////////////////////////////////////////////////////
// MeshDecompositionUTC.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	21-24 September 2005
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHDECOMPOSITION_H__INCLUDED)
#define MESHDECOMPOSITION_H__INCLUDED

class MeshContainer3d;
class MeshContainer2d;
class MeshViewSet;
class MeshPoint2d;
class MeshPoint3d;
class MeshElement;
class MeshFace;
class MeshEdge3d;
class SurfaceParametric;

#include "DPoint.h"
#include "DataVector.h"

typedef DataVector<MeshFace*> DecompFaceList;
typedef DataVector<MeshPoint3d*> DecompNodeList;
typedef DataVector<MeshEdge3d**> DecompEdgeArrayList;

/**
 * Surface mesh decompotision by cutting for 3D mesh generation
 */
class MeshDecompositionUTC
{
public:
	struct CutNode{
		CutNode() : color(-1), visited(false) {}
		MeshPoint3d* mpoint;
		DataVector<int> edges_out;
		DataVector<int> edges_in;
		int	color;
		double dist;
		double penalty;
		double min_flow_penalty;
		bool visited;
	};
	struct CutEdge{
		CutEdge() : tag(0) {}
		int node0;
		int node1;
		MeshEdge3d* medge;
		double length;
		double orientation;
		int crossing;
		int tag;
		double penalty;
	};
public:
	MeshDecompositionUTC(const string& fname);
	~MeshDecompositionUTC();
public:
	void storeOpenCutMesh(int part_id, int part_pct, MeshPoint3d** part_points, 
		int part_fct, MeshFace** part_faces, bool with_cut_faces);
	void storeCutFaces(DecompFaceList* cut_faces);
	void storeCutContours(DecompEdgeArrayList& contours);
	void storeClosedCutMesh(int part_id, 
		int part_pct, MeshPoint3d** part_points, int part_fct, MeshFace** part_faces, 
		int cut_pct, MeshPoint2d** cut_points, int cut_fct, MeshElement** cut_faces,
		SurfaceParametric* surface, bool orientation_valid);
	void storeInterfaceCutMesh(int part_id, int cut_pct, MeshPoint2d** cut_points, 
		int cut_fct, MeshElement** cut_faces, SurfaceParametric* surface);
	void storeInterfaceCutContours(int part_id, int cut_pct, MeshPoint2d** cut_points);
	void markPointsOrientation();
	DecompFaceList* selectCutTriangles();
	DecompNodeList* selectPoints(DecompFaceList* cut_faces);
	DecompEdgeArrayList findContour3D(DecompNodeList* node_list, DecompFaceList* cut_faces);
	MeshContainer2d* createCutMesh(DecompEdgeArrayList& contours);
//	bool storeCutSurfaceMeshes(MeshContainer2d* cut_mesh);
	static MeshContainer3d* readIncidenceMeshDescription(const string& fname);
	bool run();
	double countNodePenalty(CutNode& cn);
	double countEdgePenalty(CutEdge& ce);
	/// Returns the "screenshot" of this mesh for visualization
	MeshViewSet* getViewSet(MeshViewSet* set = nullptr);
private:
	MeshContainer3d* m_mesh;
	DPoint3d m_plane_center;
	DVector3d m_plane_normal;
	ofstream m_log_file;
	string m_cut_mesh_name;
	double m_epsilon;
	double m_steep_edge_threshold;
	int m_visualization;
	double m_penalty_node_dist;
	double m_penalty_node_class;
	double m_penalty_edge_length;
	double m_penalty_edge_orient;
	double m_penalty_edge_bad_orient;
	double m_penalty_edge_cross;
	bool m_orientation_switch;
};

#endif // !defined(MESHDECOMPOSITIONUTC_H__INCLUDED)
