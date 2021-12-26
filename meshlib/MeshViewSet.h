/////////////////////////////////////////////////////////////////////////////
// MeshViewSet.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHVIEWSET_H__INCLUDED)
#define MESHVIEWSET_H__INCLUDED

#include <memory>

#include "DataVector.h"
#include "DPoint.h"
#include "DRect.h"
#include "MeshData.h"
#include "DataMatrix.h"
#include "DataHashTable.h"

class MeshBlock;
class MeshPoint2d;
class MeshPoint3d;
class MeshEdge2d;
class MeshEdge3d;
class MeshElement;
class MeshFace;
class SurfaceParametric;
class MeshContainer2d;
class MeshContainer3d;
class MeshViewSet;

enum ShowMode {SHOW_THESAME, SHOW_DOMAIN, SHOW_MESH_2D, SHOW_MESH_2DQ, SHOW_MESH_SURF, SHOW_MESH_SURF_QUALITY,
				SHOW_MESH_3D, SHOW_MESH_3DQ, SHOW_MESH_CONTROL, SHOW_MESH_SURF_METRIC_GRADATION};

class MeshViewExt
{
public:
	virtual ~MeshViewExt() {}
public:
	virtual void showViewSet(const string& desc, MeshViewSet* set, int t) = 0;
	virtual void showViewSetNoReset(const string& desc, MeshViewSet* set, int t){
		showViewSet(desc, set,t); // by default, but can be changed if possible
	}
};

struct MeshViewPointData{
	MeshViewPointData(const FPoint3d& _pt, char _border = 0, int _id = 0, bool _numbered = false)
		: pt(_pt), part(0), border(_border), id(_id), numbered(_numbered), hidden(0) {}
	FPoint3d pt;
	int part;
	char border;
	int id;
	bool numbered;
	int hidden; // 0 - visible, 1 - clipped, 2 - not visible
};

struct MeshViewBlockData;

struct MeshViewEdgeData{
	MeshViewEdgeData(const FPoint3d& _pt1, const FPoint3d& _pt2, char _border = 0) 
		: pt1(_pt1), pt2(_pt2), part(0), border(_border), hidden(0) {}
	FPoint3d pt1, pt2;
	DataVector<std::shared_ptr<MeshViewBlockData>> adjacent_blocks;
	int part;
	char border;
	int hidden; // 0 - visible, 1 - clipped, 2 - not visible
};
struct MeshViewFaceData{
	MeshViewFaceData(int ct) : pts(ct), wpts(ct), indices(3*(ct-2)), part(0), hidden(0) { assert( ct > 2 ); }
	void countNormal() { normal = (pts[1]-pts[0]).crossProduct(pts[2]-pts[0]).normalized(); }
	DataVector<FPoint3d> pts;
	DataVector<float> wpts;
	DataVector<int> indices;
	FVector3d normal;
	int part;
	int area_id;
	float quality;
	int hidden; // 0 - visible, 1 - clipped, 2 - not visible
};

struct MeshViewBlockData{
	MeshViewBlockData(unsigned int pct, unsigned int fct) 
		: part(0), pts(pct), normals(fct), adjacent(fct), indices(fct*3), hidden(0) {}
	int part;
	DataVector<FPoint3d> pts;
	DataVector<FVector3d> normals;
	DataVector<std::shared_ptr<MeshViewBlockData>> adjacent;
	DataVector<unsigned char> indices;
	int area_id;
	double quality;
	int hidden; // 0 - visible, 1 - clipped, 2 - not visible
};

class MeshViewSet
{
public:
	MeshViewSet(size_t pct=100, size_t ect=100, size_t fct=100, size_t bct=100) :
		m_points(pct), m_edges(ect), m_faces(fct), m_blocks(bct),
		m_hash_blocks((unsigned int)std::max(bct, (size_t)(1024*32)), nullptr), 
		m_hash_faces((unsigned int)std::max(fct, (size_t)(1024 * 32)), nullptr),
		m_hash_edges((unsigned int)std::max(ect, (size_t)(1024 * 32)), nullptr),
		m_hash_points((unsigned int)std::max(pct, (size_t)(1024 * 32)), nullptr),
		m_mesh2d(0), m_mesh3d(0), polygon_fill(FILL_AREA) { }
private:
	MeshViewSet(const MeshViewSet& ) = delete;
public:
	struct LabelInfo {
		LabelInfo(const FPoint3d& _pt = FPoint3d::zero, const string& _label = "x") : pt(_pt), label(_label) {}
		FPoint3d pt;
		string label;
	};
	struct ClipPlane {
		ClipPlane(const FVector3d& _vn = FVector3d::zero, float _d = 0.0f) : vn(_vn), d(_d) {}
		bool clipped(const FPoint3d& pt) const {
			return pt.x*vn.x + pt.y*vn.y + pt.z*vn.z + d > 0.0;
		}
		FVector3d vn;
		float d;
	};
public:
	void prepareFreePlace(size_t pct=0, size_t ect=0, size_t fct=0, size_t bct=0){
		if(pct) m_points.prepare(pct);
		if(ect) m_edges.prepare(ect);
		if(fct) m_faces.prepare(fct);
		if(bct) m_blocks.prepare(bct);
	}
	// 2d
	void addPoint(const MeshPoint2d* point, std::shared_ptr<const SurfaceParametric> surface, int part = 0);
	void addEdge(const MeshEdge2d* edge, std::shared_ptr<const SurfaceParametric> surface);
	void addElement(const MeshElement* element, std::shared_ptr<const SurfaceParametric> surface,
		bool proper_orientation = true, int id = -2);
	void addElementWithEdges(const MeshElement* element, std::shared_ptr<const SurfaceParametric> surface,
		bool proper_orientation = true, int id = -2);
	// 3d
	void addPoint(const MeshPoint3d* point, int part = 0);
	void addPoint(const FPoint3d& point, int part = 0);
	void addPoint(const FPoint3d& point, int part, int id);
	void addEdge(const MeshEdge3d* edge, int id = -2);
	void addEdge(const FPoint3d& p0, const FPoint3d& p1, int id = -2);
	void addFace(const FPoint3d& p0, const FPoint3d& p1, const FPoint3d& p2, int id = -2);
	void addFace(const MeshFace* face, int id = -2, double shrink = 1.0, bool proper_orientation = true);
	void addFaceWithEdges(const MeshFace* face, int id = -2, double shrink = 1.0, bool proper_orientation = true);
	void addFaceWithEdgesAndPoints(const MeshFace* face, int id = -2, double shrink = 1.0, bool proper_orientation = true);
	void addFace(const MeshFace* face, const FVector3d& shift, int id = -2);
	void addBlock(const MeshBlock* block, int id = -2);
	void addBlockWithEdges(const MeshBlock* block, int id = -2);
	void addTetra(const FPoint3d& point, float dx, float dy, float dz, double q = 1.0);
	void addEmptyBlockWithEdges(const MeshBlock* block, int id = -2);
	void addEdges(const MeshBlock* block);
	void addEdges(const MeshFace* face);
	void addPolygon(const DataVector<DPoint3d> & poly, int id = -2);
	void addPolygonConvex(const DataVector<DPoint3d> & poly, int id = -2);
	void addInfo( const string & sname, int ivalue) { m_info.add(sname); m_info.add(to_string(ivalue)); }
	void addInfo( const string & sname, size_t stvalue) { m_info.add(sname); m_info.add(to_string(stvalue)); }
	void addInfo( const string & sname, double dvalue) { m_info.add(sname); m_info.add(to_string(dvalue)); }
	void addInfo( const string & sname, const string & svalue) { m_info.add(sname); m_info.add(svalue); }
	DBox getBoundingBox() const;
	void clearLabels();
	void addLabel(const FPoint3d& pt, const string& label);
	void log();
	void setMesh(const MeshContainer2d* mesh) { m_mesh2d = mesh; m_mesh3d = 0; }
	void setMesh(const MeshContainer3d* mesh) { m_mesh2d = 0; m_mesh3d = mesh; }
	bool storeMatlabFile(const string& fname, bool use_mesh, const DPoint3d& clip, double quality_clip) const;
	bool storeMatlabFile(const string& fpath, const string& fname = "mesh_view.m") const;
	static FPoint3d transProject(float trans_matrix[], const FPoint3d& pt);
	bool storeEPSFile(float trans_matrix[], float model_matrix[], int view_mode = 15, const string& fname = "mesh_view.eps") const;
	void setHiddenByClipPlane(const MeshViewSet::ClipPlane& cp);
public:
	/// polygon fill modes
	enum PolygonFill { FILL_AREA, FILL_QUALITY, FILL_NODES, FILL_LGRAY };
	void setPolygonFillMode(PolygonFill mode) { polygon_fill = mode; }
	PolygonFill getPolygonFillMode() const { return polygon_fill; }
public:
	static void showViewSet(const string& desc, MeshViewSet* set, int t = -1);
	static void showViewSetNoReset(const string& desc, MeshViewSet* set, int t = -1);
	static void showDebugMesh(const string& desc, const MeshContainer2d* mesh, 
		const MeshElement* el1, const MeshElement* el2 = nullptr, double radius = 2.0, int t = -1);
	static void showDebugMesh(const string& desc, const MeshContainer2d* mesh, 
		const MeshPoint2d* pt1, const MeshPoint2d* pt2 = nullptr, double radius = 2.0, int t = -1);
	static void showDebugMesh(const string& desc, const MeshContainer3d* mesh, 
		const MeshBlock* el1, const MeshBlock* el2 = nullptr, double radius = 2.0, int t = -1);
	static void showDebugMesh(const string& desc, const MeshContainer3d* mesh, 
		const MeshPoint3d* pt1, const MeshPoint3d* pt2 = nullptr, double radius = 2.0, int t = -1);
	static void showDebugMesh(const MeshContainer3d* boundary, const MeshContainer3d* mesh, 
		int phase, int t = -1);
public:
	DataVector<std::shared_ptr<MeshViewPointData>> m_points;
	DataVector<std::shared_ptr<MeshViewEdgeData>> m_edges;
	DataVector<std::shared_ptr<MeshViewFaceData>> m_faces;
	DataVector<std::shared_ptr<MeshViewBlockData>> m_blocks;
	DataVector<LabelInfo> m_labels;
	DataVector<string> m_info;
private:
	/// hash table for duplicate checking
	DataHashTableKeyValue<const MeshBlock*, std::shared_ptr<MeshViewBlockData>> m_hash_blocks;
	DataHashTableKeyValue<const void*, std::shared_ptr<MeshViewFaceData>> m_hash_faces;
	DataHashTableKeyValue<const void*, std::shared_ptr<MeshViewEdgeData>> m_hash_edges;
	DataHashTableKeyValue<const void*, std::shared_ptr<MeshViewPointData>> m_hash_points;
	/// surface mesh reference
	const MeshContainer2d* m_mesh2d;
	/// volume mesh reference
	const MeshContainer3d* m_mesh3d;
	/// polygon fill mode
	PolygonFill polygon_fill;
public:
	/// face/block shrink for visualization
	static float param_shrink;
	/// whether any visualization should be prepared/shown
	static int param_show_visualization;
	/// reference to the visualization module
	static MeshViewExt *m_view;
};

#endif // !defined(MESHVIEWSET_H__INCLUDED)
