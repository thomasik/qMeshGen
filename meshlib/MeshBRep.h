/////////////////////////////////////////////////////////////////////////////
// MeshBRep.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHBREP_H__INCLUDED)
#define MESHBREP_H__INCLUDED

#include <memory>

#include "DPoint.h"
#include "DMetric2d.h"
#include "DMetric3d.h"
#include "DEquation.h"
#include "DataVector.h"
#include "DataHashTable.h"
#include "Curve2dParametric.h"
#include "Curve3dParametric.h"
#include "SurfaceParametric.h"
#include "DataList.h"

class MeshContainer2d;
class MeshContainer3d;
class MeshPoint2d;
class MeshPoint3d;
class MeshBoundaryCondition;
class MeshingException;
class ControlSpace2d;
class ControlSpace3d;
class MeshDomainSurface;
class MeshDomainVolume;
class MeshFace;

/**
 * Class responsible for (parsing), validating and storing the description of the
 * domain, mesh, etc in BRep (boundary representation) format.
 */
class MeshBRep
{
public:
	MeshBRep() : hash_curves(nullptr), hash_surfaces(nullptr), hash_points(nullptr), 
		hash_faces(nullptr), hash_vertices2d(nullptr), m_ctable(nullptr) {}
	// virtual destructor
	virtual ~MeshBRep();
public:
	/// Parse model description data from file
	virtual bool parseFile(const string& fname) = 0;
	/// Check validity of data
	virtual bool validate();
	/// Create boundary mesh from parsed description
	virtual MeshContainer3d* createDomainMesh();
	/// Clear all data
	virtual void clear();
protected:
	/// Data of input points
	struct PointItem{
		PointItem(int _id) : id(_id), sid(-1), cid(-1), bid(-1) {}
		int id;
		int sid;
		int cid;
		int bid;
		DPoint3d pt;
		string label;
	};
	/// Data of input edges
	struct EdgeItem{
		EdgeItem(int _sid, int _vid0, int _vid1)
			: sid(_sid), vid0(_vid0), vid1(_vid1), fixed(true), cid(-1), t0(0.0), t1(1.0) {}
		int sid;
		int vid0, vid1;
		bool fixed;
		int cid;
		double t0, t1;
		string tag;
	};
	/// Face-point
	struct FacePoint{
		FacePoint(int _PId = 0) : pid(_PId), cid(-1), params_available(false) {}
		int pid;
		int cid;
		DPoint2d params; // params for face point
		bool params_available;
	};
	/// Data of face description
	struct FaceAreaItem{
		FaceAreaItem(int _sub_id = 1, bool _closed = true, bool _filled = true) 
			: sub_id(_sub_id), closed(_closed), filled(_filled) { }
		int sub_id;
		bool closed;
		bool filled;
		DataVector<FacePoint> fpoints;
	};
	/// Data of face description
	struct FaceItem{
		FaceItem(int _id, int _sid = 0) 
			: id(abs(_id)), sid(_sid) { }
		int id;
		int sid;
		DataVector<FaceAreaItem> areas;
		string tag;
	};
	/// Data of block description
	struct BlockItem{
		BlockItem(int _id) : id(_id) {}
		int id;
		DataVector<int> faces;
		DataVector<bool> oriented;
		//DVector3d prism_vector;
		DataCompoundList<DPoint3d> mesh_vertices;
		DataCompoundList<int> mesh_elements;

	};
	/// Data of control description
	struct Control2dItem{
		enum ControlType { CS_UNKNOWN = 0, CS_POINT = 1, CS_SEGMENT = 2, CS_CURVE = 3, CS_DOMAIN = 4};
		Control2dItem(int _fid, int _control_type) : fid(_fid), control_type(_control_type), multi_nr(0) {}
//		bool parse(istream& is, int id);
		bool parseCDS(const char* lx_str, const char* ly_str, const char* angle_str,
			DEquationConstTable *ctable);
		bool parseDomain(const char* domain1_str, const char* domain2_str,
			DEquationConstTable *ctable);
		void adjust();
		bool full_analytic() const;
		DPoint2d pt[2];
		int fid;
		int control_type;
		int multi_nr;
		int curve_id;
		bool directional;
		ControlDataMatrix2d data;
		double radius;
		DEquation equations[5];
		string str_equations[5];
	};
	/// Data of control description
	struct Control3dItem{
		enum ControlType { CS_UNKNOWN = 0, CS_POINT = 1, CS_SEGMENT = 2, 
			CS_CURVE = 3, CS_DOMAIN = 4, CS_TRIANGLE = 6, CS_INTERNAL = 7 };
		Control3dItem(int _bid, int _control_type) : bid(_bid), control_type(_control_type), multi_nr(0) {}
//		bool parse(istream& is, int id);
		bool parseCDS(const char* lx_str, const char* ly_str, const char* lz_str, 
			const char* ax_str, const char* ay_str, const char* az_str,
			DEquationConstTable *ctable);
		bool parseDomain(const char* domain1_str, const char* domain2_str, 
			DEquationConstTable *ctable);
		void adjust();
		bool full_analytic() const;		
		DataVector<DPoint3d> pts;
		int bid;
		int control_type;
		int multi_nr;
		int curve_id;
		bool directional;
		ControlDataMatrix3d data;
		double radius;
		DEquation equations[8];	// 0-2: lengths, 3-5: angles, 6-7: domains
		string str_equations[8];
		DBox bbox;
	};
protected:
	/// Create curves from description
	bool createCurves();
	/// Create surfaces from description
	bool createSurfaces();
	/// Create vertices 3d from description
	MeshContainer3d* createPoints3d();
	/// Create faces from description
	bool createFaces();
	/// Create the i-th face
	MeshDomainSurface* createFace(FaceItem* item, SurfaceConstPtr surface);
	/// Prepare vertices for face
	bool prepareFaceVertices(FaceItem* item, SurfaceConstPtr surface);
	/// Prepare (special, i.e. curvilinear) edges for face
	bool prepareFaceEdges(FaceItem* item);
	/// Create 2d area-boundary-mesh description for face in parametric space
	MeshContainer2d* createFaceBoundaryMesh2d(FaceItem* item);
	/// Prepare (user) control space for face
	std::shared_ptr<ControlSpace2d> prepareUserCS(FaceItem* item, MeshContainer2d* boundary);
	/// Create domain-surface for face 3d representation
	MeshDomainSurface* createFace3D(MeshContainer2d* boundary);
	/// Create blocks from description
	bool createBlocks(MeshContainer3d* mesh);
	/// Create mesh blocks from description
	MeshContainer3d* createMeshBlocks();
	/// Create surface-mesh from description
	MeshContainer3d* createSurfaceMesh();
	/// Prepare (user) control space for block
	std::shared_ptr<ControlSpace3d> prepareUserCS(BlockItem* item, MeshDomainVolume* volume);
	/// Clear unused surfaces, curves and points
	bool clearUnusedData(MeshContainer3d* mesh);
	/// Set border-tags for model faces, edges and vertices
	void markBoundaryTags(MeshContainer3d* mesh);
protected:
	bool points_2d_only;
	// parse lists
	DataVector<std::shared_ptr<PointItem>> point_list;
	DataVector<std::shared_ptr<PointItem>> freepoint_list;
	DataVector<std::shared_ptr<EdgeItem>> edge_list;
	DataVector< std::shared_ptr<Curve2dParametric> > curve_list;
	DataVector< std::shared_ptr<SurfaceParametric> > surface_list;
	DataVector<std::shared_ptr<FaceItem>> face_list;
	DataVector<std::shared_ptr<BlockItem>> block_list;
	DataVector<std::shared_ptr<Control2dItem>> control2d_list;
	DataVector<std::shared_ptr<Control3dItem>> control3d_list;
	DataVector<MeshPoint3d*> smesh_point_list;
	DataVector<MeshFace*> smesh_face_list;
	DataVector< std::shared_ptr<SurfaceParametric> > smesh_surface_list;
	DataVector<std::shared_ptr<Curve3dParametric> > smesh_curve_list;
	// hash tables of created objects
	DataHashTableKeyValue<int, std::shared_ptr<Curve2dParametric> >	*hash_curves;
	DataHashTableKeyValue<int, std::shared_ptr<SurfaceParametric> >	*hash_surfaces;
	DataHashTableKeyValue<int, MeshPoint3d*>			*hash_points;
	DataHashTableKeyValue<int, MeshDomainSurface*>		*hash_faces;
	DataHashTableKeyValue<int, MeshPoint2d*>			*hash_vertices2d;
	// array of dynamic boundary-condition data
	DataVector<std::shared_ptr<MeshBoundaryCondition>> mbc_list;
protected:
	DEquationConstTable *m_ctable;
};

#endif // !defined(MESHBREP_H__INCLUDED)
