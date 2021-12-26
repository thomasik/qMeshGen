/////////////////////////////////////////////////////////////////////////////
// MeshContainer3dSurface.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2010-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHCONTAINER3DSURFACE_H__INCLUDED)
#define MESHCONTAINER3DSURFACE_H__INCLUDED

#include <functional>
#include <memory>

#include "DRect.h"
#include "DataContainer.h"
#include "TagExtended.h"
#include "MeshPoint3d.h"
#include "MeshFace.h"
#include "IteratorEdge3d.h"
#include "ControlSpace3d.h"
#include "DataVector.h"
#include "DataPtrMatrix.h"
#include "MeshEdge3d.h"

#include "SurfaceParametric.h"
#include "Curve3dParametric.h"

class ControlSpace3d;
class MeshViewSet;
class MeshDomainVolume;
class Metric3dContext;

enum BorderStage { BORDER_UNKNOWN = 0, BORDER_TOPOLOGY = 1, BORDER_KNOWN = 2 };

/**
 * This class implements a container of mesh points and faces for 2.5D discretization.
 */
class MeshContainer3dSurface
{
public:
	/// Standard constructor
	MeshContainer3dSurface(int part_size );
	/// Standard destructor
	virtual ~MeshContainer3dSurface();
public:
	/// for-each face
	void forEachFace( std::function <void(MeshFace* face)> f );
	/// for-each point
	void forEachPoint( std::function <void(MeshPoint3d* point)> f );
	//// surface is just being deleted
	bool isDeleted() const { return m_points == nullptr; }
	/// show surface mesh with marking faces with sharp edges
	void showFacesWithSharpEdges(const string& label) const;
	/// show surface mesh with marking "inside-ratio" of local surfaces for vertices
	//void showLocalSurfacesQuality(Metric3dContext& mc, const string& label) const;
	/// show surface mesh with marking local surfaces
	void showLocalSurfacesForMesh(const string& label) const;
	/// check surface/surface-set counters for mesh points
	void checkSurfaceSetCounters(Metric3dContext& mc, bool log_data = true) const;
	/// check local surface approximation quality
	void checkLocalSurfaces(Metric3dContext& mc) const;
	/// Store surface mesh to .XML file
	bool storeXML(const string& fname, const string& description = "") const;
	/// Store surface mesh to .OFF file
	bool storeOFF(const string& fname) const;
	/// Convert all non-triangular faces to triangles
	void convertPolysToTriangles();
	/// Checks validity of local surface params
	bool checkLocalSurfaceParams() const;
	/// Checks validity
	bool isValid(bool allow_degenerate_faces = false) const;
	/// Clears all boundary flags from surface entities
	void clearBoundaryFlags();
	/// Sets boundary flags (edges and points) for sharp edges (f - max value for scalar product of normal vectors)
	void setBoundarySharpEdges(double fmax, TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Sets boundary flags (edges and points) for "feature" edges i.e. singe or multiply-connected
	void setBoundaryFeatureEdges();
	/// Sets boundary flags (edges and points) according to the given tag_type
	void setBoundaryTagEdges(TagExtended::TagType tag_type);
	/// Remove all tags
	void removeAllTags(TagExtended::TagType tag_type);
	/// Returns the "screenshot" of this mesh for visualization
	MeshViewSet* getViewSet(MeshViewSet* set = nullptr, 
		TagExtended::TagType color_tag_type = TagExtended::TAG_NONE) const;
	/// Returns the "screenshot" of subarea of this mesh for visualization during debug
	MeshViewSet* getDebugViewSetTopological(MeshFace* face, int radius = 2, 
		TagExtended::TagType color_tag_type = TagExtended::TAG_NONE) const;
	/// Returns the "screenshot" of subarea of this mesh for visualization during debug
	MeshViewSet* getDebugViewSetTopological(const MeshPoint3d* pt, int radius = 2, 
		TagExtended::TagType color_tag_type = TagExtended::TAG_NONE) const;
	/// Returns the "screenshot" of subarea of this mesh for visualization during debug
	MeshViewSet* getDebugViewSet(const MeshPoint3d* pt1, const MeshPoint3d* pt2 = 0, double radius = 2.0, 
		TagExtended::TagType color_tag_type = TagExtended::TAG_NONE) const;
	/// Returns the number of mesh faces with the given number of edges (i.e. triangles or quads)
	int getFacesCount(int edge_count) const;
	/// Sets a reference to the class describing the control space for mesh faces in this container
	void setControlSpace(CS3dPtr space);
	/// Switches the position of two mesh points in the container (their indices)
	void switchMeshPoints(int i, int j){ m_points->switchDataItems(i, j); }
	/// Removes all points and faces (references AND objects)
	void deleteAll();
	/// Returns the number of mesh points in the container
	int getPointsCount() const { return m_points->countInt(); }
	/// Returns the number of mesh faces in the container
	int getFacesCount() const { return m_faces->countInt(); }
	/// Returns the point (a pointer to it) at the given index within the container
	MeshPoint3d* getPointAt(int index) const { return m_points->getDataAt(index); }
	/// Returns the face (a pointer to it) at the given index within the container
	MeshFace* getFaceAt(int index) const { return m_faces->getDataAt(index); }
	/// Deletes all mesh points (both references and objects) from this container
	void deleteAllMeshPoints(){ m_points->deleteAll(); }
	/// Deletes all mesh faces (both references and objects) from this container
	void deleteAllMeshFaces(){ m_faces->deleteAll(); }
	/// Removes the reference to the point at the given index from the container
	MeshPoint3d* removeMeshPoint(int index){ return m_points->removeDataItem(index); }
	/// Removes the reference to the point
	MeshPoint3d* removeMeshPoint(MeshPoint3d* point){ return m_points->removeDataItem(point->getIndex()); }
	/// Removes the reference to the face at the given index from the container (additional ln(n) updates if heap-order is active)
	MeshFace* removeMeshFace(int index){ return m_faces->removeDataItem(index); }
	/// Removes the reference to the given face from the container
	MeshFace* removeMeshFace(MeshFace* face){ return removeMeshFace(face->getIndex()); }
	/// Adds the new reference to point into this container (no checks for duplicates)
	int addMeshPoint(MeshPoint3d *point){ return m_points->addDataItem(point); }
	/// Adds the new reference to face into this container (no checks for duplicates, additional ln(n) updates if heap-order is active)
	int addMeshFace(MeshFace *face){ return m_faces->addDataItem(face); }
	/// Returns the bounding box, containing all points and faces from this container
	const DBox& getBoundingBox( bool force_update = false );
	/// Returns the bounding box diameter, containing all points and faces from this container
	double getBoundingBoxDiameter( bool force_update = false );
	/// Returns the currently set reference to the control space
	CS3dPtr getControlSpace() const { return m_control; }
	/// Returns the iterator for browsing all edges in this mesh
	IteratorEdge3d getFirstEdge3d() const { return IteratorEdge3d(this); }
	/// Insert domain volume
	void addDomainVolume(std::shared_ptr<MeshDomainVolume> mdv){ m_volumes.add(mdv); }
	/// Removes all local surface and curve definitions and counters
	void clearLocalShapes();
	/// Search for local reparameterization surfaces (for whole mesh, or only faces with the given tag)
	//int identifyLocalSurfaces(Metric3dContext& mc, double tolerance = 0.5, 
	//	TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Search for local reparameterization curves (for boundary contours)
	//int identifyLocalCurves(Metric3dContext& mc, double tolerance = 0.5, 
	//	TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value1 = 0, int tag_value2 = 0);
	/// Search for local reparameterization curve for the given chain of edges
	Curve3dConstPtr approximateLocalCurveForChain( Metric3dContext& mc, SurfaceConstPtr surface, 
		const DataVector<MeshEdge3d*> chain_edges, const DataVector<MeshPoint3d*> chain_points, double tolerance );
	/// Search for local reparameterization curve for the given (boundary) edge, within the local surface patch
	Curve3dConstPtr approximateLocalCurve(Metric3dContext& mc, MeshEdge3d* edge, 
		SurfaceConstPtr surface, DataHashTable<MeshEdge3d*> * visited_edges = nullptr, double tolerance = 0.2);
	/// Search for local reparameterization surface for the given face, returns surface and number of ascribed faces
	SurfaceConstPtr approximateLocalSurface(Metric3dContext& mc, const DataVector<MeshFace*> &faces, 
		int & ascribed_faces, double tolerance = 0.2, int top_layers = 3, double metric_radius = 1.0,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	/// Search for local reparameterization surface for the given face, returns surface and number of ascribed faces
	SurfaceConstPtr approximateLocalSurfaceViaNormals(Metric3dContext& mc, const DataVector<MeshFace*> &faces, 
		int & ascribed_faces, double sp_min = 0.5, double tolerance = 0.2, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0);
	static void cleanAndPackPFLists( int l_outside,
			DataVector<int> & p_list, DataVector<int> & p_layers, 
			DataVector<int> & f_list, DataVector<int> & f_layers,
			DataVector<DPoint3d>& lpoints, DataVector<DPoint2d>& lparams, 
			DataVector<double>& laquality);
	void fixLocalSurfaceOrientation( SurfacePtr surface, 
		const DataVector<int> & f_list, const DataVector<int> & pref,
		const DataVector<DPoint3d>& lpoints, DataVector<DPoint2d>& lparams);
	void createLocalSurfaceDomain( Metric3dContext& mc, 
		SurfacePtr surface, const DBox& sbox, double tolerance,
		const DataVector<int> & f_list, const DataVector<int> & pref,
		const DataVector<DPoint3d>& lpoints, const DataVector<DPoint2d>& lparams,
		const DataVector<double>& laquality);
	int ascribeLocalSurface( Metric3dContext& mc, 
		SurfaceConstPtr surface, const DBox& sbox,
		const DataVector<int> & f_list,
		const DataVector<int> & p_list,
		const DataVector<DPoint2d>& lparams,
		const DataVector<double>& laquality);
	int createContoursForLocalSurface(Metric3dContext& mc, 
		SurfaceConstPtr surface, double tolerance,
		const DataVector<int> & f_list);
	void fixLocalSurfaceOrientationForInvertedFaces(
		SurfaceConstPtr surface,
		DataVector<int> & f_list,
		const DataVector<int> & pref,
		const DataVector<DPoint2d>& lparams);
	/// Tries to move inner points to fit theirs ascribed local surfaces
	//int moveInnerPointsToLocalShape(Metric3dContext& mc, TagExtended::TagType tag_type = TagExtended::TAG_NONE, 
	//	int tag_value = 0, int forbid_tag_value = 0);
	/// Tries to move border points to fit theirs ascribed local surfaces
	//int moveBorderPointsToLocalShape(Metric3dContext& mc, TagExtended::TagType tag_type = TagExtended::TAG_NONE, 
	//	int tag_value = 0, int forbid_tag_value = 0);
	/// Check average length of inner edges (with metric)
	double checkInnerEdgesLength(Metric3dContext& mc, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Check average length of border edges (with metric)
	double checkBorderEdgesLength(Metric3dContext& mc, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value1 = 1, int tag_value2 = 1);
	/// Gather set of nodes through topological neighbourhood, with face-face adjacency, no border-crossing
	bool gatherLayeredVerticesTopological(DataVector<int> & f_layers_count, 
				DataVector<int> & p_list, DataVector<int> & p_layers, 
				DataVector<int> & f_list, DataVector<int> & f_layers,
				int layers = -1, TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0) const;
	/// Gather set of nodes through geometrical neighbourhood, with point-point adjacency
	bool gatherLayeredVerticesGeometrical(Metric3dContext& mc, 
				const DPoint3d & dpt, DataVector<int> & f_layers_count, 
				DataVector<int> & p_list, DataVector<int> & p_layers, 
				DataVector<int> & f_list, DataVector<int> & f_layers,
				double radius, TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0) const;
	/// Gather set of nodes through neighbourhood, via normals, with point-point adjacency
	double gatherLayeredVerticesViaNormals(SurfaceConstPtr base_surface, double sp_min,
				DataVector<int> & f_layers_count, 
				DataVector<int> & p_list, DataVector<int> & p_layers, 
				DataVector<int> & f_list, DataVector<int> & f_layers,
				TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0) const;
	/// Gather set of nodes through neighbourhood (point-point & face-face adjacency) with surface fitting
	bool gatherLayeredVerticesForSurface(Metric3dContext& mc, 
				SurfaceConstPtr surface, DataVector<int> & f_layers_count, 
				DataVector<int> & p_list, DataVector<int> & p_layers, 
				DataVector<int> & f_list, DataVector<int> & f_layers,
				double tolerance, TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0) const;
	/// Gather set of nodes through topological neighbourhood, with face-face adjacency
	bool gatherLayeredVerticesTopological(MeshFace* start_face,	DataVector<DPoint3d> & points, 
		int layers = -1, bool cross_borders = false, int min_points = 0,
		DataVector<MeshFace*> * layer_faces = nullptr,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0) const;
	/// Gahter first layer of faces around the given point
	void gatherFirstLayerTopologicalForPoint(const MeshPoint3d* point, DataVector<int> & f_layers_count, 
				DataVector<int> & p_list, DataVector<int> & p_layers, 
				DataVector<int> & f_list, DataVector<int> & f_layers,
				TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0) const;
	/// Gahter first layer of faces for the given face (actually it's only this face(s) and its vertices)
	void gatherFirstLayerTopological(const DataVector<MeshFace*> &faces,
				DataVector<int> & f_layers_count, 
				DataVector<int> & p_list, DataVector<int> & p_layers, 
				DataVector<int> & f_list, DataVector<int> & f_layers,
				TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0) const;
	/// Gahter chain of boundary edges forming border contour, not crossing "corner"-boundary-points
	static MeshPoint3d* gatherBorderContourChain(MeshEdge3d* start_edge, MeshPoint3d* start_point, 
				DataVector<MeshEdge3d*> & edges, int layers = -1, bool only_without_curve = true);
	/// Gahter chain of boundary edges forming border contour, not crossing "corner"-boundary-points
	static MeshPoint3d* gatherBorderContourChain(MeshEdge3d* start_edge, 
				DataVector<MeshEdge3d*> & chain_edges, int layers = -1, bool only_without_curve = true);
	/// Gahter chain of boundary edges forming border contour, not crossing "corner"-boundary-points
	static MeshPoint3d* gatherBorderContourChainForSurface(MeshEdge3d* start_edge, MeshPoint3d* start_point, 
				SurfaceConstPtr surface, DataVector<MeshEdge3d*> & edges);
	/// Gahter chain of boundary edges forming border contour, not crossing "corner"-boundary-points
	static MeshPoint3d* gatherBorderContourChainForSurface(MeshEdge3d* start_edge, SurfaceConstPtr surface,
				DataVector<MeshEdge3d*> & chain_edges );
	/// Get info about border identification stage
	BorderStage getBorderStage() const { return m_border_stage; }
	/// Set info about border identification stage
	void setBorderStage(BorderStage bs) { m_border_stage = bs; }
	/// Add local surface
	void addLocalSurface(SurfaceConstPtr surface) { m_local_surfaces.add(surface); }
	/// Add local curve
	void addLocalCurve(Curve3dConstPtr curve) { m_local_curves.add(curve); }
	/// Any local shape (surface or curve) ?
	bool hasAnyLocalShape() const { return m_local_surfaces.notEmpty() || m_local_curves.notEmpty(); }
public:
	/// Show faces/points in the local sub-domain - just for test
	void showDebugLocalFaces(const string& caption, DataVector<int> & p_list, DataVector<int> & p_layers, 
			DataVector<int> & f_list, DataVector<int> & f_layers, int l_outside, DataVector<int> * pref = nullptr ) const;
protected:
	/// Set of mesh points
	DataContainer<MeshPoint3d>	*m_points;
	/// Set of mesh faces
	DataContainer<MeshFace>		*m_faces;
	/// Factor of size-increase for memory management
	int	m_part_size;
	/// The reference to the control space (which governs the size and shape of created mesh-faces)
	CS3dPtr m_control;
	/// Vector of volume blocks, for storing material id
	DataVector<std::shared_ptr<MeshDomainVolume>> m_volumes;
	/// Container for local parameterization surfaces
	DataVector<SurfaceConstPtr> m_local_surfaces;
	/// Container for local parameterization curves
	DataVector<Curve3dConstPtr> m_local_curves;
	/// Stage of border-identification 
	BorderStage m_border_stage;
	DBox m_box;
	double m_box_diameter;
public:
	enum SurfaceDomainType { SDOMAIN_HULL = 0, SDOMAIN_FACES = 1,  SDOMAIN_FACES_PLANAR = 2 };
	static int param_surface_domain_type;
};

#endif // !defined(MESHCONTAINER2D_H__INCLUDED)
