/////////////////////////////////////////////////////////////////////////////
// MeshDomainVolume.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHDOMAINVOLUME_H__INCLUDED)
#define MESHDOMAINVOLUME_H__INCLUDED

#include <memory>

#include "MeshBlock.h"
#include "DRect.h"
#include "MeshData.h"
#include "DataVector.h"
#include "ControlSpace3d.h"
#include "TagExtended.h"

class MeshFace;
class MeshContainer2d;
class MeshContainer3d;
class MeshContainer3dSurface;
class MeshPoint3d;
class MeshViewSet;

/**
 * This class implements the volume-block describing the boundary of the domain to discretize.
 */
class MeshDomainVolume : public MeshBlock  
{
public:
	/// Standard constructor
	MeshDomainVolume(DataVector<MeshFace*> &new_faces, 
		DataVector<MeshPoint3d*> &new_points, 
		DataVector<bool> &orientation);
	/// Standard constructor
	MeshDomainVolume(DataVector<MeshFace*> &new_faces, 
		DataVector<MeshPoint3d*> &new_points);
	/// Constructor for only 3d mesh
	MeshDomainVolume(MeshContainer3d* mesh3d);
	/// Constructor for only 3d surface mesh
	MeshDomainVolume(MeshContainer3dSurface* surface_mesh);
	/// Standard (empty) constructor
	MeshDomainVolume();
	/// Standard destructor
	virtual ~MeshDomainVolume();
public:
	/// Returns the "screenshot" of this domain-volume for visualization (edges of surface_mesh)
	MeshViewSet* getViewSet(MeshViewSet* set = nullptr) const;
	/// Transforms volume mesh into surface mesh (removes tetrahedra and inner entities)
	MeshContainer3dSurface* cutSurfaceMeshFromVolumeMesh();
	/// Creates surface mesh from volume mesh (copy boundary faces and entities)
	MeshContainer3dSurface* copySurfaceMeshFromVolumeMesh();
	/// classify free-point for this volume-block
	bool offerFreePoint(std::shared_ptr<MeshPoint3d> fpoint);
	/// clear some data for faster all-delete process
	virtual void preDeleteAll();
	/// Checks if control space is consistent with close boundary edges (and adjusts if necessary)
	bool checkControlForCloseBoundaryFaces(Metric3dContext& mc) const;
	/// Checks if control space is consistent with boundary edges (and adjusts if necessary)
	bool checkControlAtBoundary(Metric3dContext& mc) const;
	/// Calls several tetrahedra-smoothing procedures ("steps" times) for the discretization of this volume
	bool smoothenVolumeMesh(int steps = 1, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1, 
		int method = MeshData::SM3_LAPLACE_MIXED | MeshData::SM3_SWAP_COMPLETE);
	/// Stores existing surface/3D mesh into Amira text file (points + incidency)
	void storeAmiraMesh(const string& fname, int index = 0) const;
	/// Stores surface mesh points into text file
	void storeSurfacePoints(const string& fname, int index = 0) const;
	/// Stores existing volume mesh into text file with Abaqus format
	bool storeAbaqus(const string& fname, const string& name = "MESH3D") const;
	/// Stores existing surface mesh into text file (points + incidency)
	void storeSurfaceMeshTxt(const string& fname, int index = 0) const;
	/// Performs the Delaunay-tetrahedralization of boundary nodes (with recovery of boundary faces/edges)
	int discretizeUsingBoundary(Metric3dContext& mc);
	/// Discretizes this domain volume into tetrahedra on the basis of boundary faces (and control space)
	int createTetrahedralMesh();
	/// Creates the initial control space
	CS3dPtr createInitialControlSpace();
	/// Prepares the 3D boundary meshes for faces and the initial ACS
	int prepareBoundaryMesh();
	/// Removes the created discretization of this volume
	void clearDiscretization();
	/// Returns the bounding cubicoid enclosing the whole domain-block
	virtual DBox getBoundingBox() const;
	/// Returns the type of this object
	virtual ElementType getType() const { return BLOCK_DOMAIN; }
	/// Returns the reference to the associated surface discretization (if available)
	MeshContainer3dSurface* getSurfaceMesh() const { return m_surface_mesh; }
	/// Sets the reference to the associated surface discretization
	void setSurfaceMesh(MeshContainer3dSurface* mesh);
	/// Returns the reference to the associated discretization (if created)
	MeshContainer3d* getMesh() const { return m_mesh; }
	/// Removes and returns the reference to the associated discretization (if created)
	MeshContainer3d* removeMesh();
	/// Sets the reference to the associated discretization
	void setMesh(MeshContainer3d* mesh);
	/// Return the current control space for this domain block
	CS3dPtr getControlSpace() const { return m_control_space; }
	/// Clears current control space for this domain block
	void clearControlSpace();
	/// Return the user control space for this domain block
	CS3dPtr getUserControlSpace() const { return m_user_space; }
	/// Sets the user control space for this block
	void setUserControlSpace(CS3dPtr s) { m_user_space = s; }
	/// Sets the control space for this block
	void setControlSpace(CS3dPtr s) { m_control_space = s; }
	/// Not used in this derived class
	virtual double getVolume(Metric3dContext& /* mc */, bool /* local_metric */ = true) const { assert(false); return 0.0; }
	/// Not used in this derived class
	virtual double getVolumeNoMetric() const { assert(false); return 0.0; }
	/// Create block-faces connections
	virtual void attachToFaces() { assert(false); }
	/// Removes block-faces connections
	virtual void detachFromFaces() { assert(false); }
	/// Returns the set of free-points
	std::shared_ptr<DataVector<std::shared_ptr<MeshPoint3d>>> getFreePoints() { return m_freepoints; }
	/// Remove and join the volume mesh from the given mdv with this one
	bool combineVolumeMeshFrom(MeshDomainVolume* mdv);
protected:
	/// Pointer to the surface mesh (discretization) for the whole surface of this volume
	MeshContainer3dSurface* m_surface_mesh;
	/// Pointer to the volume mesh (discretization) for this domain_volume
	MeshContainer3d* m_mesh;
	/// Pointer to the user control space
	CS3dPtr m_user_space;
	/// Pointer to the control space for this domain volume
	CS3dPtr m_control_space;
	/// FreePoints
	std::shared_ptr<DataVector<std::shared_ptr<MeshPoint3d>>> m_freepoints
		= std::make_shared<DataVector<std::shared_ptr<MeshPoint3d>>>();
public:
	static double param_min_volume_for_simple_convex;
};


#endif // !defined(MESHDOMAINVOLUME_H__INCLUDED)
