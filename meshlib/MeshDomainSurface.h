/////////////////////////////////////////////////////////////////////////////
// MeshDomainSurface.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHDOMAINSURFACE_H__INCLUDED)
#define MESHDOMAINSURFACE_H__INCLUDED

#include <memory>

#include "common.h"
#include "MeshData.h"
#include "DPoint.h"
#include "MeshFace.h"
#include "DRect.h"
#include "DataVector.h"
#include "TagExtended.h"
#include "SurfaceParametric.h"
#include "ControlSpace2d.h"
#include "DataHashTable.h"

class MeshContainer2d;
class Metric2dContext;
class MeshDomainEdge3d;
class MeshDomainVolume;
class MeshViewSet;
class MeshContainer3d;
class MeshContainer3dSurface;

/**
 * This class implements the surface-patch describing the boundary of the domain to discretize.
 */
class MeshDomainSurface : public MeshFace  
{
public:
	/// Standard constructor
	MeshDomainSurface(int ct, MeshPoint3d** new_points, int ect, MeshEdge3d** new_edges, 
		SurfaceConstPtr surface, MeshContainer2d* boundary);
	/// Standard constructor
	MeshDomainSurface(DataVector<MeshPoint3d*> &new_points, DataVector<MeshDomainEdge3d*> &new_edges, 
		SurfaceConstPtr surface, MeshContainer2d* boundary);
	/// Standard destructor
	virtual ~MeshDomainSurface();
public:
	DVector3d getSurfaceNormalForMidEdge(const MeshPoint3d* mp0, const MeshPoint3d* mp1) const override;
	static void autoSetBoundaryTags(MeshContainer3d* mesh);
	/// Creates and returns a copy of this face (+ whole connectivity)
	virtual MeshFace* clone() const { assert(false); return nullptr; }
	/// classify free-point for this domain-surface
	bool offerFreePoint(std::shared_ptr<MeshPoint3d> fpoint);
	/// Returns the "screenshot" of this domain-surface for visualization
	MeshViewSet* getViewSet(MeshViewSet* set = nullptr) const;
	/// Removes link to this edge (fo delete-all phase) - returns true if last one
	virtual bool removeEdgeLink(MeshEdge3d* edge);
	/// Stores statistics about the discretization of this domain-surface into the given text-file
	void assessQuality(ofstream& file) const;
	/// Creates an EPS-image of this surface-patch (will be changed)
	void storeEPS(int id = 0) const;
	/// Calls several smoothing-quad procedures ("steps" times) for the discretization of this surface
	bool smoothenQuads(int steps = 1, int method = MeshData::SM_LAPLACE);
	/// Calls the specifiad conversion procedure which transforms the triangulation of this surface into quadrilateral mesh
	int convertToQuads(int method, int max_ct = -1);
	/// Calls several smoothing-triangle procedures ("steps" times) for the discretization of this surface
	bool smoothen(int steps = 1,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 0, 
		int method = MeshData::SM_LAPLACE_MIXED | MeshData::SM_TOP_SWAP 
		| MeshData::SM_DEL_SWAP_COMPLETE);
	/// Performs the Delaunay-triangulation of boundary nodes (with recovery of boundary edges)
	int triangulateBoundary(Metric2dContext& mc);
	/// Performs the triangulation of this surface on the basis of boundary nodes and control space
	int triangulate();
	/// Creates the initial boundary-discretization
	int createBoundaryMesh();
	/// Removes the created discretization of this surface
	void clearDiscretization();
	/// Clear calculated control space for this surface
	void clearControlSpace();
	/// Returns the bounding cubicoid containing the whole domain-surface (and discretization)
	virtual DBox getBoundingBox() const;
	/// Returns the reference to the associated boundary description (discretized)
	MeshContainer2d* getBoundary() const { return m_boundary; }
	/// Returns the reference to the boundary mesh (discretized)
	MeshContainer2d* getBoundaryMesh() const { return m_boundary_mesh; }
	/// Returns the reference to the associated discretization (if created)
	MeshContainer2d* getMesh() const { return m_mesh; }
	/// Returns the type of this object
	virtual ElementType getType() const { return FACE_DOMAIN; }
	/// Return the underlying parametric surface for this domain patch
	std::shared_ptr<const SurfaceParametric> getBaseSurface() const { return m_surface; }
	/// Sets the surface as underlying for this patch
	void setBaseSurface(std::shared_ptr<const SurfaceParametric> s) { m_surface = s; }
	/// Return the user control space for this domain patch
	std::shared_ptr<ControlSpace2d> getUserControlSpace() const { return m_user_space; }
	/// Sets the user control space for this patch
	void setUserControlSpace(std::shared_ptr<ControlSpace2d> s) { m_user_space = s; }
	/// Prepare 3d mesh entites from the 2d mesh on the surface
	bool addToSurfaceMesh(MeshContainer3dSurface * smesh, 
		DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> & hbpoints, MeshDomainVolume* mdv);
	/// Returns the number of edges for this domain surface (!!! different from the number of points !!!)
	virtual int getEdgeCount() const { return edge_count; }
protected:
	/// Auxiliary structure for polyline-area-grid creation
	struct PolySegment {
		PolySegment(const DPoint2d& _pt0 = DPoint2d::zero, 
			const DPoint2d& _pt1 = DPoint2d::zero, int la = -2, int ra = -2) 
			: pt0(_pt0), pt1(_pt1), left_area(la), right_area(ra), cross(0.0), incident(0) {}
		DPoint2d pt0, pt1;
		/// material info
		int left_area, right_area;
		/// value of X or Y for crossing the current scanning line
		double cross;
		/// number of segment-vertices lying directly on the current scanning line
		int incident;
	};
	/// Automatically creates the plane surface if only boundary segments are available
	bool createPlaneSurface();
	/// Automatically creates complete boundary description (if it was given by a set of points only)
	bool createSimpleBoundary();
protected:
	/// Number of edges (can be different from the count of points for domains with inner-boundary-edges)
	int edge_count;
	/// Pointer to the user control space
	std::shared_ptr<ControlSpace2d> m_user_space;
	/// Pointer to the base surface
	std::shared_ptr<const SurfaceParametric> m_surface;
	/// Pointer to the boundary description (in the parameter space: edges, areas)
	MeshContainer2d* m_boundary;
	/// Pointer to the boundary mesh (in the parameter space: discrete edges, areas)
	MeshContainer2d* m_boundary_mesh;
	/// Pointer to the surface mesh (discretization) for this domain_face: discrete edges, elements
	MeshContainer2d* m_mesh;
	/// Most recently used quality criterion for assessing (visualization) of the quality of elements for discretization
	int m_last_quality_mode;
	/// FreePoints
	DataVector<std::shared_ptr<MeshPoint3d>> m_freepoints;
};

#endif // !defined(MESHDOMAINSURFACE_H__INCLUDED)
