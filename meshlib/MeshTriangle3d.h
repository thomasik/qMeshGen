/////////////////////////////////////////////////////////////////////////////
// MeshTriangle3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHTRIANGLE3D_H__INCLUDED)
#define MESHTRIANGLE3D_H__INCLUDED

#include <memory>

#include "MeshFace.h"
#include "MeshPoint3d.h"

class MeshBlock;
class SurfaceParametric;

/**
 * This class implements a three-dimensional, triangular mesh face
 *	with several required methods.
 */
class MeshTriangle3d : public MeshFace  
{
public:
	/// Standard constructor
	MeshTriangle3d(MeshPoint3d *p1, MeshPoint3d *p2, MeshPoint3d *p3, MeshBlock* bounded_block = nullptr);
	/// Standard destructor
	virtual ~MeshTriangle3d();
public:
	// check, whether this face is valid with respect to the given surface
	virtual bool valid(std::shared_ptr<const SurfaceParametric> surface ) const override;
	/// Split face into (mesh) triangles
	virtual bool splitToTriangles( DataVector< MeshFace* > & split_faces ) const override;
	/// Split face into triangles
	virtual bool splitToTriangles( DataVector< DTriangle3d > & triangles ) const override;
	/// returns shape quality for a face
	virtual double getShapeQuality(Metric3dContext& mc) const override;
	/// returns shape quality for a face
	virtual double getShapeQuality() const override;
	/// remove point, return resulting face (the same, or a new one)
	virtual MeshFace* removePoint(const MeshPoint3d* /* point */) override { return nullptr; }
	/// Creates and returns a copy of this face (+ whole connectivity)
	virtual MeshFace* clone() const override;
	/// Returns the type of this object
	virtual ElementType getType() const override { return FACE_TRIANGLE; }
public:
	/// Area of a triangular face
	double area() const;
	/// Area of a triangular face
	double area(Metric3dContext& mc) const;
	/// Alpha quality coefficient of a triangular face
	double alphaQuality() const;
	/// Swaps the i-th edge of this triangle with the adjacent one (can be conditionally only)
	MeshEdge3d* swapWithNeighbour(Metric3dContext& mc, SurfaceConstPtr surface, int i, bool improve_only, bool local_metric = false, 
		TagExtended::TagType side_tag = TagExtended::TAG_NONE);
private:
	/// Three-element array of mesh edges
	MeshEdge3d*		triangle_edges[3];
	/// Three-element array of vertices (mesh points)
	MeshPoint3d*	triangle_points[3];
};

#endif // !defined(MESHTRIANGLE3D_H__INCLUDED)
