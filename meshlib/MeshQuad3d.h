/////////////////////////////////////////////////////////////////////////////
// MeshQuad3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHQUAD3D_H__INCLUDED)
#define MESHQUAD3D_H__INCLUDED

#include <memory>

#include "MeshFace.h"
#include "MeshPoint3d.h"

/**
 * This class implements a three-dimensional, quadrilateral mesh face
 *	with several required methods.
 */
class MeshQuad3d : public MeshFace  
{
public:
	/// Standard constructor
	MeshQuad3d(MeshPoint3d *p1, MeshPoint3d *p2, MeshPoint3d *p3, MeshPoint3d *p4, MeshBlock* bounded_block = nullptr);
	/// Standard destructor
	virtual ~MeshQuad3d();
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
	virtual MeshFace* removePoint(const MeshPoint3d* point) override;
	/// Creates and returns a copy of this face (+ whole connectivity)
	virtual MeshFace* clone() const override;
	/// Returns the type of this object
	virtual ElementType getType() const override { return FACE_QUAD; }
private:
	/// Four-element array of mesh edges
	MeshEdge3d*		quad_edges[4];
	/// Four-element array of vertices (mesh points)
	MeshPoint3d*	quad_points[4];
};

#endif // !defined(MESHQUAD3D_H__INCLUDED)
