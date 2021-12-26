/////////////////////////////////////////////////////////////////////////////
// MeshPoly3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHPOLY3D_H__INCLUDED)
#define MESHPOLY3D_H__INCLUDED

#include "MeshFace.h"
#include "MeshPoint3d.h"

/**
 * This class implements a three-dimensional, polygonial mesh face
 *	with several required methods.
 */
class MeshPoly3d : public MeshFace  
{
public:
	/// Standard constructor
	MeshPoly3d(const DataVector<MeshPoint3d*> & poly_points, MeshBlock* bounded_block = nullptr);
	/// Standard constructor
	MeshPoly3d(int poly_count, MeshPoint3d** poly_points, MeshBlock* bounded_block = nullptr);
	/// Standard destructor
	virtual ~MeshPoly3d();
public:
	// check, whether this face is valid with respect to the given surface
	virtual bool valid( SurfaceConstPtr surface ) const override;
	/// Split face into triangles
	virtual bool splitToTriangles( DataVector< DTriangle3d > & triangles ) const override;
	/// Split face into triangles
	virtual bool splitToTriangles( DataVector<MeshFace*> & split_faces ) const override;
	/// returns shape quality for a face
	virtual double getShapeQuality(Metric3dContext& mc) const override;
	/// returns shape quality for a face
	virtual double getShapeQuality() const override;
	/// remove point, return resulting face (the same, or a new one)
	virtual MeshFace* removePoint(const MeshPoint3d* point) override;
	/// Creates and returns a copy of this face (+ whole connectivity)
	virtual MeshFace* clone() const override;
	/// Returns the type of this object
	virtual ElementType getType() const override { return FACE_POLY; }
	/// Returns the vector normal to this face if possible
//	virtual bool checkAndGetNormalVector(DVector3d& vn) const;
	/// Returns the vector normal to this face
//	virtual DVector3d getNormalVector() const;
	/// Returns the description of this face for visualization
	virtual std::shared_ptr<MeshViewFaceData> getViewData(double shrink_ratio = 1.0, bool proper_orientation = true) const override;
private:
	static double shapeQuality( const DataVector<DPoint3d> & poly );
};

#endif // !defined(MESHPOLY3D_H__INCLUDED)
