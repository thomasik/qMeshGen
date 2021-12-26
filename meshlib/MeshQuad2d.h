/////////////////////////////////////////////////////////////////////////////
// MeshQuad2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHQUAD_H__INCLUDED)
#define MESHQUAD_H__INCLUDED

#include "common.h"
#include "MeshElement.h"

class MeshPoint2d;
class MeshEdge2d;
class Metric2dContext;

/**
 * This class implements a two-dimensional, quadrilateral mesh element
 *	with several required methods.
 */
class MeshQuad2d : public MeshElement  
{
public:
	/// Standard constructor
	MeshQuad2d(MeshPoint2d *p1, MeshPoint2d *p2, MeshPoint2d *p3, MeshPoint2d *p4);
	/// Standard destructor
	virtual ~MeshQuad2d();
public:
	/// clear some data for faster all-delete process
	virtual void preDeleteAll();
	/// Counts (and buffers) the quality of this quad
	virtual double countQuality(Metric2dContext& mc, bool mind_area = false, int criterion = -1);
	/// Minimum area for quad (in unitary space)
	double getMinArea() const { return 0.5; } 
	/// Returns the inner angle at the i-th vertex (different metric may be used)
	virtual double getAngle(Metric2dContext& mc, int i, bool local_metric = true) const;
	/// Returns the quality of this quad for the "alpha" criterion (different metric may be used)
	virtual double getAlphaQuality(Metric2dContext& mc, bool local_metric = true) const;
	/// Returns the area of this quad(different metric may be used)
	virtual double getArea(Metric2dContext& mc, bool local_metric = true) const;
	/// Returns the area of the triangle
	virtual double getAreaNoMetric() const;
	/// Stores the GRD textual description of quad into the given file
	virtual void exportToGRD(FILE* file) const;
	/// Checks whether this quad contains the given point
	virtual bool isPointInside(const DPoint2d& mpt) const;
	/// Checks wheter the quad is inverted (crossing edges or wrong orientation)
	virtual bool isInverted() const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return ELEMENT_MESH_QUAD; }
private:
	/// Four-element array of mesh edges
	MeshEdge2d*	quad_edges[4];
	/// Four-element array of vertices (mesh points)
	MeshPoint2d*	quad_points[4];
};

#endif // !defined(MESHQUAD_H__INCLUDED)
