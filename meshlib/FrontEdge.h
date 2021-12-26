/////////////////////////////////////////////////////////////////////////////
// FrontEdge.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(FRONTEDGE_H__INCLUDED)
#define FRONTEDGE_H__INCLUDED

#include "MeshData.h"
#include "Metric2dContext.h"

class MeshEdge2d;
class MeshPoint2d;

#define ANGLE_THRESHOLD 1.707 // 1-cos(0.75*PI)
#define ANGLE_DIFFERENCE 0.134 // 1-cos(PI/6)
#define ANGLE_1 0.076 // 1-cos(PI/8)
#define ANGLE_2 0.134 // 1-cos(PI/6)

/**
 * This class represents an edge of "running front" used during triangle-to-quad conversion.
 * Each edge has pointers to the neighbouring edges (thus forming a closed cycle) and is
 * clasified according to its length, layer number and angles between this edge and its neighbours.
 */
class FrontEdge : public IndexedObject
{
public:
	enum EdgeSide {EDGE_LEFT = 0, EDGE_RIGHT = 1, EDGE_UPPER = 2 };
public:
	/// Standard constructor
	FrontEdge(MeshEdge2d* edge = nullptr, MeshPoint2d* point = nullptr, int lev = 0);
	/// Standard destructor
	~FrontEdge();
public:
	/// clear some data for faster all-delete process
	void preDeleteAll() {} // nothing here acutally
	/// Classifies this edge (angles, length, level) and returns the pointer to it
	FrontEdge* classify(Metric2dContext& mc);
	/// Returns the average lenght of this edge and of the right- or left-adjacent edge
	//double getSideAverageLength(Metric2dContext& mc, int side) const;
	/// Returns type of the edge (Type denotes whether left or right edge can be directyly used for forming quad)
	short getType() const { return type; }
	/// Sets arbitrary type for the edge
	void setType(short t){ type = t; }
	/// Returns the level of this edge
	int getLevel() const { return level; }
	/// Increases level for the level
	void incLevel(int dl = 1) { level += dl; }
	/// Sets arbitrary level for the level
	void setLevel(int l){ level = l; }
	/// Return metric gradation ratio
	double getMetricGradation() const { return metric_gradation; }
	/// Return metric length (squared) of the edge
	double getMetricLengthSquared() const { return length_squared; }
	/// Returns the pointer to the left front edge
	FrontEdge* getSideFrontEdge(int side) const { return fedges[side]; }
	/// Sets a link to the left edge
	void setSideFrontEdge(FrontEdge* v, int side);
	/// Checks whether the left or right edge (its angle) is small enough
	bool isSideFrontEdgeNear(int side) const { return angles[side] < ANGLE_THRESHOLD; }
	/// Returns the angle between this and the left or right edge
	double getSideAngleCos(int side) const { return angles[side]; }
	/// Sets the index of point between this and the left edge (0 or 1, used for mark the direction of edge)
	void setLeftIndex(int i){ index_left = i; }
	/// Returns the index of point between this and the left edge (0 or 1, used for mark the direction of edge)
	int getLeftIndex() const { return index_left; }
	/// Returns the pointer to the mesh-edge which forms this front-edge
	MeshEdge2d* getEdge() const { return m_edge; }
	/// Clear front edge (temporarily)
	void clear();
	/// Initialized the front edge by setting the assigned mesh-edge and the mesh-point (on the left side of edge)
	void setEdge(MeshEdge2d* edge, MeshPoint2d* point);
	/// Compares the rank of this front-edge with the given one (required by the class DataContainer)
	short compareTo(const FrontEdge* item) const;
	/// Description
	friend ostream& operator<<(ostream& os, const FrontEdge& fe);
protected:
	/// Pointer to the assigned mesh-edge
	MeshEdge2d *m_edge;
	/// Pointer to the front-edge adjacent at the left and right sides
	FrontEdge *fedges[2];
	/// Index of the point (within the mesh-edge) between this and the left front-edge
	int		index_left;
	/// The angle (cosinus: 0-4) between this and the left and right edge (stored for efficiency)
	double	angles[2];
	/// The number of the "front layer" (starting from 0 and increasing)
	int		level;
	/// Type of the front-edge (based on the agnles between this and neigbouring front-edges)
	short	type;
	/// Metric gradation for this front edge
	double  metric_gradation;
	/// Metric length of edge
	double length_squared;
};

#endif // !defined(FRONTEDGE_H__INCLUDED)
