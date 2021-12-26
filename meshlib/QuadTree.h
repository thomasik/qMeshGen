/////////////////////////////////////////////////////////////////////////////
// QuadTree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(QUADTREE_H__INCLUDED)
#define QUADTREE_H__INCLUDED

#include "DPoint.h"
#include "MeshTriangle2d.h"

/**
 * This class implements a quad-tree structure for finding of containing triangle.
 */
class QuadTree  
{
public:
	/// Standard constructor
	QuadTree(const DRect& rect);
	/// Standard destructor
	~QuadTree();
public:
	/// Returns the triangle nearest the given point (within the quad tree)
	MeshTriangle2d* getNearestTriangle(const DPoint2d& point, MeshTriangle2d* last_triangle = nullptr) const;
	/// Removes the given reference to the triangle from the quad tree
	void removeTriangleLink(const MeshTriangle2d* triangle);
	/// Removes the given reference to the triangle from the quad tree
	void removeTriangleLink(const MeshTriangle2d* triangle, const DPoint2d& pt);
	/// Adds the given reference to the triangle from the quad tree
	void insertTriangleLink(MeshTriangle2d* triangle);
	/// Adds the given reference to the triangle from the quad tree
	void insertTriangleLink(MeshTriangle2d* triangle, const DPoint2d& pt);
	/// Return current value of max level
	int getMaxLevel() const { return max_level; }
protected:
	/// Number of triangles in this area
	int			count;
	/// Coordinates of the middle point for this area
	DPoint2d	middle;
	/// The distance from the currently nearest triangle from the middle-point
	double		dist;
	/// The size of rectangular area (as the distance dx, dy from the center)
	DPoint2d	size;
	/// Preferred starting triangle for this area (nearest the middle)
	MeshTriangle2d* triangle;
	/// Array of sub-trees [0 - 1 / 2 - 3]
	QuadTree*	parts[4];
	/// Maximum level (for stats)
	int max_level;
public:
	/// Threshold-value, determining when the leaf should be split into four parts
	static int param_qtree_threshold;
};

#endif // !defined(QUADTREE_H__INCLUDED)
