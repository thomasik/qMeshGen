/////////////////////////////////////////////////////////////////////////////
// OctTree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(OCTTREE_H__INCLUDED)
#define OCTTREE_H__INCLUDED

#include "DPoint.h"
#include "DRect.h"
#include "MeshTetrahedron.h"
#include "Metric3dContext.h"

/**
 * This class implements a Oct-tree structure for finding of containing tetrahedron.
 */
class OctTree  
{
public:
	/// Standard constructor
	OctTree(const DBox& box);
	/// Standard destructor
	~OctTree();
public:
	/// Returns the tetrahedron nearest the given point (within the Oct tree)
	MeshTetrahedron* getNearestTetrahedron(const DPoint3d& pt, MeshTetrahedron* last_tetrahedron = nullptr) const;
	/// Removes the given reference to the tetrahedron from the Oct tree
	void removeTetrahedronLink(const MeshTetrahedron* tetrahedron);
	/// Removes the given reference to the tetrahedron from the Oct tree
	void removeTetrahedronLink(const MeshTetrahedron* tetrahedron, const DPoint3d& pt);
	/// Adds the given reference to the tetrahedron into the Oct tree
	void insertTetrahedronLink(MeshTetrahedron* tetrahedron);
	/// Adds the given reference to the tetrahedron into the Oct tree
	void insertTetrahedronLink(MeshTetrahedron* tetrahedron, const DPoint3d& pt);
	/// Return current value of max level
	int getMaxLevel() const { return max_level; }
protected:
	/// Number of tetrahedrons in this area
	int			count;
	/// Coordinates of the middle point for this area
	DPoint3d	middle;
	/// The distance from the currently nearest tetrahedron from the middle-point
	double		dist;
	/// The size of hexahedral area (as the distance dx, dy, dz from the center)
	DPoint3d	size;
	/// Preferred starting tetrahedron for this area (nearest the middle)
	MeshTetrahedron* tetrahedron;
	/// Array of sub-trees down:[2 - 3 / 0 - 1] / up:[6 - 7 / 4 - 5]
	OctTree*	parts[8];
	/// Maximum level (for stats)
	int max_level;
public:
	/// Threshold-value, determining when the leaf should be split into eight parts
	static int param_octtree_threshold;
};

#endif // !defined(OCTTREE_H__INCLUDED)



#if !defined(OCTTREEMESHPOINTS_H__INCLUDED)
#define OCTTREEMESHPOINTS_H__INCLUDED

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"

class MeshPoint3d;

/**
 * This class implements a Octree structure for storing mesh points
 */
class OctTreeMeshPoints
{
public:
	/// Standard constructor
	OctTreeMeshPoints(const DBox& box, OctTreeMeshPoints* father = nullptr);
	/// Standard destructor
	~OctTreeMeshPoints();
public:
	/// Returns the tetrahedron nearest the given point (within the Oct tree)
	//MeshTetrahedron* getNearestTetrahedron(const DPoint3d& pt, MeshTetrahedron* last_tetrahedron = nullptr) const;
	/// Removes the given reference from the Octree
	void removeMeshPointLink(MeshPoint3d* point);
	/// Adds the given reference into the Octree
	void insertMeshPointLink(MeshPoint3d* point);
	/// Whether there is any point in the proximity of points in tree
	bool anyMeshPointInProximity(Metric3dContext& mc, const DPoint3d& pt, double max_metric_distance = SQRT2) const;
	/// Whether there is any point in the proximity of points in leaf
	bool anyMeshPointInProximityForLeaf(Metric3dContext& mc, const DMPoint3d& dpt, double max_metric_distance2) const;
	/// Whether there is any point in the proximity of points in near leaves
	bool anyMeshPointInProximityForNearLeaves(Metric3dContext& mc, const DPoint3d& pt, const DMPoint3d& dpt, 
		double h_max, double max_metric_distance2) const;
	/// Find leaf for given point
	OctTreeMeshPoints* findLeaf(const DPoint3d& pt);
	/// Find leaf for given point (const)
	const OctTreeMeshPoints* findLeaf(const DPoint3d& pt) const;
protected:
	/// Number of points in this box
	DataVector<MeshPoint3d*> m_points;
	/// Coordinates of the middle point for this area
	DPoint3d	m_middle;
	/// The size of hexahedral area (as the distance dx, dy, dz from the center)
	DVector3d	m_size;
	/// Which coordinate is split
	int m_split;
	/// Array of sub-trees
	OctTreeMeshPoints*	m_parts[2];
	/// Father
	OctTreeMeshPoints*	m_father;
};

#endif // !defined(OCTTREEMESHPOINTS_H__INCLUDED)
