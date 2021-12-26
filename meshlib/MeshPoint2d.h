/////////////////////////////////////////////////////////////////////////////
// MeshPoint2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHPOINT_H__INCLUDED)
#define MESHPOINT_H__INCLUDED

#include <memory>

#include "MeshData.h"
#include "DPoint.h"
#include "TagBorder.h"
#include "TagExtended.h"

class MeshEdge2d;
struct MeshViewPointData;
class SurfaceParametric;
class Metric2dContext;

/**
 * This class implements a mesh point in 2D with its 
 *	coordinates, incident mesh edges and additional data.
 */
class MeshPoint2d : public IndexedObject, public TagBorder, public TagExtended
{
public:
	/// Standard constructor
	MeshPoint2d(double x, double y);
	/// Standard constructor
	MeshPoint2d(const DPoint2d& pt);
	/// Copying constructor
	MeshPoint2d(const MeshPoint2d& point);
	/// Copying constructor
	MeshPoint2d(const MeshPoint2d* point);
	/// Standard destructor
	~MeshPoint2d();
public:
	/// Counts curvature for this point (based on its neighbours and Hessian approximation)
//	bool countHesjanCurvature(int ct, DPoint3d& curvature, int method = 0) const;
	/// Counts curvature for this point (based on its neighbours and Hessian approximation)
//	bool countHesjanCurvatureWithBorder(int ct, DPoint3d& curvature, int method = 0) const;
	/// Counts Hesjan structure for this point (based on its neighbours)
//	DHesjan countHesjan(int ct, bool relative = false, double ratio = 1.0, int method = 0);
	/// clear some data for faster all-delete process
	void preDeleteAll();
	/// Returns the description of this point for visualization
	std::shared_ptr<MeshViewPointData> getViewData(std::shared_ptr<const SurfaceParametric> surface = nullptr) const;
	/// Returns the number of elements incident to this point, with the given number of edges
	int getElementsCount(int edges_ct) const;
	/// Adds a reference to a mesh edge into the adjacency array (no duplication test!)
	void addEdge(MeshEdge2d* edge) { edges.add(edge); }
	/// Returns the reference to a mesh edge at the given index in the adjacency array
	MeshEdge2d* getEdge(int i) const { return edges.get(i); }
	/// Returns the reference to the mesh edge, which joins this point with the given one (nullptr if doesn't exist)
	MeshEdge2d* getEdgeToPoint(const MeshPoint2d* point) const;
	/// Removes the reference to the mesh edge from the adjacency array
	void removeEdge(MeshEdge2d* const edge) { edges.remove(edge); }
	/// Returns the coordinates of this point
	const DPoint2d& getCoordinates() const { return coord; }
	/// Returns the real space coordinates of this point
	const DPoint2d getRealCoordinates(Metric2dContext& mc) const;
	/// Returns the metric coordinates of this point
	const DMPoint2d getMetricCoordinates(Metric2dContext& mc);
	/// Returns the rank of this point (count of adjacent edges)
	int getRank() const { return (int)edges.countInt(); }
	/// Set the 2D coordinates for this mesh point
	void setCoordinates(const DPoint2d& pt) { coord = pt; }
	/// Stores the GRD textual description of this point into the given file
	void exportToGRD(FILE* file) const { fprintf(file, "%8G\t%8G\n", coord.x, coord.y); }
public:
	/// Method of assessing the quality of triangle
	static int param_max_node_count;
private:
	/// The coordinates of this point
	DPoint2d	coord;
	/// adjacent edges
	DataVector<MeshEdge2d*> edges;
};

#endif // !defined(MESHPOINT_H__INCLUDED)
