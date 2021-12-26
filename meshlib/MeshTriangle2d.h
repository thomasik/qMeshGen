/////////////////////////////////////////////////////////////////////////////
// MeshTriangle2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHTRIANGLE_H__INCLUDED)
#define MESHTRIANGLE_H__INCLUDED

#include "common.h"
#include "MeshData.h"
#include "MeshElement.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "DPoint.h"

/**
 * This class implements a two-dimensional, triangular mesh element
 *	with several required methods.
 */
class MeshTriangle2d : public MeshElement  
{
public:
	/// Standard constructor
	MeshTriangle2d(MeshPoint2d *p1, MeshPoint2d *p2, MeshPoint2d *p3);
	/// Standard destructor
	virtual ~MeshTriangle2d();
public:
	/// Returns the third vertex
	MeshPoint2d* getOtherPoint(const MeshPoint2d* p0, const MeshPoint2d* p1) const;
	/// Returns the coordinates of the outer circle center for this triangle
	DMPoint2d getOuterCircleCenter(Metric2dContext& mc, bool local_metric = true);
	/// Swaps the i-th edge of this triangle with the adjacent one (check already done before)
	bool swapWithNeighbourNoCheck(int i);
	/// Swaps the i-th edge of this triangle with the adjacent one (can be conditionally only)
	bool swapWithNeighbour(Metric2dContext& mc, int i, bool improve_only, bool local_metric = false, 
		TagExtended::TagType side_tag = TagExtended::TAG_NONE);
	/// Swaps the given edge of this triangle with the adjacent one (can be conditionally only)
	bool swapWithNeighbour(Metric2dContext& mc, const MeshEdge2d* edge, bool improve_only, bool local_metric = false);
	/// Checks whether the swapPIng of the i-th edge of this triangle with the adjacent one improves the overal quality
	bool shouldSwap(Metric2dContext& mc, int i, bool local_metric = false) const;
	/// Returns the reference to neighbouring triangle in the direction given by the point
	MeshTriangle2d* getNeighbourInDirection(const DPoint2d& dpt, bool mind_border = false) const;
	/// Checks whether the i-th edge of the triangle belongs to the boundary
	bool isBorder(int i) const{ return edges[i]->isBorder(); }
	/// Returns the weight for the i-th vertex in this triangle
	//double getVertexWeight(int i) const { return getPoint(i)->getWeight(); }
	/// Checks whether the given mesh point is withing the outer circle of this triangle (different metric may be used)
	bool isPointInOuterCircle(Metric2dContext& mc, const MeshPoint2d *point, bool local_metric = true) const 
		{ return isPointInOuterCircle(mc, point->getCoordinates(), local_metric); }
	/// Checks whether the given point is withing the outer circle of this triangle (different metric may be used)
	bool isPointInOuterCircle(Metric2dContext& mc, const DPoint2d& point, bool local_metric = true) const;
	/// Returns the length of the radius of the inner circle for this triangle (different metric may be used)
	double getInnerCircleRadius(Metric2dContext& mc, bool local_metric = true) const;
	/// Checks whether the swapPIng of common edge for two triangles (given be 4 points) is possible and advantagous
	bool checkSwapCriterion(const DMPoint2d& dpt1, const DMPoint2d& dpt2, 
		const DMPoint2d& dpt3, const DMPoint2d& dpt4) const;
	/// Returns the triangle containing the given point, starting from this one
	MeshTriangle2d* findTriangleByNeighbours(const DPoint2d& mpt, bool mind_border = false);
	/// Returns the triangle containing the given point, starting from this one
	MeshTriangle2d* findTriangleByNeighbours(MeshPoint2d* point, bool mind_border = false){
		return findTriangleByNeighbours(point->getCoordinates(), mind_border); }
public:
	/// clear some data for faster all-delete process
	virtual void preDeleteAll();
	/// add this element as a lump of mass (respecting control space)
	static void addForInertialCenterAdaptiveSplitLongestEdge(Metric2dContext& mc, 
		const DPoint2d& pt0, const DPoint2d& pt1, const DPoint2d& pt2, 
		DPoint2d& inertial_center, double& total_mass, int level = 0);
	/// add this element as a lump of mass (respecting control space)
	virtual void addForInertialCenter(Metric2dContext& mc, DPoint2d& inertial_center, double& total_mass) const;
	/// add this element as a lump of mass (respecting control space)
	static void addForInertialMomentsAdaptiveSplitLongestEdge(Metric2dContext& mc, 
		const DPoint2d& pt0, const DPoint2d& pt1, const DPoint2d& pt2, 
		const DPoint2d& inertial_center, DMatrix2d& inertial_moments, int level = 0);
	/// add this element as a lump of mass (respecting control space)
	virtual void addForInertialMoments(Metric2dContext& mc, const DPoint2d& inertial_center, DMatrix2d& inertial_moments) const;
	/// Returns the area of the triangle
	virtual double getArea(Metric2dContext& mc, bool local_metric = true) const;
	/// Returns the area of the triangle
	virtual double getAreaNoMetric() const;
	/// Returns the inner angle at the i-th vertex of this element
	virtual double getAngle(Metric2dContext& mc, int i, bool local_metric = true) const;
	/// Counts (and stores) the quality of the triangle (for global or explicite criterion)
	virtual double countQuality(Metric2dContext& mc, bool mind_area = false, int criterion = -1);
	/// Sets and returns the metric difference of this element
	virtual double countMetricDiffQuality(Metric2dContext& mc); 
	/// Returns the metric difference of this element
	virtual double countMetricDiff(Metric2dContext& mc) const;
	/// Returns "mean ratio" coefficient quality
	virtual double getMeanRatio(Metric2dContext& /*mc*/, bool /*ext_metric*/ = false) const;
	/// Returns the quality of the triangle for the "alpha" criterion
	double getAlphaQuality(Metric2dContext& mc, bool local_metric = true) const;
	/// Stores the GRD textual description of triangle into the given file
	virtual void exportToGRD(FILE* file) const;
	/// Checks whether this triangle contains the given point
	virtual bool isPointInside(const DPoint2d& mpt) const;
	/// Checks whether this triangle is inverted (wrongly orientated)
	virtual bool isInverted() const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return ELEMENT_MESH_TRIANGLE; }
	/// Returns the optimal area of triangle ( |edge| = 1.0 )
	double getOptimumArea() const { return SQRT3/4.0; }
	/// Returns the optimal radius length (squared) of the circumscribed circle  ( |edge| = 1.0 )
	double getOptimumCircumradius2() const { return 1.0/3.0; }
public:
	/// Test if triangle is OK
	bool isOK() const { return (getPoint(0) != getPoint(1)) && (getPoint(0) != getPoint(2)) &&
		(getPoint(1) != getPoint(2)); }
public:
	/// Method of assessing the quality of triangle
	static int param_quality_criterion;
	/// Method of metric selection for assessing the quality of triangle
	static int param_quality_metric_selection;
	/// Types of counters
//	enum COUNTER_TYPE { METRIC_QUALITY_1 = 0, METRIC_QUALITY_4 = 1, TRIANGLE_SEARCHING = 2 };
	/// Counter of metric calls for triangle quality assesment 
	/// (0 for one and 1 for four points mode and 2 for containing triangle searching)
//	static unsigned int m_counters[3];
	/// Get and clear counter
//	static unsigned int getAndClearCounter(int ctype) { unsigned int x = m_counters[ctype]; m_counters[ctype] = 0; return x; }
private:
	/// Three-element array of mesh edges
	MeshEdge2d*	triangle_edges[3];
	/// Three-element array of vertices (mesh points)
	MeshPoint2d*	triangle_points[3];
};

#endif // !defined(MESHTRIANGLE_H__INCLUDED)
