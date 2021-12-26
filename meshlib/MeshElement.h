/////////////////////////////////////////////////////////////////////////////
// MeshElement.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHELEMENT_H__INCLUDED)
#define MESHELEMENT_H__INCLUDED

#include <memory>

#include "MeshData.h"
#include "DataVector.h"
#include "TagExtended.h"
#include "DPoint.h"

class MeshPoint2d;
class MeshEdge2d;
class SurfaceParametric;
class Metric2dContext;
struct MeshViewFaceData;
class DMatrix2d;

/**
 * This class is a base, abstract class for a two-dimensional mesh element
 *	such as triangles or quads with some required methods.
 */
class MeshElement : public IndexedObject, public TagExtended
{
public:
	/// Standard constructor
	MeshElement(int _count, MeshEdge2d**_edges = nullptr, MeshPoint2d** _points = nullptr);
	/// Standard destructor
	virtual ~MeshElement() {}
public:
	/// Set tag for all adjacent edges
	void setTagForEdges(TagExtended::TagType type, int t = 1) const;
	/// Returns true if the sequence of the two given points complies to the orientation of element
	bool properOrientation(const MeshPoint2d *mp1, const MeshPoint2d*mp2) const;
	/// Returns i-th vertex of this element
	MeshPoint2d* getPoint(int i) const { return points[i]; }
	/// Get real angles (projected onto plane)
	DataVector<double> getRealAngles(std::shared_ptr<const SurfaceParametric> surface) const;
	/// Return representative points (using param_point_distance)
	void getRepresentativePoints(DataVector<DPoint2d> & pts) const;
	/// Return distance (squared) point-triangle (using param_point_distance)
	double distance2(const DPoint2d& pt) const;
	/// Check whether given point is too close to one of vertices (using metric space)
	MeshPoint2d* closeToVertex(Metric2dContext& mc, MeshPoint2d* point) const;
	/// Return other edge adjacent to this point
	MeshEdge2d* getOtherEdge(const MeshEdge2d* edge, const MeshPoint2d* point) const;
	/// Returns number of edges (e.g. 3 for triangles)
	int getEdgeCount() const { return count; }
	/// Returns the reference to the element incident through the i-th edge (can be nullptr)
	MeshElement* getNeighbour(int i) const;
	/// Returns the reference to the i-th edge
	MeshEdge2d* getEdge(int i) const { return edges[i]; }
	/// Returns the reference to the edge after the given one (anticlockwise)
	MeshEdge2d* getNextEdge(const MeshEdge2d* edge) const;
	/// Returns the reference to the edge before the given one (anticlockwise)
	MeshEdge2d* getPrevEdge(const MeshEdge2d* edge) const;
	/// Returns the difference between the minimum and maximum ID of vertices of this element
	int getIdSpan() const;
	/// Returns the minimum area for this element (in unitary space)
	double getMinArea() const { return 0.25; }
	/// Returns one of vertices of this element with index > k and the least rank
	MeshPoint2d* getMinVertexBiggerThen(int k) const;
	/// Returns the vertex on the edge incident to the given point, with index > k (if exists)
//	MeshPoint2d* getInnerVertexBiggerThen(const MeshPoint2d* point, int k) const;
public:
	/// clear some data for faster all-delete process
	virtual void preDeleteAll() {}
	/// add this element as a lump of mass (respecting control space)
	virtual void addForInertialCenter(Metric2dContext& mc, DPoint2d& inertial_center, double& total_mass) const;
	/// add this element as a lump of mass (respecting control space)
	virtual void addForInertialMoments(Metric2dContext& mc, const DPoint2d& inertial_center, DMatrix2d& inertial_moments) const;
	/// Sets and returns the quality of this element (with globally or explicitely selected criterion)
	virtual double countQuality(Metric2dContext& mc, bool mind_area = false, int criterion = -1) = 0;
	/// Sets and returns the metric difference of this element
	virtual double countMetricDiffQuality(Metric2dContext& /*mc*/) { return quality = 1.0; }
	/// Returns "mean ratio" coefficient quality
	virtual double getMeanRatio(Metric2dContext& /*mc*/, bool /*ext_metric*/) const { assert(false); return 1.0; }
	/// Returns the metric difference of this element
	virtual double countMetricDiff(Metric2dContext& /*mc*/) const { return LARGE_NUMBER; }
	/// Returns the i-th inner angle in this element (can be calculated in modified metric)
	virtual double getAngle(Metric2dContext& mc, int i, bool local_metric = true) const = 0;
	/// Returns the alpha-quality for this element (can be calculated in modified metric)
	virtual double getAlphaQuality(Metric2dContext& mc, bool local_metric = true) const = 0;
	/// Returns the size of this element (can be calculated in modified metric)
	virtual double getArea(Metric2dContext& mc, bool local_metric = true) const = 0;
	/// Returns the size of this element (with non-altered coordinates)
	virtual double getAreaNoMetric() const = 0;
	/// Checks whether the element is filled (for regular mesh elements it is alwyas the case)
	virtual bool isFilled() const { return true; }
	/// Stores the description of this elements into the given file (in text GRD format)
	virtual void exportToGRD(FILE* file) const = 0;
	/// Checks whether the given point is within the element
	virtual bool isPointInside(const DPoint2d& mpt) const = 0;
	/// Checks whether the given point is within the element
	bool isPointInside(MeshPoint2d *point) const;
	/// Checks whether the element is inverted (i.e. bad orientation or invalid geometry)
	virtual bool isInverted() const = 0;
	/// Returns the description of this point for visualization
	std::shared_ptr<MeshViewFaceData> getViewData(double shrink_ratio = 1.0, 
		std::shared_ptr<const SurfaceParametric> surface = nullptr, bool proper_orientation = true) const;
	/// Returns the type of object 
	virtual ElementType getType() const { return ELEMENT_UNKNOWN; }
	/// Compares this element with the given one (required if the heap-ordering of container is used)
	short compareTo(MeshElement* item);
	/// Returns the current quality of element (bufferred or recalculated)
	double getQuality(){ assert(quality >= 0.0); return quality; }
	/// Sets the quality value for this element (buffered for efficiency)
	void setQuality(double q){ quality = q; }
	/// Returns the area-ID of element
	int getAreaID() const { return area_id; }
	/// Sets the area-ID of element
	void setAreaID(int id){ area_id = id; }
	/// Returns the vertex before the given one (anticlockwise)
	MeshPoint2d* getPrevPoint(const MeshPoint2d* point) const;
	/// Returns the vertex after the given one (anticlockwise)
	MeshPoint2d* getNextPoint(const MeshPoint2d* point) const;
	/// Switches one of vertices to other mesh point (with adequate update of edges)
	void switchPointsWithEdges(const MeshPoint2d* point1, MeshPoint2d* point2);
	/// Returns the index of the given point as a vertex of this element
	int getPointIndex(const MeshPoint2d* point) const;
	/// Switches one of vertices to other mesh point 
	void switchPoints(const MeshPoint2d* point1, MeshPoint2d* point2);
	/// Switches one of edges to other mesh point 
	void switchEdges(const MeshEdge2d* edge1, MeshEdge2d* edge2);
	/// Returns the middle point for this element
	virtual DPoint2d getMiddlePoint() const;
protected:
	/// Buffered quality coefficient (invalid if < 0.0)
	double	quality;
	/// Area identifier for this element
	int		area_id;
	/// Count of edges
	int			count;
	/// Array of references to edges
	MeshEdge2d**	edges;
	/// Array of references to vertices
	MeshPoint2d**	points;
public:
	static double last_area;
	static int param_point_distance;
};

#endif // !defined(MESHELEMENT_H__INCLUDED)
