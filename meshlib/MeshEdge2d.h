/////////////////////////////////////////////////////////////////////////////
// MeshEdge2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef	MESHEDGE_H_INCLUDED
#define MESHEDGE_H_INCLUDED

#include <memory>

//#include "common.h"
//#include "MeshData.h"
#include "DPoint.h"
//#include "Curve2dParametric.h"
#include "TagBorder.h"
#include "TagExtended.h"
#include "DataVector.h"

class MeshPoint2d;
class MeshElement;
struct MeshViewEdgeData;
class Metric2dContext;
class SurfaceParametric;
class DRect;
class ControlSpace2d;
class ControlSpace2dAdaptive;

/**
 * This class implements a mesh edge in 2D and several operations.
 */
class MeshEdge2d : public TagBorder, public TagExtended
{
public:
	/// Orientation of adjacent area (for directed edges)
	enum AreaOrientation {AREA_LEFT = 0, AREA_RIGHT = 1};
	/// Auxiliary marking-type (mostly for running front)
	enum EnumTagType {TAG_LEFT = 1, TAG_RIGHT = 2, TAG_UNCHECKED = 4};
public:
	/// Standard constructor
	MeshEdge2d(MeshPoint2d *p1, MeshPoint2d *p2, MeshElement* e1 = nullptr, MeshElement* e2 = nullptr);
	/// Standard destructor
	virtual ~MeshEdge2d();
public:
	/// Returns the parameter ksi for a given point (numerical approximation)
	virtual double getParameter(const DPoint2d& pt) const;
	/// Replaces one element link with another
	bool switchElementLink(const MeshElement* old_element, MeshElement* new_element);
	/// Removes link to this point (fo delete-all phase) - returns true if last one
	bool removePointLink(MeshPoint2d* point);
	/// Returns curvature of edge shape on plane
	virtual double getPlanarCurvature(double /* ksi */) const { return 0.0; }
	/// Returns curvature of edge shape on non-planar surface
	virtual double getNonPlanarCurvature(std::shared_ptr<const SurfaceParametric> surface, double ksi, double* cdt_len = nullptr) const;
	/// Update given ACS for curvature of edge contour
	virtual bool updateACSwithCurvature(std::shared_ptr<ControlSpace2d> /* space */, 
		std::shared_ptr<const SurfaceParametric> /* surface */, double /* d2_threshold */ ) const {
		return false; }
	/// Create new edge with similar geometry
	virtual MeshEdge2d* cloneGeometric(MeshPoint2d *p0, MeshPoint2d *p1, double ksi0 = 0.0, double ksi1 = 1.0, MeshElement* e1 = nullptr, MeshElement* e2 = nullptr) const;
	/// numerical approximation -> calculate parameter for point on surface+curve
	virtual DPoint2d surfaceParameters(std::shared_ptr<const SurfaceParametric> surface, const DPoint3d& point, double & ksi, bool same_direction) const;
	/// Returns the type of this object
	virtual ElementType getType() const { return EDGE_SIMPLE; }
	/// Includes this edge into bounding rectangle
	virtual void addToBoundingRect(DRect& rect) const;
	/// Checks whether this edge is incident to some point
	bool incidentTo(const MeshPoint2d* point) const { return points[0] == point || points[1] == point; }
	/// Returns the description of this point for visualization
	std::shared_ptr<MeshViewEdgeData> getViewData(std::shared_ptr<const SurfaceParametric> surface = nullptr) const;
	/// Returns the coordinates of point within this edge, described by parameter ksi [0,1]
	virtual DPoint2d getPoint(double ksi) const;
	/// Returns the length of a subsegment of this edge (in 3D space)
	double getLength(std::shared_ptr<const SurfaceParametric> surface, double ksi0 = 0.0, double ksi1 = 1.0) const;
	/// Returns the length of a subsegment of this edge (current metric can be used)
	virtual double getLength(Metric2dContext& mc, double ksi0, double ksi1, bool local_metric = true) const;
	/// Returns the length of a subsegment of this edge (current metric can be used)
	virtual double getLengthMax(Metric2dContext& mc, double ksi0, double ksi1) const;
	/// Returns the length of a subsegment of this edge (current metric can be used)
	virtual double checkAndGetLength(Metric2dContext& mc, double ksi0, double& ksi1, double max_len = 1.0, bool local_metric = true) const;
	/// Returns the length of this edge (current metric can be used)
	virtual double getLength(Metric2dContext& mc, bool local_metric = true) const;
	/// Returns the metric length quality of this edge
	double getLengthQuality(Metric2dContext& mc, bool ext_metric) const;
	/// Returns the length of this edge (for variable metric)
	double getRequiredLength(Metric2dContext& mc, double ksi0, double &ksi1, double req_length) const;
	/// Returns the length of this edge (for variable metric)
	double getLengthMetricAdapted(Metric2dContext& mc, double ksi0, double ksi1, int level = 0) const;
	/// Returns the length of this edge (for variable metric)
	double getLengthMetricAdapted(Metric2dContext& mc) const;
	/// Returns the weight based on the type of adjacent elements (i.e. triangles or quads)
	double getElementsWeight() const;
	/// Returns the element (or nullptr) adjacent to this edge from the left side (considering the given vertex as a beginning of edge)
	MeshElement* getMeshElement(const MeshPoint2d* point) const;
	/// Switches the given vertex in the edge with the other one
	void switchPoints(const MeshPoint2d* point1, MeshPoint2d* point2);

	/// Checks wheter the specified side of edge is tagged (given point denotes the beginning of the edge)
	bool isSideTagged(const MeshPoint2d *point) const { 
		return hasIntFlag(TagExtended::TAG_SIDE_EDGE, (point == points[0]) ? TAG_LEFT : TAG_RIGHT); }
	/// Checks wheter the edge is tagged (on either side)
	bool isSideTagged() const { 
		return hasIntFlag(TagExtended::TAG_SIDE_EDGE, TAG_LEFT | TAG_RIGHT); }
	/// Marks (sets tag on) the specified side of this edge (given point denotes the beginning of the edge)
	void setSideTag(const MeshPoint2d *point) { 
		setIntFlag(TagExtended::TAG_SIDE_EDGE, (point == points[0]) ? TAG_LEFT : TAG_RIGHT); }
	/// Marks this edge (on both sides)
	void setSideTag(){ 
		setIntFlag(TagExtended::TAG_SIDE_EDGE, TAG_LEFT | TAG_RIGHT); }
	/// Unmarks (clears tag on) the specified side of this edge (given point denotes the beginning of the edge)
	void clearSideTag(const MeshPoint2d *point){ 
		clearIntFlag(TagExtended::TAG_SIDE_EDGE, (point == points[0]) ? TAG_LEFT : TAG_RIGHT); }
	/// Unmarks (clears tag on) the side of this edge specified by index
	void clearSideTag(int index){ 
		clearIntFlag(TagExtended::TAG_SIDE_EDGE, (index == 0) ? TAG_LEFT : TAG_RIGHT); }
	/// Clears marks on this edge (on both sides)
	void clearSideTag(){ 
		clearIntFlag(TagExtended::TAG_SIDE_EDGE, TAG_LEFT | TAG_RIGHT); }
	/// Sets additional "check" flag for this edge (side is unimportant)
	void setSideUnchecked(){ 
		setIntFlag(TagExtended::TAG_SIDE_EDGE, TAG_UNCHECKED); }
	/// Clears additional "check" flag for this edge (side is unimportant)
	void clearSideUnchecked(){ 
		clearIntFlag(TagExtended::TAG_SIDE_EDGE, TAG_UNCHECKED); }
	/// Checks wheter the additional "check" flag for this edge is set (side is unimportant)
	bool isSideUnchecked() const { 
		return hasIntFlag(TagExtended::TAG_SIDE_EDGE, TAG_UNCHECKED); }
	/// Returns the approximation of this edge via polyline (array of points)
	virtual void getPolyLine(DataVector<DPoint2d> & polyline) const;
	/// Returns the approximation of this edge via polyline (array of points)
	virtual void getPolyLine(DataVector<DPoint3d> & polyline, 
		std::shared_ptr<const SurfaceParametric> surface) const;
	/// Returns the index (0 or 1) of this point as a vertex in this directed edge
	int getPointIndex(const MeshPoint2d* point) const;
	/// Sets the area-identifier of the adjacent element at the given side
	void setIncidentAreaID(int side, int id);
	/// Returns the direction of this edge (relative to the given vertex)
	int getDirection(const MeshPoint2d* point) const;
	/// Checks wheter any element is adjacent to this edge
	bool isBounded() const { return elements[0] || elements[1]; }
	/// Returns the area-identifier of the element at the given side
	int getIncidentAreaID(int side) const;
	/// Changes the direction of the edge
	void addDirection(int dir, const MeshPoint2d* point);
	/// Sets the direction of the edge
	void setDirection(int dir, const MeshPoint2d* point);
	/// Checks whether the edge is directed
	bool isDirected() const { return isBorder(TagBorder::OUTER); }
	/// Adds the reference to the adjacent element at the given side (given point denotes the beginning of the edge)
	void addElementLink(MeshElement* element, const MeshPoint2d *p1);
	/// Removes the reference to the given element from the list of adjacent elements of this edge
	bool removeElementLink(const MeshElement* element);
	/// Returns on of the adjacent edge vertics
	MeshPoint2d* getMeshPoint(int i) const { return points[i]; }
	/// Returns one of the adjacent elements
	MeshElement* getMeshElement(int i) const { return elements[i]; }
	/// Returns the other vertex of this edge
	MeshPoint2d* getOtherPoint(const MeshPoint2d* pt) const { assert(points[0] == pt || points[1] == pt); return (points[0]==pt)?points[1]:points[0]; }
	/// Returns the other adjacent element of this edge
	MeshElement* getOtherElement(const MeshElement* el) const { return (elements[0]==el)?elements[1]:elements[0]; }
protected:
	/// Reverses the directionality of the edge
	virtual void switchSide();
protected:
	/// Vertices of this edge (array of references to mesh points)
	MeshPoint2d*	points[2];
	/// Adjacent elements (array of references to mesh elements)
	MeshElement*	elements[2];
};

#endif // MESHEDGE_H_INCLUDED
