 /////////////////////////////////////////////////////////////////////////////
// DRect.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DRECT_H__INCLUDED)
#define DRECT_H__INCLUDED

#include "DPoint.h"

class Metric3dContext;
class MeshViewSet;

/**
 * This class implements a rectangle in a two-dimensional space
 *  plus some basic operations
 */
class DRect  
{
public:
	/// Standard constructor (creates empty - invalid rectangle)
	DRect() : x0(0.0), x1(0.0), y0(0.0), y1(0.0), valid(false) {}
	/// Standard constructor (creates valid rectangle)
	DRect(double _x0, double _x1, double _y0, double _y1) 
		: x0(_x0), x1(_x1), y0(_y0), y1(_y1), valid(true) {}
	/// Sets all coordinates (and makes it valid)
	void set(double _x0, double _x1, double _y0, double _y1) {
		x0 = _x0; x1 = _x1; y0 = _y0; y1 = _y1; valid = true;
	}
	/// Sets X coordinate of the x0 edge of the rectangle (and makes it valid)
	void setX0(double _x0){ x0 = _x0; valid = true; }
	/// Sets Y coordinate of the y1 edge of the rectangle (and makes it valid)
	void setY1(double _y1){ y1 = _y1; valid = true; }
	/// Sets X coordinate of the x1 edge of the rectangle (and makes it valid)
	void setX1(double _x1){ x1 = _x1; valid = true; }
	/// Sets Y coordinate of the y0 edge of the rectangle (and makes it valid)
	void setY0(double _y0){ y0 = _y0; valid = true; }
	/// Returns x0 upper vertex of the rectangle as a 2D-point
	DPoint2d getX0Y1() const { return DPoint2d(x0, y1); }
	/// Returns x0 y0 vertex of the rectangle as a 2D-point
	DPoint2d getX0Y0() const { return DPoint2d(x0, y0); }
	/// Returns x1 upper vertex of the rectangle as a 2D-point
	DPoint2d getX1Y1() const { return DPoint2d(x1, y1); }
	/// Returns x1 y0 vertex of the rectangle as a 2D-point
	DPoint2d getX1Y0() const { return DPoint2d(x1, y0); }
	/// Returns middle 2D-point of the rectangle
	DPoint2d getMiddlePoint() const { return DPoint2d(0.5*(x0+x1), 0.5*(y0+y1)); }
	/// Returns the diameter of this rectangle
	double getDiameter() const { return sqrt(sqr(x1-x0)+sqr(y1-y0)); }
	/// Returns the squared diameter of this rectangle
	double getDiameter2() const { return sqr(x1-x0)+sqr(y1-y0); }
	/// Returns the width of the rectangle
	double getDX() const { return x1 - x0; }
	/// Returns the height of the rectangle
	double getDY() const { return y1 - y0; }
	/// Returns the point prescribed by parameters (s,t) with range 0..1
	DPoint2d getPoint(double s, double t) const { return DPoint2d(x0+(x1-x0)*s, y0+(y1-y0)*t); }
	/// Increases (if needed) the rectangle as a union of this and the given rectangle
	void addRect(const DRect& rct);
	/// Increases (if needed) the rectangle as a union of this and the given point
	void addPoint(const DPoint2d& point);
	/// Checks, whether this rectangle contains the given 2D-point
	bool contains(const DPoint2d& p) const { return (p.x >= x0) && (p.x <= x1) && (p.y >= y0) && (p.y <= y1); }
	/// Increases size of the rectangle by the given factor (equally in all directions)
	void inflate(double factor);
	/// Returns rectangle with dncreased size by the given factor (equally in all directions)
	DRect inflated(double factor) const;
	/// Transpose rectangle by given vector
	void transpose(const DVector2d& v) { x0 += v.x; x1 += v.x; y0 += v.y; y1 += v.y; }
	Axis getLongestAxis() const;
	double getMiddleValue(Axis axis) const;
	void split(Axis axis, double value, DRect& box0, DRect& box1) const;
	void splitAndUpdate(Axis axis, double value, bool lower);
	/// Adjusts the cooridnates of the given point, so it shoud be in this rect
	DPoint2d fitInPoint(const DPoint2d& pt) const;
	/// Clears the values of rectangle (and makes it invalid)
	void clear() {x0 = x1 = y1 = y0 = 0.0; valid = false; }
	/// Calculates two crossing points for rectangle and line 
	bool calculateCrossingPoints(const DPoint2d& line_pt, const DVector2d& line_vt, DPoint2d& pt0, DPoint2d& pt1) const;
	/// ostream
	friend ostream& operator<<(ostream& os, const DRect& box){
		return os << '[' << box.x0 << ',' << box.y0 << "]x[" << box.x1 << ',' << box.y1 << ']';
	}

	/// Returns the width of the box
	double getDLen(Axis axis) const;
	double getD0(Axis axis) const;
	double getD1(Axis axis) const;
public:
	/// X-coordinate of the x0 edge
	double x0;
	/// X-coordinate of the x1 edge
	double x1;
	/// Y-coordinate of y0 edge
	double y0;
	/// Y-coordinate of the y1 edge
	double y1;
	/// Marks validity of this rectangle (whether it was already initiated)
	bool valid;
};

/**
 * This class implements a cubicoid (mostly used as a bounding-box)
 *  plus some basic operations
 */
class DBox  
{
public:
	/// Standard constructor (creates empty - invalid box)
	DBox() : x0(0.0),x1(0.0),y0(0.0),y1(0.0),z0(0.0),z1(0.0),valid(false) {}
	/// Standard constructor (creates valid box)
	DBox(double _x0, double _x1, double _y0, double _y1, double _z0, double _z1) 
		: x0(_x0), x1(_x1),y0(_y0),y1(_y1),z0(_z0),z1(_z1),valid(true) {}
	/// Sets all coordinates (and makes it valid)
	void set(double _x0, double _x1, double _y0, double _y1, double _z0, double _z1) {
		x0 = _x0; x1 = _x1; y0 = _y0; y1 = _y1; z0 = _z0; z1 = _z1; valid = true; }
	/// Sets X coordinate of the x0 face of the box (and makes it valid)
	void setX0(double _x0){ x0 = _x0; valid = true; }
	/// Sets X coordinate of the x1 face of the box (and makes it valid)
	void setX1(double _x1){ x1 = _x1; valid = true; }
	/// Sets Y coordinate of the near face of the box (and makes it valid)
	void setY0(double _y0){ y0 = _y0; valid = true; }
	/// Sets Y coordinate of the far face of the box (and makes it valid)
	void setY1(double _y1){ y1 = _y1; valid = true; }
	/// Sets Z coordinate of the z0 face of the box (and makes it valid)
	void setZ0(double _z0){ z0 = _z0; valid = true; }
	/// Sets Z coordinate of the upper face of the box (and makes it valid)
	void setZ1(double _z1){ z1 = _z1; valid = true; }
	/// Returns the width of the box
	double getDX() const { return x1 - x0; }
	/// Returns the depth of the box
	double getDY() const { return y1 - y0; }
	/// Returns the height of the box
	double getDZ() const { return z1 - z0; }
	/// Returns the diameter of this box
	double getDiameter() const { return sqrt(sqr(x1-x0)+sqr(z1-z0)+sqr(y1-y0)); }
	/// Returns the diameter of this box
	DVector3d getDiameterVec() const { return DVector3d(x1 - x0, y1 - y0, z1 - z0); }
	/// Returns the volume of this box
	double getVolume() const { return (x1-x0)*(y1-y0)*(z1-z0); }
	/// Increases (if needed) the box as a union of this and the given box
	void addBox(const DBox& box);
	/// Increases (if needed) the box as a union of this box and the given point
	void addPoint(const DPoint3d& point);
	/// Increases (if needed) the box as a union of this box and the given point
	void addPoint(const FPoint3d& point);
	/// Increases (if needed) the box as a union of this box and the given point
	void addPoint(const DMPoint3d& point);
	/// Checks, whether this rectangle contains the given 2D-point
	bool contains(const DPoint3d& p) const { return (p.x >= x0) && (p.x <= x1) && 
		(p.y >= y0) && (p.y <= y1) && (p.z >= z0) && (p.z <= z1); }
	/// Checks, whether this rectangle contains the given 2D-point
	bool containsEps(const DPoint3d& p, double eps = 1e-5) const {
		double epsx = (x1 - x0)*eps;
		double epsy = (y1 - y0)*eps;
		double epsz = (z1 - z0)*eps;
		return (p.x - x0 > -epsx) && (x1 - p.x > -epsx) &&
			(p.y - y0 > -epsy) && (y1 - p.y > -epsy) && 
			(p.z - z0 > -epsz) && (z1 - p.z > -epsz);
	}
	/// Returns middle 3D-point of the box
	DPoint3d getMiddlePoint() const { return DPoint3d(0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1)); }
	/// Returns corner of this box
	DPoint3d getMinPoint() const { return DPoint3d(x0, y0, z0); }
	/// Returns corner of this box
	DPoint3d getMaxPoint() const { return DPoint3d(x1, y1, z1); }
	/// Increases size of the box by the given factor (proportionally in all directions)
	void inflate(double factor);
	/// Increases size of the box by the given length (equally in all directions)
	void grow(double offset);
	/// Increases size of the box by the given length (in metric, in all directions)
	void growDM(Metric3dContext& mc, double mlen);
	/// Increases X-size of the box by the given length
	void growDX(double dx) { x0 -= dx; x1 += dx; }
	/// Increases Y-size of the box by the given length
	void growDY(double dy) { y0 -= dy; y1 += dy; }
	/// Increases Z-size of the box by the given length
	void growDZ(double dz) { z0 -= dz; z1 += dz; }
	/// Clears the values of box (and makes it invalid)
	void clear() {x0 = x1 = y0 = y1 = z0 = z1 = 0.0; valid = false; }
	/// Adjusts the cooridnates of the given point, so it shoud be in this box
	DPoint3d fitInPoint(const DPoint3d& pt) const;
	DPoint3d fitInPoint(const DPoint3d& pt, const DVector3d& dv) const;
	/// Draw edges of this box
	MeshViewSet* draw( MeshViewSet* set, int id = -2) const; 
	MeshViewSet* drawXYZ(MeshViewSet* set, int id = -2) const;
	/// Stores the coordinates of the box into the stream
	friend ostream& operator<<(ostream& os, const DBox& box);
	/// Loads the coordinates of the box from the stream
	friend istream& operator>>(istream& is, DBox& box);

	Axis getLongestAxis() const;
	double getMiddleValue(Axis axis) const;
	DPoint3d getRandomPoint() const;
	DPoint3d getRandomPoint(double margin) const;
	DPoint3d getRandomPoint(double margin_x, double margin_y, double margin_z) const;

	void split(Axis axis, double value, DBox& box0, DBox& box1) const;
	void splitAndUpdate(Axis axis, double value, bool lower);
	DBox splitOct(const DPoint3d& split_pt, int v) const;
	void splitOctAndUpdate(const DPoint3d& split_pt, int v);
	DPoint3d getVertex(int v) const;
	DPoint3d getMiddlePointForFace(int fi) const;
	DPoint3d getMiddlePointForEdge(int ei) const;

	/// Returns the width of the box
	double getDLen(Axis axis) const;
	double getD0(Axis axis) const;
	double getD1(Axis axis) const;

	void forRegularGrid(int grid_size, const std::function<void(const DPoint3d&pt)>& fg) const;
public:
	//-- static
	static OctVertexWhich edge_to_vertex[12][2];
	static OctFaceWhich edge_to_face[12][2];
	static OctEdgeWhich vertex_to_edge[8][3];
	static OctFaceWhich vertex_to_face[8][3];
	static OctVertexWhich face_to_vertex[6][4]; // for each side of leaf -> 4 sub-leaves
	static OctEdgeWhich face_to_edge[6][4];
	static OctVertexWhich opposite_vertex[8];
public:
	/// X-coordinate of the x0 face
	double x0;
	/// X-coordinate of the x1 face
	double x1;
	/// Y-coordinate of the near face
	double y0;
	/// Y-coordinate of the far face
	double y1;
	/// Z-coordinate of the z0 face
	double z0;
	/// Z-coordinate of the upper face
	double z1;
	/// Marks validity of this box (whether it was already initiated)
	bool valid;
};

#endif // !defined(DRECT_H__INCLUDED)
