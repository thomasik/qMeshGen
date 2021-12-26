/////////////////////////////////////////////////////////////////////////////
// DOrientedBox.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2014-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DORIENTEDBOX_H__INCLUDED)
#define DORIENTEDBOX_H__INCLUDED

#include "DPoint.h"
#include "DMatrix.h"

class Metric3dContext;
class MeshViewSet;

/**
 * This class implements an oriented cubicoid (mostly used as a bounding-box)
 *  plus some basic operations
 */
class DOrientedBox  
{
public:
	/// Standard constructor (creates empty - invalid box)
	DOrientedBox();
	/// Standard constructor (creates empty - invalid box)
	DOrientedBox(const DPoint3d& pt0, const DVector3d& e0, const DVector3d& e1);
	/// Standard constructor (creates valid box)
	DOrientedBox(double _x0, double _x1, double _y0, double _y1, double _z0, double _z1);
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
	/// Returns the volume of this box
	double getVolume() const { return (x1-x0)*(y1-y0)*(z1-z0)*m_rot.det(); }
	/// Increases (if needed) the box as a union of this and the given box
	void addOrientedBox(const DOrientedBox& box);
	/// Increases (if needed) the box as a union of this box and the given point
	void addPoint(const DPoint3d& point);
	/// Increases (if needed) the box as a union of this box and the given point
	void addPoint(const FPoint3d& point);
	/// Increases (if needed) the box as a union of this box and the given point
	void addPoint(const DMPoint3d& point);
	/// Checks, whether this rectangle contains the given 2D-point
	bool contains(const DPoint3d& p) const;
	/// Returns middle 3D-point of the box
	DPoint3d getMiddlePoint() const;
	/// Returns corner of this box
	//DPoint3d getMinPoint() const { return DPoint3d(x0, y0, z0); }
	/// Returns corner of this box
	//DPoint3d getMaxPoint() const { return DPoint3d(x1, y1, z1); }
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
	/// Stores the coordinates of the box into the stream
	friend ostream& operator<<(ostream& os, const DOrientedBox& box);
	/// Loads the coordinates of the box from the stream
	friend istream& operator>>(istream& is, DOrientedBox& box);
	/// Draw edges of this oriented box
	MeshViewSet* draw( MeshViewSet* set, int id = -2) const; 
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
	/// Marks validity of this box (whether it was already initialised)
	bool valid;
	/* Orientation data */
	// Translation vector
	DVector3d m_trans;
	// Rotation matrix
	DMatrix3d m_rot, m_rot_inv;
};

#endif // !defined(DORIENTEDBOX_H__INCLUDED)
