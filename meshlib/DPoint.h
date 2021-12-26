/////////////////////////////////////////////////////////////////////////////
// DPoint.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DPOINT_H__INCLUDED
#define DPOINT_H__INCLUDED

#include "DVector.h"
#include "DataVector.h"

/**
 * This class implements a point in a two-dimensional space
 *  plus some basic operations.
 */
class DPoint2d
{
public:
	/// Standard constructor
	DPoint2d(double _x, double _y) : x(_x), y(_y) {}
	/// linear
	DPoint2d(const DPoint2d& a, const DPoint2d& b, double t = 0.5) 
		: x(a.x * (1.0-t) + b.x * t), y(a.y * (1.0-t) + b.y * t) {} 
	/// average
	DPoint2d(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c) 
		: x( (a.x  + b.x  + c.x) * (1.0/3.0) ), y ((a.y + b.y  + c.y) * (1.0/3.0)) {} 
	/// Standard (empty) constructor
	DPoint2d() : x(0.0), y(0.0) {}
public:
	/// Translates point by vector
	const DPoint2d	operator+(const DVector2d& v) const { return DPoint2d(x+v.x, y+v.y); }
	/// Translates point by vector
	const DPoint2d	operator-(const DVector2d& v) const { return DPoint2d(x-v.x, y-v.y); }
	/// Translates the point by the given vector
	DPoint2d&	operator+=(const DVector2d& v){ x+=v.x; y+=v.y; return *this; }
	/// Translates the point by the given vector
	DPoint2d&	operator-=(const DVector2d& v){ x-=v.x; y-=v.y; return *this; }
	/// Returns the difference of two vectors
	const DVector2d	operator-(const DPoint2d& pt) const { return DVector2d(x-pt.x, y-pt.y); }
	/// Transforms to fixed vector (from point [0,0])
	const DVector2d   fixedVector() const { return DVector2d(x, y); }

	/// Divide point by number
	DPoint2d&	operator/=(double d){ x/=d; y/=d; return *this; }
	/// Returns the point scaled by the given factor
	DPoint2d	operator*(double d) const { return DPoint2d(d*x, d*y); }
	/// Multiply point by number
	DPoint2d&	operator*=(double d){ x*=d; y*=d; return *this; }
	/// Set to coordinates
	DPoint2d& set(const DPoint2d& pt, double f = 1.0) {
		x = pt.x * f; y = pt.y * f; return *this; 
	}
	/// Add to coordinates
	DPoint2d& add(const DPoint2d& pt, double f) {
		x += pt.x * f; y += pt.y * f; return *this; 
	}
	/// Add to coordinates
	DPoint2d& add(const DPoint2d& pt) {
		x += pt.x; y += pt.y; return *this; 
	}

	// "less-than" operator, for convex-hull algorithm
    bool operator<(const DPoint2d &pt) const {
            return x < pt.x || (x == pt.x && y < pt.y);
    }

	double& operator[](Axis axis) { return *(&x + static_cast<int>(axis)); }
	const double& operator[](Axis axis) const { return *(&x + static_cast<int>(axis)); }

	/// Returns the distance between two points
	double distance(const DPoint2d& a) const { return sqrt(sqr(a.x - x) + sqr(a.y - y)); }
	/// Returns the squared distance between two points
	double distance2(const DPoint2d& a) const { return sqr(a.x - x) + sqr(a.y - y); }

	/// Average of two points
	static const DPoint2d average(const DPoint2d& a, const DPoint2d& b) {
		return DPoint2d((a.x+b.x)*0.5, (a.y+b.y)*0.5);
	}
	/// Average of three points
	static const DPoint2d average(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c) {
		return DPoint2d((a.x+b.x+c.x)*(1.0/3.0), (a.y+b.y+c.y)*(1.0/3.0));
	}
	/// Average of four points
	static const DPoint2d average(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c, const DPoint2d& d) {
		return DPoint2d((a.x+b.x+c.x+d.x)*0.25, (a.y+b.y+c.y+d.y)*0.25);
	}

	/// orientation check
	static bool properOrientation(const DataVector<DPoint2d> & points);
	/// clear edges with "zero" length (less than threshold)
	static int clearZeroEdges(DataVector<DPoint2d> & points, double threshold = VERY_SMALL_NUMBER);
	/// clear edges with "zero" length (less than threshold)
	static int clearZeroEdges(DataVector<DPoint2d> & points, DataVector<int> &trans_tab, double threshold = VERY_SMALL_NUMBER);

	bool storeSimple(ostream& os, char end_char = '\n') const;
	bool readSimple(istream& is);

	/// Stores the coordinates of the point into the stream
	friend ostream& operator<<(ostream& os, const DPoint2d& pt);
	/// Loads the coordinates of the point from the stream
	friend istream& operator>>(istream& is, DPoint2d& pt);

public:
	/// x-coordinate of the point
	double x; 
	/// y-coordinate of the point
	double y; 
	// zero point
	static const DPoint2d zero;
	// inf point
	static const DPoint2d inf;
};

/**
 * This class implements a point in a two-dimensional metric space
 *  plus some basic operations.
 */
class DMPoint2d
{
public:
	/// Standard constructor
	DMPoint2d(double _x, double _y) : x(_x), y(_y) {}
	/// linear
	DMPoint2d(const DMPoint2d& a, const DMPoint2d& b, double t) 
		: x(a.x * (1.0-t) + b.x * t), y(a.y * (1.0-t) + b.y * t) {} 
	/// Standard (empty) constructor
	DMPoint2d() : x(0.0), y(0.0) {}
public:
	/// Translates point by vector
	const DMPoint2d	operator+(const DMVector2d& v) const { return DMPoint2d(x+v.x, y+v.y); }
	/// Translates point by vector
	const DMPoint2d	operator-(const DMVector2d& v) const { return DMPoint2d(x-v.x, y-v.y); }
	/// Translates the point by the given vector
	DMPoint2d&	operator+=(const DMVector2d& v){ x+=v.x; y+=v.y; return *this; }
	/// Translates the point by the given vector
	DMPoint2d&	operator-=(const DMVector2d& v){ x-=v.x; y-=v.y; return *this; }
	/// Returns the difference of two vectors
	const DMVector2d	operator-(const DMPoint2d& pt) const { return DMVector2d(x-pt.x, y-pt.y); }
	/// Transforms to fixed vector (from point [0,0])
	const DMVector2d   fixedVector() const { return DMVector2d(x, y); }

	/// Divide point by number
	DMPoint2d&	operator/=(double d){ x/=d; y/=d; return *this; }
	/// Multiply point by number
	DMPoint2d&	operator*=(double d){ x*=d; y*=d; return *this; }
	/// Add to coordinates
	DMPoint2d& add(const DMPoint2d& pt, double f) {
		x += pt.x * f; y += pt.y * f; return *this; 
	}

	/// Returns the distance between two points
	double distance(const DMPoint2d& a) const { return sqrt(sqr(a.x - x) + sqr(a.y - y)); }
	/// Returns the squared distance between two points
	double distance2(const DMPoint2d& a) const { return sqr(a.x - x) + sqr(a.y - y); }

	/// Stores the coordinates of the point into the stream
	friend ostream& operator<<(ostream& os, const DMPoint2d& pt);
	/// Loads the coordinates of the point from the stream
	friend istream& operator>>(istream& is, DMPoint2d& pt);

public:
	/// x-coordinate of the point
	double x; 
	/// y-coordinate of the point
	double y; 
	// zero point
	static const DMPoint2d zero;
};

/**
 * This class implements a point in a three-dimensional space
 *  plus some basic operations
 */
class DPoint3d
{
public:
	/// Standard constructor
	DPoint3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	/// linear
	DPoint3d(const DPoint3d& a, const DPoint3d& b, double t) 
		: x(a.x * (1.0-t) + b.x * t), y(a.y * (1.0-t) + b.y * t), z(a.z * (1.0-t) + b.z * t) {} 
	/// Standard (empty) constructor
	DPoint3d() : x(0.0), y(0.0), z(0.0) {}
	/// Standard constructor
	DPoint3d(const DPoint2d& point, double _z) : x(point.x), y(point.y), z(_z) {}
public:
	/// Translates the point by the given vector
	const DPoint3d	operator+(const DVector3d& v) const { return DPoint3d(x+v.x, y+v.y, z+v.z); }
	/// Translates the point by the given vector
	const DPoint3d	operator-(const DVector3d& v) const { return DPoint3d(x-v.x, y-v.y, z-v.z); }
	/// Translates the point by the given vector
	DPoint3d&	operator+=(const DVector3d& v){ x+=v.x; y+=v.y; z+=v.z; return *this; }
	/// Returns the difference of two vectors
	const DVector3d	operator-(const DPoint3d& pt) const { return DVector3d(x-pt.x, y-pt.y, z-pt.z); }
	/// Translates the point by the given vector
	DPoint3d&	operator-=(const DVector3d& v){ x-=v.x; y-=v.y; z-=v.z; return *this; }
	/// Returns the point scaled by the given factor
	DPoint3d	operator*(double d) const { return DPoint3d(d*x, d*y, d*z); }
	/// Transforms to fixed vector (from [0,0,0])
	const DVector3d   fixedVector() const { return DVector3d(x, y, z); }

	bool operator==( const DPoint3d& pt ) const { return x == pt.x && y == pt.y && z == pt.z; }
	bool operator!=( const DPoint3d& pt ) const { return x != pt.x || y != pt.y || z != pt.z; }

	double& operator[](Axis axis) { return *(&x + static_cast<int>(axis)); }
	const double& operator[](Axis axis) const { return *(&x + static_cast<int>(axis)); }

	/// Divide point by number
	DPoint3d&	operator/=(double d){ x/=d; y/=d; z/=d; return *this; }
	/// Multiply point by number
	DPoint3d&	operator*=(double d){ x*=d; y*=d; z*=d; return *this; }
	/// Set to coordinates
	DPoint3d& set(const DPoint3d& pt, double f = 1.0) {
		x = pt.x * f; y = pt.y * f; z = pt.z * f; return *this; 
	}
	/// Add to coordinates
	DPoint3d& add(const DPoint3d& pt, double f) {
		x += pt.x * f; y += pt.y * f; z += pt.z * f; return *this; 
	}
	/// Add to coordinates
	DPoint3d& add(const DPoint3d& pt) {
		x += pt.x; y += pt.y; z += pt.z; return *this; 
	}
	/// Average of two points
	static const DPoint3d average(const DPoint3d& a, const DPoint3d& b) {
		return DPoint3d((a.x+b.x)*0.5, (a.y+b.y)*0.5, (a.z+b.z)*0.5);
	}
	/// Average of three points
	static const DPoint3d average(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c) {
		return DPoint3d((a.x+b.x+c.x)*(1.0/3.0), (a.y+b.y+c.y)*(1.0/3.0), (a.z+b.z+c.z)*(1.0/3.0));
	}
	/// Average of four points
	static const DPoint3d average(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d) {
		return DPoint3d((a.x+b.x+c.x+d.x)*0.25, (a.y+b.y+c.y+d.y)*0.25, (a.z+b.z+c.z+d.z)*0.25);
	}
	/// Average of five points
	static const DPoint3d average(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d, const DPoint3d& e) {
		return DPoint3d((a.x+b.x+c.x+d.x+e.x)*0.2, (a.y+b.y+c.y+d.y+e.y)*0.2, (a.z+b.z+c.z+d.z+e.z)*0.2);
	}

	/// Returns the distance between two points
	double distance(const DPoint3d& a) const { return sqrt(sqr(a.x - x) + sqr(a.y - y) + sqr(a.z - z)); }
	/// Returns the squared distance between two points
	double distance2(const DPoint3d& a) const { return sqr(a.x - x) + sqr(a.y - y) + sqr(a.z - z); }

	OctVertexWhich vertexDirection(const DPoint3d& pt) const;

	bool storeSimple(ostream& os, char end_char = '\n') const;
	bool readSimple(istream& is);

	/// Stores the coordinates of the point into the stream
	friend ostream& operator<<(ostream& os, const DPoint3d& pt);
	/// Loads the coordinates of the point from the stream
	friend istream& operator>>(istream& is, DPoint3d& pt);

public:
	/// x-coordinate of the point
	double x;
	/// y-coordinate of the point
	double y;
	/// z-coordinate of the point
	double z;
	// zero point
	static const DPoint3d zero;
	// +inf point
	static const DPoint3d inf;
};

class FPoint3d
{
public:
	FPoint3d(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) : x(_x), y(_y), z(_z) {}
	FPoint3d(const FPoint3d& a, const FPoint3d& b, float t) 
		: x(a.x * (1.f-t) + b.x * t), y(a.y * (1.f-t) + b.y * t), z(a.z * (1.f-t) + b.z * t) {} 
	FPoint3d(const DPoint3d& pt) : x((float)pt.x), y((float)pt.y), z((float)pt.z) {}
public:
	const FVector3d	operator-(const FPoint3d& pt) const { return FVector3d(x-pt.x, y-pt.y, z-pt.z); }
	const FPoint3d	operator+(const FVector3d& v) const { return FPoint3d(x+v.x, y+v.y, z+v.z); }
	FPoint3d&	operator+=(const FVector3d& v){ x+=v.x; y+=v.y; z+=v.z; return *this; }
	FPoint3d& add(const FPoint3d& pt, float f) {
		x += pt.x * f; y += pt.y * f; z += pt.z * f; return *this; 
	}
	friend ostream& operator<<(ostream& os, const FPoint3d& pt) {
		return os << pt.x << " " << pt.y << " " << pt.z;
	}
public:
	float x,y,z;
	static const FPoint3d zero;
};

/**
 * This class implements a point in a three-dimensional metric space
 *  plus some basic operations
 */
class DMPoint3d
{
public:
	/// Standard constructor
	DMPoint3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	/// linear
	DMPoint3d(const DMPoint3d& a, const DMPoint3d& b, double t) 
		: x(a.x * (1.0-t) + b.x * t), y(a.y * (1.0-t) + b.y * t), z(a.z * (1.0-t) + b.z * t) {} 
	/// Standard (empty) constructor
	DMPoint3d() : x(0.0), y(0.0), z(0.0) {}
	/// Standard constructor
	DMPoint3d(const DMPoint2d& pt, double _z) : x(pt.x), y(pt.y), z(_z) {}
public:
	/// Translates the point by the given vector
	const DMPoint3d	operator+(const DMVector3d& v) const { return DMPoint3d(x+v.x, y+v.y, z+v.z); }
	/// Translates the point by the given vector
	const DMPoint3d	operator-(const DMVector3d& v) const { return DMPoint3d(x-v.x, y-v.y, z-v.z); }
	/// Translates the point by the given vector
	DMPoint3d&	operator+=(const DMVector3d& v){ x+=v.x; y+=v.y; z+=v.z; return *this; }
	/// Returns the difference of two vectors
	DMVector3d	operator-(const DMPoint3d& pt) const { return DMVector3d(x-pt.x, y-pt.y, z-pt.z); }

	/// Average of two points
	static const DMPoint3d average(const DMPoint3d& a, const DMPoint3d& b) {
		return DMPoint3d((a.x+b.x)*0.5, (a.y+b.y)*0.5, (a.z+b.z)*0.5);
	}
	/// Average of three points
	static const DMPoint3d average(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c) {
		return DMPoint3d((a.x+b.x+c.x)*(1.0/3.0), (a.y+b.y+c.y)*(1.0/3.0), (a.z+b.z+c.z)*(1.0/3.0));
	}

	/// simple conversion
	const DPoint3d toRealSpace() const { return DPoint3d(x, y, z); }

	/// Divide point by number
	DMPoint3d&	operator/=(double d){ x/=d; y/=d; z/=d; return *this; }
	/// Add to coordinates
	DMPoint3d& add(const DMPoint3d& pt) {
		x += pt.x; y += pt.y; z += pt.z; return *this; 
	}
	/// Add to coordinates
	DMPoint3d& add(const DMPoint3d& pt, double f) {
		x += pt.x * f; y += pt.y * f; z += pt.z * f; return *this; 
	}
	/// Returns the distance between two points
	double distance(const DMPoint3d& a) const { return sqrt(sqr(a.x - x) + sqr(a.y - y) + sqr(a.z - z)); }
	/// Returns the squared distance between two points
	double distance2(const DMPoint3d& a) const { return sqr(a.x - x) + sqr(a.y - y) + sqr(a.z - z); }

	/// Stores the coordinates of the point into the stream
	friend ostream& operator<<(ostream& os, const DMPoint3d& pt);
	/// Loads the coordinates of the point from the stream
	friend istream& operator>>(istream& is, DMPoint3d& pt);

public:
	/// x-coordinate of the point
	double x;
	/// y-coordinate of the point
	double y;
	/// z-coordinate of the point
	double z;
	// zero point
	static const DMPoint3d zero;
};

#endif // !defined(DPOINT_H__INCLUDED)
