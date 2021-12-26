/////////////////////////////////////////////////////////////////////////////
// DVector.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DVECTOR_H__INCLUDED
#define DVECTOR_H__INCLUDED

class DPoint2d;
class DMPoint2d;
class FPoint3d;
class DPoint3d;
class DMPoint3d;

#include "common.h"
#include <iostream>

/**
 * This class implements a vector in a two-dimensional space
 *  plus some basic operations.
 */
class DVector2d
{
public:
	/// Standard constructor
	DVector2d(double _x, double _y) : x(_x), y(_y) {}
	/// Standard constructor
	DVector2d(const DPoint2d& v0, const DPoint2d& v1);
	/// Standard (emvy) constructor
	DVector2d() : x(0.0), y(0.0) {}
public:
	/// Returns the sum of two vectors
	const DVector2d	operator+(const DVector2d& v) const { return DVector2d(x+v.x, y+v.y); }
	/// Adds two vectors
	DVector2d&	operator+=(const DVector2d& v){ x+=v.x; y+=v.y; return *this; }
	/// Returns the difference of two vectors
	const DVector2d	operator-(const DVector2d& v) const { return DVector2d(x-v.x, y-v.y); }
	/// Substracts two vectors
	DVector2d&	operator-=(const DVector2d& v){ x-=v.x; y-=v.y; return *this; }
	/// Returns the vector scaled by the given factor
	const DVector2d	operator*(double d) const { return DVector2d(d*x, d*y); }
	/// Returns the vector scaled by the given factor
	const DVector2d	operator/(double d) const { return DVector2d(x/d, y/d); }
	/// Scales up the vector by the given factor
	DVector2d&	operator*=(double d){ x*=d, y*=d; return *this; }
	/// Scales down the vector by the given factor
	DVector2d&	operator/=(double d){ assert(d!=0.0); x/=d; y/=d; return *this; }

	double& operator[](Axis axis) { return *(&x + static_cast<int>(axis)); }
	const double& operator[](Axis axis) const { return *(&x + static_cast<int>(axis)); }

	/// Returns the length of the vector
	double length() const { return sqrt(x*x + y*y); }
	/// Returns the squared length of the vector
	double length2() const { return x*x + y*y; }

	/// Rotates the vector counterclockwise by the given angle (expressed by sin and cos values)
	void turn(double sinAlpha, double cosAlpha);
	/// Returns the vector rotated counterclockwise by the given angle (expressed by sin and cos values)
	const DVector2d turned(double sinAlpha, double cosAlpha) const;
	/// Returns orthogonal vector
	const DVector2d orthogonal() const { return DVector2d(-y, x); }
	/// Returns orthonormal vector
	const DVector2d orthonormal() const { return orthogonal().normalized(); }
	/// Returns the cross product of vectors |a| and |b|
	double crossProduct(const DVector2d& vb) const { return x * vb.y - vb.x * y; }
	/// Returns the scalar product of vectors |a| and |b|
	double scalarProductNormalized(const DVector2d& vb) const;
	/// Returns the scalar product of vectors va and vb
	double scalarProduct(const DVector2d& vb) const { return x * vb.x + y * vb.y; }
	/// Returns the angle between vectors |a| and |b|
	double getAngle(const DVector2d& vb) const;
	/// Returns the cosinus angle between vectors |a| and |b|
	double getAngleCos(const DVector2d& vb) const;
	/// Reverse vector
	DVector2d& reverse() { x = -x; y = -y; return *this; }
	/// Normalizes the vector (makes its length equal to 1)
	DVector2d& normalize(const double& nlen = 1.0);
	/// Return normalized vector
	const DVector2d normalized(const double& nlen = 1.0) const;

	/// Returns the cosinus angle between vectors |ab| and |ac|
	static double angle(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c);
	/// is zero?
	bool isZero() const { return (x == 0.0) && (y == 0.0); }

	/// Stores the coordinates of the v into the stream
	friend std::ostream& operator<<(std::ostream& os, const DVector2d& v);
	/// Loads the coordinates of the v from the stream
	friend std::istream& operator>>(std::istream& is, DVector2d& v);

public:
	/// x-coordinate of the vector
	double x; 
	/// y-coordinate of the vector
	double y; 
	// zero vector
	static const DVector2d zero;
	// ox vector
	static const DVector2d v_ox;
	// oy vector
	static const DVector2d v_oy;
};

/**
 * This class implements a vector in a two-dimensional metric space
 *  plus some basic operations.
 */
class DMVector2d
{
public:
	/// Standard constructor
	DMVector2d(double _x, double _y) : x(_x), y(_y) {}
	/// Standard constructor
	DMVector2d(const DMPoint2d& v0, const DMPoint2d& v1);
	/// Standard (emvy) constructor
	DMVector2d() : x(0.0), y(0.0) {}
public:
	/// Returns the sum of two vectors
	const DMVector2d	operator+(const DMVector2d& v) const { return DMVector2d(x+v.x, y+v.y); }
	/// Adds two vectors
	DMVector2d&	operator+=(const DMVector2d& v){ x+=v.x; y+=v.y; return *this; }
	/// Returns the difference of two vectors
	const DMVector2d	operator-(const DMVector2d& v) const { return DMVector2d(x-v.x, y-v.y); }
	/// Substracts two vectors
	DMVector2d&	operator-=(const DMVector2d& v){ x-=v.x; y-=v.y; return *this; }
	/// Returns the vector scaled by the given factor
	const DMVector2d	operator*(double d) const { return DMVector2d(d*x, d*y); }
	/// Returns the vector scaled by the given factor
	const DMVector2d	operator/(double d) const { return DMVector2d(x/d, y/d); }
	/// Scales up the vector by the given factor
	DMVector2d&	operator*=(double d){ x*=d, y*=d; return *this; }
	/// Scales down the vector by the given factor
	DMVector2d&	operator/=(double d){ assert(d!=0.0); x/=d; y/=d; return *this; }

	/// Returns the length of the vector
	double length() const { return sqrt(x*x + y*y); }
	/// Returns the squared length of the vector
	double length2() const { return x*x + y*y; }

	/// Rotates the vector counterclockwise by the given angle (expressed by sin and cos values)
	void turn(double sinAlpha, double cosAlpha);
	/// Returns the vector rotated counterclockwise by the given angle (expressed by sin and cos values)
	const DMVector2d turned(double sinAlpha, double cosAlpha) const;
	/// Returns the cross product of vectors |a| and |b|
	double crossProduct(const DMVector2d& v) const { return x*v.y - v.x*y; }
	/// Returns the scalar product of vectors |a| and |b|
	double scalarProductNormalized(const DMVector2d& v) const;
	/// Returns the scalar product of vectors va and vb
	double scalarProduct(const DMVector2d& vb) const { return x * vb.x + y * vb.y; }
	/// Returns the angle between vectors |a| and |b|
	double getAngle(const DMVector2d& v) const;
	/// Returns the cosinus angle between vectors |a| and |b|
	double getAngleCos(const DMVector2d& v) const;
	/// Normalizes the vector (makes its length equal to 1)
	DMVector2d& normalize();
	/// Return normalized vector
	const DMVector2d normalized() const;
	/// is zero?
	bool isZero() const { return (x == 0.0) && (y == 0.0); }

	/// Returns the cosinus angle between vectors |ab| and |ac|
	static double angleCos(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c);
	/// Returns the angle between vectors |ab| and |ac|
	static double angle(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c);

	/// Stores the coordinates of the v into the stream
	friend std::ostream& operator<<(std::ostream& os, const DMVector2d& v);
	/// Loads the coordinates of the v from the stream
	friend std::istream& operator>>(std::istream& is, DMVector2d& v);

public:
	/// x-coordinate of the vector
	double x; 
	/// y-coordinate of the vector
	double y; 
};

/**
 * This class implements a vector in a three-dimensional space
 *  plus some basic operations
 */
class DVector3d
{
public:
	/// Standard constructor
	DVector3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	/// Standard constructor
	DVector3d(const DPoint3d& pt0, const DPoint3d& pt1);
	/// Standard (emvy) constructor
	DVector3d() : x(0.0), y(0.0), z(0.0) {}
	/// Standard constructor
	DVector3d(const DVector2d& v, double _z) : x(v.x), y(v.y), z(_z) {}
public:
	/// Add to coordinates
	DVector3d& add(const DVector3d& v, double f);
	/// Returns the sum of two vectors
	const DVector3d	operator+(const DVector3d& v) const { return DVector3d(x+v.x, y+v.y, z+v.z); }
	/// Adds two vectors
	DVector3d&	operator+=(const DVector3d& v){ x+=v.x; y+=v.y; z+=v.z; return *this; }
	/// Returns the difference of two vectors
	const DVector3d	operator-(const DVector3d& v) const { return DVector3d(x-v.x, y-v.y, z-v.z); }
	/// Returns the vector scaled by the given factor
	const DVector3d	operator*(double d) const { return DVector3d(d*x, d*y, d*z); }
	/// Returns the vector scaled by the given factor
	const DVector3d	operator/(double d) const { assert(d!=0.0); return DVector3d(x/d, y/d, z/d); }
	/// Scales up the vector by the given factor
	DVector3d&	operator*=(double d){ x*=d, y*=d, z*=d; return *this; }
	/// Scales down the vector by the given factor
	DVector3d&	operator/=(double d){ assert(d!=0.0); x/=d; y/=d; z/=d; return *this; }
	/// Returns negation of the vector
	const DVector3d operator-() const { return DVector3d(-x, -y, -z); }
	/// array operator
	double& operator[](int i) { assert( i >= 0 && i <= 2 ); return (i==0)?x:((i==1)?y:z); }
	const double& operator[](int i) const { assert( i >= 0 && i <= 2 ); return (i==0)?x:((i==1)?y:z); }
	double& operator[](Axis axis) { return *(&x + static_cast<int>(axis)); }
	const double& operator[](Axis axis) const { return *(&x + static_cast<int>(axis)); }

	DVector3d rotatedAroundAxis(Axis axis, double sin_alpha, double cos_alpha) const;

	/// is zero?
	bool isZero() const { return (x == 0.0) && (y == 0.0) && (z == 0.0); }
	/// Reverse vector
	DVector3d& reverse() { x = -x; y = -y; z = -z; return *this; }
	/// Normalizes the vector (makes its length equal to 1)
	DVector3d& normalize(const double& nlen = 1.0);
	/// Returnes the normalized vector (length equal to 1)
	const DVector3d normalized(const double& nlen = 1.0) const;
	/// is normalized (with epsilon)
	bool isNormal( const double &nlen = 1.0, const double & eps = 1e-5 ) const {
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, std::abs( nlen - length2() ) << " xxx " << eps);
		return std::abs( nlen - length2() ) < eps; 
	}
	/// Returns the length of the vector
	double length() const { return sqrt(x*x + y*y + z*z); }
	/// Returns the squared length of the vector
	double length2() const { return x*x + y*y + z*z; }
	/// Returns the cross product of vectors va and vb
	const DVector3d crossProduct(const DVector3d& vb) const;
	/// Returns the normalized scalar product of vectors va and vb
	double scalarProductNormalized(const DVector3d& b) const;
	/// Returns the scalar product of vectors va and vb
	double scalarProduct(const DVector3d& vb) const { return x * vb.x + y * vb.y + z * vb.z; }
	/// Returns the angle between vectors |a| and |b|
	double getAngle(const DVector3d& b) const;
	/// Counts two orthogonal vectors 
	void orthogonalVectors(DVector3d& v1, DVector3d& v2) const;
	/// Counts two orthonormal vectors
	void orthonormalVectors(DVector3d& v1, DVector3d& v2) const;

	/// Returns the cross product of vectors |ab| and |ac|
	static const DVector3d crossProduct(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c);

	/// Stores the coordinates of the v into the stream
	friend ostream& operator<<(ostream& os, const DVector3d& v);
	/// Loads the coordinates of the v from the stream
	friend istream& operator>>(istream& is, DVector3d& v);
public:
	/// x-coordinate of the v (vector)
	double x;
	/// y-coordinate of the v (vector)
	double y;
	/// z-coordinate of the v (vector)
	double z;
	// zero vector
	static const DVector3d zero;
	// ox vector
	static const DVector3d v_ox;
	static const DVector3d v_oy;
	static const DVector3d v_oz;
};

class FVector3d
{
public:
	FVector3d(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) : x(_x), y(_y), z(_z) {}
	FVector3d(const DVector3d& vt) : x((float)vt.x), y((float)vt.y), z((float)vt.z) {}
public:
	const FVector3d	operator*(float d) const { return FVector3d(d*x, d*y, d*z); }
	float length() const { return sqrtf(x*x + y*y + z*z); }
	void invert() { x = -x; y = -y; z = -z; }
	bool isZero() const { return (x == 0.0) && (y == 0.0) && (z == 0.0); }
	const FVector3d normalized() const;
	const FVector3d crossProduct(const FVector3d& vb) const {
		return FVector3d(y*vb.z - vb.y*z, vb.x*z - x*vb.z, x*vb.y - vb.x*y);
	}
	/// Returns the cross product of vectors |ab| and |ac|
	static const FVector3d crossProduct(const FPoint3d& a, const FPoint3d& b, const FPoint3d& c);
	float scalarProduct(const FVector3d& vb) const { return x * vb.x + y * vb.y + z * vb.z; }
public:
	float x,y,z;
	// zero vector
	static const FVector3d zero;
};

/**
 * This class implements a vector in a three-dimensional space
 *  plus some basic operations
 */
class DMVector3d
{
public:
	/// Standard constructor
	DMVector3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	/// Standard constructor
	DMVector3d(const DMPoint3d& pt0, const DMPoint3d& pt1);
	/// Standard (emvy) constructor
	DMVector3d() : x(0.0), y(0.0), z(0.0) {}
	/// Standard constructor
	DMVector3d(const DMVector2d& v, double _z) : x(v.x), y(v.y), z(_z) {}
public:
	/// Returns the sum of two vectors
	const DMVector3d	operator+(const DMVector3d& v) const { return DMVector3d(x+v.x, y+v.y, z+v.z); }
	/// Adds two vectors
	DMVector3d&	operator+=(const DMVector3d& v){ x+=v.x; y+=v.y; z+=v.z; return *this; }
	/// Returns the difference of two vectors
	const DMVector3d	operator-(const DMVector3d& v) const { return DMVector3d(x-v.x, y-v.y, z-v.z); }
	/// Returns the vector scaled by the given factor
	const DMVector3d	operator*(double d) const { return DMVector3d(d*x, d*y, d*z); }
	/// Returns the vector scaled by the given factor
	const DMVector3d	operator/(double d) const { assert(d!=0.0); return DMVector3d(x/d, y/d, z/d); }
	/// Scales up the vector by the given factor
	DMVector3d&	operator*=(double d){ x*=d, y*=d, z*=d; return *this; }
	/// Scales down the vector by the given factor
	DMVector3d&	operator/=(double d){ assert(d!=0.0); x/=d; y/=d; z/=d; return *this; }

	/// is zero?
	bool isZero() const { return (x == 0.0) && (y == 0.0) && (z == 0.0); }
	/// Normalizes the vector (makes its length equal to 1)
	DMVector3d& normalize();
	/// Returnes the normalized vector (length equal to 1)
	const DMVector3d normalized() const;
	/// Returns the length of the vector
	double length() const { return sqrt(x*x + y*y + z*z); }
	/// Returns the squared length of the vector
	double length2() const { return x*x + y*y + z*z; }
	/// Returns the cross product of vectors va and vb
	const DMVector3d crossProduct(const DMVector3d& vb) const;
	/// Returns the normalized scalar product of vectors va and vb
	double scalarProductNormalized(const DMVector3d& b) const;
	/// Returns the scalar product of vectors va and vb
	double scalarProduct(const DMVector3d& vb) const { return x * vb.x + y * vb.y + z * vb.z; }
	/// Returns the angle between vectors |a| and |b|
	double getAngle(const DMVector3d& b) const;

	/// Returns the cross product of vectors |ab| and |ac|
	static const DMVector3d crossProduct(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c);

	/// Stores the coordinates of the v into the stream
	friend ostream& operator<<(ostream& os, const DMVector3d& v);
	/// Loads the coordinates of the v from the stream
	friend istream& operator>>(istream& is, DMVector3d& v);

public:
	/// x-coordinate of the v (vector)
	double x;
	/// y-coordinate of the v (vector)
	double y;
	/// z-coordinate of the v (vector)
	double z;
};

#endif // !defined(DVECTOR_H__INCLUDED)
