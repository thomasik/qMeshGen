/////////////////////////////////////////////////////////////////////////////
// DVector.cpp
// This class implements a vector in a two-dimensional space
//  plus some basic operations.
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DVector.h"
#include "DPoint.h"
#include "DEquation.h"

const DVector2d DVector2d::zero(0.0, 0.0);
const DVector2d DVector2d::v_ox(1.0, 0.0);
const DVector2d DVector2d::v_oy(0.0, 1.0);

DVector2d::DVector2d(const DPoint2d& v0, const DPoint2d& v1) : x(v1.x-v0.x), y(v1.y-v0.y) {}

/**
 * \param sinAlpha sine of rotation \a angle
 * \param cosAlpha cosine of rotation \a angle
 * \remark attributes ::x and ::y are changed
 * \return \c void
 */
void DVector2d::turn(double sinAlpha, double cosAlpha)
{
	double tx = x;
	x = tx * cosAlpha - y * sinAlpha;
	y = tx * sinAlpha + y * cosAlpha;
}

/**
 * \param sinAlpha sine of rotation \a angle
 * \param cosAlpha cosine of rotation \a angle
 * \return \c rotated point
 */
const DVector2d DVector2d::turned(double sinAlpha, double cosAlpha) const
{
	return DVector2d(x * cosAlpha - y * sinAlpha,
		x * sinAlpha + y * cosAlpha);
}

/**
 * \param v second vector
 * \remark both vectors are normalized before calculating the scalar product
 * \return \c double from range [-1.0, 1.0]
 */
double DVector2d::scalarProductNormalized(const DVector2d &v) const
{
	double len = sqrt(this->length2() * v.length2());
	double value = (x * v.x + y * v.y) / len;
	if(value < -1.0) return -1.0;
	else if(value > 1.0) return 1.0;
	else return value;
}

/**
 * \param v the second vector
 * \remark The angle |a|->|b| in radians, 2-dimensional version
 * \return \c double from range [0, 2PI]
 */
double DVector2d::getAngle(const DVector2d &v) const
{
	double alpha = acos(this->scalarProductNormalized(v));
	return (this->crossProduct(v) < 0) ? (2*PI - alpha) : alpha;
}

/**
 * \param b the second vector
 * \remark The adjusted cosinus of angle |a|->|b|, 2-dimensional version
 * \return \c double from range [0, 4]
 */
double DVector2d::getAngleCos(const DVector2d &b) const
{
	double cos_alpha = this->scalarProductNormalized(b);
	if(this->crossProduct(b) < 0){
		return 3+cos_alpha; // [2,4)
	}else{
		return 1-cos_alpha; // [0,2)
	}
}

DVector2d& DVector2d::normalize(const double& nlen)
{
	double len_ratio = nlen / length();
	x *= len_ratio;
	y *= len_ratio;
	return *this;
}

const DVector2d DVector2d::normalized(const double& nlen) const
{
	double len_ratio = nlen / length();
	return DVector2d(x*len_ratio, y*len_ratio);
}

/// Returns the cosinus angle between vectors |ab| and |ac|
double DVector2d::angle(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c)
{
	const DVector2d ab = b - a;
	const DVector2d ac = c - a;

	double alpha = acos(ab.scalarProductNormalized(ac));
	return (ab.crossProduct(ac) < 0) ? (2*PI - alpha) : alpha;
}

ostream& operator<<(ostream& os, const DVector2d& v){
	string str_x, str_y;
	DEquation::doubleToString(v.x, DEquation::v_length, str_x);
	DEquation::doubleToString(v.y, DEquation::v_length, str_y);
	os << "<x>" << str_x << "</x> <y>" << str_y << "</y>";
	return os;
}

istream& operator>>(istream& is, DVector2d& v){
	string str_x, str_y;
	string tmp;
	char z;
	getline(is, tmp, '>'); is >> z; // <x>
	getline(is, str_x, '<');
	getline(is, tmp, '>'); is >> z; // </x>
	getline(is, tmp, '>'); is >> z; // <y>
	getline(is, str_y, '<');
	getline(is, tmp, '>'); is >> z; // </y>
	DEquation eq;
	if(eq.parse(str_x.c_str()))	v.x = eq.getValue(0.0);
	if(eq.parse(str_y.c_str()))	v.y = eq.getValue(0.0);
	return is;
}

DMVector2d::DMVector2d(const DMPoint2d& v0, const DMPoint2d& v1) : x(v1.x-v0.x), y(v1.y-v0.y) {}

/**
 * \param sinAlpha sine of rotation \a angle
 * \param cosAlpha cosine of rotation \a angle
 * \remark attributes ::x and ::y are changed
 * \return \c void
 */
void DMVector2d::turn(double sinAlpha, double cosAlpha)
{
	double tx = x;
	x = tx * cosAlpha - y * sinAlpha;
	y = tx * sinAlpha + y * cosAlpha;
}

/**
 * \param sinAlpha sine of rotation \a angle
 * \param cosAlpha cosine of rotation \a angle
 * \return \c rotated point
 */
const DMVector2d DMVector2d::turned(double sinAlpha, double cosAlpha) const
{
	return DMVector2d(x * cosAlpha - y * sinAlpha,
		x * sinAlpha + y * cosAlpha);
}

/**
 * \param v second vector
 * \remark both vectors are normalized before calculating the scalar product
 * \return \c double from range [-1.0, 1.0]
 */
double DMVector2d::scalarProductNormalized(const DMVector2d &v) const
{
	double len = sqrt(this->length2() * v.length2());
	double value = (x * v.x + y * v.y) / len;
	if(value < -1.0) return -1.0;
	else if(value > 1.0) return 1.0;
	else return value;
}

/**
 * \param v the second vector
 * \remark The angle |a|->|b| in radians, 2-dimensional version
 * \return \c double from range [0, 2PI]
 */
double DMVector2d::getAngle(const DMVector2d &v) const
{
	double alpha = acos(this->scalarProductNormalized(v));
	return (this->crossProduct(v) < 0) ? (2*PI - alpha) : alpha;
}

/**
 * \param b the second vector
 * \remark The adjusted cosinus of angle |a|->|b|, 2-dimensional version
 * \return \c double from range [0, 4]
 */
double DMVector2d::getAngleCos(const DMVector2d &b) const
{
	double cos_alpha = this->scalarProductNormalized(b);
	if(this->crossProduct(b) < 0){
		return 3+cos_alpha; // [2,4)
	}else{
		return 1-cos_alpha; // [0,2)
	}
}

DMVector2d& DMVector2d::normalize()
{
	double len_ratio = 1.0 / length();
	x *= len_ratio;
	y *= len_ratio;
	return *this;
}

const FVector3d FVector3d::normalized() const
{
	float len_ratio = 1.0f / length();
	return FVector3d(x*len_ratio, y*len_ratio, z*len_ratio);
}

const DMVector2d DMVector2d::normalized() const
{
	double len_ratio = 1.0 / length();
	return DMVector2d(x*len_ratio, y*len_ratio);
}

/// Returns the cosinus angle between vectors |ab| and |ac|
double DMVector2d::angleCos(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c)
{
	const DMVector2d ab = b - a;
	const DMVector2d ac = c - a;
	double cos_alpha = ab.scalarProductNormalized(ac);
	if(ab.crossProduct(ac) < 0){
		return 3+cos_alpha; // [2,4)
	}else{
		return 1-cos_alpha; // [0,2)
	}
}

/// Returns the cangle between vectors |ab| and |ac|
double DMVector2d::angle(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c)
{
	const DMVector2d ab = b - a;
	const DMVector2d ac = c - a;
	double alpha = acos(ab.scalarProductNormalized(ac));
	return (ab.crossProduct(ac) < 0) ? (2*PI - alpha) : alpha;
}

ostream& operator<<(ostream& os, const DMVector2d& v){
	string str_x, str_y;
	DEquation::doubleToString(v.x, DEquation::v_length, str_x);
	DEquation::doubleToString(v.y, DEquation::v_length, str_y);
	os << "<x>" << str_x << "</x> <y>" << str_y << "</y>";
	return os;
}

istream& operator>>(istream& is, DMVector2d& v){
	string str_x, str_y;
	string tmp;
	char z;
	getline(is, tmp, '>'); is >> z; // <x>
	getline(is, str_x, '<');
	getline(is, tmp, '>'); is >> z; // </x>
	getline(is, tmp, '>'); is >> z; // <y>
	getline(is, str_y, '<');
	getline(is, tmp, '>'); is >> z; // </y>
	DEquation eq;
	if(eq.parse(str_x.c_str()))	v.x = eq.getValue(0.0);
	if(eq.parse(str_y.c_str()))	v.y = eq.getValue(0.0);
	return is;
}

const FVector3d FVector3d::crossProduct(const FPoint3d& a, const FPoint3d& b, const FPoint3d& c) {
	return (b-a).crossProduct(c-a); 
}

const DVector3d DVector3d::zero(0.0, 0.0, 0.0);
const DVector3d DVector3d::v_ox(1.0, 0.0, 0.0);
const DVector3d DVector3d::v_oy(0.0, 1.0, 0.0);
const DVector3d DVector3d::v_oz(0.0, 0.0, 1.0);

DVector3d::DVector3d(const DPoint3d& pt0, const DPoint3d& pt1) 
	: x(pt1.x-pt0.x), y(pt1.y-pt0.y), z(pt1.z-pt0.z) 
{}

//////////////////////////////////////////////////////////////////////
// Wylicza iloczyn skalarny (cosinus k¹ta) dwóch wektorów 
double DVector3d::scalarProductNormalized(const DVector3d &b) const
{
	double len = this->length() * b.length();
	double value = (x * b.x + y * b.y + z * b.z) / len;
	if(value < -1.0){
		return -1.0;
	}else if(value > 1.0){
		return 1.0;
	}else
		return value;
}

//////////////////////////////////////////////////////////////////////
// cross product of two given vectors
const DVector3d DVector3d::crossProduct(const DVector3d &vb) const
{
	return DVector3d(y*vb.z - vb.y*z, vb.x*z - x*vb.z, x*vb.y - vb.x*y);
}

/// Add to coordinates
DVector3d& DVector3d::add(const DVector3d& v, double f) {
	x += v.x * f; y += v.y * f; z += v.z * f; return *this; 
}

DVector3d DVector3d::rotatedAroundAxis(Axis axis, double sa, double ca) const
{
	switch (axis) {
		case Axis::X:
			return DVector3d(x, y*ca - z*sa,  y*sa + z*ca);
		case Axis::Y:
			return DVector3d(x*ca + z*sa, y, -x*sa + z*ca);
		case Axis::Z:
			return DVector3d(x*ca - y*sa,  x*sa + y*ca, z);
		default:
			assert(false);
			return DVector3d::zero;
	}
}

DVector3d& DVector3d::normalize(const double& nlen)
{
	double len = length();
	assert( len > 0.0 );
	double len_ratio = nlen / len;
	x *= len_ratio;
	y *= len_ratio;
	z *= len_ratio;
	return *this;
}

const DVector3d DVector3d::normalized(const double& nlen) const
{
	double len = length();
	assert( len > 0.0 );
	double len_ratio = nlen / len;
	return DVector3d(x*len_ratio, y*len_ratio, z*len_ratio);
}

//////////////////////////////////////////////////////////////////////
// Zwraca wartoœæ k¹ta o pomiêdzy |a|->|b|.
double DVector3d::getAngle(const DVector3d &b) const
{
	return acos(this->scalarProductNormalized(b));
}

void DVector3d::orthogonalVectors(DVector3d& v1, DVector3d& v2) const
{
	if(abs(x) < abs(y) && abs(x) < abs(z)){
		v1.x = 0.0;
		v1.y = z;
		v1.z = -y;
	}else if(abs(y) < abs(z)){
		v1.x = z;
		v1.y = 0.0;
		v1.z = -x;
	}else{
		v1.x = y;
		v1.y = -x;
		v1.z = 0.0;
	}

	v2 = this->crossProduct(v1);
}

void DVector3d::orthonormalVectors(DVector3d& v1, DVector3d& v2) const
{
	orthogonalVectors(v1, v2);
	v1.normalize();
	v2.normalize();
}

ostream& operator<<(ostream& os, const DVector3d& v){
	string str_x, str_y, str_z;
	DEquation::doubleToString(v.x, DEquation::v_length, str_x);
	DEquation::doubleToString(v.y, DEquation::v_length, str_y);
	DEquation::doubleToString(v.z, DEquation::v_length, str_z);
	os << "<x>" << str_x << "</x> <y>" << str_y << "</y> <z>" << str_z  << "</z>";
	return os;
}

istream& operator>>(istream& is, DVector3d& v){
	string str_x, str_y, str_z;
	string tmp;
	char z;
	getline(is, tmp, '>'); is >> z; // <x>
	getline(is, str_x, '<');
	getline(is, tmp, '>'); is >> z; // </x>
	getline(is, tmp, '>'); is >> z; // <y>
	getline(is, str_y, '<');
	getline(is, tmp, '>'); is >> z; // </y>
	getline(is, tmp, '>'); is >> z; // <z>
	getline(is, str_z, '<');
	getline(is, tmp, '>'); is >> z; // </z>
	DEquation eq;
	if(eq.parse(str_x.c_str()))	v.x = eq.getValue(0.0);
	if(eq.parse(str_y.c_str()))	v.y = eq.getValue(0.0);
	if(eq.parse(str_z.c_str()))	v.z = eq.getValue(0.0);
	return is;
}

const FVector3d FVector3d::zero(0.0f, 0.0f, 0.0f);

DMVector3d::DMVector3d(const DMPoint3d& pt0, const DMPoint3d& pt1) 
	: x(pt1.x-pt0.x), y(pt1.y-pt0.y), z(pt1.z-pt0.z) 
{}

//////////////////////////////////////////////////////////////////////
// Wylicza iloczyn skalarny (cosinus k¹ta) dwóch wektorów 
double DMVector3d::scalarProductNormalized(const DMVector3d &b) const
{
	double len = this->length() * b.length();
	double value = (x * b.x + y * b.y + z * b.z) / len;
	if(value < -1.0){
		return -1.0;
	}else if(value > 1.0){
		return 1.0;
	}else
		return value;
}

//////////////////////////////////////////////////////////////////////
// Wylicza iloczyn wektorowy dwóch podanych wektorów
const DMVector3d DMVector3d::crossProduct(const DMVector3d &vb) const
{
	return DMVector3d(y*vb.z - vb.y*z, vb.x*z - x*vb.z, x*vb.y - vb.x*y);
}

DMVector3d& DMVector3d::normalize()
{
	double len_ratio = 1.0 / length();
	x *= len_ratio;
	y *= len_ratio;
	z *= len_ratio;
	return *this;
}

const DMVector3d DMVector3d::normalized() const
{
	double len_ratio = 1.0 / length();
	return DMVector3d(x*len_ratio, y*len_ratio, z*len_ratio);
}

//////////////////////////////////////////////////////////////////////
// Zwraca wartoœæ k¹ta o pomiêdzy |a|->|b|. [0-2PI)
double DMVector3d::getAngle(const DMVector3d &b) const
{
	return acos(this->scalarProductNormalized(b));
}

/// Returns the cross product of vectors |ab| and |ac|
const DVector3d DVector3d::crossProduct(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c) 
{
	return (b-a).crossProduct(c-a); 
}

/// Returns the cross product of vectors |ab| and |ac|
const DMVector3d DMVector3d::crossProduct(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c) 
{
	return (b-a).crossProduct(c-a); 
}

ostream& operator<<(ostream& os, const DMVector3d& v){
	string str_x, str_y, str_z;
	DEquation::doubleToString(v.x, DEquation::v_length, str_x);
	DEquation::doubleToString(v.y, DEquation::v_length, str_y);
	DEquation::doubleToString(v.z, DEquation::v_length, str_z);
	os << "<x>" << str_x << "</x> <y>" << str_y << "</y> <z>" << str_z  << "</z>";
	return os;
}

istream& operator>>(istream& is, DMVector3d& v){
	string str_x, str_y, str_z;
	string tmp;
	char z;
	getline(is, tmp, '>'); is >> z; // <x>
	getline(is, str_x, '<');
	getline(is, tmp, '>'); is >> z; // </x>
	getline(is, tmp, '>'); is >> z; // <y>
	getline(is, str_y, '<');
	getline(is, tmp, '>'); is >> z; // </y>
	getline(is, tmp, '>'); is >> z; // <z>
	getline(is, str_z, '<');
	getline(is, tmp, '>'); is >> z; // </z>
	DEquation eq;
	if(eq.parse(str_x.c_str()))	v.x = eq.getValue(0.0);
	if(eq.parse(str_y.c_str()))	v.y = eq.getValue(0.0);
	if(eq.parse(str_z.c_str()))	v.z = eq.getValue(0.0);
	return is;
}
