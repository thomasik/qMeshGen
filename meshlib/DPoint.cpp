/////////////////////////////////////////////////////////////////////////////
// DPoint.cpp
// This class implements a point in a two-dimensional space
//  plus some basic operations.
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include <cfloat>

#include "common.h"

#include "DPoint.h"
#include "DEquation.h"

const DPoint2d DPoint2d::zero(0.0, 0.0);
const DPoint2d DPoint2d::inf(DBL_MAX, DBL_MAX);

/// clear edges with "zero" length (less than threshold)
int DPoint2d::clearZeroEdges(DataVector<DPoint2d> & points, double threshold)
{
	size_t pct = points.countInt();
	double T2 = threshold*threshold;
	for(size_t i = 0; (pct > 2) && (i < pct); ){
		double dist2 = points[i].distance2(points[(i+1)%pct]);
		if(dist2 < T2){
			points.removeOrderedAt(i);
			--pct;
		}else
			++i;
	}
	return (int)pct;
}

/// clear edges with "zero" length (less than threshold)
int DPoint2d::clearZeroEdges(DataVector<DPoint2d> & points, DataVector<int> &trans_tab, double threshold)
{
	size_t pct = points.countInt();
	trans_tab.clear();
	trans_tab.prepare(pct);
	for(size_t i = 0; i < pct; i++) trans_tab.add((int)i);

	double T2 = threshold*threshold;
	for(size_t i = 0; (pct > 2) && (i < pct); ){
		double dist2 = points[i].distance2(points[(i+1)%pct]);
		if(dist2 < T2){
			points.removeOrderedAt(i);
			trans_tab.removeOrderedAt(i);
			--pct;
		}else
			++i;
	}
	return (int)pct;
}

// Sum over the edges, (x2-x1)(y2+y1). 
// If the result is positive the curve is clockwise, if it's negative the curve is counter-clockwise. 
// (The result is twice the enclosed area, with a +/- convention.)
bool DPoint2d::properOrientation(const DataVector<DPoint2d> & points)
{
	size_t pct = points.countInt();
	if(pct < 3) return false;

	double d = 0.0;

	for(size_t i = 0; i < pct; i++){
		const DPoint2d& p0 = points[i];
		const DPoint2d& p1 = points[(i + 1) % pct];

		d += (p1.x-p0.x)*(p1.y+p0.y);
	}
	return d < 0.0;
}

bool DPoint2d::storeSimple(ostream& os, char end_char) const {
	os << x << '\t' << y << end_char;
	return os.good();
}

bool DPoint2d::readSimple(istream& is) {
	is >> x >> y;
	return is.good();
}

ostream& operator<<(ostream& os, const DPoint2d& pt){
	string str_x, str_y;
	DEquation::doubleToString(pt.x, DEquation::v_length, str_x);
	DEquation::doubleToString(pt.y, DEquation::v_length, str_y);
	os << "<u>" << str_x << "</u> <v>" << str_y << "</v>";
	return os;
}

istream& operator>>(istream& is, DPoint2d& pt){
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
	if(eq.parse(str_x.c_str()))	pt.x = eq.getValue(0.0);
	if(eq.parse(str_y.c_str()))	pt.y = eq.getValue(0.0);
	return is;
}

const DMPoint2d DMPoint2d::zero(0.0, 0.0);

ostream& operator<<(ostream& os, const DMPoint2d& pt){
	string str_x, str_y;
	DEquation::doubleToString(pt.x, DEquation::v_length, str_x);
	DEquation::doubleToString(pt.y, DEquation::v_length, str_y);
	os << "<u>" << str_x << "</u> <v>" << str_y << "</v>";
	return os;
}

istream& operator>>(istream& is, DMPoint2d& pt){
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
	if(eq.parse(str_x.c_str()))	pt.x = eq.getValue(0.0);
	if(eq.parse(str_y.c_str()))	pt.y = eq.getValue(0.0);
	return is;
}

const DPoint3d DPoint3d::zero(0.0, 0.0, 0.0);
const DPoint3d DPoint3d::inf(DBL_MAX, DBL_MAX, DBL_MAX);
const FPoint3d FPoint3d::zero(0.0, 0.0, 0.0);

bool DPoint3d::storeSimple(ostream& os, char end_char) const {
	os << x << '\t' << y << '\t' << z << end_char;
	return os.good();
}
bool DPoint3d::readSimple(istream& is) {
	is >> x >> y >> z;
	return is.good();
}

ostream& operator<<(ostream& os, const DPoint3d& pt){
	string str_x, str_y, str_z;
	DEquation::doubleToString(pt.x, DEquation::v_length, str_x);
	DEquation::doubleToString(pt.y, DEquation::v_length, str_y);
	DEquation::doubleToString(pt.z, DEquation::v_length, str_z);
	os << "<x>" << str_x << "</x> <y>" << str_y << "</y> <z>" << str_z  << "</z>";
	return os;
}

istream& operator>>(istream& is, DPoint3d& pt){
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
	if(eq.parse(str_x.c_str()))	pt.x = eq.getValue(0.0);
	if(eq.parse(str_y.c_str()))	pt.y = eq.getValue(0.0);
	if(eq.parse(str_z.c_str()))	pt.z = eq.getValue(0.0);
	return is;
}

const DMPoint3d DMPoint3d::zero(0.0, 0.0, 0.0);

OctVertexWhich DPoint3d::vertexDirection(const DPoint3d & pt) const
{
	if (pt.x < x)
		if (pt.y < y)
			return (pt.z < z) ? VX3D_LSW : VX3D_HSW;
		else
			return (pt.z < z) ? VX3D_LNW : VX3D_HNW;
	else if (pt.y < y)
		return (pt.z < z) ? VX3D_LSE : VX3D_HSE;
	else
		return (pt.z < z) ? VX3D_LNE : VX3D_HNE;
}
