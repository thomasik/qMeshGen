/////////////////////////////////////////////////////////////////////////////
// DOrientedBox.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2014-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "DOrientedBox.h"
#include "DEquation.h"
#include "MeshViewSet.h"
#include "Metric3dContext.h"

/// Standard constructor (creates empty - invalid box)
DOrientedBox::DOrientedBox() 
	: x0(0.0),x1(0.0),y0(0.0),y1(0.0),z0(0.0),z1(0.0),valid(false) 
{
	m_rot.setIdentity(); 
	m_rot_inv.setIdentity();
}

/// Standard constructor (creates empty - invalid box)
DOrientedBox::DOrientedBox(const DPoint3d& pt0, const DVector3d& e0, const DVector3d& e1) 
	: m_trans(pt0.x, pt0.y, pt0.z), m_rot(e0, e1, e0.crossProduct(e1) ),
		x0(0.0),x1(0.0),y0(0.0),y1(0.0),z0(0.0),z1(0.0),valid(false) 
{
	m_rot_inv = m_rot.inverse();
}

/// Standard constructor (creates valid box)
DOrientedBox::DOrientedBox(double _x0, double _x1, double _y0, double _y1, double _z0, double _z1) 
	: x0(_x0), x1(_x1),y0(_y0),y1(_y1),z0(_z0),z1(_z1),valid(true) 
{ 
	m_rot.setIdentity(); 
	m_rot_inv.setIdentity();
}

/////////////////////////////////////////////////////////////////////////////
// Metoda powoduje powiêkszenie danego prosy1ad³oœcianu tak, aby obejmowa³ zadany box.
void DOrientedBox::addOrientedBox(const DOrientedBox& box)
{
	assert(false); // to be updated with m_rot/m_trans
	if(!box.valid) return;
	if(!valid){
		x0	= box.x0;
		y0	= box.y0;
		z0	= box.z0;
		x1	= box.x1;
		y1	= box.y1;
		z1	= box.z1;
		valid	= box.valid;
	}else{
		x0	= std::min(x0,	box.x0);
		x1	= std::max(x1,	box.x1);
		y0	= std::min(y0,	box.y0);
		y1	= std::max(y1,	box.y1);
		z0	= std::min(z0,	box.z0);
		z1	= std::max(z1,	box.z1);
	}
}

/////////////////////////////////////////////////////////////////////////////
// Metoda powoduje powiêkszenie danego prostok¹ta tak, aby obejmowa³ zadany punkt.
void DOrientedBox::addPoint(const DPoint3d& point)
{
	const DPoint3d op = m_rot_inv * (point - m_trans);

	if(!valid){
		x0	= x1	= op.x;
		y0	= y1	= op.y;
		z1	= z0	= op.z;
		valid		= true;
	}else{
		x0	= std::min(x0,	op.x);
		x1	= std::max(x1,	op.x);
		y0	= std::min(y0,	op.y);
		y1	= std::max(y1,	op.y);
		z0	= std::min(z0,	op.z);
		z1	= std::max(z1,	op.z);
	}
}

void DOrientedBox::addPoint(const FPoint3d& point)
{
	addPoint( DPoint3d( point.x, point.y, point.z) );
}

/////////////////////////////////////////////////////////////////////////////
// Metoda powoduje powiêkszenie danego prostok¹ta tak, aby obejmowa³ zadany punkt.
void DOrientedBox::addPoint(const DMPoint3d& point)
{
	addPoint( DPoint3d( point.x, point.y, point.z) );
}


/// Checks, whether this rectangle contains the given 2D-point
bool DOrientedBox::contains(const DPoint3d& p) const 
{ 
	const DPoint3d op = m_rot_inv * (p - m_trans);
	return (op.x >= x0) && (op.x <= x1) && (op.y >= y0) && (op.y <= y1) && (op.z >= z0) && (op.z <= z1); 
}

/// Returns middle 3D-point of the box
DPoint3d DOrientedBox::getMiddlePoint() const 
{ 
	DPoint3d rel_p(0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1)); 
	return m_rot * rel_p + m_trans;
}

/// Increases size of the box by the given length (in metric, in all directions)
void DOrientedBox::growDM(Metric3dContext& mc, double mlen)
{
	
	growDX( mlen / mc.transformRStoMS( m_rot.column(0) ).length() );
	growDY( mlen / mc.transformRStoMS( m_rot.column(1) ).length() );
	growDZ( mlen / mc.transformRStoMS( m_rot.column(2) ).length() );
}

void DOrientedBox::inflate(double factor)
{
	double dx = factor * (x1 - x0);
	x0 -= dx;
	x1 += dx;

	double dz = factor * (z1 - z0);
	z0 -= dz;
	z1 += dz;

	double dy = factor * (y1 - y0);
	y0 -= dy;
	y1 += dy;
}

void DOrientedBox::grow(double offset)
{
	x0 -= offset;
	x1 += offset;
	z0 -= offset;
	z1 += offset;
	y0 -= offset;
	y1 += offset;
}

/////////////////////////////////////////////////////////////////////////////
// Method adjusts the coordinates of the given point, so it should be inside this box
DPoint3d DOrientedBox::fitInPoint(const DPoint3d& pt) const
{
	assert(false); // to be updated with m_rot/m_trans
	DPoint3d p = pt;

	if(p.x < x0) p.x = x0; else if(p.x > x1) p.x = x1;
	if(p.z < z0) p.z = z0; else if(p.z > z1) p.z = z1;
	if(p.y < y0) p.y = y0; else if(p.y > y1) p.y = y1;

	return p;
}

ostream& operator<<(ostream& os, const DOrientedBox& box){
	assert(false); // to be updated with m_rot/m_trans
	os << "[" << box.x0 << "," << box.x1 << ",";
	os <<        box.y0 << "," << box.y1 << ",";
	os <<        box.z0 << "," << box.z1 << "]";
	return os;
}

istream& operator>>(istream& is, DOrientedBox& box){
	assert(false); // to be updated with m_rot/m_trans
	string str_x0, str_y0, str_z0;
	string str_x1, str_y1, str_z1;
	char z;
	is >> ws >> z; // '['
	getline(is, str_x0, ',');
	getline(is, str_x1, ',');
	getline(is, str_y0, ',');
	getline(is, str_y1, ',');
	getline(is, str_z0, ',');
	getline(is, str_z1, ']');
	DEquation eq;
	if(eq.parse(str_x0.c_str()))	box.x0 = eq.getValue(0.0);
	if(eq.parse(str_x1.c_str()))	box.x1 = eq.getValue(0.0);
	if(eq.parse(str_y0.c_str()))	box.y0 = eq.getValue(0.0);
	if(eq.parse(str_y1.c_str()))	box.y1 = eq.getValue(0.0);
	if(eq.parse(str_z0.c_str()))	box.z0 = eq.getValue(0.0);
	if(eq.parse(str_z1.c_str()))	box.z1 = eq.getValue(0.0);
	return is;
}

MeshViewSet*  DOrientedBox::draw( MeshViewSet* set, int id) const
{
	if(!valid) return set;
	assert( set != nullptr );

	const DPoint3d p0 = m_rot * DPoint3d(x0, y0, z0) + m_trans;
	const DPoint3d p1 = m_rot * DPoint3d(x1, y0, z0) + m_trans;
	const DPoint3d p2 = m_rot * DPoint3d(x0, y1, z0) + m_trans;
	const DPoint3d p3 = m_rot * DPoint3d(x1, y1, z0) + m_trans;
	const DPoint3d p4 = m_rot * DPoint3d(x0, y0, z1) + m_trans;
	const DPoint3d p5 = m_rot * DPoint3d(x1, y0, z1) + m_trans;
	const DPoint3d p6 = m_rot * DPoint3d(x0, y1, z1) + m_trans;
	const DPoint3d p7 = m_rot * DPoint3d(x1, y1, z1) + m_trans;

	set->addEdge( p0, p1, id );
	set->addEdge( p2, p3, id );
	set->addEdge( p0, p2, id );
	set->addEdge( p1, p3, id );

	set->addEdge( p4, p5, id );
	set->addEdge( p6, p7, id );
	set->addEdge( p4, p6, id );
	set->addEdge( p5, p7, id );

	set->addEdge( p0, p4, id );
	set->addEdge( p1, p5, id );
	set->addEdge( p2, p6, id );
	set->addEdge( p3, p7, id );

	return set;
}
