/////////////////////////////////////////////////////////////////////////////
// DRect.cpp
// Klasa implementuj¹ca prostok¹t okreœlony w uk³adzie wspó³rzêdnych
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////

#include "DRect.h"
#include "DEquation.h"
#include "MeshViewSet.h"
#include "Metric3dContext.h"

/////////////////////////////////////////////////////////////////////////////
// Metoda powoduje powiêkszenie danego prostok¹ta tak, aby obejmowa³ zadany prostok¹t rct.
void DRect::addRect(const DRect & rct)
{
	if(!rct.valid) return;
	if(!valid){
		x0	= rct.x0;
		x1	= rct.x1;
		y0	= rct.y0;
		y1	= rct.y1;
		valid	= rct.valid;
	}else{
		x0	= std::min(x0,	rct.x0);
		x1	= std::max(x1,	rct.x1);
		y0	= std::min(y0,	rct.y0);
		y1	= std::max(y1,	rct.y1);
	}
}

/////////////////////////////////////////////////////////////////////////////
// Metoda powoduje powiêkszenie danego prostok¹ta tak, aby obejmowa³ zadany punkt.
void DRect::addPoint(const DPoint2d& point)
{
	if(!valid){
		x0 = x1 = point.x;
		y1	= y0 = point.y;
		valid = true;
	}else{
		x0 = std::min(x0,	point.x);
		y1		= std::max(y1,		point.y);
		x1	= std::max(x1,	point.x);
		y0	= std::min(y0,	point.y);
	}
}

void DRect::inflate(double factor)
{
	double dx = factor * (x1-x0);
	x0 -= dx;
	x1 += dx;

	double dy = factor * (y1-y0);
	y0 -= dy;
	y1 += dy;
}

DRect DRect::inflated(double factor) const
{
	double dx = factor * (x1-x0);
	double dy = factor * (y1-y0);
	return DRect(x0-dx, y1+dy, x1+dx, y0-dy);
}

/////////////////////////////////////////////////////////////////////////////
// Method adjusts the coordinates of the given point, so it should be inside this rectangle
DPoint2d DRect::fitInPoint(const DPoint2d& pt) const
{
	DPoint2d p = pt;

	if(p.x < x0)	p.x = x0;	else if(p.x > x1) p.x = x1;
	if(p.y < y0) p.y = y0;	else if(p.y > y1) p.y = y1;

	return p;
}

/// Calculates two crossing points for rectangle and line 
bool DRect::calculateCrossingPoints(const DPoint2d& line_pt, const DVector2d& line_vt, DPoint2d& pt0, DPoint2d& pt1) const
{
	if(std::abs(line_vt.x) < SMALL_NUMBER){
		// line OY
		assert(line_pt.x >= x0 && line_pt.x <= x1);
		if(line_pt.x < x0 || line_pt.x > x1) return false;
		pt0.x = pt1.x = line_pt.x;
		pt0.y = y0;
		pt1.y = y1;
	}else if(std::abs(line_vt.y) < SMALL_NUMBER){
		// line OX
		assert(line_pt.y >= y0 && line_pt.y <= y1); 
		if(line_pt.y < y0 || line_pt.y > y1) return false; // outside the box
		pt0.y = pt1.y = line_pt.y;
		pt0.x = x0;
		pt1.x = x1;
	}else{
		// x0, x1, y0, y1
		double t[4] = { (x0   - line_pt.x) / line_vt.x, 
						(x1  - line_pt.x) / line_vt.x,
						(y0 - line_pt.y) / line_vt.y,
						(y1    - line_pt.y) / line_vt.y};
		// sort
		int ind[4] = {0, 1, 2, 3};
		for(int i = 0; i < 3; i++){
			// find min
			int i_min = i;
			for(int j = i+1; j < 4; j++)
				if(t[ind[j]] < t[ind[i_min]]) i_min = j;
			// switch
			int tmp = ind[i]; ind[i] = ind[i_min]; ind[i_min] = tmp;
		}
		// check crossing conditions
		if((ind[0] + ind[1] == 1) ||	// both two first vertical
			(ind[0] + ind[1] == 5))		// both two first horizontal
			return false;
		// select two inner t-values
		pt0 = line_pt + line_vt * t[ind[1]];
		pt1 = line_pt + line_vt * t[ind[2]];
	}
	return true;
}

void DRect::splitAndUpdate(Axis axis, double value, bool lower)
{
	switch (axis) {
	case Axis::X:
		assert(value > x0 && value < x1);
		if (lower) x1 = value; else x0 = value;
		break;
	case Axis::Y:
		assert(value > y0 && value < y1);
		if (lower) y1 = value; else y0 = value;
		break;
	case Axis::Z:
		assert(false);
		break;
	}
}

void DRect::split(Axis axis, double value, DRect & box0, DRect & box1) const
{
	switch (axis) {
	case Axis::X:
		assert(value > x0 && value < x1);
		box0.set(x0, value, y0, y1);
		box1.set(value, x1, y0, y1);
		break;
	case Axis::Y:
		assert(value > y0 && value < y1);
		box0.set(x0, x1, y0, value);
		box1.set(x0, x1, value, y1);
		break;
	case Axis::Z:
		assert(false);
		break;
	}
}

Axis DRect::getLongestAxis() const
{
	return ((x1 - x0) > (y1 - y0)) ? Axis::X : Axis::Y;
}

double DRect::getMiddleValue(Axis axis) const
{
	switch (axis) {
	case Axis::X: return (x0 + x1)*0.5;
	case Axis::Y: return (y0 + y1)*0.5;
	default: assert(false); return 0.0;
	}
}

double DRect::getDLen(Axis axis) const
{
	switch (axis) {
	case Axis::X: return x1 - x0;
	case Axis::Y: return y1 - y0;
	default: assert(false); return 0.0;
	}
}

double DRect::getD0(Axis axis) const
{
	switch (axis) {
	case Axis::X: return x0;
	case Axis::Y: return y0;
	default: assert(false); return 0.0;
	}
}

double DRect::getD1(Axis axis) const
{
	switch (axis) {
	case Axis::X: return x1;
	case Axis::Y: return y1;
	default: assert(false); return 0.0;
	}
}

OctVertexWhich DBox::edge_to_vertex[12][2] = {
	{ VX3D_LSW, VX3D_LSE }, // ED_LS
	{ VX3D_LNW, VX3D_LNE }, // ED_LN
	{ VX3D_LSW, VX3D_LNW }, // ED_LW
	{ VX3D_LSE, VX3D_LNE }, // ED_LE
	{ VX3D_LSW, VX3D_HSW }, // ED_SW
	{ VX3D_LSE, VX3D_HSE }, // ED_SE
	{ VX3D_LNW, VX3D_HNW }, // ED_NW
	{ VX3D_LNE, VX3D_HNE }, // ED_NE
	{ VX3D_HSW, VX3D_HSE }, // ED_HS
	{ VX3D_HNW, VX3D_HNE }, // ED_HN
	{ VX3D_HSW, VX3D_HNW }, // ED_HW
	{ VX3D_HSE, VX3D_HNE }  // ED_HE
};

OctFaceWhich DBox::edge_to_face[12][2] = {
	{ FC3D_LOW, FC3D_SOUTH }, // ED3D_LS
	{ FC3D_LOW, FC3D_NORTH }, // ED3D_LN
	{ FC3D_LOW, FC3D_WEST }, // ED3D_LW
	{ FC3D_LOW, FC3D_EAST }, // ED3D_LE
	{ FC3D_SOUTH, FC3D_WEST }, // ED3D_SW
	{ FC3D_SOUTH, FC3D_EAST }, // ED3D_SE
	{ FC3D_NORTH, FC3D_WEST }, // ED3D_NW
	{ FC3D_NORTH, FC3D_EAST }, // ED3D_NE
	{ FC3D_HIGH, FC3D_SOUTH }, // ED3D_HS
	{ FC3D_HIGH, FC3D_NORTH }, // ED3D_HN
	{ FC3D_HIGH, FC3D_WEST }, // ED3D_HW
	{ FC3D_HIGH, FC3D_EAST }  // ED3D_HE
};

OctFaceWhich DBox::vertex_to_face[8][3] = {
	{ FC3D_LOW, FC3D_SOUTH, FC3D_WEST }, // VX3D_LSW
	{ FC3D_LOW, FC3D_SOUTH, FC3D_EAST }, // VX3D_LSE
	{ FC3D_LOW, FC3D_NORTH, FC3D_WEST }, // VX3D_LNW
	{ FC3D_LOW, FC3D_NORTH, FC3D_EAST }, // VX3D_LNE
	{ FC3D_HIGH, FC3D_SOUTH, FC3D_WEST }, // VX3D_HSW
	{ FC3D_HIGH, FC3D_SOUTH, FC3D_EAST }, // VX3D_HSE
	{ FC3D_HIGH, FC3D_NORTH, FC3D_WEST }, // VX3D_HNW
	{ FC3D_HIGH, FC3D_NORTH, FC3D_EAST }  // VX3D_HNE
};

OctEdgeWhich DBox::vertex_to_edge[8][3] = {
	{ ED3D_LS, ED3D_LW, ED3D_SW }, // VX3D_LSW
	{ ED3D_LS, ED3D_LE, ED3D_SE }, // VX3D_LSE
	{ ED3D_LN, ED3D_LW, ED3D_NW }, // VX3D_LNW
	{ ED3D_LN, ED3D_LE, ED3D_NE }, // VX3D_LNE
	{ ED3D_HS, ED3D_HW, ED3D_SW }, // VX3D_HSW
	{ ED3D_HS, ED3D_HE, ED3D_SE }, // VX3D_HSE
	{ ED3D_HN, ED3D_HW, ED3D_NW }, // VX3D_HNW
	{ ED3D_HN, ED3D_HE, ED3D_NE }  // VX3D_HNE
};

OctVertexWhich DBox::face_to_vertex[6][4] = { // for each side of leaf -> 4 sub-leaves
	{ VX3D_LSW, VX3D_LSE, VX3D_HSW, VX3D_HSE }, // FC3D_SOUTH
	{ VX3D_LNW, VX3D_LNE, VX3D_HNW, VX3D_HNE }, // FC3D_NORTH
	{ VX3D_LSW, VX3D_LNW, VX3D_HSW, VX3D_HNW }, // FC3D_WEST
	{ VX3D_LSE, VX3D_LNE, VX3D_HSE, VX3D_HNE }, // FC3D_EAST
	{ VX3D_LSW, VX3D_LNW, VX3D_LSE, VX3D_LNE }, // FC3D_LOW
	{ VX3D_HSW, VX3D_HNW, VX3D_HSE, VX3D_HNE }  // FC3D_HIGH
};

OctEdgeWhich DBox::face_to_edge[6][4] = { // for each side of leaf -> 4 edges
	{ ED3D_LS, ED3D_SE, ED3D_HS, ED3D_SW }, // FC3D_SOUTH
	{ ED3D_LN, ED3D_NE, ED3D_HN, ED3D_NW }, // FC3D_NORTH
	{ ED3D_LW, ED3D_SW, ED3D_HW, ED3D_NW }, // FC3D_WEST
	{ ED3D_LE, ED3D_SE, ED3D_HE, ED3D_NE }, // FC3D_EAST
	{ ED3D_LS, ED3D_LE, ED3D_LN, ED3D_LW }, // FC3D_LOW
	{ ED3D_HS, ED3D_HE, ED3D_HN, ED3D_HW }  // FC3D_HIGH
};

OctVertexWhich DBox::opposite_vertex[8] = {
	VX3D_HNE, VX3D_HNW, VX3D_HSE, VX3D_HNW,
	VX3D_LNE, VX3D_LNW, VX3D_LSE, VX3D_LNW
};

/////////////////////////////////////////////////////////////////////////////
// Metoda powoduje powiêkszenie danego prosy1ad³oœcianu tak, aby obejmowa³ zadany prosy1ad³oœcian box.
void DBox::addBox(const DBox& box)
{
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
void DBox::addPoint(const DPoint3d& point)
{
	if(!valid){
		x0	= x1	= point.x;
		y0	= y1	= point.y;
		z1	= z0	= point.z;
		valid		= true;
	}else{
		x0	= std::min(x0,	point.x);
		x1	= std::max(x1,	point.x);
		y0	= std::min(y0,	point.y);
		y1	= std::max(y1,	point.y);
		z0	= std::min(z0,	point.z);
		z1	= std::max(z1,	point.z);
	}
}

void DBox::addPoint(const FPoint3d& point)
{
	addPoint( DPoint3d( point.x, point.y, point.z) );
}

/////////////////////////////////////////////////////////////////////////////
// Metoda powoduje powiêkszenie danego prostok¹ta tak, aby obejmowa³ zadany punkt.
void DBox::addPoint(const DMPoint3d& point)
{
	if(!valid){
		x0	= x1	= point.x;
		y0	= y1	= point.y;
		z1	= z0	= point.z;
		valid		= true;
	}else{
		x0	= std::min(x0,	point.x);
		x1	= std::max(x1,	point.x);
		y0	= std::min(y0,	point.y);
		y1	= std::max(y1,	point.y);
		z0	= std::min(z0,	point.z);
		z1	= std::max(z1,	point.z);
	}
}

void DBox::inflate(double factor)
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

void DBox::grow(double offset)
{
	x0 -= offset;
	x1 += offset;
	z0 -= offset;
	z1 += offset;
	y0 -= offset;
	y1 += offset;
}

/// Increases size of the box by the given length (in metric, in all directions)
void DBox::growDM(Metric3dContext& mc, double mlen)
{
	growDX( mlen / mc.transformRStoMS( DVector3d(1.0, 0.0, 0.0) ).length() );
	growDY( mlen / mc.transformRStoMS( DVector3d(0.0, 1.0, 0.0) ).length() );
	growDZ( mlen / mc.transformRStoMS( DVector3d(0.0, 0.0, 1.0) ).length() );
}

/////////////////////////////////////////////////////////////////////////////
// Method adjusts the coordinates of the given point, so it should be inside this box
DPoint3d DBox::fitInPoint(const DPoint3d& pt) const
{
	DPoint3d p = pt;

	if(p.x < x0) p.x = x0; else if(p.x > x1) p.x = x1;
	if(p.z < z0) p.z = z0; else if(p.z > z1) p.z = z1;
	if(p.y < y0) p.y = y0; else if(p.y > y1) p.y = y1;

	return p;
}

DPoint3d DBox::fitInPoint(const DPoint3d& pt, const DVector3d& dv) const
{
	DPoint3d cpt = pt + dv;
	if (!contains(cpt)) {
		// fit
		const DPoint3d fpt = fitInPoint(cpt);
		double tx = (abs(cpt.x - pt.x) > mesh_data.relative_small_number) ? ((fpt.x - pt.x) / (cpt.x - pt.x)) : 1.0;
		double ty = (abs(cpt.y - pt.y) > mesh_data.relative_small_number) ? ((fpt.y - pt.y) / (cpt.y - pt.y)) : 1.0;
		double tz = (abs(cpt.z - pt.z) > mesh_data.relative_small_number) ? ((fpt.z - pt.z) / (cpt.z - pt.z)) : 1.0;
		cpt = fitInPoint(pt + dv * std::min(std::min(tx, ty), tz));
	}
	return cpt;
}

Axis DBox::getLongestAxis() const
{
	double dx = x1 - x0;
	double dy = y1 - y0;
	double dz = z1 - z0;
	if (dx > dy && dx > dz) return Axis::X;
	else return (dy > dz) ? Axis::Y : Axis::Z;
}

double DBox::getMiddleValue(Axis axis) const
{
	switch (axis) {
	case Axis::X: return (x0 + x1)*0.5;
	case Axis::Y: return (y0 + y1)*0.5;
	case Axis::Z: return (z0 + z1)*0.5;
	}
	assert(false);
	return 0.0;
}

DPoint3d DBox::getRandomPoint() const
{
	static RandomGen rg(1);
	return DPoint3d(
		x0 + rg.doub() * (x1-x0),
		y0 + rg.doub() * (y1-y0),
		z0 + rg.doub() * (z1-z0)
		);
}

DPoint3d DBox::getRandomPoint(double margin) const
{
	static RandomGen rg(1);
	double mm = 2 * margin;
	double dx = x1 - x0;
	double dy = y1 - y0;
	double dz = z1 - z0;
	return DPoint3d(
		(dx > mm) ? (x0 + margin + rg.doub() * (dx - mm)) : (x0 + 0.5 * dx),
		(dy > mm) ? (y0 + margin + rg.doub() * (dy - mm)) : (y0 + 0.5 * dy),
		(dz > mm) ? (z0 + margin + rg.doub() * (dz - mm)) : (z0 + 0.5 * dz)
	);
}

DPoint3d DBox::getRandomPoint(double margin_x, double margin_y, double margin_z) const
{
	static RandomGen rg(1);
	double mmx = 2 * margin_x;
	double mmy = 2 * margin_y;
	double mmz = 2 * margin_z;
	double dx = x1 - x0;
	double dy = y1 - y0;
	double dz = z1 - z0;
	return DPoint3d(
		(dx > mmx) ? (x0 + margin_x + rg.doub() * (dx - mmx)) : (x0 + 0.5 * dx),
		(dy > mmy) ? (y0 + margin_y + rg.doub() * (dy - mmy)) : (y0 + 0.5 * dy),
		(dz > mmz) ? (z0 + margin_z + rg.doub() * (dz - mmz)) : (z0 + 0.5 * dz)
	);
}

void DBox::split(Axis axis, double value, DBox & box0, DBox & box1) const
{
	switch (axis) {
	case Axis::X:
		//assert(value > x0 && value < x1);
		box0.set(x0, value, y0, y1, z0, z1);
		box1.set(value, x1, y0, y1, z0, z1);
		break;
	case Axis::Y:
		//assert(value > y0 && value < y1);
		box0.set(x0, x1, y0, value, z0, z1);
		box1.set(x0, x1, value, y1, z0, z1);
		break;
	case Axis::Z:
		//assert(value > z0 && value < z1);
		box0.set(x0, x1, y0, y1, z0, value);
		box1.set(x0, x1, y0, y1, value, z1);
		break;
	}
}

void DBox::splitAndUpdate(Axis axis, double value, bool lower)
{
	switch (axis) {
	case Axis::X:
		//assert(value > x0 && value < x1);
		if (lower) x1 = value; else x0 = value;
		break;
	case Axis::Y:
		//assert(value > y0 && value < y1);
		if (lower) y1 = value; else y0 = value;
		break;
	case Axis::Z:
		//assert(value > z0 && value < z1);
		if (lower) z1 = value; else z0 = value;
		break;
	}
}

DBox DBox::splitOct(const DPoint3d & spt, int v) const
{
	switch (v) {
	case VX3D_LSW:
		return DBox(x0, spt.x, y0, spt.y, z0, spt.z);
	case VX3D_LSE:
		return DBox(spt.x, x1, y0, spt.y, z0, spt.z);
	case VX3D_LNW:
		return DBox(x0, spt.x, spt.y, y1, z0, spt.z);
	case VX3D_LNE:
		return DBox(spt.x, x1, spt.y, y1, z0, spt.z);
	case VX3D_HSW:
		return DBox(x0, spt.x, y0, spt.y, spt.z, z1);
	case VX3D_HSE:
		return DBox(spt.x, x1, y0, spt.y, spt.z, z1);
	case VX3D_HNW:
		return DBox(x0, spt.x, spt.y, y1, spt.z, z1);
	case VX3D_HNE:
		return DBox(spt.x, x1, spt.y, y1, spt.z, z1);
	default:
		assert(false);
	}
	return DBox();
}

void DBox::splitOctAndUpdate(const DPoint3d & spt, int v)
{
	switch (v) {
	case VX3D_LSW:
		x1 = spt.x, y1 = spt.y, z1 = spt.z; break;
	case VX3D_LSE:
		x0 = spt.x, y1 = spt.y, z1 = spt.z; break;
	case VX3D_LNW:
		x1 = spt.x, y0 = spt.y, z1 = spt.z; break;
	case VX3D_LNE:
		x0 = spt.x, y0 = spt.y, z1 = spt.z; break;
	case VX3D_HSW:
		x1 = spt.x, y1 = spt.y, z0 = spt.z; break;
	case VX3D_HSE:
		x0 = spt.x, y1 = spt.y, z0 = spt.z; break;
	case VX3D_HNW:
		x1 = spt.x, y0 = spt.y, z0 = spt.z; break;
	case VX3D_HNE:
		x0 = spt.x, y0 = spt.y, z0 = spt.z; break;
	default:
		assert(false);
	}
}

DPoint3d DBox::getVertex(int v) const
{
	switch (v) {
	case VX3D_LSW:	return DPoint3d(x0, y0, z0);
	case VX3D_LSE:	return DPoint3d(x1, y0, z0);
	case VX3D_LNW:	return DPoint3d(x0, y1, z0);
	case VX3D_LNE:	return DPoint3d(x1, y1, z0);
	case VX3D_HSW:	return DPoint3d(x0, y0, z1);
	case VX3D_HSE:	return DPoint3d(x1, y0, z1);
	case VX3D_HNW:	return DPoint3d(x0, y1, z1);
	case VX3D_HNE:	return DPoint3d(x1, y1, z1);
	default:
		assert(false);
	}
	return DPoint3d::zero;
}

DPoint3d DBox::getMiddlePointForFace(int fi) const
{
	switch (fi) {
	case FC3D_SOUTH:	return DPoint3d((x0 + x1) / 2, y0, (z0 + z1) / 2);
	case FC3D_NORTH:	return DPoint3d((x0 + x1) / 2, y1, (z0 + z1) / 2);
	case FC3D_WEST:		return DPoint3d(x0, (y0 + y1) / 2, (z0 + z1) / 2);
	case FC3D_EAST:		return DPoint3d(x1, (y0 + y1) / 2, (z0 + z1) / 2);
	case FC3D_LOW:		return DPoint3d((x0 + x1) / 2, (y0 + y1) / 2, z0);
	case FC3D_HIGH:		return DPoint3d((x0 + x1) / 2, (y0 + y1) / 2, z1);
	default:
		assert(false);
	}
	return DPoint3d::zero;
}

DPoint3d DBox::getMiddlePointForEdge(int ei) const
{
	switch (ei) {
	case ED3D_LS:	return DPoint3d((x0 + x1) / 2, y0, z0);
	case ED3D_LN:	return DPoint3d((x0 + x1) / 2, y1, z0);
	case ED3D_LW:	return DPoint3d(x0, (y0 + y1) / 2, z0);
	case ED3D_LE:	return DPoint3d(x1, (y0 + y1) / 2, z0);
	case ED3D_SW:	return DPoint3d(x0, y0, (z0 + z1) / 2);
	case ED3D_SE:	return DPoint3d(x1, y0, (z0 + z1) / 2);
	case ED3D_NW:	return DPoint3d(x0, y1, (z0 + z1) / 2);
	case ED3D_NE:	return DPoint3d(x1, y1, (z0 + z1) / 2);
	case ED3D_HS:	return DPoint3d((x0 + x1) / 2, y0, z1);
	case ED3D_HN:	return DPoint3d((x0 + x1) / 2, y1, z1);
	case ED3D_HW:	return DPoint3d(x0, (y0 + y1) / 2, z1);
	case ED3D_HE:	return DPoint3d(x1, (y0 + y1) / 2, z1);
	default:
		assert(false);
	}
	return DPoint3d::zero;
}

double DBox::getDLen(Axis axis) const
{
	switch (axis) {
	case Axis::X: return x1 - x0;
	case Axis::Y: return y1 - y0;
	case Axis::Z: return z1 - z0;
	default: assert(false); return 0.0;
	}
}

double DBox::getD0(Axis axis) const
{
	switch (axis) {
	case Axis::X: return x0;
	case Axis::Y: return y0;
	case Axis::Z: return z0;
	default: assert(false); return 0.0;
	}
}

double DBox::getD1(Axis axis) const
{
	switch (axis) {
	case Axis::X: return x1;
	case Axis::Y: return y1;
	case Axis::Z: return z1;
	default: assert(false); return 0.0;
	}
}

void DBox::forRegularGrid(int grid_size, const std::function<void(const DPoint3d&pt)>& fg) const {
	double fx = 1.0 / (double)grid_size;
	double dx = getDX() * fx;
	double dy = getDY() * fx;
	double dz = getDZ() * fx;

	DPoint3d pt(x0 + 0.5*dx, y0, z0);
	for (; pt.x < x1; pt.x += dx) {
		for (pt.y = y0 + 0.5*dy; pt.y < y1; pt.y += dy) {
			for (pt.z = z0 + 0.5*dz; pt.z < z1; pt.z += dz) {
				fg(pt);
			}
		}
	}
}

ostream& operator<<(ostream& os, const DBox& box){
	os << "[" << box.x0 << "," << box.x1 << ",";
	os <<        box.y0 << "," << box.y1 << ",";
	os <<        box.z0 << "," << box.z1 << "]";
	return os;
}

istream& operator>>(istream& is, DBox& box){
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

MeshViewSet*  DBox::draw( MeshViewSet* set, int id) const
{
	if(!valid) return set;
	assert( set != nullptr );

	const DPoint3d p0(x0, y0, z0);
	const DPoint3d p1(x1, y0, z0);
	const DPoint3d p2(x0, y1, z0);
	const DPoint3d p3(x1, y1, z0);
	const DPoint3d p4(x0, y0, z1);
	const DPoint3d p5(x1, y0, z1);
	const DPoint3d p6(x0, y1, z1);
	const DPoint3d p7(x1, y1, z1);

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

MeshViewSet*  DBox::drawXYZ(MeshViewSet* set, int id) const
{
	if (!valid) return set;
	assert(set != nullptr);
	draw(set, id);

	set->addLabel(DPoint3d(x0, y0, z0), "0");
	set->addLabel(DPoint3d(x1, y0, z0), "x");
	set->addLabel(DPoint3d(x0, y1, z0), "y");
	set->addLabel(DPoint3d(x0, y0, z1), "z");

	return set;
}