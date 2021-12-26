/////////////////////////////////////////////////////////////////////////////
// DLine.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DLine.h"
#include "DMatrix.h"

/// Standard constructor
DLine2d::DLine2d(const DPoint2d& pt, const DVector2d& vt) 
	: m_pt(pt), m_vt(vt.normalized()), m_vn(vt.y, -vt.x) { }

/// Standard constructor
DLine2d::DLine2d(const DPoint2d& p1, const DPoint2d& p2) 
	: m_pt(p1), m_vt((p2-p1).normalized()), m_vn(m_vt.y, -m_vt.x)  {}

const DPoint2d DLine2d::crossPoint(const DPoint2d &p1, const DPoint2d &p2, const DPoint2d &p3, const DPoint2d &p4)
{
	return crossPoint(p1, p2-p1, p3, p4-p3);
}

/// set line attributes
void DLine2d::setPtVn(const DPoint2d& pt, const DVector2d& vn)
{
	m_pt = pt;
	m_vn = vn;
	m_vt = m_vn.orthogonal();
}

/// set line attributes
void DLine2d::setPtVt(const DPoint2d& pt, const DVector2d& vt)
{
	m_pt = pt;
	m_vt = vt;
	m_vn = m_vt.orthonormal();
}

const DPoint2d DLine2d::crossPoint(const DPoint2d& p1, const DVector2d& v1, const DPoint2d& p2, const DVector2d& v2)
{
	double delta = v2.x * v1.y - v1.x * v2.y;
	//if(abs(delta)<0.000000001) return false;
	assert(delta != 0.0);

	double t = (v2.x * (p2.y-p1.y) - (p2.x-p1.x) * v2.y) / delta;
	return p1 + v1*t;;
}

/// Returns distance (squared) from line to point
double DLine2d::distanceToPoint2(const DPoint2d& line_p, const DVector2d& line_v, const DPoint2d& pt)
{
	double d = distanceToPoint(line_p, line_v, pt);
	return d*d;
}

/// Returns distance from segment to point
double DLine2d::distanceToPoint(const DPoint2d& line_p, const DVector2d& line_v, const DPoint2d& pt)
{
	const DVector2d vn = DVector2d(line_v.y, -line_v.x).normalized();
	double c = vn.x * line_p.x + vn.y * line_p.y;

	return vn.x * pt.x + vn.y * pt.y - c;
}

/// Returns distance (squared) from line to point
double DLine2d::distanceToPoint2(const DPoint2d& _pt) const
{
	double d = distanceToPoint(_pt);
	return d*d;
}

/// Returns distance from segment to point
double DLine2d::distanceToPoint(const DPoint2d& pt) const
{
	return m_vn.x * (pt.x - m_pt.x) + m_vn.y * (pt.y - m_pt.y);
}

/// Returns point on the line for the given parameter
DPoint2d DLine2d::getPoint(double t) const
{
	return m_pt + m_vt * t;
}

/// Returns parameter for point on the line, nearest to pt
double DLine2d::paramForPoint(const DPoint2d& pt) const
{
	return m_vt.scalarProduct(pt - m_pt);
}

bool DLine2d::projectToLine(const DPoint2d& point, double& t, double& z) const
{
	assert( !m_vn.isZero() );
	DMatrix2d A;
	A.setColumn(0, m_vt);
	A.setColumn(1, m_vn);
	DVector2d res;
	if(!A.solve(point - m_pt, res)) return false;
	t = res.x;
	z = res.y;
	return true;
}

ostream& operator<<(ostream& os, const DLine2d & line)
{
	return os << line.m_pt << " " << line.m_vt;
}

istream& operator>>(istream& is, DLine2d & line)
{
	is >> ws >> line.m_pt >> ws >> line.m_vt;
	line.m_vn = line.m_vt.orthogonal();

	return is;
}

/// Standard constructor
DLine3d::DLine3d(const DPoint3d& _line_p, const DVector3d& _line_v) 
	: m_pt(_line_p), m_vt(_line_v) { assert( m_vt.isNormal() ); // should be normalised!
}

/// Standard constructor
DLine3d::DLine3d(const DPoint3d& p1, const DPoint3d& p2) 
	: m_pt(p1), m_vt((p2-p1).normalized()) {}

/// Returns distance (squared) from line to point
double DLine3d::distanceToPoint2(const DPoint3d& line_p, const DVector3d& line_v, const DPoint3d& pt, 
	bool line_v_normalized)
{
	double tq = line_v.scalarProduct(pt - line_p);
	if (line_v_normalized) {
		assert(line_v.isNormal());
	} else {
		tq /= line_v.scalarProduct(line_v);
	}
	return (pt - (line_p + line_v * tq)).length2();
}

/// Returns distance from segment to point
double DLine3d::distanceToPoint(const DPoint3d& line_p, const DVector3d& line_v, const DPoint3d& pt, 
	bool line_v_normalized)
{
	return sqrt(distanceToPoint2(line_p, line_v, pt, line_v_normalized));
}

/// Returns distance (squared) from line to point
double DLine3d::distanceToPoint2(const DPoint3d& pt) const
{
	return distanceToPoint2(m_pt, m_vt, pt, true);
}

/// Returns distance from segment to point
double DLine3d::distanceToPoint(const DPoint3d& pt) const
{
	return distanceToPoint(m_pt, m_vt, pt, true);
}

/// Returns parameter for point on the line, nearest to pt
double DLine3d::paramForPoint(const DPoint3d& pt) const
{
	return m_vt.scalarProduct(pt - m_pt);
}

/// Returns point on the line for the given parameter
DPoint3d DLine3d::getPoint(double t) const
{
	return m_pt + m_vt * t;
}

/// set line attributes
void DLine3d::setPtVt(const DPoint3d& pt, const DVector3d& vt)
{
	m_pt = pt;
	assert( vt.isNormal() ); // should be normalised!
	m_vt = vt;
}

ostream& operator<<(ostream& os, const DLine3d & line)
{
	return os << line.m_pt << " " << line.m_vt;
}

istream& operator>>(istream& is, DLine3d & line)
{
	is >> ws >> line.m_pt >> ws >> line.m_vt;

	return is;
}

