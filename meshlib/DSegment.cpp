/////////////////////////////////////////////////////////////////////////////
// DSegment.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DSegment.h"
#include "GeometricPredicates.h"

/// Returns distance (squared) from segment to point
double DSegment2d::distanceToPoint2(const DPoint2d& p1, const DPoint2d& p2, const DPoint2d& pt)
{
	const DVector2d dv = p2 - p1;
	const DVector2d vt = pt - p1;
	double t = dv.scalarProduct(vt);

	if (t <= 0) return vt.length2();

	double d = dv.length2();
	if (t >= d) return (pt - p2).length2();

	// closest point is interior to segment
	return vt.length2() - t * t / d;
}

/// Returns distance from segment to point
double DSegment2d::distanceToPoint(const DPoint2d& p1, const DPoint2d& p2, const DPoint2d& pt)
{
	return sqrt(DSegment2d::distanceToPoint2(p1, p2, pt));
}

/// Returns squared distance between two segments
double DSegment3d::distance2ToSegment(
		const DPoint3d& spt0_a, const DPoint3d& spt0_b, 
		const DPoint3d& spt1_a, const DPoint3d& spt1_b)
{
	const DVector3d v0(spt0_a, spt0_b);
	const DVector3d v1(spt1_a, spt1_b);
	const DVector3d u(spt1_a, spt0_a);

	double a = v0.scalarProduct(v0);
	double b = v0.scalarProduct(v1);
	double c = v1.scalarProduct(v1);
	double d = v0.scalarProduct(u);
	double e = v1.scalarProduct(u);
	double det = a*c - b*b;

	double eps = SMALL_NUMBER * sqrt(std::max(a, c)); // relative to length of segments

	double sn = 0;
	double sd = det;
	double tn = e;
	double td = c;

	if(det > eps){ // if not (near) parallel - (else the default values as an arbitrary choice)
		// find parameter values of closest points on segment lines
		sn = b*e - c*d;
		tn = a*e - b*d;
	}

	// check s
	if(sn < 0.0){
		sn = 0.0;
		tn = e;
	}else if(sn > det){
		sn = det;
		tn = e+b;
	}else{
		td = det;
	}

	// check t
	if(tn < 0.0){
		tn = 0.0;
		if(-d < 0.0){
			sn = 0.0;
		}else if(-d > a){
			sn = sd;
		}else{
			sn = -d;
			sd = a;
		}
	}else if(tn > td){
		tn = td;
		if((-d+b) < 0.0){
			sn = 0.0;
		}else if((-d+b) > a){
			sn = sd;
		}else{
			sn = -d+b;
			sd = a;
		}
	}

	// parameters of nearest points
	double s = sn / sd;
	double t = tn / td;

	const DPoint3d pt_s = spt0_a + v0*s;
	const DPoint3d pt_t = spt1_a + v1*t;

	return pt_s.distance2(pt_t);
}

/// Returns squared distance between two segments
double DMSegment3d::distance2ToSegment(
		const DMPoint3d& spt0_a, const DMPoint3d& spt0_b, 
		const DMPoint3d& spt1_a, const DMPoint3d& spt1_b)
{
	const DMVector3d v0(spt0_a, spt0_b);
	const DMVector3d v1(spt1_a, spt1_b);
	const DMVector3d u(spt1_a, spt0_a);

	double a = v0.scalarProduct(v0);
	double b = v0.scalarProduct(v1);
	double c = v1.scalarProduct(v1);
	double d = v0.scalarProduct(u);
	double e = v1.scalarProduct(u);
	double det = a*c - b*b;

	double eps = SMALL_NUMBER * sqrt(std::max(a, c)); // relative to length of segments

	double sn = 0;
	double sd = det;
	double tn = e;
	double td = c;

	if(det > eps){ // if not (near) parallel - (else the default values as an arbitrary choice)
		// find parameter values of closest points on segment lines
		sn = b*e - c*d;
		tn = a*e - b*d;
	}

	// check s
	if(sn < 0.0){
		sn = 0.0;
		tn = e;
	}else if(sn > det){
		sn = det;
		tn = e+b;
	}else{
		td = det;
	}

	// check t
	if(tn < 0.0){
		tn = 0.0;
		if(-d < 0.0){
			sn = 0.0;
		}else if(-d > a){
			sn = sd;
		}else{
			sn = -d;
			sd = a;
		}
	}else if(tn > td){
		tn = td;
		if((-d+b) < 0.0){
			sn = 0.0;
		}else if((-d+b) > a){
			sn = sd;
		}else{
			sn = -d+b;
			sd = a;
		}
	}

	// parameters of nearest points
	double s = sn / sd;
	double t = tn / td;

	const DMPoint3d pt_s = spt0_a + v0*s;
	const DMPoint3d pt_t = spt1_a + v1*t;

	return pt_s.distance2(pt_t);
}


bool DSegment2d::crossingSegment(const DSegment2d& other) const
{
	double det1 = GeometricPredicates::orient2d(pt_a, pt_b, other.pt_a);
	double det2 = GeometricPredicates::orient2d(pt_a, pt_b, other.pt_b);
	if(det1 * det2 > 0.0) return false; // on same side
	det1 = GeometricPredicates::orient2d(other.pt_a, other.pt_b, pt_a);
	det2 = GeometricPredicates::orient2d(other.pt_a, other.pt_b, pt_b);
	return det1 * det2 <= 0.0;
}
