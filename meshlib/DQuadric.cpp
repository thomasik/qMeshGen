/////////////////////////////////////////////////////////////////////////////
// DQuadric.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DQuadric.h"
#include "DVectorN.h"
#include "DPoint.h"

/// Standard constructor
DQuadric::DQuadric() : vq(0.0) { }

DQuadric::DQuadric(const DVectorN<N> & v) : vq(v) { }

double DQuadric::distance(const DPoint3d& pt) const
{
	const double v[N] = { 1,  pt.x, pt.y, pt.z,
		pt.x * pt.x, pt.y * pt.y, pt.z * pt.z,
		pt.x * pt.y, pt.x * pt.z, pt.y * pt.z };

	double dist = 0.0;
	for(int i = 0; i < N; i++) dist += v[i] * vq[i];
	return dist;
}

/// solve Q( pt + t * vt) == 0, with possible two solutions (or some other cases)
bool DQuadric::solve(const DPoint3d& pt, const DVector3d& vt, double& t1, double &t2) const
{
	double a =  vq[4] * vt.x * vt.x + vq[7] * vt.x * vt.y + vq[5] * vt.y * vt.y + 
				vq[8] * vt.x * vt.z + vq[9] * vt.y * vt.z + vq[6] * vt.z * vt.z;
	double b =  vq[1] * vt.x + 2 * vq[4] * pt.x * vt.x + vq[7] * vt.x * pt.y + 
				vq[2] * vt.y + vq[7] * pt.x * vt.y + 2 * vq[5] * pt.y * vt.y + 
				vq[8] * vt.x * pt.z + vq[9] * vt.y * pt.z + vq[3] * vt.z + 
				vq[8] * pt.x * vt.z + vq[9] * pt.y * vt.z + 2 * vq[6] * pt.z * vt.z;
	double c =  vq[0] + vq[1] * pt.x + vq[4] * pt.x * pt.x + vq[2] * pt.y + 
				vq[7] * pt.x * pt.y + vq[5] * pt.y * pt.y + vq[3] * pt.z + 
				vq[8] * pt.x * pt.z + vq[9] * pt.y * pt.z + vq[6] * pt.z * pt.z;

	const double EPS = SMALL_NUMBER * std::max(std::max(abs(a), abs(b)), abs(c));

	if(abs(a) < EPS){ // b x + c = 0
		if(abs(b) < EPS) return false;
		t1 = t2 = -c / b;
		return true;
	}else{ 
		// quadratic
		double delta = b*b - 4*a*c;
		if(delta < -EPS) return false; // no real solutions
		if(delta < EPS){ // delta close to 0
			t1 = t2 = -b / (2*a);
		}else{
			delta = sqrt(delta);
			t1 = (-b-delta) / (2*a);
			t2 = (-b+delta) / (2*a);
		}
	}
	return true;
}
