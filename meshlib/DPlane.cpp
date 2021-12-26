/////////////////////////////////////////////////////////////////////////////
// DPlane.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DPlane.h"
#include "DVector.h"
#include "DPoint.h"
#include "DMatrix.h"

DPlane::DPlane() {}

DPlane::DPlane(const DPoint3d& pt, const DVector3d& e0, const DVector3d& e1)
	: p0(pt), e0(e0), e1(e1)
{
	countNormal();
}

DPlane::DPlane(const DPoint3d& pt, const DVector3d& vn)
	: p0(pt), vn(vn)
{
	vn.orthonormalVectors(e0, e1);
}

void DPlane::countNormal()
{
	vn = e0.crossProduct(e1).normalized();
}

DPoint3d DPlane::projectToSpace(const DPoint2d& pt) const
{
	return p0 + e0*pt.x + e1*pt.y;
}

bool DPlane::switchOrientation()
{
	e0 = -e0;
	vn = -vn;
	return true;
}

DPoint2d DPlane::projectToPlane(const DPoint3d& point) const
{
	DMatrix3d A;
	A.setColumn(0, e0);
	A.setColumn(1, e1);
	A.setColumn(2, vn);
	DVector3d res;
	if(A.solve(point - p0, res)){
		return DPoint2d(res.x, res.y);
	}else{
		// old method
		double Wxy = e0.x * e1.y - e0.y * e1.x;
		double Wxz = e0.x * e1.z - e0.z * e1.x;
		double Wyz = e0.y * e1.z - e0.z * e1.y;
		double Ws, Wt;
		DPoint2d result;

		DVector3d point_diff = point - p0;

		if(abs(Wxy) > abs(Wxz)){
			if(abs(Wxy) > abs(Wyz)){
				//xy
				Ws = point_diff.x * e1.y - point_diff.y * e1.x;
				Wt = e0.x * point_diff.y - e0.y * point_diff.x;
				result.x = Ws / Wxy;
				result.y = Wt / Wxy;
			}else{
				//yz
				Ws = point_diff.y * e1.z - point_diff.z * e1.y;
				Wt = e0.y * point_diff.z - e0.z * point_diff.y;
				result.x = Ws / Wyz;
				result.y = Wt / Wyz;
			}
		}else{
			if(abs(Wxz) > abs(Wyz)){
				//xz
				Ws = point_diff.x * e1.z- point_diff.z* e1.x;
				Wt = e0.x * point_diff.z- e0.z* point_diff.x;
				result.x = Ws / Wxz;
				result.y = Wt / Wxz;
			}else{
				//yz
				Ws = point_diff.y * e1.z - point_diff.z * e1.y;
				Wt = e0.y * point_diff.z - e0.z * point_diff.y;
				result.x = Ws / Wyz;
				result.y = Wt / Wyz;
			}
		}
		return result;
	}
}

bool DPlane::projectToPlane(const DPoint3d& point, DPoint2d& param, double& z) const
{
	assert( !vn.isZero() );
	DMatrix3d A;
	A.setColumn(0, e0);
	A.setColumn(1, e1);
	A.setColumn(2, vn);
	DVector3d res;
	if(!A.solve(point - p0, res)) return false;
	param.x = res.x;
	param.y = res.y;
	z = res.z;
	return true;
}

ostream& operator<<(ostream& os, const DPlane & plane)
{
	return os << plane.p0 << " " << plane.e0 << " " << plane.e1;
}

istream& operator>>(istream& is, DPlane & plane)
{
	is >> ws >> plane.p0 >> ws >> plane.e0 >> ws >> plane.e1;
	plane.countNormal();

	return is;
}

DOrientedBox DPlane::getOrientedBox() const
{
	return DOrientedBox( p0, e0, e1 );
}

DOrientedBox DPlane::getOrientedBox( const DataVector<DPoint3d>& points ) const
{
	DOrientedBox obox( p0, e0, e1 );
	for(size_t i = 0; i < points.countInt(); i++) {
		obox.addPoint( points[i] );
	}
	return obox;
}

DOrientedBox DPlane::getOrientedBoxOpt( const DataVector<DPoint3d>& points ) const
{
	const DVector3d v1 = e0.normalized();
	const DVector3d v2 = e1.normalized();
	DOrientedBox obox( p0, (v2+v1)*0.5, (v2-v1)*0.5 );
	for(size_t i = 0; i < points.countInt(); i++) {
		obox.addPoint( points[i] );
	}
	return obox;
}
