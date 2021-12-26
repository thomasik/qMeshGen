/////////////////////////////////////////////////////////////////////////////
// DQuad.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "DQuad.h"
#include "common.h"
#include "DTriangle.h"

/**
 * @param b is a reference to the second point of the quad,
 * @param c is a refenence to the third point of the quad,
 * @param d is a refenence to the fourth point of the quad,
 *
 * The coefficient is equal to:
 *  - 1 for equiangular quads,
 *  - 0 for degenerated quads, 
 *  - less than 0 for inverted quads(ie. oriented clockwise or skewed)
 *
 * \remark This is a "2-dimensional" version - it doesn't use any surface curvature information!
 */
double DQuad2d::alphaQuality(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c, const DPoint2d& d)
{
	const DPoint2d* pt[4];
	pt[0] = &a;
	pt[1] = &b;
	pt[2] = &c;
	pt[3] = &d;
	double alpha[4];
	int i, j;
	for(i = 0; i < 4; i++){
		double aq = DTriangle2d::alphaQuality(*pt[i], *pt[(i+1)%4], *pt[(i+2)%4]);
		// Sorting of alpha-values of 4 triangles
		for(j = i-1; j >= 0; j--){
			if(alpha[j] > aq) alpha[j+1] = alpha[j];
			else break;
		}
		alpha[j+1] = aq;
	}

	double denominator = alpha[2] * alpha[3];
	if(abs(denominator) > SMALL_NUMBER){
		return alpha[0] * alpha[1] / denominator;
	}else{
		return 0.0;
	}
}

/// Whether is valid (not inverted or self-crossing)
bool DQuad2d::valid(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c, const DPoint2d& d)
{
	double area1 = DTriangle2d::det(a, b, c);
	double area2 = DTriangle2d::det(c, d, a);

	if(area1 * area2 < 0.0){
		area1 = DTriangle2d::det(b, c, d);
		area2 = DTriangle2d::det(d, a, b);
	}

	if(area1 * area2 < 0.0) return false;
	else if(area1 < 0.0 && area2 < 0.0) return false;
	else return true;

}


/**
 * @param b is a reference to the second point of the quad,
 * @param c is a refenence to the third point of the quad,
 * @param d is a refenence to the fourth point of the quad,
 *
 * The coefficient is equal to:
 *  - 1 for equiangular quads,
 *  - 0 for degenerated quads, 
 *  - less than 0 for inverted quads(ie. oriented clockwise or skewed)
 *
 * \remark This is a "2-dimensional" version - it doesn't use any surface curvature information!
 */
double DQuad2d::alphaQuality(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c, const DMPoint2d& d)
{
	const DMPoint2d* pt[4];
	pt[0] = &a;
	pt[1] = &b;
	pt[2] = &c;
	pt[3] = &d;
	double alpha[4];
	int i, j;
	for(i = 0; i < 4; i++){
		double aq = DTriangle2d::alphaQuality(*pt[i], *pt[(i+1)%4], *pt[(i+2)%4]);
		// Sorting of alpha-values of 4 triangles
		for(j = i-1; j >= 0; j--){
			if(alpha[j] > aq) alpha[j+1] = alpha[j];
			else break;
		}
		alpha[j+1] = aq;
	}

	double denominator = alpha[2] * alpha[3];
	if(abs(denominator) > SMALL_NUMBER){
		return alpha[0] * alpha[1] / denominator;
	}else{
		return 0.0;
	}
}

/**
 * @param b is a reference to the second point of the quad,
 * @param c is a refenence to the third point of the quad,
 * @param d is a refenence to the fourth point of the quad,
 *
 * The coefficient is equal to:
 *  - 1 for equiangular quads,
 *  - 0 for degenerated quads, 
 *  - less than 0 for inverted quads(ie. oriented clockwise or skewed)
 *
 * \remark This is a "3-dimensional" version
 */
double DQuad3d::alphaQuality(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d)
{
	const DPoint3d* pt[4];
	pt[0] = &a;
	pt[1] = &b;
	pt[2] = &c;
	pt[3] = &d;
	double alpha[4];
	int i, j;
	for(i = 0; i < 4; i++){
		double aq = DTriangle3d::alphaQuality(*pt[i], *pt[(i+1)%4], *pt[(i+2)%4]);
		// Sorting of alpha-values of 4 triangles
		for(j = i-1; j >= 0; j--){
			if(alpha[j] > aq) alpha[j+1] = alpha[j];
			else break;
		}
		alpha[j+1] = aq;
	}

	double denominator = alpha[2] * alpha[3];
	if(abs(denominator) > SMALL_NUMBER){
		return alpha[0] * alpha[1] / denominator;
	}else{
		return 0.0;
	}
}

double DQuad3d::shapeQuality(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d)
{
	const DVector3d vn1 = DVector3d::crossProduct(a, b, c);
	const DVector3d vn2 = DVector3d::crossProduct(a, c, d);
	const DVector3d vn3 = DVector3d::crossProduct(a, b, d);
	const DVector3d vn4 = DVector3d::crossProduct(b, c, d);

	double s12 = (vn1.isZero() || vn2.isZero()) ? 0.0 : vn1.scalarProductNormalized(vn2);
	double s34 = (vn3.isZero() || vn4.isZero()) ? 0.0 : vn3.scalarProductNormalized(vn4);

	const double S_THR = 0.1;

	bool s12_pos = (s12 > S_THR);
	bool s34_pos = (s34 > S_THR);

	//if(s12 < 0.0 || s34 < 0.0){
	//	MeshViewSet* set = new MeshViewSet();
	//	set->addFace(this);
	//	set->addPoint(dpt0, 0, 0);
	//	set->addPoint(dpt1, 0, 1);
	//	set->addPoint(dpt2, 0, 2);
	//	set->addPoint(dpt3, 0, 3);
	//	ostringstream ostr;
	//	ostr << "Quad3d (#" << getIndex() << "): s12= " << s12 << " s34= " << s34;
	//	SHOW_MESH(ostr.str(), set);
	//}

	if(s12_pos && s34_pos) return 1.0; // convex
	else if(s12_pos || s34_pos) return 0.5; // concave
	else if(s12 == 0.0 || s34 == 0.0) return 0.0; // degenerate
	else return -1.0; // self-intersections or non-planar
}
