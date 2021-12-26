/////////////////////////////////////////////////////////////////////////////
// DTetrahedron.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "DTetrahedron.h"
#include "common.h"
#include "DMatrix.h"

double DTetrahedron::insphere(const DPoint3d& pt, const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d)
{
  double aex = a.x - pt.x;
  double bex = b.x - pt.x;
  double cex = c.x - pt.x;
  double dex = d.x - pt.x;
  double aey = a.y - pt.y;
  double bey = b.y - pt.y;
  double cey = c.y - pt.y;
  double dey = d.y - pt.y;
  double aez = a.z - pt.z;
  double bez = b.z - pt.z;
  double cez = c.z - pt.z;
  double dez = d.z - pt.z;

  double ab = aex * bey - bex * aey;
  double bc = bex * cey - cex * bey;
  double cd = cex * dey - dex * cey;
  double da = dex * aey - aex * dey;

  double ac = aex * cey - cex * aey;
  double bd = bex * dey - dex * bey;

  double abc = aez * bc - bez * ac + cez * ab;
  double bcd = bez * cd - cez * bd + dez * bc;
  double cda = cez * da + dez * ac + aez * cd;
  double dab = dez * ab + aez * bd + bez * da;

  double alift = aex * aex + aey * aey + aez * aez;
  double blift = bex * bex + bey * bey + bez * bez;
  double clift = cex * cex + cey * cey + cez * cez;
  double dlift = dex * dex + dey * dey + dez * dez;

  return (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
}

double DTetrahedron::insphere(const DMPoint3d& pt, const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d)
{
	return DTetrahedron::insphere(
		pt.toRealSpace(),
		a.toRealSpace(), b.toRealSpace(),
		c.toRealSpace(), d.toRealSpace());
}

double DTetrahedron::outerSphereRadius(const DPoint3d &a, const DPoint3d &b, const DPoint3d &c, const DPoint3d &d)
{
	const DVector3d ad = a-d;
	const DVector3d bd = b-d;
	const DVector3d cd = c-d;
	double denominator = 2.0 * DMatrix3d::det33(ad.x, ad.y, ad.z, bd.x, bd.y, bd.z, cd.x, cd.y, cd.z);
	assert(	abs(denominator) > SMALL_NUMBER);
	DVector3d pt = bd.crossProduct(cd) * ad.length2() + 
		cd.crossProduct(ad) * bd.length2() + 
		ad.crossProduct(bd) * cd.length2();
	return (pt / denominator).length(); 
}

double DTetrahedron::outerSphereRadius(const DMPoint3d &a, const DMPoint3d &b, const DMPoint3d &c, const DMPoint3d &d)
{
	const DMVector3d ad = a-d;
	const DMVector3d bd = b-d;
	const DMVector3d cd = c-d;
	double denominator = 2.0 * DMatrix3d::det33(ad.x, ad.y, ad.z, bd.x, bd.y, bd.z, cd.x, cd.y, cd.z);
	assert(	abs(denominator) > SMALL_NUMBER);
	DMVector3d pt = bd.crossProduct(cd) * ad.length2() + 
		cd.crossProduct(ad) * bd.length2() + 
		ad.crossProduct(bd) * cd.length2();
	return (pt / denominator).length(); 
}

double DTetrahedron::outerSphereRadius2(const DPoint3d &a, const DPoint3d &b, const DPoint3d &c, const DPoint3d &d)
{
	const DVector3d ad = a-d;
	const DVector3d bd = b-d;
	const DVector3d cd = c-d;
	double denominator = 2.0 * DMatrix3d::det33(ad.x, ad.y, ad.z, bd.x, bd.y, bd.z, cd.x, cd.y, cd.z);
	assert(	abs(denominator) > SMALL_NUMBER);
	DVector3d pt = bd.crossProduct(cd) * ad.length2() + 
		cd.crossProduct(ad) * bd.length2() + 
		ad.crossProduct(bd) * cd.length2();
	return (pt / denominator).length2(); 
}

double DTetrahedron::outerSphereRadius2(const DMPoint3d &a, const DMPoint3d &b, const DMPoint3d &c, const DMPoint3d &d)
{
	const DMVector3d ad = a-d;
	const DMVector3d bd = b-d;
	const DMVector3d cd = c-d;
	double denominator = 2.0 * DMatrix3d::det33(ad.x, ad.y, ad.z, bd.x, bd.y, bd.z, cd.x, cd.y, cd.z);
//	assert(	abs(denominator) > SMALL_NUMBER);
	assert(	denominator != 0.0);
	DMVector3d pt = bd.crossProduct(cd) * ad.length2() + 
		cd.crossProduct(ad) * bd.length2() + 
		ad.crossProduct(bd) * cd.length2();
	return (pt / denominator).length2(); 
}

DPoint3d DTetrahedron::outerSphereCenter(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d)
{
	const DVector3d ad = a-d;
	const DVector3d bd = b-d;
	const DVector3d cd = c-d;
	double denominator = 2.0 * DMatrix3d::det33(ad.x, ad.y, ad.z, bd.x, bd.y, bd.z, cd.x, cd.y, cd.z);
//	assert(	abs(denominator) > SMALL_NUMBER);
	assert(	denominator != 0.0);
	DVector3d vt = bd.crossProduct(cd) * ad.length2() + 
		cd.crossProduct(ad) * bd.length2() + 
		ad.crossProduct(bd) * cd.length2();
	return d + vt / denominator; 
}

DMPoint3d DTetrahedron::outerSphereCenter(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d)
{
	const DMVector3d ad = a-d;
	const DMVector3d bd = b-d;
	const DMVector3d cd = c-d;
	double denominator = 2.0 * DMatrix3d::det33(ad.x, ad.y, ad.z, bd.x, bd.y, bd.z, cd.x, cd.y, cd.z);
//	assert(	abs(denominator) > SMALL_NUMBER);
	assert(	denominator != 0.0);
	DMVector3d vt = bd.crossProduct(cd) * ad.length2() + 
		cd.crossProduct(ad) * bd.length2() + 
		ad.crossProduct(bd) * cd.length2();
	return d + vt / denominator; 
}

double DTetrahedron::aspectRatio(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d)
{
	double V = DTetrahedron::volume(a, b, c, d);
	if(V < 0.0) return -1.0;

	//MeshData::QUALITY3D_ASPECT_RATIO
	const double f = 6.0 * sqrt(6.0);	// scaling factor
	double h1 = a.distance2(b);
	double h2 = b.distance2(c);
	double h3 = c.distance2(a);
	double h4 = a.distance2(d);
	double h5 = b.distance2(d);
	double h6 = c.distance2(d);
	double hM = sqrt(std::max(std::max(h1,h2),std::max(std::max(h3,h4), std::max(h5,h6))));

	double S = 
		DTriangle3d::area(a, b, c) + 
		DTriangle3d::area(a, b, d) + 
		DTriangle3d::area(b, c, d) + 
		DTriangle3d::area(c, a, d);

	return f * V / (hM * S);
}

double DTetrahedron::aspectRatio(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d)
{
	double V = DTetrahedron::volume(a, b, c, d);
	if(V < -METRIC_SMALL_NUMBER) return -1.0;
	else if(V < METRIC_SMALL_NUMBER) return 0.0;

	//MeshData::QUALITY3D_ASPECT_RATIO
	const double f = 6.0 * sqrt(6.0);	// scaling factor
	double h1 = a.distance2(b);
	double h2 = b.distance2(c);
	double h3 = c.distance2(a);
	double h4 = a.distance2(d);
	double h5 = b.distance2(d);
	double h6 = c.distance2(d);
	double hM = sqrt(std::max(std::max(h1,h2),std::max(std::max(h3,h4), std::max(h5,h6))));

	double S = 
		DMTriangle3d::area(a, b, c) + 
		DMTriangle3d::area(a, b, d) + 
		DMTriangle3d::area(b, c, d) + 
		DMTriangle3d::area(c, a, d);

	return f * V / (hM * S);
}
