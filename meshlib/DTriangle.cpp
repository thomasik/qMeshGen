/////////////////////////////////////////////////////////////////////////////
// DTriangle.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "DTriangle.h"
#include "common.h"


/// Whether contains the given point
bool DTriangle2d::containsPoint(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c, 
								const DPoint2d& pt, double eps)
{
	if(DTriangle2d::det(a, b, c) > 0.0){
		return  DTriangle2d::det(a, b, pt) >= -eps &&
				DTriangle2d::det(b, c, pt) >= -eps &&
				DTriangle2d::det(c, a, pt) >= -eps;
	}else{
		return  DTriangle2d::det(a, b, pt) <= eps &&
				DTriangle2d::det(b, c, pt) <= eps &&
				DTriangle2d::det(c, a, pt) <= eps;
	}
}

/// Check if point is within the circumscribed circle for a triangle given by three points
double DTriangle2d::inSphereCheck(const DPoint2d& pt) const
{
	double ax = pt_a.x - pt.x;	double ay = pt_a.y - pt.y;
	double bx = pt_b.x - pt.x;	double by = pt_b.y - pt.y;
	double cx = pt_c.x - pt.x;	double cy = pt_c.y - pt.y;

	double aa = ax*ax + ay*ay;
	double bb = bx*bx + by*by;
	double cc = cx*cx + cy*cy;

	return ax*by*cc + ay*bb*cx + aa*bx*cy - cx*by*aa - bx*ay*cc - ax*cy*bb;
}

/// Check if point is within the circumscribed circle for a triangle given by three points
double DTriangle2d::inSphereCheck(const DPoint2d& pt, const DPoint2d& a, const DPoint2d& b, const DPoint2d& c)
{
	double ax = a.x - pt.x;	double ay = a.y - pt.y;
	double bx = b.x - pt.x;	double by = b.y - pt.y;
	double cx = c.x - pt.x;	double cy = c.y - pt.y;

	double aa = ax*ax + ay*ay;
	double bb = bx*bx + by*by;
	double cc = cx*cx + cy*cy;

	return ax*by*cc + ay*bb*cx + aa*bx*cy - cx*by*aa - bx*ay*cc - ax*cy*bb;
}

/// Check if point is within the circumscribed circle for a triangle given by three points
double DTriangle2d::inSphereCheck(const DMPoint2d& pt, const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c)
{
	double ax = a.x - pt.x;	double ay = a.y - pt.y;
	double bx = b.x - pt.x;	double by = b.y - pt.y;
	double cx = c.x - pt.x;	double cy = c.y - pt.y;

	double aa = ax*ax + ay*ay;
	double bb = bx*bx + by*by;
	double cc = cx*cx + cy*cy;

	return ax*by*cc + ay*bb*cx + aa*bx*cy - cx*by*aa - bx*ay*cc - ax*cy*bb;
}

double DTriangle2d::innerCircleRadius(const DPoint2d &a, const DPoint2d &b, const DPoint2d &c)
{
	double A = a.distance(b);
	double B = b.distance(c);
	double C = a.distance(c);
	double G = A + B + C;
	double D = a.x*(b.y - c.y) - b.x*(a.y - c.y) + c.x*(a.y - b.y);

	return D / G;
}

double DTriangle2d::innerCircleRadius(const DMPoint2d &a, const DMPoint2d &b, const DMPoint2d &c)
{
	double A = a.distance(b);
	double B = b.distance(c);
	double C = a.distance(c);
	double G = A + B + C;
	double D = a.x*(b.y - c.y) - b.x*(a.y - c.y) + c.x*(a.y - b.y);

	return D / G;
}

double DTriangle2d::outerCircleRadius(const DPoint2d &a, const DPoint2d &b, const DPoint2d &c)
{
	double nominator = 0.5 * sqrt(a.distance2(b) * a.distance2(c) * b.distance2(c));
	double denominator = (a.x-c.x)*(b.y-c.y) - (b.x-c.x)*(a.y-c.y);

	return nominator / denominator;
}

double DTriangle2d::outerCircleRadius(const DMPoint2d &a, const DMPoint2d &b, const DMPoint2d &c)
{
	double nominator = 0.5 * sqrt(a.distance2(b) * a.distance2(c) * b.distance2(c));
	double denominator = (a.x-c.x)*(b.y-c.y) - (b.x-c.x)*(a.y-c.y);

	return nominator / denominator;
}

double DTriangle2d::outerCircleRadius2(const DPoint2d &a, const DPoint2d &b, const DPoint2d &c)
{
	double nominator = 0.25 * a.distance2(b) * a.distance2(c) * b.distance2(c);
	double denominator = (a.x-c.x)*(b.y-c.y) - (b.x-c.x)*(a.y-c.y);

	return nominator / (denominator*denominator);
}

double DTriangle2d::outerCircleRadius2(const DMPoint2d &a, const DMPoint2d &b, const DMPoint2d &c)
{
	double nominator = 0.25 * a.distance2(b) * a.distance2(c) * b.distance2(c);
	double denominator = (a.x-c.x)*(b.y-c.y) - (b.x-c.x)*(a.y-c.y);

	return nominator / (denominator*denominator);
}

/**
 * @param b is a reference to the second point of the triangle,
 * @param c is a refenence to the third point of the triangle,
 *
 * The coefficient is equal to:
 *  - 1 for equilateral triangles,
 *  - 0 for degenerated triangles, 
 *  - less than 0 for inverted triangles (ie. oriented clockwise)
 *
 * \remark This is a "2-dimensional" version - it doesn't use any surface curvature information!
 */
double DTriangle2d::alphaQuality(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c)
{
	double denominator = a.distance2(b) + b.distance2(c) + a.distance2(c);
	return 2*SQRT3 * (b-a).crossProduct(c-a) / denominator;
}

/**
 * @param b is a reference to the second point of the triangle,
 * @param c is a refenence to the third point of the triangle,
 *
 * The coefficient is equal to:
 *  - 1 for equilateral triangles,
 *  - 0 for degenerated triangles, 
 *  - less than 0 for inverted triangles (ie. oriented clockwise)
 *
 * \remark This is a "2-dimensional" version - it doesn't use any surface curvature information!
 */
double DTriangle2d::alphaQuality(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c)
{
	double denominator = a.distance2(b) + b.distance2(c) + a.distance2(c);
	return 2*SQRT3 * (b-a).crossProduct(c-a) / denominator;
}

const DPoint2d DTriangle2d::outerCircleCenter(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c)
{
	double ax = a.x - c.x;	double ay = a.y - c.y;
	double bx = b.x - c.x;	double by = b.y - c.y;
	double denominator = 2.0 * (ax*by - ay*bx);
	if(abs(denominator) < SMALL_NUMBER) return DPoint2d::zero;
	denominator = 1.0/denominator;
	
	return DPoint2d(
		c.x - (ay*(bx*bx+by*by) - by*(ax*ax+ay*ay)) * denominator,
		c.y + (ax*(bx*bx+by*by) - bx*(ax*ax+ay*ay)) * denominator);
}

const DMPoint2d DTriangle2d::outerCircleCenter(const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c)
{
	double ax = a.x - c.x;	double ay = a.y - c.y;
	double bx = b.x - c.x;	double by = b.y - c.y;
	double denominator = 2.0 * (ax*by - ay*bx);
	if(abs(denominator) < SMALL_NUMBER) return DMPoint2d::zero;
	denominator = 1.0/denominator;
	
	return DMPoint2d(
		c.x - (ay*(bx*bx+by*by) - by*(ax*ax+ay*ay)) * denominator,
		c.y + (ax*(bx*bx+by*by) - bx*(ax*ax+ay*ay)) * denominator);
}

/// Distance (squared) from triangle to point
double DTriangle2d::distance2ToPoint(const DPoint2d& pt, const DPoint2d& a, const DPoint2d& b, const DPoint2d& c)
{
	return distance2ToPoint(pt, a, b-a, c-a);
}

/// Distance (squared) from triangle to point
double DTriangle2d::distance2ToPoint(const DPoint2d& pt) const
{
	return distance2ToPoint(pt, pt_a, pt_b, pt_c);
}

/// Distance (squared) from triangle to point
double DTriangle2d::distance2ToPoint(const DPoint2d& pt, const DPoint2d& a, const DVector2d& v0, const DVector2d& v1)
{
    const DVector2d d  = a-pt;

	double A00 = v0.length2();
	double A01 = v0.scalarProduct(v1);
	double A11 = v1.length2();
	double B0 = d.scalarProduct(v0); 
	double B1 = d.scalarProduct(v1); 
	double C = d.length2();
	double det_A = abs(A00*A11-A01*A01);
    double s = A01*B1-A11*B0;
    double t = A01*B0-A00*B1;
    double sqr_dist = 0.0;

	if(s+t <= det_A){
		if(s < 0.0){
			if(t < 0.0){  // region 4
                if(B0 < 0.0){
                    t = 0.0;
                    if(-B0 >= A00){
                        s = 1.0; 
						sqr_dist = A00+2*B0+C;
                    }else{
                        s = -B0/A00;
                        sqr_dist = B0*s+C;
                    }
                }else{
                    s = 0.0;
                    if(B1 >= 0.0){
                        t = 0.0;
                        sqr_dist = C;
                    }else if(-B1 >= A11){
                        t = 1.0;
                        sqr_dist = A11+2*B1+C;
                    }else{
                        t = -B1/A11;
                        sqr_dist = B1*t+C;
                    }
                }
			}else{  // region 3
                s = 0.0;
                if(B1 >= 0.0){
                    t = 0.0;
                    sqr_dist = C;
                }else if(-B1 >= A11){
                    t = 1.0;
                    sqr_dist = A11+2*B1+C;
                }else{
                    t = -B1/A11;
                    sqr_dist = B1*t+C;
                }
            }
		}else if(t < 0.0){  // region 5
            t = 0.0;
            if(B0 >= 0.0){
                s = 0.0;
                sqr_dist = C;
            }else if(-B0 >= A00){
                s = 1.0;
                sqr_dist = A00+2*B0+C;
            }else{
                s = -B0/A00;
                sqr_dist = B0*s+C;
            }
		}else{  // region 0
            // minimum at interior point
            // double inv_det = 1.0/det_A;
            // s *= inv_det;
            // t *= inv_det;
            // sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
			sqr_dist = 0.0;
     }
    }else{
        double tmp0, tmp1, numer, denom;

		if(s < 0.0){  // region 2
            tmp0 = A01 + B0;
            tmp1 = A11 + B1;
            if(tmp1 > tmp0){
                numer = tmp1 - tmp0;
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqr_dist = A00+2*B0+C;
                }else{
                    s = numer/denom;
                    t = 1.0 - s;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }else{
                s = 0.0;
                if(tmp1 <= 0.0){
                    t = 1.0;
                    sqr_dist = A11+2*B1+C;
                }else if(B1 >= 0.0){
                    t = 0.0;
                    sqr_dist = C;
                }else{
                    t = -B1/A11;
                    sqr_dist = B1*t+C;
                }
            }
		}else if(t < 0.0){  // region 6
            tmp0 = A01 + B1;
            tmp1 = A00 + B0;
            if(tmp1 > tmp0){
                numer = tmp1 - tmp0;
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    t = 1.0;
                    s = 0.0;
                    sqr_dist = A11+2*B1+C;
                }else{
                    t = numer/denom;
                    s = 1.0 - t;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }else{
                t = 0.0;
                if(tmp1 <= 0.0){
                    s = 1.0;
                    sqr_dist = A00+2*B0+C;
                }else if(B0 >= 0.0){
                    s = 0.0;
                    sqr_dist = C;
                }else{
                    s = -B0/A00;
                    sqr_dist = B0*s+C;
                }
            }
		}else{  // region 1
            numer = A11 + B1 - A01 - B0;
            if(numer <= 0.0){
                s = 0.0;
                t = 1.0;
                sqr_dist = A11+2*B1+C;
            }else{
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqr_dist = A00+2*B0+C;
                }else{
                    s = numer/denom;
                    t = 1.0 - s;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }
        }
    }

    return abs(sqr_dist);
}

/// Distance (squared) from triangle to point
double DTriangle2d::distance2ToPoint(const DMPoint2d& pt, const DMPoint2d& a, const DMPoint2d& b, const DMPoint2d& c)
{
	return distance2ToPoint(pt, a, b-a, c-a);
}

/// Distance (squared) from triangle to point
double DTriangle2d::distance2ToPoint(const DMPoint2d& pt, const DMPoint2d& a, const DMVector2d& v0, const DMVector2d& v1)
{
    const DMVector2d d  = a-pt;

	double A00 = v0.length2();
	double A01 = v0.scalarProduct(v1);
	double A11 = v1.length2();
	double B0 = d.scalarProduct(v0); 
	double B1 = d.scalarProduct(v1); 
	double C = d.length2();
	double det_A = abs(A00*A11-A01*A01);
    double s = A01*B1-A11*B0;
    double t = A01*B0-A00*B1;
    double sqr_dist = 0.0;

	if(s+t <= det_A){
		if(s < 0.0){
			if(t < 0.0){  // region 4
                if(B0 < 0.0){
                    t = 0.0;
                    if(-B0 >= A00){
                        s = 1.0; 
						sqr_dist = A00+2*B0+C;
                    }else{
                        s = -B0/A00;
                        sqr_dist = B0*s+C;
                    }
                }else{
                    s = 0.0;
                    if(B1 >= 0.0){
                        t = 0.0;
                        sqr_dist = C;
                    }else if(-B1 >= A11){
                        t = 1.0;
                        sqr_dist = A11+2*B1+C;
                    }else{
                        t = -B1/A11;
                        sqr_dist = B1*t+C;
                    }
                }
			}else{  // region 3
                s = 0.0;
                if(B1 >= 0.0){
                    t = 0.0;
                    sqr_dist = C;
                }else if(-B1 >= A11){
                    t = 1.0;
                    sqr_dist = A11+2*B1+C;
                }else{
                    t = -B1/A11;
                    sqr_dist = B1*t+C;
                }
            }
		}else if(t < 0.0){  // region 5
            t = 0.0;
            if(B0 >= 0.0){
                s = 0.0;
                sqr_dist = C;
            }else if(-B0 >= A00){
                s = 1.0;
                sqr_dist = A00+2*B0+C;
            }else{
                s = -B0/A00;
                sqr_dist = B0*s+C;
            }
		}else{  // region 0
            // minimum at interior point
            // double inv_det = 1.0/det_A;
            // s *= inv_det;
            // t *= inv_det;
            // sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
			sqr_dist = 0.0;
     }
    }else{
        double tmp0, tmp1, numer, denom;

		if(s < 0.0){  // region 2
            tmp0 = A01 + B0;
            tmp1 = A11 + B1;
            if(tmp1 > tmp0){
                numer = tmp1 - tmp0;
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqr_dist = A00+2*B0+C;
                }else{
                    s = numer/denom;
                    t = 1.0 - s;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }else{
                s = 0.0;
                if(tmp1 <= 0.0){
                    t = 1.0;
                    sqr_dist = A11+2*B1+C;
                }else if(B1 >= 0.0){
                    t = 0.0;
                    sqr_dist = C;
                }else{
                    t = -B1/A11;
                    sqr_dist = B1*t+C;
                }
            }
		}else if(t < 0.0){  // region 6
            tmp0 = A01 + B1;
            tmp1 = A00 + B0;
            if(tmp1 > tmp0){
                numer = tmp1 - tmp0;
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    t = 1.0;
                    s = 0.0;
                    sqr_dist = A11+2*B1+C;
                }else{
                    t = numer/denom;
                    s = 1.0 - t;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }else{
                t = 0.0;
                if(tmp1 <= 0.0){
                    s = 1.0;
                    sqr_dist = A00+2*B0+C;
                }else if(B0 >= 0.0){
                    s = 0.0;
                    sqr_dist = C;
                }else{
                    s = -B0/A00;
                    sqr_dist = B0*s+C;
                }
            }
		}else{  // region 1
            numer = A11 + B1 - A01 - B0;
            if(numer <= 0.0){
                s = 0.0;
                t = 1.0;
                sqr_dist = A11+2*B1+C;
            }else{
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqr_dist = A00+2*B0+C;
                }else{
                    s = numer/denom;
                    t = 1.0 - s;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }
        }
    }

    return abs(sqr_dist);
}
/**
 * @param b is a reference to the second point of the triangle,
 * @param c is a refenence to the third point of the triangle,
 *
 * The coefficient is equal to:
 *  - 1 for equilateral triangles,
 *  - 0 for degenerated triangles, 
 *  - less than 0 for inverted triangles (ie. oriented clockwise)
 *
 * \remark This is a "3-dimensional" version
 */
double DTriangle3d::alphaQuality(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c)
{
	double denominator = a.distance2(b) + b.distance2(c) + a.distance2(c);
	return 2*SQRT3 * (b-a).crossProduct(c-a).length() / denominator;
}

double DMTriangle3d::alphaQuality(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c)
{
	double denominator = a.distance2(b) + b.distance2(c) + a.distance2(c);
	return 2*SQRT3 * (b-a).crossProduct(c-a).length() / denominator;
}

double DTriangle3d::orient3d(const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, const DPoint3d& d)
{
  double adx = a.x - d.x;
  double bdx = b.x - d.x;
  double cdx = c.x - d.x;
  double ady = a.y - d.y;
  double bdy = b.y - d.y;
  double cdy = c.y - d.y;
  double adz = a.z - d.z;
  double bdz = b.z - d.z;
  double cdz = c.z - d.z;

  double bdxcdy = bdx * cdy;
  double cdxbdy = cdx * bdy;

  double cdxady = cdx * ady;
  double adxcdy = adx * cdy;

  double adxbdy = adx * bdy;
  double bdxady = bdx * ady;

  return adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
}

double DTriangle3d::orient3d(const DPoint3d& d) const
{
	return orient3d(pt_a, pt_b, pt_c, d);
}

double DMTriangle3d::orient3d(const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c, const DMPoint3d& d)
{
  double adx = a.x - d.x;
  double bdx = b.x - d.x;
  double cdx = c.x - d.x;
  double ady = a.y - d.y;
  double bdy = b.y - d.y;
  double cdy = c.y - d.y;
  double adz = a.z - d.z;
  double bdz = b.z - d.z;
  double cdz = c.z - d.z;

  double bdxcdy = bdx * cdy;
  double cdxbdy = cdx * bdy;

  double cdxady = cdx * ady;
  double adxcdy = adx * cdy;

  double adxbdy = adx * bdy;
  double bdxady = bdx * ady;

  return adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
}

double DMTriangle3d::orient3d(const DMPoint3d& d) const
{
	return orient3d(pt_a, pt_b, pt_c, d);
}

/// Whether the point lies within the triangle after planar projection
bool DTriangle3d::containsProjectedPoint(const DPoint3d& pt, const DPoint3d& a, const DPoint3d& b, const DPoint3d& c, double eps )
{
	return containsProjectedPoint( pt, a, b-a, c-a, eps );
}
	
/// Whether the point lies within the triangle after planar projection
bool DTriangle3d::containsProjectedPoint(const DPoint3d& pt, const DPoint3d& a, const DVector3d& ab, const DVector3d& ac, double eps )
{
    const DVector3d d  = a-pt;

	double A00 = ab.length2();
	double A01 = ab.scalarProduct(ac);
	double A11 = ac.length2();
	double B0 = d.scalarProduct(ab); 
	double B1 = d.scalarProduct(ac); 
	double C = d.length2();
	double det_A = abs(A00*A11-A01*A01);
    double s = A01*B1-A11*B0;
    double t = A01*B0-A00*B1;

	double rel_eps = eps * det_A;
	return ( s >= -rel_eps && t >= -rel_eps && (s+t <= (det_A+rel_eps)) );
}

/// Distance (squared) from triangle to point
double DTriangle3d::distance2ToPoint(const DPoint3d& pt, const DPoint3d& a, const DPoint3d& b, const DPoint3d& c)
{
	return distance2ToPoint(pt, a, b-a, c-a);
}

/// Distance (squared) from triangle to point
double DTriangle3d::distance2ToPointOrthogonal(const DPoint3d& pt, const DPoint3d& a, const DPoint3d& b, const DPoint3d& c)
{
	const DVector3d v0 = b-a;
	const DVector3d v1 = c-a;
    const DVector3d d  = a-pt;

	double A00 = v0.length2();
	double A01 = v0.scalarProduct(v1);
	double A11 = v1.length2();
	double B0 = d.scalarProduct(v0); 
	double B1 = d.scalarProduct(v1); 
	double C = d.length2();
	double det_A = abs(A00*A11-A01*A01);
    double s = A01*B1-A11*B0;
    double t = A01*B0-A00*B1;

	if( s >= 0.0 && t >= 0.0 && (s+t <= det_A) ){
		double inv_det = 1.0/det_A;
        s *= inv_det;
        t *= inv_det;
        return s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
	}else
		return -1.0;
}

/// Distance (squared) from triangle to point
double DTriangle3d::distance2ToPoint(const DPoint3d& pt, const DPoint3d& a, const DVector3d& v0, const DVector3d& v1)
{
    const DVector3d d  = a-pt;

	double A00 = v0.length2();
	double A01 = v0.scalarProduct(v1);
	double A11 = v1.length2();
	double B0 = d.scalarProduct(v0); 
	double B1 = d.scalarProduct(v1); 
	double C = d.length2();
	double det_A = abs(A00*A11-A01*A01);
    double s = A01*B1-A11*B0;
    double t = A01*B0-A00*B1;
    double sqr_dist = 0.0;

	if(s+t <= det_A){
		if(s < 0.0){
			if(t < 0.0){  // region 4
                if(B0 < 0.0){
                    t = 0.0;
                    if(-B0 >= A00){
                        s = 1.0; 
						sqr_dist = A00+2*B0+C;
                    }else{
                        s = -B0/A00;
                        sqr_dist = B0*s+C;
                    }
                }else{
                    s = 0.0;
                    if(B1 >= 0.0){
                        t = 0.0;
                        sqr_dist = C;
                    }else if(-B1 >= A11){
                        t = 1.0;
                        sqr_dist = A11+2*B1+C;
                    }else{
                        t = -B1/A11;
                        sqr_dist = B1*t+C;
                    }
                }
			}else{  // region 3
                s = 0.0;
                if(B1 >= 0.0){
                    t = 0.0;
                    sqr_dist = C;
                }else if(-B1 >= A11){
                    t = 1.0;
                    sqr_dist = A11+2*B1+C;
                }else{
                    t = -B1/A11;
                    sqr_dist = B1*t+C;
                }
            }
		}else if(t < 0.0){  // region 5
            t = 0.0;
            if(B0 >= 0.0){
                s = 0.0;
                sqr_dist = C;
            }else if(-B0 >= A00){
                s = 1.0;
                sqr_dist = A00+2*B0+C;
            }else{
                s = -B0/A00;
                sqr_dist = B0*s+C;
            }
		}else{  // region 0
            // minimum at interior point
            double inv_det = 1.0/det_A;
            s *= inv_det;
            t *= inv_det;
            sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
     }
    }else{
        double tmp0, tmp1, numer, denom;

		if(s < 0.0){  // region 2
            tmp0 = A01 + B0;
            tmp1 = A11 + B1;
            if(tmp1 > tmp0){
                numer = tmp1 - tmp0;
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqr_dist = A00+2*B0+C;
                }else{
                    s = numer/denom;
                    t = 1.0 - s;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }else{
                s = 0.0;
                if(tmp1 <= 0.0){
                    t = 1.0;
                    sqr_dist = A11+2*B1+C;
                }else if(B1 >= 0.0){
                    t = 0.0;
                    sqr_dist = C;
                }else{
                    t = -B1/A11;
                    sqr_dist = B1*t+C;
                }
            }
		}else if(t < 0.0){  // region 6
            tmp0 = A01 + B1;
            tmp1 = A00 + B0;
            if(tmp1 > tmp0){
                numer = tmp1 - tmp0;
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    t = 1.0;
                    s = 0.0;
                    sqr_dist = A11+2*B1+C;
                }else{
                    t = numer/denom;
                    s = 1.0 - t;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }else{
                t = 0.0;
                if(tmp1 <= 0.0){
                    s = 1.0;
                    sqr_dist = A00+2*B0+C;
                }else if(B0 >= 0.0){
                    s = 0.0;
                    sqr_dist = C;
                }else{
                    s = -B0/A00;
                    sqr_dist = B0*s+C;
                }
            }
		}else{  // region 1
            numer = A11 + B1 - A01 - B0;
            if(numer <= 0.0){
                s = 0.0;
                t = 1.0;
                sqr_dist = A11+2*B1+C;
            }else{
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqr_dist = A00+2*B0+C;
                }else{
                    s = numer/denom;
                    t = 1.0 - s;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }
        }
    }

    return abs(sqr_dist);
}

/// Distance (squared) from triangle to point
double DMTriangle3d::distance2ToPoint(const DMPoint3d& pt, const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c)
{
	const DMVector3d v0 = b-a;
	const DMVector3d v1 = c-a;
    const DMVector3d d  = a-pt;

	double A00 = v0.length2();
	double A01 = v0.scalarProduct(v1);
	double A11 = v1.length2();
	double B0 = d.scalarProduct(v0); 
	double B1 = d.scalarProduct(v1); 
	double C = d.length2();
	double det_A = abs(A00*A11-A01*A01);
    double s = A01*B1-A11*B0;
    double t = A01*B0-A00*B1;
    double sqr_dist;

	if(s+t <= det_A){
		if(s < 0.0){
			if(t < 0.0){  // region 4
                if(B0 < 0.0){
                    t = 0.0;
                    if(-B0 >= A00){
                        s = 1.0; 
						sqr_dist = A00+2*B0+C;
                    }else{
                        s = -B0/A00;
                        sqr_dist = B0*s+C;
                    }
                }else{
                    s = 0.0;
                    if(B1 >= 0.0){
                        t = 0.0;
                        sqr_dist = C;
                    }else if(-B1 >= A11){
                        t = 1.0;
                        sqr_dist = A11+2*B1+C;
                    }else{
                        t = -B1/A11;
                        sqr_dist = B1*t+C;
                    }
                }
			}else{  // region 3
                s = 0.0;
                if(B1 >= 0.0){
                    t = 0.0;
                    sqr_dist = C;
                }else if(-B1 >= A11){
                    t = 1.0;
                    sqr_dist = A11+2*B1+C;
                }else{
                    t = -B1/A11;
                    sqr_dist = B1*t+C;
                }
            }
		}else if(t < 0.0){  // region 5
            t = 0.0;
            if(B0 >= 0.0){
                s = 0.0;
                sqr_dist = C;
            }else if(-B0 >= A00){
                s = 1.0;
                sqr_dist = A00+2*B0+C;
            }else{
                s = -B0/A00;
                sqr_dist = B0*s+C;
            }
		}else{  // region 0
            // minimum at interior point
            double inv_det = 1.0/det_A;
            s *= inv_det;
            t *= inv_det;
            sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
        }
    }else{
        double tmp0, tmp1, numer, denom;

		if(s < 0.0){  // region 2
            tmp0 = A01 + B0;
            tmp1 = A11 + B1;
            if(tmp1 > tmp0){
                numer = tmp1 - tmp0;
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqr_dist = A00+2*B0+C;
                }else{
                    s = numer/denom;
                    t = 1.0 - s;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }else{
                s = 0.0;
                if(tmp1 <= 0.0){
                    t = 1.0;
                    sqr_dist = A11+2*B1+C;
                }else if(B1 >= 0.0){
                    t = 0.0;
                    sqr_dist = C;
                }else{
                    t = -B1/A11;
                    sqr_dist = B1*t+C;
                }
            }
		}else if(t < 0.0){  // region 6
            tmp0 = A01 + B1;
            tmp1 = A00 + B0;
            if(tmp1 > tmp0){
                numer = tmp1 - tmp0;
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    t = 1.0;
                    s = 0.0;
                    sqr_dist = A11+2*B1+C;
                }else{
                    t = numer/denom;
                    s = 1.0 - t;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }else{
                t = 0.0;
                if(tmp1 <= 0.0){
                    s = 1.0;
                    sqr_dist = A00+2*B0+C;
                }else if(B0 >= 0.0){
                    s = 0.0;
                    sqr_dist = C;
                }else{
                    s = -B0/A00;
                    sqr_dist = B0*s+C;
                }
            }
		}else{  // region 1
            numer = A11 + B1 - A01 - B0;
            if(numer <= 0.0){
                s = 0.0;
                t = 1.0;
                sqr_dist = A11+2*B1+C;
            }else{
                denom = A00-2*A01+A11;
                if(numer >= denom){
                    s = 1.0;
                    t = 0.0;
                    sqr_dist = A00+2*B0+C;
                }else{
                    s = numer/denom;
                    t = 1.0 - s;
                    sqr_dist = s*(A00*s+A01*t+2*B0) + t*(A01*s+A11*t+2*B1)+C;
                }
            }
        }
    }

    return abs(sqr_dist);
}

/// Distance (squared) from triangle to point
double DTriangle3d::distance2ToPoint(const DPoint3d& pt) const
{
	return distance2ToPoint(pt, pt_a, pt_b, pt_c);
}

/// Distance (squared) from triangle to point
double DMTriangle3d::distance2ToPoint(const DMPoint3d& pt) const
{
	return distance2ToPoint(pt, pt_a, pt_b, pt_c);
}

/// Check if segment |p1p2| is crossing triangle |p3p4p5|
bool DTriangle3d::crossSegment(const DPoint3d& sp1, const DPoint3d& sp2, 
					const DPoint3d& a, const DPoint3d& b, const DPoint3d& c)
{
	double v01 = DTriangle3d::orient3d(a, b, c, sp1);
	//if(abs(v01) < eps) return false; // one end of the segment too close to the triangle
	double v02 = DTriangle3d::orient3d(a, b, c, sp2);
	//if(abs(v02) < eps) return false; // other end of the segment too close to the triangle

	if(v01 * v02 >= 0.0) return false; // both vertices of segment are on the same side of the face

	double v1 = DTriangle3d::orient3d(a, b, sp1, sp2);
	double v2 = DTriangle3d::orient3d(b, c, sp1, sp2);
	double v3 = DTriangle3d::orient3d(c, a, sp1, sp2);

	return ((v1 >= 0.0 && v2 >= 0.0 && v3 >= 0.0) || // segment must be on the same "side" of all edges of face
		(v1 <= 0.0 && v2 <= 0.0 && v3 <= 0.0));
}

/// Check if segment |p1p2| is crossing triangle |p3p4p5|
bool DTriangle3d::crossSegment(const DPoint3d& sp1, const DPoint3d& sp2) const
{
	return crossSegment(sp1, sp2, pt_a, pt_b, pt_c);
}

/// Check if segment |p1p2| is crossing triangle |p3p4p5|
bool DMTriangle3d::crossSegment(const DMPoint3d& sp1, const DMPoint3d& sp2, 
					const DMPoint3d& a, const DMPoint3d& b, const DMPoint3d& c)
{
	double v01 = DMTriangle3d::orient3d(a, b, c, sp1);
	//if(abs(v01) < eps) return false; // one end of the segment too close to the triangle
	double v02 = DMTriangle3d::orient3d(a, b, c, sp2);
	//if(abs(v02) < eps) return false; // other end of the segment too close to the triangle

	if(v01 * v02 >= 0.0) return false; // both vertices of segment are on the same side of the face

	double v1 = DMTriangle3d::orient3d(a, b, sp1, sp2);
	double v2 = DMTriangle3d::orient3d(b, c, sp1, sp2);
	double v3 = DMTriangle3d::orient3d(c, a, sp1, sp2);

	return ((v1 >= 0.0 && v2 >= 0.0 && v3 >= 0.0) || // segment must be on the same "side" of all edges of face
		(v1 <= 0.0 && v2 <= 0.0 && v3 <= 0.0));
}

/// Check if segment |p1p2| is crossing triangle |p3p4p5|
bool DMTriangle3d::crossSegment(const DMPoint3d& sp1, const DMPoint3d& sp2) const
{
	return crossSegment(sp1, sp2, pt_a, pt_b, pt_c);
}

/// Return normal vector
DVector3d DTriangle3d::normalVector() const
{
	return DVector3d::crossProduct(pt_a, pt_b, pt_c).normalized();
}

/// Return normal vector
DMVector3d DMTriangle3d::normalVector() const
{
	return DMVector3d::crossProduct(pt_a, pt_b, pt_c).normalized();
}
