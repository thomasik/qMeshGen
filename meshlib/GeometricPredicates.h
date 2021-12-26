/**
 *  Routines for Arbitrary Precision Floating-point Arithmetic
 *  and Fast Robust Geometric Predicates
 *  (predicates.c)
 *  
 *  May 18, 1996
 *
 *  Placed in the public domain by
 *  Jonathan Richard Shewchuk
 *  School of Computer Science
 *  Carnegie Mellon University
 *  5000 Forbes Avenue
 *  PIttsburgh, Pennsylvania  15213-3891
 *  jrs@cs.cmu.edu
*/

#if !defined(GEOMETRICPREDICATES_H__INCLUDED)
#define GEOMETRICPREDICATES_H__INCLUDED

#include "DPoint.h"

class GeometricPredicates  
{
public:
	struct InnerData{
		double splitter;     ///< = 2^ceiling(p / 2) + 1.  Used to split floats in half.
		double epsilon;      ///< = 2^(-p).  Used to estimate roundoff errors.
		/* A set of coefficients used to calculate maximum roundoff errors.          */
		double resulterrbound;
		double ccwerrboundA, ccwerrboundB, ccwerrboundC;
		double o3derrboundA, o3derrboundB, o3derrboundC;
		double iccerrboundA, iccerrboundB, iccerrboundC;
		double isperrboundA, isperrboundB, isperrboundC;
	public:
		/// Constructor - Initialize the variables used for exact arithmetic.
		InnerData();
	};
private:
	/// Adaptive exact 2D incircle test.
	static double incircle_adapt(const DPoint2d& pa, const DPoint2d& pb, const DPoint2d& pc, const DPoint2d& pd, double permanent);
	/// Adaptive exact 3D incircle test.
	static double insphere_adapt(const DPoint3d& pa, const DPoint3d& pb, const DPoint3d& pc, const DPoint3d& pd, const DPoint3d& pe, double permanent);
	/// Exact 3D insphere test. Robust.
	static double insphere_exact(const DPoint3d& pa, const DPoint3d& pb, const DPoint3d& pc, const DPoint3d& pd, const DPoint3d& pe);
	/// Adaptive exact 2D orientation test.
	static double orient2d_adapt(const DPoint2d& pa, const DPoint2d& pb, const DPoint2d&, double detsum);
	/// Adaptive exact 3D orientation test.
	static double orient3d_adapt(const DPoint3d& pa, const DPoint3d& pb, const DPoint3d& pc, const DPoint3d& pd, double permanent);
public:
	/// Adaptive exact 3D insphere test.  Robust.
	static double insphere(const DPoint3d& pa, const DPoint3d& pb, const DPoint3d& pc, const DPoint3d& pd, const DPoint3d& pe);
	/// Approximate 3D insphere test.  Nonrobust.
	static double insphere_fast(const DPoint3d& pa, const DPoint3d& pb, const DPoint3d& pc, const DPoint3d& pd, const DPoint3d& pe);
	/// Adaptive exact 2D incircle test.  Robust.
	static double incircle(const DPoint2d& pa, const DPoint2d& pb, const DPoint2d& pc, const DPoint2d& pd);
	/// Approximate 2D incircle test.  Nonrobust.
	static double incircle_fast(const DPoint2d& pa, const DPoint2d& pb, const DPoint2d& pc, const DPoint2d& pd);
	/// Adaptive exact tetrahedra area (may be negative).  Robust.
	static double tetrahedra_volume(const DPoint3d& pa, const DPoint3d& pb, const DPoint3d& pc, const DPoint3d& pd){
		return (1.0/6.0)*orient3d(pa, pb, pc, pd);
	}
	/// Adaptive exact 3D orientation test.  Robust.
	static double orient3d(const DPoint3d& pa, const DPoint3d& pb, const DPoint3d& pc, const DPoint3d& pd);
	/// Approximate 3D orientation test.  Nonrobust.
	static double orient3d_fast(const DPoint3d& pa, const DPoint3d& pb, const DPoint3d& pc, const DPoint3d& pd);
	/// Adaptive exact triangle area (may be negative).  Robust.
	static double triangle_area(const DPoint2d& pa, const DPoint2d& pb, const DPoint2d& pc){
		return 0.5*orient2d_fast(pa, pb, pc);
	}
	/// Adaptive exact 2D orientation test.  Robust.
	static double orient2d(const DPoint2d& pa, const DPoint2d& pb, const DPoint2d& pc);
	/// Approximate 2D orientation test.  Nonrobust.
	static double orient2d_fast(const DPoint2d& pa, const DPoint2d& pb, const DPoint2d& pc);
	/// Produce a one-word estimate of an expansion's value.
	static double estimate(int elen, double* e);
	/// Compress an expansion.
	static int compress(int elen, double* e, double* h);
	/// Multiply an expansion by a scalar, eliminating zero components from the output expansion.
	static int scale_expansion_zeroelim(int elen, double* e, double b, double* h);
	/// Multiply an expansion by a scalar.
	static int scale_expansion(int elen, double* e, double b, double* h);
	/// Sum two expansions, eliminating zero components from the output expansion.
	static int linear_expansion_sum_zeroelim(int elen, double* e, int flen, double* f, double* h);
	/// Sum two expansions.
	static int linear_expansion_sum(int elen, double* e, int flen, double* f, double* h);
	/// Sum two expansions, eliminating zero components from the output expansion.
	static int fast_expansion_sum_zeroelim(int elen, double* e, int flen, double* f, double* h);
	/// Sum two expansions.
	static int fast_expansion_sum(int elen, double* e, int flen, double* f, double* h);
	/// Sum two expansions, eliminating zero components from the output expansion.
	static int expansion_sum_zeroelim2(int elen, double* e, int flen, double* f, double* h);
	/// Sum two expansions, eliminating zero components from the output expansion.
	static int expansion_sum_zeroelim1(int elen, double* e, int flen, double* f, double* h);
	/// Sum two expansions.
	static int expansion_sum(int elen, double*  e, int flen, double* f, double* h);
	/// Add a scalar to an expansion, eliminating zero components from the output expansion.
	static int grow_expansion_zeroelim(int elen, double* e, double b, double* h);
	/// Add a scalar to an expansion.
	static int grow_expansion(int elen, double* e, double b, double* h);
private:
	static InnerData inner;
};

#endif // !defined(GEOMETRICPREDICATES_H__INCLUDED)
