/////////////////////////////////////////////////////////////////////////////
// Curve2dBSpline.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVEBSPLINE_H__INCLUDED)
#define CURVEBSPLINE_H__INCLUDED

#include "DPoint.h"
#include "Curve2dParametric.h"
#include "DataVector.h"

class DMetric2d;

/**
 * This class implements a B-spline curve in 2D
 *	with several required methods.
 */
class Curve2dBSpline : public Curve2dParametric
{
public:
	/// Standard constructor
	Curve2dBSpline(int ct, const DPoint2d points[], bool opened = true);
	/// Standard constructor
	Curve2dBSpline(const DataVector<DPoint2d> & points, bool opened = true);
	/// Standard constructor
	Curve2dBSpline() : Curve2dParametric(0.0, 1.0), m_opened(true) {}
	/// Copying constructor
	Curve2dBSpline(const Curve2dBSpline& spline);
public:
	virtual string getSimpleDescription() const override { return m_opened ? "bspline2d-opened" : "bspline2d-closed"; }
	/// for storing in XML file as a name of the curve-type
	virtual string typenameXML() const { return "bspline"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return CURVE_BSPLINE; }
	/// Returns the point of this curve for the given parameter
	virtual DPoint2d getPoint(double t) const;
	/// Returns a parameter t for the given point 2D (t > ts), in range
	virtual double getParameterInRange( const DPoint2d& pt, double ts, double t_min, double t_max ) const;
	/// Returns the derivatives (dx/dt, dy/dt) of the curve for the given parameter
	virtual DVector2d getDerivative(double t) const;
	/// Returns the second derivatives (d2x/dt2, d2y/dt2) of the curve for the given parameter
	virtual DVector2d getSecondDerivative(double t) const;
	/// Returns the third derivatives (d3x/dt3, d3y/dt3) of the curve for the given parameter
	virtual DVector2d getThirdDerivative(double t) const;
	/// Checks wheterh this curve is valid (i.e. is properly initialized)
	virtual bool isValid() const { return m_points.countInt() > 1; }
	/// Stores the parameters of this curve into the stream
	friend ostream& operator<<(ostream& os, const Curve2dBSpline* spline);
	/// Restores the parameters of this curve from the stream
	friend istream& operator>>(istream& is, Curve2dBSpline* spline);
public:
	/// Returns maximum t-parameter for this curve
	double getMaxParameter() const { return m_opened? (double)(m_points.countInt()-1) : (double)m_points.countInt(); }
public:
	/// fit to point-cloud using direct fit, return max-dist
	double fitToPoints(const DataVector<DPoint2d> & points, 
		const DPoint2d& ptA, const DPoint2d& ptB, const DMetric2d& dm);
	/// adaptive subroutinge - fit to point-cloud using direct fit, return number of inserted points
	int fitToPointsAdaptive(const DataVector<DPoint2d> & points, int ind, const DMetric2d& dm);
	/// fit to point-cloud using direct fit, return max-dist
	void fitThroughPoints(const DataVector<DPoint2d> & points, 
		const DPoint2d& ptA, const DPoint2d& ptB, const DMetric2d& dm);
	/// adaptive subroutinge - fit to point-cloud using direct fit, return number of inserted points
	int fitThroughPointsAdaptive(const DataVector<DPoint2d> & points, int ind, const DMetric2d& dm);
private:
	/// max distance
	double maxDistanceForSubSegment(const DataVector<DPoint2d>& points, int ind, const DMetric2d& dm) const;
	/// max distance
	double maxDistanceForSpline(const DataVector<DPoint2d>& points, const DMetric2d& dm) const;
private:
	/// Shape functions coefficients for B-Spline
	static double BSplineMatrix[4][4];
	/// First derivative shape functions coefficients for B-Spline
	static double BDerivativeSplineMatrix[3][4];
	/// Second derivative shape functions coefficients for B-Spline
	static double BDDerivativeSplineMatrix[2][4];
	/// Third derivative shape functions coefficients for B-Spline
	static double BDDDerivativeSplineMatrix[4];
	/// Make friends! 
	friend class SurfaceBSplineCylindrical;	// - since it's using the same equations anyway...
	friend class SurfaceBSplinePlanar;		// - since it's using the same equations anyway...
	friend class Curve3dBSpline;			// - since it's using the same equations anyway...
private:
	/// Count functional knots for this B-spline
	bool countNodes();
	/// Returns the number of B-Spline segment and sets t as a parameter inside this segment
	int countSplineIndices(double& t) const;
private:
	/// Geometric nodes [0..pts_ct+1]
	DataVector<DPoint2d> m_nodes;
	/// Geometric points [1..pts_ct]
	DataVector<DPoint2d> m_points;
	/// Whether the B-spline is opened
	bool m_opened;
};

#endif // !defined(CURVEBSPLINE_H__INCLUDED)
