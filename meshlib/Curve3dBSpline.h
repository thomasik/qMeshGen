/////////////////////////////////////////////////////////////////////////////
// Curve3dBSpline.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CURVE3DBSPLINE_H__INCLUDED)
#define CURVE3DBSPLINE_H__INCLUDED

#include "DPoint.h"
#include "Curve3dParametric.h"
#include "DataVector.h"

class DMetric3d;
class Metric3dContext;

/**
 * This class implements a B-spline curve in 3D
 *	with several required methods.
 */
class Curve3dBSpline : public Curve3dParametric
{
public:
	/// Standard constructor
	Curve3dBSpline(int ct, const DPoint3d points[], bool opened = true );
	/// Standard constructor
	Curve3dBSpline(const DataVector<DPoint3d> & points, bool opened = true );
	/// Copying constructor
	Curve3dBSpline(const Curve3dBSpline& spline);
public:
	virtual string getSimpleDescription() const override { return m_opened ? "bspline3d-opened" : "bspline3d-closed"; }
	/// for storing in XML file as a name of the curve-type
	virtual string typenameXML() const { return "bspline3d-curve"; }
	/// Store XML description to stream
	virtual ostream& storeXML(ostream& os, const string& prefix = "") const;
	/// Returns the object-specific type
	virtual ElementType getType() const { return CURVE_3D_BSPLINE; }
	/// Returns the point of this curve for the given parameter
	virtual DPoint3d getPoint(double t) const;
	/// Returns a parameter t for the given point 3D (t > ts)
	virtual double getParameterInRange(const DPoint3d& pt, double ts, double t_min, double t_max) const override;
	/// Returns a parameter t for the given point 3D (t > ts)
	virtual double getParameter(const DPoint3d& pt, double ts) const override;
	/// Returns the derivatives (dx/dt, dy/dt) of the curve for the given parameter
	virtual DVector3d getDerivative(double t) const;
	/// Returns the second derivatives (d2x/dt2, d2y/dt2) of the curve for the given parameter
	virtual DVector3d getSecondDerivative(double t) const;
	/// Returns the third derivatives (d3x/dt3, d3y/dt3) of the curve for the given parameter
	virtual DVector3d getThirdDerivative(double t) const;
	/// Checks wheterh this curve is valid (i.e. is properly initialized)
	virtual bool isValid() const { return m_points.countInt() > 1; }
	/// Stores the parameters of this curve into the stream
	friend ostream& operator<<(ostream& os, const Curve3dBSpline* spline);
	/// Restores the parameters of this curve from the stream
	friend istream& operator>>(istream& is, Curve3dBSpline* spline);
public:
	/// Returns maximum t-parameter for this curve
	double getMaxParameter() const { return m_opened? (double)(m_points.countInt()-1) : (double)(m_points.countInt()); }
public:
	/// fit to point-cloud using direct fit
	static Curve3dPtr fitToPointSequence(const DataVector<DPoint3d> & points, 
		Metric3dContext & mc, const double& tolerance );
	/// fit to point-cloud using direct fit
	void fitToPointSequenceAdaptive(const DataVector<DPoint3d> & points, 
		Metric3dContext & mc, const double& tolerance );
	/// fit to point-cloud using direct fit
	static Curve3dPtr fitToPoints(const DataVector<DPoint3d> & points, 
		const DPoint3d& ptA, const DPoint3d& ptB);
	/// adaptive subroutine - fit to point-cloud using direct fit, return number of inserted points
	int fitToPointsAdaptive(const DataVector<DPoint3d> & points, int ind);
	/// fit to point-cloud using direct fit
	static Curve3dPtr fitThroughPoints(const DataVector<DPoint3d> & points, 
		const DPoint3d& ptA, const DPoint3d& ptB);
	/// adaptive subroutine - fit to point-cloud using direct fit, return number of inserted points
	int fitThroughPointsAdaptive(const DataVector<DPoint3d> & points, int ind);
private:
	/// max distance
	double maxDistanceForSubSegment(const DataVector<DPoint3d>& points, int ind, const DMetric3d& dm) const;
	/// max distance
	double maxDistanceForSpline(const DataVector<DPoint3d>& points, const DMetric3d& dm) const;
private:
	/// Count functional knots for this B-spline
	bool countNodes();
	/// Returns the number of B-Spline segment and sets t as a parameter inside this segment
	int countSplineIndices(double& t) const;
private:
	/// Geometric nodes [0..pts_ct+1]
	DataVector<DPoint3d> m_nodes;
	/// Geometric points [1..pts_ct]
	DataVector<DPoint3d> m_points;
	/// Whether the B-spline is opened
	bool m_opened;
};

#endif // !defined(CURVE3DBSPLINE_H__INCLUDED)
