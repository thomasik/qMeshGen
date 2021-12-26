/////////////////////////////////////////////////////////////////////////////
// DMetric2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DMETRIC2D_H__INCLUDED)
#define DMETRIC2D_H__INCLUDED

#include <memory>

#include "DPoint.h"	
#include "DMatrix.h"
#include "DataVector.h"

class SurfaceParametric;
class Curve2dParametric;

class ControlDataStretch2d
{
public:
	/// Constructor
	ControlDataStretch2d(const double _lx, const double _ly, const double _angle = 0.0)
		: lx(_lx), ly(_ly), angle(_angle) {}
	/// Constructor
	ControlDataStretch2d() : lx(0.0), ly(0.0), angle(0.0) {}
	/// Stores the coordinates of the data into the stream
	friend ostream& operator<<(ostream& os, const ControlDataStretch2d& data);
	/// Loads the coordinates of the data from the stream
	friend istream& operator>>(istream& is, ControlDataStretch2d& data);

	/// stretching
	double lx, ly;
	/// angle
	double angle;
};

class ControlDataMatrix2d
{
public:
	/// Constructor
	ControlDataMatrix2d(const double _m11, const double _m12, const double _m22) 
		: m11(_m11), m12(_m12), m22(_m22) { }
	/// Constructor
	ControlDataMatrix2d(const ControlDataMatrix2d& data) : m11(data.m11), m12(data.m12), m22(data.m22) {}
	/// Constructor
	ControlDataMatrix2d() : m11(0.0), m12(0.0), m22(0.0) {}
	/// Constructor
	ControlDataMatrix2d(const DPoint2d& pt0, const DPoint2d& pt1, const DPoint2d& pt2);
	/// Constructor (from eigensystem)
	ControlDataMatrix2d(const DVector2d& e0, const DVector2d& e1, double d[]);
	/// Clear data (set to 0.0)
	void reset() {m11 = m12 = m22 = 0.0; }
	/// Set to identity
	void setIdentity() { m11 = m22 = 1.0; m12 = 0.0; }
	/// Count minimum of two metrics (using simultaneous reduction, returns true if any changes made)
	bool setMinimum(const ControlDataMatrix2d& m);
	/// Count inverse matrix
	ControlDataMatrix2d inverse() const {
		double f = 1.0 / (m11 * m22 - m12*m12);
		return ControlDataMatrix2d(f*m22, -f*m12, f*m11);
	}
	/// Count minimum eigenvalue
	double minEigenvalue() const;
	/// Count eigenvalues only
	bool eigenvalues(double eigenvalues[]) const;
	/// Count eigenvalues and eigenvectors
	bool eigensystem(DMatrix2d& eigenvectors, double eigenvalues[]) const;
	/// Calculate metric from simplex
	static ControlDataMatrix2d countMetric(const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2,
		DVector3d& e0, DVector3d& e1);
	/// Calculate metric tensor from simplex
	static ControlDataMatrix2d countMetricTensor(const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2,
		DVector3d& e0, DVector3d& e1);
	/// Calculate metric from simplex
	static ControlDataMatrix2d countMetric(const DPoint2d& pt0, const DPoint2d& pt1, const DPoint2d& pt2);
	/// Calculate metric tensor from simplex
	static ControlDataMatrix2d countMetricTensor(const DPoint2d& pt0, const DPoint2d& pt1, const DPoint2d& pt2);
	/// Transformation matrix to metric tensor
	ControlDataMatrix2d transformationToTensor() const { return this->squared(); }
	/// Metric tensor to transformation matrix
	ControlDataMatrix2d tensorToTransformation() const;
	/// Count det
	double det() const { return m11*m22 - m12*m12; }
	/// Count difference of two metrics
	double countDifferenceRR(const ControlDataMatrix2d& data) const;
	/// Multiplies matrix and vector
	DPoint2d operator*(const DPoint2d& pt) const {
		return DPoint2d(m11*pt.x+m12*pt.y, m12*pt.x+m22*pt.y);
	}
	/// Multiplies matrix and vector
	DVector2d operator*(const DVector2d& v) const {
		return DVector2d(m11*v.x+m12*v.y, m12*v.x+m22*v.y);
	}
	/// Adds two matrices
	ControlDataMatrix2d&	operator+=(const ControlDataMatrix2d& m){ 
		m11+=m.m11; m12+=m.m12; m22+=m.m22; return *this; }
	/// Adds two matrices
	ControlDataMatrix2d	operator+(const ControlDataMatrix2d& m){ 
		return ControlDataMatrix2d(m11+m.m11, m12+m.m12, m22+m.m22); }
	/// Scales down the control matrix by the given factor
	ControlDataMatrix2d&	operator/=(double d){ assert(d!=0.0); m11/=d; m12/=d; m22/=d; return *this; }
	/// Returns the control matrix scaled by the given factor
	ControlDataMatrix2d	operator*(double d) const { return ControlDataMatrix2d(d*m11, d*m12, d*m22); }
	/// Returns squared control matrix
	ControlDataMatrix2d	squared() const { 
		return ControlDataMatrix2d(m11*m11 + m12*m12, m11*m12 + m12*m22, m12*m12 + m22*m22); }
	/// Returns product of two control matrices
	DMatrix2d	operator*(const ControlDataMatrix2d& cdm) const { 
		return DMatrix2d(m11*cdm.m11 + m12*cdm.m12, m11*cdm.m12 + m12*cdm.m22, 
			m12*cdm.m11 + m22*cdm.m12, m12*cdm.m12 + m22*cdm.m22); }
	/// Scales the control matrix by the given factor
	ControlDataMatrix2d&	operator*=(double d){ m11*=d; m12*=d; m22*=d; return *this; }
	/// Returns the control matrix scaled by the given factor
	ControlDataMatrix2d	operator/(double d) const { return ControlDataMatrix2d(m11/d, m12/d, m22/d); }

	/// matrix elements
	double m11, m12, m22;
	/// identity
	static const ControlDataMatrix2d IDENTITY;

	friend std::ostream& operator<<(std::ostream& os, const ControlDataMatrix2d& cdm) {
		return os << "[" << cdm.m11 << ", " << cdm.m12 << ", " << cdm.m22 << "]";
	}
};

double diffKdValue(const ControlDataMatrix2d& cdm1, const ControlDataMatrix2d& cdm2);

typedef ControlDataMatrix2d CDM2d;
/**
 * This class is responsible for Riemman-metric transformations
 *  and managing the transformation matrix.
 */
class DMetric2d  
{
public:
	/// Standard constructor (parameteric matrix)
	
	DMetric2d(std::shared_ptr<const SurfaceParametric> surface, const DPoint2d& pt);
	DMetric2d(const SurfaceParametric* surface, const DPoint2d& pt);
	/// Standard constructor (sizing matrix)
	DMetric2d(const ControlDataMatrix2d& scdm);
	/// Standard constructor (sizing matrix + parametric matrix)
	DMetric2d(const ControlDataMatrix2d& scdm, const ControlDataMatrix2d& _pcdm);
	/// Standard constructor (both)
	DMetric2d(const ControlDataMatrix2d& scdm, 
		const SurfaceParametric* surface, const DPoint2d& pt);
	DMetric2d(const ControlDataMatrix2d& scdm,
		std::shared_ptr<const SurfaceParametric> surface, const DPoint2d& pt);
	/// Standard constructor (empty, identity)
	DMetric2d();
public:
	/// Transforms explicite data (width, height, angle) into metric transformation tensor
	static ControlDataMatrix2d stretchToMatrix(const ControlDataStretch2d& stretch);
	/// Transforms explicite data (width, height, angle) into metric transformation tensor
	static ControlDataMatrix2d stretchToMatrixWithAdjust(const ControlDataStretch2d& stretch);
	/// Transforms metric transformation tensor into explicite data (width, height, angle)
	static ControlDataStretch2d matrixToStretch(const ControlDataMatrix2d& matrix);
	/// Calculates sqrt from symmetric matrix form (i.e. MsT*Ms = M)
	static ControlDataMatrix2d matrixSquareRoot(const ControlDataMatrix2d& cdm);
	/// Sets new data for the transformation matrix
	void setData(const ControlDataMatrix2d& s_cdm, const ControlDataMatrix2d& p_cdm);
	/// Calculates coordinates of the 2D point via reverse transformation (using current metric)
	DMPoint2d transformPStoMS(const DPoint2d& pt) const;
	/// Calculates coordinates of the 2D point via reverse transformation (using current metric)
	DPoint2d transformPStoRS(const DPoint2d& pt) const;
	/// Calculates coordinates of the 2D point via reverse transformation (using current metric)
	DPoint2d transformRStoPS(const DPoint2d& pt) const;
	/// Calculates coordinates of the 2D point via direct transformation (using current metric)
	DPoint2d transformMStoPS(const DMPoint2d& pt) const;
	/// Calculates coordinates of the 2D vector via reverse transformation (using current metric)
	DMVector2d transformPStoMS(const DVector2d& v) const;
	/// Calculates coordinates of the 2D vector via reverse transformation (using current metric)
	DVector2d transformPStoRS(const DVector2d& v) const;
	/// Calculates coordinates of the 2D vector via reverse transformation (using current metric)
	DVector2d transformRStoPS(const DVector2d& v) const;
	/// Calculates coordinates of the 2D vector via direct transformation (using current metric)
	DVector2d transformMStoPS(const DMVector2d& v) const;
	/// Resets the transformation matrix (sets it to identity)
	void setToIdentity();
	/// Returns metric tensor for current metric
	ControlDataMatrix2d getMetricTensor() const {
		DMatrix2d mm = m * m;
		return ControlDataMatrix2d(mm.m[0][0], m.m[0][1], m.m[1][1]);
	}
	double detPStoMS() const { return mr.det(); }
public:
	/// Adjust metric lengths according to min/max/stretch ratio
	static void adjustLengths(double & lx, double & ly);
protected:
	/// Matrix for direct transformation
	DMatrix2d m;
	/// Matrix for reverse transformation
	DMatrix2d mr;
	/// Data for parametric transformation
	ControlDataMatrix2d pcdm;
};

class ControlDataExtMatrix2d
{
public:
	ControlDataExtMatrix2d(double r1, const ControlDataMatrix2d& d1,
			double r2, const ControlDataMatrix2d& d2)  
		: radius1(r1), radius2(r2), data1(d1), data2(d2) {}
	ControlDataExtMatrix2d(double r1, const ControlDataMatrix2d& d1)  
		: radius1(r1), radius2(0.0), data1(d1), data2(d1) {}
	virtual ~ControlDataExtMatrix2d() {}
public:
	virtual bool isPointWithin(const DPoint2d& pt) const = 0;
	virtual ControlDataMatrix2d getControlDataMatrix(const DPoint2d& pt) const = 0;
	double totalRadius() const { return radius1 + radius2; }
	double getRadius1() const { return radius1; }
	double getRadius2() const { return radius2; }
	const ControlDataMatrix2d& getInnerData() const { return data1; }
	const ControlDataMatrix2d& getOuterData() const { return data2; }
protected:
	double radius1, radius2;	// in real space
	ControlDataMatrix2d data1, data2;
};

class ControlDataExtMatrix2dRadial : public ControlDataExtMatrix2d
{
public:
	ControlDataExtMatrix2dRadial(const DPoint2d& pt, double r1, const ControlDataMatrix2d& d1,
		double r2, const ControlDataMatrix2d& d2, std::shared_ptr<const SurfaceParametric> surface);
	ControlDataExtMatrix2dRadial(const DPoint2d& pt, double r1, const ControlDataMatrix2d& d1,
		std::shared_ptr<const SurfaceParametric> surface);
public:
	bool isPointWithin(const DPoint2d& pt) const;
	ControlDataMatrix2d getControlDataMatrix(const DPoint2d& pt) const;
	const DPoint2d& getMiddle() const { return middle; }
protected:
	DPoint2d middle;			// in parametric space
	DMetric2d dmp;
};

class ControlDataExtMatrix2dSegment : public ControlDataExtMatrix2d
{
public:
	ControlDataExtMatrix2dSegment(const DPoint2d& pt, const DVector2d& dv, 
		double r1, const ControlDataMatrix2d& d1,
		double r2, const ControlDataMatrix2d& d2, 
		std::shared_ptr<const SurfaceParametric> surface);
	ControlDataExtMatrix2dSegment(const DPoint2d& pt, const DVector2d& dv, 
		double r1, const ControlDataMatrix2d& d1, 
		std::shared_ptr<const SurfaceParametric> surface);
public:
	bool isPointWithin(const DPoint2d& pt) const;
	ControlDataMatrix2d getControlDataMatrix(const DPoint2d& pt) const;
	double getLocalDistance(const DPoint2d& pt) const;
	DPoint2d getMiddle() const { return pt0+dv*0.5; }
	const DVector2d& getDv() const { return dv; }
	const DPoint2d& getStart() const { return pt0; }
protected:
	DPoint2d pt0;
	DVector2d dv;
	DVector2d nv;
	double w, dr;
};

class ControlDataExtMatrix2dCurve : public ControlDataExtMatrix2d
{
public:
	ControlDataExtMatrix2dCurve(std::shared_ptr<const Curve2dParametric> fig, double t0, double t1,
		double r1, const ControlDataMatrix2d& d1,
		double r2, const ControlDataMatrix2d& d2, std::shared_ptr<const SurfaceParametric> surface);
	ControlDataExtMatrix2dCurve(std::shared_ptr<const Curve2dParametric> fig, double t0, double t1,
		double r1, const ControlDataMatrix2d& d1, std::shared_ptr<const SurfaceParametric> surface);
public:
	bool isPointWithin(const DPoint2d& pt) const;
	ControlDataMatrix2d getControlDataMatrix(const DPoint2d& pt) const;
protected:
//	Curve2dConstPtr curve;
//	double t0, t1;
public:
	DataVector<std::shared_ptr<ControlDataExtMatrix2dSegment>> poly_control;
};

#endif // !defined(DMETRIC2D_H__INCLUDED)
