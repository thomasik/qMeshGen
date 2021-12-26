/////////////////////////////////////////////////////////////////////////////
// DMetric3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DMETRIC3D_H__INCLUDED)
#define DMETRIC3D_H__INCLUDED

#include <memory>

#include "DPoint.h"	
#include "DRect.h"
#include "DataVector.h"
#include "DMatrix.h"
#include "DataHashTable.h"

class ControlDataMatrix2d;
class ControlNode3d;
class SurfaceParametric;

class ControlDataStretch3d
{
public:
	/// Constructor
	ControlDataStretch3d(const double _lx, const double _ly, const double _lz, 
		const double _ax= 0.0, const double _ay = 0.0 , const double _az = 0.0) 
		: lx(_lx), ly(_ly), lz(_lz), ax(_ax), ay(_ay), az(_az) {}
	/// Constructor
	ControlDataStretch3d(const double _lxyz) 
		: lx(_lxyz), ly(_lxyz), lz(_lxyz), ax(0.0), ay(0.0), az(0.0) {}
	/// Constructor
	ControlDataStretch3d() : lx(0.0), ly(0.0), lz(0.0), ax(0.0), ay(0.0), az(0.0) {}
	/// Stores the coordinates of the data into the stream
	friend ostream& operator<<(ostream& os, const ControlDataStretch3d& data);
	/// Loads the coordinates of the data from the stream
	friend istream& operator>>(istream& is, ControlDataStretch3d& data);

	/// stretching
	double lx, ly, lz;
	/// angle
	double ax, ay, az;
};

class ControlDataMatrix3d
{
public:
	/// Constructor
	ControlDataMatrix3d(const double _dxx, const double _dyy, const double _dzz,
			const double _dxy, const double _dxz, const double _dyz) 
		: dxx(_dxx), dyy(_dyy), dzz(_dzz), dxy(_dxy), dxz(_dxz), dyz(_dyz) {}
	/// Constructor
	explicit ControlDataMatrix3d(const double d)
		: dxx(d), dyy(d), dzz(d), dxy(0.0), dxz(0.0), dyz(0.0) {
		assert(d >= 0.0); 
	}
	/// Constructor (copying)
	ControlDataMatrix3d(const ControlDataMatrix3d& data) 
		: dxx(data.dxx), dyy(data.dyy), dzz(data.dzz),
			dxy(data.dxy), dxz(data.dxz), dyz(data.dyz) {}
	/// Constructor (empty)
	ControlDataMatrix3d() : dxx(0.0), dyy(0.0), dzz(0.0), dxy(0.0), dxz(0.0), dyz(0.0) {}
	/// Constructor (from eigensystem)
	ControlDataMatrix3d(const DVector3d& e0, const DVector3d& e1, const DVector3d& e2, double d[]);
	/// If zero
	bool isZero() const { return dxx == 0.0 && dyy == 0.0 && dzz == 0.0 && dxy == 0.0 && dxz == 0.0 && dyz == 0.0; }
	/// Clear data (set to 0.0)
	void reset() { dxx = dyy = dzz = dxy = dxz = dyz = 0.0; }
	/// Set to Identity matrix
	void setIdentity() { dxx = dyy = dzz = 1.0; dxy = dxz = dyz = 0.0; }
	/// Is sufficiently close to isotropic
	bool isotropic(double eps = SMALL_NUMBER) const;
	/// Check for Identity matrix
	bool isIdentity() const { return (dxx == 1.0) && (dyy == 1.0) && (dzz == 1.0) && (dxy == 0.0) && (dxz == 0.0) && (dyz == 0.0); }
	/// Count minimum of two metrics
	bool setMinimum(const ControlDataMatrix3d& m, double * d_max = nullptr);
	/// Set from eigensystem
	void setEigensystem(DVector3d e[], double d[]);
	/// Set from eigensystem
	void setEigensystem(const DMatrix3d& e, double d[]);
	/// Count minimum of two metrics
	double countDifferenceRR(const ControlDataMatrix3d& m) const;
	/// Count inverse matrix
	ControlDataMatrix3d inverse() const {
		double f = 1.0 / det();
		return ControlDataMatrix3d(
			f*(dyy*dzz-dyz*dyz), f*(dxx*dzz-dxz*dxz), f*(dxx*dyy-dxy*dxy),
			f*(dxz*dyz-dxy*dzz), f*(dxy*dyz-dxz*dyy), f*(dxz*dxy-dxx*dyz));
	}
	/// Count eigenvector in column 'ie' for eigenvalue d
	void countEigenvector(double d, int ie, DMatrix3d& eigenvectors) const;
	/// Count eigenvalues and eigenvectors
	bool eigensystem(DMatrix3d& eigenvectors, double eigenvalues[]) const;
	/// Count minimum eigenvalue
	double minEigenvalue() const;
	/// Count minimum eigenvalue
	double maxEigenvalue() const;
	/// Count eigenvalues only
	bool eigenvalues(double eigenvalues[]) const;
	/// Calculate metric from simplex
	static ControlDataMatrix3d countMetric(const DPoint3d& pt0, const DPoint3d& pt1, 
		const DPoint3d& pt2, const DPoint3d& pt3);
	/// Calculate metric tensor from simplex
	static ControlDataMatrix3d countMetricTensor(const DPoint3d& pt0, const DPoint3d& pt1, 
		const DPoint3d& pt2, const DPoint3d& pt3);
	/// Transformation matrix to metric tensor
	ControlDataMatrix3d transformationToTensor() const;
	/// Metric tensor to transformation matrix
	ControlDataMatrix3d tensorToTransformation() const;
	/// Count determinant
	double det() const {
		return dxx*dyy*dzz + 2*dxy*dxz*dyz - dyy*dxz*dxz - dzz*dxy*dxy - dxx*dyz*dyz;
	}
	/// Adds two matrices
	ControlDataMatrix3d&	operator+=(const ControlDataMatrix3d& m){ 
		dxx+=m.dxx; dyy+=m.dyy; dzz+=m.dzz; 
		dxy+=m.dxy; dxz+=m.dxz; dyz+=m.dyz; 
		return *this; 
	}
	/// Adds two matrices
	ControlDataMatrix3d	operator+(const ControlDataMatrix3d& m) const { 
		return ControlDataMatrix3d(dxx+m.dxx, dyy+m.dyy, dzz+m.dzz, dxy+m.dxy, dxz+m.dxz, dyz+m.dyz); 
	}
	/// Scales down the control matrix by the given factor
	ControlDataMatrix3d&	operator/=(double d){ 
		assert(d!=0.0); dxx/=d; dyy/=d; dzz/=d; dxy/=d; dxz/=d; dyz/=d; 
		return *this; 
	}
	bool operator==(const ControlDataMatrix3d& m) const {
		return dxx == m.dxx && dyy == m.dyy && dzz == m.dzz &&
			dxy == m.dxy && dxz == m.dxz && dyz == m.dyz;
	}
	/// Scales the control matrix by the given factor
	ControlDataMatrix3d&	operator*=(double d){ 
		dxx*=d; dyy*=d; dzz*=d; dxy*=d; dxz*=d; dyz*=d; 
		return *this; 
	}
	/// Returns the control matrix scaled by the given factor
	ControlDataMatrix3d	operator*(double d) const { 
		return ControlDataMatrix3d(dxx*d, dyy*d, dzz*d, dxy*d, dxz*d, dyz*d); 
	}
	/// Product of two matrices
	DMatrix3d operator*(const ControlDataMatrix3d& m) const {
		return DMatrix3d(
			dxx*m.dxx + dxy*m.dxy + dxz*m.dxz,
			dxx*m.dxy + dxy*m.dyy + dxz*m.dyz,
			dxx*m.dxz + dxy*m.dyz + dxz*m.dzz,

			dxy*m.dxx + dyy*m.dxy + dyz*m.dxz,
			dxy*m.dxy + dyy*m.dyy + dyz*m.dyz,
			dxy*m.dxz + dyy*m.dyz + dyz*m.dzz,

			dxz*m.dxx + dyz*m.dxy + dzz*m.dxz,
			dxz*m.dxy + dyz*m.dyy + dzz*m.dyz,
			dxz*m.dxz + dyz*m.dyz + dzz*m.dzz);
	}
	/// Product of matrix and vector
	DPoint3d operator*(const DPoint3d& pt) const {
		return DPoint3d(dxx * pt.x + dxy * pt.y + dxz * pt.z, 
			dxy * pt.x + dyy * pt.y + dyz * pt.z,
			dxz * pt.x + dyz * pt.y + dzz * pt.z);
	}
	DMPoint3d multiplyRStoMS(const DPoint3d& pt) const {
		return DMPoint3d(dxx * pt.x + dxy * pt.y + dxz * pt.z, 
			dxy * pt.x + dyy * pt.y + dyz * pt.z,
			dxz * pt.x + dyz * pt.y + dzz * pt.z);
	}
	DPoint3d multiplyMStoRS(const DMPoint3d& pt) const {
		return DPoint3d(dxx * pt.x + dxy * pt.y + dxz * pt.z, 
			dxy * pt.x + dyy * pt.y + dyz * pt.z,
			dxz * pt.x + dyz * pt.y + dzz * pt.z);
	}
	/// Product of matrix and vector
	DVector3d operator*(const DVector3d& v) const {
		return DVector3d(dxx * v.x + dxy * v.y + dxz * v.z, 
			dxy * v.x + dyy * v.y + dyz * v.z,
			dxz * v.x + dyz * v.y + dzz * v.z);
	}
	DMVector3d multiplyRStoMS(const DVector3d& v) const {
		return DMVector3d(dxx * v.x + dxy * v.y + dxz * v.z, 
			dxy * v.x + dyy * v.y + dyz * v.z,
			dxz * v.x + dyz * v.y + dzz * v.z);
	}
	DVector3d multiplyMStoRS(const DMVector3d& v) const {
		return DVector3d(dxx * v.x + dxy * v.y + dxz * v.z, 
			dxy * v.x + dyy * v.y + dyz * v.z,
			dxz * v.x + dyz * v.y + dzz * v.z);
	}
	/// Returns the control matrix scaled by the given factor
	ControlDataMatrix3d	operator/(double d) const { 
		return ControlDataMatrix3d(dxx/d, dyy/d, dzz/d, dxy/d, dxz/d, dyz/d); 
	}

	DRect prepareGraphics(DataVector<DPoint2d> & pts0, DataVector<DPoint2d> & pts1, DataVector<bool> & front) const;
	static void fixEigensystem(DMatrix3d& e, double d[]);
	static void switchEigenvectors(DMatrix3d& e, double d[], int k0, int k1);
	static void drawInterpolationEllipsoids(const ControlDataMatrix3d& cdm1, const ControlDataMatrix3d& cdm2);
	static void drawIntersectionEllipsoids(const ControlDataMatrix3d& cdm1, const ControlDataMatrix3d& cdm2);

	bool storeSimple(ostream& os, char end_char = '\n') const;
	bool readSimple(istream& is);

	double getMaxHRatioIsotropic(const ControlDataMatrix3d & cdm1, double * h_ratio0_ret = nullptr, double * h_ratio1_ret = nullptr) const;
	double getMaxHRatioDirect(const ControlDataMatrix3d & cdm1, double * h_ratio0_ret = nullptr, double * h_ratio1_ret = nullptr) const;
	double getMaxHRatioReduction(const ControlDataMatrix3d & cdm1, double * h_ratio0_ret = nullptr, double * h_ratio1_ret = nullptr) const;

	/// Log
	friend ostream& operator<<(ostream& os, const ControlDataMatrix3d& cdm){
		return os << '[' << cdm.dxx << ',' << cdm.dyy << ',' << cdm.dzz << ',' <<
			cdm.dxy << ',' << cdm.dxz << ',' << cdm.dyz << ']';
	}
	friend istream& operator>>(istream& is, ControlDataMatrix3d& cdm);
	/// Second derivatives
	double dxx, dyy, dzz;
	double dxy, dxz, dyz;
	/// identity control data matrix
	static ControlDataMatrix3d identity;
	/// temporary counter
	static unsigned int m_intersection_counter;
};

double diffKdValue(const ControlDataMatrix3d& cdm1, const ControlDataMatrix3d& cdm2);

typedef ControlDataStretch3d CDS3d;
typedef ControlDataMatrix3d CDM3d;

/**
 * This class is responsible for Riemman-metric transformations
 *  and managing the transformation matrix.
 */
class DMetric3d  
{
public:
	/// Standard constructor
	DMetric3d() { setToIdentity(); }
	DMetric3d(const ControlDataMatrix3d& cdm) { setData(cdm); }
public:
	/// Transforms 2D metric transformation tensor into 3D
	static bool projectCDMto3D(const ControlDataMatrix2d& cdm2d, 
		std::shared_ptr<const SurfaceParametric> surface, const DPoint2d& pt2d,
		ControlDataMatrix3d& cdm, DVector3d * cdm_e0 = 0, DVector3d * cdm_e1 = 0);
	/// Transforms explicite data into metric transformation tensor
	static ControlDataMatrix3d stretchToMatrix(const ControlDataStretch3d& stretch);
	/// Transforms explicite data into metric transformation tensor
	static ControlDataMatrix3d stretchToMatrixWithAdjust(const ControlDataStretch3d& stretch);
	/// Sets new data for the transformation matrix
//	void setData(const ControlDataStretch3d& data);
	/// Sets new data for the transformation matrix
	void setData(const ControlDataMatrix3d& data);
	/// Calculates coordinates of the 3D point via reverse transformation (using current metric)
	DMPoint3d transformRStoMS(const DPoint3d& pt) const;
	/// Calculates coordinates of the 3D point via direct transformation (using current metric)
	DPoint3d transformMStoRS(const DMPoint3d& pt) const;
	/// Calculates coordinates of the 3D vector via reverse transformation (using current metric)
	DMVector3d transformRStoMS(const DVector3d& pt) const;
	/// Calculates coordinates of the 3D vector via direct transformation (using current metric)
	DVector3d transformMStoRS(const DMVector3d& pt) const;
	/// Resets the transformation matrix (sets it to identity)
	void setToIdentity();
	/// Returns maximumem length
	double maxLength() const { return m_matrix.maxEigenvalue(); }
	/// Returns minimumem length
	double minLength() const { return m_matrix.minEigenvalue(); }

	static double getMetricGradation(const CDM3d& cdm0, const CDM3d& cdm1, const DVector3d& dv,
		double * h_ratio0 = nullptr, double * h_ratio1 = nullptr, double diff_eps = 1e-2);
	static bool_pair adjustMetricGradation(ControlNode3d * cn0, ControlNode3d * cn1);
	static bool adjustMetricGradationNb(ControlNode3d* cn0);
	static bool adjustMetricGradationIterativeNb(ControlNode3d * cn0, int step);
public:
	/// Adjust metric lengths according to min/max/stretch ratio
	static void adjustLengths(double & lx, double & ly, double & lz);
protected:
	/// Matrix for direct transformation
	ControlDataMatrix3d m_matrix;
	/// Matrix for reverse transformation
	ControlDataMatrix3d m_reverse;
};

class ControlDataExtMatrix3d
{
public:
	ControlDataExtMatrix3d(double r1, const ControlDataMatrix3d& d1,
			double r2, const ControlDataMatrix3d& d2)  
		: radius1(r1), radius2(r2), data1(d1), data2(d2) {}
	ControlDataExtMatrix3d(double r1, const ControlDataMatrix3d& d1)  
		: radius1(r1), radius2(0.0), data1(d1), data2(d1) {}
	virtual ~ControlDataExtMatrix3d() {} // empty virtual destructor
public:
	virtual bool isPointWithin(const DPoint3d& pt) const = 0;
	virtual ControlDataMatrix3d getControlDataMatrix3d(const DPoint3d& pt) const = 0;
	double totalRadius() const { return radius1 + radius2; }
	double getRadius1() const { return radius1; }
	double getRadius2() const { return radius2; }
	ControlDataMatrix3d getInnerData() const { return data1; }
	ControlDataMatrix3d getOuterData() const { return data2; }
protected:
	double radius1, radius2;	// in real space
	ControlDataMatrix3d data1, data2;
};

class ControlDataExtMatrix3dSphere : public ControlDataExtMatrix3d
{
public:
	ControlDataExtMatrix3dSphere(const DPoint3d& pt, double r1, const ControlDataMatrix3d& d1,
		double r2, const ControlDataMatrix3d& d2);
	ControlDataExtMatrix3dSphere(const DPoint3d& pt, double r1, const ControlDataMatrix3d& d1);
public:
	bool isPointWithin(const DPoint3d& pt) const;
	ControlDataMatrix3d getControlDataMatrix3d(const DPoint3d& pt) const;
	const DPoint3d& getMiddle() const { return middle; }
protected:
	DPoint3d middle;
};

class ControlDataExtMatrix3dSegment : public ControlDataExtMatrix3d
{
public:
	ControlDataExtMatrix3dSegment(const DPoint3d& pt, const DVector3d& dv, 
		double r1, const ControlDataMatrix3d& d1,
		double r2, const ControlDataMatrix3d& d2);
	ControlDataExtMatrix3dSegment(const DPoint3d& pt, const DVector3d& dv, 
		double r1, const ControlDataMatrix3d& d1);
public:
	bool isPointWithin(const DPoint3d& pt) const;
	ControlDataMatrix3d getControlDataMatrix3d(const DPoint3d& pt) const;
	DPoint3d getMiddle() const { return pt0+dv[0]*0.5; }
	const DVector3d& getDv() const { return dv[0]; }
	const DPoint3d& getStart() const { return pt0; }
protected:
	DPoint3d pt0;
	DVector3d dv[3];
	double w, dr;
	DMatrix3d dm;
};

class ControlDataExtMatrix3dTriangle : public ControlDataExtMatrix3d
{
public:
	ControlDataExtMatrix3dTriangle(const DPoint3d& pt, 
		const DVector3d& dv1, const DVector3d& dv2,
		double r1, const ControlDataMatrix3d& d1,
		double r2, const ControlDataMatrix3d& d2);
	ControlDataExtMatrix3dTriangle(const DPoint3d& pt, 
		const DVector3d& dv1, const DVector3d& dv2,
		double r1, const ControlDataMatrix3d& d1);
public:
	bool isPointWithin(const DPoint3d& pt) const;
	ControlDataMatrix3d getControlDataMatrix3d(const DPoint3d& pt) const;
	DPoint3d getMiddle() const { return pt0 + dv[0]*0.5 + dv[1]*0.5; }
	const DVector3d& getDv(int i) const { return dv[i]; }
	const DPoint3d& getStart() const { return pt0; }
protected:
	DPoint3d pt0;
	DVector3d dv[3];
	double w, dr;
	DMatrix3d dm;
};

#endif // !defined(DMETRIC3D_H__INCLUDED)
