/////////////////////////////////////////////////////////////////////////////
// DMatrix.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(DMATRIX_H__INCLUDED)
#define DMATRIX_H__INCLUDED

#include "DPoint.h"
#include "DVector.h"

class DMatrix2d
{
public:
	DMatrix2d(double m11 = 0.0, double m12 = 0.0, double m21 = 0.0, double m22 = 0.0){
		m[0][0] = m11; m[0][1] = m12; m[1][0] = m21; m[1][1] = m22;
	}
	double& operator()(int i, int j) { return m[i][j]; }
	const double& operator()(int i, int j) const { return m[i][j]; }
	double det() const { return m[0][0]*m[1][1] - m[0][1]*m[1][0]; }
	/// Count solution for linear equation
	bool solve(const DVector2d& b, DVector2d& x) const;
	void setIdentity() { m[0][0] = m[1][1] = 1.0; m[0][1] = m[1][0] = 0.0; }
	DMatrix2d inverse() const{
		double f = 1.0 / (m[0][0] * m[1][1] - m[0][1]*m[1][0]);
		return DMatrix2d(f*m[1][1], -f*m[0][1], -f*m[1][0], f*m[0][0]);
	}
	DMatrix2d transposed() const{
		return DMatrix2d(m[0][0], m[1][0], m[0][1], m[1][1]);
	}
	DMatrix2d operator*(const DMatrix2d& dm) const {
		return DMatrix2d(m[0][0]*dm.m[0][0]+m[0][1]*dm.m[1][0], 
			m[0][0]*dm.m[0][1]+m[0][1]*dm.m[1][1], 
			m[1][0]*dm.m[0][0]+m[1][1]*dm.m[1][0], 
			m[1][0]*dm.m[0][1]+m[1][1]*dm.m[1][1]);
	}
	DMatrix2d operator+(const DMatrix2d& dm) const {
		return DMatrix2d(m[0][0] + dm.m[0][0], m[0][1] + dm.m[0][1], 
			m[1][0] + dm.m[1][0], m[1][1] + dm.m[1][1]);
	}
	DPoint2d operator*(const DPoint2d& pt) const {
		return DPoint2d(m[0][0] * pt.x + m[0][1] * pt.y, 
			m[1][0] * pt.x + m[1][1] * pt.y);
	}
	DPoint2d multiplyMStoPS(const DMPoint2d& pt) const {
		return DPoint2d(m[0][0] * pt.x + m[0][1] * pt.y, 
			m[1][0] * pt.x + m[1][1] * pt.y);
	}
	DMPoint2d multiplyPStoMS(const DPoint2d& pt) const {
		return DMPoint2d(m[0][0] * pt.x + m[0][1] * pt.y, 
			m[1][0] * pt.x + m[1][1] * pt.y);
	}
	DVector2d operator*(const DVector2d& v) const {
		return DVector2d(m[0][0] * v.x + m[0][1] * v.y, 
			m[1][0] * v.x + m[1][1] * v.y);
	}
	DMVector2d multiplyPStoMS(const DVector2d& v) const {
		return DMVector2d(m[0][0] * v.x + m[0][1] * v.y, 
			m[1][0] * v.x + m[1][1] * v.y);
	}
	DVector2d multiplyMStoPS(const DMVector2d& v) const {
		return DVector2d(m[0][0] * v.x + m[0][1] * v.y, 
			m[1][0] * v.x + m[1][1] * v.y);
	}
	bool eigensystem(DMatrix2d& eigenvectors, double eigenvalues[]) const;
	/// Get i-th column
	DVector2d column(int i) const { assert(i == 0 || i == 1); return DVector2d(m[0][i], m[1][i]); }
	/// Set i-th column of this matrix from the given vector
	void setColumn(int i, const DVector2d& v) { m[0][i] = v.x; m[1][i] = v.y; }
	/// Returns i-th row
	DVector2d row(int i) const { assert(i == 0 || i == 1); return DVector2d(m[i][0], m[i][1]); }
public:
	double m[2][2];
};

class DMatrix3d
{
public:
	/// Standard constructor
	DMatrix3d(double m11 = 0.0, double m12 = 0.0, double m13 = 0.0,
			double m21 = 0.0, double m22 = 0.0, double m23 = 0.0,
			double m31 = 0.0, double m32 = 0.0, double m33 = 0.0)
	{
		m[0][0] = m11; m[0][1] = m12; m[0][2] = m13;
		m[1][0] = m21; m[1][1] = m22; m[1][2] = m23;
		m[2][0] = m31; m[2][1] = m32; m[2][2] = m33;
	}
	/// Standard constructor
	DMatrix3d(const DVector3d& v0, const DVector3d& v1, const DVector3d& v2)
	{
		m[0][0] = v0.x; m[1][0] = v0.y; m[2][0] = v0.z; 
		m[0][1] = v1.x; m[1][1] = v1.y; m[2][1] = v1.z; 
		m[0][2] = v2.x; m[1][2] = v2.y; m[2][2] = v2.z; 
	}
	double& operator()(int i, int j) { return m[i][j]; }
	const double& operator()(int i, int j) const { return m[i][j]; }
	/// Count determinant
	double det() const {
		return m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] -
			m[0][2]*m[1][1]*m[2][0] - m[0][1]*m[1][0]*m[2][2] - m[1][2]*m[2][1]*m[0][0];
	}
	/// Count solution for linear equation
	bool solve(const DVector3d& b, DVector3d& x) const;
	/// Counts the determinant of 3x2 matrix (third column equals "1")
	static double det32(double a1, double b1,
		double a2, double b2, double a3, double b3);
	/// Counts the determinant of 3x3 matrix
	static double det33(double a1, double b1, double c1,
		double a2, double b2, double c2, 
		double a3, double b3, double c3);
	/// Counts the determinant of 4x3 matrix (fourth column equals "1")
	static double det43(double a1, double b1, double c1,
		double a2, double b2, double c2, 
		double a3, double b3, double c3, 
		double a4, double b4, double c4);
	/// Counts the determinant of 4x4 matrix
	static double det44(double a1, double b1, double c1, double d1,
		double a2, double b2, double c2, double d2,
		double a3, double b3, double c3, double d3,
		double a4, double b4, double c4, double d4);
	/// Set as identity matrix
	void setIdentity(double d = 1.0) {
		m[0][0] = m[1][1] = m[2][2] = d;
		m[0][1] = m[0][2] = m[1][0] = m[1][2] = m[2][0] = m[2][1] = 0.0;
	}
	/// Set as a rotation matrix around vector v
	void setRotationAroundNormalizedVector(const DVector3d& vn, double _sin, double _cos){
		// rotation around vector (x,y,z), normalized
		assert(abs(1.0 - vn.length()) < SMALL_NUMBER);
		// | xx(1-c)+c	xy(1-c)-zs  xz(1-c)+ys	|
		// | yx(1-c)+zs	yy(1-c)+c   yz(1-c)-xs	|
		// | xz(1-c)-ys	yz(1-c)+xs  zz(1-c)+c	|
		double c = 1.0 - _cos;
		m[0][0] = vn.x * vn.x * c + _cos;
		m[0][1] = vn.x * vn.y * c - vn.z * _sin;
		m[0][2] = vn.x * vn.z * c + vn.y * _sin;
		m[1][0] = vn.y * vn.x * c + vn.z * _sin;
		m[1][1] = vn.y * vn.y * c + _cos;
		m[1][2] = vn.y * vn.z * c - vn.x * _sin;
		m[2][0] = vn.z * vn.x * c - vn.y * _sin;
		m[2][1] = vn.z * vn.y * c + vn.x * _sin;
		m[2][2] = vn.z * vn.z * c + _cos;
	}
	/// Count eigenvector in column 'ie' for eigenvalue d
	void countEigenvector(double d, int ie, DMatrix3d& eigenvectors);
	/// Count two orthonormal eigenvectors in columns i0 and i1 for the third eigenvector 'ie'
	static void countOrthonormalEigenvectors(int ie, int i0, int i1, DMatrix3d& eigenvectors, bool normalize = false);
	/// Returns inverse matrix
	DMatrix3d inverse() const {
		double f = 1.0 / det();
		return DMatrix3d(
			f*(m[1][1]*m[2][2]-m[1][2]*m[2][1]),
			f*(m[0][2]*m[2][1]-m[0][1]*m[2][2]),
			f*(m[0][1]*m[1][2]-m[0][2]*m[1][1]),
			f*(m[1][2]*m[2][0]-m[1][0]*m[2][2]),
			f*(m[0][0]*m[2][2]-m[0][2]*m[2][0]),
			f*(m[0][2]*m[1][0]-m[0][0]*m[1][2]),
			f*(m[1][0]*m[2][1]-m[1][1]*m[2][0]),
			f*(m[0][1]*m[2][0]-m[0][0]*m[2][1]),
			f*(m[0][0]*m[1][1]-m[0][1]*m[1][0]));
	}
	/// Returns transposed matrix
	DMatrix3d transposed() const {
		return DMatrix3d(m[0][0], m[1][0], m[2][0], 
			m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2]);
	}
	/// Return product of two matrices
	DMatrix3d operator*(const DMatrix3d& dm) const {
		return DMatrix3d(
			m[0][0]*dm.m[0][0]+m[0][1]*dm.m[1][0]+m[0][2]*dm.m[2][0], 
			m[0][0]*dm.m[0][1]+m[0][1]*dm.m[1][1]+m[0][2]*dm.m[2][1], 
			m[0][0]*dm.m[0][2]+m[0][1]*dm.m[1][2]+m[0][2]*dm.m[2][2], 
			m[1][0]*dm.m[0][0]+m[1][1]*dm.m[1][0]+m[1][2]*dm.m[2][0], 
			m[1][0]*dm.m[0][1]+m[1][1]*dm.m[1][1]+m[1][2]*dm.m[2][1], 
			m[1][0]*dm.m[0][2]+m[1][1]*dm.m[1][2]+m[1][2]*dm.m[2][2], 
			m[2][0]*dm.m[0][0]+m[2][1]*dm.m[1][0]+m[2][2]*dm.m[2][0], 
			m[2][0]*dm.m[0][1]+m[2][1]*dm.m[1][1]+m[2][2]*dm.m[2][1], 
			m[2][0]*dm.m[0][2]+m[2][1]*dm.m[1][2]+m[2][2]*dm.m[2][2]);
	}
	/// Return product of matrix and number
	DMatrix3d operator*(double d) const {
		return DMatrix3d(
			m[0][0] * d, m[0][1] * d, m[0][2] * d, 
			m[1][0] * d, m[1][1] * d, m[1][2] * d, 
			m[2][0] * d, m[2][1] * d, m[2][2] * d);
	}
	/// Returns product of matrix and point
	DPoint3d operator*(const DPoint3d& pt) const {
		return DPoint3d(
			m[0][0] * pt.x + m[0][1] * pt.y + m[0][2] * pt.z, 
			m[1][0] * pt.x + m[1][1] * pt.y + m[1][2] * pt.z,
			m[2][0] * pt.x + m[2][1] * pt.y + m[2][2] * pt.z);
	}
	/// Returns product of matrix and vector
	DVector3d operator*(const DVector3d& v) const {
		return DVector3d(
			m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z, 
			m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
			m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z);
	}
	/// Solves eigensystem for three-dimensional matrix (assumes all eigenvalues to be positive)
	bool eigensystem(DMatrix3d& eigenvectors, double eigenvalues[]);
	/// Get i-th column
	DVector3d column(int i) const { assert(i >= 0 && i <= 2); return DVector3d(m[0][i], m[1][i], m[2][i]); }
	/// Set i-th column of this matrix from the given vector
	void setColumn(int i, const DVector3d& v) { m[0][i] = v.x; m[1][i] = v.y; m[2][i] = v.z; }
public:
	/// Solves cubic equation (assuming x to be real and positive)
	static int countCubicRoots(double coeff[4], double x[3]);
	/// Log
	friend ostream& operator<<(ostream& os, const DMatrix3d& dm){
		return os << "{{" << dm.m[0][0] << ',' << dm.m[0][1] << ',' << dm.m[0][2] << "},{" <<
			dm.m[1][0] << ',' << dm.m[1][1] << ',' << dm.m[1][2] << "},{" <<
			dm.m[2][0] << ',' << dm.m[2][1] << ',' << dm.m[2][2] << "}}";
	}
public:
	/// matrix data
	double m[3][3];
};

#endif // !defined(DMATRIX_H__INCLUDED)
