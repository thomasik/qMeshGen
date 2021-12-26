// DMatrix.cpp: implementation of the DMetric2d class.
//
//////////////////////////////////////////////////////////////////////

#include "DMatrix.h"

/// Count solution for linear equation
bool DMatrix2d::solve(const DVector2d& b, DVector2d& vx) const
{
	double w = det();
	double wx = b.x*m[1][1] - m[0][1]*b.y;
	double wy = m[0][0]*b.y	- b.x*m[1][0];
	if(w == 0.0) return false; // more epsilon-like test?
	vx.x = wx/w;
	vx.y = wy/w;
	return true;
}

bool DMatrix2d::eigensystem(DMatrix2d& eigenvectors, double eigenvalues[]) const
{
	if(abs(m[0][1])+abs(m[1][0]) < SMALL_NUMBER){
		eigenvectors.m[0][0] = eigenvectors.m[1][1] = 1.0;
		eigenvectors.m[0][1] = eigenvectors.m[1][0] = 0.0;
		eigenvalues[0] = m[0][0];
		eigenvalues[1] = m[1][1];
		return true;
	}
	double delta = sqr(m[0][0] - m[1][1]) + 4.0*m[0][1]*m[1][0];
	if(delta < 0.0) return false;
	delta = sqrt(delta);
	eigenvalues[0] = 0.5 * (m[0][0] + m[1][1] + delta);
	eigenvalues[1] = 0.5 * (m[0][0] + m[1][1] - delta);

	for(int i = 0; i < 2; i++){
		double a11 = m[0][0]-eigenvalues[i];
		double a22 = m[1][1]-eigenvalues[i];
		// select "better" equation
		if(abs(a11)+abs(m[0][1]) > abs(m[1][0])+abs(a22)){
			eigenvectors.m[1][i] = 1.0;
			eigenvectors.m[0][i] = -m[0][1]/a11;
		}else{
			eigenvectors.m[1][i] = 1.0;
			eigenvectors.m[0][i] = -a22/m[1][0];
		}
	}

	return true;
}

/*******************************************************************************
 * FindCubicRoots (http://www.worldserver.com/turk/opensource/FindCubicRoots.c.txt)
 *
 *	Solve:
 *		coeff[3] * x^3 + coeff[2] * x^2 + coeff[1] * x + coeff[0] = 0
 *
 *	returns:
 *		3 - 3 real roots
 *		1 - 1 real root (2 complex conjugate)
 *******************************************************************************/
int DMatrix3d::countCubicRoots(double coeff[4], double x[3])
{
	double a1 = coeff[2] / coeff[3];
	double a2 = coeff[1] / coeff[3];
	double a3 = coeff[0] / coeff[3];

	double Q = (a1 * a1 - 3 * a2) / 9.0;
	double R = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
	double Qcubed = Q * Q * Q;
	double d = Qcubed - R * R;

	/* Three real roots */
	if (d > -0.1) {
//	if (d >= 0.0) {
		if(Q > 0.0){
			double temp = R / sqrt(Qcubed);
			if(temp > 1.0) temp = 1.0;
			else if(temp < -1.0) temp = -1.0;
			double theta = acos(temp);
			double sqrtQ = sqrt(Q);
			x[0] = -2.0 * sqrtQ * cos( theta             / 3.0) - a1 / 3.0;
			x[1] = -2.0 * sqrtQ * cos((theta + 2.0 * PI) / 3.0) - a1 / 3.0;
			x[2] = -2.0 * sqrtQ * cos((theta + 4.0 * PI) / 3.0) - a1 / 3.0;
		}else{
			x[0] = x[1] = x[2] = - a1 / 3.0;
		}
		return 3;
	}else{
	/* One real root */
		double e = pow(sqrt(-d) + abs(R), 1.0 / 3.0);
		if (R > 0.0)
			e = -e;
		x[0] = (e + Q / e) - a1 / 3.0;
		return 1;
	}
}

void DMatrix3d::countEigenvector(double d, int ie, DMatrix3d& eigenvectors)
{
	for(int i = 0; i < 3; i++) m[i][i] -= d;

	double w[3];

	for(int i = 0; i < 3; i++){
		int i0 = (i + 1) % 3;
		int i1 = (i + 2) % 3;
		w[i] = m[i0][i0]*m[i1][i1] - m[i0][i1]*m[i1][i0];
	}

	int min_i = (abs(w[0]) > abs(w[1])) ? 0 : 1;
	if(abs(w[2]) > abs(w[min_i])) min_i = 2;
//	assert(abs(w[min_i]) > SMALL_NUMBER);

	int i0 = (min_i + 1) % 3;
	int i1 = (min_i + 2) % 3;

	double w0 = m[i1][min_i]*m[i0][i1] - m[i0][min_i]*m[i1][i1];
	double w1 = m[i0][min_i]*m[i1][i0] - m[i0][i0]*m[i1][min_i];

	eigenvectors.m[min_i][ie] = 1.0;
	eigenvectors.m[i0][ie] = w0/w[min_i];
	eigenvectors.m[i1][ie] = w1/w[min_i];

	for(int i = 0; i < 3; i++) m[i][i] += d;
}

void DMatrix3d::countOrthonormalEigenvectors(int ie, int i0, int i1, DMatrix3d& e, bool normalize)
{
	int min_i = (abs(e.m[0][ie]) < abs(e.m[1][ie])) ? 0 : 1;
	if(abs(e.m[2][ie]) < abs(e.m[min_i][ie])) min_i = 2;

	e.m[min_i][i0] = 0.0;
	int j0 = (min_i+1)%3;
	int j1 = (min_i+2)%3;
	e.m[j0][i0] =  e.m[j1][ie];
	e.m[j1][i0] = -e.m[j0][ie];

	e.m[0][i1] = e.m[1][ie]*e.m[2][i0] - e.m[1][i0]*e.m[2][ie];
	e.m[1][i1] = e.m[0][i0]*e.m[2][ie] - e.m[0][ie]*e.m[2][i0];
	e.m[2][i1] = e.m[0][ie]*e.m[1][i0] - e.m[0][i0]*e.m[1][ie];

	if(normalize){
		double f = 1.0 / sqrt(sqr(e.m[0][i0])+sqr(e.m[1][i0])+sqr(e.m[2][i0]));
		for(int i = 0; i < 3; i++) e.m[i][i0] *= f;
		f = 1.0 / sqrt(sqr(e.m[0][i1])+sqr(e.m[1][i1])+sqr(e.m[2][i1]));
		for(int i = 0; i < 3; i++) e.m[i][i1] *= f;
	}
}

bool DMatrix3d::eigensystem(DMatrix3d& eigenvectors, double d[])
{
	double sf = 0.0;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			if(abs(m[i][j]) > sf) sf = abs(m[i][j]);

	double sfr = 1.0 / sf;
	DMatrix3d smat = *this * sfr;

	double coeff[4] = {
		-smat.m[0][2]*smat.m[1][1]*smat.m[2][0] + smat.m[0][1]*smat.m[1][2]*smat.m[2][0] 
			+ smat.m[0][2]*smat.m[1][0]*smat.m[2][1] - smat.m[0][0]*smat.m[1][2]*smat.m[2][1] 
			- smat.m[0][1]*smat.m[1][0]*smat.m[2][2] + smat.m[0][0]*smat.m[1][1]*smat.m[2][2],
		smat.m[0][1]*smat.m[1][0] - smat.m[0][0]*smat.m[1][1] 
			+ smat.m[0][2]*smat.m[2][0] + smat.m[1][2]*smat.m[2][1] 
			- smat.m[0][0]*smat.m[2][2] - smat.m[1][1]*smat.m[2][2],
		smat.m[0][0] + smat.m[1][1] + smat.m[2][2],
		-1.0};

	int dct = countCubicRoots(coeff, d);
	assert(dct == 3);
	if(dct != 3) return false;
	d[0] *= sf;
	d[1] *= sf;
	d[2] *= sf;

	const double SAME_EPS = 1e-3;
	bool same_01 = (abs(1.0-d[0]/d[1]) < SAME_EPS);
	bool same_02 = (abs(1.0-d[0]/d[2]) < SAME_EPS);
	bool same_12 = (abs(1.0-d[1]/d[2]) < SAME_EPS);

	if(same_01 && same_02){
		// all eignevalues equal
		eigenvectors.setIdentity();
		return true;
	}else if(same_01){
		countEigenvector(d[2], 2, eigenvectors);
		// set 0,1 as orthonormal
		countOrthonormalEigenvectors(2, 0, 1, eigenvectors, false);
	}else if(same_02){
		countEigenvector(d[1], 1, eigenvectors);
		// set 0,2 as orthonormal
		countOrthonormalEigenvectors(1, 0, 2, eigenvectors, false);
	}else if(same_12){
		countEigenvector(d[0], 0, eigenvectors);
		// set 1,2 as orthonormal
		countOrthonormalEigenvectors(0, 1, 2, eigenvectors, false);
	}else{
		// count all eigenvectors[0..2][i]
		for(int i = 0; i < 3; i++){
			countEigenvector(d[i], i, eigenvectors);
		}
	}

	return true;
}

/// Counts the determinant of 3x2 matrix (third column equals "1")
double DMatrix3d::det32(double a1, double b1, double a2, double b2, double a3, double b3)
{
	return (a1-a3)*(b2-b3) - (b1-b3)*(a2-a3);
}

/// Counts the determinant of 3x3 matrix
double DMatrix3d::det33(double a1, double b1, double c1,
		double a2, double b2, double c2, double a3, double b3, double c3)
{
	return a1*b2*c3+a2*b3*c1+a3*b1*c2-a1*b3*c2-a2*b1*c3-a3*b2*c1;
}

/// Counts the determinant of 4x3 matrix
double DMatrix3d::det43(double a1, double b1, double c1,	double a2, double b2, double c2, 
		double a3, double b3, double c3, double a4, double b4, double c4)
{
	return (a1-a4)*((b2-b4)*(c3-c4)-(c2-c4)*(b3-b4)) - 
		(b1-b4)*((a2-a4)*(c3-c4)-(c2-c4)*(a3-a4)) + 
		(c1-c4)*((a2-a4)*(b3-b4)-(b2-b4)*(a3-a4));
}

/// Counts the determinant of 4x4 matrix
double DMatrix3d::det44(double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2,
		double a3, double b3, double c3, double d3, double a4, double b4, double c4, double d4)
{
	return a1*det33(b2,c2,d2,b3,c3,d3,b4,c4,d4) - b1*det33(a2,c2,d2,a3,c3,d3,a4,c4,d4)
		+ c1*det33(a2,b2,d2,a3,b3,d3,a4,b4,d4) - d1*det33(a2,b2,c2,a3,b3,c3,a4,b4,c4);
}

/// Count solution for linear equation
bool DMatrix3d::solve(const DVector3d& b, DVector3d& vx) const
{
	double w = det();
	double wx = b.x*m[1][1]*m[2][2] + m[0][1]*m[1][2]*b.z + m[0][2]*b.y*m[2][1] -
			m[0][2]*m[1][1]*b.z - m[0][1]*b.y*m[2][2] - m[1][2]*m[2][1]*b.x;
	double wy = m[0][0]*b.y*m[2][2] + b.x*m[1][2]*m[2][0] + m[0][2]*m[1][0]*b.z -
			m[0][2]*b.y*m[2][0] - b.x*m[1][0]*m[2][2] - m[1][2]*b.z*m[0][0];
	double wz = m[0][0]*m[1][1]*b.z + m[0][1]*b.y*m[2][0] + b.x*m[1][0]*m[2][1] -
			b.x*m[1][1]*m[2][0] - m[0][1]*m[1][0]*b.z - b.y*m[2][1]*m[0][0];
	if(w == 0.0) return false; // more epsilon-like test?
	vx.x = wx/w;
	vx.y = wy/w;
	vx.z = wz/w;
	return true;
}
