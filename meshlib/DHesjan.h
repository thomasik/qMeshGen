// DHesjan.h: interface for the DHesjan class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DHESJAN_H__INCLUDED)
#define DHESJAN_H__INCLUDED

#include "DMetric2d.h"

class DHesjan  
{
public:
	DHesjan();
	virtual ~DHesjan();
public:
	void QRN(int N, int M, double A[][9], double R[][9], double* D);
	void LES(int N, int M, double A[][9], double* B, double* X);
	void countEigenvalues();
	void calculateLength();
	void qrsolv(int n, double A[][9], double *b, double *c, double *d);
	bool qrdcmp(int n, double A[][9], double *c, double *d);
	bool solveBG(int n, double A[][9], double *f, double f0, int method = 0, bool count_eigenvalues = true);
//	DTriple getTriple(){ return DTriple(angle, lx, ly); }
public:
	double Dxx, Dyy, Dxy, Dx, Dy;
	double ax[9];
	double ay[9];
	double axx[9];
	double axy[9];
	double ayy[9];
	double w1, w2;
	ControlDataStretch2d stretch;
	int point_ct;
	int	point_nrs[9];
	int point_type[9];
	bool valid, singular;
public:
	static void convertToLength(double& len_x, double& len_y);
	static double m_factor;
	static bool m_trim_min;
	static bool m_trim_max;
	static double m_min_len;
	static double m_max_len;
	static bool m_with_sqrt;
	static double bxx[9], bxy[9], byy[9];	// Wektory do uk³adu równañ dla hesjanu
	static double bx[9], by[9];				// Wektory do uk³adu równañ dla hesjanu
public:
	static double param_hesjan_factor;
	static double param_hesjan_min_ratio;
	static double param_hesjan_max_ratio;
};

#endif // !defined(DHESJAN_H__INCLUDED)
