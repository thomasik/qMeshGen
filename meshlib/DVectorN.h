/////////////////////////////////////////////////////////////////////////////
// DVectorN.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DVECTORN_H__INCLUDED
#define DVECTORN_H__INCLUDED

#include <iostream>
#include <cmath>

/// template for one-dimensional vector of double for numerical applications
template<int N>
class DVectorN
{
public:
	DVectorN() {} // default, no initialization
	DVectorN(double val) { std::fill_n(m, N, val); }
	double& operator[](int i) { return m[i]; }
	const double& operator[](int i) const { return m[i]; }
public:
	DVectorN operator*( double d ) {
		DVectorN result;
		for(int i = 0; i < N; i++) result.m[i] = m[i] * d;
		return result;
	}
	DVectorN& operator+=( const DVectorN & v ) {
		for(int i = 0; i < N; i++) 
			m[i] += v.m[i];
		return *this;
	}
	DVectorN& operator*=( double d ) {
		for(int i = 0; i < N; i++) 
			m[i] *= d;
		return *this;
	}
public:
	/// Returns the length of the vector
	double length() const { return sqrt(length2()); }
	/// Returns the squared length of the vector
	double length2() const { 
		double sum = 0.0;
		for(int i = 0; i < N; i++) sum += m[i] * m[i];
		return sum;
	}
	/// Normalizes the vector (makes its length equal to 1)
	DVectorN& normalize() {
		double len_ratio = 1.0 / length();
		for(int i = 0; i < N; i++) m[i] *= len_ratio;
		return *this;
	}
	/// Returnes the normalized vector (length equal to 1)
	const DVectorN normalized() const{
		double len_ratio = 1.0 / length();
		DVectorN result;
		for(int i = 0; i < N; i++) result.m[i] = m[i] * len_ratio;
		return result;
	}
	friend std::ostream& operator<<(std::ostream& os, const DVectorN<N> & vn)
	{
		os << vn[0];
		for(int i = 1; i < N; i++) os << " " << vn[i];
		return os;
	}
	friend std::istream& operator>>(std::istream& is, DVectorN<N> & vn)
	{
		for(int i = 0; i < N; i++) is >> std::ws >> vn[i];
		return is;
	}

private:
	double m[N];
};

class DVectorNV
{
public:
	DVectorNV(int _n) { // default, no initialization
		m = new double[n = _n];
	} 
	DVectorNV(int _n, double val) {
		m = new double[n = _n];
		for(int i = 0; i < n; i++) m[i] = val;
	}
	~DVectorNV() {
		delete[] m;
	}
	double& operator[](int i) { return m[i]; }
	const double& operator[](int i) const { return m[i]; }
public:
	/// Returns the length of the vector
	double length() const { return std::sqrt(length2()); }
	/// Returns the squared length of the vector
	double length2() const { 
		double sum = 0.0;
		for(int i = 0; i < n; i++) sum += m[i] * m[i];
		return sum;
	}
	/// Normalizes the vector (makes its length equal to 1)
	DVectorNV& normalize() {
		double len_ratio = 1.0 / length();
		for(int i = 0; i < n; i++) m[i] *= len_ratio;
		return *this;
	}
	/// Returnes the normalized vector (length equal to 1)
	const DVectorNV normalized() const{
		double len_ratio = 1.0 / length();
		DVectorNV result(n);
		for(int i = 0; i < n; i++) result.m[i] = m[i] * len_ratio;
		return result;
	}
private:
	int n;
	double* m;
};

#endif // !defined(DVECTORN_H__INCLUDED)
