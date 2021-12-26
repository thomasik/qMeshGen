/////////////////////////////////////////////////////////////////////////////
// DMatrix.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(DMATRIXN_H__INCLUDED)
#define DMATRIXN_H__INCLUDED

#include "DVectorN.h"
#include <string.h>

/// template for two-dimensional matrix of doubles for numerical applications
template<int N>
class DMatrixN
{
public:
	DMatrixN() {}
	DMatrixN(double val) {
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
				m[i][j] = val;
	}
	double& operator()(int i, int j) { return m[i][j]; }
	const double& operator()(int i, int j) const { return m[i][j]; }
public:
	// solve linear system using LU, from NumericalRecipes#3
	bool solveLU(const DVectorN<N>& b, DVectorN<N>& x) {
		int i, ii=0, imax, ip, j, k;
		double sum, big, temp;
		//double d;
		double vv[N];
		int indx[N];
		double lu[N][N];;
		double aref[N][N];

		memcpy(lu, m, N*N*sizeof(double));
		memcpy(aref, m, N*N*sizeof(double));
		
		// LU decompose
		for (i=0; i<N; i++) { // Loop over rows to get the implicit scaling information
			big=0.0; 
			for (j=0; j<N; j++)
				if ((temp=abs(lu[i][j])) > big) big=temp;
			if (big == 0.0) return false; // throw("Singular matrix in LUdcmp");
										  // No nonzero largest element.
			vv[i]=1.0/big; //Save the scaling.
		}
		for (k=0; k<N; k++) { // This is the outermost kij loop.
			big=0.0; // Initialize for the search for largest PIvot element.
			for (i=k; i<N; i++) {
				temp=vv[i]*abs(lu[i][k]);
				if (temp > big) { // Is the figure of merit for the PIvot better than
					big=temp;		// the best so far?
					imax=i;
				}
			}
			if (k != imax) { // Do we need to interchange rows?
				for (j=0;j<N;j++) { // Yes, do so...
					temp=lu[imax][j];
					lu[imax][j]=lu[k][j];
					lu[k][j]=temp;
				}
				//d = -d; // ...and change the parity of d.
				vv[imax]=vv[k]; // Also interchange the scale factor.
			}
			indx[k]=imax;
			if (lu[k][k] == 0.0) lu[k][k]=VERY_SMALL_NUMBER;
			for (i=k+1;i<N;i++) {
				temp=lu[i][k] /= lu[k][k]; // Divide by the PIvot element.
				for (j=k+1;j<N;j++) // Innermost loop: reduce remaining submatrix.
					lu[i][j] -= temp*lu[k][j];
			}
		}

		// solve
		for (i=0; i<N; i++) x[i] = b[i];
		for (i=0; i<N; i++) { 
			ip=indx[i];
			sum=x[ip];
			x[ip]=x[i];
			if (ii != 0)
				for (j=ii-1; j<i; j++) sum -= lu[i][j]*x[j];
			else if (sum != 0.0) // A nonzero element was encountered, so from now on we
				ii=i+1;			// will have to do the sums in the loop above.
			x[i]=sum;
		}
		for (i=N-1; i>=0; i--) { // Now we do the backsubstitution
			sum=x[i];
			for (j=i+1; j<N; j++) sum -= lu[i][j]*x[j];
			x[i]=sum/lu[i][i]; //Store a component of the solution vector X.
		} // All done!

		return true;
	}

	double smallestEigenvectorNormalized(DVectorN<N>& ve){
		// Jacobi method
		DMatrixN<N> v(0.0); // matrix with columns -> eigenvectors
		DVectorN<N> d; // eigenvalues
		int nrot = 0;
		const double EPS = SMALL_NUMBER;
		// ---
		DVectorN<N> b;
		DVectorN<N> z(0.0); // This vector will accumulate terms of the form tapq.
		for(int i=0; i < N; i++){
			// Initialize v to the identity matrix.
			v.m[i][i] = 1.0;
			// ... Initialize b and d to the diagonal of a.
			b[i] = d[i] = m[i][i];
		}
		for(int i = 1;i <= 50; i++){
			double sm = 0.0;
			for (int ip = 0; ip < N-1; ip++){ // Sum magnitude of off-diagonal elements
				for(int iq = ip+1; iq < N; iq++)
					sm += abs(m[ip][iq]);
			}
			if(sm == 0.0){ 
				// The normal return, which relies on quadratic convergence to
				// machine underflow.
				// -> eigsrt(d,&v);
				// ... find minimum eigenvalue
				int imin = 0;
				for(int j = 1; j < N; j++)
					if(d[j] < d[imin]) imin = j; // with abs ?
				// ... return proper eigenvector
				for(int j = 0; j < N; j++)
					ve[j] = v.m[j][imin];
				ve.normalize();
				return d[imin];
			}
			double tresh = (i < 4) ? 0.2*sm/(N*N) : 0.0;
			for(int ip = 0; ip < N-1; ip++){
				for(int iq = ip+1; iq < N; iq++){
					double g = 100.0 * abs(m[ip][iq]);
					//After four sweeps, skip the rotation if the off-diagonal element is small.
					if(i > 4 && g <= EPS*abs(d[ip]) && g <= EPS*abs(d[iq]))
						m[ip][iq] = 0.0;
					else if(abs(m[ip][iq]) > tresh){
						double h = d[iq] - d[ip];
						double t;
						if(g <= EPS*abs(h))
							t = (m[ip][iq]) / h;
						else{
							double theta = 0.5*h/(m[ip][iq]);
							t=1.0/(abs(theta)+sqrt(1.0+theta*theta));
							if(theta < 0.0) t = -t;
						}
						double c = 1.0/sqrt(1+t*t);
						double s = t*c;
						double tau = s/(1.0+c);
						h= t * m[ip][iq];
						z[ip] -= h;
						z[iq] += h;
						d[ip] -= h;
						d[iq] += h;
						m[ip][iq]=0.0;
						for(int j = 0; j < ip; j++)
							rot(*this, s, tau, j, ip, j, iq);
						for(int j = ip+1; j < iq; j++)
							rot(*this, s, tau, ip, j, j, iq);
						for(int j = iq+1; j < N; j++)
							rot(*this, s, tau, ip, j, iq, j);
						for(int j = 0; j < N; j++)
							rot(v, s, tau, j, ip, j, iq);
						++nrot;
					}
				}
			}
			for(int ip = 0; ip < N; ip++){
				b[ip] += z[ip];
				d[ip] = b[ip]; //Update d with the sum of tapq,
				z[ip] = 0.0; //and reinitialize z.
			}
		}
		// throw("Too many iterations in routine jacobi");
		return -1.0;
	}
private:
	static void rot(DMatrixN<N> &a, double s, double tau, int i, int j, int k, int l){
		double g = a.m[i][j];
		double h = a.m[k][l];
		a.m[i][j] = g-s*(h+g*tau);
		a.m[k][l] = h+s*(g-h*tau);
	}
private:
	double m[N][N];
};

// no-template version
class DMatrixNV
{
public:
	DMatrixNV(int _n) : n(_n) {
		m = new double*[n];
		for(int i = 0; i < n; i++)
			m[i] = new double[n];
	}
	DMatrixNV(int _n, double val) : n(_n) {
		m = new double*[n];
		for(int i = 0; i < n; i++)
			m[i] = new double[n];

		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				m[i][j] = val;
	}
	~DMatrixNV() {
		for(int i = 0; i < n; i++)
			delete[] m[i];
		delete[] m;
	}
	double& operator()(int i, int j) { return m[i][j]; }
	const double& operator()(int i, int j) const { return m[i][j]; }
public:
	double smallestEigenvectorNormalized(DVectorNV& ve){
		// Jacobi method
		DMatrixNV v(n, 0.0); // matrix with columns -> eigenvectors
		DVectorNV d(n); // eigenvalues
		int nrot = 0;
		const double EPS = SMALL_NUMBER;
		// ---
		DVectorNV b(n);
		DVectorNV z(n, 0.0); // This vector will accumulate terms of the form tapq.
		for(int i=0; i < n; i++){
			// Initialize v to the identity matrix.
			v.m[i][i] = 1.0;
			// ... Initialize b and d to the diagonal of a.
			b[i] = d[i] = m[i][i];
		}
		for(int i = 1;i <= 50; i++){
			double sm = 0.0;
			for (int ip = 0; ip < n-1; ip++){ // Sum magnitude of off-diagonal elements
				for(int iq = ip+1; iq < n; iq++)
					sm += abs(m[ip][iq]);
			}
			if(sm == 0.0){ 
				// The normal return, which relies on quadratic convergence to
				// machine underflow.
				// -> eigsrt(d,&v);
				// ... find minimum eigenvalue
				int imin = 0;
				for(int j = 1; j < n; j++)
					if(d[j] < d[imin]) imin = j; // with abs ?
				// ... return proper eigenvector
				for(int j = 0; j < n; j++)
					ve[j] = v.m[j][imin];
				ve.normalize();
				return d[imin];
			}
			double tresh = (i < 4) ? 0.2*sm/(n*n) : 0.0;
			for(int ip = 0; ip < n-1; ip++){
				for(int iq = ip+1; iq < n; iq++){
					double g = 100.0 * abs(m[ip][iq]);
					//After four sweeps, skip the rotation if the off-diagonal element is small.
					if(i > 4 && g <= EPS*abs(d[ip]) && g <= EPS*abs(d[iq]))
						m[ip][iq] = 0.0;
					else if(abs(m[ip][iq]) > tresh){
						double h = d[iq] - d[ip];
						double t;
						if(g <= EPS*abs(h))
							t = (m[ip][iq]) / h;
						else{
							double theta = 0.5*h/(m[ip][iq]);
							t=1.0/(abs(theta)+sqrt(1.0+theta*theta));
							if(theta < 0.0) t = -t;
						}
						double c = 1.0/sqrt(1+t*t);
						double s = t*c;
						double tau = s/(1.0+c);
						h= t * m[ip][iq];
						z[ip] -= h;
						z[iq] += h;
						d[ip] -= h;
						d[iq] += h;
						m[ip][iq]=0.0;
						for(int j = 0; j < ip; j++)
							rot(*this, s, tau, j, ip, j, iq);
						for(int j = ip+1; j < iq; j++)
							rot(*this, s, tau, ip, j, j, iq);
						for(int j = iq+1; j < n; j++)
							rot(*this, s, tau, ip, j, iq, j);
						for(int j = 0; j < n; j++)
							rot(v, s, tau, j, ip, j, iq);
						++nrot;
					}
				}
			}
			for(int ip = 0; ip < n; ip++){
				b[ip] += z[ip];
				d[ip] = b[ip]; //Update d with the sum of tapq,
				z[ip] = 0.0; //and reinitialize z.
			}
		}
		// throw("Too many iterations in routine jacobi");
		return -1.0;
	}
private:
	static void rot(DMatrixNV &a, double s, double tau, int i, int j, int k, int l){
		double g = a.m[i][j];
		double h = a.m[k][l];
		a.m[i][j] = g-s*(h+g*tau);
		a.m[k][l] = h+s*(g-h*tau);
	}
private:
	int n;
	double** m;
};

#endif // !defined(DMATRIXN_H__INCLUDED)
