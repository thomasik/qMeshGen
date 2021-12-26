// DHesjan.cpp: implementation of the DHesjan class.
//
//////////////////////////////////////////////////////////////////////

#include "DHesjan.h"
#include "common.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double	DHesjan::m_factor = 0.25;
bool	DHesjan::m_trim_min = false;
bool	DHesjan::m_trim_max = false;
double	DHesjan::m_min_len  = 1.5;
double	DHesjan::m_max_len  = 1.5;
bool	DHesjan::m_with_sqrt= false;

double DHesjan::bx[9]  = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double DHesjan::by[9]  = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double DHesjan::bxx[9] = {0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double DHesjan::byy[9] = {0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double DHesjan::bxy[9] = {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

double DHesjan::param_hesjan_factor    = 0.01;
double DHesjan::param_hesjan_min_ratio = 0.01;
double DHesjan::param_hesjan_max_ratio = 0.2;

DHesjan::DHesjan()
{
	valid = false;
	point_ct = 0;
}

DHesjan::~DHesjan()
{

}

bool DHesjan::solveBG(int n, double A[][9], double *f, double f0, int method, bool count_eigenvalues)
{
	double c[9];
	double d[9];
	int i, j;

	assert( n <= 9 );

	point_ct = n;

	if(method == 0){	// Trzy uk³ady równañ dla wspó³czynników + metoda SQR
		singular = qrdcmp(n, A, c, d);

		for(i = 0; i < n; i++){
			ax[i] = bx[i];
			ay[i] = by[i];
			axx[i] = bxx[i];
			axy[i] = bxy[i];
			ayy[i] = byy[i];
		}

		qrsolv(n, A, ax, c, d);
		qrsolv(n, A, ay, c, d);
		qrsolv(n, A, axx, c, d);
		qrsolv(n, A, axy, c, d);
		qrsolv(n, A, ayy, c, d);

		Dx = Dy = 0.0;
		Dxx = Dyy = Dxy = 0.0;
		for(i = 0; i < n; i++){
			Dx += ax[i] * (f[i] - f0);
			Dy += ay[i] * (f[i] - f0);
			Dxx += axx[i] * (f[i] - f0);
			Dyy += ayy[i] * (f[i] - f0);
			Dxy += axy[i] * (f[i] - f0);
		}
	}else if(method == 1){	// Trzy uk³ady równañ dla wspó³czynników + metoda LES
		double AA[9][9];

		singular = false;	// nie uzywana

		for(i = 0; i < n; i++){
			c[i] = bx[i];
			for(j = 0; j < n; j++) AA[i][j] = A[i][j];
		}
		LES(n, n, AA, c, ax);
		for(i = 0; i < n; i++){
			c[i] = by[i];
			for(j = 0; j < n; j++) AA[i][j] = A[i][j];
		}
		LES(n, n, AA, c, ay);
		for(i = 0; i < n; i++){
			c[i] = bxx[i];
			for(j = 0; j < n; j++) AA[i][j] = A[i][j];
		}
		LES(n, n, AA, c, axx);
		for(i = 0; i < n; i++){
			c[i] = bxy[i];
			for(j = 0; j < n; j++) AA[i][j] = A[i][j];
		}
		LES(n, n, AA, c, axy);
		for(i = 0; i < n; i++){
			c[i] = byy[i];
			for(j = 0; j < n; j++) AA[i][j] = A[i][j];
		}
		LES(n, n, AA, c, ayy);

		Dx = Dy = 0.0;
		Dxx = Dyy = Dxy = 0.0;
		for(i = 0; i < n; i++){
			Dx += ax[i] * (f[i] - f0);
			Dy += ay[i] * (f[i] - f0);
			Dxx += axx[i] * (f[i] - f0);
			Dyy += ayy[i] * (f[i] - f0);
			Dxy += axy[i] * (f[i] - f0);
		}
	}else if(method == 2){	// Orkisz + LES
		double AA[9][9];

		singular = false;	// nie uzywana

		for(i = 0; i < n; i++)
			for(j = 0; j < n; j++) 
				AA[i][j] = A[j][i];		
		for(i = 0; i < n; i++){
			AA[i][2] *= 0.5;
			AA[i][3] *= 0.5;
			c[i] = f[i] - f0;
		}
		LES(n, 5, AA, c, d);

		Dx = d[0];
		Dy = d[1];
		Dxx = d[2];
		Dyy = d[3];
		Dxy = d[4];
	}else{
		assert(false);
	}

	if(count_eigenvalues)
		countEigenvalues();

	return valid = true;
}

#define SIGN(a,b) ((b) > 0 ? abs(a) : -abs(a))

bool DHesjan::qrdcmp(int n, double A[][9], double *c, double *d)
{
	int i, j, k;
	double scale, sigma, sum, tau;
	bool sing = false;

	for(k = 0; k < n-1; k++){
		scale=0.0;
		for(i = k; i < n; i++) scale=std::max(scale, abs(A[i][k]));
		if(scale == 0.0){	// Singular case.
			sing = true;
			c[k] = d[k] = 0.0;
		}else{	// Form Q k and Q k . A.
			for(i = k; i < n; i++) A[i][k] /= scale;
			for(sum = 0.0, i = k; i < n; i++) sum += A[i][k] * A[i][k];
			sigma=SIGN(sqrt(sum),A[k][k]);
			A[k][k] += sigma;
			c[k] = sigma * A[k][k];
			d[k] = -scale * sigma;
			for(j = k+1; j < n; j++){
				for(sum = 0.0, i = k; i < n; i++) sum += A[i][k] * A[i][j];
				tau = sum/c[k];
				for(i = k; i < n; i++) A[i][j] -= tau * A[i][k];
			}
		}
	}
	d[n-1]=A[n-1][n-1];
	if(d[n-1] == 0.0) sing = true;

	return sing;
}

void DHesjan::qrsolv(int n, double A[][9], double *b, double *c, double *d)
{
	int i, j;
	double sum, tau;

	for(j = 0; j < n-1; j++){ // Form Q T. b.
		for(sum = 0.0, i = j; i < n; i++) sum += A[i][j] * b[i];
		tau = sum / c[j];
		for(i = j; i < n;i++) b[i] -= tau * A[i][j];
	}

	b[n-1] /= d[n-1];
	for(i = n-2; i >= 0; i--){
		for(sum = 0.0, j = i+1; j < n; j++) sum += A[i][j] * b[j];
		b[i] = (b[i] - sum) / d[i];
	}
}

void DHesjan::countEigenvalues()
{
	// Wylicz dwie wartoœci w³asne
	double delta = sqr(Dxx - Dyy) + 4.0*Dxy*Dxy;
	if(delta < 0.0){
		valid = false;
		return;
	}
	delta = sqrt(delta);
	w1 = 0.5 * (Dxx + Dyy + delta);
	w2 = 0.5 * (Dxx + Dyy - delta);
	// Wylicz k¹t
	if(abs(Dxy) < SMALL_NUMBER){
		stretch.angle = 0.0;
	}else{
		double x2 = -(Dxx - w1) / Dxy;
		stretch.angle = atan(x2);
	}
	calculateLength();
}

void DHesjan::calculateLength()
{
	// Przekszta³æ wartoœci w³asne na wyd³u¿enia lx, ly
	if(abs(w1) < SMALL_NUMBER){
		stretch.lx = m_max_len;
	}else{
		stretch.lx = m_factor / (m_with_sqrt ? sqrt(abs(w1)) : abs(w1));
		if(m_trim_min && stretch.lx < m_min_len) stretch.lx = m_min_len;
		if(m_trim_max && stretch.lx > m_max_len) stretch.lx = m_max_len;
	}
	if(abs(w2) < SMALL_NUMBER){
		stretch.ly = m_max_len;
	}else{
		stretch.ly = m_factor / (m_with_sqrt ? sqrt(abs(w2)) : abs(w2));
		if(m_trim_min && stretch.ly < m_min_len) stretch.ly = m_min_len;
		if(m_trim_max && stretch.ly > m_max_len) stretch.ly = m_max_len;
	}
}

void DHesjan::convertToLength(double &len_x, double &len_y)
{
	// Przekszta³æ wartoœci w³asne na wyd³u¿enia lx, ly
	if(abs(len_x) < SMALL_NUMBER){
		len_x = LARGE_NUMBER;
	}else{
		len_x = m_factor / (m_with_sqrt ? sqrt(abs(len_x)) : abs(len_x));
		if(m_trim_min && len_x < m_min_len) len_x = m_min_len;
		if(m_trim_max && len_x > m_max_len) len_x = m_max_len;
	}
	if(abs(len_y) < SMALL_NUMBER){
		len_y = LARGE_NUMBER;
	}else{
		len_y = m_factor / (m_with_sqrt ? sqrt(abs(len_y)) : abs(len_y));
		if(m_trim_min && len_y < m_min_len) len_y = m_min_len;
		if(m_trim_max && len_y > m_max_len) len_y = m_max_len;
	}
}

void DHesjan::LES(int N, int M, double A[][9], double *B, double *X)
{
	int i, j;
	
	double R[9][9], D[9], BB[9];

	QRN(N, M, A, R, D);
	
	for(j = 0; j < M; j++) {
		BB[j] = 0.0;
		for(i = 0; i < N; i++) {
			A[i][j] /= D[j];
			BB[j] += A[i][j] * B[i];
		}
	}

	for(i = 0; i < M; i++)
		B[i] = BB[i];
	
	X[M-1] = B[M-1];
	for(i = M-2; i >= 0; i--) {
		X[i] = B[i];
		for(j = M-1; j > i; j--)
			X[i] -= R[i][j] * X[j];
	}
}

void DHesjan::QRN(int N, int M, double A[][9], double R[][9], double *D)
{
	int i, j, k;
	
	D[0] = 0.0;
	for(i = 0; i < N; i++)
		D[0] += A[i][0] * A[i][0];
	R[0][0] = 1.0;
	
	for(j = 1; j < M; j++) {					// for all column q_j
		for(i = 0; i < j; i++) {				// for all column q_i : i=0..j-1
			R[j][i] = R[i][j] = 0.0;			// calculating <q_j,a_i>/||q_j||
			for(k = 0; k < N; k++)
				R[i][j] += A[k][i] * A[k][j];
			R[i][j] /= D[i];

			for(k = 0; k < N; k++)				// updating q_j
				A[k][j] -= R[i][j] * A[k][i];
		}

		D[j] = 0.0;
		for(i = 0; i < N; i++)
			D[j] += A[i][j] * A[i][j];
		R[j][j] = 1.0;			// calculating norm of the first column of matrix A
	}
}
