// DMetric3d.cpp: implementation of the DMetric2d class.
//
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include <vector>

#include "DMetric3d.h"
#include "MeshData.h"
#include "EPSFile.h"
#include "DEquation.h"
#include "ControlSpace2dAdaptive.h"
#include "ControlSpace3dAdaptive.h"

unsigned int ControlDataMatrix3d::m_intersection_counter = 0;
ControlDataMatrix3d ControlDataMatrix3d::identity(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

void DMetric3d::setToIdentity()
{
	m_matrix.setIdentity();
	m_reverse.setIdentity();
}

DPoint3d DMetric3d::transformMStoRS(const DMPoint3d &pt) const
{
	return m_matrix.multiplyMStoRS(pt);
}

DMPoint3d DMetric3d::transformRStoMS(const DPoint3d &pt) const
{
	return m_reverse.multiplyRStoMS(pt);
}

DVector3d DMetric3d::transformMStoRS(const DMVector3d &pt) const
{
	return m_matrix.multiplyMStoRS(pt);
}

DMVector3d DMetric3d::transformRStoMS(const DVector3d &pt) const
{
	return m_reverse.multiplyRStoMS(pt);
}
/*
void DMetric3d::setData(const ControlDataStretch3d& data)
{
	m_matrix = stretchToMatrix(data);
	m_reverse = m_matrix.inverse();
}
*/
void DMetric3d::setData(const ControlDataMatrix3d& data)
{
	m_matrix = data;
	m_reverse = m_matrix.inverse();
}

ControlDataMatrix3d DMetric3d::stretchToMatrix(const ControlDataStretch3d& stretch)
{
/*
	DVector3d e0 = DVector3d::v_ox;
	DVector3d e1 = DVector3d::v_oy;
	DVector3d e2 = DVector3d::v_oz;

	if (stretch.ax != 0.0) {
		double sa = sin(stretch.ax);
		double ca = cos(stretch.ax);
		e0 = e0.rotatedAroundAxis(Axis::X, sa, ca);
		e1 = e1.rotatedAroundAxis(Axis::X, sa, ca);
		e2 = e2.rotatedAroundAxis(Axis::X, sa, ca);
	}

	if (stretch.ay != 0.0) {
		double sa = sin(stretch.ay);
		double ca = cos(stretch.ay);
		e0 = e0.rotatedAroundAxis(Axis::Y, sa, ca);
		e1 = e1.rotatedAroundAxis(Axis::Y, sa, ca);
		e2 = e2.rotatedAroundAxis(Axis::Y, sa, ca);
	}

	if (stretch.az != 0.0) {
		double sa = sin(stretch.az);
		double ca = cos(stretch.az);
		e0 = e0.rotatedAroundAxis(Axis::Z, sa, ca);
		e1 = e1.rotatedAroundAxis(Axis::Z, sa, ca);
		e2 = e2.rotatedAroundAxis(Axis::Z, sa, ca);
	}
	double d[] = { stretch.lx, stretch.ly, stretch.lz };
	return CDM3d(e0, e1, e2, d);
*/

	double sin_ax = sin(stretch.ax);
	double sin2_ax = sin_ax*sin_ax;
	double sin_ay = sin(stretch.ay);
	double sin2_ay = sin_ay*sin_ay;
	double sin_az = sin(stretch.az);
	double sin2_az = sin_az*sin_az;
	double cos_ax = cos(stretch.ax);
	double cos2_ax = cos_ax*cos_ax;
	double cos_ay = cos(stretch.ay);
	double cos2_ay = cos_ay*cos_ay;
	double cos_az = cos(stretch.az);
	double cos2_az = cos_az*cos_az;

	double ayy_x = sin_ax*cos_az - cos_ax*sin_ay*sin_az;
	double ayy_y = cos_ax*cos_az + sin_ax*sin_ay*sin_az;
	double azz_x = sin_ax*sin_az + cos_ax*sin_ay*cos_az;
	double azz_y = cos_ax*sin_az - sin_ax*sin_ay*cos_az;
	double ayy2_x = ayy_x * ayy_x;
	double ayy2_y = ayy_y * ayy_y;
	double azz2_x = azz_x * azz_x;
	double azz2_y = azz_y * azz_y;

	return ControlDataMatrix3d(
		stretch.lx * cos2_ax * cos2_ay + stretch.ly * sin2_ax * cos2_ay +	//dxx
			stretch.lz * sin2_ay, 
		stretch.lx * ayy2_x + stretch.ly * ayy2_y +							// dyy
			stretch.lz * cos2_ay * sin2_az,    
		stretch.lx * azz2_x + stretch.ly * azz2_y +							// dzz
			stretch.lz * cos2_ay * cos2_az,
		stretch.lx * cos_ax * cos_ay * ayy_x -								// dxy
			stretch.ly * sin_ax * cos_ay * ayy_y + 
			stretch.lz * cos_ay * sin_ay * sin_az, 
		stretch.lx * cos_ax * cos_ay * azz_x -								// dxz
			stretch.ly * sin_ax * cos_ay * azz_y - 
			stretch.lz * cos_ay * sin_ay * cos_az, 
		stretch.lx * ayy_x * azz_x +										// dyz
			stretch.ly * ayy_y * azz_y - 
			stretch.lz * cos2_ay * sin_az * cos_az);
}

ControlDataMatrix3d DMetric3d::stretchToMatrixWithAdjust(const ControlDataStretch3d& stretch)
{
	ControlDataStretch3d cds = stretch;
	adjustLengths(cds.lx, cds.ly, cds.lz);
	return stretchToMatrix(cds);
}

ControlDataMatrix3d::ControlDataMatrix3d(const DVector3d& e0, const DVector3d& e1, 
		const DVector3d& e2, double d[])
{
	DMatrix3d e;
	e.setColumn(0, e0);
	e.setColumn(1, e1);
	e.setColumn(2, e2);

	DMatrix3d mt = e;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			mt.m[i][j] *= d[j];
	mt = mt*e.transposed();

	assert(abs(mt.m[0][1]) - abs(mt.m[1][0]) < SMALL_NUMBER);
	assert(abs(mt.m[0][2]) - abs(mt.m[2][0]) < SMALL_NUMBER);
	assert(abs(mt.m[2][1]) - abs(mt.m[1][2]) < SMALL_NUMBER);

	dxx = mt.m[0][0];
	dyy = mt.m[1][1];
	dzz = mt.m[2][2];
	dxy = mt.m[0][1];
	dxz = mt.m[0][2];
	dyz = mt.m[1][2];
}

void ControlDataMatrix3d::setEigensystem(const DMatrix3d& e, double d[])
{
	DMatrix3d mt = e;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			mt.m[i][j] *= d[j];
	mt = mt*e.transposed();

	assert(abs(mt.m[0][1]) - abs(mt.m[1][0]) < SMALL_NUMBER);
	assert(abs(mt.m[0][2]) - abs(mt.m[2][0]) < SMALL_NUMBER);
	assert(abs(mt.m[2][1]) - abs(mt.m[1][2]) < SMALL_NUMBER);

	dxx = mt.m[0][0];
	dyy = mt.m[1][1];
	dzz = mt.m[2][2];
	dxy = mt.m[0][1];
	dxz = mt.m[0][2];
	dyz = mt.m[1][2];
}

void ControlDataMatrix3d::setEigensystem(DVector3d ev[], double d[])
{
	DMatrix3d e;
	e.setColumn(0, ev[0]);
	e.setColumn(1, ev[1]);
	e.setColumn(2, ev[2]);

	DMatrix3d mt = e;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			mt.m[i][j] *= d[j];
	mt = mt*e.transposed();

	assert(abs(mt.m[0][1]) - abs(mt.m[1][0]) < SMALL_NUMBER);
	assert(abs(mt.m[0][2]) - abs(mt.m[2][0]) < SMALL_NUMBER);
	assert(abs(mt.m[2][1]) - abs(mt.m[1][2]) < SMALL_NUMBER);

	dxx = mt.m[0][0];
	dyy = mt.m[1][1];
	dzz = mt.m[2][2];
	dxy = mt.m[0][1];
	dxz = mt.m[0][2];
	dyz = mt.m[1][2];
}

void ControlDataMatrix3d::countEigenvector(double d, int ie, DMatrix3d& eigenvectors) const
{
	double w[3] = {
		(dyy-d)*(dzz-d) - dyz*dyz,
		(dxx-d)*(dzz-d) - dxz*dxz,
		(dxx-d)*(dyy-d) - dxy*dxy}; 

	int mi = (abs(w[0])>abs(w[1])) ? 0 : 1;
	if(abs(w[2]) > abs(w[mi])) mi = 2;
//	assert(abs(w[mi]) > SMALL_NUMBER);

	eigenvectors.m[mi][ie] = 1.0;
	double w0, w1;
	int itab[][2] = {{1,2},{0,2},{0,1}};

	switch(mi){
		case 0:	
			w0 = dxz*dyz - dxy*(dzz-d);
			w1 = dxy*dyz - dxz*(dyy-d);
			break;
		case 1:	
			w0 = dxz*dyz - dxy*(dzz-d);
			w1 = dxy*dxz - (dxx-d)*dyz;
			break;
		case 2:	
		default:
			w0 = dxy*dyz - dxz*(dyy-d);
			w1 = dxy*dxz - (dxx-d)*dyz;
			break;
	}

	eigenvectors.m[itab[mi][0]][ie] = w0/w[mi];
	eigenvectors.m[itab[mi][1]][ie] = w1/w[mi];

	double f = 1.0 / sqrt(sqr(eigenvectors.m[0][ie])+sqr(eigenvectors.m[1][ie])+sqr(eigenvectors.m[2][ie]));
	for(int i = 0; i < 3; i++)	eigenvectors.m[i][ie] *= f;
}

bool ControlDataMatrix3d::eigenvalues(double d[]) const
{
	// normalize coefficients -> better numerical stability
	double sf = abs(dxx);
	if(abs(dyy) > sf) sf = abs(dyy);
	if(abs(dzz) > sf) sf = abs(dzz);
	if(abs(dxy) > sf) sf = abs(dxy);
	if(abs(dxz) > sf) sf = abs(dxz);
	if(abs(dyz) > sf) sf = abs(dyz);

	double sfr = 1.0 / sf;
	double sdxx = sfr*dxx;
	double sdyy = sfr*dyy;
	double sdzz = sfr*dzz;
	double sdxy = sfr*dxy;
	double sdxz = sfr*dxz;
	double sdyz = sfr*dyz;

	double coeff[4] = {
		2*sdxy*sdyz*sdxz - sdxz*sdyy*sdxz - sdxx*sdyz*sdyz - sdxy*sdxy*sdzz + sdxx*sdyy*sdzz,
		sdxy*sdxy - sdxx*sdyy + sdxz*sdxz + sdyz*sdyz - sdxx*sdzz - sdyy*sdzz,
		sdxx + sdyy + sdzz,
		-1.0};

	int dct = DMatrix3d::countCubicRoots(coeff, d);
	assert(dct == 3);
	d[0] *= sf;
	d[1] *= sf;
	d[2] *= sf;
	return (dct == 3);
}

double ControlDataMatrix3d::minEigenvalue() const
{
	double d[3];
	eigenvalues(d);
	return std::min(std::min(d[0],d[1]),d[2]);
}

double ControlDataMatrix3d::maxEigenvalue() const
{
	double d[3];
	eigenvalues(d);
	return std::max(std::max(d[0],d[1]),d[2]);
}

bool ControlDataMatrix3d::eigensystem(DMatrix3d& e, double d[]) const
{
	bool result = eigenvalues(d);
	if(!result) return false;

	const double SAME_EPS = 1e-3;
	bool same_01 = (abs(1.0-d[0]/d[1]) < SAME_EPS);
	bool same_02 = (abs(1.0-d[0]/d[2]) < SAME_EPS);
	bool same_12 = (abs(1.0-d[1]/d[2]) < SAME_EPS);

	if(same_01 && same_02){
		// all eignevalues equal
		e.setIdentity();
		return true;
	}else if(same_01){
		countEigenvector(d[2], 2, e);
		// set 0,1 as orthonormal
		DMatrix3d::countOrthonormalEigenvectors(2, 0, 1, e, true);
	}else if(same_02){
		countEigenvector(d[1], 1, e);
		// set 0,2 as orthonormal
		DMatrix3d::countOrthonormalEigenvectors(1, 0, 2, e, true);
	}else if(same_12){
		countEigenvector(d[0], 0, e);
		// set 1,2 as orthonormal
		DMatrix3d::countOrthonormalEigenvectors(0, 1, 2, e, true);
	}else{
		// count all eigenvectors[0..2][i]
		for(int i = 0; i < 3; i++)
			countEigenvector(d[i], i, e);
	}

	// check sign of determinant
/*
	if(e.det() < 0){
		// switch columns 0 and 1 (and respecting eigenvalues),
		double d1 = d[0]; d[0] = d[1]; d[1] = d1;
		for(int i = 0; i < 3; i++){
			d1 = e.m[i][0]; e.m[i][0] = e.m[i][1]; e.m[i][1] = d1;
		}
	}
*/
	return true;
}

ControlDataMatrix3d ControlDataMatrix3d::transformationToTensor() const
{
	return ControlDataMatrix3d(dxx*dxx+dxy*dxy+dxz*dxz, dxy*dxy+dyy*dyy+dyz*dyz,
		dxz*dxz+dyz*dyz+dzz*dzz, dxx*dxy+dxy*dyy+dxz*dyz,
		dxx*dxz+dxy*dyz+dxz*dzz, dxy*dxz+dyy*dyz+dyz*dzz);
}

ControlDataMatrix3d ControlDataMatrix3d::tensorToTransformation() const
{
	DMatrix3d e;
	double d[3];
	bool success = eigensystem(e, d);
	assert(success);
	if(!success) return ControlDataMatrix3d(*this);
	assert(abs(1.0 - abs(e.det())) < 1e-5); // should be "rotation" matrix
	assert((d[0] > 0.0) && (d[1] > 0.0) && (d[2] > 0.0));
	double d0 = sqrt(d[0]);
	double d1 = sqrt(d[1]);
	double d2 = sqrt(d[2]);
	// e.Identity(d).Transpose(e)
	return ControlDataMatrix3d(
		sqr(e.m[0][0])*d0 + sqr(e.m[0][1])*d1 + sqr(e.m[0][2])*d2,
		sqr(e.m[1][0])*d0 + sqr(e.m[1][1])*d1 + sqr(e.m[1][2])*d2,
		sqr(e.m[2][0])*d0 + sqr(e.m[2][1])*d1 + sqr(e.m[2][2])*d2,
		e.m[0][0]*e.m[1][0]*d0 + e.m[0][1]*e.m[1][1]*d1 + e.m[0][2]*e.m[1][2]*d2,
		e.m[0][0]*e.m[2][0]*d0 + e.m[0][1]*e.m[2][1]*d1 + e.m[0][2]*e.m[2][2]*d2,
		e.m[1][0]*e.m[2][0]*d0 + e.m[1][1]*e.m[2][1]*d1 + e.m[1][2]*e.m[2][2]*d2);
}

bool ControlDataMatrix3d::setMinimum(const ControlDataMatrix3d& m, double * d_max)
{
	bool any_changes = false;

	if(abs(dxy) < SMALL_NUMBER && abs(m.dxy) < SMALL_NUMBER &&
		abs(dxz) < SMALL_NUMBER && abs(m.dxz) < SMALL_NUMBER &&
		abs(dyz) < SMALL_NUMBER && abs(m.dyz) < SMALL_NUMBER)
	{
		if(m.dxx < dxx){ any_changes = true; dxx = m.dxx; }
		if(m.dyy < dyy){ any_changes = true; dyy = m.dyy; }
		if(m.dzz < dzz){ any_changes = true; dzz = m.dzz; }
		return any_changes;
	}

	m_intersection_counter++;

	ControlDataMatrix3d mt11 = this->transformationToTensor();
	ControlDataMatrix3d mt22 = m.transformationToTensor();

	DMatrix3d n = mt11.inverse() * mt22;
	DMatrix3d e;
	double d[3];

	bool success = n.eigensystem(e, d);
	assert(success);
	if(!success) return false;

	for(int i = 0; i < 3; i++){
		d[i] = std::min(
			(e.m[0][i]*mt11.dxx + e.m[1][i]*mt11.dxy + e.m[2][i]*mt11.dxz)*e.m[0][i] + 
				(e.m[0][i]*mt11.dxy + e.m[1][i]*mt11.dyy + e.m[2][i]*mt11.dyz)*e.m[1][i] +
				(e.m[0][i]*mt11.dxz + e.m[1][i]*mt11.dyz + e.m[2][i]*mt11.dzz)*e.m[2][i],
			(e.m[0][i]*mt22.dxx + e.m[1][i]*mt22.dxy + e.m[2][i]*mt22.dxz)*e.m[0][i] + 
				(e.m[0][i]*mt22.dxy + e.m[1][i]*mt22.dyy + e.m[2][i]*mt22.dyz)*e.m[1][i] +
				(e.m[0][i]*mt22.dxz + e.m[1][i]*mt22.dyz + e.m[2][i]*mt22.dzz)*e.m[2][i]);
	}

	if (d_max != nullptr) {
		*d_max = sqrt(std::max(d[0], std::max(d[1], d[2])));
	}

	DMatrix3d p = e.inverse();

	DMatrix3d dm33 = p.transposed();
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			dm33.m[i][j] *= d[j];
	dm33 = dm33*p;

	assert(abs(dm33.m[0][1]) - abs(dm33.m[1][0]) < SMALL_NUMBER);
	assert(abs(dm33.m[0][2]) - abs(dm33.m[2][0]) < SMALL_NUMBER);
	assert(abs(dm33.m[2][1]) - abs(dm33.m[1][2]) < SMALL_NUMBER);
	ControlDataMatrix3d mt33(dm33.m[0][0], dm33.m[1][1], dm33.m[2][2], 
		dm33.m[0][1], dm33.m[0][2], dm33.m[1][2]);

	ControlDataMatrix3d m33 = mt33.tensorToTransformation();

#ifdef _DEBUG
	ControlDataMatrix3d mt33x = m33.transformationToTensor();
	double check_diff = mt33.countDifferenceRR(mt33x);
	assert(check_diff < 0.001);
#endif

	any_changes = this->countDifferenceRR(m33) > SMALL_NUMBER;
	if(any_changes)	*this = m33;

	return any_changes;
}

double ControlDataMatrix3d::countDifferenceRR(const ControlDataMatrix3d& data) const {
	const ControlDataMatrix3d mr = data.inverse();
	const ControlDataMatrix3d r = this->inverse();

	// || (M11 * M22.inv - I) + (M11.inv * M22 - I) ||

	double rr_xx = r.dxx*data.dxx + r.dxy*data.dxy + r.dxz*data.dxz
		+ mr.dxx*dxx + mr.dxy*dxy + mr.dxz*dxz;
	double rr_xy = r.dxx*data.dxy + r.dxy*data.dyy + r.dxz*data.dyz
		+ mr.dxx*dxy + mr.dxy*dyy + mr.dxz*dyz;
	double rr_xz = r.dxx*data.dxz + r.dxy*data.dyz + r.dxz*data.dzz
		+ mr.dxx*dxz + mr.dxy*dyz + mr.dxz*dzz;
	double rr_yy = r.dxy*data.dxy + r.dyy*data.dyy + r.dyz*data.dyz
		+ mr.dxy*dxy + mr.dyy*dyy + mr.dyz*dyz;
	double rr_yz = r.dxy*data.dxz + r.dyy*data.dyz + r.dyz*data.dzz
		+ mr.dxy*dxz + mr.dyy*dyz + mr.dyz*dzz;
	double rr_zz = r.dxz*data.dxz + r.dyz*data.dyz + r.dzz*data.dzz
		+ mr.dxz*dxz + mr.dyz*dyz + mr.dzz*dzz;

	return sqrt(sqr(rr_xx-2.0) + sqr(rr_yy-2.0) + sqr(rr_zz-2.0) 
		+ sqr(rr_xy) + sqr(rr_xz) + sqr(rr_yz));
}

void ControlDataMatrix3d::switchEigenvectors(DMatrix3d& e, double d[], int k0, int k1)
{
	double t = d[k0]; d[k0] = d[k1]; d[k1] = t;
	for(int i = 0; i < 3; i++){
		t = e.m[i][k0]; e.m[i][k0] = e.m[i][k1]; e.m[i][k1] = t;
	}
}

void ControlDataMatrix3d::fixEigensystem(DMatrix3d& e, double d[])
{
	int k = (abs(e.m[0][0]) > abs(e.m[0][1])) ? 0 : 1;
	if(abs(e.m[0][2]) > abs(e.m[0][k])) k = 2;
	if(k != 0) switchEigenvectors(e, d, 0, k);
	if(abs(e.m[1][2]) > abs(e.m[1][1]))
		switchEigenvectors(e, d, 1, 2);
}

DRect ControlDataMatrix3d::prepareGraphics(
	DataVector<DPoint2d> & pts0, DataVector<DPoint2d> & pts1, 
	DataVector<bool> & front) const
{
	const int RANGE_U = 40;
	const int RANGE_V = RANGE_U/2;
	DPoint3d p[RANGE_U][RANGE_V+1];
	DPoint2d p2[RANGE_U][RANGE_V+1];
	const double du = 2*PI/RANGE_U;
	const double dv = PI/RANGE_V;

	for(int i = 0; i < RANGE_U; i++){
		double u = (i+1) * du;
		for(int j = 0; j <= RANGE_V; j++){
			double v = j*dv;
			p[i][j].x = cos(u)*sin(v);
			p[i][j].y = sin(u)*sin(v);
			p[i][j].z = cos(v);
		}
	}

	// transform
	DMatrix3d e;
	double d[3];
	bool success = eigensystem(e, d);
	assert(success);
	if(!success) return DRect();
//	fixEigensystem(e, d);
	DMatrix3d tr;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			tr.m[j][i] = e.m[j][i] * d[i];

	tr = tr * e.transposed();

	for(int i = 0; i < RANGE_U; i++)
		for(int j = 0; j <= RANGE_V; j++)
			p[i][j] = tr * p[i][j];

	// project to 2D
	const double f = 0.2;
	DRect rect;
	for(int i = 0; i < RANGE_U; i++)
		for(int j = 0; j <= RANGE_V; j++){
			p2[i][j].x = p[i][j].x + f * p[i][j].y;
			p2[i][j].y = p[i][j].z + f * p[i][j].y;
			rect.addPoint(p2[i][j]);
		}

	// draw
	pts0.prepare(2*RANGE_U * RANGE_V);
	pts1.prepare(2*RANGE_U * RANGE_V);
	front.prepare(2*RANGE_U * RANGE_V);
	for(int i = 0; i < RANGE_U; i++){
		int i1 = (i+1)%RANGE_U;
		for(int j = 1; j <= RANGE_V; j++){
			bool is_front = ((p2[i1][j] - p2[i][j]).crossProduct(p2[i][j-1] - p2[i][j]) >= 0.0);
			front.add(is_front);
			pts0.add(p2[i][j-1]);
			pts1.add(p2[i][j]);
			if(j < RANGE_V){
				front.add(is_front);
				pts0.add(p2[i][j]);
				pts1.add(p2[i1][j]);
			}
		}
	}

	return rect;
}

void ControlDataMatrix3d::drawInterpolationEllipsoids(
	const ControlDataMatrix3d& cdm1, const ControlDataMatrix3d& cdm2)
{
	const int STEP = 6;

	DataVector<DPoint2d> pts0[STEP+1];
	DataVector<DPoint2d> pts1[STEP+1];
	DataVector<bool> front[STEP+1];
	DRect brect[STEP+1];
	
	brect[0] = cdm1.prepareGraphics(pts0[0], pts1[0], front[0]);
	brect[STEP] = cdm2.prepareGraphics(pts0[STEP], pts1[STEP], front[STEP]);
	double dt = 1.0 / STEP;
	for(int i = 1; i < STEP; i++){
		double t = i * dt;
		ControlDataMatrix3d cdm = cdm1 * (1-t) + cdm2 * t;
		brect[i] = cdm.prepareGraphics(pts0[i], pts1[i], front[i]);
	}

	DRect total_rect = brect[0];
	DVector2d dv[STEP+1];
	//dv[0] = DVector2d::zero;

	double ddx[] = {0.6, 0.9, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	for(int i = 1 ; i <= STEP; i++){
		dv[i] = dv[i-1] + DVector2d(ddx[i-1]*brect[i-1].getDX(), 0.0);
		brect[i].transpose(dv[i]);
		total_rect.addRect(brect[i]);
	}

	EPSFile epsf("metric3-interpolation.eps", total_rect.x0, total_rect.x1, total_rect.y0, total_rect.y1);
	for(int i = 0; i <= STEP; i++){
		int ct = pts0[i].countInt();
		for(int j = 0; j < ct; j++)
			if(!front[i][j]) epsf.drawLineGray(pts0[i][j] + dv[i], pts1[i][j] + dv[i], 0.8);
		if(i == 0 || i == STEP){
			for(int j = 0; j < ct; j++)
				if(front[i][j])	 epsf.drawLineGray(pts0[i][j] + dv[i], pts1[i][j] + dv[i], 0.0);
		}else{
			for(int j = 0; j < ct; j++)
				if(front[i][j])	 epsf.drawLineRGB(pts0[i][j] + dv[i], pts1[i][j] + dv[i], 0.7, 0.0, 0.0);
		}
	}

}

void ControlDataMatrix3d::drawIntersectionEllipsoids(
	const ControlDataMatrix3d& cdm1, const ControlDataMatrix3d& cdm2)
{
	DataVector<DPoint2d> pts0[3];
	DataVector<DPoint2d> pts1[3];
	DataVector<bool> front[3];

	DRect total_rect = cdm1.prepareGraphics(pts0[0], pts1[0], front[0]);
	total_rect.addRect(cdm2.prepareGraphics(pts0[1], pts1[1], front[1]));
	ControlDataMatrix3d cdm = cdm1;
	cdm.setMinimum(cdm2);
	total_rect.addRect(cdm.prepareGraphics(pts0[2], pts1[2], front[2]));

	EPSFile epsf("metric3-intersection.eps", total_rect.x0, total_rect.x1, total_rect.y0, total_rect.y1);
	double gray[] = {0.6, 0.6, 0.0};
	for(int i = 0; i < 3; i++){
		int ct = pts0[i].countInt();
		for(int j = 0; j < ct; j++)
			if(!front[i][j]) epsf.drawLineGray(pts0[i][j], pts1[i][j], 0.8);
		if(i == 2){
			for(int j = 0; j < ct; j++)
				if(front[i][j])	 epsf.drawLineRGB(pts0[i][j], pts1[i][j], 0.7, 0.0, 0.0);
		}else{
			for(int j = 0; j < ct; j++)
				if(front[i][j])	 epsf.drawLineGray(pts0[i][j], pts1[i][j], gray[i]);
		}
	}

}

ControlDataMatrix3d ControlDataMatrix3d::countMetric(const DPoint3d& pt0, const DPoint3d& pt1, 
		const DPoint3d& pt2, const DPoint3d& pt3)
{
	return countMetricTensor(pt0, pt1, pt2, pt3).tensorToTransformation();
}

#define SIGN(a,b) ((b) > 0 ? abs(a) : -abs(a))

ControlDataMatrix3d ControlDataMatrix3d::countMetricTensor(const DPoint3d& pt0, const DPoint3d& pt1, 
		const DPoint3d& pt2, const DPoint3d& pt3)
{
	DVector3d vAB = pt1-pt0;
	DVector3d vAC = pt2-pt0;
	DVector3d vAD = pt3-pt0;
	DVector3d vBC = pt2-pt1;
	DVector3d vBD = pt3-pt1;
	DVector3d vCD = pt3-pt2;

	const int N = 6;

	double A[N][N] = {
		{sqr(vAB.x), sqr(vAB.y), sqr(vAB.z), 2*vAB.x*vAB.y, 2*vAB.x*vAB.z, 2*vAB.y*vAB.z},
		{sqr(vAC.x), sqr(vAC.y), sqr(vAC.z), 2*vAC.x*vAC.y, 2*vAC.x*vAC.z, 2*vAC.y*vAC.z},
		{sqr(vAD.x), sqr(vAD.y), sqr(vAD.z), 2*vAD.x*vAD.y, 2*vAD.x*vAD.z, 2*vAD.y*vAD.z},
		{sqr(vBC.x), sqr(vBC.y), sqr(vBC.z), 2*vBC.x*vBC.y, 2*vBC.x*vBC.z, 2*vBC.y*vBC.z},
		{sqr(vBD.x), sqr(vBD.y), sqr(vBD.z), 2*vBD.x*vBD.y, 2*vBD.x*vBD.z, 2*vBD.y*vBD.z},
		{sqr(vCD.x), sqr(vCD.y), sqr(vCD.z), 2*vCD.x*vCD.y, 2*vCD.x*vCD.z, 2*vCD.y*vCD.z}
	};

	double b[N] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; 
	double c[N], d[N];

	// qrdcmp
	bool sing = false;
	for(int k = 0; k < N-1; k++){
		double scale=0.0;
		for(int i = k; i < N; i++) scale = std::max(scale, abs(A[i][k]));
		if(scale == 0.0){	// Singular case.
			sing = true;
			c[k] = d[k] = 0.0;
		}else{	// Form Q k and Q k . A.
			for(int i = k; i < N; i++) A[i][k] /= scale;
			double sum = 0.0;
			for(int i = k; i < N; i++) sum += A[i][k] * A[i][k];
			double sigma=SIGN(sqrt(sum), A[k][k]);
			A[k][k] += sigma;
			c[k] = sigma * A[k][k];
			d[k] = -scale * sigma;
			for(int j = k+1; j < N; j++){
				sum = 0.0;
				for(int i = k; i < N; i++) sum += A[i][k] * A[i][j];
				double tau = sum/c[k];
				for(int i = k; i < N; i++) A[i][j] -= tau * A[i][k];
			}
		}
	}
	d[N-1]=A[N-1][N-1];
	if(d[N-1] == 0.0) sing = true;

	if(sing) LOG4CPLUS_WARN(MeshLog::logger_console, "metric3d from simplex - singular matrix");

	// qrsolve
	for(int j = 0; j < N-1; j++){ // Form Q T. b.
		double sum = 0.0;
		for(int i = j; i < N; i++) sum += A[i][j] * b[i];
		double tau = sum / c[j];
		for(int i = j; i < N;i++) b[i] -= tau * A[i][j];
	}

	b[N-1] /= d[N-1];
	for(int i = N-2; i >= 0; i--){
		double sum = 0.0;
		for(int j = i+1; j < N; j++) sum += A[i][j] * b[j];
		b[i] = (b[i] - sum) / d[i];
	}

	ControlDataMatrix3d cdmm(b[0], b[1], b[2], b[3], b[4], b[5]);
	return cdmm.inverse();
}

ControlDataExtMatrix3dSphere::ControlDataExtMatrix3dSphere(const DPoint3d& pt, 
		double r1, const ControlDataMatrix3d& d1,
		double r2, const ControlDataMatrix3d& d2)
	: ControlDataExtMatrix3d(r1, d1, r2, d2), middle(pt)
{
	assert(r1+r2 > 0.0);
}

ControlDataExtMatrix3dSphere::ControlDataExtMatrix3dSphere(const DPoint3d& pt, 
		double r1, const ControlDataMatrix3d& d1)
	: ControlDataExtMatrix3d(r1, d1), middle(pt)
{
	assert(r1 > 0.0);
}

bool ControlDataExtMatrix3dSphere::isPointWithin(const DPoint3d& pt) const 
{
	return middle.distance2(pt) < sqr(radius1+radius2);
}

ControlDataMatrix3d ControlDataExtMatrix3dSphere::getControlDataMatrix3d(const DPoint3d& pt) const 
{
	double dist = middle.distance(pt);
	if(dist <= radius1 || radius2 == 0.0) 
		return data1;
	else if(dist <= radius1+radius2){
		double t = (dist-radius1)/radius2;
		return data1*(1.0-t) + data2*t;
	}else return data2;
}

ControlDataExtMatrix3dSegment::ControlDataExtMatrix3dSegment(
		const DPoint3d& _pt, const DVector3d& _dv, 
		double r1, const ControlDataMatrix3d& d1,
		double r2, const ControlDataMatrix3d& d2)
	: ControlDataExtMatrix3d(r1, d1, r2, d2), pt0(_pt)
{
	assert(r1+r2 > 0.0);
	dv[0] = _dv;
	dv[0].orthonormalVectors(dv[1], dv[2]);
	dv[1] *= (r1+r2);
	dv[2] *= (r1+r2);
	dr = r1 / (r1+r2);
	dm.setColumn(0, dv[0]);
	dm.setColumn(1, dv[1]);
	dm.setColumn(2, dv[2]);
	w = dm.det();
}

ControlDataExtMatrix3dSegment::ControlDataExtMatrix3dSegment(
		const DPoint3d& _pt, const DVector3d& _dv, 
		double r1, const ControlDataMatrix3d& d1)
	: ControlDataExtMatrix3d(r1, d1), pt0(_pt)
{
	assert(r1 > 0.0);
	dv[0] = _dv;
	dv[0].orthonormalVectors(dv[1], dv[2]);
	if(r1 > 0.0) dv[1] *= r1;
	if(r1 > 0.0) dv[2] *= r1;
	dr = 1.0;
	dm.setColumn(0, dv[0]);
	dm.setColumn(1, dv[1]);
	dm.setColumn(2, dv[2]);
	w = dm.det();
}

bool ControlDataExtMatrix3dSegment::isPointWithin(const DPoint3d& pt) const 
{
	const DVector3d dpt = pt - pt0;
	bool within_cylinder = true;
	DMatrix3d m = dm;
	double t0[3] = {0.0, -1.0, -1.0};
	for(int i = 0; i < 3; i++){
		m.setColumn(i, dpt);
		double t = m.det() / w;
		m.setColumn(i, dv[i]);
		if(t < t0[i] || t > 1.0){
			within_cylinder = false;
			break;
		}
	}
	if(within_cylinder) return true;
	// additional spherical endings
	if(pt.distance2(pt0) <= sqr(totalRadius())) return true;
	if(pt.distance2(pt0+dv[0]) <= sqr(totalRadius())) return true;
	return false;
}

ControlDataMatrix3d ControlDataExtMatrix3dSegment::getControlDataMatrix3d(const DPoint3d& pt) const 
{
	const DVector3d dpt = pt - pt0;
	double t[3];
	DMatrix3d m = dm;
	for(int i = 0; i < 3; i++){
		m.setColumn(i, dpt);
		t[i] = m.det() / w;
		m.setColumn(i, dv[i]);
	}
	double tr = 0.0;
	if(t[0] < 0.0 || t[0] > 1.0 || t[1] < -1.0 || t[1] > 1.0 || t[2] < -1.0 || t[2] > 1.0){
		tr = pt.distance(pt0) / (totalRadius());
		if(tr > 1.0) tr = pt.distance(pt0+dv[0]) / (totalRadius());
	}else{
		tr = sqrt(t[1]*t[1] + t[2]*t[2]);
	}
	if(tr < dr || dr == 1.0) return data1;
	else if(tr <= 1.0){
		double s = (tr-dr)/(1.0-dr);
		return data1*(1.0-s) + data2*s;
	}else return data2;
}

ControlDataExtMatrix3dTriangle::ControlDataExtMatrix3dTriangle(
		const DPoint3d& _pt, 
		const DVector3d& _dv1, const DVector3d& _dv2,
		double r1, const ControlDataMatrix3d& d1,
		double r2, const ControlDataMatrix3d& d2)
	: ControlDataExtMatrix3d(r1, d1, r2, d2), pt0(_pt)
{
	assert(r1+r2 > 0.0);
	dv[0] = _dv1;
	dv[1] = _dv2;
	dv[2] = _dv1.crossProduct(_dv2).normalized() * (r1+r2);
	dr = r1 / (r1+r2);
	dm.setColumn(0, dv[0]);
	dm.setColumn(1, dv[1]);
	dm.setColumn(2, dv[2]);
	w = dm.det();
}

ControlDataExtMatrix3dTriangle::ControlDataExtMatrix3dTriangle(
		const DPoint3d& _pt, 
		const DVector3d& _dv1, const DVector3d& _dv2,
		double r1, const ControlDataMatrix3d& d1)
	: ControlDataExtMatrix3d(r1, d1), pt0(_pt)
{
	assert(r1 >= 0.0);
	dv[0] = _dv1;
	dv[1] = _dv2;
	dv[2] = _dv1.crossProduct(_dv2).normalized();
	if(r1 > 0)  dv[2] *= r1;
	dr = 1.0;
	dm.setColumn(0, dv[0]);
	dm.setColumn(1, dv[1]);
	dm.setColumn(2, dv[2]);
	w = dm.det();
}

bool ControlDataExtMatrix3dTriangle::isPointWithin(const DPoint3d& pt) const 
{
	const DVector3d dpt = pt - pt0;
	bool within_prism = true;
	DMatrix3d m = dm;
	double t[3];
	double t0[3] = {0.0, 0.0, -1.0};
	for(int i = 0; i < 3; i++){
		m.setColumn(i, dpt);
		t[i] = m.det() / w;
		m.setColumn(i, dv[i]);
		if(t[i] < t0[i] || t[i] > 1.0){
			within_prism = false;
			break;
		}
	}
	if(t[0] + t[1] > 1.0) within_prism = false;

	return within_prism;
//	if(within_prism) return true;
	// additional spherical endings
//	if(pt.distance2(pt0) <= sqr(totalRadius())) return true;
//	if(pt.distance2(pt0+dv[0]) <= sqr(totalRadius())) return true;
//	return false;
}

ControlDataMatrix3d ControlDataExtMatrix3dTriangle::getControlDataMatrix3d(const DPoint3d& pt) const 
{
	const DVector3d dpt = pt - pt0;
	double t[3];
	DMatrix3d m = dm;
	for(int i = 0; i < 3; i++){
		m.setColumn(i, dpt);
		t[i] = m.det() / w;
		m.setColumn(i, dv[i]);
	}
	if(t[0] < 0.0 || t[0] > 1.0 || t[1] < 0.0 || t[1] > 1.0 || t[0]+t[1] > 1.0)
		return data2; // outside the triangle

	double tr = abs(t[2]);
	if(tr > 1.0) return data2;
	else if(tr < dr) return data1;
	else{
		double s = (tr-dr)/(1.0-dr);
		return data1*(1.0-s) + data2*s;
	}
}

ostream& operator<<(ostream& os, const ControlDataStretch3d& data){
	string str_lx, str_ly, str_lz;
	string str_ax, str_ay, str_az;
	DEquation::doubleToString(data.lx, DEquation::v_length, str_lx);
	DEquation::doubleToString(data.ly, DEquation::v_length, str_ly);
	DEquation::doubleToString(data.lz, DEquation::v_length, str_lz);
	DEquation::doubleToString(data.ax, DEquation::v_value, str_ax);
	DEquation::doubleToString(data.ay, DEquation::v_value, str_ay);
	DEquation::doubleToString(data.az, DEquation::v_value, str_az);
	os << "[" << str_lx << "," << str_ly << "," << str_lz << "," 
		<< str_ax << "," << str_ay << "," << str_az << "]";
	return os;
}

istream& operator>>(istream& is, ControlDataStretch3d& data){
	string str_lx, str_ly, str_lz;
	string str_ax, str_ay, str_az;
	char z;
	is >> ws >> z; // '['
	getline(is, str_lx, ',');
	getline(is, str_ly, ',');
	getline(is, str_lz, ',');
	getline(is, str_ax, ',');
	getline(is, str_ay, ',');
	getline(is, str_az, ']');
	DEquation eq;
	if(eq.parse(str_lx.c_str()))	data.lx = eq.getValue(0.0);
	if(eq.parse(str_ly.c_str()))	data.ly = eq.getValue(0.0);
	if(eq.parse(str_lz.c_str()))	data.lz = eq.getValue(0.0);
	if(eq.parse(str_ax.c_str()))	data.ax = eq.getValue(0.0);
	if(eq.parse(str_ay.c_str()))	data.ay = eq.getValue(0.0);
	if(eq.parse(str_az.c_str()))	data.az = eq.getValue(0.0);
	return is;
}

bool DMetric3d::adjustMetricGradationNb( ControlNode3d * cn0)
{
	double dr = ControlSpace2dAdaptive::param_gradation_ratio;
	assert(dr >= 1.0);
	if (dr < 1.0) dr = 1.0;

	double max_diff = 0.0;

	std::shared_ptr<ControlNode3d> cn1;
	for (auto it = cn0->nbs->iterator(); it.valid(); it.moveNext()) {
		cn1 = it.item();
		double diff = cn0->control_data.countDifferenceRR(cn1->control_data);
		if (diff > max_diff) max_diff = diff;
	}
	if (max_diff < 0.01) {
		cn0->max_gradation_ratio = 1.0;
		return false;
	}

	struct CnData {
		CnData(ControlNode3d* _cn) : cn1(_cn) {}
	public:
		double checkGradation(ControlNode3d* cn0, double req_gratio , bool & cn0_changed, int phase) {
			const DMetric3d dm0(cn0->control_data);
			const DMetric3d dm1(cn1->control_data);
			DVector3d dv = cn1->coord - cn0->coord;
			double lm0 = dm0.transformRStoMS(dv).length();
			double lm1 = dm1.transformRStoMS(dv).length();
			if (phase == 0 && (lm0 < SMALL_NUMBER || lm1 < SMALL_NUMBER)) {
				too_close = true;
				cn0_changed |= cn0->control_data.setMinimum(cn1->control_data);
				if (cn1->control_data.setMinimum(cn0->control_data)) cn1->wi = -1;
				return 1.0;
			}
			else if (phase == 1 && too_close) {
				if (cn1->control_data.setMinimum(cn0->control_data)) cn1->wi = -1;
				return 1.0;
			}
			else if (too_close) return 1.0;

			double lm_ave = 0.5*(lm0 + lm1);
			assert(phase >= 0 && phase < 3);
			req_hratio[phase] = (req_gratio - 1.0) * lm_ave + 1.0;
			return DMetric3d::getMetricGradation(cn0->control_data, cn1->control_data, dv,
				&hratio0[phase], &hratio1[phase], 0.0);
		}
		bool intersectCn0(ControlNode3d* cn0) {
			return !too_close && (req_hratio[0] < hratio1[0]) && 
				cn0->control_data.setMinimum(cn1->control_data * req_hratio[0]);
		}
		void intersectCn1(ControlNode3d* cn0) {
			if(!too_close && (req_hratio[1] < hratio0[1]) && 
				cn1->control_data.setMinimum(cn0->control_data * req_hratio[1]))
				cn1->wi = -1;
		}
	public:
		ControlNode3d* cn1 = nullptr;
		bool too_close = false;
		double req_hratio[3];
		double hratio0[3], hratio1[3];
	};

	DataVector<std::shared_ptr<CnData>> cn_vec;
	bool cn0_changed = false;

	double req_gratio = ControlSpace2dAdaptive::param_gradation_ratio;
	static const double GRADATION_EPS = 1e-3;

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "= " << cn0->control_data);

	double max_oryg_gratio = 0.0;
	for (auto it = cn0->nbs->iterator(); it.valid(); it.moveNext()) {
		cn1 = it.item();
		auto cnd = std::make_shared<CnData>(cn1.get());
		cn_vec.add(cnd);
		double gratio = cnd->checkGradation(cn0, req_gratio, cn0_changed, 0);
		if (gratio > max_oryg_gratio) max_oryg_gratio = gratio;
	}
	if ((max_oryg_gratio - req_gratio) < GRADATION_EPS)
		return cn0_changed;

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "== " << cn0->control_data);

	cn_vec.forEach([&cn0,&cn0_changed](auto cnd){
		cn0_changed |= cnd->intersectCn0(cn0);
	});

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "=== " << cn0->control_data);

	double max_intersect_gratio = 0.0;
	cn_vec.forEach([&](auto cnd) {
		cnd->checkGradation(cn0, req_gratio, cn0_changed, 1); // update req_hratio/hratio
		cnd->intersectCn1(cn0);
		double gratio = cnd->checkGradation(cn0, req_gratio, cn0_changed, 2);
		if (gratio > max_intersect_gratio) max_intersect_gratio = gratio;
	});
	if (abs(max_intersect_gratio - req_gratio) > GRADATION_EPS) {
		// correction
	}
	
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "==== " << cn0->control_data);

	return cn0_changed;
}

bool_pair DMetric3d::adjustMetricGradation(ControlNode3d * cn0, ControlNode3d * cn1)
{
	CDM3d& cdm0 = cn0->control_data;
	CDM3d& cdm1 = cn1->control_data;

	double diff = cdm0.countDifferenceRR(cdm1);
	if (diff < 0.01) return {false, false}; //std::make_pair(false,false);

	const DMetric3d dm0(cdm0);
	const DMetric3d dm1(cdm1);
	DVector3d dv = cn1->coord - cn0->coord;
	double lm0 = dm0.transformRStoMS(dv).length();
	double lm1 = dm1.transformRStoMS(dv).length();
	if (lm0 < SMALL_NUMBER || lm1 < SMALL_NUMBER)
		return std::make_pair(
			cdm0.setMinimum(cdm1), 
			cdm1.setMinimum(cdm0));

	double lm_ave = 0.5*(lm0 + lm1);

	CDM3d oryg_cdm0 = cdm0;
	CDM3d oryg_cdm1 = cdm1;
	auto result = make_pair(false, false);

	double req_gratio = ControlSpace2dAdaptive::param_gradation_ratio;
	double req_hratio = (req_gratio - 1.0) * lm_ave + 1.0;
	static const double GRADATION_EPS = 1e-3;

	double oryg_hratio0, oryg_hratio1;
	double oryg_gratio = DMetric3d::getMetricGradation(cdm0, cdm1, dv, 
		&oryg_hratio0, &oryg_hratio1, 0.0);
	if ((oryg_gratio - req_gratio) < GRADATION_EPS) return result;

	if (req_hratio < oryg_hratio1)
		result.first = cdm0.setMinimum(oryg_cdm1 * req_hratio);
	if (req_hratio < oryg_hratio0)
		result.second = cdm1.setMinimum(oryg_cdm0 * req_hratio);
	assert(result.first || result.second);

	static string rec_indent = "";

	double intersect_hratio0, intersect_hratio1;
	double intersect_gratio = DMetric3d::getMetricGradation(cdm0, cdm1, dv, 
		&intersect_hratio0, &intersect_hratio1, 0.0);
	if (abs(intersect_gratio - req_gratio) > GRADATION_EPS) {
		double intersect_hratio01 = std::max(intersect_hratio0, intersect_hratio1);

		// interpolate (intersect|oryg)

		const DMetric3d intersect_dm0(cdm0);
		const DMetric3d intersect_dm1(cdm1);
		double intersect_lm0 = intersect_dm0.transformRStoMS(dv).length();
		double intersect_lm1 = intersect_dm1.transformRStoMS(dv).length();
		double intersect_lm_ave = 0.5*(intersect_lm0 + intersect_lm1);

		if (result.first && result.second) {
			// change both
			double k = (req_gratio - 1.0) * intersect_lm_ave / (intersect_hratio01 - 1.0);

			cdm0 *= k;
			cdm1 *= k;
			intersect_gratio = DMetric3d::getMetricGradation(cdm0, cdm1, dv, nullptr, nullptr, 0.0);
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, rec_indent << "check_gr(2x) -> " << intersect_gratio);
		}
		else { // change only one
			double a = 2.0 * intersect_hratio01;
			double mb = (req_gratio - 1.0) * (result.first ? intersect_lm1 : intersect_lm0) + 2.0; // minus
			double mc = (req_gratio - 1.0) * (result.first ? intersect_lm0 : intersect_lm1); // minus
			double delta = mb*mb + 4.0 * a * mc;
			double delta_sqrt = sqrt(delta);

			double k0 = (mb + delta_sqrt) / (2.0*a);
			//double k1 = (mb - delta_sqrt) / (2.0*a);

			if(result.first)  cdm0 *= k0; 
			else cdm1 *= k0;

			double res_hratio0, res_hratio1;
			intersect_gratio = DMetric3d::getMetricGradation(cdm0, cdm1, dv, &res_hratio0, &res_hratio1, 0.0);
			if (abs(intersect_gratio - req_gratio) > GRADATION_EPS) {
				//LOG4CPLUS_INFO(MeshLog::logger_mesh, rec_indent << "check_gr(1x) -> recursive -> " << intersect_gratio);
				rec_indent += " ";
				assert(rec_indent.length() < 10);
				auto rec_result = adjustMetricGradation(cn0, cn1);
				rec_indent.pop_back();
				result.first  |= rec_result.first;
				result.second |= rec_result.second;
				intersect_gratio = DMetric3d::getMetricGradation(cdm0, cdm1, dv, &res_hratio0, &res_hratio1, 0.0);
			}
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, rec_indent << "check_gr(1x) -> " << intersect_gratio);
		}
	}
	else {
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, rec_indent << "check_gr(no) -> " << intersect_gratio);
	}

		//for (int i = 0; i < 10; i++) {
		//	cdm0 = oryg_cdm0;
		//	cdm1 = oryg_cdm1;
		//	result = 0;
		//	double k_next = (k1 - k0)*(req_g_ratio - gr0) / (gr1 - gr0) + k0;
		//	if (i % 2 == 0) {
		//		k = 2 * k_next - k;
		//		if (k <= k0 || k >= k1) k = 0.5*(k0 + k1);
		//	}
		//	else k = k_next;
		//	//assert(k >= k0 && k <= k1);
		//	k_h_ratio = k * h_ratio01 + (1.0 - k) * req_h_ratio;
		//	if ((k_h_ratio < h_ratio1) && cdm0.setMinimum(oryg_cdm1 * k_h_ratio)) result += 1;
		//	if ((k_h_ratio < h_ratio0) && cdm1.setMinimum(oryg_cdm0 * k_h_ratio)) result += 2;
		//	last_gr = DMetric3d::getMetricGradation(cdm0, cdm1, dv);
		//	stat_gr.add(last_gr);
		//	stat_k.add(k);
		//	double dgr = (last_gr - req_g_ratio);
		//	if (abs(dgr) < GRADATION_EPS) break;
		//	if (dgr < 0) {
		//		k0 = k; gr0 = last_gr;
		//	}
		//	else {
		//		k1 = k; gr1 = last_gr;
		//	}
		//}
		//MESHLOGSTAT << "\t" << stat_gr.countInt();
		//for (int i = stat_gr.countInt() - 1; i >= 0; i--)
		//	MESHLOGSTAT << "\t" << stat_gr[i] << "[" << stat_k[i] << "]";
		//MESHLOGSTAT << endl;

	cn0->setGradationUnknown();
	cn1->setGradationUnknown();

	return result;
}

bool DMetric3d::adjustMetricGradationIterativeNb(ControlNode3d * cn0, int step)
{
	CDM3d& cdm0 = cn0->control_data;

	double max_diff = 0.0;
	std::shared_ptr<ControlNode3d> cn1;
	int cn1_count = 0;
	for (auto it = cn0->nbs->iterator(); it.valid(); it.moveNext()) {
		cn1 = it.item();
		cn1_count++;
		double diff = cn0->control_data.countDifferenceRR(cn1->control_data);
		if (diff > max_diff) max_diff = diff;
	}
	if (max_diff < 0.01) {
		cn0->max_gradation_ratio = 1.0;
		return false;
	}

	const DMetric3d dm0(cdm0);

	double req_g_ratio = ControlSpace2dAdaptive::param_gradation_ratio;

	struct CnData
	{
	public:
		CnData(ControlNode3d* cn, const DVector3d& _dv) 
			: cn1(cn), cdm1(cn->control_data), oryg_cdm1(cdm1), dv(_dv) {}
		double computeGradationK1(const CDM3d& cdm0) {
			double gr = DMetric3d::getMetricGradation(cdm0, cdm1, dv, &h_ratio0, &h_ratio1);
			h_ratio01 = std::max(h_ratio0, h_ratio1);
			return gr;
		}
		double adjustGradationK0(const CDM3d& cdm0) {
			CDM3d tmp_cdm0 = cdm0;
			CDM3d tmp_cdm1 = cdm1;
			if (req_h_ratio < h_ratio1) tmp_cdm0.setMinimum(cdm1 * req_h_ratio);
			if (req_h_ratio < h_ratio0) tmp_cdm1.setMinimum(cdm0 * req_h_ratio);
			return DMetric3d::getMetricGradation(tmp_cdm0, tmp_cdm1, dv);
		}
		bool adjustGradationDirect(CDM3d& cdm0, const CDM3d& oryg_cdm0, int _step) {
			bool res0 = cdm0.setMinimum(oryg_cdm1 * req_h_ratio);
			result1 = cdm1.setMinimum(oryg_cdm0 * req_h_ratio);
			return res0;
		}
		void computeKHRatio(double k) {
			k_h_ratio = k * h_ratio01 + (1.0 - k) * req_h_ratio;
		}
		bool adjustGradationKForCN0(CDM3d& cdm0) {
			return (k_h_ratio < h_ratio1) && cdm0.setMinimum(oryg_cdm1 * k_h_ratio);
		}
		double adjustGradationKForCN1(CDM3d& cdm0, const CDM3d& oryg_cdm0) {
			cdm1 = oryg_cdm1;
			result1 = (k_h_ratio < h_ratio0) && cdm1.setMinimum(oryg_cdm0 * k_h_ratio);
			return DMetric3d::getMetricGradation(cdm0, cdm1, dv);
		}
	public:
		ControlNode3d* cn1;
		CDM3d& cdm1;
		CDM3d oryg_cdm1;
		DVector3d dv;
		double req_h_ratio = 0.0;
		double h_ratio0;
		double h_ratio1;
		double h_ratio01;
		double k_h_ratio;
		bool result1 = false;
	};

	CDM3d oryg_cdm0 = cdm0;

	DataVector<std::shared_ptr<CnData>> cn_data_list(cn1_count);

	for (auto it = cn0->nbs->iterator(); it.valid(); it.moveNext()) {
		cn1 = it.item();
		const DMetric3d dm1(cn1->control_data);
		auto cnd = std::make_shared<CnData>(cn1.get(), cn1->coord - cn0->coord);
		double lm0 = dm0.transformRStoMS(cnd->dv).length();
		double lm1 = dm1.transformRStoMS(cnd->dv).length();
		double lm_ave = 0.5*(lm0 + lm1);
		cnd->req_h_ratio = (req_g_ratio - 1.0) * lm_ave + 1.0;
		cn_data_list.add(cnd);
	}

	bool result0 = false;

	if (true) { // iterative
		double max_gr1 = 0.0;
		cn_data_list.forEach([&cdm0,&max_gr1](auto& cnd) {	
			double gr = cnd->computeGradationK1(cdm0); 
			if (gr > max_gr1) max_gr1 = gr;
		});
		if (max_gr1 < 1.01 * req_g_ratio) return false;

		double max_gr0 = 0.0;
		cn_data_list.forEach([&cdm0, &max_gr0](auto& cnd) {
			double gr = cnd->adjustGradationK0(cdm0);
			if (gr > max_gr0) max_gr0 = gr;
		});

		double last_gr = max_gr0;
		double k = 0.0, k0 = 0.0, k1 = 1.0;

		for (int i = 0; i < 6; i++) {
			cdm0 = oryg_cdm0;
			result0 = false;

			double k_next = (k1 - k0)*(req_g_ratio - max_gr0) / (max_gr1 - max_gr0) + k0;
			if (i % 2 == 0) {
				k = 2 * k_next - k;
				if (k <= k0 || k >= k1) k = 0.5*(k0 + k1);
			}
			else k = k_next;

			// all cmd1 -> cmd0
			cn_data_list.forEach([&result0,&cdm0,k](auto& cnd) {
				cnd->computeKHRatio(k);
				result0 |= cnd->adjustGradationKForCN0(cdm0);
			});
			// cmd0 -> all cmd1
			last_gr = 0.0;
			cn_data_list.forEach([&](auto& cnd) {
				double gr = cnd->adjustGradationKForCN1(cdm0, oryg_cdm0);
				if (gr > last_gr) last_gr = gr;
			});

			double dgr = (last_gr - req_g_ratio);
			if (abs(dgr) < SMALL_NUMBER) break;
			if (dgr < 0) {
				k0 = k; max_gr0 = last_gr;
			}
			else {
				k1 = k; max_gr1 = last_gr;
			}
		}
	}
	else {
		cn_data_list.forEach([&](auto& cnd) {
			result0 |= cnd->adjustGradationDirect(cdm0, oryg_cdm0, step);
		});
	}

	bool result_any = result0;
	cn_data_list.forEach([&result_any](auto& cnd) {
		if (cnd->result1) {
			cnd->cn1->setGradationUnknown();
			result_any = true;
		}
	});

	if (result_any) {
		cn0->setGradationUnknown();
	}

	return result_any;
}

double ControlDataMatrix3d::getMaxHRatioIsotropic(
	const ControlDataMatrix3d & cdm1,
	double * h_ratio0_ret, double * h_ratio1_ret) const
{
	const CDM3d & cdm0 = *this;
	assert(cdm0.isotropic());
	assert(cdm1.isotropic());

	double h_ratio0, h_ratio1;
	double max_h_ratio0 = 0.0;
	double max_h_ratio1 = 0.0;
	double d0[] = { cdm0.dxx, cdm0.dyy, cdm0.dzz };
	double d1[] = { cdm1.dxx, cdm1.dyy, cdm1.dzz };
	for (int i = 0; i < 3; i++) {
		if (d1[i] >= d0[i]) {
			h_ratio0 = d1[i] / d0[i];
			if (h_ratio0 > max_h_ratio0) max_h_ratio0 = h_ratio0;
		}
		else {
			h_ratio1 = d0[i] / d1[i];
			if (h_ratio1 > max_h_ratio1) max_h_ratio1 = h_ratio1;
		}
	}

	if (h_ratio0_ret != nullptr) *h_ratio0_ret = max_h_ratio0;
	if (h_ratio1_ret != nullptr) *h_ratio1_ret = max_h_ratio1;

	return std::max(max_h_ratio0, max_h_ratio1);
}

double DMetric3d::getMetricGradation(const CDM3d & cdm0, const CDM3d & cdm1,
	const DVector3d & dv, double * h_ratio0_ret, double * h_ratio1_ret,
	double diff_eps)
{
	if (diff_eps > 0.0 && cdm0.countDifferenceRR(cdm1) < diff_eps) 
		return 1.0;

	const DMetric3d dm0(cdm0);
	const DMetric3d dm1(cdm1);
	double lm0 = dm0.transformRStoMS(dv).length();
	if (lm0 < VERY_SMALL_NUMBER) return LARGE_NUMBER;
	double lm1 = dm1.transformRStoMS(dv).length();
	if (lm1 < VERY_SMALL_NUMBER) return LARGE_NUMBER;
	double lm_ave = 0.5*(lm0 + lm1);

//	double max_h_ratio_red = cdm0.getMaxHRatioReduction(cdm1, h_ratio0_ret, h_ratio1_ret);
	double max_h_ratio = cdm0.getMaxHRatioDirect(cdm1, h_ratio0_ret, h_ratio1_ret);

	double g_ratio = (max_h_ratio - 1.0) / lm_ave + 1.0;

	return g_ratio;
}

double ControlDataMatrix3d::getMaxHRatioDirect(
	const ControlDataMatrix3d & cdm1,
	double * h_ratio0_ret, double * h_ratio1_ret) const
{
	const CDM3d & cdm0 = *this;

	if (cdm0.isotropic() && cdm1.isotropic())
		return getMaxHRatioIsotropic(cdm1, h_ratio0_ret, h_ratio1_ret);

	DMatrix3d e0, e1;
	double d0[3], d1[3];
	bool res = cdm0.eigensystem(e0, d0) && cdm1.eigensystem(e1, d1);
	assert(res == true);
	if (!res) return LARGE_NUMBER;

	static int counter = 0;
	++counter;

	double h_ratio0, h_ratio1;
	double max_h_ratio0 = 0.0;
	double max_h_ratio1 = 0.0;
	for (int i = 0; i < 3; i++) {
		DVector3d ev = e0.column(i);
		//double h0x = (cdm0 * ev).length();
		//assert(abs((cdm0 * ev).length() - d0[i]) < METRIC_SMALL_NUMBER);
		double h1 = (cdm1 * ev).length();
		if (h1 >=  d0[i]) {
			h_ratio0 = h1 / d0[i];
			if (h_ratio0 > max_h_ratio0) max_h_ratio0 = h_ratio0;
		}
		else {
			h_ratio1 = d0[i] / h1;
			if (h_ratio1 > max_h_ratio1) max_h_ratio1 = h_ratio1;
		}

		ev = e1.column(i);
		double h0 = (cdm0 * ev).length();
		//double h1x = (cdm1 * ev).length();
		//assert(abs((cdm1 * ev).length() - d1[i]) < METRIC_SMALL_NUMBER);
		if (h0 >= d1[i]) {
			h_ratio1 = h0 / d1[i];
			if (h_ratio1 > max_h_ratio1) max_h_ratio1 = h_ratio1;
		}
		else {
			h_ratio0 = d1[i] / h0;
			if (h_ratio0 > max_h_ratio0) max_h_ratio0 = h_ratio0;
		}
	}

	if (h_ratio0_ret != nullptr) *h_ratio0_ret = max_h_ratio0;
	if (h_ratio1_ret != nullptr) *h_ratio1_ret = max_h_ratio1;

	return std::max(max_h_ratio0, max_h_ratio1);
}

inline bool ControlDataMatrix3d::isotropic(double eps) const
{
	return abs(dxy) < eps && abs(dxz) < eps && abs(dyz) < eps;
}

double ControlDataMatrix3d::getMaxHRatioReduction(
	const ControlDataMatrix3d & cdm1,
	double * h_ratio0_ret, double * h_ratio1_ret) const
{
	const CDM3d & cdm0 = *this;
	if (cdm0.isotropic() && cdm1.isotropic())
		return getMaxHRatioIsotropic(cdm1, h_ratio0_ret, h_ratio1_ret);

	ControlDataMatrix3d mt00 = cdm0.transformationToTensor();
	ControlDataMatrix3d mt11 = cdm1.transformationToTensor();

	DMatrix3d n = mt00.inverse() * mt11;
	DMatrix3d e;
	double d[3];

	bool success = n.eigensystem(e, d);
	assert(success);
	if (!success) return false;

	double h_ratio0, h_ratio1;
	double max_h_ratio0 = 0.0;
	double max_h_ratio1 = 0.0;
	for (int i = 0; i < 3; i++) {
		double d0 = 
			(e.m[0][i] * mt00.dxx + e.m[1][i] * mt00.dxy + e.m[2][i] * mt00.dxz)*e.m[0][i] +
			(e.m[0][i] * mt00.dxy + e.m[1][i] * mt00.dyy + e.m[2][i] * mt00.dyz)*e.m[1][i] +
			(e.m[0][i] * mt00.dxz + e.m[1][i] * mt00.dyz + e.m[2][i] * mt00.dzz)*e.m[2][i];
		double d1 = 
			(e.m[0][i] * mt11.dxx + e.m[1][i] * mt11.dxy + e.m[2][i] * mt11.dxz)*e.m[0][i] +
			(e.m[0][i] * mt11.dxy + e.m[1][i] * mt11.dyy + e.m[2][i] * mt11.dyz)*e.m[1][i] +
			(e.m[0][i] * mt11.dxz + e.m[1][i] * mt11.dyz + e.m[2][i] * mt11.dzz)*e.m[2][i];
		if (d1 >= d0) {
			h_ratio0 = sqrt(d1 / d0);
			if (h_ratio0 > max_h_ratio0) max_h_ratio0 = h_ratio0;
		}
		else {
			h_ratio1 = sqrt(d0 / d1);
			if (h_ratio1 > max_h_ratio1) max_h_ratio1 = h_ratio1;
		}
	//d[i] = std::min(d0, d1);
	}

//	DMatrix3d p = e.inverse();
//
//	DMatrix3d dm33 = p.transposed();
//	for (int i = 0; i < 3; i++)
//		for (int j = 0; j < 3; j++)
//			dm33.m[i][j] *= d[j];
//	dm33 = dm33*p;
//
//	assert(abs(dm33.m[0][1]) - abs(dm33.m[1][0]) < SMALL_NUMBER);
//	assert(abs(dm33.m[0][2]) - abs(dm33.m[2][0]) < SMALL_NUMBER);
//	assert(abs(dm33.m[2][1]) - abs(dm33.m[1][2]) < SMALL_NUMBER);
//	ControlDataMatrix3d mt33(dm33.m[0][0], dm33.m[1][1], dm33.m[2][2],
//		dm33.m[0][1], dm33.m[0][2], dm33.m[1][2]);
//
//	ControlDataMatrix3d m33 = mt33.tensorToTransformation();
//
//#ifdef _DEBUG
//	ControlDataMatrix3d mt33x = m33.transformationToTensor();
//	double check_diff = mt33.countDifferenceRR(mt33x);
//	assert(check_diff < 0.001);
//#endif
//
//	any_changes = this->countDifferenceRR(m33) > SMALL_NUMBER;
//	if (any_changes)	*this = m33;
//
//	return any_changes;

	if (h_ratio0_ret != nullptr) *h_ratio0_ret = max_h_ratio0;
	if (h_ratio1_ret != nullptr) *h_ratio1_ret = max_h_ratio1;

	return std::max(max_h_ratio0, max_h_ratio1);
}

void DMetric3d::adjustLengths(double & lx, double & ly, double & lz)
{
	double p = mesh_data.getModelDiameter();
	double max_len = p * ControlSpace3dAdaptive::param_max_diameter_ratio;
	double min_len = std::min(p * ControlSpace2dAdaptive::param_min_diameter_ratio, ControlSpace2dAdaptive::param_min_length);

	if(lx > max_len) lx = max_len; 
	else if(lx < min_len) lx = min_len; 

	if(ly > max_len) ly = max_len; 
	else if(ly < min_len) ly = min_len; 

	if(lz > max_len) lz = max_len; 
	else if(lz < min_len) lz = min_len; 

	if(lx < ly && lx < lz){ // lx is minimum
		if((ly / lx) > ControlSpace2dAdaptive::param_stretch_max_ratio)
			ly = lx * ControlSpace2dAdaptive::param_stretch_max_ratio;
		if((lz / lx) > ControlSpace2dAdaptive::param_stretch_max_ratio)
			lz = lx * ControlSpace2dAdaptive::param_stretch_max_ratio;
	}else if(ly < lz){ // ly is minimum
		if((lx / ly) > ControlSpace2dAdaptive::param_stretch_max_ratio)
			lx = ly * ControlSpace2dAdaptive::param_stretch_max_ratio;
		if((lz / ly) > ControlSpace2dAdaptive::param_stretch_max_ratio)
			lz = ly * ControlSpace2dAdaptive::param_stretch_max_ratio;
	}else{ // lz is minimum
		if((lx / lz) > ControlSpace2dAdaptive::param_stretch_max_ratio)
			lx = lz * ControlSpace2dAdaptive::param_stretch_max_ratio;
		if((ly / lz) > ControlSpace2dAdaptive::param_stretch_max_ratio)
			ly = lz * ControlSpace2dAdaptive::param_stretch_max_ratio;
	}
}

bool ControlDataMatrix3d::storeSimple(ostream& os, char end_char) const {
	os << dxx << '\t' << dyy << '\t' << dzz 
		<< dxy << '\t' << dxz << '\t' << dyz 
		<< end_char;
	return os.good();
}
bool ControlDataMatrix3d::readSimple(istream& is) {
	is >> dxx >> dyy >> dzz >> dxy >> dxz >> dyz;
	return is.good();
}

istream& operator>>(istream& is, ControlDataMatrix3d& cdm)
{
	string str_xx, str_yy, str_zz;
	string str_xy, str_xz, str_yz;
	char z;
	is >> ws >> z; // '['
	getline(is, str_xx, ',');
	getline(is, str_yy, ',');
	getline(is, str_zz, ',');
	getline(is, str_xy, ',');
	getline(is, str_xz, ',');
	getline(is, str_yz, ']');
	DEquation eq;
	if(eq.parse(str_xx.c_str()))	cdm.dxx = eq.getValue(0.0);
	if(eq.parse(str_yy.c_str()))	cdm.dyy = eq.getValue(0.0);
	if(eq.parse(str_zz.c_str()))	cdm.dzz = eq.getValue(0.0);
	if(eq.parse(str_xy.c_str()))	cdm.dxy = eq.getValue(0.0);
	if(eq.parse(str_xz.c_str()))	cdm.dxz = eq.getValue(0.0);
	if(eq.parse(str_yz.c_str()))	cdm.dyz = eq.getValue(0.0);
	return is;
}

/// Transforms 2D metric transformation tensor into 3D
bool DMetric3d::projectCDMto3D(const ControlDataMatrix2d& cdm2d, SurfaceConstPtr surface, const DPoint2d& pt2d, 
							   ControlDataMatrix3d& cdm, DVector3d * cdm_e0, DVector3d * cdm_e1)
{
	DMatrix2d e;
	double d[3];
	bool success = cdm2d.eigensystem(e, d);
	assert(success); if(!success) return false;
	// 
	const DVector3d pu = surface->getDerivative(DEquation::deriv_du, pt2d);
	const DVector3d pv = surface->getDerivative(DEquation::deriv_dv, pt2d);
	const DVector2d pt_u = e.column(0);
	const DVector2d pt_v = e.column(1);
	const DVector3d e0 = DVector3d(pt_u.x * pu.x + pt_u.y * pv.x,
			pt_u.x * pu.y + pt_u.y * pv.y, pt_u.x * pu.z + pt_u.y * pv.z).normalized();
	const DVector3d e1 = DVector3d(pt_v.x * pu.x + pt_v.y * pv.x,
			pt_v.x * pu.y + pt_v.y * pv.y, pt_v.x * pu.z + pt_v.y * pv.z);
	// after surface transformation -> e0 and e1 are not necessarily orthogonal
	const DVector3d e2 = e0.crossProduct(e1).normalized();
	const DVector3d e1x = e2.crossProduct(e0).normalized(); 
	// rescale d[1] ??????
	d[2] = ControlSpace2dAdaptive::param_stretch_max_ratio * std::min(d[0], d[1]);
	cdm = ControlDataMatrix3d(e0, e1x, e2, d);

	if(cdm_e0) *cdm_e0 = e0;
	if(cdm_e1) *cdm_e1 = e1x;

	return true;
}

double diffKdValue(const ControlDataMatrix3d& cdm1, const ControlDataMatrix3d& cdm2) {
	return cdm1.countDifferenceRR(cdm2);
}

