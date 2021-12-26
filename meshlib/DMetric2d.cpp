// DMetric2d.cpp: implementation of the DMetric2d class.
//
//////////////////////////////////////////////////////////////////////

#include "DMetric2d.h"

#include "DMetric3d.h"
#include "DEquation.h"
//#include "MeshData.h"
#include "SurfaceParametric.h"
#include "SurfacePlane.h"
#include "ControlSpace2dAdaptive.h"
#include "ControlSpace3dAdaptive.h"
#include "Curve2dParametric.h"

const ControlDataMatrix2d ControlDataMatrix2d::IDENTITY(1.0, 0.0, 1.0);

DMetric2d::DMetric2d(const ControlDataMatrix2d& scdm)
{
	setData(scdm, ControlDataMatrix2d::IDENTITY);
}

/// Standard constructor (sizing matrix + parametric matrix)
DMetric2d::DMetric2d(const ControlDataMatrix2d& scdm, const ControlDataMatrix2d& _pcdm)
{
	setData(scdm, _pcdm);
}

DMetric2d::DMetric2d(std::shared_ptr<const SurfaceParametric> surface, const DPoint2d& pt)
	: DMetric2d(surface.get(), pt) { }

DMetric2d::DMetric2d(const SurfaceParametric* surface, const DPoint2d& pt)
{
	if(!surface){
		setToIdentity();
		return;
	}
	double g_ratio;
	setData(ControlDataMatrix2d::IDENTITY, surface->countParameterizationMatrix(pt, g_ratio));
}

DMetric2d::DMetric2d(const ControlDataMatrix2d& scdm,
	std::shared_ptr<const SurfaceParametric> surface, const DPoint2d& pt)
	: DMetric2d(scdm, surface.get(), pt) {}

DMetric2d::DMetric2d(const ControlDataMatrix2d& scdm, const SurfaceParametric* surface, const DPoint2d& pt)
{
	if(surface){
		double g_ratio;
		setData(scdm, surface->countParameterizationMatrix(pt, g_ratio));
	}else{
		setData(scdm, ControlDataMatrix2d::IDENTITY);
	}
}

DMetric2d::DMetric2d()
{
	m.setIdentity();
	mr.setIdentity();
	pcdm.setIdentity();
}

void DMetric2d::setToIdentity()
{
	m.setIdentity();
	mr.setIdentity();
	pcdm.setIdentity();
}

DPoint2d DMetric2d::transformMStoPS(const DMPoint2d &pt) const
{
	return m.multiplyMStoPS(pt);
}

DPoint2d DMetric2d::transformPStoRS(const DPoint2d &pt) const
{
	return pcdm*pt;
}

DPoint2d DMetric2d::transformRStoPS(const DPoint2d &pt) const
{
	return pcdm.inverse()*pt;
}

DMPoint2d DMetric2d::transformPStoMS(const DPoint2d &pt) const
{
	return mr.multiplyPStoMS(pt);
}

DVector2d DMetric2d::transformMStoPS(const DMVector2d &v) const
{
	return m.multiplyMStoPS(v);
}

DVector2d DMetric2d::transformPStoRS(const DVector2d &v) const
{
	return pcdm*v;
}

DVector2d DMetric2d::transformRStoPS(const DVector2d &v) const
{
	return pcdm.inverse()*v;
}

DMVector2d DMetric2d::transformPStoMS(const DVector2d &v) const
{
	return mr.multiplyPStoMS(v);
}

void DMetric2d::setData(const ControlDataMatrix2d& s_cdm, const ControlDataMatrix2d& p_cdm)
{
	pcdm = p_cdm;
	m = s_cdm * p_cdm.inverse();
	mr = m.inverse();
}

ControlDataMatrix2d::ControlDataMatrix2d(const DPoint2d& ptA, const DPoint2d& ptB, const DPoint2d& ptC)
{
	const DVector2d AB = ptB-ptA;
	const DVector2d AC = ptC-ptA;
	const DVector2d BC = ptC-ptB;
	DMatrix3d a(sqr(AB.x), sqr(AB.y), 2*AB.x*AB.y,
		sqr(AC.x), sqr(AC.y), 2*AC.x*AC.y, sqr(BC.x), sqr(BC.y), 2*BC.x*BC.y);

	double w = a.det();
	assert(abs(w) > SMALL_NUMBER);

	const DVector3d b(1.0, 1.0, 1.0);

	DMatrix3d a11 = a; 
	a11.setColumn(0, b);
	double w11 = a11.det();
	DMatrix3d a22 = a; 
	a22.setColumn(1, b);
	double w22 = a22.det();
	DMatrix3d a12 = a; 
	a12.setColumn(2, b);
	double w12 = a12.det();

	*this = DMetric2d::matrixSquareRoot(ControlDataMatrix2d(w11/w,w12/w,w22/w)).inverse();
}

ControlDataMatrix2d::ControlDataMatrix2d(const DVector2d& e0, const DVector2d& e1, double d[])
{
	DMatrix2d e;
	e.setColumn(0, e0);
	e.setColumn(1, e1);

	DMatrix2d et = e.transposed();
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			et.m[i][j] *= d[j];
	et = et*e;

	assert(abs(et.m[0][1]) - abs(et.m[1][0]) < SMALL_NUMBER);

	m11 = et.m[0][0];
	m22 = et.m[1][1];
	m12 = et.m[0][1];
}

ControlDataMatrix2d ControlDataMatrix2d::countMetric(const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2,
		DVector3d& e0, DVector3d& e1)
{
	return DMetric2d::matrixSquareRoot(countMetricTensor(pt0, pt1, pt2, e0, e1));
}

ControlDataMatrix2d ControlDataMatrix2d::countMetricTensor(const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2,
		DVector3d& e0, DVector3d& e1)
{
	e0 = (pt1 - pt0).normalized();
	// find orthogonal base ABxAE
	const DVector3d e2 = e0.crossProduct(pt2 - pt0);
	e1 = e2.crossProduct(e0).normalized();
	SurfacePlane plane(pt0, e0, e1);

	return countMetricTensor(DPoint2d(0.0, 0.0), plane.getParameters(pt1), plane.getParameters(pt2));
}

ControlDataMatrix2d ControlDataMatrix2d::countMetric(const DPoint2d& pt0, const DPoint2d& pt1, const DPoint2d& pt2)
{
	return DMetric2d::matrixSquareRoot(countMetricTensor(pt0, pt1, pt2));
}

ControlDataMatrix2d ControlDataMatrix2d::countMetricTensor(const DPoint2d& pt0, const DPoint2d& pt1, const DPoint2d& pt2)
{
	DVector2d vAB = pt1-pt0;
	DVector2d vAC = pt2-pt0;
	DVector2d vBC = pt2-pt1;

	double wmax = std::max(abs(vAB.x), abs(vAB.y));
	wmax = std::max(wmax, std::max(abs(vAC.x), abs(vAC.y)));
	wmax = std::max(wmax, std::max(abs(vBC.x), abs(vBC.y)));

	vAB /= wmax;
	vAC /= wmax;
	vBC /= wmax;

	DMatrix3d dm(	sqr(vAB.x), 2*vAB.x*vAB.y, sqr(vAB.y),
					sqr(vAC.x), 2*vAC.x*vAC.y, sqr(vAC.y),
					sqr(vBC.x), 2*vBC.x*vBC.y, sqr(vBC.y));
	
	double w = dm.det();
	assert(abs(w) > SMALL_NUMBER);
	DVector3d col_unit(1.0, 1.0, 1.0);
	// 11
	DVector3d col11 = dm.column(0);
	dm.setColumn(0, col_unit);
	double w11 = dm.det();
	dm.setColumn(0, col11);
	// 12
	DVector3d col12 = dm.column(1);
	dm.setColumn(1, col_unit);
	double w12 = dm.det();
	dm.setColumn(1, col12);
	// 22
	//DVector3d col22 = dm.column(2);
	dm.setColumn(2, col_unit);
	double w22 = dm.det();

	ControlDataMatrix2d cdmm(w11/(w*wmax*wmax), w12/(w*wmax*wmax), w22/(w*wmax*wmax));
	return cdmm.inverse();
}

ControlDataMatrix2d DMetric2d::matrixSquareRoot(const ControlDataMatrix2d& cdm)
{
	if(abs(cdm.m12) < mesh_data.relative_small_number)
		return ControlDataMatrix2d(sqrt(cdm.m11), cdm.m12, sqrt(cdm.m22));

	// -> stretch
	double delta = sqrt(sqr(cdm.m11 - cdm.m22) + 4.0*sqr(cdm.m12));
	double lx, ly;
//	if(abs(cdm.m11) > abs(cdm.m22)){
		lx = 0.5 * (cdm.m11 + cdm.m22 + delta);
		ly = 0.5 * (cdm.m11 + cdm.m22 - delta);
//	}else{
//		lx = 0.5 * (cdm.m11 + cdm.m22 - delta),
//		ly = 0.5 * (cdm.m11 + cdm.m22 + delta);
//	}

	double t_sincos = cdm.m12 / (lx-ly);
	double t_w = (ly-lx)*(ly+lx);
	double t_sin2 = (cdm.m11*ly-cdm.m22*lx)/t_w;
	double t_cos2 = (cdm.m22*ly-cdm.m11*lx)/t_w;
	// sqrt
	assert(lx >= 0.0);
	assert(ly >= 0.0);
	lx = sqrt(lx);
	ly = sqrt(ly);

	return ControlDataMatrix2d(
		lx * t_cos2 + ly * t_sin2,
		(lx - ly) * t_sincos,
		lx * t_sin2 + ly * t_cos2);
}

ControlDataMatrix2d DMetric2d::stretchToMatrix(const ControlDataStretch2d& stretch)
{
	double sin_dir = sin(stretch.angle);
	double cos_dir = cos(stretch.angle);
	double sin2 = sqr(sin_dir);
	double cos2 = sqr(cos_dir);
	double dxy = (stretch.lx - stretch.ly) * sin_dir * cos_dir;
	return ControlDataMatrix2d(stretch.lx * cos2 + stretch.ly * sin2, 
		dxy, stretch.lx * sin2 + stretch.ly * cos2);
}

ControlDataMatrix2d DMetric2d::stretchToMatrixWithAdjust(const ControlDataStretch2d& stretch)
{
	ControlDataStretch2d cds = stretch;
	adjustLengths(cds.lx, cds.ly);
	return stretchToMatrix(cds);
}

ControlDataStretch2d DMetric2d::matrixToStretch(const ControlDataMatrix2d& m)
{
	if(abs(m.m12) < mesh_data.relative_small_number)
		return ControlDataStretch2d(abs(m.m11), abs(m.m22), 0.0);
	double delta = sqrt(sqr(m.m11 - m.m22) + 4.0*sqr(m.m12));
	ControlDataStretch2d data;
	if(abs(m.m11) > abs(m.m22)){
		data.lx = 0.5 * (m.m11 + m.m22 + delta);
		data.ly = 0.5 * (m.m11 + m.m22 - delta);
	}else{
		data.lx = 0.5 * (m.m11 + m.m22 - delta),
		data.ly = 0.5 * (m.m11 + m.m22 + delta);
	}
	data.angle = atan((data.lx - m.m11) / m.m12);

	return data;
}

ostream& operator<<(ostream& os, const ControlDataStretch2d& data){
	string str_lx, str_ly, str_angle;
	DEquation::doubleToString(data.lx, DEquation::v_length, str_lx);
	DEquation::doubleToString(data.ly, DEquation::v_length, str_ly);
	DEquation::doubleToString(data.angle, DEquation::v_value, str_angle);
	os << "[" << str_lx << "," << str_ly << "," << str_angle << "]";
	return os;
}

istream& operator>>(istream& is, ControlDataStretch2d& data){
	string str_lx, str_ly, str_angle;
	char z;
	is >> ws >> z; // '['
	getline(is, str_lx, ',');
	getline(is, str_ly, ',');
	getline(is, str_angle, ']');
	DEquation eq;
	if(eq.parse(str_lx.c_str()))	data.lx = eq.getValue(0.0);
	if(eq.parse(str_ly.c_str()))	data.ly = eq.getValue(0.0);
	if(eq.parse(str_angle.c_str()))	data.angle = eq.getValue(0.0);
	return is;
}

double ControlDataMatrix2d::countDifferenceRR(const ControlDataMatrix2d& data) const {
	const ControlDataMatrix2d mr = data.inverse();
	const ControlDataMatrix2d r = this->inverse();

	// || (M11 * M22.inv - I) + (M11.inv * M22 - I) ||

	double rr11 = r.m11*data.m11 + r.m12*data.m12 + mr.m11*m11 + mr.m12*m12;
	double rr12 = r.m11*data.m12 + r.m12*data.m22 + mr.m11*m12 + mr.m12*m22;
	double rr21 = r.m12*data.m11 + r.m22*data.m12 + mr.m12*m11 + mr.m22*m12;
	double rr22 = r.m12*data.m12 + r.m22*data.m22 + mr.m12*m12 + mr.m22*m22;

	return sqrt(sqr(rr11-2.0) + sqr(rr12) + sqr(rr21) + sqr(rr22-2.0));
}

ControlDataMatrix2d ControlDataMatrix2d::tensorToTransformation() const
{
	DMatrix2d e;
	double d[2];
	bool success = eigensystem(e, d);
	assert(success);
	if(!success) return *this;
	assert(abs(1.0 - e.det()) < SMALL_NUMBER); // should be "rotation" matrix
	assert((d[0] > 0.0) && (d[1] > 0.0));
	double d0 = sqrt(d[0]);
	double d1 = sqrt(d[1]);
	// e.Identity(d).Transpose(e)
	// e[0][0] == e[1][1], e[0][1] == -e[1][0]
	double e0 = e.m[0][0];
	double e1 = e.m[1][0];
	return ControlDataMatrix2d(e0*e0*d0 + e1*e1*d1, (d0-d1)*e0*e1, e0*e0*d1 + e1*e1*d0);
}

bool ControlDataMatrix2d::setMinimum(const ControlDataMatrix2d& m) // using simultaneous reduction
{
	bool any_changes = false;

	if(abs(m12) < SMALL_NUMBER && abs(m.m12) < SMALL_NUMBER){
		if(m.m11 < m11){ any_changes = true; m11 = m.m11; }
		if(m.m22 < m22){ any_changes = true; m22 = m.m22; }
		return any_changes;
	}

	ControlDataMatrix2d mt11 = this->transformationToTensor();
	ControlDataMatrix2d mt22 = m.transformationToTensor();

	DMatrix2d n = mt11.inverse() * mt22;
	DMatrix2d e;
	double d[2];

	bool success = n.eigensystem(e, d);
	assert(success);
	if(!success) return false;

	for(int i = 0; i < 2; i++){
		d[i] = std::min(
			(e.m[0][i]*mt11.m11 + e.m[1][i]*mt11.m12)*e.m[0][i] + 
				(e.m[0][i]*mt11.m12 + e.m[1][i]*mt11.m22)*e.m[1][i],
			(e.m[0][i]*mt22.m11 + e.m[1][i]*mt22.m12)*e.m[0][i] + 
				(e.m[0][i]*mt22.m12 + e.m[1][i]*mt22.m22)*e.m[1][i]);
	}

	DMatrix2d p = e.inverse();

	DMatrix2d dm33 = p.transposed() * DMatrix2d(d[0], 0.0, 0.0, d[1]) * p;

	assert(abs(dm33.m[0][1]) - abs(dm33.m[1][0]) < SMALL_NUMBER);
	ControlDataMatrix2d mt33(dm33.m[0][0], dm33.m[0][1], dm33.m[1][1]);

	ControlDataMatrix2d m33 = mt33.tensorToTransformation();

	any_changes = this->countDifferenceRR(m33) > SMALL_NUMBER;
	if(any_changes)	*this = m33;

	return any_changes;
}

ControlDataExtMatrix2dRadial::ControlDataExtMatrix2dRadial(const DPoint2d& pt, 
		double r1, const ControlDataMatrix2d& d1,
		double r2, const ControlDataMatrix2d& d2, SurfaceConstPtr surface)
	: ControlDataExtMatrix2d(r1, d1, r2, d2), middle(pt), dmp(surface, pt)
{
	assert(r1+r2 > 0.0);
}

ControlDataExtMatrix2dRadial::ControlDataExtMatrix2dRadial(const DPoint2d& pt, 
		double r1, const ControlDataMatrix2d& d1, SurfaceConstPtr surface)
	: ControlDataExtMatrix2d(r1, d1), middle(pt), dmp(surface, pt)
{
	assert(r1 > 0.0);
}

bool ControlDataExtMatrix2dRadial::isPointWithin(const DPoint2d& pt) const 
{
	return dmp.transformPStoRS(middle-pt).length2() < sqr(radius1+radius2);
}

ControlDataMatrix2d ControlDataExtMatrix2dRadial::getControlDataMatrix(const DPoint2d& pt) const 
{
	double dist = dmp.transformPStoRS(middle-pt).length();
	if(dist <= radius1) 
		return data1;
	else if(dist <= radius1+radius2){
		double t = (dist-radius1)/radius2;
		return data1*(1.0-t) + data2*t;
	}else return data2;
}

ControlDataExtMatrix2dSegment::ControlDataExtMatrix2dSegment(
		const DPoint2d& _pt, const DVector2d& _dv, 
		double r1, const ControlDataMatrix2d& d1,
		double r2, const ControlDataMatrix2d& d2, SurfaceConstPtr surface)
	: ControlDataExtMatrix2d(r1, d1, r2, d2), pt0(_pt), dv(_dv)
{
	assert(r1+r2 > 0.0);
	if(surface){
		const DMetric2d dmp(surface, pt0+dv*0.5);
		const DVector2d dv_mp = dmp.transformPStoRS(dv);
		const DVector2d nv_mp(dv_mp.y, -dv_mp.x);
		nv = dmp.transformRStoPS(nv_mp * ((r1+r2)/nv_mp.length()));
	}else{
		nv = DVector2d(dv.y, -dv.x);
		nv *= (r1+r2)/nv.length();
	}
	w = dv.x * nv.y - dv.y * nv.x;
	assert(w != 0.0);
	dr = r1 / (r1+r2);
}

ControlDataExtMatrix2dSegment::ControlDataExtMatrix2dSegment(
		const DPoint2d& _pt, const DVector2d& _dv, 
		double r1, const ControlDataMatrix2d& d1, SurfaceConstPtr surface)
	: ControlDataExtMatrix2d(r1, d1), pt0(_pt), dv(_dv)
{
	assert(r1 > 0.0);
	if(surface){
		const DMetric2d dmp(surface, pt0+dv*0.5);
		const DVector2d dv_mp = dmp.transformPStoRS(dv);
		const DVector2d nv_mp(dv_mp.y, -dv_mp.x);
		nv = dmp.transformRStoPS(nv_mp * (r1/nv_mp.length()));
	}else{
		nv = DVector2d(dv.y, -dv.x);
		nv *= r1/nv.length();
	}
	w = dv.x * nv.y - dv.y * nv.x;
	assert(w != 0.0);
	dr = 1.0;
}

double ControlDataExtMatrix2dSegment::getLocalDistance(const DPoint2d& pt) const
{
	const DVector2d dp = pt-pt0;
	double u = (dp.x * nv.y - dp.y * nv.x) / w;
	if(u < 0.0 || u > 1.0) return LARGE_NUMBER; // orthogonal line doesn't cross the segment

	return (dv.x * dp.y - dv.y * dp.x) / w; // should be [-1,1] if withing the ControlExtSegment
}

bool ControlDataExtMatrix2dSegment::isPointWithin(const DPoint2d& pt) const 
{
	double v = getLocalDistance(pt);
	return (v >= -1.0 && v <= 1.0);
}

ControlDataMatrix2d ControlDataExtMatrix2dSegment::getControlDataMatrix(const DPoint2d& pt) const 
{
	double v = getLocalDistance(pt);

	if(v < dr) return data1;
	else if(v < 1.0){
		double t = (v-dr)/(1.0-dr);
		return data1*(1.0-t) + data2*t;
	}else return data2;
}

ControlDataExtMatrix2dCurve::ControlDataExtMatrix2dCurve(
		Curve2dConstPtr fig, double _t0, double _t1, 
		double r1, const ControlDataMatrix2d& d1,
		double r2, const ControlDataMatrix2d& d2, SurfaceConstPtr surface)
	: ControlDataExtMatrix2d(r1, d1, r2, d2)//, curve(fig), t0(_t0), t1(_t1)
{
	assert(fig);
	DataVector<double> nodes;
	fig->getPolyLineInRange(_t0, _t1, nodes);
	for(size_t i = 1; i < nodes.countInt(); i++){
		const DPoint2d pt0 = fig->getPoint(nodes[i-1]);
		const DPoint2d pt1 = fig->getPoint(nodes[i]);
		poly_control.add(std::make_shared<ControlDataExtMatrix2dSegment>(pt0, pt1-pt0, r1, d1, r2, d2, surface));
	}
}

ControlDataExtMatrix2dCurve::ControlDataExtMatrix2dCurve(
		Curve2dConstPtr fig, double _t0, double _t1, 
		double r1, const ControlDataMatrix2d& d1, SurfaceConstPtr surface)
	: ControlDataExtMatrix2d(r1, d1)//, curve(fig), t0(_t0), t1(_t1)
{
	assert(fig);
	DataVector<double> nodes;
	fig->getPolyLineInRange(_t0, _t1, nodes);
	for(size_t i = 1; i < nodes.countInt(); i++){
		const DPoint2d pt0 = fig->getPoint(nodes[i-1]);
		const DPoint2d pt1 = fig->getPoint(nodes[i]);
		poly_control.add(std::make_shared< ControlDataExtMatrix2dSegment>(pt0, pt1-pt0, r1, d1, surface));
	}
}

bool ControlDataExtMatrix2dCurve::isPointWithin(const DPoint2d& pt) const 
{
	for(size_t i = 0; i < poly_control.countInt(); i++)
		if(poly_control[i]->isPointWithin(pt)) return true;

	return false;
}

ControlDataMatrix2d ControlDataExtMatrix2dCurve::getControlDataMatrix(const DPoint2d& pt) const 
{
	double min_v = LARGE_NUMBER;
	size_t min_i = 0;
	for(size_t i = 0; i < poly_control.countInt(); i++){
		double v = abs(poly_control[i]->getLocalDistance(pt));
		if(v < min_v){
			min_v = v;
			min_i = i;
		}
	}

	if(min_v > 1.0) return data2;
	else return poly_control[min_i]->getControlDataMatrix(pt);
}

bool ControlDataMatrix2d::eigensystem(DMatrix2d& eigenvectors, double eigenval[]) const
{
	if(abs(m12) < SMALL_NUMBER){
		eigenvectors.m[0][0] = eigenvectors.m[1][1] = 1.0;
		eigenvectors.m[0][1] = eigenvectors.m[1][0] = 0.0;
		eigenval[0] = m11;
		eigenval[1] = m22;
		return true;
	}
	double delta = sqr(m11 - m22) + 4.0*m12*m12;
	if(delta < 0.0) return false;
	delta = sqrt(delta);
	eigenval[0] = 0.5 * (m11 + m22 + delta);
	eigenval[1] = 0.5 * (m11 + m22 - delta);

	double a11 = m11-eigenval[0];
	double a22 = m22-eigenval[0];
	// select "better" equation
	if(abs(a11) > abs(a22)){
		eigenvectors.m[1][0] = 1.0;
		eigenvectors.m[0][0] = -m12/a11;
	}else{
		eigenvectors.m[1][0] = 1.0;
		eigenvectors.m[0][0] = -a22/m12;
	}
	// normalize
	double len = sqrt(sqr(eigenvectors.m[0][0]) + sqr(eigenvectors.m[1][0]));
	if(eigenvectors.m[0][0] < 0.0) len = -len;
	eigenvectors.m[0][0] /= len;
	eigenvectors.m[1][0] /= len;

	// second eigenvector - must be orthogonal
	eigenvectors.m[0][1] = -eigenvectors.m[1][0];
	eigenvectors.m[1][1] =  eigenvectors.m[0][0];

	return true;
}

bool ControlDataMatrix2d::eigenvalues(double eigenval[]) const
{
	if(abs(m12) < SMALL_NUMBER){
		eigenval[0] = m11;
		eigenval[1] = m22;
		return true;
	}
	double delta = sqr(m11 - m22) + 4.0*m12*m12;
	if(delta < 0.0) return false;
	delta = sqrt(delta);
	eigenval[0] = 0.5 * (m11 + m22 + delta);
	eigenval[1] = 0.5 * (m11 + m22 - delta);
	return true;
}

double ControlDataMatrix2d::minEigenvalue() const
{
	if(abs(m12) < SMALL_NUMBER)
		return std::min(m11, m22);
	double delta = sqr(m11 - m22) + 4.0*m12*m12;
	assert(delta >= 0.0);
	if(delta < 0.0) delta = 0.0;
	return 0.5 * (m11 + m22 - sqrt(delta));
}

void DMetric2d::adjustLengths(double & lx, double & ly)
{
	double p = mesh_data.getModelDiameter();
	double max_len = p * std::min(ControlSpace2dAdaptive::param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
	double min_len = std::max(p * ControlSpace2dAdaptive::param_min_diameter_ratio, ControlSpace2dAdaptive::param_min_length);

	if(lx > max_len) lx = max_len; 
	else if(lx < min_len) lx = min_len; 

	if(ly > max_len) ly = max_len; 
	else if(ly < min_len) ly = min_len; 

	if(ly > lx){
		if((ly / lx) > ControlSpace2dAdaptive::param_stretch_max_ratio)
			ly = lx * ControlSpace2dAdaptive::param_stretch_max_ratio;
	}else{
		if((lx / ly) > ControlSpace2dAdaptive::param_stretch_max_ratio)
			lx = ly * ControlSpace2dAdaptive::param_stretch_max_ratio;
	}	
}

double diffKdValue(const ControlDataMatrix2d& cdm1, const ControlDataMatrix2d& cdm2) {
	return cdm1.countDifferenceRR(cdm2);
}
