// ControlSpace3dAdaptive.cpp: implementation of the ControlSpace3dAdaptive class.
//
//////////////////////////////////////////////////////////////////////

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "ControlSpace3dAdaptive.h"
#include "ControlSpace2dAdaptive.h"
#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshFace.h"
#include "SurfacePlane.h"
#include "Metric3dContext.h"
#include "DataHashTable.h"
#include "Curve3dParametric.h"
#include "DPlane.h"

double ControlSpace3dAdaptive::param_max_diameter_ratio;

/*
void ControlSpace3dAdaptive::addControlSegment(
	const DPoint3d &pt1, const DPoint3d &pt2, 
	const ControlDataStretch3d& data1, double r1, 
	const ControlDataStretch3d& data2, double r2)
{
	double seg_len = pt1.distance(pt1);
	if(seg_len < mesh_data.relative_small_number) return;

	const DPoint3d dv = pt2-pt1;
	double dt = std::min(0.01, 2.0*std::min(std::min(data1.lx, data1.ly), data1.lz) / seg_len);

	addSegmentControlPoint(pt1, dv, dt, data1, r1);
	if(r2 > 0.0){
		dt = std::min(0.01, 2.0*std::min(std::min(data2.lx, data2.ly), data2.lz) / seg_len);
		addSegmentControlPoint(pt1, dv, dt, data2, r1+r2);
	}

// TODO
//	if(r1+r2 > 0.0)
//		m_ext_source_points.add(new ControlDataExtMatrixSegment3d(pt1, dv,
//			r1, DMetric3d::stretchToMatrix(data1), 
//			r2, DMetric3d::stretchToMatrix(data2)));
}
*/

/*
void ControlSpace3dAdaptive::addControlSegment(
	const Curve2dParametric *curve, double t0, double t1, 
	const ControlDataStretch2d& data1, double r1, 
	const ControlDataStretch2d& data2, double r2)
{
	DataVector<double> polyline(50);
	curve->getPolyLine(t0, t1, polyline);

	int ct = polyline.countInt();
	assert(ct > 1);
	DPoint2d pt0 = curve->getPoint(polyline.get(0));
	for(int i = 1; i < polyline.countInt(); i++){
		DPoint2d pt1 = curve->getPoint(polyline.get(i));
		addControlSegment(pt0, pt1, data1, r1, data2, r2);
		pt0 = pt1;
	}
}
*/

double ControlSpace3dAdaptive::getMaxMetricLen(double ratio) const
{
	return m_box.getDiameter() * ((ratio > 0) ? ratio : param_max_diameter_ratio);
}

void ControlSpace3dAdaptive::setMaxMetric(double ratio)
{
	double p = m_box.getDiameter();
	double max_len = p * ((ratio > 0) ? ratio : param_max_diameter_ratio);

	ControlDataMatrix3d max_data(max_len);
	if (getControlNodesCount() > 0) {
		forEachControlNode([&](ControlNode3d& node) {
			node.control_data = max_data;
			node.w = -1.0;
		});
		m_initialized = 1;
	}
	else
		setGlobalMetric(max_data);
}

void ControlSpace3dAdaptive::setGlobalMetric(const ControlDataMatrix3d& cdm)
{
	assert(getControlNodesCount() > 0);

	forEachControlNode([&](ControlNode3d& node) {
		node.control_data = cdm;
		node.w = -1.0;
	});
	m_initialized = 1;
}

// no adaptation, just set min value
int ControlSpace3dAdaptive::setMinControlValue(const ControlDataExtMatrix3d& data)
{
	// 1. check all grid vertices (***)
	int completely_new_count = 0;
	forEachControlNode([&](ControlNode3d& qv) {
		if (data.isPointWithin(qv.coord)) {
			if (qv.w < 0.0)	// node already initialized
				qv.control_data.setMinimum(data.getControlDataMatrix3d(qv.coord));
			else {
				qv.control_data = data.getControlDataMatrix3d(qv.coord);
				qv.w = -1.0;
				++completely_new_count;
			}
		}
	});

	return completely_new_count;
}

bool ControlSpace3dAdaptive::applyAsMinimum(CS3dPtr space)
{
	assert(m_initialized > 0 && m_initialized < 10);

	bool any_changes = false;
	// 1: adjust structure of this space with other space
	space->forEachControlNode([&](ControlNode3d& qv) {
		any_changes |= setMinControl(qv.coord, qv.control_data, false);
	});
	// 2: calculate minimum metric for own nodes
	forEachControlNode([&](ControlNode3d& qv) {
		any_changes |= qv.control_data.setMinimum(space->getMetricAtPoint(qv.coord));
	});

	return any_changes;
}

bool ControlSpace3dAdaptive::applyAsMinimum(CS2dPtr space)
{
	assert(m_initialized > 0 && m_initialized < 10);
	SurfaceConstPtr surface = space->getBaseSurface();
	ControlDataMatrix3d cdm;
	bool any_changes = false;
	// 1: adjust structure of this space with other space
	space->forEachControlNode([&](ControlNode2d& qv) {
		if (DMetric3d::projectCDMto3D(qv.control_data, surface, qv.coord, cdm))
			any_changes |= setMinControl(surface->getPoint(qv.coord), cdm, false);
	});

	return any_changes;
}

double ControlSpace3dAdaptive::compareMetricWith(CS3dPtr space)
{
	assert(m_initialized > 0);
	assert(!space->isAdaptive() || space->getAsAdaptive()->initializationState() > 0);
	double max_diff = 0.0;
	// 1: check metric for each quad-vertex of this space with other space
	forEachControlNode([&](const ControlNode3d & qv) {
		double diff = qv.control_data.countDifferenceRR(space->getMetricAtPoint(qv.coord));
		if (diff > max_diff) max_diff = diff;
	});
	// 2: ... and the other way round
	space->forEachControlNode([&](const ControlNode3d & qv) {
		double diff = qv.control_data.countDifferenceRR(getMetricAtPoint(qv.coord));
		if (diff > max_diff) max_diff = diff;
	});

	return max_diff;
}

/*
bool ControlSpace3dAdaptive::setMinControl(const ControlDataExtMatrix2dRadial& data)
{
	assert(m_initialized > 0 && m_initialized < 10);

	double r = data.totalRadius();
	if(r <= 0.0){
		return setMinControl(data.getMiddle(), data.getInnerData());
	}

#ifdef STORE_CONTROL_EPS
	storeEPS("control-min-control-before");
#endif

	bool any_changes = false;

	DMetric2d dm(base_surface, data.getMiddle());
	if(data.getRadius2() > 0.0){
		const ControlDataStretch2d cds = DMetric2d::matrixToStretch(data.getOuterData());
		double min_len = 2*std::min(cds.lx, cds.ly);
		const double da = 2*PI/param_radial_parts;
		double a = 0.0;
		for(int i = 0; i < param_radial_parts; i++, a+=da)
			any_changes |= setRadialControlPoint(dm, data.getMiddle(), data.totalRadius(), 
								data.getOuterData(), a, a+da, min_len);
	}

#ifdef STORE_CONTROL_EPS
//	storeEPS("control-min-control-middle");
#endif

	if(data.getRadius1() > 0.0){
		const ControlDataStretch2d cds = DMetric2d::matrixToStretch(data.getInnerData());
		double min_len = 2*std::min(cds.lx, cds.ly);
		const double da = 2*PI/param_radial_parts;
		double a = 0.0;
		for(int i = 0; i < param_radial_parts; i++, a+=da)
			any_changes |= setRadialControlPoint(dm, data.getMiddle(), data.getRadius1(), 
								data.getInnerData(), a, a+da, min_len);
	}

	any_changes |= (setMinControlValue(data) > 0);

#ifdef STORE_CONTROL_EPS
//	storeEPS("control-min-control-after");
#endif
	return any_changes;
}
*/
/*
bool ControlSpace2dAdaptive::setMinControl(const ControlDataExtMatrix2dSegment& data)
{
	assert(m_initialized > 0 && m_initialized < 10);

	DMetric2d dmp(base_surface, data.getMiddle());
	double seg_len = dmp.transform(data.getDv()).length();
	if(seg_len < mesh_data.relative_small_number) return false;

	const ControlDataStretch2d cds1 = DMetric2d::matrixToStretch(data.getInnerData());
	double dt1 = std::min(0.1, 2.0*std::min(cds1.lx, cds1.ly) / seg_len);


	double r = data.totalRadius();
	if(r <= 0.0){
		return setSegmentMinControl(data.getStart(), data.getDv(), dt1, 
			data.getInnerData(), 0.0);
	}

#ifdef STORE_CONTROL_EPS
//	storeEPS("control-min-control-before");
#endif

	bool any_changes = false;

	if(data.getRadius2() > 0.0){
		const ControlDataStretch2d cds2 = DMetric2d::matrixToStretch(data.getOuterData());
		double dt2 = std::min(0.1, 2.0*std::min(cds2.lx, cds2.ly) / seg_len);
		any_changes |= setSegmentMinControl(data.getStart(), data.getDv(), dt2, 
							data.getOuterData(), data.totalRadius());
	}

#ifdef STORE_CONTROL_EPS
//	storeEPS("control-min-control-middle");
#endif

	if(data.getRadius1() > 0.0){
		any_changes |= setSegmentMinControl(data.getStart(), data.getDv(), dt1, 
							data.getInnerData(), data.getRadius1());
	}

	any_changes |= (setMinControlValue(data) > 0);

#ifdef STORE_CONTROL_EPS
//	storeEPS("control-min-control-after");
#endif
	return any_changes;
}
*/

/*
bool ControlSpace2dAdaptive::setRadialControlPoint(const DMetric2d& dmp, 
	const DPoint2d& pt, double r, const ControlDataMatrix2d& cdm, 
	double a0, double a1, double min_len)
{
	assert(m_initialized > 0 && m_initialized < 10);
	DPoint2d v0(1.0, 0.0);
	v0.turn(sin(a0), cos(a0));
	DPoint2d v1(1.0, 0.0);
	v1.turn(sin(a1), cos(a1));
	v0 *= (r/dmp.transform(v0).length());
	v1 *= (r/dmp.transform(v1).length());
	double dist = dmp.transform(v0-v1).length2();
	if(dist > min_len){
		bool any_changes = setRadialControlPoint(dmp, pt, r, cdm, a0, 0.5*(a0+a1), min_len);
		any_changes |= setRadialControlPoint(dmp, pt, r, cdm, 0.5*(a0+a1), a1, min_len);
		return any_changes;
	}else{
		return setMinControl(pt+v0, cdm, false);
	}
}
*/

/*
void ControlSpace2dAdaptive::addRadialControlPoint(const DMetric2d& dmp, 
	const DPoint2d& pt, double r, const ControlDataMatrix2d& cdm, 
	double a0, double a1, double min_len)
{
	assert(m_initialized == 0);
	DPoint2d v0(1.0, 0.0);
	v0.turn(sin(a0), cos(a0));
	DPoint2d v1(1.0, 0.0);
	v1.turn(sin(a1), cos(a1));
	v0 *= (r/dmp.transform(v0).length());
	v1 *= (r/dmp.transform(v1).length());
	double dist = dmp.transform(v0-v1).length2();
	if(dist > min_len){
		addRadialControlPoint(dmp, pt, r, cdm, a0, 0.5*(a0+a1), min_len);
		addRadialControlPoint(dmp, pt, r, cdm, 0.5*(a0+a1), a1, min_len);
	}else{
		addControlNode(ControlNode2d(pt+v0, cdm));
	}
}
*/

/*
void ControlSpace3dAdaptive::addSegmentControlPoint(
	const DPoint3d& pt0, const DPoint3d& dv, double dt,
	const ControlDataStretch3d& data, double r)
{
	assert(m_initialized == 0);
	const ControlDataMatrix3d cdm = DMetric3d::stretchToMatrix(data);

	if(r > 0.0){
		// TODO -> radial
		const DPoint2d nv(dv_mp.y, -dv_mp.x);
		const DPoint2d nv = dmp.transformReverse(nv_mp * (r/nv_mp.length()));

		for(double t = 0.0; t <= 1.0; t += dt){
			const DPoint2d pt = pt0 + dv * t;
			addControlNode(ControlNode2d(pt+nv, cdm));
			addControlNode(ControlNode2d(pt-nv, cdm));
		}
	}else{
		for(double t = 0.0; t <= 1.0; t += dt)
			addControlNode(ControlNode3d(pt0 + dv * t, cdm));
	}
}
*/

/*
bool ControlSpace3dAdaptive::setSegmentMinControl(
	const DPoint3d& pt0, const DPoint3d& dv, double dt,
	const ControlDataMatrix3d& cdm, double r)
{
	assert(m_initialized > 0 && m_initialized < 10);
	if(r > 0.0){
		// TODO -> radial
		const DPoint2d dv_mp = dmp.transform(dv);
		const DPoint2d nv_mp(dv_mp.y, -dv_mp.x);
		double rn = r/nv_mp.length();
		const DPoint2d nv = dmp.transformReverse(nv_mp * rn);
		const DPoint2d nv_r = dmp.transformReverse(nv_mp * -rn);

		bool any_changes = false;
		for(double t = 0.0; t <= 1.0; t += dt){
			const DPoint2d pt = pt0 + dv * t;
			any_changes |= setMinControl(pt+nv, cdm, false);
			any_changes |= setMinControl(pt+nv_r, cdm, false);
		}
		return any_changes;
	}else{
		bool any_changes = false;
		for(double t = 0.0; t <= 1.0; t += dt)
			any_changes |= setMinControl(pt0 + dv * t, cdm, false);
		return any_changes;
	}
}
*/

/*
void ControlSpace3dAdaptive::addControlPoint(const DPoint3d& pt, 
	const ControlDataStretch3d& data1, double r1, 
	const ControlDataStretch3d& data2, double r2)
{
	assert(m_initialized == 0);
	const ControlNode3d qv(pt, DMetric3d::stretchToMatrix(data1));

	double r = r1 + r2;
	if(r <= 0.0){
		addControlNode(qv);
		return;
	}

	if(r1 > 0.0){
		double min_len = 2*std::min(std::min(data1.lx, data1.ly), data1.lz);
		const ControlDataMatrix3d cdm = DMetric3d::stretchToMatrix(data1);
		const double da = 2*PI/param_radial_parts;
		double a = 0.0;
		// TODO -> spherical
		for(int i = 0; i < param_radial_parts; i++, a+=da)
			addRadialControlPoint(dm, pt, r1, cdm, a, a+da, min_len);
	}
	
	if(r2 > 0.0){
		double min_len = 2*std::min(std::min(data2.lx, data2.ly), data2.lz);
		const ControlDataMatrix3d cdm = DMetric3d::stretchToMatrix(data2);
		const double da = 2*PI/param_radial_parts;
		double a = 0.0;
		// TODO -> spherical
		for(int i = 0; i < param_radial_parts; i++, a+=da)
			addRadialControlPoint(dm, pt, r1+r2, cdm, a, a+da, min_len);
	}
	
//	m_ext_source_points.add(new ControlDataExtMatrixSphere(pt, r1, 
//		qv.control_data, r2, DMetric3d::stretchToMatrix(data2)));
}
*/

void ControlSpace3dAdaptive::addSimpleBoundaryControlData(
	DataVector<MeshFace*> & bfaces, DataVector<MeshPoint3d*> & bpoints)
{
	assert(m_initialized > 0 && m_initialized < 10);
	double stretch_max_ratio = std::min(
		ControlSpace2dAdaptive::param_contour_stretch_max_ratio,
		ControlSpace2dAdaptive::param_stretch_max_ratio);

	int pct = (int)bpoints.countInt();
	DataHashTableKeyValue<MeshPoint3d*, int> bindex(2*pct, nullptr);
	for(int i = 0; i < pct; i++)
		bindex.insert(bpoints[i], i);

	DataVector<ControlDataMatrix3d> cdm_ave(pct, ControlDataMatrix3d());
	DataVector<int> cdm_counter(pct, 0);

	int fct = (int)bfaces.countInt();
	for(int i = 0; i < fct; i++){
		MeshFace* face = bfaces[i];
		if(face->getEdgeCount() != 3) continue;
		const DPoint3d& pt0 = face->getPoint(0)->getCoordinates();
		const DPoint3d& pt1 = face->getPoint(1)->getCoordinates();
		const DPoint3d& pt2 = face->getPoint(2)->getCoordinates();
		//  - 2D metric from simplex (three points)
		const DVector3d vAB = pt1 - pt0;
		// find orthogonal base ABxAE
		const DVector3d vAD = vAB.crossProduct(pt2 - pt0);
		const DVector3d vAE = vAD.crossProduct(vAB);
		SurfacePlane plane(pt0, vAB.normalized(), vAE.normalized());
		const DPoint2d pt0x(0.0, 0.0); // = plane.getParameters(pt0)
		ControlDataMatrix2d cdm2d = ControlDataMatrix2d::countMetric(pt0x, 
			plane.getParameters(pt1), plane.getParameters(pt2));
		//  - plus orthogonal third direction and 3D metric
		//  - construct 3D metric from 2D metric + third direction
		DMatrix2d e;
		double d[3];
		bool result = cdm2d.eigensystem(e, d);
		assert(result);
		if(!result) return;
		const DVector3d v1 = plane.getPoint(pt0x+e.column(0)) - pt0;
		const DVector3d v2 = plane.getPoint(pt0x+e.column(1)) - pt0;
		d[2] = stretch_max_ratio * std::min(d[0], d[1]);
		// eigensystem -> 3d metric
		ControlDataMatrix3d cdm(v1, v2, v1.crossProduct(v2).normalized(), d);

		for(int j = 0; j < 3; j++){
			int ind = bindex.getValue(face->getPoint(j), -1);
			assert(ind >= 0);
			cdm_ave[ind] += cdm;
			cdm_counter[ind]++;
		}
	}

	for(int i = 0; i < pct; i++){
		if(cdm_counter[i] > 0)
			setMinControl(bpoints[i]->getCoordinates(), cdm_ave[i] / cdm_counter[i]);
	}
}

bool ControlSpace3dAdaptive::addSimpleBoundaryControlData(MeshContainer3dSurface * surface_mesh)
{
	assert(m_initialized > 0 && m_initialized < 10);
	double stretch_max_ratio = std::min(
		ControlSpace2dAdaptive::param_contour_stretch_max_ratio,
		ControlSpace2dAdaptive::param_stretch_max_ratio);

	int pct = surface_mesh->getPointsCount();

	DataVector<ControlDataMatrix3d> cdm_ave(pct, ControlDataMatrix3d());
	DataVector<int> cdm_counter(pct, 0);

	int fct = surface_mesh->getFacesCount();
	for(int i = 0; i < fct; i++){
		MeshFace* face = surface_mesh->getFaceAt(i);
		if(face->getEdgeCount() != 3) continue;
		const DPoint3d& pt0 = face->getPoint(0)->getCoordinates();
		const DPoint3d& pt1 = face->getPoint(1)->getCoordinates();
		const DPoint3d& pt2 = face->getPoint(2)->getCoordinates();
		//  - 2D metric from simplex (three points)
		const DVector3d vAB = pt1 - pt0;
		// find orthogonal base ABxAE
		const DVector3d vAD = vAB.crossProduct(pt2 - pt0);
		const DVector3d vAE = vAD.crossProduct(vAB);
		SurfacePlane plane(pt0, vAB.normalized(), vAE.normalized());
		const DPoint2d pt0x(0.0, 0.0); // = plane.getParameters(pt0)
		ControlDataMatrix2d cdm2d = ControlDataMatrix2d::countMetric(pt0x, 
			plane.getParameters(pt1), plane.getParameters(pt2));
		//  - plus orthogonal third direction and 3D metric
		//  - construct 3D metric from 2D metric + third direction
		DMatrix2d e;
		double d[3];
		bool result = cdm2d.eigensystem(e, d);
		assert(result);
		if(!result) continue;
		const DVector3d v1 = plane.getPoint(pt0x+e.column(0)) - pt0;
		const DVector3d v2 = plane.getPoint(pt0x+e.column(1)) - pt0;
		d[2] = stretch_max_ratio * std::min(d[0], d[1]);
		// eigensystem -> 3d metric
		ControlDataMatrix3d cdm(v1, v2, v1.crossProduct(v2).normalized(), d);

		for(int j = 0; j < 3; j++){
			int ind = face->getPoint(j)->getIndex();
			cdm_ave[ind] += cdm;
			cdm_counter[ind]++;
		}
	}

	bool any_change = false;
	for(int i = 0; i < pct; i++){
		if(cdm_counter[i] > 0)
			any_change |= setMinControl(surface_mesh->getPointAt(i)->getCoordinates(), cdm_ave[i] / cdm_counter[i]);
	}

	return any_change;
}

bool ControlSpace3dAdaptive::addSimpleMeshControlData(MeshContainer3d * mesh)
{
	assert(m_initialized > 0 && m_initialized < 10);

	int pct = mesh->getPointsCount();

	DataVector<ControlDataMatrix3d> cdm_ave(pct, ControlDataMatrix3d());
	DataVector<int> cdm_counter(pct, 0);

	int bct = mesh->getBlocksCount();
	for(int i = 0; i < bct; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		int bpct = block->getPointCount();
		if(bpct != 4) continue;
		ControlDataMatrix3d cdm_simplex = ControlDataMatrix3d::countMetric(
			block->getPoint(0)->getCoordinates(), 
			block->getPoint(1)->getCoordinates(), 
			block->getPoint(2)->getCoordinates(), 
			block->getPoint(3)->getCoordinates());

		double d[3];
		DMatrix3d e;
		if(!cdm_simplex.eigensystem(e, d)) continue;
		double min_d = std::min( std::min(d[0], d[1]), d[2]);
		double max_d = std::max( std::max(d[0], d[1]), d[2]);
		double anisotropy_ratio = max_d / min_d;
		if(anisotropy_ratio > ControlSpace2dAdaptive::param_stretch_max_ratio){
			double min_dx = max_d / ControlSpace2dAdaptive::param_stretch_max_ratio;
			for(int j = 0; j < 3; j++)
				d[j] = std::max(d[j], min_dx);
			cdm_simplex.setEigensystem(e, d);
		}

		for(int j = 0; j < bpct; j++){
			int ind = block->getPoint(j)->getIndex();
			cdm_ave[ind] += cdm_simplex;
			cdm_counter[ind]++;
		}
	}

	bool any_change = false;
	for(int i = 0; i < pct; i++){
		if(cdm_counter[i] > 0)
			any_change |= setMinControl(mesh->getPointAt(i)->getCoordinates(), cdm_ave[i] / cdm_counter[i]);
	}

	return any_change;
}

bool ControlSpace3dAdaptive::updateForBoundarySegment(const MeshEdge3d* edge)
{
	double len = 1.2*edge->getLength();
	if(len < ControlSpace2dAdaptive::param_min_length) 
		len = ControlSpace2dAdaptive::param_min_length;

	const DPoint3d& pt0 = edge->getMeshPoint(0)->getCoordinates();
	const DPoint3d& pt1 = edge->getMeshPoint(1)->getCoordinates();

	const DVector3d v01 = pt1-pt0;
	const DVector3d e0 = v01.normalized();
	DVector3d e1, e2;
	e0.orthonormalVectors(e1, e2); // orthonormal
	double d[3] = {
		len, 
		len * ControlSpace2dAdaptive::param_stretch_max_ratio, 
		len * ControlSpace2dAdaptive::param_stretch_max_ratio};
	ControlDataMatrix3d cdm(e0, e1, e2, d);	// from eigensystem

	bool result = setMinControl(pt0, cdm);
	result |= setMinControl(pt1, cdm);
	result |= setMinControl(DPoint3d(pt0, pt1, 0.5), cdm);
	return result;

//	ControlDataExtMatrix3dSegment data(pt0, v01, d[0], cdm);
//	ControlDataExtMatrix3dSegment data(pt0, v01, 0.0, cdm);
//	return setMinControl(data);
}

void ControlSpace3dAdaptive::addSimpleBoundaryControlDataIsotroPIc(
	DataVector<MeshFace*> & /* bfaces */, DataVector<MeshPoint3d*> & /* bpoints */)
{
	assert(false);
/*
	assert(m_initialized > 0 && m_initialized < 10);
	// boundary_mesh -- 3D mesh containing boundary nodes
	int bct = boundary_mesh->getBlocksCount();
	for(int i = 0; i < bct; i++){
		MeshBlock* block = boundary_mesh->getBlockAt(i);
		int fct = block->getFaceCount();
		for(int j = 0; j < fct; j++){
			MeshFace* face = block->getFace(j);
			if(face->getEdgeCount() != 3) continue;
			const DPoint3d& pt0 = face->getPoint(0)->getCoordinates();
			const DPoint3d& pt1 = face->getPoint(1)->getCoordinates();
			const DPoint3d& pt2 = face->getPoint(2)->getCoordinates();
			double len = (pt0.distance(pt1) + pt1.distance(pt2) + pt2.distance(pt0))*(1.0/3.0);
			setMinControl(DPoint3d::average(pt0, pt1, pt2), ControlDataMatrix3d(len, len, len, 0.0, 0.0, 0.0));
			//addControlNode(node);
		}
	}
*/
}

int ControlSpace3dAdaptive::smoothenMetricForNodes(ControlNode3d * cn1, ControlNode3d * cn2, const DVector3d& dv)
{
	int result = 0;
	double dr = ControlSpace2dAdaptive::param_gradation_ratio;
	assert(dr >= 1.0);
	if(dr < 1.0) dr = 1.0;
	double a_max = 2*(dr-1)/(dr+1);

	double diff = cn1->control_data.countDifferenceRR(cn2->control_data);
	if(diff < 0.01) return 0;

	double dlen = dv.length();
	//assert(dlen > VERY_SMALL_NUMBER);
	if (dlen < VERY_SMALL_NUMBER) {
		if (cn1->control_data.setMinimum(cn2->control_data))
			result += 1;
		if (cn2->control_data.setMinimum(cn1->control_data))
			result += 2;
		return result;
	}
	const DVector3d dvn = dv * (1.0/dlen);
	// count required length of element along the given line 
	//	- according to metric at each node
	double len_m[2] = {(cn1->control_data * dvn).length(), (cn2->control_data * dvn).length()};
	int vmax = (len_m[1]>len_m[0])?1:0;
	int vmin = 1-vmax;
	// count max ratio:
	DMatrix3d e1, e2;
	double d1[3], d2[3];
	cn1->control_data.eigensystem(e1, d1);
	cn2->control_data.eigensystem(e2, d2);
	const ControlDataMatrix3d m1inv = cn1->control_data.inverse();
	const ControlDataMatrix3d m2inv = cn2->control_data.inverse();
	double real_max_ratio_1 = std::max( std::max((m2inv * (e1.column(0) * d1[0])).length(), 
		(m2inv * (e1.column(1) * d1[1])).length()), (m2inv * (e1.column(2) * d1[2])).length() );
	double real_max_ratio_2 = std::max( std::max((m1inv * (e2.column(0) * d2[0])).length(), 
		(m1inv * (e2.column(1) * d2[1])).length()), (m1inv * (e2.column(2) * d2[2])).length() );
	// check ratio of metric-lengths
	double a_m = (len_m[vmax]-len_m[vmin])/dlen;
	if(a_m > a_max)
		len_m[vmax] = (a_m=a_max) * dlen + len_m[vmin];
	// count max diff
	double a1 = len_m[vmin]/(1.0-0.5*a_m);
	double max_diff_ratio = 1.0-(1.0-dr)*dlen/a1;
	double real_dr[2] = { dr, dr };

	if(real_max_ratio_1 > max_diff_ratio){
		if(cn1->control_data.setMinimum(cn2->control_data * max_diff_ratio)) 
			result += 1;
	}else real_dr[0] = 1.0 - a1*(1.0-real_max_ratio_1)/dlen;

	if(real_max_ratio_2 > max_diff_ratio){
		if(cn2->control_data.setMinimum(cn1->control_data * max_diff_ratio)) 
			result += 2;
	}else real_dr[1] = 1.0 - a1*(1.0-real_max_ratio_2)/dlen;

	if(cn1->gradationUnknown() || real_dr[0] > cn1->max_gradation_ratio) 
		cn1->max_gradation_ratio = real_dr[0];
	if(cn2->gradationUnknown() || real_dr[1] > cn2->max_gradation_ratio) 
		cn2->max_gradation_ratio = real_dr[1];

	return result;
}

bool ControlSpace3dAdaptive::adjustMetricForNodesLeft(ControlDataMatrix3d& cdm0,
	const ControlDataMatrix3d& cdm1, const DVector3d& dv)
{
	double lm1 = (cdm1.inverse() * dv).length();
	double s1 = (ControlSpace2dAdaptive::param_gradation_ratio - 1.0) * lm1 + 1.0;
	return cdm0.setMinimum(cdm1 * s1);
}

bool ControlSpace3dAdaptive::adjustMetricForNodesLeft(ControlDataMatrix3d& cdm0,
	const ControlDataMatrix3d& cdm1, const DVector3d& dv,
	double & h0_max, double h1_min)
{
	double lm1 = (cdm1.inverse() * dv).length();
	double s1 = (ControlSpace2dAdaptive::param_gradation_ratio - 1.0) * lm1 + 1.0;
	return (h1_min * s1 < h0_max) && cdm0.setMinimum(cdm1 * s1, &h0_max);
}

bool ControlSpace3dAdaptive::smoothenMetricForNodesLeft(ControlNode3d * cn1, 
	const ControlNode3d * cn2, const DVector3d& dv)
{
	double dr = ControlSpace2dAdaptive::param_gradation_ratio;
	assert(dr >= 1.0);
	if (dr < 1.0) dr = 1.0;
	double a_max = 2 * (dr - 1) / (dr + 1);

	double diff = cn1->control_data.countDifferenceRR(cn2->control_data);
	if (diff < 0.01) return false;

	double dlen = dv.length();
	//assert(dlen > VERY_SMALL_NUMBER);
	if (dlen < VERY_SMALL_NUMBER)
		return cn1->control_data.setMinimum(cn2->control_data);

	const DVector3d dvn = dv * (1.0 / dlen);
	// count required length of element along the given line 
	//	- according to metric at each node
	double len_m[2] = { (cn1->control_data * dvn).length(), (cn2->control_data * dvn).length() };
	int vmax = (len_m[1]>len_m[0]) ? 1 : 0;
	int vmin = 1 - vmax;
	// count max ratio:
	DMatrix3d e1;
	double d1[3];
	cn1->control_data.eigensystem(e1, d1);
	const ControlDataMatrix3d m2inv = cn2->control_data.inverse();
	double real_max_ratio_1 = std::max(std::max(
		(m2inv * (e1.column(0) * d1[0])).length(),
		(m2inv * (e1.column(1) * d1[1])).length()), 
		(m2inv * (e1.column(2) * d1[2])).length());

	// check ratio of metric-lengths
	double a_m = (len_m[vmax] - len_m[vmin]) / dlen;
	if (a_m > a_max)
		len_m[vmax] = (a_m = a_max) * dlen + len_m[vmin];
	// count max diff
	double a1 = len_m[vmin] / (1.0 - 0.5*a_m);
	double max_diff_ratio = 1.0 - (1.0 - dr)*dlen / a1;
	double real_dr = dr;

	bool result = false;
	if (real_max_ratio_1 > max_diff_ratio) {
		if (cn1->control_data.setMinimum(cn2->control_data * max_diff_ratio))
			result = true;
	}
	else real_dr = 1.0 - a1 * (1.0 - real_max_ratio_1) / dlen;

	if (cn1->gradationUnknown() || real_dr > cn1->max_gradation_ratio)
		cn1->max_gradation_ratio = real_dr;

	return result;
}

/// Smoothen variance of metric within the control space (with specified gradation ratio)

bool ControlSpace3dAdaptive::smoothenWithGradation(double gradation) {
	double tmp = ControlSpace2dAdaptive::param_gradation_ratio;
	ControlSpace2dAdaptive::param_gradation_ratio = gradation;
	bool result = smoothen();
	ControlSpace2dAdaptive::param_gradation_ratio = tmp; return result;
}

bool ControlSpace3dAdaptive::logMetricAlongSegment(const string& fname,
		const DPoint3d& pt0, const DPoint3d& pt1, 
		const DVector3d& vn, int n)
{
	ofstream ofs(fname.c_str());

	ofs << "Lp\tvn\tvr\tvt" << endl;
	double dt = 1.0 / (double)n;

	const DVector3d vn_norm = vn.normalized();
	const DVector3d vr = pt1-pt0;
	const DVector3d vr_norm = vr.normalized();
	const DVector3d vt_norm = vn.crossProduct(vr_norm).normalized();

	Metric3dContext mc(this);
	int counter = 0;
	for(double t = 0.0; t <= 1.0; t += dt){
		const DPoint3d pt = pt0 + vr * t;
		mc.countMetricAtPoint(pt);
		ofs << ++counter << '\t';
		ofs << (1.0 / mc.transformRStoMS(vn_norm).length()) << '\t';
		ofs << (1.0 / mc.transformRStoMS(vr_norm).length()) << '\t';
		ofs << (1.0 / mc.transformRStoMS(vt_norm).length()) << endl;
	}
	return true;
}

bool ControlSpace3dAdaptive::testMetric(const string& test_name)
{
	const DPoint3d middle = m_box.getMiddlePoint();
	const DVector2d dv(m_box.getDiameter(), 0.0);
	const DVector3d vn(0.0, 0.0, 1.0);

	double dvx = m_box.getDX()/2;
	double dvy = m_box.getDY()/2;

	double alpha = 0.0;
	const int N = 10;
	for(int i = 0; i < N; i++, alpha += PI/N){
		double sin_alpha = sin(alpha);
		double cos_alpha = cos(alpha);
		DVector2d dvt = dv.turned(sin_alpha, cos_alpha);
		if(abs(dvt.x) > dvx) dvt *= (dvx / abs(dvt.x));
		if(abs(dvt.y) > dvy) dvt *= (dvy / abs(dvt.y));
		const DVector3d dv3(dvt.x, dvt.y, 0.0);

		const DPoint3d pt0 = m_box.fitInPoint(middle - dv3);
		const DPoint3d pt1 = m_box.fitInPoint(middle + dv3);

		ostringstream oss;
		oss << "metric-test-" << test_name << "-" << i << ".txt";

		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Metric test [" << oss.str() << " ], pt0=" << pt0 << ", pt1=" << pt1);

		logMetricAlongSegment(oss.str(), pt0, pt1, vn);
	}

	return true;
}

/// Metric gradation tests

double ControlSpace3dAdaptive::calculateMaxGradationRatioViaNodes() const 
{
	int cn_ct = getControlNodesCount();
	DataVector<const ControlNode3d*> cnodes(cn_ct);
	forEachControlNode([&cnodes](const ControlNode3d& cn) { cnodes.add(&cn); });

	//CDM3d cdm2 = DMetric3d::stretchToMatrix(ControlDataStretch3d(2.0));
	//CDM3d cdm4 = DMetric3d::stretchToMatrix(ControlDataStretch3d(4.0));
	//CDM3d cdm6 = DMetric3d::stretchToMatrix(ControlDataStretch3d(6.0));
	//DVector3d dv(3.0, 0.0, 0.0);
	//double dr24 = DMetric3d::getMetricGradation(cdm2, cdm4, dv);
	//double dr26 = DMetric3d::getMetricGradation(cdm2, cdm6, dv);


	double max_dr = 0.0;
	for (int i = 0; i < cn_ct - 1; i++) {
		for (int j = i + 1; j < cn_ct; j++) {
			double dr = DMetric3d::getMetricGradation(cnodes[i]->control_data, cnodes[j]->control_data,
				cnodes[j]->coord - cnodes[i]->coord);
//			LOG4CPLUS_INFO(MeshLog::logger_mesh, "i=" << i << " j=" << j << "\tdr=" << dr);
			if (dr > max_dr) max_dr = dr;
		}
	}
	return max_dr; 
}

double ControlSpace3dAdaptive::calculateMaxGradationRatioRandom(int probe_count) const 
{ 
	double max_dr = 0.0;
	for (int i = 0; i < probe_count; i++) {
		DPoint3d pt1 = m_box.getRandomPoint();
		DPoint3d pt2 = m_box.getRandomPoint();
		double dr = DMetric3d::getMetricGradation(
			getMetricAtPoint(pt1), getMetricAtPoint(pt2), pt2 - pt1);
		if (dr > max_dr) max_dr = dr;
	}
	return max_dr;
}

double ControlSpace3dAdaptive::calculateMaxGradationRatioRegular(int csize) const 
{ 
	int csize2 = csize*csize;
	CDM3d* cache = new CDM3d[csize2 * csize];

	double fx = 1.0 / (double)csize;
	double cdx = m_box.getDX() * fx;
	double cdy = m_box.getDY() * fx;
	double cdz = m_box.getDZ() * fx;

	DVector3d dvx(cdx, 0.0, 0.0);
	DVector3d dvy(0.0, cdy, 0.0);
	DVector3d dvz(0.0, 0.0, cdz);

	DPoint3d pt(m_box.x0 + 0.5*cdx, m_box.y0, m_box.z0);
	int ix = 0, iy, iz;
	for (; pt.x < m_box.x1; pt.x += cdx, ix++) {
		for (iy = 0, pt.y = m_box.y0 + 0.5*cdy; pt.y < m_box.y1; pt.y += cdy, iy++) {
			for (iz = 0, pt.z = m_box.z0 + 0.5*cdz; pt.z < m_box.z1; pt.z += cdz, iz++) {
				cache[ix*csize2 + iy * csize + iz] = getMetricAtPoint(pt);
			}
		}
	}

	double dr, max_dr = 0.0;
	for (int i = 0; i < csize; i++) {
		for (int j = 0; j < csize; j++) {
			for (int k = 1; k < csize; k++) {
				int k0 = k - 1;
				dr = DMetric3d::getMetricGradation(
					cache[i*csize2 + j * csize + k], cache[i*csize2 + j * csize + k0], dvz);
				if (dr > max_dr) max_dr = dr;
				dr = DMetric3d::getMetricGradation(
					cache[i*csize2 + k * csize + j], cache[i*csize2 + k0 * csize + j], dvy);
				if (dr > max_dr) max_dr = dr;
				dr = DMetric3d::getMetricGradation(
					cache[k*csize2 + i * csize + j], cache[k0*csize2 + i * csize + j], dvx);
				if (dr > max_dr) max_dr = dr;
			}
		}
	}
	return max_dr;
}

/*
bool ControlSpace3dAdaptive::updateForBoundarySegment(const MeshEdge3d* edge)
{
	const DPoint3d& pt0 = edge->getMeshPoint(0)->getCoordinates();
	const DPoint3d& pt1 = edge->getMeshPoint(1)->getCoordinates();

	const DVector3d v01 = pt1-pt0;
	double len = 1.2*v01.length();
	const DVector3d e0 = v01.normalized();
	DVector3d e1,e2;
	e0.orthonormalVectors(e1, e2);
	double d[3] = {len, len*param_stretch_max_ratio, len*param_stretch_max_ratio};
	ControlDataMatrix3d cdm(e0, e1, e2, d);	// from eigensystem

	ControlDataExtMatrix3dSegment data(pt0, v01, d[0], cdm);
	return setMinControl(data);
}
*/

bool ControlSpace3dAdaptive::setMinControl(const ControlDataExtMatrix3dSphere& data)
{
	assert(m_initialized > 0 && m_initialized < 10);

	double r = data.totalRadius();
	if(r <= 0.0)
		return setMinControl(data.getMiddle(), data.getInnerData(), true);

	bool any_changes = false;
	if(data.getRadius2() > 0.0){
		any_changes |= setSphereMinControl(data.getMiddle(), data.getOuterData(), data.totalRadius());
	}

	if(data.getRadius1() > 0.0){
		any_changes |= setSphereMinControl(data.getMiddle(), data.getInnerData(), data.getRadius1());
	}

	any_changes |= (setMinControlValue(data) > 0);

	return any_changes;
}

bool ControlSpace3dAdaptive::setMinControl(const ControlDataExtMatrix3dSegment& data)
{
	assert(m_initialized > 0 && m_initialized < 10);

	double seg_len = data.getDv().length();
	if(seg_len < mesh_data.relative_small_number) return false;

	double dt1 = std::min(0.25, 2.0*data.getInnerData().minEigenvalue() / seg_len);


	double r = data.totalRadius();
	if(r <= 0.0)
		return setSegmentMinControl(data.getStart(), data.getDv(), dt1, 
			data.getInnerData(), 0.0);

	bool any_changes = false;
	if(data.getRadius2() > 0.0){
		double dt2 = std::min(0.25, 2.0*data.getOuterData().minEigenvalue() / seg_len);
		any_changes |= setSegmentMinControl(data.getStart(), data.getDv(), dt2, 
							data.getOuterData(), data.totalRadius());
	}

	if(data.getRadius1() > 0.0){
		any_changes |= setSegmentMinControl(data.getStart(), data.getDv(), dt1, 
							data.getInnerData(), data.getRadius1());
	}

	any_changes |= (setMinControlValue(data) > 0);

	return any_changes;
}

bool ControlSpace3dAdaptive::setMinControl(const ControlDataExtMatrix3dTriangle& data)
{
	assert(m_initialized > 0 && m_initialized < 10);

	double seg_len[2];
	double dt[2];
	double inner_min = data.getInnerData().minEigenvalue();
	for(int i = 0; i < 2; i++){
		seg_len[i] = data.getDv(i).length();
		if(seg_len[i] < mesh_data.relative_small_number) return false;
		dt[i] = std::min(0.1, 2.0*inner_min / seg_len[i]);
	}

	double r = data.totalRadius();
	if(r <= 0.0)
		return setTriangleMinControl(data.getStart(), 
			data.getDv(0), data.getDv(1), dt[0], dt[1],
			data.getInnerData(), 0.0);

	bool any_changes = false;
	if(data.getRadius2() > 0.0){
		double dt_outer[2];
		double outer_min = data.getOuterData().minEigenvalue();
		for(int i = 0; i < 2; i++)
			dt_outer[i] = std::min(0.1, 2.0*outer_min / seg_len[i]);
		any_changes |= setTriangleMinControl(data.getStart(), 
			data.getDv(0), data.getDv(1), dt_outer[0], dt_outer[1],
			data.getOuterData(), data.totalRadius());
	}

	if(data.getRadius1() > 0.0){
		any_changes |= setTriangleMinControl(data.getStart(), 
			data.getDv(0), data.getDv(1), dt[0], dt[1],
			data.getInnerData(), data.getRadius1());
	}

	any_changes |= (setMinControlValue(data) > 0);

	return any_changes;
}

bool ControlSpace3dAdaptive::setSphereMinControl(
	const DPoint3d& pt, const ControlDataMatrix3d& cdm, double r)
{
	assert(m_initialized > 0 && m_initialized < 10);
	if(r > 0.0){
		DVector3d nvx(r, 0.0, 0.0);
		DVector3d nvy(0.0, r, 0.0);
		DVector3d nvz(0.0, 0.0, r);

		double min_len = cdm.minEigenvalue();

		const double HALF_PI = 0.5 * PI;
		double arc_len = HALF_PI * r; // one-fourth of the full circle
		double arc_ratio = arc_len / min_len;
		int nt = 1 + (int)arc_ratio;

		bool any_changes = false;

		// -> big circle
		double dalpha = HALF_PI / nt;
		for(int j = 0; j < nt; j++){
			double sa = sin(j * dalpha);
			double ca = cos(j * dalpha);
			any_changes |= setMinControlTranslated(pt, nvx *  sa + nvy *  ca, cdm, false);
			any_changes |= setMinControlTranslated(pt, nvx * -sa + nvy *  ca, cdm, false);
			any_changes |= setMinControlTranslated(pt, nvx *  sa + nvy * -ca, cdm, false);
			any_changes |= setMinControlTranslated(pt, nvx * -sa + nvy * -ca, cdm, false);
		}
		// -> smaller ones (smaller circles, and lower/higher)
		for(int i = 1; i < nt; i++){
			double sratio = cos(i * dalpha) * arc_ratio;
			int snt = 1 + (int)sratio;
			double sdalpha = HALF_PI / snt;
			for(int j = 0; j < snt; j++){
				double sa = sin(j * sdalpha);
				double ca = cos(j * sdalpha);
				DVector3d snvz = nvz * sin(i * dalpha);
				any_changes |= setMinControlTranslated(pt, nvx *  sa + nvy *  ca + snvz, cdm, false);
				any_changes |= setMinControlTranslated(pt, nvx * -sa + nvy *  ca + snvz, cdm, false);
				any_changes |= setMinControlTranslated(pt, nvx *  sa + nvy * -ca + snvz, cdm, false);
				any_changes |= setMinControlTranslated(pt, nvx * -sa + nvy * -ca + snvz, cdm, false);
				any_changes |= setMinControlTranslated(pt, nvx *  sa + nvy *  ca - snvz, cdm, false);
				any_changes |= setMinControlTranslated(pt, nvx * -sa + nvy *  ca - snvz, cdm, false);
				any_changes |= setMinControlTranslated(pt, nvx *  sa + nvy * -ca - snvz, cdm, false);
				any_changes |= setMinControlTranslated(pt, nvx * -sa + nvy * -ca - snvz, cdm, false);
			}
		}
		// -> poles
		any_changes |= setMinControlTranslated(pt,  nvz, cdm, false);
		any_changes |= setMinControlTranslated(pt, -nvz, cdm, false);
		// -> finish
		return any_changes;
	}else{
		return setMinControl(pt, cdm, false);
	}
}

bool ControlSpace3dAdaptive::setSegmentMinControl(
	const DPoint3d& pt0, const DVector3d& dv, double dt,
	const ControlDataMatrix3d& cdm, double r)
{
	assert(m_initialized > 0 && m_initialized < 10);
	if(r > 0.0){
		DVector3d nv1, nv2;
		dv.orthonormalVectors(nv1, nv2);
		nv1 *= r;
		nv2 *= r;

		const double TWO_PI = 2.0 * PI;
		double dx = dt * dv.length(); //
		double dtr = dx / r;

		DataVector<double> sin_tr(2 + (int)(TWO_PI / dtr));
		DataVector<double> cos_tr(2 + (int)(TWO_PI / dtr));
		for(double tr = 0.0; tr < TWO_PI; tr += dtr){
			sin_tr.add(sin(tr));
			cos_tr.add(cos(tr));
		}

		bool any_changes = false;
		for(double t = 0.0; t <= 1.0; t += dt){
			const DPoint3d pt = m_box.fitInPoint(pt0, dv * t);
			any_changes |= setMinControl(pt, cdm, false);
			// circle
			for(int i = 0; i < sin_tr.countInt(); i++)
				any_changes |= setMinControlTranslated(pt, nv1 * sin_tr[i] + nv2 * cos_tr[i], cdm, false);
		}
		for(int k = 0; k < 3; k++){
			const DPoint3d pt = pt0 + dv * (k*0.5);
			// disk
			for(double s = dtr; s < 1.0; s += dtr)
				for(int i = 0; i < sin_tr.countInt(); i++)
					any_changes |= setMinControlTranslated(pt, nv1 * (s * sin_tr[i]) + nv2 * (s * cos_tr[i]), cdm, false);
		}
		return any_changes;
	}else{
		bool any_changes = false;
		for(double t = 0.0; t <= 1.0; t += dt)
			any_changes |= setMinControlTranslated(pt0, dv * t, cdm, true);
		return any_changes;
	}
}

bool ControlSpace3dAdaptive::setTriangleMinControl(
	const DPoint3d& pt0, const DVector3d& dv0, const DVector3d& dv1, 
	double dt0, double dt1, const ControlDataMatrix3d& cdm, double r)
{
	assert(m_initialized > 0 && m_initialized < 10);
	bool any_changes = false;
	if(r > 0.0){
		const DVector3d dn = dv0.crossProduct(dv1).normalized() * r;
		const DVector3d dnr = dn * -1.0;
		for(double t0 = 0.0; t0 <= 1.0; t0 += dt0){
			for(double t1 = 0.0; t1 < 1.0-t0; t1 += dt1){
				const DPoint3d pt = pt0 + dv0 * t0 + dv1 * t1;
				if(m_box.contains(pt)){
					any_changes |= setMinControlTranslated(pt, dn,  cdm, false);
					any_changes |= setMinControlTranslated(pt, dnr, cdm, false);
				}
			}
		}
	}else{
		for(double t0 = 0.0; t0 <= 1.0; t0 += dt0){
			for(double t1 = 0.0; t1 < 1.0-t0; t1 += dt1){
				const DPoint3d pt = pt0 + dv0 * t0 + dv1 * t1;
				if(m_box.contains(pt))
					any_changes |= setMinControl(pt, cdm, true);
			}
		}
	}
	return any_changes;
}

bool ControlSpace3dAdaptive::setMinControlTranslated(const DPoint3d& pt, const DVector3d& dv, 
		const ControlDataMatrix3d& cdm, bool min_value_set)
{
	return setMinControl(m_box.fitInPoint(pt, dv), cdm, min_value_set);
}

/// Invoke (test) function in regular grid over this space
void ControlSpace3dAdaptive::forRegularGrid(int grid_size, const std::function<void(const DPoint3d&pt)>& fg) {
	m_box.forRegularGrid(grid_size, fg);
}

/// Refines control space basing on other control space, within given rectangular area
bool ControlSpace3dAdaptive::setMinControlRect
	(Metric3dContext& mc, CS3dPtr space,
	const DPoint3d& pt00, const DVector3d& dv0, const DVector3d& dv1, 
	const DVector3d& dvn, double r)
{
	const DPoint3d& middle = pt00 + (dv0 + dv1)*0.5;

	// apply metric in the middle point
	mc.countMetricAtPoint(middle);
	const DVector3d xdvn = dvn / mc.transformRStoMS(dvn).length();
	bool any_change = setMinControlPointVicinity(space, middle, xdvn, r);

	// check dimensions of quadtree cell
	double len0 = mc.transformRStoMS(dv0).length2();
	double len1 = mc.transformRStoMS(dv1).length2();
	//if(std::max(len0, len1) < 0.25) return any_change; // sqrt(len) < 0.5
	if(std::max(len0, len1) < 1.0) return any_change; // sqrt(len) < 1.0
	// need further refinement
	if(len0 > len1){
		const DVector3d xdv0 = dv0 * 0.5;
		any_change |= setMinControlRect(mc, space, pt00, xdv0, dv1, dvn, r);
		any_change |= setMinControlRect(mc, space, pt00+xdv0, xdv0, dv1, dvn, r);
	}else{
		const DVector3d xdv1 = dv1 * 0.5;
		any_change |= setMinControlRect(mc, space, pt00, dv0, xdv1, dvn, r);
		any_change |= setMinControlRect(mc, space, pt00+xdv1, dv0, xdv1, dvn, r);
	}

	return any_change;
}

/// Refines control space basing on other control space, along given segment
bool ControlSpace3dAdaptive::setMinControlLine(Metric3dContext& mc, CS3dPtr space,
	const DPoint3d& pt0, const DPoint3d& pt1, 
	const DVector3d& dvn, double r)
{
	double t0 = 0.0;
	double dt = 0.01;
	DPoint3d xpt0 = pt0;
	const DVector3d dv = pt1 - pt0;
	bool last_time = false;
	bool any_change = false;
	while(true)
	{
		mc.countMetricAtPoint(xpt0);
		// calculate proper length of normal vector
		const DVector3d xdvn = dvn / mc.transformRStoMS(dvn).length();
		// insert new metric source
		any_change |= setMinControlPointVicinity(space, xpt0, xdvn, r);

		if(last_time) break;
		// calculate next point
		double t1 = t0 + dt;
		DPoint3d xpt1 = pt0 + dv * t1;
		// normalize for unit (metric) segment along selected line
		dt /= mc.transformRStoMS(xpt1 - xpt0).length();
		// calculate proper point
		t0 += 0.5*dt; // move by half of "unit" segment length
		if(t0 >= 1.0){ t0 = 1.0; last_time = true; }
		xpt0 = pt0 + dv * t0;
	}

	return any_change;
}

bool ControlSpace3dAdaptive::setMinControlPointVicinity(CS3dPtr space,
	const DPoint3d& pt, const DVector3d& dvn, double r, double dt)
{
	bool any_change = false;
	// insert new metric source
	// ... at point
	any_change |= setMinControl(pt, space->getMetricAtPoint(pt));
	// ... at width (in both directions)
	DPoint3d xpt = pt + dvn * r;
	any_change |= setMinControl(xpt, space->getMetricAtPoint(xpt));
	xpt = pt - dvn * r;
	any_change |= setMinControl(xpt, space->getMetricAtPoint(xpt));
	// ... and in between
	for(double t = dt; t < r; t += dt){
		xpt = pt + dvn * t;
		any_change |= setMinControl(xpt, space->getMetricAtPoint(xpt));
		xpt = pt - dvn * t;
		any_change |= setMinControl(xpt, space->getMetricAtPoint(xpt));
	}
	return any_change;
}

/// Calculates control space for the given point using curvature of the given surface
bool ControlSpace3dAdaptive::calculateControlFromSurfaceCurvature(const DPoint3d& pt, 
		SurfaceConstPtr surface, double model_diameter, 
		ControlDataMatrix3d& surf_cdm, const DPoint2d * pt2d, SurfaceCurvature * curvature2d)
{
	assert(surface);
	double g_ratio;
	DPoint2d cpt2d = pt2d ? *pt2d : surface->getParameters(pt);
	const SurfaceCurvature curvature = surface->getCurvature(cpt2d, g_ratio);
	if( curvature2d ) *curvature2d = curvature;
	if(g_ratio <= MIN_PARAM_GRATIO) return false;


	double max_len = model_diameter * std::min(ControlSpace2dAdaptive::param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
	double min_len = std::max(model_diameter * ControlSpace2dAdaptive::param_min_diameter_ratio, ControlSpace2dAdaptive::param_min_length);

	const ControlDataStretch2d data = ControlSpace2dAdaptive::adjustCurvatureData(curvature, 
				ControlSpace2dAdaptive::param_curvature_ratio, model_diameter, min_len, max_len, 
				ControlSpace2dAdaptive::param_stretch_max_ratio);
	// Przeliczenie (lx, ly, alfa) -> (dxx, dxy, dyy) -> 3D
	return DMetric3d::projectCDMto3D(DMetric2d::stretchToMatrix(data), surface, cpt2d, surf_cdm);
}

// Refines control space at the given point with the given extended metric
	/// Refines control space at the given point using curvature of the given surface
bool ControlSpace3dAdaptive::setMinControlFromSurfaceCurvature(const DPoint3d& pt, 
		SurfaceConstPtr surface, double model_diameter, bool min_value_set,
		ControlDataMatrix3d* surf_cdm, SurfaceCurvature *curvature2d )
{
	ControlDataMatrix3d cdm;
	if( calculateControlFromSurfaceCurvature(pt, surface, model_diameter, cdm, nullptr, curvature2d) ) {
		if(surf_cdm) *surf_cdm = cdm;
		return setMinControl(pt, cdm, min_value_set);
	}else
		return false;
}

// Refines control space at the given point with the given extended metric
	/// Refines control space at the given point using curvature of the given surface
bool ControlSpace3dAdaptive::setMinControlFromSurfaceCurvature(MeshFace* face, 
		SurfaceConstPtr surface, double model_diameter, bool /* min_value_set */,
		ControlDataMatrix3d* surf_cdm, SurfaceCurvature *curvature2d )
{
	ControlDataMatrix3d cdm;
	const DPoint3d fmiddle = face->getMiddlePoint();
	if( calculateControlFromSurfaceCurvature( fmiddle, surface, model_diameter, cdm, nullptr, curvature2d) ) {
		if(surf_cdm) *surf_cdm = cdm;

		DataVector< DTriangle3d > triangles;
		bool any_change = false;
		if( face->splitToTriangles( triangles ) ) {
			for(int i = 0; i < triangles.countInt(); i++) {
				const DTriangle3d& dtri = triangles[i];
				any_change |= setMinControl( ControlDataExtMatrix3dTriangle( dtri.pt_a, 
					dtri.pt_b - dtri.pt_a, dtri.pt_c - dtri.pt_a, 0.0, cdm ) );
			}
		}

		return any_change;
	}else
		return false;
}

/// Refines control space at the given point using curvature of the given curve
bool ControlSpace3dAdaptive::setMinControlFromCurveCurvature(const DPoint3d& pt, Curve3dConstPtr curve, 
															 double model_diameter, bool min_value_set)
{
	assert(curve);

	assert(false); // should have given "near_t" parameter somehow...
	double t = curve->getParameter(pt, 0.0);
	double c = curve->getCurvature(t);

	if(c < mesh_data.relative_small_number) return false;
	double max_len = model_diameter * std::min(ControlSpace2dAdaptive::param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
	double min_len = std::max(model_diameter * ControlSpace2dAdaptive::param_min_diameter_ratio, ControlSpace2dAdaptive::param_min_length);

	double len = ControlSpace2dAdaptive::param_curvature_ratio / c;
	if(len > max_len) len = max_len;
	else if(len < min_len) len = min_len; 

	DVector3d e0 = curve->getDerivative(t);
	DVector3d e1, e2;

	e0.orthogonalVectors(e1, e2);
	double d[] = { len, 
		len * ControlSpace2dAdaptive::param_stretch_max_ratio,
		len * ControlSpace2dAdaptive::param_stretch_max_ratio };

	return setMinControl(pt, ControlDataMatrix3d(e0, e1, e2, d), min_value_set);
}
