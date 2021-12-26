// ControlSpace2dAdaptive.cpp: implementation of the ControlSpace2dAdaptive class.
//
//////////////////////////////////////////////////////////////////////

#include <algorithm>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "ControlSpace2dAdaptive.h"
#include "SurfaceParametric.h"
#include "Curve2dParametric.h"
#include "Curve2dSegment.h"
#include "DEquation.h"
#include "MeshData.h"
#include "MeshContainer2d.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "MeshEdge2dCurve.h"
#include "EPSFile.h"
#include "ControlSpace3d.h"
#include "ControlSpace3dAdaptive.h"
#include "Metric2dContext.h"

#include "MeshViewSet.h"

#ifdef _DEBUG
//#define STORE_CONTROL_EPS
#endif


double ControlSpace2dAdaptive::param_curvature_ratio = 0.15;
double ControlSpace2dAdaptive::param_max_diameter_ratio = 0.5;
double ControlSpace2dAdaptive::param_inner_boundary_ratio = 1.0;
double ControlSpace2dAdaptive::param_min_diameter_ratio = 0.0001;
double ControlSpace2dAdaptive::param_min_length = 0.0;
double ControlSpace2dAdaptive::param_stretch_max_ratio = 2.0;
int ControlSpace2dAdaptive::param_use_surface_curvature = 1;
int ControlSpace2dAdaptive::param_use_contour_curvature = 1;
double ControlSpace2dAdaptive::param_threshold_diff = 0.2;

/// Number of probe points for radial metric definition in adaptive control space
int ControlSpace2dAdaptive::param_radial_parts = 8;

double ControlSpace2dAdaptive::param_contour_stretch_max_ratio = 2.0;
double ControlSpace2dAdaptive::param_gradation_ratio = 2.0;

/// starting X*Y resolution for quadtree control/uniform probing, ect.
int ControlSpace2dAdaptive::param_control_nxy = 100;

void ControlSpace2dAdaptive::addControlSegment(
	const DPoint2d &pt1, const DPoint2d &pt2, 
	const ControlDataMatrix2d& data, double r)
{
	DMetric2d dmp(base_surface, DPoint2d::average(pt1, pt2));
	double seg_len = dmp.transformPStoRS(pt2-pt1).length();
	if(seg_len < mesh_data.relative_small_number) return;

	const DVector2d dv = pt2-pt1;
	double dt = std::min(0.01, 2.0 * data.minEigenvalue() / seg_len);

	addSegmentControlPoint(pt1, dv, dt, data, r);

	if(r > 0.0)
		m_ext_source_points.add(
			std::make_shared<ControlDataExtMatrix2dSegment>(pt1, dv,
			r, data, base_surface));
}

void ControlSpace2dAdaptive::addControlSegment(
	Curve2dConstPtr curve, double t0, double t1, 
	const ControlDataMatrix2d& data, double r)
{
	DataVector<double> polyline(50);
	curve->getPolyLineInRange(t0, t1, polyline);

	int ct = (int)polyline.countInt();
	assert(ct > 1);
	DPoint2d pt0 = curve->getPoint(polyline.get(0));
	for(int i = 1; i < ct; i++){
		DPoint2d pt1 = curve->getPoint(polyline.get(i));
		addControlSegment(pt0, pt1, data, r);
		pt0 = pt1;
	}
}

bool ControlSpace2dAdaptive::setContourCurvatureControlData(const MeshContainer2d* boundary)
{
	double p = m_box.getDiameter();
	double max_len = p * std::min(param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
	double min_len = std::max(p * param_min_diameter_ratio, param_min_length);

	bool any_point = false;
	if(m_initialized == 0){
		ControlDataMatrix2d max_data(max_len, 0.0, max_len);
		forEachControlNode([&](ControlNode2d& node) {
			node.control_data = max_data;
			node.w = -1.0;
		});
		m_initialized = 1;
	}

	bool non_planar = true;	// should be: (surface->getType() != SURFACE_PLANE) or (non-standard parameterization);
	const double rf = 1.5;

	double stretch_max_ratio = std::min(
		ControlSpace2dAdaptive::param_contour_stretch_max_ratio,
		ControlSpace2dAdaptive::param_stretch_max_ratio);



	int pct = boundary->getPointsCount();
	for(int i = 0; i < pct; i++){

		MeshPoint2d* point = boundary->getPointAt(i);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0)
				continue;
			if(non_planar){
				double dksi = 0.01;
				for(double ksi = 0.0; ksi <= 1.0; ){
					double c_dlen;
					double c = edge->getNonPlanarCurvature(base_surface, ksi, &c_dlen);

					if(c < mesh_data.relative_small_number){
						ksi += dksi; continue;
					}
					double len = param_curvature_ratio / c;
					if(len > max_len) len = max_len;
					else if(len < min_len) len = min_len; 
					double c_dksi = 0.2*(len / c_dlen);

					DPoint2d pt = (ksi < 1.0) ? edge->getPoint(ksi) : edge->getPoint(ksi - 0.1*c_dksi);
					const DVector2d dpt = edge->getPoint(std::min(ksi+c_dksi, 1.0)) - pt;

					DMetric2d dmp(base_surface, pt);
					double angle = dmp.transformPStoRS(DVector2d(len, 0.0)).getAngle(
						dmp.transformPStoRS(dpt));
					ControlDataStretch2d cds(len, max_len, angle);
					if(max_len/len > stretch_max_ratio)
						cds.ly = len*stretch_max_ratio;
					const ControlDataMatrix2d contour_cdm = DMetric2d::stretchToMatrix(cds);
					const ControlDataMatrix2d local_cdm = getMetricAtPoint(pt);
					double diff = local_cdm.countDifferenceRR(contour_cdm);
					if(diff > param_threshold_diff){
						ControlDataMatrix2d min_cdm = contour_cdm;
						min_cdm.setMinimum(local_cdm);
						diff = local_cdm.countDifferenceRR(min_cdm);
						if(diff > param_threshold_diff){
//							ControlDataExtMatrix2dSegment data(pt, dpt, rf*len, 
//								contour_cdm, rf*max_len, max_data, base_surface);
							//ControlDataExtMatrix2dRadial data(pt, rf*len, contour_cdm, base_surface);
							//setMinControl(data);
							setMinControl(pt, contour_cdm);
							any_point = true;
						}
					}
					ksi += c_dksi;
				}
			}else if(edge->getType() != EDGE_SIMPLE){
				double dksi = 0.01;
				for(double ksi = 0.0; ksi <= 1.0; ksi += dksi){
					double c = edge->getPlanarCurvature(ksi);

					if(c < mesh_data.relative_small_number)
						continue;
					double len = param_curvature_ratio / c;
					if(len > max_len) len = max_len;
					else if(len < min_len) len = min_len; 

					DPoint2d pt = (ksi < 1.0) ? edge->getPoint(ksi) : edge->getPoint(ksi - 0.5*dksi);
					const DVector2d dv = edge->getPoint(std::min(ksi+dksi, 1.0)) - pt;

					DMetric2d dmp(base_surface, pt);
					double angle = dmp.transformPStoRS(DVector2d(len, 0.0)).getAngle(
						dmp.transformPStoRS(dv));
					ControlDataStretch2d cds(len, max_len, angle);
					if(max_len/len > stretch_max_ratio)
						cds.ly = len*stretch_max_ratio;
					ControlDataExtMatrix2dSegment data(pt, dv, rf*len, 
						DMetric2d::stretchToMatrix(cds), base_surface);

					setMinControl(data);
					any_point = true;
				}
			}
		}
	}

#ifdef STORE_CONTROL_EPS
	storeEPS("control-contour-curvature");
#endif

	return any_point;
}

void ControlSpace2dAdaptive::setMaxMetric(double ratio)
{
	const DMetric2d dmp(base_surface, m_box.getMiddlePoint());
	double p = dmp.transformPStoRS(m_box.getX1Y1() - m_box.getX0Y0()).length();
//	double p = m_box.getDiameter();
	double max_len = p * (ratio > 0 ? ratio : 
		std::min(param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio));

//	LOG4CPLUS_INFO(MeshLog::logger_console, "CS3A::diameter", p);
//	LOG4CPLUS_INFO(MeshLog::logger_console, "CS3A::max_len", max_len);

	ControlDataMatrix2d max_data(max_len, 0.0, max_len);
	forEachControlNode([&](ControlNode2d& node) {
		node.control_data = max_data;
		node.w = -1.0;
	});
	m_initialized = 1;
}

// no adaptation, just set min value
int ControlSpace2dAdaptive::setMinControlValue(const ControlDataExtMatrix2d& data)
{
	// 1. check all grid vertices (***)
	int completely_new_count = 0;
	forEachControlNode([&](ControlNode2d& qv) {
		if (data.isPointWithin(qv.coord))
			if (qv.w < 0.0)	// node already initialized
				qv.control_data.setMinimum(data.getControlDataMatrix(qv.coord));
			else {
				qv.control_data = data.getControlDataMatrix(qv.coord);
				qv.w = -1.0;
				++completely_new_count;
			}
	});
	return completely_new_count;
}

bool ControlSpace2dAdaptive::applyAsMinimum(CS2dPtr space)
{
	assert(m_initialized);

	bool any_changes = false;
	// 1: adjust structure of this space with other space
	space->forEachControlNode([&](ControlNode2d& qv) {
		any_changes |= setMinControl(qv.coord, qv.control_data, false);
	});
	// 2: calculate minimum metric for own nodes
	forEachControlNode([&](ControlNode2d& qv) {
		any_changes |= qv.control_data.setMinimum(space->getMetricAtPoint(qv.coord));
	});

#ifdef STORE_CONTROL_EPS
	storeEPS("control-minimum-destination");
	if(space->isAdaptive()){
		((CS2dPtr)space)->storeEPS("control-minimum-source");
	}
#endif

#ifdef STORE_CONTROL_EPS
	storeEPS("control-minimum-result");
#endif
	return any_changes;
}

bool ControlSpace2dAdaptive::applyAsMinimum(CS3dPtr space)
{
	assert(m_initialized);
	bool any_changes = false;

	//showMetricLength("applyAsMinimum - start");
	//storeEPS("cs1");

	if(base_surface->getType() == SURFACE_PLANE){

		//MeshViewSet* set = new MeshViewSet();
		//set->addEdge(base_surface->getPoint(m_box.getX0Y0()), base_surface->getPoint(m_box.getX0Y1()), 1);
		//set->addEdge(base_surface->getPoint(m_box.getX1Y0()), base_surface->getPoint(m_box.getX1Y1()), 1);
		//set->addEdge(base_surface->getPoint(m_box.getX0Y0()), base_surface->getPoint(m_box.getX1Y0()), 1);
		//set->addEdge(base_surface->getPoint(m_box.getX0Y1()), base_surface->getPoint(m_box.getX1Y1()), 1);

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "1: adjust structure of this space with other space");

		// 1: adjust structure of this space with other space
		space->forEachControlNode([&](const ControlNode3d & qv) {
			const DPoint2d pt2d = base_surface->getParameters(qv.coord);
			if (m_box.contains(pt2d)) {
				const DPoint3d pt3d = base_surface->getPoint(pt2d);
				double dv_len = (qv.coord - pt3d).length();
				//double cs_len = qv.control_data.minEigenvalue();
				bool close_to_surface = dv_len < mesh_data.relative_small_number;
				if (!close_to_surface) {
					double max_v = qv.control_data.maxEigenvalue();
					close_to_surface = (dv_len < max_v);
				}
				//set->addPoint(qv.coord, close_to_surface ? (
				//	(dv_len < mesh_data.relative_small_number)?2:3) : 4);
				if (close_to_surface) {
					auto cdm3d = space->getMetricAtPoint(pt3d);
					const ControlDataMatrix2d cdm = base_surface->projectTransformationTensor(
						pt2d, cdm3d);
					any_changes |= setMinControl(pt2d, cdm, false);
				}
			}
			else {
				//set->addPoint(qv.coord, 1);
			}
		});

		//SHOW_MESH("cs3d cnodes", set);
	}

	//storeEPS("cs2");
	//showMetricLength("applyAsMinimum - cs3d cnodes into this");

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "2: calculate minimum metric for own nodes");

	// 2: calculate minimum metric for own nodes
	forEachControlNode([&](ControlNode2d & qv) {
		const DPoint3d spt = base_surface->getPoint(qv.coord);
		if (space->containsPoint(spt)) {
			const ControlDataMatrix3d cdm3d = space->getMetricAtPoint(spt);
			const ControlDataMatrix2d cdm = base_surface->projectTransformationTensor(
				qv.coord, cdm3d);
			// apply
			any_changes |= qv.control_data.setMinimum(cdm);
		}
	});

	//showMetricLength("applyAsMinimum - min of own cnodes");

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "3: adapt with respect to cs3d");
	// 3: adapt with respect to cs3d
	any_changes |= adaptToControlSpace3d(space);
#ifdef STORE_CONTROL_EPS
	storeEPS("control-minimum-result");
#endif

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "4: done");
	//showMetricLength("applyAsMinimum - adapt");

	return any_changes;
}

double ControlSpace2dAdaptive::compareMetricWith(CS2dPtr space)
{
	assert(m_initialized);
	assert(!space->isAdaptive() || space->getAsAdaptive()->initializationState() > 0);
	double max_diff = 0.0;

	// 1: check metric for each quad-vertex of this space with other space
	forEachControlNode([&](const ControlNode2d & qv) {
		double diff = qv.control_data.countDifferenceRR(space->getMetricAtPoint(qv.coord));
		if (diff > max_diff) max_diff = diff;
	});
	// 2: ... and the other way round
	space->forEachControlNode([&](const ControlNode2d & qv) {
		double diff = qv.control_data.countDifferenceRR(getMetricAtPoint(qv.coord));
		if (diff > max_diff) max_diff = diff;
	});

	return max_diff;
}

bool ControlSpace2dAdaptive::setMinControl(const ControlDataExtMatrix2dRadial& data)
{
	assert(m_initialized > 0);

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
		double min_len = 2*data.getOuterData().minEigenvalue();
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
		double min_len = 2*data.getInnerData().minEigenvalue();
		const double da = 2*PI/param_radial_parts;
		double a = 0.0;
		for(int i = 0; i < param_radial_parts; i++, a+=da){
			any_changes |= setRadialControlPoint(dm, data.getMiddle(), data.getRadius1(), 
								data.getInnerData(), a, a+da, min_len);
		}
	}

	any_changes |= (setMinControlValue(data) > 0);

#ifdef STORE_CONTROL_EPS
//	storeEPS("control-min-control-after");
#endif
	return any_changes;
}

bool ControlSpace2dAdaptive::setMinControl(const ControlDataExtMatrix2dSegment& data)
{
	assert(m_initialized > 0);

	DMetric2d dmp(base_surface, data.getMiddle());
	double seg_len = dmp.transformPStoRS(data.getDv()).length();
	if(seg_len < mesh_data.relative_small_number) return false;

	double dt1 = std::min(0.1, 2.0*data.getInnerData().minEigenvalue() / seg_len);


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
		double dt2 = std::min(0.1, 2.0*data.getOuterData().minEigenvalue() / seg_len);
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

bool ControlSpace2dAdaptive::setMinControl(const ControlDataExtMatrix2dCurve& data)
{
	bool any_changes = false;
	for(size_t i = 0; i < data.poly_control.countInt(); i++)
		any_changes |= setMinControl(*(data.poly_control[i]));

	return any_changes;
}

bool ControlSpace2dAdaptive::setMinControlTranslated(const DPoint2d& pt, const DVector2d& dv, 
		const ControlDataMatrix2d& cdm, bool min_value_set)
{
	DPoint2d cpt = pt+dv;
	if(!m_box.contains(cpt)){
		// fit
		const DPoint2d fpt = m_box.fitInPoint(cpt);
		double tx = (abs(cpt.x-pt.x) > mesh_data.relative_small_number) ? ((fpt.x-pt.x)/(cpt.x-pt.x)) : 1.0;
		double ty = (abs(cpt.y-pt.y) > mesh_data.relative_small_number) ? ((fpt.y-pt.y)/(cpt.y-pt.y)) : 1.0;
		cpt = pt + dv * std::min(tx, ty);
	}
	return setMinControl(cpt, cdm, min_value_set);
}

bool ControlSpace2dAdaptive::setRadialControlPoint(const DMetric2d& dmp, 
	const DPoint2d& pt, double r, const ControlDataMatrix2d& cdm, 
	double a0, double a1, double min_len)
{
	assert(m_initialized > 0);
	DVector2d v0(1.0, 0.0);
	v0.turn(sin(a0), cos(a0));
	DVector2d v1(1.0, 0.0);
	v1.turn(sin(a1), cos(a1));
	v0 *= (r/dmp.transformPStoRS(v0).length());
	v1 *= (r/dmp.transformPStoRS(v1).length());
	double dist = dmp.transformPStoRS(v0-v1).length2();
	if(dist > min_len){
		bool any_changes = setRadialControlPoint(dmp, pt, r, cdm, a0, 0.5*(a0+a1), min_len);
		any_changes |= setRadialControlPoint(dmp, pt, r, cdm, 0.5*(a0+a1), a1, min_len);
		return any_changes;
	}else{
		return setMinControlTranslated(pt, v0, cdm, false);
	}
}

void ControlSpace2dAdaptive::addControlNodeTranslated(const ControlNode2d& node, const DVector2d& dv)
{
	ControlNode2d cnode = node;
	cnode.coord += dv;
	if(!m_box.contains(cnode.coord)){
		// fit
		const DPoint2d fpt = m_box.fitInPoint(cnode.coord);
		double tx = (abs(cnode.coord.x-node.coord.x) > mesh_data.relative_small_number) 
			? (fpt.x-node.coord.x)/(cnode.coord.x-node.coord.x) : 1.0;
		double ty = (abs(cnode.coord.y-node.coord.y) > mesh_data.relative_small_number) 
			? (fpt.y-node.coord.y)/(cnode.coord.y-node.coord.y) : 1.0;

		cnode.coord = node.coord + dv * std::min(tx, ty);
	}
	addControlNode(cnode);
}

void ControlSpace2dAdaptive::addRadialControlPoint(const DMetric2d& dmp, 
	const DPoint2d& pt, double r, const ControlDataMatrix2d& cdm, 
	double a0, double a1, double min_len)
{
	assert(m_initialized == 0);
	DVector2d v0(1.0, 0.0);
	v0.turn(sin(a0), cos(a0));
	DVector2d v1(1.0, 0.0);
	v1.turn(sin(a1), cos(a1));
	v0 *= (r/dmp.transformPStoRS(v0).length());
	v1 *= (r/dmp.transformPStoRS(v1).length());
	double dist = dmp.transformPStoRS(v0-v1).length2();
	if(dist > min_len){
		addRadialControlPoint(dmp, pt, r, cdm, a0, 0.5*(a0+a1), min_len);
		addRadialControlPoint(dmp, pt, r, cdm, 0.5*(a0+a1), a1, min_len);
	}else{
		addControlNodeTranslated(ControlNode2d(pt, cdm), v0);
	}
}

void ControlSpace2dAdaptive::addSegmentControlPoint(
	const DPoint2d& pt0, const DVector2d& dv, double dt,
	const ControlDataMatrix2d& data, double r)
{
	assert(m_initialized == 0);

	ControlNode2d cn(pt0, data);
	if(r > 0.0){
		const DMetric2d dmp(base_surface, pt0+dv*0.5);
		const DVector2d dv_mp = dmp.transformPStoRS(dv);
		const DVector2d nv_mp(dv_mp.y, -dv_mp.x);
		double rn = r/nv_mp.length();
		const DVector2d nv = dmp.transformRStoPS(nv_mp * rn);
		const DVector2d nv_r = dmp.transformRStoPS(nv_mp * -rn);

		for(double t = 0.0; t <= 1.0; t += dt){
			cn.coord = pt0 + dv * t;
			addControlNodeTranslated(cn, nv);
			addControlNode(cn);
			addControlNodeTranslated(cn, nv_r);
		}
		const DPoint2d pts[3] = {pt0, pt0 + dv * 0.5, pt0 + dv};
		for(double t = dt; t < 1.0; t += dt){
			const DVector2d nvt = nv * t;
			const DVector2d nvrt = nv_r * t;
			for(int i = 0; i < 3; i++){
				cn.coord = pts[i];
				addControlNodeTranslated(cn, nvt);
				addControlNodeTranslated(cn, nvrt);
			}
		}
	}else{		
		for(double t = 0.0; t <= 1.0; t += dt)
			addControlNodeTranslated(cn, dv * t);

	}
}

bool ControlSpace2dAdaptive::setSegmentMinControl(
	const DPoint2d& pt0, const DVector2d& dv, double dt,
	const ControlDataMatrix2d& cdm, double r)
{
	assert(m_initialized > 0);
	if(r > 0.0){
		const DMetric2d dmp(base_surface, pt0+dv*0.5);
		const DVector2d dv_mp = dmp.transformPStoRS(dv);
		const DVector2d nv_mp(dv_mp.y, -dv_mp.x);
		double rn = r/nv_mp.length();
		const DVector2d nv = dmp.transformRStoPS(nv_mp * rn);
		const DVector2d nv_r = dmp.transformRStoPS(nv_mp * -rn);

		bool any_changes = false;
		for(double t = 0.0; t <= 1.0; t += dt){
			const DPoint2d pt = pt0 + dv * t;
			any_changes |= setMinControlTranslated(pt, nv, cdm, false);
			any_changes |= setMinControl(pt, cdm, false);
			any_changes |= setMinControlTranslated(pt, nv_r, cdm, false);
		}
		for(int i = 0; i < 3; i++){
			const DPoint2d pt = pt0 + dv * (i*0.5);
			for(double t = dt; t < 1.0; t += dt){
				any_changes |= setMinControlTranslated(pt, nv * t, cdm, false);
				any_changes |= setMinControlTranslated(pt, nv_r * t, cdm, false);
			}
		}
		return any_changes;
	}else{
		bool any_changes = false;
		for(double t = 0.0; t <= 1.0; t += dt)
			any_changes |= setMinControlTranslated(pt0, dv * t, cdm, false);
		return any_changes;
	}
}

void ControlSpace2dAdaptive::addControlPoint(const DPoint2d& pt, 
	const ControlDataMatrix2d& data, double r)
{
	assert(m_initialized == 0);
	const ControlNode2d qv(pt, data);

	if(r <= 0.0){
		addControlNode(qv);
		return;
	}

	DMetric2d dm(base_surface, qv.coord);

	if(r > 0.0){
		double min_len = 2 * data.minEigenvalue();
		const double da = 2*PI/param_radial_parts;
		double a = 0.0;
		for(int i = 0; i < param_radial_parts; i++, a+=da)
			addRadialControlPoint(dm, pt, r, data, a, a+da, min_len);
	}

	m_ext_source_points.add(
		std::make_shared<ControlDataExtMatrix2dRadial>(pt, r, data, base_surface));
}

ControlDataStretch2d ControlSpace2dAdaptive::adjustCurvatureData(
	const SurfaceCurvature & c, double c_ratio, double /* model_diameter */, 
	double min_len, double max_len, double max_stretch_ratio)
{
	// c -> (x,y,z) = (c_lx, c_ly, angle)
	assert(c.c1 > -mesh_data.relative_small_number);
	assert(c.c2 > -mesh_data.relative_small_number);
	double c1 = (c.c1 > mesh_data.relative_small_number) ? c.c1 : mesh_data.relative_small_number;
	double c2 = (c.c2 > mesh_data.relative_small_number) ? c.c2 : mesh_data.relative_small_number;
	
	ControlDataStretch2d data(c_ratio / c1, c_ratio /c2, c.angle);

	if(data.lx > max_len) data.lx = max_len;	
	if(data.ly > max_len) data.ly = max_len;	
	if(data.lx < min_len) data.lx = min_len; 
	if(data.ly < min_len) data.ly = min_len; 

	if(data.lx / data.ly > max_stretch_ratio) 
		data.lx = data.ly * max_stretch_ratio;
	else if(data.ly / data.lx > max_stretch_ratio) 
		data.ly = data.lx * max_stretch_ratio;

	return data;
}

bool ControlSpace2dAdaptive::updateForBoundarySegment(const MeshEdge2d* edge)
{
	double len = 1.2*edge->getLength(base_surface);
	if(len < param_min_length) len = param_min_length;

	const DPoint2d& pt0 = edge->getMeshPoint(0)->getCoordinates();
	const DPoint2d& pt1 = edge->getMeshPoint(1)->getCoordinates();

	const DVector2d v01 = pt1-pt0;
	const DVector2d du = v01.normalized();
	const DVector2d dv(du.y, -du.x); // orthogonal
	double d[2] = {len, len*param_stretch_max_ratio};
	ControlDataMatrix2d cdm(du, dv, d);	// from eigensystem

	ControlDataExtMatrix2dSegment data(pt0, v01, d[0], cdm, base_surface);
	return setMinControl(data);
}

bool ControlSpace2dAdaptive::updateForBoundaryTriangle(const MeshEdge2d* edge, const MeshPoint2d* point)
{
	const DPoint2d& pt0 = edge->getMeshPoint(0)->getCoordinates();
	const DPoint2d& pt1 = edge->getMeshPoint(1)->getCoordinates();
	const DPoint2d& pt2 = point->getCoordinates();
	const DVector2d du = pt1-pt0;
	const DVector2d dv(du.y, -du.x);
	const DVector2d dpt02 = pt2-pt0;

	double W = du.x*dv.y-du.y*dv.x;
	double t = (dpt02.x*dv.y-dpt02.y*dv.x)/W;
	if(t < 0.0 || t > 1.0) return false;
	double s = (du.x*dpt02.y-du.y*dpt02.x)/W;


	DMetric2d dmp(base_surface, pt0);
	double h = dmp.transformPStoRS(dv*s).length();
	if(h < param_min_length) h = param_min_length;
	double d[2] = { h*param_stretch_max_ratio, h };
	ControlDataMatrix2d cdm(du.normalized(), dv.normalized(), d);

	ControlDataExtMatrix2dRadial data(pt2, h, cdm, base_surface);
	return setMinControl(data);
}

bool ControlSpace2dAdaptive::updateForBoundaryShape(const MeshEdge2d* edge)
{
	double p = m_box.getDiameter();
	double max_len = p * std::min(param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
	double min_len = std::max(p * param_min_diameter_ratio, param_min_length);
	const DPoint2d pt0 = edge->getPoint(0.5);

//	double len = edge->getLength(base_surface);
//	ControlDataExtMatrix2dRadial data(pt0, len, 
//		ControlDataMatrix2d(len/2, 0, len/2), base_surface);

	double ksi[3] = {0.25, 0.5, 0.75 };
	double c = 0.0;
	for(int i = 0; i < 3; i++){
		double c0 = edge->getNonPlanarCurvature(base_surface, ksi[i]);
		if(c0 > c) c = c0;
	}
	assert(c > mesh_data.relative_small_number);
	if(c < mesh_data.relative_small_number) return false;
	double len = param_curvature_ratio / c;
	if(len > max_len) len = max_len;
	else if(len < min_len) len = min_len; 

	ControlDataExtMatrix2dRadial data(pt0, len, 
		ControlDataMatrix2d(len, 0.0, len), base_surface);

	return setMinControl(data);
}

double ControlSpace2dAdaptive::getLocalResolution(const DPoint2d& /* pt */) const {
	return std::min(m_box.getDX(), m_box.getDY());
}

int ControlSpace2dAdaptive::smoothenMetricForNodes(ControlNode2d * cn1, ControlNode2d * cn2, const DVector2d& dv)
{
	double dr = ControlSpace2dAdaptive::param_gradation_ratio;
	assert(dr >= 1.0);
	if(dr < 1.0) dr = 1.0;
	double a_max = 2*(dr-1)/(dr+1);

	double diff = cn1->control_data.countDifferenceRR(cn2->control_data);
	if(diff < 0.01) return 0;

	double dlen = dv.length();
	const DVector2d dvn = dv.normalized();
	// count required length of element along the given line 
	//	- according to metric at each node
	double len_m[] = {(cn1->control_data * dvn).length(), (cn2->control_data * dvn).length()};
	int vmax = (len_m[1]>len_m[0])?1:0;
	int vmin = 1-vmax;
	// check ratio of metric-lengths
	double a_m = (len_m[vmax]-len_m[vmin])/dlen;
	if(a_m > a_max)
		len_m[vmax] = (a_m=a_max) * dlen + len_m[vmin];
	// count max ratio:
	DMatrix2d e1, e2;
	double d1[2], d2[2];
	cn1->control_data.eigensystem(e1, d1);
	cn2->control_data.eigensystem(e2, d2);
	const ControlDataMatrix2d m1inv = cn1->control_data.inverse();
	const ControlDataMatrix2d m2inv = cn2->control_data.inverse();
	double real_max_ratio_1 = std::max( (m2inv * (e1.column(0) * d1[0])).length(), (m2inv * (e1.column(1) * d1[1])).length() );
	double real_max_ratio_2 = std::max( (m1inv * (e2.column(0) * d2[0])).length(), (m1inv * (e2.column(1) * d2[1])).length() );
	// count max acceptable diff
	double a1 = len_m[vmin]/(1.0-0.5*a_m);
	double max_diff_ratio = 1.0-(1.0-dr)*dlen/a1;
	int result = 0;
	double real_dr[2] = { dr, dr };

	if(real_max_ratio_1 > max_diff_ratio){
		if(cn1->control_data.setMinimum(cn2->control_data * max_diff_ratio)) 
			result += 1;
	}else real_dr[0] = 1.0 - a1*(1.0-real_max_ratio_1)/dlen;

	if(real_max_ratio_2 > max_diff_ratio){
		if(cn2->control_data.setMinimum(cn1->control_data * max_diff_ratio)) 
			result += 2;
	}else real_dr[1] = 1.0 - a1*(1.0-real_max_ratio_2)/dlen;

	if(real_dr[0] > cn1->max_gradation_ratio) 
		cn1->max_gradation_ratio = real_dr[0];
	if(real_dr[1] > cn2->max_gradation_ratio) 
		cn2->max_gradation_ratio = real_dr[1];

	return result;
}

void ControlSpace2dAdaptive::storeSizingEPS(const MeshContainer2d* mesh, const char* fname, int id, double vlen)
{
	ostringstream filename;
	filename << fname << "-size-" << id << ".eps";

	const DRect& box = getBoundingRect();
	EPSFile eps(filename.str(), box.x0, box.x1, box.y0, box.y1);

	const DVector2d du(vlen, 0.0);
	const DVector2d dv(0.0, vlen);

	// control nodes
	forEachControlNode([&](const ControlNode2d& cn) {
		const DVector2d vx = cn.control_data * du;
		const DVector2d vy = cn.control_data * dv;
		eps.drawLine(cn.coord - vx, cn.coord + vx);
		eps.drawLine(cn.coord - vy, cn.coord + vy);
	});

	// control bounding box
	eps.drawLine(box.getX0Y0(), box.getX1Y0());
	eps.drawLine(box.getX0Y1(),	box.getX1Y1());
	eps.drawLine(box.getX0Y0(), box.getX0Y1());
	eps.drawLine(box.getX1Y0(), box.getX1Y1());

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "StoreSizingEPS: bounding-box: " << box);

	// mesh boundary
	if(mesh){
		for(int i = 0; i < mesh->getPointsCount(); i++){
			const MeshPoint2d* mp = mesh->getPointAt(i);
			if(!mp->isBorder()) continue;
			for(int j = 0; j < mp->getRank(); j++){
				const MeshEdge2d* edge = mp->getEdge(j);
				if(edge->isBorder() && edge->getPointIndex(mp) == 0) // only border and only once
					eps.drawLine(mp->getCoordinates(), edge->getMeshPoint(1)->getCoordinates(), true);
			}
		}
	}
}

/// Refines control space basing on other control space, along given segment
bool ControlSpace2dAdaptive::setMinControlLine(Metric2dContext& mc, const CS2dPtr space, 
	const DPoint2d& pt0, const DPoint2d& pt1, 
	const DVector2d& dvn, double r)
{
	double t0 = 0.0;
	double dt = 0.01;
	DPoint2d xpt0 = pt0;
	const DVector2d dv = pt1 - pt0;
	bool last_time = false;
	bool any_change = false;
	while(true)
	{
		mc.countMetricAtPoint(xpt0);
		// calculate proper length of normal vector
		const DVector2d xdvn = dvn / mc.transformPStoMS(dvn).length();
		// insert new metric source
		any_change |= setMinControlPointVicinity(space, xpt0, xdvn, r);

		if(last_time) break;
		// calculate next point
		double t1 = t0 + dt;
		DPoint2d xpt1 = pt0 + dv * t1;
		// normalize for unit (metric) segment along selected line
		dt /= mc.transformPStoMS(xpt1 - xpt0).length();
		// calculate proper point
		t0 += 0.5*dt; // move by half of "unit" segment length
		if(t0 >= 1.0){ t0 = 1.0; last_time = true; }
		xpt0 = pt0 + dv * t0;
	}

	return any_change;
}

bool ControlSpace2dAdaptive::setMinControlPointVicinity(const CS2dPtr space, 
	const DPoint2d& pt, const DVector2d& dvn, double r, double dt)
{
	bool any_change = false;
	// insert new metric source
	// ... at point
	any_change |= setMinControl(pt, space->getMetricAtPoint(pt));
	// ... at width (in both directions)
	DPoint2d xpt = pt + dvn * r;
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

/// Smoothen variance of metric within the control space (with specified gradation ratio)
bool ControlSpace2dAdaptive::smoothenWithGradation(double gradation) 
{ 
	double tmp = param_gradation_ratio;   
	param_gradation_ratio = gradation;

	bool result = smoothen();

	param_gradation_ratio = tmp; 
	return result; 
}

void ControlSpace2dAdaptive::showMetricLength(const string& label) const
{
	static const int RES = 30;
	double dx = m_box.getDX() / (RES - 1);
	double dy = m_box.getDY() / (RES - 1);

	MeshViewSet* set = new MeshViewSet();

	DPoint2d pt;
	double t;

	double min_l = -1.0, max_l;

	for (int iy = 0; iy < RES; iy++) {
		t = (double)iy / (RES - 1);
		pt.y = m_box.y0 * (1.0 - t) + m_box.y1 * t;
		for (int ix = 0; ix < RES; ix++) {
			t = (double)ix / (RES - 1);
			pt.x = m_box.x0 * (1.0 - t) + m_box.x1 * t;

			auto cdm = getMetricAtPoint(pt);
			double mind = cdm.minEigenvalue();

			if (min_l < 0.0) min_l = max_l = mind;
			else {
				if (mind < min_l) min_l = mind;
				if (mind > max_l) max_l = mind;
			}
		}
	}

	for (int iy = 0; iy < RES; iy++) {
		t = (double)iy / (RES - 1);
		pt.y = m_box.y0 * (1.0 - t) + m_box.y1 * t;
		for (int ix = 0; ix < RES; ix++) {
			t = (double)ix / (RES - 1);
			pt.x = m_box.x0 * (1.0 - t) + m_box.x1 * t;

			auto cdm = getMetricAtPoint(pt);
			double mind = cdm.minEigenvalue();

			double q = (mind - min_l) / (max_l - min_l);
			set->addTetra(FPoint3d((float)pt.x, (float)pt.y, 0.0f), 
				(float)(dx * 0.5), (float)(dy * 0.5), (float)(dx * 0.5), q);
		}
	}

	set->setPolygonFillMode(MeshViewSet::FILL_QUALITY);
	set->addInfo("min_min_d", min_l);
	set->addInfo("max_min_d", max_l);
	SHOW_MESH( ((label!="") ? label : "metric min length"), set);
}