// MeshTetrahedron.cpp: implementation of the MeshTetrahedron class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshTetrahedron.h"
#include "MeshTriangle3d.h"
#include "MeshData.h"
#include "ControlSpace3d.h"
#include "MeshViewSet.h"
#include "Metric3dContext.h"
#include "MeshGenerator3d.h"
#include "MeshEdge3d.h"
#include "DTriangle.h"
#include "DTetrahedron.h"
#include "GeometricPredicates.h"

int MeshTetrahedron::param_quality_criterion = MeshData::QUALITY3D_CIRCLE_SPACE_AND_EDGES;
int MeshTetrahedron::param_quality_metric_selection = MeshData::QM_MIDDLE;
const double MeshTetrahedron::opt_metric_height = sqrt(6.0)/3.0;
unsigned int MeshTetrahedron::traverse_counter = 0;

/// Available sequences of vertices
const unsigned char MeshTetrahedron::PCONF[4][4] = {{1, 3, 2, 0}, {2, 3, 0, 1}, {0, 3, 1, 2}, {0, 1, 2, 3}};
/// Edge to vertices relation
const unsigned char MeshTetrahedron::ECONF[6][2] = {{1,3},{3,2},{2,1},{1,0},{3,0},{2,0}};

MeshTetrahedron::MeshTetrahedron(MeshPoint3d *p1, MeshPoint3d *p2, MeshPoint3d *p3, MeshPoint3d *p4) :
	MeshBlock(4, 4, face_array, point_array)
{
	points[0] = p1;
	points[1] = p2;
	points[2] = p3;
	points[3] = p4;

	double volume = DTriangle3d::orient3d(
		p1->getCoordinates(), p2->getCoordinates(), 
		p3->getCoordinates(), p4->getCoordinates());
	if(volume < 0.0){
		static int case_count = 0;
		++case_count;
		if (case_count <= 10) {
			LOG4CPLUS_DEBUG(MeshLog::logger_console, "Tetrahedron() -> orient3d: " << volume);
			if(case_count == 10)
				LOG4CPLUS_DEBUG(MeshLog::logger_console, "skipping (tracing) next negative orient3d cases...");
		}
		else{
			LOG4CPLUS_TRACE(MeshLog::logger_mesh, "Tetrahedron() -> orient3d: " << volume);
		}
	}

	for(int i = 0; i < 4; i++){
		faces[i] = points[PCONF[i][0]]->getFaceToPoints(points[PCONF[i][1]], points[PCONF[i][2]]);
		if(!faces[i]){	// Nowa krawêdŸ
			faces[i] = new MeshTriangle3d(points[PCONF[i][0]], points[PCONF[i][1]], points[PCONF[i][2]], this);
		}else{			
			faces[i]->setBlockLink(this, points[PCONF[i][0]], points[PCONF[i][1]]);
		}
	}
}

MeshTetrahedron::~MeshTetrahedron()
{
	if(faces[0]){
		// if the face is connected to this block only, it should be removed...
		if(faces[0]->clearBlockLink(this)) delete faces[0];
		if(faces[1]->clearBlockLink(this)) delete faces[1];
		if(faces[2]->clearBlockLink(this)) delete faces[2];
		if(faces[3]->clearBlockLink(this)) delete faces[3];
	}
}

/// clear some data for faster all-delete process
void MeshTetrahedron::preDeleteAll()
{
	faces[0] = 0; // block normal adjacency update
}

MeshTetrahedron* MeshTetrahedron::findTetrahedronByNeighbours(const DPoint3d &mpt, bool mind_border)
{
	MeshTetrahedron* tetrahedron = this;
	MeshTetrahedron* prev_tetrahedron = nullptr;

	traverse_counter = 0;
//	int second_counter = 0;

	while(true){
		// Is point inside the tetrahedron?
		if(!tetrahedron){
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "MeshTetrahedron:find empty...");
		 	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "FT: error-found-empty");
			return nullptr;
		}
		if(tetrahedron->isPointInside(mpt)) break;
		// If not ...
		if(++traverse_counter > 50000){
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"MeshTetrahedron: error findTetrahedronByNeighbours - 50000 steps exceeded.");
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "MeshTetrahedron:find error findTetrahedronByNeighbours - 50000 steps exceeded.");
		 	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "FT: error-50000-limit-reached");
			return nullptr;
		}
		// Select the tetrahedron in the direction of the search point
		MeshTetrahedron* t = tetrahedron->getNeighbourInFaceDirection(mpt, mind_border);
/*
		if(true){
			MeshViewSet* set = new MeshViewSet;
			set->addBlockWithEdges(tetrahedron, 0);
			if(t) set->addEmptyBlockWithEdges(t, 1);
			set->addPoint(mc.transformMStoRS(mpt), 2);
			SHOW_MESH("findTetrahedron", set);
		}
*/
		if(prev_tetrahedron == t){
			t = tetrahedron->getNeighbourInFaceDirection(mpt, mind_border, true);
//			++second_counter;
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "MeshTetrahedron:find second guess...");
		}
		if(prev_tetrahedron == t){
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "MeshTetrahedron:find loop...");
		 	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "FT: error-loop");
			return nullptr;
		}
		prev_tetrahedron = tetrahedron;
		tetrahedron = t;
	}
 	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "FT: " << counter);
/*
 	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "FindTetrahedron-F: " << counter;
 	if(second_counter) LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "\t" << second_counter;
 	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;
*/
	return tetrahedron;
}

/// Not used currently (worse then the FaceDirection version)
MeshTetrahedron* MeshTetrahedron::getNeighbourInBlockDirection(const DPoint3d &point, bool mind_border, bool second_best)
{
	int best = -1;
	int best2 = -1;
	double best_angle = -1.0;
	double best_angle2 = -1.0;
	DVector3d dv = (point - getMiddlePoint()).normalized();
	for(int i = 0; i < 4; i++){
		// Wyznacz cosinus k¹ta pomiêdzy wektorem normalnym do j-tej 
		//	œciany a wektorem kierunkowym
		DVector3d nv = faces[i]->getNormalVector();
		if(faces[i]->getBlockIndex(this) == MeshFace::BLOCK_UP) nv *= -1.0;
		double cos_angle = dv.scalarProduct(nv);
		if(cos_angle > best_angle){
			best2 = best;
			best = i;
			best_angle2 = best_angle;
			best_angle = cos_angle;
		}else if(cos_angle > best_angle2){
			best2 = i;
			best_angle2 = cos_angle;
		}
	}
	if(second_best) best = best2;

	if(mind_border && faces[best]->isBorder()) return nullptr;
	return (MeshTetrahedron*)(faces[best]->getOtherBlock(this));
}

MeshTetrahedron* MeshTetrahedron::getNeighbourInFaceDirection(const DPoint3d &mpt, bool mind_border, bool second_best)
{
	int best = -1;
	int best2 = -1;
	double best_angle = -1.0;
	double best_angle2 = -1.0;
	for(int i = 0; i < 4; i++){
		// Wyznacz cosinus k¹ta pomiêdzy wektorem normalnym do j-tej 
		//	œciany a wektorem kierunkowym
		const DPoint3d& dpt0 = faces[i]->getPoint(0)->getCoordinates();
		const DPoint3d& dpt1 = faces[i]->getPoint(1)->getCoordinates();
		const DPoint3d& dpt2 = faces[i]->getPoint(2)->getCoordinates();

		const DPoint3d dmiddle = DPoint3d::average(dpt0, dpt1, dpt2);

		DVector3d dv = (mpt - dmiddle).normalized();
		DVector3d nv = (dpt1 - dpt0).crossProduct(dpt2 - dpt0).normalized();
		if(faces[i]->getBlockIndex(this) == MeshFace::BLOCK_UP) nv *= -1.0;
		double cos_angle = dv.scalarProduct(nv);
		if(cos_angle > best_angle){
			best2 = best;
			best = i;
			best_angle2 = best_angle;
			best_angle = cos_angle;
		}else if(cos_angle > best_angle2){
			best2 = i;
			best_angle2 = cos_angle;
		}
	}
	if(second_best) best = best2;

	if(mind_border && faces[best]->isBorder()) return nullptr;
	return (MeshTetrahedron*)(faces[best]->getOtherBlock(this));
}

bool MeshTetrahedron::isPointInside(const DPoint3d &pt0) const
{
	const DPoint3d& pt1 = points[0]->getCoordinates();
	const DPoint3d& pt2 = points[1]->getCoordinates();
	const DPoint3d& pt3 = points[2]->getCoordinates();
	const DPoint3d& pt4 = points[3]->getCoordinates();

	return DTriangle3d::orient3d(pt1, pt2, pt3, pt0) >= 0.0 &&
		DTriangle3d::orient3d(pt1, pt4, pt2, pt0) >= 0.0 &&
		DTriangle3d::orient3d(pt2, pt4, pt3, pt0) >= 0.0 &&
		DTriangle3d::orient3d(pt3, pt4, pt1, pt0) >= 0.0;
}

bool MeshTetrahedron::isPointInOuterSphere(Metric3dContext& mc, const DMPoint3d& mpt, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint3d pt1 = points[0]->getMetricCoordinates(mc);
	const DMPoint3d pt2 = points[1]->getMetricCoordinates(mc);
	const DMPoint3d pt3 = points[2]->getMetricCoordinates(mc);
	const DMPoint3d pt4 = points[3]->getMetricCoordinates(mc);

	return DTetrahedron::insphere(mpt, pt1, pt2, pt3, pt4) >= 0.0;
}

//////////////////////////////////////////////////////////////////////
// Wylicza jakoœæ czworoœcianu wg. aktualnie przyjêtego kryterium
double MeshTetrahedron::countQuality(Metric3dContext& mc, int criterion)
{
	mc.countMetricAtPoints(points);
	const DMPoint3d p1 = points[0]->getMetricCoordinates(mc);
	const DMPoint3d p2 = points[1]->getMetricCoordinates(mc);
	const DMPoint3d p3 = points[2]->getMetricCoordinates(mc);
	const DMPoint3d p4 = points[3]->getMetricCoordinates(mc);

	int quality_criterion = criterion;
	if(criterion == -1) quality_criterion = param_quality_criterion;

	switch(quality_criterion){
	case MeshData::QUALITY3D_SPACE:
		{
			double volume = DTetrahedron::volume(p1, p2, p3, p4);
			double wanted_volume = getOptimumVolume();
			return quality = wanted_volume / volume;
		}
	case MeshData::QUALITY3D_CIRCLE_SPACE:
		{
			double wanted_radius = getOptimumCircumradius2();
			double r2 = DTetrahedron::outerSphereRadius2(p1, p2, p3, p4);
			return quality = wanted_radius / r2;
		}
	case MeshData::QUALITY3D_CIRCLE_SPACE_AND_EDGES:
		{
			double wanted_radius = getOptimumCircumradius2();
			double r2 = DTetrahedron::outerSphereRadius2(p1, p2, p3, p4);
			quality = wanted_radius / r2;
			if(quality > MeshGenerator3d::param_quality_threshold){
				// additionally check metric length of edges
				if( getEdge(0)->getLength(mc, true) > 1.5 ||
					getEdge(1)->getLength(mc, true) > 1.5 || 
					getEdge(2)->getLength(mc, true) > 1.5 || 
					getEdge(3)->getLength(mc, true) > 1.5 || 
					getEdge(4)->getLength(mc, true) > 1.5 || 
					getEdge(5)->getLength(mc, true) > 1.5)
				{	// needs refining
					quality = 0.95 * MeshGenerator3d::param_quality_threshold;
				}
			}
			return quality;
		}
/*	case MeshData::QUALITY3D_ASPECT_RATIO:
		{
			const double a = 6.0 * sqrt(6.0);	// scaling factor
			double h1 = p1.distance2(p2);
			double h2 = p2.distance2(p3);
			double h3 = p3.distance2(p1);
			double h4 = p1.distance2(p4);
			double h5 = p2.distance2(p4);
			double h6 = p3.distance2(p4);
			double hM = sqrt(std::max(std::max(h1,h2),std::max(std::max(h3,h4), std::max(h5,h6))));
			double S = p1.triangleArea(p2, p3) + p1.triangleArea(p2, p4) + 
				p2.triangleArea(p3, p4) + p3.triangleArea(p1, p4);

			return quality = a * volume / (hM * S);
		}
	case MeshData::QUALITY3D_RADIUS_RATIO:
		{
			double r_outer = DMPoint3d::countOuterSphereRadius(p1, p2, p3, p4);
			double S = p1.triangleArea(p2, p3) + p1.triangleArea(p2, p4) + 
				p2.triangleArea(p3, p4) + p3.triangleArea(p1, p4);
			double r_inner = 3.0 * volume / S;

			return quality = 3.0 * r_inner / r_outer;
		}
*/
	case MeshData::QUALITY3D_MEAN_RATIO:
		{
			double L = p1.distance2(p2) + p2.distance2(p3) + p3.distance2(p1) +
				p1.distance2(p4) + p2.distance2(p4) + p3.distance2(p4);
			double volume = DTetrahedron::volume(p1, p2, p3, p4);
			return quality = 12.0 * pow(9.0 * volume * volume, 1.0/3.0) / L;
		}
	}

	LOG4CPLUS_WARN(MeshLog::logger_console, 
		"Unknown quality_criterion in MeshTetrahedron::countQuality()");
	return 0.0;
}

MeshPoint3d* MeshTetrahedron::getOtherPoint(const MeshPoint3d *pt0, const MeshPoint3d *pt1, const MeshPoint3d *pt2) const
{
	for(int i = 0; i < 4; i++)
		if(points[i] != pt0 && points[i] != pt1 && points[i] != pt2)
			return points[i];
	return nullptr;
}

bool MeshTetrahedron::getOtherPoints(const MeshPoint3d* pt0, const MeshPoint3d* pt1, MeshPoint3d* pts[]) const
{
	int j = 0;
	for(int i = 0; i < 4; i++){
		if(points[i] != pt0 && points[i] != pt1){
			pts[j++] = points[i];
			if(j == 2) return true;
		}
	}
	assert(false);
	return false;
}


double MeshTetrahedron::getOuterSphereRadius(Metric3dContext& mc, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint3d pt1 = points[0]->getMetricCoordinates(mc);
	const DMPoint3d pt2 = points[1]->getMetricCoordinates(mc);
	const DMPoint3d pt3 = points[2]->getMetricCoordinates(mc);
	const DMPoint3d pt4 = points[3]->getMetricCoordinates(mc);

	return DTetrahedron::outerSphereRadius(pt1, pt2, pt3, pt4);
}

DMPoint3d MeshTetrahedron::getOuterSphereCenter(Metric3dContext& mc, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint3d pt1 = points[0]->getMetricCoordinates(mc);
	const DMPoint3d pt2 = points[1]->getMetricCoordinates(mc);
	const DMPoint3d pt3 = points[2]->getMetricCoordinates(mc);
	const DMPoint3d pt4 = points[3]->getMetricCoordinates(mc);

	return DTetrahedron::outerSphereCenter(pt1, pt2, pt3, pt4);
}

void MeshTetrahedron::getFaceIndices(int i, int& i1, int& i2, int& i3) const {
	i1 = PCONF[i][0];
	i2 = PCONF[i][1];
	i3 = PCONF[i][2];
}

/// Returns face opposite the given point
MeshFace* MeshTetrahedron::getOppositeFace(const MeshPoint3d* mp) const
{
	for(int i = 0; i < 4; i++)
		if(points[i] == mp) return faces[i];
	assert(false);
	return nullptr;
}

/// Returns point opposite the given face
MeshPoint3d* MeshTetrahedron::getOppositePoint(const MeshFace* mf) const
{
	for(int i = 0; i < 4; i++)
		if(faces[i] == mf) return points[i];
	assert(false);
	return nullptr;
}

void MeshTetrahedron::getEdgeIndices(int i, int& i1, int& i2) const {
	i1 = ECONF[i][0];
	i2 = ECONF[i][1];
}

MeshEdge3d* MeshTetrahedron::getEdge(int i) const
{
	return points[ECONF[i][0]]->getEdgeToPoint(points[ECONF[i][1]]);
}

bool MeshTetrahedron::isInverted() const
{
	double det_adaptive = GeometricPredicates::orient3d(
		points[0]->getCoordinates(), 
		points[1]->getCoordinates(), 
		points[2]->getCoordinates(), 
		points[3]->getCoordinates());
//#ifdef _DEBUG
//	double det_simple = DTriangle3d::orient3d(
//		points[0]->getCoordinates(), 
//		points[1]->getCoordinates(), 
//		points[2]->getCoordinates(), 
//		points[3]->getCoordinates());
//	if(det_simple * det_adaptive <= 0.0){
//		LOG4CPLUS_WARN(MeshLog::logger_console, ("Tetrahedron::inverted - simple", det_simple);
//		LOG4CPLUS_WARN(MeshLog::logger_console, ("Tetrahedron::inverted - adaptive", det_adaptive);
//	}
//#endif // _DEBUG
	return det_adaptive <= 0.0;
}

std::shared_ptr<MeshViewBlockData> MeshTetrahedron::getViewData(double shrink_ratio) const
{
	auto data = std::make_shared<MeshViewBlockData>(4, 4);
	data->area_id = area_id;
	data->quality = quality;
	const DPoint3d middle = getMiddlePoint();

	for(int i = 0; i < 4; i++)
		data->pts.add( middle + (points[i]->getCoordinates() - middle) * shrink_ratio );

	//------------
	for(int i = 0; i < 12; i++)
		data->indices.add( PCONF[i/3][i%3] );

	for(int i = 0; i < 4; i++)
		data->normals.add( FVector3d::crossProduct(data->pts[PCONF[i][0]], 
			data->pts[PCONF[i][1]], data->pts[PCONF[i][2]]).normalized() );

//	data->adjacent[0] = data->adjacent[1] = data->adjacent[2] = data->adjacent[3] = nullptr;

	return data;
}

double MeshTetrahedron::getVolume(Metric3dContext& mc, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint3d pt1 = getPoint(0)->getMetricCoordinates(mc);
	const DMPoint3d pt2 = getPoint(1)->getMetricCoordinates(mc);
	const DMPoint3d pt3 = getPoint(2)->getMetricCoordinates(mc);
	const DMPoint3d pt4 = getPoint(3)->getMetricCoordinates(mc);

	return DTetrahedron::volume(pt1, pt2, pt3, pt4);
}

double MeshTetrahedron::getVolumeNoMetric() const
{
	const DPoint3d& pt1 = getPoint(0)->getCoordinates();
	const DPoint3d& pt2 = getPoint(1)->getCoordinates();
	const DPoint3d& pt3 = getPoint(2)->getCoordinates();
	const DPoint3d& pt4 = getPoint(3)->getCoordinates();

	return GeometricPredicates::tetrahedra_volume(pt1, pt2, pt3, pt4);
}

/// Calculates the minimum dihedral angle
double MeshTetrahedron::getMinDihedralAngleSin() const
{
	if(isInverted()) return 0.0;

	DVector3d dns[4];
	for(int i = 0; i < 4; i++){
		dns[i] = DVector3d::crossProduct(
			points[PCONF[i][0]]->getCoordinates(), 
			points[PCONF[i][1]]->getCoordinates(), 
			points[PCONF[i][2]]->getCoordinates()).normalized();
		dns[i].reverse();
	}
	double min_angle_sc = 1;  // 1 -> 180o, 0 -> 90o, -1 -> 0o
	for(int i = 0; i < 3; i++){
		for(int j = i+1; j < 4; j++){
			double angle_sc = dns[i].scalarProduct(dns[j]);
			if(angle_sc < min_angle_sc) min_angle_sc = angle_sc;
		}
	}

	return 1.0 + min_angle_sc;
}

/// Returns "mean ratio" coefficient quality
double MeshTetrahedron::getMeanRatio(Metric3dContext& mc, bool ext_metric) const
{
	int qms = MeshTetrahedron::param_quality_metric_selection;
	MeshTetrahedron::param_quality_metric_selection = 
		ext_metric ? MeshData::QM_VERTICES_AVE : MeshData::QM_MIDDLE;
	mc.countMetricAtPoints(points);
	MeshTetrahedron::param_quality_metric_selection = qms;

	const DMPoint3d p1 = getPoint(0)->getMetricCoordinates(mc);
	const DMPoint3d p2 = getPoint(1)->getMetricCoordinates(mc);
	const DMPoint3d p3 = getPoint(2)->getMetricCoordinates(mc);
	const DMPoint3d p4 = getPoint(3)->getMetricCoordinates(mc);

	double vol = DTetrahedron::volume(p1, p2, p3, p4);

	return sgn(vol) * 12.0 * pow(9.0 * sqr(vol), (1.0/3.0)) / 
		(p1.distance2(p2) + p1.distance2(p3) + p1.distance2(p4) + 
		 p2.distance2(p3) + p2.distance2(p4) + p3.distance2(p4));
}

double MeshTetrahedron::countMetricDiff(Metric3dContext& mc) const
{
	ControlDataMatrix3d cdmm_simplex = 
		ControlDataMatrix3d::countMetricTensor(
			points[0]->getCoordinates(), points[1]->getCoordinates(), 
			points[2]->getCoordinates(), points[3]->getCoordinates());

	ControlDataMatrix3d cdmm_space = 
		mc.getControlSpace()->getMetricAtPoint(getMiddlePoint()).transformationToTensor();
	return cdmm_simplex.countDifferenceRR(cdmm_space);
}


double MeshTetrahedron::countMetricDiffQuality(Metric3dContext& mc)
{
	return quality = (1.0 - 0.1 * countMetricDiff(mc));
}

// half (12) of all permutations (i.e. even permutations of 0-1-2-3)
bool MeshTetrahedron::properOrientation(const MeshPoint3d* mpt0, 
		const MeshPoint3d* mpt1, const MeshPoint3d* mpt2) const
{
	if(mpt0 == points[0]){
		return 	(mpt1 == points[1] && mpt2 == points[2]) ||
				(mpt1 == points[2] && mpt2 == points[3]) ||
				(mpt1 == points[3] && mpt2 == points[1]);
	}else if(mpt0 == points[1]){
		return 	(mpt1 == points[2] && mpt2 == points[0]) ||
				(mpt1 == points[0] && mpt2 == points[3]) ||
				(mpt1 == points[3] && mpt2 == points[2]);
	}else if(mpt0 == points[2]){
		return 	(mpt1 == points[0] && mpt2 == points[1]) ||
				(mpt1 == points[1] && mpt2 == points[3]) ||
				(mpt1 == points[3] && mpt2 == points[0]);
	}else if(mpt0 == points[3]){
		return 	(mpt1 == points[0] && mpt2 == points[2]) ||
				(mpt1 == points[2] && mpt2 == points[1]) ||
				(mpt1 == points[1] && mpt2 == points[0]);
	}else return false;
}

/// add this element as a lump of mass (respecting control space)
void MeshTetrahedron::addForInertialCenter(Metric3dContext& mc, DPoint3d& inertial_center, double& total_mass) const
{
	addForInertialCenterAdaptiveSplitLongestEdge(mc, 
		getPoint(0)->getCoordinates(), 
		getPoint(1)->getCoordinates(),
		getPoint(2)->getCoordinates(),
		getPoint(3)->getCoordinates(),
		inertial_center, total_mass);
}

/// add this element as a lump of mass (respecting control space)
void MeshTetrahedron::addForInertialCenterAdaptiveSplitLongestEdge(Metric3dContext& mc, 
	const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, const DPoint3d& pt3,
	DPoint3d& inertial_center, double& total_mass, int level)
{
	const DPoint3d middle = DPoint3d::average(pt0, pt1, pt2, pt3);
	mc.countMetricAtPoint(middle);
	const DMPoint3d mpt[4] = {
		mc.transformRStoMS(pt0),	mc.transformRStoMS(pt1),
		mc.transformRStoMS(pt2),	mc.transformRStoMS(pt3) };

	// check recursion depth
	bool last_step = (level <= 10);

	double len[6];
	int longest_i = 0;

	if(!last_step){ // check max length
		// calculate lengths
		for(int i = 0; i < 6; i++)
			len[i] = mpt[ECONF[i][0]].distance2(mpt[ECONF[i][1]]);
		// find longest
		for(int i = 1; i < 6; i++)
			if(len[i] > len[longest_i])
				longest_i = i;
		// check threshold
		last_step = (len[longest_i] < 4.0);  // len^2 !!!
	}

	// calculate
	if(last_step){
		// ====> min metric
		double volume = DTetrahedron::volume(mpt[0], mpt[1], mpt[2], mpt[3]);
		total_mass += volume;
		inertial_center.add(middle, volume);
	}else{
		// {{1,3},{3,2},{2,1},{1,0},{3,0},{2,0}};
		if(longest_i == 0){
			const DPoint3d pt = DPoint3d::average(pt1, pt3);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt0, pt1, pt2, pt, 
				inertial_center, total_mass, level+1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt, pt2, pt0, pt3,
				inertial_center, total_mass, level+1);
		}else if(longest_i == 1){
			const DPoint3d pt = DPoint3d::average(pt2, pt3);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt0, pt1, pt2, pt, 
				inertial_center, total_mass, level+1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt, pt0, pt1, pt3,
				inertial_center, total_mass, level+1);
		}else if(longest_i == 2){
			const DPoint3d pt = DPoint3d::average(pt1, pt2);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt0, pt1, pt, pt3, 
				inertial_center, total_mass, level+1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt, pt2, pt0, pt3,
				inertial_center, total_mass, level+1);
		}else if(longest_i == 3){
			const DPoint3d pt = DPoint3d::average(pt0, pt1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt2, pt0, pt, pt3, 
				inertial_center, total_mass, level+1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt, pt1, pt2, pt3,
				inertial_center, total_mass, level+1);
		}else if(longest_i == 4){
			const DPoint3d pt = DPoint3d::average(pt0, pt3);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt0, pt1, pt2, pt, 
				inertial_center, total_mass, level+1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt, pt1, pt2, pt3,
				inertial_center, total_mass, level+1);
		}else if(longest_i == 5){
			const DPoint3d pt = DPoint3d::average(pt0, pt2);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt1, pt2, pt, pt3, 
				inertial_center, total_mass, level+1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt, pt0, pt1, pt3,
				inertial_center, total_mass, level+1);
		}
	}
}

/// add this element as a lump of mass (respecting control space)
void MeshTetrahedron::addForInertialMoments(Metric3dContext& mc, const DPoint3d& inertial_center, DMatrix3d& inertial_moments) const
{
	addForInertialMomentsAdaptiveSplitLongestEdge(mc, 
		getPoint(0)->getCoordinates(), 
		getPoint(1)->getCoordinates(),
		getPoint(2)->getCoordinates(),
		getPoint(3)->getCoordinates(),
		inertial_center, inertial_moments);
}

/// add this element as a lump of mass (respecting control space)
void MeshTetrahedron::addForInertialMomentsAdaptiveSplitLongestEdge(Metric3dContext& mc, 
	const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, const DPoint3d& pt3, 
	const DPoint3d& inertial_center, DMatrix3d& inertial_moments, int level)
{
	const DPoint3d middle = DPoint3d::average(pt0, pt1, pt2, pt3);
	mc.countMetricAtPoint(middle);
	const DMPoint3d mpt[4] = {
		mc.transformRStoMS(pt0),	mc.transformRStoMS(pt1),
		mc.transformRStoMS(pt2),	mc.transformRStoMS(pt3) };

	// check recursion depth
	bool last_step = (level <= 10);

	double len[6];
	int longest_i = 0;

	if(!last_step){ // check max length
		// calculate lengths
		for(int i = 0; i < 6; i++)
			len[i] = mpt[ECONF[i][0]].distance2(mpt[ECONF[i][1]]);
		// find longest
		for(int i = 1; i < 6; i++)
			if(len[i] > len[longest_i]) 
				longest_i = i;
		// check threshold
		last_step = (len[longest_i] < 4.0);  // len^2 !!!
	}

	// calculate
	if(last_step){
		// ====> min metric
		double volume = DTetrahedron::volume(mpt[0], mpt[1], mpt[2], mpt[3]);
		const DVector3d dv = middle - inertial_center;
		inertial_moments.m[0][0] += volume * (dv.y*dv.y + dv.z*dv.z); // Ixx
		inertial_moments.m[1][1] += volume * (dv.x*dv.x + dv.z*dv.z); // Iyy
		inertial_moments.m[2][2] += volume * (dv.y*dv.y + dv.x*dv.x); // Izz
		inertial_moments.m[0][1] += -volume * (dv.x*dv.y); // Ixy
		inertial_moments.m[1][0] += -volume * (dv.x*dv.y); // Iyx = Ixy
		inertial_moments.m[0][2] += -volume * (dv.x*dv.z); // Ixz
		inertial_moments.m[2][0] += -volume * (dv.x*dv.z); // Izx = Ixz
		inertial_moments.m[1][2] += -volume * (dv.y*dv.z); // Iyz
		inertial_moments.m[2][1] += -volume * (dv.y*dv.z); // Izy = Iyz
	}else{
		// {{1,3},{3,2},{2,1},{1,0},{3,0},{2,0}};
		if(longest_i == 0){
			const DPoint3d pt = DPoint3d::average(pt1, pt3);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt0, pt1, pt2, pt, 
				inertial_center, inertial_moments, level+1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt, pt2, pt0, pt3,
				inertial_center, inertial_moments, level+1);
		}else if(longest_i == 1){
			const DPoint3d pt = DPoint3d::average(pt2, pt3);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt0, pt1, pt2, pt, 
				inertial_center, inertial_moments, level+1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt, pt0, pt1, pt3,
				inertial_center, inertial_moments, level+1);
		}else if(longest_i == 2){
			const DPoint3d pt = DPoint3d::average(pt1, pt2);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt0, pt1, pt, pt3, 
				inertial_center, inertial_moments, level+1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt, pt2, pt0, pt3,
				inertial_center, inertial_moments, level+1);
		}else if(longest_i == 3){
			const DPoint3d pt = DPoint3d::average(pt0, pt1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt2, pt0, pt, pt3, 
				inertial_center, inertial_moments, level+1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt, pt1, pt2, pt3,
				inertial_center, inertial_moments, level+1);
		}else if(longest_i == 4){
			const DPoint3d pt = DPoint3d::average(pt0, pt3);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt0, pt1, pt2, pt, 
				inertial_center, inertial_moments, level+1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt, pt1, pt2, pt3,
				inertial_center, inertial_moments, level+1);
		}else if(longest_i == 5){
			const DPoint3d pt = DPoint3d::average(pt0, pt2);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt1, pt2, pt, pt3, 
				inertial_center, inertial_moments, level+1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt, pt0, pt1, pt3,
				inertial_center, inertial_moments, level+1);
		}
	}
}

bool MeshTetrahedron::checkSwitchPointsWithFaces(const MeshPoint3d *point1, MeshPoint3d *point2) const
{
	int pi = -1;
	while(++pi < 4) if(points[pi] == point1) break;
	if(pi >= 4) return false; // doesn't contain the first point

	MeshPoint3d* xpoints[4] = {
		points[0], points[1], points[2], points[3] };
	xpoints[pi] = point2;

	for(int i = 0; i < 4; i++){
		assert(faces[i]->incidentToBlock(this));
		// if affected
		if(!faces[i]->incidentToPoint(point1)) continue;
		// check new face
		MeshFace* f = xpoints[PCONF[i][0]]->getFaceToPoints(xpoints[PCONF[i][1]], xpoints[PCONF[i][2]]);
		if(f){
			if(f->isBoundedBothSides()) return false; // can't have third block incident to a face ...
			if(f->incidentToBlock(this)) return false; // already incident ...
		}
	}

	return true;
}

void MeshTetrahedron::switchPointsWithFaces(const MeshPoint3d *point1, MeshPoint3d *point2)
{
	int pi = -1;
	while(++pi < 4) if(points[pi] == point1) break;
	assert(pi < 4);

	points[pi] = point2;

	for(int i = 0; i < 4; i++){
		assert(faces[i]->incidentToBlock(this));
		// if affected
		if(!faces[i]->incidentToPoint(point1)) continue;
		// remove old reference
		if(faces[i]->removeBlockLink(this)) delete faces[i];
		// create new reference (and face if necessary)
		faces[i] = points[PCONF[i][0]]->getFaceToPoints(points[PCONF[i][1]], points[PCONF[i][2]]);
		if(!faces[i]){	// new one
			faces[i] = new MeshTriangle3d(points[PCONF[i][0]], points[PCONF[i][1]], points[PCONF[i][2]], this);
		}else{
			faces[i]->setBlockLink(this, points[PCONF[i][0]], points[PCONF[i][1]]);
		}
	}
}

void MeshTetrahedron::switchPointsWithFacesBoundary(const MeshPoint3d *point1, MeshPoint3d *point2)
{
	int pi = -1;
	while(++pi < 4) if(points[pi] == point1) break;
	assert(pi < 4);

	points[pi] = point2;

	for(int i = 0; i < 4; i++){
		// if affected
		if(!faces[i]->incidentToPoint(point1)) continue;
		// remove old reference
		MeshFace* old_face = faces[i];
		// create new reference (and face if necessary)
		faces[i] = points[PCONF[i][0]]->getFaceToPoints(points[PCONF[i][1]], points[PCONF[i][2]]);
		if(!faces[i]){	// new one
			if(old_face->isBoundedBothSides()){
				old_face->removeBlockLink(this);
				faces[i] = new MeshTriangle3d(points[PCONF[i][0]], points[PCONF[i][1]], points[PCONF[i][2]], this);
				if(old_face->isBorder()){
					faces[i]->setBorder(old_face->getBorderFlags());
					faces[i]->copyAllTags(old_face);
				}
				for(int j = 0; j < 3; j++){
					MeshEdge3d* old_edge = old_face->getEdge(j);
					if(old_edge->isBorder() && old_edge->incidentTo(point1)){
						MeshEdge3d* new_edge = old_edge->getOtherPoint(point1)->getEdgeToPoint(point2);
						assert(new_edge);
						new_edge->setBorderFlags(old_edge->getBorderFlags());
					}
				}
			}else{
				faces[i] = old_face;
				faces[i]->switchPointsWithEdgesBoundary(point1, point2);
			}
		}else{
			if(old_face->removeBlockLinkDisregardBoundary(this)){
				old_face->detachFromEdgesDisregardBoundary();
				delete old_face;
			}
			faces[i]->setBlockLink(this, points[PCONF[i][0]], points[PCONF[i][1]]);
		}
	}
}

MeshBlock * MeshTetrahedron::clone() const
{
	MeshTetrahedron* block = new MeshTetrahedron(points[0], points[1], points[2], points[3]);
	block->copyAllTags(this);
	return block;
}

/// Create face-edges connections
void MeshTetrahedron::attachToFaces()
{
	for(int i = 0; i < 4; i++){
		faces[i] = points[PCONF[i][0]]->getFaceToPoints(points[PCONF[i][1]], points[PCONF[i][2]]);
		if(!faces[i]){	// New face
			faces[i] = new MeshTriangle3d(points[PCONF[i][0]], points[PCONF[i][1]], points[PCONF[i][2]], this);
		}else{			
			faces[i]->setBlockLink(this, points[PCONF[i][0]], points[PCONF[i][1]]);
		}
	}
}

/// Removes face-edges connections
void MeshTetrahedron::detachFromFaces()
{
	// if the face is connected to this block only, it should be removed...
	for(int i = 0; i < 4; i++){
		if(faces[i]){
			if(faces[i]->clearBlockLink(this)) 
				delete faces[i];
			faces[i] = nullptr;
		}
	}
}
