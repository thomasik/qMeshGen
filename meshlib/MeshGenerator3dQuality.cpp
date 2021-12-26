// MeshGenerator3dQuality.cpp: implementation of the MeshGenerator3dQuality class.
//
//////////////////////////////////////////////////////////////////////

//#define USE_OPENMP_HERE

#include <iomanip>
#include <algorithm>

#include "MeshGenerator3dQuality.h"

#include "MeshContainer3d.h"
#include "MeshDomainVolume.h"
#include "OctTree.h"
#include "MeshTetrahedron.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshFace.h"
#include "MeshBlock.h"
#include "DataVector.h"
#include "MeshData.h"
#include "MeshLog.h"
#include "MeshBoundaryCondition.h"
#include "ControlSpace3d.h"
#include "MeshGenerator3d.h"
#include "MeshGenerator3dAdapt.h"
#include "MeshGenerator3dDelaunayBoundary.h"
#include "DataStatistics.h"
#include "DPoint.h"
#include "MeshViewSet.h"
#include "DataStatistics.h"
#include "GeometricPredicates.h"

int MeshGenerator3dQuality::param_swap_criterion = MeshData::SWAP3_MIN_VOLUME;

void MeshGenerator3dQuality::testTetrahedralMesh(const MeshContainer3d* mesh)
{
	assert(mesh);
	int bct = mesh->getBlocksCount();
	for(int i = 0; i < bct; i++){
		assert(mesh->getBlockAt(i));
		assert(mesh->getBlockAt(i)->getType() == BLOCK_TETRA);
		// .. some additional tests if required
	}
}

bool MeshGenerator3dQuality::tryMovingPoint(Metric3dContext& mc, 
		MeshPoint3d *point, const DPoint3d& new_pt,
		bool no_boundary)
{
	LOG_ASSERT(!no_boundary || !point->isBorder());
	// MeshingException("MG3:tryMovingPoint moving boundary-node not allowed"), false);

	const DPoint3d old_pt = point->getCoordinates();
	const DVector3d vector = new_pt - old_pt;
	int rank = point->getRank();

	// gather blocks incident to this vertex
	DataVector<MeshBlock*> blocks(2*rank);
	point->adjacentBlocks(blocks);
	int block_count = blocks.countInt();

	double factor = 1.0;
	for(int k = 0; k < 10; k++){	// 10 tries
		point->setCoordinates(old_pt + vector*factor);
		// Check for validity
		bool fault = false;
		for(int j = 0; j < block_count; j++){
			const MeshBlock* block = blocks.get(j);
			fault = (block->getVolume(mc, false) < METRIC_SMALL_NUMBER);
			//fault = block->isInverted();
			if(fault) break;
		}
		if(!fault){
//			LOG4CPLUS_INFO(MeshLog::logger_mesh, "tryMovingPoint: success for factor=" << factor);
			return true;
		}
		else factor /= 2.0;
	}
	// 10x fault -> cancel movement
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "tryMovingPoint: failure for end-factor=" << factor);
	point->setCoordinates(old_pt);
	return false;
}

bool MeshGenerator3dQuality::movePointByLaplace(Metric3dContext& mc, MeshPoint3d *point)
{
	if(point->isBorder()) return false;
#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
	point->setActive();
#endif

	DPoint3d new_pt, old_pt = point->getCoordinates();
	int rank = point->getRank();
	double weight = 0.0;
	for(int j = 0; j < rank; j++){
		MeshEdge3d* edge = point->getEdge(j);
		MeshPoint3d* other_point = edge->getOtherPoint(point);
		#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
		if(other_point->isActive()){
			point->clearActive();
			return false;
		}
		#endif
		new_pt.add(other_point->getCoordinates());// * w;
		weight += 1.0; //w;
	}
	new_pt /= weight;

	bool res = tryMovingPoint(mc, point, new_pt);
#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
	point->clearActive();
#endif
	return res;
}

bool MeshGenerator3dQuality::moveBoundaryPointByLaplace(Metric3dContext& mc, MeshPoint3d *point)
{
	

	if(!point->isBorder()) return false;
	if(point->isBorder(TagBorder::CORNER | TagBorder::FIXED)) return false; // corner and/or fixed points are immovable
	bool is_ridge = point->isBorder(TagBorder::RIDGE); // ridge-point moved only with respect to other ridge-points

#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
	point->setActive();
#endif

	DPoint3d new_pt, old_pt = point->getCoordinates();
	int rank = point->getRank();
	double weight = 0.0;
	mc.countMetricAtPoint( old_pt );
	for(int j = 0; j < rank; j++){
		MeshEdge3d* edge = point->getEdge(j);
		if(!edge->isBorder()) continue; // select only boundary edges
		if(is_ridge && !edge->isBorder(TagBorder::RIDGE)) continue;
		MeshPoint3d* other_point = edge->getOtherPoint(point);
		#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
		if(other_point->isActive()){
			point->clearActive();
			return false;
		}
		#endif
		new_pt.add(other_point->getCoordinates());// * w;
		weight += 1.0; //w;
	}
	new_pt /= weight;

	bool res = false;
	assert(false); // need to be rewritten, as in Generator3dSurface::moveByLaplace...
	//bool res = tryMovingPoint(mc, point, point->getCoordinatesWithLocalShape(mc, new_pt), false);
#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
	point->clearActive();
#endif
	return res;
}

bool MeshGenerator3dQuality::smoothenBlocks(MeshContainer3d* boundary, int steps, 
		TagExtended::TagType tag_type, int tag_value, int method)
{
	START_CLOCK("MG3d::smoothenBlocks");
	int block_count = boundary->getBlocksCount();
	for(int i = 0; i < block_count; i++){
		MeshDomainVolume* domain_volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		assert(domain_volume && (domain_volume->getType() == BLOCK_DOMAIN));
		domain_volume->smoothenVolumeMesh(steps, tag_type, tag_value, method);
	}
	STOP_CLOCK("MG3d::smoothenBlocks");
	return true;
}

double MeshGenerator3dQuality::calculateCavityMinAngleSc(const DPoint3d & dpt, 
		DataVector<DVector3d> tri_fvec[3], DataVector<DPoint3d> tri_fpts[3], 
		const DataVector<DVector3d> & normals, DataVector<int> & min_i)
{
	int bct = normals.countInt();
	min_i.clear();
	min_i.add(0);

	double min_min_angle_sc = 2.0;
	for(int i = 0; i < bct; i++){
		double min_angle_sc = 2.0;
		if( GeometricPredicates::orient3d(tri_fpts[0][i], tri_fpts[1][i], tri_fpts[2][i], dpt) >= 0.0){
			min_angle_sc = 0.0;
		}else{
			for(int k = 0; k < 3; k++){ // calculate local min angle (per face)
				DVector3d fvn = (dpt - tri_fpts[k][i]).crossProduct(tri_fvec[k][i]).normalized();
				double angle_sc = normals[i].scalarProduct(fvn);
				if(angle_sc < min_angle_sc) min_angle_sc = angle_sc;
			}
		}
		if(i > 0){ // update active set and global min angle
			if( abs(min_angle_sc - min_min_angle_sc) < SMALL_NUMBER ) 
				min_i.add(i);
			else if(min_angle_sc < min_min_angle_sc){
				min_i.clear();
				min_i.add(i);
				min_min_angle_sc = min_angle_sc;
			}
		}else{
			min_min_angle_sc = min_angle_sc;
		}
	}
	return min_min_angle_sc;
}

/// Auxiliary function for calcuating optimum vector for optimization
DVector3d MeshGenerator3dQuality::calculateCavityOptVec(const DPoint3d & dpt,
	const DataVector<DPoint3d> & opt_vertices, const DataVector<int> & min_i)
{
	int mct = min_i.countInt();
	assert(mct > 0);
	DVector3d opt_vec;
	for(int i = 0; i < mct; i++){
		DVector3d vec = opt_vertices[min_i[i]] - dpt;
		opt_vec += vec;
	}
	if(mct > 1 ) opt_vec /= mct;

	return opt_vec;
}

/// Tries to improve the mesh locally by moving the given point using optimization-based smoothing
double MeshGenerator3dQuality::movePointOpt(Metric3dContext& mc, MeshPoint3d *point, double threshold)
{
	DataVector<MeshBlock*> blocks(point->getRank());
	if(!point->adjacentBlocks(blocks))
		return -1.0;

	DPoint3d dpt = point->getCoordinates();
	DPoint3d dpt_initial = dpt;
	mc.countMetricAtPoint(dpt);

	int bct = blocks.countInt();
	DataVector<DVector3d> normals(bct);
	DataVector<DPoint3d> opt_vertices(bct);
	DataVector<DVector3d> tri_fvec[3] = {
		DataVector<DVector3d>(bct), DataVector<DVector3d>(bct), DataVector<DVector3d>(bct) };
	DataVector<DPoint3d> tri_fpts[3] = {
		DataVector<DPoint3d>(bct), DataVector<DPoint3d>(bct), DataVector<DPoint3d>(bct) };
	DataVector<int> min_i(bct);

	if(threshold < 1.0){
		double real_min_angle_sc = 2.0;
		// calculate real min angle and check with threshold
		for(int i = 0; i < bct; i++){
			MeshTetrahedron* tetra = (MeshTetrahedron*)blocks[i];
			if(tetra->isInverted()){
				real_min_angle_sc = 0.0; break;
			}
			double a = tetra->getMinDihedralAngleSin();
			if(a < real_min_angle_sc) real_min_angle_sc = a;
		}
		if(real_min_angle_sc >= (threshold+1)) 
			return real_min_angle_sc;
	}

	for(int i = 0; i < bct; i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*)blocks[i];
		// calculate optimum position of this node for this tetrahedron
		MeshFace* face = tetra->getOppositeFace(point);
		DVector3d vn = face->getNormalVector();
		if(face->getBlockIndex(tetra) == 1){
			tri_fvec[0].add( (face->getPoint(1)->getCoordinates() - face->getPoint(0)->getCoordinates()).normalized() );
			tri_fvec[1].add( (face->getPoint(2)->getCoordinates() - face->getPoint(1)->getCoordinates()).normalized() );
			tri_fvec[2].add( (face->getPoint(0)->getCoordinates() - face->getPoint(2)->getCoordinates()).normalized() );
			tri_fpts[0].add( face->getPoint(1)->getCoordinates() );
			tri_fpts[1].add( face->getPoint(2)->getCoordinates() );
			tri_fpts[2].add( face->getPoint(0)->getCoordinates() );
		}else{
			vn.reverse();
			tri_fvec[0].add( (face->getPoint(2)->getCoordinates() - face->getPoint(0)->getCoordinates()).normalized() );
			tri_fvec[1].add( (face->getPoint(1)->getCoordinates() - face->getPoint(2)->getCoordinates()).normalized() );
			tri_fvec[2].add( (face->getPoint(0)->getCoordinates() - face->getPoint(1)->getCoordinates()).normalized() );
			tri_fpts[0].add( face->getPoint(2)->getCoordinates() );
			tri_fpts[1].add( face->getPoint(1)->getCoordinates() );
			tri_fpts[2].add( face->getPoint(0)->getCoordinates() );
		}

		assert( (dpt - face->getMiddlePoint()).scalarProduct(vn) > 0);
		normals.add(vn);
		double len = mc.transformRStoMS(vn).length();
		const DPoint3d opt_pt = face->getMiddlePoint() + vn / len;
		opt_vertices.add(opt_pt);
	}

	double min_min_angle_sc = MeshGenerator3dQuality::calculateCavityMinAngleSc(dpt, tri_fvec, tri_fpts, normals, min_i);
	if(min_min_angle_sc > 0.0){
		for(int i = 0; i < bct; i++){
			if(blocks[i]->isInverted()){
				min_min_angle_sc = 0.0; break;
			}
		}
	}

	//if(true){
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "==================================================");
	//	double real_min = 2.0;
	//	for(int i = 0; i < bct; i++){
	//		double min_a = ((MeshTetrahedron*)blocks[i])->getMinDihedralAngleSin();
	//		if(min_a < real_min) real_min = min_a;
	//	}
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "min_angle_sc = " << (1.0+min_min_angle_sc)
	//		<< ", min_angle_sc for whole tetrahedron = " << real_min << endl;
	//}

	// ... select step along the opt_vec, (1) until next change of active-set
	// ...                                (2) with direct sampling

	DVector3d opt_vec;
	for(int i_step = 0; i_step < 2; i_step++){
		DPoint3d dpt_oryg = dpt;
		opt_vec = MeshGenerator3dQuality::calculateCavityOptVec(dpt, opt_vertices, min_i);

		double dt = 0.1;
		DataVector<int> temp_min_i(bct);
		for(int i = 0; i < 3; i++){
			double best_t = 0.0;
			for(double t = dt; t <= 1.0; t += dt){
				DPoint3d t_dpt = dpt + opt_vec * t;
				double a = MeshGenerator3dQuality::calculateCavityMinAngleSc(t_dpt, tri_fvec, tri_fpts, normals, temp_min_i);
				if(a > 0.0){
					for(int j = 0; j < bct; j++){
						if(blocks[j]->isInverted()){
							a = 0.0; break;
						}
					}
				}
				if(a > min_min_angle_sc){
					min_min_angle_sc = a; 
					best_t = t;
				}
				if(!min_i.contains(temp_min_i)) // change of active-set
					break;
			}
			dpt += opt_vec * best_t;
			opt_vec *= dt;
		}

		//if(true){
		//	point->setCoordinates(dpt);
		//	double real_min = 2.0;
		//	for(int i = 0; i < bct; i++){
		//		double min_a = ((MeshTetrahedron*)blocks[i])->getMinDihedralAngleSin();
		//		if(min_a < real_min) real_min = min_a;
		//	}
		//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "min_angle_sc = " << (1.0+min_min_angle_sc)
		//		<< ", min_angle_sc for whole tetrahedron = " << real_min << endl;
		//}

		if(min_i.addAllIfNew(temp_min_i) == 0) break;
		if(dpt.distance2(dpt_oryg) < SMALL_NUMBER) break;
	}

	point->setCoordinates(dpt);

	//if(true){
	//	movePointUsingAngleOptimization(point);
	//	double real_min = 2.0;
	//	for(int i = 0; i < bct; i++){
	//		double min_a = ((MeshTetrahedron*)blocks[i])->getMinDihedralAngleSin();
	//		if(min_a < real_min) real_min = min_a;
	//	}
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "min_angle_sc for whole tetrahedron (+opt) = " << real_min);
	//}

	//if(true){
	//	point->setCoordinates(dpt_initial);
	//	movePointUsingAngleOptimization(point);
	//	double real_min = 2.0;
	//	for(int i = 0; i < bct; i++){
	//		double min_a = ((MeshTetrahedron*)blocks[i])->getMinDihedralAngleSin();
	//		if(min_a < real_min) real_min = min_a;
	//	}
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "min_angle_sc for whole tetrahedron (opt only) = " << real_min);
	//}

	//if(true){
	//	MeshViewSet* set = new MeshViewSet;
	//	set->addPoint( dpt_initial );
	//	set->addPoint( dpt );
	//	for(int i = 0; i < min_i.countInt(); i++){
	//		set->addEdge( dpt, opt_vertices[min_i[i]] );
	//	}
	//	set->addEdge(dpt, dpt+opt_vec, 1);
	//	set->addPoint( point, 2 );
	//	SHOW_MESH("point optimization", set);
	//}

	// int i_step = 0;
	// Compute f(x0) and A(x0)
	// while( (xi != x*) && (a > MIN_STEP) && (i_step < MAX_ITER) && ( abs(A(xi)-A(xi-1)) > MIN_IMP))
	// {
	//    Compute the gradients gi
	//    Compute search direction si
	//    Compute a
	//    while (STEP_NOT_ACCEPTED && (a > MIN_STEP))
	//    {
	//       Compute xi+1 = a si
	//       Compute f(xi+1) and A(xi+1)
	//       Test for step acceptance
	//       a = a/2
	//    }
	//    i_step++;
	// }

	double real_min_angle_sc = 2.0;
	// calculate real min angle and check with threshold
	for(int i = 0; i < bct; i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*)blocks[i];
		if(tetra->isInverted()){
			point->setCoordinates(dpt_initial);
			return 0.0;
		}
		double a = tetra->getMinDihedralAngleSin();
		if(a < real_min_angle_sc) real_min_angle_sc = a;
	}
	return real_min_angle_sc;
}

/// Calls the optimization-smoothing procedure for all points in the given mesh
bool MeshGenerator3dQuality::smoothenOptMixed(Metric3dContext& mc, MeshContainer3d* mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int pct = mesh->getPointsCount();
	if(pct < 1) return false;

	const double THRESHOLD_1 = 0.2; // ~ sin(10o)
	const double THRESHOLD_2 = 0.1; // ~ sin(5o)

	START_CLOCK("MG3dQ::smoothenOpt");

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if(tag_type != TagExtended::TAG_NONE && !point->checkIntTag(tag_type, tag_value))
			continue;
		if(point->isBorder()) continue;

		MeshGenerator3dQuality::movePointByLaplaceForVariableMetric(mc, point);
		double qt = MeshGenerator3dQuality::movePointOpt(mc, point, THRESHOLD_1);
		if(qt < THRESHOLD_2) 
			MeshGenerator3dQuality::movePointUsingAngleOptimization(point);
	}

	STOP_CLOCK("MG3dQ::smoothenOpt");

	return true;
}

/// Calls the optimization-smoothing procedure for all points in the given mesh
bool MeshGenerator3dQuality::smoothenOpt(Metric3dContext& mc, MeshContainer3d* mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int pct = mesh->getPointsCount();
	if(pct < 1) return false;

	START_CLOCK("MG3dQ::smoothenOpt");

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if(tag_type != TagExtended::TAG_NONE && !point->checkIntTag(tag_type, tag_value))
			continue;
		if(point->isBorder()) continue;
		MeshGenerator3dQuality::movePointOpt(mc, point);
	}

	STOP_CLOCK("MG3dQ::smoothenOpt");

	return true;
}

bool MeshGenerator3dQuality::smoothenLaplace(Metric3dContext& mc, MeshContainer3d* mesh, 
		TagExtended::TagType tag_type, int tag_value, bool with_boundary)
{
//	if(!mesh || mesh->getInnerEdgesCount() > 0) return false;
	if(!mesh) return false;
	int pct = mesh->getPointsCount();
	if(pct < 1) return true;

//	SHOW_STEP_BREAKABLE(2, "* Wyg³adzanie Laplace'a.", 0.0);

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "MG3d::smoothenLaplace");
	START_CLOCK("MG3d::smoothenLaplace");

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if(tag_value != TagExtended::TAG_NONE && !point->checkIntTag(tag_type, tag_value))
			continue;
		if(point->isBorder()){
			if(with_boundary)
				moveBoundaryPointByLaplace(mc, point);
		}else
			movePointByLaplace(mc, point);
	}
	STOP_CLOCK("MG3d::smoothenLaplace");

	return true;
}

bool MeshGenerator3dQuality::smoothenSwap(Metric3dContext& mc, MeshContainer3d* mesh,
			TagExtended::TagType /* tag_type */, int /* tag_value */)
{
//	if(!mesh || mesh->getInnerEdgesCount() > 0) return false;
	if(!mesh) return false;
	int bct = mesh->getBlocksCount();
	if(bct < 1) return true;
	int pct = mesh->getPointsCount();

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "MG3d::smoothenSwap");
	START_CLOCK("MG3d::smoothenSwap");

	// swap23
	int swap23_count = 0;
	for(int i = 0; i < bct; i++){
		MeshTetrahedron* tetrahedron = (MeshTetrahedron*)mesh->getBlockAt(i);
		if(tetrahedron->getType() != BLOCK_TETRA) continue;
		for(int j = 0; j < 4; j++){
			MeshFace* face = tetrahedron->getFace(j);
			if(face->isBorder()) continue;
			MeshTetrahedron* tetrahedron_2 = (MeshTetrahedron*)face->getOtherBlock(tetrahedron);
			if(tetrahedron_2->getType() != BLOCK_TETRA) continue;
			mc.countMetricAtPoint(face->getMiddlePoint());
			// swap if it helps
			if(MeshGenerator3d::swap23(mc, mesh, face, nullptr, param_swap_criterion)){
				++swap23_count;
				// blocks were removed/inserted -> update count, skip this tetrahedron
				bct = mesh->getBlocksCount();
				break;
			}
		}
	}
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
		"Tetrahedral smoothing - " << swap23_count << " 23-swaps");

	// swap32
	int swap32_count = 0;
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(edge->isBorder()) continue;
			if(edge->getPointIndex(point) == 0){
				mc.countMetricAtPoint(edge->getPoint(0.5));
				// swap if it helps
				if(MeshGenerator3d::swap32(mc, mesh, edge, nullptr, param_swap_criterion)) ++swap32_count;
			}
		}
	}
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Tetrahedral smoothing - " << swap32_count << " 32-swaps");

	STOP_CLOCK("MG3d::smoothenSwap");

	return true;
}

bool MeshGenerator3dQuality::smoothenSwapParallel(Metric3dContext& mc, MeshContainer3d* mesh)
{
//	if(!mesh || mesh->getInnerEdgesCount() > 0) return false;
	if(!mesh) return false;
	int pct = mesh->getPointsCount();
	if(pct < 1) return true;

	int swap23_count = 0;
	int swap32_count = 0;

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "MG3d::smoothenSwapParallel");
	START_CLOCK("MG3d::smoothenSwapParallel");

	const int IS_USED = 1;
	const int IS_DONE = 2;
	int pct_left = pct;
	int np = 1;
	int ip = 0;
	int i = 0;
#pragma omp parallel shared(pct, pct_left) private(np)
	{
#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
		np = omp_get_num_threads();
		ip = omp_get_thread_num();
		i = ip*pct/np;
#endif //_OPENMP
		while(pct_left > 0){

			if(pct_left <= np && ip > 0) break; // leave last np points to master...

			// select next free point
			MeshPoint3d* selected_point = nullptr;
			for(int j = 0; !selected_point && j < pct; j++){
				MeshPoint3d* point = mesh->getPointAt((i+j)%pct);

				#pragma omp critical (tag_mesh)
				{
					if(!point->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED | IS_DONE)){
						point->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
						selected_point = point;
	//					LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Thread " << ip << '/' << np << " - checking point " 
	//						<< selected_point->getIndex() << endl;
					}
				}
			}
			if(!selected_point) continue; // apparently some nodes were being used, maybe pct_left become 0 in the meantime...

			int change_count = 0;
			int collision_count = 0;
			int attempt = 0;
			do{
				++attempt;
				collision_count = change_count = 0;
				for(int ei = 0; ei < selected_point->getRank(); ei++){
					MeshEdge3d* edge = selected_point->getEdge(ei);
					MeshPoint3d* other_point = edge->getOtherPoint(selected_point);

					bool collision_leave = false;
					#pragma omp critical (tag_mesh)
					{
						if(other_point->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED)){
							++collision_count; 
							collision_leave = true; 
						}
						else other_point->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
					}

					if(collision_leave) continue;

					if(!edge->isBorder() && 
						!edge->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_DONE) && 
						edge->getFaceCount() == 3)
					{ // try gather points for swap32
						MeshPoint3d* points[3] = {
							edge->getFaceAt(0)->getOtherPoint(selected_point, other_point), 
							edge->getFaceAt(1)->getOtherPoint(selected_point, other_point), 
							edge->getFaceAt(2)->getOtherPoint(selected_point, other_point)};

						#pragma omp critical (tag_mesh)
						{
							if(points[0]->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED) || 
								points[1]->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED) ||
								points[2]->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED))
							{ 
								points[0] = nullptr; 
							}else{
								points[0]->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
								points[1]->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
								points[2]->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
							}
						}
						if(points[0]){ // try swap32
							const DPoint3d middle = DPoint3d::average(
								selected_point->getCoordinates(),
								other_point->getCoordinates(),
								points[0]->getCoordinates(),
								points[1]->getCoordinates(),
								points[2]->getCoordinates());
							mc.countMetricAtPoint(middle);
							MeshTetrahedron* tetrahedra[2];
							if(MeshGenerator3d::swap32(mc, mesh, edge, tetrahedra, param_swap_criterion)){
								++swap32_count;
								++change_count;
								tetrahedra[0]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP_PAR);
								tetrahedra[1]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP_PAR);
								other_point->clearIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
								for(int k = 0; k < 3; k++)
									points[k]->clearIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
								continue;
							}else{
								edge->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_DONE);
								for(int k = 0; k < 3; k++)
									points[k]->clearIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
							}
						}else ++collision_count;
					}

					for(int fi = 0; fi < edge->getFaceCount(); fi++){
						MeshFace* face = edge->getFaceAt(fi);
						if(face->isBorder() || face->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_DONE)) continue;
						MeshTetrahedron* tetrahedra[3] = {(MeshTetrahedron*)face->getBlock(0),
							(MeshTetrahedron*)face->getBlock(1), nullptr};

						MeshPoint3d* points[3] = { face->getOtherPoint(selected_point, other_point),
							tetrahedra[0]->getOppositePoint(face), tetrahedra[1]->getOppositePoint(face) };

						bool collision_leave = false;
						#pragma omp critical (tag_mesh) // gather points
						{
							if(points[0]->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED) || 
								points[1]->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED) || 
								points[2]->hasIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED)) 
							{ 
								++collision_count; collision_leave = true; 
							}else{
								points[0]->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
								points[1]->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
								points[2]->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
							}
						}
						if(collision_leave) continue;
						// try swap23
						const DPoint3d middle = DPoint3d::average(
							selected_point->getCoordinates(),
							other_point->getCoordinates(),
							points[0]->getCoordinates(),
							points[1]->getCoordinates(),
							points[2]->getCoordinates());
						mc.countMetricAtPoint(middle);
						if(MeshGenerator3d::swap23(mc, mesh, face, tetrahedra, param_swap_criterion)){
							++swap23_count;
							++change_count;
							for(int k = 0; k < 3; k++){
								tetrahedra[k]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP_PAR);
								points[k]->clearIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
							}
							other_point->clearIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
							continue;
						}else{
							face->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_DONE);
							for(int k = 0; k < 3; k++)
								points[k]->clearIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
						}
					}
					other_point->clearIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
				}
			}while(change_count > 0 || (collision_count > 0 && attempt <= 10));
			if(attempt <= 10) selected_point->setIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_DONE);
			selected_point->clearIntFlag(TagExtended::TAG_MG3D_SM_SWAP_PAR, IS_USED);
		}
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Tetrahedral smoothing (parallel) - " << swap23_count << " 23-swaps");
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Tetrahedral smoothing (parallel) - " << swap32_count << " 32-swaps");
	STOP_CLOCK("MG3d::smoothenSwapParallel");

	mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP_PAR);


	return true;
}

#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
bool MeshGenerator3dQuality::smoothenSwapComplete(Metric3dContext& mc, MeshContainer3d* mesh)
{
	if(!mesh) return false;
	int pct = mesh->getPointsCount();
	if(pct < 4) return false;

	LOG4CPLUS_INFO(MeshLog::logger_console, "MG3d::smoothenSwapComplete");
	START_CLOCK("MG3d::smoothenSwapComplete");

	int * active_points  = new int[pct];
	int * points_counter = new int[pct];
	const int MAX_COUNTER = 100;

	#pragma omp parallel for
	for(int i = 0; i < pct; i++){
		active_points[i] = i;
		points_counter[i] = 0;
	}

	int s32_checked = 0;
	int s32_swapped = 0;
	int s23_checked = 0;
	int s23_swapped = 0;
	int pts_checked = 0;

	int pct_left = pct;

#pragma omp parallel shared(pct, pct_left, active_points, points_counter, mesh, mc) \
					reduction(+:s32_checked) reduction(+:s32_swapped) \
					reduction(+:s23_checked) reduction(+:s23_swapped) \
					reduction(+:pts_checked)
	{
		int thr_id = omp_get_thread_num();
		int thr_ct = omp_get_num_threads();
		int i = thr_id * pct / thr_ct;
		MeshPoint3d* points[5];
		MeshTetrahedron** tetrahedra = new MeshTetrahedron*[3];
		Metric3dContext mc_local(mesh->getControlSpace());

		while(true){
			int selected = -1;

			#pragma omp critical (swap3d_complete)
			{
				if(pct_left > 0 && (pct_left > thr_ct || thr_id == 0)){
					i = (i+1) % pct_left;
					selected = active_points[i];
					active_points[i] = active_points[--pct_left];
				}
			}

			if(selected < 0) break;

			points[0] = mesh->getPointAt(selected);

			#pragma omp critical (swap3d_complete)
			{
				if(points[0]->isActive()) points[0] = nullptr;
				else points[0]->setActive();
			}

			if(points[0]) points[0]->setIntTag(TagExtended::TAG_MG3D_SM_SWAP);
			else continue;

			++pts_checked;
			for(int j = 0; j < points[0]->getRank(); j++){
				// check edge
				MeshEdge3d* edge = points[0]->getEdge(j);
				if(edge->isBorder()) continue;
				points[1] = edge->getOtherPoint(points[0]);

				#pragma omp critical (swap3d_complete)
				{
					if(points[1]->isActive()) points[1] = nullptr;
					else points[1]->setActive();
				}

				if(!points[1]) continue;

				bool swap_ok = false;
				points[2] = nullptr;
				if(edge->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP) && 
					edge->getFaceCount() == 3)
				{ // if tagged -> already checked
					points[2] = edge->getFaceAt(0)->getOtherPoint(points[0], points[1]); 
					points[3] = edge->getFaceAt(1)->getOtherPoint(points[0], points[1]);
					points[4] = edge->getFaceAt(2)->getOtherPoint(points[0], points[1]);

					#pragma omp critical (swap3d_complete)
					{
						if(points[2]->isActive() || points[3]->isActive() || points[4]->isActive()){
							points[2] = nullptr;
						}else{
							points[2]->setActive();
							points[3]->setActive();
							points[4]->setActive();
						}
					}

					if(!points[2]){
						points[1]->clearActive();
						continue;
					}

					++s32_checked;
					const DPoint3d middle = DPoint3d::average( // metric set for 5 points (don't depend on swap)
						points[0]->getCoordinates(), 
						points[1]->getCoordinates(), 
						points[2]->getCoordinates(),
						points[3]->getCoordinates(), 
						points[4]->getCoordinates());
					mc_local.countMetricAtPoint(middle);
					// try swap
					
					swap_ok = MeshGenerator3d::swap32(mc_local, mesh, edge, tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
				}
				if(swap_ok){
					++s32_swapped;
					j = 0; 
					tetrahedra[0]->clearTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
					tetrahedra[1]->clearTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
					for(int m = 0; m < 5; m++){
						if(points[m]->nonZeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){
							int id = points[m]->getIndex();

							#pragma omp critical (swap3d_complete)
							{
								if(++points_counter[id] < MAX_COUNTER){
									active_points[pct_left++] = id;
									points[m]->clearTag(TagExtended::TAG_MG3D_SM_SWAP);
								}
							}
						}
					}
					points[1]->clearActive();
					points[2]->clearActive();
					points[3]->clearActive();
					points[4]->clearActive();
					continue; 
				}else{
					edge->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
					if(points[2]){
						points[2]->clearActive();
						points[3]->clearActive();
						points[4]->clearActive();
					}
				}

				// check faces
				for(int k = 0; k < edge->getFaceCount(); k++){
					MeshFace* face = edge->getFaceAt(k);
					if(face->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){ // if tagged -> already checked
						tetrahedra[0] = (MeshTetrahedron*)face->getBlock(0);
						tetrahedra[1] = (MeshTetrahedron*)face->getBlock(1);
						points[2] = face->getOtherPoint(points[0], points[1]);
						points[3] = tetrahedra[0]->getOppositePoint(face); 
						points[4] = tetrahedra[1]->getOppositePoint(face);

						#pragma omp critical (swap3d_complete)
						{
							if(points[2]->isActive() || points[3]->isActive() || points[4]->isActive()){
								points[2] = nullptr;
							}else{
								points[2]->setActive();
								points[3]->setActive();
								points[4]->setActive();
							}
						}

						if(!points[2]) continue;

						++s23_checked;
						const DPoint3d middle = DPoint3d::average(
							points[0]->getCoordinates(),
							points[1]->getCoordinates(), 
							points[2]->getCoordinates(),
							points[3]->getCoordinates(),
							points[4]->getCoordinates());
						mc_local.countMetricAtPoint(middle);
						swap_ok = MeshGenerator3d::swap23(mc_local, mesh, face, tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
					}else
						points[2] = nullptr;

					if(swap_ok){ // start again from first edge
						++s23_swapped;
						j = 0; 
						for(int m = 0; m < 3; m++) 
							tetrahedra[m]->clearTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
						for(int m = 0; m < 5; m++)
							if(points[m]->nonZeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){
								int id = points[m]->getIndex();

								#pragma omp critical (swap3d_complete)
								{
									if(++points_counter[id] < MAX_COUNTER){
										active_points[pct_left++] = id;
										points[m]->clearTag(TagExtended::TAG_MG3D_SM_SWAP);
									}
								}
						}
						points[2]->clearActive();
						points[3]->clearActive();
						points[4]->clearActive();
						break; 
					}else{
						face->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
						if(points[2]){
							points[2]->clearActive();
							points[3]->clearActive();
							points[4]->clearActive();
						}
					}
				}

				points[1]->clearActive();
			}
			points[0]->clearActive();
		}
		delete[] tetrahedra;
	}

	STOP_CLOCK("MG3d::smoothenSwapComplete");

	delete[] active_points;
	delete[] points_counter;
	mesh->clearAllTags(TagExtended::TAG_MG3D_SM_SWAP);  // clear all tags for points, edges and faces

	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Tetrahedra smoothing complete, pts: " << pts_checked << '/' << pct 
		<< ", swap23: " << s23_checked << '/' << s23_swapped 
		<< ", swap32: " << s32_checked << '/' << s32_swapped << endl;

	return true;
}
#else
bool MeshGenerator3dQuality::smoothenSwapComplete(Metric3dContext& mc, MeshContainer3d* mesh, 
		TagExtended::TagType tag_type, int tag_value, bool with_boundary)
{
	if(!mesh) return false;
	int pct = mesh->getPointsCount();
	if(pct < 4) return false;

//	LOG4CPLUS_INFO(MeshLog::logger_console, "MG3d::smoothenSwapComplete");
	START_CLOCK("MG3d::smoothenSwapComplete");

	DataVector<int> active_points(pct, 0);
	DataVector<int> points_counter(pct, 0);
	const int MAX_COUNTER = 100;

	for(int i = 0; i < pct ; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if(tag_value == TagExtended::TAG_NONE || point->checkIntTag(tag_type, tag_value))
			active_points.add(i);
	}

	MeshTetrahedron** tetrahedra = new MeshTetrahedron*[3];
	MeshPoint3d* points[5];

	int s32_checked = 0;
	int s32_swapped = 0;
	int s23_checked = 0;
	int s23_swapped = 0;
	int s22_checked = 0;
	int s22_swapped = 0;
	int pts_checked = 0;

	while(active_points.countInt() > 0){
		points[0] = mesh->getPointAt(active_points.removeLast());
		points[0]->setIntTag(TagExtended::TAG_MG3D_SM_SWAP);
		++pts_checked;
		for(int j = 0; j < points[0]->getRank(); j++){
			// check edge
			MeshEdge3d* edge = points[0]->getEdge(j);
			if(edge->isBorder()){ 
				if(with_boundary){
					points[1] = edge->getOtherPoint(points[0]);
					bool swap_ok = false;
					if(edge->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP) && 
						edge->getFaceCount() == 3)
					{ // if tagged -> already checked
						++s22_checked;
						points[2] = edge->getFaceAt(0)->getOtherPoint(points[0], points[1]); 
						points[3] = edge->getFaceAt(1)->getOtherPoint(points[0], points[1]);
						points[4] = edge->getFaceAt(2)->getOtherPoint(points[0], points[1]);
						const DPoint3d middle = DPoint3d::average( // metric set for 5 points (don't depend on swap)
							points[0]->getCoordinates(), 
							points[1]->getCoordinates(), 
							points[2]->getCoordinates(),
							points[3]->getCoordinates(), 
							points[4]->getCoordinates());
						mc.countMetricAtPoint(middle);
						// try swap
						swap_ok = MeshGenerator3d::swap22(mc, mesh, edge, tetrahedra) != nullptr;
					}
					if(swap_ok){
						++s22_swapped;
						j = 0; 
						for(int m = 0; m < 2; m++) tetrahedra[m]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
						for(int m = 0; m < 5; m++) 
							if(points[m]->nonZeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){
								int id = points[m]->getIndex();
								if(++points_counter[id] < MAX_COUNTER){
									if(tag_value == TagExtended::TAG_NONE || points[m]->checkIntTag(tag_type, tag_value))
										active_points.add(id);
									points[m]->setZeroTag(TagExtended::TAG_MG3D_SM_SWAP);
								}
							}
						continue; 
					}else edge->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
				}
				continue;
			}
			points[1] = edge->getOtherPoint(points[0]);
			bool swap_ok = false;
			if(edge->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP) && 
				edge->getFaceCount() == 3)
			{ // if tagged -> already checked
				++s32_checked;
				points[2] = edge->getFaceAt(0)->getOtherPoint(points[0], points[1]); 
				points[3] = edge->getFaceAt(1)->getOtherPoint(points[0], points[1]);
				points[4] = edge->getFaceAt(2)->getOtherPoint(points[0], points[1]);
				const DPoint3d middle = DPoint3d::average( // metric set for 5 points (don't depend on swap)
					points[0]->getCoordinates(), 
					points[1]->getCoordinates(), 
					points[2]->getCoordinates(),
					points[3]->getCoordinates(), 
					points[4]->getCoordinates());
				mc.countMetricAtPoint(middle);
				// try swap
				swap_ok = MeshGenerator3d::swap32(mc, mesh, edge, tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
			}
			if(swap_ok){
				++s32_swapped;
				j = 0; 
				for(int m = 0; m < 2; m++) tetrahedra[m]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
				for(int m = 0; m < 5; m++) 
					if(points[m]->nonZeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){
						int id = points[m]->getIndex();
						if(++points_counter[id] < MAX_COUNTER){
							if(tag_value == TagExtended::TAG_NONE || points[m]->checkIntTag(tag_type, tag_value))
								active_points.add(id);
							points[m]->setZeroTag(TagExtended::TAG_MG3D_SM_SWAP);
						}
					}
				continue; 
			}else edge->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
			// check faces
			for(int k = 0; k < edge->getFaceCount(); k++){
				MeshFace* face = edge->getFaceAt(k);
				if(face->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){ // if tagged -> already checked
					++s23_checked;
					tetrahedra[0] = (MeshTetrahedron*)face->getBlock(0);
					tetrahedra[1] = (MeshTetrahedron*)face->getBlock(1);
					points[2] = face->getOtherPoint(points[0], points[1]);
					points[3] = tetrahedra[0]->getOppositePoint(face); 
					points[4] = tetrahedra[1]->getOppositePoint(face);
					const DPoint3d middle = DPoint3d::average(
						points[0]->getCoordinates(),
						points[1]->getCoordinates(), 
						points[2]->getCoordinates(),
						points[3]->getCoordinates(),
						points[4]->getCoordinates());
					mc.countMetricAtPoint(middle);
					swap_ok = MeshGenerator3d::swap23(mc, mesh, face, tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
				}
				if(swap_ok){ // start again from first edge
					++s23_swapped;
					j = 0; 
					for(int m = 0; m < 3; m++) tetrahedra[m]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
					for(int m = 0; m < 5; m++)
						if(points[m]->nonZeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){
							int id = points[m]->getIndex();
							if(++points_counter[id] < MAX_COUNTER){
								if(tag_value == TagExtended::TAG_NONE || points[m]->checkIntTag(tag_type, tag_value))
									active_points.add(id);
								points[m]->setZeroTag(TagExtended::TAG_MG3D_SM_SWAP);
							}
					}
					break; 
				}else face->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
			}
		}
	}

	delete[] tetrahedra;

	STOP_CLOCK("MG3d::smoothenSwapComplete");

	mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);  // clear all tags for points, edges and faces

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
		"Tetrahedra smoothing complete, pts: " << pts_checked << '/' << pct 
		<< ", swap23: " << s23_checked << '/' << s23_swapped 
		<< ", swap32: " << s32_checked << '/' << s32_swapped
		<< ", swap22: " << s22_checked << '/' << s22_swapped);

	return true;
}
#endif //_OPENMP

bool MeshGenerator3dQuality::smoothenLaplaceMixed(Metric3dContext& mc, MeshContainer3d *mesh, 
		TagExtended::TagType tag_type, int tag_value, bool with_boundary)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return true;

//	LOG4CPLUS_INFO(MeshLog::logger_console, "MG3d::smoothenLaplaceMixed");
	START_CLOCK("MG3d::smoothenLaplaceMixed");

#pragma omp parallel for
	for(int i = 0; i < count; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if(tag_value != TagExtended::TAG_NONE && !point->checkIntTag(tag_type, tag_value))
			continue;

		double mg = mc.getMetricGradation(point->getCoordinates());
		if(mg > 1.5){
			if(point->isBorder()){
				if(with_boundary) moveBoundaryPointByLaplaceForVariableMetric(mc, point);
			}else
				movePointByLaplaceForVariableMetric(mc, point);
		}else{ 
			if(point->isBorder()){
				if(with_boundary) moveBoundaryPointByLaplace(mc, point);
			}else
				movePointByLaplace(mc, point);
		}
	}

	STOP_CLOCK("MG3d::smoothenLaplaceMixed");
	return true;
}

//bool MeshGenerator3dQuality::movePointByLaplaceForVariableMetric(
//	Metric3dContext& mc, MeshPoint3d *point)
//{
//	if(point->isBorder()) return false;
//#ifdef _OPENMP
//	point->setActive();
//#endif
//
//	DMVector3d total_mdv(0.0, 0.0, 0.0);
//	int rank = point->getRank();
//	double total_weight = 0.0;
//	for(int j = 0; j < rank; j++){
//		const MeshEdge3d* edge = point->getEdge(j);
//		MeshPoint3d* other_point = edge->getOtherPoint(point);
//		#ifdef _OPENMP
//		if(other_point->isActive()){
//			point->clearActive();
//			return false;
//		}
//		#endif
//		mc.countMetricAtPoint(edge->getPoint(0.5));
//		const DMVector3d mdv = other_point->getMetricCoordinates(mc) - point->getMetricCoordinates(mc);
//		total_mdv += mdv;
//		total_weight += 1.0;
//	}
//	total_mdv /= total_weight;
//	mc.countMetricAtPoint(point->getCoordinates());
//	bool res = tryMovingPoint(point, point->getCoordinates() + mc.transformMStoRS(total_mdv));
//#ifdef _OPENMP
//	point->clearActive();
//#endif
//	return res;
//}

// -> this version: better shape quality, smoothen edges, worse metric-match
bool MeshGenerator3dQuality::movePointByLaplaceForVariableMetric(
		Metric3dContext& mc, MeshPoint3d *point)
{
	if(point->isBorder()) return false;
#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
	point->setActive();
#endif

	const DPoint3d& dpt = point->getCoordinates();
	int rank = point->getRank();
	DataVector< DVector3d > dvs(rank);
	DataVector< double > lens(rank);
	double ave_len = 0.0;
	for(int j = 0; j < rank; j++){
		MeshEdge3d* edge = point->getEdge(j);
		MeshPoint3d* other_point = edge->getOtherPoint(point);
		#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
		if(other_point->isActive()){
			point->clearActive();
			return false;
		}
		#endif
		DPoint3d other_dpt = other_point->getCoordinates();
		mc.countMetricAtPoints(point, other_point);
		
		double len = mc.transformRStoMS( other_dpt - dpt).length();
		lens.add( len );
		dvs.add( other_dpt - dpt );
		ave_len += len;
	}
	ave_len *= (0.9 / rank); // calculate average length of edges and make it a little shorter

	DVector3d total_dv;
	double total_w = 0.0;
	for(int j = 0; j < rank; j++){
		double len_inv = 1.0 / lens[j];
		double ratio = (lens[j] - ave_len) * len_inv;
		DVector3d dv = dvs[j] * ratio;;
		total_dv += dv * len_inv;
		total_w += len_inv;
	}
	total_dv /= total_w;
	mc.countMetricAtPoint(point->getCoordinates());

	bool res =  tryMovingPoint(mc, point, dpt + total_dv);
#if defined(USE_OPENMP_HERE) && defined(_OPENMP)
	point->clearActive();
#endif
	return res;
}

bool MeshGenerator3dQuality::moveBoundaryPointByLaplaceForVariableMetric(
		Metric3dContext& mc, MeshPoint3d *point)
{
	if(!point->isBorder()) return false;
	if(point->isBorder(TagBorder::CORNER | TagBorder::FIXED)) return false; // corner and/or fixed points are immovable

	bool is_ridge = point->isBorder(TagBorder::RIDGE); // ridge-point moved only with respect to other ridge-points

	DPoint3d dpt = point->getCoordinates();
	mc.countMetricAtPoint( dpt );
	int rank = point->getRank();
	DataVector< DVector3d > dvs(rank);
	DataVector< double > lens(rank);
	double ave_len = 0.0;
	for(int j = 0; j < rank; j++){
		MeshEdge3d* edge = point->getEdge(j);
		if(!edge->isBorder()) continue; // select only boundary edges
		if(is_ridge && !edge->isBorder(TagBorder::RIDGE)) continue;
		MeshPoint3d* other_point = edge->getOtherPoint(point);
		DPoint3d other_sdpt = other_point->getCoordinates();
		mc.countMetricAtPoints(point, other_point);
		
		double len = mc.transformRStoMS( other_sdpt - dpt ).length();
		lens.add( len );
		dvs.add( other_sdpt - dpt );
		ave_len += len;
	}

	rank = lens.countInt();
	ave_len *= (0.9 / rank); // calculate average length of edges and make it a little shorter

	DVector3d total_dv;
	double total_w = 0.0;
	for(int j = 0; j < rank; j++){
		double len_inv = 1.0 / lens[j];
		double ratio = (lens[j] - ave_len) * len_inv;
		DVector3d dv = dvs[j] * ratio;;
		total_dv += dv * len_inv;
		total_w += len_inv;
	}
	total_dv /= total_w;
	mc.countMetricAtPoint(dpt);

	assert(false);
	// need to be rewritten, as in Generator3dSurface ...
	//return tryMovingPoint(mc, point, point->getCoordinatesWithLocalShape(mc, dpt + total_dv), false);
	return false;
}

bool MeshGenerator3dQuality::smoothen(Metric3dContext& mc, MeshContainer3d* mesh, 
	int steps, TagExtended::TagType tag_type, int tag_value, bool with_boundary, 
	int method)
{
	DataVector<MeshContainer3d::AdditionalBoundaryNode>& inserted_boundary_points = 
		mesh->getAdditionalBoundaryNodes();

//	SHOW_MESH("smoothen-start", mesh->getViewSet());

	for(int i = 0; i < steps; i++){
		if((method & MeshData::SM3_SWAP) != 0){
			MeshGenerator3dQuality::smoothenSwap(mc, mesh, tag_type, tag_value); 
//			SHOW_MESH("after smoothenSwap", mesh->getViewSet());
		}
		if((method & MeshData::SM3_SWAP_COMPLETE) != 0){
			MeshGenerator3dQuality::smoothenSwapComplete(mc, mesh, tag_type, tag_value, with_boundary);
//			SHOW_MESH("after smoothenSwapComplete", mesh->getViewSet());
		}
		if((method & MeshData::SM3_LAPLACE) != 0){
			MeshGenerator3dQuality::smoothenLaplace(mc, mesh, tag_type, tag_value, with_boundary); 
//			SHOW_MESH("after smoothenLaplace", mesh->getViewSet());
		}
		if((method & MeshData::SM3_LAPLACE_MIXED) != 0){
			MeshGenerator3dQuality::smoothenLaplaceMixed(mc, mesh, tag_type, tag_value, with_boundary); 
//			SHOW_MESH("after smoothenLaplaceMixed", mesh->getViewSet());
		}
		if((method & MeshData::SM3_OPT_MIXED) != 0){
			MeshGenerator3dQuality::smoothenOptMixed(mc, mesh, tag_type, tag_value);
//			SHOW_MESH("after smoothenOptMixed", mesh->getViewSet());
		}
		if((method & MeshData::SM3_OPT_SIMPLE) != 0){
			MeshGenerator3dQuality::smoothenOpt(mc, mesh, tag_type, tag_value);
//			SHOW_MESH("after smoothenOpt", mesh->getViewSet());
		}
		if(inserted_boundary_points.countInt() > 0){
			MeshGenerator3dDelaunayBoundary::fixAdditionalBoundaryNodes(mc, mesh, inserted_boundary_points);
			LOG4CPLUS_DEBUG(MeshLog::logger_console, 
				"Additional boundary nodes (final-smoothing) count: " << inserted_boundary_points.countInt());
//			SHOW_MESH("after fixAdditionalBoundaryNodes", mesh->getViewSet());
		}
		if((method & MeshData::SM3_BORDER_PRUNE) != 0){
			MeshGenerator3dQuality::removeBoundarySlivers(mc, mesh, 0.2, tag_type, tag_value);
			SHOW_MESH("after removeBoundarySlivers", mesh->getViewSet());
		}
	}

	return true;
}

/// Attempts to improve the worst tetrahedra by cavity retriangulation
bool MeshGenerator3dQuality::improveTetrahedron(Metric3dContext& mc, MeshContainer3d* mesh, MeshTetrahedron* tetra)
{
	if(tetra->getType() != BLOCK_TETRA) return false; // tetrahedra only

	int inner_point_count = 0;
	for(int j = 0; j < tetra->getPointCount(); j++){
		if(!tetra->getPoint(j)->isBorder())
			inner_point_count++;
	}
	if(inner_point_count > 0) return false;

	DataVector<MeshBlock*> blocks;
	blocks.add(tetra);
	MeshEdge3d* edge = nullptr;
	for(int j = 0; j < tetra->getEdgeCount(); j++){
		edge = tetra->getEdge(j);
		if(edge->isBorder()) continue;
		edge->adjacentBlocks(blocks);
		break;
	}
	if(!edge) return false;

	// insert node at the middle of the edge
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Insert new node for tetrahedron improvement");

	MeshPoint3d* mid_point = new MeshPoint3d( edge->getPoint(0.5) );
	mesh->addMeshPoint(mid_point);
	MeshPoint3d * ep0 = edge->getMeshPoint(0);
	MeshPoint3d * ep1 = edge->getMeshPoint(1);

	for(int j = 0; j < blocks.countInt(); j++){
		MeshTetrahedron * t = (MeshTetrahedron*)blocks[j];
		mesh->removeMeshBlock(t);

		MeshTetrahedron * tn = new MeshTetrahedron(
			(t->getPoint(0) == ep0) ? mid_point : t->getPoint(0),
			(t->getPoint(1) == ep0) ? mid_point : t->getPoint(1),
			(t->getPoint(2) == ep0) ? mid_point : t->getPoint(2),
			(t->getPoint(3) == ep0) ? mid_point : t->getPoint(3));

		mesh->addMeshBlock(tn);
		tn->setAreaID( t->getAreaID() );

		tn = new MeshTetrahedron(
			(t->getPoint(0) == ep1) ? mid_point : t->getPoint(0),
			(t->getPoint(1) == ep1) ? mid_point : t->getPoint(1),
			(t->getPoint(2) == ep1) ? mid_point : t->getPoint(2),
			(t->getPoint(3) == ep1) ? mid_point : t->getPoint(3));

		mesh->addMeshBlock(tn);
		tn->setAreaID( t->getAreaID() );
	}

	movePointUsingSimpleOptimization(mc, mid_point);

	return true;
}

/// Attempts to improve tetrahedra with point relocation using simple optimization
bool MeshGenerator3dQuality::optimizeNearBoundarySimple(Metric3dContext& mc, MeshContainer3d* mesh)
{
	DataStatistics stat;

	int pct = mesh->getPointsCount();
	for(int pi = 0; pi < pct; pi++){
		MeshPoint3d* point = mesh->getPointAt(pi);
		if(point->isBorder()) continue;
		int border_count = 0;
		for(int j = 0; j < point->getRank(); j++)
			if(point->getEdge(j)->getOtherPoint(point)->isBorder())
				border_count++;

		if(border_count == 0) continue;
		// gather adjacent blocks
		DataVector<MeshBlock*> blocks;
		if(!point->adjacentBlocks(blocks)) continue;

		DBox ubox;
		int bct = blocks.countInt();
		for(int j = 0; j < bct; j++){
			int bpct = blocks[j]->getPointCount();
			for(int k = 0; k < bpct; k++)
				ubox.addPoint(blocks[j]->getPoint(k)->getCoordinates());
		}

		static const int MAX_RES = 15;
		static const int MAX_SIZE = 4 * MAX_RES * MAX_RES * MAX_RES;
		static const int REFILL_SIZE = 100;

		double dx = ubox.getDX() / MAX_RES;
		double dy = ubox.getDY() / MAX_RES;
		double dz = ubox.getDZ() / MAX_RES;

		DataVector<PointNode> cells(MAX_SIZE);
		PointNode cell(DPoint3d(0.0, 0.0, ubox.z0 + 0.5*dz));
		for(int i = 0; i < MAX_RES; i++){
			cell.coord.y = ubox.y0 + 0.5*dy;
			for(int j = 0; j < MAX_RES; j++){
				cell.coord.x = ubox.x0 + 0.5*dx;
				for(int k = 0; k < MAX_RES; k++){
					cells.add(cell);
					cell.coord.x += dx;
				}
				cell.coord.y += dy;
			}
			cell.coord.z += dz;
		}
		cells.add(PointNode(point->getCoordinates()));
		
		dx *= 0.5;
		dy *= 0.5;
		dz *= 0.5;
		double dr = std::max(std::max(dx, dy), dz);

		double best_min_quality = 1.0;
		for(int j = 0; j < bct; j++){
			double q = blocks[j]->getMeanRatio(mc);
			if(q < best_min_quality) best_min_quality = q;
		}

		ostringstream log_line;
		log_line << fixed << setprecision(3) << best_min_quality << " ";
		double last_q = best_min_quality;

		int max_steps = 3;
		while(max_steps--){
			DPoint3d dpt = point->getCoordinates();
			// calculate quality
			for(int i = 0; i < cells.countInt(); i++){
				PointNode & pn = cells[i];
				point->setCoordinates(pn.coord);
				double min_q = 1.0;
				for(int j = 0; j < bct; j++){
					double q = blocks[j]->getMeanRatio(mc);
					if(q < min_q) min_q = q;
					if(min_q <= 0.0) break;
				}
				if(min_q > 0) pn.quality = min_q;
				else cells.removeAt(i--);
			}

			// - if no hope...
			if(cells.empty()){
				point->setCoordinates(dpt);
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "---");
				break;
			}

			if(max_steps){
				// - replace cells with lower resolution
				// -- clear low-quality cells if necessary
				if(14*cells.countInt() > REFILL_SIZE){
					std::sort( cells.begin(), cells.end() );
					while(14*cells.countInt() > REFILL_SIZE){
						// - remove last (which is lowest, after sort)
						cells.removeLast();
					}
				}

				log_line << fixed << setprecision(3) << setw(6) << (cells[0].quality - last_q) << " ";
				last_q = cells[0].quality;

				// -- split (it's easier to do from the end
				static const double DRATIO = 0.67;
				dr *= DRATIO;
				dx *= DRATIO;
				dy *= DRATIO;
				dz *= DRATIO;
				for(int i = cells.countInt()-1; i >=0; i--){
					DPoint3d pt = cells[i].coord;
					cells.removeAt(i);
					cells.add(PointNode(DPoint3d(pt.x + dx, pt.y + dy, pt.z + dz)));
					cells.add(PointNode(DPoint3d(pt.x - dx, pt.y + dy, pt.z + dz)));
					cells.add(PointNode(DPoint3d(pt.x + dx, pt.y - dy, pt.z + dz)));
					cells.add(PointNode(DPoint3d(pt.x - dx, pt.y - dy, pt.z + dz)));
					cells.add(PointNode(DPoint3d(pt.x + dx, pt.y + dy, pt.z - dz)));
					cells.add(PointNode(DPoint3d(pt.x - dx, pt.y + dy, pt.z - dz)));
					cells.add(PointNode(DPoint3d(pt.x + dx, pt.y - dy, pt.z - dz)));
					cells.add(PointNode(DPoint3d(pt.x - dx, pt.y - dy, pt.z - dz)));
					cells.add(PointNode(DPoint3d(pt.x + dx, pt.y, pt.z)));
					cells.add(PointNode(DPoint3d(pt.x - dx, pt.y, pt.z)));
					cells.add(PointNode(DPoint3d(pt.x, pt.y + dy, pt.z)));
					cells.add(PointNode(DPoint3d(pt.x, pt.y - dy, pt.z)));
					cells.add(PointNode(DPoint3d(pt.x, pt.y, pt.z + dz)));
					cells.add(PointNode(DPoint3d(pt.x, pt.y, pt.z - dz)));
				}
			}
		}
		
		int best_i = 0;
		for(int j = 1; j < cells.countInt(); j++){
			if(cells[j].quality > cells[best_i].quality) best_i = j;
		}

		point->setCoordinates(cells[best_i].coord);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
			log_line.str() << fixed << "* " << setprecision(3) << cells[best_i].quality);

		stat.add( cells[best_i].quality - best_min_quality);

	}

	if(stat.calculate()){
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Average increase: " << stat.average());
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Maximum increase: " << stat.maximum());
	}

	return false;
}

/// Attempts to improve tetrahedra with point relocation using simple optimization
bool MeshGenerator3dQuality::movePointUsingSimpleOptimization(Metric3dContext& mc, MeshPoint3d* point)
{
	if(point->isBorder()) return false;

	// gather adjacent blocks
	DataVector<MeshBlock*> blocks;
	if(!point->adjacentBlocks(blocks)) return false;

	DBox ubox;
	int bct = blocks.countInt();
	for(int j = 0; j < bct; j++){
		int bpct = blocks[j]->getPointCount();
		for(int k = 0; k < bpct; k++)
			ubox.addPoint(blocks[j]->getPoint(k)->getCoordinates());
	}

	static const int MAX_RES = 15;
	static const int MAX_SIZE = 4 * MAX_RES * MAX_RES * MAX_RES;

	double dx = ubox.getDX() / MAX_RES;
	double dy = ubox.getDY() / MAX_RES;
	double dz = ubox.getDZ() / MAX_RES;

	DataVector<PointNode> cells(MAX_SIZE);
	PointNode cell(DPoint3d(0.0, 0.0, ubox.z0 + 0.5*dz));
	for(int i = 0; i < MAX_RES; i++){
		cell.coord.y = ubox.y0 + 0.5*dy;
		for(int j = 0; j < MAX_RES; j++){
			cell.coord.x = ubox.x0 + 0.5*dx;
			for(int k = 0; k < MAX_RES; k++){
				cells.add(cell);
				cell.coord.x += dx;
			}
			cell.coord.y += dy;
		}
		cell.coord.z += dz;
	}
	cells.add(PointNode(point->getCoordinates()));
	
	dx *= 0.5;
	dy *= 0.5;
	dz *= 0.5;
	double dr = std::max(std::max(dx, dy), dz);

	double best_min_quality = 1.0;
	for(int j = 0; j < bct; j++){
		double q = blocks[j]->getMeanRatio(mc);
		if(q < best_min_quality) best_min_quality = q;
	}

	ostringstream log_line;
	log_line << fixed << setprecision(3) << best_min_quality << " ";
	double last_q = best_min_quality;

	int max_steps = 3;
	while(max_steps--){
		DPoint3d dpt = point->getCoordinates();
		// calculate quality
		for(int i = 0; i < cells.countInt(); i++){
			PointNode & pn = cells[i];
			point->setCoordinates(pn.coord);
			double min_q = 1.0;
			for(int j = 0; j < bct; j++){
				double q = blocks[j]->getMeanRatio(mc);
				if(q < min_q) min_q = q;
				if(min_q <= 0.0) break;
			}
			if(min_q > 0) pn.quality = min_q;
			else cells.removeAt(i--);
		}

		// - if no hope...
		if(cells.empty()){
			point->setCoordinates(dpt);
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, log_line.str() << "---");
			break;
		}

		if(max_steps){
			// - replace cells with lower resolution
			// -- clear low-quality cells if necessary
			if(14*cells.countInt() > MAX_SIZE){
				std::sort( cells.begin(), cells.end() );
				while(14*cells.countInt() > MAX_SIZE){
					// - remove last (which is lowest, after sort)
					cells.removeLast();
				}
			}

			log_line << fixed << setprecision(3) << setw(6) << (cells[0].quality - last_q) << " ";
			last_q = cells[0].quality;

			// -- split (it's easier to do from the end
			static const double DRATIO = 0.67;
			dr *= DRATIO;
			dx *= DRATIO;
			dy *= DRATIO;
			dz *= DRATIO;
			for(int i = cells.countInt()-1; i >=0; i--){
				DPoint3d pt = cells[i].coord;
				cells.removeAt(i);
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y + dy, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y + dy, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y - dy, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y - dy, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y + dy, pt.z - dz)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y + dy, pt.z - dz)));
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y - dy, pt.z - dz)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y - dy, pt.z - dz)));
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y, pt.z)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y, pt.z)));
				cells.add(PointNode(DPoint3d(pt.x, pt.y + dy, pt.z)));
				cells.add(PointNode(DPoint3d(pt.x, pt.y - dy, pt.z)));
				cells.add(PointNode(DPoint3d(pt.x, pt.y, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x, pt.y, pt.z - dz)));
			}
		}
	}
	
	int best_i = 0;
	for(int j = 1; j < cells.countInt(); j++){
		if(cells[j].quality > cells[best_i].quality) best_i = j;
	}

	point->setCoordinates(cells[best_i].coord);
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
		log_line.str() << fixed << "* " << setprecision(3) << cells[best_i].quality);

	return true;
}

/// Attempts to improve tetrahedra with point relocation using simple optimization
bool MeshGenerator3dQuality::movePointUsingAngleOptimization(MeshPoint3d* point)
{
	if(point->isBorder()) return false;

	// gather adjacent blocks
	DataVector<MeshBlock*> blocks;
	if(!point->adjacentBlocks(blocks)) return false;

	DBox ubox;
	int bct = blocks.countInt();
	for(int j = 0; j < bct; j++){
		int bpct = blocks[j]->getPointCount();
		for(int k = 0; k < bpct; k++)
			ubox.addPoint(blocks[j]->getPoint(k)->getCoordinates());
	}

	static const int MAX_RES = 10;
	static const int MAX_SIZE = 4 * MAX_RES * MAX_RES * MAX_RES;

	double dx = ubox.getDX() / MAX_RES;
	double dy = ubox.getDY() / MAX_RES;
	double dz = ubox.getDZ() / MAX_RES;

	DataVector<PointNode> cells(MAX_SIZE);
	PointNode cell(DPoint3d(0.0, 0.0, ubox.z0 + 0.5*dz));
	for(int i = 0; i < MAX_RES; i++){
		cell.coord.y = ubox.y0 + 0.5*dy;
		for(int j = 0; j < MAX_RES; j++){
			cell.coord.x = ubox.x0 + 0.5*dx;
			for(int k = 0; k < MAX_RES; k++){
				cells.add(cell);
				cell.coord.x += dx;
			}
			cell.coord.y += dy;
		}
		cell.coord.z += dz;
	}
	cells.add(PointNode(point->getCoordinates()));
	
	dx *= 0.5;
	dy *= 0.5;
	dz *= 0.5;
	double dr = std::max(std::max(dx, dy), dz);

	double best_min_quality = PI;
	for(int j = 0; j < bct; j++){
		double q = ((MeshTetrahedron*)blocks[j])->getMinDihedralAngleSin();
		if(q < best_min_quality) best_min_quality = q;
	}

//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, fixed << setprecision(3) << best_min_quality << " ";
	double last_q = best_min_quality;

	int max_steps = 3;
	while(max_steps--){
		DPoint3d dpt = point->getCoordinates();
		// calculate quality
		for(int i = 0; i < cells.countInt(); i++){
			PointNode & pn = cells[i];
			point->setCoordinates(pn.coord);
			double min_q = 1.0;
			for(int j = 0; j < bct; j++){
				double q = ((MeshTetrahedron*)blocks[j])->getMinDihedralAngleSin();
				if(q < min_q) min_q = q;
				if(min_q <= 0.0) break;
			}
			if(min_q > 0) pn.quality = min_q;
			else cells.removeAt(i--);
		}

		// - if no hope...
		if(cells.empty()){
			point->setCoordinates(dpt);
//			LOG4CPLUS_INFO(MeshLog::logger_mesh, "---");
			break;
		}

		if(max_steps){
			// - replace cells with lower resolution
			// -- clear low-quality cells if necessary
			if(14*cells.countInt() > MAX_SIZE){
				std::sort( cells.begin(), cells.end() );
//				while(14*cells.countInt() > MAX_SIZE){
				while(cells.countInt() > 10){
					// - remove last (which is lowest, after sort)
					cells.removeLast();
				}
			}

//			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, fixed << setprecision(3) << setw(6) << (cells[0].quality - last_q) << " ";
			last_q = cells[0].quality;

			// -- split (it's easier to do from the end
			static const double DRATIO = 0.67;
			dr *= DRATIO;
			dx *= DRATIO;
			dy *= DRATIO;
			dz *= DRATIO;
			for(int i = cells.countInt()-1; i >=0; i--){
				DPoint3d pt = cells[i].coord;
				cells.removeAt(i);
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y + dy, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y + dy, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y - dy, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y - dy, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y + dy, pt.z - dz)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y + dy, pt.z - dz)));
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y - dy, pt.z - dz)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y - dy, pt.z - dz)));
				cells.add(PointNode(DPoint3d(pt.x + dx, pt.y, pt.z)));
				cells.add(PointNode(DPoint3d(pt.x - dx, pt.y, pt.z)));
				cells.add(PointNode(DPoint3d(pt.x, pt.y + dy, pt.z)));
				cells.add(PointNode(DPoint3d(pt.x, pt.y - dy, pt.z)));
				cells.add(PointNode(DPoint3d(pt.x, pt.y, pt.z + dz)));
				cells.add(PointNode(DPoint3d(pt.x, pt.y, pt.z - dz)));
			}
		}
	}
	
	int best_i = 0;
	for(int j = 1; j < cells.countInt(); j++){
		if(cells[j].quality > cells[best_i].quality) best_i = j;
	}

	point->setCoordinates(cells[best_i].coord);
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, fixed << "* " << setprecision(3) << cells[best_i].quality);

	return true;
}

/// Atempts to improve tetrahedra with quality below the given threshold, using simple optimization
int MeshGenerator3dQuality::optimizeBadTetrahedra(Metric3dContext& mc, MeshContainer3d* mesh, double threshold)
{
	if(!mesh) return 0;

	START_CLOCK("MeshGenerator3dQuality::optimizeBadTetrahedra");

	int counter = 0;
	for(int i = 0; i < mesh->getBlocksCount(); i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*)mesh->getBlockAt(i);
		if(tetra->getType() != BLOCK_TETRA) continue; // tetrahedra only

		//  - 3D shape quality in metric
		double qt = tetra->getMeanRatio(mc);
		if(qt >= threshold) continue;

		//LOG4CPLUS_INFO(MeshLog::logger_console, "Mean ratio, init", qt);
		//if(true){
		//	MeshViewSet* set = new MeshViewSet;
		//	set->addBlockWithEdges(tetra);
		//	SHOW_MESH("Bad tetrahedron", set);
		//}

		for(int j = 0; j < tetra->getPointCount(); j++){
			MeshPoint3d* point = tetra->getPoint(j);
			if(point->isBorder()) continue;
			movePointUsingSimpleOptimization(mc, point);
		}
		qt = tetra->getMeanRatio(mc);

		//LOG4CPLUS_INFO(MeshLog::logger_console, "Mean ratio, opt vertices", qt);
		//if(true){
		//	MeshViewSet* set = new MeshViewSet;
		//	set->addBlockWithEdges(tetra);
		//	SHOW_MESH_NORESET("Bad tetrahedron, opt vertices", set);
		//}
/*
		if(qt >= threshold) continue;

		// TODO -> try to collapse (from shortest non-border edge)
		//		-> if collapse + smoothing successfull, (then maybe additiona flipPIng) and OK
		//		-> else split longest non-border edge
		//		->	then smoothen + flipPIng

		// for now, split first non-border edge...
		DataVector<MeshBlock*> blocks;
		blocks.add(tetra);
		MeshEdge3d* edge = nullptr;
		for(int j = 0; j < tetra->getEdgeCount(); j++){
			edge = tetra->getEdge(j);
			if(edge->isBorder()) continue;
			edge->adjacentBlocks(blocks);
			break;
		}
		if(!edge) continue;

		double min_qt = qt;
		for(int j = 0; j < blocks.countInt(); j++){
			if(blocks[j] == tetra) continue;
			double local_qt = blocks[j]->getMeanRatio(mc);
			if(local_qt < min_qt) min_qt = local_qt;
		}
		LOG4CPLUS_INFO(MeshLog::logger_console, "Min mean ratio before", min_qt);

		if(true){
			MeshViewSet* set = new MeshViewSet;
			for(int j = 0; j < blocks.countInt(); j++){
				set->addBlockWithEdges(blocks[j]);
			}
			SHOW_MESH_NORESET("Bad tetrahedron + neighbours, for split", set);
		}
	
		// insert node at the middle of the edge
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Insert new node for tetrahedron improvement");

		MeshPoint3d* mid_point = new MeshPoint3d( edge->getPoint(0.5) );
		mesh->addMeshPoint(mid_point);
		MeshPoint3d * ep0 = edge->getMeshPoint(0);
		MeshPoint3d * ep1 = edge->getMeshPoint(1);

		DataVector<MeshBlock*> new_blocks(2*blocks.countInt());

		for(int j = 0; j < blocks.countInt(); j++){
			MeshTetrahedron * t = (MeshTetrahedron*)blocks[j];
			MeshPoint3d* tpoints[4] = {
				t->getPoint(0), t->getPoint(1), t->getPoint(2), t->getPoint(3) };
			int area_id = t->getAreaID();
			delete mesh->removeMeshBlock(t);

			MeshTetrahedron * tn = new MeshTetrahedron(
				(tpoints[0] == ep0) ? mid_point : tpoints[0],
				(tpoints[1] == ep0) ? mid_point : tpoints[1],
				(tpoints[2] == ep0) ? mid_point : tpoints[2],
				(tpoints[3] == ep0) ? mid_point : tpoints[3]);

			new_blocks.add(tn);
			mesh->addMeshBlock(tn);
			tn->setAreaID( area_id );

			tn = new MeshTetrahedron(
				(tpoints[0] == ep1) ? mid_point : tpoints[0],
				(tpoints[1] == ep1) ? mid_point : tpoints[1],
				(tpoints[2] == ep1) ? mid_point : tpoints[2],
				(tpoints[3] == ep1) ? mid_point : tpoints[3]);

			new_blocks.add(tn);
			mesh->addMeshBlock(tn);
			tn->setAreaID( area_id );
		}

		movePointUsingSimpleOptimization(mc, mid_point);

		min_qt = 1.0;
		for(int j = 0; j < new_blocks.countInt(); j++){
			double local_qt = new_blocks[j]->getMeanRatio(mc);
			if(local_qt < min_qt) min_qt = local_qt;
		}
		LOG4CPLUS_INFO(MeshLog::logger_console, "Min mean ratio after", min_qt);

		if(true){
			MeshViewSet* set = new MeshViewSet;
			for(int j = 0; j < new_blocks.countInt(); j++){
				set->addBlockWithEdges(new_blocks[j]);
			}
			set->addPoint(mid_point);
			SHOW_MESH_NORESET("Bad tetrahedron + neighbours, after split", set);
		}

		// -> plus face flips ?
*/
		counter++;
	}
	STOP_CLOCK("MeshGenerator3dQuality::optimizeBadTetrahedra");
	return counter;
}

/// Insert new point (using Delaunay) for refinement of the given tetrahedra
bool MeshGenerator3dQuality::insertPointForTetrahedronRefinement(Metric3dContext& mc, 
		MeshContainer3d* mesh, MeshTetrahedron* tetrahedron)
{
	mc.countMetricAtPoint(tetrahedron->getMiddlePoint());

	MeshPoint3d* points[] = {
		tetrahedron->getPoint(0), 
		tetrahedron->getPoint(1),
		tetrahedron->getPoint(2), 
		tetrahedron->getPoint(3)};

	MeshTetrahedron* main_tetrahedron = tetrahedron;
	const DMPoint3d dnew_metric = tetrahedron->getOuterSphereCenter(mc, false);
	DPoint3d dnew = mc.transformMStoRS(dnew_metric);

	if(!tetrahedron->isPointInside(dnew)){
		tetrahedron = tetrahedron->findTetrahedronByNeighbours(dnew, true);
		if(!tetrahedron) return false;
	}
	for(int i = 0; i < 4; i++){
		if(dnew_metric.distance2(tetrahedron->getPoint(i)->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER)
			return false;
	}

	mc.countMetricAtPoint(dnew);
	if(main_tetrahedron != tetrahedron){
		const DMPoint3d dnew_second_metric = mc.transformRStoMS(dnew);
		assert(tetrahedron->isPointInOuterSphere(mc, dnew_second_metric, false)); // since dnew is inside the tetrahedron !!!
		if(!main_tetrahedron->isPointInOuterSphere(mc, dnew_second_metric, false))
			return false;

		for(int i = 0; i < 4; i++){
			if(dnew_second_metric.distance2(tetrahedron->getPoint(i)->getMetricCoordinates(mc)) < 0.25)
				return false;
		}
	}

	MeshPoint3d* p3 = new MeshPoint3d(dnew);
	if(!MeshGenerator3d::addPointToTriangulation(mc, mesh, p3, tetrahedron, main_tetrahedron, true)){
		delete p3;
		return false;
	}

	return true;
}

/// Attempts to improve tetrahedra with angles below the given threshold, using local modifications
int MeshGenerator3dQuality::optimizeTetrahedraForMinAngleMod(Metric3dContext& mc, 
			MeshContainer3d* mesh, double threshold)
{
	if(!mesh) return 0;

	START_CLOCK("MeshGenerator3dQuality::optimizeTetrahedraForMinAngleMod");

	int counter = 0;
	for(int i = 0; i < mesh->getBlocksCount(); i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*)mesh->getBlockAt(i);
		if(tetra->getType() != BLOCK_TETRA) continue; // tetrahedra only

		double min_angle = tetra->getMinDihedralAngleSin();
		if(min_angle >= threshold) continue;

		// 1. Try collapse one of edges
		//int ect = tetra->getEdgeCount();
		//for(int j = 0; j < ect; j++){
		//	MeshEdge3d* edge = tetra->getEdge(j);
		//	if(edge->isBorder()) continue;
		//	MeshPoint3d* point = MeshGenerator3dAdapt::collapseInnerEdge(mc, mesh, edge);
		//	if(point){
		//		DataVector<MeshBlock*> blocks(50);
		//		point->adjacentBlocks(blocks);
		//		for(int k = 0; k < blocks.countInt(); k++){
		//			MeshTetrahedron* t = (MeshTetrahedron*)blocks[k];
		//			if(t->getMinDihedralAngleSin() < threshold){
		//				if(MeshGenerator3dQuality::insertPointForTetrahedronRefinement(mc, mesh, t))
		//					break;
		//			}
		//		}
		//		break;
		//	}
		//}

		// 2. Then insert new point using Delaunay?
		MeshGenerator3dQuality::insertPointForTetrahedronRefinement(mc, mesh, tetra);

		counter++;
	}
	STOP_CLOCK("MeshGenerator3dQuality::optimizeTetrahedraForMinAngleMod");
	return counter;
}

/// Attempts to improve tetrahedra with angles below the given threshold, using simple optimization
int MeshGenerator3dQuality::optimizeTetrahedraForMinAngle(MeshContainer3d* mesh, double threshold)
{
	if(!mesh) return 0;

	START_CLOCK("MeshGenerator3dQuality::optimizeTetrahedraForMinAngle");

	int counter = 0;
	for(int i = 0; i < mesh->getBlocksCount(); i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*)mesh->getBlockAt(i);
		if(tetra->getType() != BLOCK_TETRA) continue; // tetrahedra only

		double min_angle = tetra->getMinDihedralAngleSin();
		if(min_angle >= threshold) continue;

		for(int j = 0; j < 4; j++){
			MeshPoint3d* point = tetra->getPoint(j);
			if(point->isBorder()) continue;
			movePointUsingAngleOptimization(point);
		}

		counter++;
	}
	STOP_CLOCK("MeshGenerator3dQuality::optimizeTetrahedraForMinAngle");
	return counter;
}

/// Generate (and log) statistics about mesh
void MeshGenerator3dQuality::statMesh(Metric3dContext& mc, MeshContainer3d* mesh, const string& caption)
{
	DataStatistics mean_ratio_stats;
	DataStatistics edge_stats;

	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "============ " << caption << " ============");

	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "NT\t" << mesh->getBlocksCount());
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "NP\t" << mesh->getPointsCount());

	mesh->statMeanRatio(mc, mean_ratio_stats, false);
	mesh->statMetricEdgeLength(mc, edge_stats, false);
	if(mean_ratio_stats.calculate()){
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricMeanRatio-ave\t" << mean_ratio_stats.average());
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricMeanRatio-min\t" << mean_ratio_stats.minimum());
	}
	if(edge_stats.calculate()){
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricEdgeLength-cnt\t" << edge_stats.countInt());
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricEdgeLength-ave\t" << edge_stats.average());
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricEdgeLength-max\t" << edge_stats.maximum());
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricEdgeLength-min\t" << edge_stats.minimum());
	}
	MeshGenerator3dQuality::statMinDihedralAngles(mesh);
}

/// Removes slivers on boundary surface
int MeshGenerator3dQuality::removeBoundarySlivers(Metric3dContext& mc, MeshContainer3d* mesh, 
	double mean_ratio_threshold, TagExtended::TagType tag_type, int tag_value)
{
	int bct = mesh->getBlocksCount();
	int removed_count = 0;
	for(int i = 0; i < bct; i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*)mesh->getBlockAt(i);
		assert(tetra->getType() == BLOCK_TETRA);
		if(tag_type != TagExtended::TAG_NONE && 
			!tetra->checkIntTagForAnyPoint(tag_type, tag_value)) continue;
		int tfct = tetra->getFaceCount(); assert(tfct == 4);
		DataVector<MeshFace*> bfaces(4);
		DataVector<MeshFace*> ifaces(4);
		for(int j = 0; j < tfct; j++){
			MeshFace* face = tetra->getFace(j);
			if(face->isBorder()) 
				bfaces.add(face);
			else
				ifaces.add(face);
		}
		int bfct = bfaces.countInt();
		if(bfct < 2) continue;
		bool valid = true;
		for(int j = 0; valid && j < bfct; j++)
			if(bfaces[j]->isBorder(TagBorder::FIXED)) valid = false;
		if(!valid) continue;

		DataVector<MeshEdge3d*> bedges(6);
		for(int j = 0; valid && j < bfct-1; j++){
			for(int k = j+1; k < bfct; k++){
				MeshEdge3d* edge = bfaces[j]->getCommonEdge(bfaces[k]);
				if(!edge || edge->isBorder( TagBorder::RIDGE | TagBorder::FIXED) ){
					valid = false; 
					break;
				}
				assert(edge->getFaceCount() == 2);
				bedges.add(edge);
			}
		}
		if(!valid) continue;

		int ifct = ifaces.countInt();
		for(int j = 0; valid && j < ifct-1; j++){
			for(int k = j+1; k < ifct; k++){
				MeshEdge3d* edge = ifaces[j]->getCommonEdge(ifaces[k]);
				if(!edge || edge->isBorder() ){
					valid = false; 
					break;
				}
			}
		}
		if(!valid) continue;

		double qt = tetra->countQuality(mc, MeshData::QUALITY3D_MEAN_RATIO);
		if(qt > mean_ratio_threshold) continue;


		//if(bfct > 2){
		//	double vol = tetra->getVolume(mc);
		//	ostringstream ostr;
		//	ostr << "Tetrahedron, vol=" << vol << ", qt=" << qt;
		//	MeshViewSet* set = new MeshViewSet;
		//	set->addBlockWithEdges(tetra, 2);
		//	for(int j = 0; j < 4; j++){
		//		MeshBlock* block = tetra->getNeighbour(j);
		//		if(block) set->addBlock(block);
		//	}
		//	SHOW_MESH(ostr.str(), set);
		//}

		char bf = bfaces[0]->getBorderFlags();

		for(int j = 0; j < ifct; j++){
			ifaces[j]->setBorderFlagsWithEdges(bf);
			ifaces[j]->copyAllTags(bfaces[0]);
		}
		for(int j = 0; j < bfct; j++){
			bfaces[j]->clearBorder();
		}

		for(int j = 0; j < bedges.countInt(); j++){
			bedges[j]->clearBorder();
		}

		MeshPoint3d* points[4] = {
			tetra->getPoint(0), tetra->getPoint(1),
			tetra->getPoint(2), tetra->getPoint(3)
		};

		delete mesh->removeMeshBlock(tetra);

		for(int j = 0; j < 4; j++){
			if(points[j]->getRank() == 0) 
				delete mesh->removeMeshPoint(points[j]);
		}

		i--;
		bct--;
		removed_count++;
	}


	return removed_count;
}

/// Prints statistical information about minimum inner angles (dihedral)
double MeshGenerator3dQuality::statMinDihedralAngles(MeshContainer3d* mesh)
{
	DataStatistics angle_stats;
	mesh->statMinDihedralAngles(angle_stats);
	if(angle_stats.calculate()){
		double angles_min = angle_stats.minimum() * 180.0 / PI;
		int angles_05 = angle_stats.getDataCountBottom(PI/36);
		int angles_10 = angle_stats.getDataCountBottom(PI/18);
		int angles_15 = angle_stats.getDataCountBottom(PI/12);
		int angles_20 = angle_stats.getDataCountBottom(PI/9);
		LOG4CPLUS_INFO(MeshLog::logger_console, "Angles (min) \t" << angles_min);
		LOG4CPLUS_INFO(MeshLog::logger_console, "------ (<05) \t" << angles_05);
		LOG4CPLUS_INFO(MeshLog::logger_console, "------ (<10) \t" << angles_10);
		LOG4CPLUS_INFO(MeshLog::logger_console, "------ (<15) \t" << angles_15);
		LOG4CPLUS_INFO(MeshLog::logger_console, "------ (<20) \t" << angles_20);
		int ct = angle_stats.countInt();
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Angles-below05o\t" << angles_05 << "\t(" << (100.0*angles_05 / ct) << "%)");
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Angles-below10o\t" << angles_10 << "\t(" << (100.0*angles_10 / ct) << "%)");
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Angles-below15o\t" << angles_15 << "\t(" << (100.0*angles_15 / ct) << "%)");
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Angles-below20o\t" << angles_20 << "\t(" << (100.0*angles_20 / ct) << "%)");

		return angles_min;
	}
	return -1.0;
}
