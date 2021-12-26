// MeshGenerator3d.cpp: implementation of the MeshGenerator3d class.
//
//////////////////////////////////////////////////////////////////////

//#define USE_OPENMP_HERE

#include <iomanip>
#include <algorithm>

#include "MeshGenerator3d.h"

#include "MeshContainer3d.h"
#include "MeshDomainVolume.h"
#include "MeshDomainSurface.h"
#include "MeshContainer2d.h"
#include "OctTree.h"
#include "MeshTetrahedron.h"
#include "MeshPoint2d.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshFace.h"
#include "MeshTriangle3d.h"
#include "MeshBlock.h"
#include "DataVector.h"
#include "DataList.h"
#include "DataHeapVector.h"
#include "MeshData.h"
#include "MeshLog.h"
#include "MeshBoundaryCondition.h"
#include "ControlSpace3d.h"
#include "ControlSpace3dOctree.h"
#include "ControlSpace2dIdentity.h"
#include "MeshGenerator2d.h"
#include "MeshGenerator3dQuality.h"
#include "MeshGenerator3dDelaunayBoundary.h"
#include "MeshViewSet.h"
#include "SurfacePlane.h"
#include "IteratorEdge3d.h"
#include "IteratorFace.h"
#include "ControlSpace2dAdaptive.h"
#include "GeometricPredicates.h"
#include "DTetrahedron.h"
#include "DTriangle.h"

#include "MeshDomainEdge3d.h"

int MeshGenerator3d::param_triangulation_type = MeshData::TRIANGULATE_CIRCLE;
double MeshGenerator3d::param_quality_improvement = 0.1;
double MeshGenerator3d::param_quality_threshold = 0.57;

/// Volume threshold during 3D generation (retriangulation, swapPIng, etc.)
double MeshGenerator3d::param_volume_threshold = 1e-5;

int MeshGenerator3d::param_max_auto_retriangulation_count = 1;

// Method used for boundary triangulation
// 0 - mixed, 1 - Delaunay only
int MeshGenerator3d::param_boundary_triangulation_method = 0;

//#define SHOW_DEL_SWAP
#define LOG_DEL_SWAP


#define IMPACT_SPHERE_FACTOR 2.0

//#define LOG_TRIANGULATION_CAVITY_STATS

int MeshGenerator3d::prepareBoundaryMesh(MeshContainer3d *boundary)
{
	START_CLOCK("MG3d::PrepareBoundary");
	int block_count = boundary->getBlocksCount();
	boundary->clearDiscretization3d();
	int finished_count = 0;
	for(int i = 0; i < block_count; i++){
		MeshDomainVolume* domain_volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		assert(domain_volume && (domain_volume->getType() == BLOCK_DOMAIN));
		assert(domain_volume->getMesh() == nullptr);
		int prepare_res = domain_volume->prepareBoundaryMesh();
		if(prepare_res > 0) ++finished_count;
		LOG_ASSERT(prepare_res > 0);
		//, MeshingException("Error preparing 3D boundary mesh"), 0);
	}
	STOP_CLOCK("MG3d::PrepareBoundary");
	return finished_count;
}

int MeshGenerator3d::prepareBoundarySurfaceMesh(MeshContainer3d *boundary)
{
	LOG_ASSERT(boundary != nullptr);
	START_CLOCK("MG3d::PrepareBoundarySurface");
	boundary->clearDiscretization3d();
	int total_fct = 0;
	// count
	for(IteratorFace it = boundary->getFirstFace(); it.isValid(); it.nextFace()){
		MeshDomainSurface* mds = (MeshDomainSurface*)it.getFace();
		MeshContainer2d* mds_mesh = mds->getMesh();
		if(mds_mesh)
			total_fct += mds_mesh->getElementsCount();
	}
/*
	// prepare
	MeshContainer3dSurface* mesh_surf = new MeshContainer3dSurface(total_fct);
	// fill
	for(IteratorFace it = boundary->getFirstFace(); it.isValid(); it.nextFace()){
		MeshDomainSurface* mds = (MeshDomainSurface*)it.getFace();
		MeshContainer2d* mds_mesh = mds->getMesh();
		if(!mds_mesh) continue;
		int pct = mds_mesh->getPointsCount();
		for(int i = 0; i < pct; i++){
			MeshPoint2d* mp = mds_mesh->getPointAt(i);
			MeshPoint3d* mp3d = mp->getPtrTag(TagExtended::TAG_
		}
	}

	int block_count = boundary->getBlocksCount();
	int finished_count = 0;
	for(int i = 0; i < block_count; i++){
		MeshDomainVolume* domain_volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		assert(domain_volume && (domain_volume->getType() == BLOCK_DOMAIN));
		assert(domain_volume->getMesh() == nullptr);
		if(domain_volume->prepareBoundaryMesh()) ++finished_count;
		else LOG4CPLUS_ASSERT(MeshLog::logger_mesh, false, MeshingException("Error preparing 3D boundary mesh"), 0);
	}
*/
	STOP_CLOCK("MG3d::PrepareBoundarySurface");
	return total_fct;
}

int MeshGenerator3d::createTetrahedra(MeshContainer3d *boundary)
{
	START_CLOCK("MG3d::CreateTetrahedralMesh");
	int block_count = boundary->getBlocksCount();
	int finished_count = 0;
	for(int i = 0; i < block_count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		assert(volume && (volume->getType() == BLOCK_DOMAIN));
		if(volume->createTetrahedralMesh())
			++finished_count;
	}
	STOP_CLOCK("MG3d::CreateTetrahedralMesh");
	return finished_count;
}

bool MeshGenerator3d::addPointToTriangulation(Metric3dContext& mc, MeshContainer3d *mesh, 
		MeshPoint3d* point, MeshTetrahedron *containing_tetrahedron,
		MeshTetrahedron* improved_tetrahedron, bool improving_only)
{
	assert(mesh->isValid());

	#pragma omp critical (mesh)
	{
		mesh->addMeshPoint(point);
	}

	bool result = false;

	switch(param_triangulation_type){
	case MeshData::TRIANGULATE_CIRCLE:
		result = addPointToTriangulationBySphere(mc, mesh, point, 
			containing_tetrahedron, improved_tetrahedron, improving_only);
		break;
	case MeshData::TRIANGULATE_ANGLE:
		// containing_tetrahedron == improved_tetrahedron, anyway...
		result = addPointToTriangulationBySwap(mc, mesh, point, containing_tetrahedron);
		break;
	default:
		assert(false);
	}

	if(!result){
		assert(point->getRank() == 0);
		#pragma omp critical (mesh)
		{
			mesh->removeMeshPoint(point->getIndex());
		}
	}

//	assert(mesh->isValid());

	return result;
}

#if defined(USE_OPENMP_HERE) && defined(_OPENMP)

bool MeshGenerator3d::addPointToTriangulationBySphere(Metric3dContext& mc, MeshContainer3d *mesh, 
		MeshPoint3d *point, MeshTetrahedron *containing_tetrahedron,
		MeshTetrahedron* improved_tetrahedron, bool improving_only)
{
	// Find and mark the triangle containing this point
	const DPoint3d dpt = point->getCoordinates();
	if(!containing_tetrahedron){
		MeshTetrahedron* start_tetrahedron = mesh->getNearTetrahedron(dpt);
		assert(start_tetrahedron);
		containing_tetrahedron = start_tetrahedron->findTetrahedronByNeighbours(dpt);
	}
	if(!containing_tetrahedron){
		LOG4CPLUS_WARN(MeshLog::logger_console, "switching to linear search for containing tetrahedron");
		#pragma omp critical (mesh)
		{
			int bct = mesh->getBlocksCount();
			for(int k = 0; k < bct; k++){
				MeshTetrahedron* tetrahedron = (MeshTetrahedron*)mesh->getBlockAt(k);
				if(tetrahedron->isPointInside(dpt)){
					containing_tetrahedron = tetrahedron;
					break;
				}
			}
		}
	}
	if(!containing_tetrahedron){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "containing tetrahedron not found");
		return false;
	}
	containing_tetrahedron->setIntTag(TagExtended::TAG_TRI_INSPHERE);
	DataCompoundList<MeshTetrahedron*> cavity_tetra;
	cavity_tetra.append(containing_tetrahedron);
	int area_id = containing_tetrahedron->getAreaID();

	// Find all neighbours, having this point in the circumscribed sphere
	for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
		for(int j = 0; j < 4; j++){	// For all neighbours
			MeshTetrahedron* tetrahedron = (MeshTetrahedron*)it.item()->getNeighbour(j);
			if(tetrahedron){
				if(!tetrahedron->availableTag(TagExtended::TAG_TRI_INSPHERE)){
					bool allowed = mc.withinImpactRadius(tetrahedron->getMiddlePoint());
					if(mesh->getConstrainingPhase() != MeshData::CONSTRAIN_NONE){
						allowed = !it.item()->isBorder(j);
					}
					if(allowed && mesh->getConstrainingPhase() != MeshData::CONSTRAIN_NONE){
					//if(allowed && mesh->getConstrainingPhase() == MeshData::CONSTRAIN_ACTIVE){
						for(int k = 0; k < 4; k++){
							if(tetrahedron->isBorder(k)){
								MeshTetrahedron* t = (MeshTetrahedron*)tetrahedron->getNeighbour(k);
								if(t && t->availableTag(TagExtended::TAG_TRI_INSPHERE))
								{	// if is incident to other tetrahedron of core, but through boundary face
									//LOG4CPLUS_INFO(MeshLog::logger_mesh, "** stranded boundary face while gathering core - fixing.");
									allowed = false;
									break;
								}
							}
						}
						// check for stranded boundary edges:
						for(int k = 0; allowed && (k < 6); k++){
							MeshEdge3d* edge = tetrahedron->getEdge(k);
							// * boundary edge
							if(edge->isBorder()){
								int edge_rank = edge->getFaceCount();
								allowed = false;
								// * all incident faces are normal
								// * and all neighbouring blocks (with this) would be tagged
								for(int l = 0; l < edge_rank; l++){
									MeshFace* face = edge->getFaceAt(l);
									if(face->isBorder()){
										allowed = true;	// this edge is safe - there is incident boundary face
										break;
									}else{
										MeshBlock* block = face->getBlock(0);
										if(block != tetrahedron && !block->availableTag(TagExtended::TAG_TRI_INSPHERE)){
											allowed = true; // this edge is safe - there is other not-tagged block
											break;
										}
										block = face->getBlock(1);
										if(block != tetrahedron && !block->availableTag(TagExtended::TAG_TRI_INSPHERE)){
											allowed = true; // this edge is safe - there is other not-tagged block
											break;
										}
									}
								}
//								if(!allowed){
//									LOG4CPLUS_INFO(MeshLog::logger_mesh, "** stranded boundary edge while gathering core - fixing.");
//								}
							}
						}
					}
#ifdef LOG_TRIANGULATION_CAVITY_STATS
					++checked_count;
#endif
					if(allowed && tetrahedron->isPointInOuterSphere(mc, point, false)){
						tetrahedron->setIntTag(TagExtended::TAG_TRI_INSPHERE);
						// Add to cavity-set
						cavity_tetra.append(tetrahedron);
					}
				}
			}
		}
	}

//	bool refine_mesh = (MeshGenerator3d::param_refined_ICD_ratio > 0.0) && 
//		(mesh->getConstrainingPhase() != MeshData::CONSTRAIN_DONE);

#ifdef LOG_TRIANGULATION_CAVITY_STATS
	int cavity_count = cavity_tetra.countInt();
#endif

	// Set faces of the core
	DataCompoundList<TFace> faces;
	double min_volume = 0.001;
	int step = 0;
	do{
		bool inconsistent_core;	
		do{
			faces.clear();
			inconsistent_core = false;

			for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
				if(!it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE)) continue;
				for(int j = 0; j < 4; j++){
					MeshTetrahedron* tetrahedron = (MeshTetrahedron*)it.item()->getNeighbour(j);
					if(!tetrahedron || !tetrahedron->availableTag(TagExtended::TAG_TRI_INSPHERE)){
						// Boundary face or existing not-tagged neighbour
						TFace tf;
						tf.face = it.item()->getFace(j);
						tf.face_index = tf.face->getBlockIndex(it.item());
						// Check orientation of the face
						MeshPoint3d* fourth_point = it.item()->getPoint(j);
						assert(!tf.face->incidentToPoint(fourth_point));
						if(tf.face->sameSide(fourth_point->getCoordinates(), 
							point->getCoordinates()))
						{
							double vol = DTetrahedron::volume(
								tf.face->getPoint(0)->getMetricCoordinates(mc),
								tf.face->getPoint(1)->getMetricCoordinates(mc),
								tf.face->getPoint(2)->getMetricCoordinates(mc),
								point->getMetricCoordinates(mc));
							if(abs(vol) < min_volume){
								LOG4CPLUS_WARN(MeshLog::logger_console, ("addPointBySphere:skipPIng face - small volume", vol);
								it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
								inconsistent_core = true;	// mark global repeat
								break;
							}else
								faces.append(tf);
						}else{
							//LOG4CPLUS_INFO(MeshLog::logger_mesh, "** core face visibility failed - fixing.");
							it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
							inconsistent_core = true;	// mark global repeat
							break;
						}
					}
				}
			}
	//		LOG4CPLUS_INFO(MeshLog::logger_console, " --> cavity faces", faces.countInt());
		}while(inconsistent_core);

		if(improving_only){
			// check proximity
			for (auto itf = faces.iterator(); itf.valid(); itf.moveNext()) {
				TFace& tf = itf.item();
				for(int i = 0; i < 3; i++){
					MeshPoint3d* fpoint = tf.face->getPoint(i);
					double dist2 = fpoint->getMetricCoordinates(mc).distance2(
						point->getMetricCoordinates(mc));
					if(dist2 < 0.25){
						//LOG4CPLUS_WARN(MeshLog::logger_console, "Skip improve-point near to other point", sqrt(dist2));
						for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
							if (it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE))
								it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
						}
						return false;
					}
				}
				// check distance from boundary faces (if any)
				if(tf.face->isBorder()){
					double dist2 = DTriangle3d::distance2ToPoint(
						point->getMetricCoordinates(mc),
						tf.face->getPoint(0)->getMetricCoordinates(mc),
						tf.face->getPoint(1)->getMetricCoordinates(mc),
						tf.face->getPoint(2)->getMetricCoordinates(mc));
					if(dist2 < 0.04){
						//LOG4CPLUS_WARN(MeshLog::logger_console, "Skip improve-point near to boundary face", sqrt(dist2));
						for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
							if (it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE))
								it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
						}
						return false;
					}
				}
			}
		}

		if(improving_only || ++step > 2) break;
		if(faces.empty()){
			LOG4CPLUS_WARN(MeshLog::logger_console, ("addPointBySphere, step", step);
			min_volume *= 0.01;
			for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
				it.item()->setIntTag(TagExtended::TAG_TRI_INSPHERE);
			}
		}

	}while(faces.empty());

//	LOG4CPLUS_INFO(MeshLog::logger_console, "==============================");
//	LOG4CPLUS_INFO(MeshLog::logger_console, " --> cavity tetrahedra", cavity_tetra.countInt());
//	LOG4CPLUS_INFO(MeshLog::logger_console, " --> cavity faces", faces.countInt());

	if(improved_tetrahedron){
		bool improved_tetrahedron_in_cavity = false;
		for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
			if(it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE) &&
				(it.item() == improved_tetrahedron))
			{
				improved_tetrahedron_in_cavity = true;
				break;
			}
		}
		if(!improved_tetrahedron_in_cavity){
			LOG4CPLUS_WARN(MeshLog::logger_console, "improved tetrahedron not in cavity");
			for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
				if(it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE))
					it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
			return false;
		}
	}

#ifdef LOG_TRIANGULATION_CAVITY_STATS
	int final_cavity_count = 0;
#endif

	// Remove obsolete tetrahedra
	#pragma omp critical (mesh)
	{
		for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
			if(it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE)){
				delete mesh->removeMeshTetrahedron(it.item());
#ifdef LOG_TRIANGULATION_CAVITY_STATS
				++final_cavity_count;
#endif
			}
		}
	}

//	SHOW_STEP_PT_BREAKABLE(3, "Pusta wnêka.", point->getCoordinates(), false);
	
	if(faces.countInt() < 3){
		LOG4CPLUS_WARN(MeshLog::logger_console, "too few faces in cavity", faces.countInt());
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "too few faces in cavity");

//		MeshViewSet *set = new MeshViewSet;
//		set->addEmptyBlockWithEdges(containing_tetrahedron, 0);
//		set->addPoint(point);
//		SHOW_MESH("null cavity?", set);

		if(improving_only) return false;
		else return addPointToTriangulationBySwap(mc, mesh, point, containing_tetrahedron, false);
	}

#ifdef _DEBUG
	for (auto itf = faces.iterator(); itf.valid(); itf.moveNext()) {
		assert(!(itf.item().face->getBlock(0) && itf.item().face->getBlock(1)));
	}
#endif // _DEBUG
	
	// Create new tetrahedra
	for (auto itf = faces.iterator(); itf.valid(); itf.moveNext()) {
		TFace& tf = itf.item();
		if(tf.face_index == 0){
			tf.tetrahedron = new MeshTetrahedron(tf.face->getPoint(0), 
				tf.face->getPoint(1), tf.face->getPoint(2), point);
		}else{
			tf.tetrahedron = new MeshTetrahedron(tf.face->getPoint(0), 
				tf.face->getPoint(2), tf.face->getPoint(1), point);
		}
		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE){
			tf.tetrahedron->setAreaID(area_id);
			tf.tetrahedron->countQuality(mc);
		}
		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_ACTIVE){
//		if(true){
			double vol = tf.tetrahedron->getVolume(mc, false);
			if(vol < METRIC_SMALL_NUMBER){
				LOG4CPLUS_WARN(MeshLog::logger_console, "addPoint - creating very small volume", vol);
			}
		}
		#pragma omp critical (mesh)
		{
			mesh->addMeshTetrahedron(tf.tetrahedron);
		}
	}

#ifdef LOG_TRIANGULATION_CAVITY_STATS
	if(mesh->getConstrainingPhase() == 0){ // boundary
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "stat-traverse: " << mesh->getPointsCount() << '\t';
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, MeshTetrahedron::traverse_counter << '\t';
		MeshTetrahedron::traverse_counter = 0;
		LOG4CPLUS_INFO(MeshLog::logger_mesh, mesh->getMaxSearchTreeLevel());
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "stat-cavity: " << mesh->getPointsCount() << '\t' << mesh->getConstrainingPhase() << '\t';
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, (mc.getCallCounter() - previous_cs_calls) << '\t';
	previous_cs_calls = mc.getCallCounter();
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, checked_count << '\t' << cavity_count << '\t';
	LOG4CPLUS_INFO(MeshLog::logger_mesh, final_cavity_count << '\t' << faces.countInt());
#endif

	return true;
}
#else
bool MeshGenerator3d::addPointToTriangulationBySphere(Metric3dContext& mc, MeshContainer3d *mesh, 
		MeshPoint3d *point, MeshTetrahedron *containing_tetrahedron,
		MeshTetrahedron* improved_tetrahedron, bool improving_only)
{
	// Find and mark the triangle containing this point
	const DPoint3d dpt = point->getCoordinates();
	if(!containing_tetrahedron){
		MeshTetrahedron* start_tetrahedron = mesh->getNearTetrahedron(dpt);
		assert(start_tetrahedron);
		containing_tetrahedron = start_tetrahedron->findTetrahedronByNeighbours(dpt);
	}
	if(!containing_tetrahedron){
		LOG4CPLUS_WARN(MeshLog::logger_console, "switching to linear search for containing tetrahedron");
		#pragma omp critical (mesh)
		{
			int bct = mesh->getBlocksCount();
			for(int k = 0; k < bct; k++){
				MeshTetrahedron* tetrahedron = (MeshTetrahedron*)mesh->getBlockAt(k);
				if(tetrahedron->isPointInside(dpt)){
					containing_tetrahedron = tetrahedron;
					break;
				}
			}
		}
	}
	if(!containing_tetrahedron){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "containing tetrahedron not found");
		return false;
	}
	containing_tetrahedron->setIntTag(TagExtended::TAG_TRI_INSPHERE);
	DataCompoundList<MeshTetrahedron*> cavity_tetra;
	cavity_tetra.append(containing_tetrahedron);
	int area_id = containing_tetrahedron->getAreaID();

#ifdef LOG_TRIANGULATION_CAVITY_STATS
	int checked_count = 0;
	static unsigned int previous_cs_calls = mc.getCallCounter();
#endif

	// Find all neighbours, having this point in the circumscribed sphere
	for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
		for(int j = 0; j < 4; j++){	// For all neighbours
			MeshTetrahedron* tetrahedron = (MeshTetrahedron*)it.item()->getNeighbour(j);
			if(tetrahedron){
				if(!tetrahedron->availableTag(TagExtended::TAG_TRI_INSPHERE)){
					bool allowed = true;
					if(mesh->getConstrainingPhase() != MeshData::CONSTRAIN_NONE){
						allowed = !it.item()->isBorder(j);
					}
					if(allowed && mesh->getConstrainingPhase() != MeshData::CONSTRAIN_NONE){
					//if(allowed && mesh->getConstrainingPhase() == MeshData::CONSTRAIN_ACTIVE){
						for(int k = 0; k < 4; k++){
							if(tetrahedron->isBorder(k)){
								MeshTetrahedron* t = (MeshTetrahedron*)tetrahedron->getNeighbour(k);
								if(t && t->availableTag(TagExtended::TAG_TRI_INSPHERE))
								{	// if is incident to other tetrahedron of core, but through boundary face
									//LOG4CPLUS_INFO(MeshLog::logger_mesh, "** stranded boundary face while gathering core - fixing.");
									allowed = false;
									break;
								}
							}
						}
						// check for stranded boundary edges:
						for(int k = 0; allowed && (k < 6); k++){
							MeshEdge3d* edge = tetrahedron->getEdge(k);
							// * boundary edge
							if(edge->isBorder()){
								int edge_rank = edge->getFaceCount();
								allowed = false;
								// * all incident faces are normal
								// * and all neighbouring blocks (with this) would be tagged
								for(int l = 0; l < edge_rank; l++){
									MeshFace* face = edge->getFaceAt(l);
									if(face->isBorder()){
										allowed = true;	// this edge is safe - there is incident boundary face
										break;
									}else{
										MeshBlock* block = face->getBlock(0);
										if(block != tetrahedron && !block->availableTag(TagExtended::TAG_TRI_INSPHERE)){
											allowed = true; // this edge is safe - there is other not-tagged block
											break;
										}
										block = face->getBlock(1);
										if(block != tetrahedron && !block->availableTag(TagExtended::TAG_TRI_INSPHERE)){
											allowed = true; // this edge is safe - there is other not-tagged block
											break;
										}
									}
								}
//								if(!allowed){
//									LOG4CPLUS_INFO(MeshLog::logger_mesh, "** stranded boundary edge while gathering core - fixing.");
//								}
							}
						}
					}
#ifdef LOG_TRIANGULATION_CAVITY_STATS
					++checked_count;
#endif
					if(allowed && tetrahedron->isPointInOuterSphere(mc, point, false)){
						tetrahedron->setIntTag(TagExtended::TAG_TRI_INSPHERE);
						// Add to cavity-set
						cavity_tetra.append(tetrahedron);
					}
				}
			}
		}
	}

//	bool refine_mesh = (MeshGenerator3d::param_refined_ICD_ratio > 0.0) && 
//		(mesh->getConstrainingPhase() != MeshData::CONSTRAIN_DONE);

#ifdef LOG_TRIANGULATION_CAVITY_STATS
	int cavity_count = cavity_tetra.countInt();
#endif

	// Set faces of the core
	DataCompoundList<TFace> faces;
	double min_volume = 0.001;
	int step = 0;
	do{
		bool inconsistent_core;	
		do{
			faces.clear();
			inconsistent_core = false;

			for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
				if(!it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE)) continue;
				for(int j = 0; j < 4; j++){
					MeshTetrahedron* tetrahedron = (MeshTetrahedron*)it.item()->getNeighbour(j);
					if(!tetrahedron || !tetrahedron->availableTag(TagExtended::TAG_TRI_INSPHERE)){
						// Boundary face or existing not-tagged neighbour
						TFace tf;
						tf.face = it.item()->getFace(j);
						tf.face_index = tf.face->getBlockIndex(it.item());
						// Check orientation of the face
						MeshPoint3d* fourth_point = it.item()->getPoint(j);
						assert(!tf.face->incidentToPoint(fourth_point));
						if(tf.face->sameSide(fourth_point->getCoordinates(), 
							point->getCoordinates()))
						{
							double vol = DTetrahedron::volume(
								tf.face->getPoint(0)->getMetricCoordinates(mc),
								tf.face->getPoint(1)->getMetricCoordinates(mc),
								tf.face->getPoint(2)->getMetricCoordinates(mc),
								point->getMetricCoordinates(mc));
							if(abs(vol) < min_volume){
								// LOG4CPLUS_WARN(MeshLog::logger_console, ("addPointBySphere:skipPIng face - small volume", vol);
								it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
								inconsistent_core = true;	// mark global repeat
								break;
							}else
								faces.append(tf);
						}else{
							//LOG4CPLUS_INFO(MeshLog::logger_mesh, "** core face visibility failed - fixing.");
							it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
							inconsistent_core = true;	// mark global repeat
							break;
						}
					}
				}
			}
	//		LOG4CPLUS_INFO(MeshLog::logger_console, " --> cavity faces", faces.countInt());
		}while(inconsistent_core);

		if(improving_only){
			// check proximity
			for (auto itf = faces.iterator(); itf.valid(); itf.moveNext()) {
				TFace& tf = itf.item();
				for(int i = 0; i < 3; i++){
					MeshPoint3d* fpoint = tf.face->getPoint(i);
					double dist2 = fpoint->getMetricCoordinates(mc).distance2(
						point->getMetricCoordinates(mc));
					if(dist2 < 0.25){
						//LOG4CPLUS_WARN(MeshLog::logger_console, "Skip improve-point near to other point", sqrt(dist2));
						for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
							if (it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE))
								it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
						}
						return false;
					}
				}
				// check distance from boundary faces (if any)
				if(tf.face->isBorder()){
					double dist2 = DMTriangle3d::distance2ToPoint(
						point->getMetricCoordinates(mc),
						tf.face->getPoint(0)->getMetricCoordinates(mc),
						tf.face->getPoint(1)->getMetricCoordinates(mc),
						tf.face->getPoint(2)->getMetricCoordinates(mc));
					if(dist2 < 0.04){
						//LOG4CPLUS_WARN(MeshLog::logger_console, "Skip improve-point near to boundary face", sqrt(dist2));
						for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
							if (it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE))
								it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
						}
						return false;
					}
				}
			}
		}

		if(improving_only || ++step > 2) break;
		if(faces.empty()){
			LOG4CPLUS_WARN(MeshLog::logger_console, "addPointBySphere, step " << step);
			min_volume *= 0.01;
			for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
				it.item()->setIntTag(TagExtended::TAG_TRI_INSPHERE);
			}
		}

	}while(faces.empty());

//	LOG4CPLUS_INFO(MeshLog::logger_console, "==============================");
//	LOG4CPLUS_INFO(MeshLog::logger_console, " --> cavity tetrahedra", cavity_tetra.countInt());
//	LOG4CPLUS_INFO(MeshLog::logger_console, " --> cavity faces", faces.countInt());

	if(improved_tetrahedron){
		bool improved_tetrahedron_in_cavity = false;
		for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
			if(it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE) &&
				(it.item() == improved_tetrahedron))
			{
				improved_tetrahedron_in_cavity = true;
				break;
			}
		}
		if(!improved_tetrahedron_in_cavity){
			//LOG4CPLUS_WARN(MeshLog::logger_console, "improved tetrahedron not in cavity");
			for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
				if (it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE))
					it.item()->removeTag(TagExtended::TAG_TRI_INSPHERE);
			}
			return false;
		}
	}

#ifdef LOG_TRIANGULATION_CAVITY_STATS
	int final_cavity_count = 0;
#endif

	// Remove obsolete tetrahedra
	for (auto it = cavity_tetra.iterator(); it.valid(); it.moveNext()) {
		if (it.item()->availableTag(TagExtended::TAG_TRI_INSPHERE)) {
			delete mesh->removeMeshTetrahedron(it.item());
#ifdef LOG_TRIANGULATION_CAVITY_STATS
			++final_cavity_count;
#endif
		}
	}

//	SHOW_STEP_PT_BREAKABLE(3, "Pusta wnêka.", point->getCoordinates(), false);
	
	if(faces.countInt() < 3){
		LOG4CPLUS_WARN(MeshLog::logger_console, "too few faces in cavity -> " << faces.countInt());
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "too few faces in cavity");

//		MeshViewSet *set = new MeshViewSet;
//		set->addEmptyBlockWithEdges(containing_tetrahedron, 0);
//		set->addPoint(point);
//		SHOW_MESH("null cavity?", set);

		if(improving_only) return false;
		else return addPointToTriangulationBySwap(mc, mesh, point, containing_tetrahedron, false);
	}

#ifdef _DEBUG
	for (auto itf = faces.iterator(); itf.valid(); itf.moveNext()) {
		assert(!(itf.item().face->getBlock(0) && itf.item().face->getBlock(1)));
	}
#endif // _DEBUG
	
	// Create new tetrahedra
	for (auto itf = faces.iterator(); itf.valid(); itf.moveNext()) {
		TFace& tf = itf.item();
		if(tf.face_index == 0){
			tf.tetrahedron = new MeshTetrahedron(tf.face->getPoint(0), 
				tf.face->getPoint(1), tf.face->getPoint(2), point);
		}else{
			tf.tetrahedron = new MeshTetrahedron(tf.face->getPoint(0), 
				tf.face->getPoint(2), tf.face->getPoint(1), point);
		}
		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE){
			tf.tetrahedron->setAreaID(area_id);
			tf.tetrahedron->countQuality(mc);
		}
		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_ACTIVE){
//		if(true){
			double vol = tf.tetrahedron->getVolume(mc, false);
			if(vol < METRIC_SMALL_NUMBER){
				LOG4CPLUS_WARN(MeshLog::logger_console, "addPoint - creating very small volume, vol=" << vol);
			}
		}
		mesh->addMeshTetrahedron(tf.tetrahedron);
	}

#ifdef LOG_TRIANGULATION_CAVITY_STATS
	if(mesh->getConstrainingPhase() == 0){ // boundary
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "stat-traverse: " << mesh->getPointsCount() << '\t';
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, MeshTetrahedron::traverse_counter << '\t';
		MeshTetrahedron::traverse_counter = 0;
		LOG4CPLUS_INFO(MeshLog::logger_mesh, mesh->getMaxSearchTreeLevel());
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "stat-cavity: " << mesh->getPointsCount() << '\t' << mesh->getConstrainingPhase() << '\t';
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, (mc.getCallCounter() - previous_cs_calls) << '\t';
	previous_cs_calls = mc.getCallCounter();
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, checked_count << '\t' << cavity_count << '\t';
	LOG4CPLUS_INFO(MeshLog::logger_mesh, final_cavity_count << '\t' << faces.countInt());
#endif

	return true;
}
#endif

#if defined(USE_OPENMP_HERE) && defined(_OPENMP)

bool MeshGenerator3d::addPointToTriangulationBySwap(Metric3dContext& mc, MeshContainer3d *mesh, 
							MeshPoint3d *point, MeshTetrahedron *containing_tetrahedron, 
							bool allow_swap)
{
	const DPoint3d& pt = point->getCoordinates();
	if(!containing_tetrahedron){
		MeshTetrahedron* start_tetrahedron = mesh->getNearTetrahedron(pt);
		assert(start_tetrahedron);
		containing_tetrahedron = start_tetrahedron->findTetrahedronByNeighbours(pt);
	}
	if(!containing_tetrahedron){
		LOG4CPLUS_WARN(MeshLog::logger_console, "switching to linear search for containing tetrahedron");
		#pragma omp critical (mesh)
		{
			int bct = mesh->getBlocksCount();
			for(int k = 0; k < bct; k++){
				MeshTetrahedron* tetrahedron = (MeshTetrahedron*)mesh->getBlockAt(k);
				if(tetrahedron->isPointInside(pt)){
					containing_tetrahedron = tetrahedron;
					break;
				}
			}
		}
	}
	if(!containing_tetrahedron) return false;

	// split tetrahedron into four blocks
	DataCompoundList<MeshTetrahedron*> new_tetrahedra;
	int area_id1 = containing_tetrahedron->getAreaID();
	MeshTetrahedron *tetrahedron = containing_tetrahedron;


	MeshPoint3d* points[4] = {
		tetrahedron->getPoint(0), tetrahedron->getPoint(1),
		tetrahedron->getPoint(2), tetrahedron->getPoint(3)};

	//mc.countMetricAtPoint(pt);

	// check for edge-collision
	for(int i = 0; i < 6; i++){
		int i1, i2;
		tetrahedron->getEdgeIndices(i, i1, i2);
		MeshEdge3d* edge = tetrahedron->getEdge(i);
		if(!edge->isBorder()){
			const DMPoint3d dpt1 = points[i1]->getMetricCoordinates(mc);
			const DMPoint3d dpt2 = points[i2]->getMetricCoordinates(mc);
			const DMPoint3d dpt0 = point->getMetricCoordinates(mc);
			double d12 = dpt1.distance(dpt2);
			double d01 = dpt0.distance(dpt1);
			double d02 = dpt0.distance(dpt2);
			double h = 2*METRIC_SMALL_NUMBER;
			if(d01 < d12 && d02 < d12){
				double p = 0.5 * (d01 + d02 + d12);
				double s = sqrt(p * (p - d01) * (p - d02) * (p - d12));
				h = (2 * s) / d12;
			}
			if(h < METRIC_SMALL_NUMBER){
				// point being inserted "at" the edge
				MeshPoint3d* mpt1 = points[i1];
				MeshPoint3d* mpt2 = points[i2];
#ifdef LOG_DEL_SWAP
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "Del-by-swap: inserting at edge");
#endif
#ifdef SHOW_DEL_SWAP
				MeshView::showDebugMesh("edge-collision", mesh, mpt1, mpt2);
#endif
				while(true){
					int fct = edge->getFaceCount();
					MeshFace* face = edge->getFaceAt(0);
					assert(face->isBounded() && !face->isBorder());
					tetrahedron = (MeshTetrahedron*)face->getBlock(0);
					if(!tetrahedron) tetrahedron = (MeshTetrahedron*)face->getBlock(1);
					assert(tetrahedron);
					assert(!tetrahedron->isInverted());
					MeshPoint3d *tpoints[4] = { tetrahedron->getPoint(0), tetrahedron->getPoint(1),
						tetrahedron->getPoint(2), tetrahedron->getPoint(3)};

					#pragma omp critical (mesh)
					{
						delete mesh->removeMeshTetrahedron(tetrahedron);

						MeshPoint3d *epoints[2] = { mpt1, mpt2 };
						for(int j = 0; j < 2; j++){
							mesh->addMeshTetrahedron(tetrahedron = new MeshTetrahedron(
								(tpoints[0] == epoints[j]) ? point : tpoints[0],
								(tpoints[1] == epoints[j]) ? point : tpoints[1],
								(tpoints[2] == epoints[j]) ? point : tpoints[2],
								(tpoints[3] == epoints[j]) ? point : tpoints[3]));
							assert(!tetrahedron->isInverted());
							new_tetrahedra.append(tetrahedron);
						}
					}

					if(fct == 2) break;	// it was last tetrahedron to split
				}
				break; // the for loop
			}
		}else{
			// point can't be inserted at the border edge
		}
	}
	if(new_tetrahedra.countInt() == 0){
		// check for face-collision
		for(int i = 0; i < 4; i++){
			int i1, i2, i3;
			tetrahedron->getFaceIndices(i, i1, i2, i3);
			MeshFace* face = tetrahedron->getFace(i);
			if(!face->isBorder() && 
				DTriangle3d::distance2ToPoint(
					point->getMetricCoordinates(mc), 
					points[i1]->getMetricCoordinates(mc), 
					points[i2]->getMetricCoordinates(mc),
					points[i3]->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER)
			{
#ifdef LOG_DEL_SWAP
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "Del-by-swap: inserting at face");
#endif
/*
				if(true){
					MeshViewSet *set = new MeshViewSet;
					set->addFace(face, 0);
					set->addPoint(point, 1);
					SHOW_MESH("face-collision", set);
				}
*/
				MeshTetrahedron* tetra[2] = { tetrahedron, 
					(MeshTetrahedron*)face->getOtherBlock(tetrahedron) };
/*
				if(true){
					MeshViewSet *set = new MeshViewSet;
					set->addFace(face, 0);
					set->addPoint(point, 1);
					set->addEmptyBlockWithEdges(tetra[0]);
					set->addEmptyBlockWithEdges(tetra[1]);
					SHOW_MESH("face-collision", set);
				}
*/
				MeshPoint3d* fpoints[3] = { face->getPoint(0),
					face->getPoint(1), face->getPoint(2) };
				for(int ti = 0; ti < 2; ti++){
					MeshPoint3d* tpoints[4] = { 
						tetra[ti]->getPoint(0), tetra[ti]->getPoint(1), 
						tetra[ti]->getPoint(2), tetra[ti]->getPoint(3) };
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "** tetra volume = " << tetra[ti]->getVolumeNoMetric());

					#pragma omp critical (mesh)
					{
						delete mesh->removeMeshTetrahedron(tetra[ti]);
						for(int fi = 0; fi < 3; fi++){
							mesh->addMeshTetrahedron(tetrahedron = new MeshTetrahedron(
								(tpoints[0]==fpoints[fi]) ? point : tpoints[0],
								(tpoints[1]==fpoints[fi]) ? point : tpoints[1],
								(tpoints[2]==fpoints[fi]) ? point : tpoints[2],
								(tpoints[3]==fpoints[fi]) ? point : tpoints[3]));
							assert(!tetrahedron->isInverted());
							LOG4CPLUS_INFO(MeshLog::logger_mesh, "*** tetra volume = " << tetrahedron->getVolumeNoMetric());
							new_tetrahedra.append(tetrahedron);
						}
					}
				}
				break; // the for loop
			}
		}
	}
	if(new_tetrahedra.countInt() == 0){
		// insert point within the tetrahedron
		MeshPoint3d *tpoints[4] = { tetrahedron->getPoint(0), tetrahedron->getPoint(1),
			tetrahedron->getPoint(2), tetrahedron->getPoint(3)};

		#pragma omp critical (mesh)
		{
			delete mesh->removeMeshTetrahedron(tetrahedron);

			for(int j = 0; j < 4; j++){
				mesh->addMeshTetrahedron(tetrahedron = new MeshTetrahedron(
					(j == 0) ? point : tpoints[0],
					(j == 1) ? point : tpoints[1],
					(j == 2) ? point : tpoints[2],
					(j == 3) ? point : tpoints[3]));
				assert(!tetrahedron->isInverted());
				new_tetrahedra.append(tetrahedron);
			}
		}

#ifdef LOG_DEL_SWAP
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Del-by-swap: inserting within the block");
#endif
/*
		if(true){
			MeshViewSet *set = new MeshViewSet;
			for(int i = 0; i < 4; i++){
				set->addEmptyBlockWithEdges(new_tetrahedrons.get(i), 0);
				set->addPoint(tpoints[i], 0);
			}
			set->addPoint(point, 1);
			SHOW_MESH("tetrahedra-split", set);
		}
*/
#ifdef SHOW_DEL_SWAP
		MeshViewSet *set = new MeshViewSet;
		for(int i = 0; i < 4; i++)
			set->addBlockWithEdges(new_tetrahedrons.get(i), 0);
		set->addPoint(point);
		SHOW_MESH("tetrahedra-split", set);
#endif
	}
#ifdef LOG_DEL_SWAP
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Del-by-swap: new tetrahedra=" << new_tetrahedra.countInt());
#endif
	assert(new_tetrahedra.countInt() > 0);
	if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE){
		for (auto it = new_tetrahedra.iterator(); it.valid(); it.moveNext()) {
			it.item()->setAreaID(area_id1);
		}
	}

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Del-by-swap: split, mesh valid " << mesh->isValid());

	// swap
	if(allow_swap) 
		iterativeTetrahedraSwap(mc, mesh, point, 5, 100, true, 
			mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE);

	return true;
}
#else
bool MeshGenerator3d::addPointToTriangulationBySwap(Metric3dContext& mc, MeshContainer3d *mesh, 
							MeshPoint3d *point, MeshTetrahedron *containing_tetrahedron, 
							bool allow_swap)
{
	const DPoint3d& pt = point->getCoordinates();
	if(!containing_tetrahedron){
		MeshTetrahedron* start_tetrahedron = mesh->getNearTetrahedron(pt);
		assert(start_tetrahedron);
		containing_tetrahedron = start_tetrahedron->findTetrahedronByNeighbours(pt);
	}
	if(!containing_tetrahedron){
		LOG4CPLUS_WARN(MeshLog::logger_console, "switching to linear search for containing tetrahedron");
		int bct = mesh->getBlocksCount();
		for(int k = 0; k < bct; k++){
			MeshTetrahedron* tetrahedron = (MeshTetrahedron*)mesh->getBlockAt(k);
			if(tetrahedron->isPointInside(pt)){
				containing_tetrahedron = tetrahedron;
				break;
			}
		}
	}
	if(!containing_tetrahedron) return false;

	// split tetrahedron into four blocks
	DataCompoundList<MeshTetrahedron*> new_tetrahedra;
	int area_id1 = containing_tetrahedron->getAreaID();
	MeshTetrahedron *tetrahedron = containing_tetrahedron;


	MeshPoint3d* points[4] = {
		tetrahedron->getPoint(0), tetrahedron->getPoint(1),
		tetrahedron->getPoint(2), tetrahedron->getPoint(3)};

	//mc.countMetricAtPoint(pt);

	// check for edge-collision
	for(int i = 0; i < 6; i++){
		int i1, i2;
		tetrahedron->getEdgeIndices(i, i1, i2);
		MeshEdge3d* edge = tetrahedron->getEdge(i);
		if(!edge->isBorder()){
			const DMPoint3d dpt1 = points[i1]->getMetricCoordinates(mc);
			const DMPoint3d dpt2 = points[i2]->getMetricCoordinates(mc);
			const DMPoint3d dpt0 = point->getMetricCoordinates(mc);
			double d12 = dpt1.distance(dpt2);
			double d01 = dpt0.distance(dpt1);
			double d02 = dpt0.distance(dpt2);
			double h = 2*METRIC_SMALL_NUMBER;
			if(d01 < d12 && d02 < d12){
				double p = 0.5 * (d01 + d02 + d12);
				double s = sqrt(p * (p - d01) * (p - d02) * (p - d12));
				h = (2 * s) / d12;
			}
			if(h < METRIC_SMALL_NUMBER){
				// point being inserted "at" the edge
				MeshPoint3d* mpt1 = points[i1];
				MeshPoint3d* mpt2 = points[i2];
#ifdef LOG_DEL_SWAP
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Del-by-swap: inserting at edge");
#endif
#ifdef SHOW_DEL_SWAP
				MeshView::showDebugMesh("edge-collision", mesh, mpt1, mpt2);
#endif
				while(true){
					int fct = edge->getFaceCount();
					MeshFace* face = edge->getFaceAt(0);
					assert(face->isBounded() && !face->isBorder());
					tetrahedron = (MeshTetrahedron*)face->getBlock(0);
					if(!tetrahedron) tetrahedron = (MeshTetrahedron*)face->getBlock(1);
					assert(tetrahedron);
					assert(!tetrahedron->isInverted());
					MeshPoint3d *tpoints[4] = { tetrahedron->getPoint(0), tetrahedron->getPoint(1),
						tetrahedron->getPoint(2), tetrahedron->getPoint(3)};
					delete mesh->removeMeshTetrahedron(tetrahedron);

					MeshPoint3d *epoints[2] = { mpt1, mpt2 };
					for(int j = 0; j < 2; j++){
						mesh->addMeshTetrahedron(tetrahedron = new MeshTetrahedron(
							(tpoints[0] == epoints[j]) ? point : tpoints[0],
							(tpoints[1] == epoints[j]) ? point : tpoints[1],
							(tpoints[2] == epoints[j]) ? point : tpoints[2],
							(tpoints[3] == epoints[j]) ? point : tpoints[3]));
						assert(!tetrahedron->isInverted());
						new_tetrahedra.append(tetrahedron);
					}

					if(fct == 2) break;	// it was last tetrahedron to split
				}
				break; // the for loop
			}
		}else{
			// point can't be inserted at the border edge
		}
	}
	if(new_tetrahedra.countInt() == 0){
		// check for face-collision
		for(int i = 0; i < 4; i++){
			MeshFace* face = tetrahedron->getFace(i);
			int i1, i2, i3;
			tetrahedron->getFaceIndices(i, i1, i2, i3);
			bool on_face = DMTriangle3d::distance2ToPoint(
					point->getMetricCoordinates(mc), 
					points[i1]->getMetricCoordinates(mc), 
					points[i2]->getMetricCoordinates(mc),
					points[i3]->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER;
			if(!on_face){
				MeshPoint3d* fpoint = tetrahedron->getPoint(i); // opposite point
				on_face = GeometricPredicates::orient3d(
					((points[0] == fpoint) ? point : points[0])->getCoordinates(),
					((points[1] == fpoint) ? point : points[1])->getCoordinates(),
					((points[2] == fpoint) ? point : points[2])->getCoordinates(),
					((points[3] == fpoint) ? point : points[3])->getCoordinates()) < mesh_data.relative_small_number;
			}
			if(on_face){
				if(face->isBorder()) return false;
#ifdef LOG_DEL_SWAP
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Del-by-swap: inserting at face");
#endif
/*
				if(true){
					MeshViewSet *set = new MeshViewSet;
					set->addFace(face, 0);
					set->addPoint(point, 1);
					SHOW_MESH("face-collision", set);
				}
*/
				MeshTetrahedron* tetra[2] = { tetrahedron, 
					(MeshTetrahedron*)face->getOtherBlock(tetrahedron) };
/*
				if(true){
					MeshViewSet *set = new MeshViewSet;
					set->addFace(face, 0);
					set->addPoint(point, 1);
					set->addEmptyBlockWithEdges(tetra[0]);
					set->addEmptyBlockWithEdges(tetra[1]);
					SHOW_MESH("face-collision", set);
				}
*/
				MeshPoint3d* fpoints[3] = { face->getPoint(0),
					face->getPoint(1), face->getPoint(2) };
				for(int ti = 0; ti < 2; ti++){
					MeshPoint3d* tpoints[4] = { 
						tetra[ti]->getPoint(0), tetra[ti]->getPoint(1), 
						tetra[ti]->getPoint(2), tetra[ti]->getPoint(3) };
					LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
						"** tetra volume = " << tetra[ti]->getVolumeNoMetric());
					delete mesh->removeMeshTetrahedron(tetra[ti]);
					for(int fi = 0; fi < 3; fi++){
						mesh->addMeshTetrahedron(tetrahedron = new MeshTetrahedron(
							(tpoints[0]==fpoints[fi]) ? point : tpoints[0],
							(tpoints[1]==fpoints[fi]) ? point : tpoints[1],
							(tpoints[2]==fpoints[fi]) ? point : tpoints[2],
							(tpoints[3]==fpoints[fi]) ? point : tpoints[3]));
						assert(!tetrahedron->isInverted());
						LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
							"*** tetra volume = " << tetrahedron->getVolumeNoMetric());
						new_tetrahedra.append(tetrahedron);
					}
				}
				break; // the for loop
			}
		}
	}
	if(new_tetrahedra.countInt() == 0){
		// insert point within the tetrahedron
		delete mesh->removeMeshTetrahedron(tetrahedron);

		for(int j = 0; j < 4; j++){
			mesh->addMeshTetrahedron(tetrahedron = new MeshTetrahedron(
				(j == 0) ? point : points[0],
				(j == 1) ? point : points[1],
				(j == 2) ? point : points[2],
				(j == 3) ? point : points[3]));
			//if(tetrahedron->isInverted()){
			//	MeshViewSet* set = new MeshViewSet;
			//	set->addEdge(points[0]->getCoordinates(), points[3]->getCoordinates());
			//	set->addEdge(points[1]->getCoordinates(), points[3]->getCoordinates());
			//	set->addEdge(points[2]->getCoordinates(), points[3]->getCoordinates());
			//	set->addEdge(points[0]->getCoordinates(), points[2]->getCoordinates());
			//	set->addEdge(points[1]->getCoordinates(), points[2]->getCoordinates());
			//	set->addEdge(points[0]->getCoordinates(), points[1]->getCoordinates());
			//	set->addPoint(point);
			//	SHOW_MESH("inverted tetrahedron", set);
			//}
			assert(!tetrahedron->isInverted());
			new_tetrahedra.append(tetrahedron);
		}

#ifdef LOG_DEL_SWAP
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Del-by-swap: inserting within the block");
#endif
/*
		if(true){
			MeshViewSet *set = new MeshViewSet;
			for(int i = 0; i < 4; i++){
				set->addEmptyBlockWithEdges(new_tetrahedrons.get(i), 0);
				set->addPoint(tpoints[i], 0);
			}
			set->addPoint(point, 1);
			SHOW_MESH("tetrahedra-split", set);
		}
*/
#ifdef SHOW_DEL_SWAP
		MeshViewSet *set = new MeshViewSet;
		for(int i = 0; i < 4; i++)
			set->addBlockWithEdges(new_tetrahedrons.get(i), 0);
		set->addPoint(point);
		SHOW_MESH("tetrahedra-split", set);
#endif
	}
#ifdef LOG_DEL_SWAP
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Del-by-swap: new tetrahedra=" << new_tetrahedra.countInt());
#endif
	assert(new_tetrahedra.countInt() > 0);
	if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE){
		for (auto it = new_tetrahedra.iterator(); it.valid(); it.moveNext()) {
			it.item()->setAreaID(area_id1);
		}
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Del-by-swap: split, mesh valid " << mesh->isValid());

	// swap
	if(allow_swap) 
		iterativeTetrahedraSwap(mc, mesh, point, 5, 100, true, 
			mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE);

	return true;
}
#endif // USE_OPENMP_HERE && _OPENMP 

#if defined(USE_OPENMP_HERE) && defined(_OPENMP)

void MeshGenerator3d::iterativeTetrahedraSwap(Metric3dContext& mc, 
		MeshContainer3d* mesh, MeshPoint3d *init_point, 
		int max_layers, int max_swaps, bool local_metric, bool recalculate_quality)
{
	int mpct = mesh->getPointsCount();

	DataVector<MeshPoint3d*> active_points(mpct);
	assert(init_point);
	init_point->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); 
	active_points.add(init_point);
	MeshTetrahedron* tetrahedra[3];
	MeshPoint3d* points[5];

	int s32_checked = 0;
	int s32_swapped = 0;
	int s23_checked = 0;
	int s23_swapped = 0;

	int swap_count[3] = {-100, -100, -100};
	int layer = 0;

	Metric3dContext local_mc(mc);

	while(swap_count[0] != 0){
		swap_count[0] = 0;
		for(int i = 0; i < active_points.countInt(); i++){
			points[0] = active_points[i];
			for(int j = 0; j < points[0]->getRank(); j++){
				// check edge
				MeshEdge3d* edge = points[0]->getEdge(j);
				if(edge->isBorder()) continue;
				points[1] = edge->getOtherPoint(points[0]);
				bool swap_ok = false;
				if(edge->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP) && 
					edge->getFaceCount() == 3)
				{ // if tagged -> already checked
					++s32_checked;
					points[2] = edge->getFaceAt(0)->getOtherPoint(points[0], points[1]); 
					points[3] = edge->getFaceAt(1)->getOtherPoint(points[0], points[1]);
					points[4] = edge->getFaceAt(2)->getOtherPoint(points[0], points[1]);
					if(local_metric){
						const DPoint3d middle = DPoint3d::average( // metric set for 5 points (don't depend on swap)
							points[0]->getCoordinates(),
							points[1]->getCoordinates(), 
							points[2]->getCoordinates(), 
							points[3]->getCoordinates(), 
							points[4]->getCoordinates());
						local_mc.countMetricAtPoint(middle);
					}
					swap_ok = swap32(local_mc, mesh, edge, tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
				}
				if(swap_ok){
					++swap_count[0];
					++s32_swapped;
					if(s32_swapped + s23_swapped > max_swaps){
						#pragma omp critical (mesh)
						{
							mesh->clearAllTags(TagExtended::TAG_MG3D_SM_SWAP);
						}
						return;
					}
					j = 0; 
					for(int m = 0; m < 2; m++){
						tetrahedra[m]->clearTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
						if(recalculate_quality){
							tetrahedra[m]->countQuality(mc);
							if(mesh->isHeapOrder()) 
								mesh->updateBlockPosition(tetrahedra[m]);
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
						if(local_metric){
							const DPoint3d middle = DPoint3d::average(
								points[0]->getCoordinates(),
								points[1]->getCoordinates(), 
								points[2]->getCoordinates(),
								points[3]->getCoordinates(), 
								points[4]->getCoordinates());
							mc.countMetricAtPoint(middle);
						}
						swap_ok = swap23(mc, mesh, face, tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
					}
					if(swap_ok){ // start again from first edge
						++swap_count[0];
						++s23_swapped;
						if(s32_swapped + s23_swapped > max_swaps){
							mesh->clearAllTags(TagExtended::TAG_MG3D_SM_SWAP);
							return;
						}
						j = 0; 
						for(int m = 0; m < 3; m++){
							tetrahedra[m]->clearTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
							if(recalculate_quality){
								tetrahedra[m]->countQuality(mc);
								if(mesh->isHeapOrder()) 
									mesh->updateBlockPosition(tetrahedra[m]);
							}
						}
						break; 
					}else face->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
				}
			}
		}
		if(swap_count[0] == 0){ // expand set of active points (by neighbors)
			if(++layer > max_layers){
				mesh->clearAllTags(TagExtended::TAG_MG3D_SM_SWAP);
				return;
			}
			int act = active_points.countInt();
			int new_points_count = 0;
			for(int i = 0; i < act; i++){
				MeshPoint3d* point = active_points[i];
				for(int j = 0; j < point->getRank(); j++){
					MeshPoint3d* other_point = point->getEdge(j)->getOtherPoint(point);
					if(other_point->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP) && 
						(other_point->availableTag(TagExtended::TAG_BOUNDARY_POINT)))
					{
						other_point->setIntTag(TagExtended::TAG_MG3D_SM_SWAP);
						active_points.add(other_point);
						swap_count[0] = -layer;
						++new_points_count;
					}
				}
			}
			if(new_points_count == 0){
				mesh->clearAllTags(TagExtended::TAG_MG3D_SM_SWAP);
				return;
			}
			assert(swap_count[0] != 0);
		}
		if(swap_count[2] == swap_count[1] && swap_count[2] == swap_count[0]){
			LOG4CPLUS_INFO(MeshLog::logger_console, " ==== [restore] loop check exit =====");
			break;
		}
		swap_count[2] = swap_count[1];
		swap_count[1] = swap_count[0];
	}

	mesh->clearAllTags(TagExtended::TAG_MG3D_SM_SWAP);
}
#else
void MeshGenerator3d::iterativeTetrahedraSwap(Metric3dContext& mc, 
		MeshContainer3d* mesh, MeshPoint3d *init_point, 
		int max_layers, int max_swaps, bool local_metric, bool recalculate_quality)
{
	int mpct = mesh->getPointsCount();

	DataVector<MeshPoint3d*> active_points(mpct);
	assert(init_point);
	init_point->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); 
	active_points.add(init_point);
	MeshTetrahedron* tetrahedra[3];
	MeshPoint3d* points[5];

	int s32_checked = 0;
	int s32_swapped = 0;
	int s23_checked = 0;
	int s23_swapped = 0;

	int swap_count[3] = {-100, -100, -100};
	int layer = 0;

	while(swap_count[0] != 0){
		swap_count[0] = 0;
		for(size_t i = 0; i < active_points.countInt(); i++){
			points[0] = active_points[i];
			for(int j = 0; j < points[0]->getRank(); j++){
				// check edge
				MeshEdge3d* edge = points[0]->getEdge(j);
				if(edge->isBorder()) continue;
				points[1] = edge->getOtherPoint(points[0]);
				bool swap_ok = false;
				if(edge->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP) && 
					edge->getFaceCount() == 3)
				{ // if tagged -> already checked
					++s32_checked;
					points[2] = edge->getFaceAt(0)->getOtherPoint(points[0], points[1]); 
					points[3] = edge->getFaceAt(1)->getOtherPoint(points[0], points[1]);
					points[4] = edge->getFaceAt(2)->getOtherPoint(points[0], points[1]);
					if(local_metric){
						const DPoint3d middle = DPoint3d::average( // metric set for 5 points (don't depend on swap)
							points[0]->getCoordinates(),
							points[1]->getCoordinates(), 
							points[2]->getCoordinates(), 
							points[3]->getCoordinates(), 
							points[4]->getCoordinates());
						mc.countMetricAtPoint(middle);
					}
					swap_ok = swap32(mc, mesh, edge, tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
				}
				if(swap_ok){
					++swap_count[0];
					++s32_swapped;
					if(s32_swapped + s23_swapped > max_swaps){
						mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);
						return;
					}
					j = 0; 
					for(int m = 0; m < 2; m++){
						tetrahedra[m]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
						if(recalculate_quality){
							tetrahedra[m]->countQuality(mc);
							if(mesh->isHeapOrder()) 
								mesh->updateBlockPosition(tetrahedra[m]);
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
						if(local_metric){
							const DPoint3d middle = DPoint3d::average(
								points[0]->getCoordinates(),
								points[1]->getCoordinates(), 
								points[2]->getCoordinates(),
								points[3]->getCoordinates(), 
								points[4]->getCoordinates());
							mc.countMetricAtPoint(middle);
						}
						swap_ok = swap23(mc, mesh, face, tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
					}
					if(swap_ok){ // start again from first edge
						++swap_count[0];
						++s23_swapped;
						if(s32_swapped + s23_swapped > max_swaps){
							mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);
							return;
						}
						j = 0; 
						for(int m = 0; m < 3; m++){
							tetrahedra[m]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
							if(recalculate_quality){
								tetrahedra[m]->countQuality(mc);
								if(mesh->isHeapOrder()) 
									mesh->updateBlockPosition(tetrahedra[m]);
							}
						}
						break; 
					}else face->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
				}
			}
		}
		if(swap_count[0] == 0){ // expand set of active points (by neighbors)
			if(++layer > max_layers){
				mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);
				return;
			}
			size_t act = active_points.countInt();
			int new_points_count = 0;
			for(size_t i = 0; i < act; i++){
				MeshPoint3d* point = active_points[i];
				for(int j = 0; j < point->getRank(); j++){
					MeshPoint3d* other_point = point->getEdge(j)->getOtherPoint(point);
					if(other_point->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP) && 
						(other_point->availableTag(TagExtended::TAG_BOUNDARY_POINT)))
					{
						other_point->setIntTag(TagExtended::TAG_MG3D_SM_SWAP);
						active_points.add(other_point);
						swap_count[0] = -layer;
						++new_points_count;
					}
				}
			}
			if(new_points_count == 0){
				mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);
				return;
			}
			assert(swap_count[0] != 0);
		}
		if(swap_count[2] == swap_count[1] && swap_count[2] == swap_count[0]){
			LOG4CPLUS_DEBUG(MeshLog::logger_console, " ==== [restore] loop check exit =====");
			break;
		}
		swap_count[2] = swap_count[1];
		swap_count[1] = swap_count[0];
	}

	mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);
}
#endif

bool MeshGenerator3d::swap44(Metric3dContext& mc, MeshContainer3d* mesh, MeshFace *face1, 
		MeshEdge3d *edge, int /* phase */)
{
	assert(face1->isBorder() == false);
	assert(edge->getFaceCount() == 4);

//	assert(mesh->isValid());

	MeshTetrahedron* tetrahedrons[4] = {(MeshTetrahedron*)(face1->getBlock(0)), 
		nullptr, nullptr, (MeshTetrahedron*)(face1->getBlock(1))};
	MeshFace* faces[4] = {face1, nullptr, nullptr, nullptr};

	int i, fct = edge->getFaceCount();
	for(i = 0; i < fct; i++){
		MeshFace *face = edge->getFaceAt(i);
		if(face != face1){
			if(face->incidentToBlock(tetrahedrons[0])){
				faces[1] = face;
				tetrahedrons[1] = (MeshTetrahedron*)(face->getOtherBlock(tetrahedrons[0]));
			}else if(face->incidentToBlock(tetrahedrons[3])){
				faces[3] = face;
				tetrahedrons[2] = (MeshTetrahedron*)(face->getOtherBlock(tetrahedrons[3]));
			}else{
				faces[2] = face;
			}
		}
	}
	assert(tetrahedrons[1] != nullptr);
	assert(tetrahedrons[2] != nullptr);

#ifdef SHOW_DEL_SWAP
	if(phase > 1){
		MeshViewSet *set = new MeshViewSet;
		for(i = 0; i < 4; i++){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "block " << tetrahedrons[i]->getVolumeNoMetric());
			set->addBlockWithEdges(tetrahedrons[i]);
		}
		SHOW_MESH("swap44 (before)", set);
	}
#endif

	MeshPoint3d* edge_pt0 = edge->getMeshPoint(0);
	MeshPoint3d* edge_pt1 = edge->getMeshPoint(1);
	MeshPoint3d* points[4];

	for(i = 0; i < 4; i++){
		if(faces[i]->isBorder()) return false;
		points[i] = faces[i]->getOtherPoint(edge_pt0, edge_pt1);
	}

/*
	if(phase == 4){
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "=xx==> swap44: edge [" << edge_pt0->getID() 
			<< " == " << edge_pt1->getID() << "], other four: ["
			<< points[0]->getID() << ','
			<< points[1]->getID() << ','
			<< points[2]->getID() << ','
			<< points[3]->getID() << "]" << endl;
	}
*/

	// check (just in case...)
	const DMPoint3d dpt0 = points[0]->getMetricCoordinates(mc);
	const DMPoint3d dpt1 = points[1]->getMetricCoordinates(mc);
	const DMPoint3d dpt2 = points[2]->getMetricCoordinates(mc);
	const DMPoint3d dpt3 = points[3]->getMetricCoordinates(mc);
	const DMPoint3d dpte0 = edge_pt0->getMetricCoordinates(mc);
	const DMPoint3d dpte1 = edge_pt1->getMetricCoordinates(mc);

	double v0 = DMTriangle3d::orient3d(dpt0, dpt2, dpt1, dpte0);
	double v1 = DMTriangle3d::orient3d(dpt0, dpt1, dpt2, dpte1);
	double v2 = DMTriangle3d::orient3d(dpt0, dpt3, dpt2, dpte0);
	double v3 = DMTriangle3d::orient3d(dpt0, dpt2, dpt3, dpte1);

	if(abs(v0) < METRIC_SMALL_NUMBER || abs(v1) < METRIC_SMALL_NUMBER ||
		abs(v2) < METRIC_SMALL_NUMBER || abs(v3) < METRIC_SMALL_NUMBER)
	{
		LOG4CPLUS_WARN(MeshLog::logger_console, "Swap44 impossible -- would create too small tetrahedron.");
		return false;
	}

//	static int run_counter = 0;

/*
	if(phase == 4){
		LOG4CPLUS_INFO(MeshLog::logger_console, "swap44 before, mesh_valid", mesh->isValid());
		LOG4CPLUS_INFO(MeshLog::logger_console, "=== counter", ++run_counter);
		LOG4CPLUS_INFO(MeshLog::logger_console, "=== v0", v0);
		LOG4CPLUS_INFO(MeshLog::logger_console, "=== v1", v1);
		LOG4CPLUS_INFO(MeshLog::logger_console, "=== v2", v2);
		LOG4CPLUS_INFO(MeshLog::logger_console, "=== v3", v3);
	}
*/
	if(v0 > 0.0 && v1 > 0.0 && v2 > 0.0 && v3 > 0.0){
		// OK
	}else if(v0 < 0.0 && v1 < 0.0 && v2 < 0.0 && v3 < 0.0){
		// swap
		MeshPoint3d* mpt = edge_pt0; edge_pt0 = edge_pt1; edge_pt1 = mpt;
	}else{
		LOG4CPLUS_WARN(MeshLog::logger_console, "Swap44 impossible.");
		return false;
	}

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Starting swap44 - Mesh invalid, leaving."); return false;
		}
#endif

	int material_id_0 = tetrahedrons[0]->getAreaID();
	int material_id_1 = tetrahedrons[3]->getAreaID();
	for(i = 0; i < 4; i++){
		delete mesh->removeMeshTetrahedron(tetrahedrons[i]);
	}
	tetrahedrons[0] = new MeshTetrahedron(points[0], points[2], points[1], edge_pt0);
	tetrahedrons[1] = new MeshTetrahedron(points[0], points[1], points[2], edge_pt1);
	tetrahedrons[2] = new MeshTetrahedron(points[0], points[3], points[2], edge_pt0);
	tetrahedrons[3] = new MeshTetrahedron(points[0], points[2], points[3], edge_pt1);
	for(i = 0; i < 4; i++){
		tetrahedrons[i]->setAreaID((i<2)?material_id_0:material_id_1);
		mesh->addMeshTetrahedron(tetrahedrons[i]);
	}
	assert(edge_pt0->getEdgeToPoint(edge_pt1) == nullptr);
	assert(points[0]->getEdgeToPoint(points[2]) != nullptr);

#ifdef SHOW_DEL_SWAP
	if(phase > 1){
		MeshViewSet *set = new MeshViewSet;
		for(i = 0; i < 4; i++){
			double vol = tetrahedrons[i]->getVolumeNoMetric();
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "block " << vol);
			if(vol >= 0.0)
				set->addEmptyBlockWithEdges(tetrahedrons[i]);
			else set->addBlockWithEdges(tetrahedrons[i]);
		}
		SHOW_MESH("swap44 (after)", set);
	}
#endif

//	assert(mesh->isValid());
//	if(phase == 4) LOG4CPLUS_INFO(MeshLog::logger_console, "swap44 after, mesh_valid", mesh->isValid()); 

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Done swap44 - Mesh invalid."); 
		}
#endif

	return true;
}

MeshEdge3d* MeshGenerator3d::swap22(Metric3dContext& mc, MeshContainer3d* mesh, MeshEdge3d *edge, MeshTetrahedron** created_tetrahedrons)
{
	if(!edge->isBorder() || edge->isBorder(TagBorder::RIDGE | TagBorder::FIXED)) return nullptr;
	if(edge->getFaceCount() != 3) return nullptr; // 3 faces and 2 blocks is the only working scenario

//	assert(mesh->isValid());

	//for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
	//	MeshFace* face = it.getFace();
	//	if(face->isBoundedBothSides() && face->isBorder()){
	//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Face " << face->getPoint(0)->getIndex()
	//			<< "x" << face->getPoint(1)->getIndex()
	//			<< "x" << face->getPoint(2)->getIndex() << endl;
	//		SHOW_MESH("swap22-start, mis-labelled face", mesh->getDebugViewSet(face->getBlock(0)));
	//	}
	//}

	MeshFace* faces[3] = {nullptr, nullptr, nullptr};

	for(int i = 0; i < 3; i++){
		MeshFace* face = edge->getFaceAt(i);
		if(face->isBoundedBothSides()){
			assert(!face->isBorder());
			faces[2] = face;
		}else{
			assert(face->isBorder());
			if(faces[0]) faces[1] = face;
			else faces[0] = face;
		}
	}

	MeshTetrahedron* tetrahedrons[2] = {
		(MeshTetrahedron*)faces[0]->getOtherBlock(nullptr),
		(MeshTetrahedron*)faces[1]->getOtherBlock(nullptr)};

#ifdef SHOW_DEL_SWAP
	if(true){
		MeshViewSet *set = new MeshViewSet;
		for(int i = 0; i < 2; i++){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "block " << tetrahedrons[i]->getVolumeNoMetric());
			set->addBlockWithEdges(tetrahedrons[i]);
		}
		SHOW_MESH("swap22 (before)", set);
	}
#endif

	MeshPoint3d* epts[4] = {edge->getMeshPoint(0), edge->getMeshPoint(1), nullptr, nullptr };
	epts[2] = faces[0]->getOtherPoint(epts[0], epts[1]);
	epts[3] = faces[1]->getOtherPoint(epts[0], epts[1]);

	if(epts[2]->getEdgeToPoint(epts[3])){
		//SHOW_MESH("crossing edge", mesh->getDebugViewSet(epts[0], epts[1]));
		return nullptr;
	}

	// check quality condition
	MeshPoint3d* tp0[4] = {
		tetrahedrons[0]->getPoint(0), tetrahedrons[0]->getPoint(1),
		tetrahedrons[0]->getPoint(2), tetrahedrons[0]->getPoint(3)};
	MeshPoint3d* tp1[4] = {
		tetrahedrons[1]->getPoint(0), tetrahedrons[1]->getPoint(1),
		tetrahedrons[1]->getPoint(2), tetrahedrons[1]->getPoint(3)};

	DMPoint3d dmp0[4] = {
		tp0[0]->getMetricCoordinates(mc), tp0[1]->getMetricCoordinates(mc), 
		tp0[2]->getMetricCoordinates(mc), tp0[3]->getMetricCoordinates(mc)};
	DMPoint3d dmp1[4] = {
		tp1[0]->getMetricCoordinates(mc), tp1[1]->getMetricCoordinates(mc), 
		tp1[2]->getMetricCoordinates(mc), tp1[3]->getMetricCoordinates(mc)};

	double q0_pre = DTetrahedron::aspectRatio( dmp0[0], dmp0[1], dmp0[2], dmp0[3] );
	double q1_pre = DTetrahedron::aspectRatio( dmp1[0], dmp1[1], dmp1[2], dmp1[3] );

	double q0_after = DTetrahedron::aspectRatio( 
		tp0[0]==epts[0] ? epts[3]->getMetricCoordinates(mc) : dmp0[0], 
		tp0[1]==epts[0] ? epts[3]->getMetricCoordinates(mc) : dmp0[1], 
		tp0[2]==epts[0] ? epts[3]->getMetricCoordinates(mc) : dmp0[2], 
		tp0[3]==epts[0] ? epts[3]->getMetricCoordinates(mc) : dmp0[3] );
	double q1_after = DTetrahedron::aspectRatio( 
		tp1[0]==epts[1] ? epts[2]->getMetricCoordinates(mc) : dmp1[0], 
		tp1[1]==epts[1] ? epts[2]->getMetricCoordinates(mc) : dmp1[1], 
		tp1[2]==epts[1] ? epts[2]->getMetricCoordinates(mc) : dmp1[2], 
		tp1[3]==epts[1] ? epts[2]->getMetricCoordinates(mc) : dmp1[3] );

	if(std::min(q0_after, q1_after) <= std::min(q0_pre, q1_pre)) return nullptr; // current min quality is better than after swap

	char border_edge = edge->getBorderFlags();
	edge->clearBorder();
	char border_face0 = faces[0]->getBorderFlags();
	faces[0]->clearBorder();
	char border_face1 = faces[1]->getBorderFlags();
	faces[1]->clearBorder();
	TagExtended tag_face0, tag_face1;
	tag_face0.copyAllTags(faces[0]);
	tag_face1.copyAllTags(faces[1]);

	MeshPoint3d* ept4 = tetrahedrons[0]->getOppositePoint(faces[0]); assert(ept4);
	char border_edge_02 = epts[0]->getEdgeToPoint(epts[2])->getBorderFlags();
	char border_edge_12 = epts[1]->getEdgeToPoint(epts[2])->getBorderFlags();
	char border_edge_04 = epts[0]->getEdgeToPoint(ept4)->getBorderFlags();
	char border_edge_14 = epts[1]->getEdgeToPoint(ept4)->getBorderFlags();

	MeshPoint3d* tmp_point = new MeshPoint3d(epts[0]->getCoordinates());
	tetrahedrons[1]->switchPointsWithFaces(epts[1], tmp_point);
	tetrahedrons[0]->switchPointsWithFaces(epts[0], epts[3]);
	tetrahedrons[1]->switchPointsWithFaces(tmp_point, epts[2]);

	edge = epts[2]->getEdgeToPoint(epts[3]); assert(edge);
	edge->setBorder(border_edge);
	faces[0] = edge->getFaceToPoint(epts[1]); assert(faces[0]);
	faces[0]->setBorder(border_face0);
	faces[0]->copyAllTags(&tag_face0);
	faces[1] = edge->getFaceToPoint(epts[0]); assert(faces[1]);
	faces[1]->setBorder(border_face1);
	faces[1]->copyAllTags(&tag_face1);

	epts[0]->getEdgeToPoint(epts[2])->setBorder(border_edge_02);
	epts[1]->getEdgeToPoint(epts[2])->setBorder(border_edge_12);
	epts[0]->getEdgeToPoint(ept4)->setBorder(border_edge_04);
	epts[1]->getEdgeToPoint(ept4)->setBorder(border_edge_14);

	//for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
	//	MeshFace* face = it.getFace();
	//	if(!face->isBoundedBothSides()){
	//		if(!face->isBorder()){
	//			MeshViewSet* set = new MeshViewSet;
	//			set->addBlockWithEdges(tetrahedrons[0]);
	//			set->addBlockWithEdges(tetrahedrons[1]);
	//			set->addFace(face);
	//			for(int i = 0; i < 4; i++)
	//				set->addPoint(epts[i]);
	//			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Non-border face: " << face->getPoint(0)->getIndex()
	//				<< "x" << face->getPoint(1)->getIndex()
	//				<< "x" << face->getPoint(2)->getIndex() << endl;
	//			SHOW_MESH("Problematic case 22", set);
	//		}
	//		assert(face->isBorder());
	//	}
	//}
	assert(tmp_point->getRank() == 0);
	delete tmp_point;

	if(created_tetrahedrons){
		created_tetrahedrons[0] = tetrahedrons[0];
		created_tetrahedrons[1] = tetrahedrons[1];
	}

#ifdef SHOW_DEL_SWAP
	if(true){
		MeshViewSet *set = new MeshViewSet;
		for(i = 0; i < 2; i++){
			double vol = tetrahedrons[i]->getVolumeNoMetric();
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "block " << vol);
			if(vol >= 0.0)
				set->addEmptyBlockWithEdges(tetrahedrons[i]);
			else set->addBlockWithEdges(tetrahedrons[i]);
		}
		SHOW_MESH("swap22 (after)", set);
	}
#endif

//	assert(mesh->isValid());

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Done swap22 - Mesh invalid."); 
		}
#endif

	//for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
	//	MeshFace* face = it.getFace();
	//	if(face->isBoundedBothSides()){
	//		assert(!face->isBorder());
	//	}
	//}

	return edge;
}

bool MeshGenerator3d::swap23(Metric3dContext& mc, MeshContainer3d *mesh, MeshFace *face, 
							 MeshTetrahedron** created_tetrahedrons, int quality_condition)
{
	if(face->isBorder()) return false;

//	assert(mesh->isValid());

	// points of the common face
	MeshPoint3d* mpt0 = face->getPoint(0);
	MeshPoint3d* mpt1 = face->getPoint(1);
	MeshPoint3d* mpt2 = face->getPoint(2);

	MeshTetrahedron* tetrahedron0 = (MeshTetrahedron*)face->getBlock(0);
	MeshTetrahedron* tetrahedron1 = (MeshTetrahedron*)face->getBlock(1);
	assert(tetrahedron0 != nullptr);
	assert(tetrahedron1 != nullptr);
	assert(tetrahedron0->getType() == BLOCK_TETRA);
	assert(tetrahedron1->getType() == BLOCK_TETRA);
	if(tetrahedron0->getType() != BLOCK_TETRA) return false;
	if(tetrahedron1->getType() != BLOCK_TETRA) return false;

	// points of the recovered edge
	MeshPoint3d* mpt3 = tetrahedron0->getOppositePoint(face);
	MeshPoint3d* mpt4 = tetrahedron1->getOppositePoint(face);
	assert(mpt3 && mpt4);

	// check if swap is possible (e.g. whether it would create inverted tetra)
	const DMPoint3d dpt0 = mpt0->getMetricCoordinates(mc);
	const DMPoint3d dpt1 = mpt1->getMetricCoordinates(mc);
	const DMPoint3d dpt2 = mpt2->getMetricCoordinates(mc);
	const DMPoint3d dpt3 = mpt3->getMetricCoordinates(mc);
	const DMPoint3d dpt4 = mpt4->getMetricCoordinates(mc);

	// quality
	// (pt0, pt1, pt2, pt3)  -> (pt3, pt0, pt1, pt4)
	// (pt0, pt1, pt2, pt4)  -> (pt3, pt1, pt2, pt4)
	//						 -> (pt3, pt2, pt0, pt4)
	switch(quality_condition){
	case MeshData::SWAP3_ALWAYS:
		{
			if(DMTriangle3d::orient3d(dpt3, dpt4, dpt2, dpt1) < param_volume_threshold) return false;
			if(DMTriangle3d::orient3d(dpt3, dpt4, dpt1, dpt0) < param_volume_threshold) return false;
			if(DMTriangle3d::orient3d(dpt3, dpt4, dpt0, dpt2) < param_volume_threshold) return false;
			break;
		}
	case MeshData::SWAP3_MIN_VOLUME:
		{
			double vol_10 = DMTriangle3d::orient3d(dpt3, dpt4, dpt2, dpt1);
			if(vol_10 <= 0.0) return false;
			double vol_11 = DMTriangle3d::orient3d(dpt3, dpt4, dpt1, dpt0);
			if(vol_11 <= 0.0) return false;
			double vol_12 = DMTriangle3d::orient3d(dpt3, dpt4, dpt0, dpt2);
			if(vol_12 <= 0.0) return false;

			double vol_00 = DMTriangle3d::orient3d(dpt0, dpt1, dpt2, dpt3);
			double vol_01 = DMTriangle3d::orient3d(dpt0, dpt2, dpt1, dpt4);

			double min_volume_before = std::min(abs(vol_00), abs(vol_01));
			double min_volume_after = std::min(std::min(vol_10, vol_11), vol_12);
			if(min_volume_before > min_volume_after) // cancel swap
				return false;
			break;
		}
	case MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY:
		{
			double q_00 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt2, dpt3);
			double q_01 = DTetrahedron::aspectRatio(dpt0, dpt2, dpt1, dpt4);
			assert(q_00 >= 0.0);
			assert(q_01 >= 0.0);
			double q_10 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt2, dpt1);
			double q_11 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt1, dpt0);
			double q_12 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt0, dpt2);

			double min_quality_before = std::min(q_00, q_01);
			double min_quality_after = std::min(std::min(q_10, q_11), q_12);
			if(min_quality_before > min_quality_after) // cancel swap
				return false;
		}
		break;
	case MeshData::SWAP3_AVE_TETRAHEDRON_QUALITY:
		{
			double q_00 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt2, dpt3);
			double q_01 = DTetrahedron::aspectRatio(dpt0, dpt2, dpt1, dpt4);
			assert(q_00 > 0.0);
			assert(q_01 > 0.0);
			double q_10 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt2, dpt1);
			if(q_10 <= 0.0) return false;
			double q_11 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt1, dpt0);
			if(q_11 <= 0.0) return false;
			double q_12 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt0, dpt2);
			if(q_12 <= 0.0) return false;

			double ave_quality_before = 0.5*(q_00 + q_01);
			double ave_quality_after = (1.0/3.0) * (q_10 + q_11 + q_12);
			if(ave_quality_before > ave_quality_after) // cancel swap
				return false;
		}
		break;
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console, "MG3d::swap23 -> unknown quality condition");
		break;
	}

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Starting swap23 - Mesh invalid, leaving."); return false;
		}
#endif

	int id = tetrahedron0->getAreaID();

	#pragma omp critical (mesh)
	{
		// remove 2 tetrahedrons
		delete mesh->removeMeshTetrahedron(tetrahedron0);
		delete mesh->removeMeshTetrahedron(tetrahedron1);

		// create 3 tetrahedrons
		tetrahedron0 = new MeshTetrahedron(mpt0, mpt4, mpt3, mpt1);
		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE)
			tetrahedron0->setAreaID(id);
		mesh->addMeshTetrahedron(tetrahedron0);
		if(created_tetrahedrons) created_tetrahedrons[0] = tetrahedron0;

		tetrahedron0 = new MeshTetrahedron(mpt1, mpt4, mpt3, mpt2);
		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE)
			tetrahedron0->setAreaID(id);
		mesh->addMeshTetrahedron(tetrahedron0);
		if(created_tetrahedrons) created_tetrahedrons[1] = tetrahedron0;

		tetrahedron0 = new MeshTetrahedron(mpt2, mpt4, mpt3, mpt0);
		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE)
			tetrahedron0->setAreaID(id);
		mesh->addMeshTetrahedron(tetrahedron0);
		if(created_tetrahedrons) created_tetrahedrons[2] = tetrahedron0;
	}

	assert(mpt3->getEdgeToPoint(mpt4) != nullptr);

//	assert(mesh->isValid());

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Done swap23 - Mesh invalid."); 
		}
#endif

	return true;
}

bool MeshGenerator3d::swap23Advantageous(Metric3dContext& mc, MeshFace *face, int quality_condition)
{
	if(face->isBorder()) return false;

	// points of the common face
	MeshPoint3d* mpt0 = face->getPoint(0);
	MeshPoint3d* mpt1 = face->getPoint(1);
	MeshPoint3d* mpt2 = face->getPoint(2);

	MeshTetrahedron* tetrahedron0 = (MeshTetrahedron*)face->getBlock(0);
	MeshTetrahedron* tetrahedron1 = (MeshTetrahedron*)face->getBlock(1);
	assert(tetrahedron0 != nullptr);
	assert(tetrahedron1 != nullptr);
	assert(tetrahedron0->getType() == BLOCK_TETRA);
	assert(tetrahedron1->getType() == BLOCK_TETRA);

	// points of the recovered edge
	MeshPoint3d* mpt3 = tetrahedron0->getOppositePoint(face);
	MeshPoint3d* mpt4 = tetrahedron1->getOppositePoint(face);
	assert(mpt3 && mpt4);

	const DMPoint3d dpt0 = mpt0->getMetricCoordinates(mc);
	const DMPoint3d dpt1 = mpt1->getMetricCoordinates(mc);
	const DMPoint3d dpt2 = mpt2->getMetricCoordinates(mc);
	const DMPoint3d dpt3 = mpt3->getMetricCoordinates(mc);
	const DMPoint3d dpt4 = mpt4->getMetricCoordinates(mc);

	switch(quality_condition){
	case MeshData::SWAP3_ALWAYS:
		{
			if(DMTriangle3d::orient3d(dpt3, dpt4, dpt2, dpt1) < param_volume_threshold) return false;
			if(DMTriangle3d::orient3d(dpt3, dpt4, dpt1, dpt0) < param_volume_threshold) return false;
			return (DMTriangle3d::orient3d(dpt3, dpt4, dpt0, dpt2) >= param_volume_threshold);
		}
	case MeshData::SWAP3_MIN_VOLUME:
		{
			double vol_10 = DTetrahedron::volume(dpt3, dpt4, dpt2, dpt1);
			if(vol_10 <= 0.0) return false;
			double vol_11 = DTetrahedron::volume(dpt3, dpt4, dpt1, dpt0);
			if(vol_11 <= 0.0) return false;
			double vol_12 = DTetrahedron::volume(dpt3, dpt4, dpt0, dpt2);
			if(vol_12 <= 0.0) return false;

			double vol_00 = DMTriangle3d::orient3d(dpt0, dpt1, dpt2, dpt3);
			double vol_01 = DMTriangle3d::orient3d(dpt0, dpt2, dpt1, dpt4);

			double min_volume_before = std::min(abs(vol_00), abs(vol_01));
			double min_volume_after = std::min(std::min(vol_10, vol_11), vol_12);
			return(min_volume_after > min_volume_before);
		}
	case MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY:
		{
			double q_00 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt2, dpt3);
			double q_01 = DTetrahedron::aspectRatio(dpt0, dpt2, dpt1, dpt4);
			assert(q_00 > 0.0);
			assert(q_01 > 0.0);
			double q_10 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt2, dpt1);
			double q_11 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt1, dpt0);
			double q_12 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt0, dpt2);

			double min_quality_before = std::min(q_00, q_01);
			double min_quality_after = std::min(std::min(q_10, q_11), q_12);
			return (min_quality_after > min_quality_before);
		}
	case MeshData::SWAP3_AVE_TETRAHEDRON_QUALITY:
		{
			double q_00 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt2, dpt3);
			double q_01 = DTetrahedron::aspectRatio(dpt0, dpt2, dpt1, dpt4);
			assert(q_00 > 0.0);
			assert(q_01 > 0.0);
			double q_10 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt2, dpt1);
			if(q_10 <= 0.0) return false;
			double q_11 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt1, dpt0);
			if(q_11 <= 0.0) return false;
			double q_12 = DTetrahedron::aspectRatio(dpt3, dpt4, dpt0, dpt2);
			if(q_12 <= 0.0) return false;

			double ave_quality_before = 0.5*(q_00 + q_01);
			double ave_quality_after = (q_10 + q_11 + q_12) / 3.0;
			return (ave_quality_after > ave_quality_before);
		}
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MG3d::swap23Advantageous -> unknown quality condition");
		break;
	}

	return false;
}

bool MeshGenerator3d::swap23possible(Metric3dContext& mc, MeshFace *face)
{
	if(face->isBorder()) return false;

	// points of the common face
	MeshPoint3d* mpt0 = face->getPoint(0);
	MeshPoint3d* mpt1 = face->getPoint(1);
	MeshPoint3d* mpt2 = face->getPoint(2);

	MeshTetrahedron* tetrahedron0 = (MeshTetrahedron*)face->getBlock(0);
	MeshTetrahedron* tetrahedron1 = (MeshTetrahedron*)face->getBlock(1);
	assert(tetrahedron0 != nullptr);
	assert(tetrahedron1 != nullptr);
	assert(tetrahedron0->getType() == BLOCK_TETRA);
	assert(tetrahedron1->getType() == BLOCK_TETRA);

	// points of the recovered edge
	MeshPoint3d* mpt3 = tetrahedron0->getOppositePoint(face);
	MeshPoint3d* mpt4 = tetrahedron1->getOppositePoint(face);
	assert(mpt3 && mpt4);

	// check if swap is possible (e.g. whether it would create inverted tetra)
	// check if swap is possible (e.g. whether it would create inverted tetra)
	const DMPoint3d dpt0 = mpt0->getMetricCoordinates(mc);
	const DMPoint3d dpt1 = mpt1->getMetricCoordinates(mc);
	const DMPoint3d dpt2 = mpt2->getMetricCoordinates(mc);
	const DMPoint3d dpt3 = mpt3->getMetricCoordinates(mc);
	const DMPoint3d dpt4 = mpt4->getMetricCoordinates(mc);

	return DMTriangle3d::orient3d(dpt3, dpt4, dpt1, dpt2) >= param_volume_threshold &&
		DMTriangle3d::orient3d(dpt3, dpt4, dpt2, dpt0) >= param_volume_threshold &&
		DMTriangle3d::orient3d(dpt3, dpt4, dpt0, dpt1) >= param_volume_threshold;
}

bool MeshGenerator3d::swap32(Metric3dContext& mc, MeshContainer3d *mesh, MeshEdge3d *edge, MeshTetrahedron** created_tetrahedrons, int quality_condition)
{
	if(edge->isBorder()) return false;

//	assert(mesh->isValid());

	if(edge->getFaceCount() != 3) 
		return false;

	MeshFace* face = edge->getFaceAt(0);
	assert(!face->isBorder());

	// points of this edge
	MeshPoint3d* mpt0 = edge->getMeshPoint(0);
	MeshPoint3d* mpt1 = edge->getMeshPoint(1);
	if(!face->properOrientation(mpt0, mpt1)){
		mpt0 = edge->getMeshPoint(1);
		mpt1 = edge->getMeshPoint(0);
	}

	// tetrahedrons adjacent to this edge
	MeshPoint3d* mpt2 = face->getOtherPoint(mpt0, mpt1);
	MeshTetrahedron* tetrahedron0 = (MeshTetrahedron*)face->getBlock(0);

	face = tetrahedron0->getIncidentFace(face, edge);
	assert(!face->isBorder());
	MeshPoint3d* mpt3 = face->getOtherPoint(mpt0, mpt1);
	MeshTetrahedron* tetrahedron1 = (MeshTetrahedron*)face->getOtherBlock(tetrahedron0);

	face = tetrahedron1->getIncidentFace(face, edge);
	assert(!face->isBorder());
	MeshPoint3d* mpt4 = face->getOtherPoint(mpt0, mpt1);
	MeshTetrahedron* tetrahedron2 = (MeshTetrahedron*)face->getOtherBlock(tetrahedron1);

	assert(tetrahedron0->getType() == BLOCK_TETRA);
	assert(tetrahedron1->getType() == BLOCK_TETRA);
	assert(tetrahedron2->getType() == BLOCK_TETRA);
	if(tetrahedron0->getType() != BLOCK_TETRA) return false;
	if(tetrahedron1->getType() != BLOCK_TETRA) return false;
	if(tetrahedron2->getType() != BLOCK_TETRA) return false;

	// check if swap is possible (e.g. whether it would create inverted tetra)
	const DMPoint3d dpt0 = mpt0->getMetricCoordinates(mc);
	const DMPoint3d dpt1 = mpt1->getMetricCoordinates(mc);
	const DMPoint3d dpt2 = mpt2->getMetricCoordinates(mc);
	const DMPoint3d dpt3 = mpt3->getMetricCoordinates(mc);
	const DMPoint3d dpt4 = mpt4->getMetricCoordinates(mc);

	switch(quality_condition){
	case MeshData::SWAP3_ALWAYS:
		{
			if(DMTriangle3d::orient3d(dpt2, dpt3, dpt4, dpt1) < param_volume_threshold) return false;
			if(DMTriangle3d::orient3d(dpt2, dpt4, dpt3, dpt0) < param_volume_threshold) return false;
			break;
		}
	case MeshData::SWAP3_MIN_VOLUME:
		{
			double vol_10 = DMTriangle3d::orient3d(dpt2, dpt3, dpt4, dpt1);
			if(vol_10 <= 0.0) return false;
			double vol_11 = DMTriangle3d::orient3d(dpt2, dpt4, dpt3, dpt0);
			if(vol_11 <= 0.0) return false;

			double vol[3] = {
				DMTriangle3d::orient3d(dpt0, dpt1, dpt2, dpt3),
				DMTriangle3d::orient3d(dpt0, dpt1, dpt3, dpt4),
				DMTriangle3d::orient3d(dpt0, dpt1, dpt4, dpt2)};

			double min_volume_before = std::min(std::min(vol[0], vol[1]), vol[2]);
			double min_volume_after = std::min(vol_10, vol_11);
			if(min_volume_before > min_volume_after) // cancel swap
				return false;
		}
		break;
	case MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY:
		{
			double q_00 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt2, dpt3);
			double q_01 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt3, dpt4);
			double q_02 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt4, dpt2);
			assert(q_00 >= 0.0);
			assert(q_01 >= 0.0);
			assert(q_02 >= 0.0);
			double q_10 = DTetrahedron::aspectRatio(dpt2, dpt3, dpt4, dpt1);
			double q_11 = DTetrahedron::aspectRatio(dpt2, dpt4, dpt3, dpt0);

			double min_quality_before = std::min(std::min(q_00, q_01), q_02);
			double min_quality_after = std::min(q_10, q_11);
			if(min_quality_before > min_quality_after) // cancel swap
				return false;
		}
		break;
	case MeshData::SWAP3_AVE_TETRAHEDRON_QUALITY:
		{
			double q_00 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt2, dpt3);
			double q_01 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt3, dpt4);
			double q_02 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt4, dpt2);
			assert(q_00 >= 0.0);
			assert(q_01 >= 0.0);
			assert(q_02 >= 0.0);
			double q_10 = DTetrahedron::aspectRatio(dpt2, dpt3, dpt4, dpt1);
			if(q_10 <= 0.0) return false;
			double q_11 = DTetrahedron::aspectRatio(dpt2, dpt4, dpt3, dpt0);
			if(q_11 <= 0.0) return false;

			double ave_quality_before = (1.0/3.0) * (q_00 + q_01 + q_02);
			double ave_quality_after = 0.5 * (q_10 + q_11);
			if(ave_quality_before > ave_quality_after) // cancel swap
				return false;
		}
		break;
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MG3d::swap32 -> unknown quality condition");
		break;
	}

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Starting swap32 - Mesh invalid, leaving."); return false;
		}
#endif

	int id = tetrahedron0->getAreaID();

	#pragma omp critical (mesh)
	{
		// remove 3 tetrahedrons
		delete mesh->removeMeshTetrahedron(tetrahedron0);
		delete mesh->removeMeshTetrahedron(tetrahedron1);
		delete mesh->removeMeshTetrahedron(tetrahedron2);

		// create 2 tetrahedrons
		tetrahedron0 = new MeshTetrahedron(mpt2, mpt4, mpt3, mpt0);
		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE)
			tetrahedron0->setAreaID(id);
		mesh->addMeshTetrahedron(tetrahedron0);

		tetrahedron1 = new MeshTetrahedron(mpt2, mpt3, mpt4, mpt1);
		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE)
			tetrahedron1->setAreaID(id);
		mesh->addMeshTetrahedron(tetrahedron1);
	}

	if(created_tetrahedrons){
		created_tetrahedrons[0] = tetrahedron0;
		created_tetrahedrons[1] = tetrahedron1;
	}

	assert(mpt2->getFaceToPoints(mpt3, mpt4) != nullptr);

//	assert(mesh->isValid());

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Done swap32 - Mesh invalid."); 
		}
#endif

	return true;
}

bool MeshGenerator3d::swap32possible(Metric3dContext& mc, MeshEdge3d *edge)
{
	if(edge->isBorder() || edge->getFaceCount() != 3) return false;

	MeshFace* face = edge->getFaceAt(0);
	assert(!face->isBorder());

	// points of this edge
	MeshPoint3d* mpt0 = edge->getMeshPoint(0);
	MeshPoint3d* mpt1 = edge->getMeshPoint(1);
	if(!face->properOrientation(mpt0, mpt1)){
		mpt0 = edge->getMeshPoint(1);
		mpt1 = edge->getMeshPoint(0);
	}

	// tetrahedrons adjacent to this edge
	MeshPoint3d* mpt2 = face->getOtherPoint(mpt0, mpt1);
	MeshTetrahedron* tetrahedron0 = (MeshTetrahedron*)face->getBlock(0);

	face = tetrahedron0->getIncidentFace(face, edge);
	assert(!face->isBorder());
	MeshPoint3d* mpt3 = face->getOtherPoint(mpt0, mpt1);
	MeshTetrahedron* tetrahedron1 = (MeshTetrahedron*)face->getOtherBlock(tetrahedron0);

	face = tetrahedron1->getIncidentFace(face, edge);
	assert(!face->isBorder());
	MeshPoint3d* mpt4 = face->getOtherPoint(mpt0, mpt1);
	MeshTetrahedron* tetrahedron2 = (MeshTetrahedron*)face->getOtherBlock(tetrahedron1);

	assert(tetrahedron0->getType() == BLOCK_TETRA);
	assert(tetrahedron1->getType() == BLOCK_TETRA);
	assert(tetrahedron2->getType() == BLOCK_TETRA);
	if(tetrahedron0->getType() != BLOCK_TETRA) return false;
	if(tetrahedron1->getType() != BLOCK_TETRA) return false;
	if(tetrahedron2->getType() != BLOCK_TETRA) return false;

	// check if swap is possible (e.g. whether it would create inverted tetra)
	const DMPoint3d dpt0 = mpt0->getMetricCoordinates(mc);
	const DMPoint3d dpt1 = mpt1->getMetricCoordinates(mc);
	const DMPoint3d dpt2 = mpt2->getMetricCoordinates(mc);
	const DMPoint3d dpt3 = mpt3->getMetricCoordinates(mc);
	const DMPoint3d dpt4 = mpt4->getMetricCoordinates(mc);

	return (DMTriangle3d::orient3d(dpt2, dpt3, dpt4, dpt1) >= param_volume_threshold) &&
			(DMTriangle3d::orient3d(dpt2, dpt4, dpt3, dpt0) >= param_volume_threshold);
}

bool MeshGenerator3d::swap32Advantageous(Metric3dContext& mc, MeshEdge3d *edge, int quality_condition)
{
	if(edge->isBorder()) return false;

//	assert(mesh->isValid());

	if(edge->getFaceCount() != 3) return false;

	MeshFace* face = edge->getFaceAt(0);
	assert(!face->isBorder());

	// points of this edge
	MeshPoint3d* mpt0 = edge->getMeshPoint(0);
	MeshPoint3d* mpt1 = edge->getMeshPoint(1);
	if(!face->properOrientation(mpt0, mpt1)){
		mpt0 = edge->getMeshPoint(1);
		mpt1 = edge->getMeshPoint(0);
	}

	// tetrahedrons adjacent to this edge
	MeshPoint3d* mpt2 = face->getOtherPoint(mpt0, mpt1);
	MeshTetrahedron* tetrahedron0 = (MeshTetrahedron*)face->getBlock(0);

	face = tetrahedron0->getIncidentFace(face, edge);
	assert(!face->isBorder());
	MeshPoint3d* mpt3 = face->getOtherPoint(mpt0, mpt1);
	MeshTetrahedron* tetrahedron1 = (MeshTetrahedron*)face->getOtherBlock(tetrahedron0);

	face = tetrahedron1->getIncidentFace(face, edge);
	assert(!face->isBorder());
	MeshPoint3d* mpt4 = face->getOtherPoint(mpt0, mpt1);
	MeshTetrahedron* tetrahedron2 = (MeshTetrahedron*)face->getOtherBlock(tetrahedron1);

	assert(tetrahedron0->getType() == BLOCK_TETRA);
	assert(tetrahedron1->getType() == BLOCK_TETRA);
	assert(tetrahedron2->getType() == BLOCK_TETRA);
	if(tetrahedron0->getType() != BLOCK_TETRA) return false;
	if(tetrahedron1->getType() != BLOCK_TETRA) return false;
	if(tetrahedron2->getType() != BLOCK_TETRA) return false;

	// check if swap is possible (e.g. whether it would create inverted tetra)
	const DMPoint3d dpt0 = mpt0->getMetricCoordinates(mc);
	const DMPoint3d dpt1 = mpt1->getMetricCoordinates(mc);
	const DMPoint3d dpt2 = mpt2->getMetricCoordinates(mc);
	const DMPoint3d dpt3 = mpt3->getMetricCoordinates(mc);
	const DMPoint3d dpt4 = mpt4->getMetricCoordinates(mc);

	switch(quality_condition){
	case MeshData::SWAP3_ALWAYS:
		{
			return (
				DMTriangle3d::orient3d(dpt2, dpt3, dpt4, dpt1) >= param_volume_threshold &&
				DMTriangle3d::orient3d(dpt2, dpt4, dpt3, dpt0) >= param_volume_threshold);
		}
	case MeshData::SWAP3_MIN_VOLUME:
		{
			double vol_10 = DMTriangle3d::orient3d(dpt2, dpt3, dpt4, dpt1);
			if(vol_10 <= 0.0) return false;
			double vol_11 = DMTriangle3d::orient3d(dpt2, dpt4, dpt3, dpt0);
			if(vol_11 <= 0.0) return false;

			double vol[3] = {
				DMTriangle3d::orient3d(dpt0, dpt1, dpt2, dpt3),
				DMTriangle3d::orient3d(dpt0, dpt1, dpt3, dpt4),
				DMTriangle3d::orient3d(dpt0, dpt1, dpt4, dpt2)};

			double min_volume_before = std::min(std::min(vol[0], vol[1]), vol[2]);
			double min_volume_after = std::min(vol_10, vol_11);
			return (min_volume_after > min_volume_before);
		}
	case MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY:
		{
			double q_00 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt2, dpt3);
			double q_01 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt3, dpt4);
			double q_02 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt4, dpt2);
			assert(q_00 >= 0.0);
			assert(q_01 >= 0.0);
			assert(q_02 >= 0.0);
			double q_10 = DTetrahedron::aspectRatio(dpt2, dpt3, dpt4, dpt1);
			double q_11 = DTetrahedron::aspectRatio(dpt2, dpt4, dpt3, dpt0);

			double min_quality_before = std::min(std::min(q_00, q_01), q_02);
			double min_quality_after = std::min(q_10, q_11);
			return (min_quality_after > min_quality_before);
		}
	case MeshData::SWAP3_AVE_TETRAHEDRON_QUALITY:
		{
			double q_00 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt2, dpt3);
			double q_01 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt3, dpt4);
			double q_02 = DTetrahedron::aspectRatio(dpt0, dpt1, dpt4, dpt2);
			assert(q_00 >= 0.0);
			assert(q_01 >= 0.0);
			assert(q_02 >= 0.0);
			double q_10 = DTetrahedron::aspectRatio(dpt2, dpt3, dpt4, dpt1);
			if(q_10 <= 0.0) return false;
			double q_11 = DTetrahedron::aspectRatio(dpt2, dpt4, dpt3, dpt0);
			if(q_11 <= 0.0) return false;

			double ave_quality_before = (1.0/3.0) * (q_00 + q_01 + q_02);
			double ave_quality_after = 0.5 * (q_10 + q_11);
			return (ave_quality_after > ave_quality_before);
		}
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MG3d::swap32 -> unknown quality condition");
		break;
	}

	return false;
}

#if defined(USE_OPENMP_HERE) && defined(_OPENMP)

int MeshGenerator3d::addInnerNodes(Metric3dContext& mc, MeshContainer3d *mesh)
{
	if(!mesh) return 0;

	if(!MeshGenerator2d::param_triangulate_with_inner_nodes) return 0;

//	if(!mesh || mesh->getInnerEdgesCount() > 0) return 0;
	mesh->clearSearchTree();

	START_CLOCK("MG3d::addInnerNodes");
	
	int border_count = mesh->getPointsCount();

	// Calculate quality of the available tetrahedra
	int count = mesh->getBlocksCount();
	if(count < 1) return 0;

	#pragma omp parallel for shared(mesh,count) firstprivate(mc)
	for(int i = 0; i < count; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		block->countQuality(mc);
	}

	//START_CLOCK("MG3d::initialGuess");
	//double ev_nt_count = 0.0;

	//#pragma omp parallel for shared(mesh,count) firstprivate(mc) reduction(+:ev_nt_count)
	//for(int i = 0; i < count; i++){
	//	ev_nt_count += mesh->getBlockAt(i)->getVolume(mc, false);
	//}

	//ev_nt_count *= 12.0/SQRT2;
	//STOP_CLOCK("MG3d::initialGuess");
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Initial guess NT3d: " << ev_nt_count);

	mesh->setHeapOrder(true);
	int new_ct = 0;
	int steps = 0;
	const DBox bounding_box = mesh->getBoundingBox();

	int count_pt_outside_box = 0;
	int count_pt_no_tetra_found = 0;
	int count_pt_too_near_other_pt = 0;
	int count_pt_in_other_tetra = 0;
	int count_pt_too_diff_metrics = 0;
	int count_pt_too_near_other_pt_2 = 0;

//	clock_t clock_start = clock();
//	int insert_circum_count = 0;
//	int insert_edge_count = 0;

//	int last_nt_ev_step = -1;

	Metric3dContext** pmc;

#pragma omp parallel shared(mesh, pmc, bounding_box) firstprivate(mc) \
					reduction(+:count_pt_outside_box) reduction(+:count_pt_no_tetra_found) \
					reduction(+:count_pt_too_near_other_pt) reduction(+:count_pt_in_other_tetra) \
					reduction(+:count_pt_too_diff_metrics) reduction(+:count_pt_too_near_other_pt_2) \
					reduction(+:steps) reduction(+:new_ct)
	{
		int thr_id = omp_get_thread_num();
		int thr_ct = omp_get_num_threads();

		#pragma omp single
		pmc = new Metric3dContext*[thr_ct];

		pmc[thr_id] = &mc; // pointer to thread-local variable

		#pragma omp barrier

		while(true){
			// Worst tetrahedron (with account for minimum volume) is located in the root
			MeshTetrahedron *tetrahedron = nullptr;
			MeshTetrahedron* main_tetrahedron = nullptr;
			int active_ispheres = 0;

			#pragma omp critical (mesh)
			{
				int pct = mesh->getPointsCount();
				for(int i = 0; i < pct; i++){
					tetrahedron = (MeshTetrahedron*)mesh->getBlockAt(i);
					double tq = tetrahedron->getQuality();
					if(tq > param_quality_threshold) break; // actually it's only heap, not sorted, but let's better break
					// quality is low enough
					// ... check collision with "impact sphere" of other threads
					DPoint3d middle = tetrahedron->getMiddlePoint();
					mc.countMetricAtPoint(middle);
					for(int j = 0; j < thr_ct; j++){
						if(j == thr_id) continue; // don't check with itself
						if(pmc[j]->withinImpactRadius(middle)){
							// collision!
							tetrahedron = nullptr;
							break;
						}
					}
					if(tetrahedron){ // accepted
						mc.setRadius(IMPACT_SPHERE_FACTOR * tetrahedron->getOuterSphereRadius(mc, false));
						main_tetrahedron = tetrahedron;
						break;
					}
				}

				if(!tetrahedron){ // if not found, check whether to finish
					for(int j = 0; j < thr_ct; j++)
						if(! pmc[j]->noRadius()) ++active_ispheres;
				}
			}

			if(!tetrahedron)
				if(active_ispheres) continue; // try again
				else break;  // finish

			//int current_nt_ev_step = (int)(20 * tq);
			//if(current_nt_ev_step > last_nt_ev_step && tq >= 0.01){
			//	int nt3d = mesh->getBlocksCount();
			//	if(true){
			//		double nt_expected = nt3d * (param_quality_threshold) / tq;
			//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "NT3d= " << setw(8) << nt3d << ", qt= " << fixed << setprecision(2) << tq 
			//			<< ", predicted= " << setw(8) << (int)nt_expected << endl;
			//	}

			//	last_nt_ev_step = current_nt_ev_step;
			//}

			++steps;

			MeshPoint3d* points[] = {
				tetrahedron->getPoint(0), 
				tetrahedron->getPoint(1),
				tetrahedron->getPoint(2), 
				tetrahedron->getPoint(3)};
			DPoint3d dnew;
			if(tq < param_quality_improvement){
				//case MeshData::IMPROVE_LONGEST_EDGE:
				// Find longest edge
				int imax = -1;
				double d2, dmax = 0.0;
				for(int i = 0; i < 6; i++){
					if(tetrahedron->getEdge(i)->isBorder()) continue;
					int i1, i2;
					tetrahedron->getEdgeIndices(i, i1, i2);
					if((d2 = points[i1]->getMetricCoordinates(mc).distance2(
						points[i2]->getMetricCoordinates(mc))) > dmax){
						dmax = d2;
						imax = i;
					}
				}
				// If edge is too short -> cancel
				if(dmax < 4.0){
					tetrahedron->setQuality(mesh_data.relative_infinity);

					#pragma omp critical (mesh)
					{
						mesh->updateBlockPosition(tetrahedron);
						mc.setNoRadius();
					}
					continue;
				}
				dnew = tetrahedron->getEdge(imax)->getPoint(0.5);
	
				#pragma omp critical (mesh)
				{
					mc.countMetricAtPoint(dnew);
					for(int j = 0; j < thr_ct; j++){
						if(j == thr_id) continue; // don't check with itself
						if(mc.withinImpactSphere(dnew)){
							// collision!
							mc.setNoRadius();
							break;
						}
					}
					if(! mc.noRadius()) // no collision, so update
						mc.setRadius(IMPACT_SPHERE_FACTOR * tetrahedron->getOuterSphereRadius(mc, false));
				}

				if(mc.noRadius())
					continue; // select another tetrahedron to improve

	//			++insert_edge_count;
			}else{
				//case MeshData::IMPROVE_OUTER_CIRCLE:
				const DMPoint3d dnew_metric = tetrahedron->getOuterSphereCenter(mc, false);
				dnew = mc.transformMStoRS(dnew_metric);
				if(!bounding_box.contains(dnew)){
					main_tetrahedron->setQuality(mesh_data.relative_infinity);
					#pragma omp critical (mesh)
					{
						mesh->updateBlockPosition(main_tetrahedron);
						mc.setNoRadius();
					}
					++count_pt_outside_box;
					continue;
				}
				if(!tetrahedron->isPointInside(dnew)){
					#pragma omp critical (mesh)
					{
						tetrahedron = tetrahedron->findTetrahedronByNeighbours(dnew, true);
					}
					if(!tetrahedron){
						main_tetrahedron->setQuality(mesh_data.relative_infinity);
						#pragma omp critical (mesh)
						{
							mesh->updateBlockPosition(main_tetrahedron);
							mc.setNoRadius();
						}
						++count_pt_no_tetra_found;
						continue;
					}
				}
				int i;
				for(i = 0; i < 4; i++){
					if(dnew_metric.distance2(
						tetrahedron->getPoint(i)->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER){
	//						SHOW_STEP_PT(2, "Nowy punkt zbyt blisko innego.", dnew);
						main_tetrahedron->setQuality(mesh_data.relative_infinity);
						#pragma omp critical (mesh)
						{
							mesh->updateBlockPosition(main_tetrahedron);
							mc.setNoRadius();
						}
						++count_pt_too_near_other_pt;
						break;
					}
				}
				if(i < 4) continue; // jeœli powy¿sza pêtla zosta³a przerwana

				#pragma omp critical (mesh)
				{
					mc.countMetricAtPoint(dnew);
					for(int j = 0; j < thr_ct; j++){
						if(j == thr_id) continue; // don't check with itself
						if(pmc[j]->withinImpactRadius(dnew)){
							// collision!
							mc.setNoRadius();
							break;
						}
					}
					if(! mc.noRadius()) // no collision, so update
						mc.setRadius(IMPACT_SPHERE_FACTOR * tetrahedron->getOuterSphereRadius(mc, false));
				}
				if(mc.noRadius())
					continue; // select another tetrahedron to improve

				if(main_tetrahedron != tetrahedron){
					const DMPoint3d dnew_second_metric = mc.transformRStoMS(dnew);
					assert(tetrahedron->isPointInOuterSphere(mc, dnew_second_metric, false)); // since dnew is inside the tetrahedron !!!
					if(!main_tetrahedron->isPointInOuterSphere(mc, dnew_second_metric, false)){
	//						SHOW_STEP_PT(2, "Zbyt ró¿ne metryki dla wstawienia punktu.", dnew);
						main_tetrahedron->setQuality(mesh_data.relative_infinity);
						#pragma omp critical (mesh)
						{
							mesh->updateBlockPosition(main_tetrahedron);
							mc.setNoRadius();
						}
						++count_pt_too_diff_metrics;
						continue;
					}
					for(i = 0; i < 4; i++){
						if(dnew_second_metric.distance2(
							tetrahedron->getPoint(i)->getMetricCoordinates(mc)) < 0.25)
						{
	//							SHOW_STEP_PT(2, "Nowy punkt zbyt blisko innego.", dnew);
							main_tetrahedron->setQuality(mesh_data.relative_infinity);
							#pragma omp critical (mesh)
							{
								mesh->updateBlockPosition(main_tetrahedron);
								mc.setNoRadius();
							}
							++count_pt_too_near_other_pt_2;
							break;
						}
					}
					if(i < 4) continue; // jeœli powy¿sza pêtla zosta³a przerwana
					++count_pt_in_other_tetra;
				}
	//			++insert_circum_count;
			}
			MeshPoint3d* p3 = new MeshPoint3d(dnew);
			if(!addPointToTriangulation(mc, mesh, p3, tetrahedron, main_tetrahedron, true)){
				delete p3;
				main_tetrahedron->setQuality(mesh_data.relative_infinity);
				#pragma omp critical (mesh)
				{
					mesh->updateBlockPosition(main_tetrahedron);
					mc.setNoRadius();
				}
				continue;
			}

			++new_ct;
	//		if(new_ct % 1000 == 0) LOG4CPLUS_INFO(MeshLog::logger_mesh, "Adding 3d inner nodes - " << new_ct);
	//		count = border_count + new_ct;
	//		SHOW_STATUS(count);

			#pragma omp critical (mesh)
			{
				mc.setNoRadius();
			}
		}
	}

	delete[] pmc;

//	double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
//	LOG4CPLUS_INFO(MeshLog::logger_console, "inner-time", sec);
//	LOG4CPLUS_INFO(MeshLog::logger_console, "inner-circum", ((100.0 * insert_circum_count) / (insert_circum_count + insert_edge_count)));

	LOG4CPLUS_INFO(MeshLog::logger_mesh, new_ct << " inner nodes inserted.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_in_other_tetra << " pts in other tetra.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Cancelled insertions: ";
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, (count_pt_outside_box + count_pt_no_tetra_found + count_pt_too_near_other_pt +
		count_pt_too_diff_metrics + count_pt_too_near_other_pt_2) << endl;
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_outside_box << " pts outside box.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_no_tetra_found << " pts with no containing tetra found.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_too_near_other_pt << " pts too near other pt.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_too_diff_metrics << " pts inserted in too different metrics.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_too_near_other_pt_2 << " pts too near other_pt (if other tetra).");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "final block count = " << mesh->getBlocksCount());

//	SHOW_STATUS_END(count);

	STOP_CLOCK("MG3d::addInnerNodes");

//#ifdef STAT_COUNT
//	MeshRepository::addMark("MeshGenerator2d::addInnerNodes - end.");
//#endif
		
	mesh->setHeapOrder(false);

	DataVector<MeshContainer3d::AdditionalBoundaryNode>& inserted_boundary_points = mesh->getAdditionalBoundaryNodes();
	if(inserted_boundary_points.countInt() > 0){
		MeshGenerator3dDelaunayBoundary::fixAdditionalBoundaryNodes(mc, mesh, inserted_boundary_points);
		LOG4CPLUS_INFO(MeshLog::logger_console, "Additional boundary nodes (final-inner)", inserted_boundary_points.countInt());
	}

	return new_ct;
}

#else

int MeshGenerator3d::addInnerNodes(Metric3dContext& mc, MeshContainer3d *mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return 0;

	if(!MeshGenerator2d::param_triangulate_with_inner_nodes) return 0;

//	if(!mesh || mesh->getInnerEdgesCount() > 0) return 0;
	mesh->clearSearchTree();

	START_CLOCK("MG3d::addInnerNodes");
	
	int border_count = mesh->getPointsCount();

	// Wyznacz jakoœæ dostêpnych trójk¹tów
	int count = mesh->getBlocksCount();
	if(count < 1) return 0;

	for(int i = 0; i < count; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		if(tag_type == TagExtended::TAG_NONE || block->checkIntTagForAnyPoint(tag_type, tag_value)){
			block->countQuality(mc);
		}else{
			block->setQuality(mesh_data.relative_infinity);
		}
	}

	//START_CLOCK("MG3d::initialGuess");
	//double ev_nt_count = 0.0;

	//#pragma omp parallel for shared(mesh,count) firstprivate(mc) reduction(+:ev_nt_count)
	//for(int i = 0; i < count; i++){
	//	ev_nt_count += mesh->getBlockAt(i)->getVolume(mc, false);
	//}

	//ev_nt_count *= 12.0/SQRT2;
	//STOP_CLOCK("MG3d::initialGuess");
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Initial guess NT3d: " << ev_nt_count);

	mesh->setHeapOrder(true);
	int new_ct = 0;
	int steps = 0;
	const DBox bounding_box = mesh->getBoundingBox();

	int count_pt_outside_box = 0;
	int count_pt_no_tetra_found = 0;
	int count_pt_too_near_other_pt = 0;
	int count_pt_in_other_tetra = 0;
	int count_pt_too_diff_metrics = 0;
	int count_pt_too_near_other_pt_2 = 0;

//	clock_t clock_start = clock();
//	int insert_circum_count = 0;
//	int insert_edge_count = 0;

//	int last_nt_ev_step = -1;

	while(true){
		// Worst tetrahedron (with account for minimum volume) is located in the root
		MeshTetrahedron *tetrahedron = (MeshTetrahedron*)mesh->getBlockAt(0);
		MeshTetrahedron* main_tetrahedron = tetrahedron;

//		char text[100];
//		sprintf(text, "Najgorszy trójk¹t (jakoœæ = %0.2f)", triangle->getQuality());
//		SHOW_STEP_PT_BREAKABLE(2, text, triangle->getMiddlePoint(), true);
		
		// Is the quality low enough ?
		double tq = tetrahedron->getQuality();
		if(tq > param_quality_threshold) break;

		//const int MAX_NT_COUNT = 0;
		const int MAX_NT_COUNT = 10000000;
		if (MAX_NT_COUNT > 0 && mesh->getBlocksCount() >= MAX_NT_COUNT) {
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "*STOP* Reached max number of NT  : " << fixed << MAX_NT_COUNT);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " ** current (skip) tetra quality : " << fixed << tq);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " ** threshold for tetra quality  : " << fixed << param_quality_threshold);
			break;
		}

		//int current_nt_ev_step = (int)(20 * tq);
		//if(current_nt_ev_step > last_nt_ev_step && tq >= 0.01){
		//	int nt3d = mesh->getBlocksCount();
		//	if(true){
		//		double nt_expected = nt3d * (param_quality_threshold) / tq;
		//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "NT3d= " << setw(8) << nt3d << ", qt= " << fixed << setprecision(2) << tq 
		//			<< ", predicted= " << setw(8) << (int)nt_expected << endl;
		//	}

		//	last_nt_ev_step = current_nt_ev_step;
		//}

		++steps;

		//const int MAX_STEP_COUNT = 0;
		const int MAX_STEP_COUNT = 10000000;
		if (MAX_STEP_COUNT > 0 && steps >= MAX_STEP_COUNT) {
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "*STOP* Reached max number of steps  : " << fixed << MAX_STEP_COUNT);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " ** current (skip) tetra quality : " << fixed << tq);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " ** threshold for tetra quality  : " << fixed << param_quality_threshold);
			break;
		}


		if (steps % 100000 == 0) {
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
				"inserting inner nodes: step=" << steps << "\t NP=" << mesh->getPointsCount());
		}

		mc.countMetricAtPoint(tetrahedron->getMiddlePoint());

		// IMPACT-SPHERE-RADIUS
		//double Rcavity = main_tetrahedron->getOuterSphereRadius(mc, false);
		//double qq = ((int)(tq * 100 + 0.5)) * 0.01;
		//LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl << qq << "\t" << Rcavity << "\t";

		MeshPoint3d* points[] = {
			tetrahedron->getPoint(0), 
			tetrahedron->getPoint(1),
			tetrahedron->getPoint(2), 
			tetrahedron->getPoint(3)};
		DPoint3d dnew;
		if(tq < param_quality_improvement){
			//case MeshData::IMPROVE_LONGEST_EDGE:
			// Find longest edge
			int imax = -1;
			double d2, dmax = 0.0;
			for(int i = 0; i < 6; i++){
				if(tetrahedron->getEdge(i)->isBorder()) continue;
				int i1, i2;
				tetrahedron->getEdgeIndices(i, i1, i2);
				if((d2 = points[i1]->getMetricCoordinates(mc).distance2(
					points[i2]->getMetricCoordinates(mc))) > dmax){
					dmax = d2;
					imax = i;
				}
			}
			// If edge is too short -> cancel
			if(dmax < 4.0){
				tetrahedron->setQuality(mesh_data.relative_infinity);
				mesh->updateBlockPosition(tetrahedron);
				continue;
			}
			dnew = tetrahedron->getEdge(imax)->getPoint(0.5);
			mc.countMetricAtPoint(dnew);
//			++insert_edge_count;
		}else{
			//case MeshData::IMPROVE_OUTER_CIRCLE:
			const DMPoint3d dnew_metric = tetrahedron->getOuterSphereCenter(mc, false);
			dnew = mc.transformMStoRS(dnew_metric);
			if(!bounding_box.contains(dnew)){
				main_tetrahedron->setQuality(mesh_data.relative_infinity);
				mesh->updateBlockPosition(main_tetrahedron);
				++count_pt_outside_box;
				continue;
			}
			if(!tetrahedron->isPointInside(dnew)){
				tetrahedron = tetrahedron->findTetrahedronByNeighbours(dnew, true);
				if(!tetrahedron){
					main_tetrahedron->setQuality(mesh_data.relative_infinity);
					mesh->updateBlockPosition(main_tetrahedron);
					++count_pt_no_tetra_found;
					continue;
				}
			}
			int i;
			for(i = 0; i < 4; i++){
				if(dnew_metric.distance2(
					tetrahedron->getPoint(i)->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER){
//						SHOW_STEP_PT(2, "Nowy punkt zbyt blisko innego.", dnew);
					main_tetrahedron->setQuality(mesh_data.relative_infinity);
					mesh->updateBlockPosition(main_tetrahedron);
					++count_pt_too_near_other_pt;
					break;
				}
			}
			if(i < 4) continue; // jeœli powy¿sza pêtla zosta³a przerwana

			mc.countMetricAtPoint(dnew);
			if(main_tetrahedron != tetrahedron){
				const DMPoint3d dnew_second_metric = mc.transformRStoMS(dnew);
				assert(tetrahedron->isPointInOuterSphere(mc, dnew_second_metric, false)); // since dnew is inside the tetrahedron !!!
				if(!main_tetrahedron->isPointInOuterSphere(mc, dnew_second_metric, false)){
//						SHOW_STEP_PT(2, "Zbyt ró¿ne metryki dla wstawienia punktu.", dnew);
					main_tetrahedron->setQuality(mesh_data.relative_infinity);
					mesh->updateBlockPosition(main_tetrahedron);
					++count_pt_too_diff_metrics;
					continue;
				}
				for(i = 0; i < 4; i++){
					if(dnew_second_metric.distance2(
						tetrahedron->getPoint(i)->getMetricCoordinates(mc)) < 0.25)
					{
//							SHOW_STEP_PT(2, "Nowy punkt zbyt blisko innego.", dnew);
						main_tetrahedron->setQuality(mesh_data.relative_infinity);
						mesh->updateBlockPosition(main_tetrahedron);
						++count_pt_too_near_other_pt_2;
						break;
					}
				}
				if(i < 4) continue; // jeœli powy¿sza pêtla zosta³a przerwana
				++count_pt_in_other_tetra;
			}
//			++insert_circum_count;
		}
		MeshPoint3d* p3 = new MeshPoint3d(dnew);
		if(!addPointToTriangulation(mc, mesh, p3, tetrahedron, main_tetrahedron, true)){
			delete p3;
			main_tetrahedron->setQuality(mesh_data.relative_infinity);
			mesh->updateBlockPosition(main_tetrahedron);
			continue;
		}

		if(tag_type != TagExtended::TAG_NONE)
			p3->setIntTag(tag_type, tag_value);

		++new_ct;
//		if(new_ct % 1000 == 0) LOG4CPLUS_INFO(MeshLog::logger_mesh, "Adding 3d inner nodes - " << new_ct);

		//static const int MAX_INNER_PCT = 1000000;
		//if (new_ct >= MAX_INNER_PCT) {
		//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "inserting inner nodes - reached max " << new_ct << "new nodes: STOP");
		//	break;
		//}

		count = border_count + new_ct;
//		SHOW_STATUS(count);
	}

//	double sec = (clock() - clock_start)/(double)CLOCKS_PER_SEC;
//	LOG4CPLUS_INFO(MeshLog::logger_console, "inner-time", sec);
//	LOG4CPLUS_INFO(MeshLog::logger_console, "inner-circum", ((100.0 * insert_circum_count) / (insert_circum_count + insert_edge_count)));

	LOG4CPLUS_INFO(MeshLog::logger_mesh, new_ct << " inner nodes inserted.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_in_other_tetra << " pts in other tetra.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"Cancelled insertions: " 
		<< (count_pt_outside_box + count_pt_no_tetra_found + count_pt_too_near_other_pt +
			count_pt_too_diff_metrics + count_pt_too_near_other_pt_2));
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_outside_box << " pts outside box.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_no_tetra_found << " pts with no containing tetra found.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_too_near_other_pt << " pts too near other pt.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_too_diff_metrics << " pts inserted in too different metrics.");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "  - " << count_pt_too_near_other_pt_2 << " pts too near other_pt (if other tetra).");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "final block count = " << mesh->getBlocksCount());

//	SHOW_STATUS_END(count);

	STOP_CLOCK("MG3d::addInnerNodes");

//#ifdef STAT_COUNT
//	MeshRepository::addMark("MeshGenerator2d::addInnerNodes - end.");
//#endif
		
	mesh->setHeapOrder(false);

	DataVector<MeshContainer3d::AdditionalBoundaryNode>& inserted_boundary_points = mesh->getAdditionalBoundaryNodes();
	if(inserted_boundary_points.countInt() > 0){
		MeshGenerator3dDelaunayBoundary::fixAdditionalBoundaryNodes(mc, mesh, inserted_boundary_points);
		LOG4CPLUS_INFO(MeshLog::logger_console, 
			"Additional boundary nodes (final-inner): " << inserted_boundary_points.countInt());
	}

	return new_ct;
}

#endif // USE_OPENMP_HERE && _OPENMP

bool MeshGenerator3d::autoTriangulate(MeshContainer3d* boundary, int sm_count)
{
	boundary->clearDiscretization();
	int pcount = boundary->getPointsCount();
	int bct = boundary->getBlocksCount();
	//if(pcount < 1) return false;
	int loop_counter = 0;
	int invalid_count = 0;

	for(int i = 0; i < bct; i++){
		MeshDomainVolume* mdv = ((MeshDomainVolume*)boundary->getBlockAt(i));
		mdv->clearControlSpace();
		mdv->createInitialControlSpace();
	}

	do{
		++loop_counter;
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Auto-control-adjustment [3D]: loop #" << loop_counter);
		// -> GEN-1D + GEN-2D = surfaces
		START_CLOCK("MG3d::AutoGen-2D");
		MeshGenerator2d::autoTriangulate(boundary, sm_count, false, true);
		STOP_CLOCK("MG3d::AutoGen-2D");
		invalid_count = 0;
		// -> GEN-3D = volumes
		START_CLOCK("MG3d::AutoGen-3D-Boundary");
		for(int i = 0; i < bct; i++){
			MeshDomainVolume* domain_volume = (MeshDomainVolume*)boundary->getBlockAt(i);
			if(domain_volume->getMesh() == nullptr || 
			   domain_volume->getMesh()->getDiscretizationState() == 0)
			{
				domain_volume->prepareBoundaryMesh();
				CS3dPtr control_space = domain_volume->getControlSpace();
				Metric3dContext mc(control_space);
				if(control_space->isAdaptive() && loop_counter < param_max_auto_retriangulation_count){
					//LOG4CPLUS_INFO(MeshLog::logger_console, "CS3d - pre-boundary");
					//((CS3dPtr)control_space)->logDescription();
					if(!domain_volume->checkControlAtBoundary(mc)) {
						++invalid_count;
						LOG4CPLUS_DEBUG(MeshLog::logger_console, "Re-discretization of boundary surface required.");
					}else{
						domain_volume->discretizeUsingBoundary(mc);
						if(!domain_volume->checkControlForCloseBoundaryFaces(mc)){
							++invalid_count;
							LOG4CPLUS_DEBUG(MeshLog::logger_console, "Re-discretization of boundary surface required.");
						}
					}
					//LOG4CPLUS_INFO(MeshLog::logger_console, "CS3d - post-boundary");
					//((CS3dPtr)control_space)->logDescription();
				}else 
					domain_volume->discretizeUsingBoundary(mc);
			}
		}
		STOP_CLOCK("MG3d::AutoGen-3D-Boundary");
	}while(invalid_count > 0);

	START_CLOCK("MG3d::AutoGen-3D-Inner-total");
	for(int i = 0; i < bct; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		MeshContainer3d* mesh = volume->getMesh();
		LOG_ASSERT(mesh != nullptr);
		// MeshingException("Empty mesh"), false);
		CS3dPtr cs = mesh->getControlSpace();
		if(cs->isAdaptive())
			cs->getAsAdaptive()->compact();
		Metric3dContext mc(cs);
/*
		if(MeshGenerator2d::param_mesh_decomposition == 11 || 
				MeshGenerator2d::param_mesh_decomposition == 12){ 
			// wd decomposition
			MeshGenerator3d::createDecompositionMesh(mesh);
			MeshContainer3d* second_mesh = mesh->splitByBlocks(2);

			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NT", mesh->getBlocksCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NT", second_mesh->getBlocksCount());

//			SHOW_MESH("decomposition mesh [0]", mesh->getViewSet());
//			SHOW_MESH("decomposition mesh [1]", second_mesh->getViewSet());
			// continue meshing
			START_CLOCK("DCMP3d - mesh [0] final refinement");
			MeshGenerator3d::addInnerNodes(mc, mesh);
			MeshGenerator3dQuality::smoothen(mc, mesh, sm_count);
			STOP_CLOCK("DCMP3d - mesh [0] final refinement");

			START_CLOCK("DCMP3d - mesh [1] final refinement");
			MeshGenerator3d::addInnerNodes(mc, second_mesh);
			MeshGenerator3dQuality::smoothen(mc, second_mesh, sm_count);
			STOP_CLOCK("DCMP3d - mesh [1] final refinement");

			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NT", mesh->getBlocksCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NT", second_mesh->getBlocksCount());

			mesh->logQuality();
			second_mesh->logQuality();

//			SHOW_MESH("final mesh [0]", mesh->getViewSet());
//			SHOW_MESH("final mesh [1]", second_mesh->getViewSet());

//			MeshViewSet* vset = mesh->getViewSet();
//			vset = second_mesh->getViewSet(vset);
//			SHOW_MESH("decomposition mesh [0+1]", vset);
			delete second_mesh;
		}else{
*/
			MeshGenerator3d::addInnerNodes(mc, mesh);
			MeshGenerator3dQuality::smoothen(mc, mesh, sm_count, TagExtended::TAG_NONE, 1, true,
				MeshData::SM3_OPT_MIXED | MeshData::SM3_SWAP_COMPLETE | MeshData::SM3_BORDER_PRUNE);
//		}
		mesh->setDiscretizationState(2);
	}
	STOP_CLOCK("MG3d::AutoGen-3D-Inner-total");

	return true;
}

bool MeshGenerator3d::autoTriangulateRemesh(MeshContainer3d* boundary, int sm_count)
{
	return false;
}

bool MeshGenerator3d::createDecompositionMesh(MeshContainer3d* mesh)
{
	// Calculate quality of present triangles + inertia center + some axis

//	SHOW_MESH("coarse mesh", mesh->getViewSet());

	int count = mesh->getBlocksCount();
	if(count < 1) return 0;

	CS3dPtr main_space = mesh->getControlSpace();
	assert(main_space);
	Metric3dContext mc(main_space);

	double total_mass = 0;
	DPoint3d inertial_center;

	START_CLOCK("DCMP3d - NT3 evaulation");

	for(int i = 0; i < count; i++)
		mesh->getBlockAt(i)->addForInertialCenter(mc, inertial_center, total_mass);
	inertial_center /= total_mass;

	STOP_CLOCK("DCMP3d - NT3 evaulation");

	int ev_count = (int) (total_mass * (12.0/SQRT2));
	LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP evaluation NT3 count: " << ev_count);

	DMatrix3d inertial_moments;
	for(int i = 0; i < count; i++)
		mesh->getBlockAt(i)->addForInertialMoments(mc, inertial_center, inertial_moments);

	LOG4CPLUS_DEBUG(MeshLog::logger_console, 
		"** inertial center 3d = " << inertial_center);

	DMatrix3d v; // eigenvectors
	double d[3]; // eigenvalues;
	inertial_moments.eigensystem(v, d);

	for(int i = 0; i < 3; i++){
		LOG4CPLUS_DEBUG(MeshLog::logger_console, 
			"** vector-" << i << " = " << v.column(i).normalized());
		LOG4CPLUS_DEBUG(MeshLog::logger_console, 
			"**  value-" << i << " = " << d[i]);
	}

/*
	double ev_nt_count = 0.0;
	for(int i = 0; i < count; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		block->countQuality(mc);
		double vol = block->getVolume(mc, false);
		ev_nt_count += vol;
		inertia_center += block->getMiddlePoint() * vol;
	}

	// normalize inertia center (divide by total mass)
	inertia_center /= ev_nt_count;
*/

	// calculate center and axis of inertia
	const DBox box = mesh->getBoundingBox();
	// ... for now simple plane OXY, crossing the inertia center
	const DPoint3d dcmp_pt00(box.x0, box.y0, inertial_center.z);
	const DPoint3d dcmp_pt01(box.x1, box.y0, inertial_center.z);
	const DPoint3d dcmp_pt10(box.x0, box.y1, inertial_center.z);
	const DPoint3d dcmp_pt11(box.x1, box.y1, inertial_center.z);
	// ... for now simple plane OYZ, crossing the inertia center
//	const DPoint3d dcmp_pt00(inertial_center.x, box.y0, box.z0);
//	const DPoint3d dcmp_pt01(inertial_center.x, box.y1, box.z0);
//	const DPoint3d dcmp_pt10(inertial_center.x, box.y0, box.z1);
//	const DPoint3d dcmp_pt11(inertial_center.x, box.y1, box.z1);
		
	auto decomposition_space = std::make_shared<ControlSpace3dOctree>(box);

	// create special control space
	START_CLOCK("DCMP3d - special CS creation");
	
	if(MeshGenerator2d::param_mesh_decomposition == 11){
		// "1 - special control space as created copy"

		// set max
		decomposition_space->getAsAdaptive()->setMaxMetric(1.0);

		const DVector3d dvn = DVector3d::crossProduct(dcmp_pt00, dcmp_pt01, dcmp_pt10); // normal vector

		// copy metric along border lines
		const double r = MeshGenerator2d::param_mesh_decomposition_width;
		decomposition_space->getAsAdaptive()->setMinControlLine(mc, main_space, dcmp_pt00, dcmp_pt01, dvn, r);
		decomposition_space->getAsAdaptive()->setMinControlLine(mc, main_space, dcmp_pt00, dcmp_pt10, dvn, r);
		decomposition_space->getAsAdaptive()->setMinControlLine(mc, main_space, dcmp_pt11, dcmp_pt01, dvn, r);
		decomposition_space->getAsAdaptive()->setMinControlLine(mc, main_space, dcmp_pt11, dcmp_pt10, dvn, r);

		// copy metric within rectangular area (recursive quadtree)
		decomposition_space->getAsAdaptive()->setMinControlRect(mc, main_space,
			dcmp_pt00, dcmp_pt01 - dcmp_pt00, dcmp_pt10 - dcmp_pt00, dvn, r);

		// smoothen
		if(MeshGenerator2d::param_mesh_decomposition_gradation > 0)
			decomposition_space->getAsAdaptive()->smoothenWithGradation(MeshGenerator2d::param_mesh_decomposition_gradation);

	}else if(MeshGenerator2d::param_mesh_decomposition == 12){
		// "12 - special control space as filtered layer"
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Not implemented yet.");
	}

	// switch spaces
	mesh->setControlSpace(decomposition_space);

	STOP_CLOCK("DCMP3d - special CS creation");

	START_CLOCK("DCMP3d - special mesh refinement");
	Metric3dContext dcmp_mc(decomposition_space);
	MeshGenerator3d::addInnerNodes(dcmp_mc, mesh);

	MeshGenerator3dQuality::smoothenLaplaceMixed(dcmp_mc, mesh);
	MeshGenerator3dQuality::smoothenSwap(dcmp_mc, mesh);
	MeshGenerator3dQuality::smoothenLaplaceMixed(dcmp_mc, mesh);

	// mark elements
	count = mesh->getBlocksCount();
	for(int i = 0; i < count; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		int det_neg = 0;
		int det_pos = 0;
		for(int j = 0; j < block->getPointCount(); j++){
			double det = DTriangle3d::orient3d(dcmp_pt00, dcmp_pt01, dcmp_pt10, block->getPoint(j)->getCoordinates());
			if(det >= 0) ++det_pos;
			else ++det_neg;
		}
//		if(det_neg == 0) element->setAreaID(2);
//		else if(det_pos == 0) element->setAreaID(3);
//		else element->setAreaID(4);
		block->setAreaID((det_pos > det_neg) ? 2 : 6);
	}

	// reduce cut-size
	for(int i = 0; i < count; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		if(block->getPointCount() != 4) continue;
		int area_id = block->getAreaID();
		int sim_count[2] = {0,0};
		for(int j = 0; j < 4; j++){
			MeshBlock* other_block = block->getNeighbour(j);
			if(!other_block) continue;
			sim_count[(other_block->getAreaID() == area_id) ? 0 : 1]++;
		}
		if(sim_count[1] > sim_count[0]) 
			block->setAreaID((2+6) - area_id); // switch
	}

	STOP_CLOCK("DCMP3d - special mesh refinement");

	// count interface faces
	if(true){
		int interface_count = 0;
		for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
			MeshFace* face = it.getFace();
			MeshBlock* block0 = face->getBlock(0);
			MeshBlock* block1 = face->getBlock(1);
			if(block0 && block1 && (block0->getAreaID() != block1->getAreaID()))
				++interface_count;
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_console, 
			"Interface NF (after reduce): count=" << interface_count);
	}

	LOG4CPLUS_INFO(MeshLog::logger_console, 
		"Decomposition mesh: NT=" << count);
//	SHOW_MESH("decomposition mesh", mesh->getViewSet());

	// and switch space back
	// ... decomposition_space should be removed automatically
	mesh->setControlSpace(main_space); 

	return true;
}

MeshPoint3d * MeshGenerator3d::splitEdgeSimple(MeshContainer3d * mesh, MeshEdge3d * edge, const DPoint3d & pt)
{
	assert(!edge->isBorder()); // for now, TODO later for boundary-edge case
	assert(!edge->hasLocalCurve());

	MeshPoint3d* p0 = edge->getMeshPoint(0);
	MeshPoint3d* p1 = edge->getMeshPoint(1);

	//assert(!p0->hasLocalSurface());
	//assert(!p1->hasLocalSurface());

	MeshPoint3d* point = new MeshPoint3d(pt);
	//point->setBaseNormal(((p0->getBaseNormal() + p1->getBaseNormal()) * 0.5).normalized());

	mesh->addMeshPoint(point);
	// set border info
	//char bflags = edge->getBorderFlags();
	//point->setBorder(bflags); // ie. ridge or normal boundary

	// rearange the mesh (i.e. "split" blocks)
	DataVector<MeshBlock*> blocks;
	if (!edge->adjacentBlocks(blocks)) return nullptr;

	blocks.forEach([&](MeshBlock* block) {
		MeshBlock* block_clone = block->clone();
		block->switchPointsWithFaces(p1, point);
		block_clone->switchPointsWithFaces(p0, point);
		mesh->addMeshBlock(block_clone);
	});

	//if (bflags != 0) {
	//	delete edge;
	//	edge = p0->getEdgeToPoint(point);
	//	assert(edge);
	//	edge->setBorderFlags(bflags);
	//	edge = p1->getEdgeToPoint(point);
	//	assert(edge);
	//	edge->setBorderFlags(bflags);
	//}

	if (false) {
		SHOW_MESH("After split", mesh->getDebugViewSet(point));
	}

	return point;
}

MeshPoint3d * MeshGenerator3d::splitFaceSimple(MeshContainer3d * mesh, MeshFace * face, const DPoint3d & pt)
{
	assert(!face->isBorder()); // for now, TODO later for boundary-face case
	assert(!face->hasLocalSurface());

	MeshPoint3d* p0 = face->getPoint(0);
	MeshPoint3d* p1 = face->getPoint(1);
	MeshPoint3d* p2 = face->getPoint(2);

	//assert(!p0->hasLocalSurface());
	//assert(!p1->hasLocalSurface());

	MeshPoint3d* point = new MeshPoint3d(pt);
	//point->setBaseNormal(((p0->getBaseNormal() + p1->getBaseNormal()) * 0.5).normalized());

	mesh->addMeshPoint(point);
	// set border info
	//char bflags = face->getBorderFlags();
	//point->setBorder(bflags); // ie. ridge or normal boundary

	// rearange the mesh (i.e. "split" blocks)
	MeshBlock* blocks[2] = { face->getBlock(0), face->getBlock(1) };

	for(MeshBlock* block : blocks){
		MeshBlock* block_clone0 = block->clone();
		MeshBlock* block_clone1 = block->clone();
		block->switchPointsWithFaces(p0, point);
		block_clone0->switchPointsWithFaces(p1, point);
		block_clone1->switchPointsWithFaces(p2, point);
		mesh->addMeshBlock(block_clone0);
		mesh->addMeshBlock(block_clone1);
	}

	if (false) {
		SHOW_MESH("After split", mesh->getDebugViewSet(point));
	}

	return point;
}

/// Create triangulation of a set of points
MeshContainer3d* MeshGenerator3d::triangulatePoints(
	std::shared_ptr<DataVector<std::shared_ptr<MeshPoint3d>>> vec)
{
	if(!vec) return nullptr;

	DBox bbox;
	int pct = vec->countInt();
	for(int i = 0; i < pct; i++)
		bbox.addPoint(vec->get(i)->getCoordinates());

	// initial mesh
	auto cs = std::make_shared<ControlSpace3dIdentity>(bbox.getDiameter() * SMALL_NUMBER);
	MeshContainer3d* mesh = MeshGenerator3dDelaunayBoundary::createInitialMesh(pct, bbox, cs);
	if(!mesh) return nullptr;

	srand(0);

	Metric3dContext mc(cs);
	mc.countMetricAtPoint(bbox.getMiddlePoint());

	OctTreeMeshPoints octree_points(bbox);

	DataVector<int> sequence(pct);
	for(int i = 0; i < pct; i++)
		sequence.add(i);
	std::random_shuffle( sequence.begin(), sequence.end() );

	// * insert new points into the mesh
	for(int i = 0; i < pct; i++){
		// Create new point as a clone
		auto boundary_point = vec->get(sequence[i]);
		MeshPoint3d *point = new MeshPoint3d(*boundary_point);

		// Retriangulate
		if(!MeshGenerator3d::addPointToTriangulation(mc, mesh, point)){
			delete point;
//			delete mesh; 
//			return nullptr;
		}
	}

	return mesh;
}
