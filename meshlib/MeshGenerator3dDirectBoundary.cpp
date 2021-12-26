
#include <algorithm>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "MeshGenerator3dDirectBoundary.h"

#include "DPoint.h"
#include "DTetrahedron.h"
#include "MeshContainer3d.h"
#include "MeshEdge3d.h"
#include "MeshFace.h"
#include "MeshTriangle3d.h"
#include "MeshTetrahedron.h"
#include "Metric3dContext.h"
#include "ControlSpace3dAdaptive.h"
#include "MeshDomainVolume.h"
#include "MeshGenerator3d.h"
#include "MeshViewSet.h"
#include "DataHashTable.h"
#include "FrontFace.h"
#include "DataContainer.h"
#include "MeshContainer3dSurface.h"

class MeshBoundaryCondition;

/// Generate boundary constrained tetrahedral mesh
//MeshContainer3d* MeshGenerator3dDirectBoundary::createBoundaryConstrainedMesh(
//	Metric3dContext& mc, const DataVector<MeshFace*> & bfaces, 
//	const DataVector<MeshPoint3d*> & bpoints, MeshDomainVolume* mdv)
//{
//	LOG4CPLUS_INFO(MeshLog::logger_console, "Gen3d. Triangulating boundary (Direct)...");
//	START_CLOCK("MG3d::triangulateBoundary");
//
//	DPoint3d center;
//	double min_volume = minStarTetrahedraVolume(mc, bfaces, bpoints, mdv, center);
//	MeshContainer3d* mesh = nullptr;
//
//	const double MIN_VOL_WITH_SIZING = MeshGenerator3d::param_quality_threshold*MeshTetrahedron::getOptimumVolume();
//	const double MIN_VOL_WITHOUT_SIZING = 0.001;
//
//	if(min_volume >= MIN_VOL_WITH_SIZING)
//		mesh = createSimpleConvexBoundaryConstrainedMesh(mc, bfaces, bpoints, mdv, center);
////	if(!mesh) mesh = createFrontalBoundaryConstrainedMesh(mc, bfaces, bpoints, mdv);
//	if(!mesh && min_volume >= MIN_VOL_WITHOUT_SIZING)
//		mesh = createSimpleConvexBoundaryConstrainedMesh(mc, bfaces, bpoints, mdv, center);
//
//	STOP_CLOCK("MG3d::triangulateBoundary");
//
//	assert(!mesh || mesh->isValid());
//	return mesh;
//}

/// Check "star-shape" property of boundary mesh - by calculating minimum volume
double MeshGenerator3dDirectBoundary::minStarTetrahedraVolume(Metric3dContext& mc, 
	MeshContainer3dSurface* surface_mesh, MeshDomainVolume* mdv, DPoint3d& center)
{
	START_CLOCK("MG3dDB::minStarTetrahedraVolume");
	center = calculateInertialCenterFromSurface( surface_mesh );
	// TODO? -> check other points for better min_volume?

	auto freepoints = mdv->getFreePoints();
	if(freepoints && !freepoints->empty()){
		// select center as the freepoint closest to the oryginal center
		double min_dist2 = center.distance2(freepoints->get(0)->getCoordinates());
		for(size_t i = 1; i < freepoints->countInt(); i++){
			double dist2 = center.distance2(freepoints->get(i)->getCoordinates());
			if(dist2 < min_dist2){
				freepoints->switchData(0, i);
				min_dist2 = dist2;
			}
		}
		center = freepoints->get(0)->getCoordinates();			
	}

	// check potential tetrahedra volume and orientation using mdv reference
	double min_vol = 0.0;
	double ev_nt_count = 0.0;
	int fct = surface_mesh->getFacesCount();
	for(int i = 0; i < fct; i++){
		MeshFace* face = surface_mesh->getFaceAt(i);
		mc.countMetricAtPoint(face->getMiddlePoint());
		double vol = DTetrahedron::volume(
			face->getPoint(0)->getMetricCoordinates(mc),
			face->getPoint(1)->getMetricCoordinates(mc),
			face->getPoint(2)->getMetricCoordinates(mc),
			mc.transformRStoMS(center));
		if(face->getBlock(1) == mdv) vol = -vol;
		if(i == 0 || vol < min_vol) min_vol = vol;
		ev_nt_count += vol;
	}
	ev_nt_count *= 12.0/SQRT2;
	STOP_CLOCK("MG3dDB::minStarTetrahedraVolume");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Initial guess (minStarTetrahedraVolume): " << ev_nt_count);

	return min_vol;
}

/**
 * Create simple boundary constrained mesh for convex (star shaped) surface mesh
 */
MeshContainer3d* MeshGenerator3dDirectBoundary::createSimpleConvexBoundaryConstrainedMesh(
			Metric3dContext& mc, MeshContainer3dSurface* surface_mesh, 
			MeshDomainVolume* mdv, const DPoint3d& center)
{
	START_CLOCK("MG3dDB::createSimpleConvexBoundaryConstrainedMesh");
	// create mesh
	int pct = surface_mesh->getPointsCount();
	int fct = surface_mesh->getFacesCount();
	MeshContainer3d* mesh = new MeshContainer3d(fct);
	mesh->setConstrainingPhase(MeshData::CONSTRAIN_NONE);
	// -> add points
	DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> points_ref(2*pct, nullptr);
	for(int i = 0; i < pct; i++){
		MeshPoint3d* bpoint = surface_mesh->getPointAt(i);
		MeshPoint3d* new_point = new MeshPoint3d(*bpoint);
		points_ref.insert(bpoint, new_point);
		mesh->addMeshPoint(new_point);
		new_point->setBorder();
		new_point->setPtrTag(TagExtended::TAG_BOUNDARY_POINT, bpoint);
	}

	std::shared_ptr<OctTree> oct_tree;
	auto freepoints = mdv->getFreePoints();
	size_t fcount = freepoints ? freepoints->countInt() : 0;
	if(fcount && OctTree::param_octtree_threshold){
		DBox box = surface_mesh->getBoundingBox();
		oct_tree = std::make_shared<OctTree>(box);
		mesh->setSearchTree(oct_tree);
	}


	//  + center point
	MeshPoint3d* mpc = new MeshPoint3d(center);
	mesh->addMeshPoint(mpc);
	// -> add tetrahedra
	for(int i = 0; i < fct; i++){
		MeshFace* face = surface_mesh->getFaceAt(i);
		MeshPoint3d *bpts[3], *mpts[3];
		for(int j = 0; j < 3; j++){
			bpts[j] = face->getPoint(j);
			mpts[j] = points_ref.getValue(bpts[j], nullptr);
			assert(mpts[j] != nullptr);
		}
		MeshTetrahedron* block = nullptr;
		if(face->getBlock(0) == mdv){ // proper orientation
			block = new MeshTetrahedron(mpts[0], mpts[1], mpts[2], mpc);
		}else{
			block = new MeshTetrahedron(mpts[0], mpts[2], mpts[1], mpc);
		}
		block->setAreaID(mdv->getAreaID());
		mesh->addMeshTetrahedron(block);
		MeshFace* mface = block->getOppositeFace(mpc);
		mface->setBorder();
		for (int j = 0; j < 3; j++) {
			MeshEdge3d* bedge = bpts[j]->getEdgeToPoint(bpts[(j + 1) % 3]);
			assert(bedge != nullptr);
			MeshEdge3d* medge = mpts[j]->getEdgeToPoint(mpts[(j + 1) % 3]);
			assert(medge != nullptr);
			medge->copyBorderFlagsFrom(bedge);
			medge->setBorderFlags(TagBorder::OUTER);
			MeshBoundaryCondition* bc = (MeshBoundaryCondition*)bedge->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
			if (bc) medge->setPtrTag(TagExtended::TAG_BOUNDARY_COND, bc);
		}
		MeshBoundaryCondition* bc = (MeshBoundaryCondition*) face->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
		if(bc) mface->setPtrTag(TagExtended::TAG_BOUNDARY_COND, bc);
	}
	// final update
	mesh->setControlSpace(surface_mesh->getControlSpace());

	STOP_CLOCK("MG3dDB::createSimpleConvexBoundaryConstrainedMesh");

	// + freepoints
	if(!fcount){
		mesh->setConstrainingPhase(MeshData::CONSTRAIN_DONE);
		return mesh;
	}

	START_CLOCK("MG3dDB::createSimpleConvexBoundaryConstrainedMesh-freepoints");
	// freepoints[0] -> special case
	mpc->setBorder();
	mpc->setPtrTag(TagExtended::TAG_BOUNDARY_POINT, freepoints->get(0).get());

	fcount = freepoints->countInt() - 1;
	if(fcount > 0){
		DataVector<size_t> sequence(fcount);
		for(size_t i = 0; i < fcount; i++)
			sequence.add(1+i);
		std::random_shuffle( sequence.begin(), sequence.end() );


		for(size_t i = 0; i < fcount; i++){
			// Select one of points in random
			MeshPoint3d* new_point = new MeshPoint3d(*freepoints->get(sequence[i]));

			if(MeshGenerator3d::addPointToTriangulation(mc, mesh, new_point)){
				new_point->setBorder();
				new_point->setPtrTag(TagExtended::TAG_BOUNDARY_POINT, freepoints->get(sequence[i]).get());
			}else{
				// maybe try harder, searching for the containing tetrahedron...
				MeshTetrahedron* containing_tetrahedron = nullptr;
				for(int k = 0; k < mesh->getBlocksCount(); k++){
					MeshTetrahedron* tetrahedron = (MeshTetrahedron*)mesh->getBlockAt(k);
					if(tetrahedron->isPointInside(new_point->getCoordinates())){
						containing_tetrahedron = tetrahedron;
						break;
					}
				}
				if(containing_tetrahedron && 
					MeshGenerator3d::addPointToTriangulation(mc, mesh, new_point, containing_tetrahedron))
				{
					LOG4CPLUS_WARN(MeshLog::logger_console, "fixed no-containing-tetra for freepoint");
					new_point->setBorder();
					new_point->setPtrTag(TagExtended::TAG_BOUNDARY_POINT, freepoints->get(sequence[i]).get());
				}else{
					delete new_point;
				}
			}
		}
	}
	STOP_CLOCK("MG3dDB::createSimpleConvexBoundaryConstrainedMesh-freepoints");

	// finish
	mesh->setConstrainingPhase(MeshData::CONSTRAIN_DONE);
	return mesh;
}

/// Create boundary constrained mesh using frontal approach
MeshContainer3d* MeshGenerator3dDirectBoundary::createFrontalBoundaryConstrainedMesh(
		Metric3dContext& mc, MeshContainer3dSurface* surface_mesh, const MeshDomainVolume* mdv)
{
	return nullptr; // not ready yet -> still problem with proper detection of collisions

	// create mesh
	int pct = surface_mesh->getPointsCount();
	int fct = surface_mesh->getFacesCount();
	MeshContainer3d* mesh = new MeshContainer3d(fct);
	// -> create & add points
	DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> points_ref(2*pct, nullptr);
	for(int i = 0; i < pct; i++){
		MeshPoint3d* bpoint = surface_mesh->getPointAt(i);
		MeshPoint3d* new_point = new MeshPoint3d(*bpoint);
		points_ref.insert(bpoint, new_point);
		mesh->addMeshPoint(new_point);
		new_point->setBorder();
		new_point->setPtrTag(TagExtended::TAG_BOUNDARY_POINT, bpoint);
	}
	// -> create faces
	DataContainer<FrontFace> front_faces(fct);

	for(int i = 0; i < fct; i++){
		MeshFace* bface = surface_mesh->getFaceAt(i);
		MeshPoint3d* pts[3] = {
			points_ref.getValue(bface->getPoint(0), nullptr),
			points_ref.getValue(bface->getPoint(1), nullptr),
			points_ref.getValue(bface->getPoint(2), nullptr) 
		};
		// create face
		MeshFace* mface = nullptr;
		if(bface->getBlock(0) == mdv)
			mface = new MeshTriangle3d(pts[0], pts[1], pts[2]);
		else
			mface = new MeshTriangle3d(pts[0], pts[2], pts[1]);
		mface->setBorderWithEdges();
		MeshBoundaryCondition* bc = (MeshBoundaryCondition*) bface->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
		if(bc) mface->setPtrTag(TagExtended::TAG_BOUNDARY_COND, bc);
		// create front face
		front_faces.addDataItem(new FrontFace(mface, true));
		if(bface->getBlock(0) == mdv && bface->getBlock(1) == mdv) // face with the same block on both sides
			front_faces.addDataItem(new FrontFace(mface, false));
	}

	// -> gather and classify faces as front-faces
	for(int i = 0; i < front_faces.countInt(); i++)
		front_faces.getDataAt(i)->init(mc);
	front_faces.setHeapOrder();

	// check direct-greedy procedure with heap
	for(int i_face = 0; i_face < front_faces.countInt(); i_face++){

		if(true){
			// visualization
			MeshViewSet* set = new MeshViewSet;
			// ... boundary edges
			// ... vertices
			for(int i = 0; i < pct; i++){
				MeshPoint3d* bpoint = surface_mesh->getPointAt(i);
				int rank = bpoint->getRank();
				set->addPoint(bpoint->getCoordinates(), 0, i);
				for(int j = 0; j < rank; j++){
					MeshEdge3d* edge = bpoint->getEdge(j);
					if(edge->getPointIndex(bpoint) == 0)
						set->addEdge(edge, 2);
				}
			}
			// ... front-faces
			for(int i = 0; i < front_faces.countInt(); i++){
				FrontFace* ff = front_faces.getDataAt(i);
				set->addFace(ff->getFace(), 1, 0.95, ff->isOriented());
			}
			SHOW_MESH("frontal boundary conforming mesh creation", set);
		}

		FrontFace* fface = front_faces.getDataAt(i_face); // find first valid candidature
		// check angle - if too large, continue;
		int min_i = fface->getMinAngleIndex();
		if(fface->getAngle(min_i) > PI) continue;
		// identify vertices (in proper order)
		MeshFace* face = fface->getFace();
		MeshEdge3d* edge = face->getEdge(min_i);
		MeshPoint3d* mpoints[4];
		mpoints[2] = face->getOtherPoint(
			edge->getMeshPoint(0), edge->getMeshPoint(1));
		mpoints[3] = fface->getIncidentFace(min_i)->getOtherPoint(
			edge->getMeshPoint(0), edge->getMeshPoint(1));
		if(fface->isOriented()){
			mpoints[0] = face->getNextPoint(mpoints[2]);
			mpoints[1] = face->getPrevPoint(mpoints[2]);
		}else{
			mpoints[0] = face->getPrevPoint(mpoints[2]);
			mpoints[1] = face->getNextPoint(mpoints[2]);
		};
		// check side faces - if collision, continue
		// ... identify side front-faces
		FrontFace* other_fface = fface->getNeighbor(min_i);
		FrontFace* ff12 = nullptr; // side faces of this_fface
		FrontFace* ff02 = nullptr;
		FrontFace* ff13 = nullptr; // side faces of other_fface
		FrontFace* ff03 = nullptr;
		for(int j = 0; j < 3; j++){
			FrontFace* ff = fface->getNeighbor(j);
			if(ff != other_fface){
				if(ff->getPoint(0) == mpoints[0] || ff->getPoint(1) == mpoints[0] || ff->getPoint(2) == mpoints[0])
					ff02 = ff;
				else ff12 = ff;
			}
			ff = other_fface->getNeighbor(j);
			if(ff != fface){
				if(ff->getPoint(0) == mpoints[0] || ff->getPoint(1) == mpoints[0] || ff->getPoint(2) == mpoints[0])
					ff03 = ff;
				else ff13 = ff;
			}
		}
		assert(ff02 && ff12 && ff03 && ff13);
		mc.countMetricAtPoint(DPoint3d::average(
			mpoints[0]->getCoordinates(),
			mpoints[1]->getCoordinates(),
			mpoints[2]->getCoordinates()));
		// volume of the future tetrahedron - can't be negative (or too small)
		double vol0123 = DTetrahedron::volume(
			mpoints[0]->getMetricCoordinates(mc),
			mpoints[1]->getMetricCoordinates(mc),
			mpoints[2]->getMetricCoordinates(mc),
			mpoints[3]->getMetricCoordinates(mc));
		if(vol0123 < METRIC_SMALL_NUMBER) continue;

		if(true){
			MeshViewSet* set = new MeshViewSet;
			// ... boundary edges
			for(int i = 0; i < pct; i++){
				MeshPoint3d* bpoint = surface_mesh->getPointAt(i);
				int rank = bpoint->getRank();
				for(int j = 0; j < rank; j++){
					MeshEdge3d* tedge = bpoint->getEdge(j);
					if(tedge->getPointIndex(bpoint) == 0)
						set->addEdge(tedge, 2);
				}
			}
			// ... points
			for(int i = 0; i < 4; i++)
				set->addPoint(mpoints[i]->getCoordinates(), 0, mpoints[i]->getIndex());
			// ... faces
			set->addFace(fface->getFace(), 0, 1.0, fface->isOriented());
			set->addFace(other_fface->getFace(), 1, 1.0, other_fface->isOriented());
			// ... tetrahedra
			for(int i = 0; i < mesh->getBlocksCount(); i++)
				set->addBlock(mesh->getBlockAt(i), 3);
			SHOW_MESH("front-faces", set);
		}

		if(ff12 != ff13){ // if no common side-face
			// --> if (angle >= PI) skip test
			if(fface->getAngle(ff12) < PI){
				double vol12 = DTetrahedron::volume(
					mpoints[1]->getMetricCoordinates(mc),
					mpoints[2]->getMetricCoordinates(mc),
					mpoints[3]->getMetricCoordinates(mc),
					ff12->getFace()->getOtherPoint(mpoints[1], mpoints[2])->getMetricCoordinates(mc));
				if(vol12 < METRIC_SMALL_NUMBER) continue;
			}
			if(other_fface->getAngle(ff13) < PI){
				double vol13 = DTetrahedron::volume(
					mpoints[1]->getMetricCoordinates(mc),
					mpoints[2]->getMetricCoordinates(mc),
					mpoints[3]->getMetricCoordinates(mc),
					ff13->getFace()->getOtherPoint(mpoints[1], mpoints[3])->getMetricCoordinates(mc));
				if(vol13 < METRIC_SMALL_NUMBER) continue;
			}
		}
		// similar for ff02 and ff03
		if(ff02 != ff03){ // if no common side-face
			if(fface->getAngle(ff02) < PI){
				double vol02 = DTetrahedron::volume(
					mpoints[0]->getMetricCoordinates(mc),
					mpoints[3]->getMetricCoordinates(mc),
					mpoints[2]->getMetricCoordinates(mc),
					ff02->getFace()->getOtherPoint(mpoints[0], mpoints[2])->getMetricCoordinates(mc));
				if(vol02 < METRIC_SMALL_NUMBER) continue;
			}
			if(fface->getAngle(ff03) < PI){
				double vol03 = DTetrahedron::volume(
					mpoints[0]->getMetricCoordinates(mc),
					mpoints[3]->getMetricCoordinates(mc),
					mpoints[2]->getMetricCoordinates(mc),
					ff03->getFace()->getOtherPoint(mpoints[0], mpoints[3])->getMetricCoordinates(mc));
				if(vol03 < METRIC_SMALL_NUMBER) continue;
			}
		}

		if(mpoints[2]->getEdgeToPoint(mpoints[3]) == nullptr){
			// ... check for crossing edges ...
			// TODO current solution doesn't work
			bool ok = true;
			for(int i = 0; i < front_faces.countInt(); i++){
				FrontFace* ff = front_faces.getDataAt(i);
				if(i == i_face || ff == other_fface ||
					ff == ff12 || ff == ff13 || ff == ff02 || ff == ff03) continue;
				MeshFace* mf = ff->getFace();
				if( DMTriangle3d::crossSegment(
						mpoints[2]->getMetricCoordinates(mc),
						mpoints[3]->getMetricCoordinates(mc),
						mf->getPoint(0)->getMetricCoordinates(mc),
						mf->getPoint(1)->getMetricCoordinates(mc),
						mf->getPoint(2)->getMetricCoordinates(mc)))
				{
					ok = false;
					break;
				}
			}
			if(!ok) continue;
		}

		// check ok -> create tetrahedron
		MeshTetrahedron* tetra = new MeshTetrahedron(
			mpoints[0], mpoints[1], mpoints[2], mpoints[3]);
		tetra->setAreaID(mdv->getAreaID());
		mesh->addMeshBlock(tetra);

		// remove obsolete front faces
		for(int i = 0; i < 4; i++){
			FrontFace* ff = nullptr;
			MeshFace* tface = tetra->getFace(i);
			if(tface->getBlock(0) == tetra && (ff = (FrontFace*)tface->getPtrTag(TagExtended::TAG_FRONT_1))){
				delete front_faces.removeDataItem(ff->getIndex());
			}else if(tface->getBlock(1) == tetra && (ff = (FrontFace*)tface->getPtrTag(TagExtended::TAG_FRONT_2))){
				delete front_faces.removeDataItem(ff->getIndex());
			}
		}

		// create new front faces & update neighbors
		DataVector<FrontFace*> new_ffaces;
		for(int i = 0; i < 4; i++){
			MeshFace* tface = tetra->getFace(i);
			if(tface->isBorder() || tface->getOtherBlock(tetra) != nullptr) continue;
			new_ffaces.add(new FrontFace(tface, tface->getBlock(1) == tetra));
		}

		if(true){
			// visualization
			MeshViewSet* set = new MeshViewSet;
			// ... new front faces
			for(size_t i = 0; i < new_ffaces.countInt(); i++){
				FrontFace* ff = new_ffaces[i];
				set->addFace(ff->getFace(), 0, 0.95, ff->isOriented());
			}
			// ... tetrahedra
			set->addBlock(tetra, 3);
			SHOW_MESH("new front faces", set);
		}

		// update connectivity for new front-faces
		size_t new_ff_ct = new_ffaces.countInt();
		for(size_t i = 0; i < new_ff_ct; i++)
			new_ffaces[i]->init(mc);

		// gather neighbors of the new front-faces
		for(size_t i = 0; i < new_ff_ct; i++)
			new_ffaces[i]->gatherNeighbors(new_ffaces);

		// update connectivity for old front-faces adjacent to new front-faces
		for(size_t i = new_ff_ct; i < new_ffaces.countInt(); i++){
			front_faces.removeDataItem(new_ffaces[i]->getIndex());
			new_ffaces[i]->init(mc);
		}

		// update heap
		for(size_t i = 0; i < new_ffaces.countInt(); i++)
			front_faces.addDataItem(new_ffaces[i]);

		// should be OK...

		i_face = 0; // start checking front-face heap from the beginning
	}

	// for now
	//delete mesh;

	if(front_faces.countInt() == 0){
		return mesh;
	}else{
		front_faces.deleteAll();
		delete mesh;
		return nullptr;
	}
}

/// calculate inertial center from the surface mesh (weighted)
DPoint3d MeshGenerator3dDirectBoundary::calculateInertialCenterFromSurface(
	MeshContainer3dSurface* surface_mesh)
{
	DPoint3d center;

	int pct = surface_mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		const DPoint3d& pt = surface_mesh->getPointAt(i)->getCoordinates();
		center.add(pt);
	}
	center /= pct;

	return center;
}

