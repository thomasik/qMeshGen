// MeshSpecialRoutinesUTC.cpp: implementation of the MeshSpecialRoutinesUTC class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshData.h"
#include "MeshSpecialRoutinesUTC.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshBlock.h"
#include "MeshTetrahedron.h"
#include "MeshTriangle3d.h"
#include "MeshEdge3d.h"
#include "MeshViewSet.h"
#include "ControlSpace3dOctree.h"
#include "MeshGenerator3d.h"
#include "MeshGenerator3dAdapt.h"
#include "MeshGenerator3dQuality.h"
#include "MeshSplit3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshGenerator3dSurface.h"
#include "MeshDomainVolume.h"
#include "DLeastSquaresFitting.h"
#include "SurfacePlane.h"
#include "GeometricPredicates.h"
#include "DPlane.h"
#include "MeshLog.h"

/// marks elements within the crack (described by box), return number of marked elements
int MeshSpecialRoutinesUTC::markCrack(MeshContainer3d* mesh, const DBox& box)
{
	int pct = mesh->getPointsCount();
	DataVector<bool> marked_points(pct, false);
	for(int i = 0; i < pct; i++)
		if(box.contains(mesh->getPointAt(i)->getCoordinates()))
			marked_points[i] = true;

	int bct = mesh->getBlocksCount();
	int marked_count = 0;
	for(int i = 0; i < bct; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		int bpct = block->getPointCount();
		for(int j = 0; j < bpct; j++){
			if(marked_points[block->getPoint(j)->getIndex()]){
				block->setAreaID(2);
				block->setIntTag(TagExtended::TAG_CUT_ADAPT, 1);
				marked_count++;
				break;
			}
		}
	}

	return marked_count;
}

/// retriangulate given mesh with next-step-CS and smoothing of crack surfaces
void MeshSpecialRoutinesUTC::remeshCrack(Metric3dContext& mc, MeshContainer3d* mesh, int mode)
{
	// load / prepare next-step-CS
	//  -> refine ahead
	//  -> smoothe crack surfaces
	//  -> coarsen behind

	if(mode == 0)
		MeshGenerator3dAdapt::remeshWithLocalCutRegions(mc, mesh);
	else
		MeshGenerator3dAdapt::remeshWithLocalTransformations(mc, mesh);

/*
	int pct = mesh->getPointsCount();
	int bct = mesh->getBlocksCount();
	MeshViewSet* set = new MeshViewSet(0, 2*bct, 0, 0);
	DataVector<double> point_len(pct, 0.0);
	DataVector<int> point_w(pct, 0);
	DataHashTable<MeshEdge3d*> visited_edges(4*bct, nullptr);
	DataVector<MeshFace*> crack_faces(3*bct);

	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
		MeshFace* face = it.getFace();
		if(face->availableTag(TagExtended::TAG_ADAPT_SURF)){
			crack_faces.add(face);
			int fect = face->getEdgeCount();
			for(int i = 0; i < fect; i++){
				MeshEdge3d* edge = face->getEdge(i);
				if(visited_edges.insert(edge)){
					double len = edge->getLength(mc, true);
					for(int j = 0; j < 2; j++){
						MeshPoint3d* point = edge->getMeshPoint(j);
						int ind = point->getIndex();
						point_len[ind] += len;
						point_w[ind]++;
					}
				}
			}
		}
	}
	for(int i = 0; i < pct; i++)
		if(point_w[i] > 0) point_len[i] /= point_w[i];
	DataHashTable<MeshPoint3d*> hpoints(4*crack_faces.countInt(), nullptr);
	for(int i = 0; i < crack_faces.countInt(); i++){
		MeshFace* face = crack_faces[i];
		double max_len = 0.0;
		int fpct = face->getPointCount();
		for(int j = 0; j < fpct; j++){
			double len = point_len[face->getPoint(j)->getIndex()];
			if(len > max_len) max_len = len;
		}
		if(max_len < 0.5){
			face->setIntTag(TagExtended::TAG_ADAPT_SURF, 2);
			for(int j = 0; j < fpct; j++)
				hpoints.insert(face->getPoint(j));
			set->addFaceWithEdges(face, 1, MeshViewSet::param_shrink, face->getBlock(0) != nullptr);
		}else{
			face->setIntTag(TagExtended::TAG_ADAPT_SURF, 1);
			set->addEdges(face);
		}
	}
	SHOW_MESH("crack faces", set);

	DataVector<MeshPoint3d*> crack_points;
	hpoints.getValues(crack_points);
	MeshContainer3d* crack_surface_mesh = new MeshContainer3d(crack_points.countInt());
	const DVector3d DV(0.07, 0.0, 0.0);
	DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> hmpoints(2*crack_points.countInt(), nullptr);
	for(int i = 0; i < crack_points.countInt(); i++){
		// duplicate poins (with slight translation for debug)
		MeshPoint3d* point = new MeshPoint3d(* crack_points[i]);
		point->setCoordinates( point->getCoordinates() + DV );
		crack_surface_mesh->addMeshPoint(point);
		hmpoints.insert(crack_points[i], point);
	}

	set = new MeshViewSet(0, crack_faces.countInt(), 0, 0);

	for(int i = 0; i < crack_faces.countInt(); i++){
		// create duplicate faces with adjusted orientation
		MeshFace* face = crack_faces[i];
		assert(face->getPointCount() == 3);
		if(face->getIntTag(TagExtended::TAG_ADAPT_SURF, 0) != 2) continue;
		MeshFace* mface = (face->getBlock(0) != nullptr) 
			? new MeshTriangle3d(
					hmpoints.getValue(face->getPoint(0), nullptr),
					hmpoints.getValue(face->getPoint(1), nullptr),
					hmpoints.getValue(face->getPoint(2), nullptr))
			: new MeshTriangle3d(
					hmpoints.getValue(face->getPoint(0), nullptr),
					hmpoints.getValue(face->getPoint(2), nullptr),
					hmpoints.getValue(face->getPoint(1), nullptr));
		mface->setBorder();
		set->addFace(mface);
	}

	for(IteratorEdge3d it = crack_surface_mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
		if(it.getEdge()->getFaceCount() == 1)
			it.getEdge()->setBorder();
	}

	SHOW_MESH("Duplicated crack surface", set);

	delete crack_surface_mesh;
*/
}

/// prepare (analytically, for now) control space for next simulation step
CS3dAPtr MeshSpecialRoutinesUTC::prepareNextCS(MeshContainer3d* mesh,
		const DataCompoundList<DTriangle3d>& crack, const DVector3d& dv)
{
	const double DMAX = 0.3;
	const double DMIN = 0.003;

//#define NEXT_CS_WITH_CURVATURE

#ifdef NEXT_CS_WITH_CURVATURE
	// 1. Re-create source control space (sCS) - directly from surface mesh elements
	MeshContainer3dSurface* mesh_surface = MeshGenerator3dSurface::copySurfaceMeshFromVolumeMesh(mesh);
	mesh_surface->clearBoundaryFlags();
	mesh_surface->setBoundaryFeatureEdges();
	mesh_surface->setBoundarySharpEdges(MeshGenerator3dSurface::param_sharp_edge_threshold);
	CS3dPtr source_cs = MeshGenerator3dSurface::createACSfromMeshFaces(mesh_surface);
	// 2. Use this sCS to identify parameterization surfaces
	Metric3dContext source_mc(source_cs);
	mesh_surface->identifyLocalSurfaces(source_mc, 0.2);
	// 3. Create adapted control space 	- using:
	//		3a. global maximum
	//		3b. curvature of parameterization surfaces
	//		3c. distance between ridge vertices?
	CS3dPtr space = new ControlSpace3dOctree(source_cs->getBoundingBox());
	space->setGlobalMetric(DMetric3d::stretchToMatrix(ControlDataStretch3d(DMAX)));

	MeshGenerator3dSurface::updateACSwithLocalCurvature(space, mesh_surface);
	space->smoothen();

#else // ! NEXT_CS_WITH_CURVATURE

	DBox box;
	for(int i = 0; i < mesh->getPointsCount(); i++)
		box.addPoint(mesh->getPointAt(i)->getCoordinates());
	box.inflate(ControlSpace2d::param_inflate_box_factor);
	auto space = std::make_shared<ControlSpace3dOctree>(box);

	space->setGlobalMetric(DMetric3d::stretchToMatrix(ControlDataStretch3d(DMAX)));

#endif // NEXT_CS_WITH_CURVATURE

	// 4. Insert additional discrete sources
	for (auto it = crack.iterator(); it.valid(); it.moveNext()) {
		const DTriangle3d& triangle = it.item();
		space->setMinControl(ControlDataExtMatrix3dTriangle(
				triangle.pt_a + dv, triangle.pt_b - triangle.pt_a, triangle.pt_c - triangle.pt_a, 2*DMIN, 
			DMetric3d::stretchToMatrix(ControlDataStretch3d(DMIN))));
	}
	
	// ... postprocess
	space->smoothen();
	space->compact();
	return space;
}

/// identify set of 3d triangles approximating the surface of crack
int MeshSpecialRoutinesUTC::identifyCrackPlanes(Metric3dContext& mc, 
	const DataVector<MeshBlock*> & crack_blocks, double tolerance,
	DataVector<std::shared_ptr<DTriangle3d>> & crack_planes, 
	DataVector<int> & crack_indices)
{
	int crack_ct = crack_blocks.countInt();
	DataVector<DPoint3d>  crack_centers(crack_ct);
	for(int i = 0; i < crack_ct; i++){
		crack_centers.add( crack_blocks[i]->getMiddlePoint() );
	}

	DPlane plane;
	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(
		crack_centers, plane);
	mc.countMetricAtPoint(plane.p0);
	double metric_f = mc.transformRStoMS(plane.vn).length();
	double metric_max_dist = plane_max_dist * metric_f;
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Crack surface (plane), max_dist = " << metric_max_dist);
	if(metric_max_dist <= tolerance){
		DRect brect; // ... assuming rectangular area, split into two triangles
		for(int i = 0; i < crack_ct; i++){
			DPoint2d pt = plane.projectToPlane(crack_centers[i]);
			brect.addPoint(pt);
		}
		DPoint2d ptA = brect.getX0Y1();
		DPoint2d ptB = brect.getX1Y0();
		auto t0 = std::make_shared<DTriangle3d>( plane.projectToSpace(ptA), plane.projectToSpace(brect.getX0Y0()), plane.projectToSpace(ptB) );
		auto t1 = std::make_shared<DTriangle3d>( plane.projectToSpace(ptB), plane.projectToSpace(brect.getX1Y1()), plane.projectToSpace(ptA) );
		crack_planes.add(t0);
		crack_planes.add(t1);
		for(int i = 0; i < crack_ct; i++){
			DPoint2d pt = plane.projectToPlane(crack_centers[i]);
			crack_indices[i] = (GeometricPredicates::orient2d(ptA, ptB, pt) < 0) ? 0 : 1;
		}

		//MeshViewSet* set = new MeshViewSet;
		//set->addEdge(t0->pt_a, t0->pt_b);
		//set->addEdge(t0->pt_b, t0->pt_c);
		//set->addEdge(t0->pt_c, t0->pt_a);
		//set->addEdge(t1->pt_a, t1->pt_b, 2);
		//set->addEdge(t1->pt_b, t1->pt_c, 2);
		//set->addEdge(t1->pt_c, t1->pt_a, 2);
		//for(int i = 0; i < crack_ct; i++){
		//	if(crack_indices[i] == 0) set->addPoint(crack_centers[i]);
		//	else set->addPoint(crack_centers[i], 2);
		//}
		//SHOW_MESH("crack surface points", set);
		return 2;
	}else{ 
		return 0; // ... temporary, TODO -> identify successive triangles / patches
	}

//		if(metric_max_dist <= tolerance){
//			surface = new SurfacePlane(plane_pt, plane_e0, plane_e1);
//		}

	return crack_planes.countInt();
}

/// remove elements within the crack (described by set of (triangle) planes), return number of removed elements
int MeshSpecialRoutinesUTC::createCrack(Metric3dContext & mc, MeshContainer3d* mesh, const DataCompoundList<DTriangle3d>& crack)
{
	DataVector<MeshBlock*> crack_blocks(100);
	int crack_ct = MeshSpecialRoutinesUTC::gatherCrackBlocks(mesh, crack, crack_blocks);

	// ...

	// 1. create set of planes/triangles approximating the centers of crack_blocks
	DataVector<std::shared_ptr<DTriangle3d>> crack_planes(50);
	DataVector<int> crack_indices(crack_ct, 0);
	int planes_ct = MeshSpecialRoutinesUTC::identifyCrackPlanes(mc, crack_blocks, 0.5, crack_planes, crack_indices);

	// 2. create offset-planes and ascribe offset-planes to crack boundary surfaces
	DataHashTableKeyValue<MeshFace*,int> hash_faces(6*crack_ct, nullptr);
	for(int i = 0; i < planes_ct; i++){
		DPoint3d ptA = crack_planes[i]->pt_a;
		DVector3d e0 = crack_planes[i]->pt_b - ptA;
		DVector3d dn = e0.crossProduct(crack_planes[i]->pt_c - ptA).normalized();
		DVector3d e1 = dn.crossProduct(e0);
		double D =  -(dn.x * ptA.x + dn.y * ptA.y + dn.z * ptA.z);
		//SurfacePlane plane(pt_x, e0, e1);

		//MeshViewSet* sets[3] = {
		//	new MeshViewSet, new MeshViewSet, new MeshViewSet 
		//};

		double sdist[2] = {0.0, 0.0};
		DataVector<MeshFace*> up_faces(crack_ct);
		DataVector<MeshFace*> down_faces(crack_ct);
		DataVector<MeshFace*> * sfaces[2] = { &up_faces, &down_faces };
		DataVector<MeshFace*> other_faces(crack_ct);

		const double e = mesh_data.relative_small_number;
		for(int j = 0; j < crack_ct; j++){
			MeshBlock* block = crack_blocks[j];
			if(crack_indices[j] != i) continue; // only for current plane
			int bfct = block->getFaceCount();
			for(int k = 0; k < bfct; k++){
				MeshFace* face = block->getFace(k);
				if(face->isBorder() || hash_faces.contains(face)) continue;
				MeshBlock* other_block = face->getOtherBlock(block);
				if(other_block->availableTag(TagExtended::TAG_ADAPT_CRACK)) 
					continue;
				// ok
				double dist[3];
				for(int m = 0; m < 3; m++){
					DPoint3d pt = face->getPoint(m)->getCoordinates();
					dist[m] = dn.x * pt.x + dn.y * pt.y + dn.z * pt.z + D;
				}
				if(dist[0] > -e && dist[1] > -e && dist[2] > -e){
					// up-plane
					//sets[0]->addFace(face);
					hash_faces.insert(face, 2*i);
					sfaces[0]->add(face);
					sdist[0] += (dist[0] + dist[1] + dist[2]);
				}else if(dist[0] < e && dist[1] < e && dist[2] < e){
					// down-plane
					//sets[1]->addFace(face);
					hash_faces.insert(face, 2*i+1);
					sfaces[1]->add(face);
					sdist[1] += (dist[0] + dist[1] + dist[2]);
				}else{
					// nothing
					//sets[2]->addFace(face);
					hash_faces.insert(face, -1);
					other_faces.add(face);
				}
			}
		}

		while(other_faces.notEmpty()){
			DataVector<MeshFace*> face_set(100);
			face_set.add(other_faces.removeLast());
			int k0 = 0, k1 = 1;
			int other_tag = -1;
			// gather adjacent crack-faces
			for(int j = 0; j < 4; j++){
				for(int k = k0; k < k1; k++){
					MeshFace* face = face_set[k];
					for(int m = 0; m < 3; m++){
						MeshEdge3d* edge = face->getEdge(m);
						int efct = edge->getFaceCount();
						for(int l = 0; l < efct; l++){
							MeshFace* eface = edge->getFaceAt(l);
							if(eface == face) continue;
							int tag = hash_faces.getValue(eface, -2);
							if(tag == -2) continue;
							else if(tag == -1){
								size_t id = other_faces.find(eface);
								if(id != std::string::npos){
									other_faces.removeAt(id);
									face_set.add(eface);
								}
							}else{
								if(other_tag == -1) other_tag = tag;
								else if(other_tag != tag) other_tag = -2; // 
							}
						}
					}
				}
				k0 = k1;
				k1 = face_set.countInt();
			}
			// check and decide
			if(other_tag >= 0){ // whole set of -1 faces should be actually part of up/down surface
				int m = other_tag % 2;
				for(int j = 0; j < k1; j++){
					sfaces[m]->add(face_set[j]);
					hash_faces.setValue(face_set[j], other_tag);
					for(int k = 0; k < 3; k++){
						DPoint3d pt = face_set[j]->getPoint(k)->getCoordinates();
						sdist[m] += dn.x * pt.x + dn.y * pt.y + dn.z * pt.z + D;
					}
				}
			}else{
				//for(int j = 0; j < k1; j++){
				//	face_set[j]->setBorderWithEdgesAndPoints(TagBorder::FIXED);
				//	//sets[2]->addFace(face_set[j]);
				//}
			}
		}

		for(int m = 0; m < 2; m++){
			int sct = sfaces[m]->countInt();
			if(sct > 0){
				double ave_dist = sdist[m] / (3*sct);
				SurfacePtr plane(new SurfacePlane(ptA + dn*ave_dist, e0, e1));
				mesh->addLocalSurface(plane);
				// ... faces
				DataHashTable<MeshPoint3d*> hash_points(6*sct, nullptr);
				for(int j = 0; j < sct; j++){
					MeshFace* face = sfaces[m]->get(j);
					//sets[m]->addFace(face);
					face->setLocalSurface(plane);
					for(int k = 0; k < 3; k++)
						hash_points.insert(face->getPoint(k));
				}
				// ... points
				int spct = hash_points.countInt();
				DataVector<MeshPoint3d*> spoints(spct);
				hash_points.getValues(spoints);
				DataHashTableKeyValue< SurfaceParametricSet*, SurfaceParametricSet* > extended_sets(100, nullptr);

				// TODO_CURVES: after removing SurfaceSet, implement LOCAL_CURVE for boundary edges with >1 adjacent faces
				LOG4CPLUS_INFO(MeshLog::logger_console, 
					"TODO: implement local_curve for boundary conoturs of the crack...");

				for(int j = 0; j < spct; j++){
					MeshPoint3d* point = spoints[j];
					SurfaceConstPtr psurface = point->getLocalValidSurface();
					if(!psurface) point->setLocalSurface(plane, plane->getParameters( point->getCoordinates() ), 1.0 );
				}
			}
		}

		//for(int j = 0; j < 3; j++){
		//	sets[j]->addEdge(crack_planes[i]->pt_a, crack_planes[i]->pt_b);
		//	sets[j]->addEdge(crack_planes[i]->pt_b, crack_planes[i]->pt_c);
		//	sets[j]->addEdge(crack_planes[i]->pt_c, crack_planes[i]->pt_a);
		//	SHOW_MESH("test face-plane assign", sets[j]);
		//}
	}

	// Remove crack blocks
	for(int i = 0; i < crack_ct; i++){
		MeshBlock * block = crack_blocks[i];
		int bfct = block->getFaceCount();
		for(int j = 0; j < bfct; j++)
			block->getFace(j)->setBorderWithEdgesAndPoints( TagBorder::NONE );
		delete mesh->removeMeshBlock(block);
	}
	// ... remove stranded points (if any)
	int pct = mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if(point->getRank() == 0){
			delete mesh->removeMeshPoint(point);
			i--;
			pct--;
		}
	}

//	MeshViewSet* set = new MeshViewSet(0, 0, 1000, 0);
	// ... fix boundary marking ?
	DataSimpleList<MeshFace*> crack_faces;
	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
		MeshFace* face = it.getFace();
		if(face->isBorder() || face->isBoundedBothSides()) 
			continue;
		crack_faces.append(face);
		face->setBorderWithEdgesAndPoints();
		face->setIntTag(TagExtended::TAG_ADAPT_CRACK, 1);
		face->setIntTagForPoints(TagExtended::TAG_ADAPT_CRACK, 1);
		face->setIntTagForPoints(TagExtended::TAG_ADAPT_VOLUME, 1);
//		set->addFace(face, -2, MeshViewSet::param_shrink, face->getBlock(0) != nullptr);
	}
//	SHOW_MESH("crack", set);

	DataSimpleList<MeshEdge3d*> ridge_edges;
	while(crack_faces.notEmpty()){
		MeshFace* face = crack_faces.removeFirst();
		if(!face->hasLocalSurface())
			face->setBorderWithEdgesAndPoints( TagBorder::OUTER | TagBorder::FIXED);
		int fect = face->getEdgeCount();
		for(int i = 0; i < fect; i++){
			MeshEdge3d* edge = face->getEdge(i);
			int efct = edge->getFaceCount();
			for(int j = 0; j < efct; j++){
				MeshFace* f = edge->getFaceAt(j);
				if(!f->isBorder() || f == face) continue;
				if(!f->availableTag(TagExtended::TAG_ADAPT_CRACK)){
					edge->setBorderFlags( TagBorder::RIDGE );
					edge->getMeshPoint(0)->setBorderFlags( TagBorder::RIDGE );
					edge->getMeshPoint(1)->setBorderFlags( TagBorder::RIDGE );
					ridge_edges.append(edge);
					break;
				}
			}
		}
	}

	while(ridge_edges.notEmpty()){
		MeshEdge3d* edge = ridge_edges.removeFirst();
		for(int i = 0; i < 2; i++){
			MeshPoint3d* point = edge->getMeshPoint(i);
			if(point->isBorder( TagBorder::CORNER )) continue;
			int rct = 0;
			int ect = point->getRank();
			for(int j = 0; j < ect; j++){
				if(point->getEdge(j)->isBorder( TagBorder::RIDGE )) rct++;
			}
			if(rct > 2) point->setBorderFlags( TagBorder::CORNER );
		}
	}

//	SHOW_MESH("createCrack", mesh->getViewSet());
	return crack_ct;
}

int MeshSpecialRoutinesUTC::gatherCrackBlocks(MeshContainer3d* mesh, 
		const DataCompoundList<DTriangle3d>& crack, DataVector<MeshBlock*> & crack_blocks)
{
	int pct = mesh->getPointsCount();
	int bct = mesh->getBlocksCount();

	DataVector<bool> crack_mark_points(pct, false);

	for(int i = 0; i < bct; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		assert(block->getType() == BLOCK_TETRA);
		const DPoint3d& pt0 = block->getPoint(0)->getCoordinates();
		const DPoint3d& pt1 = block->getPoint(1)->getCoordinates();
		const DPoint3d& pt2 = block->getPoint(2)->getCoordinates();
		const DPoint3d& pt3 = block->getPoint(3)->getCoordinates();
		for (auto it = crack.iterator(); it.valid(); it.moveNext()) {
			const DTriangle3d& triangle = it.item();
			//bool up = false, down = false;
			//if(triangle.orient3d(pt0) >= 0.0) up = true; else down = true;
			//if(triangle.orient3d(pt1) >= 0.0) up = true; else down = true;
			//if(triangle.orient3d(pt2) >= 0.0) up = true; else down = true;
			//if(triangle.orient3d(pt3) >= 0.0) up = true; else down = true;
			//if(up && down){
			if( triangle.crossSegment(pt0, pt1) || triangle.crossSegment(pt0, pt2) ||
				triangle.crossSegment(pt0, pt3) || triangle.crossSegment(pt1, pt2) ||
				triangle.crossSegment(pt1, pt3) || triangle.crossSegment(pt2, pt3))
			{
				crack_blocks.add(block);
				block->setIntTag(TagExtended::TAG_ADAPT_CRACK);
				for(int j = 0; j < block->getPointCount(); j++)
					crack_mark_points[ block->getPoint(j)->getIndex() ];
				break;
			}
		}
	}
	LOG4CPLUS_INFO(MeshLog::logger_console, "cut-removed blocks: " << crack_blocks.countInt());

	// ... postprocess - avoid single edge/vertex connection
	int last_removed_count;
	do{
		last_removed_count = crack_blocks.countInt();
		// ... check edges (fix single edge connections)
		DataSimpleList<MeshEdge3d*> bad_edges;
		DataHashTable<MeshEdge3d*> hash_edges(last_removed_count * 8, nullptr);
		// ----> gather bad edges
		for(int i = 0; i < last_removed_count; i++){
			MeshBlock* block = crack_blocks[i];
			int bect = block->getEdgeCount();
			for(int j = 0; j < bect; j++){
				MeshEdge3d* edge = block->getEdge(j);
				if(!hash_edges.insert(edge)) continue;
				int efct = edge->getFaceCount();
				int half_bounded_count = 0;
				for(int k = 0; k < efct; k++){
					MeshFace* face = edge->getFaceAt(k);
					int non_crack_adjacent = 0;
					MeshBlock* b = face->getBlock(0);
					if(b && !b->availableTag(TagExtended::TAG_ADAPT_CRACK)) non_crack_adjacent++;
					b = face->getBlock(1);
					if(b && !b->availableTag(TagExtended::TAG_ADAPT_CRACK)) non_crack_adjacent++;

					if(non_crack_adjacent == 1) half_bounded_count++;
				}
				if(half_bounded_count > 2){
					//MeshViewSet* set = new MeshViewSet;
					//DataVector<MeshBlock*> eblocks(edge->getFaceCount());
					//edge->adjacentBlocks(eblocks);
					//for(int m = 0; m < eblocks.countInt(); m++){
					//	set->addBlockWithEdges(eblocks[m],
					//		crack_mark[eblocks[m]->getIndex()] ? 2 : 0);
					//}
					//set->addPoint(edge->getMeshPoint(0));
					//set->addPoint(edge->getMeshPoint(1));
					//SHOW_MESH("contact edges", set);
					bad_edges.append(edge);
				}
			}
		}
		if(bad_edges.notEmpty())
			LOG4CPLUS_WARN(MeshLog::logger_console, "createCrack - contact edges: " << bad_edges.countInt());
		// ----> mark _all_ blocks adjacent to bad edges
		while(bad_edges.notEmpty()){
			MeshEdge3d* edge = bad_edges.removeFirst();
			DataVector<MeshBlock*> blocks;
			if(edge->adjacentBlocks(blocks)){
				for(int i = 0; i < blocks.countInt(); i++){
					MeshBlock* block = blocks[i];
					if(!block->availableTag(TagExtended::TAG_ADAPT_CRACK)){
						crack_blocks.add(block);
						for(int j = 0; j < block->getPointCount(); j++)
							crack_mark_points[ block->getPoint(j)->getIndex() ];
						block->setIntTag(TagExtended::TAG_ADAPT_CRACK);
					}
				}
			}
		}
		// ... check nodes (remove stranded, fix single vertex connections)
		for(int i = 0; i < pct; i++){
			if(!crack_mark_points[i]) continue;
			// check for single vertex connection
			MeshPoint3d* point = mesh->getPointAt(i);
			DataVector<MeshBlock*> blocks(100);
			point->adjacentBlocks(blocks);
			DataVector<MeshBlock*> pblocks(blocks.countInt());
			for(int j = 0; j < blocks.countInt(); j++){
				if(blocks[j]->availableTag(TagExtended::TAG_ADAPT_CRACK)){
					pblocks.add(blocks[j]);
				}
			}
			int pbct = pblocks.countInt();
			assert(pbct > 0);
			DataVector<int> check_list(pbct);
			DataVector<bool> check_mark(pbct, false);
			check_list.add(0); // select first crack-block from the adjacency list
			check_mark[0] = true;
			for(int j = 0; check_list.countInt() < pbct && j < check_list.countInt(); j++){
				int k = check_list[j];
				MeshBlock* block = pblocks[k];
				int bfct = block->getFaceCount();
				for(int m = 0; m < bfct; m++){
					MeshBlock* mblock = block->getNeighbour(m);
					if(!mblock || !mblock->availableTag(TagExtended::TAG_ADAPT_CRACK)) continue;
					size_t mid = pblocks.find(mblock);
					if(mid != std::string::npos && mid > 0 && !check_mark[mid]){
						assert(mid <= INT_MAX);
						check_list.add((int)mid);
						check_mark[mid] = true;
					}
				}
			}
			if(check_list.countInt() != pbct){ // found single vertex connection
				LOG4CPLUS_WARN(MeshLog::logger_console, "createCrack - found single contact vertex");
				for(int j = 0; j < blocks.countInt(); j++){
					if(!blocks[j]->availableTag(TagExtended::TAG_ADAPT_CRACK)){
						MeshBlock* block = blocks[j];
						crack_blocks.add(block);
						for(int k = 0; k < block->getPointCount(); k++)
							crack_mark_points[ block->getPoint(k)->getIndex() ];
						block->setIntTag(TagExtended::TAG_ADAPT_CRACK);
					}
				}
			}
		}
	}while(last_removed_count != crack_blocks.countInt());

	return last_removed_count;
}

/// store binary file with 3d tetrahedral mesh
	//FileBinRax=Fopen(Filename,"rb");
	//fread(NbNoeud,sizeof(int),1,FileBinRax);
	//fread(NbEle,sizeof(int),1,FileBinRax);
	//fread(Coord,sizeof(float),3*(*NbNoeud),FileBinRax);
	//fread(Numele,sizeof(int),4*(*NbEle),FileBinRax);
	//fread(NbNoeSkin,sizeof(int),1,FileBinRax);
	//fclose(FileBinRax);
bool MeshSpecialRoutinesUTC::storeBinFile(const string& fname, MeshContainer3d* mesh)
{
	int pct = mesh->getPointsCount();
	int bct = mesh->getBlocksCount();
	if(pct == 0 || bct == 0) return false;

	ofstream f(fname, ostream::binary);
	if (!f) return false;

	f.write((const char*)&pct, sizeof(int));	// number of points
	f.write((const char*)&bct, sizeof(int));	// number of blocks (tetrahedra)
	for(int i = 0; i < pct; i++){
		const DPoint3d& pt = mesh->getPointAt(i)->getCoordinates();
		float coord[3] = { (float)pt.x, (float)pt.y, (float)pt.z };
		f.write((const char*)coord, 3*sizeof(float));
	}
	const int rct[] = {1, 0, 2, 3}; // for switching orientation of tetrahedra (sequence of vertices)
	for(int i = 0; i < bct; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		assert(block->getType() == BLOCK_TETRA);
		int indices[4] = { 
			block->getPoint(rct[0])->getIndex(),
			block->getPoint(rct[1])->getIndex(),
			block->getPoint(rct[2])->getIndex(),
			block->getPoint(rct[3])->getIndex() };
		f.write((const char*)indices, 4*sizeof(int));
	}
	int d = 0;
	f.write((const char*)&d, sizeof(int));
	f.close();

	return true;
}

/// load binary file with 3d tetrahedral mesh
	//FileBinRax=Fopen(Filename,"rb");
	//fread(NbNoeud,sizeof(int),1,FileBinRax);
	//fread(NbEle,sizeof(int),1,FileBinRax);
	//fread(Coord,sizeof(float),3*(*NbNoeud),FileBinRax);
	//fread(Numele,sizeof(int),4*(*NbEle),FileBinRax);
	//fread(NbNoeSkin,sizeof(int),1,FileBinRax);
	//fclose(FileBinRax);
MeshContainer3d* MeshSpecialRoutinesUTC::readBinFile(const string& fname)
{
	int pct, bct;
	ifstream f(fname, ifstream::binary);
	if( !f ) return nullptr;

	f.read((char*)&pct, sizeof(int));	// number of points
	f.read((char*)&bct, sizeof(int));	// number of blocks (tetrahedra)
	if(pct == 0 || bct == 0) return nullptr;

	// create mesh...
	MeshContainer3d* mesh = new MeshContainer3d(pct);
	// ... mesh points ...
	for(int i = 0; i < pct; i++){
		float coord[3];
		f.read((char*)coord, 3*sizeof(float));
		mesh->addMeshPoint(new MeshPoint3d((double)coord[0], (double)coord[1], (double)coord[2]));
	}
	const int rct[] = {1, 0, 2, 3}; // for switching orientation of tetrahedra (sequence of vertices)
	for(int i = 0; i < bct; i++){
		int indices[4];
		f.read((char*)indices, 4*sizeof(int));
		MeshBlock* block = new MeshTetrahedron(
			mesh->getPointAt(indices[rct[0]]),
			mesh->getPointAt(indices[rct[1]]),
			mesh->getPointAt(indices[rct[2]]),
			mesh->getPointAt(indices[rct[3]]));
		block->setAreaID(1);
		mesh->addMeshBlock(block);
	}
//	int d = 0;
//	fread(&d, sizeof(int), 1, f);
	f.close();

	mesh->setConstrainingPhase(MeshData::CONSTRAIN_DONE);
	mesh->setDiscretizationState(2);

	mesh->markOuterBoundary();
	mesh->markSharpBoundaryEdges();

	assert(mesh->isValid());
	return mesh;
}
