// MeshGenerator3dDelaunayBoundary.cpp: implementation of the MeshGenerator3d class.
//
//////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <algorithm>

#include "MeshGenerator3dDelaunayBoundary.h"

#include "MeshGenerator3d.h"
#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshDomainVolume.h"
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
#include "MeshGenerator2d.h"
#include "MeshGenerator3dQuality.h"
#include "MeshViewSet.h"
#include "SurfacePlane.h"
#include "IteratorEdge3d.h"
#include "IteratorFace.h"
#include "ControlSpace2dAdaptive.h"
#include "GeometricPredicates.h"
#include "DTetrahedron.h"
#include "DTriangle.h"
#include "DSegment.h"

int MeshGenerator3dDelaunayBoundary::param_boundary_nodes_for_recovery = 1;
double MeshGenerator3dDelaunayBoundary::param_initial_mesh_size_ratio = 0.5;

/// Volume threshold during 3D boundary recovery for checking line-face crossing etc.
double MeshGenerator3dDelaunayBoundary::param_recovery_volume_threshold = 1e-5;

/// Special metric ratio for boundary recovery 
//    (for avoiding edges in undesired directions)
double MeshGenerator3dDelaunayBoundary::param_special_metric_ratio = 1e-2;

double MeshGenerator3dDelaunayBoundary::param_refined_ICD_ratio = 0.15;

//#define LOG_RECOVERY_PIPES
#ifdef _DEBUG
//#define CONSTRAINING_MESH_VALIDITY_CHECK
#endif
//#define DEBUG_RECOVERY_PIPES

#define SPECIAL_METRIC_SWAP_LAYERS 5


/// Generate boundary constrained tetrahedral mesh
MeshContainer3d* MeshGenerator3dDelaunayBoundary::createBoundaryConstrainedMesh(
	Metric3dContext& mc, MeshContainer3dSurface* surface_mesh, MeshDomainVolume* mdv)
{
	LOG4CPLUS_INFO(MeshLog::logger_console, "Gen3d. Triangulating boundary (Delaunay)...");
	START_CLOCK("MG3d::triangulateBoundary");

	auto freepoints = mdv->getFreePoints();
	int fcount = freepoints ? freepoints->countInt() : 0;
	int bcount = surface_mesh->getPointsCount();

	// bounding box
	DBox bbox = mdv->getBoundingBox();

	LOG4CPLUS_INFO(MeshLog::logger_mesh, " initial-mesh");
	// initial mesh
	MeshContainer3d* mesh = createInitialMesh(bcount, bbox, mdv->getControlSpace());
	if(!mesh) return nullptr;

	LOG4CPLUS_INFO(MeshLog::logger_mesh, " initial-mesh done");

	srand(0);

	bool refine_mesh = (MeshGenerator3dDelaunayBoundary::param_refined_ICD_ratio > 0.0);
	int refined_count = 0;

	OctTreeMeshPoints octree_points(bbox);
	if(refine_mesh){
		for(int i = 0; i < bcount; i++)
			octree_points.insertMeshPointLink(surface_mesh->getPointAt(i));
		for(int i = 0; i < fcount; i++)
			octree_points.insertMeshPointLink(freepoints->get(i).get());
	}

	// Total number of points = boundary + free
	int points_count = bcount + fcount;

	DataVector<int> sequence(points_count);
	for(int i = 0; i < points_count; i++)
		sequence.add(i);
	std::random_shuffle( sequence.begin(), sequence.end() );

	// * insert new points into the mesh
	for(int i = 0; i < points_count; i++){
		// Create new point as a clone
		MeshPoint3d *boundary_point = 
			sequence[i] < bcount ? surface_mesh->getPointAt(sequence[i]) : freepoints->get(sequence[i] - bcount).get();
		MeshPoint3d *point = new MeshPoint3d(*boundary_point);
		//point->setID(i);

		// Retriangulate
		mc.countMetricAtPoint(point->getCoordinates());
		if(MeshGenerator3d::addPointToTriangulation(mc, mesh, point)){
			// mark point as boundary
			point->setPtrTag(TagExtended::TAG_BOUNDARY_POINT, boundary_point);
			point->setBorder();
			void * mbc = boundary_point->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
			if(mbc) point->setPtrTag(TagExtended::TAG_BOUNDARY_COND, mbc);
			if(refine_mesh){
				octree_points.removeMeshPointLink(boundary_point);
				// check tetrahedra around inserted point and insert additional "refinement" nodes if necessary
				// judge with "geometric quality in metric" + not-too-close
				refined_count += refineMeshForBoundaryTriangulation(mc, mesh, point, octree_points);
			}
		}else{
//			hard_points.add(point);
			delete point;
			delete mesh; 
			return nullptr;
		}
	}
	STOP_CLOCK("MG3d::triangulateBoundary");

	assert(mesh->isValid());

	if(refine_mesh){
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "NP-boundary\t" << bcount);
		if(fcount > 0)	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "NP-free\t" << fcount);
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "NP-refined\t" << refined_count);
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "Gen3d. Constraining...");
	START_CLOCK("MG3d::constrainToBoundary");
	bool success = constrainToBorder(mc, mesh, mdv, surface_mesh);
	mesh->clearSearchTree();
	if(success) MeshGenerator3dQuality::smoothenSwapComplete(mc, mesh);
	STOP_CLOCK("MG3d::constrainToBoundary");

	if(success){
		return mesh;
	}else{
		delete mesh; 
		return nullptr;
	}
}

MeshContainer3d* MeshGenerator3dDelaunayBoundary::createInitialMesh(int points_count, 
						const DBox& bbox, CS3dPtr cs)
{
	if(points_count < 1) return nullptr;

	MeshContainer3d* mesh = new MeshContainer3d(points_count);
	mesh->setControlSpace(cs);

	mesh->setConstrainingPhase(MeshData::CONSTRAIN_NONE);

	// Initial mesh: 5 tetrahedra (making cube)
	DBox box = bbox;
	double offset = param_initial_mesh_size_ratio * box.getDiameter();
	box.grow(offset);

	MeshPoint3d * vertices[8];
	mesh->addMeshPoint(vertices[0] = new MeshPoint3d( box.x0,	box.y0,	box.z0));
	mesh->addMeshPoint(vertices[1] = new MeshPoint3d( box.x1,	box.y0,	box.z0));
	mesh->addMeshPoint(vertices[2] = new MeshPoint3d( box.x1,	box.y0,	box.z1));
	mesh->addMeshPoint(vertices[3] = new MeshPoint3d( box.x0,	box.y0,	box.z1));
	mesh->addMeshPoint(vertices[4] = new MeshPoint3d( box.x0,	box.y1,	box.z0));
	mesh->addMeshPoint(vertices[5] = new MeshPoint3d( box.x1,	box.y1,	box.z0));
	mesh->addMeshPoint(vertices[6] = new MeshPoint3d( box.x1,	box.y1,	box.z1));
	mesh->addMeshPoint(vertices[7] = new MeshPoint3d( box.x0,	box.y1,	box.z1));

	for(int i = 0; i < 8; i++) {
		vertices[i]->setBorder();
		vertices[i]->setIntTag(TagExtended::TAG_OUTERHULL_POINT);
	}

	// Drzewo wyszukiwania czworoœcianów w pobli¿u wstawianego punktu
	if(OctTree::param_octtree_threshold){
		mesh->setSearchTree(std::make_shared<OctTree>(box));
	}

	// Five initial tetrahedra 
	//	(with automatic creation of neighbourhood-connections and incidence tables)
	int conf[5][4] = {{0, 1, 2, 5}, {0, 2, 3, 7}, {0, 5, 7, 4}, {2, 5, 6, 7}, {0, 5, 2, 7}};

	for(int i = 0; i < 5; i++)
		mesh->addMeshTetrahedron(new MeshTetrahedron(
				vertices[conf[i][0]],vertices[conf[i][1]],
				vertices[conf[i][2]],vertices[conf[i][3]]));

	// mark boundary faces (and edges)
	vertices[0]->getFaceToPoints(vertices[1], vertices[2])->setBorderWithEdges(); //--	
	vertices[0]->getFaceToPoints(vertices[2], vertices[3])->setBorderWithEdges();
	vertices[0]->getFaceToPoints(vertices[1], vertices[5])->setBorderWithEdges(); //--
	vertices[0]->getFaceToPoints(vertices[5], vertices[4])->setBorderWithEdges();
	vertices[2]->getFaceToPoints(vertices[5], vertices[1])->setBorderWithEdges(); //--
	vertices[2]->getFaceToPoints(vertices[6], vertices[5])->setBorderWithEdges();
	vertices[2]->getFaceToPoints(vertices[7], vertices[3])->setBorderWithEdges(); //--
	vertices[2]->getFaceToPoints(vertices[6], vertices[7])->setBorderWithEdges();
	vertices[0]->getFaceToPoints(vertices[7], vertices[4])->setBorderWithEdges(); //--
	vertices[0]->getFaceToPoints(vertices[3], vertices[7])->setBorderWithEdges();
	vertices[5]->getFaceToPoints(vertices[4], vertices[7])->setBorderWithEdges(); //--
	vertices[5]->getFaceToPoints(vertices[7], vertices[6])->setBorderWithEdges();

	return mesh;
}

bool MeshGenerator3dDelaunayBoundary::constrainToBorder(Metric3dContext& mc, 
		MeshContainer3d* mesh, const MeshDomainVolume* mdv, 
		MeshContainer3dSurface* surface_mesh)
{
//	SHOW_STEP_BREAKABLE(1, "Constraining of the mesh and the recovery of the boundary.", false);
	const int phase_count = param_boundary_nodes_for_recovery ? 5 : 4;

	int pct = mesh->getPointsCount();
	assert(pct > 3);

	mesh->setConstrainingPhase(MeshData::CONSTRAIN_ACTIVE);

	DataVector<MissingEdge> missing_edges(pct);
	DataVector<MissingFace> missing_faces(pct);

	int stat_edges[5][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
	int stat_missing_edges[5] = {0, 0, 0, 0, 0};
	int stat_faces[5][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
	int stat_missing_faces[5] = {0, 0, 0, 0, 0};

	int spct = surface_mesh->getPointsCount();
	int sfct = surface_mesh->getFacesCount();

	// back-trace hash array for points
	DataHashTableKeyValue<MeshPoint3d*,MeshPoint3d*> bmpoints(2*spct, nullptr);
	for(int i = 0; i < pct; i++){
		MeshPoint3d* mp = mesh->getPointAt(i);
		MeshPoint3d* bp = (MeshPoint3d*)mp->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
		if(bp){
		   	bmpoints.insert(bp, mp);
			mp->setBorder();
		}
	}
	// check edges only once
	DataHashTable<MeshEdge3d*> unique_edges(5*sfct, nullptr);

	// Check whether all boundary edges and faces are correctly represented
	int mdv_id = mdv->getAreaID();

	for(int i = 0; i < sfct; i++){
		MeshFace* b_face = surface_mesh->getFaceAt(i);
		MeshPoint3d* b_pts[3] = {
			b_face->getPoint(0), b_face->getPoint(1), b_face->getPoint(2) 
		};
		MeshPoint3d* m_pts[3] = {
			bmpoints.getValue(b_pts[0], nullptr),
			bmpoints.getValue(b_pts[1], nullptr),
			bmpoints.getValue(b_pts[2], nullptr)
		};
		assert(m_pts[0] && m_pts[1] && m_pts[2]); // have to be there
		// check face
		MeshFace* m_face = m_pts[0]->getFaceToPoints(m_pts[1], m_pts[2]);

		if(m_face){
			m_face->setBorder();
			m_face->setPtrTag(TagExtended::TAG_BOUNDARY_FACE, b_face);
			// Ascribe proper "material id" for blocks on both side of the face
			//  value -1 denotes block to remove ...
			assert(m_face->getBlock(0)!=nullptr && m_face->getBlock(1)!=nullptr);
			bool oriented = m_face->properOrientation( m_pts[0], m_pts[1]);
			MeshBlock* block = b_face->getBlock(0);
			m_face->getBlock(oriented ? 0 : 1)->setAreaID((block == mdv) ? mdv_id : -1);
			block = b_face->getBlock(1);
			m_face->getBlock(oriented ? 1 : 0)->setAreaID((block == mdv) ? mdv_id : -1);

			MeshBoundaryCondition* boundary_condition = 
				(MeshBoundaryCondition*)b_face->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
			if(boundary_condition)
				m_face->setPtrTag(TagExtended::TAG_BOUNDARY_COND, boundary_condition);
		}else{
			missing_faces.add(MissingFace(b_face, m_pts[0], m_pts[1], m_pts[2]));
		}
		// check edges
		for(int j = 0; j < 3; j++){
			MeshEdge3d* b_edge = b_face->getEdge(j);
			if(!unique_edges.insert(b_edge)) continue;
			MeshEdge3d* m_edge = m_pts[j]->getEdgeToPoint(m_pts[(j+1)%3]);
			if(m_edge){
				m_edge->copyBorderFlagsFrom(b_edge);
				m_edge->setBorderFlags(TagBorder::OUTER); // additionally
				MeshBoundaryCondition* boundary_condition =
					(MeshBoundaryCondition*)b_edge->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
				if(boundary_condition)
					m_edge->setPtrTag(TagExtended::TAG_BOUNDARY_COND, boundary_condition);
			}else{
				missing_edges.add(MissingEdge(b_edge, m_pts[j], m_pts[(j+1)%3]));
			}
		}
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "-----------------------");
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
		"Missing edges = " << missing_edges.countInt() << "   total = " << unique_edges.countInt());
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
		"Missing faces = " << missing_faces.countInt() << "   total = " << sfct);
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Total nodes count = " << pct);

	DataVector<MeshContainer3d::AdditionalBoundaryNode>& inserted_boundary_points = 
		mesh->getAdditionalBoundaryNodes();
	inserted_boundary_points.clear();

	const char* labels[] = { "phase-0", "phase-1", "phase-2", 
		"phase-3", "phase-4", "phase-5", "phase-6" };

//	MeshViewSet::showDebugMesh(boundary, mesh, 0);

	for(int phase = 0; phase < phase_count && 
		(!missing_edges.empty() || !missing_faces.empty()); phase++)
	{

		if (false) {
			// view - all edges (missing marked) + missing faces
			MeshViewSet* set = new MeshViewSet;
			for (auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()) {
				MeshEdge3d* edge = it.getEdge();
				if (edge->isBorder()) set->addEdge(edge, -1);
			}
			missing_edges.forEach([set](const MissingEdge & me) {
				set->addEdge(me.edge, 2);
			});
			missing_faces.forEach([set](const MissingFace & mf) {
				set->addFace(mf.face);
			});
			set->addInfo("phase", to_string(phase));
			SHOW_MESH("missing edges & faces", set);
		}

		// view - missing faces cavities?

		START_CLOCK(labels[phase]);

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Starting recovery - Mesh invalid, leaving. Phase:", phase); return false;
		}
#endif
//		MeshView::showDebugMesh(boundary, mesh, phase);
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Recovery phase = " << phase);

		stat_missing_edges[phase] = (int)missing_edges.countInt();
		stat_missing_faces[phase] = (int)missing_faces.countInt();
		bool any_change = true;
		while(any_change){
			LOG4CPLUS_DEBUG(MeshLog::logger_console, "Missing edges count=" << missing_edges.countInt());
			any_change = false;
			// for phase==3 and inserting boundary edges/faces
			DataVector<MissingEdge> new_missing_edges;
			//	Restore missing edges
			for(int i = 0; i < missing_edges.countInt(); i++){
				MissingEdge& me = missing_edges[i];
				// --- try to recover
				//mc.countMetricAtPoint(m_pt0->getCoordinates());
				int result = recoverBoundaryEdge(mc, mesh, me, missing_faces, phase);

				if(!result && !me.edge->availableTag(TagExtended::TAG_FIXED)){
					// check if the boundary set of faces+edges cannot be adapted ...
					result = matchBoundaryEdge(mc, me, bmpoints);
				}

				//if(!result && (phase == 2 || phase == 3)){
				//	mc.setSpecialMetric(me.getSpecialMetric());
				//	result = recoverBoundaryEdgeWithSpecialMetricSwap(mc, mesh, m_pt0, m_pt1,
				//		(phase==2) ? SPECIAL_METRIC_SWAP_LAYERS : 2*SPECIAL_METRIC_SWAP_LAYERS);
				//	LOG4CPLUS_DEBUG(MeshLog::logger_console, "Special metric swap for edge", result);
				//	if(!result){ 
				//		//mc.countMetricAtPoint(m_pt0->getCoordinates());
				//		result = recoverBoundaryEdge(mc, mesh, m_pt0, m_pt1, missing_faces_data, phase);
				//	}
				//}

				if(!result && phase == 4){ // insert point in the middle of missing edge
					// insert
//					mc.setSpecialMetric(me.getSpecialMetric());
//					result = recoverBoundaryEdgeWithSpecialMetricSwap(mc, mesh, m_pt0, m_pt1);
//					if(!result){
					MeshGenerator3d::iterativeTetrahedraSwap(mc, mesh, 
								me.m_pt0, SPECIAL_METRIC_SWAP_LAYERS);
/*
					if(true){
						MeshViewSet *set = new MeshViewSet;
						int lct = mesh->getBlocksCount();
						for(int li = 0; li < lct; li++){
							MeshBlock* block = mesh->getBlockAt(li);
							int fct = block->getFaceCount();
							for(int lj = 0; lj < fct; lj++){
								MeshFace* face = block->getFace(lj);
								if(face->isBorder() && face->getBlockIndex(block) == 0 &&
									face->getPoint(0)->getID() > 0 &&
									face->getPoint(1)->getID() > 0 &&
									face->getPoint(2)->getID() > 0)
								{
									if(block->getAreaID() >= 0)
										set->addFace(face, 0);
									else{
										set->addFace(face->getPoint(1)->getCoordinates(),
											face->getPoint(0)->getCoordinates(),
											face->getPoint(2)->getCoordinates(), 0);
									}
								}
							}
						}
//						for(int li = 0; li < missing_faces.countInt(); li++){
//							set->addFace(missing_faces[li].face, 1);
//						}
						set->addPoint(m_pt0, 2);
						set->addPoint(m_pt1, 2);
						SHOW_MESH("inserting inner node", set);
					}
*/
					insertBoundaryPointForEdgeRecovery(mc, mesh, me, missing_faces,
						new_missing_edges, inserted_boundary_points, phase);
					// remove old edge
					missing_edges.removeAt(i--);
					any_change = true;
//					}
				}
				if(result){
					assert(result < 5);
					//if(phase > 1) any_change = true;
					++stat_edges[phase][result];
					MeshEdge3d* m_edge = me.m_pt0->getEdgeToPoint(me.m_pt1);
					assert(m_edge);
					if (me.edge) m_edge->copyBorderFlagsFrom(me.edge);
					m_edge->setBorderFlags(TagBorder::OUTER);

					MeshBoundaryCondition* bc = nullptr;
					if(me.edge) bc = (MeshBoundaryCondition*)me.edge->getPtrTag(
										TagExtended::TAG_BOUNDARY_COND);
					if(bc) m_edge->setPtrTag(TagExtended::TAG_BOUNDARY_COND, bc);
					// remove
					missing_edges.removeAt(i--);	// decreases count and index
				}
			}
			if(new_missing_edges.countInt() > 0){
				for(int i = 0; i < new_missing_edges.countInt(); i++)
					missing_edges.add(new_missing_edges[i]);
				new_missing_edges.clear();
			}
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "----------- Edges ------------");
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Missing edges = " << missing_edges.countInt() 
				<< "   total = " << total_edges << endl;
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Total nodes count = " << mesh->getPointsCount());
#endif
			assert(mesh->isValid());

			LOG4CPLUS_DEBUG(MeshLog::logger_console, "Missing faces count= " << missing_faces.countInt());
			//	Restore missing faces
			DataVector<MissingFace> new_missing_faces;
			for(int i = 0; i < missing_faces.countInt(); i++){
				MissingFace mf = missing_faces[i];
				// --- try to recover
				mc.countMetricAtPoint(mf.m_pt[0]->getCoordinates());
				int result = recoverBoundaryFace(mc, mesh, missing_faces, 
								mf.m_pt[0], mf.m_pt[1], mf.m_pt[2], phase);

				//if(!result && phase == 3){
				//	mc.setSpecialMetric(mf.getSpecialMetric());
				//	result = recoverBoundaryFaceWithSpecialMetricSwap(mc, mesh, 
				//				m_pt1, m_pt2, m_pt3,
				//		SPECIAL_METRIC_SWAP_LAYERS);
				//	LOG4CPLUS_DEBUG(MeshLog::logger_console, "Special metric swap for edge", result);
				//}

				if(!result && phase == 4 && 
					mf.m_pt[0]->getEdgeToPoint(mf.m_pt[1]) &&
					mf.m_pt[1]->getEdgeToPoint(mf.m_pt[2]) && 
					mf.m_pt[2]->getEdgeToPoint(mf.m_pt[0]))
				{ 
					MeshGenerator3d::iterativeTetrahedraSwap(mc, mesh, 
							mf.m_pt[0], SPECIAL_METRIC_SWAP_LAYERS);
					// insert point in the middle of missing face (if all edges are ok)
					insertBoundaryPointForFaceRecovery(mc, mesh, mf, missing_edges, 
						missing_faces, new_missing_faces, inserted_boundary_points, phase);
					missing_faces.removeAt(i--);
					any_change = true;
				}
				if(result){
					assert(result < 4);
					//if(phase > 1) any_change = true;
					++stat_faces[phase][result];
					// update
					MeshFace* m_face = mf.m_pt[0]->getFaceToPoints(mf.m_pt[1], mf.m_pt[2]);
					assert(m_face != nullptr);
					assert(m_face->isBorder() == false);
					m_face->setBorderWithEdges();
					for (int i = 0; i < 3; i++) {
						MeshEdge3d* m_edge = m_face->getEdge(i);
						MeshPoint3d* bp0 = (MeshPoint3d*)(m_edge->getMeshPoint(0))->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
						MeshPoint3d* bp1 = (MeshPoint3d*)(m_edge->getMeshPoint(1))->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
						//assert(bp0 && bp1);
						if (bp0 && bp1) {
							MeshEdge3d* b_edge = bp0->getEdgeToPoint(bp1);
							if (b_edge) m_edge->copyBorderFlagsFrom(b_edge);
						}
						m_edge->setBorderFlags(TagBorder::OUTER);
					}
					m_face->setPtrTag(TagExtended::TAG_BOUNDARY_FACE, mf.face);
					// Ascribe proper "material id" for blocks on both side of the face
					//  value -1 denotes block to remove ...
					assert(m_face->getBlock(0)!=nullptr && m_face->getBlock(1)!=nullptr);
					bool oriented = m_face->properOrientation(mf.m_pt[0], mf.m_pt[1]);
					MeshBlock* block = mf.face->getBlock(0);
					m_face->getBlock(oriented ? 0 : 1)->setAreaID((block==mdv) ? mdv_id : -1);
					block = mf.face->getBlock(1);
					m_face->getBlock(oriented ? 1 : 0)->setAreaID((block==mdv) ? mdv_id : -1);

					MeshBoundaryCondition* boundary_condition = 
						(MeshBoundaryCondition*)mf.face->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
					if(boundary_condition){
						m_face->setPtrTag(TagExtended::TAG_BOUNDARY_COND, boundary_condition);
					}
					// remove
					missing_faces.removeAt(i--);
				}
			}
			if(new_missing_faces.countInt() > 0){
				for(int i = 0; i < new_missing_faces.countInt(); i++)
					missing_faces.add(new_missing_faces[i]);
				new_missing_faces.clear();
			}
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "----------- Faces ------------");
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Missing faces = " << missing_faces.countInt() 
				<< "   total = " << total_faces << endl;
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Total nodes count = " << mesh->getPointsCount());
#endif
		}

		STOP_CLOCK(labels[phase]);
	}

	// Set area-IDs for all tetrahedra ...

	bool success = (missing_edges.countInt() + missing_faces.countInt() == 0);

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Finished recovery - Mesh invalid, leaving. Success:", success); return false;
		}
#endif

	if(success){
		int mesh_count = mesh->getBlocksCount();
		// Tetrahedra marked as valid
		DataVector<MeshTetrahedron*> approved(mesh_count);

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "----------- Cleaning ------------");
		ostringstream log_file;
		log_file << mesh_count << " blocks, ";

		for(int i = 0; i < mesh_count; i++){
			MeshTetrahedron* tetrahedron  = (MeshTetrahedron*)mesh->getBlockAt(i);
			if(tetrahedron->getAreaID() > -1){	// Does contain material ?
				tetrahedron->setIntTag(TagExtended::TAG_MG3D);
				approved.add(tetrahedron);
			}
		}
/*
		MeshViewSet *set = new MeshViewSet(0, 0, 0, approved_ct);
		for(int i = 0; i < approved_ct; i++){
			set->addBlockWithEdges(approved[i]);
		}
		SHOW_MESH("marked blocks", set);
*/
		log_file << approved.countInt() << " at boundary with material, ";

		// For each marked tetrahedron add its neighbours ...
		for(int i = 0; i < approved.countInt(); i++){
			for(int j = 0; j < 4; j++){
				MeshTetrahedron* tetrahedron = (MeshTetrahedron*) approved[i]->getNeighbour(j);
				if(tetrahedron && 
					!tetrahedron->availableTag(TagExtended::TAG_MG3D) && 
					!approved[i]->isBorder(j))
				{
					tetrahedron->setIntTag(TagExtended::TAG_MG3D);
					approved.add(tetrahedron);
					tetrahedron->setAreaID(approved[i]->getAreaID());
				}
			}
		}

		log_file << approved.countInt() << " valid (total).";
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, log_file.str());

//		MeshViewSet* set = new MeshViewSet(0, 6* mesh->getBlocksCount(), 0, mesh->getBlocksCount());
		// Remove unmarked tetrahedra
		for(int i = 0; i < mesh->getBlocksCount(); ){
			MeshTetrahedron *tetrahedron = (MeshTetrahedron*)mesh->getBlockAt(i);
			if(tetrahedron->availableTag(TagExtended::TAG_MG3D)){
				tetrahedron->removeTag(TagExtended::TAG_MG3D);
				++i;
//				set->addBlockWithEdges(tetrahedron);
			}else{
	//			SHOW_STEP_PT(3, "Usuwany niepotrzebny trójk¹t.", tetrahedron->getMiddlePoint()); 
				delete mesh->removeMeshTetrahedron(tetrahedron);
				// While removing an element, the total number of blocks decreases and instead the
				//	removed element some other is inserted (unless it was already the last one)
				//  - so we DON'T decrease i.
			}
		}

//		SHOW_MESH("boundary constrained mesh", set);

		// remove auxiliary points
		int nodes_unconnected = 0;
		int nodes_negative_ID = 0;
		int non_boundary = 0;
		pct = mesh->getPointsCount();
		for(int i = 0; i < pct; ){	
			MeshPoint3d* point = mesh->getPointAt(i);
			if(point->getRank() < 1){
				delete mesh->removeMeshPoint(i);
				++nodes_unconnected;
				--pct;
			}else if(point->availableTag(TagExtended::TAG_OUTERHULL_POINT)){
				// remove boundary edges of the initial triangualation, which are obsolete now ...
				while(point->getRank() > 0){
					MeshEdge3d* edge = point->getEdge(0);
					edge->clearBorder();
					int j = edge->getFaceCount();
					if(j == 0) delete edge;
					else while(j--) delete edge->getFaceAt(0);
						// edge should be deleted automatically, by last removed face
				}
				delete mesh->removeMeshPoint(i);
				++nodes_negative_ID;
				--pct;
			}else{
				if(!point->isBorder()){
					++non_boundary;
					MeshGenerator3dQuality::movePointByLaplace(mc, point);
				}
				++i;
			}
		}

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Removed obsolete entities - Mesh invalid, leaving."); return false;
		}
#endif
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Removed " << nodes_unconnected << " unconnected nodes.");
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Removed " << nodes_negative_ID << " initial-mesh nodes.");

		if(inserted_boundary_points.countInt() > 0){
			fixAdditionalBoundaryNodes(mc, mesh, inserted_boundary_points);
			LOG4CPLUS_DEBUG(MeshLog::logger_console, 
				"Additional boundary nodes (final): " << inserted_boundary_points.countInt());
		}

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Non-boundary nodes: " << non_boundary);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Additional boundary nodes: " << inserted_boundary_points.countInt());

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
			"Constraining done, " << mesh->getBlocksCount() << " blocks and "
			<< mesh->getPointsCount() << " points left. ");
	}

	mesh->setConstrainingPhase(MeshData::CONSTRAIN_DONE);
	
	// show stats
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "|---------------------------------------------------------|");
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "|   |  E  |  F  | e-0 | e-t | e-p0| e-pt| f-0 | f-t | f-pt|");
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
		"| - |" 
		<< setw(4) << unique_edges.countInt() << " |"  
		<< setw(4) << sfct << " |");
	for(int i = 0; i < phase_count; i++){
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
			"| " << i << " |"
			<< setw(4) << stat_missing_edges[i] << " |"
			<< setw(4) << stat_missing_faces[i] << " |");

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh_stat, "Missing-edges-" << labels[i] << "\t" << stat_missing_edges[i]);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh_stat, "Missing-faces-" << labels[i] << "\t" << stat_missing_faces[i]);

		ostringstream log_line;
		for(int j = 1; j <= 4; j++){
			log_line << setw(4) << stat_edges[i][j] << " |";
		}
		for(int j = 1; j <= 3; j++){
			log_line << setw(4) << stat_faces[i][j] << " |";
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, log_line.str());
	}
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
		"| = |" 
		<< setw(4) << missing_edges.countInt() << " |"
		<< setw(4) << missing_faces.countInt() << " |");

//	SHOW_MESH("Volume mesh after fix", mesh->getViewSet());
	//MeshViewSet* set = new MeshViewSet(0, 0, 0, mesh->getBlocksCount());
	//for(int i = 0; i < mesh->getBlocksCount(); i++)
	//	set->addBlock(mesh->getBlockAt(i));
	//SHOW_MESH("Volume mesh", set);

	return success;
}

bool MeshGenerator3dDelaunayBoundary::fixAdditionalBoundaryNodes(Metric3dContext& mc, MeshContainer3d* mesh,
		DataVector<MeshContainer3d::AdditionalBoundaryNode> & inserted_boundary_points)
{
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "fixAdditionalBoundaryNodes - start");
	int i = inserted_boundary_points.countInt();
	while(--i >= 0){
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "fix extra b-nodes start, count= " << i);
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "mesh " << (mesh->isValid() ? "valid" : "invalid"));
		MeshContainer3d::AdditionalBoundaryNode& bn = inserted_boundary_points[i];
//		restoreFromSpecialMetricSwap(mc, mesh, bn.point);
		mc.countMetricAtPoint(bn.point->getCoordinates());
		if(bn.other_points.countInt() == 3){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "fix-3 -> TO BE DONE CORRECTLY!");
			// recover "face-case"
			MeshPoint3d* points[] = {bn.other_points[0], bn.other_points[1], bn.other_points[2]};
			MeshFace* faces[] = {bn.point->getFaceToPoints(points[0], points[1]),
				bn.point->getFaceToPoints(points[1], points[2]),
				bn.point->getFaceToPoints(points[2], points[0])};
			if(!faces[0] || !faces[1] || !faces[2]){
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"fixAdditionalBoundaryNodes: can't fix face - subfaces?");
				continue;
			}
			MeshTetrahedron* tetrahedra[] = {(MeshTetrahedron*)faces[0]->getBlock(0),
				(MeshTetrahedron*)faces[0]->getBlock(1)};
			assert(tetrahedra[0] || tetrahedra[1]);
			//count min len
			const DMPoint3d dpt = bn.point->getMetricCoordinates(mc);
			double min_len2 = dpt.distance2(
						bn.point->getEdge(0)->getOtherPoint(bn.point)->getMetricCoordinates(mc));
			int rank = bn.point->getRank();
			for(int j = 1; j < rank; j++){
				double len2 = dpt.distance2(
						bn.point->getEdge(j)->getOtherPoint(bn.point)->getMetricCoordinates(mc));
				if(len2 < min_len2) min_len2 = len2;
			}
			const DMVector3d dnv = DMVector3d::crossProduct(
				points[0]->getMetricCoordinates(mc), 
				points[1]->getMetricCoordinates(mc), 
				points[2]->getMetricCoordinates(mc)).normalized() * (sqrt(min_len2) * 0.5);
			const DVector3d nv = mc.transformMStoRS(dnv);
			DPoint3d moved_pt = bn.point->getCoordinates() + nv;
			double vol0 = DTriangle3d::orient3d(
				points[0]->getCoordinates(), 
				points[1]->getCoordinates(), 
				points[2]->getCoordinates(), 
				moved_pt);
			if(tetrahedra[0] && tetrahedra[1]){ // two-side				
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "fix-3 two-side -> TO BE DONE CORRECTLY!");
				// orientation check
				double vol1 = DTriangle3d::orient3d(
					points[0]->getCoordinates(), 
					points[1]->getCoordinates(), 
					points[2]->getCoordinates(), 
					tetrahedra[0]->getOppositePoint(faces[0])->getCoordinates());
				if(vol0 * vol1 < 0.0){ // different sides
					MeshTetrahedron* t = tetrahedra[0];
					tetrahedra[0] = tetrahedra[1];
					tetrahedra[1] = t;
				}
				MeshPoint3d* new_point = new MeshPoint3d(*bn.point);
				DataVector<MeshBlock*> blocks;
				bn.point->adjacentBlocks(blocks);
				for(int j = 0; j < blocks.countInt(); j++){
					MeshBlock* block = blocks[j];
					vol1 = DTriangle3d::orient3d(
						points[0]->getCoordinates(), 
						points[1]->getCoordinates(), 
						points[2]->getCoordinates(), 
						block->getMiddlePoint());
					if(vol0 * vol1 < 0.0) // if different sides -> switch block to new point
						block->switchPointsWithFaces(bn.point, new_point);
				}
				MeshGenerator3dQuality::tryMovingPoint(mc, bn.point, moved_pt, false);
				mesh->addMeshPoint(new_point);
				moved_pt = bn.point->getCoordinates() - nv;
				MeshGenerator3dQuality::tryMovingPoint(mc, new_point, moved_pt, false);
				bn.point->clearBorder();
				for(int j = 0; j < 3; j++){
					assert(bn.point->getEdgeToPoint(points[j]));
					bn.point->getEdgeToPoint(points[j])->clearBorder();
					faces[j]->clearBorder();
				}
				MeshTetrahedron* new_tetrahedra[2] = {
					new MeshTetrahedron(points[0], points[1], points[2], bn.point),
					new MeshTetrahedron(points[2], points[1], points[0], new_point)};
				for(int j = 0; j < 2; j++){
					new_tetrahedra[j]->setAreaID(tetrahedra[j]->getAreaID());
					mesh->addMeshBlock(new_tetrahedra[j]);
				}
				MeshFace* new_face = new_tetrahedra[0]->getOppositeFace(bn.point);
				new_face->setBorder(TagBorder::INNER); // inner boundary
				MeshGenerator3dQuality::movePointByLaplace(mc, bn.point);
				MeshGenerator3dQuality::movePointByLaplace(mc, new_point);
//				LOG4CPLUS_INFO(MeshLog::logger_console, "mesh valid **", mesh->isValid());
//				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "--> additional boundary (face-inner) node " << 
//						bn.point->getIndex() << " fixed." << endl;
				inserted_boundary_points.removeAt(i);
			}else{
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "fix-3 one-side -> TO BE DONE CORRECTLY!");
				// one-side
				MeshTetrahedron* tetrahedron = tetrahedra[0] ? tetrahedra[0] : tetrahedra[1];
				double vol1 = DTriangle3d::orient3d(
					points[0]->getCoordinates(), 
					points[1]->getCoordinates(), 
					points[2]->getCoordinates(), 
					tetrahedron->getOppositePoint(faces[0])->getCoordinates());
				if(vol0 * vol1 < 0.0) // switch orientation
					moved_pt = bn.point->getCoordinates() - nv;
				MeshBoundaryCondition * condition = 
					(MeshBoundaryCondition*)faces[0]->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
				void * bface_ptr = faces[0]->getPtrTag(TagExtended::TAG_BOUNDARY_FACE);
				if(MeshGenerator3dQuality::tryMovingPoint(mc, bn.point, moved_pt, false)){
					bn.point->clearBorder();
					for(int j = 0; j < 3; j++){
						assert(bn.point->getEdgeToPoint(points[j]));
						bn.point->getEdgeToPoint(points[j])->clearBorder();
						faces[j]->clearBorder();
					}
					MeshTetrahedron* new_tetrahedron = new MeshTetrahedron(
						points[0], points[1], points[2], bn.point);
					assert(new_tetrahedron);
					new_tetrahedron->setAreaID(tetrahedron->getAreaID());
					mesh->addMeshBlock(new_tetrahedron);
					MeshFace* new_face = new_tetrahedron->getOppositeFace(bn.point);
					new_face->setBorderWithEdges();
					if(condition)
						new_face->setPtrTag(TagExtended::TAG_BOUNDARY_COND, condition);
					if(bface_ptr)
						new_face->setPtrTag(TagExtended::TAG_BOUNDARY_FACE, bface_ptr);
					MeshGenerator3dQuality::movePointByLaplace(mc, bn.point);
					LOG4CPLUS_DEBUG(MeshLog::logger_console, "one-side-3 mesh valid *** " << mesh->isValid());
//					LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "--> additional boundary (face) node " << 
//						bn.point->getIndex() << " fixed." << endl;
					inserted_boundary_points.removeAt(i);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "fixAdditionalBoundaryNodes: failed moving point");
				}
			}
		}else if(bn.other_points.countInt() == 4){
			// edge [0]-[1], [2] and [3] -> adjacent faces
			MeshPoint3d* points[] = {bn.other_points[0], bn.other_points[2], 
				bn.other_points[1], bn.other_points[3]};
			MeshFace* faces[] = {bn.point->getFaceToPoints(points[0], points[1]),
				bn.point->getFaceToPoints(points[1], points[2]),
				bn.point->getFaceToPoints(points[2], points[3]),
				bn.point->getFaceToPoints(points[3], points[0])};
			if(!faces[0] || !faces[1] || !faces[2] || !faces[3]){
				LOG4CPLUS_WARN(MeshLog::logger_console, "fixAdditionalBoundaryNodes: can't fix face - subfaces?");
				continue;
			}
			// recover "edge-case"
			MeshTetrahedron* tetrahedra[] = {(MeshTetrahedron*)faces[0]->getBlock(0),
				(MeshTetrahedron*)faces[0]->getBlock(1)};
			if(!tetrahedra[0]){
				tetrahedra[0] = tetrahedra[1];
				tetrahedra[1] = nullptr;
			}
			assert(tetrahedra[0]);
/*			
			if(true){
				MeshViewSet *set = new MeshViewSet;
				for(int i = 0; i < 4; i++){
					set->addFace(faces[i], 1);
					set->addPoint(points[i], 1);
					set->addEmptyBlockWithEdges(faces[i]->getOtherBlock(nullptr));
				}
				set->addPoint(bn.point, 1);
				SHOW_MESH("fix boundary edge", set);
			}
*/
			// check orientation
			if(!tetrahedra[0]->properOrientation(points[0], points[1], bn.point)){
				// switch points[]
				points[1] = bn.other_points[3];
				points[3] = bn.other_points[2];
				// switch faces[]
				MeshFace* tf = faces[0]; faces[0] = faces[3]; faces[3] = tf;
				tf = faces[1]; faces[1] = faces[2]; faces[2] = tf;
				// update tetrahedra (since face[0] has changed)
				tetrahedra[0] = (MeshTetrahedron*)faces[0]->getBlock(0);
				tetrahedra[1] = (MeshTetrahedron*)faces[0]->getBlock(1);
				if(!tetrahedra[0]){
					tetrahedra[0] = tetrahedra[1];
					tetrahedra[1] = nullptr;
				}
				assert(tetrahedra[0]);
				// should be OK now...
				assert(tetrahedra[0]->properOrientation(points[0], points[1], bn.point));
			}
			//count min len
			const DMPoint3d dpt = bn.point->getMetricCoordinates(mc);
			double min_len2 = dpt.distance2(
					bn.point->getEdge(0)->getOtherPoint(bn.point)->getMetricCoordinates(mc));
			int rank = bn.point->getRank();
			for(int j = 1; j < rank; j++){
				double len2 = dpt.distance2(
						bn.point->getEdge(j)->getOtherPoint(bn.point)->getMetricCoordinates(mc));
				if(len2 < min_len2) min_len2 = len2;
			}
//			LOG4CPLUS_INFO(MeshLog::logger_console, "fix extra point, min_len", sqrt(min_len2));
			DVector3d nv = DVector3d::crossProduct(
				points[0]->getCoordinates(), 
				points[2]->getCoordinates(), 
				points[1]->getCoordinates());
			nv *= (sqrt(min_len2) * 0.5) / mc.transformRStoMS(nv).length();
			DPoint3d moved_pt = bn.point->getCoordinates() + nv;

/*
			double vvv0 = GeometricPredicates::orient3d(
				points[0]->getCoordinates(), 
				points[1]->getCoordinates(), 
				points[2]->getCoordinates(), 
				moved_pt);
*/
			if(tetrahedra[0] && tetrahedra[1]){ // two-side
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "fix-4 two-side -> TO BE DONE CORRECTLY!");
				MeshPoint3d* new_point = new MeshPoint3d(*bn.point);
				DataVector<MeshBlock*> blocks;
				bn.point->adjacentBlocks(blocks);
				for(int j = 0; j < blocks.countInt(); j++){
					MeshBlock* block = blocks[j];
					double vol1 = DTriangle3d::orient3d(
						points[0]->getCoordinates(), 
						points[1]->getCoordinates(), 
						points[2]->getCoordinates(), 
						block->getMiddlePoint());
					if(vol1 < 0.0) // if different side -> switch block to new point
						block->switchPointsWithFaces(bn.point, new_point);
				}
				MeshGenerator3dQuality::tryMovingPoint(mc, bn.point, moved_pt, false);
				mesh->addMeshPoint(new_point);
				moved_pt = bn.point->getCoordinates() - nv;
				MeshGenerator3dQuality::tryMovingPoint(mc, new_point, moved_pt, false);
				bn.point->clearBorder();
				for(int j = 0; j < 4; j++){
					assert(bn.point->getEdgeToPoint(points[j]));
					bn.point->getEdgeToPoint(points[j])->clearBorder();
					faces[j]->clearBorder();
				}
				MeshTetrahedron* new_tetrahedra[4] = {
					new MeshTetrahedron(points[0], points[1], points[2], bn.point),
					new MeshTetrahedron(points[2], points[1], points[0], new_point),
					new MeshTetrahedron(points[2], points[3], points[0], bn.point),
					new MeshTetrahedron(points[0], points[3], points[2], new_point)};
				for(int j = 0; j < 4; j++){
					new_tetrahedra[j]->setAreaID(tetrahedra[j%2]->getAreaID());
					mesh->addMeshBlock(new_tetrahedra[j]);
				}
				for(int j = 0; j < 2; j++){
					MeshFace* new_face = new_tetrahedra[2*j]->getOppositeFace(bn.point);
					new_face->setBorder(TagBorder::INNER); // inner boundary
				}
				points[0]->getEdgeToPoint(points[2])->setBorder(TagBorder::INNER); // inner boundary
				MeshGenerator3dQuality::movePointByLaplace(mc, bn.point);
				MeshGenerator3dQuality::movePointByLaplace(mc, new_point);
				LOG4CPLUS_DEBUG(MeshLog::logger_console, "two-side-4 mesh valid **** " << mesh->isValid());
//					LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "--> additional boundary (edge-inner) node " << 
//							bn.point->getIndex() << " fixed." << endl;
				inserted_boundary_points.removeAt(i);
			}else{ // one-side
				MeshBoundaryCondition* condition = 
					(MeshBoundaryCondition*)faces[0]->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
				void * bface_ptr = faces[0]->getPtrTag(TagExtended::TAG_BOUNDARY_FACE);
				assert(mesh->isValid());
//					double vv0x = GeometricPredicates::orient3d(
//						points[0]->getCoordinates(), points[1]->getCoordinates(), 
//						points[2]->getCoordinates(), bn.point->getCoordinates());
				if(MeshGenerator3dQuality::tryMovingPoint(mc, bn.point, moved_pt, false)){
					LOG4CPLUS_DEBUG(MeshLog::logger_console, "mesh valid xx " << mesh->isValid());
					assert(mesh->isValid());

					double vv0 = GeometricPredicates::orient3d(
						points[0]->getCoordinates(), points[1]->getCoordinates(), 
						points[2]->getCoordinates(), bn.point->getCoordinates());
					LOG4CPLUS_DEBUG(MeshLog::logger_console, "Volume T[0] (no metric): " << vv0);
					if(vv0 <= 0.0) continue;
					//if(abs(vv0) < METRIC_SMALL_NUMBER) continue;
					double vv1 = GeometricPredicates::orient3d(
						points[2]->getCoordinates(), points[3]->getCoordinates(), 
						points[0]->getCoordinates(), bn.point->getCoordinates());
					LOG4CPLUS_DEBUG(MeshLog::logger_console, "Volume T[1] (no metric): " << vv1);
					if(vv1 <= 0.0) continue;
					//if(abs(vv1) < METRIC_SMALL_NUMBER) continue;

/*
					if(true){
						MeshViewSet *set = new MeshViewSet;
						for(int i = 0; i < 4; i++){
							set->addFace(faces[i], 1);
							set->addPoint(points[i], 1);
							LOG4CPLUS_INFO(MeshLog::logger_mesh, "point " << i << " -> id = " << points[i]->getID());
							set->addEmptyBlockWithEdges(faces[i]->getOtherBlock(nullptr));
						}
						set->addPoint(bn.point, 1);
						set->addPoint(moved_pt, 2);
						SHOW_MESH("fix boundary edge", set);
					}
*/
					bn.point->clearBorder();
					for(int j = 0; j < 4; j++){
						assert(bn.point->getEdgeToPoint(points[j]));
						bn.point->getEdgeToPoint(points[j])->clearBorder();
						faces[j]->clearBorder();
					}
					MeshTetrahedron* new_tetrahedra[] = {
						new MeshTetrahedron(points[0], points[1], points[2], bn.point),
						new MeshTetrahedron(points[2], points[3], points[0], bn.point)};
					for(int j = 0; j < 2; j++){
						new_tetrahedra[j]->setAreaID(tetrahedra[0]->getAreaID());
						mesh->addMeshBlock(new_tetrahedra[j]);
						MeshFace* new_face = new_tetrahedra[j]->getOppositeFace(bn.point);
						new_face->setBorder();
						if(condition)
							new_face->setPtrTag(TagExtended::TAG_BOUNDARY_COND, condition);
						if(bface_ptr)
							new_face->setPtrTag(TagExtended::TAG_BOUNDARY_FACE, bface_ptr);
					}

					MeshEdge3d* m_edge = points[0]->getEdgeToPoint(points[2]);
					MeshPoint3d* bp0 = (MeshPoint3d*)points[0]->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
					MeshPoint3d* bp1 = (MeshPoint3d*)points[2]->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
					assert(bp0 && bp1); // not sure
					MeshEdge3d* b_edge = bp0->getEdgeToPoint(bp1);
					if (b_edge) m_edge->copyBorderFlagsFrom(b_edge);
					m_edge->setBorderFlags(TagBorder::OUTER);
					MeshGenerator3dQuality::movePointByLaplace(mc, bn.point);
					assert(mesh->isValid());
					LOG4CPLUS_DEBUG(MeshLog::logger_console, "mesh valid ***** " << mesh->isValid());
					inserted_boundary_points.removeAt(i);
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, 
						"fixAdditionalBoundaryNodes: one-side planar 4-case"
						<< " - error tryMovingPoint");
				}
			}
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "fixAdditionalBoundaryNodes: n-case " << i);
		}
	}

//	LOG4CPLUS_INFO(MeshLog::logger_console, "fix extra b-nodes end, mesh valid", mesh->isValid());

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "fixAdditionalBoundaryNodes - pre-finish");
	assert(mesh->isValid());
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "fixAdditionalBoundaryNodes - finish");
	return inserted_boundary_points.empty();
}

bool MeshGenerator3dDelaunayBoundary::insertBoundaryPointForFaceRecovery(
		Metric3dContext& mc, MeshContainer3d* mesh, MissingFace& mf, 
		DataVector<MissingEdge> & missing_edges, 
		const DataVector<MissingFace> & missing_faces, DataVector<MissingFace> & new_missing_faces, 
		DataVector<MeshContainer3d::AdditionalBoundaryNode> & inserted_boundary_points, int phase)
{
	if(mf.depth >= 10) return false;

#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "InsertBoundaryPointForFaceRecovery: depth = " << mf.depth << ", " 
		<< m_pt1->getIndex() << "-" << m_pt2->getIndex() << "-" << m_pt3->getIndex() << endl;
#endif

	// - new mesh point
	MeshPoint3d* points[3] = {mf.m_pt[0], mf.m_pt[1], mf.m_pt[2]};
	MeshEdge3d* edges[3] = { mf.m_pt[0]->getEdgeToPoint(mf.m_pt[1]), 
		mf.m_pt[1]->getEdgeToPoint(mf.m_pt[2]), mf.m_pt[2]->getEdgeToPoint(mf.m_pt[0]) };
	DataVector<MeshEdge3d*> cross_edges;
	gatherCrossingEdgesForFace(mc, edges, points, cross_edges);
	assert(cross_edges.countInt() > 0);
	DPoint3d new_dpt;
	MeshTetrahedron* containing_tetrahedron = nullptr;
	if(cross_edges.countInt() > 0){
		// crossing point
		MeshPoint3d* m_pt4 = cross_edges[0]->getMeshPoint(0);
		MeshPoint3d* m_pt5 = cross_edges[0]->getMeshPoint(1);
		double vol0 = DMTriangle3d::orient3d(
			mf.m_pt[0]->getMetricCoordinates(mc),
			mf.m_pt[1]->getMetricCoordinates(mc),
			mf.m_pt[2]->getMetricCoordinates(mc),
			m_pt4->getMetricCoordinates(mc));
		double vol1 = DMTriangle3d::orient3d(
			mf.m_pt[0]->getMetricCoordinates(mc),
			mf.m_pt[2]->getMetricCoordinates(mc),
			mf.m_pt[1]->getMetricCoordinates(mc),
			m_pt5->getMetricCoordinates(mc));

		new_dpt.add(m_pt4->getCoordinates(), vol1);
		new_dpt.add(m_pt5->getCoordinates(), vol0);
		new_dpt /= (vol0+vol1);

//		MeshViewSet *set = new MeshViewSet;
//		set->addFace(points[0]->getCoordinates(), points[1]->getCoordinates(), 
//			points[2]->getCoordinates(), 1);
//		set->addEdge(cross_edges[0]);
//		set->addPoint(new_dpt, 2);
//		SHOW_MESH("new crossing point at face", set);

		containing_tetrahedron = (MeshTetrahedron*)cross_edges[0]->getFaceAt(0)->getBlock(0);
		assert(containing_tetrahedron);
	}else{
		// middle point of the face
		new_dpt = DPoint3d::average(
			mf.m_pt[0]->getCoordinates(), 
			mf.m_pt[1]->getCoordinates(), 
			mf.m_pt[2]->getCoordinates());
		MeshTetrahedron* near_block = (MeshTetrahedron*)mf.m_pt[0]->getEdge(0)->getFaceAt(0)->getBlock(0);
		containing_tetrahedron = near_block->findTetrahedronByNeighbours(new_dpt);
	}
	mc.countMetricAtPoint(new_dpt);
	MeshPoint3d* new_mpt = new MeshPoint3d(new_dpt);
	int ins_result = MeshGenerator3d::addPointToTriangulation(mc, mesh, 
							new_mpt, containing_tetrahedron);
	assert(ins_result > 0);
	if(!ins_result) return false;
	// -store for further back-recovery
	MeshContainer3d::AdditionalBoundaryNode bn(new_mpt, 3);
	bn.other_points.add(mf.m_pt[0]);
	bn.other_points.add(mf.m_pt[1]);
	bn.other_points.add(mf.m_pt[2]);
	inserted_boundary_points.add(bn);
	new_mpt->setBorder();
	// update missing_edges:
	// - remove old edge
	// check new edges
	for(int j = 0; j < 3; j++){
		MissingEdge fme(nullptr, mf.m_pt[j], new_mpt, mf.depth + 1);
		if(recoverBoundaryEdge(mc, mesh, fme, missing_faces, phase)){
			MeshEdge3d* m_edge = mf.m_pt[j]->getEdgeToPoint(new_mpt);
			assert(m_edge);
			m_edge->setBorder();
		}else{
			missing_edges.add(fme);
		}
	}
	// set new faces to be checked
	for(int j = 0; j < 3; j++){
		new_missing_faces.add(MissingFace(mf.face, points[j], points[(j+1)%3], new_mpt, mf.depth+1));
	}
	return true;
}

bool MeshGenerator3dDelaunayBoundary::insertBoundaryPointForEdgeRecovery(
		Metric3dContext& mc, MeshContainer3d* mesh,	MissingEdge& me, 
		DataVector<MissingFace> & missing_faces, DataVector<MissingEdge> & new_missing_edges, 
		DataVector<MeshContainer3d::AdditionalBoundaryNode> & inserted_boundary_points, int phase)
{
	if(me.depth >= 10) return false;

#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "InsertBoundaryPointForEdgeRecovery: depth = " << me.depth << ", " <<
		m_pt0->getIndex() << "-" << m_pt1->getIndex() << endl;
#endif

	// - new mesh point
	DPoint3d new_dpt;
	MeshTetrahedron* containing_tetrahedron = nullptr;
	PipeNode * node = nextPipeFromPoint(mc, me.m_pt0, me.m_pt1);	// build recovery-pipe
	if(node){
		new_dpt = node->pt;
		double dist0 = me.m_pt0->getCoordinates().distance(new_dpt);
		double dist1 = me.m_pt1->getCoordinates().distance(new_dpt);
		double t = dist0 / (dist0+dist1);
		if(t < 0.2){
			new_dpt = DPoint3d(me.m_pt0->getCoordinates(), me.m_pt1->getCoordinates(), 0.2);
			MeshTetrahedron* near_block = (MeshTetrahedron*)me.m_pt0->getEdge(0)->getFaceAt(0)->getBlock(0);
			containing_tetrahedron = near_block->findTetrahedronByNeighbours(new_dpt);
		}else if(t > 0.8){
			new_dpt = DPoint3d(me.m_pt0->getCoordinates(), me.m_pt1->getCoordinates(), 0.8);
			MeshTetrahedron* near_block = (MeshTetrahedron*)me.m_pt0->getEdge(0)->getFaceAt(0)->getBlock(0);
			containing_tetrahedron = near_block->findTetrahedronByNeighbours(new_dpt);
		}else{
			if(node->element_type == PIPE_BLOCK_FACE){ // face
				containing_tetrahedron = (MeshTetrahedron*)((MeshFace*)node->element2)->getBlock(0);
			}else{ // edge
				containing_tetrahedron = (MeshTetrahedron*)((MeshEdge3d*)node->element2)->getFaceAt(0)->getBlock(0);
			}
		}
	}else{
		new_dpt = DPoint3d::average(me.m_pt0->getCoordinates(), me.m_pt1->getCoordinates());
		MeshTetrahedron* near_block = (MeshTetrahedron*)me.m_pt0->getEdge(0)->getFaceAt(0)->getBlock(0);
		containing_tetrahedron = near_block->findTetrahedronByNeighbours(new_dpt);
	}

	MeshPoint3d* new_mpt = new MeshPoint3d(new_dpt);

	mc.countMetricAtPoint(new_dpt);
	int ins_result = MeshGenerator3d::addPointToTriangulation(mc, mesh, 
							new_mpt, containing_tetrahedron);
	assert(ins_result > 0);
	if(!ins_result) return false;

//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "=xx==> Boundary point: " << m_pt0->getID() << " == " << 
//		new_mpt->getID() << " == " << m_pt1->getID() << endl;

	new_mpt->setBorder();
	// -store for further back-recovery
	MeshContainer3d::AdditionalBoundaryNode bn(new_mpt, 4);
	bn.other_points.add(me.m_pt0);
	bn.other_points.add(me.m_pt1);
	// - add two new edges
	MeshBoundaryCondition* condition = nullptr;
	if(me.edge) condition = (MeshBoundaryCondition*)me.edge->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
	char bf = TagBorder::OUTER;
	if (me.edge) {
		MeshEdge3d* be = (MeshEdge3d*)me.edge->getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if (be) bf |= be->getBorderFlags();
	}
	assert(bf > 0);
	MissingEdge me0(me.edge, me.m_pt0, new_mpt, me.depth + 1);
	if(recoverBoundaryEdge(mc, mesh, me0, missing_faces, phase)){
		MeshEdge3d* m_edge = me.m_pt0->getEdgeToPoint(new_mpt);
		assert(m_edge);
		m_edge->setBorder(bf);
		if(condition) m_edge->setPtrTag(TagExtended::TAG_BOUNDARY_COND, condition);
	}else{
		new_missing_edges.add(me0);
	}

	MissingEdge me1(me.edge, me.m_pt1, new_mpt, me.depth + 1);
	if(recoverBoundaryEdge(mc, mesh, me1, missing_faces, phase)){
		MeshEdge3d* m_edge = me.m_pt1->getEdgeToPoint(new_mpt);
		assert(m_edge);
		m_edge->setBorder(bf);
		if(condition) m_edge->setPtrTag(TagExtended::TAG_BOUNDARY_COND, condition);
	}else{
		new_missing_edges.add(me1);
	}
	// update missing faces ...
	for(int j = 0; j < missing_faces.countInt(); j++){
		MissingFace mf = missing_faces[j];
		int k = 0;
		while(k < 3){
			if( (mf.m_pt[k] == me.m_pt0 && mf.m_pt[(k+1)%3] == me.m_pt1) ||
				(mf.m_pt[(k+1)%3] == me.m_pt0 && mf.m_pt[k] == me.m_pt1)) 
			{
				break;
			}else k++;
		}
		if(k < 3){ // i.e. if the missing face is incident to split-for-recovery-edge
			// remove such face from missing
			missing_faces.removeAt(j--);
			// create and add two new
			// check for middle edge
			bn.other_points.add(mf.m_pt[(k+2)%3]);
			MissingEdge fme(nullptr, mf.m_pt[(k + 2) % 3], new_mpt, me.depth + 1);
			if(recoverBoundaryEdge(mc, mesh, fme, missing_faces, phase)){
				MeshEdge3d* m_edge = mf.m_pt[(k+2)%3]->getEdgeToPoint(new_mpt);
				assert(m_edge);
				m_edge->setBorder();
			}else{
				new_missing_edges.add(fme);
			}
			if(mf.m_pt[k] == me.m_pt0){
				missing_faces.add(MissingFace(mf.face, 
						me.m_pt0, new_mpt, mf.m_pt[(k+2)%3], me.depth+1));
				missing_faces.add(MissingFace(mf.face, 
						new_mpt, me.m_pt1, mf.m_pt[(k+2)%3], me.depth+1));
			}else{
				missing_faces.add(MissingFace(mf.face, 
						me.m_pt0, mf.m_pt[(k+2)%3], new_mpt, me.depth+1));
				missing_faces.add(MissingFace(mf.face, 
						new_mpt, mf.m_pt[(k+2)%3], me.m_pt1, me.depth+1));
			}
		}
	}
	inserted_boundary_points.add(bn);
	return true;
}

MeshGenerator3dDelaunayBoundary::PipeNode* MeshGenerator3dDelaunayBoundary::buildPipe(Metric3dContext& mc, 
		MeshPoint3d* end_point, PipeNode* last_node, MeshPoint3d* start_point)
{
	// build recovery-pipe
	PipeNode * pipe = nullptr, *Pipe_end = nullptr;
	PipeType last_Pipe = (last_node)?(last_node->element_type):PIPE_START_POINT;

	DataVector<void*> loop_check;

	mc.countMetricAtPoint(DPoint3d::average(start_point->getCoordinates(), 
											end_point->getCoordinates()));

	while(true){
		PipeNode * node = nullptr;		
		switch(last_Pipe){
		case PIPE_START_POINT:
			node = nextPipeFromPoint(mc, start_point, end_point);
			break;
		case PIPE_FACE_POINT:
		case PIPE_BLOCK_POINT:
			if(last_node && !pipe){
				node = nextPipeFromPoint(mc, (MeshPoint3d*)(last_node->element2), end_point);
			}else{
				assert(Pipe_end != nullptr);
				node = nextPipeFromPoint(mc, (MeshPoint3d*)(Pipe_end->element2), end_point);
			}
			break;
		case PIPE_FACE_EDGE:
		case PIPE_BLOCK_EDGE:
			if(last_node && !pipe){
				node = nextPipeFromEdge(mc, (MeshEdge3d*)(last_node->element2), 
							last_node->pt, end_point);
			}else{
				assert(Pipe_end != nullptr);
				node = nextPipeFromEdge(mc, (MeshEdge3d*)(Pipe_end->element2), 
							Pipe_end->pt, end_point);
			}
			break;
		case PIPE_BLOCK_FACE:
			if(last_node && !pipe){
				node = nextPipeFromFace(mc, (MeshBlock*)(last_node->element1), 
					(MeshFace*)(last_node->element2), last_node->pt, end_point);
			}else{
				assert(Pipe_end != nullptr);
				node = nextPipeFromFace(mc, (MeshBlock*)(Pipe_end->element1), 
					(MeshFace*)(Pipe_end->element2), Pipe_end->pt, end_point);
			}
			break;
		}

		if(node){
			if(!loop_check.addIfNew(node->element2)){ // if already there
				LOG4CPLUS_WARN(MeshLog::logger_console, "PIPE-RECOVERY - loop encountered.");
				while(pipe){
					PipeNode* temp = pipe->next;
					delete pipe;
					pipe = temp;
				}
				return nullptr;
			}
#ifdef LOG_RECOVERY_PIPES
			switch(node->element_type){
			case PIPE_FACE_POINT:
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "  -[face]->point.");
				break;
			case PIPE_FACE_EDGE:
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "  -[face]->edge[" << ((MeshEdge3d*)(node->element2))->getFaceCount() << "].");
				break;
			case PIPE_BLOCK_POINT:
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "  -[block]->point.");
				break;
			case PIPE_BLOCK_EDGE:
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "  -[block]->edge[" << ((MeshEdge3d*)(node->element2))->getFaceCount() << "].");
				break;
			case PIPE_BLOCK_FACE:
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "  -[block]->face.");
				break;
			case PIPE_START_POINT:
				break;
			}
#endif
			// insert into the linked list
			last_Pipe = node->element_type;
			if(pipe){
				Pipe_end = Pipe_end->next = node;
			}else{
				pipe = Pipe_end = node;
			}
			
			// check if ending point reached
			if(node->element_type <= PIPE_BLOCK_POINT){
#ifdef LOG_RECOVERY_PIPES
				if((MeshPoint3d*)(node->element2) == end_point){
					LOG4CPLUS_INFO(MeshLog::logger_mesh, " * pipe finished successfully.");
				}else{
					LOG4CPLUS_INFO(MeshLog::logger_mesh, " * pipe finished with a wrong point ???.");
				}
#endif
				// Leave the loop
				return pipe;
			}
		}else{
			while(pipe){
				PipeNode* temp = pipe->next;
				delete pipe;
				pipe = temp;
			}
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_WARN(MeshLog::logger_console, "PIPE-RECOVERY - broken pipe.");
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " * Broken pipe (leaving ...).");
#endif
			return nullptr;
		}		
	}
}

bool MeshGenerator3dDelaunayBoundary::insertPointAtLongestAdjacentEdge(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshPoint3d *point, const DataVector<MissingFace> & missing_faces)
{
#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "insertPointAtLongestAdjacentEdge.");
#endif
	mc.countMetricAtPoint(point->getCoordinates());

	double max_len = 2.0;
	MeshEdge3d* best_edge = 0;

	int rank = point->getRank();
	int mct = missing_faces.countInt();
	for(int i = 0; i < rank; i++){
		MeshEdge3d* edge = point->getEdge(i);
		if(edge->isBorder()) continue;
		double len = edge->getLength(mc);
		if(len > max_len){
			// check for missing_faces (shouldn't be really necessary...)
			const DPoint3d dpt = edge->getPoint(0.5);
			bool safe = true;
			for(int j = 0; safe && j < mct; j++){
				const FaceNode node(
					missing_faces[j].m_pt[0]->getCoordinates(), 
					missing_faces[j].m_pt[1]->getCoordinates(),
					missing_faces[j].m_pt[2]->getCoordinates());
				if(dpt.distance2(node.middle) > node.max_dist2) continue;
				if(DTriangle3d::distance2ToPoint(dpt, node.pt, node.v0, node.v1) < node.min_dist2){
					//LOG4CPLUS_WARN(MeshLog::logger_console, "insertPointAtLongestAdjacentEdge - too near some missing face???");
					safe = false;
				}
			}
			if(safe){
				best_edge = edge;
				max_len = len;
			}
		}
	}
	if(!best_edge) return false;

#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "longest-edge= " << max_len);
#endif

#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
	if(!mesh->isValid()){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "insertPointAtLongestAdjacentEdge start - Mesh invalid, leaving."); return false;
	}
#endif

	MeshPoint3d* new_point = new MeshPoint3d(best_edge->getPoint(0.5));
	mesh->addMeshPoint(new_point);
	MeshTetrahedron* tetrahedron = (MeshTetrahedron*)best_edge->getFaceAt(0)->getBlock(0);
	if(!tetrahedron) tetrahedron = (MeshTetrahedron*)best_edge->getFaceAt(0)->getBlock(1);
	bool success = MeshGenerator3d::addPointToTriangulationBySphere(mc, mesh, 
						new_point, tetrahedron, nullptr, true);
//	bool success = addPointToTriangulationBySwap(mc, mesh, new_point, tetrahedron);
#ifdef LOG_RECOVERY_PIPES
	if(!success)
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "insertPointAtLongestAdjacentEdge - error adding point to mesh.");
#endif
#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
	if(!mesh->isValid()){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "insertPointAtLongestAdjacentEdge done - Mesh invalid, leaving. Last point=", 
			new_point->getID()); return false;
	}
#endif

	return success;
}

MeshGenerator3dDelaunayBoundary::PipeNode* MeshGenerator3dDelaunayBoundary::buildPipeWithInsertingNodes(
	Metric3dContext& mc, MeshContainer3d* mesh, 
	const DataVector<MissingFace> & missing_faces, 
	MeshPoint3d* &pt2, MeshPoint3d* &pt1, int phase)
{
	PipeNode * pipe = buildPipe(mc, pt2, nullptr, pt1);
	if(pipe) return pipe;

#ifdef LOG_RECOVERY_PIPES
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Checking mesh...");
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, (mesh->isValid() ? "OK" : "Invalid"));
#endif
	pipe = buildPipe(mc, pt1, nullptr, pt2);
	if(pipe){
#ifdef LOG_RECOVERY_PIPES
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Broken pipe found, if starting from the other side...");
#endif
		LOG4CPLUS_WARN(MeshLog::logger_console, "Broken pipe found, if starting from the other side...");
		MeshPoint3d* temp = pt1;
		pt1 = pt2;
		pt2 = temp;
		return pipe;
	}

	if(phase < 2) return nullptr;
#ifdef LOG_RECOVERY_PIPES
//		MeshView::showDebugMesh("Broken pipe... (from start)", mesh, pt1, pt2);
#endif
	// try inserting point
	const int MAX_TRIES = 100;
	for(int i = 0; i < MAX_TRIES; i++){
		bool success1 = insertPointAtLongestAdjacentEdge(mc, mesh, pt1, missing_faces);
		if(success1){
			LOG4CPLUS_WARN(MeshLog::logger_console, "Point inserted for broken pipe.");
			if(pt1->getEdgeToPoint(pt2) != nullptr) return nullptr;
			else if((pipe = buildPipe(mc, pt2, nullptr, pt1))) return pipe;
			else if((pipe = buildPipe(mc, pt1, nullptr, pt2))){
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "Broken pipe found, if starting from the other side...");
#endif
				MeshPoint3d* temp = pt1;
				pt1 = pt2;
				pt2 = temp;
				return pipe;
			}
			//MeshView::showDebugMesh("Broken pipe - after inserting (pt1)", mesh, pt1, pt2);
		}
		bool success2 = insertPointAtLongestAdjacentEdge(mc, mesh, pt2, missing_faces);
		if(success2){
			LOG4CPLUS_WARN(MeshLog::logger_console, "Point inserted for broken pipe.");
			if(pt1->getEdgeToPoint(pt2) != nullptr) return nullptr;
			else if((pipe = buildPipe(mc, pt2, nullptr, pt1))) return pipe;
			else if((pipe = buildPipe(mc, pt1, nullptr, pt2))){
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "Broken pipe found, if starting from the other side...");
#endif
				MeshPoint3d* temp = pt1;
				pt1 = pt2;
				pt2 = temp;
				return pipe;
			}
			//MeshView::showDebugMesh("Broken pipe - after inserting (pt2)", mesh, pt1, pt2);
		}
		if(!success1 && !success2) return nullptr; // no point could be inserted
	}
	return nullptr;
}

MeshViewSet* MeshGenerator3dDelaunayBoundary::createPipeViewSet(PipeNode* pipe, 
						MeshPoint3d* start_point, MeshPoint3d* last_point)
{
	MeshViewSet* set = new MeshViewSet;
	set->addEdge(start_point->getCoordinates(), last_point->getCoordinates(), 3);
	set->addPoint(start_point, 2);
	set->addPoint(last_point, 3);

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Missing edge: " << start_point->getIndex() << " - " << last_point->getIndex());
	ostringstream log_pipe;
	log_pipe << "pipe: ";

	for( ; pipe && pipe->next; pipe = pipe->next){
		switch(pipe->element_type){
		case PIPE_FACE_POINT:
		case PIPE_BLOCK_POINT:
			log_pipe << " FACE|BLOCK_POINT ";
			set->addPoint((MeshPoint3d*)pipe->element2, 3);
			break;
		case PIPE_FACE_EDGE:
		case PIPE_BLOCK_EDGE:
			{
				log_pipe << " FACE|BLOCK_EDGE ";
				MeshEdge3d* edge = (MeshEdge3d*)pipe->element2;
				set->addEdge(edge);
				int rank = edge->getFaceCount();
				DataVector<MeshBlock*> blocks;
				DataVector<MeshPoint3d*> points;
				for(int i = 0; i < rank; i++){
					blocks.addIfNew(edge->getFaceAt(i)->getBlock(0));
					blocks.addIfNew(edge->getFaceAt(i)->getBlock(1));
				}
				for(int i = 0; i < blocks.countInt(); i++){
					set->addBlockWithEdges(blocks[i]);
					int pct = blocks[i]->getPointCount();
					for(int j = 0; j < pct; j++)
						points.addIfNew(blocks[i]->getPoint(j));
				}
				for(int i = 0; i < points.countInt(); i++){
					set->addPoint(points[i]);
				}
			}
			break;
		case PIPE_BLOCK_FACE:
			log_pipe << " BLOCK_FACE ";
			set->addFace((MeshFace*)pipe->element2, 3);
			break;
		case PIPE_START_POINT:
			break;
		}
	}
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh_stat, log_pipe.str());
	return set;
}

int MeshGenerator3dDelaunayBoundary::insertPointInWorstBlock(Metric3dContext& mc, 
				MeshContainer3d* mesh, MeshPoint3d* /* start_point */, 
				PipeNode* pipe, const DataVector<FaceNode> &missing_faces)
{
	MeshTetrahedron* worst_block = nullptr;
	double min_r = -1.0;
	int block_count = 0;
	MeshPoint3d* end_point = nullptr;
	for(PipeNode* p = pipe; p; p = p->next){
		if(p->element_type == PIPE_BLOCK_POINT ||
			p->element_type == PIPE_BLOCK_EDGE ||
			p->element_type == PIPE_BLOCK_FACE)
		{
			if(p->element_type == PIPE_BLOCK_POINT)
				end_point = (MeshPoint3d*)p->element2;
			MeshTetrahedron* block = (MeshTetrahedron*)p->element1;
			++block_count;
			//const DPoint3d middle = block->getMiddlePoint();
			const DPoint3d middle = mc.transformMStoRS(block->getOuterSphereCenter(mc));
			bool acceptable = true;
			// check for collisions with faces-to-be-recovered
			for(int i = 0; i < missing_faces.countInt(); i++){
				const FaceNode& node = missing_faces.get(i);
				if(middle.distance2(node.middle) > node.max_dist2)
					continue;
				if(DTriangle3d::distance2ToPoint(middle, node.pt, node.v0, node.v1) < node.min_dist2){
					acceptable = false;
					break;
				}
			}
			if(acceptable){
				double r = block->getOuterSphereRadius(mc);
				if(r > min_r){
					min_r = r;
					worst_block = block;
				}
			}
		}else if(p->element_type == PIPE_FACE_EDGE &&
					p->next->element_type == PIPE_FACE_POINT)
		{
			end_point = (MeshPoint3d*)p->next->element2;
			MeshFace* first_face = (MeshFace*)p->element1;
			MeshFace* next_face = (MeshFace*)p->next->element1;
			MeshEdge3d* edge = (MeshEdge3d*)p->element2;
			int edge_rank = edge->getFaceCount();
			for(int i = 0; i < edge_rank; i++){
				MeshFace* face = edge->getFaceAt(i);
				if(face == first_face || face == next_face) continue;
				for(int j = 0; j < 2; j++){
					MeshTetrahedron* block = (MeshTetrahedron*)face->getBlock(j);
					if(block->isAdjacentTo(first_face) || block->isAdjacentTo(next_face)) continue;
					const DPoint3d middle = block->getMiddlePoint();
					bool acceptable = true;
					// check for collisions with faces-to-be-recovered
					for(int k = 0; k < missing_faces.countInt(); k++){
						const FaceNode& node = missing_faces.get(k);
						if(middle.distance2(node.middle) > node.max_dist2) continue;
						if(DTriangle3d::distance2ToPoint(middle, node.pt, node.v0, node.v1) < node.min_dist2){
							acceptable = false;
							break;
						}
					}
					if(acceptable){
						double r = block->getOuterSphereRadius(mc, false);
						if(r > min_r){
							min_r = r;
							worst_block = block;
						}
					}
				}
			}
		}
	}

	if(worst_block){
		// insert node
		DPoint3d new_pt = mc.transformMStoRS(worst_block->getOuterSphereCenter(mc));
		new_pt = mesh->getBoundingBox().fitInPoint(new_pt);

		MeshTetrahedron* ct = worst_block->findTetrahedronByNeighbours(new_pt);
		if(!ct){
			ct = mesh->getNearTetrahedron(new_pt);
			assert(ct);
			ct = ct->findTetrahedronByNeighbours(new_pt);
		}
		if(!ct){
#ifdef LOG_RECOVERY_PIPES
			//MeshViewSet* set = mesh->getViewSet();
			//set->addPoint(new_pt, 1);
			//SHOW_MESH("linear search?", set);
			LOG4CPLUS_WARN(MeshLog::logger_console, "Linear search for containing tetrahedron");
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " -- linear search for containing tetrahedron");
#endif
			int bct = mesh->getBlocksCount();
			for(int i = 0; i < bct; i++){
				MeshTetrahedron* block = (MeshTetrahedron*)mesh->getBlockAt(i);
				if(block->isPointInside(new_pt)){
					ct = block;
					break;
				}
			}
		}
		if(!ct){
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " -- insertion cancelled, can't find containing tetrahedron");
#endif
			return 0;
		}

		for(int i = 0; i < 4; i++){
			double dist2 = mc.transformRStoMS(new_pt - ct->getPoint(i)->getCoordinates()).length2();
			if(dist2 < 0.125){
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " -- insertion cancelled, distance to some vertex= " << sqrt(dist2));
#endif
				return 0;	// too near (len < 0.3 in metric space)
			}
		}

		MeshPoint3d* new_point = new MeshPoint3d(new_pt);
#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Inserting node to recover pipe - Mesh invalid, leaving."); return 0;
		}
#endif

/*
		{
			MeshViewSet* set = createPipeViewSet(pipe, start_point, end_point);
			set->addPoint(new_point, 4);
			SHOW_MESH("insertingPointInWorstBlock", set);
		}
*/
		bool success = MeshGenerator3d::addPointToTriangulation(mc, mesh, 
							new_point, ct, ct, true);
#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
		if(!mesh->isValid()){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Inserted node to recover pipe - Mesh invalid, leaving."); return 0;
		}
#endif
		return success ? 1 : 0;
	}else{
//		MeshViewSet* set = createPipeViewSet(pipe, start_point, end_point);
//		SHOW_MESH("insertingPointInWorstBlock - none found", set);
	}

	return 0;
}

int MeshGenerator3dDelaunayBoundary::insertPointAtLongestPipeEdges(Metric3dContext& mc, 
			MeshContainer3d* mesh, MeshPoint3d* start_point, 
			PipeNode* pipe, const DataVector<MissingFace> & missing_faces)
{
	DataVector<MeshPoint3d*> check_points;
	check_points.add(start_point);
	for( ; pipe; pipe = pipe->next){
		switch(pipe->element_type){
		case PIPE_FACE_POINT:
		case PIPE_BLOCK_POINT:
			check_points.add((MeshPoint3d*)pipe->element2);
			break;
		case PIPE_FACE_EDGE:
		case PIPE_BLOCK_EDGE:
			check_points.addIfNew(((MeshEdge3d*)pipe->element2)->getMeshPoint(0));
			check_points.addIfNew(((MeshEdge3d*)pipe->element2)->getMeshPoint(1));
			break;
		case PIPE_BLOCK_FACE:
			check_points.addIfNew(((MeshFace*)pipe->element2)->getPoint(0));
			check_points.addIfNew(((MeshFace*)pipe->element2)->getPoint(1));
			check_points.addIfNew(((MeshFace*)pipe->element2)->getPoint(2));
			break;
		case PIPE_START_POINT:
			break;
		}
	}
	int count = check_points.countInt();
	int inserted = 0;
	for(int i = 0; i < count; i++)
		if(insertPointAtLongestAdjacentEdge(mc, mesh, check_points[i], missing_faces))
			++inserted;

	return inserted;
}

int MeshGenerator3dDelaunayBoundary::recoverBoundaryEdgeWithSpecialMetricSwap(Metric3dContext& mc, 
		MeshContainer3d* mesh, MeshPoint3d *pt1, MeshPoint3d *pt2, int max_layers)
{
#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "EdgeRecovery (with special metric swap): " << 
		pt1->getIndex() << "-" << pt2->getIndex() << endl;
#endif

	int mpct = mesh->getPointsCount() - 8; // no outside points

	DataVector<MeshPoint3d*> active_points(mpct);
	pt1->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); active_points.add(pt1);
	pt2->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); active_points.add(pt2);
	MeshTetrahedron** tetrahedra = new MeshTetrahedron*[3];

	int s32_checked = 0;
	int s32_swapped = 0;
	int s23_checked = 0;
	int s23_swapped = 0;

	int swap_count[3] = {-100, -100, -100};
	int layer = 0;

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "======== SPECIAL METRIC SWAP FOR EDGE ============");
	while(swap_count[0] != 0){
		swap_count[0] = 0;
		for(int i = 0; i < active_points.countInt(); i++){
			MeshPoint3d* point = active_points[i];
			for(int j = 0; j < point->getRank(); j++){
				// check edge
				MeshEdge3d* edge = point->getEdge(j);
				if(edge->isBorder()) continue;
				bool swap_ok = false;
				if(edge->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){ // if tagged -> already checked
					++s32_checked;
					swap_ok = MeshGenerator3d::swap32(mc, mesh, edge, 
								tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
				}
				if(swap_ok){
					++swap_count[0];
					++s32_swapped;
					j = 0; 
					for(int m = 0; m < 2; m++) 
						tetrahedra[m]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
					continue; 
				}else edge->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
				// check faces
				for(int k = 0; k < edge->getFaceCount(); k++){
					MeshFace* face = edge->getFaceAt(k);
					if(face->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){ // if tagged -> already checked
						++s23_checked;
						swap_ok = MeshGenerator3d::swap23(mc, mesh, face, 
									tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
					}
					if(swap_ok){ // start again from first edge
						++swap_count[0];
						++s23_swapped;
						if(pt1->getEdgeToPoint(pt2)){ 
							LOG4CPLUS_DEBUG(MeshLog::logger_console, "**** > edge successfully recovered");
#ifdef LOG_RECOVERY_PIPES
							LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "*** fixed, swap23: " << s23_checked << '/' << s23_swapped <<
								", swap32: " << s32_checked << '/' << s32_swapped << endl;
#endif
							delete[] tetrahedra;
							mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);
							return REDUCE_TRANSFORM; 
						}
						j = 0; 
						for(int m = 0; m < 3; m++) 
							tetrahedra[m]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
						break; 
					}else face->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
				}
			}
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "----> swap_count = " << swap_count[0]);
		if((swap_count[0] == 0) && active_points.countInt() < mpct){ // expand set of active points (by neighbors)
			if(++layer > max_layers) break;
			int act = active_points.countInt();
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
					}
				}
			}
//			if(layer > 1) LOG4CPLUS_INFO(MeshLog::logger_console, "[edge special metric swap] --> next layer", layer);
//			LOG4CPLUS_DEBUG(MeshLog::logger_console, "--> active points", active_points.countInt());
			assert(swap_count[0] != 0);
		}
		if(swap_count[2] == swap_count[1] && swap_count[2] == swap_count[0]){
			LOG4CPLUS_DEBUG(MeshLog::logger_console, " ==== [edge] loop check exit =====");
			break;
		}
		swap_count[2] = swap_count[1];
		swap_count[1] = swap_count[0];
	}

//	LOG4CPLUS_DEBUG(MeshLog::logger_console, "======== END SPECIAL METRIC SWAP FOR EDGE ========> active/all", active_points.countInt() / (double)mpct);
	delete[] tetrahedra;
#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "--- failed, swap23: " << s23_checked << '/' << s23_swapped <<
			", swap32: " << s32_checked << '/' << s32_swapped << endl;
#endif
	mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);
	return REDUCE_FAILED;
}

int MeshGenerator3dDelaunayBoundary::recoverBoundaryFaceWithSpecialMetricSwap(Metric3dContext& mc, 
		MeshContainer3d* mesh, MeshPoint3d *pt1, MeshPoint3d *pt2, MeshPoint3d *pt3, int max_layers)
{
#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "FaceRecovery (with special metric swap): " << 
		pt1->getIndex() << "-" << pt2->getIndex() << "-" << pt3->getIndex() << endl;
#endif

	int mpct = mesh->getPointsCount() - 8;

	DataVector<MeshPoint3d*> active_points(mpct);
	pt1->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); active_points.add(pt1);
	pt2->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); active_points.add(pt2);
	pt3->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); active_points.add(pt3);
	MeshTetrahedron** tetrahedra = new MeshTetrahedron*[3];

	int s32_checked = 0;
	int s32_swapped = 0;
	int s23_checked = 0;
	int s23_swapped = 0;
	int layer = 0;

	int swap_count[3] = {-100, -100, -100};
	while(swap_count[0] != 0){
		swap_count[0] = 0;
		for(int i = 0; i < active_points.countInt(); i++){
			MeshPoint3d* point = active_points[i];
			for(int j = 0; j < point->getRank(); j++){
				// check edge
				MeshEdge3d* edge = point->getEdge(j);
				if(edge->isBorder()) continue;
				bool swap_ok = false;
				if(edge->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){ // if tagged -> already checked
					++s32_checked;
					swap_ok = MeshGenerator3d::swap32(mc, mesh, edge, 
								tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
				}
				if(swap_ok){ 
					++s32_swapped;
					j = 0; 
					if(pt1->getFaceToPoints(pt2, pt3)){
						delete[] tetrahedra; 
#ifdef LOG_RECOVERY_PIPES
						LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "*** fixed, swap23: " << s23_checked << '/' << s23_swapped <<
							", swap32: " << s32_checked << '/' << s32_swapped << endl;
#endif
						mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);
						return REDUCE_TRANSFORM; 
					}
					for(int m = 0; m < 2; m++) 
						tetrahedra[m]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
					++swap_count[0]; 
					continue; 
				}else edge->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
				// check faces
				for(int k = 0; k < edge->getFaceCount(); k++){
					MeshFace* face = edge->getFaceAt(k);
					if(face->zeroIntTag(TagExtended::TAG_MG3D_SM_SWAP)){ // if tagged -> already checked
						++s23_checked;
						swap_ok = MeshGenerator3d::swap23(mc, mesh, face, 
									tetrahedra, MeshData::SWAP3_MIN_TETRAHEDRON_QUALITY);
					}
					if(swap_ok){ // start again from first edge
						++s23_swapped;
						j = 0; 
						for(int m = 0; m < 3; m++) 
							tetrahedra[m]->removeTagForEdgesAndFaces(TagExtended::TAG_MG3D_SM_SWAP);
						++swap_count[0]; 
						break; 
					}else face->setIntTag(TagExtended::TAG_MG3D_SM_SWAP); // checked
				}
			}
		}
		if((swap_count[0] == 0) && active_points.countInt() < mpct){ // expand set of active points (by neighbors)
			if(++layer > max_layers) break;
			int act = active_points.countInt();
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
					}
				}
			}
			assert(swap_count[0] != 0);
		}

		if(swap_count[2] == swap_count[1] && swap_count[2] == swap_count[0]){
			LOG4CPLUS_DEBUG(MeshLog::logger_console, " ==== [face] loop check exit =====");
			break;
		}
		swap_count[2] = swap_count[1];
		swap_count[1] = swap_count[0];
	}

	delete[] tetrahedra;
#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "--- failed, swap23: " << s23_checked << '/' << s23_swapped <<
			", swap32: " << s32_checked << '/' << s32_swapped << endl;
#endif
	mesh->removeAllTags(TagExtended::TAG_MG3D_SM_SWAP);
	return REDUCE_FAILED;
}

int MeshGenerator3dDelaunayBoundary::recoverBoundaryEdge(Metric3dContext& mc, MeshContainer3d* mesh, 
	MissingEdge& me, const DataVector<MissingFace> & missing_faces, int phase)
{

//	if(phase == 4){ LOG4CPLUS_INFO(MeshLog::logger_mesh, "=xx==> Recover edge: " << pt1->getID() << " == " << pt2->getID()); }
	MeshPoint3d* pt1 = me.m_pt0;
	MeshPoint3d* pt2 = me.m_pt1;
	if(pt1->getEdgeToPoint(pt2) != nullptr)
		return REDUCE_CLEAN;

#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "EdgeRecovery: " << pt1->getIndex() << "-" << pt2->getIndex());
#endif

//	if (phase > 3) {
	if (false) {
		recoverBoundaryWithLocalCavityRemeshing(mc, mesh, pt1, pt2, missing_faces);
		if(pt1->getEdgeToPoint(pt2) != nullptr)
			return REDUCE_CLEAN;
	}
	assert(mesh->isValid());

	// build recovery-pipe
	PipeNode * pipe = buildPipeWithInsertingNodes(mc, mesh, missing_faces, pt2, pt1, phase);	

	if(!pipe)
		return (pt1->getEdgeToPoint(pt2)) ? REDUCE_CLEAN : REDUCE_FAILED;

	//if(phase > 3){
	if (false) {
		MeshViewSet* set = createPipeViewSet(pipe, pt1, pt2);
		SHOW_MESH("pipe (edge) for phase 3", set);

		// cut through ...
		assert(pipe->next != nullptr);
		assert(me.inserted_points.empty());
		while (pipe->next != nullptr) {
			MeshPoint3d* pt1_new = cutFirstPipeNode(mc, mesh, pt1, pt2, pipe);
			assert(pt1_new != nullptr);
			me.inserted_points.add(pt1_new);
			PipeNode* pipe_new = buildPipeWithInsertingNodes(mc, mesh, missing_faces, pt2, pt1, phase);
			if (!pipe_new) return REDUCE_FAILED;

			MeshViewSet* set = createPipeViewSet(pipe_new, pt1, pt2);
			SHOW_MESH("new_pipe (edge) for phase 3", set);

			// start from new point
			pt1 = pt1_new;
			// clean the pipe
			while (pipe) { PipeNode* temp = pipe->next; delete pipe; pipe = temp; }
			// ... and set new one
			pipe = pipe_new;
		}

		return REDUCE_CUT_SPLIT;
	}

#ifdef LOG_RECOVERY_PIPES
//	if(phase > 1){
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Checking mesh (before solving pipe)...");
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, (mesh->isValid() ? "OK" : "Invalid"));
//	}
#endif

	// pipe built successfully, start the recovery process
	assert(pipe->next != nullptr);
	int modification_count = 0;
	int inserted_count = 0;
	int transform_vicinity_count = 0;
	int max_invalid_pipe_count = 3;
	while(pipe->next != nullptr){	// i.e. while pipe has at least two chains
		bool no_change = true;
		PipeNode* last_node = nullptr;
		// try to recover
		for(PipeNode* node = pipe; node->next != nullptr; ){
			bool invalid_Pipe; // is set inside "reducePipeNodes"
			if(reducePipeNodes(mc, mesh, last_node, node, phase, invalid_Pipe)){
				++modification_count;
				no_change = false;
				if(!invalid_Pipe && node->next == nullptr) break;
			}else{
				last_node = node;
				node = node->next;
			}
			if(invalid_Pipe){
				// clean the pipe
				while(pipe){PipeNode* temp = pipe->next; delete pipe;pipe = temp; }
				if(--max_invalid_pipe_count < 0) return REDUCE_FAILED;
				// continue
				pipe = buildPipeWithInsertingNodes(mc, mesh, missing_faces, pt2, pt1, phase);
				if(!pipe){
#ifdef LOG_RECOVERY_PIPES
//					MeshView::showDebugMesh("Broken pipe (during transform) ... ", mesh, pt1, pt2);
#endif
					return REDUCE_FAILED;
				}
				break;
			}
		}
		if(no_change && (phase > 1)){
			// try to transform
			for(PipeNode* node = pipe; node->next != nullptr; node=node->next){
				if(transformPipeNodes(mc, mesh, node)){
					++modification_count;
					// clean the pipe
					while(pipe){PipeNode* temp = pipe->next;delete pipe;pipe = temp;}
					// continue
					pipe = buildPipeWithInsertingNodes(mc, mesh, missing_faces, pt2, pt1, phase);
					if(!pipe){
#ifdef LOG_RECOVERY_PIPES
//						MeshView::showDebugMesh("Broken pipe (during transform) ... ", mesh, pt1, pt2);
#endif
						return REDUCE_FAILED;
					}
					no_change = false;
					break;
				}
			}
			DataVector<int> old_faces(30);
			DataVector<int> old_edges(20);
			if((transform_vicinity_count < 4) && transformPipeVicinity(mc, mesh, pipe, old_faces, old_edges)){
				++transform_vicinity_count;
				// clean the pipe
				while(pipe){PipeNode* temp = pipe->next;delete pipe;pipe = temp;}
				// continue
				pipe = buildPipeWithInsertingNodes(mc, mesh, missing_faces, pt2, pt1, phase);
				if(!pipe){
#ifdef LOG_RECOVERY_PIPES
//					MeshView::showDebugMesh("Broken pipe (during transform) ... ", mesh, pt1, pt2);
#endif
					return REDUCE_FAILED;
				}
				no_change = false;
			}
		}
		// failed
/*
		if(no_change && (phase > 2) && inserted_count < 10){
			int inserted = insertPointInWorstBlock(mc, mesh, pt1, pipe, missing_faces);
			//int inserted = insertPointAtLongestPipeEdges(mc, mesh, pt1, pipe, missing_faces);
			if(inserted > 0){
				// clean the pipe
				while(pipe){
					PipeNode* temp = pipe->next;
					delete pipe;
					pipe = temp;
				}
				if(pt1->getEdgeToPoint(pt2) != nullptr){
#ifdef LOG_RECOVERY_PIPES
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "-> edge recovered after insertion.");
#endif
					return REDUCE_POINT_CLEAN;
				}
				// continue
				pipe = buildPipeWithInsertingNodes(mc, mesh, missing_faces, pt2, pt1, phase);
				if(!pipe){
#ifdef LOG_RECOVERY_PIPES
//					MeshView::showDebugMesh("Broken pipe... (after insertion)", mesh, pt1, pt2);
#endif
					return REDUCE_FAILED;
				}
				inserted_count += inserted;
				no_change = false;
			}
		}
*/
		if(no_change){
			// Clean the pipe
			while(pipe){
				PipeNode* temp = pipe->next;
				delete pipe;
				pipe = temp;
			}

#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " x - failed.");
#endif
			return REDUCE_FAILED;
		}
	}

	assert(pipe->next == nullptr);
	assert(pt1->getEdgeToPoint(pt2) != nullptr);

	while(pipe){
		PipeNode* temp = pipe->next;
		delete pipe;
		pipe = temp;
	}
	
#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh, " v - fixed.");
	assert(mesh->isValid());
#endif
	return (inserted_count==0)?REDUCE_TRANSFORM:REDUCE_POINT_TRANSFORM;

}

int MeshGenerator3dDelaunayBoundary::insertPointForFaceRecovery(Metric3dContext& mc, MeshContainer3d* mesh,
		const DPoint3d& pt, const DataVector<MissingFace> &missing_faces, MeshTetrahedron* near_block)
{
	for(int k = 0; k < missing_faces.countInt(); k++){
		const FaceNode node(
			missing_faces[k].m_pt[0]->getCoordinates(), 
			missing_faces[k].m_pt[1]->getCoordinates(),
			missing_faces[k].m_pt[2]->getCoordinates());
		if(pt.distance2(node.middle) > node.max_dist2) continue;
		if(DTriangle3d::distance2ToPoint(pt, node.pt, node.v0, node.v1) < node.min_dist2)
			return false;
	}

	if(near_block){
		near_block = near_block->findTetrahedronByNeighbours(pt);
		if(!near_block){
			near_block = mesh->getNearTetrahedron(pt);
			assert(near_block);
			near_block = near_block->findTetrahedronByNeighbours(pt);
		}
	}else{
		near_block = mesh->getNearTetrahedron(pt);
		assert(near_block);
		near_block = near_block->findTetrahedronByNeighbours(pt);
	}
	if(!near_block){
#ifdef LOG_RECOVERY_PIPES
		//MeshViewSet* set = mesh->getViewSet();
		//set->addPoint(pt, 1);
		//SHOW_MESH("linear search?", set);
		LOG4CPLUS_WARN(MeshLog::logger_console, "Linear search for containing tetrahedron");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " -- linear search for containing tetrahedron");
#endif
		int bct = mesh->getBlocksCount();
		for(int i = 0; i < bct; i++){
			MeshTetrahedron* block = (MeshTetrahedron*)mesh->getBlockAt(i);
			if(block->isPointInside(pt)){
				near_block = block;
				break;
			}
		}
	}
	if(!near_block){
#ifdef LOG_RECOVERY_PIPES
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " -- insertion cancelled, can't find containing tetrahedron");
#endif
		return 0;
	}

	const DMPoint3d dpt = mc.transformRStoMS(pt);
	for(int i = 0; i < 4; i++){
		double dist2 = dpt.distance2(near_block->getPoint(i)->getMetricCoordinates(mc));
		if(dist2 < 0.25){ //(0.5)^2
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " -- insertion cancelled, distance to some vertex= " << sqrt(dist2));
#endif
			return 0;	// too near (len < 0.25 in metric space)
		}
	}

	MeshPoint3d* point = new MeshPoint3d(pt);

	if(!MeshGenerator3d::addPointToTriangulation(mc, mesh, point, near_block, nullptr, true)){
		delete point;
		return 0;
	}
#ifdef CONSTRAINING_MESH_VALIDITY_CHECK
	if(!mesh->isValid()){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Inserted node to recover face - Mesh invalid, leaving."); return 0;
	}
#endif

#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh, " -- inserted node id " << point->getID());
#endif

	return 1;
}

int MeshGenerator3dDelaunayBoundary::recoverBoundaryFace(Metric3dContext& mc, MeshContainer3d* mesh, 
			const DataVector<MissingFace> & missing_faces, 
			MeshPoint3d *pt1, MeshPoint3d *pt2, MeshPoint3d *pt3, int phase)
{
	if(pt1->getFaceToPoints(pt2, pt3) != nullptr)
		return REDUCE_CLEAN;

#ifdef LOG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "FaceRecovery:");
#endif

//	if (phase > 3) {
	if(false){
		recoverBoundaryWithLocalCavityRemeshing(mc, mesh, pt1, pt2, pt3);
		if(pt1->getFaceToPoints(pt2, pt3) != nullptr)
			return REDUCE_CLEAN;
	}
	assert(mesh->isValid());

	MeshPoint3d* points[3] = {pt1, pt2, pt3};
	MeshEdge3d* edges[3];
	short missing_count = 0;
	for(int i = 0; i < 3; i++){
		edges[i] = points[i]->getEdgeToPoint(points[(i+1)%3]);
		if(edges[i] == nullptr) ++missing_count;
	}

	if(missing_count > 0 && phase == 0){ // try inserting single points
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " - " << missing_count << " edges missing, trying inserting points. ");
#endif
		const DMPoint3d dmiddle = DMPoint3d::average(
			points[0]->getMetricCoordinates(mc),
			points[1]->getMetricCoordinates(mc),
			points[2]->getMetricCoordinates(mc));
		DMVector3d nv = DMVector3d::crossProduct(
			points[0]->getMetricCoordinates(mc),
			points[1]->getMetricCoordinates(mc), 
			points[2]->getMetricCoordinates(mc));
		nv *= (MeshTetrahedron::opt_metric_height / nv.length()); 
		// for this face
		const DBox box = mesh->getBoundingBox();
		DPoint3d pts[2] = { mc.transformMStoRS(dmiddle+nv), // from metric to real space
			mc.transformMStoRS(dmiddle-nv)};
		int inserted = 0;
		const DPoint3d middle = mc.transformMStoRS(dmiddle);
		// insert two points (with special metric)
		for(int i = 0; i < 2; i++){
			if(!box.contains(pts[i]))
				pts[i] = DPoint3d(middle, box.fitInPoint(pts[i]), 0.75);
			inserted += insertPointForFaceRecovery(mc, mesh, pts[i], missing_faces,
				(MeshTetrahedron*)pt1->getEdge(0)->getFaceAt(0)->getBlock(0));
		}
		//mc.countMetricAtPoint(dmiddle);
		// try again recovering edges
		if(inserted > 0){
			missing_count = 0;
			for(int i = 0; i < 3; i++){
				if(edges[i]) continue;
				edges[i] = points[i]->getEdgeToPoint(points[(i+1)%3]);
				if(edges[i] == nullptr){
					DataVector<MissingFace> missing_faces_local; // empty
					MissingEdge me(nullptr, points[i], points[(i + 1) % 3]);
					recoverBoundaryEdge(mc, mesh, me, missing_faces_local, phase);
					edges[i] = points[i]->getEdgeToPoint(points[(i+1)%3]);
					if(edges[i] == nullptr){
#ifdef LOG_RECOVERY_PIPES
						LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Points: " << points[i]->getID() << ", ";
						LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, points[(i+1)%3]->getID() << ", ";
						LOG4CPLUS_INFO(MeshLog::logger_mesh, points[(i+2)%3]->getID());
						LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Lengths: " << points[i]->getMetricCoordinates(mc).distance(
							points[(i+1)%3]->getMetricCoordinates(mc)) << ", ";
						LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, points[(i+1)%3]->getMetricCoordinates(mc).distance(
							points[(i+2)%3]->getMetricCoordinates(mc)) << ", ";
						LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, points[(i+2)%3]->getMetricCoordinates(mc).distance(
							points[i]->getMetricCoordinates(mc)) << endl;
						//MeshView::showDebugMesh("still missing", mesh, points[i], points[(i+1)%3], 1.5);
#endif
						++missing_count;
					}
				}
			}
		}
#ifdef LOG_RECOVERY_PIPES
		if(missing_count == 0){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " - missing edges successfully recovered. ");
			if(pt1->getFaceToPoints(pt2, pt3) != nullptr)
				return REDUCE_POINT_CLEAN;
		}
#endif
	}
	if(missing_count > 0 && phase == 1){ // try inserting points in the vicinity
		// near faces
		// - find
		DataVector<MeshFace*> near_faces;
		for(int i = 0; i < 3; i++){
			int ect = points[i]->getRank();
			for(int j = 0; j < ect; j++){
				MeshEdge3d* edge = points[i]->getEdge(j);
				int fct = edge->getFaceCount();
				for(int k = 0; k < fct; k++){
					MeshFace* face = edge->getFaceAt(k);
					if(face->isBorder()) near_faces.addIfNew(face);
				}
			}
		}
		// - insert points
		const DBox box = mesh->getBoundingBox();
		int inserted = 0;
		for(int i = 0; i < near_faces.countInt(); i++){
			int inserted_temp = 0;
			MeshFace* face = near_faces[i];
			mc.countMetricAtPoint(face->getMiddlePoint());
			MeshPoint3d* fpts[] = {face->getPoint(0), face->getPoint(1), face->getPoint(2)};
			const DMPoint3d fdmiddle = DMPoint3d::average(
				fpts[0]->getMetricCoordinates(mc),
				fpts[1]->getMetricCoordinates(mc), 
				fpts[2]->getMetricCoordinates(mc));
			DMVector3d nv = DMVector3d::crossProduct(
				fpts[0]->getMetricCoordinates(mc),
				fpts[1]->getMetricCoordinates(mc), 
				fpts[2]->getMetricCoordinates(mc));
			nv *= (MeshTetrahedron::opt_metric_height / nv.length());
			DPoint3d pts[2] = { 
				mc.transformMStoRS(fdmiddle + nv),
				mc.transformMStoRS(fdmiddle - nv)};
			const DPoint3d fmiddle = mc.transformMStoRS(fdmiddle);
			// insert (with special metric)
			for(int j = 0; j < 2; j++){
				if(!box.contains(pts[j])) 
					pts[j] = DPoint3d(fmiddle, box.fitInPoint(pts[j]), 0.75);
				inserted_temp += insertPointForFaceRecovery(mc, mesh, pts[j], missing_faces,
					(MeshTetrahedron*)face->getBlock(0));
			}
			inserted += inserted_temp;
			if(inserted_temp > 0){
				for(int j = 0; j < 3; j++)
					if(!edges[j]) edges[j] = points[j]->getEdgeToPoint(points[(j+1)%3]);
				if(edges[0] && edges[1] && edges[2]) break;
			}
		}
		// try again recovering edges
		//mc.countMetricAtPoint(pt1->getCoordinates());
		if(inserted > 0){
			missing_count = 0;
			for(int i = 0; i < 3; i++){
				if(edges[i]) continue;
				edges[i] = points[i]->getEdgeToPoint(points[(i+1)%3]);
				if(edges[i]) continue;
				DataVector<MissingFace> missing_faces_empty; // empty
				MissingEdge me(nullptr, points[i], points[(i + 1) % 3]);
				recoverBoundaryEdge(mc, mesh, me, missing_faces_empty, phase);
				edges[i] = points[i]->getEdgeToPoint(points[(i+1)%3]);
				if(edges[i] == nullptr) ++missing_count;
			}
		}
#ifdef LOG_RECOVERY_PIPES
		if(missing_count == 0){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " - missing edges successfully recovered. ");
			if(pt1->getFaceToPoints(pt2, pt3) != nullptr)
				return REDUCE_POINT_CLEAN;
		}
#endif
	}

	if(missing_count > 0 && phase == 2){
		int total_inserted = 0;
		for(int i = 0; i < 3; i++){
			int inserted;
			do{
				inserted = 0;
				int rank = points[i]->getRank();
				mc.countMetricAtPoint(points[i]->getCoordinates());
				const DMPoint3d dpt = points[i]->getMetricCoordinates(mc);
				for(int j = 0; j < rank; j++){
					const MeshEdge3d* edge = points[i]->getEdge(j);
					if(edge->isBorder()) continue;
					MeshPoint3d* other_point = edge->getOtherPoint(points[i]);
					double dist2 = dpt.distance2(other_point->getMetricCoordinates(mc));
					if(dist2 > 4.0){
						const DMPoint3d ndpt = dpt + (other_point->getMetricCoordinates(mc)-dpt).normalized();
						const DPoint3d npt = mc.transformMStoRS(ndpt);
						// check
						int mct = missing_faces.countInt();
						bool safe = true;
						for(int k = 0; safe && k < mct; k++){
							const FaceNode node(
								missing_faces[k].m_pt[0]->getCoordinates(), 
								missing_faces[k].m_pt[1]->getCoordinates(),
								missing_faces[k].m_pt[2]->getCoordinates());
							if(npt.distance2(node.middle) > node.max_dist2) continue;
							if(DTriangle3d::distance2ToPoint(npt, node.pt, node.v0, node.v1) < node.min_dist2){
								//LOG4CPLUS_WARN(MeshLog::logger_console, "insertPointForUnitEdge - too near some missing face???");
								safe = false;
								break;
							}
						}
						if(safe){
							inserted += insertPointForFaceRecovery(mc, mesh, npt, missing_faces,
								(MeshTetrahedron*)edge->getFaceAt(0)->getBlock(0));
							total_inserted += inserted;
							if(inserted) break;
						}
					}
				}
			}while(inserted > 0);
		}
		if(total_inserted > 0){
			missing_count = 0;
			for(int i = 0; i < 3; i++){
				if(edges[i]) continue;
				edges[i] = points[i]->getEdgeToPoint(points[(i+1)%3]);
				if(edges[i]) continue;
				DataVector<MissingFace> missing_faces_empty; // empty
				MissingEdge me(nullptr, points[i], points[(i + 1) % 3]);
				recoverBoundaryEdge(mc, mesh, me, missing_faces_empty, phase);
				edges[i] = points[i]->getEdgeToPoint(points[(i+1)%3]);
				if(edges[i] == nullptr) ++missing_count;
			}
		}
	}

	if(missing_count > 0){
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " - " << missing_count << " edges still missing. ");
#endif
			return REDUCE_FAILED;
	}

	if(pt1->getFaceToPoints(pt2, pt3) != nullptr)
		return REDUCE_TRANSFORM;

	while(true){
		DataVector<MeshEdge3d*> cross_edges;
		gatherCrossingEdgesForFace(mc, edges, points, cross_edges);
#ifdef LOG_RECOVERY_PIPES
		for(int i = 0; i < cross_edges.countInt(); i++){
			MeshEdge3d* cedge = cross_edges[i];
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " found crossing edge: ";
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - " << (cedge->isBorder()?'b':'i') << '[';
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, (cedge->getMeshPoint(0)->isBorder()?'b':'i') 
				<< cedge->getMeshPoint(0)->getID() << ',';
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, (cedge->getMeshPoint(1)->isBorder()?'b':'i') 
				<< cedge->getMeshPoint(1)->getID() << ']' << endl;
		}
#endif

		// draw tetrahedrons and 3 points
//----------------------------------------------
/*
		if(phase > 1){
			MeshViewSet *set = new MeshViewSet;			
			for(int i = 0; i < cross_edges.countInt(); i++){
				set->addEdge(cross_edges[i]);
			}
			set->addFace(pt1->getCoordinates(), 
				pt2->getCoordinates(), pt3->getCoordinates());
			SHOW_MESH("face recovery", set);
		}
*/
//----------------------------------------------

		int count = cross_edges.countInt();
		bool changed = false;
		for(int i = 0; i < count; i++){
			MeshEdge3d* edge = cross_edges.get(i);
#ifdef LOG_RECOVERY_PIPES
			assert(mesh->isValid());
			MeshPoint3d* pt4 = edge->getMeshPoint(0);
			MeshPoint3d* pt5 = edge->getMeshPoint(1);
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - edge[" << (pt4->isBorder()?'b':'i');
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, pt4->getID() << ',';
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, (pt5->isBorder()?'b':'i');
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, pt5->getID() << "] -> swap 3x2 ";
#endif
			// try swap
			if(MeshGenerator3d::swap32(mc, mesh, edge)){
#ifdef LOG_RECOVERY_PIPES
				assert(mesh->isValid());
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "...ok.");
#endif
				changed = true;
			}else{
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "...failed.");
#endif
			}
		}

		if(!changed) return REDUCE_FAILED;
		
		if(pt1->getFaceToPoints(pt2, pt3) != nullptr)
			return REDUCE_TRANSFORM;
	}

	//return REDUCE_FAILED;
}

MeshPoint3d * MeshGenerator3dDelaunayBoundary::cutFirstPipeNode(
	Metric3dContext & mc, MeshContainer3d * mesh, 
	MeshPoint3d * pt1, MeshPoint3d * pt2,
	PipeNode* pipe)
{
	assert(pipe != nullptr);
	switch (pipe->element_type) {
	case PIPE_FACE_POINT: 
	case PIPE_BLOCK_POINT: 
		return (MeshPoint3d*)(pipe->element2);
	case PIPE_FACE_EDGE: 
	{
		MeshEdge3d* end_edge = (MeshEdge3d*)(pipe->element2);
		assert(!end_edge->isBorder());
		MeshPoint3d* cross_pt = MeshGenerator3d::splitEdgeSimple(mesh, end_edge, pipe->pt);
		return cross_pt;
	}
	case PIPE_BLOCK_EDGE:
	{
		MeshEdge3d* end_edge = (MeshEdge3d*)(pipe->element2);
		assert(!end_edge->isBorder());
		MeshPoint3d* cross_pt = MeshGenerator3d::splitEdgeSimple(mesh, end_edge, pipe->pt);
		return cross_pt;
	}
	case PIPE_BLOCK_FACE:
	{
		MeshFace* end_face = (MeshFace*)(pipe->element2);
		assert(!end_face->isBorder());
		MeshPoint3d* cross_pt = MeshGenerator3d::splitFaceSimple(mesh, end_face, pipe->pt);
		return cross_pt;
	}
	}
	return nullptr;
}

void MeshGenerator3dDelaunayBoundary::gatherCrossingEdgesForFace(Metric3dContext& mc, MeshEdge3d* edges[], 
		MeshPoint3d* points[], DataVector<MeshEdge3d*> & cross_edges)
{
	// set metric (just in case)
	mc.countMetricAtPoint(DPoint3d::average(
		points[0]->getCoordinates(),
		points[1]->getCoordinates(),
		points[2]->getCoordinates()));

	// directly
	for(int m = 0; m < 3; m++){
		int face_ct = edges[m]->getFaceCount();
		for(int i = 0; i < face_ct; i++){
			MeshFace* face = edges[m]->getFaceAt(i);
			for(int j = 0; j < 2; j++){
				MeshTetrahedron* block = (MeshTetrahedron*)face->getBlock(j);
				if(block == nullptr) continue;
				MeshPoint3d* pts[2];
				block->getOtherPoints(points[m], points[(m+1)%3], pts);
				//
				if(DMTriangle3d::crossSegment(
					pts[0]->getMetricCoordinates(mc), 
					pts[1]->getMetricCoordinates(mc),
					points[0]->getMetricCoordinates(mc), 
					points[1]->getMetricCoordinates(mc), 
					points[2]->getMetricCoordinates(mc)))
				{
					MeshEdge3d* cedge = pts[0]->getEdgeToPoint(pts[1]);
					cross_edges.addIfNew(cedge);
					//cross_blocks.addIfNew(block);
				}
			}
		}
	}

	// through faces incident to crossing edges
	for(int i = 0; i < cross_edges.countInt(); i++){
		MeshEdge3d* edge = cross_edges.get(i);
		MeshPoint3d* pt4 = edge->getMeshPoint(0);
		MeshPoint3d* pt5 = edge->getMeshPoint(1);
		int face_ct = edge->getFaceCount();
		for(int j = 0; j < face_ct; j++){
			MeshFace* face = edge->getFaceAt(j);
			MeshPoint3d* pt6 = face->getOtherPoint(pt4, pt5);
			if(pt6 == points[0] || pt6 == points[1] || pt6 == points[2])
				continue;
			MeshEdge3d* edge1 = pt4->getEdgeToPoint(pt6);
			if((!cross_edges.contains(edge1)) &&
				DMTriangle3d::crossSegment(
					pt4->getMetricCoordinates(mc), 
					pt6->getMetricCoordinates(mc), 
					points[0]->getMetricCoordinates(mc), 
					points[1]->getMetricCoordinates(mc), 
					points[2]->getMetricCoordinates(mc)))
			{
				cross_edges.add(edge1);
			}
			edge1 = pt5->getEdgeToPoint(pt6);
			if((!cross_edges.contains(edge1)) &&
				DMTriangle3d::crossSegment(
					pt5->getMetricCoordinates(mc), 
					pt6->getMetricCoordinates(mc), 
					points[0]->getMetricCoordinates(mc), 
					points[1]->getMetricCoordinates(mc), 
					points[2]->getMetricCoordinates(mc)))
			{
				cross_edges.add(edge1);
			}
		}
	}

	// through edges incident to vertices of crossing edges
	DataVector<MeshPoint3d*> vertices(50);
	vertices.add(points[0]);
	vertices.add(points[1]);
	vertices.add(points[2]);
	for(int i = 0; i < cross_edges.countInt(); i++){
		MeshEdge3d* edge = cross_edges.get(i);
		vertices.addIfNew(edge->getMeshPoint(0));
		vertices.addIfNew(edge->getMeshPoint(1));
	}
	for(int i = 3; i < vertices.countInt(); i++){
		MeshPoint3d* pt = vertices.get(i);
		int rank = pt->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = pt->getEdge(j);
			if(!cross_edges.contains(edge)){
				MeshPoint3d* other_pt = edge->getOtherPoint(pt);
				if(other_pt == points[0] || other_pt == points[1] || other_pt == points[2])
					continue;
				if(DMTriangle3d::crossSegment(
					pt->getMetricCoordinates(mc), 
					other_pt->getMetricCoordinates(mc),	
					points[0]->getMetricCoordinates(mc), 
					points[1]->getMetricCoordinates(mc), 
					points[2]->getMetricCoordinates(mc)))
				{
					cross_edges.add(edge);
					vertices.addIfNew(other_pt);
				}
			}
		}
	}

	// end
}

MeshGenerator3dDelaunayBoundary::PipeNode* MeshGenerator3dDelaunayBoundary::nextPipeFromPoint(Metric3dContext& mc, 
		MeshPoint3d *last_point, MeshPoint3d *end_point)
{
	int count = last_point->getRank();
	DataHeapVector<CandidateEdge> candidate_edges(count);
	//mc.countMetricAtPoint(last_point->getCoordinates());
	const DMPoint3d dpt0 = last_point->getMetricCoordinates(mc);
	const DMPoint3d dpt1 =  end_point->getMetricCoordinates(mc);

	//evaluate all edges incident to this node
	const DMVector3d dv01 = (dpt1 - dpt0).normalized();
	for(int i = 0; i < count; i++){
		CandidateEdge ce;
		ce.edge = last_point->getEdge(i);
		const DMPoint3d dpt2 = ce.edge->getOtherPoint(last_point)->getMetricCoordinates(mc);
		const DMVector3d dv02 = (dpt2 - dpt0).normalized();
		ce.direction = dv01.scalarProduct(dv02);
		candidate_edges.add(ce);
	}

	DataVector<MeshBlock*> blocks(100);
	DataVector<MeshFace*>  faces(100);

	while(candidate_edges.countInt() > 0){
		// select best edge (heap -> best is always on top)
		MeshEdge3d* edge = candidate_edges.remove().edge;
		int fct = edge->getFaceCount();
		for(int j = 0; j < fct; j++){
			MeshFace* face = edge->getFaceAt(j);
			if(faces.addIfNew(face)){
				for(int k = 0; k < 2; k++){
					MeshTetrahedron* block = (MeshTetrahedron*)face->getBlock(k);
					if(block && blocks.addIfNew(block)){
						MeshFace* cross_face = block->getOppositeFace(last_point);
						const DMPoint3d fdpts[3] = {
							cross_face->getPoint(0)->getMetricCoordinates(mc),
							cross_face->getPoint(1)->getMetricCoordinates(mc),
							cross_face->getPoint(2)->getMetricCoordinates(mc)
						};
						// check volumes (whether plane crosses segment)
						double vol0 = DMTriangle3d::orient3d(fdpts[0], fdpts[1], fdpts[2], dpt0);
						double vol1 = DMTriangle3d::orient3d(fdpts[0], fdpts[1], fdpts[2], dpt1);
						if(vol0 * vol1 > 0.0) continue; // if same side -> no crossing...
						// check volumes (whether line crosses face)
						double dvol[3];
						for(int m = 0; m < 3; m++){
							int mnext = (m+1)%3;
							dvol[m] = DMTriangle3d::orient3d(dpt0, dpt1, fdpts[m], fdpts[mnext]);
							if(abs(dvol[m]) < param_recovery_volume_threshold){
								// means 'missing edge' and some face are at the same plane -> check for crossing
								// create additional point above the plane for volume calculation
								const DMPoint3d dpt4 = dpt0 + DMVector3d::crossProduct(dpt0, fdpts[m], fdpts[mnext]);
								double vold0 = DMTriangle3d::orient3d(dpt0, fdpts[m], dpt1, dpt4);
								double vold1 = DMTriangle3d::orient3d(dpt0, dpt1, fdpts[mnext], dpt4);
								if(vold0 * vold1 > 0.0){
									// pipe: --face-->edge
									double s = vold0 / (vold0+vold1);
									const DPoint3d cross_point(
										cross_face->getPoint(m)->getCoordinates(),
										cross_face->getPoint(mnext)->getCoordinates(), 
										s);
									MeshEdge3d* cross_edge = cross_face->getEdge(m);
									MeshFace* through_face = cross_edge->getFaceToPoint(last_point);
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
									if(false){
//									if(end_point->getID() == 533){
										MeshViewSet *set = new MeshViewSet;			
										set->addPoint(last_point, 1);
										set->addPoint(end_point, 1);
										set->addPoint(cross_point, 2);
										set->addFace(through_face, -2, 0.95);
										set->addEdge(cross_edge, 0);
										SHOW_MESH("pipe from point through face to edge", set);
									}
#endif						//---------------------------
									return new PipeNode(PIPE_FACE_EDGE, through_face, cross_edge, cross_point);
								}else{
									// no crossing, skip 
								}
							}
						}
						if((dvol[0] > 0 && dvol[1] > 0 && dvol[2] > 0) ||
								(dvol[0] < 0 && dvol[1] < 0 && dvol[2] < 0)){
							// pipe: --block-->face
							DPoint3d cross_point;
							cross_point.add(cross_face->getPoint(0)->getCoordinates(), dvol[1]);
							cross_point.add(cross_face->getPoint(1)->getCoordinates(), dvol[2]);
							cross_point.add(cross_face->getPoint(2)->getCoordinates(), dvol[0]);
							cross_point /= (dvol[0] + dvol[1] + dvol[2]);
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
//								if(end_point->getID() == 66){
								if(false){
									MeshViewSet *set = new MeshViewSet;			
									set->addPoint(last_point, 1);
									set->addPoint(end_point, 1);
									set->addPoint(cross_point, 2);
									set->addFace(cross_face);
									LOG4CPLUS_INFO(MeshLog::logger_mesh, "vol= " << dvol[0] << " , " << dvol[1] << " , " << dvol[2]);
									SHOW_MESH("pipe from point through block to face", set);
								}
#endif						//---------------------------
							return new PipeNode(PIPE_BLOCK_FACE, block, cross_face, cross_point);
						}
					}
				}
			}
		}
	}

	LOG4CPLUS_WARN(MeshLog::logger_console, "Failed finding next pipe element (pipe from point).");
	return nullptr;
}

MeshGenerator3dDelaunayBoundary::PipeNode* MeshGenerator3dDelaunayBoundary::nextPipeFromEdge(Metric3dContext& mc, 
		MeshEdge3d *last_edge, const DPoint3d& cross_point, MeshPoint3d *end_point)
{
	int fcount = last_edge->getFaceCount();

	// Check faces (edge->point)
	MeshPoint3d* pt1 = last_edge->getMeshPoint(0);
	MeshPoint3d* pt2 = last_edge->getMeshPoint(1);
	//mc.countMetricAtPoint(cross_point);
	const DMPoint3d dpt1 = pt1->getMetricCoordinates(mc);
	const DMPoint3d dpt2 = pt2->getMetricCoordinates(mc);
	const DMVector3d dv12 = (dpt2 - dpt1).normalized();
	const DMPoint3d dptc = mc.transformRStoMS(cross_point);
	const DMPoint3d dpte = end_point->getMetricCoordinates(mc);
	const DMVector3d dvce = (dpte - dptc).normalized();

	DataHeapVector<CandidateFace> candidate_faces(fcount);

	for(int i = 0; i < fcount; i++){
		CandidateFace cf;
		cf.face = last_edge->getFaceAt(i);
		MeshPoint3d* pt3 = cf.face->getOtherPoint(pt1, pt2);
		// if the face contains the ending point, return it as a pipe-node
		if(pt3 == end_point)
			return new PipeNode(PIPE_FACE_POINT, cf.face, end_point, end_point->getCoordinates());
		// else check orientation and insert into candidate faces
		const DMVector3d dvc3 = (pt3->getMetricCoordinates(mc) - dptc).normalized();
		cf.direction = dvce.scalarProduct(dvc3);
		candidate_faces.add(cf);
	}

	DataVector<MeshBlock*> blocks(fcount);

	// Check faces
	while(candidate_faces.countInt() > 0){
		MeshFace* face = candidate_faces.remove().face;
		blocks.addIfNew(face->getBlock(0));
		blocks.addIfNew(face->getBlock(1));
		MeshPoint3d* pt3 = face->getOtherPoint(pt1, pt2);
		const DMPoint3d dpt3 = pt3->getMetricCoordinates(mc);
		double vol_b0 = DMTriangle3d::orient3d(dpt1, dpt2, dpte, dpt3);
		if(abs(vol_b0) < param_recovery_volume_threshold){
			// --face-->edge
			const DMPoint3d dpt4 = dpt1 + dv12.crossProduct(dpt3-dpt1).crossProduct(dv12); // make it orthogonal
			SurfacePlane plane(
				ControlDataMatrix3d::identity.multiplyMStoRS(dpt1), 
				ControlDataMatrix3d::identity.multiplyMStoRS(dpt2), 
				ControlDataMatrix3d::identity.multiplyMStoRS(dpt4));
			const DPoint2d pt2d_1(0.0, 0.0); // = plane.getParameters(dpt1);
			const DPoint2d pt2d_2(1.0, 0.0); // = plane.getParameters(dpt2);
			const DPoint2d pt2d_3 = plane.getParameters(
				ControlDataMatrix3d::identity.multiplyMStoRS(dpt3));
			const DPoint2d pt2d_e = plane.getParameters(
				ControlDataMatrix3d::identity.multiplyMStoRS(dpte));
			const DPoint2d pt2d_c = plane.getParameters(
				ControlDataMatrix3d::identity.multiplyMStoRS(dptc));
			// same side ?
			double vol0 = DTriangle2d::det(pt2d_1, pt2d_2, pt2d_e);
			double vol1 = DTriangle2d::det(pt2d_1, pt2d_2, pt2d_3);
			if(vol0 * vol1 <= 0.0) continue; // if no, skip this face
			// check first edge (pt1--pt3)
			vol0 = DTriangle2d::det(pt2d_1, pt2d_c, pt2d_e);
			vol1 = DTriangle2d::det(pt2d_c, pt2d_3, pt2d_e);
			if(vol0 * vol1 > 0.0){
				double s = vol0 / (vol0+vol1);
				const DPoint3d new_cross_point(pt1->getCoordinates(), pt3->getCoordinates(), s);
				MeshEdge3d* new_cross_edge = pt1->getEdgeToPoint(pt3);
				assert(new_cross_edge);
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
									if(false){
//									if(end_point->getID() == 533){
										MeshViewSet *set = new MeshViewSet;
										set->addPoint(cross_point, 1);
										set->addPoint(end_point, 1);
										set->addFace(face);
										set->addEdge(new_cross_edge);
										set->addPoint(new_cross_point, 2);
										SHOW_MESH("pipe from edge through face to edge (1)", set);
									}
#endif						//---------------------------
				return new PipeNode(PIPE_FACE_EDGE, face, new_cross_edge, new_cross_point);
			}
			// check second edge (pt2--pt3)
			vol0 = DTriangle2d::det(pt2d_c, pt2d_e, pt2d_2);
			if((vol0 > 0 && vol1 > 0) || (vol0 < 0 && vol1 < 0)){
				double s = vol0 / (vol0+vol1);
				const DPoint3d new_cross_point(pt2->getCoordinates(), pt3->getCoordinates(), s);
				MeshEdge3d* new_cross_edge = pt2->getEdgeToPoint(pt3);
				assert(new_cross_edge);
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
									if(false){
//									if(end_point->getID() == 533){
										MeshViewSet *set = new MeshViewSet;
										set->addPoint(cross_point, 1);
										set->addPoint(end_point, 1);
										set->addFace(face);
										set->addEdge(new_cross_edge);
										set->addPoint(new_cross_point, 2);
										SHOW_MESH("pipe from edge through face to edge (2)", set);
									}
#endif						//---------------------------
				return new PipeNode(PIPE_FACE_EDGE, face, new_cross_edge, new_cross_point);
			} // else the missing edge is at the same plane as face, but it doesn't cross it
		}
	}
	// else check blocks
	for(int i = 0; i < blocks.countInt(); i++){
		MeshTetrahedron* block = (MeshTetrahedron*)blocks[i];
		MeshPoint3d *pts34[2];
		block->getOtherPoints(pt1, pt2, pts34);
		const DMPoint3d dpt3 = pts34[0]->getMetricCoordinates(mc);
		const DMPoint3d dpt4 = pts34[1]->getMetricCoordinates(mc);
		double vol0 = DMTriangle3d::orient3d(dptc, dpte, dpt3, dpt4);
		if(abs(vol0) < param_recovery_volume_threshold){
			// segments c.e and 3.4 are co-planar, but do they cross ?
//			const DMPoint3d dptx = dptc + DTriangle3d::crossProduct(dptc, dpte, dpt3);
//			double xvol0 = DTriangle3d::orient3d(dpt3, dpt4, dptc, dptx);
//			double xvol1 = DTriangle3d::orient3d(dpt3, dpt4, dpte, dptx);
//			if(xvol0 * xvol1 < 0.0){ //ok, check second pair
//				xvol0 = DTriangle3d::orient3d(dptc, dpte, dpt3, dptx);
//				xvol1 = DTriangle3d::orient3d(dptc, dpte, dpt4, dptx);
//			}
			double xvol0 = DMTriangle3d::orient3d(dpt1, dpt2, dpt3, dpte);
			double xvol1 = DMTriangle3d::orient3d(dpt1, dpt2, dpt3, dpt4);
			if(xvol0 * xvol1 > 0.0){ //ok, check second pair
				xvol0 = DMTriangle3d::orient3d(dpt1, dpt2, dpt4, dpte);
				// xvol1 = DTriangle3d::orient3d(dpt1, dpt2, dpt4, dpt3);  // so it's sufficient to change sign
				xvol1 = -xvol1;
			}
			if(xvol0 * xvol1 > 0.0){ // ok
				// pipe: --block-->edge
				double vol_b0 = DMTriangle3d::orient3d(dpt1, dpt2, dpte, dpt3);
				double vol_b1 = DMTriangle3d::orient3d(dpt1, dpt2, dpt4, dpte);

//				if(vol_b0 * vol_b1 <= 0.0){
//					MeshViewSet *set = new MeshViewSet;
//					set->addEmptyBlockWithEdges(block, 0);
//					set->addPoint(cross_point, 1);
//					set->addPoint(end_point, 1);
//					set->addEdge(cross_point, end_point->getCoordinates(), 1);
//					SHOW_MESH("pipe from edge through block to edge", set);
//				}

				assert(vol_b0 * vol_b1 > 0.0);
				double s = vol_b0 / (vol_b0+vol_b1);
				const DPoint3d new_cross_point(pts34[0]->getCoordinates(), pts34[1]->getCoordinates(), s);
				MeshEdge3d* new_cross_edge = pts34[0]->getEdgeToPoint(pts34[1]);
				assert(new_cross_edge);
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
								if(false){
									MeshViewSet *set = new MeshViewSet;
									set->addPoint(cross_point, 1);
									set->addPoint(end_point, 1);
									//set->addFace(face);
									set->addEdge(new_cross_edge);
									set->addPoint(new_cross_point, 2);
									SHOW_MESH("pipe from edge through block to edge", set);
								}
#endif						//---------------------------
				return new PipeNode(PIPE_BLOCK_EDGE, block, new_cross_edge, new_cross_point);
			}
		}
		double vol1 = DMTriangle3d::orient3d(dptc, dpte, dpt4, dpt2);
		if(vol0 * vol1 > 0.0){
			// check third vol
			double vol2 = DMTriangle3d::orient3d(dptc, dpte, dpt2, dpt3);
			if(vol0 * vol2 > 0.0){
				// check orientation
				double volx0 = DMTriangle3d::orient3d(dpt2, dpt3, dpt4, dptc);
				double volx1 = DMTriangle3d::orient3d(dpt2, dpt3, dpt4, dpte);
				if(volx0 * volx1 < 0.0){ // if both points are at different sides
					MeshFace* new_cross_face = block->getOppositeFace(pt1);
					DPoint3d new_cross_point;
					new_cross_point.add(pt2->getCoordinates(), vol0);
					new_cross_point.add(pts34[0]->getCoordinates(), vol1);
					new_cross_point.add(pts34[1]->getCoordinates(), vol2);
					new_cross_point /= (vol0 + vol1 + vol2);
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
								if(false){
//								if(end_point->getID() == 533){
									MeshViewSet *set = new MeshViewSet;
									set->addPoint(cross_point, 1);
									set->addPoint(end_point, 1);
									set->addFace(new_cross_face);
									set->addPoint(new_cross_point, 2);
									SHOW_MESH("pipe from edge through block to face", set);
								}
#endif						//---------------------------
					return new PipeNode(PIPE_BLOCK_FACE, block, new_cross_face, new_cross_point);
				}
			}
		}
		// check second face...
		vol1 = DMTriangle3d::orient3d(dptc, dpte, dpt4, dpt1);
		if(vol0 * vol1 > 0.0){
			// check third vol
			double vol2 = DMTriangle3d::orient3d(dptc, dpte, dpt1, dpt3);
			if(vol0 * vol2 > 0.0){
				double volx0 = DMTriangle3d::orient3d(dpt1, dpt3, dpt4, dptc);
				double volx1 = DMTriangle3d::orient3d(dpt1, dpt3, dpt4, dpte);
				if(volx0 * volx1 < 0.0){ // if both points are at different sides
					MeshFace* new_cross_face = block->getOppositeFace(pt2);
					DPoint3d new_cross_point;
					new_cross_point.add(pt1->getCoordinates(), vol0);
					new_cross_point.add(pts34[0]->getCoordinates(), vol1);
					new_cross_point.add(pts34[1]->getCoordinates(), vol2);
					new_cross_point /= (vol0 + vol1 + vol2);
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
							if(false){
//							if(end_point->getID() == 533){
									MeshViewSet *set = new MeshViewSet;
									set->addPoint(cross_point, 1);
									set->addPoint(end_point, 1);
									set->addFace(new_cross_face);
									set->addPoint(new_cross_point, 2);
									LOG4CPLUS_INFO(MeshLog::logger_mesh, "vol= " << vol0 << " , " << vol1 << " , " << vol2);
									SHOW_MESH("pipe from edge through block to face", set);
							}
#endif						//---------------------------
					return new PipeNode(PIPE_BLOCK_FACE, block, new_cross_face, new_cross_point);
				}
			}
		}
	}

#ifdef _DEBUG
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
		MeshViewSet *set = new MeshViewSet;
		set->addPoint(cross_point, 1);
		set->addPoint(end_point, 1);
		set->addEdge(cross_point, end_point->getCoordinates(), 1);
		set->addPoint(pt1, 2);
		set->addPoint(pt2, 2);
		int rank = last_edge->getFaceCount();
		for(int i = 0; i < rank; i++){
			MeshFace* mface = last_edge->getFaceAt(i);
			set->addFace(mface);
			set->addPoint(mface->getOtherPoint(pt1, pt2));
		}
		SHOW_MESH("Failed finding next pipe element (pipe from edge)", set);
#endif						//---------------------------
#endif //_DEBUG

	LOG4CPLUS_WARN(MeshLog::logger_console, "Failed finding next pipe element (pipe from edge).");
	return nullptr;
}

MeshGenerator3dDelaunayBoundary::PipeNode* MeshGenerator3dDelaunayBoundary::nextPipeFromFace(
	Metric3dContext& mc, MeshBlock* last_block, 
	MeshFace *last_face, const DPoint3d &cross_point, MeshPoint3d *end_point)
{
	MeshTetrahedron* tetrahedron = (MeshTetrahedron*)last_face->getOtherBlock(last_block);
	assert(tetrahedron);
	MeshPoint3d* pts[3] = { last_face->getPoint(0), last_face->getPoint(1), last_face->getPoint(2)};
	MeshPoint3d* pt3 = tetrahedron->getOppositePoint(last_face);
	if(pt3 == end_point)
		return new PipeNode(PIPE_BLOCK_POINT, tetrahedron, end_point, end_point->getCoordinates());			

	//mc.countMetricAtPoint(cross_point);

	const DMPoint3d dpts[] = { pts[0]->getMetricCoordinates(mc), pts[1]->getMetricCoordinates(mc),
		pts[2]->getMetricCoordinates(mc) };
	const DMPoint3d  dpt3 = pt3->getMetricCoordinates(mc);
	const DMPoint3d  dptc = mc.transformRStoMS(cross_point);
	const DMPoint3d  dpte = end_point->getMetricCoordinates(mc);

	// Check blocks (--block-->edge)
	double vols[3];
	for(int i = 0; i < 3; i++){
		vols[i] = DMTriangle3d::orient3d(dptc, dpte, dpts[i], dpt3);
		if(abs(vols[i]) < param_recovery_volume_threshold){
			int in = (i+1)%3;
			double vol0 = DMTriangle3d::orient3d(dptc, dpte, dpts[i], dpts[in]);
			double vol1 = DMTriangle3d::orient3d(dptc, dpt3, dpte, dpts[in]);
			if(vol0 * vol1 > 0.0){
				double s = vol0 / (vol0+vol1);
				const DPoint3d new_cross_point(pts[i]->getCoordinates(), pt3->getCoordinates(), s);
				MeshEdge3d* new_cross_edge = pts[i]->getEdgeToPoint(pt3);
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
							if(false){
									MeshViewSet *set = new MeshViewSet;
									set->addPoint(cross_point, 1);
									set->addPoint(end_point, 1);
									set->addEdge(new_cross_edge);
									set->addBlock(tetrahedron);
									set->addPoint(new_cross_point, 2);
									LOG4CPLUS_INFO(MeshLog::logger_mesh, "vol = " << tetrahedron->getVolume(mc, false));
									SHOW_MESH("pipe from face through block to edge", set);
							}
#endif						//---------------------------
				return new PipeNode(PIPE_BLOCK_EDGE, tetrahedron, new_cross_edge, new_cross_point);	
			}
		}
	}

	// Check blocks (--block-->face)
	for(int i = 0; i < 3; i++){
		int in = (i+1)%3;
		if(abs(vols[i]) < param_recovery_volume_threshold ||
			abs(vols[in]) < param_recovery_volume_threshold) continue;
		if(vols[i] * vols[in] > 0.0) continue;
		double vol = DMTriangle3d::orient3d(dptc, dpte, dpts[in], dpts[i]);
		if(vols[i] * vol < 0.0) continue;
		MeshFace* face = tetrahedron->getOppositeFace(pts[(in+1)%3]);
		assert(face);
		DPoint3d new_cross_point;
		new_cross_point.add(pts[i]->getCoordinates(), -vols[in]);
		new_cross_point.add(pt3->getCoordinates(),     vol);
		new_cross_point.add(pts[in]->getCoordinates(), vols[i]);
		new_cross_point /= (-vols[in] + vol + vols[i]);
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
//								if(end_point->getID() == 521){
								if(false){
									MeshViewSet *set = new MeshViewSet;
									set->addPoint(cross_point, 1);
									set->addPoint(end_point, 1);
									set->addFace(face);
									set->addPoint(new_cross_point, 2);
									LOG4CPLUS_INFO(MeshLog::logger_mesh, "vol= " << -vols[in] << " , " << vol << " , " << vols[i]);
									SHOW_MESH("pipe from face through block to face", set);
								}
#endif						//---------------------------
		return new PipeNode(PIPE_BLOCK_FACE, tetrahedron, face, new_cross_point);
	}

	LOG4CPLUS_WARN(MeshLog::logger_console, "Failed finding next pipe element (pipe from face).");
/*
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "vols = [" << vols[0] << ',' << vols[1] << ',' << vols[2] << ']');

	mc.countMetricAtPoint(cross_point);
	const DPoint3d mdpts[] = { pts[0]->getMetricCoordinates(mc), pts[1]->getMetricCoordinates(mc),
		pts[2]->getMetricCoordinates(mc) };
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "lens = [" << mdpts[0].distance(mdpts[1]) << ',' 
		<< mdpts[1].distance(mdpts[2]) << ',' 
		<< mdpts[2].distance(mdpts[0]) << ']' << endl;

	const DPoint3d& mdpt3 = pt3->getMetricCoordinates(mc);
	const DPoint3d  mdptc = mc.transformRStoMS(cross_point);
	const DPoint3d  mdpte = end_point->getMetricCoordinates(mc);
	for(int i = 0; i < 3; i++)
		vols[i] = DTriangle3d::orient3d(mdptc, mdpte, mdpts[i], mdpt3);
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "vols = [" << vols[0] << ',' << vols[1] << ',' << vols[2] << ']');
*/
//#ifdef _DEBUG
#ifdef DEBUG_RECOVERY_PIPES //---------------------------
		MeshViewSet *set = new MeshViewSet;
		set->addFace(last_face);
		set->addPoint(cross_point, 1);
		set->addPoint(end_point, 1);
		set->addEdge(cross_point, end_point->getCoordinates(), 1);
		set->addEmptyBlockWithEdges(tetrahedron, 2);
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Failed finding next pipe element (pipe from face): end_point->id == " 
			<< end_point->getID() << endl;
		SHOW_MESH("Failed finding next pipe element (pipe from face)", set);
#endif						//---------------------------
//#endif //_DEBUG
	return nullptr;
}

bool MeshGenerator3dDelaunayBoundary::reducePipeNodes(Metric3dContext& mc, MeshContainer3d* mesh, PipeNode* last_node, 
							PipeNode *node, int phase, bool& invalid_Pipe)
{
	invalid_Pipe = false;
	PipeNode* node2 = node->next;
	if(node->element_type == PIPE_FACE_EDGE &&
		node2->element_type == PIPE_FACE_POINT)
	{
		MeshEdge3d* common_edge = (MeshEdge3d*)(node->element2);
		if(common_edge->getFaceCount() > 4){
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - fixing for swap 4x4 ...";
#endif
			if(reduceEdgeRankForSwap44(mc, mesh, common_edge, (MeshFace*)node->element1, (MeshFace*)node2->element1)){
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " ok.");
#endif
			}else{
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " failed.");
#endif
			}
		}
		if(common_edge->getFaceCount() == 4){
			// Case 4->4
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - swap 4x4 ... ";
#endif
			if(MeshGenerator3d::swap44(mc, mesh, (MeshFace*)node->element1, (MeshEdge3d*)node->element2, phase)){
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " ok.");
#endif
				node->next = node2->next;
				node->element2 = node2->element2;
				if(last_node){
					node->element_type = PIPE_FACE_POINT;
					assert(last_node->element_type == PIPE_FACE_EDGE ||
						last_node->element_type == PIPE_BLOCK_EDGE);
					MeshEdge3d * edge = (MeshEdge3d*)last_node->element2;
					assert(edge);
					MeshFace* new_face = ((MeshPoint3d*)(node->element2))->getFaceToPoints(
						edge->getMeshPoint(0), edge->getMeshPoint(1));
					node->element1 = new_face;
				}else{
					node->element_type = PIPE_START_POINT;
				}
				delete node2;
				return true;
			}else{
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " failed.");
#endif
				return false;
			}
		}
	}else if (node->element_type == PIPE_FACE_EDGE &&
		node2->element_type == PIPE_FACE_EDGE)
	{
		MeshEdge3d* common_edge = (MeshEdge3d*)(node->element2);
		if(common_edge->getFaceCount() > 4){
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - fixing for swap 4x4 ...";
#endif
			if(reduceEdgeRankForSwap44(mc, mesh, common_edge, (MeshFace*)node->element1, (MeshFace*)node2->element1)){
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " ok.");
#endif
			}else{
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " failed.");
#endif
			}
		}
		if(common_edge->getFaceCount() == 4){
			// Case 4->4
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - swap 4x4 ... ";
#endif
			//if(phase == 4) LOG4CPLUS_INFO(MeshLog::logger_console, " edge-4 before 4x4 (2), mesh_valid", mesh->isValid());
			if(MeshGenerator3d::swap44(mc, mesh, (MeshFace*)node->element1, (MeshEdge3d*)node->element2, phase)){
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " ok.");
#endif
				//if(phase == 4) LOG4CPLUS_INFO(MeshLog::logger_console, " edge-4 after 4x4 (2), mesh_valid", mesh->isValid());
				invalid_Pipe = true;
				return true;
			}else{
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " failed.");
#endif
				return false;
			}
		}
	}else if(node->element_type == PIPE_BLOCK_FACE){
		if(last_node == nullptr){
			MeshFace* common_face = (MeshFace*)(node->element2);
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - swap 2x3 (point-face) ...";
#endif
			if(MeshGenerator3d::swap23(mc, mesh, common_face)){
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "ok.");
#endif
				node->next = node2->next;
				node->element2 = node2->element2;
				if(node->next == nullptr){
					node->element_type = PIPE_START_POINT;
					delete node2;
				}else if(node2->element_type == PIPE_BLOCK_FACE){
					node->element_type = PIPE_BLOCK_FACE;
					MeshFace* next_face = (MeshFace*)(node2->element2);
					assert(next_face->getType() == FACE_TRIANGLE);
					node->element1 = next_face->getOtherBlock((MeshBlock*)node->next->element1);
					delete node2;
				}else{
					invalid_Pipe = true; // -- to edge, but can be either through block or face...
				}
				return true;
			}else{
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " failed.");
#endif
				return false;
			}
		}else if(node2->element_type == PIPE_BLOCK_POINT){
			MeshFace* common_face = (MeshFace*)(node->element2);
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - swap 2x3 (face-point) ...";
#endif
			if(MeshGenerator3d::swap23(mc, mesh, common_face)){
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "ok.");
#endif
				if(last_node->element_type == PIPE_BLOCK_FACE){
					node->next = node2->next;
					node->element2 = node2->element2;
					node->element_type = PIPE_BLOCK_POINT;
					MeshFace* last_face = (MeshFace*)(last_node->element2);
					MeshBlock* last_block = (MeshBlock*)(last_node->element1);
					MeshBlock* this_block = last_face->getOtherBlock(last_block);
					assert(this_block);
					node->element1 = this_block;
					delete node2;
				}else{
					invalid_Pipe = true;
				}
				return true;
			}else{
#ifdef LOG_RECOVERY_PIPES
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " failed.");
#endif
				return false;
			}
		}
	}

	return false;
}

bool MeshGenerator3dDelaunayBoundary::transformPipeVicinity(Metric3dContext& mc, MeshContainer3d* mesh, PipeNode *pipe,
			DataVector<int> &old_faces, DataVector<int> &old_edges)
{
	DataVector<MeshFace*> swap23_faces;
	DataVector<MeshEdge3d*> swap32_edges;
	for(PipeNode* node = pipe; node; node=node->next){
		switch(node->element_type){
			case PIPE_BLOCK_EDGE:
			case PIPE_BLOCK_FACE:
			case PIPE_BLOCK_POINT:
				{
					MeshTetrahedron* t = (MeshTetrahedron*)node->element1;
					for(int i = 0; i < 4; i++){
						MeshFace* face = t->getFace(i);
						int p[3] = {face->getPoint(0)->getIndex(),
							face->getPoint(1)->getIndex(), face->getPoint(2)->getIndex()};
						int chk[3][2] = {{0,1},{1,2},{0,1}};
						for(int j = 0; j < 3; j++){
							int i0 = chk[j][0];
							int i1 = chk[j][1];
							if(p[i0] > p[i1]){
								int px = p[i0]; p[i0] = p[i1]; p[i1] = px;
							}
						}
						assert((p[0] < p[1]) && (p[1] < p[2]));
						int ct = old_faces.countInt() / 3;
						bool found = false;
						for(int j = 0; j < ct; j++){
							if(old_faces[3*j] == p[0] &&
									old_faces[3*j+1] == p[1] &&
									old_faces[3*j+2] == p[2]){
								found = true;
								break;
							}
						}
						if(!found){
							old_faces.add(p[0]);
							old_faces.add(p[1]);
							old_faces.add(p[2]);
							assert(!swap23_faces.contains(face));
							if(MeshGenerator3d::swap23possible(mc, face)) swap23_faces.add(face);
						}
					}
					for(int i = 0; i < 6; i++){
						MeshEdge3d* edge = t->getEdge(i);
						int p[2] = {edge->getMeshPoint(0)->getIndex(),
							edge->getMeshPoint(1)->getIndex()};
						if(p[0] > p[1]){
							int px = p[0]; p[0] = p[1]; p[1] = px;
						}

						int ct = old_edges.countInt() / 2;
						bool found = false;
						for(int j = 0; j < ct; j++){
							if(old_edges[2*j] == p[0] && old_edges[2*j+1] == p[1]){
								found = true; break;
							}
						}
						if(!found){
							old_edges.add(p[0]);
							old_edges.add(p[1]);
							assert(!swap32_edges.contains(edge));
							if(MeshGenerator3d::swap32possible(mc, edge)) swap32_edges.add(edge);
						}
					}
					break;
				}
			default:
				break;
		}
	}

	int fct = swap23_faces.countInt();
	int ect = swap32_edges.countInt();
#ifdef DEBUG_RECOVERY_PIPES
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - vicinity transformation: " << fct << " faces for 2x3, ";
	LOG4CPLUS_INFO(MeshLog::logger_mesh, ect << " edges for 3x2.");
#endif

	int total_count = fct + ect; 
	if(total_count == 0) return false;

	int i = rand() % total_count;

	if(i < ect){
#ifdef DEBUG_RECOVERY_PIPES
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " -- swap32 for edge nr " << i);
#endif
		bool result = MeshGenerator3d::swap32(mc, mesh, swap32_edges.get(i));
		if(!result){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   
				"transformPipeVicinity: swap32 should be possible, but isn't ???");
		}
		return result;
	}else{
		i -= ect;
#ifdef DEBUG_RECOVERY_PIPES
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " -- swap23 for face nr " << i);
#endif
		bool result = MeshGenerator3d::swap23(mc, mesh, swap23_faces.get(i));
		if(!result){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   
				"transformPipeVicinity: swap23 should be possible, but isn't ???");
		}
		return result;
	}
}

bool MeshGenerator3dDelaunayBoundary::transformPipeNodes(Metric3dContext& mc, MeshContainer3d* mesh, PipeNode *node)
{
	PipeNode* node2 = node->next;
	if(node->element_type == PIPE_BLOCK_FACE &&
		node2->element_type == PIPE_BLOCK_FACE){
		MeshFace* common_face = (MeshFace*)(node->element2);
#ifdef LOG_RECOVERY_PIPES
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " - swap 2x3 (face-face) ...";
#endif
		if(MeshGenerator3d::swap23(mc, mesh, common_face)){
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "ok.");
#endif
			return true;
		}else{
#ifdef LOG_RECOVERY_PIPES
			LOG4CPLUS_INFO(MeshLog::logger_mesh, " failed.");
#endif
			return false;
		}
	}

	return false;
}

bool MeshGenerator3dDelaunayBoundary::reduceEdgeRankForSwap44(Metric3dContext& mc, MeshContainer3d* mesh, MeshEdge3d *edge, MeshFace *first_face, MeshFace *second_face)
{
	for(int i = 0; i < 2; i++){
		MeshFace* face0 = first_face;
		MeshBlock* block0 = face0->getBlock(i);
		MeshFace* face1 = block0->getIncidentFace(face0, edge);
		MeshBlock* block1 = face1->getOtherBlock(block0);
		MeshFace* face2 = block1->getIncidentFace(face1, edge);
		while(face2 != second_face){
			MeshBlock* next_block = face2->getOtherBlock(block1);
			if(MeshGenerator3d::swap23(mc, mesh, face1) && 
			   edge->getFaceCount() == 4) return true;
			block0 = block1;
			block1 = next_block;
			face0 = face1;
			face1 = face2;
			face2 = block1->getIncidentFace(face1, edge);
		}
	}
	return false;
}

bool MeshGenerator3dDelaunayBoundary::MissingFace::adjacentToPoints(MeshPoint3d* pt1, MeshPoint3d* pt2) const
{
	return ( (m_pt[0] == pt1 || m_pt[1] == pt1 || m_pt[2] == pt1) &&
			 (m_pt[0] == pt2 || m_pt[1] == pt2 || m_pt[2] == pt2));
}

const ControlDataMatrix3d& MeshGenerator3dDelaunayBoundary::MissingFace::getSpecialMetric()
{
	if(valid_metric) return special_metric;

	DVector3d e[3];
	double v[3];
	e[0] = face->getPoint(1)->getCoordinates() - face->getPoint(0)->getCoordinates();
	v[0] = v[1] = e[0].length();
	e[0] /= v[0];
	v[2] = MeshGenerator3dDelaunayBoundary::param_special_metric_ratio * v[0];
	e[1] = face->getPoint(2)->getCoordinates() - face->getPoint(1)->getCoordinates();

	e[2] = e[0].crossProduct(e[1]).normalized();
	e[1] = e[0].crossProduct(e[2]).normalized();

	special_metric.setEigensystem(e, v);
	valid_metric = true;
	return special_metric;
}

const ControlDataMatrix3d& MeshGenerator3dDelaunayBoundary::MissingEdge::getSpecialMetric()
{
	if(valid_metric) return special_metric;

	DVector3d e[3];
	double v[3];
	e[0] = edge->getMeshPoint(1)->getCoordinates() - edge->getMeshPoint(0)->getCoordinates();
	v[0] = e[0].length();
	e[0] /= v[0];
	v[1] = v[2] = MeshGenerator3dDelaunayBoundary::param_special_metric_ratio * v[0];
	// second direction
	MeshEdge3d* other_edge = edge->getMeshPoint(0)->getEdge(0);
	if(other_edge == edge){
		assert(edge->getMeshPoint(0)->getRank() > 1);
		other_edge = edge->getMeshPoint(0)->getEdge(1);
	}
	e[1] = other_edge->getMeshPoint(1)->getCoordinates() - other_edge->getMeshPoint(0)->getCoordinates();
	e[2] = e[0].crossProduct(e[1]).normalized();
	e[1] = e[0].crossProduct(e[2]).normalized();

	special_metric.setEigensystem(e, v);
	valid_metric = true;
	return special_metric;
}

int MeshGenerator3dDelaunayBoundary::refineMeshForBoundaryTriangulation(Metric3dContext& mc, 
		MeshContainer3d* mesh, MeshPoint3d* point, const MeshContainer3d* boundary_mesh, int bindex)
{
//	assert(mc.validMetricAtPoint(point)); // shouldn't be used when skipPIng is turned off
	int inserted_count = 0;

	while(true){
		DataVector<MeshBlock*> blocks;
		point->adjacentBlocks(blocks);
		// select worst block
		MeshTetrahedron* worst_block = nullptr;
		for(int i = 0; i < blocks.countInt(); i++){
			if(blocks[i]->getQuality() < 0.0) // which is default for new "uncalculated" blocks
				blocks[i]->countQuality(mc, MeshData::QUALITY3D_MEAN_RATIO);
			if(blocks[i]->getQuality() > param_refined_ICD_ratio) continue;
			if(!worst_block || blocks[i]->getQuality() < worst_block->getQuality())
				worst_block = (MeshTetrahedron*)blocks[i];
		}
		if(!worst_block) return inserted_count;
		// insert point (longest edge?)
		int edge_count = worst_block->getEdgeCount();
		double max_len = 0.0;
		MeshEdge3d* longest_edge = nullptr;
		for(int j = 0; j < edge_count; j++){
			MeshEdge3d* edge = worst_block->getEdge(j);
			double len = edge->getLength(mc);
			if(len > 9.0 && (!longest_edge || len > max_len)){
				longest_edge = edge;
				max_len = len;
			}
		}
		if(!longest_edge){
			worst_block->setQuality(1.0); continue;
		}
		DPoint3d dnew = longest_edge->getPoint(0.5);
		mc.countMetricAtPoint(dnew);
		// check other boundary points
		int bpct = boundary_mesh->getPointsCount();
		const DMPoint3d dmnew = mc.transformRStoMS(dnew);
		for(int i = bindex; i < bpct; i++){
			if(dmnew.distance2(boundary_mesh->getPointAt(i)->getMetricCoordinates(mc)) < 2.0){
				worst_block->setQuality(1.0);
				longest_edge = nullptr;
				break;
			}
		}
		if(!longest_edge) continue;
		// ok, insert into mesh
		MeshPoint3d* new_point = new MeshPoint3d(dnew);
		if(MeshGenerator3d::addPointToTriangulation(mc, mesh, new_point, worst_block, worst_block, true)){
			++inserted_count;
		}else{
			delete new_point;
			worst_block->setQuality(1.0); 
			continue;
		}
	}
}

int MeshGenerator3dDelaunayBoundary::refineMeshForBoundaryTriangulation(Metric3dContext& mc, 
			MeshContainer3d* mesh, MeshPoint3d* point, const OctTreeMeshPoints& octree_points)
{
//	if(point->getRank() < 9) return 0;

//	assert(mc.validMetricAtPoint(point)); // shouldn't be used when skipPIng is turned off
	int inserted_count = 0;

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "* refinement - rank = " << point->getRank());

	while(true){
		DataVector<MeshBlock*> blocks;
		point->adjacentBlocks(blocks);
		// prune
		for(int i = 0; i < blocks.countInt(); ){
			if(blocks[i]->getQuality() < 0.0) // which is default for new "uncalculated" blocks
				blocks[i]->countQuality(mc, MeshData::QUALITY3D_MEAN_RATIO);
			if(blocks[i]->getQuality() > param_refined_ICD_ratio) 
				blocks.removeAt(i);
			else ++i;

		}
		if(blocks.empty()) return inserted_count;
		// refine worst blocks
		bool no_change = true;
		while(blocks.countInt() > 0){
			// select worst block
			int worst = 0;
			for(int i = 1; i < blocks.countInt(); i++)
				if(blocks[i]->getQuality() < blocks[worst]->getQuality()) worst = i;
			// insert point (longest edge?)
			int edge_count = blocks[worst]->getEdgeCount();
			double max_len = 0.0;
			MeshEdge3d* longest_edge = nullptr;
			for(int j = 0; j < edge_count; j++){
				MeshEdge3d* edge = blocks[worst]->getEdge(j);
				double len = edge->getLength(mc);
				if(len > 9.0 && (!longest_edge || len > max_len)){
					longest_edge = edge;
					max_len = len;
				}
			}
			if(!longest_edge){
				blocks.removeAt(worst); continue;
			}
			DPoint3d dnew = longest_edge->getPoint(0.5);
			mc.countMetricAtPoint(dnew);
			// check other boundary points
			if(octree_points.anyMeshPointInProximity(mc, dnew)){
				blocks.removeAt(worst); continue;
			}
			// ok, insert into mesh
			MeshPoint3d* new_point = new MeshPoint3d(dnew);
			if(MeshGenerator3d::addPointToTriangulation(mc, mesh, new_point, (MeshTetrahedron*)blocks[worst], 
				(MeshTetrahedron*)blocks[worst], true))
			{
				++inserted_count;
				blocks.clear();
				no_change = false;
			}else{
				delete new_point;
				blocks.removeAt(worst);
			}
		}
		if(no_change) return inserted_count;
//		int bad_block_count = 0;
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, " * refinement block check - bad count = " << bad_block_count);
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, " * refinement point inserted - current rank = " << point->getRank());
	}
}

int MeshGenerator3dDelaunayBoundary::refineMeshForBoundaryTriangulationByEdges(Metric3dContext& mc, 
			MeshContainer3d* mesh, MeshPoint3d* point, const OctTreeMeshPoints& octree_points)
{
	const int MIN_RANK = 8;

//	assert(mc.validMetricAtPoint(point)); // shouldn't be used when skipPIng is turned off
	int inserted_count = 0;

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "* refinement - rank = " << point->getRank());

	while(point->getRank() > MIN_RANK){
		int rank = point->getRank();
		DataVector<MeshEdge3d*> edges(rank);
		DataVector<double> lengths(rank);
		mc.countMetricAtPoint(point->getCoordinates());
		const DMPoint3d dpt = point->getMetricCoordinates(mc);
		// count
		for(int i = 0; i < rank; i++){
			MeshEdge3d* edge = point->getEdge(i);
			if(edge->isBorder()) continue;
			double dist2 = dpt.distance2(edge->getOtherPoint(point)->getMetricCoordinates(mc));
			if(dist2 > 81.0){
				edges.add(edge);
				lengths.add(dist2);
			}
		}
		if(edges.empty()) return inserted_count;
		// refine longest edge
		bool no_change = true;
		while(edges.countInt() > 0){
			// select longest edge
			int longest = 0;
			for(int i = 1; i < edges.countInt(); i++)
				if(lengths[i] > lengths[longest]) longest = i;
			DPoint3d dnew = edges[longest]->getPoint(0.5);
			mc.countMetricAtPoint(dnew);
			// check other boundary points
			if(octree_points.anyMeshPointInProximity(mc, dnew)){
				edges.removeAt(longest); lengths.removeAt(longest); continue;
			}
			// ok, insert into mesh
			MeshPoint3d* new_point = new MeshPoint3d(dnew);
			// edge is not border, so face is not border, so it has to have both blocks
			MeshTetrahedron* tetrahedron = (MeshTetrahedron*)edges[longest]->getFaceAt(0)->getBlock(0);
			if(MeshGenerator3d::addPointToTriangulation(mc, mesh, new_point, tetrahedron, tetrahedron, true))
			{
				++inserted_count;
				edges.clear();
				no_change = false;
			}else{
				delete new_point;
				edges.removeAt(longest);
				lengths.removeAt(longest);
			}
		}
		if(no_change) return inserted_count;
//		int bad_block_count = 0;
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, " * refinement block check - bad count = " << bad_block_count);
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, " * refinement point inserted - current rank = " << point->getRank());
	}
	return inserted_count;
}

/// Match missing edge (and adjacent faces) by matching with border-entities (with transformations on both sides)
int MeshGenerator3dDelaunayBoundary::matchBoundaryEdge(Metric3dContext& mc, const MissingEdge& me,
		DataHashTableKeyValue<MeshPoint3d*,MeshPoint3d*> & bmpoints)
{
	assert(!me.edge->availableTag(TagExtended::TAG_FIXED));
	if(me.edge->getFaceCount() != 2) return 0;
	MeshFace* bfaces[2] = { me.edge->getFaceAt(0), me.edge->getFaceAt(1) };
	if(bfaces[0]->availableTag(TagExtended::TAG_FIXED)) return 0;
	if(bfaces[1]->availableTag(TagExtended::TAG_FIXED)) return 0;
	// other edges of bfaces[x] has to be already recovered 
	// ---> (for now only 2-face cases - later maybe extended version...)
	// -> identify two additional mpoints
	MeshPoint3d* ebpoints[2] = {
		me.edge->getMeshPoint(0), me.edge->getMeshPoint(1) 
	};
	MeshPoint3d* bpoints[2] = {
		bfaces[0]->getOtherPoint(ebpoints[0], ebpoints[1]),
		bfaces[1]->getOtherPoint(ebpoints[0], ebpoints[1])
	};
	MeshPoint3d* mpoints[2] = {
		bmpoints.getValue(bpoints[0], nullptr),
		bmpoints.getValue(bpoints[1], nullptr)
	};
	assert( mpoints[0] && mpoints[1] );
	// -> check edges
	MeshPoint3d* empoints[2] = {
		bmpoints.getValue(ebpoints[0], nullptr),
		bmpoints.getValue(ebpoints[1], nullptr)
	};
	assert(empoints[0] && empoints[1] );
	assert(empoints[0]->getEdgeToPoint(empoints[1]) == nullptr);

	MeshEdge3d* e00 = empoints[0]->getEdgeToPoint(mpoints[0]);
	if(!e00) return 0;
	MeshEdge3d* e01 = empoints[0]->getEdgeToPoint(mpoints[1]);
	if(!e01) return 0;
	MeshEdge3d* e10 = empoints[1]->getEdgeToPoint(mpoints[0]);
	if(!e10) return 0;
	MeshEdge3d* e11 = empoints[1]->getEdgeToPoint(mpoints[1]);
	if(!e11) return 0;

	double dist2 = DSegment3d::distance2ToSegment(
		empoints[0]->getMetricCoordinates(mc).toRealSpace(),
		empoints[1]->getMetricCoordinates(mc).toRealSpace(),
		mpoints[0]->getMetricCoordinates(mc).toRealSpace(),
		mpoints[1]->getMetricCoordinates(mc).toRealSpace());
	// ???
	if(dist2 > 0.01) return 0; // distance 0.1 in metric space

	return 0;
}

//#define SHOW_RECOVER_BOUNDARY_CAVITY

/// Try recover boundary edge (and adjacent faces) by retrieving local cavity and remeshing
bool MeshGenerator3dDelaunayBoundary::recoverBoundaryWithLocalCavityRemeshing(
	Metric3dContext& mc, MeshContainer3d* mesh, 
	MeshPoint3d* pt1, MeshPoint3d* pt2, 
	const DataVector<MissingFace> & missing_faces)
{
	//if(true){
	//	MeshPoint3d* bpt1 = (MeshPoint3d*) pt1->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
	//	MeshPoint3d* bpt2 = (MeshPoint3d*) pt2->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "recoverBoundaryWithLocalCavityRemeshing: ");
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, pt1->getIndex() << "(" << (bpt1 ? bpt1->getIndex() : -1) << ") -> " << pt1->getCoordinates());
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, pt2->getIndex() << "(" << (bpt2 ? bpt2->getIndex() : -1) << ") -> " << pt2->getCoordinates());
	//}

	// 1. Gather group of adjacent missing faces

	DataVector<MissingEdge> local_edges(10);
	local_edges.add(MissingEdge(nullptr, pt1, pt2));

	DataHashTable<MeshPoint3d*> cavity_fpoints(50, nullptr);
	DataVector<MeshEdge3d*> cavity_fedges(10);

	DataVector<int> local_faces(missing_faces.countInt());
	for(int k = 0; k < local_edges.countInt(); k++){
		MissingEdge me = local_edges[k];
		for(int i = 0; i < missing_faces.countInt(); i++){
			const MissingFace& mf = missing_faces[i];
			if(mf.adjacentToPoints(me.m_pt0, me.m_pt1)){
				if(local_faces.addIfNew(i)){
					cavity_fpoints.insert(mf.m_pt[0]);
					cavity_fpoints.insert(mf.m_pt[1]);
					cavity_fpoints.insert(mf.m_pt[2]);
					// check for additional missing edges in this face
					for(int m = 0; m < 3; m++){
						if( (mf.m_pt[m] != me.m_pt0) && (mf.m_pt[m] != me.m_pt1) ) { // mf.m_pt[m] is the other point
							MeshEdge3d* edge = mf.m_pt[m]->getEdgeToPoint(mf.m_pt[(m+1)%3]);
							if( edge ) cavity_fedges.addIfNew(edge);
							else local_edges.add(MissingEdge(nullptr, mf.m_pt[m], mf.m_pt[(m+1)%3]));
							edge = mf.m_pt[m]->getEdgeToPoint(mf.m_pt[(m+2)%3]);
							if( edge ) cavity_fedges.addIfNew(edge);
							else local_edges.add(MissingEdge(nullptr, mf.m_pt[m], mf.m_pt[(m+2)%3]));
							break;
						}
					}
				}
			}
		}
	}

//	if(local_edges.countInt() > 1) return false;

	// 2. Identify all crossing blocks
	DataVector<MeshBlock*> crossing_blocks[3];
	DataHashTable<MeshEdge3d*> used_edges(100, nullptr);

	mc.countMetricAtPoints(pt1, pt2);
	for(int i = 0; i < local_faces.countInt(); i++){
		const MissingFace& mf = missing_faces[local_faces[i]];
		for(int m = 0; m < 3; m++){
			int m1 = (m+1)%3;
			MeshEdge3d* edge = mf.m_pt[m]->getEdgeToPoint(mf.m_pt[m1]);
			if(!edge) continue;
			int efct = edge->getFaceCount();
			bool edge_processed = false;
			for(int j = 0; !edge_processed && (j < efct); j++){
				MeshFace* face = edge->getFaceAt(j);
				MeshPoint3d* other_pt = face->getOtherPoint(mf.m_pt[m], mf.m_pt[m1]);
				if(cavity_fpoints.contains(other_pt)){ // check for face-on-face collisions
					bool colliding = true;
					MeshEdge3d* check_edge = other_pt->getEdgeToPoint(mf.m_pt[m]);
					assert(check_edge);
					if(! used_edges.contains(check_edge)){
						for(int k = 0; colliding && (k < cavity_fedges.countInt()); k++){
							MeshPoint3d* ce_mp0 = cavity_fedges[k]->getMeshPoint(0);
							MeshPoint3d* ce_mp1 = cavity_fedges[k]->getMeshPoint(1);
							if( (ce_mp0 == other_pt && ce_mp1 == mf.m_pt[m]) ||
								(ce_mp1 == other_pt && ce_mp0 == mf.m_pt[m]) ) colliding = false; 
						}
						if(colliding){
							used_edges.insert(check_edge);
							if(!check_edge->isBorder()){
								check_edge->adjacentBlocks(crossing_blocks[0]);
								edge_processed = true;
							}
						}
					} else edge_processed = true;
					if(!colliding){ // check the other edge
						check_edge = other_pt->getEdgeToPoint(mf.m_pt[m1]);
						assert(check_edge);
						colliding = true;
						if(! used_edges.contains(check_edge)){
							for(int k = 0; colliding && (k < cavity_fedges.countInt()); k++){
								MeshPoint3d* ce_mp0 = cavity_fedges[k]->getMeshPoint(0);
								MeshPoint3d* ce_mp1 = cavity_fedges[k]->getMeshPoint(1);
								if( (ce_mp0 == other_pt && ce_mp1 == mf.m_pt[m1]) ||
									(ce_mp1 == other_pt && ce_mp0 == mf.m_pt[m1]) ) colliding = false; 
							}
							if(colliding){
								used_edges.insert(check_edge);
								if(!check_edge->isBorder()){
									check_edge->adjacentBlocks(crossing_blocks[0]);
									edge_processed = true;
								}
							}
						} else edge_processed = true;
					}
				}
			}
			if(!edge_processed){
				DMPoint3d dpts[3] = {
					mf.m_pt[0]->getMetricCoordinates(mc),
					mf.m_pt[1]->getMetricCoordinates(mc),
					mf.m_pt[2]->getMetricCoordinates(mc) };
				// check all adjacent blocks for edges crossing any missing boundary face
				DataVector<MeshBlock*> adjacent_blocks;
				if(edge->adjacentBlocks(adjacent_blocks)){
					for(int j = 0; !edge_processed && (j < adjacent_blocks.countInt()); j++){
						MeshTetrahedron* tetra = (MeshTetrahedron*) adjacent_blocks[j];
						MeshEdge3d* opposite_edge = nullptr;
						for(int k = 0; k < tetra->getEdgeCount(); k++)
							if( ! (opposite_edge=tetra->getEdge(k))->incidentTo(edge)) break;
						if( ! used_edges.contains(opposite_edge) ){
							if( DMTriangle3d::crossSegment(
									opposite_edge->getMeshPoint(0)->getMetricCoordinates(mc),
									opposite_edge->getMeshPoint(1)->getMetricCoordinates(mc),
									dpts[0], dpts[1], dpts[2]) )
							{
								used_edges.insert(opposite_edge);
								if(!opposite_edge->isBorder()){
									opposite_edge->adjacentBlocks(crossing_blocks[0]);
									edge_processed = true;
								}
							}
						}else edge_processed = true;
					}
				}
			}
		}
	}

	if(crossing_blocks[0].empty()){
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Edge recovery cavity: crossing blocks empty");
		return false;
	}

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
	if(true){
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < crossing_blocks[0].countInt(); i++){
			set->addBlockWithEdges(crossing_blocks[0][i]);
		}
		set->addPoint(pt1, 3);
		set->addPoint(pt2, 3);
		set->addEdge(pt1->getCoordinates(), pt2->getCoordinates(), 3);

		SHOW_MESH("cavity blocks", set);
	}
#endif

	// additional check of edges of cavity blocks
	for(int i = 0; i < crossing_blocks[0].countInt(); i++)
		crossing_blocks[1].add(crossing_blocks[0][i]);

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Recovery edge cavity, checking crossing blocks... ");
	for(int i = 0; i < crossing_blocks[1].countInt(); i++){
		MeshBlock* block = crossing_blocks[1][i];
		double vol = block->getVolume(mc, false);
		double q = ((MeshTetrahedron*)block)->getMeanRatio(mc, false);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Recovery edge cavity, crossing volume = " << vol << ", q = " << q);
		bool block_too_small = (vol < 0.01);
//		bool block_too_small = false;
		int ect = block->getEdgeCount();
		for(int j = 0; j < ect; j++){
			MeshEdge3d* edge = block->getEdge(j);
			if(edge->isBorder() || used_edges.contains(edge)) continue;
			MeshPoint3d* ept0 = edge->getMeshPoint(0);
			MeshPoint3d* ept1 = edge->getMeshPoint(1);
			DMSegment3d seg(ept0->getMetricCoordinates(mc), ept1->getMetricCoordinates(mc));
			for(int k = 0; k < local_edges.countInt(); k++){
				double dist2 = 0.0;
				if(!block_too_small){ // if too smal, add always, else actually check the dist2
					const MissingEdge& me = local_edges[k];
					if( me.m_pt0 == ept0 || me.m_pt0 == ept1 ||
						me.m_pt1 == ept0 || me.m_pt1 == ept1) continue; // adjacent - don't check for crossing
					dist2 = seg.distance2ToSegment(
						me.m_pt0->getMetricCoordinates(mc),
						me.m_pt1->getMetricCoordinates(mc));
				}
				if(dist2 < 0.0001){
					used_edges.insert(edge);
					edge->adjacentBlocks(crossing_blocks[1]);
				}
			}
		}
	}

	DataVector< std::shared_ptr<DataHashTable<MeshBlock*>> > hblocks;
	hblocks.add(std::make_shared<DataHashTable<MeshBlock*>>(4*crossing_blocks[0].countInt(), nullptr));
	hblocks.add(std::make_shared<DataHashTable<MeshBlock*>>(4*crossing_blocks[1].countInt(), nullptr));
	hblocks.add(std::make_shared<DataHashTable<MeshBlock*>>(4*crossing_blocks[2].countInt(), nullptr));

	for(int i = 0; i < crossing_blocks[0].countInt(); i++)
		hblocks[0]->insert(crossing_blocks[0][i]);

	for(int i = 0; i < crossing_blocks[1].countInt(); i++){
		crossing_blocks[2].add(crossing_blocks[1][i]);
		hblocks[1]->insert(crossing_blocks[1][i]);
		hblocks[2]->insert(crossing_blocks[1][i]);
	}

	// + extra additional blocks
	for(int i = 0; i < crossing_blocks[2].countInt(); i++){
		MeshBlock* block = crossing_blocks[2][i];
		int bfct = block->getFaceCount();
		for(int j = 0; j < bfct; j++){
			MeshBlock* other_block = block->getNeighbour(j);
			if(!other_block || hblocks[2]->contains(other_block)) continue;
			// ... if at least two neighbours of this extra block are in cavity, add also this one
			int obfct = other_block->getFaceCount();
			int cav_count = 0;
			for(int k = 0; k < obfct; k++){
				MeshBlock* nblock = other_block->getNeighbour(k);
				if(nblock && hblocks[2]->contains(nblock)) cav_count++;
			}
			if(cav_count > 1){
				crossing_blocks[2].add(other_block);
				hblocks[2]->insert(other_block);
			}
		}
	}

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
	if(crossing_blocks[0].countInt() < crossing_blocks[2].countInt()){
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < crossing_blocks[2].countInt(); i++){
			set->addBlockWithEdges(crossing_blocks[2][i]);
		}
		set->addPoint(pt1, 3);
		set->addPoint(pt2, 3);
		set->addEdge(pt1->getCoordinates(), pt2->getCoordinates(), 3);

		SHOW_MESH("cavity blocks extended", set);
	}
#endif

	for(int step_i = 2; step_i >= 0; step_i--){
		if(retriangulateCavity(mc, mesh, crossing_blocks[step_i], *hblocks[step_i], 
				cavity_fpoints, cavity_fedges, local_edges, local_faces, missing_faces))
			return true;
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
		"Edge recovery cavity: couldn't find star center for all set of crossing blocks");
	return false;

	//  (b) greedy frontal approach?
}

bool MeshGenerator3dDelaunayBoundary::retriangulateCavity(
		Metric3dContext& mc, MeshContainer3d* mesh, 
		const DataVector<MeshBlock*> & crossing_blocks,
		const DataHashTable<MeshBlock*> & hblocks,
		const DataHashTable<MeshPoint3d*> & cavity_fpoints,
		const DataVector<MeshEdge3d*> & cavity_fedges,
		const DataVector<MeshGenerator3dDelaunayBoundary::MissingEdge> & local_edges,
		const DataVector<int> local_faces,
		const DataVector<MeshGenerator3dDelaunayBoundary::MissingFace> & missing_faces)
{
	// 3. Create two cavities (divide faces into two sets)
	//  + additional check for extra blocks (if some "cavity" face is still crossing the missing edge)
	DataVector<MeshFace*> cfaces(10*crossing_blocks.countInt());

	for(int i = 0; i < crossing_blocks.countInt(); i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*) crossing_blocks[i];
		int fct = tetra->getFaceCount();
		for(int j = 0; j < fct; j++){
			MeshFace* face = tetra->getFace(j);
			MeshBlock* block = face->getOtherBlock(tetra);
			if(!block) cfaces.add(face);
			else if(!hblocks.contains(block)){
				if( cavity_fpoints.contains(face->getPoint(0)) ||
					cavity_fpoints.contains(face->getPoint(1)) ||
					cavity_fpoints.contains(face->getPoint(2)))
				{
					cfaces.add(face); // adjacent to base_cavity edges, definitely add
				}else{
					// check, if face is crossing any missing edge for this case
					DMPoint3d dpts[3] = {
						face->getPoint(0)->getMetricCoordinates(mc),
						face->getPoint(1)->getMetricCoordinates(mc),
						face->getPoint(2)->getMetricCoordinates(mc) };
					bool crossing = false;
					for(int k = 0; !crossing && (k < local_edges.countInt()); k++){
						crossing = DMTriangle3d::crossSegment(
								local_edges[k].m_pt0->getMetricCoordinates(mc),
								local_edges[k].m_pt1->getMetricCoordinates(mc),
								dpts[0], dpts[1], dpts[2]);
						if(crossing) {
							return false;
						}else
							cfaces.add(face);
					}
				}
			}
		}
	}

	// 3a. Orientate faces
	DataHashTableKeyValue<MeshFace*,int> hfaces(2*cfaces.countInt(), nullptr);
	for(int i = 0; i < cfaces.countInt(); i++)
		hfaces.insert(cfaces[i], i);
	DataHashTable<MeshEdge3d*> hedges(2*cavity_fedges.countInt(), nullptr);
	for(int i = 0; i < cavity_fedges.countInt(); i++)
		hedges.insert(cavity_fedges[i]);

	DataVector<bool> cavity_upper_faces(cfaces.countInt(), false);
	DataVector<int> cavity_stack(cfaces.countInt());
	cavity_upper_faces[0] = true;
	cavity_stack.add(0);
	for(int i = 0; i < cavity_stack.countInt(); i++){
		MeshFace* face = cfaces[cavity_stack[i]];
		int ect = face->getEdgeCount();
		for(int j = 0; j < ect; j++){
			MeshEdge3d* edge = face->getEdge(j);
			if(hedges.contains(edge)) continue; // without "crossing" the cavity base contour
			int efct = edge->getFaceCount();
			for(int k = 0; k < efct; k++){
				MeshFace* other_face = edge->getFaceAt(k);
				int id = hfaces.getValue(other_face, -1);
				if(id < 0 || cavity_upper_faces[id]) continue; // if not-cavity or already marked as "upper"
				// ... else
				cavity_upper_faces[id] = true;
				cavity_stack.add(id);
			}
		}
	}

	assert(cavity_stack.countInt() < cfaces.countInt()); 

	DataVector<bool> proper_orientation(cfaces.countInt(), true);
	for(int i = 0; i < cfaces.countInt(); i++){
		MeshBlock* block = cfaces[i]->getBlock(0);
		if(!block || !hblocks.contains(block)) 
			proper_orientation[i] = false;
	}

	// --> check if the cavities are "watertight"
	{
		MeshContainer3dSurface upper_test_mesh(3*cfaces.countInt());
		MeshContainer3dSurface lower_test_mesh(3*cfaces.countInt());
		DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> upper_test_hpoints(3*cfaces.countInt(), nullptr);
		DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> lower_test_hpoints(3*cfaces.countInt(), nullptr);
		for(int i = 0; i < cfaces.countInt(); i++){
			MeshFace* face = cfaces[i];
			int fpct = face->getPointCount(); assert(fpct == 3);
			DataVector<MeshPoint3d*> test_face_points(fpct);
			for(int j = 0; j < fpct; j++){
				MeshPoint3d* point = face->getPoint(j);
				MeshPoint3d* test_point = nullptr;
				if(cavity_upper_faces[i]){
					test_point = upper_test_hpoints.getValue(point, nullptr);
					if(!test_point){
						upper_test_hpoints.insert( point, test_point = new MeshPoint3d(*point) );
						upper_test_mesh.addMeshPoint(test_point);
					}
				}else{
					test_point = lower_test_hpoints.getValue(point, nullptr);
					if(!test_point){
						lower_test_hpoints.insert( point, test_point = new MeshPoint3d(*point) );
						lower_test_mesh.addMeshPoint(test_point);
					}
				}
				test_face_points.add(test_point);
			}
			// orientation is not really important here...
			MeshFace* test_face = 
				new MeshTriangle3d(test_face_points[0], test_face_points[1], test_face_points[2]);
			if(cavity_upper_faces[i])
				upper_test_mesh.addMeshFace(test_face);
			else
				lower_test_mesh.addMeshFace(test_face);
		}
		for(int i = 0; i < local_faces.countInt(); i++){
			const MissingFace& mf = missing_faces[local_faces[i]];
			DataVector<MeshPoint3d*> upper_test_face_points(3);
			DataVector<MeshPoint3d*> lower_test_face_points(3);
			for(int j = 0; j < 3; j++){
				MeshPoint3d* test_point = upper_test_hpoints.getValue(mf.m_pt[j], nullptr);
				if(!test_point){
					LOG4CPLUS_WARN(MeshLog::logger_mesh, "Edge recovery cavity: missing upper test point");
					return false;
				}
				upper_test_face_points.add(test_point);
				test_point = lower_test_hpoints.getValue(mf.m_pt[j], nullptr);
				if(!test_point){
					LOG4CPLUS_WARN(MeshLog::logger_mesh, "Edge recovery cavity: missing lower test point");
					return false;
				}
				lower_test_face_points.add(test_point);
			}
			upper_test_mesh.addMeshFace(new MeshTriangle3d(upper_test_face_points[0], 
				upper_test_face_points[1], upper_test_face_points[2]));
			lower_test_mesh.addMeshFace(new MeshTriangle3d(lower_test_face_points[0], 
				lower_test_face_points[1], lower_test_face_points[2]));
		}
		// and now, check these surface meshes
		for(IteratorEdge3d it = upper_test_mesh.getFirstEdge3d(); it.isValid(); it.nextEdge())
			if(it.getEdge()->getFaceCount() != 2){
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Edge recovery cavity: upper test face != 2");
				return false;
			}
		for(IteratorEdge3d it = lower_test_mesh.getFirstEdge3d(); it.isValid(); it.nextEdge())
			if(it.getEdge()->getFaceCount() != 2){
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Edge recovery cavity: upper test face != 2");
				return false;
			}

		// plus, check for "lost" nodes
		for(int i = 0; i < crossing_blocks.countInt(); i++){
			MeshBlock* block = crossing_blocks[i];
			for(int j = 0; j < block->getPointCount(); j++){
				MeshPoint3d* point = block->getPoint(j);
				if( upper_test_hpoints.contains(point) ) continue; // ok
				if( lower_test_hpoints.contains(point) ) continue; // ok
				// else -> there is point within the blocks, that is not on the cavity boundary ...
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Edge recovery cavity: extra stranded point for cavity");
				return false;
			}
		}
	}

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
	if(true){
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < cfaces.countInt(); i++){
			if(cavity_upper_faces[i])
				set->addFaceWithEdges(cfaces[i], -2, MeshViewSet::param_shrink, proper_orientation[i]);
		}
		set->addPoint(local_edges[0].m_pt0, 3);
		set->addPoint(local_edges[0].m_pt1, 3);
		set->addEdge(local_edges[0].m_pt0->getCoordinates(), local_edges[0].m_pt1->getCoordinates(), 3);

		SHOW_MESH("cavity upper faces", set);

		set = new MeshViewSet;
		for(int i = 0; i < cfaces.countInt(); i++){
			if(!cavity_upper_faces[i])
				set->addFaceWithEdges(cfaces[i], -2, MeshViewSet::param_shrink, proper_orientation[i]);
		}
		set->addPoint(local_edges[0].m_pt0, 3);
		set->addPoint(local_edges[0].m_pt1, 3);
		set->addEdge(local_edges[0].m_pt0->getCoordinates(), local_edges[0].m_pt1->getCoordinates(), 3);

		SHOW_MESH("cavity lower faces", set);
	}
#endif

	DataVector<DMTriangle3d> orient_faces_upper(cfaces.countInt() + 2*local_faces.countInt());
	DataVector<DMTriangle3d> orient_faces_lower(cfaces.countInt() + 2*local_faces.countInt());

	for(int i = 0; i < cfaces.countInt(); i++){
		DMTriangle3d dmt;
		if(proper_orientation[i])
			dmt = DMTriangle3d(
					cfaces[i]->getPoint(0)->getMetricCoordinates(mc),
					cfaces[i]->getPoint(1)->getMetricCoordinates(mc),
					cfaces[i]->getPoint(2)->getMetricCoordinates(mc));
		else
			dmt = DMTriangle3d(
					cfaces[i]->getPoint(1)->getMetricCoordinates(mc),
					cfaces[i]->getPoint(0)->getMetricCoordinates(mc),
					cfaces[i]->getPoint(2)->getMetricCoordinates(mc));
		if(cavity_upper_faces[i])
			orient_faces_upper.add(dmt);
		else
			orient_faces_lower.add(dmt);
	}

	// - add missing faces for upper/lower cavity (using edge info)
	DataVector<bool> local_orientation_upper(local_faces.countInt());
	DataVector<bool> local_orientation_lower(local_faces.countInt());
	for(int i = 0; i < local_faces.countInt(); i++){
		const MissingFace& mf = missing_faces[local_faces[i]];
		bool upper_done = false;
		bool lower_done = false;
		for(int j = 0; j < cfaces.countInt(); j++){
			if(cavity_upper_faces[j]){
				if(upper_done) continue;
				MeshFace* face = cfaces[j];
				for(int k = 0; k < 3; k++){
					if(face->incidentToPoints(mf.m_pt[k], mf.m_pt[(k+1)%3])){
						bool oriented_ok = proper_orientation[j];
						if(face->properOrientation(mf.m_pt[k], mf.m_pt[(k+1)%3])){
							oriented_ok = !oriented_ok;
						}
						if(oriented_ok){
							orient_faces_upper.add( 
								DMTriangle3d(
									mf.m_pt[0]->getMetricCoordinates(mc),
									mf.m_pt[1]->getMetricCoordinates(mc),
									mf.m_pt[2]->getMetricCoordinates(mc)));
						}else{
							orient_faces_upper.add( 
								DMTriangle3d(
									mf.m_pt[1]->getMetricCoordinates(mc),
									mf.m_pt[0]->getMetricCoordinates(mc),
									mf.m_pt[2]->getMetricCoordinates(mc)));
						}
						local_orientation_upper.add(oriented_ok);
						upper_done = true;
					}
				}
			}else{
				if(lower_done) continue;
				MeshFace* face = cfaces[j];
				for(int k = 0; k < 3; k++){
					if(face->incidentToPoints(mf.m_pt[k], mf.m_pt[(k+1)%3])){
						bool oriented_ok = proper_orientation[j];
						if(face->properOrientation(mf.m_pt[k], mf.m_pt[(k+1)%3])){
							oriented_ok = !oriented_ok;
						}
						if(oriented_ok){
							orient_faces_lower.add( 
								DMTriangle3d(
									mf.m_pt[0]->getMetricCoordinates(mc),
									mf.m_pt[1]->getMetricCoordinates(mc),
									mf.m_pt[2]->getMetricCoordinates(mc)));
						}else{
							orient_faces_lower.add( 
								DMTriangle3d(
									mf.m_pt[1]->getMetricCoordinates(mc),
									mf.m_pt[0]->getMetricCoordinates(mc),
									mf.m_pt[2]->getMetricCoordinates(mc)));
						}
						local_orientation_lower.add(oriented_ok);
						lower_done = true;
					}
				}
			}
		}
	}

	DMPoint3d center_upper, center_lower;
	if( findStarCenter(orient_faces_upper, center_upper) &&
		findStarCenter(orient_faces_lower, center_lower) )
	{

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
		MeshViewSet* set = new MeshViewSet;
		double min_volume = LARGE_NUMBER;
#endif

		// 4. Retriangulate cavities 
		//	- remove crossing blocks
		for(int i = 0; i < crossing_blocks.countInt(); i++)
			delete mesh->removeMeshBlock(crossing_blocks[i]);
		//	- insert new mesh points
		MeshPoint3d* mp_upper = new MeshPoint3d(mc.transformMStoRS(center_upper));
		mesh->addMeshPoint(mp_upper);
		MeshPoint3d* mp_lower = new MeshPoint3d(mc.transformMStoRS(center_lower));
		mesh->addMeshPoint(mp_lower);
		//	- insert new tetrahedra for cavity faces
		for(int i = 0; i < cfaces.countInt(); i++){
			MeshFace* face = cfaces[i];
			MeshTetrahedron* tetra = nullptr;
			if(cavity_upper_faces[i]){
				if(proper_orientation[i])
					tetra = new MeshTetrahedron(
						face->getPoint(0), face->getPoint(1), face->getPoint(2), mp_upper);
				else
					tetra = new MeshTetrahedron(
						face->getPoint(1), face->getPoint(0), face->getPoint(2), mp_upper);
			}else{
				if(proper_orientation[i])
					tetra = new MeshTetrahedron(
						face->getPoint(0), face->getPoint(1), face->getPoint(2), mp_lower);
				else
					tetra = new MeshTetrahedron(
						face->getPoint(1), face->getPoint(0), face->getPoint(2), mp_lower);
			}
			mesh->addMeshBlock(tetra);

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
			set->addBlockWithEdges(tetra);
			double v = tetra->getVolume(mc, false);
			if(v < min_volume) min_volume = v;
#endif

		}
		//	- reconstruct missing faces
		for(int i = 0; i < local_faces.countInt(); i++){
			const MissingFace& mf = missing_faces[local_faces[i]];
			MeshTetrahedron* tetra = nullptr;
			if(local_orientation_upper[i])
				tetra = new MeshTetrahedron(mf.m_pt[0], mf.m_pt[1], mf.m_pt[2], mp_upper);
			else
				tetra = new MeshTetrahedron(mf.m_pt[1], mf.m_pt[0], mf.m_pt[2], mp_upper);
			mesh->addMeshBlock(tetra);

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
			set->addBlockWithEdges(tetra);
			double v = tetra->getVolume(mc, false);
			if(v < min_volume) min_volume = v;
#endif

			if(local_orientation_lower[i])
				tetra = new MeshTetrahedron(mf.m_pt[0], mf.m_pt[1], mf.m_pt[2], mp_lower);
			else
				tetra = new MeshTetrahedron(mf.m_pt[1], mf.m_pt[0], mf.m_pt[2], mp_lower);
			mesh->addMeshBlock(tetra);

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
			set->addBlockWithEdges(tetra);
			v = tetra->getVolume(mc, false);
			if(v < min_volume) min_volume = v;
#endif
		}

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "New tetrahedra min volume", min_volume);
		SHOW_MESH("New tetrahedra", set);
#endif

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Recovery edge cavity, OK.");
		return true;

	}

	return false;
}

/// Try recover boundary facee (single) by retrieving local cavity and remeshing
bool MeshGenerator3dDelaunayBoundary::recoverBoundaryWithLocalCavityRemeshing(
	Metric3dContext& mc, MeshContainer3d* mesh, 
	MeshPoint3d* pt1, MeshPoint3d* pt2, MeshPoint3d* pt3)
{
	MeshPoint3d* fpoints[3] = { pt1, pt2, pt3 };
	MeshEdge3d* fedges[3] = {
		pt1->getEdgeToPoint(pt2),
		pt2->getEdgeToPoint(pt3),
		pt3->getEdgeToPoint(pt1)};

	if(!fedges[0] || !fedges[1] || !fedges[2]) return false;

	// 1. Identify all crossing blocks
	DataVector<MeshBlock*> crossing_blocks(100);
	DataHashTable<MeshEdge3d*> used_edges(100, nullptr);

	mc.countMetricAtPoints(pt1, pt2);
	DMPoint3d dpts[3] = {
		pt1->getMetricCoordinates(mc),
		pt2->getMetricCoordinates(mc),
		pt3->getMetricCoordinates(mc) };
	for(int m = 0; m < 3; m++){
		int m1 = (m+1)%3;
		// check all adjacent blocks for edges crossing the missing boundary face
		DataVector<MeshBlock*> adjacent_blocks;
		if(fedges[m]->adjacentBlocks(adjacent_blocks)){
			//if(true){
			//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "== adjacent blocks ==");
			//	MeshViewSet* dset = new MeshViewSet;
			//	for(int k = 0; k < adjacent_blocks.countInt(); k++){
			//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Vol-mc: " << adjacent_blocks[k]->getVolume(mc, false));
			//		dset->addBlockWithEdges(adjacent_blocks[k]);
			//	}
			//	for(int j = 0; j < 3; j++){
			//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Point id=" << fpoints[j]->getIndex();
			//		MeshPoint3d* bpoint = (MeshPoint3d*)fpoints[j]->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
			//		if(bpoint)
			//			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, ", bid=" << bpoint->getIndex();
			//		LOG4CPLUS_INFO(MeshLog::logger_mesh, ", coord= " << fpoints[j]->getCoordinates());
			//	}
			//	dset->addEdge(fedges[m], 3);
			//	dset->addEdge(fedges[(m+1)%3], 2);
			//	dset->addEdge(fedges[(m+2)%3], 2);
			//	dset->addPoint(pt1, 3);
			//	dset->addPoint(pt2, 3);
			//	dset->addPoint(pt3, 3);
			//	SHOW_MESH("xxx check", dset);
			//}
			for(int j = 0; j < adjacent_blocks.countInt(); j++){
				MeshTetrahedron* tetra = (MeshTetrahedron*) adjacent_blocks[j];
				MeshEdge3d* opposite_edge = nullptr;
				for(int k = 0; k < tetra->getEdgeCount(); k++)
					if( ! (opposite_edge=tetra->getEdge(k))->incidentTo(fedges[m])) break;
				if( ! used_edges.contains(opposite_edge) ){
					if( DMTriangle3d::crossSegment(
							opposite_edge->getMeshPoint(0)->getMetricCoordinates(mc),
							opposite_edge->getMeshPoint(1)->getMetricCoordinates(mc),
							dpts[0], dpts[1], dpts[2]) )
					{
						used_edges.insert(opposite_edge);
						if(!opposite_edge->isBorder()){
							opposite_edge->adjacentBlocks(crossing_blocks);
							break;
						}
					}
				}
			}
		}
	}

	if(crossing_blocks.empty()) return false;

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
	if(true){
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < crossing_blocks.countInt(); i++){
			set->addBlockWithEdges(crossing_blocks[i]);
		}
		set->addPoint(pt1, 3);
		set->addPoint(pt2, 3);
		set->addPoint(pt3, 3);

		SHOW_MESH("cavity blocks", set);
	}
#endif

	// additional check of edges of cavity blocks
	int cross_count_before = crossing_blocks.countInt();
	for(int i = 0; i < crossing_blocks.countInt(); i++){
		MeshBlock* block = crossing_blocks[i];
//		bool block_too_small = (block->getVolume(mc, false) < METRIC_SMALL_NUMBER);
		bool block_too_small = false;
		int ect = block->getEdgeCount();
		for(int j = 0; j < ect; j++){
			MeshEdge3d* edge = block->getEdge(j);
			if(edge->isBorder() || used_edges.contains(edge)) continue;
			MeshPoint3d* ept0 = edge->getMeshPoint(0);
			MeshPoint3d* ept1 = edge->getMeshPoint(1);
			if( DMTriangle3d::crossSegment(
					ept0->getMetricCoordinates(mc),
					ept1->getMetricCoordinates(mc),
					dpts[0], dpts[1], dpts[2]) )
			{
				used_edges.insert(edge);
				edge->adjacentBlocks(crossing_blocks);
			}
		}
	}

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
	if(cross_count_before < crossing_blocks.countInt()){
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < crossing_blocks.countInt(); i++){
			set->addBlockWithEdges(crossing_blocks[i]);
		}
		set->addPoint(pt1, 3);
		set->addPoint(pt2, 3);
		set->addPoint(pt3, 3);

		SHOW_MESH("cavity blocks extended", set);
	}
#endif

	// 3. Create two cavities (divide faces into two sets)
	DataHashTable<MeshBlock*> hblocks(4*crossing_blocks.countInt(), nullptr);
	for(int i = 0; i < crossing_blocks.countInt(); i++)
		hblocks.insert(crossing_blocks[i]);
	DataVector<MeshFace*> cfaces(10*crossing_blocks.countInt());

	for(int i = 0; i < crossing_blocks.countInt(); i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*) crossing_blocks[i];
		int fct = tetra->getFaceCount();
		for(int j = 0; j < fct; j++){
			MeshFace* face = tetra->getFace(j);
			MeshBlock* block = face->getOtherBlock(tetra);
			if(!block || !hblocks.contains(block)) cfaces.add(face);
		}
	}

	// 3a. Orientate faces
	DataHashTableKeyValue<MeshFace*,int> hfaces(2*cfaces.countInt(), nullptr);
	for(int i = 0; i < cfaces.countInt(); i++)
		hfaces.insert(cfaces[i], i);
	DataHashTable<MeshEdge3d*> hedges(6, nullptr);
	for(int i = 0; i < 3; i++)
		hedges.insert(fedges[i]);

	DataVector<bool> cavity_upper_faces(cfaces.countInt(), false);
	DataVector<int> cavity_stack(cfaces.countInt());
	cavity_upper_faces[0] = true;
	cavity_stack.add(0);
	for(int i = 0; i < cavity_stack.countInt(); i++){
		MeshFace* face = cfaces[cavity_stack[i]];
		int ect = face->getEdgeCount();
		for(int j = 0; j < ect; j++){
			MeshEdge3d* edge = face->getEdge(j);
			if(hedges.contains(edge)) continue; // without "crossing" the cavity base contour
			int efct = edge->getFaceCount();
			for(int k = 0; k < efct; k++){
				MeshFace* other_face = edge->getFaceAt(k);
				int id = hfaces.getValue(other_face, -1);
				if(id < 0 || cavity_upper_faces[id]) continue; // if not-cavity or already marked as "upper"
				// ... else
				cavity_upper_faces[id] = true;
				cavity_stack.add(id);
			}
		}
	}

	assert(cavity_stack.countInt() < cfaces.countInt()); 

	DataVector<bool> proper_orientation(cfaces.countInt(), true);
	for(int i = 0; i < cfaces.countInt(); i++){
		MeshBlock* block = cfaces[i]->getBlock(0);
		if(!block || !hblocks.contains(block)) 
			proper_orientation[i] = false;
	}

	// --> check if the cavities are "watertight"
	{
		MeshContainer3dSurface upper_test_mesh(3*cfaces.countInt());
		MeshContainer3dSurface lower_test_mesh(3*cfaces.countInt());
		DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> upper_test_hpoints(3*cfaces.countInt(), nullptr);
		DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> lower_test_hpoints(3*cfaces.countInt(), nullptr);
		for(int i = 0; i < cfaces.countInt(); i++){
			MeshFace* face = cfaces[i];
			int fpct = face->getPointCount(); assert(fpct == 3);
			DataVector<MeshPoint3d*> test_face_points(fpct);
			for(int j = 0; j < fpct; j++){
				MeshPoint3d* point = face->getPoint(j);
				MeshPoint3d* test_point = nullptr;
				if(cavity_upper_faces[i]){
					test_point = upper_test_hpoints.getValue(point, nullptr);
					if(!test_point){
						upper_test_hpoints.insert( point, test_point = new MeshPoint3d(*point) );
						upper_test_mesh.addMeshPoint(test_point);
					}
				}else{
					test_point = lower_test_hpoints.getValue(point, nullptr);
					if(!test_point){
						lower_test_hpoints.insert( point, test_point = new MeshPoint3d(*point) );
						lower_test_mesh.addMeshPoint(test_point);
					}
				}
				test_face_points.add(test_point);
			}
			// orientation is not really important here...
			MeshFace* test_face = 
				new MeshTriangle3d(test_face_points[0], test_face_points[1], test_face_points[2]);
			if(cavity_upper_faces[i])
				upper_test_mesh.addMeshFace(test_face);
			else
				lower_test_mesh.addMeshFace(test_face);
		}
		DataVector<MeshPoint3d*> upper_test_face_points(3);
		DataVector<MeshPoint3d*> lower_test_face_points(3);
		for(int j = 0; j < 3; j++){
			MeshPoint3d* test_point = upper_test_hpoints.getValue(fpoints[j], nullptr);
			if(!test_point) return false;
			upper_test_face_points.add(test_point);
			test_point = lower_test_hpoints.getValue(fpoints[j], nullptr);
			if(!test_point) return false;
			lower_test_face_points.add(test_point);
		}
		upper_test_mesh.addMeshFace(new MeshTriangle3d(upper_test_face_points[0], 
			upper_test_face_points[1], upper_test_face_points[2]));
		lower_test_mesh.addMeshFace(new MeshTriangle3d(lower_test_face_points[0], 
			lower_test_face_points[1], lower_test_face_points[2]));

		// and now, check these surface meshes
		for(IteratorEdge3d it = upper_test_mesh.getFirstEdge3d(); it.isValid(); it.nextEdge())
			if(it.getEdge()->getFaceCount() != 2) return false;
		for(IteratorEdge3d it = lower_test_mesh.getFirstEdge3d(); it.isValid(); it.nextEdge())
			if(it.getEdge()->getFaceCount() != 2) return false;

		// plus, check for "lost" nodes
		for(int i = 0; i < crossing_blocks.countInt(); i++){
			MeshBlock* block = crossing_blocks[i];
			for(int j = 0; j < block->getPointCount(); j++){
				MeshPoint3d* point = block->getPoint(j);
				if( upper_test_hpoints.contains(point) ) continue; // ok
				if( lower_test_hpoints.contains(point) ) continue; // ok
				// else -> there is point within the blocks, that is not on the cavity boundary ...
				return false;
			}
		}
	}

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
	if(true){
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < cfaces.countInt(); i++){
			if(cavity_upper_faces[i])
				set->addFaceWithEdges(cfaces[i], -2, MeshViewSet::param_shrink, proper_orientation[i]);
		}
		set->addPoint(pt1, 3);
		set->addPoint(pt2, 3);
		set->addPoint(pt3, 3);

		SHOW_MESH("cavity upper faces", set);

		set = new MeshViewSet;
		for(int i = 0; i < cfaces.countInt(); i++){
			if(!cavity_upper_faces[i])
				set->addFaceWithEdges(cfaces[i], -2, MeshViewSet::param_shrink, proper_orientation[i]);
		}
		set->addPoint(pt1, 3);
		set->addPoint(pt2, 3);
		set->addPoint(pt3, 3);

		SHOW_MESH("cavity lower faces", set);
	}
#endif

	DataVector<DMTriangle3d> orient_faces_upper(cfaces.countInt() + 1);
	DataVector<DMTriangle3d> orient_faces_lower(cfaces.countInt() + 1);

	for(int i = 0; i < cfaces.countInt(); i++){
		DMTriangle3d dmt;
		if(proper_orientation[i])
			dmt = DMTriangle3d(
					cfaces[i]->getPoint(0)->getMetricCoordinates(mc),
					cfaces[i]->getPoint(1)->getMetricCoordinates(mc),
					cfaces[i]->getPoint(2)->getMetricCoordinates(mc));
		else
			dmt = DMTriangle3d(
					cfaces[i]->getPoint(1)->getMetricCoordinates(mc),
					cfaces[i]->getPoint(0)->getMetricCoordinates(mc),
					cfaces[i]->getPoint(2)->getMetricCoordinates(mc));
		if(cavity_upper_faces[i])
			orient_faces_upper.add(dmt);
		else
			orient_faces_lower.add(dmt);
	}

	// - add missing faces for upper/lower cavity (using edge info)
	bool local_orientation_upper, local_orientation_lower;
	bool upper_done = false;
	bool lower_done = false;
	for(int j = 0; j < cfaces.countInt(); j++){
		if(cavity_upper_faces[j]){
			if(upper_done) continue;
			MeshFace* face = cfaces[j];
			for(int k = 0; k < 3; k++){
				if(face->incidentToPoints(fpoints[k], fpoints[(k+1)%3])){
					bool oriented_ok = proper_orientation[j];
					if(face->properOrientation(fpoints[k], fpoints[(k+1)%3])){
						oriented_ok = !oriented_ok;
					}
					if(oriented_ok){
						orient_faces_upper.add( 
							DMTriangle3d(
								fpoints[0]->getMetricCoordinates(mc),
								fpoints[1]->getMetricCoordinates(mc),
								fpoints[2]->getMetricCoordinates(mc)));
					}else{
						orient_faces_upper.add( 
							DMTriangle3d(
								fpoints[1]->getMetricCoordinates(mc),
								fpoints[0]->getMetricCoordinates(mc),
								fpoints[2]->getMetricCoordinates(mc)));
					}
					local_orientation_upper = oriented_ok;
					upper_done = true;
				}
			}
		}else{
			if(lower_done) continue;
			MeshFace* face = cfaces[j];
			for(int k = 0; k < 3; k++){
				if(face->incidentToPoints(fpoints[k], fpoints[(k+1)%3])){
					bool oriented_ok = proper_orientation[j];
					if(face->properOrientation(fpoints[k], fpoints[(k+1)%3])){
						oriented_ok = !oriented_ok;
					}
					if(oriented_ok){
						orient_faces_lower.add( 
							DMTriangle3d(
								fpoints[0]->getMetricCoordinates(mc),
								fpoints[1]->getMetricCoordinates(mc),
								fpoints[2]->getMetricCoordinates(mc)));
					}else{
						orient_faces_lower.add( 
							DMTriangle3d(
								fpoints[1]->getMetricCoordinates(mc),
								fpoints[0]->getMetricCoordinates(mc),
								fpoints[2]->getMetricCoordinates(mc)));
					}
					local_orientation_lower = oriented_ok;
					lower_done = true;
				}
			}
		}
	}

	DMPoint3d center_upper, center_lower;
	if( findStarCenter(orient_faces_upper, center_upper) &&
		findStarCenter(orient_faces_lower, center_lower) )
	{

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
		MeshViewSet* set = new MeshViewSet;
		double min_volume = LARGE_NUMBER;
#endif

		// 4. Retriangulate cavities 
		//	- remove crossing blocks
		for(int i = 0; i < crossing_blocks.countInt(); i++)
			delete mesh->removeMeshBlock(crossing_blocks[i]);
		//	- insert new mesh points
		MeshPoint3d* mp_upper = new MeshPoint3d(mc.transformMStoRS(center_upper));
		mesh->addMeshPoint(mp_upper);
		MeshPoint3d* mp_lower = new MeshPoint3d(mc.transformMStoRS(center_lower));
		mesh->addMeshPoint(mp_lower);
		//	- insert new tetrahedra for cavity faces
		for(int i = 0; i < cfaces.countInt(); i++){
			MeshFace* face = cfaces[i];
			MeshTetrahedron* tetra = nullptr;
			if(cavity_upper_faces[i]){
				if(proper_orientation[i])
					tetra = new MeshTetrahedron(
						face->getPoint(0), face->getPoint(1), face->getPoint(2), mp_upper);
				else
					tetra = new MeshTetrahedron(
						face->getPoint(1), face->getPoint(0), face->getPoint(2), mp_upper);
			}else{
				if(proper_orientation[i])
					tetra = new MeshTetrahedron(
						face->getPoint(0), face->getPoint(1), face->getPoint(2), mp_lower);
				else
					tetra = new MeshTetrahedron(
						face->getPoint(1), face->getPoint(0), face->getPoint(2), mp_lower);
			}
			mesh->addMeshBlock(tetra);

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
			set->addBlockWithEdges(tetra);
			double v = tetra->getVolume(mc, false);
			if(v < min_volume) min_volume = v;
#endif

		}
		//	- reconstruct missing face
		MeshTetrahedron* tetra = nullptr;
		if(local_orientation_upper)
			tetra = new MeshTetrahedron(fpoints[0], fpoints[1], fpoints[2], mp_upper);
		else
			tetra = new MeshTetrahedron(fpoints[1], fpoints[0], fpoints[2], mp_upper);
		mesh->addMeshBlock(tetra);

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
		set->addBlockWithEdges(tetra);
		double v = tetra->getVolume(mc, false);
		if(v < min_volume) min_volume = v;
#endif

		if(local_orientation_lower)
			tetra = new MeshTetrahedron(fpoints[0], fpoints[1], fpoints[2], mp_lower);
		else
			tetra = new MeshTetrahedron(fpoints[1], fpoints[0], fpoints[2], mp_lower);
		mesh->addMeshBlock(tetra);

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
		set->addBlockWithEdges(tetra);
		v = tetra->getVolume(mc, false);
		if(v < min_volume) min_volume = v;

		LOG4CPLUS_DEBUG(MeshLog::logger_console, "New tetrahedra min volume", min_volume);
		SHOW_MESH("New tetrahedra", set);
#endif

		return true;

	}

	return false;

	//  (b) greedy frontal approach?
}

/// Find good center of star-shaped cavity (if possible)
bool MeshGenerator3dDelaunayBoundary::findStarCenter(const DataVector<DMTriangle3d> & orient_faces, DMPoint3d& center)
{
	if(orient_faces.empty()) return false;

	int ofct = orient_faces.countInt();
	DataVector<double> orient_faces_det(ofct);
	for(int i = 0; i < ofct; i++)
		orient_faces_det.add(abs(orient_faces[i].det()));

	DBox ubox;
	for(int i = 0; i < ofct; i++){
		ubox.addPoint(orient_faces[i].pt_a);
		ubox.addPoint(orient_faces[i].pt_b);
		ubox.addPoint(orient_faces[i].pt_c);
	}

	static const int MAX_RES = 15;
	static const int MAX_SIZE = 4 * MAX_RES * MAX_RES * MAX_RES;

	double dx = ubox.getDX() / MAX_RES;
	double dy = ubox.getDY() / MAX_RES;
	double dz = ubox.getDZ() / MAX_RES;

	DataVector<PointNode> cells(MAX_SIZE);
	PointNode cell(DMPoint3d(0.0, 0.0, ubox.z0 + 0.5*dz));
	for(int i = 0; i < MAX_RES; i++){
		cell.pt.y = ubox.y0 + 0.5*dy;
		for(int j = 0; j < MAX_RES; j++){
			cell.pt.x = ubox.x0 + 0.5*dx;
			for(int k = 0; k < MAX_RES; k++){
				cells.add(cell);
				cell.pt.x += dx;
			}
			cell.pt.y += dy;
		}
		cell.pt.z += dz;
	}
	
	dx *= 0.5;
	dy *= 0.5;
	dz *= 0.5;
	double dr = std::max(std::max(dx, dy), dz);

	int max_steps = 10;
	while(max_steps--){
		// calculate visibility
		for(int i = 0; i < cells.countInt(); i++){
			PointNode & pn = cells[i];
			for(int j = 0; j < ofct; j++){
				double h = orient_faces[j].orient3d(cells[i].pt) / orient_faces_det[j];
				if(h < -dr){
					pn.visibility = 0; break;
				}else if (h < dr){
					double v = 0.5*(dr+h)/dr;
					if(v < pn.visibility) pn.visibility = v;
				}
			}
		}

#ifdef SHOW_RECOVER_BOUNDARY_CAVITY
		if(true){
			MeshViewSet* set = new MeshViewSet;
			for(int i = 0; i < ofct; i++){
				set->addFace(
					orient_faces[i].pt_a.toRealSpace(),
					orient_faces[i].pt_b.toRealSpace(),
					orient_faces[i].pt_c.toRealSpace() );
			}

			for(int i = 0; i < cells.countInt(); i++){
				if(cells[i].visibility > 0.0){
					if(cells[i].visibility < 1.0)
						set->addPoint(cells[i].pt.toRealSpace(), 1);
					else
						set->addPoint(cells[i].pt.toRealSpace(), 2);
				}else
					set->addPoint(cells[i].pt.toRealSpace());
			}

			SHOW_MESH("orient faces + cells", set);
		}
#endif

		// - if any 1, first calculate barycenter of all 1-cells, then select 1-cell closest to barycenter
		DMPoint3d barycenter;
		int w_count = 0;
		for(int i = 0; i < cells.countInt(); i++){
			if(cells[i].visibility == 1.0){
				barycenter.add(cells[i].pt);
				w_count++;
			}
		}

		if(w_count > 0){
			barycenter /= w_count;
			int best_i = -1;
			double min_dist2 = 0.0;
			for(int i = 0; i < cells.countInt(); i++){
				if(cells[i].visibility == 1){
					double dist2 = barycenter.distance2(cells[i].pt);
					if(best_i < 0 || dist2 < min_dist2){
						best_i = i;
						min_dist2 = dist2;
					}
				}
			}
			center = cells[best_i].pt;
			return true;
		}

		// ... else, remove all 0-cells
		for(int i = 0; i < cells.countInt(); i++){
			if(cells[i].visibility == 0.0){
				cells.removeAt(i);
				i--;
			}
		}
		// - if no hope...
		if(cells.empty()){
			LOG4CPLUS_WARN(MeshLog::logger_console, "Failed searching for star-center, all cells are 0");
			return false;
		}

		// - replace cells with lower resolution
		// -- clear low-visibility cells if necessary (TODO better!)
		if(14*cells.countInt() > MAX_SIZE){
			std::sort( cells.begin(), cells.end() );
			while(14*cells.countInt() > MAX_SIZE){
				// - remove last (which is lowest, after sort)
				cells.removeLast();
			}
		}
		// -- split (it's easier to do from the end
		static const double DRATIO = 0.67;
		dr *= DRATIO;
		dx *= DRATIO;
		dy *= DRATIO;
		dz *= DRATIO;
		for(int i = cells.countInt()-1; i >=0; i--){
			DMPoint3d pt = cells[i].pt;
			cells.removeAt(i);
			cells.add(PointNode(DMPoint3d(pt.x + dx, pt.y + dy, pt.z + dz)));
			cells.add(PointNode(DMPoint3d(pt.x - dx, pt.y + dy, pt.z + dz)));
			cells.add(PointNode(DMPoint3d(pt.x + dx, pt.y - dy, pt.z + dz)));
			cells.add(PointNode(DMPoint3d(pt.x - dx, pt.y - dy, pt.z + dz)));
			cells.add(PointNode(DMPoint3d(pt.x + dx, pt.y + dy, pt.z - dz)));
			cells.add(PointNode(DMPoint3d(pt.x - dx, pt.y + dy, pt.z - dz)));
			cells.add(PointNode(DMPoint3d(pt.x + dx, pt.y - dy, pt.z - dz)));
			cells.add(PointNode(DMPoint3d(pt.x - dx, pt.y - dy, pt.z - dz)));
			cells.add(PointNode(DMPoint3d(pt.x + dx, pt.y, pt.z)));
			cells.add(PointNode(DMPoint3d(pt.x - dx, pt.y, pt.z)));
			cells.add(PointNode(DMPoint3d(pt.x, pt.y + dy, pt.z)));
			cells.add(PointNode(DMPoint3d(pt.x, pt.y - dy, pt.z)));
			cells.add(PointNode(DMPoint3d(pt.x, pt.y, pt.z + dz)));
			cells.add(PointNode(DMPoint3d(pt.x, pt.y, pt.z - dz)));
		}
	}

	LOG4CPLUS_WARN(MeshLog::logger_console, "Failed searching for star-center, reached max-steps");
	return false;
}
