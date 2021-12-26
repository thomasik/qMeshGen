// MeshGenerator3dAdapt.cpp: implementation of the MeshGenerator3dAdapt class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshGenerator3dAdapt.h"

#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshTetrahedron.h"
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
#include "ControlSpace3dAdaptive.h"
#include "ControlSpace3dOctree.h"
#include "MeshGenerator3d.h"
#include "MeshGenerator3dSurface.h"
#include "MeshGenerator3dQuality.h"
#include "MeshGenerator3dDelaunayBoundary.h"
#include "MeshSplit3d.h"
#include "MeshViewSet.h"
#include "SurfacePlane.h"
#include "SurfacePlanarQuadric.h"
#include "SurfaceBSplinePlanar.h"
#include "IteratorEdge3d.h"
#include "IteratorFace.h"
#include "ControlSpace2dAdaptive.h"
#include "GeometricPredicates.h"
#include "DTetrahedron.h"
#include "DTriangle.h"
#include "DPlane.h"
#include "DPlanarQuadric.h"
#include "DLeastSquaresFitting.h"

// #define STAT_MESH

/// Modify tetrahedral mesh according to the given CS3d, with local operators (edge collapsing and smoothing)
bool MeshGenerator3dAdapt::remeshWithLocalTransformations(Metric3dContext& mc, MeshContainer3d* mesh)
{
	START_CLOCK("MeshGenerator3dAdapt::remeshWithLocalTransformations");

	START_CLOCK("TRANSFORM-ADAPT:mark");
	DataSimpleList< DataVector<MeshBlock*> > adapt_regions;
	int adapt_count = MeshGenerator3dAdapt::findAdaptRegions(mc, mesh, adapt_regions, 0.01, 10, 1);

	if(adapt_count == 0) return false;
	LOG4CPLUS_INFO(MeshLog::logger_console, "Adapt regions: count=" << adapt_regions.countInt());

	// mark nodes ...
	DataVector<bool> active_points(mesh->getPointsCount(), false);
	for(auto it = adapt_regions.iterator(); it.valid(); it.moveNext()){
		const DataVector<MeshBlock*> rblocks = it.item();
		int rbct = rblocks.countInt();
		for(int i = 0; i < rbct; i++){
			MeshBlock* block = rblocks[i];
			int bpct = block->getPointCount();
			for(int j = 0; j < bpct; j++){
				MeshPoint3d* point = block->getPoint(j);
				int pid = point->getIndex();
				if(!active_points[pid]){
					active_points[pid] = true;
					point->setIntTag(TagExtended::TAG_ADAPT_VOLUME, 1); // same tag for all regions
				}
			}
		}
	}
	STOP_CLOCK("TRANSFORM-ADAPT:mark");

#ifdef STAT_MESH
	MeshGenerator3dQuality::statMesh(mc, mesh, "START");
#endif

//	mesh->clearLocalSurfaces();
//	MeshGenerator3dAdapt::identifyLocalSurfaces(mc, mesh, 0.5, TagExtended::TAG_ADAPT_CRACK, 1);
//	MeshGenerator3dAdapt::identifyLocalSurfaces(mc, mesh, 0.2, TagExtended::TAG_ADAPT_VOLUME, 1);

	//if(true){
	//	MeshViewSet* set = new MeshViewSet;
	//	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
	//		MeshFace* face = it.getFace();
	//		if(face->hasLocalSurface())
	//			set->addFaceWithEdges(face);
	//	}
	//	SHOW_MESH("Local surfaces - init", set);
	//}

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
		"Start local-transform: NT=" << mesh->getBlocksCount()
		<< ", NP=" << mesh->getPointsCount());

//	SHOW_MESH("Volume mesh, init remeshing", mesh->getViewSet());

	START_CLOCK("TRANSFORM-ADAPT:collapse");
	int iter = 0;
	static const int MIN_RE = 5;
	static const int MAX_ITER = 10;
	int re = MIN_RE+1; // removed edges during remeshing
	while(re > MIN_RE && iter < MAX_ITER){
		iter++;

		int ire = MeshGenerator3dAdapt::collapseInnerEdges(mc, mesh, TagExtended::TAG_ADAPT_VOLUME, 1);
#ifdef STAT_MESH
		MeshGenerator3dQuality::statMesh(mc, mesh, "AFTER-COLLAPSE-INNER");
#endif

		int bre = MeshGenerator3dAdapt::collapseBoundaryEdges(mc, mesh, TagExtended::TAG_ADAPT_VOLUME, 1);
		re = ire + bre;
#ifdef STAT_MESH
		MeshGenerator3dQuality::statMesh(mc, mesh, "AFTER-COLLAPSE-BORDER");
#endif

		//if(true){
		//	MeshViewSet* set = new MeshViewSet;
		//	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
		//		MeshFace* face = it.getFace();
		//		if(face->hasLocalSurface())
		//			set->addFaceWithEdges(face);
		//	}
		//	SHOW_MESH("Local surfaces - bc", set);
		//}

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
			iter << ". After collapse (" << ire << "ie + " << bre << "be)"
			<< ": NT=" << mesh->getBlocksCount() << ", NP=" << mesh->getPointsCount());

//		SHOW_MESH_NORESET("Volume mesh after collapse", mesh->getViewSet());

		MeshGenerator3dQuality::smoothen(mc, mesh, 1, TagExtended::TAG_ADAPT_VOLUME, 1, true);

		int bs = MeshGenerator3dQuality::removeBoundarySlivers(mc, mesh, 0.2, TagExtended::TAG_ADAPT_VOLUME, 1);

		//if(true){
		//	MeshViewSet* set = new MeshViewSet;
		//	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
		//		MeshFace* face = it.getFace();
		//		if(face->hasLocalSurface())
		//			set->addFaceWithEdges(face);
		//	}
		//	SHOW_MESH("Local surfaces - sm1", set);
		//}

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
			iter << ". After smoothen"
			<< ": NT=" << mesh->getBlocksCount() << ", NP=" << mesh->getPointsCount());

#ifdef STAT_MESH
		MeshGenerator3dQuality::statMesh(mc, mesh, "AFTER-SM1");
#endif
		//if(true){
		//	ostringstream ostr;
		//	ostr << "Volume mesh after iteration #" << iter << ", removed edges: " << re;
		//	SHOW_MESH_NORESET(ostr.str(), mesh->getViewSet());
		//}

	}

	MeshGenerator3dQuality::smoothen(mc, mesh, 1, TagExtended::TAG_ADAPT_VOLUME, 1, true,
		MeshData::SM3_LAPLACE_MIXED | MeshData::SM3_OPT_SIMPLE | MeshData::SM3_SWAP_COMPLETE | MeshData::SM3_BORDER_PRUNE);
#ifdef STAT_MESH
	MeshGenerator3dQuality::statMesh(mc, mesh, "AFTER-SM2");
#endif
	STOP_CLOCK("TRANSFORM-ADAPT:collapse");
//	SHOW_MESH_NORESET("Volume mesh after collapse, final", mesh->getViewSet());

	//if(true){
	//	MeshViewSet* set = new MeshViewSet;
	//	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
	//		MeshFace* face = it.getFace();
	//		if(face->hasLocalSurface())
	//			set->addFaceWithEdges(face);
	//	}
	//	SHOW_MESH("Local surfaces - sm2", set);
	//}

	int ins = MeshGenerator3d::addInnerNodes(mc, mesh, TagExtended::TAG_ADAPT_VOLUME, 1);
	int bns = MeshGenerator3dAdapt::addBoundaryNodesSplit(mc, mesh, TagExtended::TAG_ADAPT_VOLUME, 1);

	//if(true){
	//	MeshViewSet* set = new MeshViewSet;
	//	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
	//		MeshFace* face = it.getFace();
	//		if(face->hasLocalSurface())
	//			set->addFaceWithEdges(face);
	//	}
	//	SHOW_MESH("Local surfaces - bns", set);
	//}

//	SHOW_MESH_NORESET("Volume mesh after split", mesh->getViewSet());

//	int ins = MeshGenerator3d::addInnerNodes(mc, mesh, TagExtended::TAG_ADAPT_VOLUME, 1);

#ifdef STAT_MESH
	MeshGenerator3dQuality::statMesh(mc, mesh, "AFTER-ADD-INNER");
#endif
	MeshGenerator3dQuality::smoothen(mc, mesh, 2, TagExtended::TAG_ADAPT_VOLUME, 1, true,
		MeshData::SM3_LAPLACE_MIXED | MeshData::SM3_OPT_SIMPLE | MeshData::SM3_SWAP_COMPLETE | MeshData::SM3_BORDER_PRUNE);
	//MeshGenerator3dQuality::smoothen(mc, mesh, 1, TagExtended::TAG_ADAPT_VOLUME, 1, true,
	//	MeshData::SM3_OPT_MIXED | MeshData::SM3_SWAP_COMPLETE | MeshData::SM3_BORDER_PRUNE);
#ifdef STAT_MESH
	MeshGenerator3dQuality::statMesh(mc, mesh, "AFTER-SM3");
#endif

	STOP_CLOCK("MeshGenerator3dAdapt::remeshWithLocalTransformations");

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
		"After split (" << bns << "bp + " << ins << "ip)"
		<< ": NT=" << mesh->getBlocksCount() << ", NP=" << mesh->getPointsCount());

	//SHOW_MESH_NORESET("Volume mesh after collapse+split, final", mesh->getViewSet());

	return true;
}

/// Modify tetrahedral mesh according to the given CS3d, with local operators (edge collapsing and smoothing)
bool MeshGenerator3dAdapt::remeshWithLocalCutRegions(Metric3dContext& mc, MeshContainer3d* mesh)
{
	START_CLOCK("MeshGenerator3dAdapt::remeshWithLocalCutRegions");

	DataSimpleList< DataVector<MeshBlock*> > adapt_regions;
	int adapt_count = MeshGenerator3dAdapt::findAdaptRegions(mc, mesh, adapt_regions, 0.01, 100, 1);

	if(adapt_count == 0) return false;

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "Adapt regions count=" << adapt_regions.countInt());

//	MeshViewSet* set = new MeshViewSet(0, 2*adapt_count, 0, adapt_count);
	int counter = 0;
	for(auto it = adapt_regions.iterator(); it.valid(); it.moveNext()){
		const DataVector<MeshBlock*> rblocks = it.item();
		counter++;
		for(int i = 0; i < rblocks.countInt(); i++){
			rblocks[i]->setIntTag(TagExtended::TAG_CUT_ADAPT, counter);
//			set->addBlockWithEdges(rblocks[i], counter);
		}
	}
//	SHOW_MESH("blocks to adapt", set);

	for(int i = 0; i < counter; i++){
		MeshContainer3d* cut_mesh = MeshSplit3d::splitMeshByBlocks(mesh, TagExtended::TAG_CUT_ADAPT, i+1);

		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Region to adapt, NT=" << cut_mesh->getBlocksCount());
		//MeshViewSet* set = new MeshViewSet;
		//for(IteratorEdge3d it = cut_mesh->getFirstEdge3d(); it.isValid(); it.nextEdge())
		//	if(it.getEdge()->isBorder())
		//		set->addEdge(it.getEdge(), -1);
		//int bct = cut_mesh->getBlocksCount();
		//for(int j = 0; j < bct; j++){
		//	MeshTetrahedron* t = (MeshTetrahedron*) cut_mesh->getBlockAt(j);
		//	int bfct = 0;
		//	for(int k = 0; k < t->getFaceCount(); k++)
		//		if(t->getFace(k)->isBorder()) bfct++;
		//	if(bfct < 2) continue;
		//	double q = t->getMeanRatio(mc);
		//	if(q < 0.2)
		//		set->addBlock(t);
		//}
		//SHOW_MESH("cut mesh", set);


		cut_mesh->setControlSpace(mesh->getControlSpace());
		MeshDomainVolume* cut_mdv = new MeshDomainVolume(cut_mesh);

		MeshContainer3dSurface* mesh_surface = cut_mdv->cutSurfaceMeshFromVolumeMesh();
//		SHOW_MESH("Surface mesh", mesh_surface->getViewSet());

		// mark elements of mesh_surface
		int fct = mesh_surface->getFacesCount();
		for(int i = 0; i < fct; i++){
			MeshFace* face = mesh_surface->getFaceAt(i);
			int face_tag = face->getIntTag(TagExtended::TAG_ADAPT_SURF, -1);
			if(face_tag < 0){
				face_tag = face->availableTag(TagExtended::TAG_CUT_3D) ? 4 : 2;
				// 4 -> adaptation forbidden
				// 2 -> adaptation allowed, but separately from "1"
				face->setIntTag(TagExtended::TAG_ADAPT_SURF, face_tag); 
			}else{
				assert(face_tag == 1); // should be "1" for crack-surface
			}
			// + mark points
			int fpct = face->getPointCount();
			for(int j = 0; j < fpct; j++)
				face->getPoint(j)->setIntFlag(TagExtended::TAG_ADAPT_SURF, face_tag);
		}

		// ... mark boundaries for suface mesh
		mesh_surface->clearBoundaryFlags();
		mesh_surface->setBoundaryFeatureEdges();
		mesh_surface->setBoundaryTagEdges(TagExtended::TAG_ADAPT_SURF);
		mesh_surface->setBoundarySharpEdges(MeshGenerator3dSurface::param_sharp_edge_threshold, 
			TagExtended::TAG_ADAPT_SURF, 2); // only for the "oryginal" surface of the model

		mesh_surface->clearLocalShapes();
		//mesh_surface->identifyLocalSurfaces(mc, 0.5, TagExtended::TAG_ADAPT_SURF, 1);
		//mesh_surface->identifyLocalSurfaces(mc, 0.2, TagExtended::TAG_ADAPT_SURF, 2);
		MeshGenerator3dSurface::param_local_shape_tolerance = 0.5;
		MeshGenerator3dSurface::identifyLocalSurfaces(mc, mesh_surface, TagExtended::TAG_ADAPT_SURF, 1);
		MeshGenerator3dSurface::param_local_shape_tolerance = 0.2;
		MeshGenerator3dSurface::identifyLocalSurfaces(mc, mesh_surface, TagExtended::TAG_ADAPT_SURF, 2);
		
		MeshGenerator3dSurface::remeshCrackSurfaceMeshWithLocalTransformations(mc, mesh_surface);

		START_CLOCK("MSR-UTC:meshing3dregion");
		if(cut_mdv->discretizeUsingBoundary(mc) > 0){
			MeshContainer3d *new_cut_mesh = cut_mdv->getMesh();
			if(new_cut_mesh){
				MeshGenerator3d::addInnerNodes(mc, new_cut_mesh);
				MeshGenerator3dQuality::smoothen(mc, new_cut_mesh, 2);

//				SHOW_MESH("Remeshed volume", new_cut_mesh->getViewSet());
				LOG4CPLUS_DEBUG(MeshLog::logger_console, 
					"Region adapted, NT=" << new_cut_mesh->getBlocksCount());
			}
		}
		STOP_CLOCK("MSR-UTC:meshing3dregion");

		MeshContainer3d* new_cut_mesh = cut_mdv->removeMesh();
		// Remove boundary-pointers from tags, since domain_volume will be deleted shortly
		// ... also copy cut-data
		for(int j = 0; j < new_cut_mesh->getPointsCount(); j++){
			MeshPoint3d* point = new_cut_mesh->getPointAt(j);
			if(!point->isBorder()){
				assert(!point->availableTag(TagExtended::TAG_BOUNDARY_POINT));
				continue;
			}
			MeshPoint3d* bpoint = (MeshPoint3d*) point->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
			if(!bpoint) continue;
			point->removeTag(TagExtended::TAG_BOUNDARY_POINT);
			if(bpoint->availableTag(TagExtended::TAG_CUT_MP_3D))
				point->setPtrTag(TagExtended::TAG_CUT_MP_3D, bpoint->getPtrTag(TagExtended::TAG_CUT_MP_3D) );
			if(bpoint->availableTag(TagExtended::TAG_CUT_ORYG_DATA))
				point->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, bpoint->getIntTag(TagExtended::TAG_CUT_ORYG_DATA) );
		}
		for(IteratorFace it = new_cut_mesh->getFirstFace(); it.isValid(); it.nextFace()){
			MeshFace* face = it.getFace();
			if(!face->isBorder()){
				if(face->availableTag(TagExtended::TAG_BOUNDARY_FACE)){
					MeshViewSet * set = new MeshViewSet(0, 6*new_cut_mesh->getBlocksCount(), 3*new_cut_mesh->getBlocksCount());
					for(IteratorFace itx = new_cut_mesh->getFirstFace(); itx.isValid(); itx.nextFace()){
						if(itx.getFace()->isBorder())
							set->addFaceWithEdges(itx.getFace(), -2, MeshViewSet::param_shrink, 
									itx.getFace()->getBlock(0) != nullptr);
					}
					SHOW_MESH("Border face consistency?", set);
					//assert(!face->availableTag(TagExtended::TAG_BOUNDARY_FACE));
					continue;
				}
			}
			MeshFace* bface = (MeshFace*) face->getPtrTag(TagExtended::TAG_BOUNDARY_FACE);
			if(!bface) continue;
			face->removeTag(TagExtended::TAG_BOUNDARY_FACE);
			if(bface->availableTag(TagExtended::TAG_CUT_MF_3D))
				face->setPtrTag(TagExtended::TAG_CUT_MF_3D, bface->getPtrTag(TagExtended::TAG_CUT_MF_3D) );
			if(bface->availableTag(TagExtended::TAG_CUT_ORYG_DATA))
				face->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, bface->getIntTag(TagExtended::TAG_CUT_ORYG_DATA) );
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Merging cut mesh");
		MeshSplit3d::mergeMeshes(mesh, new_cut_mesh);
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Merging done");
//		delete new_cut_mesh;
		delete cut_mdv;
	}

	mesh->markOuterBoundary();

	STOP_CLOCK("MeshGenerator3dAdapt::remeshWithLocalCutRegions");
	return true;
}

/*  VERSION ORIGINAL ( -> ComPIegne, 30.10.2010)
 *
/// Modify tetrahedral mesh according to the given CS3d, with local operators (edge collapsing and smoothing)
bool MeshGenerator3dAdapt::remeshWithLocalCutRegions(Metric3dContext& mc, MeshContainer3d* mesh)
{
	START_CLOCK("MeshGenerator3dAdapt::remeshWithLocalCutRegions");

	MeshGenerator3d::addInnerNodes(mc, mesh);
	MeshGenerator3dQuality::smoothen(mc, mesh, 2);

	DataSimpleList< DataVector<MeshBlock*> > adapt_regions;
	int adapt_count = MeshGenerator3dAdapt::findAdaptRegions(mc, mesh, adapt_regions, 0.01, 100, 1);

	if(adapt_count == 0) return;

	LOG4CPLUS_INFO(MeshLog::logger_console, "Adapt regions", adapt_regions.countInt());

//	MeshViewSet* set = new MeshViewSet(0, 2*adapt_count, 0, adapt_count);
	int counter = 0;
	for(auto it = adapt_regions.iterator(); it.valid(); it.moveNext()){
		const DataVector<MeshBlock*> rblocks = it.item();
		counter++;
		for(int i = 0; i < rblocks.countInt(); i++){
			rblocks[i]->setIntTag(TagExtended::TAG_CUT_ADAPT, counter);
//			set->addBlockWithEdges(rblocks[i], counter);
		}
	}
//	SHOW_MESH("blocks to adapt", set);

	for(int i = 0; i < counter; i++){
		MeshContainer3d* cut_mesh = MeshSplit3d::splitMeshByBlocks(mesh, TagExtended::TAG_CUT_ADAPT, i+1);

		LOG4CPLUS_INFO(MeshLog::logger_console, "Region to adapt, tetrahedra", cut_mesh->getBlocksCount());
		//MeshViewSet* set = new MeshViewSet;
		//for(IteratorEdge3d it = cut_mesh->getFirstEdge3d(); it.isValid(); it.nextEdge())
		//	if(it.getEdge()->isBorder())
		//		set->addEdge(it.getEdge(), -1);
		//int bct = cut_mesh->getBlocksCount();
		//for(int j = 0; j < bct; j++){
		//	MeshTetrahedron* t = (MeshTetrahedron*) cut_mesh->getBlockAt(j);
		//	int bfct = 0;
		//	for(int k = 0; k < t->getFaceCount(); k++)
		//		if(t->getFace(k)->isBorder()) bfct++;
		//	if(bfct < 2) continue;
		//	double q = t->getMeanRatio(mc);
		//	if(q < 0.2)
		//		set->addBlock(t);
		//}
		//SHOW_MESH("cut mesh", set);


		cut_mesh->setControlSpace(cs_next);
		MeshDomainVolume* cut_mdv = new MeshDomainVolume(cut_mesh);
		//m_model_mesh->addMeshBlock(cut_mdv);

//		MeshGenerator3d::remeshWithLocalTransformations(mc, cut_mesh);

		MeshContainer3dSurface* mesh_surface = cut_mdv->cutSurfaceMeshFromVolumeMesh();
//		SHOW_MESH("Surface mesh", mesh_surface->getViewSet());

		// mark elements of mesh_surface
		int fct = mesh_surface->getFacesCount();
		for(int i = 0; i < fct; i++){
			MeshFace* face = mesh_surface->getFaceAt(i);
			int face_tag = face->getIntTag(TagExtended::TAG_ADAPT_SURF, -1);
			if(face_tag < 0){
				face_tag = face->availableTag(TagExtended::TAG_CUT_3D) ? 4 : 2;
				// 4 -> adaptation forbidden
				// 2 -> adaptation allowed, but separately from "1"
				face->setIntTag(TagExtended::TAG_ADAPT_SURF, face_tag); 
			}else{
				assert(face_tag == 1); // should be "1" for crack-surface
			}
			// + mark points
			int fpct = face->getPointCount();
			for(int j = 0; j < fpct; j++)
				face->getPoint(j)->setIntFlag(TagExtended::TAG_ADAPT_SURF, face_tag);
		}

		// ... mark boundaries for suface mesh
		mesh_surface->clearBoundaryFlags();
		mesh_surface->setBoundaryFeatureEdges();
		mesh_surface->setBoundaryTagEdges(TagExtended::TAG_ADAPT_SURF);
		mesh_surface->setBoundarySharpEdges(MeshGenerator3dSurface::param_sharp_edge_threshold, 
			TagExtended::TAG_ADAPT_SURF, 2); // only for the "oryginal" surface of the model

		mesh_surface->clearLocalSurfaces();
		mesh_surface->identifyLocalSurfaces(mc, 0.5, TagExtended::TAG_ADAPT_SURF, 1);
		mesh_surface->identifyLocalSurfaces(mc, 0.2, TagExtended::TAG_ADAPT_SURF, 2);

		MeshGenerator3dSurface::remeshSurfaceMeshWithLocalTransformations(mc, mesh_surface);

		START_CLOCK("MSR-UTC:meshing3dregion");
		if(cut_mdv->discretizeUsingBoundary(mc) > 0){
			MeshContainer3d *new_cut_mesh = cut_mdv->getMesh();
			if(new_cut_mesh){
				MeshGenerator3d::addInnerNodes(mc, new_cut_mesh);
				MeshGenerator3dQuality::smoothen(mc, new_cut_mesh, 2);

//				SHOW_MESH("Remeshed volume", new_cut_mesh->getViewSet());
				LOG4CPLUS_INFO(MeshLog::logger_console, "Region adapted, tetrahedra", new_cut_mesh->getBlocksCount());
			}
		}
		STOP_CLOCK("MSR-UTC:meshing3dregion");

		MeshContainer3d* new_cut_mesh = cut_mdv->removeMesh();
		// Remove boundary-pointers from tags, since domain_volume will be deleted shortly
		// ... also copy cut-data
		for(int j = 0; j < new_cut_mesh->getPointsCount(); j++){
			MeshPoint3d* point = new_cut_mesh->getPointAt(j);
			if(!point->isBorder()){
				assert(!point->availableTag(TagExtended::TAG_BOUNDARY_POINT));
				continue;
			}
			MeshPoint3d* bpoint = (MeshPoint3d*) point->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
			if(!bpoint) continue;
			point->removeTag(TagExtended::TAG_BOUNDARY_POINT);
			if(bpoint->availableTag(TagExtended::TAG_CUT_MP_3D))
				point->setPtrTag(TagExtended::TAG_CUT_MP_3D, bpoint->getPtrTag(TagExtended::TAG_CUT_MP_3D) );
			if(bpoint->availableTag(TagExtended::TAG_CUT_ORYG_DATA))
				point->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, bpoint->getIntTag(TagExtended::TAG_CUT_ORYG_DATA) );
		}
		for(IteratorFace it = new_cut_mesh->getFirstFace(); it.isValid(); it.nextFace()){
			MeshFace* face = it.getFace();
			if(!face->isBorder()){
				if(face->availableTag(TagExtended::TAG_BOUNDARY_FACE)){
					MeshViewSet * set = new MeshViewSet(0, 6*new_cut_mesh->getBlocksCount(), 3*new_cut_mesh->getBlocksCount());
					for(IteratorFace itx = new_cut_mesh->getFirstFace(); itx.isValid(); itx.nextFace()){
						if(itx.getFace()->isBorder())
							set->addFaceWithEdges(itx.getFace(), -2, MeshViewSet::param_shrink, 
									itx.getFace()->getBlock(0) != nullptr);
					}
					SHOW_MESH("Border face consistency?", set);
					//assert(!face->availableTag(TagExtended::TAG_BOUNDARY_FACE));
					continue;
				}
			}
			MeshFace* bface = (MeshFace*) face->getPtrTag(TagExtended::TAG_BOUNDARY_FACE);
			if(!bface) continue;
			face->removeTag(TagExtended::TAG_BOUNDARY_FACE);
			if(bface->availableTag(TagExtended::TAG_CUT_MF_3D))
				face->setPtrTag(TagExtended::TAG_CUT_MF_3D, bface->getPtrTag(TagExtended::TAG_CUT_MF_3D) );
			if(bface->availableTag(TagExtended::TAG_CUT_ORYG_DATA))
				face->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, bface->getIntTag(TagExtended::TAG_CUT_ORYG_DATA) );
		}
		LOG4CPLUS_INFO(MeshLog::logger_console, "Merging cut mesh");
		MeshSplit3d::mergeMeshes(mesh, new_cut_mesh);
		LOG4CPLUS_INFO(MeshLog::logger_console, "Done");
//		delete new_cut_mesh;
		delete cut_mdv;
	}

	mesh->markOuterBoundary();

	STOP_CLOCK("MeshGenerator3dAdapt::remeshWithLocalCutRegions");
	return true;
}
*/

/// Create new adaptive control space with sizing taken from mesh blocks
CS3dAPtr MeshGenerator3dAdapt::createACSfromMeshBlocks(MeshContainer3d* mesh)
{
	if(!mesh) return nullptr;
	int pct = mesh->getPointsCount();
	if(pct == 0) return nullptr;

	DBox box;
	for(int i = 0; i < pct; i++)
		box.addPoint(mesh->getPointAt(i)->getCoordinates());

	box.inflate(ControlSpace2d::param_inflate_box_factor);
	auto space = std::make_shared<ControlSpace3dOctree>(box);
	space->setMaxMetric();
	if(space->addSimpleMeshControlData(mesh))
		space->smoothen();

	return space;
}

//#define SHOW_COLLAPSE

/// Collapse too short edges, return number of removed edges
MeshPoint3d* MeshGenerator3dAdapt::collapseInnerEdge(Metric3dContext& mc, MeshContainer3d* mesh,
		MeshEdge3d* edge, DataContainer<MeshEdge3d::ActiveEdge> * active_edges)
{
	if(!edge || edge->isBorder()) return nullptr;

	MeshPoint3d* p0 = edge->getMeshPoint(0);
	MeshPoint3d* p1 = edge->getMeshPoint(1);

	if(p0->isBorder() && p1->isBorder()) // one at most
		return nullptr;

	int fct = edge->getFaceCount();
	DataVector<MeshFace*> faces(fct);
	DataVector<MeshPoint3d*> points(fct);
	for(int i = 0; i < fct; i++){
		MeshFace* face = edge->getFaceAt(i);
		faces.add(face);
		points.add(face->getOtherPoint(p0, p1));
	}

#ifdef SHOW_COLLAPSE
	if(true){
		ostringstream ostr;
		ostr << "Selected edge to collapse ";
		ostr << ", p0 - " << (p0->isBorder() ? "border" : "inner");
		ostr << ", p1 - " << (p1->isBorder() ? "border" : "inner");
		SHOW_MESH(ostr.str(), mesh->getDebugViewSet(p0, p1, 4.0));
		LOG4CPLUS_INFO(MeshLog::logger_mesh, ostr.str());
	}
#endif

	DataVector<MeshBlock*> blocks(fct);
	edge->adjacentBlocks(blocks);

	for(int i = 0; i < blocks.countInt(); i++)
		blocks[i]->detachFromFaces();
	DPoint3d old_pt;
	DataVector<MeshBlock*> switched_blocks(10);

	if(p1->isBorder()){
		MeshPoint3d* temp_p = p1; p1 = p0; p0 = temp_p;
	}

	p1->adjacentBlocks(switched_blocks, false);

	// check extra conditions
	bool valid = true;
	DataHashTable<MeshEdge3d*> check_edges((unsigned int)(2*switched_blocks.countInt()), nullptr);
	for(int i = 0; valid && i < switched_blocks.countInt(); i++){
		MeshBlock* sblock = switched_blocks[i];
		for(int j = 0; valid && j < sblock->getEdgeCount(); j++){
			MeshEdge3d* edge = sblock->getEdge(j);
			if( edge->getMeshPoint(0)->getEdgeToPoint(p0) != nullptr &&
				edge->getMeshPoint(1)->getEdgeToPoint(p0) != nullptr &&
				!check_edges.insert(edge))
			{
				valid = false;
			}
		}
	}

	if(!valid){
		for(int i = 0; i < blocks.countInt(); i++)
			blocks[i]->attachToFaces();
		return nullptr;
	}

	if(p0->isBorder()){ // --> move towards p0
		if(active_edges){
			// remove tag from edges which will be removed
			for(int i = 0; i < p1->getRank(); i++){
				MeshEdge3d* temp_edge = p1->getEdge(i);
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					temp_edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
					delete active_edges->removeDataItem(aact->getIndex());
				}
			}
		}
		// switch all blocks incident to point
		for(int i = 0; i < switched_blocks.countInt(); i++){
			switched_blocks[i]->switchPointsWithFaces(p1, p0);
		}
		assert(p1->getRank() == 0);
	}else{ // --> move using Laplace
		if(active_edges){
			// remove tag from edges which will be removed
			for(int i = 0; i < p1->getRank(); i++){
				MeshEdge3d* temp_edge = p1->getEdge(i);
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					temp_edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
					delete active_edges->removeDataItem(aact->getIndex());
				}
			}
		}
		// switch all blocks incident to point
		for(int i = 0; i < switched_blocks.countInt(); i++){
			switched_blocks[i]->switchPointsWithFaces(p1, p0);
		}
		assert(p1->getRank() == 0);
		old_pt = p0->getCoordinates();
		MeshGenerator3dQuality::movePointByLaplaceForVariableMetric(mc, p0);
	}

	// Check if the mesh is OK
	valid = true;
	// ... check volume for switched blocks
	for(int i = 0; valid && i < switched_blocks.countInt(); i++){
//			valid = !switched_blocks[i]->isInverted();
		valid = (switched_blocks[i]->getVolume(mc) > METRIC_SMALL_NUMBER);
//			double vol = switched_blocks[i]->getVolumeNoMetric();
//			LOG4CPLUS_INFO(MeshLog::logger_mesh, "vol = " << vol);
	}
	if(valid){
		for(int i = 0; i < blocks.countInt(); i++)
			delete mesh->removeMeshBlock(blocks[i]);
		delete mesh->removeMeshPoint(p1);

		if(active_edges){
			// update active edges incident to p0
			for(int i = 0; i < p0->getRank(); i++){
				MeshEdge3d* edge = p0->getEdge(i);
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					mc.countMetricAtPoints(p0, edge->getOtherPoint(p0));
					double len = edge->getLength(mc);
					if(len != aact->len){
						aact->len = len;
						active_edges->updateDataItemPosition(aact);
					}
				}
			}
		}
		//SHOW_MESH_NORESET("After collapse", mesh->getDebugViewSet(p0, nullptr, 4.0));
		return p0;
	}else{
		//SHOW_MESH("Error after collapse?", mesh->getDebugViewSet(p0));
		// ... restore 
		if(!p0->isBorder()) p0->setCoordinates(old_pt);
		for(int i = 0; i < switched_blocks.countInt(); i++){
			switched_blocks[i]->switchPointsWithFaces(p0, p1);
		}
		for(int i = 0; i < blocks.countInt(); i++){
			blocks[i]->attachToFaces();
		}
		//SHOW_MESH_NORESET("Restored after collapse", mesh->getDebugViewSet(p0, nullptr, 4.0));
	}

	return nullptr;
}

//#define SHOW_ACTIVE_EDGES

/// Collapse too short edges (inner)
int MeshGenerator3dAdapt::collapseInnerEdges(Metric3dContext& mc, MeshContainer3d* mesh,
		TagExtended::TagType tag_type, int tag_value)
{
//	START_CLOCK("MeshGenerator3d::collapseInnerEdges");
	static const double MAX_LEN = 0.6;
	// gather edges

#ifdef SHOW_ACTIVE_EDGES
	MeshViewSet* set = new MeshViewSet;
#endif

	DataContainer<MeshEdge3d::ActiveEdge> active_edges(mesh->getPointsCount(), true);
	for(IteratorEdge3d it = mesh->getFirstEdge3d(tag_type, tag_value); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(edge->isBorder()){
			//set->addEdge(edge);
			continue;
		}
		MeshPoint3d* p0 = edge->getMeshPoint(0);
		MeshPoint3d* p1 = edge->getMeshPoint(1);

		if(p0->isBorder() && p1->isBorder()) continue;

		mc.countMetricAtPoints(p0, p1);
		double len = edge->getLength(mc);
		if(len < MAX_LEN){
			MeshEdge3d::ActiveEdge* act = new MeshEdge3d::ActiveEdge(edge, len);
			active_edges.addDataItem(act);
			edge->setPtrTag(TagExtended::TAG_COLLAPSE_3D, act);
#ifdef SHOW_ACTIVE_EDGES
			set->addEdge(edge);
#endif
		}
	}

#ifdef SHOW_ACTIVE_EDGES
	SHOW_MESH("Inner edges to collapse", set);
#endif

//	mesh->setHeapOrder(false); // should be false anyway...

//	assert(mesh->isValid());

	int removed_count = 0;
//	int counter = 0;

	while(active_edges.countInt() > 0){
//		counter++;
		MeshEdge3d::ActiveEdge* act = active_edges.removeDataItem(0);
		act->edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
		MeshEdge3d* act_edge = act->edge;
		double act_len = act->len;
		delete act;

		if(act_len > MAX_LEN || act_edge->isBorder())
			continue;

		if(MeshGenerator3dAdapt::collapseInnerEdge(mc, mesh, act_edge, &active_edges))
			removed_count++;

	}

//	STOP_CLOCK("MeshGenerator3d::collapseInnerEdges");
	return removed_count;
}

/// Collapse border edge, return resulting point or nullptr - if impossible
MeshPoint3d* MeshGenerator3dAdapt::collapseBoundaryEdge(Metric3dContext& mc, MeshContainer3d* mesh,
	MeshEdge3d* act_edge, DataContainer<MeshEdge3d::ActiveEdge> * active_edges)
{
	MeshPoint3d* p0 = act_edge->getMeshPoint(0);
	MeshPoint3d* p1 = act_edge->getMeshPoint(1);

	bool p0_corner = p0->isBorder(TagBorder::CORNER | TagBorder::FIXED);
	bool p1_corner = p1->isBorder(TagBorder::CORNER | TagBorder::FIXED);
	bool p0_ridge = p0->isBorder(TagBorder::RIDGE);
	bool p1_ridge = p1->isBorder(TagBorder::RIDGE);

	if(p0_corner && p1_corner) // one at most
		return nullptr;

	if(p1_corner){ // switch direction of collapse, if necessary
		if(!p0_corner){
			MeshPoint3d* tmp_p = p0; p0 = p1; p1 = tmp_p;
			bool tmp = p0_corner; p0_corner = p1_corner; p1_corner = tmp;
			tmp = p0_ridge; p0_ridge = p1_ridge; p1_ridge = tmp;
		}
	}else if(p1_ridge){
		if(!p0_corner && !p0_ridge){
			MeshPoint3d* tmp_p = p0; p0 = p1; p1 = tmp_p;
			bool tmp = p0_corner; p0_corner = p1_corner; p1_corner = tmp;
			tmp = p0_ridge; p0_ridge = p1_ridge; p1_ridge = tmp;
		}
	}

	char p0_border = p0->getBorderFlags();
	char p1_border = p1->getBorderFlags();

	int fct = act_edge->getFaceCount();
	DataVector<MeshFace*> faces(fct);
	DataVector<MeshPoint3d*> fpoints(fct);
	for(int i = 0; i < fct; i++){
		MeshFace* face = act_edge->getFaceAt(i);
		faces.add(face);
		fpoints.add(face->getOtherPoint(p0, p1));
	}

	DataVector<MeshBlock*> blocks(fct);
	act_edge->adjacentBlocks(blocks, false);
	for(int i = 0; i < blocks.countInt(); i++){
		int one_side_count = 0;
		int bfct = blocks[i]->getFaceCount();
		for(int j = 0; j < bfct; j++){
			MeshFace* face = blocks[i]->getFace(j);
			if(face->incidentToEdge(act_edge)) continue;
			if(face->isBoundedBothSides()) continue;
			assert(face->isBounded());
			one_side_count++;
		}
		if(one_side_count > 1)
			return nullptr;
	}


	//if(true){
	//	ostringstream ostr;
	//	ostr << "Selected bedge to collapse";
	//	ostr << ": p0 [" << p0->getIndex() << "] " << (p0_corner ? "corner" : (p0_ridge ? "ridge" : "normal"));
	//	ostr << ", p1 [" << p1->getIndex() << "] " << (p1_corner ? "corner" : (p1_ridge ? "ridge" : "normal"));
	//	SHOW_MESH(ostr.str(), mesh->getDebugViewSet(p0, p1, 4.0));
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, ostr.str());
	//}

	// store border data
	DataVector<std::shared_ptr<BorderEdge>> border_edges(100);
	DataVector<std::shared_ptr<BorderFace>> border_faces(100);
	// ... for edges
	for(int i = 0; i < p0->getRank(); i++){
		MeshEdge3d* edge = p0->getEdge(i);
		if(edge->isBorder()){
			border_edges.add(
				std::make_shared<BorderEdge>(p0, edge->getOtherPoint(p0), edge->getBorderFlags()));
			edge->clearBorder();
			if(active_edges){
				auto aact = (MeshEdge3d::ActiveEdge*)edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					delete active_edges->removeDataItem(aact->getIndex());
					edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
				}
			}
		}
	}
	for(int i = 0; i < p1->getRank(); i++){
		MeshEdge3d* edge = p1->getEdge(i);
		if(edge->isBorder()){
			border_edges.add(
				std::make_shared<BorderEdge>(p1, edge->getOtherPoint(p1), edge->getBorderFlags()));
			edge->clearBorder();
			if(active_edges){
				auto aact = (MeshEdge3d::ActiveEdge*)edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					delete active_edges->removeDataItem(aact->getIndex());
					edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
				}
			}
		}
	}
	// ... for faces
	DataVector<MeshFace*> adj_faces(200);
	p0->adjacentFaces(adj_faces);
	p1->adjacentFaces(adj_faces);
	for(int i = 0; i < adj_faces.countInt(); i++){
		MeshFace* face = adj_faces[i];
		if(face->isBorder()){
			auto bf = std::make_shared<BorderFace>(
				face->getPoint(0), face->getPoint(1), face->getPoint(2), face->getBorderFlags());
			bf->tags.copyAllTags(face);
			border_faces.add(bf);
			face->clearBorder();
		}
	}

	// ... and the point to be removed
	char border_point = p1->getBorderFlags();
	p1->clearBorder();

	for(int i = 0; i < blocks.countInt(); i++)
		blocks[i]->detachFromFaces();
	faces.clear(); // --> no longer valid!
	DataVector<MeshBlock*> switched_blocks(10);

	p1->adjacentBlocks(switched_blocks, false);

	// check extra conditions
	bool valid = true;

	// additional validity check
	DataHashTable<MeshEdge3d*> check_edges(2*switched_blocks.countInt(), nullptr);
	for(int i = 0; valid && i < switched_blocks.countInt(); i++){
		MeshBlock* sblock = switched_blocks[i];
		for(int j = 0; valid && j < sblock->getEdgeCount(); j++){
			MeshEdge3d* edge = sblock->getEdge(j);
			if( edge->getMeshPoint(0)->getEdgeToPoint(p0) != nullptr &&
				edge->getMeshPoint(1)->getEdgeToPoint(p0) != nullptr &&
				!check_edges.insert(edge))
			{
				valid = false;
			}
		}
	}

	if(!valid){
		// restore changes
		// ... reattach faces
		for(int i = 0; i < blocks.countInt(); i++)
			blocks[i]->attachToFaces();
		// ... restore border flags
		for(int i = 0; i < border_edges.countInt(); i++){
			auto be = border_edges[i];
			MeshEdge3d* edge = be->p0->getEdgeToPoint(be->p1);
			if(edge) edge->setBorder(be->border);
		}
		for(int i = 0; i < border_faces.countInt(); i++){
			auto bf = border_faces[i];
			MeshFace* face = bf->p0->getFaceToPoints(bf->p1, bf->p2);
			if(face){
				assert(!face->isBoundedBothSides());
				face->setBorder(bf->border);
				face->copyAllTags(&(bf->tags));
			}
		}
		p1->setBorder(border_point);

		return nullptr;
	}

	// switch all elements incident to point
	for(int i = 0; valid && (i < switched_blocks.countInt()); i++){
		if(switched_blocks[i]->checkSwitchPointsWithFaces(p1, p0))
			switched_blocks[i]->switchPointsWithFaces(p1, p0);
		else{
			for(int j = 0; j < i; j++) // switch back ...
				switched_blocks[j]->switchPointsWithFaces(p0, p1);
			valid = false;
		}
	}

	if(!valid){
		// restore changes
		// ... reattach faces
		for(int i = 0; i < blocks.countInt(); i++)
			blocks[i]->attachToFaces();
		// ... restore border flags
		for(int i = 0; i < border_edges.countInt(); i++){
			auto be = border_edges[i];
			auto edge = be->p0->getEdgeToPoint(be->p1);
			if(edge) edge->setBorder(be->border);
		}
		for(int i = 0; i < border_faces.countInt(); i++){
			auto bf = border_faces[i];
			auto face = bf->p0->getFaceToPoints(bf->p1, bf->p2);
			if(face){
				assert(!face->isBoundedBothSides());
				face->setBorder(bf->border);
				face->copyAllTags(&(bf->tags));
			}
		}
		p1->setBorder(border_point);

		return nullptr;
	}

	assert(p1->getRank() == 0);

	// Check if the mesh is OK

	// ... check volume for switched blocks
	for(int i = 0; valid && i < switched_blocks.countInt(); i++){
		//valid = !switched_blocks[i]->isInverted();
		valid = (switched_blocks[i]->getVolume(mc) > METRIC_SMALL_NUMBER);
	}

	for(int i = 0; valid && (i < border_faces.countInt()); i++){
		auto bf = border_faces[i];
		MeshPoint3d* bp0 = bf->p0;
		MeshPoint3d* bp1 = bf->p1;
		MeshPoint3d* bp2 = bf->p2;
		if(bp0 == p1){
			bp0 = p0; if(bp1 == p0 || bp2 == p0) continue;
		}else if(bp1 == p1){
			bp1 = p0; if(bp0 == p0 || bp2 == p0) continue;
		}else if(bf->p2 == p1){
			bp2 = p0; if(bp0 == p0 || bp1 == p0) continue;
		}
		MeshFace* face = bp0->getFaceToPoints(bp1, bp2);
		if(face && face->isBoundedBothSides()){
			int vct = 0;
			for(int j = 0; j < 2; j++){
				MeshBlock* block = face->getBlock(j);
				for(int k = 0; k < blocks.countInt(); k++){
					if(block == blocks[k]){
						vct++; break;
					}
				}
			}
			if(vct != 1) valid = false; // can not have "not-in-blocks[]" on both sides
		}
	}

	if(valid){
		for(int i = 0; i < blocks.countInt(); i++)
			delete mesh->removeMeshBlock(blocks[i]);
		delete mesh->removeMeshPoint(p1);

		for(int i = 0; i < border_edges.countInt(); i++){
			auto be = border_edges[i];
			auto edge = p0->getEdgeToPoint(be->p1); // p0 stays p0, b1 becomes p0
			if(edge) edge->setBorderFlags(be->border);
		}
		for(int i = 0; i < border_faces.countInt(); i++){
			auto bf = border_faces[i];
			if(bf->p0 == p1){
				bf->p0 = p0; if(bf->p1 == p0 || bf->p2 == p0) continue;
			}else if(bf->p1 == p1){
				bf->p1 = p0; if(bf->p0 == p0 || bf->p2 == p0) continue;
			}else if(bf->p2 == p1){
				bf->p2 = p0; if(bf->p0 == p0 || bf->p1 == p0) continue;
			}
			MeshFace* face = bf->p0->getFaceToPoints(bf->p1, bf->p2);
			if(face){
				assert(!face->isBoundedBothSides());
				face->setBorder(bf->border);
				face->copyAllTags(&(bf->tags));
			}
		}

		if((p0_border == p1_border) && !p0_corner){ // --> smoothen if both points are in the same category and it is not corner
			MeshGenerator3dQuality::moveBoundaryPointByLaplaceForVariableMetric(mc, p0);
		}

		//SHOW_MESH_NORESET("After collapse", mesh->getDebugViewSet(p0));
		return p0;
	}else{ 
		// invalid, revert...
		for(int i = 0; i < switched_blocks.countInt(); i++){
			switched_blocks[i]->switchPointsWithFaces(p0, p1);
		}
		for(int i = 0; i < blocks.countInt(); i++){
			blocks[i]->attachToFaces();
		}
		// restore border flags
		for(int i = 0; i < border_edges.countInt(); i++){
			auto be = border_edges[i];
			auto edge = be->p0->getEdgeToPoint(be->p1);
			if(edge) edge->setBorder(be->border);
		}
		for(int i = 0; i < border_faces.countInt(); i++){
			auto bf = border_faces[i];
			auto face = bf->p0->getFaceToPoints(bf->p1, bf->p2);
			if(face){
				assert(!face->isBoundedBothSides());
				face->setBorder(bf->border);
				face->copyAllTags(&(bf->tags));
			}
		}
		p1->setBorder(border_point);

		//SHOW_MESH("Restored after collapse", mesh->getDebugViewSet(p0));
	}

	return nullptr;
}

/// Collapse too short edges (boundary)
int MeshGenerator3dAdapt::collapseBoundaryEdges(Metric3dContext& mc, MeshContainer3d* mesh,
		TagExtended::TagType tag_type, int tag_value)
{
//	START_CLOCK("MeshGenerator3d::collapseBoundaryEdges");
	static const double MAX_LEN = 0.5;
	// gather edges

#ifdef SHOW_ACTIVE_EDGES
	MeshViewSet* set = new MeshViewSet;
#endif

	DataContainer<MeshEdge3d::ActiveEdge> active_edges(mesh->getPointsCount(), true);
	for(IteratorEdge3d it = mesh->getFirstEdge3d(tag_type, tag_value); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(!edge->isBorder()) continue;

		MeshPoint3d* p0 = edge->getMeshPoint(0);
		MeshPoint3d* p1 = edge->getMeshPoint(1);

		mc.countMetricAtPoints(p0, p1);
		double len = edge->getLength(mc);
		if(len < MAX_LEN){
			MeshEdge3d::ActiveEdge* act = new MeshEdge3d::ActiveEdge(edge, len);
			active_edges.addDataItem(act);
			edge->setPtrTag(TagExtended::TAG_COLLAPSE_3D, act);
#ifdef SHOW_ACTIVE_EDGES
			set->addEdge(edge);
#endif
		}
	}

#ifdef SHOW_ACTIVE_EDGES
	if(active_edges.countInt() > 0)
		SHOW_MESH("Boundary edges to collapse", set);
	else delete set;
#endif

//	assert(mesh->isValid());

	int removed_count = 0;
//	int counter = 0;

	while(active_edges.countInt() > 0){
//		counter++;
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "counter == " << counter);

		MeshEdge3d::ActiveEdge* act = active_edges.removeDataItem(0);
		act->edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
		MeshEdge3d* act_edge = act->edge;
		double act_len = act->len;
		delete act;

		if(act_len > MAX_LEN || !(act_edge->isBorder()))
			continue;

		if(MeshGenerator3dAdapt::collapseBoundaryEdge(mc, mesh, act_edge, &active_edges))
			removed_count++;
	}
//	STOP_CLOCK("MeshGenerator3dSurface::collapseBoundaryEdges");
	return removed_count;
}

/// Insert nodes at boundary, for too long edges, return number of inserted points
int MeshGenerator3dAdapt::addBoundaryNodesSplit(Metric3dContext& mc, MeshContainer3d* mesh,
	TagExtended::TagType tag_type, int tag_value)
{
	static const double MIN_LEN = 1.8;

	START_CLOCK("MeshGenerator3dAdapt::addBoundaryNodesSplit");
	// gather edges

#ifdef SHOW_ACTIVE_EDGES
	MeshViewSet* set = new MeshViewSet;
#endif

	DataContainer<MeshEdge3d::ActiveEdge> active_edges(mesh->getPointsCount(), true);
	for(IteratorEdge3d it = mesh->getFirstEdge3d(tag_type, tag_value); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(!edge->isBorder()) continue;

		MeshPoint3d* p0 = edge->getMeshPoint(0);
		MeshPoint3d* p1 = edge->getMeshPoint(1);

		mc.countMetricAtPoints(p0, p1);
		double len = edge->getLength(mc);
		if(len > MIN_LEN){
			MeshEdge3d::ActiveEdge* act = new MeshEdge3d::ActiveEdge(edge, 1.0/len);
			active_edges.addDataItem(act);
//			edge->setPtrTag(TagExtended::TAG_SPLIT_EDGE_3D, act);
#ifdef SHOW_ACTIVE_EDGES
			set->addEdge(edge);
#endif
		}
	}

#ifdef SHOW_ACTIVE_EDGES
	if(active_edges.countInt() > 0)
		SHOW_MESH("Boundary edges to split", set);
	else delete set;
#endif


	int inserted_count = 0;

	while(active_edges.countInt() > 0){
		MeshEdge3d::ActiveEdge* act = active_edges.removeDataItem(0);
//		act->edge->removeTag(TagExtended::TAG_SPLIT_EDGE_3D);
		MeshEdge3d* act_edge = act->edge;
		if(! act_edge->isBorder()){
			delete act;
			continue;
		}

		MeshPoint3d* p0 = act_edge->getMeshPoint(0);
		MeshPoint3d* p1 = act_edge->getMeshPoint(1);
		inserted_count++;

		bool show_case = false;
		if(show_case){
			ostringstream ostr;
			ostr << "Selected bedge to split, len=" << (1.0/act->len);
			ostr << endl << "#" << inserted_count;
			ostr << ": p0 [" << p0->getIndex() << "], p1 [" << p1->getIndex() << "] ";
			ostr << (act_edge->isBorder(TagBorder::RIDGE) ? "ridge" : "normal");
			SHOW_MESH(ostr.str(), mesh->getDebugViewSet(p0, p1, 4.0));
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, ostr.str());
		}

		DataVector<std::shared_ptr<BorderFace>> border_faces(100);
		char border_edge = act_edge->getBorderFlags();
		act_edge->clearBorder();

		int fct = act_edge->getFaceCount();
		for(int i = 0; i < fct; i++){
			MeshFace* face = act_edge->getFaceAt(i);
			if(face->isBorder()){
				auto bf = std::make_shared<BorderFace>(p0, p1, face->getOtherPoint(p0, p1), face->getBorderFlags());
				bf->tags.copyAllTags(face);
				border_faces.add(bf);
				face->clearBorder();
			}
		}

		DataVector<MeshBlock*> blocks(fct);
		act->edge->adjacentBlocks(blocks, false);
		delete act;

		MeshPoint3d* point = new MeshPoint3d(act_edge->getPoint(0.5));
		mesh->addMeshPoint(point);
		point->setBorder(border_edge); // ie. ridge or normal boundary
		if(tag_type != TagExtended::TAG_NONE)
			point->setIntTag(tag_type, tag_value);

		// move and duplicate tetrahedra
		for(int i = 0; i < blocks.countInt(); i++){
			MeshBlock* block = blocks[i];
			block->switchPointsWithFaces(p1, point);
			MeshPoint3d* pts[4] = {
				block->getPoint(0), block->getPoint(1), block->getPoint(2), block->getPoint(3) 
			};
			for(int j = 0; j < 4; j++) if(pts[j] == point) { pts[j] = p1; break; }
			for(int j = 0; j < 4; j++) if(pts[j] == p0) { pts[j] = point; break; }
			MeshBlock* block_x = new MeshTetrahedron(pts[0], pts[1], pts[2], pts[3]);
			block_x->setAreaID(block->getAreaID());
			mesh->addMeshBlock(block_x);
		}

		// restore border-info
		for(int i = 0; i < border_faces.countInt(); i++){
			auto bf = border_faces[i];
			auto edge = point->getEdgeToPoint(bf->p2); assert(edge);
			edge->setBorder(bf->border);
			auto f0 = edge->getFaceToPoint(bf->p0); assert(f0);
			auto f1 = edge->getFaceToPoint(bf->p1); assert(f1);
			f0->setBorder(bf->border);
			f0->copyAllTags(&(bf->tags));
			f1->setBorder(bf->border);
			f1->copyAllTags(&(bf->tags));
		}
		MeshEdge3d* edge_p0 = point->getEdgeToPoint(p0); assert(edge_p0);
		edge_p0->setBorder(border_edge);
		MeshEdge3d* edge_p1 = point->getEdgeToPoint(p1); assert(edge_p0);
		edge_p1->setBorder(border_edge);

		// project new point onto associated surface

		// swap tetrahedra in the vicinity + Laplace for the new point
		DataVector<MeshEdge3d*> test_edges(point->getRank());
		for(int i = 0; i < point->getRank(); i++){
			MeshEdge3d* edge = point->getEdge(i);
			if(!edge->isBorder()) continue;
			MeshEdge3d* edge_sw = MeshGenerator3d::swap22(mc, mesh, edge);
			if(edge_sw){ // if swap successful
				i--;
				test_edges.add(edge_sw);
			}else{
				test_edges.add(edge);
			}
		}
		MeshGenerator3dQuality::moveBoundaryPointByLaplaceForVariableMetric(mc, point);

		// add new edges to active_edges
		for(int i = 0; i < test_edges.countInt(); i++){
			MeshEdge3d* edge = test_edges[i];
			assert(edge->isBorder());

			MeshPoint3d* p0 = edge->getMeshPoint(0);
			MeshPoint3d* p1 = edge->getMeshPoint(1);

			if(tag_type != TagExtended::TAG_NONE){
				if( !p0->checkIntTag(tag_type, tag_value) &&
					!p1->checkIntTag(tag_type, tag_value)) continue;
			}

			mc.countMetricAtPoints(p0, p1);
			double len = edge->getLength(mc);
			if(len > MIN_LEN){
				MeshEdge3d::ActiveEdge* act = new MeshEdge3d::ActiveEdge(edge, 1.0/len);
				active_edges.addDataItem(act);
	//			edge->setPtrTag(TagExtended::TAG_SPLIT_EDGE_3D, act);
			}
		}

		if(show_case){
			SHOW_MESH_NORESET("After split", mesh->getDebugViewSet(point));
		}

		// ... or make outer loop with final smoothing for selected set of nodes,
		//		then repeating gathering of active edges again (for selected set of nodes)
	}

	STOP_CLOCK("MeshGenerator3dAdapt::addBoundaryNodesSplit");

	return inserted_count;
}

/// returns total number of blocks in sub-regions with blocks deviating from the given control space
///   TODO -> should take border faces into account (i.e. make separate regions in such case...)
int MeshGenerator3dAdapt::findAdaptRegions(Metric3dContext& mc, MeshContainer3d* mesh,
	DataSimpleList< DataVector<MeshBlock*> > & regions,
	double minq, int min_region_size, int layers)
{
	START_CLOCK("MeshGenerator3dAdapt::findAdaptRegions");
	// Start with metric-quality check
//	double minq_inv = 1.0 / minq;
	int bct = mesh->getBlocksCount();
	DataVector<int> qblocks(bct, 0);
	for(int i = 0; i < bct; i++){
		MeshTetrahedron* tetra = (MeshTetrahedron*)mesh->getBlockAt(i);
		if(tetra->getType() != BLOCK_TETRA) continue;
		double q = tetra->countMetricDiffQuality(mc);		
		if(q < minq) qblocks[i] = -1;
		//int ect = tetra->getEdgeCount();
		//double edge_len = 0.0;
		//for(int j = 0; j < ect; j++)
		//	edge_len += tetra->getEdge(j)->getLength(mc, true);
		//edge_len /= ect;
		//if(edge_len < minq || edge_len > minq_inv) qblocks[i] = -1;
	}
	// Check size of regions
	int region_counter = 0;
	for(int i = 0; i < bct; i++){
		if(qblocks[i] != -1) continue;
		// start new region
		qblocks[i] = ++region_counter;
		// find other members of this region by adjacency check
		DataVector<int> region_blocks(bct);
		region_blocks.add(i);
		int k = 0;
		while(k < region_blocks.countInt()){
			MeshBlock* block = mesh->getBlockAt(region_blocks[k++]);
			int bfct = block->getFaceCount();
			for(int j = 0; j < bfct; j++){
				MeshBlock* other_block = block->getNeighbour(j);
				if(!other_block) continue;
				int ind = other_block->getIndex();
				if(qblocks[ind] == -1){
					qblocks[ind] = region_counter;
					region_blocks.add(ind);
				}
			}
		}
		// check count
		if(region_blocks.countInt() < min_region_size){ // cancel region
			for(int j = 0; j < region_blocks.countInt(); j++){
				qblocks[region_blocks[j]] = 0; // clear
			}
			region_counter--;
		}
	}
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "=== Total regions before layers: " << region_counter);
	// add layers (and possibly combine regions if necessary)
	for(int i = 0; i < layers; i++){
		// mark adjacent blocks
		DataHashTable<MeshPoint3d*> visited_points(2*mesh->getPointsCount(), nullptr);
		for(int j = 0; j < bct; j++){
			int region_id = qblocks[j];
			if(region_id > 0){
				// check all blocks adjacent to vertices of this block
				MeshBlock* block = mesh->getBlockAt(j);
				int bpct = block->getPointCount();
				for(int k = 0; k < bpct; k++){
					MeshPoint3d* point = block->getPoint(k);
					if(!visited_points.insert(point)) continue;
					DataVector<MeshBlock*> pblocks;
					if(point->adjacentBlocks(pblocks)){
						for(int m = 0; m < pblocks.countInt(); m++){
							int ind = pblocks[m]->getIndex();
							if(qblocks[ind] == 0) qblocks[ind] = -region_id;
						}
					}
				}
			}
		}
		for(int j = 0; j < bct; j++){
			int region_id = qblocks[j];
			if(region_id < 0) qblocks[j] = -region_id;
		}
		// combine regions if necessary
		for(int j = 0; j < bct; j++){
			int region_id = qblocks[j];
			if(region_id == 0) continue;
			MeshBlock* block = mesh->getBlockAt(j);
			int bfct = block->getFaceCount();
			for(int k = 0; k < bfct; k++){
				MeshBlock* other_block = block->getNeighbour(k);
				if(!other_block) continue;
				int ind = other_block->getIndex();
				if(ind < j || qblocks[ind] == 0 || qblocks[ind] == region_id) continue;
				// contact with other region -> join
				int other_region_id = qblocks[ind];
				for(int m = 0; m < bct; m++)
					if(qblocks[m] == other_region_id) qblocks[m] = region_id;
			}
		}
	}
	// Prepare result
	int block_count = 0;
	for(int i = 0; i < bct; i++){
		int region_id = qblocks[i];
		if(region_id > 0){
			DataVector<MeshBlock*> rblocks(bct-i);
			rblocks.add(mesh->getBlockAt(i));
			qblocks[i] = -region_id;
			for(int j = i+1; j < bct; j++)
				if(qblocks[j] == region_id){
					rblocks.add(mesh->getBlockAt(j));
					qblocks[j] = -region_id;
				}
			block_count += rblocks.countInt();
			regions.append(rblocks);
		}
	}
	STOP_CLOCK("MeshGenerator3dAdapt::findAdaptRegions");

	return block_count;
}

// Search for local reparameterization surfaces (for whole mesh, or only points with the given tag)
int MeshGenerator3dAdapt::identifyLocalSurfaces(MeshContainer3d* mesh, double tolerance, 
		TagExtended::TagType tag_type, int tag_value)
{
	int count = 0;

	// Browse border points
	int pct = mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if(!point->isBorder()) continue;
		if(tag_type != TagExtended::TAG_NONE && !point->hasAnyIntFlags(tag_type, tag_value))
			continue;
		int pect = point->getRank();
		for(int j = 0; j < pect; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(!edge->isBorder()) continue;
			int efct = edge->getFaceCount();
			for(int k = 0; k < efct; k++){
				MeshFace* face = edge->getFaceAt(k);
				if(face->isBorder() && !face->hasLocalSurface()){
					count += MeshGenerator3dAdapt::approximateLocalSurface(mesh, face, tolerance, tag_type, tag_value);
				}
			}
		}
	}

	return count;
}

#define SHOW_LOCAL_SURFACES

// Search for local reparameterization surface for the given (starting) face 
//   and the mesh in the vicinity (without crossing ridge boundaries!)
int MeshGenerator3dAdapt::approximateLocalSurface(MeshContainer3d* mesh,
		MeshFace* face_init, double tolerance, TagExtended::TagType tag_type, int tag_value)
{
	assert(mesh); if(!mesh) return 0;
	assert(face_init); if(!face_init) return 0;

	// 1. Gather local neigborhood on boundary surface (without crossing boundary ridge)

	int pct = mesh->getPointsCount();

	DataVector<MeshFace*> faces_list(pct);
	DataVector<int> faces_src(pct);
	DataVector<int> faces_layer(pct);
	DataVector<MeshPoint3d*> faces_other_pts(pct);
	DataHashTable<MeshEdge3d*> edges_visited(pct, nullptr);
	DataHashTable<MeshFace*> faces_visited(pct, nullptr); // if visited
	DataVector<int> points_layer(pct, -1);

	DataVector< std::shared_ptr<DataVector<int>> > points_list(50);
	DataVector<double> threshold_layers(50);

	int face_0 = 0;
	int face_1 = 1;
	faces_list.add(face_init);
	faces_src.add(0);
	faces_layer.add(0);
	faces_other_pts.add(nullptr);
	faces_visited.insert(face_init);
	auto ids = std::make_shared<DataVector<int>>(3);
	points_list.add(ids);
	for(int i = 0; i < face_init->getPointCount(); i++){
		int id = face_init->getPoint(i)->getIndex();
		points_layer[id] = 0; // zero layer
		ids->add(id);
	}
	int local_points_count = 3;

	double min_dist2 = std::min( std::min(
		face_init->getPoint(0)->getCoordinates().distance2(face_init->getPoint(1)->getCoordinates()),
		face_init->getPoint(0)->getCoordinates().distance2(face_init->getPoint(2)->getCoordinates())),
		face_init->getPoint(1)->getCoordinates().distance2(face_init->getPoint(2)->getCoordinates()));

	while(face_1 > face_0){
		for(int fi = face_0; fi < face_1; fi++){
			// check face fi
			MeshFace* face = faces_list[fi];
			int fect = face->getEdgeCount();
			for(int i = 0; i < fect; i++){
				MeshEdge3d* edge = face->getEdge(i);
				if(edge->isBorder( TagBorder::RIDGE )) continue; // without crossing sharp edges
				if(!edges_visited.insert( edge )) continue; // already visited
				MeshPoint3d* ep0 = edge->getMeshPoint(0);
				MeshPoint3d* ep1 = edge->getMeshPoint(1);
				if(tag_type != TagExtended::TAG_NONE && // at least one point of the edge should be appropriately tagged
						!ep0->hasAnyIntFlags(tag_type, tag_value) && 
						!ep1->hasAnyIntFlags(tag_type, tag_value)) 
					continue;

				int efct = edge->getFaceCount();
				MeshFace* other_face = nullptr;
				for(int j = 0; j < efct; j++){
					MeshFace* f = edge->getFaceAt(j);
					if(f != face && f->isBorder()){
						if(other_face){ // more than 2 border faces for edge -> it should be marked as ridge!
							edge->setBorderFlags( TagBorder::RIDGE );
							other_face = nullptr;
							break;
						}else{
							other_face = f;
						}
					}
				}
				if(!other_face || faces_visited.contains(other_face)) continue;
				// ok
				faces_list.add(other_face);
				faces_src.add(fi); // since face -> (through edge) -> other_face
				// local/ext ?
				MeshPoint3d* other_point = other_face->getOtherPoint(ep0, ep1);
				int id = other_point->getIndex();
				int layer = points_layer[id];
				if(layer == -1){
					points_layer[id] = layer = 1+std::min( points_layer[ep0->getIndex()], points_layer[ep1->getIndex()] );
					local_points_count++;
					while(layer >= points_list.countInt()) 
						points_list.add(std::make_shared<DataVector<int>>(100));
					points_list[layer]->add(id);
				}
				
				faces_visited.insert(other_face);
				faces_layer.add(layer);
				faces_other_pts.add(other_point);

				min_dist2 = std::min(min_dist2, 
					ep0->getCoordinates().distance2( other_point->getCoordinates()));
				min_dist2 = std::min(min_dist2, 
					ep1->getCoordinates().distance2( other_point->getCoordinates()));
			}
		}
		face_0 = face_1;
		face_1 = faces_list.countInt();
	}

	// Gather all identified points into a single list, sorting by layer number
	DataVector<DPoint3d> point_matrix( local_points_count );
	int layer_count = points_list.countInt();
	DataVector<int> points_count( layer_count);
	int total_count = 0;
	for(int i = 0; i < layer_count; i++){
		auto ids = points_list[i];
		int ct = ids->countInt();
		total_count += ct;
		points_count.add(total_count);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Layer #" << i << ", total_count = " << total_count);
		for(int j = 0; j < ct; j++)
			point_matrix.add( mesh->getPointAt( ids->get(j) )->getCoordinates() );
	}
	assert( points_count[ layer_count-1 ] == local_points_count );

	double eps2 = tolerance * tolerance * min_dist2;
	// Try to fit surface, starting from the most extended number of layers until the zero-layer (single face)
	SurfacePtr surface;
	int layer_i = layer_count;
	while(!surface && (--layer_i >= 0)){
		point_matrix.leaveOnly(points_count[layer_i]);
		surface = MeshGenerator3dAdapt::fitLocalSurface(point_matrix, eps2);
	}

	assert(surface); if(!surface) return 0;

	mesh->addLocalSurface(surface);

	int fct = faces_list.countInt();
	if(layer_i < (layer_count-1)){
#ifdef SHOW_LOCAL_SURFACES
		MeshViewSet* set = new MeshViewSet;
#endif
		// extend faces for this local surface (by neighbors)
		for(int i = 0; i < fct; i++){
			if(faces_layer[i] <= layer_i){
#ifdef SHOW_LOCAL_SURFACES
				set->addFaceWithEdges(faces_list[i], 0);
#endif
				continue;
			}
			int source_id = faces_src[i];
			if(faces_layer[source_id] > layer_i) continue;
			MeshPoint3d* other_point = faces_other_pts[i];
			int p_layer = points_layer[other_point->getIndex()];
			if(p_layer > layer_i){
				// check distance
				const DPoint3d& other_pt = other_point->getCoordinates();
				const DPoint3d surf_pt = surface->getPoint(surface->getParameters(other_pt));
				double dist2 = surf_pt.distance2(other_pt);
				if(dist2 < eps2) // distance ok, include (modify layer number)
					points_layer[other_point->getIndex()] = p_layer = layer_i;
				else continue;
			}
			faces_layer[i] = layer_i;
#ifdef SHOW_LOCAL_SURFACES
			set->addFaceWithEdges(faces_list[i], 1);
#endif
		}
#ifdef SHOW_LOCAL_SURFACES
		SHOW_MESH("Set of faces for a local surface", set);
#endif
	}

	// set surface for faces
	int sface_count = 0;
	for(int i = 0; i < fct; i++){
		if(faces_layer[i] > layer_i) continue;
		MeshFace* face = faces_list[i];
		// ... else, set surface 
		SurfaceConstPtr face_surface = face->getCurrentLocalSurface();
		if(face_surface) continue; // already set with different surface, let it be
		face->setLocalSurface(surface);
		sface_count++;
	}

	// TODO? Check if this surface include any/some of other, less vast surfaces

	// set surface for points
	//  TODO_CURVES ...
	LOG4CPLUS_WARN(MeshLog::logger_console, "TODO_CURVES...");
	for(int i = 0; i < layer_count; i++){
		auto ids = points_list[i];
		int ct = ids->countInt();
		for(int j = 0; j < ct; j++){
			int id = ids->get(j);
			if((i > layer_i) && (points_layer[id] > layer_i)) continue;
			MeshPoint3d* point = mesh->getPointAt(id);
			SurfaceConstPtr psurface = point->getLocalValidSurface();
			if(!psurface) point->setLocalSurface(surface, 
				surface->getParameters( point->getCoordinates() ), 1.0 ); // TODO param should be reused, from earlier...
		}
	}

	return sface_count;
}

/// Fit surface to cloud of points
SurfacePtr MeshGenerator3dAdapt::fitLocalSurface( const DataVector<DPoint3d> & points, double eps2 )
{
	// a) plane
	DPlane plane;
	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(points, plane);
	if( (plane_max_dist*plane_max_dist) <= eps2){
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
			"Local surface (plane), max_dist ratio = " << (plane_max_dist / sqrt(eps2)));
#ifdef SHOW_LOCAL_SURFACES
		if(true){
			MeshViewSet* set = new MeshViewSet;
			for(int i = 0; i < points.countInt(); i++){
				set->addPoint(points[i]);
			}
			set->addPoint(plane.p0, 2, 0);
			set->addEdge(plane.p0, plane.p0 + plane.e0, 1);
			set->addEdge(plane.p0, plane.p0 + plane.e1, 1);
			SHOW_MESH("plane fit", set);
		}
#endif
		return std::make_shared<SurfacePlane>(plane);
	}

	DPlanarQuadric pquadric;
	double pquadric_max_dist = DLeastSquaresFitting::fitPlanarQuadric(points, plane, pquadric);
	if( (pquadric_max_dist*pquadric_max_dist) <= eps2){
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
			"Local surface (planar quadric), max_dist ratio = " << (pquadric_max_dist / sqrt(eps2)));
#ifdef SHOW_LOCAL_SURFACES
		if(true){
			MeshViewSet* set = new MeshViewSet(points.countInt(), 2*DPlanarQuadric::SKETCH_LINES, 1, 0);
			pquadric.createViewSetForPoints(set, points);
			SHOW_MESH("planar quadric fit", set);
		}
#endif
		return std::make_shared<SurfacePlanarQuadric>(pquadric);
	}


	return nullptr;
}
