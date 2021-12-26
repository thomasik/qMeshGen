// MeshSplit3d.cpp: implementation of the MeshSplit3d class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshData.h"
#include "MeshSplit3d.h"
#include "MeshContainer3d.h"
#include "MeshEdge3d.h"
#include "MeshFace.h"
#include "MeshTriangle3d.h"
#include "MeshViewSet.h"

/*
/// **** not correct yet, and 3 times longer... ***
/// Partition mesh (remove marked blocks and transfer them to other mesh)
MeshContainer3d* MeshSplit3d::splitMeshByBlocks(MeshContainer3d* mesh, TagExtended::TagType tag_type, int tag_value)
{
	START_CLOCK("MeshSplit3d::splitMeshByBlocks");

	int block_count = mesh->getBlocksCount();
	MeshContainer3d* other_mesh = new MeshContainer3d(block_count/2);
	bool heap_order = mesh->isHeapOrder();
	if(heap_order) mesh->setHeapOrder(false);

	//  1. browse blocks, mark adjacent faces/edges/points with 
	//		TAG_CUT_3D as |1 or |2 depending on split destination
	for(int i = 0; i < block_count; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		// 2 - to the cut-mesh
		// 1 - stays in the old mesh
		int cut_tag = (block->getIntTag(tag_type, tag_value-1) == tag_value) ? 2 : 1;
		// --> flag faces
		int fbct = block->getFaceCount();
		for(int j = 0; j < fbct; j++)
			block->getFace(j)->setIntFlag(TagExtended::TAG_CUT_3D, cut_tag);
		// --> flag edges
		int ebct = block->getEdgeCount();
		for(int j = 0; j < ebct; j++)
			block->getEdge(j)->setIntFlag(TagExtended::TAG_CUT_3D, cut_tag);
		// --> flag vertices
		int pbct = block->getPointCount();
		for(int j = 0; j < pbct; j++)
			block->getPoint(j)->setIntFlag(TagExtended::TAG_CUT_3D, cut_tag);
	}

	//  2. duplicate faces/edges/points with TAG_CUT_3D as |3
	//  2a. --> duplicate and move points
	int pct = mesh->getPointsCount();
	DataHashTableKeyValue<MeshPoint3d*,MeshPoint3d*> hpoints(2*pct, nullptr);
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		int cut_tag = point->getIntTag(TagExtended::TAG_CUT_3D, 0);
		// cut_tag == 0 or == 1 -> stays in the old mesh
		if(cut_tag == 2){
			// move
			other_mesh->addMeshPoint(mesh->removeMeshPoint(point));
			pct--; // since the point has been removed from the old mesh ...
			i--;
		}else{ // i.e. cut_tag == 3
			// duplicate
			MeshPoint3d* mpoint = new MeshPoint3d(*point);
			hpoints.insert(point, mpoint);
			other_mesh->addMeshPoint(mpoint);
			// info for merging
			mpoint->setPtrTag(TagExtended::TAG_CUT_MP_3D, point);
			point->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, point->isBorder() ? 1 : 0);
			// mark as border now
			point->setBorder();
			mpoint->setBorder();
		}
	}
	// 2b. --> duplicate edges
	for(IteratorEdge3d it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		int cut_tag = edge->getIntTag(TagExtended::TAG_CUT_3D, 0);
		if(cut_tag == 3){
			// create new edge (clone)
			MeshPoint3d* point0 = edge->getMeshPoint(0);
			MeshPoint3d* point1 = edge->getMeshPoint(1);
			MeshEdge3d* medge = new MeshEdge3d(
				hpoints.getValue(point0, nullptr), 
				hpoints.getValue(point1, nullptr));
			edge->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, edge->isBorder() ? 1 : 0);
			// mark as border now
			edge->setBorder();
			medge->setBorder();
		}
	}
	// 2c. --> duplicate faces
	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
		MeshFace* face = it.getFace();
		int cut_tag = face->getIntTag(TagExtended::TAG_CUT_3D, 0);
		if(cut_tag == 3){
			// create new face (clone)
			MeshPoint3d* point0 = face->getPoint(0);
			MeshPoint3d* point1 = face->getPoint(1);
			MeshPoint3d* point2 = face->getPoint(2);
			assert(face->getType() == FACE_TRIANGLE);
			MeshTriangle3d* mface = new MeshTriangle3d(
				hpoints.getValue(point0, nullptr), 
				hpoints.getValue(point1, nullptr), 
				hpoints.getValue(point2, nullptr));
			face->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, face->isBorder() ? 1 : 0);
			// mark as border now
			face->setBorder();
			mface->setBorder();
		}
	}
	// 3. move blocks
	for(int i = 0; i < block_count; ){
		MeshBlock* block = mesh->getBlockAt(i);
		if(block->getIntTag(tag_type, tag_value-1) == tag_value){
			int pbct = block->getPointCount();
			for(int j = 0; j < pbct; j++){
				MeshPoint3d* point = block->getPoint(j);
				if(point->getIntTag(TagExtended::TAG_CUT_3D, 0) == 3)
					block->switchPointsWithFaces(point, hpoints.getValue(point, nullptr));
			}
			other_mesh->addMeshBlock(mesh->removeMeshBlock(i));
			--block_count;
		}else ++i;
	}

	other_mesh->setConstrainingPhase(mesh->getConstrainingPhase());
	other_mesh->setControlSpace(mesh->getControlSpace());
	other_mesh->setDiscretizationState(mesh->getDiscretizationState());

	mesh->clearAllTags(TagExtended::TAG_CUT_3D);
	other_mesh->clearAllTags(TagExtended::TAG_CUT_3D);

	if(heap_order){
		mesh->setHeapOrder(true);
		other_mesh->setHeapOrder(true);
	}

	STOP_CLOCK("MeshSplit3d::splitMeshByBlocks");

	return other_mesh;
}
*/

/// Partition mesh (remove marked blocks and transfer them to other mesh)
MeshContainer3d* MeshSplit3d::splitMeshByBlocks(MeshContainer3d* mesh, TagExtended::TagType tag_type, int tag_value)
{
	START_CLOCK("MeshSplit3d::splitMeshByBlocks");

	int block_count = mesh->getBlocksCount();
	MeshContainer3d* other_mesh = new MeshContainer3d(block_count/2);
	bool heap_order = mesh->isHeapOrder();
	if(heap_order) mesh->setHeapOrder(false);

	int pct = mesh->getPointsCount();
	DataHashTable<MeshEdge2d*> cut_edges(3*pct, 0);
	DataHashTable<MeshPoint2d*> near_cut_points(pct, 0);

	DataVector<MeshPoint3d*> point_adjacency(pct);

	// identify split surface
	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
		MeshFace* face = it.getFace();
		MeshBlock* block0 = face->getBlock(0);
		MeshBlock* block1 = face->getBlock(1);
		if(!block0 || !block1) continue;
		bool cut_block0 = (block0->getIntTag(tag_type, tag_value-1) == tag_value);
		bool cut_block1 = (block1->getIntTag(tag_type, tag_value-1) == tag_value);
		if(cut_block0 != cut_block1){
			face->setIntTag(TagExtended::TAG_CUT_3D, 1);
			int fpct = face->getEdgeCount();
			assert(fpct == 3);
			DataVector<MeshPoint3d*> points(fpct);
			DataVector<MeshPoint3d*> mpoints(fpct);
			for(int j = 0; j < fpct; j++){
				MeshPoint3d* point = face->getPoint(j);
				points.add(point);
				if(point->availableTag(TagExtended::TAG_CUT_3D)){
					mpoints.add(point_adjacency[point->getIntTag(TagExtended::TAG_CUT_3D)]);
				}else{ // clone point
					MeshPoint3d* mpoint = new MeshPoint3d(*point);
					mpoints.add(mpoint);
					other_mesh->addMeshPoint(mpoint);
					point->setIntTag(TagExtended::TAG_CUT_3D, (int)point_adjacency.add(mpoint));
					mpoint->setIntTag(TagExtended::TAG_CUT_3D, -1);
					mpoint->setPtrTag(TagExtended::TAG_CUT_MP_3D, point);
					mpoint->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, point->isBorder() ? 1 : 0);
				}
			}
			MeshFace* mface = nullptr;
			if(fpct == 3){
				mface = new MeshTriangle3d(mpoints[0], mpoints[1], mpoints[2]);
				mface->copyAllTags(face);
				mface->setIntTag(TagExtended::TAG_CUT_3D, -1);
				mface->setPtrTag(TagExtended::TAG_CUT_MF_3D, face);
				mface->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, face->getBorderFlags());
			}else LOG4CPLUS_ERROR(MeshLog::logger_console,   "MeshSplit3d: Split face not a triangle!");
			//MeshEdge2d* medge = new MeshEdge2d(mpoints[0], mpoints[1]);
			// TODO copy edge parameters or create some new type of edge for inter-partition connections
			if(face->isBorder())
				mface->copyBorderFlagsFrom(face);
			else{
				face->setBorder();
				mface->setBorder();
			}
			for(int j = 0; j < fpct; j++){
				points[j]->setBorder();
				mpoints[j]->setBorder();
				face->getEdge(j)->setBorder();
				mface->getEdge(j)->setBorder();
			}
		}
	}

	// additional check for isolated interface edges/vertices
	DataVector<int> point_marks(pct, 0);
	for(int i = 0; i < block_count; i++){
		MeshBlock* block = mesh->getBlockAt(i);
		bool cut_block = (block->getIntTag(tag_type, tag_value-1) == tag_value);
		int bpct = block->getPointCount();
		for(int j = 0; j < bpct; j++){
			MeshPoint3d* point = block->getPoint(j);
			int cut_tag = point->getIntTag(TagExtended::TAG_CUT_3D, -1);
			if(cut_tag >= 0) continue; // already marked as part of interface
			point_marks[point->getIndex()] |= (cut_block ? 1 : 2);
		}
	}
	// ... check for points on interface (marked as 1+2=3)
	for(int i = 0; i < pct; i++)
		if(point_marks[i] == 3){ // yes, clone
			//LOG4CPLUS_INFO(MeshLog::logger_console, "Found isolated cut-interface point");
			MeshPoint3d* point = mesh->getPointAt(i);
			MeshPoint3d* mpoint = new MeshPoint3d(*point);
			other_mesh->addMeshPoint(mpoint);
			point->setIntTag(TagExtended::TAG_CUT_3D, (int)point_adjacency.add(mpoint));
			mpoint->setIntTag(TagExtended::TAG_CUT_3D, -1);
			mpoint->setPtrTag(TagExtended::TAG_CUT_MP_3D, point);
			mpoint->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, point->isBorder() ? 1 : 0);
			point->setBorder();
			mpoint->setBorder();
		}

	// move inner points (of split volume)
	DataSimpleList<FaceInfo>   border_faces_near_split;
	DataSimpleList<Edge3dInfo> border_edges_near_split;
	for(IteratorFace it = mesh->getFirstFace(); it.isValid(); it.nextFace()){
		MeshFace* face = it.getFace();
		if(face->availableTag(TagExtended::TAG_CUT_3D)) 
			continue;

		MeshBlock* block0 = face->getBlock(0);
		MeshBlock* block1 = face->getBlock(1);
		//bool is_cut_face = (block0 && (block0->getIntTag(tag_type, tag_value-1) == tag_value))
		//		|| (block1 && (block1->getIntTag(tag_type, tag_value-1) == tag_value));
		//if(!is_cut_face) continue;

		int fpct = face->getEdgeCount();
		assert(fpct == 3);
		int split_count = 0;
		for(int i = 0; i < fpct; i++)
			if(face->getPoint(i)->availableTag(TagExtended::TAG_CUT_3D)) ++split_count;
		if(face->isBorder() && split_count > 0){ // border face incident to split, but not part of it
			border_faces_near_split.insert(
				FaceInfo(face->getPoint(0), face->getPoint(1), face->getPoint(2), 
							face->getBorderFlags()));
			face->clearBorder();
			// plus edges ?
			for(int j = 0; j < fpct; j++){
				MeshEdge3d* edge = face->getEdge(j);
				if(edge->isBorder()){
					int split_fcount = 0;
					int fect = edge->getFaceCount();
					for(int k = 0; k < fect; k++) 
						if(edge->getFaceAt(k)->availableTag(TagExtended::TAG_CUT_3D)) ++split_fcount;
					if(split_fcount == 0){
						border_edges_near_split.insert(Edge3dInfo(edge->getMeshPoint(0), edge->getMeshPoint(1)));
						edge->clearBorder();
					}
				}
			}
		}
		int tag0 = block0 ? block0->getIntTag(tag_type, tag_value-1) : tag_value-1;
		int tag1 = block1 ? block1->getIntTag(tag_type, tag_value-1) : tag_value-1;
		if(tag0 == tag_value || tag1 == tag_value){
			for(int j = 0; j < fpct; j++){
				MeshPoint3d* point = face->getPoint(j);
				if(!point->availableTag(TagExtended::TAG_CUT_3D)){
					other_mesh->addMeshPoint(mesh->removeMeshPoint(point));
					point->setIntTag(TagExtended::TAG_CUT_3D,-2);
				}
			}
		}
	}

	// move blocks
	for(int i = 0; i < block_count; ){
		MeshBlock* block = mesh->getBlockAt(i);
		if(block->getIntTag(tag_type, tag_value-1) == tag_value){
			other_mesh->addMeshBlock(mesh->removeMeshBlock(i));
			--block_count;
		}else ++i;
	}

	// adjust incidency for split contour
	int split_block_count = other_mesh->getBlocksCount();
	for(int i = 0; i < split_block_count; i++){
		MeshBlock* block = other_mesh->getBlockAt(i);
		int bpct = block->getPointCount();
		for(int j = 0; j < bpct; j++){
			MeshPoint3d* point = block->getPoint(j);
			assert(point->availableTag(TagExtended::TAG_CUT_3D));
			int tag = point->getIntTag(TagExtended::TAG_CUT_3D);
			if(tag >= 0){ // i.e. split contour (but original points)
				int fct = block->getFaceCount();
				// save oryginal faces
				DataVector<MeshFace*> oryg_faces(fct);
				DataVector<TagExtended> oryg_tags(fct, TagExtended());
				for(int m = 0; m < fct; m++){
					MeshFace* face = block->getFace(m);
					oryg_faces.add(face);
					oryg_tags[m].copyAllTags(face);
				}
				block->switchPointsWithFaces(point, point_adjacency[tag]);
				// restore info
				for(int m = 0; m < fct; m++){
					MeshFace* face = block->getFace(m);
					face->copyAllTags(&(oryg_tags[m]));
					//face->setIntTag(TagExtended::TAG_CUT_3D, -1);
				}
			}
		}
	}

	// adjust border faces near split surface
	while(border_faces_near_split.notEmpty()){
		FaceInfo info = border_faces_near_split.removeFirst();
		MeshFace* face = info.pt0->getFaceToPoints(info.pt1, info.pt2);
		if(face) face->setBorder(info.btype); // this mesh
		else{ // other_mesh
			if(info.pt0->availableTag(TagExtended::TAG_CUT_3D)){
				int tag = info.pt0->getIntTag(TagExtended::TAG_CUT_3D);
				if(tag >= 0) info.pt0 = point_adjacency[tag];
			}
			if(info.pt1->availableTag(TagExtended::TAG_CUT_3D)){
				int tag = info.pt1->getIntTag(TagExtended::TAG_CUT_3D);
				if(tag >= 0) info.pt1 = point_adjacency[tag];
			}
			if(info.pt2->availableTag(TagExtended::TAG_CUT_3D)){
				int tag = info.pt2->getIntTag(TagExtended::TAG_CUT_3D);
				if(tag >= 0) info.pt2 = point_adjacency[tag];
			}
			face = info.pt0->getFaceToPoints(info.pt1, info.pt2);
			if(face) face->setBorder(info.btype);
			else LOG4CPLUS_ERROR(MeshLog::logger_console,   "MeshSplit3d - missing face?");
		}
	}
	// adjust border edges near split surface
	while(border_edges_near_split.notEmpty()){
		Edge3dInfo info = border_edges_near_split.removeFirst();
		MeshEdge3d* edge = info.pt0->getEdgeToPoint(info.pt1);
		if(edge) edge->setBorder(); // this mesh
		else{ // other mesh
			//MeshPoint3d* old_pt0 = info.pt0;
			//MeshPoint3d* old_pt1 = info.pt1;

			if(info.pt0->availableTag(TagExtended::TAG_CUT_3D)){
				int tag = info.pt0->getIntTag(TagExtended::TAG_CUT_3D);
				if(tag > 0) info.pt0 = point_adjacency[tag];
			}
			if(info.pt1->availableTag(TagExtended::TAG_CUT_3D)){
				int tag = info.pt1->getIntTag(TagExtended::TAG_CUT_3D);
				if(tag > 0) info.pt1 = point_adjacency[tag];
			}
			edge = info.pt0->getEdgeToPoint(info.pt1);
			if(edge) edge->setBorder();
			else{
				//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Edge pt0, id = " << old_pt0->getIndex());
				//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Edge pt1, id = " << old_pt1->getIndex());
				//MeshViewSet* set = new MeshViewSet;
				//set->addPoint(info.pt0->getCoordinates(), 1, (info.pt0 == old_pt0) ? 1 : 11);
				//set->addPoint(info.pt1->getCoordinates(), 1, (info.pt1 == old_pt1) ? 2 : 12);
				//set->addEdge(info.pt0->getCoordinates(), info.pt1->getCoordinates());
				//DataVector<MeshBlock*> blocks10(100);
				//DataVector<MeshBlock*> blocks11(100);
				//DataVector<MeshBlock*> blocks20(100);
				//DataVector<MeshBlock*> blocks21(100);
				//info.pt0->adjacentBlocks(blocks20);
				//info.pt1->adjacentBlocks(blocks21);
				//if(info.pt0 != old_pt0) old_pt0->adjacentBlocks(blocks10);
				//if(info.pt1 != old_pt1) old_pt1->adjacentBlocks(blocks11);

				//for(int i = 0; i < blocks10.countInt(); i++) set->addBlock(blocks10[i], 0);
				//for(int i = 0; i < blocks11.countInt(); i++) set->addBlock(blocks11[i], 0);
				//for(int i = 0; i < blocks20.countInt(); i++){
				//	int tag = (blocks20[i]->getIntTag(tag_type, tag_value-1) == tag_value) ? 2 : 1;
				//	set->addBlock(blocks20[i], tag);
				//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Block 20-" << i << " - points: ";
				//	for(int j = 0; j < blocks20[i]->getPointCount(); j++){
				//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, blocks20[i]->getPoint(j)->getIndex() << " ";
				//		set->addPoint(blocks20[i]->getPoint(j));
				//	}
				//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;
				//}
				//for(int i = 0; i < blocks21.countInt(); i++){
				//	int tag = (blocks21[i]->getIntTag(tag_type, tag_value-1) == tag_value) ? 2 : 1;
				//	set->addBlock(blocks21[i], tag);
				//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Block 21-" << i << " - points: ";
				//	for(int j = 0; j < blocks21[i]->getPointCount(); j++){
				//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, blocks21[i]->getPoint(j)->getIndex() << " ";
				//		set->addPoint(blocks21[i]->getPoint(j));
				//	}
				//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;
				//}

				//SHOW_MESH("missing edge???", set);
				LOG4CPLUS_ERROR(MeshLog::logger_console, "MeshSplit3d - missing edge?");
			}
		}
	}

	other_mesh->setConstrainingPhase(mesh->getConstrainingPhase());
	other_mesh->setControlSpace(mesh->getControlSpace());
	other_mesh->setDiscretizationState(mesh->getDiscretizationState());

	mesh->removeAllTags(TagExtended::TAG_CUT_3D);

	if(heap_order){
		mesh->setHeapOrder(true);
		other_mesh->setHeapOrder(true);
	}

	STOP_CLOCK("MeshSplit3d::splitMeshByBlocks");

	return other_mesh;
}

/// Merge two meshes into one (mesh_add will be erased!)
bool MeshSplit3d::mergeMeshes(MeshContainer3d* mesh_base, MeshContainer3d* &mesh_add)
{
	START_CLOCK("MeshSplit3d::mergeMeshes");
	// process faces
	int counter = 0;
//	DataVector<MeshFace*> border_faces(3*mesh_add->getBlocksCount());
	for(IteratorFace it = mesh_add->getFirstFace(); it.isValid(); it.nextFace()){
		MeshFace* face = it.getFace();
		++counter;
		int face_tag = face->getIntTag(TagExtended::TAG_ADAPT_SURF, -1);
		if(face->isBoundedBothSides()) continue; // cannot be part of cut-interface
		assert(face->getPointCount() == 3);
		// retrieve connection to boundary surface of this cut-mesh
		MeshFace* oryg_face = (MeshFace*)face->getPtrTag(TagExtended::TAG_CUT_MF_3D);
		if(!oryg_face) continue; // definitely not a part of cut-interface
		// ... otherwise
		assert(face->availableTag(TagExtended::TAG_CUT_ORYG_DATA));
		char btype = face->getIntTag(TagExtended::TAG_CUT_ORYG_DATA, -2);
		oryg_face->setBorder( btype );
//		if(btype >= 0) 
//			border_faces.add(oryg_face);
//		else{
//			for(int i = 0; i < face->getEdgeCount(); i++)  // first, remove border mark from all adjacent to inner faces
//				face->getEdge(i)->clearBorder();
//		}
		face->setBorderWithEdges(TagBorder::NONE); // prepare to be removed later...
	}

//	for(int i = 0; i < border_faces.countInt(); i++){ // ... and then, add border mark for all adjacent to (cut) border faces
//		MeshFace* face = border_faces[i];
//		for(int j = 0; j < face->getEdgeCount(); j++)
//			face->getEdge(j)->setBorder();
//	}

	// process blocks
	while(mesh_add->getBlocksCount() > 0){
		MeshBlock* block = mesh_add->removeMeshBlock(0);
		mesh_base->addMeshBlock(block);
		// check vertices
		int pct = block->getPointCount();
		for(int i = 0; i < pct; i++){
			MeshPoint3d* point = block->getPoint(i);
			if(!point->isBorder()) continue;
			MeshPoint3d* oryg_point = (MeshPoint3d*)point->getPtrTag(TagExtended::TAG_CUT_MP_3D);
			if(!oryg_point) continue;
			// ... otherwise
			oryg_point->setBorder(point->getIntTag(TagExtended::TAG_CUT_ORYG_DATA, TagBorder::NONE));
			block->switchPointsWithFaces(point, oryg_point);
			if(point->getRank() == 0)
				delete mesh_add->removeMeshPoint(point);
		}
	}
	// move other non-cut points
	while(mesh_add->getPointsCount() > 0){
		mesh_base->addMeshPoint(mesh_add->removeMeshPoint(0));
	}
	// clean
	delete mesh_add;
	mesh_add = nullptr;

	STOP_CLOCK("MeshSplit3d::mergeMeshes");
	return true;
}
