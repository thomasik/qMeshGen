// MeshGenerator2d.cpp: implementation of the MeshGenerator2d class.
//
//////////////////////////////////////////////////////////////////////

#include <iomanip>
#include "MeshGenerator1d.h"
#include "MeshGenerator2d.h"
#include "MeshGenerator3d.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshDomainVolume.h"
#include "MeshDomainSurface.h"
#include "MeshDomainEdge3d.h"
#include "MeshTriangle2d.h"
#include "QuadTree.h"
#include "ControlSpace2d.h"
#include "ControlSpace2dIdentity.h"
#include "ControlSpace2dMatrixUniform.h"
#include "ControlSpace2dMesh.h"
#include "ControlSpace2dProjected.h"
#include "ControlSpace2dQuadTree.h"
#include "ControlSpace2dKdTreeL.h"
#include "ControlSpace2dKdTreeV.h"
#include "SurfaceParametric.h"
#include "MeshViewSet.h"
#include "DEquation.h"
#include "MeshEdge3d.h"
#include "DataVector.h"
#include "DataList.h"
#include "IteratorEdge2d.h"
#include "IteratorEdge3d.h"
#include "DTriangle.h"
#include "OctTree.h"
#include "ControlSpace3dOctree.h"
#include "ControlSpace3dKdTreeL.h"
#include "ControlSpace3dKdTreeV.h"
#include "DLeastSquaresFitting.h"
#include "DPlane.h"
#include "SurfacePlane.h"
#include "MeshArea.h"
#include "Metric2dContext.h"

#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#pragma warning( disable : 4005 )	// disable warning: macro redefinition ...
#define START_CLOCK(text)
#define STOP_CLOCK(text)

//#define PRINT_METRIC_CALL_COUNTERS
//#define STAT_SWAP_INNER_ANGLES

int MeshGenerator2d::param_triangulate_with_inner_nodes = 1;
int MeshGenerator2d::param_triangulation_type = -123; // MeshData::TRIANGULATE_CIRCLE; <- set elsewhere
int MeshGenerator2d::param_quality_improvement = 5432; // MeshData::IMPROVE_OUTER_CIRCLE; <- set elsewhere
double MeshGenerator2d::param_quality_threshold = 0.58;
double MeshGenerator2d::param_initial_mesh_size_ratio = 0.1;
int MeshGenerator2d::param_swap_criterion = MeshData::SWAP_ANGLES;
int MeshGenerator2d::param_swap_maximum = 100000;
int MeshGenerator2d::param_max_auto_retriangulation_count = 10;
bool MeshGenerator2d::show_prediction = true;
int MeshGenerator2d::param_mesh_decomposition = 0;
double MeshGenerator2d::param_mesh_decomposition_width = 1.0;
double MeshGenerator2d::param_mesh_decomposition_gradation = 8.0;

#ifdef PRINT_METRIC_CALL_COUNTERS
	#define INC_COUNT_PT_OUTSIDE_BOX count_pt_outside_box++
	#define INC_COUNT_PT_TOO_NEAR_BOUNDARY count_pt_too_near_boundary++
	#define INC_COUNT_PT_NO_TETRA_FOUND count_pt_no_tetra_found++
	#define INC_COUNT_PT_TOO_NEAR_OTHER_PT count_pt_too_near_other_pt++
	#define INC_COUNT_PT_IN_OTHER_TETRA count_pt_in_other_tetra++
	#define INC_COUNT_PT_TOO_DIFF_METRICS count_pt_too_diff_metrics++
	unsigned int stat_retriangulation_angles_check;
	unsigned int stat_retriangulation_angles_swap;
	unsigned int stat_retriangulation_circle_check;
	unsigned int stat_retriangulation_circle_remove;
	unsigned int stat_retriangulation_new_elements;
#else
	#define INC_COUNT_PT_OUTSIDE_BOX
	#define INC_COUNT_PT_TOO_NEAR_BOUNDARY
	#define INC_COUNT_PT_NO_TETRA_FOUND
	#define INC_COUNT_PT_TOO_NEAR_OTHER_PT
	#define INC_COUNT_PT_IN_OTHER_TETRA
	#define INC_COUNT_PT_TOO_DIFF_METRICS
#endif //PRINT_METRIC_CALL_COUNTERS

#ifdef STAT_SWAP_INNER_ANGLES
	unsigned int stat_swap_inner_angles_counter = 0;
#endif

int MeshGenerator2d::prepareBoundaryMesh(MeshContainer3d *boundary)
{
	int block_count = boundary->getBlocksCount();
#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::prepareBoundaryMesh - clear previous.");
#endif

	START_CLOCK("MG2d::prepareBoundaryMesh");

	boundary->clearDiscretization2d();
#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::prepareBoundaryMesh - create.");
#endif
	int finished_count = 0;
	for(int i = 0; i < block_count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		assert(volume && (volume->getType() == BLOCK_DOMAIN));
		int face_count = volume->getFaceCount();
		for(int j = 0; j < face_count; j++){
			MeshDomainSurface* domain_surface = (MeshDomainSurface*)volume->getFace(j);
//			Metric2dContext mc(domain_surface->getBoundary()->getControlSpace());
			if(domain_surface->getMesh() == nullptr && domain_surface->createBoundaryMesh()){
				++finished_count;
			}
		}
	}

	STOP_CLOCK("MG2d::prepareBoundaryMesh");

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::prepareBoundaryMesh - end.");
#endif
	return finished_count;
}

int MeshGenerator2d::triangulateFaces(MeshContainer3d *boundary)
{
	int block_count = boundary->getBlocksCount();
	int finished_count = 0;

	START_CLOCK("MG2d::triangulateFaces");

	for(int i = 0; i < block_count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		assert(volume && (volume->getType() == BLOCK_DOMAIN));
		int face_count = volume->getFaceCount();
		for(int j = 0; j < face_count; j++){
			MeshDomainSurface* domain_surface = (MeshDomainSurface*)volume->getFace(j);
			if(domain_surface->getMesh() != nullptr && domain_surface->triangulate()){
				++finished_count;
			}
		}
	}

	STOP_CLOCK("MG2d::triangulateFaces");

	return finished_count;
}

MeshContainer2d* MeshGenerator2d::createInitialMesh(
	int pct, SurfaceConstPtr surface, const DRect& brect, CS2dPtr space)
{
	if(pct < 1) return nullptr;

	MeshContainer2d* mesh = new MeshContainer2d(pct);
	mesh->setSurface(surface);
	mesh->setConstrainingPhase(MeshData::CONSTRAIN_NONE);

	// Initial construction - two triangles (forming rectangle) containing all points
	DRect bounding_rect = brect;
	bounding_rect.inflate(param_initial_mesh_size_ratio);

	MeshPoint2d *p1, *p2, *p3, *p4;
	mesh->addMeshPoint(p1 = new MeshPoint2d(bounding_rect.x0,	bounding_rect.y0));
	mesh->addMeshPoint(p2 = new MeshPoint2d(bounding_rect.x1,	bounding_rect.y0));
	mesh->addMeshPoint(p3 = new MeshPoint2d(bounding_rect.x1,	bounding_rect.y1));
	mesh->addMeshPoint(p4 = new MeshPoint2d(bounding_rect.x0,	bounding_rect.y1));

	p1->setBorder(TagBorder::OUTER | TagBorder::CORNER | TagBorder::FIXED);	p1->setIntTag(TagExtended::TAG_OUTERHULL_POINT);
	p2->setBorder(TagBorder::OUTER | TagBorder::CORNER | TagBorder::FIXED);	p2->setIntTag(TagExtended::TAG_OUTERHULL_POINT);
	p3->setBorder(TagBorder::OUTER | TagBorder::CORNER | TagBorder::FIXED);	p3->setIntTag(TagExtended::TAG_OUTERHULL_POINT);
	p4->setBorder(TagBorder::OUTER | TagBorder::CORNER | TagBorder::FIXED);	p4->setIntTag(TagExtended::TAG_OUTERHULL_POINT);

	// QTree for triangle scanning
	QuadTree* quad_tree = nullptr;
	if(QuadTree::param_qtree_threshold){
		quad_tree = new QuadTree(bounding_rect);
		mesh->setSearchTree(quad_tree);
	}

	// Two initial triangles (plus incidency)
	mesh->addMeshTriangle(new MeshTriangle2d(p1, p2, p4));
	mesh->addMeshTriangle(new MeshTriangle2d(p2, p3, p4));

	p1->getEdgeToPoint(p2)->setBorder();
	p2->getEdgeToPoint(p3)->setBorder();
	p3->getEdgeToPoint(p4)->setBorder();
	p4->getEdgeToPoint(p1)->setBorder();

	if (!space)
		space = std::make_shared<ControlSpace2dIdentity>(
			std::min(brect.getDX(), brect.getDY()) * SMALL_NUMBER);
	mesh->setControlSpace(space);

	return mesh;
}

MeshContainer2d* MeshGenerator2d::createInitialMesh(
	const MeshContainer2d* boundary, const DRect* bounding_area)
{
	int points_count = boundary->getPointsCount();
	if(points_count < 1) return nullptr;

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::createInitialMesh - start.");
#endif

	START_CLOCK("MG2d::createInitialMesh");

	DRect bounding_rect;
	if(bounding_area){
		bounding_rect = *bounding_area;
	}else{
		bounding_rect = boundary->getBoundingRect();
	}

	MeshContainer2d* mesh = createInitialMesh(points_count, boundary->getSurface(),
		bounding_rect, boundary->getControlSpace());

	STOP_CLOCK("MG2d::createInitialMesh");

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::createInitialMesh - end.");
#endif

	return mesh;
}

bool MeshGenerator2d::addPointToTriangulation(Metric2dContext& mc, MeshContainer2d *mesh, 
		MeshPoint2d *point, MeshTriangle2d *containing_triangle,
		TagExtended::TagType tag_type, int tag_value)
{
	mesh->addMeshPoint(point);

	switch(param_triangulation_type){
	case MeshData::TRIANGULATE_CIRCLE:
		return addPointToTriangulationByCircle(mc, mesh, point, containing_triangle, tag_type, tag_value);
	case MeshData::TRIANGULATE_ANGLE:
		return addPointToTriangulationByAngles(mc, mesh, point, containing_triangle, true, tag_type, tag_value);
	default:
		assert(false);
	}

	return false;
}

bool MeshGenerator2d::addPointToTriangulationByCircle(Metric2dContext& mc, MeshContainer2d *mesh, 
		MeshPoint2d *point, MeshTriangle2d *containing_triangle,
		TagExtended::TagType tag_type, int tag_value)
{

	const DPoint2d& pt = point->getCoordinates();
	// Metric at the inserted point
	mc.countMetricAtPoint(pt);
	// Find and mark the triangle containg the point being inserted
	if(!containing_triangle){
		MeshTriangle2d* start_triangle = mesh->getNearTriangle(pt);
		assert(start_triangle);
		containing_triangle = start_triangle->findTriangleByNeighbours(point);
	}
	if(!containing_triangle){
		LOG4CPLUS_WARN(MeshLog::logger_console, "addPointToTriangulationByCircle: error finding the containing triangle");
		LOG4CPLUS_WARN(MeshLog::logger_console, " -> switching to linear search (for this point only) ...");
		int tct = mesh->getElementsCount();
		for(int i = 0; i < tct; i++){
			MeshElement* element = mesh->getElementAt(i);
			if(element->isPointInside(point)){
				containing_triangle = (MeshTriangle2d*)element;
				break;
			}
		}
	}
	if(!containing_triangle){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "addPointToTriangulationByCircle: error finding the containing triangle");
		return false;
	}

	if(containing_triangle->closeToVertex(mc, point)) return false;

	containing_triangle->setIntTag(TagExtended::TAG_TRI_INSPHERE);
	DataCompoundList<MeshTriangle2d*> cavity_triangles;
	cavity_triangles.append(mesh->removeMeshTriangle(containing_triangle));
	int area_id = containing_triangle->getAreaID();
	int new_edge_ct = 0;	// Number of cavity edges after removing of triangles

	// * Find (and remove) all neighbors with this point in circumcircle
	for (auto it = cavity_triangles.iterator(); it.valid(); it.moveNext()) {
		for(int j = 0; j < 3; j++){	// For all neighbors
			MeshTriangle2d* triangle = (MeshTriangle2d*)it.item()->getNeighbour(j);
			if(triangle){
				if(triangle->zeroIntTag(TagExtended::TAG_TRI_INSPHERE)){
					bool allowed = true;
					if(mesh->getConstrainingPhase() != MeshData::CONSTRAIN_NONE){
						allowed = !it.item()->isBorder(j);
					}
					if((tag_type != TagExtended::TAG_NONE) &&
						! triangle->checkIntTag(tag_type, tag_value))
						allowed = false;

					if(allowed){
#ifdef PRINT_METRIC_CALL_COUNTERS
						++stat_retriangulation_circle_check;
#endif //PRINT_METRIC_CALL_COUNTERS
						if(triangle->isPointInOuterCircle(mc, point, false)){
							triangle->setIntTag(TagExtended::TAG_TRI_INSPHERE);
							cavity_triangles.append(mesh->removeMeshTriangle(triangle));
						}else ++new_edge_ct;
					}else ++new_edge_ct;
				}
			}else ++new_edge_ct;
		}
	}

	// Set edges of cavity
	struct TEdges{
		MeshPoint2d *p1, *p2;		// wierzcho³ki wspólnej krawêdzi
		MeshTriangle2d* triangle;	// nowy trójk¹t
	} *edges = new TEdges[new_edge_ct];
	int ect = 0;
	for (auto it = cavity_triangles.iterator(); it.valid(); it.moveNext()) {
		for(int j = 0; j < 3; j++){
			MeshTriangle2d* triangle = (MeshTriangle2d*)it.item()->getNeighbour(j);
			if(!triangle){
				// Boundary edge of domain
				edges[ect].p1 = it.item()->getPoint(j);
				edges[ect].p2 = it.item()->getPoint((j+1)%3);
				++ect;
			}else if(triangle->zeroIntTag(TagExtended::TAG_TRI_INSPHERE)){
				// Non-removed neighboring triangle
				edges[ect].p1 = it.item()->getPoint(j);
				edges[ect].p2 = it.item()->getPoint((j+1)%3);
				++ect;
			}
		}
	}

	// Remove obsolete triangles
	for (auto it = cavity_triangles.iterator(); it.valid(); it.moveNext()) {
		delete it.item();
	}

	// Create new triangles
	for(int k = 0; k < ect; k++){
		edges[k].triangle = new MeshTriangle2d(edges[k].p1, edges[k].p2, point);
		if(tag_type != TagExtended::TAG_NONE)
			edges[k].triangle->setIntTag(tag_type, tag_value);

		if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE){
			edges[k].triangle->setAreaID(area_id);
			edges[k].triangle->countQuality(mc, true);
		}
		mesh->addMeshTriangle(edges[k].triangle);
	}
/*
	if(!mesh->validHeapOrder()){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MG2d: [3] invalid heap order");
		return false;
	}
*/
#ifdef PRINT_METRIC_CALL_COUNTERS
	stat_retriangulation_circle_remove += 1+cavity_triangles.countInt();
	stat_retriangulation_new_elements += ect;
#endif //PRINT_METRIC_CALL_COUNTERS

	// Free
	delete[] edges;
	return true;
}

bool MeshGenerator2d::addPointToTriangulationByAngles(Metric2dContext& mc, MeshContainer2d *mesh, 
		MeshPoint2d *point, MeshTriangle2d *containing_triangle, bool allow_swap,
		TagExtended::TagType tag_type, int tag_value)
{
	mc.countMetricAtPoint(point->getCoordinates());
	//ZnajdŸ i zaznacz trójk¹t który zawiera ten punkt
	if(!containing_triangle){
		MeshTriangle2d *start_triangle = mesh->getNearTriangle(point->getCoordinates());
		assert(start_triangle);
		containing_triangle = start_triangle->findTriangleByNeighbours(point);
	}
	if(!containing_triangle) return false;

//	char text[100];
//	sprintf(text, "Nowy punkt [%d] do retriangulacji wg kryterium k¹tów.", mesh->getPointsCount());
//	SHOW_STEP_PT_BREAKABLE(2, text, point->getCoordinates(), false);

	// split triangle in three (or four) sub-triangles
	mesh->removeMeshTriangle(containing_triangle);
	MeshTriangle2d **new_triangles = new MeshTriangle2d*[4];
	int new_ct = 3;
	int i = 0;
	int area_id1 = containing_triangle->getAreaID();
	MeshTriangle2d *triangle = containing_triangle;
	MeshTriangle2d *second_triangle = nullptr;


	for(int j = 0; j < 3; j++){
		MeshPoint2d* p1 = triangle->getPoint(j);
		MeshPoint2d* p2 = triangle->getPoint((j+1)%3);

		double area = DTriangle2d::area(
			point->getMetricCoordinates(mc), 
			p1->getMetricCoordinates(mc), 
			p2->getMetricCoordinates(mc));
		if(area > METRIC_SMALL_NUMBER){
			//assert(area > 0);
			new_triangles[i] = new MeshTriangle2d(point, p1, p2);
			mesh->addMeshTriangle(new_triangles[i], true);
			if(tag_type != TagExtended::TAG_NONE)
				new_triangles[i]->setIntTag(tag_type, tag_value);
			if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE)
				new_triangles[i]->setAreaID(area_id1);
			++i;
		}else{
			assert(area >= -METRIC_SMALL_NUMBER);
			// Point is located at the edge
			MeshEdge2d *edge = p1->getEdgeToPoint(p2);
			if(new_ct > 3 || !edge || edge->isBorder()){	//
				for(int l = 0; l < i; l++)
					delete new_triangles[l];
				delete[] new_triangles;
				return false;
			}
			second_triangle = (MeshTriangle2d*)triangle->getNeighbour(j);
			if(second_triangle){
				mesh->removeMeshTriangle(second_triangle);
				MeshPoint2d* p3 = second_triangle->getOtherPoint(p1, p2);
				int area_id2 = second_triangle->getAreaID();
				// Two new triangles
				mesh->addMeshTriangle(new_triangles[i] = new MeshTriangle2d(p1, p3, point), true);
				mesh->addMeshTriangle(new_triangles[i+1] = new MeshTriangle2d(point, p3, p2), true);

				assert(new_triangles[i]->getArea(mc, false) > 0);
				assert(new_triangles[i+1]->getArea(mc, false) > 0);

				if(tag_type != TagExtended::TAG_NONE){
					new_triangles[i]->setIntTag(tag_type, tag_value);
					new_triangles[i+1]->setIntTag(tag_type, tag_value);
				}
				// Area identifiers in triangles
				if(mesh->getConstrainingPhase() == MeshData::CONSTRAIN_DONE){
					new_triangles[i]->setAreaID(area_id2);
					new_triangles[i+1]->setAreaID(area_id2);
				}
				new_ct++;
				i += 2;
			}else{
				// outer boundary edge
				for(int l = 0; l < i; l++)
					delete new_triangles[l];
				delete[] new_triangles;
				return false;
			}
		}
	}

	delete triangle;
	if(second_triangle) delete second_triangle;

	DataCompoundList<MeshTriangle2d*> check_triangles;

	if(allow_swap){
		for(int j = 0; j < new_ct; j++){
			check_triangles.append(new_triangles[j]);
		}
	}

	DataCompoundList<MeshTriangle2d*> modified_triangles;

	if(allow_swap){ // mark edges
		for(i = 0; i < new_ct; i++)
			for(int j = 0; j < 3; j++)
				new_triangles[i]->getEdge(j)->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
	}

	if(mesh->isHeapOrder()){
		if(allow_swap){ // mark for further update
			for(i = 0; i < new_ct; i++){
				modified_triangles.insert(new_triangles[i]);
				new_triangles[i]->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
			}
		}else // update now
			for(i = 0; i < new_ct; i++)
				mesh->updateElementPosition(new_triangles[i]);
	}

	delete[] new_triangles;

	if(!allow_swap)	return true;

	int swapped = 0;

	QuadTree* search_tree = mesh->getSearchTree();

	// Swap edges with inner angles criterion...
	while(check_triangles.notEmpty()){
		triangle = check_triangles.removeFirst();
		// Check inner angles (and if quad is convex) for neighbors
		for(int j = 0; j < 3; j++){
			MeshEdge2d* common_edge = triangle->getEdge(j);
			if(!common_edge->availableTag(TagExtended::TAG_MG2D_SM_SWAP)) continue;
			second_triangle = (MeshTriangle2d*)triangle->getNeighbour(j);
			if(second_triangle && 
				second_triangle->getType() == ELEMENT_MESH_TRIANGLE &&
				((tag_type != TagExtended::TAG_NONE) && !second_triangle->checkIntTag(tag_type, tag_value))){
				// swap ?
#ifdef PRINT_METRIC_CALL_COUNTERS
				++stat_retriangulation_angles_check;
#endif //PRINT_METRIC_CALL_COUNTERS
				if(triangle->shouldSwap(mc, j, false)){
#ifdef PRINT_METRIC_CALL_COUNTERS
					++stat_retriangulation_angles_swap;
#endif //PRINT_METRIC_CALL_COUNTERS
					++swapped;
					if(search_tree){
						// update tree
						if(!triangle->availableTag(TagExtended::TAG_MG2D_SM_SWAP)) 
							search_tree->removeTriangleLink(triangle);
						if(!second_triangle->availableTag(TagExtended::TAG_MG2D_SM_SWAP)) 
							search_tree->removeTriangleLink(second_triangle);
					}
					triangle->swapWithNeighbourNoCheck(j);	// condition already checked
					if(mesh->isHeapOrder()){
						if(!triangle->availableTag(TagExtended::TAG_MG2D_SM_SWAP)){
							modified_triangles.insert(triangle);
							triangle->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
						}
						if(!second_triangle->availableTag(TagExtended::TAG_MG2D_SM_SWAP)){
							modified_triangles.insert(second_triangle);
							second_triangle->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
						}
					}
					// Add to list
					check_triangles.append(triangle);
					check_triangles.append(second_triangle);
					// Set tag for adjacent edges
					for(int k = 0; k < 3; k++){
						triangle->getEdge(k)->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
						second_triangle->getEdge(k)->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
					}
					for(int k = 0; k < 3; k++){ // remove tag for common edge for this pair
						if(triangle->getNeighbour(k) == second_triangle){
							triangle->getEdge(k)->removeTag(TagExtended::TAG_MG2D_SM_SWAP);
							break;
						}
					}
					break;
				}else{
					common_edge->removeTag(TagExtended::TAG_MG2D_SM_SWAP);
				}
			}else{
				common_edge->removeTag(TagExtended::TAG_MG2D_SM_SWAP);
			}
		}
		// Preventing infinite loop
		if(swapped > param_swap_maximum){  // default = 100000
			// After (param_swap_maximum swaps, stop (clear list)
			LOG4CPLUS_WARN(MeshLog::logger_console, "Inner angles retriangulation - reached max swaps");
			// clear tags
			for(auto it = check_triangles.iterator(); it.valid(); it.moveNext()){
				for(int j = 0; j < 3; j++)
					it.item()->getEdge(j)->removeTag(TagExtended::TAG_MG2D_SM_SWAP);
			}
			check_triangles.clear();
		}
	}

#ifdef PRINT_METRIC_CALL_COUNTERS
	stat_retriangulation_new_elements += modified_triangles.countInt();
#endif //PRINT_METRIC_CALL_COUNTERS

#ifdef STAT_SWAP_INNER_ANGLES
	stat_swap_inner_angles_counter += swapped;
#endif

	while(modified_triangles.notEmpty()){
		MeshTriangle2d* t = modified_triangles.removeFirst();
		t->countQuality(mc, true);
		t->removeTag(TagExtended::TAG_MG2D_SM_SWAP);
		mesh->updateElementPosition(t);
		if(search_tree)	search_tree->insertTriangleLink(t);
	}

	return true;
}

bool MeshGenerator2d::constrainToBorder(Metric2dContext& mc, const MeshContainer2d *boundary, 
		MeshContainer2d *mesh)
{

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::constrainToBorder - start.");
#endif


//	SHOW_STEP_BREAKABLE(1, "* Kontrola brzegu i identyfikacja elementów.", false);

	START_CLOCK("MG2d::constrainToBorder");

	int pct = mesh->getPointsCount();
	assert(pct > 3);

	mesh->setConstrainingPhase(MeshData::CONSTRAIN_ACTIVE);

	int boundary_count = boundary->getPointsCount();
	DataVector<MeshPoint2d*> bm_points(boundary_count, 0); // boundary->mesh points
	for(int i = 0; i < pct; i++){
		MeshPoint2d* mpt = mesh->getPointAt(i);
		MeshPoint2d* bpt = (MeshPoint2d*)mpt->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
		if(!bpt) continue;
		bm_points[bpt->getIndex()] = mpt;
	}

	// whether border edges are included in the triangulation
	for(int i = 0; i < boundary_count; i++){
		MeshPoint2d* p1 = boundary->getPointAt(i);
		MeshPoint2d* mp1 = bm_points[i];
		bm_points[i] = 0; // mark as already checked
		int rank = p1->getRank();
		mp1->setBorder();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = p1->getEdge(j);
			MeshPoint2d* p2 = edge->getOtherPoint(p1);
			MeshPoint2d* mp2 = bm_points[p2->getIndex()];
			if(!mp2) continue; // already checked (or missing?)

			MeshEdge2d* mesh_edge = mp1->getEdgeToPoint(mp2);
			int dir = edge->getDirection(p1);
			char border_flags = edge->getBorderFlags();
			int area1 = edge->getIncidentAreaID(MeshEdge2d::AREA_LEFT);
			int area2 = edge->getIncidentAreaID(MeshEdge2d::AREA_RIGHT);
			// Inner edge ?
			if(area1 > -1 && area2 > -1) border_flags = TagBorder::INNER;
			int mesh_index = (mesh_edge)?(mesh_edge->getPointIndex(mp1)):0;
			if(dir == 0){
				if(edge->getPointIndex(p1) != mesh_index){
					int temp = area1;
					area1 = area2;
					area2 = temp;
				}
			}
			if(mesh_edge){	// if exists fill info about this edge
				mesh_edge->setDirection(dir, mp1);
				mesh_edge->setBorder(border_flags);
				if(area1 > -1) mesh_edge->setIncidentAreaID(MeshEdge2d::AREA_LEFT, area1);
				if(area2 > -1) mesh_edge->setIncidentAreaID(MeshEdge2d::AREA_RIGHT, area2);
			}else{
				// Odtwórz krawêdŸ
				bool result = makeEdgeByEdgesSwapPIng(mc, mp1, mp2, dir, border_flags, 
						area1, area2);
				assert(result);
				if(!result){
					LOG4CPLUS_ERROR(MeshLog::logger_console,   "SWAPPING-FAILED");
					return false;
				}
				mesh_edge = mp1->getEdgeToPoint(mp2);
			}
			mesh_edge->setPtrTag(TagExtended::TAG_BOUNDARY_EDGE, edge); // mesh => boundary
			mesh_edge->setPtrTag(TagExtended::TAG_ME_2D_3D,		// mesh => boundary 3D
				edge->getPtrTag(TagExtended::TAG_ME_2D_3D));
		}
	}

	int mesh_count = mesh->getElementsCount();
	// Approved triangles
	DataVector<MeshTriangle2d*> approved(mesh_count);

	for(int i = 0; i < mesh_count; i++){
		MeshTriangle2d* triangle = (MeshTriangle2d*)mesh->getElementAt(i);
		if(triangle->getAreaID() > -1){	// "in"
			triangle->setIntTag(TagExtended::TAG_MG2D);
			approved.add(triangle);
		}
	}

	// Dodaj do listy zatwierdzonych trójk¹tów ich s¹siadów i s¹siadów ich s¹siadów ...
	for(int i = 0; i < approved.countInt(); i++){
		for(int j = 0; j < 3; j++){
			MeshTriangle2d* triangle = (MeshTriangle2d*) approved[i]->getNeighbour(j);
			if(triangle && !triangle->availableTag(TagExtended::TAG_MG2D) && !approved[i]->isBorder(j)){
				triangle->setIntTag(TagExtended::TAG_MG2D);
				approved.add(triangle);
				triangle->setAreaID(approved[i]->getAreaID());
			}
		}
	}
/*
	if(true){
		MeshViewSet *set = new MeshViewSet;
		SurfaceParametric* surface = mesh->getSurface();
		for(int i = 0; i < mesh->getElementsCount(); i++){
			MeshTriangle2d* vt = (MeshTriangle2d*)mesh->getElementAt(i);
			if(vt->isTagged()){
				set->addElement(vt, surface);
			}
		}
		SHOW_MESH("constraining", set);
	}
*/

	// Remove unmarked triangles
	for(int i = 0; i < mesh->getElementsCount(); ){
		MeshTriangle2d* triangle = (MeshTriangle2d*)mesh->getElementAt(i);
		if(triangle->availableTag(TagExtended::TAG_MG2D)){
			triangle->removeTag(TagExtended::TAG_MG2D);
			++i;
		}else{
//			SHOW_STEP_PT(3, "Usuwany niepotrzebny trójk¹t.", triangle->getMiddlePoint()); 
			delete mesh->removeMeshElement(i);
			// Uwaga: usuwanie powoduje, ¿e zmienia siê ogólna liczba elementów,
			//	a na miejsce usuniêtego wstawiany jest inny (o ile usuwany nie jest ostatni)
			//	- dlatego jeœli usuwamy, nie zmieniamy indeksu i
		}
	}

	// remove auxiliary points
	for(int i = 0; i < pct; ){
		MeshPoint2d* point = mesh->getPointAt(i);
		if(point->availableTag(TagExtended::TAG_OUTERHULL_POINT) || point->getRank() < 1){
			delete mesh->removeMeshPoint(i);
			--pct;
		}else ++i;
	}

	mesh->setConstrainingPhase(MeshData::CONSTRAIN_DONE);
	mesh->setDiscretizationState(1);

	STOP_CLOCK("MG2d::constrainToBorder");

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::constrainToBorder - end.");
#endif
	
	return true;
}

/*
CS2dPtr MeshGenerator2d::createUniformControlSpace(const SurfaceParametric *surface, 
	const MeshContainer2d *boundary)
{
	START_CLOCK("MG2d::createUniformControlSpace");
	int i, j, pct = boundary->getPointsCount();
	double min_len = mesh_data.relative_infinity;
	for(i = 0; i < pct; i++){
		MeshPoint2d* point = boundary->getPointAt(i);
		int rank = point->getRank();
		for(j = 0; j < rank; j++){
			MeshEdge2d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0)
				min_len = std::min(min_len, edge->getLength());
		}
	}

	DRect rect = boundary->getBoundingRect();
	rect.inflate(ControlSpace2d::param_inflate_box_factor);

	double ratio_x = rect.getDY() / min_len;
	double ratio_y = rect.getDX() / min_len;
	if(ratio_x > 400.0) ratio_x = 400.0;
	if(ratio_y > 400.0) ratio_y = 400.0;
	ControlSpace2dMatrixUniform* space = new ControlSpace2dMatrixUniform(surface, rect, 
		(int)(0.25 * ratio_x), (int)(0.25 * ratio_y));

	for(i = 0; i < pct; i++){
		MeshPoint2d* point = boundary->getPointAt(i);
		int rank = point->getRank();
		for(j = 0; j < rank; j++){
			MeshEdge2d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0){
				double len = edge->getLength(surface);
				space->addControlPoint(edge->getPoint(0.5), ControlDataStretch2d(len, len, 0.0));
			}
		}
	}

	space->interpolate();

	STOP_CLOCK("MG2d::createUniformControlSpace");

	return space;
}
*/

bool MeshGenerator2d::makeEdgeByEdgesSwapPIng(Metric2dContext& mc, MeshPoint2d *pt1, MeshPoint2d *pt2, 
		int dir, char border_type, int area1, int area2)
{
	MeshTriangle2d *triangle = nullptr;
	const int MAX_STEPS = 30;
	const DMPoint2d dpt1 = pt1->getMetricCoordinates(mc);
	const DMPoint2d dpt2 = pt2->getMetricCoordinates(mc);

	// Find triangle for first point
	int rank = pt1->getRank();
	int steps = 0;
	int index = -1;
	for(int i = 0; i < rank; i++){
		const MeshEdge2d* edge = pt1->getEdge(i);
		MeshTriangle2d *triangle1 = (MeshTriangle2d*)edge->getMeshElement(pt1);
		if(!triangle1) triangle1 = (MeshTriangle2d*)edge->getOtherElement(nullptr);
		if(triangle1 && triangle1->getEdgeCount() == 3){
			index = (triangle1->getPointIndex(pt1) + 1) % 3;
			assert(index > -1);
			const DMPoint2d dpt3 = triangle1->getPoint(index)->getMetricCoordinates(mc);
			const DMPoint2d dpt4 = triangle1->getPoint((index+1)%3)->getMetricCoordinates(mc);
			// If segments (1-2) and (3-4) cross ...
			double det13 = DTriangle2d::det(dpt1, dpt2, dpt3);
			double det14 = DTriangle2d::det(dpt1, dpt2, dpt4);
			if(abs(det13) < METRIC_SMALL_NUMBER || abs(det14) < METRIC_SMALL_NUMBER){
				MeshPoint2d* pt = (abs(det13) < abs(det14)) ? 
						triangle1->getPoint(index) : triangle1->getPoint((index+1)%3);
				const DMPoint2d dpt = pt->getMetricCoordinates(mc);
				double dist12 = dpt1.distance2(dpt2);
				if(dpt.distance2(dpt1) < dist12 && dpt.distance2(dpt2) < dist12){
					if(++steps > MAX_STEPS) return false;
					if(!pt->isBorder()){
						// Move slightly away from the edge being recovered
						const DVector2d v = (pt->getCoordinates() - pt1->getCoordinates()).turned(
							sin(PI/40.0), cos(PI/40.0)); // 4.5 degree
						tryMovingPoint(pt, pt1->getCoordinates() + v);
						i = -1; // start from the beginning
						continue;
					}
				}
			}
			double det2 = DTriangle2d::det(dpt3, dpt4, dpt1) * DTriangle2d::det(dpt3, dpt4, dpt2);
			// If they do cross ...
			if(det13 * det14 < 0 && det2 < 0){
				triangle = triangle1; break;
			}
		}
	}
//	assert(triangle);
	if(!triangle) return false;

	DataSimpleList<MeshEdge2d*> flip_edges;
	flip_edges.append(triangle->getEdge(index));

	int ct = 1;
	steps = 0;

	// Prepare list of edge to flip ...
	while(true){
		MeshEdge2d* edge = flip_edges.getLast();
		MeshTriangle2d* other_triangle = (MeshTriangle2d*)edge->getOtherElement(triangle);
		if(!other_triangle || edge->isBorder() || other_triangle->getEdgeCount() != 3){
			//assert(false);
			return false;
		}
		MeshPoint2d* pt3 = triangle->getPoint(index);
		int other_index = other_triangle->getPointIndex(pt3);
		assert(other_index > -1);
		MeshPoint2d* pt4 = other_triangle->getPoint((other_index+1)%3);
		if(pt4 == pt2) break; // found ending point
		const DMPoint2d dpt3 = pt3->getMetricCoordinates(mc);
		const DMPoint2d dpt4 = pt4->getMetricCoordinates(mc);
		// If segments (1-2) and (3-4) don't cross ...
		double det13 = DTriangle2d::det(dpt1, dpt2, dpt3);
		double det14 = DTriangle2d::det(dpt1, dpt2, dpt4);
		if(abs(det13) < METRIC_SMALL_NUMBER || abs(det14) < METRIC_SMALL_NUMBER)
		{
			MeshPoint2d* pt = (abs(det13) < abs(det14)) ? pt3 : pt4;
			if(++steps > MAX_STEPS) return false;
			if(!pt->isBorder()){
				const DVector2d v = (pt->getCoordinates() - pt1->getCoordinates()).turned(
					sin(PI/40), cos(PI/40)); // 4.5 degree
				tryMovingPoint(pt, pt1->getCoordinates() + v);
				// check this triangle again ...
				continue;
			}
		}
		double det2 = DTriangle2d::det(dpt3, dpt4, dpt1) * DTriangle2d::det(dpt3, dpt4, dpt2);
		if(det13 * det14 > 0 || det2 > 0){
			// select next edge
			other_index = (other_index+1)%3;
		}
		flip_edges.append(other_triangle->getEdge(other_index));
		triangle = other_triangle;
		index = other_index;
		++ct;
	}

	// Flip edges ...
	int step = 0;
	while(flip_edges.notEmpty()){
		if(++step > 3*ct){
			// Error recovering the edge
			//assert(false);
			return false;
		}
		MeshEdge2d* edge = flip_edges.removeFirst();
		MeshPoint2d* pt3 = edge->getMeshPoint(0);
		MeshPoint2d *pt4 = edge->getMeshPoint(1);
		triangle = (MeshTriangle2d*)edge->getMeshElement(pt3);
		MeshTriangle2d* other_triangle = (MeshTriangle2d*)edge->getMeshElement(pt4);
		pt3 = triangle->getPrevPoint(pt3); // opposite points
		pt4 = other_triangle->getPrevPoint(pt4);
		if(triangle->swapWithNeighbour(mc, edge, false)){
			// Edge after flip
			const DMPoint2d dpt3 = pt3->getMetricCoordinates(mc);
			const DMPoint2d dpt4 = pt4->getMetricCoordinates(mc);
			// If segments (1-2) and (3-4) cross ...
			if(flip_edges.notEmpty() &&
				DTriangle2d::det(dpt1, dpt2, dpt3) * DTriangle2d::det(dpt1, dpt2, dpt4) < 0 &&
				DTriangle2d::det(dpt3, dpt4, dpt1) * DTriangle2d::det(dpt3, dpt4, dpt2) < 0)
			{
				// Edge is still in the way - na koniec listy
				flip_edges.append(pt3->getEdgeToPoint(pt4));
			}
		}else{
			// Move to end
			flip_edges.append(edge);
		}
	}

	// Set properties
	if(dir > -3){
		MeshEdge2d* edge = pt1->getEdgeToPoint(pt2);
		assert(edge);
		edge->setDirection(dir, pt1);
		edge->setBorder(border_type);
		edge->setIncidentAreaID(MeshEdge2d::AREA_LEFT, area1);
		edge->setIncidentAreaID(MeshEdge2d::AREA_RIGHT, area2);
	}

	return true;
}

int MeshGenerator2d::addInnerNodes(Metric2dContext& mc, MeshContainer2d *mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return 0;

	if(!MeshGenerator2d::param_triangulate_with_inner_nodes) return 0;

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::addInnerNodes - start.");
#endif

	START_CLOCK("MG2d::addInnerNodes");

	int border_count = mesh->getPointsCount();

	// Calculate quality of present triangles
	int count = mesh->getElementsCount();
	if(count < 1) return 0;

	mesh->clearSearchTree();

//	double ev_nt_count = 0.0;
	for(int i = 0; i < count; i++){
		MeshElement* element = mesh->getElementAt(i);
		if((tag_type == TagExtended::TAG_NONE) ||
			(element->checkIntTag(tag_type, tag_value)))
		{
			element->countQuality(mc, true);
		}else{
			element->setQuality(mesh_data.relative_infinity);
		}
//		ev_nt_count += MeshElement::last_area;
	}
	
	// normalize NT evaluation (total metric area -> number of triangles)
//	ev_nt_count *= 4.0/SQRT3;
//	if(show_prediction) LOG4CPLUS_INFO(MeshLog::logger_console, "Initial guess NT2d", ev_nt_count);

	mesh->setHeapOrder(true);

	int new_ct = 0;
	int steps = 0;
	const DRect bounding_rect = mesh->getBoundingRect();

#ifdef STAT_SWAP_INNER_ANGLES
	stat_swap_inner_angles_counter = 0;
#endif

	//int last_nt_ev_step = -1;
	while(true){
		// Najgorszy trójk¹t z uwzglêdnieniem minimalnej powierzchni trójk¹tów znajduje siê w korzeniu
		MeshTriangle2d *triangle = (MeshTriangle2d*)mesh->getElementAt(0);

		// Is quality low enough
		double quality = triangle->getQuality();
		//LOG4CPLUS_INFO(MeshLog::logger_console, "Q", quality);
		if(quality > param_quality_threshold) break;

		//int current_nt_ev_step = (int)(20 * quality);
		//if(show_prediction && current_nt_ev_step > last_nt_ev_step && quality >= 0.01){
		//	int nt2d = mesh->getElementsCount();
		//	double nt_expected = nt2d * param_quality_threshold / quality;
		//	ostringstream os;
		//	os << "NT2d= " << setw(7) << nt2d << ", qt= " << fixed << setprecision(2) << quality << ", predicted= " << setw(7) << (int)nt_expected;
		//	LOG4CPLUS_INFO(MeshLog::logger_console, os.str());

		//	last_nt_ev_step = current_nt_ev_step;
		//}

		++steps;
		mc.countMetricAtPoint(triangle->getMiddlePoint());

		DMPoint2d dpoints[3];
		for(int i = 0; i < 3; i++){
			dpoints[i] = triangle->getPoint(i)->getMetricCoordinates(mc);
		}
		DPoint2d dnew;
		int qi_method = param_quality_improvement;
		if(qi_method == MeshData::IMPROVE_LONGEST_EDGE_AND_OUTER_CIRCLE)
			if(triangle->getQuality() < 0.05)
				qi_method = MeshData::IMPROVE_LONGEST_EDGE;
			else qi_method = MeshData::IMPROVE_OUTER_CIRCLE;

		switch(qi_method){
		case MeshData::IMPROVE_LONGEST_EDGE:
			{
				if(!countNewNodeLongestEdge(mc, dpoints, dnew, triangle)){
					triangle->setQuality(mesh_data.relative_infinity);
					mesh->updateElementPosition(triangle);
					continue;
				}
				if((tag_type != TagExtended::TAG_NONE) &&
					triangle && !triangle->checkIntTag(tag_type, tag_value))
				{
					triangle->setQuality(mesh_data.relative_infinity);
					mesh->updateElementPosition(triangle);
					continue;
				}
				break;
			}
		case MeshData::IMPROVE_OUTER_CIRCLE:
			{
				MeshTriangle2d* main_triangle = triangle;
				DMPoint2d dnew_metric = DTriangle2d::outerCircleCenter(dpoints[0], dpoints[1], dpoints[2]);
				dnew = mc.transformMStoPS(dnew_metric);
				if(!bounding_rect.contains(dnew)){
					INC_COUNT_PT_OUTSIDE_BOX;
					main_triangle->setQuality(mesh_data.relative_infinity);
					mesh->updateElementPosition(main_triangle);
					continue;
				}
				if(!triangle->isPointInside(dnew)){
					triangle = triangle->findTriangleByNeighbours(dnew, true);
					if(!triangle){
						INC_COUNT_PT_NO_TETRA_FOUND;
						main_triangle->setQuality(mesh_data.relative_infinity);
						mesh->updateElementPosition(main_triangle);
						continue;
					}
				}
				if((tag_type != TagExtended::TAG_NONE) &&
					triangle && !triangle->checkIntTag(tag_type, tag_value))
				{
					main_triangle->setQuality(mesh_data.relative_infinity);
					mesh->updateElementPosition(main_triangle);
					continue;
				}

				int i;
				for(i = 0; i < 3; i++){
					MeshEdge2d* edge = triangle->getEdge(i);
					if(edge->isBorder()){
						// Jeœli trójk¹t by³by zbyt "w¹ski" przy brzegu, zaniechaj wstawiania
						const DMPoint2d dp0 = triangle->getPoint(i)->getMetricCoordinates(mc);
						const DMPoint2d dp1 = triangle->getPoint((i+1)%3)->getMetricCoordinates(mc);
						double ratio = (dp0.distance2(dp1) * SQRT3 * 0.25) / 	
							DTriangle2d::area(dp0, dp1, dnew_metric);
						if(ratio > 1.5){
//							sprintf(text, "Przewidywany trójk¹t [opt / real = %G]", ratio);
//							SHOW_STEP_PT(2, text, dnew);
							INC_COUNT_PT_TOO_NEAR_BOUNDARY;
							main_triangle->setQuality(mesh_data.relative_infinity);
							mesh->updateElementPosition(main_triangle);
							break;
						}						
					}
				}
				if(i < 3) continue;	// jeœli powy¿sza pêtla zosta³a przerwana

				if(main_triangle != triangle){
					mc.countMetricAtPoint(dnew);
//					assert(triangle->isPointInOuterCircle(dnew)); // since dnew is inside the triangle !!!
					if(!main_triangle->isPointInOuterCircle(mc, dnew, false)){
//						SHOW_STEP_PT(2, "Zbyt ró¿ne metryki dla wstawienia punktu.", dnew);
						INC_COUNT_PT_TOO_DIFF_METRICS;
						if(!countNewNodeLongestEdge(mc, dpoints, dnew, main_triangle)){
							main_triangle->setQuality(mesh_data.relative_infinity);
							mesh->updateElementPosition(main_triangle);
							continue;
						}else triangle = main_triangle;
					}
					dnew_metric = mc.transformPStoMS(dnew);
					for(i = 0; i < 3; i++){
						if(dnew_metric.distance2(triangle->getPoint(i)->getMetricCoordinates(mc)) < 0.25){
//							SHOW_STEP_PT(2, "Nowy punkt zbyt blisko innego.", dnew);
							INC_COUNT_PT_TOO_NEAR_OTHER_PT;
							break;
						}
					}
					if(i < 3){
						if(!countNewNodeLongestEdge(mc, dpoints, dnew, main_triangle)){
							main_triangle->setQuality(mesh_data.relative_infinity);
							mesh->updateElementPosition(main_triangle);
							continue; // jeœli powy¿sza pêtla zosta³a przerwana
						}else triangle = main_triangle;
					}
					INC_COUNT_PT_IN_OTHER_TETRA;
				}
				break;
			}
		}

		MeshPoint2d* p3 = new MeshPoint2d(dnew);

		if(tag_type != TagExtended::TAG_NONE)
			p3->setIntTag(tag_type, tag_value);

//		LOG4CPLUS_INFO(MeshLog::logger_console, "=== step", steps);
		if(!addPointToTriangulation(mc, mesh, p3, triangle, tag_type, tag_value)){
			delete p3;
			mesh->setHeapOrder(false);
			return new_ct;
		}
//		LOG4CPLUS_INFO(MeshLog::logger_console, "inserted OK");

		if (mesh->getPointsCount() >= MeshPoint2d::param_max_node_count) {
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Maximum number of mesh nodes exceeded.");
			return 0;
		}

		++new_ct;
		count = border_count + new_ct;

#ifdef STAT_COUNT
		MeshData::StatData stats;
		if(getIncidenceInfo(mesh, stats)){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "inner\t" << count << "\t" << stats.average << "\t" << stats.maximum);
		}
#endif

	}

	STOP_CLOCK("MG2d::addInnerNodes");

	mesh->setHeapOrder(false);

#ifdef STAT_SWAP_INNER_ANGLES
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average diagonal swaps (inner): " << (double)stat_swap_inner_angles_counter/new_ct);
#endif

//	SHOW_STATUS_END(count);

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::addInnerNodes - end.");
#endif

	LOG4CPLUS_TRACE(MeshLog::logger_mesh, "MG2d::addInnerNodes (" << new_ct << ") done.");

	return new_ct;
}

bool MeshGenerator2d::countNewNodeLongestEdge(Metric2dContext& mc, const DMPoint2d dpoints[], 
	DPoint2d& dnew, MeshTriangle2d* triangle)
{
	// Find longest edge
	int imax = -1;
	double d2, dmax = 0.0;
	for(int i = 0; i < 3; i++){
		if((d2 = dpoints[i].distance2(dpoints[(i+1)%3])) > dmax){
			dmax = d2;
			imax = i;
		}
	}
	// New node in the middle
	MeshEdge2d* edge = triangle->getEdge(imax);
	// Jeœli punkt znajduje siê na krawêdzi -> anuluj wstawianie
	if(edge->isBorder()) return false;
	dnew = edge->getPoint(0.5);
	const DMPoint2d dnew_metric = mc.transformPStoMS(dnew);
	double min_dist = dnew_metric.distance2(dpoints[0]);
	d2 = dnew_metric.distance2(dpoints[1]);
	if(d2 < min_dist) min_dist = d2;
	d2 = dnew_metric.distance2(dpoints[2]);
	if(d2 < min_dist) min_dist = d2;
	if(min_dist < 0.5) return false;
	return true;
}


bool MeshGenerator2d::smoothenMetric(Metric2dContext& mc, MeshContainer2d* mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return true;
//	START_CLOCK("MG2d::smoothenMetric");

	for(int i = 0; i < count; i++){
		MeshPoint2d* point = mesh->getPointAt(i);
		if(point->isBorder()) continue;
		if((tag_type != TagExtended::TAG_NONE) &&
			! point->checkIntTag(tag_type, tag_value)) continue;

		mc.countMetricAtPoint(point->getCoordinates());
		int rank = point->getRank();
		DPoint2d ave_point;
		int ave_count = 0;
		for(int j = 0; j < rank; j++){
			MeshTriangle2d* triangle = (MeshTriangle2d*)point->getEdge(j)->getMeshElement(point);
			if(!triangle || triangle->getEdgeCount() != 3) continue;
			MeshPoint2d* point1 = triangle->getNextPoint(point);
			MeshPoint2d* point2 = triangle->getNextPoint(point1);
			const DPoint2d edge_middle = DPoint2d::average(point1->getCoordinates(), point2->getCoordinates());
			//mc.countMetricAtPoint(edge_middle);
			const DMVector2d mdv = point2->getMetricCoordinates(mc) - point1->getMetricCoordinates(mc);
			double len2 = mdv.length2();
			double h2 = 1.0 - 0.25*len2;
			const DMVector2d mnv = DMVector2d(-mdv.y, mdv.x) * sqrt(h2/len2);
			ave_point.add(mc.transformMStoPS(mc.transformPStoMS(edge_middle) + mnv));
			ave_count++;
		}
		if(ave_count > 0){
			ave_point /= ave_count;
			tryMovingPoint(point, ave_point);
		}
	}

//	STOP_CLOCK("MG2d::smoothenMetric");
	return true;
}

bool MeshGenerator2d::smoothenLaplace(Metric2dContext& mc, MeshContainer2d *mesh, bool variable_metric,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return true;

//	START_CLOCK("MG2d::smoothenLaplace");

	for(int i = 0; i < count; i++){
		MeshPoint2d* point = mesh->getPointAt(i);
		if((tag_type != TagExtended::TAG_NONE) &&
			! point->checkIntTag(tag_type, tag_value)) continue;

		if(variable_metric) movePointByLaplaceForVariableMetric(mc, point);
		else movePointByLaplace(point);
	}

//	STOP_CLOCK("MG2d::smoothenLaplace");
	return true;
}

bool MeshGenerator2d::smoothenPostCheck(Metric2dContext& mc, MeshContainer2d *mesh)
{
	if(!mesh) return false;
	int ect = mesh->getElementsCount();
	if(ect < 1) return true;

//	START_CLOCK("MG2d::smoothenPostCheck");

	for(int i = 0; i < ect; i++){
		MeshElement* element = mesh->getElementAt(i);
		if(!element) continue;
		mc.countMetricAtPoint(element->getMiddlePoint());
		double last_aq = element->getAlphaQuality(mc, false);
		if(element->isInverted() || last_aq < 0.1){
//			LOG4CPLUS_INFO(MeshLog::logger_console, "postcheck element, alpha quality", last_aq);
//			MeshView::showDebugMesh("bad element", mesh, element);
			int edge_ct = element->getEdgeCount();
			DataVector<MeshPoint2d*> free_points(edge_ct);
			for(int j = 0; j < edge_ct; j++){
				MeshPoint2d* point = element->getPoint(j);
				if(!point->isBorder()) free_points.add(point);
			}
			while(free_points.countInt() > 1){
				// select point with greatest increase of quality for the bad element
				double max_diff_aq = 0.0;
				int best_j = -1;
				for(int j = 0; j < free_points.countInt(); j++){
					MeshPoint2d* point = free_points[j];
					const DPoint2d old_pt = point->getCoordinates();
					// try moving
					movePointByLaplace(point);
					double diff_aq = element->getAlphaQuality(mc, false) - last_aq;
					if(diff_aq > max_diff_aq){
						max_diff_aq = diff_aq;
						best_j = j;
					}
					// move back
					point->setCoordinates(old_pt);
				}
				if(best_j >= 0){
					movePointByLaplace(free_points[best_j]);
					free_points.removeAt(best_j);
					last_aq += max_diff_aq;
//					LOG4CPLUS_INFO(MeshLog::logger_console, "postcheck improvement, aq", element->getAlphaQuality(mc, false));
//					MeshView::showDebugMesh("bad element", mesh, element);
				}else break;
			}
			if(free_points.countInt() == 1)
				movePointByLaplace(free_points[0]);
//			LOG4CPLUS_INFO(MeshLog::logger_console, "postcheck smoothing finished, aq", element->getAlphaQuality(mc, false));
//			MeshView::showDebugMesh("finished", mesh, element);
		}
	}

//	STOP_CLOCK("MG2d::smoothenPostCheck");
	return true;
}

bool MeshGenerator2d::smoothenLaplaceMixed(Metric2dContext& mc, MeshContainer2d *mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return true;

//	START_CLOCK("MG2d::smoothenLaplaceMixed");

	for(int i = 0; i < count; i++){
		MeshPoint2d* point = mesh->getPointAt(i);
		if((tag_type != TagExtended::TAG_NONE) &&
			! point->checkIntTag(tag_type, tag_value)) continue;

		double mg = mc.getMetricGradation(point->getCoordinates());
		if(mg > 1.5) movePointByLaplaceForVariableMetric(mc, point);
		else movePointByLaplace(point);
	}

//	STOP_CLOCK("MG2d::smoothenLaplaceMixed");
	return true;
}

bool MeshGenerator2d::tryMovingPoint(MeshPoint2d *point, const DPoint2d& new_pt)
{
	if(point->isBorder()) return false;

	const DPoint2d old_pt = point->getCoordinates();
	const DVector2d vector = new_pt - old_pt;
	int rank = point->getRank();

	double factor = 1.0;
	for(int k = 0; k < 10; k++){	// 10 tries
		point->setCoordinates(old_pt + vector*factor);
		// Check validity of adjacent mesh
		bool valid = true;
		for(int j = 0; valid && j < rank; j++){
			// Get element to the left of this edge
			const MeshElement* element = point->getEdge(j)->getMeshElement(point);
			if(!element) element = point->getEdge(j)->getOtherElement(element);
			if(!element) continue;
			valid = !element->isInverted();
		}
		if(valid) return true;
		else factor *= 0.5;
	}
	// 10x fault -> cancel movement
	point->setCoordinates(old_pt);
	return false;
}

bool MeshGenerator2d::movePointByLaplace(MeshPoint2d *point)
{
	if(point->isBorder()) return false;

	DPoint2d new_pt(0.0, 0.0);
	int rank = point->getRank();
	double weight = 0.0;
	for(int j = 0; j < rank; j++){
		const MeshEdge2d* edge = point->getEdge(j);
		const DPoint2d pt = edge->getOtherPoint(point)->getCoordinates();
		double w = edge->getElementsWeight();
		new_pt.add(pt, w);
		weight += w;
	}
	new_pt /= weight;
	return tryMovingPoint(point, new_pt);
}

//bool MeshGenerator2d::movePointByLaplaceForVariableMetric(Metric2dContext& mc, MeshPoint2d *point)
//{
//	if(point->isBorder()) return false;
////	SHOW_STEP_PT(3, "[movePointByLaplace].", point->getCoordinates());
//	DMVector2d total_mdv(0.0, 0.0);
//	int rank = point->getRank();
//	double total_weight = 0.0;
//	for(int j = 0; j < rank; j++){
//		const MeshEdge2d* edge = point->getEdge(j);
//		mc.countMetricAtPoint(edge->getPoint(0.5));
//		double w = edge->getElementsWeight();
//		const DMVector2d mdv = edge->getOtherPoint(point)->getMetricCoordinates(mc) 
//			- point->getMetricCoordinates(mc);
//		total_mdv += mdv * w;
//		total_weight += w;
//	}
//	total_mdv /= total_weight;
//	mc.countMetricAtPoint(point->getCoordinates());
//	return tryMovingPoint(point, point->getCoordinates() + mc.transformMStoPS(total_mdv));
//}

bool MeshGenerator2d::movePointByLaplaceForVariableMetric(Metric2dContext& mc, MeshPoint2d *point)
{
	if(point->isBorder()) return false;

	const DPoint2d& dpt = point->getCoordinates();
	int rank = point->getRank();
	DataVector< DVector2d > dvs(rank);
	DataVector< double > lens(rank);
	double ave_len = 0.0;
	for(int j = 0; j < rank; j++){
		MeshEdge2d* edge = point->getEdge(j);
		MeshPoint2d* other_point = edge->getOtherPoint(point);
		const DPoint2d& other_dpt = other_point->getCoordinates();
		mc.countMetricAtPoints(point, other_point);
		
		double len = mc.transformPStoMS( other_dpt - dpt).length();
		lens.add( len );
		dvs.add( other_dpt - dpt );
		ave_len += len;
	}
	ave_len *= (0.9 / rank); // calculate average length of edges and make it a little shorter

	DVector2d total_dv;
	double total_w = 0.0;
	for(int j = 0; j < rank; j++){
		double len_inv = 1.0 / lens[j];
		double ratio = (lens[j] - ave_len) * len_inv;
		DVector2d dv = dvs[j] * ratio;
		double w = point->getEdge(j)->getElementsWeight() * len_inv;
		total_dv += dv * w;
		total_w += w;
	}
	total_dv /= total_w;
	mc.countMetricAtPoint(dpt);

	return tryMovingPoint(point, dpt + total_dv);
}

bool MeshGenerator2d::smoothenTopologicalSwap(MeshContainer2d *mesh,
		TagExtended::TagType tag_type, int tag_value)
{

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::smoothenTopologicalSwap - start.");
#endif
	
//	LOG("\n*** Wyg³adzanie topologiczne [zamiana] ***\n", 0);
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return false;
//	SHOW_STEP_BREAKABLE(1, "* Zamiana krawêdzi (wyrównanie rzêdów wierzcho³ków).", 0);

	ControlSpace2dIdentity csi;
	Metric2dContext temp_mc(&csi);

//	START_CLOCK("MG2d::smoothenTopologicalSwap");

	for(int i = 0; i < count; i++){
		MeshPoint2d* point = mesh->getPointAt(i);
		if(point->isBorder()) continue;
		if((tag_type != TagExtended::TAG_NONE) &&
			! point->checkIntTag(tag_type, tag_value)) continue;

		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0 && !edge->isBorder()){
				MeshTriangle2d* triangle1 = (MeshTriangle2d*)edge->getMeshElement(0);
				MeshTriangle2d* triangle2 = (MeshTriangle2d*)edge->getMeshElement(1);
				if(triangle1->getEdgeCount() != 3 || triangle2->getEdgeCount() != 3){
					continue;
				}
				if((tag_type != TagExtended::TAG_NONE) &&
					(!(triangle1->checkIntTag(tag_type, tag_value)) ||
					 !(triangle2->checkIntTag(tag_type, tag_value))))
					 continue;
				MeshPoint2d* p1 = edge->getOtherPoint(point);
				if(p1->isBorder()) continue;
				MeshPoint2d* p2 = triangle1->getPrevPoint(point);
				if(p2->isBorder()) continue;
				MeshPoint2d* p3 = triangle2->getNextPoint(point);
				if(p3->isBorder()) continue;
				int diff = rank + p1->getRank() - p2->getRank() - p3->getRank();
				if(diff > 2){
//					char text[100];
//					sprintf(text, "KrawêdŸ do zamiany (ró¿nica %d).", diff);
//					SHOW_STEP_PT_BREAKABLE(2, text, edge->getPoint(0.5), operation_count);	
					triangle1->swapWithNeighbour(temp_mc, edge, false);
					rank--;
				}
			}
		}
	}

	for(int i = 0; i < count; i++){
		MeshPoint2d* point = mesh->getPointAt(i);
		if(point->isBorder()) continue;
		if((tag_type != TagExtended::TAG_NONE) &&
			! point->checkIntTag(tag_type, tag_value)) continue;

		int rank = point->getRank();
		if(rank < 5){
			MeshGenerator2d::removeTriangulationPoint(temp_mc, mesh, point);
			count = mesh->getPointsCount();
			i = 0;
		}
	}
//	STOP_CLOCK("MG2d::smoothenTopologicalSwap");

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::smoothenTopologicalSwap - end.");
#endif
	
	return true;
}

bool MeshGenerator2d::smoothenDelaunaySwapComplete(Metric2dContext& mc, MeshContainer2d *mesh,
		TagExtended::TagType tag_type, int tag_value) //THROWS_EXCEPTION
{

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::smoothenDelaunaySwapComplete - start.");
#endif
	
	if(!mesh) return false;
	int tri_ct = mesh->getElementsCount(3);
	if(tri_ct < 1) return false;
	int all_ct = mesh->getElementsCount();

	int total_count = 0;

//	START_CLOCK("MG2d::smoothenDelaunaySwapComplete");

	DataVector<MeshTriangle2d*> active_triangles(tri_ct);
	for(int i = 0; i < all_ct; i++){
		MeshTriangle2d* triangle = (MeshTriangle2d*)mesh->getElementAt(i);
		if(triangle->getEdgeCount() != 3) continue;
		if((tag_type != TagExtended::TAG_NONE) &&
			! triangle->checkIntTag(tag_type, tag_value)) continue;

		triangle->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
		triangle->setTagForEdges(TagExtended::TAG_MG2D_SM_SWAP);
		active_triangles.add(triangle);
	}

//	int counter = 0;
	size_t active_count[3] = { std::string::npos };
	while(active_triangles.countInt() > 0){
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Del-swap-complete, step " << ++counter << ", active " << active_triangles.countInt());
		for(size_t i = 0; i < active_triangles.countInt(); ){
			MeshTriangle2d * triangle = active_triangles[i];
			triangle->setZeroTag(TagExtended::TAG_MG2D_SM_SWAP);
			for(int j = 0; j < 3; j++){
				MeshEdge2d* edge = triangle->getEdge(j);
				if(edge->zeroIntTag(TagExtended::TAG_MG2D_SM_SWAP) || edge->isBorder()) continue;
				MeshTriangle2d* other_triangle = (MeshTriangle2d*)edge->getOtherElement(triangle);
				if(triangle->swapWithNeighbour(mc, j, true, true, TagExtended::TAG_MG2D_SM_SWAP)){
					++total_count;
					triangle->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
					if(other_triangle->zeroIntTag(TagExtended::TAG_MG2D_SM_SWAP)){
						if((tag_type == TagExtended::TAG_NONE) ||
							other_triangle->checkIntTag(tag_type, tag_value))
						{
							other_triangle->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
							active_triangles.add(other_triangle);
						}
					}
					break;
				}else edge->setZeroTag(TagExtended::TAG_MG2D_SM_SWAP);
			}
			if(triangle->zeroIntTag(TagExtended::TAG_MG2D_SM_SWAP)) 
				active_triangles.removeAt(i);
			else ++i;
		}
		active_count[0] = active_count[1];
		active_count[1] = active_count[2];
		active_count[2] = active_triangles.countInt();
		if(active_count[2] == active_count[1] || active_count[2] == active_count[0]) break;
	}

	mesh->removeAllTags(TagExtended::TAG_MG2D_SM_SWAP);
//	STOP_CLOCK("MG2d::smoothenDelaunaySwapComplete");
	LOG4CPLUS_TRACE(MeshLog::logger_mesh, "Del-swapped edges: " << total_count);

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::smoothenTopologicalSwapComplete - end.");
#endif
	
	return true;
}

bool MeshGenerator2d::smoothenDelaunaySwap(Metric2dContext& mc, MeshContainer2d *mesh,
		TagExtended::TagType tag_type, int tag_value) //THROWS_EXCEPTION
{

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::smoothenDelaunaySwap - start.");
#endif
	
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return false;
//	SHOW_STEP_BREAKABLE(1, "* Zamiana krawêdzi.", 0);

	int total_count = 0;

//	START_CLOCK("MG2d::smoothenDelaunaySwap");

	for(int i = 0; i < count; i++){
		const MeshPoint2d* point = mesh->getPointAt(i);
		if((tag_type != TagExtended::TAG_NONE) &&
			! point->checkIntTag(tag_type, tag_value)) continue;

		int rank = point->getRank();
		for(int j = 0; j < rank; ){
			const MeshEdge2d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0 && !edge->isBorder()){
				MeshTriangle2d* triangle = (MeshTriangle2d*)edge->getMeshElement(0);
				if(triangle->getEdgeCount() == 3 && triangle->swapWithNeighbour(mc, edge, true, true)){
					rank--;
					total_count++;
				}else{
					j++;
				}
			}else j++;
		}
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "Del-swapped edges: " << total_count);

//	STOP_CLOCK("MG2d::smoothenDelaunaySwap");

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::smoothenTopologicalSwap - end.");
#endif
	
	return true;
}

bool MeshGenerator2d::smoothenFaces(MeshContainer3d* boundary, int steps,
		TagExtended::TagType tag_type, int tag_value, int method)
{
	START_CLOCK("MG2d::smoothenFaces");
	int block_count = boundary->getBlocksCount();
	for(int i = 0; i < block_count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		assert(volume && (volume->getType() == BLOCK_DOMAIN));
		int face_count = volume->getFaceCount();
		for(int j = 0; j < face_count; j++){
			MeshDomainSurface* domain_surface = (MeshDomainSurface*)volume->getFace(j);
			if(domain_surface->getMesh() != nullptr)
				domain_surface->smoothen(steps, tag_type, tag_value, method);
		}
	}
	STOP_CLOCK("MG2d::smoothenFaces");

//	boundary->logMetricQuality();

	return true;
}

//////////////////////////////////////////////////////////////////////
// Removes point (and adjacent triangles) from the mesh, updates mesh
// by inserting new triangles within thus created cavity
// Elements adjacent to point have to be triangular.
// Delaunay criterion is not guaranteed.
bool MeshGenerator2d::removeTriangulationPoint(Metric2dContext& mc, MeshContainer2d* mesh, MeshPoint2d *point) 
{
	assert(point);
	assert(point->getElementsCount(4) == 0);
	int rank = point->getRank();
	if(rank > 4 || rank < 3) return false;	// for now, only these simple cases are implemented
	MeshPoint2d *points[4];
	// Remove adjacent triangles (with edges), in anti-clockwise order
	MeshEdge2d* edge = point->getEdge(0);
	MeshElement* element = edge->getMeshElement(point);
	assert(element);
	int id = element->getAreaID();
	points[0] = edge->getOtherPoint(point);
	assert(points[0] == element->getNextPoint(point));
	for(int i = 1; i < rank; i++){
		points[i] = element->getNextPoint(points[i-1]);
		edge = point->getEdgeToPoint(points[i]);
		assert(edge);
		MeshElement* last_element = element;
		element = edge->getOtherElement(last_element);
		assert(element);
		delete mesh->removeMeshElement(last_element);
	}
	delete mesh->removeMeshElement(element);
	// Remove the point
	assert(point->getRank() == 0);
	delete mesh->removeMeshPoint(point);

	// Create new triangles
	if(rank == 3){
		MeshTriangle2d* triangle = new MeshTriangle2d(points[0], points[1], points[2]);
		triangle->setAreaID(id);
		mesh->addMeshElement(triangle);
	}else if(rank == 4){
		// in case of concave quadrilateral, one must choose proper two triangles
		const DMPoint2d dpt0 = points[0]->getMetricCoordinates(mc);
		const DMPoint2d dpt1 = points[1]->getMetricCoordinates(mc);
		const DMPoint2d dpt2 = points[2]->getMetricCoordinates(mc);
		const DMPoint2d dpt3 = points[3]->getMetricCoordinates(mc);
		double area1_1 = DTriangle2d::det(dpt0, dpt1, dpt2);
		double area1_2 = DTriangle2d::det(dpt2, dpt3, dpt0);
		double area2_1 = DTriangle2d::det(dpt1, dpt2, dpt3);
		double area2_2 = DTriangle2d::det(dpt3, dpt0, dpt1);
		MeshTriangle2d *triangle1, *triangle2;
		if(std::min(area1_1, area1_2) > std::min(area2_1, area2_2)){
			assert((area1_1 > 0.0) && (area1_2 > 0.0));
			triangle1 = new MeshTriangle2d(points[0], points[1], points[2]);
			triangle2 = new MeshTriangle2d(points[2], points[3], points[0]);
		}else{
			assert((area2_1 > 0.0) && (area2_2 > 0.0));
			triangle1 = new MeshTriangle2d(points[1], points[2], points[3]);
			triangle2 = new MeshTriangle2d(points[3], points[0], points[1]);
		}
		triangle1->setAreaID(id);
		mesh->addMeshElement(triangle1);
		triangle2->setAreaID(id);
		mesh->addMeshElement(triangle2);
	}
	assert(mesh->isValid());
	return true;
}

//////////////////////////////////////////////////////////////////////
// Tries to optimize locally the quality of the mesh by swapPIng edges
// of the given triangle
bool MeshGenerator2d::optimizeTriangleBySwapPIng(Metric2dContext& mc, MeshTriangle2d *triangle)
{
	if(!triangle || triangle->getEdgeCount() != 3) return false;
	if(triangle->swapWithNeighbour(mc, 0, true)) return true;
	if(triangle->swapWithNeighbour(mc, 1, true)) return true;
	return triangle->swapWithNeighbour(mc, 2, true);
}

int MeshGenerator2d::addBoundaryNodes(Metric2dContext& mc, MeshContainer2d *mesh, MeshContainer2d* boundary_mesh) 
{
	if(!mesh || !boundary_mesh) return 0;

//----------------------------------------------
/*
	// ********** show				
	SurfaceParametric* surface_xxx = boundary_mesh->getSurface();
	boundary_mesh->setSurface(nullptr);
	if(MeshView::isRunning())
		MeshView::setMeshToUpdate(boundary_mesh);
	else
		_beginthread(MeshView::startLoopGL2D, 0, boundary_mesh);
	cout << "Press enter...";
	string text;
	cin >> text;
	boundary_mesh->setSurface(surface_xxx);
*/
//----------------------------------------------

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::addBoundaryNodes - start.");
#endif

	START_CLOCK("MG2d::addBoundaryNodes");
	
	int points_count = boundary_mesh->getPointsCount();
	srand(0);

#ifdef PRINT_METRIC_CALL_COUNTERS
	ofstream ofs("stat-boundary-mc.txt");
	ofs << "pt\t" << "Nc\t" << "Nmp\t" << "Npm\t" << "Npmc" 
		<< "\tan-ch" << "\tan-sw" << "\tcr-ch" << "\tcr-rm" << "\tnew_tr" 
		<< "\ttri_search" << "\tqt-max-level" << endl;
	mc.clearCounters();
	ofs << "0\t" << mc.getCallCounter() << "\t" 
		<< mc.getMPTransformationCounter() << "\t" 
		<< mc.getPMTransformationCounter() << "\t" 
		<< mc.getPMCTransformationCounter() << "\t"
		<< (stat_retriangulation_angles_check = 0) << "\t"
		<< (stat_retriangulation_angles_swap = 0) << "\t"
		<< (stat_retriangulation_circle_check = 0) << "\t"
		<< (stat_retriangulation_circle_remove = 0) << "\t"
		<< (stat_retriangulation_new_elements = 0) << "\t"
		<< MeshTriangle2d::getAndClearCounter(MeshTriangle2d::TRIANGLE_SEARCHING) << "\t"
		<< mesh->getMaxSearchTreeLevel() << endl;
#endif //PRINT_METRIC_CALL_COUNTERS

#ifdef STAT_SWAP_INNER_ANGLES
	stat_swap_inner_angles_counter = 0;
#endif

	DataVector<int> bpoints_to_remove;
	// * Insert boundary nodes into triangulation
	for(int i = 0; i < points_count; i++){
		// Select in random one of the remaining nodes
		int j = i + rand() % (points_count - i);
		if(j != i) boundary_mesh->switchMeshPoints(i, j);

		// Create new point as a copy of the given one
		MeshPoint2d *boundary_point = boundary_mesh->getPointAt(i);
		MeshPoint2d *point = new MeshPoint2d(boundary_point);
		point->setPtrTag(TagExtended::TAG_BOUNDARY_POINT, boundary_point);
		//boundary_point->setPtrTag(TagExtended::TAG_BOUNDARY_POINT, point);
		point->setPtrTag(TagExtended::TAG_MP_2D_3D, 
			boundary_point->getPtrTag(TagExtended::TAG_MP_2D_3D));

		// Retriangulate
		if(!MeshGenerator2d::addPointToTriangulation(mc, mesh, point)){
			if(boundary_point->getRank() == 0){ // which means, it is freepoint
				delete point;
				bpoints_to_remove.add(i);
			}else
				return 0;
		}

#ifdef STAT_COUNT
		MeshData::StatData stats;
		if(getIncidenceInfo(mesh, stats)){
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "boundary\t" << i << "\t" << stats.average 
				<< "\t" << stats.maximum << endl;
		}
#endif

#ifdef PRINT_METRIC_CALL_COUNTERS
	ofs << (i+1) << "\t" 
		<< mc.getCallCounter() << "\t" 
		<< mc.getMPTransformationCounter() << "\t" 
		<< mc.getPMTransformationCounter() << "\t" 
		<< mc.getPMCTransformationCounter() << "\t"
		<< stat_retriangulation_angles_check << "\t"
		<< stat_retriangulation_angles_swap << "\t"
		<< stat_retriangulation_circle_check << "\t"
		<< stat_retriangulation_circle_remove << "\t"
		<< stat_retriangulation_new_elements << "\t"
		<< MeshTriangle2d::getAndClearCounter(MeshTriangle2d::TRIANGLE_SEARCHING) << "\t"
		<< mesh->getMaxSearchTreeLevel() << endl;
#endif //PRINT_METRIC_CALL_COUNTERS

	}

	while(!bpoints_to_remove.empty()){
		delete boundary_mesh->removeMeshPoint(bpoints_to_remove.removeLast());
	}

#ifdef STAT_SWAP_INNER_ANGLES
	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Average diagonal swaps (boundary): " 
		<< (double)stat_swap_inner_angles_counter/points_count << endl;
#endif

	STOP_CLOCK("MG2d::addBoundaryNodes");

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2d::addBoundaryNodes - end.");
#endif
	
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Boundary-nodes:\t" << points_count);

	return points_count;
}

bool MeshGenerator2d::getIncidenceInfo(MeshContainer2d *mesh, MeshData::StatData & stats)
{
	stats.average = 0;
	stats.minimum = 1000;
	stats.maximum = -1;
	int pct = mesh->getPointsCount();
	if(pct < 1) return false;
	for(int i = 0; i < pct; i++){
		int rank = mesh->getPointAt(i)->getRank();
		if(rank < stats.minimum) stats.minimum = rank;
		if(rank > stats.maximum) stats.maximum = rank;
		stats.average += rank;
	}
	stats.average /= pct;
	return true;
}

bool MeshGenerator2d::autoTriangulate(MeshContainer3d* boundary, 
			int sm_count, bool start_clean, bool use_cs3d)
{
	if(start_clean) boundary->clearDiscretization();

	int pcount = boundary->getPointsCount();
	if(pcount < 1) return false;
	int loop_counter = 0;
	int invalid_count = 0;

	int fleft = scatterFreePoints(boundary);
	if(fleft > 0)
		LOG4CPLUS_WARN(MeshLog::logger_console, fleft << " free points left unassigned.");

	int bct = boundary->getBlocksCount();
	DataVector<bool> cs3d_changed(bct, false);

	do{
		++loop_counter;
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Auto-control-adjustment [2D]: loop #" << loop_counter);
		// GEN-1D
		for(IteratorEdge3d it = boundary->getFirstEdge3d(); it.isValid(); it.nextEdge()){
			MeshDomainEdge3d* domain_edge = (MeshDomainEdge3d*)it.getEdge();
			assert(domain_edge->getType() == EDGE_DOMAIN_3D);
			if(!domain_edge->isValidDiscretization()){
				int face_count = domain_edge->getFaceCount();
				for(int k = 0; k < face_count; k++)
					((MeshDomainSurface*)domain_edge->getFaceAt(k))->clearDiscretization();
				MeshGenerator1d::discretizeEdgeMin(domain_edge);

				if(use_cs3d){
					DataVector<MeshBlock*> blocks;
					if(domain_edge->adjacentBlocks(blocks, false)){
						for(int i = 0; i < blocks.countInt(); i++){
							MeshDomainVolume* mdv = (MeshDomainVolume*) blocks[i];
							//if(mdv->getType() != VOLUME_DOMAIN) continue;
							auto cs3d = mdv->getControlSpace();
							if (!cs3d) continue;
							auto cs3da = cs3d->getAsAdaptive();
							if(cs3da) 
								if(!cs3da->isCompacted() && domain_edge->updateACS(cs3d))
									cs3d_changed[i] = true;
						}
					}
				}
			}
		}
		if(use_cs3d){
			for(int i = 0; i < bct; i++){
				MeshDomainVolume* mdv = (MeshDomainVolume*)boundary->getBlockAt(i);
				//if(mdv->getType() != VOLUME_DOMAIN) continue;
				auto cs3d = mdv->getControlSpace();
				if (!cs3d) continue;
				auto cs3da = cs3d->getAsAdaptive();
				if (cs3da){
					if(cs3d_changed[i] && cs3da->smoothen()){
						int fct = mdv->getFaceCount();
						for(int j = 0; j < fct; j++){
							MeshDomainSurface* mds = (MeshDomainSurface*) mdv->getFace(j);
							//if(mds->getType() != FACE_DOMAIN) continue;
							auto cs2d = mds->getBoundary()->getControlSpace();
							if (!cs2d) continue;
							auto cs2da = cs2d->getAsAdaptive();
							if (cs2da && cs2da->applyAsMinimum(cs3d))
								mds->clearDiscretization();
						}
						cs3d_changed[i] = false;
					}
				}
			}
		}
		invalid_count = 0;
		// GEN-2D
		for(IteratorBoundary2d it = boundary->getFirstValidBoundary2d(); it.isValid(); it.nextValidBoundary()){
			MeshDomainSurface* domain_surface = it.getDomainSurface();
			if(domain_surface->getMesh() == nullptr || 
				domain_surface->getMesh()->getDiscretizationState() == 0)
			{
				domain_surface->createBoundaryMesh();
				CS2dPtr control_space = domain_surface->getBoundary()->getControlSpace();
				Metric2dContext mc(control_space);
				if(control_space->isAdaptive() && loop_counter < param_max_auto_retriangulation_count){
					if(!MeshGenerator2d::checkControlAtBoundary(mc, domain_surface->getBoundaryMesh())) {
						++invalid_count;
						LOG4CPLUS_DEBUG(MeshLog::logger_console, "Re-discretization of boundary contours required.");
					}else{
						domain_surface->triangulateBoundary(mc);
						if(loop_counter < 4 && !MeshGenerator2d::checkControlForCloseBoundaryEdges(mc, domain_surface->getMesh())){
							++invalid_count;
							LOG4CPLUS_DEBUG(MeshLog::logger_console, "Re-discretization of boundary contours required.");
						}
					}
				}else{
					domain_surface->triangulateBoundary(mc);
				}

				// mark adjacent blocks as invalid
				MeshDomainVolume* volume = (MeshDomainVolume*)domain_surface->getBlock(0);
				if(volume && volume->getMesh()){
					volume->getMesh()->setDiscretizationState(0);
					//if(use_cs3d && control_space->isAdaptive()){
					//	// 2d surface meshes -> control space 3d
					//	CS3dPtr cs3d = (CS3dPtr) volume->getControlSpace();
					//	if(cs3d && cs3d->isAdaptive()){
					//		if(cs3d->applyAsMinimum(control_space))
					//			cs3d_changed[volume->getIndex()] = true;
					//	}
					//}
				}
				volume = (MeshDomainVolume*)domain_surface->getBlock(1);
				if(volume && volume->getMesh()){
					volume->getMesh()->setDiscretizationState(0);
					//if(use_cs3d && control_space->isAdaptive()){
					//	// 2d surface meshes -> control space 3d
					//	CS3dPtr cs3d = (CS3dPtr) volume->getControlSpace();
					//	if(cs3d && cs3d->isAdaptive()){
					//		if(cs3d->applyAsMinimum(control_space))
					//			cs3d_changed[volume->getIndex()] = true;
					//	}
					//}
				}
			}
		}
		//if(use_cs3d){
		//	for(int i = 0; i < bct; i++){
		//		MeshDomainVolume* mdv = (MeshDomainVolume*)boundary->getBlockAt(i);
		//		//if(mdv->getType() != VOLUME_DOMAIN) continue;
		//		CS3dPtr cs3d = (CS3dPtr) mdv->getControlSpace();
		//		if(cs3d && cs3d->isAdaptive()){
		//			if(cs3d_changed[i] && cs3d->smoothen()){
		//				int fct = mdv->getFaceCount();
		//				for(int j = 0; j < fct; j++){
		//					MeshDomainSurface* mds = (MeshDomainSurface*) mdv->getFace(j);
		//					//if(mds->getType() != FACE_DOMAIN) continue;
		//					CS2dPtr cs = (CS2dPtr)mds->getBoundary()->getControlSpace();
		//					if(!cs->isAdaptive()) continue;
		//					if(cs->applyAsMinimum(cs3d)){
		//						mds->clearDiscretization();
		//						++invalid_count;
		//					}
		//				}
		//				cs3d_changed[i] = false;
		//			}
		//		}
		//	}
		//}

	}while(invalid_count > 0);

	for(IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
		MeshContainer2d* mesh = it.getMesh();
		if(param_mesh_decomposition == 1 || param_mesh_decomposition == 2){ // wd decomposition
			MeshGenerator2d::createDecompositionMesh(mesh);
			START_CLOCK("DCMP mesh split");
			MeshContainer2d* second_mesh = mesh->splitByElements(2);

			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NP: " << mesh->getPointsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NT: " << mesh->getElementsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NP: " << second_mesh->getPointsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NT: " << second_mesh->getElementsCount());

			STOP_CLOCK("DCMP mesh split");
			SHOW_MESH("decomposition mesh [0]", mesh->getViewSet());
			SHOW_MESH("decomposition mesh [1]", second_mesh->getViewSet());
			// continue meshing
			Metric2dContext mc(mesh->getControlSpace());
			START_CLOCK("DCMP mesh[0] refinement");
			MeshGenerator2d::addInnerNodes(mc, mesh);
			MeshGenerator2d::smoothen(mc, mesh, sm_count);
			STOP_CLOCK("DCMP mesh[0] refinement");
			START_CLOCK("DCMP mesh[1] refinement");
			MeshGenerator2d::addInnerNodes(mc, second_mesh);
			MeshGenerator2d::smoothen(mc, second_mesh, sm_count);
			STOP_CLOCK("DCMP mesh[1] refinement");

			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NP: " << mesh->getPointsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NT: " << mesh->getElementsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NP: " << second_mesh->getPointsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NT: " << second_mesh->getElementsCount());

			MeshViewSet* vset = mesh->getViewSet();
			vset = second_mesh->getViewSet(vset);
			SHOW_MESH("decomposition mesh [0+1]", vset);
			delete second_mesh;
		}else{
			auto cs = mesh->getControlSpace();
			//cs->showMetricLength();
			Metric2dContext mc(cs);
			MeshGenerator2d::addInnerNodes(mc, mesh);
			MeshGenerator2d::smoothen(mc, mesh, sm_count);
			//SHOW_MESH("mesh 2d-smoothen", mesh->getViewSet());
		}
		mesh->setDiscretizationState(2);
	}

	return true;
}

bool MeshGenerator2d::triangulateWithCS3d(MeshContainer3d* boundary, int sm_count)
{
	int pcount = boundary->getPointsCount();
	if (pcount < 1) return false;

	int fleft = scatterFreePoints(boundary);
	if (fleft > 0)
		LOG4CPLUS_WARN(MeshLog::logger_console, fleft << " free points left unassigned");

	int bct = boundary->getBlocksCount();
	for (int i = 0; i < bct; i++) {
		auto block = boundary->getBlockAt(i);
		assert(block->getType() == BLOCK_DOMAIN);
		auto cs3d = ((MeshDomainVolume*)block)->getControlSpace();
		assert(cs3d != nullptr);
		int bfct = block->getFaceCount();
		for (int j = 0; j < bfct; j++) {
			MeshDomainSurface* dface = ((MeshDomainSurface*)block->getFace(j));
			assert(dface->getType() == FACE_DOMAIN);
			auto face_boundary = dface->getBoundary();
			if (!face_boundary->getControlSpace()) {
				face_boundary->setControlSpace(
					std::make_shared<ControlSpace2dProjected>(cs3d, dface->getBaseSurface())
				);
			}
		}
	}

	// GEN-1D
	auto boundary_edges_3d = boundary->getEdges3d();
	for (auto edge : boundary_edges_3d) {
		MeshDomainEdge3d* domain_edge = (MeshDomainEdge3d*)edge;
		assert(domain_edge->getType() == EDGE_DOMAIN_3D);
		if (!domain_edge->isValidDiscretization()) {
			int face_count = domain_edge->getFaceCount();
			for (int k = 0; k < face_count; k++)
				((MeshDomainSurface*)domain_edge->getFaceAt(k))->clearDiscretization();
			MeshGenerator1d::discretizeEdgeMin(domain_edge);
		}
	}

	// GEN-2D
	for (auto domain_surface : boundary->getValidBoundaries2d()) {
		if (domain_surface->getMesh() == nullptr ||
			domain_surface->getMesh()->getDiscretizationState() == 0)
		{
			domain_surface->createBoundaryMesh();
			CS2dPtr control_space = domain_surface->getBoundary()->getControlSpace();
			Metric2dContext mc(control_space);
			domain_surface->triangulateBoundary(mc);
		}

		// mark adjacent blocks as invalid
		MeshDomainVolume* volume = (MeshDomainVolume*)domain_surface->getBlock(0);
		if (volume && volume->getMesh()) {
			volume->getMesh()->setDiscretizationState(0);
		}
		volume = (MeshDomainVolume*)domain_surface->getBlock(1);
		if (volume && volume->getMesh()) {
			volume->getMesh()->setDiscretizationState(0);
		}
	}

	for (auto mesh : boundary->getValidMeshes2d()) {
		auto cs = mesh->getControlSpace();
		//cs->showMetricLength();
		Metric2dContext mc(cs);
		MeshGenerator2d::addInnerNodes(mc, mesh);
		MeshGenerator2d::smoothen(mc, mesh, sm_count);
		//SHOW_MESH("mesh 2d-smoothen", mesh->getViewSet());
		mesh->setDiscretizationState(2);
	}

	return true;
}

bool MeshGenerator2d::autoTriangulateRemesh(MeshContainer3d* boundary,
	int sm_count, bool start_clean, bool use_cs3d)
{
	if (start_clean) boundary->clearDiscretization();

	int pcount = boundary->getPointsCount();
	if (pcount < 1) return false;
	int loop_counter = 0;
	int invalid_count = 0;

	int fleft = scatterFreePoints(boundary);
	if (fleft > 0)
		LOG4CPLUS_WARN(MeshLog::logger_console, fleft << " free points left unassigned");

	int bct = boundary->getBlocksCount();
	DataVector<bool> cs3d_changed(bct, false);

	do {
		++loop_counter;
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Auto-control-adjustment [2D], loop #" << loop_counter);
		// GEN-1D
		for (IteratorEdge3d it = boundary->getFirstEdge3d(); it.isValid(); it.nextEdge()) {
			MeshDomainEdge3d* domain_edge = (MeshDomainEdge3d*)it.getEdge();
			assert(domain_edge->getType() == EDGE_DOMAIN_3D);
			if (!domain_edge->isValidDiscretization()) {
				int face_count = domain_edge->getFaceCount();
				for (int k = 0; k < face_count; k++)
					((MeshDomainSurface*)domain_edge->getFaceAt(k))->clearDiscretization();
				MeshGenerator1d::discretizeEdgeMin(domain_edge);

				if (use_cs3d) {
					DataVector<MeshBlock*> blocks;
					if (domain_edge->adjacentBlocks(blocks, false)) {
						for (int i = 0; i < blocks.countInt(); i++) {
							MeshDomainVolume* mdv = (MeshDomainVolume*)blocks[i];
							//if(mdv->getType() != VOLUME_DOMAIN) continue;
							auto cs3d = mdv->getControlSpace();
							if (!cs3d) continue;
							auto cs3da = cs3d->getAsAdaptive();
							if (cs3da)
								if (!cs3da->isCompacted() && domain_edge->updateACS(cs3d))
									cs3d_changed[i] = true;
						}
					}
				}
			}
		}
		if (use_cs3d) {
			for (int i = 0; i < bct; i++) {
				MeshDomainVolume* mdv = (MeshDomainVolume*)boundary->getBlockAt(i);
				//if(mdv->getType() != VOLUME_DOMAIN) continue;
				auto cs3d = mdv->getControlSpace();
				if (!cs3d) continue;
				auto cs3da = cs3d->getAsAdaptive();
				if (cs3da) {
					if (cs3d_changed[i] && cs3da->smoothen()) {
						int fct = mdv->getFaceCount();
						for (int j = 0; j < fct; j++) {
							MeshDomainSurface* mds = (MeshDomainSurface*)mdv->getFace(j);
							//if(mds->getType() != FACE_DOMAIN) continue;
							auto cs2d = mds->getBoundary()->getControlSpace();
							if (!cs2d) continue;
							auto cs2da = cs2d->getAsAdaptive();
							if (cs2da && cs2da->applyAsMinimum(cs3d))
								mds->clearDiscretization();
						}
						cs3d_changed[i] = false;
					}
				}
			}
		}
		invalid_count = 0;
		// GEN-2D
		for (IteratorBoundary2d it = boundary->getFirstValidBoundary2d(); it.isValid(); it.nextValidBoundary()) {
			MeshDomainSurface* domain_surface = it.getDomainSurface();
			if (domain_surface->getMesh() == nullptr ||
				domain_surface->getMesh()->getDiscretizationState() == 0)
			{
				domain_surface->createBoundaryMesh();
				CS2dPtr control_space = domain_surface->getBoundary()->getControlSpace();
				Metric2dContext mc(control_space);
				if (control_space->isAdaptive() && loop_counter < param_max_auto_retriangulation_count) {
					if (!MeshGenerator2d::checkControlAtBoundary(mc, domain_surface->getBoundaryMesh())) {
						++invalid_count;
						LOG4CPLUS_DEBUG(MeshLog::logger_console, "Re-discretization of boundary contours required.");
					}
					else {
						domain_surface->triangulateBoundary(mc);
						if (loop_counter < 4 && !MeshGenerator2d::checkControlForCloseBoundaryEdges(mc, domain_surface->getMesh())) {
							++invalid_count;
							LOG4CPLUS_DEBUG(MeshLog::logger_console, "Re-discretization of boundary contours required.");
						}
					}
				}
				else {
					domain_surface->triangulateBoundary(mc);
				}

				// mark adjacent blocks as invalid
				MeshDomainVolume* volume = (MeshDomainVolume*)domain_surface->getBlock(0);
				if (volume && volume->getMesh()) {
					volume->getMesh()->setDiscretizationState(0);
					//if(use_cs3d && control_space->isAdaptive()){
					//	// 2d surface meshes -> control space 3d
					//	CS3dPtr cs3d = (CS3dPtr) volume->getControlSpace();
					//	if(cs3d && cs3d->isAdaptive()){
					//		if(cs3d->applyAsMinimum(control_space))
					//			cs3d_changed[volume->getIndex()] = true;
					//	}
					//}
				}
				volume = (MeshDomainVolume*)domain_surface->getBlock(1);
				if (volume && volume->getMesh()) {
					volume->getMesh()->setDiscretizationState(0);
					//if(use_cs3d && control_space->isAdaptive()){
					//	// 2d surface meshes -> control space 3d
					//	CS3dPtr cs3d = (CS3dPtr) volume->getControlSpace();
					//	if(cs3d && cs3d->isAdaptive()){
					//		if(cs3d->applyAsMinimum(control_space))
					//			cs3d_changed[volume->getIndex()] = true;
					//	}
					//}
				}
			}
		}
		//if(use_cs3d){
		//	for(int i = 0; i < bct; i++){
		//		MeshDomainVolume* mdv = (MeshDomainVolume*)boundary->getBlockAt(i);
		//		//if(mdv->getType() != VOLUME_DOMAIN) continue;
		//		CS3dPtr cs3d = (CS3dPtr) mdv->getControlSpace();
		//		if(cs3d && cs3d->isAdaptive()){
		//			if(cs3d_changed[i] && cs3d->smoothen()){
		//				int fct = mdv->getFaceCount();
		//				for(int j = 0; j < fct; j++){
		//					MeshDomainSurface* mds = (MeshDomainSurface*) mdv->getFace(j);
		//					//if(mds->getType() != FACE_DOMAIN) continue;
		//					CS2dPtr cs = (CS2dPtr)mds->getBoundary()->getControlSpace();
		//					if(!cs->isAdaptive()) continue;
		//					if(cs->applyAsMinimum(cs3d)){
		//						mds->clearDiscretization();
		//						++invalid_count;
		//					}
		//				}
		//				cs3d_changed[i] = false;
		//			}
		//		}
		//	}
		//}

	} while (invalid_count > 0);

	for (IteratorMesh2d it = boundary->getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()) {
		MeshContainer2d* mesh = it.getMesh();
		if (param_mesh_decomposition == 1 || param_mesh_decomposition == 2) { // wd decomposition
			MeshGenerator2d::createDecompositionMesh(mesh);
			START_CLOCK("DCMP mesh split");
			MeshContainer2d* second_mesh = mesh->splitByElements(2);

			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NP: " << mesh->getPointsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NT: " << mesh->getElementsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NP: " << second_mesh->getPointsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NT: " << second_mesh->getElementsCount());

			STOP_CLOCK("DCMP mesh split");
			SHOW_MESH("decomposition mesh [0]", mesh->getViewSet());
			SHOW_MESH("decomposition mesh [1]", second_mesh->getViewSet());
			// continue meshing
			Metric2dContext mc(mesh->getControlSpace());
			START_CLOCK("DCMP mesh[0] refinement");
			MeshGenerator2d::addInnerNodes(mc, mesh);
			MeshGenerator2d::smoothen(mc, mesh, sm_count);
			STOP_CLOCK("DCMP mesh[0] refinement");
			START_CLOCK("DCMP mesh[1] refinement");
			MeshGenerator2d::addInnerNodes(mc, second_mesh);
			MeshGenerator2d::smoothen(mc, second_mesh, sm_count);
			STOP_CLOCK("DCMP mesh[1] refinement");

			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NP: " << mesh->getPointsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [0] NT: " << mesh->getElementsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NP: " << second_mesh->getPointsCount());
			LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP [1] NT: " << second_mesh->getElementsCount());

			MeshViewSet* vset = mesh->getViewSet();
			vset = second_mesh->getViewSet(vset);
			SHOW_MESH("decomposition mesh [0+1]", vset);
			delete second_mesh;
		}
		else {
			Metric2dContext mc(mesh->getControlSpace());
			MeshGenerator2d::addInnerNodes(mc, mesh);
			MeshGenerator2d::smoothen(mc, mesh, sm_count);
		}
		mesh->setDiscretizationState(2);
	}

	return true;
}

bool MeshGenerator2d::statMetricQuality(Metric2dContext& mc, MeshContainer2d* mesh, DataStatistics & stats, int criterion)
{
	int ect = mesh->getElementsCount();
	if(ect < 1) return false;

	for(int i = 0; i < ect; i++){
		MeshTriangle2d* triangle = (MeshTriangle2d*)mesh->getElementAt(i);
		if(triangle->getEdgeCount() != 3) continue;
		stats.add(triangle->countQuality(mc, false, criterion));
	}

	return true;
}

/// Checks if control space is consistent with boundary edges (and adjusts if necessary)
bool MeshGenerator2d::checkControlAtBoundary(Metric2dContext& mc, MeshContainer2d* mesh)
{
	bool proper = true;
	assert(mesh);
	auto space = mesh->getControlSpace();
	auto acs = space->getAsAdaptive();
	if(!space) return true; // no change possible anyway

//	unsigned int boundary_segment_counter = 0;
//	unsigned int boundary_shape_counter = 0;
//	unsigned int boundary_triangle_counter = 0;

	double cr2 = sqr(ControlSpace2dAdaptive::param_curvature_ratio);
	double d2_threshold = 20*sqr(1-sqrt(1-cr2/4))/cr2;

//	int eps_counter = 0;
//	space->storeSizingEPS(mesh, "control-start");
	int pct = mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint2d* mpt = mesh->getPointAt(i);
		if(!mpt->isBorder()) continue;
		int rank = mpt->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = mpt->getEdge(j);
			if(!edge->isBorder()) continue;
			assert(!edge->availableTag(TagExtended::TAG_BOUNDARY_EDGE)); // i.e. this is "boundary-mesh" (only border edges, no triangles)
			MeshDomainEdge3d* edge3d = (MeshDomainEdge3d*)edge->getPtrTag(TagExtended::TAG_ME_2D_3D);
			assert(edge3d->getType() == EDGE_DOMAIN_3D);
			// check length (!!! sets metric at the middle of the edge)
			double len = edge->getLengthMetricAdapted(mc);
//			LOG4CPLUS_INFO(MeshLog::logger_mesh, "len = " << len);
			if(len < 0.8){
				acs->updateForBoundarySegment(edge);
			}
//			if(len < 0.5 || len > 2.5){
			if(len < 0.5 || len > 1.5){
//				LOG4CPLUS_DEBUG(MeshLog::logger_console, "CS-update, b-edge metric length", len);
//				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " *** boundary-segment " << edge->getMeshPoint(0)->getID() <<
//					"-" << edge->getMeshPoint(1)->getID() << " len=" << len << endl;
				edge3d->setValidDiscretization(proper=false);
				//++boundary_segment_counter;
			}
			// check shape
			if(edge->updateACSwithCurvature(space, mesh->getSurface(), d2_threshold)){
				edge3d->setValidDiscretization(proper=false);
				//++boundary_shape_counter;
			}
		}
	}

	if(ControlSpace2dAdaptive::param_gradation_ratio > 0.0){
		bool any_changes = acs->smoothen();
		//--
		if(any_changes){ // check segments ones more
//			space->storeSizingEPS(mesh, "control-smoothing");
			for(int i = 0; i < pct; i++){
				MeshPoint2d* mpt = mesh->getPointAt(i);
				if(!mpt->isBorder()) continue;
				int rank = mpt->getRank();
				for(int j = 0; j < rank; j++){
					MeshEdge2d* edge = mpt->getEdge(j);
					if(!edge->isBorder()) continue;
					MeshDomainEdge3d* edge3d = (MeshDomainEdge3d*)edge->getPtrTag(TagExtended::TAG_ME_2D_3D);
					// check length (!!! sets metric at the middle of the edge)
					double len = edge->getLengthMetricAdapted(mc);
		//			double len = edge->getLength(true);
//					LOG4CPLUS_INFO(MeshLog::logger_mesh, "* len = " << len);
					if(len > 1.5){
						LOG4CPLUS_DEBUG(MeshLog::logger_console, 
							"CS-update, b-edge metric length (after smoothing): " << len);
//						LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " *** boundary-segment " << edge->getMeshPoint(0)->getIndex() <<
//							"-" << edge->getMeshPoint(1)->getIndex() << " len=" << len << endl;
						edge3d->setValidDiscretization(proper=false);
						//++boundary_segment_counter;
					}
				}
			}
			//space->smoothen();
		}
		//space->logDescription();
//		SHOW_MESH("boundary-update-triangles", mesh->getViewSet());
	}

//	space->logDescription();

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "* boundary-segment : " << boundary_segment_counter);
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "* boundary-shape   : " << boundary_shape_counter);
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "* boundary-triangle: " << boundary_triangle_counter);

	return proper;
}

bool MeshGenerator2d::checkControlForCloseBoundaryEdges(Metric2dContext& mc, MeshContainer2d* mesh)
{
	auto space = mesh->getControlSpace()->getAsAdaptive();
	if(!space) return true; // no change possible

	bool proper = true;
	bool size_changes = false;

	int pct = mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint2d* mpt = mesh->getPointAt(i);
		if(!mpt->isBorder()) continue;
		int rank = mpt->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = mpt->getEdge(j);
			if(!edge->isBorder()){
				if(edge->getOtherPoint(mpt)->isBorder() && edge->getPointIndex(mpt) == 0){
					double len = edge->getLengthMetricAdapted(mc);
					if(len < 0.8){
						size_changes |= space->updateForBoundarySegment(edge);
					}
				}else continue;
			}
			MeshDomainEdge3d* edge3d = (MeshDomainEdge3d*)edge->getPtrTag(TagExtended::TAG_ME_2D_3D);
			// check length (!!! sets metric at the middle of the edge)
			double len = edge->getLengthMetricAdapted(mc);
			// check incident elements

			if(abs(len - 1.0) < 0.2){
//			if(false){
				MeshTriangle2d* triangle = (MeshTriangle2d*)edge->getMeshElement(mpt);
				if(!triangle || triangle->getType() != ELEMENT_MESH_TRIANGLE)
					continue;
				MeshPoint2d* mpt1 = edge->getOtherPoint(mpt);
				assert(mpt1->isBorder());
				MeshPoint2d* mpt2 = triangle->getOtherPoint(mpt, mpt1);
				assert(mpt2->isBorder());
				if(triangle->getPrevEdge(edge)->isBorder()) continue;
				if(triangle->getNextEdge(edge)->isBorder()) continue;
				double area = DTriangle2d::area(
					mpt->getMetricCoordinates(mc), 
					mpt1->getMetricCoordinates(mc), 
					mpt2->getMetricCoordinates(mc));
				assert(area > 0.0);
				if(area < 0.2){
					//MeshView::showDebugMesh("triangle to boundary improve check", mesh, triangle);
					if(space->updateForBoundaryTriangle(edge, mpt2)){

//						LOG4CPLUS_INFO(MeshLog::logger_mesh, "Area= " << area);
//						double len_x = mpt->getMetricCoordinates().distance(mpt1->getMetricCoordinates());
//						LOG4CPLUS_INFO(MeshLog::logger_mesh, "Edge= " << len_x);
//						MeshView::showDebugMesh("triangle to boundary improve", mesh, triangle);
//						triangle->setTag();

						//++boundary_triangle_counter;
						LOG4CPLUS_DEBUG(MeshLog::logger_console, "CS-update, b-triangle metric area: " << area);
						if(edge3d) edge3d->setValidDiscretization(proper=false);
					}
				}
			}

		}
	}

	if(ControlSpace2dAdaptive::param_gradation_ratio > 0.0 && (!proper || size_changes)){
		bool any_changes = space->smoothen();
		//--
		if(any_changes){ // check segments ones more
//			space->storeSizingEPS(mesh, "control-smoothing");
			for(int i = 0; i < pct; i++){
				MeshPoint2d* mpt = mesh->getPointAt(i);
				if(!mpt->isBorder()) continue;
				int rank = mpt->getRank();
				for(int j = 0; j < rank; j++){
					MeshEdge2d* edge = mpt->getEdge(j);
					if(!edge->isBorder()) continue;
					MeshDomainEdge3d* edge3d = (MeshDomainEdge3d*)edge->getPtrTag(TagExtended::TAG_ME_2D_3D);
					// check length (!!! sets metric at the middle of the edge)
					double len = edge->getLengthMetricAdapted(mc);
		//			double len = edge->getLength(true);
					if(len > 2.5){
//						LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " *** boundary-segment " << edge->getMeshPoint(0)->getIndex() <<
//							"-" << edge->getMeshPoint(1)->getIndex() << " len=" << len << endl;
						edge3d->setValidDiscretization(proper=false);
						//++boundary_segment_counter;
					}
				}
			}
			//space->smoothen();
		}
		//space->logDescription();
//		SHOW_MESH("boundary-update-triangles", mesh->getViewSet());
	}

	return proper;
}

CS2dAPtr MeshGenerator2d::createNewControlSpace(SurfaceConstPtr surface, const DRect& rect)
{
	// select control type
	switch(ControlSpace2d::param_control_type){
	case MeshData::CONTROL_UNIFORM:
		return std::make_shared<ControlSpace2dMatrixUniform>(surface, rect,
			ControlSpace2dMatrixUniform::param_uniform_nx, 
			ControlSpace2dMatrixUniform::param_uniform_nx);
	case MeshData::CONTROL_MESH:
		return std::make_shared<ControlSpace2dMesh>(surface, rect,
			ControlSpace2dMesh::param_interpolation_method);
	case MeshData::CONTROL_QUADTREE:
		return std::make_shared<ControlSpace2dQuadTree>(surface, rect);
	case MeshData::CONTROL_KDTREE_L:
		return std::make_shared<ControlSpace2dKdTreeL>(surface, rect);
	case MeshData::CONTROL_KDTREE_V:
		return std::make_shared<ControlSpace2dKdTreeV>(surface, rect);
	}
	// if none of the above
	assert(false);
	return nullptr;
}

bool MeshGenerator2d::smoothen(Metric2dContext& mc, MeshContainer2d* mesh, int steps, 
		TagExtended::TagType tag_type, int tag_value, int method)
{
	for(int i = 0; i < steps; i++){
		if((method & MeshData::SM_TOP_SWAP) != 0){
			MeshGenerator2d::smoothenTopologicalSwap(mesh, tag_type, tag_value);
		}
		if((method & MeshData::SM_LAPLACE) != 0){
			MeshGenerator2d::smoothenLaplace(mc, mesh, false, tag_type, tag_value); 
		}
		if((method & MeshData::SM_LAPLACE_MIXED) != 0){
			MeshGenerator2d::smoothenLaplaceMixed(mc, mesh, tag_type, tag_value); 
		}
		if((method & MeshData::SM_LAPLACE_METRIC) != 0){
			MeshGenerator2d::smoothenLaplace(mc, mesh, true, tag_type, tag_value); 
		}
		if((method & MeshData::SM_METRIC) != 0){
			MeshGenerator2d::smoothenMetric(mc, mesh, tag_type, tag_value); 
		}
		if((method & MeshData::SM_DEL_SWAP) != 0){
			MeshGenerator2d::smoothenDelaunaySwap(mc, mesh, tag_type, tag_value); 
		}
		if((method & MeshData::SM_DEL_SWAP_COMPLETE) != 0){
			MeshGenerator2d::smoothenDelaunaySwapComplete(mc, mesh, tag_type, tag_value);
		}
	}
	MeshGenerator2d::smoothenPostCheck(mc, mesh);

	return true;
}

bool MeshGenerator2d::createDecompositionMesh(MeshContainer2d* mesh)
{
	// Calculate quality of present triangles + inertia center + some axis
	int count = mesh->getElementsCount();
	if(count < 1) return 0;

//	LOG4CPLUS_INFO(MeshLog::logger_console, "Decomposition mesh (NT)", count);
//	SHOW_MESH("boundary based mesh", mesh->getViewSet());


	CS2dPtr main_space = mesh->getControlSpace();
	assert(main_space);
	Metric2dContext mc(main_space);

	double total_mass = 0;
	DPoint2d inertial_center;

	START_CLOCK("DCMP evaluation");

	for(int i = 0; i < count; i++)
		mesh->getElementAt(i)->addForInertialCenter(mc, inertial_center, total_mass);
	inertial_center /= total_mass;

	int ev_count = (int) (total_mass * (4.0/SQRT3));
	LOG4CPLUS_INFO(MeshLog::logger_console, "DCMP evaluation NT count: " << ev_count);

	STOP_CLOCK("DCMP evaluation");

	DMatrix2d inertial_moments;
	for(int i = 0; i < count; i++)
		mesh->getElementAt(i)->addForInertialMoments(mc, inertial_center, inertial_moments);

	DMatrix2d v; // eigenvectors
	double d[2]; // eigenvalues;
	inertial_moments.eigensystem(v, d);

//	cout << "** vector-0 = " << v.column(0) << endl;
//	cout << "** value-0  = " << d[0] << endl;
//	cout << "** vector-1 = " << v.column(1) << endl;
//	cout << "** value-1  = " << d[1] << endl;

	// selected line -> orthogonal to smallest
	//		for 2d -> it means simply larger one
	const DVector2d line_vt = (abs(d[0])>abs(d[1])) ? v.column(0) : v.column(1);
/*
	double ev_nt_count = 0.0;
	for(int i = 0; i < count; i++){
		MeshElement* element = mesh->getElementAt(i);
		element->countQuality(mc, true);
		ev_nt_count += MeshElement::last_area;
		inertial_center += element->getMiddlePoint() * MeshElement::last_area;
	}
	
	// normalize inertia center (divide by total mass)
	inertia_centerl /= ev_nt_count;
	// normalize NT evaluation (total metric area -> number of triangles)
	ev_nt_count *= 4.0/SQRT3;
*/

	const DRect rect = mesh->getBoundingRect();
	// calculate center and axis of inertia
	DPoint2d dcmp_pt0, dcmp_pt1;
	bool cross_ok = rect.calculateCrossingPoints(inertial_center, line_vt, dcmp_pt0, dcmp_pt1);
	if(!cross_ok) LOG4CPLUS_WARN(MeshLog::logger_console, "Something wrong with cross-rect-line.");
	// ... simple y=x+b line, crossing the inertia center
	//const DPoint2d dcmp_pt0(rect.left,  inertial_center.y);   // OX
	//const DPoint2d dcmp_pt1(rect.right, inertial_center.y);
	//const DPoint2d dcmp_pt0(inertial_center.x, rect.bottom);	 // OY
	//const DPoint2d dcmp_pt1(inertial_center.x, rect.top);
		
	CS2dAPtr decomposition_space;

	// create special control space
	if(param_mesh_decomposition == 1){
		// "1 - special control space as created copy"

		decomposition_space  = createNewControlSpace(mesh->getSurface(), rect);
		assert(decomposition_space );

		// set max
		decomposition_space->setMaxMetric(1.0);

		// set metric along the selected line
		const DVector2d dv = dcmp_pt1 - dcmp_pt0;
		const DVector2d dvn(dv.y, -dv.x); // orthogonal vector
		decomposition_space->setMinControlLine(mc, main_space, dcmp_pt0, dcmp_pt1, dvn, 
			param_mesh_decomposition_width);

		// adapt and smoothen
		decomposition_space->adaptToParameterization();
		if(param_mesh_decomposition_gradation > 0)
			decomposition_space->smoothenWithGradation(param_mesh_decomposition_gradation);

	}else if(param_mesh_decomposition == 2){
		// "2 - special control space as filtered layer"
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Not implemented yet.");
	}

	// switch spaces
	mesh->setControlSpace(decomposition_space);

	Metric2dContext dcmp_mc(decomposition_space);
	MeshGenerator2d::addInnerNodes(dcmp_mc, mesh);

	MeshGenerator2d::smoothenLaplaceMixed(dcmp_mc, mesh);
	MeshGenerator2d::smoothenTopologicalSwap(mesh);
	MeshGenerator2d::smoothenLaplaceMixed(dcmp_mc, mesh);

	// mark elements
	count = mesh->getElementsCount();
	for(int i = 0; i < count; i++){
		MeshElement* element = mesh->getElementAt(i);
		int det_neg = 0;
		int det_pos = 0;
		for(int j = 0; j < element->getEdgeCount(); j++){
			double det = DTriangle2d::det(dcmp_pt0, dcmp_pt1, element->getPoint(j)->getCoordinates());
			if(det >= 0) ++det_pos;
			else ++det_neg;
		}
//		if(det_neg == 0) element->setAreaID(2);
//		else if(det_pos == 0) element->setAreaID(3);
//		else element->setAreaID(4);
		element->setAreaID((det_pos > det_neg) ? 2 : 6);
	}

	// and switch space back
	// ... decomposition_space should be removed automatically
	mesh->setControlSpace(main_space); 

	return true;
}

/// Collapse too short edges
bool MeshGenerator2d::collapseEdges(Metric2dContext& mc, MeshContainer2d* mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	static const double MAX_LEN = 0.6;
	// gather edges
	DataContainer<ActiveEdge> active_edges(mesh->getPointsCount(), true);
	for(IteratorEdge2d it = mesh->getFirstEdge2d(); it.isValid(); it.nextEdge()){
		if(it.getEdge()->isBorder()) continue;
		MeshElement* element0 = it.getEdge()->getMeshElement(0);
		MeshElement* element1 = it.getEdge()->getMeshElement(1);
		if(element0 && !element0->checkIntTag(tag_type, tag_value)) continue;
		if(element1 && !element1->checkIntTag(tag_type, tag_value)) continue;

		double len = it.getEdge()->getLength(mc);
		if(len < MAX_LEN){
			ActiveEdge* act = new ActiveEdge(it.getEdge(), len);
			active_edges.addDataItem(act);
			it.getEdge()->setPtrTag(TagExtended::TAG_COLLAPSE_2D, act);
		}
	}

	mesh->setHeapOrder(false);

//	assert(mesh->isValid());

	int counter = 0;
	while(active_edges.countInt() > 0){
		++counter;
		ActiveEdge* act = active_edges.removeDataItem(0);
		act->edge->removeTag(TagExtended::TAG_COLLAPSE_2D);
		if(act->len > MAX_LEN){
			delete act;
			continue;
		}

		MeshPoint2d* p0 = act->edge->getMeshPoint(0);
		MeshPoint2d* p1 = act->edge->getMeshPoint(1);
		assert(!(p0->isBorder() && p1->isBorder()));

		MeshElement* e0 = act->edge->getMeshElement(0);
		MeshElement* e1 = act->edge->getMeshElement(1);
		MeshPoint2d* p2 = e0->getNextPoint(p1); assert(p2 != p0);
		MeshPoint2d* p3 = e1->getNextPoint(p0);	assert(p3 != p1);
		delete act;

		// check extra conditions
		bool valid = true;
		for(int i = 0; i < p0->getRank(); i++){
			MeshPoint2d* other_point = p0->getEdge(i)->getOtherPoint(p0);
			if(other_point != p1 && other_point != p2 && other_point != p3 &&
				other_point->getEdgeToPoint(p1) != nullptr)
			{
				valid = false;
				break;
			}
		}
		if(!valid) continue;
		for(int i = 0; i < p1->getRank(); i++){
			MeshPoint2d* other_point = p1->getEdge(i)->getOtherPoint(p1);
			if(other_point != p0 && other_point != p2 && other_point != p3 &&
				other_point->getEdgeToPoint(p0) != nullptr)
			{
				valid = false;
				break;
			}
		}
		if(!valid) continue;

		delete mesh->removeMeshElement(e0);
		delete mesh->removeMeshElement(e1);

		if(p0->isBorder()){
			// remove tag from edges which will be removed
			for(int i = 0; i < p1->getRank(); i++){
				MeshEdge2d* temp_edge = p1->getEdge(i);
				ActiveEdge* aact = (ActiveEdge*)temp_edge->getPtrTag(TagExtended::TAG_COLLAPSE_2D);
				if(aact){
					temp_edge->removeTag(TagExtended::TAG_COLLAPSE_2D);
					delete active_edges.removeDataItem(aact->getIndex());
				}
			}
			// switch all elements incident to point
			while(p1->getRank() > 0){
				MeshEdge2d* edge = p1->getEdge(0);
				MeshElement* temp_el0 = edge->getMeshElement(0);
				MeshElement* temp_el1 = edge->getMeshElement(1);
				if(temp_el0) temp_el0->switchPointsWithEdges(p1, p0);
				if(temp_el1) temp_el1->switchPointsWithEdges(p1, p0);
			}
			delete mesh->removeMeshPoint(p1);
		}else if(p1->isBorder()){
			// remove tag from edges which will be removed
			for(int i = 0; i < p0->getRank(); i++){
				MeshEdge2d* temp_edge = p0->getEdge(i);
				ActiveEdge* aact = (ActiveEdge*)temp_edge->getPtrTag(TagExtended::TAG_COLLAPSE_2D);
				if(aact){
					temp_edge->removeTag(TagExtended::TAG_COLLAPSE_2D);
					delete active_edges.removeDataItem(aact->getIndex());
				}
			}
			// switch all elements incident to point
			while(p0->getRank() > 0){
				MeshEdge2d* edge = p0->getEdge(0);
				MeshElement* temp_el0 = edge->getMeshElement(0);
				MeshElement* temp_el1 = edge->getMeshElement(1);
				if(temp_el0) temp_el0->switchPointsWithEdges(p0, p1);
				if(temp_el1) temp_el1->switchPointsWithEdges(p0, p1);
			}
			delete mesh->removeMeshPoint(p0);
			p0 = p1;
		}else{
			// remove tag from edges which will be removed
			for(int i = 0; i < p1->getRank(); i++){
				MeshEdge2d* temp_edge = p1->getEdge(i);
				ActiveEdge* aact = (ActiveEdge*)temp_edge->getPtrTag(TagExtended::TAG_COLLAPSE_2D);
				if(aact){
					temp_edge->removeTag(TagExtended::TAG_COLLAPSE_2D);
					delete active_edges.removeDataItem(aact->getIndex());
				}
			}
			// switch all elements incident to point
			while(p1->getRank() > 0){
				MeshEdge2d* edge = p1->getEdge(0);
				MeshElement* temp_el0 = edge->getMeshElement(0);
				MeshElement* temp_el1 = edge->getMeshElement(1);
				if(temp_el0) temp_el0->switchPointsWithEdges(p1, p0);
				if(temp_el1) temp_el1->switchPointsWithEdges(p1, p0);
			}
			delete mesh->removeMeshPoint(p1);
			movePointByLaplaceForVariableMetric(mc, p0);
		}
		// update active edges incident to p0
		for(int i = 0; i < p0->getRank(); i++){
			MeshEdge2d* edge = p0->getEdge(i);
			ActiveEdge* aact = (ActiveEdge*)edge->getPtrTag(TagExtended::TAG_COLLAPSE_2D);
			if(aact){
				double len = edge->getLength(mc);
				if(len != aact->len){
					if(len >= MAX_LEN){
						edge->removeTag(TagExtended::TAG_COLLAPSE_2D);
						delete active_edges.removeDataItem(aact->index);
					}else{
						aact->len = len;
						active_edges.updateDataItemPosition(aact);
					}
				}
			}else{
				double len = edge->getLength(mc);
				if(len < MAX_LEN){
					aact = new ActiveEdge(edge, len);
					active_edges.addDataItem(aact);
					edge->setPtrTag(TagExtended::TAG_COLLAPSE_2D, aact);
				}
			}
		}

		// debug check
//		assert(mesh->isValid());
		//for(int i = 0; i < active_edges.countInt(); i++){
		//	ActiveEdge* act = active_edges.getDataAt(i);
		//	bool ok = act->edge->availableTag(TagExtended::TAG_COLLAPSE_2D);
		//}
	}
	return true;
}

/// Scatter to blocks/patches/contours, returns number of points remaining...
int MeshGenerator2d::scatterFreePoints(MeshContainer3d* boundary)
{
	auto fps = boundary->getFreePoints();
	if(!fps || fps->empty()) return 0;
	//...
	int default_bid = -1;
	int dbct = boundary->getBlocksCount();
	if(dbct == 1){
		default_bid = boundary->getBlockAt(0)->getAreaID();
	}

	// remove duplicates ...
	DBox bbox;
	for(int i = 0; i < fps->countInt(); i++)
		bbox.addPoint(fps->get(i)->getCoordinates());
	bbox.inflate(0.05);
	auto cs = std::make_shared<ControlSpace3dIdentity>(boundary->getBoundingBox().getDiameter());
	Metric3dContext mc(cs);
	OctTreeMeshPoints octree(bbox);
	int duplicate_count = 0;
	int fcount = 0;
	for(int i = 0; i < fps->countInt(); i++){
		auto mpt = fps->get(i);
		if(octree.anyMeshPointInProximity(mc, mpt->getCoordinates(), SMALL_NUMBER)){
			fps->removeAt(i--);
			++duplicate_count;
		}else{
			octree.insertMeshPointLink(mpt.get());
			++fcount;
		}
	}
	LOG4CPLUS_INFO(MeshLog::logger_console, "Unique freepoints: " << fcount);
	if(duplicate_count > 0)
		LOG4CPLUS_WARN(MeshLog::logger_console, "Duplicate (removed) freepoints: " << duplicate_count);

	// domain-blocks per id ...
	DataHashTableKeyValue<int, MeshDomainVolume*> dblocks(dbct*2, -1);
	DataHashTableKeyValue<int, std::shared_ptr<DataVector<std::shared_ptr<MeshPoint3d>>>> dmpts(dbct*2, -1);
	for(int i = 0; i < dbct; i++){
		MeshDomainVolume* mdv = (MeshDomainVolume*) boundary->getBlockAt(i);
		int id = mdv->getAreaID();
		if(id > -1){
			dblocks.insert(id, mdv);
			dmpts.insert(id, 
				std::make_shared<DataVector<std::shared_ptr<MeshPoint3d>>>(dbct));
		}
	}

	int left_count = 0;
	while(fps->countInt() > 0){
		auto fpoint = fps->removeLast();
		int bid = fpoint->getIntTag(TagExtended::TAG_ID, -1);
		if(bid < 0) bid = default_bid;
		MeshDomainVolume* mdv = nullptr;
		if(bid >= 0) mdv = dblocks.getValue(bid, nullptr);
		if(!mdv || !mdv->offerFreePoint(fpoint)){
			++left_count;
		}else{
			auto vec = dmpts.getValue(bid, nullptr);
			vec->add(std::make_shared<MeshPoint3d>(*fpoint));
		}
	}

	DataVector<int> mdv_keys(dbct);
	dblocks.getKeys(mdv_keys);
	for(int i = 0; i < mdv_keys.countInt(); i++){
		MeshDomainVolume* mdv = dblocks.getValue(mdv_keys[i], nullptr);
		auto vec = dmpts.getValue(mdv_keys[i], nullptr);
		if(!vec->empty()){
			MeshContainer3d* mesh = MeshGenerator3d::triangulatePoints(vec);
			auto ucs = mdv->getUserControlSpace();
			if (!ucs || !ucs->isAdaptive()) {
				CS3dAPtr new_ucs(new ControlSpace3dOctree(mesh->getBoundingBox()));
				new_ucs->setMaxMetric();
				if (ucs) new_ucs->applyAsMinimum(ucs);
				mdv->setUserControlSpace(ucs = new_ucs);
			}
			auto ucsa = ucs->getAsAdaptive();
			assert(ucsa);
			// use edges as ucs-sources
			bool any_change = false;
			for(IteratorEdge3d it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
				MeshEdge3d* edge = it.getEdge();
				if(edge->getMeshPoint(0)->availableTag(TagExtended::TAG_OUTERHULL_POINT) ||
					edge->getMeshPoint(1)->availableTag(TagExtended::TAG_OUTERHULL_POINT))
					continue;
				any_change |= ucsa->updateForBoundarySegment(edge);
			}
			delete mesh;
			if(any_change) ucsa->smoothen();
		}
	}

	return left_count;
}

/// Create triangulation of polygon
MeshContainer2d* MeshGenerator2d::triangulatePoly( 
	const DataVector<DPoint2d> & poly, 
	SurfaceConstPtr surface)
{
	assert(surface);

	size_t pct = poly.countInt();
	if(pct < 3) return nullptr;

	assert( DPoint2d::properOrientation(poly) );
	
	// Create boundary mesh
	DataVector<MeshPoint2d*> mpoints2d(pct);
	MeshContainer2d* boundary = new MeshContainer2d((int)pct);

	DRect brect;
	for(size_t i = 0; i < pct; i++){
		brect.addPoint(poly[i]);
		MeshPoint2d* mpt = new MeshPoint2d(poly[i]);
		mpt->setBorder(TagBorder::OUTER | TagBorder::FIXED);
		mpt->setIntTag(TagExtended::TAG_ID, (int)i);
		mpoints2d.add(mpt);
		boundary->addMeshPoint(mpt);
	}
	for(size_t i = 0; i < pct; i++){
		MeshEdge2d* edge = new MeshEdge2d( boundary->getPointAt((int)i), boundary->getPointAt((int)((i+1)%pct)));
		edge->setBorder(TagBorder::OUTER | TagBorder::FIXED);
	}

	boundary->setSurface(surface);

	//if(true){
	//	MeshViewSet* set = new MeshViewSet();
	//	for(int i = 0; i < pct; i++){
	//		set->addPoint(DPoint3d(poly[i], 0.0), 0, i);
	//		set->addEdge(DPoint3d(poly[i], 0.0), DPoint3d(poly[(i+1)%pct], 0.0), 0);
	//	}
	//	SHOW_MESH("poly2d", set);
	//}

	MeshArea* area = new MeshArea(mpoints2d);
	area->setFilled();
	area->setAreaID(1);
	boundary->addMeshElement(area);

	MeshContainer2d* mesh = MeshGenerator2d::createInitialMesh(boundary, &brect);
	if(!mesh) return nullptr;

//	SHOW_MESH("initial mesh 2d", mesh->getViewSet());

	Metric2dContext mc(mesh->getControlSpace());
	MeshGenerator2d::addBoundaryNodes(mc, mesh, boundary);

//	SHOW_MESH("mesh with boundary nodes", mesh->getViewSet());

	if(!MeshGenerator2d::constrainToBorder(mc, boundary, mesh)) {
		delete boundary;
		delete mesh;
		return nullptr;
	}

//	SHOW_MESH("contrained", mesh->getViewSet());

	assert(pct == mesh->getPointsCount());
	for(int i = 0; i < pct; i++){
		MeshPoint2d* mp = mesh->getPointAt(i);
		MeshPoint2d* bp = (MeshPoint2d*) mp->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
		assert(bp);
		mp->setIntTag(TagExtended::TAG_ID, bp->getIntTag(TagExtended::TAG_ID, -1));
	}

	delete boundary;
	return mesh;
}


/// Create triangulation of planar 3d polygon
MeshContainer2d* MeshGenerator2d::triangulatePoly( const DataVector<DPoint3d> & poly )
{
	size_t pct = poly.countInt();
	if(pct < 3) return nullptr;

	DBox box;
	for(size_t i = 0; i < pct; i++)
		box.addPoint(poly[i]);


	// fit plane
	DPlane plane;
	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(poly, plane, true);
	const double FIT_ERR = 2e-2; 
	if(plane_max_dist < 0.0){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MG2d::triangulatePoly - error fitting plane");
		return nullptr;
	}else if(plane_max_dist > FIT_ERR * box.getDiameter()){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"MG2d::triangulatePoly - problem fitting plane, rel-dist=" << plane_max_dist/box.getDiameter());
	}

	// prepare 2d points
	DataVector<DPoint2d> poly2d(pct);
	for(size_t i = 0; i < pct; i++)
		poly2d.add( plane.projectToPlane( poly[i] ) );

	// create 2d triangulation
	return MeshGenerator2d::triangulatePoly(poly2d, std::make_shared<SurfacePlane>(plane) );

}
