// MeshGenerator2dQMorph.cpp: implementation of the MeshGenerator2dQMorph class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshGenerator2d.h"
#include "MeshGenerator2dQuad.h"
#include "MeshGenerator2dQMorph.h"
#include "MeshContainer2d.h"
#include "FrontEdge.h"
#include "FrontLines.h"
#include "MeshEdge2d.h"
#include "MeshPoint2d.h"
#include "MeshTriangle2d.h"
#include "MeshQuad2d.h"
#include "SurfaceParametric.h"
#include "ControlSpace2d.h"
#include "DEquation.h"
#include "DataVector.h"
#include "DTriangle.h"

int MeshGenerator2dQMorph::param_size_check = 0;

//////////////////////////////////////////////////////////////////////
// Smoothen mesh ahead of the quad merging front
void MeshGenerator2dQMorph::smoothenAhead(Metric2dContext& mc, MeshContainer2d* mesh, 
		MeshPoint2d *point)
{
	assert(mesh->isValid());

	const DMPoint2d mdpt = point->getMetricCoordinates(mc);
	while(true){
		bool special_case = false;
		if(point->getRank() <= 4 && point->getElementsCount(3) == 2){ // two triangles and two quads
			MeshEdge2d* short_edge = nullptr;
			double min_len = -1.0;
			for(int j = 0; j < point->getRank(); j++){
				MeshEdge2d* edge = point->getEdge(j);
				const MeshElement* element1 = edge->getMeshElement(0);
				const MeshElement* element2 = edge->getMeshElement(1);
				if(element1 && element1->getEdgeCount() == 3 &&
				   element2 && element2->getEdgeCount() == 3) {
					short_edge = edge;
				}else{
					double len = mdpt.distance2(edge->getOtherPoint(point)->getMetricCoordinates(mc));
					if(min_len < 0.0 || len < min_len) min_len = len;
				}
			}
			MeshPoint2d* point1 = short_edge ? short_edge->getOtherPoint(point) : nullptr;
			if(short_edge && !short_edge->isBorder() && 
				!point1->isBorder() &&
				point1->zeroIntTag(TagExtended::TAG_QUAD) && 
				!short_edge->isSideTagged() &&
				mdpt.distance2(point1->getMetricCoordinates(mc)) < 0.36*min_len)
			{
				MeshElement* triangle1 = short_edge->getMeshElement(0);
				MeshElement* triangle2 = short_edge->getMeshElement(1);
				bool allowed = true;
				const DPoint2d old_coordinates = point1->getCoordinates();
				point1->setCoordinates(point->getCoordinates());
				for(int i = 0; i < point1->getRank(); i++){
					const MeshElement* element = point1->getEdge(i)->getMeshElement(point1);
					if(element != triangle1 && element != triangle2 && element->isInverted()){
						allowed = false;
						break;
					}
				}
				point1->setCoordinates(old_coordinates);
				if(allowed){
					delete mesh->removeMeshElement(triangle1);
					delete mesh->removeMeshElement(triangle2);
//					MeshView::showDebugMesh("xxx", mesh, point1, point);
					while(point1->getRank() > 0){
						MeshElement* element = point1->getEdge(0)->getMeshElement(0);
						if(!element) element = point1->getEdge(0)->getMeshElement(1);
						assert(element);
						element->switchPointsWithEdges(point1, point);
					}
					delete mesh->removeMeshPoint(point1);
					assert(mesh->isValid());
					special_case = true;
				}
			}
		}
		for(int j = 0; j < point->getRank(); j++){
			MeshPoint2d* other_point = point->getEdge(j)->getOtherPoint(point);
			if(other_point->isBorder() || 
				other_point->nonZeroIntTag(TagExtended::TAG_QUAD) ||
				other_point->getElementsCount(4) > 0) continue;
			if(other_point->getRank() < 5){
				MeshGenerator2d::removeTriangulationPoint(mc, mesh, other_point);
				j = 0;
			}else
				MeshGenerator2dQuad::smoothenNode(mc, other_point);
		}
		for(int j = 0; j < point->getRank(); j++){
			MeshTriangle2d* triangle = (MeshTriangle2d*)point->getEdge(j)->getMeshElement(point);
			if(triangle && triangle->getEdgeCount() == 3) 
				MeshGenerator2d::optimizeTriangleBySwapPIng(mc, triangle);
		}
		assert(mesh->isValid());
		if(!special_case) break;
	}
}

//////////////////////////////////////////////////////////////////////
// Removes all triangles (and points if obsolete) within the area
// surrounded by front edges
int MeshGenerator2dQMorph::removeInnerTriangles(MeshContainer2d *mesh, MeshTriangle2d *start_triangle,
		const DataVector<MeshEdge2d*> & border_edges)
{
	if(start_triangle == nullptr) return 0;

	DataCompoundList<MeshTriangle2d*> remove_triangles;
	DataCompoundList<MeshPoint2d*> remove_points;

	// Gather triangles to be removed
	start_triangle->setIntTag(TagExtended::TAG_QUAD);
	remove_triangles.append(start_triangle);
	for (auto it = remove_triangles.iterator(); it.valid(); it.moveNext()) {
		MeshTriangle2d* triangle = it.item();
		for(int i = 0; i < 3; i++){
			MeshEdge2d* edge = triangle->getEdge(i);
			if(border_edges.contains(edge)) continue;
			MeshTriangle2d *second_triangle = (MeshTriangle2d*)triangle->getNeighbour(i);
			if(!second_triangle || second_triangle->nonZeroIntTag(TagExtended::TAG_QUAD)) continue;
			if(edge->isBorder() || second_triangle->getType() != ELEMENT_MESH_TRIANGLE){
				for (auto itt = remove_triangles.iterator(); it.valid(); it.moveNext())
					it.item()->removeTag(TagExtended::TAG_QUAD);
				return 0;
			}
			second_triangle->setIntTag(TagExtended::TAG_QUAD);
			remove_triangles.append(second_triangle);
		}
	}
	// Delete triangles
	for (auto it = remove_triangles.iterator(); it.valid(); it.moveNext()) {
		MeshTriangle2d* triangle = it.item();
		// Check points
		for(int i = 0; i < 3; i++){
			MeshPoint2d* point = triangle->getPoint(i);
			// This point is to be checked for removal later
			if(point->getRank() < 3)
				remove_points.append(point);
		}
		delete mesh->removeMeshElement(triangle);
	}
	// Delete points
	for (auto it = remove_points.iterator(); it.valid(); it.moveNext()) {
		if(it.item()->getRank() == 0)
			delete mesh->removeMeshPoint(it.item());
	}
	return remove_triangles.countInt();
}

MeshEdge2d* MeshGenerator2dQMorph::findNeededEdge(Metric2dContext& mc, MeshContainer2d *mesh, 
			FrontLines *front, MeshGenerator2dQuad::QuadConversionContext * qc, MeshEdge2d *edge, 
			MeshPoint2d *point, bool to_left, const DVector2d &dv, double average_length)
{
	int i = edge->getPointIndex(point);
	MeshEdge2d *next_edge, *current_edge = edge;

	const DMPoint2d dpt0 = point->getMetricCoordinates(mc);
	const DPoint2d ppt1 = point->getCoordinates() + dv;
	const DMPoint2d dpt1 = mc.transformPStoMS(ppt1);
	MeshPoint2d* pt2 = edge->getOtherPoint(point);
	MeshPoint2d* pt3 = nullptr;
	double det1 = DTriangle2d::det(dpt0, dpt1, pt2->getMetricCoordinates(mc));

	// Find edges between which the optimum edge (pt0 -> pt1) is found
	MeshElement *element = to_left ? edge->getMeshElement(i) : edge->getMeshElement(1-i);
	while(true){
		if(!element) return nullptr;
		//LOG4CPLUS_ASSERT(MeshLog::logger_mesh, element, MeshingException("MG2QM:findNeededEdge - invalid element"), 0);
		next_edge = to_left ? element->getPrevEdge(current_edge) : element->getNextEdge(current_edge);
		pt3 = next_edge->getOtherPoint(point);
		double det2 = DTriangle2d::det(dpt0, dpt1, pt3->getMetricCoordinates(mc));
		if(det1 * det2 <= 0) break;
		det1 = det2;
		pt2 = pt3;
		element = next_edge->getOtherElement(element);
		current_edge = next_edge;
		if(next_edge == edge) return nullptr;
	}
	MeshEdge2d* second_edge = to_left ? element->getNextEdge(current_edge) : element->getPrevEdge(current_edge);

	// If edge ahead belongs to front, allowed difference is increased
	double max_angle_cos = ANGLE_DIFFERENCE;
	if(second_edge->isSideTagged()) max_angle_cos *= 1.5;
	// Edge is between current_edge and next_edge

	// Count angles between both edges and the optimum
	double cos_angle2 = DMVector2d::angleCos(dpt0, dpt1, pt2->getMetricCoordinates(mc));
	double cos_angle3 = DMVector2d::angleCos(dpt0, dpt1, pt3->getMetricCoordinates(mc));

	double dist2 = dpt0.distance(pt2->getMetricCoordinates(mc));
	double dist3 = dpt0.distance(pt3->getMetricCoordinates(mc));

	if(cos_angle2 > 2) cos_angle2 = 4 - cos_angle2;
	if(cos_angle3 > 2) cos_angle3 = 4 - cos_angle3;

	MeshEdge2d* best_edge = nullptr;
	MeshEdge2d* second_best_edge = nullptr;
	double best_len = 0.0, second_best_len = 0.0;
	double best_angle_cos = 0.0;
	if(cos_angle2 < cos_angle3){
		if((cos_angle2 < max_angle_cos) && 
			(dist2 <= qc->qmorph_max_ratio * average_length || cos_angle2 < ANGLE_1) &&
			(dist2 >  qc->qmorph_min_ratio * average_length) || 
			current_edge->getOtherPoint(point)->nonZeroIntTag(TagExtended::TAG_QUAD))
		{
			best_edge = current_edge;
			best_len = dist2;
			best_angle_cos = cos_angle2;
		}
		if((cos_angle3 < max_angle_cos) && 
			(dist3 <= qc->qmorph_max_ratio * average_length || cos_angle3 < ANGLE_1) &&
			(dist3 >  qc->qmorph_min_ratio * average_length) || 
			next_edge->getOtherPoint(point)->nonZeroIntTag(TagExtended::TAG_QUAD)){
			if(best_edge){
				second_best_edge = next_edge;
				second_best_len = dist3;
			}else{
				best_edge = next_edge;
				best_len = dist3;
				best_angle_cos = cos_angle3;
			}
		}
	}else{
		if((cos_angle3 < max_angle_cos) && 
			(dist3 <= qc->qmorph_max_ratio * average_length || cos_angle3 < ANGLE_1) &&
			(dist3 >  qc->qmorph_min_ratio * average_length) || 
			next_edge->getOtherPoint(point)->nonZeroIntTag(TagExtended::TAG_QUAD))
		{
			best_edge = next_edge;
			best_len = dist3;
			best_angle_cos = cos_angle3;
		}
		if((cos_angle2 < max_angle_cos) && 
			(dist2 <= qc->qmorph_max_ratio * average_length || cos_angle2 < ANGLE_1) &&
			(dist2 >  qc->qmorph_min_ratio * average_length) || 
			current_edge->getOtherPoint(point)->nonZeroIntTag(TagExtended::TAG_QUAD))
		{
			if(best_edge){
				second_best_edge = current_edge;
				second_best_len = dist2;
			}else{
				best_edge = current_edge;
				best_len = dist2;
				best_angle_cos = cos_angle2;
			}
		}
	}

//	if(best_edge) SHOW_STEP_PT(3, "Best candidate for edge.", mesh->getSurface()->getPoint(best_edge->getPoint(0.5)));
//	if(second_best_edge) SHOW_STEP_PT(3, "Second-best candidate for edge.", mesh->getSurface()->getPoint(second_best_edge->getPoint(0.5)));

	if(best_edge){
		MeshElement* element1 = best_edge->getMeshElement(0);
		MeshElement* element2 = best_edge->getMeshElement(1);
		if(element1 && element1->getEdgeCount() != 3 && 
			element2 && element2->getEdgeCount() != 3){
//			SHOW_STEP_PT(0, "Critical error while selecting edge !.", mesh->getSurface()->getPoint(best_edge->getPoint(0.5)));
			return nullptr;
		}
	}

	if(second_best_edge){
		MeshElement* element1 = second_best_edge->getMeshElement(0);
		MeshElement* element2 = second_best_edge->getMeshElement(1);
		if(element1 && element1->getEdgeCount() != 3 && 
			element2 && element2->getEdgeCount() != 3){
//			SHOW_STEP_PT(0, "Critical error while selecting edge !.", mesh->getSurface()->getPoint(second_best_edge->getPoint(0.5)));
			return nullptr;
		}
	}

	// Try to join fronts
	if(best_edge && best_edge->getOtherPoint(point)->nonZeroIntTag(TagExtended::TAG_QUAD)){
		if(second_best_edge && second_best_edge->getOtherPoint(point)->zeroIntTag(TagExtended::TAG_QUAD))
			return second_best_edge;
		if(!best_edge->isSideTagged() && best_len < 1.5*average_length){
			MeshPoint2d* second_point = best_edge->getOtherPoint(point);
			// Check parity of front cycles
			if(MeshGenerator2dQuad::joinFrontConditionally(mc, front, point, second_point))
				return best_edge;
		}
		if(second_best_edge && !second_best_edge->isSideTagged() && second_best_len < 1.5*average_length){
			MeshPoint2d* second_point = second_best_edge->getOtherPoint(point);
			// Check parity of front cycles
			if(MeshGenerator2dQuad::joinFrontConditionally(mc, front, point, second_point))
				return second_best_edge;
		}
		if(best_angle_cos < ANGLE_1){
			if(best_len > average_length && !best_edge->isSideTagged()){
				// insert point
				const DPoint2d pt = best_edge->getPoint(0.5);
				MeshPoint2d* new_point = new MeshPoint2d(pt);
				mesh->addMeshPoint(new_point);
				MeshElement *element1 = best_edge->getMeshElement(point);
				MeshElement *element2 = best_edge->getOtherElement(element1);
				MeshPoint2d* point1 = element1->getPrevPoint(point);
				MeshPoint2d* point2 = element2->getNextPoint(point);
				MeshPoint2d *second_point = best_edge->getOtherPoint(point);
				int id = element1->getAreaID();
				delete mesh->removeMeshElement(element1);
				delete mesh->removeMeshElement(element2);
				LOG_ASSERT(point->getEdgeToPoint(second_point) == nullptr);
				// MeshingException("MG2QM:removeInnerTriangles - error new point"), nullptr);
				// Four new triangles ...
				MeshTriangle2d *triangle;
				mesh->addMeshElement(triangle = new MeshTriangle2d(point, point2, new_point));
				triangle->setAreaID(id);
				mesh->addMeshElement(triangle = new MeshTriangle2d(point2, second_point, new_point));
				triangle->setAreaID(id);
				mesh->addMeshElement(triangle = new MeshTriangle2d(second_point, point1, new_point));
				triangle->setAreaID(id);
				mesh->addMeshElement(triangle = new MeshTriangle2d(point1, point, new_point));
				triangle->setAreaID(id);
				//SHOW_STEP_PT_BREAKABLE(3, "New point at the optimal edge [middle of edge].", mesh->getSurface()->getPoint(pt), nullptr);
				return point->getEdgeToPoint(new_point);
			}else{
				return nullptr;
			}
		}else{
			best_edge = nullptr;
		}
	}

	if(!best_edge){
		// Reconstruct edge
		if(second_edge->isSideTagged()){
			// Insert point at the optimum line before edge ...
			MeshPoint2d* point1 = element->getPrevPoint(point);
			MeshPoint2d* point2 = element->getNextPoint(point);
			double area1 = DTriangle2d::det(dpt0, point1->getMetricCoordinates(mc),
				point2->getMetricCoordinates(mc));
			double area2 = DTriangle2d::det(dpt1, point2->getMetricCoordinates(mc),
				point1->getMetricCoordinates(mc));
//			assert(area1 * area2 > 0.0);
			double t = area1 / (area1 + area2);
			const DPoint2d cross_pt = point->getCoordinates() + dv * (0.5*t);
			if(!second_edge->isBorder() && 
				dpt0.distance(mc.transformPStoMS(cross_pt)) < qc->qmorph_min_ratio)
				return nullptr;
			MeshPoint2d* new_point = new MeshPoint2d(cross_pt);
			mesh->addMeshPoint(new_point);
			MeshGenerator2d::addPointToTriangulationByAngles(mc, mesh, new_point, (MeshTriangle2d*)element, false);
			//SHOW_STEP_PT_BREAKABLE(3, "New point at the optimal edge [before front].", mesh->getSurface()->getPoint(pt), nullptr);
			return point->getEdgeToPoint(new_point);
		}else{
			// will edge-swap be sufficient ?
			MeshPoint2d* point1 = to_left ? next_edge->getOtherPoint(point) : current_edge->getOtherPoint(point);
			MeshElement* second_element = second_edge->getOtherElement(element);
			MeshPoint2d *second_point = second_element->getPrevEdge(second_edge)->getOtherPoint(point1);
			MeshPoint2d* point2 = second_edge->getOtherPoint(point1);
			double cos_angle = DMVector2d::angleCos(dpt0, dpt1, second_point->getMetricCoordinates(mc));
			double dist = dpt0.distance(second_point->getMetricCoordinates(mc));
			if(cos_angle > 2) cos_angle = 4 - cos_angle;
			// Is edge (diagonal) OK ?
			if((cos_angle < max_angle_cos) && 
				(dist < qc->qmorph_max_ratio * average_length) &&
				(dist > qc->qmorph_min_ratio * average_length) &&
				MeshGenerator2d::makeEdgeByEdgesSwapPIng(mc, point, second_point, -3))
			{
				if(second_point->nonZeroIntTag(TagExtended::TAG_QUAD)){
					MeshEdge2d* edge_sec = point->getEdgeToPoint(second_point);
					// Merging fronts possible ?
					if(!edge_sec->isSideTagged() && edge_sec->getLength(mc, false) < 2*average_length &&
							MeshGenerator2dQuad::joinFrontConditionally(mc, front, point, second_point))
						return edge_sec;
				}else{
					// Edge is OK - swap
					//SHOW_STEP_PT_BREAKABLE(3, "Switching diagonal.", mesh->getSurface()->getPoint(dpt4), nullptr);
					return point->getEdgeToPoint(second_point);
				}
			}
			// if no, split two triangles into four through inserting new point
			//	at the optimum edge line
			// new_point = (point1, point2) x (pt0, pt1)
			double area1 = DTriangle2d::det(dpt0, point1->getMetricCoordinates(mc),
				point2->getMetricCoordinates(mc));
			double area2 = DTriangle2d::det(dpt1, point2->getMetricCoordinates(mc),
				point1->getMetricCoordinates(mc));
//			if(area1 * area2 <= 0.0)
//				SHOW_MESH("xxx", mesh->getViewSet(nullptr, true, true, true));
//			assert(area1 * area2 > 0.0);
			double t = area1 / (area1 + area2);
			const DPoint2d cross_pt = point->getCoordinates() + dv * t;
			// Check length
			dist = dpt0.distance(mc.transformPStoMS(cross_pt));
			if(dist < qc->qmorph_min_ratio * average_length ){
				if(cos_angle2 < cos_angle3){
					best_edge = current_edge;
					best_len = dist2;
					second_best_edge = next_edge;
					second_best_len = dist3;
				}else{
					best_edge = next_edge;
					best_len = dist3;
					second_best_edge = current_edge;
					second_best_len = dist2;
				}

				// Try merging fronts
				if(best_edge->getOtherPoint(point)->nonZeroIntTag(TagExtended::TAG_QUAD)){
					if(!best_edge->isSideTagged()){
						MeshPoint2d* spoint = best_edge->getOtherPoint(point);
						// check parity
						if(best_len < 2*average_length && 
								MeshGenerator2dQuad::joinFrontConditionally(mc, front, point, spoint))
							return best_edge;
					}
					if(second_best_edge->getOtherPoint(point)->nonZeroIntTag(TagExtended::TAG_QUAD)){
						if(!second_best_edge->isSideTagged()){
							MeshPoint2d* spoint = second_best_edge->getOtherPoint(point);
							// check parity
							if(second_best_len < 2*average_length &&
									MeshGenerator2dQuad::joinFrontConditionally(mc, front, point, spoint))
								return second_best_edge;
						}
					}
					return nullptr;
				}
				return best_edge;
			}
			MeshPoint2d* new_point = new MeshPoint2d(cross_pt);
			mesh->addMeshPoint(new_point);
			int id = element->getAreaID();
			delete mesh->removeMeshElement(element);
			delete mesh->removeMeshElement(second_element);
			// Four new triangles ...
			MeshTriangle2d *triangle;
			mesh->addMeshElement(triangle = new MeshTriangle2d(point, point2, new_point));
			triangle->setAreaID(id);
			mesh->addMeshElement(triangle = new MeshTriangle2d(point2, second_point, new_point));
			triangle->setAreaID(id);
			mesh->addMeshElement(triangle = new MeshTriangle2d(second_point, point1, new_point));
			triangle->setAreaID(id);
			mesh->addMeshElement(triangle = new MeshTriangle2d(point1, point, new_point));
			triangle->setAreaID(id);
			//SHOW_STEP_PT_BREAKABLE(3, "New point at the optimal edge [crossing].", mesh->getSurface()->getPoint(pt), nullptr);
			return point->getEdgeToPoint(new_point);
		}
	}else{
		//SHOW_STEP_PT_BREAKABLE(3, "Found edge close to optimal.", mesh->getSurface()->getPoint(best_edge->getPoint(0.5)), nullptr);
	}
	return best_edge;
}

//////////////////////////////////////////////////////////////////////
// Smoothen mesh by improving quality of quads adjacent to the given edge
bool MeshGenerator2dQMorph::smoothenFrontPoint(Metric2dContext& mc, MeshEdge2d *edge, MeshPoint2d *point)
{
	if(point->isBorder() || point->getElementsCount(4) > 2) return false;
	int index = edge->getPointIndex(point);
	MeshElement* left_element = edge->getMeshElement(1-index);
	MeshElement* right_element = edge->getMeshElement(index);
	if(!left_element || !right_element) return false;
	if(left_element->getEdgeCount() != 4 || right_element->getEdgeCount() != 4) return false;

	// Check orientation of side edges for these quads
	MeshPoint2d* point_leftB = left_element->getNextPoint(point);
	MeshPoint2d* point_leftA = left_element->getNextPoint(point_leftB);
	MeshPoint2d* point_rightB = right_element->getPrevPoint(point);
	MeshPoint2d* point_rightA = right_element->getPrevPoint(point_rightB);
	MeshPoint2d* point_bottom = edge->getOtherPoint(point);

	const DMPoint2d ptB_left = point_leftB->getMetricCoordinates(mc);
	const DMPoint2d ptA_left = point_leftA->getMetricCoordinates(mc);
	const DMPoint2d ptB_right = point_rightB->getMetricCoordinates(mc);
	const DMPoint2d ptA_right = point_rightA->getMetricCoordinates(mc);
	const DMPoint2d ptB_middle = point->getMetricCoordinates(mc);
	const DMPoint2d ptA_middle = point_bottom->getMetricCoordinates(mc);
	const DMVector2d vAB_left  = ptB_left - ptA_left;
	const DMVector2d vAB_right = ptB_right - ptA_right;
	double cos_angle = vAB_left.getAngleCos(vAB_right);
	if(cos_angle > 0.5 && cos_angle < 3.5) return false; // angle greater than 45 degrees -> front too curved

	// Change
	double ul = ptB_middle.distance(ptB_left);	// length of left upper edge
	double ur = ptB_middle.distance(ptB_right);	//	... right upper
	double bl = ptA_middle.distance(ptA_left);	//  ... left bottom
	double br = ptA_middle.distance(ptA_right);	//	... right bottom
	double diff_upper_length = (ul*br - bl*ur) / (br + bl);
	if(diff_upper_length < METRIC_SMALL_NUMBER) return false;

	double average_edge_length = 0.5 * (vAB_left.length() + vAB_right.length());

	DMVector2d dv = ptB_middle - ptA_middle;
	// v - orthogonal vector for adjusting upper edges
	DMVector2d v(-dv.y, dv.x);
	v *= (diff_upper_length / v.length());
	dv += v;	// turn vector dv with orthogonal vector v
	dv *= (average_edge_length / dv.length());	// scale

	return MeshGenerator2d::tryMovingPoint(point, mc.transformMStoPS(ptA_middle + dv));
}

bool MeshGenerator2dQMorph::createQMorphQuad(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines *front, 
		MeshGenerator2dQuad::QuadConversionContext* qc)
{
	qc->setMergeMethod(MeshData::QUADS_QMORPH);

	smoothenAhead(mc, mesh, qc->mp[0]);
	smoothenAhead(mc, mesh, qc->mp[1]);

	// Check side edges 
	qc->checkSideEdges();
	// Check seams
	if(checkSeams(mc, mesh, front, qc, 0)) return true;
	if(checkSeams(mc, mesh, front, qc, 1)) return true;

	// Check if the upper front edge can be used
	if(!checkQuadUpperEdge(mc, qc, 0))
		checkQuadUpperEdge(mc, qc, 1);

	if(!setQuadSideEdge(mc, mesh, front, qc, 0)){
#ifdef QUAD_DEBUG_LOG
		if(false){
			SHOW_MESH("Couldn't set left edge", 
				mesh->getViewSet(nullptr, true, true, true), 2);
		}
#endif
		return false;
	}
	if(!setQuadSideEdge(mc, mesh, front, qc, 1)){
		if(!qc->side_fedges_available[FrontEdge::EDGE_LEFT]) 
			delete qc->side_fedges[FrontEdge::EDGE_LEFT];
#ifdef QUAD_DEBUG_LOG
		if(false){
			SHOW_MESH("Couldn't set right edge", 
				mesh->getViewSet(nullptr, true, true, true), 2);
		}
#endif
		return false;
	}

	if(MeshGenerator2dQuad::checkForSingleFrontTriangle(mc, front, qc)) return true;

#ifdef QUAD_DEBUG_LOG
	if(qc->mp[2] == qc->mp[3]){
		stringstream sstr;
		sstr << "Same upper quad points: " << qc->mp[2]->getID() << endl;
		sstr << "   for " << *qc->base_fedge;
		SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true));
	}
#endif
	LOG_ASSERT(qc->mp[2] != qc->mp[3]);
	//	MeshingException("MG2Q:createQMorphQuad - same upper quad points"), false);

	if(checkForInnerFrontEdges(mc, qc, 0)){
		if(!qc->side_fedges_available[FrontEdge::EDGE_LEFT]) 
			delete qc->side_fedges[FrontEdge::EDGE_LEFT];
		if(!qc->side_fedges_available[FrontEdge::EDGE_RIGHT]) 
			delete qc->side_fedges[FrontEdge::EDGE_RIGHT];
#ifdef QUAD_DEBUG_LOG
		if(false){
			SHOW_MESH("Found inner front edges from left side", 
				mesh->getViewSet(nullptr, true, true, true), 2);
		}
#endif
		return false;
	}
	if(checkForInnerFrontEdges(mc, qc, 1)){
		if(!qc->side_fedges_available[FrontEdge::EDGE_LEFT]) 
			delete qc->side_fedges[FrontEdge::EDGE_LEFT];
		if(!qc->side_fedges_available[FrontEdge::EDGE_RIGHT]) 
			delete qc->side_fedges[FrontEdge::EDGE_RIGHT];
#ifdef QUAD_DEBUG_LOG
		if(false){
			SHOW_MESH("Found inner front edges from right side", 
				mesh->getViewSet(nullptr, true, true, true), 2);
		}
#endif
		return false;
	}

	if(!setQuadUpperEdge(mc, qc)){
#ifdef QUAD_DEBUG_LOG
		if(false){
			SHOW_MESH("Error setting upper edge", 
				mesh->getViewSet(nullptr, true, true, true), 2);
		}
#endif
		return false;
	}

	if(!formNewQuad(mc, mesh, qc)){
#ifdef QUAD_DEBUG_LOG
		if(false){
			SHOW_MESH("Error forming new quad", 
				mesh->getViewSet(nullptr, true, true, true), 2);
		}
#endif
		return false;
	}

	return MeshGenerator2dQuad::updateQuadFrontEdges(mc, front, qc);
}

bool MeshGenerator2dQMorph::formNewQuad(Metric2dContext& mc, MeshContainer2d* mesh,  
		MeshGenerator2dQuad::QuadConversionContext * qc)
{
	// Remove triangles from withing the quad
	MeshElement* start_element = qc->base_fedge->getEdge()->getMeshElement(qc->mp[0]);
	DataVector<MeshEdge2d*> border_edges(4);
	border_edges.add(qc->base_fedge->getEdge());
	border_edges.add(qc->side_fedges[FrontEdge::EDGE_LEFT]->getEdge());
	border_edges.add(qc->side_fedges[FrontEdge::EDGE_RIGHT]->getEdge());
	border_edges.add(qc->side_fedges[FrontEdge::EDGE_UPPER]->getEdge());
	int removed_triangles_count = removeInnerTriangles(mesh, (MeshTriangle2d*)start_element, border_edges);
	if(removed_triangles_count < 1) return false;

#ifdef QUAD_DEBUG_LOG
	if(false){
		stringstream sstr;
		sstr << "Removed inner triangles -> " << removed_triangles_count;
		SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true));
	}
#endif
	// New quad
	qc->base_element = new MeshQuad2d(qc->mp[0], qc->mp[1], qc->mp[2], qc->mp[3]);
	qc->base_element->setAreaID(qc->element_id);
	mesh->addMeshElement(qc->base_element);

	// Smoothen points before front
	smoothenAhead(mc, mesh, qc->mp[2]);
	smoothenAhead(mc, mesh, qc->mp[3]);

	MeshGenerator2dQuad::swapTrianglesAhead(mc, qc);

	//SHOW_MESH("Qmorph quad - new", mesh->getViewSet(nullptr, true, true, true), 1);

	mc.countMetricAtPoint(qc->base_element->getMiddlePoint());
	// Smoothen quad
	for(int i = 0; i < 2; i++){
		//if(qc->mp[3-i]->isBorder() || qc->side_fedges_available[i]) continue;
		if(qc->mp[3-i]->isBorder()) continue;
		// calculate desired length and direction
		const DMVector2d dv = (qc->mp[3-i]->getMetricCoordinates(mc) - 
			qc->mp[i]->getMetricCoordinates(mc)).normalized();
		// try moving there
		MeshGenerator2d::tryMovingPoint(qc->mp[3-i], 
			qc->mp[i]->getCoordinates() + mc.transformMStoPS(dv));
	}
	//SHOW_MESH("Qmorph quad - adjust", mesh->getViewSet(nullptr, true, true, true), 1);

	if(qc->base_element->getAlphaQuality(mc) < 0.3){
		MeshGenerator2dQuad::smoothenNode(mc, qc->mp[2]);
		MeshGenerator2dQuad::smoothenNode(mc, qc->mp[3]);
	}else{
		for(int i = 0; i < 2; i++)
			if(qc->side_fedges_available[i]){
				MeshQuad2d* side_quad = (MeshQuad2d*)qc->side_fedges[i]->getEdge()->getOtherElement(qc->base_element);
				if(side_quad && side_quad->getEdgeCount() == 4 && side_quad->getAlphaQuality(mc) < 0.3)
					MeshGenerator2dQuad::smoothenNode(mc, qc->mp[3-i]);
			}
	}

	if(param_size_check){
		const double MAX_LEN_DIFF = 0.4;
		for(int i = 0; i < 2; i++){
			if(qc->mp[3-i]->isBorder()) continue;
			const DMVector2d dv = qc->mp[3-i]->getMetricCoordinates(mc) - qc->mp[i]->getMetricCoordinates(mc);
			double len = dv.length();
			if(abs(len - 1.0) > MAX_LEN_DIFF)
				MeshGenerator2d::tryMovingPoint(qc->mp[3-i], 
					qc->mp[i]->getCoordinates() + mc.transformMStoPS(dv/len));
		}
		if(!qc->mp[2]->isBorder() || !qc->mp[3]->isBorder()){ // check upper
			const DMVector2d dv = qc->mp[3]->getMetricCoordinates(mc) - qc->mp[2]->getMetricCoordinates(mc);
			double len = dv.length();
			if(abs(len - 1.0) > MAX_LEN_DIFF){
				const DVector2d ddv = mc.transformMStoPS(dv * ((len-1.0)/len));
				if(qc->mp[2]->isBorder()){ // move mp[3] -> mp[2]
					MeshGenerator2d::tryMovingPoint(qc->mp[3], qc->mp[2]->getCoordinates() + ddv);
				}else if(qc->mp[3]->isBorder()){ // move mp[2] -> mp[3]
					MeshGenerator2d::tryMovingPoint(qc->mp[2], qc->mp[3]->getCoordinates() - ddv);
				}else{ // move both
					const DPoint2d pt3 = qc->mp[2]->getCoordinates() + ddv * 0.5;
					const DPoint2d pt2 = qc->mp[3]->getCoordinates() - ddv * 0.5;
					MeshGenerator2d::tryMovingPoint(qc->mp[3], pt3);
					MeshGenerator2d::tryMovingPoint(qc->mp[2], pt2);
				}
			}
		}
	}

	smoothenAhead(mc, mesh, qc->mp[2]);
	smoothenAhead(mc, mesh, qc->mp[3]);

	if(qc->side_fedges_available[0]){
		smoothenFrontPoint(mc, qc->side_fedges[0]->getEdge(), qc->mp[3]);
		//SHOW_MESH("Qmorph quad - side front adjust left", mesh->getViewSet(nullptr, true, true, true), 1);
	}
	if(qc->side_fedges_available[1]){
		smoothenFrontPoint(mc, qc->side_fedges[1]->getEdge(), qc->mp[2]);
		//SHOW_MESH("Qmorph quad - side front adjust right", mesh->getViewSet(nullptr, true, true, true), 1);
	}
	assert(mesh->isValid());


	for(int i = 2; i <= 3; i++){
		if(qc->mp[i]->isBorder()) continue;
		for(int j = 0; j < qc->mp[i]->getRank(); j++){
			MeshEdge2d* edge = qc->mp[i]->getEdge(j);
			if(!edge->isSideTagged()) continue; // check only front edges
			double dist2 = qc->mp[i]->getMetricCoordinates(mc).distance2(
				edge->getOtherPoint(qc->mp[i])->getMetricCoordinates(mc));
			if(dist2 < 0.25){
				MeshGenerator2dQuad::smoothenNode(mc, qc->mp[i]);
				break;
			}
		}
	}

#ifdef QUAD_DEBUG_LOG
	if(false){
		double len_base  = (qc->mp[0]->getMetricCoordinates(mc) - qc->mp[1]->getMetricCoordinates(mc)).length();
		double len_left  = (qc->mp[0]->getMetricCoordinates(mc) - qc->mp[3]->getMetricCoordinates(mc)).length();
		double len_right = (qc->mp[1]->getMetricCoordinates(mc) - qc->mp[2]->getMetricCoordinates(mc)).length();
		double len_upper = (qc->mp[2]->getMetricCoordinates(mc) - qc->mp[3]->getMetricCoordinates(mc)).length();
		stringstream sstr;
		sstr << "New qmorph quad + sm: [" << len_base << ',' << len_left << ',' << len_right 
			<< ',' << len_upper << ']';
		SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true), 1);
	}
#endif

	MeshGenerator2dQuad::swapTrianglesAhead(mc, qc);

	return true;
}

bool MeshGenerator2dQMorph::setQuadUpperEdge(Metric2dContext& mc, MeshGenerator2dQuad::QuadConversionContext * qc)
{
	MeshEdge2d* upper_edge = qc->mp[2]->getEdgeToPoint(qc->mp[3]);
	if(!upper_edge)
		if(MeshGenerator2d::makeEdgeByEdgesSwapPIng(mc, qc->mp[2], qc->mp[3], -3))
			upper_edge = qc->mp[2]->getEdgeToPoint(qc->mp[3]);
		else{
			for(int side = 0; side < 2; side++)
				if(!qc->side_fedges_available[side]) delete qc->side_fedges[side];
			return false;
		}

	assert(upper_edge);
	qc->side_fedges_available[FrontEdge::EDGE_UPPER] = upper_edge->isSideTagged(qc->mp[2]);
	if(qc->side_fedges_available[FrontEdge::EDGE_UPPER]){
		if(qc->side_fedges_available[FrontEdge::EDGE_LEFT]){
			qc->side_fedges[FrontEdge::EDGE_UPPER] = 
				qc->side_fedges[FrontEdge::EDGE_LEFT]->getSideFrontEdge(FrontEdge::EDGE_LEFT);
		}else if(qc->side_fedges_available[FrontEdge::EDGE_RIGHT]){
			qc->side_fedges[FrontEdge::EDGE_UPPER] = 
				qc->side_fedges[FrontEdge::EDGE_RIGHT]->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
		}else{
			delete qc->side_fedges[FrontEdge::EDGE_LEFT];
			delete qc->side_fedges[FrontEdge::EDGE_RIGHT];
			return false;
		}
	}else{
		qc->side_fedges[FrontEdge::EDGE_UPPER] = new FrontEdge(upper_edge, qc->mp[3], qc->base_fedge->getLevel()+1);
	}
	return true;
}

bool MeshGenerator2dQMorph::checkForInnerFrontEdges(Metric2dContext& mc, 
		MeshGenerator2dQuad::QuadConversionContext * qc, int side)
{
	if(!qc->side_fedges_available[side]) return false;

	MeshEdge2d* medge = qc->side_fedges[side]->getEdge();
	if(medge->getOtherPoint(qc->mp[3-side]) == qc->mp[2+side]) return false;

	double cos_angle = DMVector2d::angleCos(
				qc->mp[3-side]->getMetricCoordinates(mc),
				qc->mp[2+side]->getMetricCoordinates(mc), 
				qc->mp[side]->getMetricCoordinates(mc));
	if(cos_angle > 2) cos_angle = 4 - cos_angle;
	if(cos_angle - qc->side_fedges[side]->classify(mc)->getSideAngleCos(side) < METRIC_SMALL_NUMBER) 
		return false;

	return true;
}

bool MeshGenerator2dQMorph::setQuadSideEdge(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines *front, 
			MeshGenerator2dQuad::QuadConversionContext * qc, int side)
{
	if(qc->side_fedges[side]) return true; // already found

	const DMPoint2d dpt[2] = {
		qc->mp[0]->getMetricCoordinates(mc),
		qc->mp[1]->getMetricCoordinates(mc)};

	DMVector2d dmv;
	// set preferred direction
	if(qc->base_fedge->getSideAngleCos(side) < 1.9){ // < PI -> bisetion
		MeshPoint2d* other_point = qc->base_fedge->getSideFrontEdge(side)->getEdge()->getOtherPoint(qc->mp[side]);
		assert(other_point);
		const DMPoint2d dpt_other = other_point->getMetricCoordinates(mc);
		dmv = (dpt_other - dpt[side]).normalized() + 
			(dpt[1-side] - dpt[side]).normalized();
	}else{
		const DMVector2d v01 = dpt[1]-dpt[0];
		dmv = DMVector2d(-v01.y, v01.x);
	}
	// normalize in metric space
	const DVector2d dv = mc.transformMStoPS(dmv.normalized());

	// Calculate optimum point
	const DPoint2d opt_point = qc->mp[side]->getCoordinates() + dv;

	MeshEdge2d* side_edge = findNeededEdge(mc, mesh, front, qc, qc->base_fedge->getEdge(), 
		qc->mp[side], (side == FrontEdge::EDGE_LEFT), dv);

	int other_side = 1-side;
	if(!side_edge && qc->side_fedges[other_side]){
		if(qc->base_fedge->getSideFrontEdge(other_side)->
				getSideFrontEdge(other_side)->
				getSideFrontEdge(other_side)->
				getSideFrontEdge(other_side) == qc->base_fedge)
		{
			// Only two triangles surrounded by front -> to merge
			qc->side_fedges_available[side] = true;
			qc->side_fedges[side] = qc->base_fedge->getSideFrontEdge(side);
			side_edge = qc->side_fedges[side]->getEdge();
		}
	}

	if(!side_edge) // failure
		return false;

	qc->mp[3-side] = side_edge->getOtherPoint(qc->mp[side]);
	if(side_edge->isSideTagged(qc->mp[(side==FrontEdge::EDGE_LEFT)?3:1])){
		qc->side_fedges[side] = qc->base_fedge->getSideFrontEdge(side);
		LOG_ASSERT(side_edge == qc->side_fedges[side]->getEdge());
		//	MeshingException("MG2Q:setQuadSideEdge - inconsistent front data"), false);
		qc->side_fedges_available[side] = true;
	}else{
		MeshGenerator2d::tryMovingPoint(qc->mp[3-side], opt_point);
		qc->side_fedges[side] = new FrontEdge(side_edge, qc->mp[(side==FrontEdge::EDGE_LEFT)?0:2], 
			qc->base_fedge->getLevel()+1);
	}
	return true;
}

bool MeshGenerator2dQMorph::checkQuadUpperEdge(Metric2dContext& mc, MeshGenerator2dQuad::QuadConversionContext * qc, 
			int side)
{
	if(!qc->side_fedges_available[side] || qc->side_fedges_available[1-side]) return false;

	FrontEdge* fside = qc->side_fedges[side]->getSideFrontEdge(side);
	MeshEdge2d* side_side_edge = fside->getEdge();
	MeshPoint2d* mpt_otherside[2] = {qc->mp[1-side],			// lower
		side_side_edge->getOtherPoint(qc->mp[3-side])};		// upper
	if(mpt_otherside[0] == mpt_otherside[1]) return false;
	if(qc->side_fedges[side]->getSideAngleCos(side) > ANGLE_THRESHOLD) return false;

	const DMPoint2d dpt0 = mpt_otherside[0]->getMetricCoordinates(mc);
	const DMPoint2d dpt1 = mpt_otherside[1]->getMetricCoordinates(mc);
	if(dpt0.distance2(dpt1) > qc->qmorph_max_ratio2) return false;	// too long (expected length = 1)

	MeshEdge2d* edge_otherside = mpt_otherside[0]->getEdgeToPoint(mpt_otherside[1]);
	if(edge_otherside == nullptr){
		if(MeshGenerator2d::makeEdgeByEdgesSwapPIng(mc, mpt_otherside[0], mpt_otherside[1], -3))
			edge_otherside = mpt_otherside[0]->getEdgeToPoint(mpt_otherside[1]);
		else return false;
	}
	LOG_ASSERT(edge_otherside);
	// MeshingException("MG2Q:checkUpperFrontEdge - invalid otherside_edge"), false);

	// Set this edge as a quad-side edge from the other side
	qc->mp[2+side] = mpt_otherside[1]; // upper
	if(edge_otherside->isSideTagged(mpt_otherside[side])){	// here lower/upper depends on side
		qc->side_fedges[1-side] = qc->base_fedge->getSideFrontEdge(1-side);
		qc->side_fedges_available[1-side] = true;
	}else
		qc->side_fedges[1-side] = new FrontEdge(edge_otherside, mpt_otherside[1-side], 
										qc->base_fedge->getLevel()+1);

	return true;
}

MeshGenerator2dQMorph::SeamContext::SeamContext(Metric2dContext& mc, 
		MeshGenerator2dQuad::QuadConversionContext * qc, int side)
{
	edge1 = (side == FrontEdge::EDGE_LEFT) ? 
		qc->base_fedge->getSideFrontEdge(FrontEdge::EDGE_LEFT) : qc->base_fedge;
	edge2 = (side == FrontEdge::EDGE_LEFT) ? 
		qc->base_fedge : qc->base_fedge->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
	point = qc->mp[side];
	angle = qc->base_fedge->getSideAngleCos(side);
	mesh_edge1 = edge1->getEdge();
	mesh_edge2 = edge2->getEdge();
	id = mesh_edge2->getMeshElement(point)->getAreaID();
	level = edge1->getLevel();
	length1 = mesh_edge1->getLength(mc, false);
	length2 = mesh_edge2->getLength(mc, false);
	ratio = length1 / length2;
}

bool MeshGenerator2dQMorph::trySeam1(Metric2dContext& mc, MeshContainer2d* mesh, 
			FrontLines *front, SeamContext& sc)
{
	if(sc.edge1->getSideAngleCos(FrontEdge::EDGE_LEFT) < ANGLE_THRESHOLD ||
			sc.edge2->getSideAngleCos(FrontEdge::EDGE_RIGHT) < ANGLE_THRESHOLD ||
			sc.length1 > 2.0 || sc.length2 > 2.0)
		return false;

	// restore the edge between two points to-be-joined
	MeshPoint2d* mp1 = sc.mesh_edge1->getOtherPoint(sc.point);
	MeshPoint2d* mp2 = sc.mesh_edge2->getOtherPoint(sc.point);
	double dist2 = mp1->getMetricCoordinates(mc).distance2(mp2->getMetricCoordinates(mc));
//	LOG4CPLUS_DEBUG(MeshLog::logger_console, "Seam-1 dist",sqrt(dist2));
	if(dist2 > 1.0/METRIC_LENGTH_RATIO2) return false;
	bool to_right = true;
	if(mp1->isBorder())
		if(mp2->isBorder()) return false;
		else to_right = false;

	// (1) "Seam operation"
#ifdef QUAD_DEBUG_LOG
//	if(sc.point->getID() == 3069){
	if(false){
		stringstream sstr;
		sstr << "Seam-1 " << sc.point->getID() << ", len[seam-edge]=" 
			<< sqrt(dist2) << " [" << mp1->getID() << ',' << mp2->getID() << ']';
		SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true));
	}
#endif

	MeshEdge2d* edge3 = mp1->getEdgeToPoint(mp2);
	if(!edge3){
		bool result = MeshGenerator2d::makeEdgeByEdgesSwapPIng(mc, mp1, mp2, -3);
		if(!result) return false;
		edge3 = mp1->getEdgeToPoint(mp2);
	}
	LOG_ASSERT(edge3);
	// MeshingException("MG2QM:checkSeams - seam1, no edge"), false);
	if(edge3->isSideTagged()) return false;

	// remove triangles from within the seam
	int index = sc.mesh_edge2->getPointIndex(sc.point);
	MeshElement* element = sc.mesh_edge2->getMeshElement(index);

	DataVector<MeshEdge2d*> border_edges(3);
	border_edges.add(sc.edge1->getEdge());
	border_edges.add(sc.edge2->getEdge());
	border_edges.add(edge3);

	int ct = removeInnerTriangles(mesh, (MeshTriangle2d*)element, border_edges);
	LOG_ASSERT(ct > 0);
	// MeshingException("MG2QM:checkSeams - seam1, error removing triangles"), false);

	element = edge3->getMeshElement(mp1);
	MeshPoint2d* mp5 = element->getNextPoint(mp2);
	delete mesh->removeMeshElement(element);

	// join boundaries - reconstruct the topology
	DataVector<InvalidFrontEdge> invalid_list;
	MeshPoint2d* mpoint1 = to_right ? mp1 : mp2;
	MeshPoint2d* mpoint2 = to_right ? mp2 : mp1;

	sc.edge1->clear();
	sc.edge2->clear();

	// Join mpoint1 -> mpoint2
	ct = mpoint1->getRank();
	for(int i = 0; i < ct; i++){
		const MeshEdge2d* edge = mpoint1->getEdge(i);
		MeshPoint2d* other_point = edge->getOtherPoint(mpoint1);
		if(edge->isSideTagged(mpoint1)){
			FrontEdge* fedge1 = front->findFrontEdge(edge, edge->getPointIndex(mpoint1));
			if(fedge1 && fedge1 != sc.edge1 && fedge1 != sc.edge2){
				fedge1->clear();
				invalid_list.add(InvalidFrontEdge(fedge1, mpoint2, other_point));
			}
		}
		if(edge->isSideTagged(other_point)){
			FrontEdge* fedge2 = front->findFrontEdge(edge, edge->getPointIndex(other_point));
			if(fedge2 && fedge2 != sc.edge1 && fedge2 != sc.edge2){
				fedge2->clear();
				invalid_list.add(InvalidFrontEdge(fedge2, other_point, mpoint2));
			}
		}
	}

	while(mpoint1->getRank() > 0){
		const MeshEdge2d* edge = mpoint1->getEdge(0);
		MeshElement* element_0 = edge->getMeshElement(0);
		MeshElement* element_1 = edge->getMeshElement(1);

		if(element_0) element_0->switchPointsWithEdges(mpoint1, mpoint2);
		if(element_1) element_1->switchPointsWithEdges(mpoint1, mpoint2);
	}

	// Move point to the middle position and remove the other one
	if(!mpoint2->isBorder()){
		MeshGenerator2d::tryMovingPoint(mpoint2, 
			DPoint2d::average(mp1->getCoordinates(), mp2->getCoordinates()));
	}
	delete mesh->removeMeshPoint(mpoint1);

//	SHOW_DEBUG_MESH("Seam-1 - after switch", mesh, mpoint2);
//	if(sc.point->getID() == 3069)
//		MeshView::showDebugMesh("Seam-1 - after switch", mesh, mpoint2);

	// additionally check mesh
	int rank = mpoint2->getRank();
	for(int j = 0; j < rank; j++){
		const MeshElement* el = mpoint2->getEdge(j)->getMeshElement(mpoint2);
		if(!el) el = mpoint2->getEdge(j)->getOtherElement(nullptr);
		//assert(el);
		if(!el) continue;
		int n = 0;
		int nct = el->getEdgeCount();
		while(el->isInverted()){
//			SHOW_DEBUG_MESH("Seam-1 - inverted elements", mesh, mpoint2);
			MeshGenerator2dQuad::smoothenNode(mc, el->getPoint(n%nct));
			if(++n > 15) break;	// max 5 prób dla ka¿dego z trzech punktów
		}
	}

//	if(sc.point->getID() == 2274)
//		MeshView::showDebugMesh("Seam-1 - after validation", mesh, mpoint2);

	// mesh improvement
	MeshGenerator2dQuad::smoothenNode(mc, mp5);
	smoothenAhead(mc, mesh, mp5);
	smoothenAhead(mc, mesh, mpoint2);

	for(int i = 0; i < mpoint2->getRank(); i++){
		MeshTriangle2d* el = (MeshTriangle2d*)mpoint2->getEdge(i)->getMeshElement(mpoint2);
		if(!el) el = (MeshTriangle2d*)mpoint2->getEdge(i)->getOtherElement(nullptr);
		//assert(el);
		if(!el) continue;
		if(el->isInverted() && el->getEdgeCount()==3)
			for(int j = 0; j < 3; j++)
				if(el->swapWithNeighbour(mc, j, false)) break;
	}

//	SHOW_DEBUG_MESH("Seam-1 - after local improvement", mesh, mpoint2);
//	if(quad_conversion_step == 19480)
//		MeshView::showDebugMesh("Seam-1 - after local improvement", mesh, mpoint2);

	// Update front edges
	FrontEdge* left_edge  = sc.edge1->getSideFrontEdge(FrontEdge::EDGE_LEFT);
	FrontEdge* right_edge = sc.edge2->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
	front->removeFrontEdge(sc.edge1);
	front->removeFrontEdge(sc.edge2);
	bool last_quad = false;
	if(left_edge->getSideFrontEdge(FrontEdge::EDGE_LEFT) == right_edge){
		// last quad
		last_quad = true;
		front->removeFrontEdge(left_edge);
		front->removeFrontEdge(right_edge);
	}else{
		left_edge->setSideFrontEdge(right_edge, FrontEdge::EDGE_RIGHT);
		right_edge->setSideFrontEdge(left_edge, FrontEdge::EDGE_LEFT);
	}

	// update front edges -> mesh edge linkage + reclassify (with neighbours because of angles)
	DataVector<FrontEdge*> front_edges_to_update;
	while(invalid_list.countInt() > 0){
		InvalidFrontEdge node = invalid_list.removeLast();
		if(last_quad && (node.fedge == left_edge || node.fedge == right_edge)) continue;
		node.fedge->setEdge(node.mp1->getEdgeToPoint(node.mp2), node.mp1);
		front_edges_to_update.addIfNew(node.fedge);
		front_edges_to_update.addIfNew(node.fedge->getSideFrontEdge(FrontEdge::EDGE_LEFT));
		front_edges_to_update.addIfNew(node.fedge->getSideFrontEdge(FrontEdge::EDGE_RIGHT));
	}
	for(size_t i = 0; i < front_edges_to_update.countInt(); i++)
		front->updateFrontEdgePosition(front_edges_to_update[i]->classify(mc));

	LOG_ASSERT(front->isValid());
	// MeshingException("MG2Q:trySeam1, invalid front"), false);

#ifdef QUAD_DEBUG_LOG
	if(false){
		stringstream sstr;
		sstr << "Seam-1 closed " << sc.point->getID();
		SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true), 1);
	}
#endif
	return true;
}

bool MeshGenerator2dQMorph::trySeam2(Metric2dContext& mc, MeshContainer2d* mesh, 
			FrontLines *front, SeamContext& sc)
{
	// (2) "Transition seam operation"
#ifdef QUAD_DEBUG_LOG
	if(false){
//	if(sc.point->getID() == 1955){
		stringstream sstr;
		sstr << "Seam-2 " << sc.point->getID() << " ["
			<< sc.mesh_edge1->getOtherPoint(sc.point)->getID() << ','
			<< sc.mesh_edge2->getOtherPoint(sc.point)->getID() << "], len="
			<< sc.length1 << ',' << sc.length2 << ']';
		SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true));
	}
#endif

	//LOG4CPLUS_INFO(MeshLog::logger_console, "Seam 2");
	bool left_longer = sc.ratio > 1.0;
	MeshEdge2d* longer_edge  = left_longer ? sc.edge1->getEdge() : sc.edge2->getEdge();
	MeshEdge2d* shorter_edge = left_longer ? sc.edge2->getEdge() : sc.edge1->getEdge();
	MeshPoint2d* longer_point = longer_edge->getOtherPoint(sc.point);
	MeshPoint2d* short_point = shorter_edge->getOtherPoint(sc.point);
	if(longer_edge->isBorder()) return false;
			
	// Change shape of element adjacent to longer edge
	MeshElement* large_element = longer_edge->getMeshElement(left_longer ? sc.point : longer_edge->getOtherPoint(sc.point));
	if(!large_element || large_element->getEdgeCount() != 4 || large_element->getAlphaQuality(mc) < 0.1) 
		return false;
//	double large_area1 = large_element->getArea();
//	double large_alpha1 = large_element->getAlphaQuality(mc);
	MeshTriangle2d* large_triangle = (MeshTriangle2d*)longer_edge->getOtherElement(large_element);
//	double large_area2 = large_triangle->getArea();

//	SHOW_DEBUG_MESH2("Seam-2 - start", mesh, large_element, large_triangle);
//	MeshView::showDebugMesh("Seam-2", mesh, large_element, large_triangle);
//	if(quad_conversion_step == 19480)
//		MeshView::showDebugMesh("Seam-2", mesh, large_element, large_triangle);

	// Clear pointers to edge to-be-removed in further steps
	if(left_longer) sc.edge1->clear();
	else sc.edge2->clear();			

	// New point in the middle of the longer edge
	const DPoint2d middle = longer_edge->getPoint(0.5);
	MeshPoint2d* new_point = new MeshPoint2d(middle);
	mesh->addMeshPoint(new_point);
			
	large_element->switchPointsWithEdges(sc.point, new_point);
	large_triangle->switchPointsWithEdges(sc.point, new_point);

	LOG_ASSERT(!large_triangle->isInverted());
	// MeshingException("MG2QM:checkSeams - seam2, invalid large_triangle"), false);
			
	// Termporary reconstruction of triangular mesh within the seam
	MeshPoint2d* triangle_point = left_longer ? large_triangle->getNextPoint(new_point) : 
		large_triangle->getPrevPoint(new_point);
	MeshTriangle2d* new_triangle = new MeshTriangle2d(sc.point, 
				left_longer ? triangle_point : new_point, 
				left_longer ? new_point : triangle_point);
	LOG_ASSERT(!new_triangle->isInverted());
	new_triangle->setAreaID(sc.id);
	mesh->addMeshElement(new_triangle);
	MeshPoint2d* old_point = left_longer ? large_element->getPrevPoint(new_point) : 
		large_element->getNextPoint(new_point);
	MeshTriangle2d* second_triangle = new MeshTriangle2d(
				left_longer ? sc.point : new_point, 
				left_longer ? new_point : sc.point, 
				old_point);
	second_triangle->setAreaID(sc.id);
	mesh->addMeshElement(second_triangle);

//	if(sc.point->getID() == 1955)
//		MeshView::showDebugMesh("insert point", mesh, large_element);


	// Edge recovery
	MeshEdge2d* quad_edge = short_point->getEdgeToPoint(new_point);
	if(!quad_edge){
		bool result = MeshGenerator2d::makeEdgeByEdgesSwapPIng(mc, short_point, new_point, -3);
		if(!result){
			// Try to continue desPIte failure
#ifdef QUAD_DEBUG_LOG
			if(false){
				stringstream sstr;
				sstr << "Seam-2 failed edge recovery " << sc.point->getID();
				SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true), 2);
			}
#endif

			// Smoothing
			MeshGenerator2dQuad::smoothenNode(mc, old_point);	
			smoothenAhead(mc, mesh, old_point);
			MeshGenerator2dQuad::smoothenNode(mc, new_point);	
			smoothenAhead(mc, mesh, new_point);
			MeshGenerator2dQuad::smoothenNode(mc, longer_point);
			smoothenAhead(mc, mesh, longer_point);
			MeshGenerator2dQuad::smoothenNode(mc, short_point);
			smoothenAhead(mc, mesh, short_point);
//			SHOW_DEBUG_MESH("Seam-2 - failed edge swap, recovering, after smoothing", mesh, point);
			// .. update front
			MeshEdge2d* new_edge1 = new_point->getEdgeToPoint(old_point);
			MeshEdge2d* new_edge2 = old_point->getEdgeToPoint(sc.point);
			MeshEdge2d* long_edge1 = new_point->getEdgeToPoint(longer_point);
			LOG_ASSERT(new_edge1 && new_edge2);
			// MeshingException("MG2QM:checkSeams - seam2, bad edges"), false);
			FrontEdge *fedge1, *fedge2;
			if(left_longer){
				fedge1 = new FrontEdge(new_edge1, new_point, sc.level+3);
				fedge2 = new FrontEdge(new_edge2, old_point, sc.level+3);
				sc.edge2->setSideFrontEdge(fedge2, FrontEdge::EDGE_LEFT);
				fedge2->setSideFrontEdge(sc.edge2, FrontEdge::EDGE_RIGHT);
				fedge2->setSideFrontEdge(fedge1, FrontEdge::EDGE_LEFT);
				fedge1->setSideFrontEdge(fedge2, FrontEdge::EDGE_RIGHT);
				fedge1->setSideFrontEdge(sc.edge1, FrontEdge::EDGE_LEFT);
				sc.edge1->setSideFrontEdge(fedge1, FrontEdge::EDGE_RIGHT);
				sc.edge1->setEdge(long_edge1, longer_point);
			}else{
				fedge1 = new FrontEdge(new_edge1, old_point, sc.level+3);
				fedge2 = new FrontEdge(new_edge2, sc.point, sc.level+3);
				sc.edge2->setSideFrontEdge(fedge1, FrontEdge::EDGE_LEFT);
				fedge1->setSideFrontEdge(sc.edge2, FrontEdge::EDGE_RIGHT);
				fedge1->setSideFrontEdge(fedge2, FrontEdge::EDGE_LEFT);
				fedge2->setSideFrontEdge(fedge1, FrontEdge::EDGE_RIGHT);
				fedge2->setSideFrontEdge(sc.edge1, FrontEdge::EDGE_LEFT);
				sc.edge1->setSideFrontEdge(fedge2, FrontEdge::EDGE_RIGHT);
				sc.edge2->setEdge(long_edge1, new_point);
			}
			front->updateFrontEdgePosition(sc.edge1->classify(mc));
			front->updateFrontEdgePosition(sc.edge2->classify(mc));
			front->addFrontEdge(fedge1->classify(mc));
			front->addFrontEdge(fedge2->classify(mc));

			LOG_ASSERT(front->isValid());
			// MeshingException("MG2QM:checkSeams - seam2, invalid front"), false);
			return true;
		}
		quad_edge = short_point->getEdgeToPoint(new_point);
	}

	LOG_ASSERT(new_point->getEdgeToPoint(longer_point));
	// MeshingException("MG2Q:trySeam2, step1"), false); // !!!!!!

	MeshPoint2d* mpts[4] = {old_point, left_longer ? sc.point : new_point, short_point, 
		left_longer ? new_point: sc.point };
	DataVector<MeshEdge2d*> border_edges(4);
	for(int i = 0; i < 4; i++)
		border_edges.add(mpts[i]->getEdgeToPoint(mpts[(i+1)%4]));

	// Remove triangles from within the created quad
	removeInnerTriangles(mesh, second_triangle, border_edges);

	LOG_ASSERT(new_point->getEdgeToPoint(longer_point));
	// MeshingException("MG2Q:trySeam2, step2"), false); // !!!!!!

	// Create new quad
	MeshQuad2d* quad = new MeshQuad2d(mpts[0], mpts[1], mpts[2], mpts[3]);
	quad->setAreaID(sc.id);
	mesh->addMeshElement(quad);

//	if(sc.point->getID() == 1955)
//		MeshView::showDebugMesh("new quad", mesh, quad);

	LOG_ASSERT(new_point->getEdgeToPoint(longer_point));
	// MeshingException("MG2Q:trySeam2, step3"), false); // !!!!!!

	if(MeshGenerator2dQuad::smoothenNode(mc, sc.point)){
		int rank = sc.point->getRank();
		for(int i = 0; i < rank; i++){
			MeshEdge2d* edge = sc.point->getEdge(i);
			if(edge->isSideTagged() && edge != sc.mesh_edge1 && edge != sc.mesh_edge2){
				FrontEdge* fedge = front->findFrontEdge(edge, 0);
				if(fedge) front->updateFrontEdgePosition(fedge->classify(mc));
				fedge = front->findFrontEdge(edge, 1);
				if(fedge) front->updateFrontEdgePosition(fedge->classify(mc));
			}
		}
	}
	LOG_ASSERT(new_point->getEdgeToPoint(longer_point));
	// MeshingException("MG2Q:trySeam2, step3"),false); // !!!!!!


	smoothenAhead(mc, mesh, sc.point);
//	if(sc.point->getID() == 1955)
//		MeshView::showDebugMesh("smoothing 1", mesh, quad);
	MeshGenerator2dQuad::smoothenNode(mc, new_point);
	smoothenAhead(mc, mesh, new_point);
//	if(sc.point->getID() == 1955)
//		MeshView::showDebugMesh("smoothing 2", mesh, quad);
	MeshGenerator2dQuad::smoothenNode(mc, longer_point);
	smoothenAhead(mc, mesh, longer_point);
//	if(sc.point->getID() == 1955)
//		MeshView::showDebugMesh("smoothing 3", mesh, quad);
	MeshGenerator2dQuad::smoothenNode(mc, short_point);
	smoothenAhead(mc, mesh, short_point);
//	if(sc.point->getID() == 1955)
//		MeshView::showDebugMesh("smoothing 4", mesh, quad);

//	SHOW_DEBUG_MESH("Seam-2 - finished", mesh, point);

	// Re-classify front edges
	longer_edge = new_point->getEdgeToPoint(longer_point);
	LOG_ASSERT(longer_edge);
	// MeshingException("MG2Q:trySeam2, reclassify1"), false);
	if(left_longer){
		sc.edge1->setEdge(longer_edge, longer_point);
		sc.edge2->setEdge(quad_edge, new_point);
	}else{
		sc.edge1->setEdge(quad_edge, short_point);
		sc.edge2->setEdge(longer_edge, new_point);
	}
	sc.edge1->incLevel();
	sc.edge2->incLevel();
	front->updateFrontEdgePosition(sc.edge1->classify(mc));
	front->updateFrontEdgePosition(sc.edge2->classify(mc));
	FrontEdge* fedge = left_longer ? sc.edge2->getSideFrontEdge(FrontEdge::EDGE_RIGHT) : 
		sc.edge1->getSideFrontEdge(FrontEdge::EDGE_LEFT);
	front->updateFrontEdgePosition(fedge->classify(mc));

	LOG_ASSERT(front->isValid());
	// MeshingException("MG2Q:trySeam2, reclassify2"), false);

//	if(sc.point->getID() == 1955)
//		MeshView::showDebugMesh("Seam-2 done", mesh, quad);

#ifdef QUAD_DEBUG_LOG
//	if(sc.point->getID() == 1955){
	if(false){
		stringstream sstr;
		sstr << "Seam-2 closed " << sc.point->getID();
		SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true));
	}
#endif

	return true;
}

bool MeshGenerator2dQMorph::trySeam3(Metric2dContext& mc, MeshContainer2d* mesh, 
			FrontLines *front, SeamContext& sc)
{
	// (3) "Transition split operation"
	//LOG4CPLUS_INFO(MeshLog::logger_console, "Seam 3");

#ifdef QUAD_DEBUG_LOG
	if(false){
		stringstream sstr;
		sstr << "Seam-3 " << sc.point->getID() << " ["
			<< sc.mesh_edge1->getOtherPoint(sc.point)->getID() << ','
			<< sc.mesh_edge2->getOtherPoint(sc.point)->getID() << "], len="
			<< sc.length1 << ',' << sc.length2 << ']';
		SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true), 2);
	}
#endif

	bool left_longer = sc.ratio > 1.0;
	MeshEdge2d* longer_edge = left_longer ? sc.edge1->getEdge() : sc.edge2->getEdge();
	MeshEdge2d* shorter_edge = left_longer ? sc.edge2->getEdge() : sc.edge1->getEdge();

	if(longer_edge->isBorder()) return false;
		
	// New point withing quad
	MeshElement* large_element = longer_edge->getMeshElement(left_longer ? sc.point : 
		longer_edge->getOtherPoint(sc.point));
	if(!large_element || large_element->getEdgeCount() != 4 || 
			large_element->getAlphaQuality(mc) < 0.1) 
		return false;

//	SHOW_DEBUG_MESH("Seam-3 - start", mesh, large_element);
//	if(large_element->getIndex() == 162)
//		MeshView::showDebugMesh("Seam-3", mesh, large_element);
//	if(quad_conversion_step == 19480)
//		MeshView::showDebugMesh("Seam-3", mesh, large_element);

	FrontEdge* longer_side_fedge = nullptr;
	// Clear pointers to edge to-be-removed in further steps
	if(left_longer){
		longer_side_fedge = sc.edge1->getSideFrontEdge(FrontEdge::EDGE_LEFT);
		front->removeFrontEdge(sc.edge1); 
		sc.edge1 = nullptr;
	}else{
		longer_side_fedge = sc.edge2->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
		front->removeFrontEdge(sc.edge2); 
		sc.edge2 = nullptr;
	}

	// New point in the middle of the longer edge
	const DPoint2d middle = longer_edge->getPoint(0.5);
	MeshPoint2d* new_point = new MeshPoint2d(middle);
	mesh->addMeshPoint(new_point);
			
	MeshTriangle2d* large_triangle = (MeshTriangle2d*)longer_edge->getOtherElement(large_element);
	MeshPoint2d* center_point = new MeshPoint2d(large_element->getMiddlePoint());
	mesh->addMeshPoint(center_point);
	// Rebuild the mesh
	MeshPoint2d* longer_point = longer_edge->getOtherPoint(sc.point);
	large_element->switchPointsWithEdges(longer_point, center_point);
	large_triangle->switchPointsWithEdges(sc.point, new_point);
	LOG_ASSERT(!large_triangle->isInverted());
	//	MeshingException("MG2Q:trySeam3, invalid large_triangle"), false);

//	if(center_point->getID() == 137)
//		MeshView::showDebugMesh("Seam-3 xx2", mesh, large_triangle);

	MeshPoint2d* quad_point = left_longer ? large_element->getNextPoint(center_point) : 
		large_element->getPrevPoint(center_point);
	MeshQuad2d* second_quad = new MeshQuad2d(new_point, left_longer ? longer_point : center_point, 
		quad_point, left_longer ? center_point : longer_point);
	second_quad->setAreaID(sc.id);
	mesh->addMeshElement(second_quad);
	MeshPoint2d* triangle_point = left_longer ? large_triangle->getNextPoint(new_point) : 
		large_triangle->getPrevPoint(new_point);
	MeshTriangle2d* new_triangle = new MeshTriangle2d(left_longer ? new_point : sc.point, 
		left_longer ? sc.point : new_point, triangle_point);

//	if(center_point->getID() == 137)
//		MeshView::showDebugMesh("Seam-3 xx3", mesh, second_quad);

	LOG_ASSERT(!new_triangle->isInverted());
	// MeshingException("MG2Q:trySeam3, invalid new_triangle"), false);
	new_triangle->setAreaID(sc.id);
	mesh->addMeshElement(new_triangle);

	MeshTriangle2d* second_triangle = new MeshTriangle2d(
		left_longer ? sc.point : new_point, 
		left_longer ? new_point : sc.point, 
		center_point);
	second_triangle->setAreaID(sc.id);
	mesh->addMeshElement(second_triangle);

//	if(center_point->getID() == 137)
//		MeshView::showDebugMesh("Seam-3 xx4", mesh, second_quad);

	// Edge recovery
	MeshPoint2d* short_point = shorter_edge->getOtherPoint(sc.point);
	MeshEdge2d* quad_edge = short_point->getEdgeToPoint(new_point);
	if(!quad_edge){
		bool result = MeshGenerator2d::makeEdgeByEdgesSwapPIng(mc, short_point, new_point, -3);
		if(!result){
			// Try to continue desPIte the failure
#ifdef QUAD_DEBUG_LOG
			if(false){
				stringstream sstr;
				sstr << "Seam-3 failed recovery of edge: [" << short_point->getID() << ',' 
					<< new_point->getID() << ']';
				SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true), 2);
			}
#endif

			// Smoothing
			MeshGenerator2dQuad::smoothenNode(mc, center_point);
			smoothenAhead(mc, mesh, center_point);
			MeshGenerator2dQuad::smoothenNode(mc, new_point);
			smoothenAhead(mc, mesh, new_point);
			MeshGenerator2dQuad::smoothenNode(mc, longer_point);
			smoothenAhead(mc, mesh, longer_point);
			MeshGenerator2dQuad::smoothenNode(mc, short_point);
			smoothenAhead(mc, mesh, short_point);
			LOG_ASSERT(mesh->isValid());

			// .. update front
			MeshEdge2d* new_edge1 = sc.point->getEdgeToPoint(center_point);
			MeshEdge2d* new_edge2 = center_point->getEdgeToPoint(new_point);
			MeshEdge2d* long_edge1 = new_point->getEdgeToPoint(longer_point);
			LOG_ASSERT(new_edge1 && new_edge2);
			// MeshingException("MG2Q:trySeam3, bad edges"), false);
			FrontEdge *fedge1, *fedge2;
			if(left_longer){
				fedge1 = new FrontEdge(new_edge1, center_point, sc.level+3);
				fedge2 = new FrontEdge(new_edge2, new_point, sc.level+3);
				sc.edge1 = new FrontEdge(long_edge1, longer_point, sc.level+3);

				longer_side_fedge->setSideFrontEdge(sc.edge1, FrontEdge::EDGE_RIGHT);
				sc.edge1->setSideFrontEdge(longer_side_fedge, FrontEdge::EDGE_LEFT);
				sc.edge1->setSideFrontEdge(fedge2, FrontEdge::EDGE_RIGHT);
				fedge2->setSideFrontEdge(sc.edge1, FrontEdge::EDGE_LEFT);
				fedge2->setSideFrontEdge(fedge1, FrontEdge::EDGE_RIGHT);
				fedge1->setSideFrontEdge(fedge2, FrontEdge::EDGE_LEFT);
				fedge1->setSideFrontEdge(sc.edge2, FrontEdge::EDGE_RIGHT);
				sc.edge2->setSideFrontEdge(fedge1, FrontEdge::EDGE_LEFT);
				front->updateFrontEdgePosition(sc.edge2->classify(mc));
				front->addFrontEdge(sc.edge1->classify(mc));
			}else{
				fedge1 = new FrontEdge(new_edge1, sc.point, sc.level+3);
				fedge2 = new FrontEdge(new_edge2, center_point, sc.level+3);
				sc.edge2 = new FrontEdge(long_edge1, new_point, sc.level+3);

				longer_side_fedge->setSideFrontEdge(sc.edge2, FrontEdge::EDGE_LEFT);
				sc.edge2->setSideFrontEdge(longer_side_fedge, FrontEdge::EDGE_RIGHT);
				sc.edge2->setSideFrontEdge(fedge2, FrontEdge::EDGE_LEFT);
				fedge2->setSideFrontEdge(sc.edge2, FrontEdge::EDGE_RIGHT);
				fedge2->setSideFrontEdge(fedge1, FrontEdge::EDGE_LEFT);
				fedge1->setSideFrontEdge(fedge2, FrontEdge::EDGE_RIGHT);
				fedge1->setSideFrontEdge(sc.edge1, FrontEdge::EDGE_LEFT);
				sc.edge1->setSideFrontEdge(fedge1, FrontEdge::EDGE_RIGHT);
				front->updateFrontEdgePosition(sc.edge1->classify(mc));
				front->addFrontEdge(sc.edge2->classify(mc));
			}
			front->addFrontEdge(fedge1->classify(mc));
			front->addFrontEdge(fedge2->classify(mc));

			LOG_ASSERT(front->isValid());
			// MeshingException("MG2Q:trySeam3, invalid front"), false);
			return true;
		}
		quad_edge = short_point->getEdgeToPoint(new_point);
	}

//	if(center_point->getID() == 137)
//		MeshView::showDebugMesh("Seam-3 xx5", mesh, second_quad);

	MeshPoint2d* mpts[4] = {center_point, left_longer ? sc.point : new_point, 
		short_point, left_longer ? new_point: sc.point };
	DataVector<MeshEdge2d*> border_edges(4);
	for(int i = 0; i < 4; i++)
		border_edges.add(mpts[i]->getEdgeToPoint(mpts[(i+1)%4]));

	// Remove triangles from within the created quad
	MeshEdge2d* triangle_edge = sc.point->getEdgeToPoint(new_point);
	removeInnerTriangles(mesh, (MeshTriangle2d*)triangle_edge->getOtherElement(nullptr), border_edges);

	// Create new quad
	MeshQuad2d* third_quad = new MeshQuad2d(mpts[0], mpts[1], mpts[2], mpts[3]);
	third_quad->setAreaID(sc.id);
	mesh->addMeshElement(third_quad);

//	if(mpts[0]->getID() == 137 && mpts[1]->getID() == 136)
//		MeshView::showDebugMesh("Seam-3 xx6", mesh, third_quad);

	if(MeshGenerator2dQuad::smoothenNode(mc, sc.point)){
		int rank = sc.point->getRank();
		for(int i = 0; i < rank; i++){
			MeshEdge2d* edge = sc.point->getEdge(i);
			if(edge->isSideTagged() && edge != sc.mesh_edge1 && edge != sc.mesh_edge2){
				FrontEdge* fedge = front->findFrontEdge(edge, 0);
				if(fedge) front->updateFrontEdgePosition(fedge->classify(mc));
				fedge = front->findFrontEdge(edge, 1);
				if(fedge) front->updateFrontEdgePosition(fedge->classify(mc));
			}
		}
	}

//	if(mpts[0]->getID() == 137 && mpts[1]->getID() == 136)
//		MeshView::showDebugMesh("Seam-3 xx6", mesh, third_quad);

	smoothenAhead(mc, mesh, sc.point);
	MeshGenerator2dQuad::smoothenNode(mc, new_point);
	smoothenAhead(mc, mesh, new_point);
	MeshGenerator2dQuad::smoothenNode(mc, center_point);
	MeshGenerator2dQuad::smoothenNode(mc, longer_point);
	smoothenAhead(mc, mesh, longer_point);
	MeshGenerator2dQuad::smoothenNode(mc, short_point);
	smoothenAhead(mc, mesh, short_point);

//	if(large_element->getIndex() == 162)
//		MeshView::showDebugMesh("Seam-3 xx7", mesh, third_quad);

//	SHOW_DEBUG_MESH("Seam-3 - finish", mesh, point);
#ifdef QUAD_DEBUG_LOG
	if(false){
		stringstream sstr;
		sstr << "Seam-3 closed for point " << sc.point->getID();
		SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true), 1);
	}
#endif

	// Reclassify front edges
	// - remove old
	// Clear pointers to edge to-be-removed in further steps
	FrontEdge* shorter_side_fedge;
	if(left_longer){
		shorter_side_fedge = sc.edge2->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
		front->removeFrontEdge(sc.edge2); 
	}else{
		shorter_side_fedge = sc.edge1->getSideFrontEdge(FrontEdge::EDGE_LEFT);
		front->removeFrontEdge(sc.edge1); 
	}

	MeshEdge2d* new_edge1 = new_point->getEdgeToPoint(left_longer ? longer_point : short_point);
	MeshEdge2d* new_edge2 = new_point->getEdgeToPoint(left_longer ? short_point : longer_point);
	sc.edge1 = new FrontEdge(new_edge1, left_longer ? longer_point : short_point, sc.level+1);
	sc.edge2 = new FrontEdge(new_edge2, new_point, sc.level+1);

	sc.edge1->setSideFrontEdge(sc.edge2, FrontEdge::EDGE_RIGHT);
	sc.edge2->setSideFrontEdge(sc.edge1, FrontEdge::EDGE_LEFT);
	if(left_longer){
		sc.edge1->setSideFrontEdge(longer_side_fedge, FrontEdge::EDGE_LEFT);
		longer_side_fedge->setSideFrontEdge(sc.edge1, FrontEdge::EDGE_RIGHT);
		sc.edge2->setSideFrontEdge(shorter_side_fedge, FrontEdge::EDGE_RIGHT);
		shorter_side_fedge->setSideFrontEdge(sc.edge2, FrontEdge::EDGE_LEFT);
	}else{
		sc.edge1->setSideFrontEdge(shorter_side_fedge, FrontEdge::EDGE_LEFT);
		shorter_side_fedge->setSideFrontEdge(sc.edge1, FrontEdge::EDGE_RIGHT);
		sc.edge2->setSideFrontEdge(longer_side_fedge, FrontEdge::EDGE_RIGHT);
		longer_side_fedge->setSideFrontEdge(sc.edge2, FrontEdge::EDGE_LEFT);
	}

	front->addFrontEdge(sc.edge1->classify(mc));
	front->addFrontEdge(sc.edge2->classify(mc));
	front->updateFrontEdgePosition(shorter_side_fedge->classify(mc));
	front->updateFrontEdgePosition(longer_side_fedge->classify(mc));

	LOG_ASSERT(front->isValid());
	// MeshingException("MG2Q:trySeam3 - invalid front"), false);
	return true;
}

bool MeshGenerator2dQMorph::checkSeams(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines *front, 
		MeshGenerator2dQuad::QuadConversionContext * qc, int side)
{
	SeamContext sc(mc, qc, side);

	double threshold = (sc.point->isBorder() || sc.point->getElementsCount(4) > 4) ? ANGLE_1 : ANGLE_2;

	if(sc.angle < threshold){
		if(sc.ratio > 0.6 && sc.ratio < 1.7)
			return trySeam1(mc, mesh, front, sc);
		else 
			return trySeam2(mc, mesh, front, sc);
	}else if(sc.ratio < 0.4 || sc.ratio > 2.5)
		return trySeam3(mc, mesh, front, sc);
	else return false;
}
