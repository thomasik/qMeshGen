// MeshGenerator2dQuad.cpp: implementation of the MeshGenerator2dQuad class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshGenerator2dQuad.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshGenerator2d.h"
#include "MeshDomainVolume.h"
#include "MeshDomainSurface.h"
#include "FrontLines.h"
#include "MeshTriangle2d.h"
#include "MeshQuad2d.h"
#include "SurfaceParametric.h"
#include "MeshGenerator2dQMorph.h"
#include "DTriangle.h"
#include "DQuad.h"

/// Metric gradation threshold for mixed quad conversion
double MeshGenerator2dQuad::param_metric_gradation_threshold = 2.0;
/// Minimum front edge level for mixed quad conversion, where LeeLo method is allowed
int MeshGenerator2dQuad::param_leelo_minimum_level = 1;
/// Metric edge length (squared) threshold for mixed quad conversion
double MeshGenerator2dQuad::param_metric_length_threshold = 4.0;
int MeshGenerator2dQuad::param_convert_smoothing_method = 1;

bool MeshGenerator2dQuad::smoothenNode(Metric2dContext& mc, MeshPoint2d* point)
{
	if(param_convert_smoothing_method == 0)
		return MeshGenerator2d::movePointByLaplace(point);
	else if(param_convert_smoothing_method == 1)
		return MeshGenerator2d::movePointByLaplaceForVariableMetric(mc, point);
	else{
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"Unknown method for qmorph smoothing: " << param_convert_smoothing_method);
	}
	return false;
}


bool MeshGenerator2dQuad::QuadConversionContext::selectFrontEdge(const FrontLines* front)
{
	base_fedge = front->getBestFrontEdge();
	return (base_fedge != nullptr);
}

bool MeshGenerator2dQuad::QuadConversionContext::init()
{
	assert(base_fedge);
	const MeshEdge2d* medge = base_fedge->getEdge();
	int index = base_fedge->getLeftIndex();
	mp[0] = medge->getMeshPoint(index);
	mp[1] = medge->getMeshPoint(1-index);
	mp[2] = mp[3] = nullptr;
	base_element = base_fedge->getEdge()->getMeshElement(mp[0]);
	LOG_ASSERT(base_element);
	//	MeshingException("MG2QM:convert - missing base element"), false);
	LOG_ASSERT(base_element->getType() == ELEMENT_MESH_TRIANGLE);
	//	MeshingException("MG2QM:convert - invalid base element"), false);
	element_id = base_element->getAreaID();
	side_fedges[0] = side_fedges[1] = side_fedges[2] = nullptr;
	side_fedges_available[0] = side_fedges_available[1] = side_fedges_available[2] = false;

	index = base_element->getPointIndex(mp[0]);
	side_triangles[FrontEdge::EDGE_RIGHT] = (MeshTriangle2d*)base_element->getNeighbour((index+1)%3);
	side_triangles[FrontEdge::EDGE_LEFT] = (MeshTriangle2d*)base_element->getNeighbour((index+2)%3);
	side_points[FrontEdge::EDGE_UPPER] = base_element->getPoint((index+2)%3);
	side_points[FrontEdge::EDGE_LEFT] = side_triangles[FrontEdge::EDGE_LEFT] ? 
		side_triangles[FrontEdge::EDGE_LEFT]->getPrevPoint(mp[0]) : nullptr;
	side_points[FrontEdge::EDGE_RIGHT] = side_triangles[FrontEdge::EDGE_RIGHT] ? 
		side_triangles[FrontEdge::EDGE_RIGHT]->getNextPoint(mp[1]) : nullptr;

	return true;
}

void MeshGenerator2dQuad::QuadConversionContext::checkSideEdges()
{
	bool last_quad = (base_fedge->getSideFrontEdge(0)->
				getSideFrontEdge(0)->
				getSideFrontEdge(0)->
				getSideFrontEdge(0) == base_fedge);

	for(int i = 0; i < 2; i++){
		// This side - front edge already ?
		if(last_quad || base_fedge->getSideAngleCos(i) <= ANGLE_THRESHOLD){
			// Use this side front edge
			side_fedges_available[i] = true;
			side_fedges[i] = base_fedge->getSideFrontEdge(i);
			mp[3-i] = side_fedges[i]->getEdge()->getOtherPoint(mp[i]);
		}else{
			side_fedges_available[i] = false;
			side_fedges[i] = nullptr;
		}
	}
}

int MeshGenerator2dQuad::convertToQuadsMixed(Metric2dContext& mc, MeshContainer2d *mesh, 
		int merge_method, int max_ct)
{
	if(!mesh) return 0;
	int element_count = mesh->getElementsCount();
	if(element_count < 2) return 0;

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQ::convertToQuadsMixed - start.");
#endif

	START_CLOCK("MG2dQ::convertToQuadsMixed");

	mesh->setHeapOrder(false);
	int total_quad_ct = 0;

	// Remove inner nodes within edges
	//if(mesh->getInnerEdgesCount() > 0) mesh->setEdgeInnerPoints(mc, 0);

	int try_count = 3;
	while(try_count-- > 0 && mesh->getElementsCount(3) > 1){
		int point_count = mesh->getPointsCount();
		for(int i = 0; i < point_count; i++){
			MeshPoint2d* point = mesh->getPointAt(i);
			point->removeTag(TagExtended::TAG_QUAD);
			int rank = point->getRank();
			for(int j = 0; j < rank; j++)
				point->getEdge(j)->clearSideTag();
		}

		FrontLines front(element_count);
		if(gatherFrontEdges(mesh, &front) < 1) return 0;
		front.classifyAllEdges(mc);

		LOG_ASSERT(front.isValid());
		//	MeshingException("MG2Q:convertToQuadsMixed - invalid front"), 0);

		int quad_attempts = 0;
		FrontEdge* last_fedge = nullptr;
		int last_fcount = 0;

		QuadConversionContext qc;

		while(qc.selectFrontEdge(&front)){
			//++quad_conversion_step;
			LOG_ASSERT(mesh->isValid());

			if(qc.base_fedge == last_fedge){
				if(last_fcount < 2) ++last_fcount;
				else{ front.postpone(qc.base_fedge); continue; }
			}else{
				last_fedge = qc.base_fedge;
				last_fcount = 0;
			}
			if(max_ct > -1 && total_quad_ct >= max_ct) return total_quad_ct;
			if(quad_attempts++ > front.countInt()){
				if(!qc.lowerQMorphRatio()) break;
				else quad_attempts = 0;
			}
			// Create quad basing on this front edge
			bool success = qc.init();
			if(!success) return 0;
			mc.countMetricAtPoint(DPoint2d::average(qc.mp[0]->getCoordinates(), qc.mp[1]->getCoordinates()));
			const DMVector2d vt = qc.mp[1]->getMetricCoordinates(mc) - qc.mp[0]->getMetricCoordinates(mc);
			const DMVector2d vn(-vt.y, vt.x);
			const DPoint2d middle = qc.mp[0]->getCoordinates() + mc.transformMStoPS((vt+vn.normalized()) * 0.5);
			mc.countMetricAtPoint(middle);

			int prev_type = qc.base_fedge->getType();
			if(qc.base_fedge->classify(mc)->getType() - prev_type > 1){ // if fedge is worst than previously expected
				front.updateFrontEdgePosition(qc.base_fedge);
				continue;
			}

#ifdef QUAD_DEBUG_LOG
//			if(false){
			if(total_quad_ct % 10 == 0){
//			if(total_quad_ct > 20){
//				int fct = front.countInt();
//				for(int i = 0; i < fct; i++)
//					LOG4CPLUS_INFO(MeshLog::logger_mesh, i << " -> " << *front.getFrontEdge(i));

				stringstream sstr;
				sstr << "Quad " << total_quad_ct << ",  " << *(qc.base_fedge);
//				SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true));
				SHOW_MESH(sstr.str(), mesh->getViewSet(nullptr, true, true, true), 1);
			}
#endif

//			if(total_quad_ct == 809 && quad_attempts == 1){
//				SHOW_MESH("???", mesh->getViewSet(nullptr, true, true, true));
//			}
//			if(MeshPoint2d::max_id_counter >= 3078)
//				LOG4CPLUS_DEBUG(MeshLog::logger_console, "Max node count >= 3078, total_quad_ct", total_quad_ct);
//			if(MeshPoint2d::max_id_counter >= 3078 || total_quad_ct >= 2020)
//				MeshView::showDebugMesh("check front edge", mesh, qc.mp[0], qc.mp[1]);

			switch(merge_method){
			case MeshData::QUADS_QMORPH:
				success = MeshGenerator2dQMorph::createQMorphQuad(mc, mesh, &front, &qc);
				break;
			case MeshData::QUADS_LEELO:
				success = createLeeloQuad(mc, mesh, &front, &qc);
				break;
			case MeshData::QUADS_MIXED:
				if((qc.base_fedge->getLevel() < param_leelo_minimum_level ||
					qc.base_fedge->getSideFrontEdge(FrontEdge::EDGE_LEFT)->getLevel() < param_leelo_minimum_level ||
					qc.base_fedge->getSideFrontEdge(FrontEdge::EDGE_RIGHT)->getLevel() < param_leelo_minimum_level ||
					qc.base_fedge->getMetricGradation() < param_metric_gradation_threshold) &&
					qc.base_fedge->getMetricLengthSquared() < param_metric_length_threshold){
					success = MeshGenerator2dQMorph::createQMorphQuad(mc, mesh, &front, &qc);
				}else{
					qc.element_id += 2;
					success = createLeeloQuad(mc, mesh, &front, &qc);
				}
				break;
			default:
				assert(false);
			}

			LOG_ASSERT(mesh->isValid());

//			if(MeshPoint2d::max_id_counter >= 3078 || total_quad_ct >= 2020)
//				MeshView::showDebugMesh("check new quad", mesh, qc.mp[0], qc.mp[1]);

			if(success){
				//++quad_ct;
				++total_quad_ct;
				quad_attempts = 0;	// Reset failure counter
			}else{
				qc.base_fedge->incLevel(3);
				qc.base_fedge->setType(3);
				front.updateFrontEdgePosition(qc.base_fedge);
				continue;
			}

//			SHOW_STATUS(total_quad_ct);

#ifdef STAT_COUNT
			MeshData::StatData stats;
			if(MeshGenerator2d::getIncidenceInfo(mesh, stats)){
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "qmorph\t" << total_quad_ct<< "\t" << stats.average << "\t" << stats.maximum);
			}
#endif
		}
	}

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQ::convertToQuadsMixed - end.");
#endif
	STOP_CLOCK("MG2dQ::convertToQuadsMixed");
	LOG_ASSERT(mesh->isValid());
	//	MeshingException("MG2Q:convertToQuadsMixed - invalid mesh after complete conversion"), 0);

	int tct = mesh->getElementsCount(3);
	if(tct > 1){
#ifdef STAT_COUNT
		MeshRepository::addMark("MeshGenerator2dQMorph::convertToQuads - additional LeeLo.");
#endif
		total_quad_ct += MeshGenerator2dQuad::convertToQuadsMixed(mc, mesh, MeshData::QUADS_LEELO);
#ifdef STAT_COUNT
		MeshRepository::addMark("MeshGenerator2dQMorph::convertToQuads - additional LeeLo end.");
#endif
	}

//	SHOW_STEP(0, "Conversion finished successfully!.");
//	SHOW_STATUS_END(total_quad_ct);

#ifdef QUAD_DEBUG_LOG
	SHOW_MESH("Conversion finished successfully!", mesh->getViewSet(nullptr, true, true, true), 1);
#endif

	return total_quad_ct;
}

bool MeshGenerator2dQuad::checkForSingleFrontTriangle(Metric2dContext& mc, FrontLines *front, 
		MeshGenerator2dQuad::QuadConversionContext * qc)
{
	if(qc->merge_method == MeshData::QUADS_QMORPH){
		// Upper edge
		if(qc->mp[2] != qc->mp[3]) return false;
		if(!qc->side_fedges_available[0] ||!qc->side_fedges_available[1]) return false;
	}else if(qc->merge_method == MeshData::QUADS_LEELO){
		if(qc->side_triangles[FrontEdge::EDGE_RIGHT] &&
			!qc->base_element->getNextEdge(qc->base_fedge->getEdge())->isSideTagged())
			return false;
		if(qc->side_triangles[FrontEdge::EDGE_LEFT] &&
			!qc->base_element->getPrevEdge(qc->base_fedge->getEdge())->isSideTagged())
			return false;
		// else
		qc->side_fedges[FrontEdge::EDGE_RIGHT] = 
					qc->base_fedge->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
		qc->side_fedges[FrontEdge::EDGE_LEFT] = 
					qc->base_fedge->getSideFrontEdge(FrontEdge::EDGE_LEFT);
	}else{
		assert(false);
		return false;
	}

	FrontEdge* next_fleft  = qc->side_fedges[FrontEdge::EDGE_LEFT]->getSideFrontEdge(FrontEdge::EDGE_LEFT);
	FrontEdge* next_fright = qc->side_fedges[FrontEdge::EDGE_RIGHT]->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
	if(next_fleft != qc->side_fedges[FrontEdge::EDGE_RIGHT]){
		// Skip triangle in this front (not ending this front)
		next_fleft->setSideFrontEdge(next_fright, FrontEdge::EDGE_RIGHT);
		next_fright->setSideFrontEdge(next_fleft, FrontEdge::EDGE_LEFT);
		front->updateFrontEdgePosition(next_fleft->classify(mc));
		front->updateFrontEdgePosition(next_fright->classify(mc));
	}

	front->removeFrontEdge(qc->base_fedge);
	front->removeFrontEdge(qc->side_fedges[FrontEdge::EDGE_LEFT]);
	front->removeFrontEdge(qc->side_fedges[FrontEdge::EDGE_RIGHT]);
	LOG_ASSERT(front->isValid());
	//	MeshingException("MG2Q:checkForSingleFrontTriangle - error skipPIng triangle"), false);
	return true;
}

MeshTriangle2d* MeshGenerator2dQuad::selectSideTriangle(Metric2dContext& mc, MeshContainer2d* mesh, 
		FrontLines *front, QuadConversionContext * qc, int side)
{
	if(qc->side_points[side]->nonZeroIntTag(TagExtended::TAG_QUAD) && 
		!joinFrontConditionally(mc, front, qc->mp[side], qc->side_points[side]))
	{
		assert(qc->side_triangles[side]);
		const DPoint2d dpt = qc->side_triangles[side]->getMiddlePoint();
		MeshPoint2d *mpoint = new MeshPoint2d(dpt);
		mesh->addMeshPoint(mpoint);
		int id = qc->side_triangles[side]->getAreaID();
		delete mesh->removeMeshElement(qc->side_triangles[side]);
		MeshTriangle2d* t1 = (side == FrontEdge::EDGE_LEFT) ? 
			new MeshTriangle2d(qc->side_points[0], qc->mp[0], mpoint) : 
			new MeshTriangle2d(qc->mp[1], qc->side_points[1], mpoint);
		t1->setAreaID(id);
		mesh->addMeshElement(t1);
		MeshTriangle2d* t2 = (side == FrontEdge::EDGE_LEFT) ? 
			new MeshTriangle2d(qc->side_points[2], qc->side_points[0], mpoint) : 
			new MeshTriangle2d(qc->side_points[1], qc->side_points[2], mpoint);
		t2->setAreaID(id);
		mesh->addMeshElement(t2);
		MeshTriangle2d* t3 = (side == FrontEdge::EDGE_LEFT) ? 
			new MeshTriangle2d(qc->mp[0], qc->side_points[2], mpoint) :
			new MeshTriangle2d(qc->side_points[2], qc->mp[1], mpoint);
		t3->setAreaID(id);
		mesh->addMeshElement(t3);

		qc->side_points[side] = mpoint;
		qc->side_triangles[side] = t3;
	}
	return qc->side_triangles[side];
}

MeshTriangle2d* MeshGenerator2dQuad::selectBestSideTriangle(Metric2dContext& mc, 
		FrontLines *front, QuadConversionContext * qc)
{
	// Select one of two triangles
	bool left_allowed = true, right_allowed = true;
	if(qc->side_points[2]->nonZeroIntTag(TagExtended::TAG_QUAD)){
		if(!joinFrontConditionally(mc, front, qc->mp[0], qc->side_points[2])) 
			right_allowed = false;
		if(!joinFrontConditionally(mc, front, qc->mp[1], qc->side_points[2])) 
			left_allowed = false;
	}

	if(qc->side_points[FrontEdge::EDGE_LEFT]->nonZeroIntTag(TagExtended::TAG_QUAD) && 
			!joinFrontConditionally(mc, front, qc->mp[0], qc->side_points[FrontEdge::EDGE_LEFT]))
		left_allowed = false;
	if(qc->side_points[FrontEdge::EDGE_RIGHT]->nonZeroIntTag(TagExtended::TAG_QUAD) && 
			!joinFrontConditionally(mc, front, qc->mp[1], qc->side_points[FrontEdge::EDGE_RIGHT]))
		right_allowed = false;

	if(left_allowed && right_allowed){
		// Select better one
		const DMPoint2d dpt0 = qc->mp[0]->getMetricCoordinates(mc);
		const DMPoint2d dpt1 = qc->mp[1]->getMetricCoordinates(mc);
		const DMPoint2d dpt2 = qc->side_points[2]->getMetricCoordinates(mc);
		const DMPoint2d dpt_left  = qc->side_points[FrontEdge::EDGE_LEFT]->getMetricCoordinates(mc);
		const DMPoint2d dpt_right = qc->side_points[FrontEdge::EDGE_RIGHT]->getMetricCoordinates(mc);
		if(DQuad2d::alphaQuality(dpt0, dpt1, dpt2, dpt_left) > DQuad2d::alphaQuality(dpt0, dpt1, dpt_right, dpt2))
			return qc->side_triangles[FrontEdge::EDGE_LEFT];
		else return qc->side_triangles[FrontEdge::EDGE_RIGHT];
	}else if(left_allowed)
		return qc->side_triangles[FrontEdge::EDGE_LEFT];
	else if(right_allowed)
		return qc->side_triangles[FrontEdge::EDGE_RIGHT];
	else return nullptr;
}

bool MeshGenerator2dQuad::updateQuadFrontEdges(Metric2dContext& mc, FrontLines *front, 
		QuadConversionContext * qc)
{
	for(int i = 0; i < 3; i++){
		assert(qc->side_fedges[i]);
		qc->side_fedges[i]->incLevel();
	}
	// Update categories of edges
	if(qc->side_fedges_available[FrontEdge::EDGE_LEFT] && 
		qc->side_fedges_available[FrontEdge::EDGE_RIGHT] && 
		qc->side_fedges_available[FrontEdge::EDGE_UPPER])
	{
		FrontEdge* next_fleft  = qc->side_fedges[FrontEdge::EDGE_LEFT]->getSideFrontEdge(FrontEdge::EDGE_LEFT);
		FrontEdge* next_fright = qc->side_fedges[FrontEdge::EDGE_RIGHT]->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
		if(next_fleft != qc->side_fedges[FrontEdge::EDGE_UPPER]){
			next_fright =qc->side_fedges[FrontEdge::EDGE_UPPER]->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
			next_fleft->setSideFrontEdge(next_fright, FrontEdge::EDGE_RIGHT);
			next_fright->setSideFrontEdge(next_fleft, FrontEdge::EDGE_LEFT);
			front->updateFrontEdgePosition(next_fleft->classify(mc));
			front->updateFrontEdgePosition(next_fright->classify(mc));
		}else if(next_fright != qc->side_fedges[FrontEdge::EDGE_UPPER]){
			next_fleft = qc->side_fedges[FrontEdge::EDGE_UPPER]->getSideFrontEdge(FrontEdge::EDGE_LEFT);
			next_fleft->setSideFrontEdge(next_fright, FrontEdge::EDGE_RIGHT);
			next_fright->setSideFrontEdge(next_fleft, FrontEdge::EDGE_LEFT);
			front->updateFrontEdgePosition(next_fleft->classify(mc));
			front->updateFrontEdgePosition(next_fright->classify(mc));
		}
	}

	// Side edges
	for(int side = 0; side < 2; side++){
		if(qc->side_fedges_available[side]){
			if(!qc->side_fedges_available[FrontEdge::EDGE_UPPER]){
				FrontEdge* next_fside = qc->side_fedges[side]->getSideFrontEdge(side);
				next_fside->setSideFrontEdge(qc->side_fedges[FrontEdge::EDGE_UPPER], 1-side);
				qc->side_fedges[FrontEdge::EDGE_UPPER]->setSideFrontEdge(next_fside, side);
				front->updateFrontEdgePosition(next_fside->classify(mc));
			}
			front->removeFrontEdge(qc->side_fedges[side]);
		}else{
			FrontEdge* next_fside = qc->base_fedge->getSideFrontEdge(side);
			next_fside->setSideFrontEdge(qc->side_fedges[side], 1-side);
			qc->side_fedges[side]->setSideFrontEdge(next_fside, side);
			if(qc->side_fedges_available[FrontEdge::EDGE_UPPER]){
				FrontEdge* next_foside = qc->side_fedges[FrontEdge::EDGE_UPPER]->getSideFrontEdge(1-side);
				qc->side_fedges[side]->setSideFrontEdge(next_foside, 1-side);
				next_foside->setSideFrontEdge(qc->side_fedges[side], side);
				front->updateFrontEdgePosition(next_foside->classify(mc));
			}else{
				qc->side_fedges[side]->setSideFrontEdge(qc->side_fedges[FrontEdge::EDGE_UPPER], 1-side);
				qc->side_fedges[FrontEdge::EDGE_UPPER]->setSideFrontEdge(qc->side_fedges[side], side);
			}
			front->addFrontEdge(qc->side_fedges[side]->classify(mc));
			front->updateFrontEdgePosition(next_fside->classify(mc));
		}
	}

	// Upper edge
	if(qc->side_fedges_available[FrontEdge::EDGE_UPPER])
		front->removeFrontEdge(qc->side_fedges[FrontEdge::EDGE_UPPER]);
	else
		front->addFrontEdge(qc->side_fedges[FrontEdge::EDGE_UPPER]->classify(mc));

	// Bottom edge
	front->removeFrontEdge(qc->base_fedge);

	LOG_ASSERT(front->isValid());
	// MeshingException("MG2Q:updateQuadFrontLines - invalid front update"), false);
	return true;
}

bool MeshGenerator2dQuad::prepareQuadFrontEdges(FrontLines *front, QuadConversionContext * qc)
{
	qc->side_fedges[FrontEdge::EDGE_LEFT]  = qc->base_fedge->getSideFrontEdge(FrontEdge::EDGE_LEFT);
	qc->side_fedges[FrontEdge::EDGE_RIGHT] = qc->base_fedge->getSideFrontEdge(FrontEdge::EDGE_RIGHT);

	MeshEdge2d* edge = qc->mp[0]->getEdgeToPoint(qc->mp[3]);
	LOG_ASSERT(edge);
	// MeshingException("MG2Q:convertToQuad - invalid edge1"), false);
	qc->side_fedges_available[FrontEdge::EDGE_LEFT] = edge->isSideTagged(qc->mp[3]);
	if(!qc->side_fedges_available[FrontEdge::EDGE_LEFT])
		qc->side_fedges[FrontEdge::EDGE_LEFT] = new FrontEdge(edge, qc->mp[0], qc->base_fedge->getLevel());

	edge = qc->mp[1]->getEdgeToPoint(qc->mp[2]);
	LOG_ASSERT(edge);
	// MeshingException("MG2Q:convertToQuad - invalid edge2"), false);
	qc->side_fedges_available[FrontEdge::EDGE_RIGHT] = edge->isSideTagged(qc->mp[1]);
	if(!qc->side_fedges_available[FrontEdge::EDGE_RIGHT])
		qc->side_fedges[FrontEdge::EDGE_RIGHT] = new FrontEdge(edge, qc->mp[2], qc->base_fedge->getLevel());

	edge = qc->mp[2]->getEdgeToPoint(qc->mp[3]);
	LOG_ASSERT(edge);
	// MeshingException("MG2Q:convertToQuad - invalid edge3"), false);
	qc->side_fedges_available[FrontEdge::EDGE_UPPER] = edge->isSideTagged(qc->mp[2]);
	if(!qc->side_fedges_available[FrontEdge::EDGE_UPPER]) 
		qc->side_fedges[FrontEdge::EDGE_UPPER] = new FrontEdge(edge, qc->mp[3], qc->base_fedge->getLevel());
	else {
		qc->side_fedges[FrontEdge::EDGE_UPPER] = 
			front->findFrontEdge(edge, edge->getPointIndex(qc->mp[2]));
		LOG_ASSERT(qc->side_fedges[FrontEdge::EDGE_UPPER]);
		// MeshingException("MG2Q:convertToQuad - invalid upper fedge"), false);
	}
	return true;
}

bool MeshGenerator2dQuad::createLeeloQuad(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines * front, 
		QuadConversionContext * qc)
{
	qc->setMergeMethod(MeshData::QUADS_LEELO);

//	SHOW_MESH("Leelo quad - start", mesh->getViewSet(nullptr, true, true, true), 1);
	LOG_ASSERT(mesh->isValid());

	if(checkForSingleFrontTriangle(mc, front, qc)) return true;

	MeshTriangle2d* second_triangle = nullptr;

//	if(qc->mp[0]->getID() == 561 && qc->mp[1]->getID() == 475)
//		MeshView::showDebugMesh("leelo", mesh, qc->base_element);

	if(!qc->side_triangles[FrontEdge::EDGE_RIGHT] ||
		qc->base_element->getNextEdge(qc->base_fedge->getEdge())->isSideTagged())
	{
		second_triangle = selectSideTriangle(mc, mesh, front, qc, FrontEdge::EDGE_LEFT);
	}else if(!qc->side_triangles[FrontEdge::EDGE_LEFT] ||
		qc->base_element->getPrevEdge(qc->base_fedge->getEdge())->isSideTagged())
	{
		second_triangle = selectSideTriangle(mc, mesh, front, qc, FrontEdge::EDGE_RIGHT);
	}

	if(!second_triangle) second_triangle = selectBestSideTriangle(mc, front, qc);
	if(!second_triangle) return false;

	// Merge two triangles into quad
	if(second_triangle == qc->side_triangles[FrontEdge::EDGE_LEFT]){
		qc->mp[2] = qc->side_points[FrontEdge::EDGE_UPPER];
		qc->mp[3] = qc->side_points[FrontEdge::EDGE_LEFT];
	}else{
		qc->mp[2] = qc->side_points[FrontEdge::EDGE_RIGHT];
		qc->mp[3] = qc->side_points[FrontEdge::EDGE_UPPER];
	}

	delete mesh->removeMeshElement(qc->base_element);
	delete mesh->removeMeshElement(second_triangle);
	qc->base_element = new MeshQuad2d(qc->mp[0], qc->mp[1], qc->mp[2], qc->mp[3]);
	mesh->addMeshElement(qc->base_element);
	qc->base_element->setAreaID(qc->element_id);

//	if(qc->mp[0]->getID() == 561 && qc->mp[1]->getID() == 475)
//		MeshView::showDebugMesh("leelo", mesh, qc->base_element);

	if(!prepareQuadFrontEdges(front, qc)) return false;

	swapTrianglesAhead(mc, qc);

	return updateQuadFrontEdges(mc, front, qc);
}

//////////////////////////////////////////////////////////////////////
// Gather all front edges
int MeshGenerator2dQuad::gatherFrontEdges(MeshContainer2d *mesh, FrontLines *front)
{
	int fct = 0;
	int element_count = mesh->getElementsCount();
	// Gather cycles of front edges
	for(int i = 0; i < element_count; i++){
		const MeshElement* element = mesh->getElementAt(i);
		if(element->getEdgeCount() != 3) continue;
		for(int j = 0; j < 3; j++){
			const MeshEdge2d* edge = element->getEdge(j);
			// if boundary or quad-triangle edge...
			if(edge->isBorder() || edge->getOtherElement(element)->getEdgeCount() != 3){
				const MeshPoint2d* pt1 = element->getPoint(j);
				if(!edge->isSideTagged(pt1)){	// If edge wasn't tagged yet ...
					fct += gatherFrontEdges(front, element, j);
				}
			}
		}
	}
	return fct;
}

//////////////////////////////////////////////////////////////////////
// Gather set of edges forming closed part of front
//	(linked by left-right neighbourhood)
int MeshGenerator2dQuad::gatherFrontEdges(FrontLines *front, const MeshElement *element, int index)
{
	MeshEdge2d* mesh_edge = element->getEdge(index);
	MeshPoint2d* pt_left = element->getPoint(index);
	MeshPoint2d* pt_right = element->getPoint((index + 1) % 3);

	FrontEdge* start_front_edge = new FrontEdge(mesh_edge, pt_left);
	front->addFrontEdge(start_front_edge);
	FrontEdge* prev_front_edge = start_front_edge;

	const MeshElement* current_element = element;
	MeshEdge2d* current_edge = mesh_edge;
	MeshPoint2d* current_point = pt_right;
	int fct = 1;
	int loops = 0;
	while(true){
		++loops;
		// Find proper edge to right
		MeshEdge2d* next_edge = current_element->getNextEdge(current_edge);
		MeshElement* next_element = next_edge->getOtherElement(current_element);
		if(next_edge->isBorder() || next_element->getEdgeCount() != 3){
			if(!next_edge->isSideTagged(current_point)){
				MeshPoint2d* next_point = next_edge->getOtherPoint(current_point);
//				SHOW_STEP_PT_BREAKABLE(3, "Front edge found.", 
//					(current_point->getCoordinates() + next_point->getCoordinates()) * 0.5, fct);
				// Connect with incident front edges
				FrontEdge* item = new FrontEdge(next_edge, current_point);
				item->setSideFrontEdge(prev_front_edge, FrontEdge::EDGE_LEFT);
				prev_front_edge->setSideFrontEdge(item, FrontEdge::EDGE_RIGHT);
				// Add to collection of front edges
				front->addFrontEdge(item);
				++fct;
				loops = 0;
				prev_front_edge = item;
				// Initiate searching for next front edge
				//  - same element, new "rotation" point
				current_edge = next_edge;
				current_point = next_point;
			}else{
				// Reached beginning
				prev_front_edge->setSideFrontEdge(start_front_edge, FrontEdge::EDGE_RIGHT);
				start_front_edge->setSideFrontEdge(prev_front_edge, FrontEdge::EDGE_LEFT);
				break;
			}
		}else{
			// Continue searching for current front edge
			// Same point of "rotation" but new element
			current_element = next_element;
			current_edge = next_edge;
			LOG_ASSERT(loops < 100000);
			// MeshingException("MG2Q:convertToQuad - infinite loop?"), 0);
		}
	}
	return fct;
}

bool MeshGenerator2dQuad::joinFrontConditionally(Metric2dContext& mc, FrontLines *front, 
		MeshPoint2d *point1, MeshPoint2d *point2)
{
	// Whether edge exists (or can be recovered)
	MeshEdge2d* edge = point1->getEdgeToPoint(point2);
	if(!edge){
		if(!MeshGenerator2d::makeEdgeByEdgesSwapPIng(mc, point1, point2, -3)) return false;
		edge = point1->getEdgeToPoint(point2);
		LOG_ASSERT(edge);
	}
	if(edge->isSideTagged()) return true;

	// Check parity of cycles created after adding this edge to front
	// Find starting edges ...
	MeshEdge2d* left_edge1 = edge;
	while(!left_edge1->isSideTagged(left_edge1->getOtherPoint(point1))){
		const MeshElement* element = left_edge1->getMeshElement(point1);
		if(!element) return false;
		left_edge1 = element->getPrevEdge(left_edge1);
	}
	MeshEdge2d* left_edge2 = edge;
	while(!left_edge2->isSideTagged(left_edge2->getOtherPoint(point2))){
		const MeshElement* element = left_edge2->getMeshElement(point2);
		if(!element) return false;
		left_edge2 = element->getPrevEdge(left_edge2);
	}
	// ... and corresponding front edges
	FrontEdge* left_fedge1 = front->findFrontEdge(left_edge1, 1 - left_edge1->getPointIndex(point1));
	FrontEdge* left_fedge2 = front->findFrontEdge(left_edge2, 1 - left_edge2->getPointIndex(point2));
	if(!left_fedge1 || !left_fedge2)
		return false;

	// Check parity (traversing front edges)
	int level = left_fedge1->getLevel();
	int lct = 0, rct = 0;
	FrontEdge* fedge = left_fedge1;
	MeshPoint2d* left_point = nullptr;
	// Additional condition (fedge->getLeft() == left_fedge1) to skip case, 
	//	when reaching point from "other front side" (mostly when there are "free boundary lines")
	while(!(left_point == point1 && fedge->getSideFrontEdge(FrontEdge::EDGE_LEFT) == left_fedge1) && 
		left_point != point2)
	{
		fedge = fedge->getSideFrontEdge(FrontEdge::EDGE_LEFT);
		left_point = fedge->getEdge()->getMeshPoint(fedge->getLeftIndex());
		++lct;
	}
	bool one_cycle = (left_point == point2);
	fedge = left_fedge2;
	left_point = nullptr;
	while(left_point != point1 && 
		!(left_point == point2 && fedge->getSideFrontEdge(FrontEdge::EDGE_LEFT) == left_fedge2))
	{
		fedge = fedge->getSideFrontEdge(FrontEdge::EDGE_LEFT);
		left_point = fedge->getEdge()->getMeshPoint(fedge->getLeftIndex());
		++rct;
	}

	/*
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "front edges = " << front->countInt());
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "one_cycle = " << (one_cycle ? "true" : "false"));
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "lct = " << lct);
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "rct = " << rct);
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "point1 id=" << point1->getID());
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "point2 id=" << point2->getID());
	//SHOW_MESH("joining fronts", mesh->getViewSet());
	*/
	if(one_cycle){
		if(front->countInt() % 2 == 0){
			if(lct % 2 == 1 || rct % 2 == 1) return false;	// Jeden cykl dzielony nieparzyœcie ...
		}else{
			if(lct % 2 == 1 && rct % 2 == 1) return false;	// Jeden cykl dzielony nieparzyœcie ...
		}
		if((lct == 2 || rct == 2) && point1->isBorder() && point2->isBorder()) return false;
	}

	// Parity check OK - insert this edge to front
	FrontEdge* right_fedge1 = left_fedge1->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
	FrontEdge* right_fedge2 = left_fedge2->getSideFrontEdge(FrontEdge::EDGE_RIGHT);
	// Add edge at one side
	FrontEdge* fedge1 = new FrontEdge(edge, point1, level);
	left_fedge1->setSideFrontEdge(fedge1, FrontEdge::EDGE_RIGHT);
	fedge1->setSideFrontEdge(left_fedge1, FrontEdge::EDGE_LEFT);
	right_fedge2->setSideFrontEdge(fedge1, FrontEdge::EDGE_LEFT);
	fedge1->setSideFrontEdge(right_fedge2, FrontEdge::EDGE_RIGHT);
	// Add edge at the other side
	FrontEdge* fedge2 = new FrontEdge(edge, point2, level);
	left_fedge2->setSideFrontEdge(fedge2, FrontEdge::EDGE_RIGHT);
	fedge2->setSideFrontEdge(left_fedge2, FrontEdge::EDGE_LEFT);
	right_fedge1->setSideFrontEdge(fedge2, FrontEdge::EDGE_LEFT);
	fedge2->setSideFrontEdge(right_fedge1, FrontEdge::EDGE_RIGHT);

	// Update front edges
	front->addFrontEdge(fedge1->classify(mc));
	front->addFrontEdge(fedge2->classify(mc));
	front->updateFrontEdgePosition(left_fedge1->classify(mc));
	front->updateFrontEdgePosition(right_fedge1->classify(mc));
	front->updateFrontEdgePosition(left_fedge2->classify(mc));
	front->updateFrontEdgePosition(right_fedge2->classify(mc));

	return true;
}

int MeshGenerator2dQuad::convertFacesToQuads(MeshContainer3d *boundary, int method, int max_ct)
{
	START_CLOCK("MG2dQuad::convertFacesToQuads");
	int i,j,block_count = boundary->getBlocksCount();
	int finished_count = 0;
	for(i = 0; i < block_count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		assert(volume && (volume->getType() == BLOCK_DOMAIN));
		int face_count = volume->getFaceCount();
		for(j = 0; j < face_count; j++){
			MeshDomainSurface* domain_surface = (MeshDomainSurface*)volume->getFace(j);
			if(domain_surface->getMesh() != nullptr && domain_surface->convertToQuads(method, max_ct)){
				++finished_count;
			}
		}
	}
	STOP_CLOCK("MG2dQuad::convertFacesToQuads");
	return finished_count;
}

bool MeshGenerator2dQuad::smoothenFaces(MeshContainer3d *boundary, int steps, int method)
{
	START_CLOCK("MG2dQuad::smoothenFaces");
	int block_count = boundary->getBlocksCount();
	for(int i = 0; i < block_count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)boundary->getBlockAt(i);
		assert(volume && (volume->getType() == BLOCK_DOMAIN));
		int face_count = volume->getFaceCount();
		for(int j = 0; j < face_count; j++){
			MeshDomainSurface* domain_surface = (MeshDomainSurface*)volume->getFace(j);
			if(domain_surface->getMesh() != nullptr)
				domain_surface->smoothenQuads(steps, method);
		}
	}
	STOP_CLOCK("MG2dQuad::smoothenFaces");
	return true;
}

//////////////////////////////////////////////////////////////////////
// Funkcja ulepsza topologicznie lokaln¹ jakoœæ siatki dla przypadku 
//	szczególnego - poprzez usuniêcie zbêdnej krawêdzi
bool MeshGenerator2dQuad::improveQuadsByEdgeElimination(MeshContainer2d *mesh)
{
//	SHOW_STEP_BREAKABLE(1, "* Quad improvement - edge elimination.", 0);

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsByEdgeElimination - start.");
#endif

	int pct = mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point1 = mesh->getPointAt(i);
		int rank1 = point1->getRank();
		if(!point1->isBorder() && rank1 == 3 && point1->getElementsCount(3) == 0){
			MeshEdge2d* edge;
			MeshPoint2d* point2;
			int found = -1;
			for(int j = 0; j < 3; j++){
				edge = point1->getEdge(j);
				point2 = edge->getOtherPoint(point1);
				if(point2->isBorder()){
					found = -2;	// shouldn't be done near boundary ...
					break;
				}else if(point2->getRank() == 3 && point2->getElementsCount(3) == 0){
					if(found != -2 ) found = j;
				}
			}
			if(found >= 0){
				edge = point1->getEdge(found);
				point2 = edge->getOtherPoint(point1);
				MeshElement* element1 = edge->getMeshElement(point1);
				MeshElement* element2 = edge->getMeshElement(point2);
				if(element1->getEdgeCount() != 4 || element2->getEdgeCount() != 4) continue;
				MeshElement* element3 = element1->getPrevEdge(edge)->getOtherElement(element1);
				MeshElement* element4 = element1->getNextEdge(edge)->getOtherElement(element1);
				MeshPoint2d* point11 = element1->getPrevPoint(point1);
				MeshPoint2d* point12 = element1->getNextPoint(point2);
				MeshPoint2d* point21 = element2->getPrevPoint(point2);
				MeshPoint2d* point22 = element2->getNextPoint(point1);
				DPoint2d pt = edge->getPoint(0.5);
				//SHOW_STEP_PT(2, "Edge to remove.", mesh->getSurface()->getPoint(pt));
				delete mesh->removeMeshElement(element1);
				delete mesh->removeMeshElement(element2);
				if(point11->getRank() + point21->getRank() <= point12->getRank() + point22->getRank()){
					// Nowa krawêdŸ - 11->21
					element3->switchPointsWithEdges(point1, point21);
					element4->switchPointsWithEdges(point2, point11);
				}else{
					// Nowa krawêdŸ - 12->22
					element3->switchPointsWithEdges(point1, point12);
					element4->switchPointsWithEdges(point2, point22);
				}
				delete mesh->removeMeshPoint(point1->getIndex());
				delete mesh->removeMeshPoint(point2->getIndex());
				//SHOW_STEP_PT(2, "Edge removed successfully.", mesh->getSurface()->getPoint(pt));
				pct -= 2;
				continue;
			}
		}
	}

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsByEdgeElimination - start.");
#endif

	return true;
}

//////////////////////////////////////////////////////////////////////
// Funkcja ulepsza topologicznie lokaln¹ jakoœæ siatki dla przypadku 
//	szczególnego - poprzez zamianê krawêdzi granicz¹cej z dwoma
//	czworok¹tami
bool MeshGenerator2dQuad::improveQuadsByDiagonalSwapPIng(Metric2dContext& mc, MeshContainer2d *mesh)
{
//	SHOW_STEP_BREAKABLE(1, "* Smoothing quads - switching diagonal.", false);

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsByDiagonalSwapPIng - start.");
#endif

	int pct = mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point1 = mesh->getPointAt(i);
		int rank1 = point1->getRank();
		for(int j = 0; j < rank1; j++){
			MeshEdge2d* edge = point1->getEdge(j);
			MeshPoint2d* point2 = edge->getOtherPoint(point1);
			if(!point1->isBorder() && !point2->isBorder() && edge->getPointIndex(point1) == 0){
				int rank2 = point2->getRank();
				int n1 = rank1 + rank2;
				if(n1 > 8 && rank1 > 2 && rank2 > 2){
					MeshElement* element1 = edge->getMeshElement(point1);
					MeshElement* element2 = edge->getMeshElement(point2);
					if(element1->getEdgeCount() == 4 &&
						element2->getEdgeCount() == 4){
						MeshPoint2d* point11 = element1->getPrevPoint(point1);
						MeshPoint2d* point12 = element1->getPrevPoint(point11);
						MeshPoint2d* point21 = element2->getPrevPoint(point2);
						MeshPoint2d* point22 = element2->getPrevPoint(point21);
						int n2 = point11->getRank() + point21->getRank();
						int n3 = point12->getRank() + point22->getRank();
						if(point11->isBorder()) ++n2;
						if(point21->isBorder()) ++n2;
						if(point12->isBorder()) ++n3;
						if(point22->isBorder()) ++n3;
						DPoint2d pt = edge->getPoint(0.5);
						if(n1 > n2 + 2 && n3 >= n2){
							// switch diagonal 1-2 z 11-21
//							SHOW_STEP_PT(2, "Switching diagonal.", mesh->getSurface()->getPoint(pt));
							element1->switchPointsWithEdges(point1, point21);
							element2->switchPointsWithEdges(point2, point11);
							int counter = 5;
							do{
								MeshGenerator2dQuad::smoothenNode(mc, point1);
								MeshGenerator2dQuad::smoothenNode(mc, point11);
								MeshGenerator2dQuad::smoothenNode(mc, point12);
								MeshGenerator2dQuad::smoothenNode(mc, point2);
								MeshGenerator2dQuad::smoothenNode(mc, point21);
								MeshGenerator2dQuad::smoothenNode(mc, point22);
							}while(--counter && (element1->isInverted() || element2->isInverted()));
							if(element1->isInverted() || element2->isInverted()){
//								SHOW_STEP_PT(2, "Error switching diagonal.", mesh->getSurface()->getPoint(pt));
								element1->switchPointsWithEdges(point21, point1);
								element2->switchPointsWithEdges(point11, point2);
								counter = 5;
								do{
									MeshGenerator2dQuad::smoothenNode(mc, point1);
									MeshGenerator2dQuad::smoothenNode(mc, point11);
									MeshGenerator2dQuad::smoothenNode(mc, point12);
									MeshGenerator2dQuad::smoothenNode(mc, point2);
									MeshGenerator2dQuad::smoothenNode(mc, point21);
									MeshGenerator2dQuad::smoothenNode(mc, point22);
								}while(--counter && (element1->isInverted() || element2->isInverted()));
								assert(mesh->isValid());
							}
//							SHOW_STEP_PT(2, "After diagonal-switch.", mesh->getSurface()->getPoint(pt));
							break;
						}else if(n1 > n3 + 2 && n2 >= n3){
							// zamieñ przek¹tne 1-2 z 12-22
//							SHOW_STEP_PT(2, "Diagonal-switch.", mesh->getSurface()->getPoint(pt));
							element1->switchPointsWithEdges(point2, point22);
							element2->switchPointsWithEdges(point1, point12);
							int counter = 5;
							do{
								MeshGenerator2dQuad::smoothenNode(mc, point1);
								MeshGenerator2dQuad::smoothenNode(mc, point11);
								MeshGenerator2dQuad::smoothenNode(mc, point12);
								MeshGenerator2dQuad::smoothenNode(mc, point2);
								MeshGenerator2dQuad::smoothenNode(mc, point21);
								MeshGenerator2dQuad::smoothenNode(mc, point22);
							}while(--counter && (element1->isInverted() || element2->isInverted()));
							if(element1->isInverted() || element2->isInverted()){
//								SHOW_STEP_PT(2, "Error Diagonal-switch.", mesh->getSurface()->getPoint(pt));
								element1->switchPointsWithEdges(point22, point2);
								element2->switchPointsWithEdges(point12, point1);
								counter = 5;
								do{
									MeshGenerator2dQuad::smoothenNode(mc, point1);
									MeshGenerator2dQuad::smoothenNode(mc, point11);
									MeshGenerator2dQuad::smoothenNode(mc, point12);
									MeshGenerator2dQuad::smoothenNode(mc, point2);
									MeshGenerator2dQuad::smoothenNode(mc, point21);
									MeshGenerator2dQuad::smoothenNode(mc, point22);
								}while(--counter && (element1->isInverted() || element2->isInverted()));
								assert(mesh->isValid());
							}
//							SHOW_STEP_PT(2, "After Diagonal-switch.", mesh->getSurface()->getPoint(pt));
							break;
						}
					}
				}
			}
		}
	}

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsByDiagonalSwapPIng - end.");
#endif

	return true;
}

//////////////////////////////////////////////////////////////////////
// Funkcja ulepsza topologicznie lokaln¹ jakoœæ siatki dla przypadku 
//	szczególnego - poprzez usuniêcie zbêdnego czworok¹ta
bool MeshGenerator2dQuad::improveQuadsByElementElimination(Metric2dContext& mc, MeshContainer2d *mesh)
{
	//SHOW_STEP_BREAKABLE(1, "* Quad improvement - quad elimination.", 0);

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsByElementElimination - start.");
#endif

	int ect = mesh->getElementsCount();
	for(int i = 0; i < ect; i++){
		MeshQuad2d* quad = (MeshQuad2d*)mesh->getElementAt(i);
		if(quad->getEdgeCount() == 4){
			int j;
			MeshPoint2d *mp1, *mp2;
			for(j = 0; j < 2; j++){
				mp1 = quad->getPoint(j);
				mp2 = quad->getPoint(j+2);
				if(!mp1->isBorder() && mp1->getRank() == 3 && 
					!mp2->isBorder() && mp2->getRank() == 3 &&
					mp1->getElementsCount(3) == 0 && mp2->getElementsCount(3) == 0) break;
			}
			if(j < 2){
				// remove quad
				//SHOW_STEP_PT(2, "Quad to remove.", mesh->getSurface()->getPoint(quad->getMiddlePoint()));
				MeshElement* element1 = quad->getEdge((j+3)%4)->getOtherElement(quad);
				MeshElement* element2 = quad->getEdge(j)->getOtherElement(quad);
				delete mesh->removeMeshElement(quad);
				element1->switchPointsWithEdges(mp1, mp2);
				element2->switchPointsWithEdges(mp1, mp2);
				delete mesh->removeMeshPoint(mp1->getIndex());
				MeshGenerator2dQuad::smoothenNode(mc, mp2);
				//SHOW_STEP_PT(2, "Quad removed.", mesh->getSurface()->getPoint(mp2->getCoordinates()));
				--ect; --i;
				continue;
			}
			if(quad->getAlphaQuality(mc) <= 0.0){
				MeshPoint2d* mpts[4] = {quad->getPoint(0), quad->getPoint(1), 
					quad->getPoint(2), quad->getPoint(3)};
				// select proper two:
				int ind = 0;
				double vol012 = DTriangle2d::det(mpts[0]->getCoordinates(), 
					mpts[1]->getCoordinates(), mpts[2]->getCoordinates());
				double vol023 = DTriangle2d::det(mpts[0]->getCoordinates(), 
					mpts[2]->getCoordinates(), mpts[3]->getCoordinates());
				if(vol012 * vol023 > 0.0) ind = 1;
				MeshPoint2d* mpoint0 = mpts[1+ind];
				MeshPoint2d* mpoint1 = mpts[(3+ind)%4];
				if(mpoint0->isBorder() && mpoint1->isBorder()) continue;
				if(!mpoint0->isBorder() && !mpoint1->isBorder()) continue;
				if(!mpoint0->isBorder()){
					MeshPoint2d* mpoint2 = mpoint0; mpoint0 = mpoint1; mpoint1 = mpoint2;
				}
				if(!mpoint0->getEdgeToPoint(quad->getNextPoint(mpoint0))->isBorder() &&
					!mpoint0->getEdgeToPoint(quad->getPrevPoint(mpoint0))->isBorder())
					continue;
				//MeshView::showDebugMesh("bad quad for fixing", mesh, quad);
				delete mesh->removeMeshElement(quad);
				while(mpoint1->getRank() > 0){
					MeshEdge2d* edge = mpoint1->getEdge(0);
					assert(!edge->isBorder());
					MeshElement* element = edge->getMeshElement(0);
					if(element) element->switchPointsWithEdges(mpoint1, mpoint0);
					element = edge->getMeshElement(1);
					if(element) element->switchPointsWithEdges(mpoint1, mpoint0);
				}
				delete mesh->removeMeshPoint(mpoint1->getIndex());
				MeshGenerator2dQuad::smoothenNode(mc, mpoint0);
				--ect; --i;
				//MeshView::showDebugMesh("bad quad fixed", mesh, mpoint0);
			}
		}
	}

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsByElementElimination - end.");
#endif

	return true;
}

//////////////////////////////////////////////////////////////////////
// Funkcja ulepsza topologicznie lokaln¹ jakoœæ siatki dla przypadku 
//	szczególnego - poprzez usuniêcie zbêdnego wierzcho³ka
bool MeshGenerator2dQuad::improveQuadsByNodeElimination(MeshContainer2d *mesh)
{
	//SHOW_STEP_BREAKABLE(1, "* Quad improvement - node elimination.", 0);

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsByNodeElimination - start.");
#endif

	int pct = mesh->getPointsCount();
	int i = 0;
	while(i < pct){
		MeshPoint2d* point = mesh->getPointAt(i);
		if(!point->isBorder() && 
			point->getRank() == 2 &&
			point->getElementsCount(3) == 0){
			// Usun jeden z incydentnych czworok¹tów
			//SHOW_STEP_PT(2, "Node to remove.", mesh->getSurface()->getPoint(point->getCoordinates()));
			MeshEdge2d* edge = point->getEdge(0);
			MeshElement* quad1 = edge->getMeshElement(0);
			MeshElement* quad2 = edge->getMeshElement(1);
			MeshPoint2d* point2 = quad1->getNextPoint(quad1->getNextPoint(point));
			delete mesh->removeMeshElement(quad1);
			quad2->switchPointsWithEdges(point, point2);
			delete mesh->removeMeshPoint(point->getIndex());
			//SHOW_STEP_PT(2, "Merged quad.", mesh->getSurface()->getPoint(quad2->getMiddlePoint()));
			--pct;
		}else{
			++i;
		}
	}

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsByNodeElimination - end.");
#endif

	return true;
}

void MeshGenerator2dQuad::improveQuads(Metric2dContext& mc, MeshContainer2d* mesh, 
									   int steps, int method)
{
	for(int i = 0; i < steps; i++){
		if((method & MeshData::SM_TOP_SWAP) != 0){
			MeshGenerator2dQuad::improveQuadsTopological(mc, mesh);
		}
		if((method & MeshData::SM_LAPLACE) != 0){
			MeshGenerator2d::smoothenLaplace(mc, mesh, false);
		}
		if((method & MeshData::SM_LAPLACE_MIXED) != 0){
			MeshGenerator2d::smoothenLaplaceMixed(mc, mesh);
		}
	}
	MeshGenerator2d::smoothenPostCheck(mc, mesh);
}

void MeshGenerator2dQuad::improveQuadsTopological(Metric2dContext& mc, 
		MeshContainer2d* mesh, int steps)
{
	for(int i = 0; i < steps; i++){
		assert(mesh->isValid());
		improveQuadsAtBoundary(mc, mesh);
		assert(mesh->isValid());
		improveQuadsByNodeElimination(mesh);
		assert(mesh->isValid());
		improveQuadsByEdgeElimination(mesh);
		assert(mesh->isValid());
		improveQuadsByElementElimination(mc, mesh);
		assert(mesh->isValid());
		improveQuadsByDiagonalSwapPIng(mc, mesh);
		assert(mesh->isValid());
	}
}

int MeshGenerator2dQuad::improveQuadsAtBoundary(Metric2dContext& mc, MeshContainer2d *mesh)
{
//	assert(mesh->isValid());

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsAtBoundary - start.");
#endif

	int changes_count = 0;
	int ect = mesh->getElementsCount();
	for(int i = 0; i < ect; i++){
		MeshQuad2d* quad = (MeshQuad2d*)mesh->getElementAt(i);
		if(quad->getEdgeCount() == 4){
			int j;
			for(j = 0; j < 4; j++){
				if(quad->getEdge(j)->isBorder() && quad->getEdge((j+3)%4)->isBorder()){
					double angle = DVector2d::angle(
						quad->getPoint(j)->getCoordinates(),
						quad->getPoint((j+1)%4)->getCoordinates(),
						quad->getPoint((j+3)%4)->getCoordinates());
					if(angle > (5*PI/6)){
						MeshEdge2d* edge = quad->getEdge((j+2)%4);
						MeshPoint2d* pt0 = quad->getPoint(j);
						MeshPoint2d* pt1 = quad->getPoint((j+1)%4);
						MeshPoint2d* pt2 = quad->getPoint((j+2)%4);
						MeshPoint2d* pt3 = quad->getPoint((j+3)%4);
						bool left_quad = true;
						if(edge->isBorder()) {
							edge = quad->getEdge((j+1)%4);
							left_quad = false;
						}else{
							MeshEdge2d* second_edge = quad->getEdge((j+1)%4);
							if(!second_edge->isBorder()){
								if(pt1->getRank() > pt3->getRank()){
									edge = second_edge;
									left_quad = false;
								}
							}
						}
						if(edge->isBorder()) continue;
						// Insert new node 
						MeshQuad2d *quad2 = (MeshQuad2d*)edge->getOtherElement(quad);
						assert(quad2);
						if(quad2->getEdgeCount() != 4) continue;
						bool boundary_quad = false;
						for(int k = 0; k < 4; k++) 
							if(quad2->getEdge(k)->isBorder()) 
								boundary_quad = true;
//						LOG4CPLUS_INFO(MeshLog::logger_console, "Quad-at-boundary improvement");
						if(boundary_quad){
//							SHOW_STEP_PT(2, "Quad-at-boundary improvement.", mesh->getSurface()->getPoint(edge->getPoint(0.5)));
							DPoint2d dp6 = edge->getPoint(0.5);
							DPoint2d dp7 = quad2->getMiddlePoint();
							if(!quad2->isPointInside(dp7)) {
//								SHOW_STEP_PT(1, "Really bad quad !!!.", mesh->getSurface()->getPoint(dp7));
								break;
							}
							MeshPoint2d* pt6 = new MeshPoint2d(dp6);
							mesh->addMeshPoint(pt6);
							MeshPoint2d* pt7 = new MeshPoint2d(dp7);
							mesh->addMeshPoint(pt7);
							MeshPoint2d* pt5 = quad2->getPrevPoint(left_quad ? pt3 : pt2);
							MeshPoint2d* pt4 = quad2->getPrevPoint(pt5);
							int area_id = quad->getAreaID();
							delete mesh->removeMeshElement(quad);
							delete mesh->removeMeshElement(quad2);
							MeshQuad2d *new_quads[4];
							if(left_quad){
								new_quads[0] = new MeshQuad2d(pt0, pt6, pt7, pt3);
								new_quads[1] = new MeshQuad2d(pt0, pt1, pt2, pt6);
								new_quads[2] = new MeshQuad2d(pt2, pt4, pt7, pt6);
								new_quads[3] = new MeshQuad2d(pt3, pt7, pt4, pt5);
							}else{
								new_quads[0] = new MeshQuad2d(pt0, pt6, pt2, pt3);
								new_quads[1] = new MeshQuad2d(pt0, pt1, pt7, pt6);
								new_quads[2] = new MeshQuad2d(pt2, pt6, pt7, pt5);
								new_quads[3] = new MeshQuad2d(pt1, pt4, pt5, pt7);
							}
							DataVector<MeshPoint2d*> affected_points(16);
							for(int k = 0; k < 4; k++){
								mesh->addMeshElement(new_quads[k]);
								new_quads[k]->setAreaID(area_id);
								for(int ip = 0; ip < 4; ip++)
									affected_points.addIfNew(new_quads[k]->getPoint(ip));
							}
							for(size_t k = 0; k < affected_points.countInt(); k++)
								smoothenNode(mc, affected_points[k]);
//							SHOW_STEP_PT(2, "Quad-at-boundary improvement - completed.", mesh->getSurface()->getPoint(dp6));
						}else{
//							SHOW_STEP_PT(2, "Quad-at-boundary improvement [special case].", mesh->getSurface()->getPoint(edge->getPoint(0.5)));
							DPoint2d dp4 = edge->getPoint(0.5);
							MeshPoint2d* pt4 = new MeshPoint2d(dp4);
							mesh->addMeshPoint(pt4);
							MeshPoint2d* pt5 = quad2->getPrevPoint(left_quad ? pt3 : pt2);
							MeshPoint2d* pt6 = quad2->getPrevPoint(pt5);
							int area_id = quad->getAreaID();
							delete mesh->removeMeshElement(quad);
							delete mesh->removeMeshElement(quad2);
							MeshQuad2d *new_quads[3];
							if(left_quad){
								new_quads[0] = new MeshQuad2d(pt0, pt4, pt5, pt3);
								new_quads[1] = new MeshQuad2d(pt0, pt1, pt2, pt4);
								new_quads[2] = new MeshQuad2d(pt4, pt2, pt6, pt5);
							}else{
								new_quads[0] = new MeshQuad2d(pt0, pt1, pt6, pt4);
								new_quads[1] = new MeshQuad2d(pt0, pt4, pt2, pt3);
								new_quads[2] = new MeshQuad2d(pt2, pt4, pt6, pt5);
							}
							DataVector<MeshPoint2d*> affected_points(12);
							for(int k = 0; k < 3; k++){
								mesh->addMeshElement(new_quads[k]);
								new_quads[k]->setAreaID(area_id);
								for(int ip = 0; ip < 4; ip++)
									affected_points.addIfNew(new_quads[k]->getPoint(ip));
							}
							for(size_t k = 0; k < affected_points.countInt(); k++)
								smoothenNode(mc, affected_points[k]);

//							SHOW_STEP_PT(2, "Quad-at-boundary improvement - completed.", mesh->getSurface()->getPoint(dp4));
						}
//						assert(mesh->isValid());
						++changes_count;
						break;
					}
				}
			}
		}
	}

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::improveQuadsAtBoundary - end.");
#endif

	return changes_count;
}

int MeshGenerator2dQuad::splitToQuads(Metric2dContext& /* mc */, MeshContainer2d* /* mesh*/)
{
	int quad_ct = 0;

/*
	mesh->setHeapOrder(false);

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::splitToQuads - start.");
#endif

	// Ustalenie pojedynczych wêz³ów wewnêtrznych dla wszystkich elementów
	mesh->setEdgeInnerPoints(mc, 1, true);

	// Podzia³ kolejnych elementów
	int i, j, element_count = mesh->getElementsCount();
	for(i = 0; i < element_count; i++){
		MeshElement* element1 = mesh->getElementAt(i);
		int edges_ct = element1->getEdgeCount();
		MeshEdge2d* prev_edge = element1->getEdge(edges_ct - 1);
		MeshPoint2d* prev_inner_point = prev_edge->getInnerMeshPoint(0);
		// Nowy punkt wewn¹trz elementu
		DPoint2d middle_coord = element1->getMiddlePoint();
		if(edges_ct == 4 && element1->getAlphaQuality(mc) < 0.0){
			// Na wypadek wklês³ych czworok¹tów
			//  (œrodek jednej z przek¹tnych)
			middle_coord = (element1->getPoint(0)->getCoordinates() + 
					element1->getPoint(2)->getCoordinates()) * 0.5;
			if(!element1->isPointInside(middle_coord))
				middle_coord = (element1->getPoint(1)->getCoordinates() + 
						element1->getPoint(3)->getCoordinates()) * 0.5;
		}
		MeshPoint2d* middle_point = new MeshPoint2d(middle_coord);
		mesh->addMeshPoint(middle_point);
		for(j = 0; j < edges_ct; j++){
			MeshEdge2d* edge = element1->getEdge(j);
			MeshPoint2d* inner_point = edge->getInnerMeshPoint(0); //Punkt na œrodku tej krawêdzi
			inner_point->removeEdgeIfIncident(edge);
			MeshPoint2d* point = element1->getPoint(j);
			MeshQuad2d* quad = new MeshQuad2d(middle_point, prev_inner_point, point, inner_point);
			quad->setAreaID(element1->getAreaID());
			mesh->addMeshElement(quad);
			if(prev_edge->isBorder()){
				// Ustal now¹ krawêdŸ (po³ówkê) jako brzegow¹
				MeshEdge2d* new_edge = prev_inner_point->getEdgeToPoint(point);
				new_edge->setBorderType(prev_edge->getBorderType());
				Curve2dParametric* shape = prev_edge->getShape();
				if(shape){
					double t0 = prev_edge->getShapeParameter(0, point);
					double t1 = prev_edge->getShapeParameter(1, point);
					t1 = 0.5 * (t0 + t1);
					if(new_edge->getPointIndex(point) == 0){
						new_edge->setShape(shape, t0, t1);
					}else{
						new_edge->setShape(shape, t1, t0);
					}
				}
			}
			if(edge->isBorder()){
				inner_point->setBorder();
				// Ustal now¹ krawêdŸ (po³ówkê) jako brzegow¹
				MeshEdge2d* new_edge = point->getEdgeToPoint(inner_point);
				new_edge->setBorderType(edge->getBorderType());
				Curve2dParametric* shape = edge->getShape();
				if(shape){
					double t0 = edge->getShapeParameter(0, point);
					double t1 = edge->getShapeParameter(1, point);
					t1 = 0.5 * (t0 + t1);
					if(new_edge->getPointIndex(point) == 0){
						new_edge->setShape(shape, t0, t1);
					}else{
						new_edge->setShape(shape, t1, t0);
					}
				}
			}
			prev_edge = edge;
			prev_inner_point = inner_point;
		}
		for(j = 0; j < edges_ct; j++){
			MeshEdge2d* edge = element1->getEdge(j);
			if(edge->isBorder() && edge->getOtherElement(element1) == nullptr){
				// umo¿liwia usuniêcie niepotrzebnych ju¿ krawêdzi
				//	poprzez usuniêcie znacznika brzegu
				edge->setBorderType(-2);
			}
		}
		delete mesh->removeMeshElement(element1);
		quad_ct += edges_ct;
	}

	mesh->setEdgeInnerPoints(mc, 0);

#ifdef STAT_COUNT
	MeshRepository::addMark("MeshGenerator2dQuad::splitToQuads - end.");
#endif
*/

	return quad_ct;
}

/// Improve mesh before front
void MeshGenerator2dQuad::swapTrianglesAhead(Metric2dContext& mc, QuadConversionContext * qc)
{
	while(true){
		// Up
		bool swapped_triangle = MeshGenerator2d::optimizeTriangleBySwapPIng(mc, 
			(MeshTriangle2d*)qc->side_fedges[FrontEdge::EDGE_UPPER]->getEdge()->getMeshElement(qc->mp[3]));
		// Left
		swapped_triangle |= MeshGenerator2d::optimizeTriangleBySwapPIng(mc, 
			(MeshTriangle2d*)qc->side_fedges[FrontEdge::EDGE_LEFT]->getEdge()->getMeshElement(qc->mp[0]));
		// Right
		swapped_triangle |= MeshGenerator2d::optimizeTriangleBySwapPIng(mc, 
			(MeshTriangle2d*)qc->side_fedges[FrontEdge::EDGE_RIGHT]->getEdge()->getMeshElement(qc->mp[2]));
		if(!swapped_triangle) return;
	}
}
