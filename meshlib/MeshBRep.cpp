// MeshBRep.cpp: implementation of the base MeshBRep class.
//
//////////////////////////////////////////////////////////////////////

#include <memory>

#include "MeshBRep.h"
#include "MeshData.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "MeshEdge2dCurve.h"
#include "MeshPoint3d.h"
#include "MeshTriangle2d.h"
#include "MeshQuad2d.h"
#include "MeshDomainSurface.h"
#include "MeshDomainVolume.h"
#include "MeshArea.h"
#include "MeshTetrahedron.h"
#include "Curve2dAnalytic.h"
#include "Curve2dCircle.h"
#include "Curve2dBSpline.h"
#include "SurfaceParametric.h"
#include "SurfaceTranslated.h"
#include "SurfacePlane.h"
#include "SurfaceAnalytic.h"
#include "SurfaceBSplinePlanar.h"
#include "SurfaceMulti.h"
#include "MeshEdge3d.h"
#include "MeshDomainEdge3d.h"
#include "ControlSpace2dMatrixUniform.h"
#include "ControlSpace2dMesh.h"
#include "ControlSpace2dQuadTree.h"
#include "ControlSpace2dAnalytic.h"
#include "ControlSpace2dIdentity.h"
#include "ControlSpace3dAnalytic.h"
#include "ControlSpace3dMatrixUniform.h"
#include "ControlSpace3dOctree.h"
#include "ControlSpace3dInternal.h"
#include "DataHashTable.h"
#include "MeshBoundaryCondition.h"
#include "MeshViewSet.h"
#include "DPlane.h"
#include "DLeastSquaresFitting.h"
#include "MeshContainer3dSurface.h"
#include "MeshGenerator3dSurface.h"
#include "DataVector.h"

MeshBRep::~MeshBRep()
{
	if(hash_curves) delete hash_curves;
	if(hash_surfaces) delete hash_surfaces;
	if(hash_points) delete hash_points;
	if(hash_faces) delete hash_faces;
	if(hash_vertices2d) delete hash_vertices2d;
	if(m_ctable) delete m_ctable;
}

void MeshBRep::clear()
{
	point_list.clear();
	edge_list.clear();
	curve_list.clear();
	surface_list.clear();
	face_list.clear();
	block_list.clear();
	control2d_list.clear();
	control3d_list.clear();

	if(hash_curves){
		delete hash_curves;
		hash_curves = nullptr;
	}
	if(hash_surfaces){
		delete hash_surfaces;
		hash_surfaces = nullptr;
	}
	if(hash_points){
		delete hash_points;
		hash_points = nullptr;
	}
	if(hash_faces){
		delete hash_faces;
		hash_faces = nullptr;
	}
	if(hash_vertices2d){
		delete hash_vertices2d;
		hash_vertices2d = nullptr;
	}
}

bool MeshBRep::validate()
{
	// minimum point and face count
	bool valid = point_list.notEmpty() && face_list.notEmpty();
	// unless there is a block with discrete mesh
	if(!valid) 
		valid = block_list.notEmpty() && block_list.get(0)->mesh_vertices.notEmpty();
	// unless there is a surface mesh
	if(!valid)
		valid = smesh_point_list.notEmpty() && smesh_face_list.notEmpty();

	// ... finally
	if(!valid){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "At least one point and face, or mesh - required");
		return false;
	}

	// fix for 2D domains
	if( point_list.notEmpty() && block_list.countInt() == 0){
		if(surface_list.countInt() == 0 && points_2d_only){
			// prepare for 2D geometries
			auto plane = std::make_shared<SurfacePlane>(DPoint3d::zero,
					DVector3d::v_ox, DVector3d::v_oy);
			plane->setIntTag(TagExtended::TAG_ID, 0);
			plane->setFixed();
			surface_list.add(plane);
			// create one block with all faces
			auto item = std::make_shared<BlockItem>(1);
			block_list.add(item);
			int fct = face_list.countInt();
			for(int i = 0; i < fct; i++){
				auto& fitem = face_list[i];
				fitem->sid = 0;
				item->faces.addIfNew(std::abs(fitem->id));
				item->oriented.add(true);
			}
			int p_ct = point_list.countInt();
			for(int i = 0; i < p_ct; i++){
				auto& point_item = point_list[i];
				if(point_item->cid >= 0){
					point_item->sid = 0;
				}
			}
		}else{
			// create one block with all faces
			auto item = std::make_shared<BlockItem>(1);
			block_list.add(item);
			int fct = face_list.countInt();
			for(int i = 0; i < fct; i++){
				auto& fitem = face_list.get(i);
				//fitem->sid = 0;
				item->faces.addIfNew(std::abs(fitem->id));
				item->oriented.add(true);
			}
		}
	}

	return true;
}

/// Create surface-mesh from description
MeshContainer3d* MeshBRep::createSurfaceMesh()
{
	int pct = smesh_point_list.countInt();
	int fct = smesh_face_list.countInt();
	if(pct < 3 || fct < 1) return nullptr;

	MeshContainer3dSurface* surface_mesh = new MeshContainer3dSurface(pct);

	for(int i = 0; i < pct; i++)
		surface_mesh->addMeshPoint(smesh_point_list[i]);

	DataHashTableKeyValue<int, std::shared_ptr<MeshDomainVolume>> hblocks(10, -1);
	for(int i = 0; i < fct; i++){
		MeshFace* face = smesh_face_list[i];
		int bid0 = face->getIntTag(TagExtended::TAG_BLOCK_0, -1);
		if(bid0 > -1) {
			auto block = hblocks.getValue(bid0, nullptr);
			if(!block){
				block = std::make_shared<MeshDomainVolume>();
				surface_mesh->addDomainVolume(block);
				block->setAreaID(bid0);
				hblocks.insert(bid0, block);
			}
			face->setBlockLink(block.get(), 0);
		}
		int bid1 = face->getIntTag(TagExtended::TAG_BLOCK_1, -1);
		if(bid1 > -1) {
			auto block = hblocks.getValue(bid1, nullptr);
			if(!block){
				block = std::make_shared<MeshDomainVolume>();
				surface_mesh->addDomainVolume(block);
				block->setAreaID(bid1);
				hblocks.insert(bid1, block);
			}
			face->setBlockLink(block.get(), 1);
		}
		surface_mesh->addMeshFace(face);
	}

	for(int i = 0; i < smesh_surface_list.countInt(); i++)
		surface_mesh->addLocalSurface( smesh_surface_list[i] );
	for(int i = 0; i < smesh_curve_list.countInt(); i++)
		surface_mesh->addLocalCurve( smesh_curve_list[i] );

	return new MeshContainer3d(10, nullptr, surface_mesh);
}

MeshContainer3d* MeshBRep::createDomainMesh()
{
	// ... create curves
	if(!createCurves()) return nullptr;

	// ... create surfaces
	if(!createSurfaces()) return nullptr;

	// ... create points
	MeshContainer3d* mesh = createPoints3d();
	if(!mesh){
		mesh = createMeshBlocks();
		if(!mesh) mesh = createSurfaceMesh();

		return mesh;
	}

	// ... faces
	if(!createFaces()) return nullptr;

	// ... blocks
	if(!createBlocks(mesh)){
		delete mesh;
		mesh = nullptr;
	}

	clearUnusedData(mesh);

	if(mesh){

		markBoundaryTags(mesh);

		// info
		ostringstream text;
		text << "Data OK: " ;
		if(hash_points && hash_points->countInt()) text << hash_points->countInt() << " pts ";
		if(hash_points && hash_points->countInt() != (unsigned)mesh->getPointsCount())
			text << '(' << mesh->getPointsCount() << " used) ";
		if(hash_curves && hash_curves->countInt()) text << hash_curves->countInt() << " curves ";
		if(hash_surfaces && hash_surfaces->countInt()) text << hash_surfaces->countInt() << " surfaces ";
		if(hash_faces && hash_faces->countInt()) text << hash_faces->countInt()<< " faces ";
		if(block_list.countInt()) text << block_list.countInt() << " blocks";
		LOG4CPLUS_INFO(MeshLog::logger_console, text.str());

		DBox bounding_box = mesh->getBoundingBox();
		mesh_data.setModelDiameter(bounding_box.getDiameter());

		auto mbc_vector = std::make_shared<DataVector<std::shared_ptr<MeshBoundaryCondition>>>(
			mbc_list.countInt());
		mbc_list.forEach([&](auto& mbc) { mbc_vector->add(mbc); });
		mesh->setBoundaryConditions(mbc_vector);

		MeshDomainSurface::autoSetBoundaryTags(mesh);
	}

	return mesh;
}

bool MeshBRep::Control2dItem::full_analytic() const
{
	// constant and positive domain mark
	return equations[3].isConstant() && equations[3].getValue(0.0) > 0.0 &&
		equations[4].isConstant() && equations[4].getValue(0.0) > 0.0;
}

bool MeshBRep::Control3dItem::full_analytic() const
{
	// constant and positive domain mark
	return equations[6].isConstant() && equations[6].getValue(0.0) > 0.0 &&
		equations[7].isConstant() && equations[7].getValue(0.0) > 0.0;
}

bool MeshBRep::Control2dItem::parseCDS(
		const char* lx_str, const char* ly_str, const char* angle_str,
		DEquationConstTable *ctable)
{
	str_equations[0] = lx_str ? lx_str : "";
	str_equations[1] = ly_str ? ly_str : "";
	str_equations[2] = angle_str ? angle_str : "0";

	return  equations[0].parse(str_equations[0], ctable) &&
			equations[1].parse(str_equations[1], ctable) &&
			equations[2].parse(str_equations[2], ctable);
}

bool MeshBRep::Control2dItem::parseDomain(
		const char* domain1_str, const char* domain2_str,
		DEquationConstTable *ctable)
{
	str_equations[3] = domain1_str ? domain1_str : "1";
	str_equations[4] = domain2_str ? domain2_str : "1";

	return  equations[3].parse(str_equations[3], ctable) && 
			equations[4].parse(str_equations[4], ctable);
}

/*
bool MeshBRep::Control2dItem::parse(istream& is, int id)
{
	if(id == CS_DOMAIN){
		is >> str_equations[0] >> ws >> str_equations[1] >> ws >> str_equations[2];
		if(!is) return false;
		is >> ws >> str_equations[3];
		if(!is) str_equations[3] = "D=1";
		else is >> ws >> str_equations[4];
		if(!is) str_equations[4] = "D=1";
		bool valid = true;
		str_equations[0][0] = str_equations[0][1] = str_equations[0][2] = ' ';	// empty "LX="
		str_equations[1][0] = str_equations[1][1] = str_equations[1][2] = ' ';	// empty "LY="
		str_equations[2][0] = str_equations[2][1] = ' ';			// empty "A="
		str_equations[3][0] = str_equations[3][1] = ' ';		// empty "D="
		str_equations[4][0] = str_equations[4][1] = ' ';		// empty "D="
		valid = equations[0].parse(str_equations[0].c_str()+3);
		if(valid) valid = equations[1].parse(str_equations[1].c_str()+3);
		if(valid) valid = equations[2].parse(str_equations[2].c_str()+2);
		if(valid) valid = equations[3].parse(str_equations[3].c_str()+2);
		if(valid) valid = equations[4].parse(str_equations[4].c_str()+2);
		return valid;
	}
	// whether directional (i.e. metric alignment)
	char z = is.peek();
	directional = (z == 'd');
	if(directional) is >> z;

	// discrete data
	switch(id){
	case CS_POINT:	// one point
		is >> pt[0]; 
		break;
	case CS_SEGMENT:	// two points
		is >> pt[0] >> ws >> pt[1];
		break;
	case CS_CURVE:	// curve_id + [t1,t2] (as pt1)
		is >> curve_id >> pt[0];
		break;
	}
	if(!is) return false;

	is >> ws >> cds;
	if(!is) return false;

	string text;
	is >> ws >> text;
	if(!is || !DEquation::stringToDouble(text, DEquation::v_length, radius)) 
		radius = 0.0;

	return true;
}
*/

void MeshBRep::Control2dItem::adjust()
{
	ControlDataStretch2d cds(
		equations[0].getValue(0.0),
		equations[1].getValue(0.0),
		equations[2].getValue(0.0));

	// Adjust metric alignement for "constant direction" sources
	if(directional && control_type == CS_SEGMENT){ // line segment
		const DVector2d ev0 = (pt[1]-pt[0]).normalized();
		const DVector2d ev1 = ev0.orthogonal();
		double ed[] = { cds.lx, cds.ly };
		DMetric2d::adjustLengths(ed[0], ed[1]);
		data = ControlDataMatrix2d(ev0, ev1, ed);
	}else
		data = DMetric2d::stretchToMatrixWithAdjust(cds);
}

bool MeshBRep::Control3dItem::parseCDS(
		const char* lx_str, const char* ly_str, const char* lz_str, 
		const char* ax_str, const char* ay_str, const char* az_str,
		DEquationConstTable *ctable)
{
	str_equations[0] = lx_str ? lx_str : "";
	str_equations[1] = ly_str ? ly_str : "";
	str_equations[2] = lz_str ? lz_str : "";
	str_equations[3] = ax_str ? ax_str : "0";
	str_equations[4] = ay_str ? ay_str : "0";
	str_equations[5] = az_str ? az_str : "0";

	return  equations[0].parse(str_equations[0], ctable) &&
			equations[1].parse(str_equations[1], ctable) &&
			equations[2].parse(str_equations[2], ctable) &&
			equations[3].parse(str_equations[3], ctable) &&
			equations[4].parse(str_equations[4], ctable) &&
			equations[5].parse(str_equations[5], ctable);
}

bool MeshBRep::Control3dItem::parseDomain(
		const char* domain1_str, const char* domain2_str,
		DEquationConstTable *ctable)
{
	str_equations[6] = domain1_str ? domain1_str : "1";
	str_equations[7] = domain2_str ? domain2_str : "1";

	return  equations[6].parse(str_equations[6], ctable) && 
			equations[7].parse(str_equations[7], ctable);
}

/*
bool MeshBRep::Control3dItem::parse(istream& is, int id)
{
	if(id == CS_DOMAIN || id == CS_MULTI_DOMAIN){
		is >> str_equations[0] >> ws >> str_equations[1] >> ws >> str_equations[2];
		is >> ws >> str_equations[3] >> ws >> str_equations[4] >> ws >> str_equations[5];
		if(!is) return false;
		is >> ws >> str_equations[6];
		if(!is) str_equations[6] = "D=1";
		else is >> ws >> str_equations[7];
		if(!is) str_equations[7] = "D=1";
		bool valid = true;
		int clear_count[8] = {3, 3, 3, 3, 3, 3, 2, 2};
		for(int i = 0; valid && i < 8; i++){
			for(int j = 0; j < clear_count[i]; j++)
				str_equations[i][j] = ' ';	// clear "LX=" etc.
			valid = equations[i].parse(str_equations[i].c_str()+clear_count[i]);
		}
		return valid;
	}

	// whether directional (i.e. metric alignment)
	char z = is.peek();
	directional = (z == 'd');
	if(directional) is >> z;

	// discrete data
	switch(id){
	case CS_POINT:	// one point
		is >> pt[0]; 
		break;
	case CS_SEGMENT:	// two points
		is >> pt[0] >> ws >> pt[1];
		break;
	case CS_CURVE:	// curve_id + [t1,t2] (as pt1) // + surface_id or direct formulation
		LOG4CPLUS_WARN(MeshLog::logger_console, "Curvilinear metric source for 3D not implemented yet!");
		//is >> curve_id >> pt1;
		//break;
		return false;
	case CS_TRIANGLE: // three points (triangle)
		is >> pt[0] >> ws >> pt[1] >> ws >> pt[2];
		break;
	}
	if(!is) return false;

	is >> ws >> cds;
	if(!is) return false;

	string text;
	is >> ws >> text;
	if(!is || !DEquation::stringToDouble(text, DEquation::v_length, radius)) 
		radius = 0.0;

	return true;
}
*/

void MeshBRep::Control3dItem::adjust()
{
	ControlDataStretch3d cds(
		equations[0].getValue(0.0),
		equations[1].getValue(0.0),
		equations[2].getValue(0.0),
		equations[3].getValue(0.0),
		equations[4].getValue(0.0),
		equations[5].getValue(0.0));

	// Adjust metric alignement for "constant direction" sources
	if(directional && control_type == CS_SEGMENT){ // line segment
		assert(pts.countInt() > 1);
		const DVector3d ev0 = (pts[1]-pts[0]).normalized();
		DVector3d ev1, ev2;
		ev0.orthonormalVectors(ev1, ev2);
		double ed[3] = { cds.lx, cds.ly, cds.lz };
		DMetric3d::adjustLengths(ed[0], ed[1], ed[2]);
		if(ed[1] != ed[2])
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"Different lengths in segment3d metric source for arbitrary directions");
		data = ControlDataMatrix3d(ev0, ev1, ev2, ed);
	}else if(directional && control_type == CS_TRIANGLE){ // triangle face
		const DVector3d ev1 = (pts[1]-pts[0]).normalized();
		DVector3d ev2 = (pts[2]-pts[0]).normalized();
		const DVector3d ev0 = ev1.crossProduct(ev2).normalized();
		ev2 = ev0.crossProduct(ev1).normalized();
		double ed[3] = { cds.lx, cds.ly, cds.lz };
		DMetric3d::adjustLengths(ed[0], ed[1], ed[2]);
		if(ed[1] != ed[2])
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"Different lengths in triangle3d metric source for arbitrary directions");
		data = ControlDataMatrix3d(ev0, ev1, ev2, ed);
	}else
		data = DMetric3d::stretchToMatrixWithAdjust(cds);
}

bool MeshBRep::createCurves()
{
	// prepare hash table
	if(hash_curves) delete hash_curves;
	int c_ct = curve_list.countInt();
	if(c_ct > 0) hash_curves = new DataHashTableKeyValue<int, Curve2dPtr>((unsigned int)(2*c_ct), -1);
	else hash_curves = nullptr;

	// create curves
	for(int i = 0; i < c_ct; i++){
		auto curve = curve_list[i];
		assert(curve->availableTag(TagExtended::TAG_ID));
		int cid = curve->getIntTag(TagExtended::TAG_ID);
		hash_curves->insert(cid, curve);
	}
	return true;
}

bool MeshBRep::createSurfaces()
{
	// prepare hash table
	if(hash_surfaces) delete hash_surfaces;
	int s_ct = surface_list.countInt();
	if(s_ct > 0) hash_surfaces = new DataHashTableKeyValue<int, SurfacePtr>((unsigned int)(2*s_ct), -1);
	else hash_surfaces = nullptr;

	// create surfaces
	for(int i = 0; i < s_ct; i++){
		auto surface = surface_list[i];
		assert(surface->availableTag(TagExtended::TAG_ID));
		int sid = surface->getIntTag(TagExtended::TAG_ID);
		hash_surfaces->insert(sid, surface);
		////-----
		//DPoint2d pmin(0.0, 0.0);
		//DPoint2d pmax(1.0, 1.0);
		//for(int k = 0; k < corr_ct; k++){
		//	const DPoint2d& cparam = surface_item->correction_centers[k];
		//	if(cparam.x < pmin.x) pmin.x = cparam.x;
		//	if(cparam.y < pmin.y) pmin.y = cparam.y;
		//	if(cparam.x > pmax.x) pmax.x = cparam.x;
		//	if(cparam.y > pmax.y) pmax.y = cparam.y;
		//}
		//MeshViewSet* set = new MeshViewSet;
		//DPoint2d param;
		//double dd = 0.01;
		//for(param.y = pmin.y; param.y <= pmax.y; param.y += dd){
		//	param.x = pmin.x;
		//	DPoint3d pt0 = surf->getPoint(param);
		//	for(param.x = pmin.x+dd; param.x <= pmax.x; param.x += dd){
		//		DPoint3d pt1 = surf->getPoint(param);
		//		set->addEdge(pt0, pt1);
		//		pt0 = pt1;
		//	}
		//}
		//for(int k = 0; k < corr_ct; k++)
		//	set->addPoint(surf->getPoint(surface_item->correction_centers[k]));
		//SHOW_MESH("surface+correction", set);
		////-----
	}
	return true;
}

MeshContainer3d* MeshBRep::createPoints3d()
{
	// prepare hash table
	if(hash_points) delete hash_points;
	int p_ct = point_list.countInt();
	if(p_ct > 2) hash_points= new DataHashTableKeyValue<int, MeshPoint3d*>((unsigned int)(2*p_ct), -1);
	else{
		hash_points = nullptr;
		return nullptr;
	}

	// create mesh
	MeshContainer3d* mesh = new MeshContainer3d(p_ct);

	for(int i = 0; i < p_ct; i++){
		auto& point_item = point_list[i];
		if(point_item->cid >= 0){
			if(!hash_surfaces || !hash_surfaces->contains(point_item->sid)){
				LOG4CPLUS_ERROR(MeshLog::logger_console, 
					"Using undefined surface id: " << point_item->sid);
				delete mesh;
				return nullptr;
			}
			if(!hash_curves || !hash_curves->contains(point_item->cid)){
				LOG4CPLUS_ERROR(MeshLog::logger_console, 
					"Using undefined curve id: " << point_item->cid);
				delete mesh;
				return nullptr;
			}
			// -> calculate
			point_item->pt = hash_surfaces->getValue(point_item->sid)->getPoint(
				hash_curves->getValue(point_item->cid)->getPoint(point_item->pt.x));
		}else if(point_item->sid >= 0){
			if(!hash_surfaces || !hash_surfaces->contains(point_item->sid)){
				LOG4CPLUS_ERROR(MeshLog::logger_console, 
					"Using undefined surface id: " << point_item->sid);
				delete mesh;
				return nullptr;
			}
			// -> calculate
			point_item->pt = hash_surfaces->getValue(point_item->sid)->getPoint(
				DPoint2d(point_item->pt.x, point_item->pt.y));
		}
		MeshPoint3d* point = new MeshPoint3d(point_item->pt);
		hash_points->insert(point_item->id, point);
		point->setIntTag(TagExtended::TAG_ID, point_item->id);
		mesh->addMeshPoint(point);
	}

	DBox box = mesh->getBoundingBox();

	int fp_ct = freepoint_list.countInt();
	if(fp_ct > 0){
		auto fp_list = std::make_shared<DataVector<std::shared_ptr<MeshPoint3d>>>(fp_ct);
		for(int i = 0; i < fp_ct; i++){
			auto point = std::make_shared<MeshPoint3d>(freepoint_list[i]->pt);
			int bid = freepoint_list[i]->bid;
			if(bid > -1) point->setIntTag(TagExtended::TAG_BID, bid);
			box.addPoint(freepoint_list[i]->pt);
			// boundary conditions for free-point
			if(freepoint_list[i]->label != ""){
				auto bc = std::make_shared<MeshBoundaryCondition>(freepoint_list[i]->label);
				mbc_list.add(bc);
				point->setPtrTag(TagExtended::TAG_BOUNDARY_COND, bc.get());
			}
			fp_list->add(point);
		}
		mesh->setFreePoints(fp_list);
	}

	// set initial diameter-value, for now only from vertices (and freepoints)
	mesh_data.setModelDiameter(box.getDiameter());

	return mesh;
}

bool MeshBRep::createFaces()
{
	// prepare hash table
	if(hash_faces) delete hash_faces;
	int f_ct = face_list.countInt();
	if(f_ct > 0) hash_faces = new DataHashTableKeyValue<int, MeshDomainSurface*>(2*f_ct, -1);
	else{
		hash_faces = nullptr;
		return false;
	}

	// create faces
	double max_dist = 0.0;
	for(int i = 0; i < f_ct; i++){
		auto& face_item = face_list[i];
		//---
		if(face_item->sid > 0 && (!hash_surfaces || !hash_surfaces->contains(face_item->sid)))
		{
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"Face with non-existing surface: " << face_item->sid);
			return false;
		}
		//---
		SurfacePtr surface;
		if(hash_surfaces) surface = hash_surfaces->getValue(face_item->sid, nullptr);
		if(!surface){
			DataVector<DPoint3d> all_points(face_item->areas.countInt());
			for(int j = 0; j < face_item->areas.countInt(); j++)
				for(int k = 0; k < face_item->areas[j].fpoints.countInt(); k++)
					all_points.add(hash_points->getValue(face_item->areas[j].fpoints[k].pid, nullptr)->getCoordinates());
			DPlane plane;
			double dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(all_points, plane, true);
			if(dist < 0.0){
				LOG4CPLUS_ERROR(MeshLog::logger_console, "Error creating plane for face: " << i);
				return false;
			}
			if(dist > max_dist) max_dist = dist;
			surface = std::make_shared<SurfacePlane>(plane);
			surface->setFixed();
		}
		//---
		auto mds = createFace(face_item.get(), surface);
		if(!mds) return false;
		hash_faces->insert(face_item->id, mds);

		// boundary conditions for face
		if(face_item->tag != ""){
			auto bc = std::make_shared<MeshBoundaryCondition>(face_item->tag);
			mbc_list.add(bc);
			mds->setPtrTag(TagExtended::TAG_BOUNDARY_COND, bc.get());
		}
	}
	if(max_dist > mesh_data.relative_small_number)
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Maximum diff for created plane from vertices: " << max_dist);
	return true;
}

/// Create the i-th face
MeshDomainSurface* MeshBRep::createFace(FaceItem* item, SurfaceConstPtr surface)
{
	// prepare 2d vertices for face
	if(!prepareFaceVertices(item, surface)) return nullptr;
	// prepare edges (non-linear or special) for face
	if(!prepareFaceEdges(item)) return nullptr;

	// create 2d mesh with boundary entities for face
	MeshContainer2d* boundary = createFaceBoundaryMesh2d(item);
	if(boundary) boundary->setSurface(surface);
	else return nullptr;

	// prepare control space
	CS2dPtr space = prepareUserCS(item, boundary);

	// create/gather/link 3d entities
	MeshDomainSurface* mds = createFace3D(boundary);
	if(mds){
		mds->setIntTag(TagExtended::TAG_ID, item->id);
		mds->setUserControlSpace(space);
	}

	return mds;
}

/// Prepare vertices for face
bool MeshBRep::prepareFaceVertices(FaceItem* item, SurfaceConstPtr surface)
{
	// count multi-faces for this fid
	int m_ct = item->areas.countInt();
	int multi_total = 0;
	for(int j = 0; j < m_ct; j++){
		multi_total += item->areas[j].fpoints.countInt();
	}

	// prepare hash table
	if(hash_vertices2d) delete hash_vertices2d;
	hash_vertices2d = new DataHashTableKeyValue<int, MeshPoint2d*>((unsigned int)(2*multi_total), -1);

	// calculate, create and insert vertices into hash table
	for(int j = 0; j < m_ct; j++){
		FaceAreaItem& area = item->areas[j];

		int mf_ct = area.fpoints.countInt();
		for(int k = 0; k < mf_ct; k++){
			FacePoint& fpt = area.fpoints[k];
			MeshPoint2d* mpt2d = hash_vertices2d->getValue(fpt.pid, nullptr);
			if(!mpt2d){
				MeshPoint3d* mpt = hash_points->getValue(fpt.pid, nullptr);
				if(fpt.params_available){
					if(fpt.cid >= 0) // from curve
						fpt.params = hash_curves->getValue(fpt.cid)->getPoint(fpt.params.x);
				}else{
					// from surface (not very efficient usually, unless it's a plane...)
					fpt.params = surface->getParameters(mpt->getCoordinates());
				}
				// check
				double dist = surface->getPoint(fpt.params).distance(mpt->getCoordinates());
				if(dist > mesh_data.relative_small_number){
					LOG4CPLUS_WARN(MeshLog::logger_mesh, 
						"Inconsistent surface point coordinates for face #" <<  item->id
						<< ", point #" <<  fpt.pid << ", distance = " <<  dist);
				}
				// create
				mpt2d = new MeshPoint2d(fpt.params);
				mpt2d->setBorder(TagBorder::OUTER | TagBorder::CORNER | TagBorder::FIXED);
				mpt2d->setIntTag(TagExtended::TAG_ID, fpt.pid);
				hash_vertices2d->insert(fpt.pid, mpt2d);
				mpt2d->setPtrTag(TagExtended::TAG_MP_2D_3D, mpt);
			}else{
				fpt.cid = -1;
				fpt.params = mpt2d->getCoordinates();
			}
		}
	}
	return true;
}

/// Prepare (special, i.e. curvilinear) edges for face
bool MeshBRep::prepareFaceEdges(FaceItem* item)
{
	if(!hash_vertices2d) return false;

	for(int i = 0; i < edge_list.countInt(); i++){
		auto& edge_item = edge_list[i];
		if(edge_item->sid != item->sid && edge_item->sid >= 0) continue;
		MeshPoint2d* mpt0 = hash_vertices2d->getValue(edge_item->vid0, nullptr);
		MeshPoint2d* mpt1 = hash_vertices2d->getValue(edge_item->vid1, nullptr);
		if(!mpt0 || !mpt1) continue;
		MeshEdge2d* medge = mpt0->getEdgeToPoint(mpt1);
		if(!medge){
			// create following the prescription
			if(edge_item->cid >= 0){
				auto shape = hash_curves->getValue(edge_item->cid, nullptr);
				if(!shape){
					LOG4CPLUS_ERROR(MeshLog::logger_console, 
						"Nonlinear edge with non-existing curve id: " << edge_item->cid);
					return false;
				}
				medge = new MeshEdge2dCurve(mpt0, mpt1, shape, 
					edge_item->t0, edge_item->t1);
			}else{
				medge = new MeshEdge2d(mpt0, mpt1);
			}
		}else{
			// check with the prescription ?
		}
		if (edge_item->fixed) {
			medge->setIntTag(TagExtended::TAG_FIXED);
			medge->setBorderFlags(TagBorder::FIXED);
		}
		// boundary conditions for edge
		if(edge_item->tag != ""){
			auto bc = std::make_shared<MeshBoundaryCondition>(edge_item->tag);
			mbc_list.add(bc);
			medge->setPtrTag(TagExtended::TAG_BOUNDARY_COND, bc.get());
		}
	}
	return true;
}

/// Create 2d area-boundary-mesh description for face in parametric space
MeshContainer2d* MeshBRep::createFaceBoundaryMesh2d(FaceItem* item)
{
	if(!hash_vertices2d || hash_vertices2d->countInt() < 3) return nullptr;

	// mesh entities from all sub-faces
	MeshContainer2d* boundary = new MeshContainer2d((int)hash_vertices2d->countInt());

	DataHashTable<MeshPoint2d*> temp_hash_points2d((unsigned int)(2*hash_vertices2d->countInt()), nullptr);
	for(int j = 0; j < item->areas.countInt(); j++){
		FaceAreaItem& area = item->areas[j];
		int mf_ct = area.fpoints.countInt();
		DataVector<DPoint2d> temp_params(mf_ct);
		DataVector<MeshPoint2d*> temp_mesh_points2d(mf_ct);
		// process points
		for(int k = 0; k < mf_ct; k++){
			FacePoint& fpt = area.fpoints[k];
			MeshPoint2d* mpt2d = hash_vertices2d->getValue(fpt.pid, nullptr);
			LOG_ASSERT(mpt2d); // has to be already created
			if(temp_hash_points2d.insert(mpt2d)) // if not already there
				boundary->addMeshPoint(mpt2d);
			temp_params.add(fpt.params);
			temp_mesh_points2d.add(mpt2d);
		}
		// check orientation
		bool wrong_orientation = false;
		if(area.closed) wrong_orientation = MeshArea::getOrientation(temp_mesh_points2d) != area.filled;
		if(wrong_orientation){
			LOG4CPLUS_WARN(MeshLog::logger_mesh, 
				item->id << "(" << j << "): Inconsistent orientation");
			// Reverse points-sequence
			temp_mesh_points2d.reverse();
			if(MeshArea::getOrientation(temp_mesh_points2d) != area.filled){
				MeshViewSet* set = new MeshViewSet;
				for(int k = 0; k < temp_mesh_points2d.countInt(); k++){
					const DPoint2d& pt = temp_mesh_points2d[k]->getCoordinates();
					set->addPoint(DPoint3d(pt.x, pt.y, 0.0), 0, 
						temp_mesh_points2d[k]->getIntTag(TagExtended::TAG_ID, 0));
				}
				SHOW_MESH("Inconsistency problem", set);
			}
			assert(MeshArea::getOrientation(temp_mesh_points2d) == area.filled);
		}

		// process edges
		for(int k = 0; k < mf_ct; k++){
			int next_k = (k+1) % mf_ct;
			if(!area.closed && next_k == 0) break; // Don't add last segment for free-lines
			MeshEdge2d* edge = temp_mesh_points2d[k]->getEdgeToPoint(temp_mesh_points2d[next_k]);
			if(!edge) edge = new MeshEdge2d(temp_mesh_points2d[k], temp_mesh_points2d[next_k]);
			if(!area.closed)
				edge->setBorder(TagBorder::INNER); // "inner boundary"
		}

		if(area.closed){
			MeshArea* marea = new MeshArea(temp_mesh_points2d, area.filled);
			marea->setAreaID(area.filled ? area.sub_id : -1);
			boundary->addMeshElement(marea);
		}
	}
	return boundary;
}

/// Prepare user CS for face
CS2dPtr MeshBRep::prepareUserCS(FaceItem* item, MeshContainer2d* boundary)
{
	// control data
	int c_ct = control2d_list.countInt();
	if(c_ct == 0) return nullptr;
	int control_count = 0;
	for(int j = 0; j < c_ct; j++){
		if(control2d_list[j]->fid == item->id) ++control_count;
	}
	if(control_count == 0) return nullptr;

	DataVector<std::shared_ptr<const Control2dItem>> analytic;
	DataVector<std::shared_ptr<const Control2dItem>> part_analytic;
	DataVector<std::shared_ptr<const Control2dItem>> full_analytic;
	DataVector<std::shared_ptr<Control2dItem>> discrete;
	DataVector<std::shared_ptr<const Control2dItem>> multi_part_analytic;
	DataVector<std::shared_ptr<const Control2dItem>> multi_full_analytic;

	for(int j = 0; j < c_ct; j++){
		auto control_item = control2d_list[j];
		if(control_item->fid != item->id) continue;

		switch(control_item->control_type){
			case Control2dItem::CS_POINT:
			case Control2dItem::CS_SEGMENT:
			case Control2dItem::CS_CURVE:
				discrete.add(control_item); break;
			case Control2dItem::CS_DOMAIN:
				if(control_item->multi_nr == 0){ // i.e. no multi-domain
					analytic.add(control_item);
					if(control_item->full_analytic())
						full_analytic.add(control_item);
					else
						part_analytic.add(control_item);
				}else{
					if(control_item->full_analytic())
						multi_full_analytic.add(control_item);
					else
						multi_part_analytic.add(control_item);
				}
				break;
		}
	}

	if(multi_full_analytic.countInt() > 0 || multi_part_analytic.countInt() > 0){
		if(multi_full_analytic.countInt() != 1){
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"Only one full multi-analytic equation allowed!");
			return nullptr;
		}
		if(analytic.countInt() > 0 || discrete.countInt() > 0){
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"Multi-analytic equations can't be mixed with other sources!");
			return nullptr;
		}

		// one full analytic CS + several domain analytic CS
		auto aspace = std::make_shared<ControlSpace2dAnalytic>(
			boundary->getSurface(),
			multi_full_analytic.get(0)->str_equations[0].c_str(), 
			multi_full_analytic.get(0)->str_equations[1].c_str(), 
			multi_full_analytic.get(0)->str_equations[2].c_str());

		for(int j = multi_part_analytic.countInt()-1; j >=0 ; j--){
			aspace->addCaseAsFirst(
				multi_part_analytic.get(j)->str_equations[3].c_str(),
				multi_part_analytic.get(j)->str_equations[4].c_str(),
				multi_part_analytic.get(j)->str_equations[0].c_str(), 
				multi_part_analytic.get(j)->str_equations[1].c_str(), 
				multi_part_analytic.get(j)->str_equations[2].c_str());
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
				"multi-part as first: " 
				<< multi_part_analytic.get(j)->str_equations[3] << " and " 
				<< multi_part_analytic.get(j)->str_equations[4]);
		}
		return aspace;
	}else if(full_analytic.countInt() > 0 && part_analytic.countInt() == 0 && discrete.countInt() == 0){
		// one full analytic CS only
		return std::make_shared<ControlSpace2dAnalytic>(
			boundary->getSurface(),
			full_analytic.get(0)->str_equations[0].c_str(), 
			full_analytic.get(0)->str_equations[1].c_str(), 
			full_analytic.get(0)->str_equations[2].c_str());
	}

	// else: analytic + discrete
	DRect rect = boundary->getBoundingRect();
	rect.inflate(ControlSpace2d::param_inflate_box_factor);
	CS2dAPtr space;

	switch(ControlSpace2d::param_control_type){
	case MeshData::CONTROL_UNIFORM:
		space = std::make_shared<ControlSpace2dMatrixUniform>(boundary->getSurface(), rect,
			ControlSpace2dMatrixUniform::param_uniform_nx, 
			ControlSpace2dMatrixUniform::param_uniform_nx);
		break;
	case MeshData::CONTROL_MESH:
		space = std::make_shared<ControlSpace2dMesh>(boundary->getSurface(), rect,
			ControlSpace2dMesh::param_interpolation_method);
		break;
	case MeshData::CONTROL_QUADTREE:
		space = std::make_shared<ControlSpace2dQuadTree>(boundary->getSurface(), rect);
		break;
	default:
		LOG_ASSERT(false);
		break;
	}
	LOG_ASSERT(space);

	// 1. initialize the default CS-nodes with domain-analytic sources
	for(int m = 0; m < analytic.countInt(); m++){
		auto control_item = analytic[m];
		space->forEachControlNode([&](ControlNode2d& node) {
			if (node.w >= 0.0 && // if not already initialized
				control_item->equations[3].getValue(node.coord.x, node.coord.y) >= 0.0) // if within the equation domain
			{ 
				ControlDataStretch2d cds(
					control_item->equations[0].getValue(node.coord.x, node.coord.y),
					control_item->equations[1].getValue(node.coord.x, node.coord.y),
					control_item->equations[2].getValue(node.coord.x, node.coord.y));
				if (node.w == 0.0) {
					node.control_data = DMetric2d::stretchToMatrixWithAdjust(cds);
					node.w = -1.0;
				}
				else {
					node.control_data /= node.w;
					node.control_data.setMinimum(DMetric2d::stretchToMatrixWithAdjust(cds));
					node.w = -1.0;
				}
			}
		});
	}

	// 2. insert discrete sources
	if(full_analytic.countInt() > 0){
		// ... updating the existing CS
		space->interpolate();
		for(int m = 0; m < discrete.countInt(); m++){
			auto control_item = discrete[m];
			control_item->adjust();
			switch(control_item->control_type){
			case Control2dItem::CS_POINT:
				space->setMinControl(
					ControlDataExtMatrix2dRadial(control_item->pt[0], 
						control_item->radius, control_item->data, 
						boundary->getSurface()));
				break;
			case Control2dItem::CS_SEGMENT:
				space->setMinControl(
					ControlDataExtMatrix2dSegment(control_item->pt[0], 
						control_item->pt[1] - control_item->pt[0], 
						control_item->radius, control_item->data, 
						boundary->getSurface()));
				break;
			case Control2dItem::CS_CURVE:
				{
					assert(hash_curves->contains(control_item->curve_id));
					auto curve = hash_curves->getValue(control_item->curve_id);
					curve->setIntTag( TagExtended::TAG_ACTIVE);
					ControlDataExtMatrix2dCurve cs_emc( curve, 
							control_item->pt[0].x, control_item->pt[0].y, 
							control_item->radius, control_item->data,
							boundary->getSurface());
					space->setMinControl(cs_emc);
					break;
				}
			default:
				LOG_ASSERT(false);
			}
		}
	}else{
		// ... initializing the CS
		for(int m = 0; m < discrete.countInt(); m++){
			auto control_item = discrete[m];
			control_item->adjust();
			switch(control_item->control_type){
			case Control2dItem::CS_POINT:
				space->addControlPoint(control_item->pt[0], 
					control_item->data, control_item->radius);
				break;
			case Control2dItem::CS_SEGMENT:
				space->addControlSegment(control_item->pt[0], control_item->pt[1], 
					control_item->data, control_item->radius);
				break;
			case Control2dItem::CS_CURVE:
				{
					auto shape = hash_curves->getValue(control_item->curve_id);
					LOG_ASSERT(shape);
					shape->setIntTag(TagExtended::TAG_ACTIVE);
					space->addControlSegment(shape, 
						control_item->pt[0].x, control_item->pt[0].y, 
						control_item->data, control_item->radius);
					break;
				}
			default:
				LOG_ASSERT(false);
			}
		}
		space->interpolate();
	}

	// 3. final touch -> adjust newly added control nodes for domain-analytic sources
	for(int m = 0; m < analytic.countInt(); m++){
		auto control_item = analytic[m];
		space->forEachControlNode([&](ControlNode2d& node) {
			if (control_item->equations[3].getValue(node.coord.x, node.coord.y) >= 0.0) { // within the equation domain
				ControlDataStretch2d cds(
					control_item->equations[0].getValue(node.coord.x, node.coord.y),
					control_item->equations[1].getValue(node.coord.x, node.coord.y),
					control_item->equations[2].getValue(node.coord.x, node.coord.y));
				assert(node.w < 0.0);
				node.control_data.setMinimum(DMetric2d::stretchToMatrixWithAdjust(cds));
			}
		});
	}

	return space;
}

/// Create domain-surface for face 3d representation
MeshDomainSurface* MeshBRep::createFace3D(MeshContainer2d* boundary)
{
	int pct = boundary->getPointsCount();

	// gather points
	DataVector<MeshPoint3d*> points(pct);
	for(int i = 0; i < pct; i++){
		MeshPoint2d* mpt2d = boundary->getPointAt(i);
		int id = mpt2d->getIntTag(TagExtended::TAG_ID);
		MeshPoint3d* mpt = hash_points->getValue(id, nullptr);
		assert(mpt);
		points.add(mpt);
		mpt2d->setPtrTag(TagExtended::TAG_MP_2D_3D, mpt);
	}

	// gather edges
	DataVector<MeshDomainEdge3d*> edges(2*pct);
	for(int i = 0; i < pct; i++){
		// ... first point
		MeshPoint2d* mpt2d0 = boundary->getPointAt(i);
		MeshPoint3d* mpt0 = points[i];
		assert(mpt2d0->getPtrTag(TagExtended::TAG_MP_2D_3D) == mpt0);

		int rank = mpt2d0->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* me2d = mpt2d0->getEdge(j);
			if(me2d->getMeshPoint(0) != mpt2d0) continue;
			// ... second point
			MeshPoint2d* mpt2d1 = me2d->getMeshPoint(1);
			MeshPoint3d* mpt1 = points[mpt2d1->getIndex()];
			assert(mpt2d1->getPtrTag(TagExtended::TAG_MP_2D_3D) == mpt1);
			// ... check edge
			MeshDomainEdge3d* me = (MeshDomainEdge3d*) mpt0->getEdgeToPoint(mpt1);
			if(me){
				assert(me->getType() == EDGE_DOMAIN_3D);
				edges.add(me);
			}else{
				edges.add(me = new MeshDomainEdge3d(mpt0, mpt1));
			}
			me2d->setPtrTag(TagExtended::TAG_ME_2D_3D, me);
			if(me2d->availableTag(TagExtended::TAG_FIXED))
				me->setIntTag(TagExtended::TAG_FIXED);
			if(me2d->availableTag(TagExtended::TAG_BOUNDARY_COND))
				me->setPtrTag(TagExtended::TAG_BOUNDARY_COND,
					me2d->getPtrTag(TagExtended::TAG_BOUNDARY_COND));
		}
	}

	// create domain surface
	return new MeshDomainSurface(points, edges, boundary->getSurface(), boundary);
}

/// Create mesh blocks from description
MeshContainer3d* MeshBRep::createMeshBlocks()
{
	int b_ct = block_list.countInt();
	if(b_ct == 0) return nullptr;

	MeshContainer3d* mesh = new MeshContainer3d(b_ct);
	DBox box;

	for(int i = 0; i < b_ct; i++){
		auto& item = block_list[i];
		if(item->mesh_vertices.empty()) continue;

		int pct = item->mesh_vertices.countInt();
		MeshContainer3d* mesh3d = new MeshContainer3d(pct);
		for (auto itp = item->mesh_vertices.iterator(); itp.valid(); itp.moveNext()) {
			mesh3d->addMeshPoint(new MeshPoint3d(itp.item()));
		}
		for (auto ite = item->mesh_elements.iterator(); ite.valid(); ite.moveNext()) {
			int vt0 = ite.item(); ite.moveNext();
			int vt1 = ite.item(); ite.moveNext();
			int vt2 = ite.item(); ite.moveNext();
			int vt3 = ite.item(); 
			MeshBlock* block = new MeshTetrahedron(
				mesh3d->getPointAt(vt0),
				mesh3d->getPointAt(vt1),
				mesh3d->getPointAt(vt2),
				mesh3d->getPointAt(vt3));
			mesh3d->addMeshBlock(block);
			block->setAreaID(item->id);
		}
		mesh3d->markOuterBoundary();
		MeshDomainVolume* mdv = new MeshDomainVolume(mesh3d);
		mesh->addMeshBlock(mdv);
		box.addBox(mesh3d->getBoundingBox());
	}

	// set initial diameter-value, for now only from vertices (and freepoints)
	mesh_data.setModelDiameter(box.getDiameter());

	return mesh;
}

bool MeshBRep::createBlocks(MeshContainer3d* mesh)
{
	if(!hash_faces) return false;
	int b_ct = block_list.countInt();

	for(int i = 0; i < b_ct; i++){
		auto& item = block_list[i];

		// -> gather faces
		int bfct = item->faces.countInt();
		assert(bfct > 0);
		DataVector<MeshFace*> mdv_faces(bfct);
		DataVector<bool> mdv_orientation(bfct);
		// -> ... and points
		int total_bpct = 0;
		for(int j = 0; j < bfct; j++){
			MeshDomainSurface* mds = hash_faces->getValue(item->faces[j], nullptr);
			if(!mds){
				LOG4CPLUS_ERROR(MeshLog::logger_console, 
					"Error for block id: " << item->id 
					<< ", used non-existing face id: " << item->faces[j]); 
				return false;
			}
			total_bpct += mds->getPointCount();
		}
		DataHashTable<MeshPoint3d*> hash_mdv_points(total_bpct, nullptr);

		DataVector<MeshPoint3d*> mdv_points(total_bpct);

		int ifct = item->faces.countInt();
		for(int k = 0; k < ifct; k++){
			MeshFace* mv_face = hash_faces->getValue(item->faces[k], nullptr);
			assert(mv_face);

			mdv_faces.add(mv_face);
			mdv_orientation.add(item->oriented[k]);

			int vct = mv_face->getPointCount();
			for(int m = 0; m < vct; m++){
				MeshPoint3d* point = mv_face->getPoint(m);
				if(hash_mdv_points.insert(point)){
					mdv_points.add(point);
				}
			}
		}
		//---
		MeshDomainVolume* volume = new MeshDomainVolume(mdv_faces, mdv_points, mdv_orientation);
		volume->setAreaID(item->id);
/*
		if(block_item->prism_vector.length() > 0.0){
			// block-prism definition from a face with vector
			if(fct > 1) LOG4CPLUS_WARN(MeshLog::logger_console, "Prism block with more than one face.");
			// create additional faces.
			int fpct = marked_points.countInt();
			MeshDomainSurface* face = (MeshDomainSurface*)faces[0];
			// + opposite face (with full copy of points)
			SurfaceParametric* surface = new SurfaceTranslated(
				face->getBaseSurface(), block_item->prism_vector);
			surface->incRefCount(); 
			//---
			//MeshContainer2d* boundary = new MeshContainer2d(multi_total);
			//DataVector<MeshPoint3d*> points(multi_total);
			//DataVector<MeshDomainEdge3d*> edges(multi_total);
			//DataVector<MeshPoint2d*> mesh_points2d(multi_total);
			for(int j = 0; j < fpct; j++){
				MeshPoint3d* pt = new MeshPoint3d(
					marked_points[j]->getCoordinates() + block_item->prism_vector);
				marked_points.add(pt);
			}
			//MeshDomainSurface* mds_opp = new MeshDomainSurface(points, edges, surface, boundary);
			//mds_opp->setUserControlSpace(space);
			//faces.add(mds_opp);
			//valid_orientation.add(!valid_orientation[0]);
		}else{
			// simple block definition from a set of faces
			volume = new MeshDomainVolume(faces, valid_orientation, marked_points);
		}
*/
		// control data
		volume->setUserControlSpace(prepareUserCS(item.get(), volume));

		mesh->addMeshBlock(volume);
	}

	return true;
}

/// Prepare (user) control space for block
CS3dPtr MeshBRep::prepareUserCS(BlockItem* item, MeshDomainVolume* volume)
{
	int c_ct = control3d_list.countInt();
	if(c_ct == 0) return nullptr;
	
	int control_count = 0;
	for(int j = 0; j < c_ct; j++)
		if(control3d_list.get(j)->bid == item->id) ++control_count;
	if(control_count == 0) return nullptr;

	DataVector<std::shared_ptr<const Control3dItem>> analytic;
	DataVector<std::shared_ptr<const Control3dItem>> part_analytic;
	DataVector<std::shared_ptr<const Control3dItem>> full_analytic;
	DataVector<std::shared_ptr<Control3dItem>> discrete;
	DataVector<std::shared_ptr<const Control3dItem>> multi_part_analytic;
	DataVector<std::shared_ptr<const Control3dItem>> multi_full_analytic;

	for(int j = 0; j < c_ct; j++){
		auto control_item = control3d_list[j];
		if(control_item->bid != item->id) continue;

		switch(control_item->control_type){
			case Control3dItem::CS_INTERNAL:
				return std::make_shared<ControlSpace3dInternal>();
			case Control3dItem::CS_POINT:
			case Control3dItem::CS_SEGMENT:
			case Control3dItem::CS_CURVE:
			case Control3dItem::CS_TRIANGLE:
				discrete.add(control_item); break;
			case Control3dItem::CS_DOMAIN:
				if(control_item->multi_nr == 0){ // i.e. no multi-domain
					analytic.add(control_item);
					if(control_item->full_analytic())
						full_analytic.add(control_item);
					else
						part_analytic.add(control_item);
				}else{
					if(control_item->full_analytic())
						multi_full_analytic.add(control_item);
					else
						multi_part_analytic.add(control_item);
				}
				break;
		}
	}

	if(multi_full_analytic.countInt() > 0 || multi_part_analytic.countInt() > 0){
		if(multi_full_analytic.countInt() != 1){
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"Only one full multi-analytic equation (3D) allowed!");
			return nullptr;
		}
		if(analytic.countInt() > 0 || discrete.countInt() > 0){
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"Multi-analytic equations (3D) can't be mixed with other sources!");
			return nullptr;
		}

		// one full analytic CS + several domain analytic CS
		auto aspace = std::make_shared<ControlSpace3dAnalytic>(
			multi_full_analytic.get(0)->str_equations[0].c_str(),	// lx,ly,lz
			multi_full_analytic.get(0)->str_equations[1].c_str(), 
			multi_full_analytic.get(0)->str_equations[2].c_str(), 
			multi_full_analytic.get(0)->str_equations[3].c_str(),	// ax,ay,az
			multi_full_analytic.get(0)->str_equations[4].c_str(),
			multi_full_analytic.get(0)->str_equations[5].c_str(),
			m_ctable);
		for(int j = multi_part_analytic.countInt()-1; j >=0 ; j--){
			aspace->addCaseAsFirst(
				multi_part_analytic.get(j)->str_equations[6].c_str(),
				multi_part_analytic.get(j)->str_equations[7].c_str(),
				multi_part_analytic.get(j)->str_equations[0].c_str(), 
				multi_part_analytic.get(j)->str_equations[1].c_str(), 
				multi_part_analytic.get(j)->str_equations[2].c_str(),
				multi_part_analytic.get(j)->str_equations[3].c_str(),
				multi_part_analytic.get(j)->str_equations[4].c_str(),
				multi_part_analytic.get(j)->str_equations[5].c_str(),
				m_ctable);
			LOG4CPLUS_INFO(MeshLog::logger_console, 
				"multi-part as first: " 
					<< multi_part_analytic.get(j)->str_equations[6] 
					<< " and " 
					<< multi_part_analytic.get(j)->str_equations[7]);
		}
		return aspace;
	}else if(full_analytic.countInt() > 0 && part_analytic.countInt() == 0 && discrete.countInt() == 0){
		// one full analytic CS only
		return std::make_shared<ControlSpace3dAnalytic>(
				full_analytic.get(0)->str_equations[0].c_str(), 
				full_analytic.get(0)->str_equations[1].c_str(), 
				full_analytic.get(0)->str_equations[2].c_str(), 
				full_analytic.get(0)->str_equations[3].c_str(), 
				full_analytic.get(0)->str_equations[4].c_str(), 
				full_analytic.get(0)->str_equations[5].c_str(),
				m_ctable);
	}

	// else: analytic + discrete
	START_CLOCK("MeshBRep::prepareUserCS");
	DBox box = volume->getBoundingBox();
	for(int m = 0; m < analytic.countInt(); m++){
		auto& control_item = analytic[m];
		if(control_item->bbox.valid){
			box.addBox(control_item->bbox);
		}
	}

	box.inflate(ControlSpace2d::param_inflate_box_factor);
	CS3dAPtr space;

	switch(ControlSpace3d::param_control_type){
	case MeshData::CONTROL_UNIFORM_3D:
		space = std::make_shared<ControlSpace3dMatrixUniform>(box,
			ControlSpace3dMatrixUniform::param_uniform_nx, 
			ControlSpace3dMatrixUniform::param_uniform_nx, 
			ControlSpace3dMatrixUniform::param_uniform_nx);
		break;
	case MeshData::CONTROL_OCTREE_3D:
		space = std::make_shared<ControlSpace3dOctree>(box);
		break;
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console, 
			"Unknown CONTROL_TYPE_3D!");
		return nullptr;
	}
	LOG_ASSERT(space);

	// 1. initialize the default CS-nodes with domain-analytic sources
	for(int m = 0; m < analytic.countInt(); m++){
		auto& control_item = analytic[m];
		space->forEachControlNode([&](ControlNode3d& node) {
			if (node.w >= 0.0 && // not yet initialized
				control_item->equations[6].getValue(node.coord.x, node.coord.y, node.coord.z) >= 0.0 && // within the equation domain
				control_item->equations[7].getValue(node.coord.x, node.coord.y, node.coord.z) >= 0.0)
			{
				ControlDataStretch3d cds(
					control_item->equations[0].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[1].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[2].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[3].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[4].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[5].getValue(node.coord.x, node.coord.y, node.coord.z));
				if (node.w == 0.0) {
					node.control_data = DMetric3d::stretchToMatrixWithAdjust(cds);
					node.w = -1.0;
				}
				else {
					node.control_data /= node.w;
					node.control_data.setMinimum(DMetric3d::stretchToMatrixWithAdjust(cds));
					node.w = -1.0;
				}
			}
		});
	}

	// 2. insert discrete sources
	if(full_analytic.countInt() > 0){
		// ... updating the existing CS
		space->interpolate();
		for(int m = 0; m < discrete.countInt(); m++){
			auto& control_item = discrete[m];
			control_item->adjust();
			switch(control_item->control_type){
			case Control3dItem::CS_POINT: // in point (with spherical vicinity)
				for(int i = 0; i < control_item->pts.countInt(); i++)
					space->setMinControl(ControlDataExtMatrix3dSphere(
						control_item->pts[i], control_item->radius, control_item->data));
				break;
			case Control3dItem::CS_SEGMENT: // along segment (with cylindrical vicinity)
				space->setMinControl(ControlDataExtMatrix3dSegment(
					control_item->pts[0], control_item->pts[1] - control_item->pts[0], 
					control_item->radius, control_item->data));
				break;
			case Control3dItem::CS_CURVE: // along curve (with cylindrical vicinity)
				LOG_ASSERT(false); // TODO
				break;
			case Control3dItem::CS_TRIANGLE: // for triangle ( with prism vicinity)
				space->setMinControl(ControlDataExtMatrix3dTriangle(
					control_item->pts[0], 
					control_item->pts[1] - control_item->pts[0], 
					control_item->pts[2] - control_item->pts[0], 
					control_item->radius, control_item->data));
				break;
			default:
				LOG_ASSERT(false);
			}
		}
		space->smoothen();
		space->compact();
	}else{ // discrete yes, but no full analytic
		// ... initializing the CS
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"User control space with discrete/partial sources, but without full analytic - no ready yet");
		assert(false); // to be done
		space->interpolate();
	}

	// 3. final touch -> adjust newly added control nodes for domain-analytic sources
	for(int m = 0; m < analytic.countInt(); m++){
		auto& control_item = analytic[m];
		space->forEachControlNode([&](ControlNode3d & node) {
			if (control_item->equations[6].getValue(node.coord.x, node.coord.y, node.coord.z) >= 0.0 &&
				control_item->equations[7].getValue(node.coord.x, node.coord.y, node.coord.z) >= 0.0)
				// within the equation domain
			{
				ControlDataStretch3d cds(
					control_item->equations[0].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[1].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[2].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[3].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[4].getValue(node.coord.x, node.coord.y, node.coord.z),
					control_item->equations[5].getValue(node.coord.x, node.coord.y, node.coord.z));
				assert(node.w < 0.0);
				node.control_data.setMinimum(DMetric3d::stretchToMatrixWithAdjust(cds));
			}
		});
	}

	STOP_CLOCK("MeshBRep::prepareUserCS");

	return space;
}

/// Clear unused surfaces, curves and points
bool MeshBRep::clearUnusedData(MeshContainer3d* mesh)
{
	if(hash_surfaces && hash_surfaces->countInt() > 0){
		DataVector<SurfacePtr> surfaces(hash_surfaces->countInt());
		hash_surfaces->getValues(surfaces);
		for(int i = 0; i < surfaces.countInt(); i++){
			if(surfaces[i].use_count() == 1){
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"Unused base surface, id: " << surfaces[i]->getIntTag(TagExtended::TAG_ID));
			}
		}
	}

	if(hash_curves && hash_curves->countInt() > 0){
		DataVector<Curve2dPtr> curves(hash_curves->countInt());
		hash_curves->getValues(curves);
		for(int i = 0; i < curves.countInt(); i++){
			if(curves[i].use_count() == 1){
				if(!curves[i]->availableTag(TagExtended::TAG_ACTIVE)){
					LOG4CPLUS_WARN(MeshLog::logger_console, 
						"Unused base curve, id: " << curves[i]->getIntTag(TagExtended::TAG_ID));
				}
			}
		}
	}

	if(!mesh) return false;

	for(int i = 0; i < mesh->getPointsCount(); )
		if(mesh->getPointAt(i)->getRank() == 0)
			mesh->removeMeshPoint(i);
		else ++i;

	return true;
}

void MeshBRep::markBoundaryTags(MeshContainer3d * mesh)
{
	for (auto it = mesh->getFirstFace(); it.isValid(); it.nextFace()) {
		MeshFace* face = it.getFace();
		face->setBorderFlags(face->isBoundedBothSides() ? TagBorder::INNER : TagBorder::OUTER);
	}
	for (auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()) {
		MeshEdge3d* edge = it.getEdge();
		bool any_outer_face = false;
		int fct = edge->getFaceCount();
		for (int i = 0; !any_outer_face && i < fct; i++) {
			any_outer_face = edge->getFaceAt(i)->isBorder(TagBorder::OUTER);
		}
		edge->setBorderFlags(any_outer_face ? TagBorder::OUTER : TagBorder::INNER);
	}
	int pct = mesh->getPointsCount();
	for (int i = 0; i < pct; i++) {
		MeshPoint3d* point = mesh->getPointAt(i);
		int rank = point->getRank();
		bool any_outer_edge = false;
		for (int j = 0; !any_outer_edge && j < rank; j++) {
			MeshEdge3d* edge = point->getEdge(j);
			any_outer_edge = edge->isBorder(TagBorder::OUTER);
		}
		point->setBorderFlags(any_outer_edge ? TagBorder::OUTER : TagBorder::INNER);
	}
}
