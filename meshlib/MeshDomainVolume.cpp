// MeshDomainVolume.cpp: implementation of the MeshDomainVolume class.
//
//////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "MeshDomainVolume.h"
#include "MeshDomainSurface.h"

#include "MeshEdge2d.h"
#include "MeshFace.h"
#include "MeshPoint2d.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "SurfaceParametric.h"
#include "MeshTriangle3d.h"
#include "MeshQuad3d.h"
#include "MeshGenerator2d.h"
#include "MeshGenerator3dSurface.h"
#include "MeshGenerator3d.h"
#include "MeshGenerator3dAdapt.h"
#include "MeshGenerator3dDelaunayBoundary.h"
#include "MeshGenerator3dDirectBoundary.h"
#include "MeshGenerator3dQuality.h"
#include "MeshBoundaryCondition.h"
#include "ControlSpace3d.h"
#include "ControlSpace2d.h"
#include "ControlSpace2dAdaptive.h"
#include "ControlSpace3dMatrixUniform.h"
#include "ControlSpace3dOctree.h"
#include "ControlSpace3dKdTreeL.h"
#include "ControlSpace3dKdTreeV.h"
#include "OctTree.h"
#include "DataHashTable.h"
#include "IteratorEdge3d.h"

#include "MeshViewSet.h"

double MeshDomainVolume::param_min_volume_for_simple_convex = 1.0;

MeshDomainVolume::MeshDomainVolume(DataVector<MeshFace*> &new_faces, 
						DataVector<MeshPoint3d*> &new_points, 
						DataVector<bool> &orientation)
	: MeshBlock((int)new_faces.countInt(), (int)new_points.countInt()), 
		m_surface_mesh(nullptr), m_mesh(nullptr), 
		m_user_space(nullptr), m_control_space(nullptr)
{
	if(face_count > 0){
		faces = new MeshFace*[face_count];
		for(int i = 0; i < face_count; i++){
			faces[i] = new_faces[i]; 
			faces[i]->setBlockLink(this, 
				orientation[i] ? MeshFace::BLOCK_UP : MeshFace::BLOCK_DOWN); 
		}
	}
	if(point_count > 0){
		points = new MeshPoint3d*[point_count];
		for(int i = 0; i < point_count; i++)
			points[i] = new_points[i];
	}
}

MeshDomainVolume::MeshDomainVolume(DataVector<MeshFace*> &new_faces, 
						DataVector<MeshPoint3d*> &new_points)
	: MeshBlock((int)new_faces.countInt(), (int)new_points.countInt()), 
		m_surface_mesh(nullptr), m_mesh(nullptr), 
		m_user_space(nullptr), m_control_space(nullptr)
{
	if(face_count > 0){
		faces = new MeshFace*[face_count];
		for(int i = 0; i < face_count; i++){
			faces[i] = new_faces[i]; 
			faces[i]->setBlockLink(this, MeshFace::BLOCK_UP); 
		}
	}
	if(point_count > 0){
		points = new MeshPoint3d*[point_count];
		for(int i = 0; i < point_count; i++)
			points[i] = new_points[i];
	}
}

MeshDomainVolume::MeshDomainVolume()
	: MeshBlock(0, 0), m_surface_mesh(nullptr), m_mesh(nullptr), 
		m_user_space(nullptr), m_control_space(nullptr)
{
}

MeshDomainVolume::MeshDomainVolume(MeshContainer3d* mesh3d)
	: MeshBlock(0, 0), m_surface_mesh(nullptr), m_mesh(mesh3d), 
		m_user_space(nullptr), m_control_space(nullptr)
{
	if(mesh3d && mesh3d->getBlocksCount() > 0){
		setAreaID(mesh3d->getBlockAt(0)->getAreaID());
		//mesh3d->markBoundary();
		// // recreate (identify) surface mesh
		//copySurfaceMeshFromVolumeMesh();
		m_control_space = mesh3d->getControlSpace();
		if(!m_control_space){
			m_control_space = MeshGenerator3dAdapt::createACSfromMeshBlocks(mesh3d);
			mesh3d->setControlSpace(m_control_space);
		}
	}
}

MeshDomainVolume::MeshDomainVolume(MeshContainer3dSurface* surface_mesh)
	: MeshBlock(0, 0), m_surface_mesh(surface_mesh), m_mesh(nullptr), 
		m_user_space(nullptr), m_control_space(nullptr)
{
	if(surface_mesh){
		int fct = surface_mesh->getFacesCount();
		for(int i = 0; i < fct; i++)
			surface_mesh->getFaceAt(i)->setBlockLink(this, 0);
		m_control_space = surface_mesh->getControlSpace();
		if(!m_control_space){
			m_control_space = MeshGenerator3dSurface::createACSfromMeshFaces(surface_mesh);
			surface_mesh->setControlSpace(m_control_space);
		}
	}
}

MeshDomainVolume::~MeshDomainVolume()
{
	if(faces) delete faces;
	if(points) delete points;

	if(m_surface_mesh && !m_surface_mesh->isDeleted() ) delete m_surface_mesh;
	if(m_mesh) delete m_mesh;

	m_freepoints->forEach([](auto dv) {
		dv->preDeleteAll(); });
}

/// clear some data for faster all-delete process
void MeshDomainVolume::preDeleteAll()
{
	// nothing here ...
}

DBox MeshDomainVolume::getBoundingBox() const
{
	if(m_mesh) 
		return m_mesh->getBoundingBox();
	else if(m_surface_mesh)
		return m_surface_mesh->getBoundingBox();
	else{
		DBox box;
		for(int i = 0; i < face_count; i++){
			box.addBox(faces[i]->getBoundingBox());
		}
		return box;
	}
}

void MeshDomainVolume::clearDiscretization()
{
	if(m_mesh){
		delete m_mesh;
		m_mesh = nullptr;
	}
	if(m_surface_mesh){
		delete m_surface_mesh;
		m_surface_mesh = nullptr;
	}
}

CS3dPtr MeshDomainVolume::createInitialControlSpace()
{
	m_control_space.reset();

	if(m_user_space && !m_user_space->isAdaptive())
		return m_control_space = m_user_space;

	if(m_surface_mesh){
		CS3dPtr surf_space = m_surface_mesh->getControlSpace();
		if(surf_space)
			return m_control_space = surf_space;
	}

	const DBox box = getBoundingBox();
	CS3dAPtr space;

	switch(ControlSpace3d::param_control_type){
	case MeshData::CONTROL_UNIFORM:
		space = std::make_shared<ControlSpace3dMatrixUniform>(box, 
			ControlSpace3dMatrixUniform::param_uniform_nx, 
			ControlSpace3dMatrixUniform::param_uniform_nx, 
			ControlSpace3dMatrixUniform::param_uniform_nx);
		break;
	case MeshData::CONTROL_OCTREE_3D:
		space = std::make_shared<ControlSpace3dOctree>(box);
		break;
	case MeshData::CONTROL_KDTREE_3D_L:
		space = std::make_shared<ControlSpace3dKdTreeL>(box);
		break;
	case MeshData::CONTROL_KDTREE_3D_V:
		space = std::make_shared<ControlSpace3dKdTreeV>(box);
		break;
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Unknown CONTROL_TYPE!");
		break;
	}
	assert(space);

	space->setMaxMetric();

	if(m_user_space){
		space->applyAsMinimum(m_user_space);
		space->smoothen();
	}else if(m_surface_mesh){
		space->addSimpleBoundaryControlData(m_surface_mesh);
		space->smoothen();
	}

	return m_control_space = space;
}

int MeshDomainVolume::prepareBoundaryMesh()
{
	if( m_surface_mesh == nullptr ){
		if( face_count == 0 ) {
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "No data for preparing boundary mesh in MeshDomainVolume");
			return 0;
		}else{
			int pct_2d = 0;
			for(int i = 0; i < face_count; i++){
				MeshDomainSurface* face = (MeshDomainSurface*)faces[i];
				assert(face->getType() == FACE_DOMAIN);
				MeshContainer2d* mesh2d = face->getMesh();
				if( mesh2d != nullptr ) pct_2d += mesh2d->getPointsCount();
			}
			m_surface_mesh = new MeshContainer3dSurface( pct_2d );
			DataHashTableKeyValue<MeshPoint3d*,MeshPoint3d*> hbpoints( pct_2d, nullptr );
			for(int i = 0; i < face_count; i++){
				MeshDomainSurface* face = (MeshDomainSurface*)faces[i];
				if(!face->addToSurfaceMesh( m_surface_mesh, hbpoints, this )){
					LOG4CPLUS_ERROR(MeshLog::logger_console, 
						"Cannot prepare boundary mesh3d for domain face: " << i);
					return 0;
				}
			}
			m_surface_mesh->setControlSpace(getControlSpace());
//			SHOW_MESH("MDV::prepareBoundaryMesh", m_surface_mesh->getViewSet());
		}
	}

	assert(m_surface_mesh->getControlSpace());

	int spct = m_surface_mesh->getPointsCount();
	if( spct == 0 ){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "No data for preparing boundary mesh in MeshDomainVolume");
		return 0;
	}

	DataVector<std::shared_ptr<ControlNode3d>> control3d_nodes(spct);

	START_CLOCK("3D ControlSpace2d - calculate metric from boundary");

	double stretch_max_ratio = std::min(
		ControlSpace2dAdaptive::param_contour_stretch_max_ratio,
		ControlSpace2dAdaptive::param_stretch_max_ratio);

	if(face_count > 0){
		// retrieve control metric from faces
		for(int i = 0; i < face_count; i++){
			MeshContainer2d* face_mesh = ((MeshDomainSurface*)faces[i])->getMesh();
			if(!face_mesh) continue;
			CS2dPtr control = face_mesh->getControlSpace();
			SurfaceConstPtr surface = face_mesh->getSurface();
			assert(control);
			int pct = face_mesh->getPointsCount();
			for(int j = 0; j < pct; j++){
				const DPoint2d dpt = face_mesh->getPointAt(j)->getCoordinates();
				ControlDataMatrix2d cdm = control->getMetricAtPoint(dpt);
				DMatrix2d e;
				double d[3];
				bool success = cdm.eigensystem(e, d);
				assert(success); if(!success) continue;
				// 
				const DVector3d pu = surface->getDerivative(DEquation::deriv_du, dpt);
				const DVector3d pv = surface->getDerivative(DEquation::deriv_dv, dpt);
				const DVector2d pt_u = e.column(0);
				const DVector2d pt_v = e.column(1);
				const DVector3d e0 = DVector3d(pt_u.x * pu.x + pt_u.y * pv.x,
					pt_u.x * pu.y + pt_u.y * pv.y, pt_u.x * pu.z + pt_u.y * pv.z).normalized();
				const DVector3d e1 = DVector3d(pt_v.x * pu.x + pt_v.y * pv.x,
					pt_v.x * pu.y + pt_v.y * pv.y, pt_v.x * pu.z + pt_v.y * pv.z);
				// after surface transformation -> e0 and e1 are not necessarily orthogonal
				const DVector3d e2 = e0.crossProduct(e1).normalized();
				const DVector3d e1x = e2.crossProduct(e0).normalized(); 
				// rescale d[1] ??????
				d[2] = stretch_max_ratio * std::min(d[0],d[1]);
				control3d_nodes.add(std::make_shared<ControlNode3d>(
					surface->getPoint(dpt), ControlDataMatrix3d(e0, e1x, e2, d)));
			}
		}
	}

	STOP_CLOCK("3D ControlSpace2d - calculate metric from boundary");

/*
	DataSimpleList<MeshPoint3d*> marked_points;

	int i_face = 0;
	int i_point = 0;

	for(int i = 0; i < face_count; i++){
		MeshDomainSurface* dface = (MeshDomainSurface*)faces[i];
		MeshContainer2d* face_mesh = dface->getMesh();
		SurfaceParametric* face_surface = dface->getBaseSurface();
		MeshBoundaryCondition* boundary_condition = (MeshBoundaryCondition*)dface->getPtrTag(TagExtended::TAG_BOUNDARY_COND);

		int mfct = face_mesh->getElementsCount();
		DataVector<MeshPoint3d*> point2d_links(face_mesh->getPointsCount(), nullptr); // 3D points created for this surface mesh
		for(int j = 0; j < mfct; j++){
			MeshElement* element = face_mesh->getElementAt(j);
			assert(element);
			int ect = element->getEdgeCount();
			MeshPoint3d* temp_points[4];
			MeshPoint3d* temp_bound_points[4];
			for(int k = 0; k < ect; k++){
				MeshPoint2d* point = element->getPoint(k);
				MeshPoint3d* point3d = (MeshPoint3d*)point->getPtrTag(TagExtended::TAG_MP_2D_3D);
				if(point3d){
					temp_bound_points[k] = point3d;
					temp_points[k] = (MeshPoint3d*)point3d->getPtrTag(TagExtended::TAG_MDV);
					if(!temp_points[k]){
						temp_points[k] = new MeshPoint3d(*point3d);;
						point3d->setPtrTag(TagExtended::TAG_MDV, temp_points[k]);
						marked_points.insert(point3d);
						temp_points[k]->setBorder();
						m_boundary_mesh->addMeshPoint(temp_points[k]);
//						mesh_points[i_point] = temp_points[k];
						++i_point;
					}
				}else{
					temp_bound_points[k] = nullptr;
					temp_points[k] = point2d_links[point->getIndex()];
					if(!temp_points[k]){
						point2d_links[point->getIndex()] = temp_points[k] = 
							new MeshPoint3d(face_surface->getPoint(point->getCoordinates()));
						temp_points[k]->setBorder();
						m_boundary_mesh->addMeshPoint(temp_points[k]);
//						mesh_points[i_point] = temp_points[k];
						++i_point;
					}
				}
			}
			MeshVolume* volume = (MeshVolume*)dface->getBlock(0);
			if(!volume) volume = (MeshVolume*)dface->getBlock(1);
			assert(volume != nullptr);
			MeshFace* mface = nullptr;
			switch(ect){
			case 3:
				if(volume->validFaceOrientation(dface))
					mface = new MeshTriangle3d(temp_points[0], temp_points[1], temp_points[2]); 
				else mface = new MeshTriangle3d(temp_points[2], temp_points[1], temp_points[0]); 
				break;
			case 4:
				if(volume->validFaceOrientation(dface))
					mface = new MeshQuad3d(temp_points[0], temp_points[1], temp_points[2], temp_points[3]); 
				else mface = new MeshQuad3d(temp_points[3], temp_points[2], temp_points[1], temp_points[0]); 
				break;
			default:
				LOG4CPLUS_WARN(MeshLog::logger_console, "MeshDomainVolume::createBoundaryMesh - Strange element");
			}
			//mesh_faces[i_face] = face;
			//mesh_faces[i_face]->setIntTag(TagExtended::TAG_MDV, element->getAreaID());
			if(boundary_condition){
				boundary_condition->incRefCount();
				mface->setPtrTag(TagExtended::TAG_BOUNDARY_COND, boundary_condition);
			}
			for(int k = 0; k < ect; k++){
				int k_next = (k+1)%ect;
				if(temp_bound_points[k] && temp_bound_points[k_next]){	// boundary condition only for boundary edges
					MeshEdge3d* bound_edge = temp_bound_points[k]->getEdgeToPoint(temp_bound_points[k_next]);
					if(!bound_edge) continue;
					MeshBoundaryCondition* condition = 
						(MeshBoundaryCondition*)bound_edge->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
					if(!condition) continue;
					MeshEdge3d *temp_edge = temp_points[k]->getEdgeToPoint(temp_points[k_next]);
					assert(temp_edge);
					if(temp_edge->availableTag(TagExtended::TAG_BOUNDARY_COND)) continue;	// already set
					temp_edge->setPtrTag(TagExtended::TAG_BOUNDARY_COND, condition);
				}
			}
			MeshVolume* dvol0 = (MeshVolume*)dface->getBlock(0);
			MeshVolume* dvol1 = (MeshVolume*)dface->getBlock(1);
			MeshVolume* local_vol0 = dvol0 ? volumePairs.getValue(dvol0, nullptr) : nullptr;
			MeshVolume* local_vol1 = dvol1 ? volumePairs.getValue(dvol1, nullptr) : nullptr;
			bool valid_orient = volume->validFaceOrientation(dface);
			if(local_vol0){
				int ind = volumeIndex.getValue(local_vol0, -1);
				assert(ind >= 0);
				new_mesh_faces[ind].add(mface);
				new_mesh_faces_orientation[ind].add(valid_orient);
			}
			if(local_vol1){
				int ind = volumeIndex.getValue(local_vol1, -1);
				assert(ind >= 0);
				new_mesh_faces[ind].add(mface);
				new_mesh_faces_orientation[ind].add(!valid_orient);
			}
			mface->setBorder();
			++i_face;
		}
	}
	LOG4CPLUS_INFO(MeshLog::logger_console, "Domain-volume faces", i_face);

	while(marked_points.notEmpty())
		marked_points.removeFirst()->removeTag(TagExtended::TAG_MDV);
	*/

	START_CLOCK("MG3d - create initial ACS");

	if(!m_control_space) createInitialControlSpace();
	auto space = m_control_space->getAsAdaptive();

	if(space && !space->isCompacted() && !control3d_nodes.empty()){
		const DBox& box = space->getBoundingBox();
		for(int i = 0; i < control3d_nodes.countInt(); i++){
			auto node = control3d_nodes.get(i);
			space->setMinControl(box.fitInPoint(node->coord), node->control_data);
		}
		//else space->addSimpleBoundaryControlDataIsotroPIc(m_boundary_mesh);
		space->smoothen();
	}

	STOP_CLOCK("MG3d - create initial ACS");

	if(space) space->markInsideNodes(m_surface_mesh);

	if (false) {
		MeshViewSet* set = new MeshViewSet;
		int re = 0;
		int fe = 0;
		int oe = 0;
		int ie = 0;
		int te = 0;
		for (auto it = m_surface_mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()) {
			MeshEdge3d* edge = it.getEdge();
			if (!edge->isBorder()) continue;
			set->addEdge(edge);
			te++;
			if (edge->isBorder(TagBorder::OUTER)) oe++;
			if (edge->isBorder(TagBorder::INNER)) ie++;
			if (edge->isBorder(TagBorder::RIDGE)) re++;
			if (edge->isBorder(TagBorder::FIXED)) fe++;
		}
		set->addInfo("total edges", te);
		set->addInfo("ridge edges", re);
		set->addInfo("fixed edges", fe);
		set->addInfo("inner edges", ie);
		set->addInfo("outer edges", oe);
		SHOW_MESH("Boundary edges after prepare-surface-mesh", set);
	}


	return m_surface_mesh->getFacesCount();
}


int MeshDomainVolume::createTetrahedralMesh()
{
	if(!m_control_space) return 0;
	Metric3dContext mc(m_control_space);
	// create mesh of boundary faces only
	int boundary_count = discretizeUsingBoundary(mc);
	if(boundary_count < 1){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "ERROR TETRAHEDRALIZING BOUNDARY");
		return 0;
	}
//	m_last_quality_mode = MeshData::QVIEW_NONE;

	// create inner mesh
	LOG4CPLUS_INFO(MeshLog::logger_console, "Gen3d. Inserting inner nodes...");
	int inner_count = MeshGenerator3d::addInnerNodes(mc, m_mesh);

	LOG4CPLUS_INFO(MeshLog::logger_console, 
		"Total node count " << (boundary_count+inner_count)
		<< " (" << boundary_count << "b+" << inner_count << "i)");

	return boundary_count+inner_count;
}

int MeshDomainVolume::discretizeUsingBoundary(Metric3dContext& mc)
{

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Boundary faces: " << m_surface_mesh->getFacesCount());
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Boundary vertices: " << m_surface_mesh->getPointsCount());

//	m_last_quality_mode = MeshData::QVIEW_NONE;

	if(m_mesh) delete m_mesh;

	if(MeshGenerator3d::param_boundary_triangulation_method == 0){ // mixed
		DPoint3d center;
		double min_volume = MeshGenerator3dDirectBoundary::minStarTetrahedraVolume(
								mc, m_surface_mesh, this, center);

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "star min vol: " << min_volume);

		const double MIN_VOL_WITH_SIZING = 
			MeshGenerator3d::param_quality_threshold * MeshTetrahedron::getOptimumVolume() *
			MeshDomainVolume::param_min_volume_for_simple_convex;
		const double MIN_VOL_WITHOUT_SIZING = 0.001;

		if(min_volume >= MIN_VOL_WITH_SIZING ||
			(!m_freepoints->empty() && (min_volume >= MIN_VOL_WITHOUT_SIZING)))
		{
		//if(min_volume >= MIN_VOL_WITHOUT_SIZING)
		//{
			m_mesh = MeshGenerator3dDirectBoundary::createSimpleConvexBoundaryConstrainedMesh(
				mc, m_surface_mesh, this, center);
		}

	//	if(!m_mesh) m_mesh = createFrontalBoundaryConstrainedMesh(mc, bfaces, bpoints, this);

		if(!m_mesh)
			m_mesh = MeshGenerator3dDelaunayBoundary::createBoundaryConstrainedMesh(
				mc, m_surface_mesh, this);

		if(!m_mesh && min_volume >= MIN_VOL_WITHOUT_SIZING)
			m_mesh = MeshGenerator3dDirectBoundary::createSimpleConvexBoundaryConstrainedMesh(
							mc, m_surface_mesh, this, center);

	}else if(MeshGenerator3d::param_boundary_triangulation_method == 1){ // Delaunay only
		m_mesh = MeshGenerator3dDelaunayBoundary::createBoundaryConstrainedMesh(
			mc, m_surface_mesh, this);
	}

	assert(m_mesh && m_mesh->getControlSpace());

	auto cs = m_mesh->getControlSpace();
	if (cs)
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "discUsingBoundary: cs-" << cs->getType());
	else
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "discUsingBoundary: cs-empty");

	if (true) {
		MeshViewSet* set = new MeshViewSet;
		int re = 0;
		int fe = 0;
		int oe = 0;
		int ie = 0;
		int te = 0;
		for (auto it = m_mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()) {
			MeshEdge3d* edge = it.getEdge();
			if (!edge->isBorder()) continue;
			set->addEdge(edge);
			te++;
			if (edge->isBorder(TagBorder::OUTER)) oe++;
			if (edge->isBorder(TagBorder::INNER)) ie++;
			if (edge->isBorder(TagBorder::RIDGE)) re++;
			if (edge->isBorder(TagBorder::FIXED)) fe++;
		}
		set->addInfo("total edges", te);
		set->addInfo("ridge edges", re);
		set->addInfo("fixed edges", fe);
		set->addInfo("inner edges", ie);
		set->addInfo("outer edges", oe);
		SHOW_MESH("Boundary edges after contraining", set);
	}

	return m_mesh ? m_mesh->getPointsCount() : 0;
}

/// Stores surface mesh points into text file
void MeshDomainVolume::storeSurfacePoints(const string& fname, int index) const
{
	ostringstream fnamep;
	fnamep << fname << '-' << index << "-p.txt";
	ofstream filep(fnamep.str().c_str());

	if( m_surface_mesh == nullptr ){
		filep << "empty mesh" << endl;
	}else{
		// POINTS
		int pct = m_surface_mesh->getPointsCount();
		filep <<pct << endl;
		for(int i = 0; i < pct; i++){
			const DPoint3d& pt = m_surface_mesh->getPointAt(i)->getCoordinates();
			filep << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
		}
	}
}

void MeshDomainVolume::storeSurfaceMeshTxt(const string &fname, int idx) const
{
	if(m_control_space && m_control_space->getType() == MeshData::CONTROL_OCTREE_3D){
		ostringstream fname_cs;
		fname_cs << fname << '-' << idx << "-cs.txt";
		((ControlSpace3dOctree*)(m_control_space.get()))->storeTXT(fname_cs.str().c_str());
	}

	ostringstream fname1,fname2,fname3;
	fname1 << fname << '-' << idx << "-p.txt";
	fname2 << fname << '-' << idx << "-i.txt";
	fname3 << fname << '-' << idx << "-e.txt";
	ofstream file1(fname1.str().c_str());
	ofstream file2(fname2.str().c_str());
	ofstream file3(fname3.str().c_str());

	if( m_surface_mesh == nullptr ){
		file1 << "empty mesh" << endl;
		file2 << "empty mesh" << endl;
		file3 << "empty mesh" << endl;
	}else{
		// POINTS
		int pct = m_surface_mesh->getPointsCount();
		file1 << pct << endl;
		for(int i = 0; i < pct; i++){
			const DPoint3d& pt = m_surface_mesh->getPointAt(i)->getCoordinates();
			file1 << i << "\t" << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
		}
		// ELEMENTS
		int fct = m_surface_mesh->getFacesCount();
		file2 << fct << endl;
		for(int j = 0; j < fct; j++){
			MeshFace* face = m_surface_mesh->getFaceAt(j);
			int ect = face->getEdgeCount();
			file2 << j << '\t' << ect;
			for(int l = 0; l < ect; l++)
				file2 << '\t' << face->getPoint(l)->getIndex();
			file2 << endl;
		}
		// EDGES
		DataCompoundList<MeshEdge3d*> all_edges;
		for(IteratorEdge3d it = m_surface_mesh->getFirstEdge3d(); it.isValid(); it.nextEdge())
			all_edges.append( it.getEdge() );

		file3 << all_edges.countInt() << endl;
		int i = 0;
		for(auto it = all_edges.iterator(); it.valid(); it.moveNext()){
			file3 << (i++);
			file3 << '\t' << it.item()->getMeshPoint(0)->getIndex();
			file3 << '\t' << it.item()->getMeshPoint(1)->getIndex();
			file3 << endl;
		}
	}
}

void MeshDomainVolume::setMesh(MeshContainer3d *mesh)
{
	if(m_mesh) delete m_mesh;
	m_mesh = mesh;
}

/// Removes and returns the reference to the associated discretization (if created)
MeshContainer3d* MeshDomainVolume::removeMesh()
{
	MeshContainer3d* mesh = m_mesh;
	m_mesh = nullptr;
	return mesh;
}

void MeshDomainVolume::setSurfaceMesh(MeshContainer3dSurface *mesh)
{
	if(m_surface_mesh) delete m_surface_mesh;
	m_surface_mesh = mesh;
}

void MeshDomainVolume::clearControlSpace()
{
	m_control_space.reset();
}

bool MeshDomainVolume::smoothenVolumeMesh(int steps, 
		TagExtended::TagType tag_type, int tag_value, int method)
{
	if(!m_mesh) return false;

	Metric3dContext mc(m_mesh->getControlSpace());

	bool result = MeshGenerator3dQuality::smoothen(mc, m_mesh, steps, tag_type, tag_value, false, method);

	LOG_ASSERT(m_mesh->isValid());

	return result;
}

void MeshDomainVolume::storeAmiraMesh(const string & /* fname */, int /* idx */) const
{
/* 
	// TODO need update for multi-block volumes

	if(m_boundary_mesh){
		ostringstream fstr_name;
		fstr_name << fname << "-surface-" << idx << ".amiramesh";
		ofstream file(fstr_name.str().c_str());

		file << "# AmiraMesh ASCII 1.0\n\n";
		int pct = m_boundary_mesh->getPointsCount();
		file << "define Nodes " << pct << endl;
		int bct = m_boundary_mesh->getBlocksCount();
		int tri_count = 0;
		int quad_count = 0;
		for(int i = 0; i < bct; i++){
			MeshBlock* block = m_boundary_mesh->getBlockAt(i);
			tri_count += block->getFaceCount(3);
			quad_count += block->getFaceCount(4);
		}
		if(tri_count > 0)
			file << "define Triangles " << tri_count << endl;
		if(quad_count > 0)
			file << "define Quads " << quad_count << endl;

		file << endl << "Nodes { float[3] Coordinates } = @1" << endl;
		if(tri_count > 0)
			file << "Triangles { int[3] Nodes } = @2" << endl;
		if(quad_count > 0)
			file << "Quads { int[4] Nodes } = @3" << endl;

		file << endl << "@1" << endl;
		for(int i = 0; i < pct; i++){
			const DPoint3d pt = m_boundary_mesh->getPointAt(i)->getCoordinates();
			file << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
		}
		if(tri_count > 0){
			file << endl << "@2" << endl;
			for(int i = 0; i < bct; i++){
				MeshBlock* block = m_boundary_mesh->getBlockAt(i);
				int fct = block->getFaceCount();
				for(int j = 0; j < fct; j++){
					MeshFace* face = block->getFace(j);
					if(face->getBlock(0) == block ||
						(face->getBlock(0) == nullptr && face->getBlock(1) == block)){
						int ect = face->getEdgeCount();
						if(ect == 3)
							file << (1+face->getPoint(0)->getIndex()) <<
								'\t' << (1+face->getPoint(1)->getIndex()) <<
								'\t' << (1+face->getPoint(2)->getIndex()) << endl;
					}
				}
			}
		}
		if(quad_count > 0){
			file << endl << "@3" << endl;
			for(int i = 0; i < bct; i++){
				MeshBlock* block = m_boundary_mesh->getBlockAt(i);
				int fct = block->getFaceCount();
				for(int j = 0; j < fct; j++){
					MeshFace* face = block->getFace(j);
					if(face->getBlock(0) == block ||
						(face->getBlock(0) == nullptr && face->getBlock(1) == block)){
						int ect = face->getEdgeCount();
						if(ect == 4)
							file << (1+face->getPoint(0)->getIndex()) <<
								'\t' << (1+face->getPoint(1)->getIndex()) <<
								'\t' << (1+face->getPoint(2)->getIndex()) <<
								'\t' << (1+face->getPoint(3)->getIndex()) << endl;
					}
				}
			}
		}
	}
	if(m_boundary_mesh){
		ostringstream fstr_name;
		fstr_name << fname << "-surface-" << idx << ".stl";
		ofstream file(fstr_name.str().c_str());

		file << "solid\n";
		int bct = m_boundary_mesh->getBlocksCount();
		for(int i = 0; i < bct; i++){
			MeshBlock* block = m_boundary_mesh->getBlockAt(i);
			int fct = block->getFaceCount();
			for(int j = 0; j < fct; j++){
				MeshFace* face = block->getFace(j);
				if(face->getBlock(0) == block ||
					(face->getBlock(0) == nullptr && face->getBlock(1) == block)){
					int ect = face->getEdgeCount();
					const DVector3d nv = face->getNormalVector();
					file << "  facet normal " << nv.x << " " << nv.y << " " << nv.z << endl;
					file << "    outer loop\n";
					for(int k = 0; k < ect; k++){
						const DPoint3d pt = face->getPoint(k)->getCoordinates();
						file << "      vertex " << pt.x << " " << pt.y << " " << pt.z << endl;
					}
					file << "    endloop\n";
					file << "  endfacet\n";
				}
			}
		}
	}
	if(m_mesh){
		ostringstream fstr_name;
		fstr_name << fname << "-volume-" << idx << ".amiramesh";
		ofstream file(fstr_name.str().c_str());

		file << "# AmiraMesh ASCII 1.0\n\n";
		int pct = m_mesh->getPointsCount();
		file << "define Nodes " << pct << endl;
		int bct = m_mesh->getBlocksCount();
		file << "define Tetrahedra " << bct << endl;

		file << endl << "Nodes { float[3] Coordinates } = @1" << endl;
		file << "Tetrahedra { int[4] Nodes } = @2" << endl;

		file << endl << "@1" << endl;
		for(int i = 0; i < pct; i++){
			const DPoint3d pt = m_mesh->getPointAt(i)->getCoordinates();
			file << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
		}
		file << endl << "@2" << endl;
		for(int i = 0; i < bct; i++){
			MeshBlock* block = m_mesh->getBlockAt(i);
			int bpct = block->getPointCount();
			assert(bpct == 4);
			for(int j = 0; j < bpct; j++)
				file << (1+block->getPoint(j)->getIndex()) << '\t';
			file << endl;
		}
	}
*/
}

bool MeshDomainVolume::checkControlForCloseBoundaryFaces(Metric3dContext& mc) const
{
//	if(!m_mesh || m_freepoints.empty()) return true; // only if there are freepoints ...
/*
	if(!m_mesh) return true;
	CS3dPtr space = (CS3dPtr)m_mesh->getControlSpace();
	if(!space->isAdaptive()) return true; // no change possible

	bool size_changes = false;

	int pct = m_mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint3d* mpt = m_mesh->getPointAt(i);
		if(!mpt->isBorder()) continue;
		int rank = mpt->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = mpt->getEdge(j);
			if(!edge->isBorder()){
				if(edge->getOtherPoint(mpt)->isBorder() && edge->getPointIndex(mpt) == 0){
					double len = edge->getLengthMetricAdapted(mc);
					if(len < 0.8){
						size_changes |= space->updateForBoundarySegment(edge);
					}
				}else continue;
			}
		}
	}

	if(ControlSpace2dAdaptive::param_gradation_ratio > 0.0 && size_changes){
		if(space->smoothen())
			return checkControlAtBoundary(mc);
	}
*/
	return true;
}

bool MeshDomainVolume::checkControlAtBoundary(Metric3dContext& mc) const
{
	bool proper = true;
	
	// retrieve control metric from faces
	for(int i = 0; i < face_count; i++){
		MeshContainer2d* face_mesh = ((MeshDomainSurface*)faces[i])->getMesh();
		if(!face_mesh){
			//LOG4CPLUS_ERROR(MeshLog::logger_console,   "MDV:checkControlAtBoundary: missing face mesh");
			continue;
		}
		CS2dPtr control = face_mesh->getControlSpace();
		if(!control){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "MDV:checkControlAtBoundary: missing face control space");
			continue;
		}
		if(!control->isAdaptive()){
			LOG4CPLUS_WARN(MeshLog::logger_console, "MDV:checkControlAtBoundary: non-adaptive face control space");
			continue;
		}
		SurfaceConstPtr surface = face_mesh->getSurface();

		bool face_proper = true;
		for(IteratorEdge2d it = face_mesh->getFirstEdge2d(); it.isValid(); it.nextEdge()){
			const MeshPoint2d* mpt0 = it.getEdge()->getMeshPoint(0);
			const MeshPoint2d* mpt1 = it.getEdge()->getMeshPoint(1);
			MeshPoint3d* mpt3d = (MeshPoint3d*)mpt0->getPtrTag(TagExtended::TAG_MP_2D_3D);
			const DPoint3d dpt0 = mpt3d ? mpt3d->getCoordinates() : surface->getPoint(mpt0->getCoordinates());
			mpt3d = (MeshPoint3d*)mpt1->getPtrTag(TagExtended::TAG_MP_2D_3D);
			const DPoint3d dpt1 = mpt3d ? mpt3d->getCoordinates() : surface->getPoint(mpt1->getCoordinates());
			mc.countMetricAtPoint(DPoint3d::average(dpt0, dpt1));
			double mlen = mc.transformRStoMS(dpt1-dpt0).length();
			if(mlen < 0.5){
				LOG4CPLUS_WARN(MeshLog::logger_console, "Low metric length 3d: " << mlen);
			}else if(mlen > 2.5){
				// 
				double len = dpt0.distance(dpt1) / mlen;
				const DPoint2d& pt0 = mpt0->getCoordinates();
				const DPoint2d& pt1 = mpt1->getCoordinates();

				const DVector2d v01 = pt1-pt0;
				const DVector2d du = v01.normalized();
				const DVector2d dv(du.y, -du.x); // orthogonal
				double d[2] = {len, len * ControlSpace2dAdaptive::param_stretch_max_ratio};
				ControlDataMatrix2d cdm(du, dv, d);	// from eigensystem
				if(control->getAsAdaptive()->setMinControl(
					ControlDataExtMatrix2dSegment(pt0, v01, d[0], cdm, surface)))
					face_proper = proper = false;
			}
		}

		if(!face_proper) ((MeshDomainSurface*)faces[i])->clearDiscretization(); // needs remeshing
	}

	return proper;
}

/// classify free-point for this volume-block
bool MeshDomainVolume::offerFreePoint(std::shared_ptr<MeshPoint3d> fpoint)
{
	const double TOL2 = mesh_data.relative_small_number * mesh_data.relative_small_number;
	// check vertices
	int pct = getPointCount();
	for(int i = 0; i < pct; i++){
		MeshPoint3d * vpt = getPoint(i);
		if(vpt->getCoordinates().distance2(fpoint->getCoordinates()) < TOL2){
			MeshBoundaryCondition* mbc = (MeshBoundaryCondition*) fpoint->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
			if(mbc) vpt->setPtrTag(TagExtended::TAG_BOUNDARY_COND, mbc);
			return true;
		}
	}
	// check faces
	int fct = getFaceCount();
	for(int i = 0; i < fct; i++){
		MeshDomainSurface* mds = (MeshDomainSurface*) getFace(i);
		assert(mds->getType() == FACE_DOMAIN);
		if(mds->offerFreePoint(fpoint)) return true;
	}
	m_freepoints->add(fpoint);
	return true;
}

/// Transforms volume mesh into surface mesh (removes tetrahedra and inner entities)
MeshContainer3dSurface* MeshDomainVolume::cutSurfaceMeshFromVolumeMesh()
{
	if(!m_mesh) return nullptr;
	int pct = m_mesh->getPointsCount();
	int bct = m_mesh->getBlocksCount();
	if(bct == 0) return nullptr;

	START_CLOCK("MC3dS:cutSurfaceMeshFromVolumeMesh");

	int spct = 0;
	for(int i = 0; i < pct; i++){
		if(m_mesh->getPointAt(i)->isBorder()) spct++;
	}

	if(m_surface_mesh) delete m_surface_mesh;
	m_surface_mesh = new MeshContainer3dSurface(spct);

	DataHashTableKeyValue<int, MeshDomainVolume*> hash_mdv(100, -1);

	int total_mdv_count = 0;
	// replace tetrahedral blocks with the dummy domain_volume
	DataHashTable<MeshFace*> hfaces(2*bct, nullptr);
	while(m_mesh->getBlocksCount() > 0){
		MeshBlock* block = m_mesh->getBlockAt(0);
		int area_id = block->getAreaID();
		MeshDomainVolume* mdv = hash_mdv.getValue(area_id, nullptr);
		if(!mdv){
			total_mdv_count++;
			if(total_mdv_count == 1){ // for first "material block" use this one,
				mdv = this;
			}else{
				// else, create additional "dummy" volumes
				auto smdv = std::make_shared<MeshDomainVolume>();
				mdv = smdv.get();
				m_surface_mesh->addDomainVolume(smdv);
			}
			mdv->setAreaID(area_id);
			hash_mdv.insert(area_id, mdv);
		}
		int fct = block->getFaceCount();
		DataVector<int> indices(fct, -1);
		DataVector<MeshFace*> faces(fct);
		for(int i = 0; i < fct; i++){
			MeshFace* face = block->getFace(i);
			faces.add(face);
			if(face->isBorder()){
				indices[i] = face->getBlockIndex(block);
				if(hfaces.insert(face)){
					m_surface_mesh->addMeshFace(face);
				}
			}
		}
		delete m_mesh->removeMeshBlock(0);
		for(int i = 0; i < fct; i++){
			if(indices[i] >= 0)
				faces[i]->setBlockLink(mdv, indices[i]);
		}
	}

	// remove obsolete points and move the right ones
	while(m_mesh->getPointsCount() > 0){
		MeshPoint3d* point = m_mesh->removeMeshPoint(0);
		if(point->isBorder()){
			m_surface_mesh->addMeshPoint(point);
		}else{
			delete point;
		}
	}

	m_surface_mesh->setControlSpace(m_mesh->getControlSpace());

	delete m_mesh;
	m_mesh = nullptr;

	STOP_CLOCK("MC3dS:cutSurfaceMeshFromVolumeMesh");
	return m_surface_mesh;
}

/// Creates surface mesh from volume mesh (copy boundary faces and entities)
MeshContainer3dSurface* MeshDomainVolume::copySurfaceMeshFromVolumeMesh()
{
	if(!m_mesh) return nullptr;

	if(m_surface_mesh) delete m_surface_mesh;
	m_surface_mesh = MeshGenerator3dSurface::copySurfaceMeshFromVolumeMesh(m_mesh, this);

	return m_surface_mesh;
}

/// Returns the "screenshot" of this domain-surface for visualization
MeshViewSet* MeshDomainVolume::getViewSet(MeshViewSet* set) const
{
	if(m_surface_mesh){
		if(!set) set = new MeshViewSet(0, 3*m_surface_mesh->getPointsCount(), 0, 0);
		else set->prepareFreePlace(0, 3*m_surface_mesh->getPointsCount(), 0, 0);
		for(IteratorEdge3d it = m_surface_mesh->getFirstEdge3d(); it.isValid(); it.nextEdge())
			set->addEdge(it.getEdge(), -1);
	}else if(m_mesh){
		if(!set) set = new MeshViewSet(0, m_mesh->getPointsCount(), 0, 0);
		else set->prepareFreePlace(0, m_mesh->getPointsCount(), 0, 0);
		for(IteratorEdge3d it = m_mesh->getFirstEdge3d(); it.isValid(); it.nextEdge())
			if(it.getEdge()->isBorder())
				set->addEdge(it.getEdge(), -1);
	}

	return set;
}

/// Stores existing volume mesh into text file with Abaqus format
bool MeshDomainVolume::storeAbaqus(const string& fname, const string& name) const
{
	if(!m_mesh) return false;

	ofstream file(fname.c_str());

	if(!file){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Can't open file for writing: " << fname);
		return false;
	}

	file << "*Heading" << endl;
	file << "** Job name: " << fname << " Model name: " << name << endl;
	file << "** Generated by: QMeshGen " << mesh_data.version() << endl;
	file << "*Preprint, echo=NO, model=NO, history=NO, contact=NO" << endl;
	file << "**" << endl;
	file << "** PARTS" << endl;
	file << "**" << endl;
	file << "*Part, name=" << name << endl;

	// *Node
	// nr wezla, wsp x, wsp y, wsp z
	// ...
	file << "*Node" << endl;
	int pct = m_mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		const DPoint3d& dpt = m_mesh->getPointAt(i)->getCoordinates();
		file << setw(7) << (1+i) << ", " << fixed << setprecision(10) << dpt.x << ", " << dpt.y << ", " << dpt.z << endl;
	}

	// *Element, type=C3D4 (lub inny, w zaleznosci od rodzaju)
	// nr elementu, wezel1, wezel2, wezel3, wezel4
	// ...
	file << "*Element, type=C3D4" << endl;
	const int rct[] = {1, 0, 2, 3}; // for switching orientation of tetrahedra (sequence of vertices)
	int bct = m_mesh->getBlocksCount();
	for(int i = 0; i < bct; i++){
		MeshBlock* block = m_mesh->getBlockAt(i);
		assert(block->getPointCount() == 4);
		file << setw(6) << (1+i);
		for(int j = 0; j < block->getPointCount(); j++){
			file << ", " << (1+block->getPoint(rct[j])->getIndex());
		}
		file << endl;
	}

	// *Elset, elset=nazwa (elementy nalezace do jednego zianra)
	// element1, element2...
	DataHashTable<int> visited_grains(2*bct, -1);
	for(int i = 0; i < bct; i++){
		int area_id = m_mesh->getBlockAt(i)->getAreaID();
		if(visited_grains.contains(area_id))
			continue;
		visited_grains.insert(area_id);
		file << "*Elset, elset=G" << area_id << endl;
		file << setw(7) << (1+i);
		int counter = 0;
		for(int j = i+1; j < bct; j++){
			if(m_mesh->getBlockAt(j)->getAreaID() != area_id)
				continue;
			counter = (counter+1)%16;
			file << (counter ? ", " : "\n") << setw(7) << (1+j);
		}
		file << endl;
	}
	
	// *Nset, nset=nazwa (wezly siatki nalezace do scian zewn)
	// wezel1, wezel2, wezel3...
	DVector3d vns[6] = {
		DVector3d(-1.0, 0.0, 0.0),
		DVector3d( 1.0, 0.0, 0.0),
		DVector3d( 0.0,-1.0, 0.0),
		DVector3d( 0.0, 1.0, 0.0),
		DVector3d( 0.0, 0.0,-1.0),
		DVector3d( 0.0, 0.0, 1.0)};
	DataVector<bool> visited_points(pct, false);
	for(int k = 0; k < 6; k++){
		file << "*Nset, nset=N" << (1+k);
		int counter = -1;
		for(IteratorFace it = m_mesh->getFirstFace(); it.isValid(); it.nextFace()){
			MeshFace* face = it.getFace();
			if(face->isBoundedBothSides()) 
				continue;
			DVector3d fvn = face->getNormalVector();
			if(face->getBlock(0) == nullptr) fvn.reverse();
			if(vns[k].scalarProduct( fvn ) < 0.9) 
				continue;
			for(int i = 0; i < face->getPointCount(); i++){
				int pid = face->getPoint(i)->getIndex();
				if(visited_points[pid]) continue;
				counter = (counter+1)%16;
				file << (counter ? ", " : "\n") << setw(7) << (1+pid);
				visited_points[pid] = true;
			}
		}
		file << endl;
	}

	return true;
}

/// Remove and join the volume mesh from the given mdv with this one
bool MeshDomainVolume::combineVolumeMeshFrom(MeshDomainVolume* mdv)
{
	if(!m_mesh){
		m_mesh = mdv->removeMesh();
		return m_mesh != nullptr;
	}
	// else combine
	MeshContainer3d* other_mesh = mdv->removeMesh();

	int pct = m_mesh->getPointsCount();
	assert(pct > 0);
	DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> hpoints(2*pct, nullptr);
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = m_mesh->getPointAt(i);
		if(!point->isBorder()) continue;
		MeshPoint3d* bpoint = (MeshPoint3d*)point->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
		assert(bpoint);
		hpoints.insert(bpoint, point);
	}

	// move blocks
	while(other_mesh->getBlocksCount() > 0){
		MeshBlock* block = other_mesh->removeMeshBlock(0);
		m_mesh->addMeshBlock(block);
		for(int i = 0; i < block->getPointCount(); i++){
			MeshPoint3d* point = block->getPoint(i);
			if(!point->isBorder()) continue;
			MeshPoint3d* bpoint = (MeshPoint3d*)point->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
			assert(bpoint);
			if(hpoints.contains(bpoint))
				block->switchPointsWithFacesBoundary(point, hpoints.getValue(bpoint, nullptr));
		}
	}

	// move points
	while(other_mesh->getPointsCount() > 0){
		MeshPoint3d* point = other_mesh->removeMeshPoint(0);
		if(point->isBorder()){
			MeshPoint3d* bpoint = (MeshPoint3d*)point->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
			assert(bpoint);
			if(hpoints.contains(bpoint)){
				while(point->getRank() > 0){
					MeshEdge3d* edge = point->getEdge(0);
					if(edge->isBorder()) edge->clearBorder();
					if(edge->getFaceCount() == 0) 
						delete edge;
					else
						delete edge->getFaceAt(0);
				}
				delete point;
			}else 
				m_mesh->addMeshPoint(point);
		}else
			m_mesh->addMeshPoint(point);
	}

//	SHOW_MESH("Mesh after combining", m_mesh->getViewSet());

	delete other_mesh;

	return true;
}
