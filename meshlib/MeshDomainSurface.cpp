// MeshDomainSurface.cpp: implementation of the MeshDomainSurface class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshDomainSurface.h"
#include "MeshDomainVolume.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshDomainEdge3d.h"
#include "SurfaceParametric.h"
#include "MeshContainer2d.h"
#include "SurfacePlane.h"
#include "MeshArea.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "MeshQuad2d.h"
#include "QuadTree.h"
#include "MeshGenerator2d.h"
#include "MeshGenerator2dQuad.h"
#include "MeshGenerator2dQMorph.h"
#include "ControlSpace2d.h"
#include "ControlSpace2dMatrixUniform.h"
#include "EPSFile.h"
#include "DEquation.h"
#include "DTriangle.h"
#include "DQuad.h"
#include "MeshViewSet.h"
#include "MeshTriangle3d.h"
#include "MeshQuad3d.h"
#include "IteratorFace.h"
#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshGenerator3dSurface.h"

//#define PRINT_METRIC_QUALITY

MeshDomainSurface::MeshDomainSurface(int ct, MeshPoint3d** new_points, int ect, MeshEdge3d** new_edges, 
				SurfaceConstPtr surface, MeshContainer2d* boundary)
	: MeshFace(ct, new_edges, new_points), edge_count(ect), m_user_space(nullptr), m_surface(surface),
		m_boundary(boundary), m_boundary_mesh(nullptr), m_mesh(nullptr), 
		m_last_quality_mode(MeshData::QVIEW_NONE)
{
	for(int i = 0; i < edge_count; i++)
		if(edges[i]) edges[i]->addFaceLink(this);

	if(!surface && count > 2) createPlaneSurface();
	if(!boundary && count > 2) createSimpleBoundary();
}

MeshDomainSurface::MeshDomainSurface(DataVector<MeshPoint3d*> &new_points, 
			DataVector<MeshDomainEdge3d*> &new_edges, 
			SurfaceConstPtr surface, MeshContainer2d* boundary)
	: MeshFace((int)new_points.countInt(), nullptr, nullptr), 
		edge_count((int)new_edges.countInt()), m_user_space(nullptr), 
		m_surface(surface), m_boundary(boundary), m_boundary_mesh(nullptr), m_mesh(nullptr), 
		m_last_quality_mode(MeshData::QVIEW_NONE)
{
	points = new MeshPoint3d*[count];
	for(int i = 0; i < count; i++) points[i] = new_points[i];
	edges = new MeshEdge3d*[edge_count];
	for(int i = 0; i < edge_count; i++){
		edges[i] = new_edges[i];
		if(edges[i]) edges[i]->addFaceLink(this);
	}

	if(!surface && count > 2) createPlaneSurface();
	if(!boundary && count > 2) createSimpleBoundary();
}

//MeshDomainSurface::MeshDomainSurface(DataVector<MeshPoint3d*> &disc_points, DataVector<MeshFace*> &disc_faces, MeshDomainVolume* mdv)
//	: MeshFace(0, nullptr, nullptr), edge_count(0), m_user_space(nullptr), 
//		m_surface(nullptr), m_boundary(nullptr), m_boundary_mesh(nullptr), m_mesh(nullptr), 
//		m_mesh3d_faces(nullptr), m_mesh3d_points(nullptr), m_mesh3d_points_new(nullptr), 
//		m_last_quality_mode(MeshData::QVIEW_NONE)
//{
//	// discrete points -> mesh3d points
//	int pct = disc_points.countInt();
//	m_mesh3d_points = new DataVector<MeshPoint3d*>(pct);
//	m_mesh3d_points_new = new DataPtrVector<MeshPoint3d>(pct);
//	for(int i = 0; i < pct; i++){
//		disc_points[i]->setBorder();
//		m_mesh3d_points->add(disc_points[i]);
//		m_mesh3d_points_new->add(disc_points[i]);
//	}
//
//	// discrete faces -> mesh3d faces
//	int fct = disc_faces.countInt();
//	m_mesh3d_faces = new DataVector<MeshFace*>(fct);
//	for(int i = 0; i < fct; i++){
//		disc_faces[i]->setBorder();
//		disc_faces[i]->setBlockLink(mdv, MeshFace::BLOCK_UP);
//		m_mesh3d_faces->add(disc_faces[i]);
//	}
//}

MeshDomainSurface::~MeshDomainSurface()
{
	if(count > 0){
		if(edges[0]){
			for(int i = 0; i < edge_count; i++){	
				// Remove edge if connected with this face only.
				if(edges[i] && edges[i]->removeFaceLink(this)) delete edges[i];
			}
		}
		if(points) delete[] points;
		if(edges) delete[] edges;
	}
	if(m_boundary) delete m_boundary;
	if(m_boundary_mesh) delete m_boundary_mesh;
	if(m_mesh) delete m_mesh;
	for(int i = 0; i < m_freepoints.countInt(); i++)
		m_freepoints[i]->preDeleteAll();
}

/// Removes link to this edge (fo delete-all phase) - returns true if last one
bool MeshDomainSurface::removeEdgeLink(MeshEdge3d* edge)
{
	bool no_more_edges = true;
	for(int i = 0; i < edge_count; i++){
		if(edges[i] == edge) edges[i] = 0;
		else if(edges[i]) no_more_edges = false;
	}

	return no_more_edges;
}

bool MeshDomainSurface::createPlaneSurface()
{
	if(count < 3) return false;
	if(m_surface){
		m_surface.reset();
		if(m_boundary){
			delete m_boundary;
			m_boundary = nullptr;
		}
	}
	const DVector3d v01 = (points[1]->getCoordinates() - points[0]->getCoordinates()).normalized();
	const DVector3d v02 = points[2]->getCoordinates() - points[0]->getCoordinates();
	const DVector3d v03 = v01.crossProduct(v02);
	const DVector3d v04 = v03.crossProduct(v01).normalized();
	m_surface = std::make_shared<SurfacePlane>(points[0]->getCoordinates(), v01, v04);

	return true;
}

bool MeshDomainSurface::createSimpleBoundary()
{
	if(count < 3) return false;
	int i;

	DPoint2d* params = new DPoint2d[count];
	for(i = 0; i < count; i++)
		params[i] = m_surface->getParameters(points[i]->getCoordinates());

	//create MeshArea as boundary
	m_boundary = new MeshContainer2d(count);
	MeshPoint2d** mesh_points = new MeshPoint2d*[count];
	for(i = 0; i < count; i++){
		mesh_points[i] = new MeshPoint2d(params[i]);
		mesh_points[i]->setPtrTag(TagExtended::TAG_MP_2D_3D, points[i]);
		mesh_points[i]->setBorder();
		m_boundary->addMeshPoint(mesh_points[i]);
	}
	for(i = 0; i < edge_count; i++){
		MeshEdge2d* edge = new MeshEdge2d(mesh_points[i], mesh_points[(i+1)%count]);
		MeshEdge3d* edge3d = points[i]->getEdgeToPoint(points[(i+1)%count]);
		edge->setPtrTag(TagExtended::TAG_ME_2D_3D, edge3d);
	}
	MeshArea* area = new MeshArea(count, mesh_points);
	m_boundary->addMeshElement(area);

	delete[] params;
	
	return true;
}

void MeshDomainSurface::clearControlSpace()
{
	if(m_boundary) m_boundary->setControlSpace(nullptr);
	if(m_mesh) m_mesh->setControlSpace(nullptr);
	if(m_boundary_mesh) m_boundary_mesh->setControlSpace(nullptr);
}

void MeshDomainSurface::clearDiscretization()
{
	if(m_mesh){
		delete m_mesh;
		m_mesh = nullptr;
	}
	if(m_boundary_mesh){
		delete m_boundary_mesh;
		m_boundary_mesh = nullptr;
	}
}

int MeshDomainSurface::createBoundaryMesh()
{
	if(!m_boundary) return 0; 
	if(m_boundary_mesh) delete m_boundary_mesh;

	m_last_quality_mode = MeshData::QVIEW_NONE;

	int bcount = m_boundary->getPointsCount();
	int fcount = m_freepoints.countInt();
	int point_count = bcount + fcount;
	m_boundary_mesh = new MeshContainer2d(point_count);

	MeshPoint2d** pts = new MeshPoint2d*[point_count];

	// boundary points
	for(int i = 0; i < bcount; i++){
		MeshPoint2d* point = m_boundary->getPointAt(i);
		pts[i] = new MeshPoint2d(*point);
		pts[i]->setBorder();
		pts[i]->setPtrTag(TagExtended::TAG_MP_2D_3D, point->getPtrTag(TagExtended::TAG_MP_2D_3D));
		m_boundary_mesh->addMeshPoint(pts[i]);
	}

	// boundary - free points
	for(size_t i = 0; i < fcount; i++){
		const auto& point = m_freepoints[i];
		DPoint2d param = m_surface->getParameters(point->getCoordinates());
		pts[i+bcount] = new MeshPoint2d(param);
		pts[i+bcount]->setBorder();
		pts[i+bcount]->setPtrTag(TagExtended::TAG_MP_2D_3D, point.get());
		if(point->availableTag(TagExtended::TAG_ID))
			pts[i+bcount]->setIntTag(TagExtended::TAG_ID, point->getIntTag(TagExtended::TAG_ID, 0));
		m_boundary_mesh->addMeshPoint(pts[i+bcount]);
	}

	int area_count = m_boundary->getElementsCount();
	MeshArea** areas = new MeshArea*[area_count];

	for(int i = 0; i < area_count; i++){
		areas[i] = new MeshArea(0, nullptr);	// only dummies
		MeshElement* element = m_boundary->getElementAt(i);
		areas[i]->setAreaID(element->getAreaID());
		areas[i]->setFilled(element->isFilled());
		m_boundary_mesh->addMeshElement(areas[i]);
	}

	double max_diff = 0.0;

	for(int i = 0; i < bcount; i++){
		MeshPoint2d* point = m_boundary->getPointAt(i);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0){
				MeshArea* area_left = (MeshArea*)edge->getMeshElement(point);				
				MeshArea* area_right = (MeshArea*)edge->getOtherElement(area_left);
				if(area_left) area_left = areas[area_left->getIndex()];
				if(area_right) area_right = areas[area_right->getIndex()];
				MeshPoint2d* last_point = pts[i];
				MeshDomainEdge3d* edge3d = (MeshDomainEdge3d*)edge->getPtrTag(TagExtended::TAG_ME_2D_3D);
				if(edge3d == nullptr){
					// degenerated edge (joins the same points)
					MeshPoint2d* new_point = pts[edge->getOtherPoint(point)->getIndex()];
					MeshEdge2d* new_edge = edge->cloneGeometric(last_point, new_point, 0.0, 1.0, area_left, area_right);
					new_edge->addDirection(1, last_point);
					new_edge->copyBorderFlagsFrom(edge);
					continue;
				}
				assert(edge3d->getType() == EDGE_DOMAIN_3D);
				bool same_direction = (edge3d->getPointIndex((MeshPoint3d*)point->getPtrTag(TagExtended::TAG_MP_2D_3D)) == 0);
				// CHANGE-17.03.2010 / 9.09.2010
				double last_ksi = same_direction ? 0.0 : 1.0; // 9.09.2010 -> changed because of sPIral.xml
				// double last_ksi = 0.0;
				auto& disc_points = edge3d->getDiscretization();

				//if(true){
				//	MeshViewSet* set = new MeshViewSet;
				//	MeshPoint3d* mp = edge3d->getMeshPoint(0);
				//	set->addPoint(mp->getCoordinates(), 0, mp->getIntTag(TagExtended::TAG_ID));
				//	mp = edge3d->getMeshPoint(1);
				//	set->addPoint(mp->getCoordinates(), 0, mp->getIntTag(TagExtended::TAG_ID));
				//	for(int m = 0; m < disc_points.countInt(); m++){
				//		MeshPoint3d * point3d = same_direction ? disc_points[m] : disc_points[disc_points.countInt()-m-1];
				//		set->addPoint(point3d->getCoordinates(), 1, m);
				//	}
				//	SHOW_MESH("edge", set);
				//}

				for(int m = 0; m < disc_points.countInt(); m++){
					const auto& point3d = same_direction ? disc_points[m] : disc_points[disc_points.countInt()-m-1];
					const DPoint3d& pt3 = point3d->getCoordinates();
					double ksi = point3d->getDoubleTag(TagExtended::TAG_DISCRETIZATION_PARAM);
					DPoint2d pt = edge->getPoint(ksi); // stored parameter for patch used in generation1d (not necessarily this patch!)
					if(pt3.distance(m_surface->getPoint(pt)) > mesh_data.relative_small_number){
						ksi = last_ksi;
						pt = edge->surfaceParameters(m_surface, pt3, ksi, same_direction);
					}
					assert(ksi != last_ksi);

					// check max diff
					double diff = pt3.distance(m_surface->getPoint(pt));
					if(diff > max_diff) max_diff = diff;

					// create mesh point 2d
					MeshPoint2d* new_point = new MeshPoint2d(pt);
					new_point->setBorder();
					new_point->setPtrTag(TagExtended::TAG_MP_2D_3D, point3d.get());
					m_boundary_mesh->addMeshPoint(new_point);
					MeshEdge2d* new_edge = edge->cloneGeometric(last_point, new_point, last_ksi, ksi, area_left, area_right);
					assert(last_point->getCoordinates().distance(
						new_point->getCoordinates()) > mesh_data.relative_small_number);

					new_edge->addDirection(1, last_point);
					new_edge->copyBorderFlagsFrom(edge);
					new_edge->setPtrTag(TagExtended::TAG_ME_2D_3D, edge3d);

					last_ksi = ksi;
					last_point = new_point;
				}
				MeshPoint2d* new_point = pts[edge->getOtherPoint(point)->getIndex()];
				// CHANGE-17.03.2010
				MeshEdge2d* new_edge = edge->cloneGeometric(last_point, new_point, last_ksi, // 9.09.2010 -> changed because of sPIral.xml
					same_direction ? 1.0 : 0.0, area_left, area_right);
//				MeshEdge2d* new_edge = edge->cloneGeometric(last_point, new_point, last_ksi, 1.0,
//					area_left, area_right);
				assert(last_point->getCoordinates().distance(
					new_point->getCoordinates()) > mesh_data.relative_small_number);

				new_edge->addDirection(1, last_point);
				new_edge->copyBorderFlagsFrom(edge);
				new_edge->setPtrTag(TagExtended::TAG_ME_2D_3D, edge3d);

				if (m_boundary_mesh->getPointsCount() >= MeshPoint2d::param_max_node_count) {
					LOG4CPLUS_ERROR(MeshLog::logger_console, "Maximum number of mesh nodes exceeded.");
					return 0;
				}
			}
		}
	}

	for(int i = 0; i < bcount; i++){
		MeshPoint2d* point = m_boundary_mesh->getPointAt(i);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = point->getEdge(j);
			if(edge->getPtrTag(TagExtended::TAG_ME_2D_3D) == nullptr){
				DPoint2d middle = edge->getPoint(0.5);
				MeshPoint2d* other_point = edge->getOtherPoint(point);
				delete edge;
				int other_rank = other_point->getRank();
				for(int n = 0; n < other_rank; n++){
					MeshEdge2d* other_edge = other_point->getEdge(0);
					other_point->removeEdge(other_edge);
					other_edge->switchPoints(other_point, point);
					point->addEdge(other_edge);
				}
				delete m_boundary_mesh->removeMeshPoint(other_point->getIndex());
				point->setCoordinates(middle);
				--point_count;
				break;
			}
		}
	}

	delete[] pts;
	delete[] areas;

	m_boundary_mesh->setSurface(m_surface);
	CS2dPtr space = m_boundary->getControlSpace();
	if(!space){
		// Prepare control space (introductory factors)
		MeshDomainVolume* dvolume0 = (MeshDomainVolume*)getBlock(0);
		MeshDomainVolume* dvolume1 = (MeshDomainVolume*)getBlock(1);
		CS3dPtr ucs_0 =
			dvolume0 ? dvolume0->getUserControlSpace() : nullptr;
		CS3dPtr ucs_1 =
			(dvolume1 && (dvolume0 != dvolume1)) ? dvolume1->getUserControlSpace() : nullptr;
		m_boundary->createControlSpace(m_user_space, ucs_0, ucs_1);
		space = m_boundary->getControlSpace();
		assert(space);
	}
	m_boundary_mesh->setControlSpace(space);

	if(max_diff > mesh_data.relative_small_number){
		string text;
		DEquation::doubleToString(max_diff, DEquation::v_length, text);
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "max edge diff: " << text.c_str());
	}

	auto space_a = space->getAsAdaptive();
	if(space_a)	space_a->markInsideNodes(m_boundary_mesh);

	return 1;
}

int MeshDomainSurface::triangulate() 
{
	if(m_boundary_mesh == nullptr) return 0;
	// create boundary mesh
	Metric2dContext mc(m_boundary_mesh->getControlSpace());
	int ct = triangulateBoundary(mc);
	if(ct < 1){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "ERROR TRIANGULATING BOUNDARY");
		return 0;
	}

	m_last_quality_mode = MeshData::QVIEW_NONE;

	// create inner mesh
	ct += MeshGenerator2d::addInnerNodes(mc, m_mesh);

	return ct;
}

int MeshDomainSurface::triangulateBoundary(Metric2dContext& mc)
{
	if(m_boundary_mesh == nullptr) return 0;

	m_last_quality_mode = MeshData::QVIEW_NONE;

	//CS2dPtr control = MeshGenerator2d::createUniformControlSpace(m_surface, boundary_mesh);
	if(m_mesh) delete m_mesh;

	m_mesh = MeshGenerator2d::createInitialMesh(m_boundary_mesh);
	if(m_mesh == nullptr) return 0;
//	SHOW_STEP_BREAKABLE(0, "* Triangulacja Delaunay'a - pocz¹tkowa siatka.",true);

	if(!MeshGenerator2d::addBoundaryNodes(mc, m_mesh, m_boundary_mesh)){
//		SHOW_STEP(0, "* B³¹d wstawiania wêz³ów brzegowych.");
		delete m_mesh;
		return 0;
	}

	if(!MeshGenerator2d::constrainToBorder(mc, m_boundary_mesh, m_mesh)){
//		SHOW_STEP(0, "* B³¹d odzyskiwania brzegu.");
		delete m_mesh;
		return 0;
	}

////=========================================
//	if(true){
//		int pct = m_mesh->getPointsCount();
//		MeshViewSet* set = new MeshViewSet(pct, 3*pct, 0, 0);
//		for(IteratorEdge2d it = m_mesh->getFirstEdge2d(); it.isValid(); it.nextEdge()){
//			if(it.getEdge()->isBorder())
//				set->addEdge(it.getEdge(), m_surface);
//		}
//		set = ((CS2dPtr)m_mesh->getControlSpace())->getViewSet(set);
//		SHOW_MESH("xxx", set);
//	}
//
////=========================================


	int points_count = m_mesh->getPointsCount();
//	LOG("\t** total points = %d\n\n", points_count);

	LOG4CPLUS_TRACE(MeshLog::logger_mesh, "triangulateBoundary (" << points_count << ") done.");

	return points_count;
}

bool MeshDomainSurface::smoothen(int steps,
		TagExtended::TagType tag_type, int tag_value, 
		int method)
{
	if(!m_mesh) return false;
	m_last_quality_mode = MeshData::QVIEW_NONE;
	Metric2dContext mc(m_mesh->getControlSpace());

//	SHOW_MESH("Start smothing", m_mesh->getViewSet());

	START_CLOCK("MDS::smoothen");

	bool result = MeshGenerator2d::smoothen(mc, m_mesh, steps, tag_type, tag_value, method);

	STOP_CLOCK("MDS::smoothen");

	assert(m_mesh->isValid());

	return result;
}

int MeshDomainSurface::convertToQuads(int method, int max_ct)
{
	m_last_quality_mode = MeshData::QVIEW_NONE;

	int result = 0;

	Metric2dContext mc(m_mesh->getControlSpace());

	switch(method){
	case MeshData::QUADS_LEELO:
	case MeshData::QUADS_QMORPH:
	case MeshData::QUADS_MIXED:
		result = MeshGenerator2dQuad::convertToQuadsMixed(mc, m_mesh, method, max_ct);
		break;
	case MeshData::QUADS_SPLIT:
		result = MeshGenerator2dQuad::splitToQuads(mc, m_mesh);
		break;
	default:
		assert(false);
		return 0;
	}

	LOG_ASSERT(m_mesh->isValid());

	return result;
}

bool MeshDomainSurface::smoothenQuads(int steps, int method)
{
	if(!m_mesh) return false;
	m_last_quality_mode = MeshData::QVIEW_NONE;

	Metric2dContext mc(m_mesh->getControlSpace());
	MeshGenerator2dQuad::improveQuads(mc, m_mesh, steps, method);

#ifdef QUAD_DEBUG_LOG
	SHOW_MESH("Quad smoothing", m_mesh->getViewSet(nullptr, true, true, true), 10);
#endif

	LOG_ASSERT(m_mesh->isValid());

	return true;
}

void MeshDomainSurface::storeEPS(int id) const
{
	// TODO change into fstream
	ostringstream fname;
	fname << "mesh-" << id << ".eps";

	DRect box = m_mesh->getBoundingRect();
	EPSFile eps(fname.str(), box.x0, box.x1, box.y0, box.y1);

	int pct = m_mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = m_mesh->getPointAt(i);
		DPoint2d pt1 = point->getCoordinates();
		int rank = point->getRank();
		if(rank > 1){
			for(int j = 0; j < rank; j++){
				MeshEdge2d* edge = point->getEdge(j);
				if(edge->getPointIndex(point) == 0){
					eps.drawLine(pt1, edge->getMeshPoint(1)->getCoordinates());
				}
			}
		}
	}
}

void MeshDomainSurface::assessQuality(ofstream &file) const
{
	assert(m_mesh);
	int i,j;
	int ect = m_mesh->getElementsCount();
	int pct = m_mesh->getPointsCount();
	int tct = m_mesh->getElementsCount(3);
	int qct = ect - tct;

	file << "N = ";
	file.width(6);
	file << pct << ", E = ";
	file.width(6);
	file << ect << " [T = ";
	file.width(6);
	file << tct << " Q = ";
	file.width(6);
	file << qct << "]" << endl;

	Metric2dContext mc(m_mesh->getControlSpace());

	if(pct > 0){
		file << endl << " Overall mesh quality " << endl;
		file << "--------------------------------------" << endl;

		file << endl << " Alpha-quality (in metric)" << endl;

		const int alpha_hist_range = 22;
		const int alpha_hist_factor = 20;
		int alpha_histogram[alpha_hist_range];
		MeshData::StatData data(1.0, 0.0, 1.0);
		double n_power = 1.0 / (double)ect;
		double x = 1.0;
		for(i = 0; i < alpha_hist_range; i++)
			alpha_histogram[i] = 0;

		for(i = 0; i < ect; i++){
			MeshElement* element = m_mesh->getElementAt(i);
			double q = element->getAlphaQuality(mc);
			if(q < data.minimum) data.minimum = q;
			if(q> data.maximum) data.maximum = q;
			if(data.minimum > 0.0){
				x *= q;
				if(x < 0.1){
					data.average *= pow(x, n_power);
					x = 1.0;
				}
			}
			int index = (int)(alpha_hist_factor*q);
			if(index < 0) index = 0;
			if(index >= alpha_hist_range) index = alpha_hist_range - 1;
			alpha_histogram[index]++;
		}
		data.average *= pow(x, n_power);

		file << " * Alpha coefficient [min, max, average]" << endl;
		file << "    " << data.minimum << "    " << data.maximum << "    " << data.average << endl;
		file << " * Alpha coefficient [histogram]" << endl;
		for(i = 0; i < alpha_hist_range; i++){
			if(alpha_histogram[i] > 0){
				file.width(4);
				file << ((i+0.5) / alpha_hist_factor);
				file << "   " << alpha_histogram[i] << endl;
			}
		}
		
		file << endl << " Alpha-quality (3D)" << endl;

		//const int alpha_hist_range = 110;
		//const int alpha_hist_factor = 100;
		//int alpha_histogram[alpha_hist_range];
		//MeshData::StatData data;
		data.minimum = 1.0;
		data.maximum = 0.0;
		data.average = 1.0;
		//double n_power = 1.0 / (double)ect;
		x = 1.0;
		for(i = 0; i < alpha_hist_range; i++)
			alpha_histogram[i] = 0;

		for(i = 0; i < ect; i++){
			MeshElement* element = m_mesh->getElementAt(i);
			const DPoint2d middle = element->getMiddlePoint();
			const DPoint3d middle3D = m_surface->getPoint(middle);
			const DVector3d vs = m_surface->getDerivative(DEquation::deriv_ds, middle);
			const DVector3d vt = m_surface->getDerivative(DEquation::deriv_dt, middle);

			const SurfacePlane plane(middle3D, vs, vt);
			const DPoint2d& pt0 = plane.getParameters(m_surface->getPoint(element->getPoint(0)->getCoordinates()));
			const DPoint2d& pt1 = plane.getParameters(m_surface->getPoint(element->getPoint(1)->getCoordinates()));
			const DPoint2d& pt2 = plane.getParameters(m_surface->getPoint(element->getPoint(2)->getCoordinates()));
			double q;
			if(element->getEdgeCount() == 3){
				q = DTriangle2d::alphaQuality(pt0, pt1, pt2);
			}else{
				const DPoint2d& pt3 = plane.getParameters(m_surface->getPoint(element->getPoint(3)->getCoordinates()));
				q = DQuad2d::alphaQuality(pt0, pt1, pt2, pt3);
			}
			if(q < data.minimum) data.minimum = q;
			if(q> data.maximum) data.maximum = q;
			if(data.minimum > 0.0){
				x *= q;
				if(x < 0.1){
					data.average *= pow(x, n_power);
					x = 1.0;
				}
			}
			int index = (int)(alpha_hist_factor*q);
			if(index < 0) index = 0;
			if(index >= alpha_hist_range) index = alpha_hist_range - 1;
			alpha_histogram[index]++;
		}
		data.average *= pow(x, n_power);

		file << " * Alpha coefficient [min, max, average]" << endl;
		file << "    " << data.minimum << "    " << data.maximum << "    " << data.average << endl;
		file << " * Alpha coefficient [histogram]" << endl;
		for(i = 0; i < alpha_hist_range; i++){
			if(alpha_histogram[i] > 0){
				file.width(4);
				file << ((i+0.5) / alpha_hist_factor);
				file << "   " << alpha_histogram[i] << endl;
			}
		}
		
		file << endl << " Edge length (in metric)" << endl;

		const int hist_range = 30;
		int histogram[hist_range];
		for(i = 0; i < hist_range; i++)
			histogram[i] = 0;

		for(i = 0; i < pct; i++){
			MeshPoint2d* point = m_mesh->getPointAt(i);
			int rank = point->getRank();
			if(rank > 1){
				for(j = 0; j < rank; j++){
					MeshEdge2d* edge = point->getEdge(j);
					if(edge->getPointIndex(point) == 0){
						double len = edge->getLength(mc);
						int index = (int)(10*len);
						if(index < 0) index = 0;
						if(index >= hist_range) index = hist_range - 1;
						histogram[index]++;
					}
				}
			}
		}
		for(i = 0; i < hist_range; i++){
			if(histogram[i] > 0){
				file.width(4);
				file << (0.1 * i);
				file << "   " << histogram[i] << endl;
			}
		}
	}

	if(tct > 0){
		file << endl << " Triangular mesh quality " << endl;
		file << "--------------------------------------" << endl;
		double min_angle = PI;
		double max_angle = 0.0;
		for(i=0; i < tct; i++){
			MeshTriangle2d* triangle = (MeshTriangle2d*)m_mesh->getElementAt(i);
			if(triangle->getEdgeCount() != 3) continue;
			for(j = 0; j < 3; j++){
				double angle = triangle->getAngle(mc, j);
				if(angle < min_angle) min_angle = angle;
				if(angle > max_angle) max_angle = angle;
			}
			// min
		}
		file << " * minimum angle: " << (180.0 / PI * min_angle) << endl;
		file << " * maximum angle: " << (180.0 / PI * max_angle) << endl;
	}

	if(qct > 0){
		file << endl << " Quadrilateral mesh quality " << endl;
		file << "--------------------------------------" << endl;
		file << " Angles (in metric)" << endl;

		const int hist_range = 180;
		const int hist_factor = 20;
		int histogram_min[hist_range], histogram_max[hist_range];
		for(i = 0; i < hist_range; i++)
			histogram_min[i] = histogram_max[i] = 0;

		for(i = 0; i < qct; i++){
			MeshQuad2d* quad = (MeshQuad2d*)m_mesh->getElementAt(i);
			if(quad->getEdgeCount() != 4) continue;
			double min_angle = PI;
			double max_angle = 0.0;
			for(j = 0; j < 4; j++){
				double angle = quad->getAngle(mc, j);
				if(angle < min_angle) min_angle = angle;
				if(angle > max_angle) max_angle = angle;
			}
			// min
			int index = (int)(hist_factor*min_angle);
			if(index < 0) index = 0;
			if(index >= hist_range) index = hist_range - 1;
			histogram_min[index]++;
			// max
			index = (int)(hist_factor*max_angle);
			if(index < 0) index = 0;
			if(index >= hist_range) index = hist_range - 1;
			histogram_max[index]++;
		}
		file << " * minimum angle" << endl;
		for(i = 0; i < hist_range; i++){
			if(histogram_min[i] > 0){
				file.width(4);
				file << (180.0 * (i+0.5) / PI / hist_factor);
				file << "   " << histogram_min[i] << endl;
			}
		}
		file << " * maximum angle" << endl;
		for(i = 0; i < hist_range; i++){
			if(histogram_max[i] > 0){
				file.width(4);
				file << (180.0 * (i+0.5) / PI / hist_factor);
				file << "   " << histogram_max[i] << endl;
			}
		}

		if(m_surface){
			file << endl << " Angles (in 3D)" << endl;

			for(i = 0; i < hist_range; i++)
				histogram_min[i] = histogram_max[i] = 0;

			for(i = 0; i < qct; i++){
				MeshQuad2d* quad = (MeshQuad2d*)m_mesh->getElementAt(i);
				if(quad->getEdgeCount() != 4) continue;
				DPoint3d pts[4];
				for(j = 0; j < 4; j++)
					pts[j] = m_surface->getPoint(quad->getPoint(j)->getCoordinates());
				double min_angle = PI;
				double max_angle = 0.0;
				for(j = 0; j < 4; j++){
					DVector3d v0 = pts[(j+1)%4] - pts[j];
					DVector3d v1 = pts[(j+3)%4] - pts[j];
					double angle = v0.getAngle(v1);
					if(angle < min_angle) min_angle = angle;
					if(angle > max_angle) max_angle = angle;
				}
				// min
				int index = (int)(hist_factor*min_angle);
				if(index < 0) index = 0;
				if(index >= hist_range) index = hist_range - 1;
				histogram_min[index]++;
				// max
				index = (int)(hist_factor*max_angle);
				if(index < 0) index = 0;
				if(index >= hist_range) index = hist_range - 1;
				histogram_max[index]++;
			}
			file << " * minimum angle" << endl;
			for(i = 0; i < hist_range; i++){
				if(histogram_min[i] > 0){
					file.width(4);
					file << (180.0 * (i+0.5) / PI / hist_factor);
					file << "   " << histogram_min[i] << endl;
				}
			}
			file << " * maximum angle" << endl;
			for(i = 0; i < hist_range; i++){
				if(histogram_max[i] > 0){
					file.width(4);
					file << (180.0 * (i+0.5) / PI / hist_factor);
					file << "   " << histogram_max[i] << endl;
				}
			}
		}
	}
}

DBox MeshDomainSurface::getBoundingBox() const 
{ 
	DBox box;
	if(m_mesh){
		int pct = m_mesh->getPointsCount();
		for(int i = 0; i < pct; i++)
			box.addPoint(m_surface->getPoint(m_mesh->getPointAt(i)->getCoordinates()));
	}else if(m_boundary_mesh){
		int pct = m_boundary_mesh->getPointsCount();
		for(int i = 0; i < pct; i++)
			box.addPoint(m_surface->getPoint(m_boundary_mesh->getPointAt(i)->getCoordinates()));
	}else if(m_boundary){
		int pct = m_boundary->getPointsCount();
		for(int i = 0; i < pct; i++)
			box.addPoint(m_surface->getPoint(m_boundary->getPointAt(i)->getCoordinates()));
	}
	return box;
}

/// Returns the "screenshot" of this domain-surface for visualization
MeshViewSet* MeshDomainSurface::getViewSet(MeshViewSet* set) const
{
	if(!m_boundary) return set;

	int pct = m_boundary->getPointsCount();
	if(pct < 1) return set;

	if(set){
		set->prepareFreePlace(pct, 10*edge_count, 0, 0);
	}else{
		set = new MeshViewSet(pct, 10*edge_count, 0, 0);
	}

	for(int i = 0; i < pct; i++){
		const MeshPoint2d * mpt = m_boundary->getPointAt(i);
		const MeshPoint3d* mpt3d = (MeshPoint3d*)mpt->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
		if(mpt3d) set->addPoint(mpt3d);
	}

	m_freepoints.forEach([&set](auto fp) {
		set->addPoint(fp->getCoordinates(), 2); });

	DataVector<PolySegment> segments;

	for(IteratorEdge2d it = m_boundary->getFirstEdge2d(); it.isValid(); it.nextEdge()){
		// segments for crossing tests
		// each segment -> two points and whether in-material-state-switching
		DataVector<DPoint2d> polyline;
		MeshEdge2d* edge = it.getEdge();
		edge->getPolyLine(polyline);
		MeshElement* el0 = it.getEdge()->getMeshElement(0);
		MeshElement* el1 = it.getEdge()->getMeshElement(1);
		int la = el0 ? el0->getAreaID() : -2;
		int ra = el1 ? el1->getAreaID() : -2;
		for(int i = 1; i < polyline.countInt(); i++){
			segments.add(PolySegment(polyline[i-1], polyline[i], la, ra));
		}

		MeshDomainEdge3d* domain_edge3d = (MeshDomainEdge3d*) edge->getPtrTag(TagExtended::TAG_ME_2D_3D);
		if(domain_edge3d){
			auto fps = domain_edge3d->getFreePoints();
			if(fps){
				fps->forEach([&set](auto fp) {
					set->addPoint(fp->getCoordinates(), 3); });
			}
		}

		// boundary polyline to draw
		DataVector<DPoint3d> poly_draw;
		it.getEdge()->getPolyLine(poly_draw, m_surface);
		for(int i = 1; i < poly_draw.countInt(); i++){
			set->addEdge(poly_draw[i-1], poly_draw[i], 1);
		}
	}

	DRect rect = m_boundary->getBoundingRect();
	rect.inflate(0.05);

	const int GRID_COUNT = 9;
	const double dy = rect.getDY() / GRID_COUNT;
	const double dx = rect.getDX() / GRID_COUNT;

	// horizontal
	for(int i = 0; i <= GRID_COUNT; i++){
		double y = rect.y0 + i * dy;
		// find segments crossing this line
		DataVector<PolySegment> crossing_segments(segments.countInt());
		for(int j = 0; j < segments.countInt(); j++){
			PolySegment& poly = segments[j];
			poly.incident = 0;
			if(abs(poly.pt0.y - y) < SMALL_NUMBER) poly.incident++;
			if(abs(poly.pt1.y - y) < SMALL_NUMBER) poly.incident++;
			if(poly.incident > 0){	
				// near line
				poly.cross = std::min(poly.pt0.x, poly.pt1.x);
				crossing_segments.add(poly);
			}else if((poly.pt0.y - y) * (poly.pt1.y - y) < 0.0){ 
				// opposite sides
				poly.cross = poly.pt0.x + (y - poly.pt0.y) * (poly.pt1.x - poly.pt0.x) / (poly.pt1.y - poly.pt0.y);
				crossing_segments.add(poly);
			}
		}
		// scan segments by X
		DataSimpleList<int> areas;
		areas.insert(-1); // starting from outside
		double last_x;
		int last_half = 0;
		while(!crossing_segments.empty()){
			// find segment with minimum cross-value
			int min_j = 0;
			PolySegment poly = crossing_segments[0];
			for(int j = 1; j < crossing_segments.countInt(); j++){
				if(crossing_segments[j].cross < poly.cross){
					min_j = j;
					poly = crossing_segments[j];
				}
			}
			// poly-edge orientation
			int left_area = -2;
			int right_area = -2;
			if(poly.pt1.y > poly.pt0.y){
				left_area = poly.left_area;
				right_area = poly.right_area;
			}else{
				left_area = poly.right_area;
				right_area = poly.left_area;
			}
			// edge crossing
			bool crossing = false;
			if(poly.incident == 0){
				// clear crossing
				crossing = true;
			}else if(poly.incident == 1){
				// one vertex on scanning-line
				double dy0 = (poly.pt0.y - y);
				double dy1 = (poly.pt1.y - y);
				int half = sgn((abs(dy0) > abs(dy1)) ? dy0 : dy1);
				if(last_half == 0){
					last_half = half;
				}else if(last_half * half < 0){
					last_half = 0;
					crossing = true;
				}
			}// else nothing to do for two vertices on scanning-line

			if(crossing){
				if(left_area == -2 && right_area == -2){
					// nothing
				}else if(left_area > -2 && right_area > -2){
					// known area from both sides
					// remove left, insert right
					if(areas.empty()){ // something's wrong
						crossing_segments.clear();
						continue;
					}
					areas.removeFirst();
					if(areas.empty()){ // something's wrong
						crossing_segments.clear();
						continue;
					}
					areas.insert(right_area);
				}else if(left_area > -2){
					// known left, unknown right -> "closing area"
					if(areas.empty()){ // something's wrong
						crossing_segments.clear();
						continue;
					}
					areas.removeFirst();
					if(areas.empty()){ // something's wrong
						crossing_segments.clear();
						continue;
					}
					right_area = areas.getFirst();
				}else{
					// known right, unknown left -> "opening area"
					left_area = areas.getFirst();
					areas.insert(right_area);
				}
				// check drawing-line
				if((left_area > -1) != (right_area > -1)){
					if(right_area == -1){ // closing
						// end of segment - polyline to draw
						DataVector<DPoint3d> poly_draw;
						m_surface->getPolyLine(poly_draw, DPoint2d(last_x, y), DPoint2d(poly.cross, y));
						for(int k = 1; k < poly_draw.countInt(); k++){
							set->addEdge(poly_draw[k-1], poly_draw[k], -1);
						}
					}else{
						// beginning of segment
						last_x = poly.cross;
					}
				}
			}

			// remove segment from the crossing-set
			crossing_segments.removeAt(min_j);
		}
	}

	// vertical
	for(int i = 0; i <= GRID_COUNT; i++){
		double x = rect.x0 + i * dx;
		// find segments crossing this line
		DataVector<PolySegment> crossing_segments(segments.countInt());
		for(int j = 0; j < segments.countInt(); j++){
			PolySegment& poly = segments[j];
			poly.incident = 0;
			if(abs(poly.pt0.x - x) < SMALL_NUMBER) poly.incident++;
			if(abs(poly.pt1.x - x) < SMALL_NUMBER) poly.incident++;
			if(poly.incident > 0){	
				// near line
				poly.cross = std::min(poly.pt0.y, poly.pt1.y);
				crossing_segments.add(poly);
			}else if((poly.pt0.x - x) * (poly.pt1.x - x) < 0.0){ 
				// opposite sides
				poly.cross = poly.pt0.y + (x - poly.pt0.x) * (poly.pt1.y - poly.pt0.y) / (poly.pt1.x - poly.pt0.x);
				crossing_segments.add(poly);
			}
		}
		// scan segments by Y
		DataSimpleList<int> areas;
		areas.insert(-1); // starting from outside
		double last_y;
		int last_half = 0;
		while(!crossing_segments.empty()){
			// find segment with minimum cross-value
			int min_j = 0;
			PolySegment poly = crossing_segments[0];
			for(int j = 1; j < crossing_segments.countInt(); j++){
				if(crossing_segments[j].cross < poly.cross){
					min_j = j;
					poly = crossing_segments[j];
				}
			}
			// poly-edge orientation
			int left_area = -2;
			int right_area = -2;
			if(poly.pt1.x < poly.pt0.x){
				left_area = poly.left_area;
				right_area = poly.right_area;
			}else{
				left_area = poly.right_area;
				right_area = poly.left_area;
			}
			// edge crossing
			bool crossing = false;
			if(poly.incident == 0){
				// clear crossing
				crossing = true;
			}else if(poly.incident == 1){
				// one vertex on scanning-line
				double dx0 = (poly.pt0.x - x);
				double dx1 = (poly.pt1.x - x);
				int half = sgn((abs(dx0) > abs(dx1)) ? dx0 : dx1);
				if(last_half == 0){
					last_half = half;
				}else if(last_half * half < 0){
					last_half = 0;
					crossing = true;
				}
			}// else nothing to do for two vertices on scanning-line

			if(crossing){
				if(left_area == -2 && right_area == -2){
					// nothing
				}else if(left_area > -2 && right_area > -2){
					// known area from both sides
					// remove left, insert right
					if(areas.empty()){ // something's wrong
						crossing_segments.clear();
						continue;
					}
					areas.removeFirst();
					areas.insert(right_area);
				}else if(left_area > -2){
					// known left, unknown right -> "closing area"
					if(areas.empty()){ // something's wrong
						crossing_segments.clear();
						continue;
					}
					areas.removeFirst();
					right_area = areas.empty() ? -1 : areas.getFirst();
					//right_area = areas.getFirst();
				}else{
					// known right, unknown left -> "opening area"
					left_area = areas.getFirst();
					areas.insert(right_area);
				}
				if((left_area > -1) != (right_area > -1)){
					if(right_area == -1){ // closing
						// end of segment - polyline to draw
						DataVector<DPoint3d> poly_draw;
						m_surface->getPolyLine(poly_draw, DPoint2d(x, last_y), DPoint2d(x, poly.cross));
						for(int k = 1; k < poly_draw.countInt(); k++){
							set->addEdge(poly_draw[k-1], poly_draw[k], -1);
						}
					}else{
						// beginning of segment
						last_y = poly.cross;
					}
				}
			}

			// remove segment from the crossing-set
			crossing_segments.removeAt(min_j);
		}
	}

	return set;
}

/// Prepare 3d mesh entites from the 2d mesh on the surface
bool MeshDomainSurface::addToSurfaceMesh(
		MeshContainer3dSurface * smesh, 
		DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> & hbpoints,
		MeshDomainVolume* mdv)
{
	if(!m_surface || !m_mesh) return false;

	assert( smesh != nullptr );

	// create points
	int pct = m_mesh->getPointsCount();
	DataVector<MeshPoint3d*> points_3d(pct);

	smesh->addLocalSurface( m_surface );

	for(int i = 0; i < pct; i++){
		MeshPoint2d* mp = m_mesh->getPointAt(i);
		MeshPoint3d* mp3d = (MeshPoint3d*)mp->getPtrTag(TagExtended::TAG_MP_2D_3D);
		MeshPoint3d* mp3d_add = nullptr;
		if(mp3d != nullptr){
			mp3d_add = hbpoints.getValue( mp3d, nullptr );
			if( mp3d_add == nullptr ){
				mp3d_add = new MeshPoint3d( *mp3d );
				mp3d_add->copyBorderFlagsFrom( mp3d );
				hbpoints.insert( mp3d, mp3d_add );
				smesh->addMeshPoint( mp3d_add );
			}
		}else{
			mp3d_add = new MeshPoint3d(m_surface, mp->getCoordinates() );
			smesh->addMeshPoint( mp3d_add );
		}
		points_3d.add( mp3d_add );
	}

	// create faces
	int fct = m_mesh->getElementsCount();
	MeshBoundaryCondition* boundary_condition = 
		(MeshBoundaryCondition*)getPtrTag(TagExtended::TAG_BOUNDARY_COND);
	bool inverted = ( getBlock(1) == mdv );
	for(int i = 0; i < fct; i++){
		MeshElement* me = m_mesh->getElementAt(i);
		int edge_ct = me->getEdgeCount();
		MeshFace* mface = nullptr;
		MeshPoint2d* me_pts[4] = { nullptr };
		if(edge_ct == 3){
			me_pts[0] = me->getPoint(0);
			me_pts[1] = me->getPoint(inverted ? 2 : 1);
			me_pts[2] = me->getPoint(inverted ? 1 : 2);
			mface = new MeshTriangle3d(
				points_3d[ me_pts[0]->getIndex() ],
				points_3d[ me_pts[1]->getIndex() ],
				points_3d[ me_pts[2]->getIndex() ]);
		}else if(edge_ct == 4){
			me_pts[0] = me->getPoint(0);
			me_pts[1] = me->getPoint(inverted ? 3 : 1);
			me_pts[2] = me->getPoint(2);
			me_pts[3] = me->getPoint(inverted ? 1 : 3);
			mface = new MeshQuad3d(
				points_3d[ me_pts[0]->getIndex() ],
				points_3d[ me_pts[1]->getIndex() ],
				points_3d[ me_pts[2]->getIndex() ],
				points_3d[ me_pts[3]->getIndex() ]);
		}else{
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Incorrect mesh element, #edges: " << edge_ct);
			continue;
		}
		for (int j = 0; j < edge_ct; j++) {
			MeshEdge2d* edge_2d = me_pts[j]->getEdgeToPoint(me_pts[(j + 1) % edge_ct]);
			assert(edge_2d != nullptr);
			MeshEdge3d* edge_3d_b = (MeshEdge3d*)edge_2d->getPtrTag(TAG_ME_2D_3D);
			if (edge_3d_b != nullptr) {
				char bf = edge_3d_b->getBorderFlags();
				MeshEdge3d* edge_3d = mface->getPoint(j)->getEdgeToPoint(mface->getPoint((j+1)%edge_ct));
				assert(edge_3d);
				edge_3d->setBorderFlags(bf);
			}
		}

		mface->setBlockLink(blocks[0], inverted?1:0);
		mface->setBlockLink(blocks[1], inverted?0:1);
		if(boundary_condition)
			mface->setPtrTag(TagExtended::TAG_BOUNDARY_COND, boundary_condition);
		smesh->addMeshFace(mface);
	}

	return true;
}

/// classify free-point for this domain-surface
bool MeshDomainSurface::offerFreePoint(std::shared_ptr<MeshPoint3d> fpoint)
{
	const double TOL2 = mesh_data.relative_small_number * mesh_data.relative_small_number;
	// lying on the surface ?
	const DPoint3d& fpt = fpoint->getCoordinates();
	const DPoint2d param = m_surface->getParameters(fpt);
	DPoint3d rpt = m_surface->getPoint(param);
	if(rpt.distance2(fpt) > TOL2) return false;
	// ... yes!

	// check edges
	for(IteratorEdge2d it = m_boundary->getFirstEdge2d(); it.isValid(); it.nextEdge()){
		MeshEdge2d* edge = it.getEdge();
		double ksi = edge->getParameter(param);
		if(ksi < 0.0 || ksi > 1.0) continue;
		rpt = m_surface->getPoint(edge->getPoint(ksi));
		if(rpt.distance2(fpt) < TOL2){
			MeshDomainEdge3d* mde = (MeshDomainEdge3d*) edge->getPtrTag(TagExtended::TAG_ME_2D_3D);
			if(mde){
				mde->addFreePoint(fpoint);
				return true;
			}
		}
	}

	// TODO -> check if within the boundary...
	m_freepoints.add(fpoint);
	return true;
}

DVector3d MeshDomainSurface::getSurfaceNormalForMidEdge(const MeshPoint3d * mp0, const MeshPoint3d * mp1) const
{
	assert(m_boundary != nullptr);
	assert(m_surface);

	for (auto it = m_boundary->getFirstEdge2d(); it.isValid(); it.nextEdge()) {
		MeshEdge2d* edge = it.getEdge();
		MeshPoint2d* mp0_2d = edge->getMeshPoint(0);
		MeshPoint2d* mp1_2d = edge->getMeshPoint(1);
		MeshPoint3d* mp0_3d = (MeshPoint3d*)mp0_2d->getPtrTag(TagExtended::TAG_MP_2D_3D);
		MeshPoint3d* mp1_3d = (MeshPoint3d*)mp1_2d->getPtrTag(TagExtended::TAG_MP_2D_3D);
		if ((mp0_3d == mp0 && mp1_3d == mp1) || (mp0_3d == mp1 && mp1_3d == mp0)) {
			DPoint2d param = edge->getPoint(0.5);
			DVector3d vn = m_surface->getNormalVector(param);
			if (false) {
				MeshViewSet *set = getViewSet();
				DPoint3d pt3d = m_surface->getPoint(param);
				set->addPoint(pt3d, 3);
				set->addEdge(pt3d, pt3d + vn, 2);
				SHOW_MESH("getSurfaceNormalForMidEdge", set);
			}
			return vn;
		}
	}
	assert(false);
	return DVector3d::v_ox;
}

void MeshDomainSurface::autoSetBoundaryTags(MeshContainer3d * mesh)
{
	if (mesh == nullptr) return;

	for (auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()) {
		MeshEdge3d* edge = it.getEdge();
		if (!edge->isBorder()) continue;
		int fct = edge->getFaceCount();
		bool is_ridge = (fct != 2);
		if (!is_ridge) {
			MeshFace* f0 = edge->getFaceAt(0);
			MeshFace* f1 = edge->getFaceAt(1);
			DVector3d vn0 = f0->getSurfaceNormalForMidEdge(
				edge->getMeshPoint(0), edge->getMeshPoint(1));
			DVector3d vn1 = f1->getSurfaceNormalForMidEdge(
				edge->getMeshPoint(0), edge->getMeshPoint(1));

			double sp = vn0.scalarProduct(vn1);
			is_ridge = (sp < MeshGenerator3dSurface::param_sharp_edge_threshold);
		}
		if (is_ridge) {
			edge->setBorder(TagBorder::OUTER | TagBorder::RIDGE);
			edge->getMeshPoint(0)->setBorder(TagBorder::OUTER | TagBorder::RIDGE);
			edge->getMeshPoint(1)->setBorder(TagBorder::OUTER | TagBorder::RIDGE);
		}
	}
	int pct = mesh->getPointsCount();
	for (int i = 0; i < pct; i++) {
		MeshPoint3d* point = mesh->getPointAt(i);
		if (!point->isBorder() || point->isBorder(TagBorder::CORNER)) continue;
		int bect = point->getBorderEdgesCount();
		if (bect == 0) continue;
		if (!point->isBorder()) point->setBorderFlags(TagBorder::OUTER | TagBorder::RIDGE);
		if (bect != 2) point->setBorderFlags(TagBorder::CORNER);
	}
}

