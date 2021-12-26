/////////////////////////////////////////////////////////////////////////////
// ControlSpace2dMesh.cpp
// This class implements the control space (sizing and stretching info)
//	in the form of the triangular background mesh
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "ControlSpace2dMesh.h"
#include "ControlSpace2dIdentity.h"
#include "MeshContainer2d.h"
#include "MeshTriangle2d.h"
#include "MeshGenerator2d.h"
#include "MeshData.h"
#include "SurfaceParametric.h"
#include "EPSFile.h"
#include "Curve2dSegment.h"
#include "DTriangle.h"
#include "ControlSpace3dAdaptive.h"

/// Whether the calculated weights for the metric value inside the triangle should be squared
int ControlSpace2dMesh::param_weights_squared = 1;

/// Whether inner nodes should be inserted during preparation of this control mesh: 
///		-1 - unlimited, 0 - none, n - maximum n nodes
int ControlSpace2dMesh::param_inner_nodes;

/// Interpolation method for control mesh
int ControlSpace2dMesh::param_interpolation_method = MeshData::CONTROL_TRIANGLE;

/// Threshold for metric difference for assesing mesh-control interpolation type
double ControlSpace2dMesh::param_mixed_threshold = 0.1;

ControlSpace2dMesh::ControlSpace2dMesh(SurfaceConstPtr surface, const DRect& box, int approx) 
	: ControlSpace2dAdaptive(surface, box), m_control_mesh(nullptr), 
		m_control_points(nullptr), m_approximation_method(approx),
		m_nodes(100)
{
	LOG4CPLUS_DEBUG(MeshLog::logger_console, "Creating mesh control space, interpolation type " << approx);
}

ControlSpace2dMesh::~ControlSpace2dMesh()
{
	if(m_control_mesh){
		delete m_control_mesh;
	}
	if(m_control_points){
		delete m_control_points;
	}
}

/// Returns number of control nodes in adaptive control structure
int ControlSpace2dMesh::getControlNodesCount() const
{
	return m_nodes.countInt();
}

/// Invoke function for all control nodes of this space (read-only)
void ControlSpace2dMesh::forEachControlNode(const std::function<void(const ControlNode2d & node)>& fg) const {
	for (int i = getControlNodesCount() - 1; i >= 0; i--)
		fg(m_nodes[i]);
}

/// Invoke function for all control nodes of this space
void ControlSpace2dMesh::forEachControlNode(const std::function<void(ControlNode2d & node)>& fg) {
	for (int i = getControlNodesCount() - 1; i >= 0; i--)
		fg(m_nodes[i]);
}

/////////////////////////////////////////////////////////////////
// Creating the background mesh for the given set of points (triangulation, ?optimization, etc...)
bool ControlSpace2dMesh::interpolate()
{
	if(!m_control_points) return false;
	int pct = m_control_points->getPointsCount();
	if(pct < 1) return false;

	// triangulate
	if(m_control_mesh){
		delete m_control_mesh;
	}

	// -> initial mesh
	m_control_mesh = MeshGenerator2d::createInitialMesh(m_control_points, &m_box);
	CS2dPtr space(new ControlSpace2dIdentity(m_min_length));
	m_control_mesh->setControlSpace(space);
	Metric2dContext mc(space);

//----------------------------------------------
//	if(MeshView::isRunning())
//		MeshView::setMeshToUpdate(m_control_mesh);
//	else
//		_beginthread(MeshView::startLoopGL2D, 0, m_control_mesh);
//	cout << "Press enter...";
//	string text;
//	cin >> text;
//----------------------------------------------
	
	START_CLOCK("Control mesh interpolation");

	// -> add points
	srand(0);
	do{
		MeshPoint2d *point = m_control_points->removeMeshPoint(rand() % pct);
		MeshTriangle2d* start_triangle = m_control_mesh->getNearTriangle(point->getCoordinates());
		assert(start_triangle);
		MeshTriangle2d* containing_triangle = start_triangle->findTriangleByNeighbours(point);
		if(!containing_triangle){
			LOG4CPLUS_WARN(MeshLog::logger_console, 
				"Error finding the containing triangle -> switching to linear search (for this point only) ...");
			int tct = m_control_mesh->getElementsCount();
			for(int i = 0; i < tct; i++){
				MeshTriangle2d* triangle = (MeshTriangle2d*)m_control_mesh->getElementAt(i);
				if(triangle->isPointInside(point->getCoordinates())){
					containing_triangle = triangle;
					break;
				}
			}
		}
		if(!containing_triangle){
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"Error creating (triangulating) the background control mesh.");
			return false;
		}
		bool valid_for_insertion = true;
		for(int i = 0; i < 3; i++){
			MeshPoint2d* point1 = containing_triangle->getPoint(i);
			if(point->getMetricCoordinates(mc).distance2(
				point1->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER){
				LOG4CPLUS_DEBUG(MeshLog::logger_console, "Control points collision - calculating minimum");
				ControlNode2d* data1 = (ControlNode2d*)point1->getPtrTag(TagExtended::TAG_ACS_METRIC);
				if(data1){
					ControlNode2d* data = (ControlNode2d*)point->getPtrTag(TagExtended::TAG_ACS_METRIC);
					if(data){
						data1->control_data.setMinimum(data->control_data);
						data->control_data = data1->control_data;
						delete point;
					}else{
						delete point;
					}
				}else{
					point1->setPtrTag(TagExtended::TAG_ACS_METRIC, point->getPtrTag(TagExtended::TAG_ACS_METRIC));
					delete point;
				}
				valid_for_insertion = false;
			}
		}
		if(valid_for_insertion){
			if(!MeshGenerator2d::addPointToTriangulation(mc, m_control_mesh, point, containing_triangle)){
				LOG4CPLUS_ERROR(MeshLog::logger_console, 
					"Error creating (triangulating) the background control mesh.");
				return false;
			}
		}
//----------------------------------------------
//		MeshView::setMeshToUpdate(m_control_mesh);
//		cout << "Press enter...";
//		cin >> text;
//----------------------------------------------
	}while(--pct > 0);
//----------------------------------------------
//		cout << "Finished. Press enter...";
//		cin >> text;
//----------------------------------------------

	// -> clear points container (should be empty by now)
	assert(m_control_points->getPointsCount() == 0);
	delete m_control_points;
	m_control_points = nullptr;

	if(param_inner_nodes){
//----------------------------------------------
//		if(MeshView::isRunning())
//			MeshView::setMeshToUpdate(m_control_mesh);
//		else
//			_beginthread(MeshView::startLoopGL2D, 0, m_control_mesh);
//		cout << "Press enter...";
//		string text;
//		cin >> text;
//----------------------------------------------
		//------ find shortest edge;
		double shortest_length = mesh_data.relative_infinity;
		int ipct = m_control_mesh->getPointsCount();
		for(int i = 0; i < ipct; i++){
			MeshPoint2d* point = m_control_mesh->getPointAt(i);
			ControlNode2d* data = (ControlNode2d*)point->getPtrTag(TagExtended::TAG_ACS_METRIC);
			if(data){
				ControlDataStretch2d s = DMetric2d::matrixToStretch(data->control_data);
				double len = std::min(s.lx, s.ly);
				if(len < shortest_length) shortest_length = len;
			}
		}
		//------
		int quality_criterion = MeshTriangle2d::param_quality_criterion;	// temporarily change criterion
		MeshTriangle2d::param_quality_criterion = MeshData::QUALITY_QUOTIENT_AREAS;
		//------
		int count = addInnerNodes(mc, param_inner_nodes);
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Inner nodes for control mesh " << count);
		//------
		MeshTriangle2d::param_quality_criterion = quality_criterion;	// restore default
//----------------------------------------------
//		MeshView::setMeshToUpdate(m_control_mesh);
//		cout << "Press enter...";
//		cin >> text;
//----------------------------------------------
	}

	if(m_approximation_method == MeshData::CONTROL_MIXED)
		checkControlMeshForInterpolationType();

	STOP_CLOCK("Control mesh interpolation");

	m_initialized = 1;
	return true;
}

/// Adds new information (stretching and lengths) for some point within the domain
void ControlSpace2dMesh::addControlNode(const ControlNode2d& node)
{
	assert(m_initialized == 0);
	assert(node.control_data.det() >= 0.0);

	// Create points container, if not yet ready
	if(!m_control_points){
		m_control_points = new MeshContainer2d(100);
		m_min_length = node.control_data.minEigenvalue();
	}else{
		m_min_length = std::min(m_min_length, node.control_data.minEigenvalue());
	}

	// New point with data
	MeshPoint2d* mesh_point = new MeshPoint2d(node.coord);
	int i = m_nodes.add(node);
	mesh_point->setPtrTag(TagExtended::TAG_ACS_METRIC, &(m_nodes[i]));
	m_control_points->addMeshPoint(mesh_point);
}

double ControlSpace2dMesh::getMetricGradationRatio(const DPoint2d& pt_fit) const
{
	MeshTriangle2d* near_triangle = m_control_mesh->getNearTriangle(pt_fit);
	assert(near_triangle != nullptr);
	MeshTriangle2d* containing_triangle = near_triangle->findTriangleByNeighbours(pt_fit);

	if(!containing_triangle){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"setCurrentMetricAtPoint: error finding the containing triangle "
			<< "switching to linear triangle search (for this point only) ...");
		int tct = m_control_mesh->getElementsCount();
		for(int i = 0; i < tct; i++){
			MeshTriangle2d* triangle = (MeshTriangle2d*)m_control_mesh->getElementAt(i);
			if(triangle->isPointInside(pt_fit)){
				containing_triangle = triangle;
				break;
			}
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_console, " -> trying to find nearest point? ...");
		int pct = m_control_mesh->getPointsCount();
		double min_distance2 = LARGE_NUMBER;
		ControlNode2d* data = 0;
		for(int i = 0; i < pct; i++){
			MeshPoint2d* point = m_control_mesh->getPointAt(i);
			double distance2 = pt_fit.distance2(point->getCoordinates());
			if(distance2 < min_distance2){ 
				void* v_data = point->getPtrTag(TagExtended::TAG_ACS_METRIC);
				if(v_data){
					data = (ControlNode2d*)v_data;
					min_distance2 = distance2;
				}
			}
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_console, 
			"Best point found at the distance: " << sqrt(min_distance2));
		return data ? data->max_gradation_ratio : 1;
	}
	LOG_ASSERT(containing_triangle != nullptr);

	double tw = 0;
	double ave = 0;
	for(int i = 0; i < 3; i++){
		const MeshPoint2d* mpt = containing_triangle->getPoint(i);
		const ControlNode2d* node = (ControlNode2d*)mpt->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(node == 0) continue;
		double dist = node->coord.distance2(mpt->getCoordinates());
		if(dist < mesh_data.relative_small_number)
			return node->max_gradation_ratio;
		double w = 1.0/dist;
		ave += node->max_gradation_ratio * w;
		tw += w;
	}
	if(tw > 0.0) return ave / tw;
	else return 1.0;
}

/////////////////////////////////////////////////////////////////////
// Funkcja ustala úredniπ metrykÍ dla zadanego punktu
ControlDataMatrix2d ControlSpace2dMesh::getMetricAtPoint(const DPoint2d& pt_fit, int method, ControlDataMatrix2d * pcdm) const
{
	MeshTriangle2d* near_triangle = m_control_mesh->getNearTriangle(pt_fit);
	LOG_ASSERT(near_triangle != nullptr);
	MeshTriangle2d* containing_triangle = near_triangle->findTriangleByNeighbours(pt_fit);

	if(!containing_triangle){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"setCurrentMetricAtPoint: error finding the containing triangle"
				<< " -> switching to linear triangle search (for this point only) ...");
		int tct = m_control_mesh->getElementsCount();
		for(int i = 0; i < tct; i++){
			MeshTriangle2d* triangle = (MeshTriangle2d*)m_control_mesh->getElementAt(i);
			if(triangle->isPointInside(pt_fit)){
				containing_triangle = triangle;
				break;
			}
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_console, 
			" -> trying to find nearest point? ...");
		int pct = m_control_mesh->getPointsCount();
		double min_distance2 = mesh_data.relative_infinity;
		ControlNode2d* data = 0;
		for(int i = 0; i < pct; i++){
			MeshPoint2d* point = m_control_mesh->getPointAt(i);
			double distance2 = pt_fit.distance2(point->getCoordinates());
			if(distance2 < min_distance2){ 
				void* v_data = point->getPtrTag(TagExtended::TAG_ACS_METRIC);
				if(v_data){
					data = (ControlNode2d*)v_data;
					min_distance2 = distance2;
				}
			}
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Best point found at the distance: " << sqrt(min_distance2));
		return data->control_data;
	}
	LOG_ASSERT(containing_triangle != nullptr);

	if(pcdm){
		MeshPoint2d* mpt = containing_triangle->getPoint(0);
		ControlNode2d* node = (ControlNode2d*)mpt->getPtrTag(TagExtended::TAG_ACS_METRIC);
		double dist2 = pt_fit.distance2(mpt->getCoordinates());
		double nearest_dist2 = dist2;
		if(!node){
			node = (ControlNode2d*)containing_triangle->getPoint(1)->getPtrTag(TagExtended::TAG_ACS_METRIC);
		}else{
			dist2 = pt_fit.distance2(containing_triangle->getPoint(1)->getCoordinates());
			if(dist2 > nearest_dist2){
				nearest_dist2 = dist2;
			}
		}
		if(!node){
			node = (ControlNode2d*)containing_triangle->getPoint(2)->getPtrTag(TagExtended::TAG_ACS_METRIC);
		}else{
			dist2 = pt_fit.distance2(containing_triangle->getPoint(2)->getCoordinates());
		}
		if(node){
			*pcdm = node->param_data;
		}else{
			double g_ratio;
			node->param_data = base_surface->countParameterizationMatrix(node->coord, g_ratio);
		}
	}

	if(method < 0) method = m_approximation_method;
	if(method == MeshData::CONTROL_MIXED)
		method = containing_triangle->getIntTag(TagExtended::TAG_ACS_METRIC);

	switch(method){
	case MeshData::CONTROL_SIMPLE:
		return getMetricSimple(pt_fit, containing_triangle);
	case MeshData::CONTROL_TRIANGLE:
		return getMetricInterpolate(pt_fit, containing_triangle, false, param_weights_squared==1);
	case MeshData::CONTROL_VORONOI:
		return getMetricInterpolate(pt_fit, containing_triangle, true, param_weights_squared==1);
//	case MeshData::CONTROL_NEM_LAYERED:
//	     return getMetricInterpolateNEM(pt_fit, containing_triangle);
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console, 
			"Unknown approximation method (" << m_approximation_method << ") for control mesh! ");
		return getMetricSimple(pt_fit, containing_triangle);
	}
}

/// Return interpolated value of extended tag data from control nodes at some given point;
double ControlSpace2dMesh::interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const
{
	const DPoint2d pt_fit = m_box.fitInPoint(pt);

	MeshTriangle2d* near_triangle = m_control_mesh->getNearTriangle(pt_fit);
	assert(near_triangle != nullptr);
	MeshTriangle2d* containing_triangle = near_triangle->findTriangleByNeighbours(pt_fit);

	if(!containing_triangle){
		int tct = m_control_mesh->getElementsCount();
		for(int i = 0; i < tct; i++){
			MeshTriangle2d* triangle = (MeshTriangle2d*)m_control_mesh->getElementAt(i);
			if(triangle->isPointInside(pt_fit)){
				containing_triangle = triangle;
				break;
			}
		}
		int pct = m_control_mesh->getPointsCount();
		double min_distance2 = mesh_data.relative_infinity;
		ControlNode2d* data = 0;
		for(int i = 0; i < pct; i++){
			MeshPoint2d* point = m_control_mesh->getPointAt(i);
			double distance2 = pt_fit.distance2(point->getCoordinates());
			if(distance2 < min_distance2){ 
				void* v_data = point->getPtrTag(TagExtended::TAG_ACS_METRIC);
				if(v_data){
					data = (ControlNode2d*)v_data;
					min_distance2 = distance2;
				}
			}
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_console, 
			"Best point found at the distance: " << sqrt(min_distance2));
		return data->getDoubleTag(type);
	}
	LOG_ASSERT(containing_triangle != nullptr);

	// use MeshData::CONTROL_TRIANGLE
	Metric2dContext mc(m_control_mesh->getControlSpace());
	const DMPoint2d dpt = mc.transformPStoMS(pt);
	const DPoint2d pt_m = mc.transformPStoRS(pt);

	double total_weight = 0.0;
	double ave = 0.0;

	for(int i = 0; i < 3; i++){
		ControlNode2d* node = (ControlNode2d*)containing_triangle->getPoint(i)->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(node == 0) continue;
		double dist = pt_m.distance(mc.transformPStoRS(node->coord));
		if(dist < mesh_data.relative_small_number)
			return node->getDoubleTag(type);
		double w = 1.0 / (dist*dist);
		total_weight += w;
		ave += node->getDoubleTag(type) * w;
	}
	if(total_weight > 0.0)
		ave /= total_weight;

	return ave;
}

ControlDataMatrix2d ControlSpace2dMesh::getMetricSimple(const DPoint2d& /* pt */, MeshTriangle2d* triangle) const
{
	LOG_ASSERT(triangle != nullptr);
	for(int i = 0; i < 3; i++){
		MeshPoint2d* mpt = triangle->getPoint(i);
		ControlNode2d* data = (ControlNode2d*)mpt->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(data) return data->control_data;
	}
	LOG4CPLUS_ERROR(MeshLog::logger_console, 
		"Simple approximation for control mesh - no valid vertex!");
	return ControlDataMatrix2d();
}

ControlDataMatrix2d ControlSpace2dMesh::getMetricInterpolate(const DPoint2d& pt, 
	MeshTriangle2d* triangle, bool use_voronoi_points, bool square_weight) const
{
	LOG_ASSERT(triangle != nullptr);
	DataVector<MeshPoint2d*> voronoi_points;

	voronoi_points.add(triangle->getPoint(0));
	voronoi_points.add(triangle->getPoint(1));
	voronoi_points.add(triangle->getPoint(2));

	Metric2dContext mc(m_control_mesh->getControlSpace());
	const DMPoint2d dpt = mc.transformPStoMS(pt);

	if(use_voronoi_points){
		DataVector<MeshTriangle2d*> marked_triangles;
		marked_triangles.add(triangle);

		// Find all outer-circle-containing triangles for this points
		//	and theirs vertices
		for(size_t i = 0; i < marked_triangles.countInt(); i++){
			MeshTriangle2d* mtriangle = marked_triangles[i];
			for(int j = 0; j < 3; j++){
				MeshTriangle2d* triangle1 = (MeshTriangle2d*)mtriangle->getNeighbour(j);
				if(!triangle1) continue;
				bool in_circle = (DTriangle2d::inSphereCheck(dpt, 
					triangle1->getPoint(0)->getMetricCoordinates(mc), 
					triangle1->getPoint(1)->getMetricCoordinates(mc), 
					triangle1->getPoint(2)->getMetricCoordinates(mc)) > 0.0);
				if(in_circle){
					if(marked_triangles.addIfNew(triangle1)){
						voronoi_points.addIfNew(triangle1->getPoint(0));
						voronoi_points.addIfNew(triangle1->getPoint(1));
						voronoi_points.addIfNew(triangle1->getPoint(2));
					}
				}
			}
		}
	}

	double total_weight = 0.0;
	int vpts_ct = (int)voronoi_points.countInt();

	ControlDataMatrix2d data;
	const DPoint2d pt_m = mc.transformPStoRS(pt);
	for(int i = 0; i < vpts_ct; i++){
		ControlNode2d* node = (ControlNode2d*)voronoi_points.get(i)->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(node == 0) continue;
		double dist = pt_m.distance(mc.transformPStoRS(node->coord));
		if(dist < mesh_data.relative_small_number)
			return node->control_data;
		double w = 1.0 / (square_weight?(dist*dist):dist);
		total_weight += w;
		data += node->control_data * w;
	}
	if(total_weight > 0.0)
		data /= total_weight;
	else
		LOG4CPLUS_ERROR(MeshLog::logger_console, 
			"setCurrentMetricInterpolate: none of points has valid metric.");

	return data;
}

/*
ControlDataMatrix2d ControlSpace2dMesh::getMetricInterpolateNEM(const DPoint2d& pt, MeshTriangle2d* triangle)
{
	assert(triangle);
	DataVector<MeshPoint2d*> triangle_points;

	triangle_points.add(triangle->getPoint(0));
	triangle_points.add(triangle->getPoint(1));
	triangle_points.add(triangle->getPoint(2));

	ControlDataMatrix2d data1;
	ControlDataMatrix2d data2;
	int count1=0, count2=0;

	for(int i = 0; i < 3; i++){
		MeshPoint2d* point = triangle->getPoint(i);
		ControlDataExt *ext_data = (ControlDataExt*)point->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(ext_data){
			double dist = pt.distance(point->getCoordinates());
			if(dist < ext_data->radius1){
				if(count1++)
					data1.counstd::minimum(ext_data->data1);
				else
					data1 = ext_data->data1;
			}else if(!count1 && dist < (ext_data->radius1+ext_data->radius2)){
				if(count2++)
					data2.counstd::minimum(ext_data->getControlDataMatrix(dist));
				else
					data2 = ext_data->getControlDataMatrix(dist);
			}
		}
	}

	if(count1){
//		LOG4CPLUS_INFO(MeshLog::logger_console, "1 zone");
		return data1;
	}else if(count2){
//		LOG4CPLUS_INFO(MeshLog::logger_console, "2 zone");
		return data2;
	}

	DataVector<MeshPoint2d*> layer_points;
	for(int i = 0; i < 3; i++)
		layer_points.add(triangle_points[i]);

	// Find two layers of nodes around the containing triangle
	// 1st and 2nd layer:
	for(int k = 0; k < 2; k++){
		int lpcount = layer_points.countInt();
		for(int i =3*k; i < lpcount; i++){ // 0 for 1st layer, 3 for 2nd layer
			MeshPoint2d* point = layer_points[i];
			int rank = point->getRank();
			for(int j = 0; j < rank; j++){
			    MeshPoint2d* other_point = point->getEdge(j)->getOtherPoint(point);
				layer_points.addIfNew(other_point);
			}
		}
	}
 	// .........
	for(int i = 0; i < vpts_ct; i++){
		ControlDataExt* ext_data = (ControlDataExt*)voronoi_points.get(i)->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(ext_data == 0) continue;
		double dist = pt.distance(voronoi_points.get(i)->getCoordinates());
		if(dist < mesh_data.relative_small_number)
			return ext_data->getControlDataMatrix(dist);
		double w = 1.0 / (square_weight?(dist*dist):dist);
		total_weight += w;
		data1 += ext_data->getControlDataMatrix(dist) * w;
	}
	if(total_weight > 0.0)
		data1 /= total_weight;
	else
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "setCurrentMetricInterpolate: none of points has valid metric.");
	return data1;
}
*/

int ControlSpace2dMesh::addInnerNodes(Metric2dContext& mc, int max_ct)
{
	LOG_ASSERT(m_control_mesh != nullptr);
	if(!m_control_mesh) return 0;

	mc.countMetricAtPoint(DPoint2d(0.0, 0.0));	// is Identity anyhow

	int border_count = m_control_mesh->getPointsCount();
	
	// Wyznacz jakoúÊ dostÍpnych trÛjkπtÛw
	int count = m_control_mesh->getElementsCount();
	if(count < 1) return 0;

	for(int i = 0; i < count; i++) 
		m_control_mesh->getElementAt(i)->countQuality(mc, true);

	m_control_mesh->setHeapOrder(true);
	int new_ct = 0;
	int steps = 0;
	double threshold = 0.4;

	while(true){
		// Najgorszy trÛjkπt z uwzglÍdnieniem minimalnej powierzchni trÛjkπtÛw znajduje siÍ w korzeniu
		MeshTriangle2d *triangle = (MeshTriangle2d*)m_control_mesh->getElementAt(0);

		// Czy jakoúÊ jest wystarczajπco niska
		if(triangle->getQuality() > threshold) break;

		++steps;

		DPoint2d dpoints[3];
		for(int i = 0; i < 3; i++)
			dpoints[i] = triangle->getPoint(i)->getCoordinates();
		// Znajdü najd≥uøszπ krawÍdü
		int imax = -1;
		double d2, dmax = 0.0;
		for(int i = 0; i < 3; i++){
			if((d2 = dpoints[i].distance2(dpoints[(i+1)%3])) > dmax){
				dmax = d2;
				imax = i;
			}
		}
		// Nowy punkt w úrodku tej krawÍdzi
		MeshEdge2d* edge = triangle->getEdge(imax);
		MeshPoint2d* bpt1 = edge->getMeshPoint(0);
		MeshPoint2d* bpt2 = edge->getMeshPoint(1);
		bool was_border = false;
		// Jeúli punkt znajduje siÍ na krawÍdzi -> anuluj wstawianie
		if(edge->isBorder()){
			edge->clearBorder();
			was_border = true;
		}
		DPoint2d dnew = edge->getPoint(0.5);

		MeshPoint2d* p3 = new MeshPoint2d(dnew);
		int method = m_approximation_method;
		if(method == MeshData::CONTROL_MIXED) method = MeshData::CONTROL_VORONOI;
		ControlDataMatrix2d data_matrix = getMetricAtPoint(dnew, method);
		if(!MeshGenerator2d::addPointToTriangulation(mc, m_control_mesh, p3, triangle)){
			delete p3;
			return new_ct;
		}

		size_t i = m_nodes.add(ControlNode2d(dnew, data_matrix));
		p3->setPtrTag(TagExtended::TAG_ACS_METRIC, &(m_nodes[i]));
		if(was_border){
			MeshEdge2d* bedge = p3->getEdgeToPoint(bpt1);
			if(bedge) bedge->setBorder();
			bedge = p3->getEdgeToPoint(bpt2);
			if(bedge) bedge->setBorder();
		}

		++new_ct;
		count = border_count + new_ct;
		if(max_ct > 0 && new_ct >= max_ct)
			return new_ct;
	}
		
	return new_ct;
}

void ControlSpace2dMesh::checkControlMeshForInterpolationType()
{
	int tct = m_control_mesh->getElementsCount();
	int counter[3] = {0,0,0};
	for(int i = 0; i < tct; i++){
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "------------------------");
		MeshTriangle2d* triangle = (MeshTriangle2d*)m_control_mesh->getElementAt(i);
		DPoint2d middle = triangle->getMiddlePoint();
		int tm = getRequiredInterpolationType(middle, triangle, param_mixed_threshold);
		if(tm == 2)	{ triangle->setIntTag(TagExtended::TAG_ACS_METRIC, 2); counter[2]++; continue; }
		tm = getRequiredInterpolationType(triangle->getEdge(0)->getPoint(0.5), triangle, param_mixed_threshold);
		if(tm == 2)	{ triangle->setIntTag(TagExtended::TAG_ACS_METRIC, 2); counter[2]++; continue; }
		tm = getRequiredInterpolationType(triangle->getEdge(1)->getPoint(0.5), triangle, param_mixed_threshold);
		if(tm == 2)	{ triangle->setIntTag(TagExtended::TAG_ACS_METRIC, 2); counter[2]++; continue; }
		tm = getRequiredInterpolationType(triangle->getEdge(2)->getPoint(0.5), triangle, param_mixed_threshold);

		if(tm == 2)	{ triangle->setIntTag(TagExtended::TAG_ACS_METRIC, 2); counter[2]++; continue; }

		ControlDataMatrix2d m = getMetricInterpolate(middle, triangle, true, param_weights_squared==1);
		double max_diff = 0.0;
		for(int j = 0; j < 3; j++){
			ControlNode2d* m0 = (ControlNode2d*)triangle->getPoint(j)->getPtrTag(TagExtended::TAG_ACS_METRIC);
			if(m0){
				double diff = m.countDifferenceRR(m0->control_data);
				if(diff > max_diff) max_diff = diff;
			}
		}
		if(max_diff < param_mixed_threshold)
			tm = 0;

//		if(tm == 0) tm = 1;
		triangle->setIntTag(TagExtended::TAG_ACS_METRIC, tm); 
		counter[tm]++;
	}
	LOG4CPLUS_INFO(MeshLog::logger_console, "simple interpolation: " << counter[0]);
	LOG4CPLUS_INFO(MeshLog::logger_console, "linear interpolation: " << counter[1]);
	LOG4CPLUS_INFO(MeshLog::logger_console, "voronoi interpolation: " << counter[2]);
}

int ControlSpace2dMesh::getRequiredInterpolationType(const DPoint2d& pt, MeshTriangle2d* triangle, double threshold)
{
	ControlDataMatrix2d mv = getMetricInterpolate(pt, triangle, true, param_weights_squared==1);
	ControlDataMatrix2d m = getMetricInterpolate(pt, triangle, false, param_weights_squared==1);
	if(mv.countDifferenceRR(m) < threshold) return MeshData::CONTROL_TRIANGLE;
	else return MeshData::CONTROL_VORONOI;
}

void ControlSpace2dMesh::storeEPS(const char* name, int id)
{
	if(!m_control_mesh) return;
	m_control_mesh->storeEPS(name, id);
}

int ControlSpace2dMesh::getTrianglesCount() const { 
	return m_control_mesh->getElementsCount(); 
}

bool ControlSpace2dMesh::countRegularGridPoints(int count, DataVector<DPoint2d> & points) const
{
	if(!base_surface) return false;

	// length x
	Curve2dSegment middle_xline(DPoint2d(m_box.x0, (m_box.y0+m_box.y1)*0.5), 
		DPoint2d(m_box.x1, (m_box.y0+m_box.y1)*0.5));
	double real_dx = middle_xline.getLengthOnSurface(0.0, 1.0, base_surface);
	// length y
	Curve2dSegment middle_yline(DPoint2d((m_box.x0+m_box.x1)*0.5, m_box.y0), 
		DPoint2d((m_box.x0+m_box.x1)*0.5, m_box.y1));
	double real_dy = middle_yline.getLengthOnSurface(0.0, 1.0, base_surface);

	double ratio = real_dx / real_dy;	// adjust m_nx:m_ny ratio to real lengths of domain
	int nx = std::max(2,(int)(sqrt(ratio*count)));
	int ny = std::max(2, (int)(nx/ratio));
	nx = std::min(nx, count/2);
	ny = std::min(ny, count/2);

	double dx = m_box.getDX()  / (nx-1);
	double dy = m_box.getDY() / (ny-1);

	for(int j = 0; j < ny; j++){
		double y = m_box.y0 + j*dy;
		if(j%2 == 0)
			for(int i = 0; i < nx; i++)
				points.add(DPoint2d(m_box.x0 + i*dx, y));
		else
			for(int i = 1; i < nx; i++)
				points.add(DPoint2d(m_box.x0 - 0.5*dx + i*dx, y));
	}
	return true;
}

void ControlSpace2dMesh::setSurfaceCurvatureControlData()
{
	DataVector<DPoint2d> points;
	countRegularGridPoints(param_control_nxy, points);

	double p = m_box.getDiameter();
	double max_len = p * std::min(param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
	double min_len = std::max(p * param_min_diameter_ratio, param_min_length);

	// count for regular initial set of points
	for(size_t i = 0; i < points.countInt(); i++){
		double g_ratio;
		const SurfaceCurvature curvature = base_surface->getCurvature(points[i], g_ratio);
		if(g_ratio > MIN_PARAM_GRATIO){
			const ControlDataStretch2d data = adjustCurvatureData(curvature, param_curvature_ratio,
				p, min_len, max_len, param_stretch_max_ratio);
			addControlPoint(points[i], DMetric2d::stretchToMatrix(data));
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Improper parameterization " << g_ratio);
		}
	}

	interpolate();

	Metric2dContext mc(m_control_mesh->getControlSpace());

	// adapt
	int ect = m_control_mesh->getElementsCount();
	for(int j = 0; j < 10; j++){
		bool no_change = true;
		int last_pct = m_control_mesh->getPointsCount();
		for(int i = 0; i < ect; i++){
			MeshTriangle2d* triangle = (MeshTriangle2d*)m_control_mesh->getElementAt(i);
			assert(triangle->getEdgeCount() == 3);
			if(triangle->getPoint(0)->isBorder() || triangle->getPoint(1)->isBorder() || 
				triangle->getPoint(2)->isBorder()) continue; // skip boundary ones
			if(triangle->getPoint(0)->getIndex() >= last_pct || 
				triangle->getPoint(1)->getIndex() >= last_pct || 
				triangle->getPoint(2)->getIndex() >= last_pct) continue; // limit number of passes
			const DPoint2d middle = triangle->getMiddlePoint();

			// compare 
			double g_ratio;
			const SurfaceCurvature curvature = base_surface->getCurvature(middle, g_ratio);
			if(g_ratio > MIN_PARAM_GRATIO){
				const ControlDataStretch2d data = adjustCurvatureData(curvature, param_curvature_ratio,
					p, min_len, max_len, param_stretch_max_ratio);
				//double min_dist = middle.distance(triangle->getPoint(0)->getCoordinates());
				//min_dist = std::min(min_dist, middle.distance(triangle->getPoint(1)->getCoordinates()));
				//min_dist = std::min(min_dist, middle.distance(triangle->getPoint(2)->getCoordinates()));
				//if(min_dist < std::min(data.lx, data.ly)) continue;
				const ControlDataMatrix2d cdm = DMetric2d::stretchToMatrix(data);
				const ControlDataMatrix2d local_cdm = getMetricAtPoint(middle);
				double diff = local_cdm.countDifferenceRR(cdm);
				if(diff > param_threshold_diff){
					// update min control value
					double dmin = cdm.minEigenvalue();
					ControlSpace2dIdentity* space = 
						(ControlSpace2dIdentity*)(m_control_mesh->getControlSpace().get());
					if(dmin < space->getMinLength()) space->setMinLength(dmin);
					// insert new point
					MeshPoint2d* mesh_point = new MeshPoint2d(middle);
					bool result = MeshGenerator2d::addPointToTriangulation(mc, 
						m_control_mesh, mesh_point, triangle);
					if(result == false){ 
						delete mesh_point; 
						LOG4CPLUS_INFO(MeshLog::logger_console, 
							"error inserting point into triangulation");
					}else{
						size_t id = m_nodes.add(ControlNode2d(middle, cdm));
						mesh_point->setPtrTag(TagExtended::TAG_ACS_METRIC, &(m_nodes[id]));
						ect = m_control_mesh->getElementsCount();
						no_change = false;
					}
				}
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "Improper parameterization " << g_ratio);
			}
		}
		if(no_change) break;
	}
}

MeshViewSet* ControlSpace2dMesh::getViewSet(MeshViewSet* set, bool with_surface, TagExtended::TagType /* tag */) const 
{ 
	if(m_control_mesh != nullptr)
		return m_control_mesh->getViewSet(set, with_surface);
	else 
		return nullptr; 
}

/// Log basic information about this control space
void ControlSpace2dMesh::logDescription() const
{
	LOG4CPLUS_INFO(MeshLog::logger_mesh,
		"CS-Mesh, control nodes = " << getControlNodesCount() <<
		", vertices = " << (m_control_mesh?m_control_mesh->getPointsCount():0) << 
		", triangles = " << (m_control_mesh?m_control_mesh->getElementsCount():0));
}

/// Refines control space at the given point
bool ControlSpace2dMesh::setMinControl(const DPoint2d& pt, const ControlDataMatrix2d& cdm, bool min_value_set)
{
	if(!m_box.contains(pt)) return false;
	assert(m_initialized > 0);
	LOG_ASSERT(m_control_mesh != nullptr);
	MeshTriangle2d* near_triangle = m_control_mesh->getNearTriangle(pt);
	LOG_ASSERT(near_triangle != nullptr);
	MeshTriangle2d* containing_triangle = near_triangle->findTriangleByNeighbours(pt);
	if(!containing_triangle){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"ControlSpace2dMesh: error finding the containing triangle"
				<< " -> switching to linear triangle search (for this point only) ...");
		int tct = m_control_mesh->getElementsCount();
		for(int i = 0; i < tct; i++){
			MeshTriangle2d* triangle = (MeshTriangle2d*)m_control_mesh->getElementAt(i);
			if(triangle->isPointInside(pt)){
				containing_triangle = triangle;
				break;
			}
		}
		if(!containing_triangle) return false;
	}
	const ControlDataMatrix2d local_cdm = getMetricAtPoint(pt);
	double diff = local_cdm.countDifferenceRR(cdm);
	bool any_changes = false;
	if(diff > param_threshold_diff){
		// split may be unnecessary, if metric is greater than already set for this leaf...
		ControlDataMatrix2d min_cdm = cdm;
		min_cdm.setMinimum(local_cdm);
		if(local_cdm.countDifferenceRR(min_cdm) < param_threshold_diff){
			// different from current metric, but introduces no changes - so skip
			return false;
		}

		double dmin = cdm.minEigenvalue();
		ControlSpace2dIdentity* space = 
			(ControlSpace2dIdentity*)(m_control_mesh->getControlSpace().get());
		if(dmin < space->getMinLength()) space->setMinLength(dmin);

		DMetric2d dmp(cdm, base_surface, pt);
		double min_len = dmp.transformPStoMS(pt - containing_triangle->getPoint(0)->getCoordinates()).length2();
		min_len = std::min(min_len, dmp.transformPStoMS(pt - containing_triangle->getPoint(1)->getCoordinates()).length2());
		min_len = std::min(min_len, dmp.transformPStoMS(pt - containing_triangle->getPoint(2)->getCoordinates()).length2());
		if(min_len > METRIC_LENGTH_RATIO2){
			MeshPoint2d* mesh_point = insertNodeAndRefine(pt, containing_triangle);
			LOG_ASSERT(mesh_point);
			if(!mesh_point)	LOG4CPLUS_WARN(MeshLog::logger_console, "error inserting point into triangulation");
			else{
				double g_ratio;
				const ControlDataMatrix2d pdm = base_surface->countParameterizationMatrix(pt, g_ratio);
				size_t id = m_nodes.add(ControlNode2d(pt, cdm, pdm));
				mesh_point->setPtrTag(TagExtended::TAG_ACS_METRIC, &(m_nodes[id]));
				any_changes = true;
			}
		}
		if(min_value_set && !any_changes){
			for(int i = 0; i < 3; i++){
				MeshPoint2d* mpt = containing_triangle->getPoint(i);
				ControlNode2d* node = (ControlNode2d*)mpt->getPtrTag(TagExtended::TAG_ACS_METRIC);
				if(node){
					any_changes |= node->control_data.setMinimum(cdm);
				}else{
					size_t id = m_nodes.add(ControlNode2d(pt, cdm));
					mpt->setPtrTag(TagExtended::TAG_ACS_METRIC, &(m_nodes[id]));
					any_changes = true;
				}
			}
		}
	}

	return any_changes;
}

MeshPoint2d* ControlSpace2dMesh::insertNodeAndRefine(const DPoint2d& point, MeshTriangle2d* containing_triangle)
{
	// create and insert
	MeshPoint2d* mesh_point = new MeshPoint2d(point);
	Metric2dContext mc(m_control_mesh->getControlSpace());
	bool result = MeshGenerator2d::addPointToTriangulation(mc, m_control_mesh, mesh_point, containing_triangle);
	if(result == false){ delete mesh_point; return nullptr; }

	// refine
	int method = m_approximation_method;
	if(method == MeshData::CONTROL_MIXED) method = MeshData::CONTROL_VORONOI;
	const int MAX_STEPS = 100;
	const double QUALITY_THRESHOLD = 0.4;
	for(int k = 0; k < MAX_STEPS; k++){
		// select worst incident triangle
		MeshTriangle2d* worst_triangle = nullptr;
		double worst_quality = 1.0;
		int rank = mesh_point->getRank();
		for(int i = 0; i < rank; i++){	// check all triangles incident to this point
			MeshTriangle2d* triangle = (MeshTriangle2d*)mesh_point->getEdge(i)->getMeshElement(mesh_point);
			double q = triangle->getAlphaQuality(mc, false);
			if(!worst_triangle || q < worst_quality){
				worst_triangle = triangle;
				worst_quality = q;
			}
		}
		if(worst_quality >= QUALITY_THRESHOLD) break;
		// calculate coordinates of the new point
		// no metric necessary here
		DPoint2d new_point = DTriangle2d::outerCircleCenter(
			worst_triangle->getPoint(0)->getCoordinates(), 
			worst_triangle->getPoint(1)->getCoordinates(), 
			worst_triangle->getPoint(2)->getCoordinates());
		if(!m_box.contains(new_point)){ // count as middle of the edge
			bool found = false;
			for(int j = 0; j < 3; j++){
				MeshEdge2d* edge = worst_triangle->getEdge(j);
				if(edge->incidentTo(mesh_point)){
					if(edge->getLength(mc, false) >= 2.0){
						new_point = edge->getPoint(0.5);
						found = true;
						break;
					}
				}
			}
			if(!found) return mesh_point;
		}
		// insert to mesh and set metric
		MeshPoint2d* mpoint = new MeshPoint2d(new_point);
		const ControlDataMatrix2d data_matrix = getMetricAtPoint(new_point, method);
		MeshGenerator2d::addPointToTriangulation(mc, m_control_mesh, mpoint, 
			worst_triangle->findTriangleByNeighbours(mpoint));
		double g_ratio;
		const ControlDataMatrix2d pdm = base_surface->countParameterizationMatrix(new_point, g_ratio);
		size_t id = m_nodes.add(ControlNode2d(new_point, data_matrix, pdm));
		mpoint->setPtrTag(TagExtended::TAG_ACS_METRIC, &(m_nodes[id]));
	}
	// return
	return mesh_point;
}

/// Refines control space to parameterization variance
void ControlSpace2dMesh::adaptToParameterization()
{
	LOG_ASSERT(m_control_mesh);
	assert(m_initialized == 1);
	int pct = m_control_mesh->getPointsCount();
	// Count for existing points
	for(int i = 0; i < pct; i++){
		const MeshPoint2d* mpt = m_control_mesh->getPointAt(i);
		ControlNode2d* node = (ControlNode2d*)mpt->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(!node) continue;
		double g_ratio;
		node->param_data = base_surface->countParameterizationMatrix(mpt->getCoordinates(), g_ratio);
		if(g_ratio < MIN_PARAM_GRATIO){
			LOG4CPLUS_WARN(MeshLog::logger_console, "Improper parameterization " << g_ratio);
		}
	}
	Metric2dContext mc(m_control_mesh->getControlSpace());
	// Check difference
	for(int i = 0; i < pct; i++){
		const MeshPoint2d* mpt = m_control_mesh->getPointAt(i);
		ControlNode2d* node_0 = (ControlNode2d*)mpt->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(!node_0) continue;
		int rank = mpt->getRank();
		for(int j = 0; j < rank; j++){
			const MeshEdge2d* edge = mpt->getEdge(j);
			if(edge->isBorder() || (edge->getPointIndex(mpt) == 1)) continue;
			const MeshPoint2d* other_mpt = edge->getMeshPoint(1);
			ControlNode2d* node_1 = (ControlNode2d*)other_mpt->getPtrTag(TagExtended::TAG_ACS_METRIC);
			if(!node_1) continue;
			double diff = node_0->param_data.countDifferenceRR(node_1->param_data);
			if(diff <= param_threshold_diff) continue;
			const DPoint2d pt_middle = DPoint2d::average(node_0->coord, node_1->coord);
			double g_ratio;
			const ControlDataMatrix2d pcdm = base_surface->countParameterizationMatrix(pt_middle, g_ratio);
			if(g_ratio < MIN_PARAM_GRATIO){
				LOG4CPLUS_WARN(MeshLog::logger_console, "Improper parameterization " << g_ratio);
				continue;
			}

			//MeshView::showDebugMesh("mesh control space - before", m_control_mesh, mpt, other_mpt);

			MeshPoint2d* mesh_point = new MeshPoint2d(pt_middle);

			int method = m_approximation_method;
			if(method == MeshData::CONTROL_MIXED) method = MeshData::CONTROL_VORONOI;
			ControlDataMatrix2d cdm = getMetricAtPoint(pt_middle, method);
			if(MeshGenerator2d::addPointToTriangulation(mc, m_control_mesh, mesh_point, 
					(MeshTriangle2d*)(edge->getMeshElement(0)?edge->getMeshElement(0):edge->getMeshElement(1))))
			{
				//MeshView::showDebugMesh("mesh control space - after", m_control_mesh, mpt, other_mpt);
				size_t id = m_nodes.add(ControlNode2d(pt_middle, cdm));
				mesh_point->setPtrTag(TagExtended::TAG_ACS_METRIC, &(m_nodes[id]));
				m_nodes[id].param_data = pcdm;
//				any_changes = true;
				// count parameterization matrix
				// update control variables
				rank = mpt->getRank();	// j=0; ??
				++pct;
			}else{
				delete mesh_point;
				LOG4CPLUS_WARN(MeshLog::logger_console, "ControlMesh:: Error inserting new node");
				assert(m_control_mesh->isValid());
			}
		}
	}
	m_initialized = 2;
}

/// Smoothen variance of metric within the control space
bool ControlSpace2dMesh::smoothen()
{
	assert(m_initialized == 2); // with parameterization matrix
	// reset gradation ratio in nodes
	for(int i = 0; i < m_control_mesh->getPointsCount(); i++){
		ControlNode2d *node = (ControlNode2d*)m_control_mesh->getPointAt(i)->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(node) node->max_gradation_ratio = 1.0;
	}
	int count = 0;
	// TODO
	int pct = m_control_mesh->getPointsCount();
	DataVector<MeshPoint2d*> active_points(pct);
	for(int i = 0; i < pct; i++){
		MeshPoint2d* mpt = m_control_mesh->getPointAt(i);
		mpt->setIntTag(TagExtended::TAG_ACS);
		active_points.add(mpt);
	}
	while(active_points.countInt() > 0){
		assert(isValid());
		MeshPoint2d* mpt = active_points.removeAt(0);
		assert(mpt->nonZeroIntTag(TagExtended::TAG_ACS));
		mpt->removeTag(TagExtended::TAG_ACS);
		ControlNode2d *node0 = (ControlNode2d*)mpt->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(node0 == nullptr) continue;
		int rank = mpt->getRank();
		for(int i = 0; i < rank; i++){
			MeshEdge2d* edge = mpt->getEdge(i);
			if(edge->getPointIndex(mpt) != 0) continue;
			// check
			MeshPoint2d* mpt_1 = edge->getMeshPoint(1);
			ControlNode2d *node1 = (ControlNode2d*)mpt_1->getPtrTag(TagExtended::TAG_ACS_METRIC);
			if(node1 == nullptr) continue;
//-------------------------
			int result = smoothenMetricForNodes(node0, node1, 
				(node1->coord - node0->coord).normalized() * edge->getLength(base_surface));
			if((result & 1) != 0){
				++count;
				if(mpt->zeroIntTag(TagExtended::TAG_ACS)){
					mpt->setIntTag(TagExtended::TAG_ACS);
					active_points.add(mpt);
				}
			}
			if((result & 2) != 0){
				++count;
				if(mpt_1->zeroIntTag(TagExtended::TAG_ACS)){
					mpt_1->setIntTag(TagExtended::TAG_ACS);
					active_points.add(mpt_1);
				}
			}
		}
	}
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "ControlSpace2dMesh smoothing - " << count << " modifications.");
	return count > 0;
}

ControlDataMatrix2d ControlSpace2dMesh::getMetricAndParameterizationAtPoint(
	const DPoint2d& pt, ControlDataMatrix2d& p_cdm) const
{
	assert(m_initialized>1);
	return getMetricAtPoint(pt, -1, &p_cdm);
}

bool ControlSpace2dMesh::isValid() const
{
	MeshContainer2d* mesh = m_control_mesh;
	if(!mesh) mesh = m_control_points;
	if(!mesh) return true;
	int pct = mesh->getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint2d* mpt = mesh->getPointAt(i);
		ControlNode2d *node = (ControlNode2d*)mpt->getPtrTag(TagExtended::TAG_ACS_METRIC);
		if(node && node->control_data.det() < 0){
			return false;
		}
	}
	return true;
}

void ControlSpace2dMesh::setMaxMetric(double ratio)
{
	if(m_control_points)
		ControlSpace2dAdaptive::setMaxMetric();
	else{
		DataVector<DPoint2d> points;
		countRegularGridPoints(param_control_nxy, points);

		double p = m_box.getDiameter();
		if(ratio <= 0) ratio = std::min(param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
		double max_len = p * ratio;
		ControlDataMatrix2d max_data(max_len, 0.0, max_len);

		for(size_t i = 0; i < points.countInt(); i++)
			addControlNode(ControlNode2d(points.get(i), max_data));
	}

	if(!m_control_mesh) interpolate();
	m_initialized = 1;
}

