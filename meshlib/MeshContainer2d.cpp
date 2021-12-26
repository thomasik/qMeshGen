/////////////////////////////////////////////////////////////////////////////
// MeshContainer2d.cpp
// Klasa odpowiedzialna za przechowywanie punktów i elementów siatki
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////

#include "MeshContainer2d.h"
#include "MeshEdge2d.h"
#include "MeshEdge3d.h"
#include "DataContainer.h"
#include "Curve2dParametric.h"
#include "MeshArea.h"
#include "SurfaceParametric.h"
#include "ControlSpace2dIdentity.h"
#include "ControlSpace2dMatrixUniform.h"
#include "ControlSpace2dMesh.h"
#include "ControlSpace2dQuadTree.h"
#include "ControlSpace3d.h"
#include "QuadTree.h"
#include "SurfacePlane.h"
#include "EPSFile.h"
#include "MeshViewSet.h"
#include "MeshBoundaryCondition.h"
#include "MeshGenerator2d.h"
#include "IteratorEdge2d.h"
#include "DataHashTable.h"
#include "DTriangle.h"
#include "DQuad.h"

/////////////////////////////////////////////////////////////////////////////
// Standardowy konstruktor (part_size) okreœla mno¿nik rozmiaru tablic, które
//	zwiêkszane s¹ w razie potrzeby krokowo o t¹ w³aœnie wartoœæ
MeshContainer2d::MeshContainer2d(int part_size) :
	m_part_size(part_size), m_last_triangle(nullptr), m_quad_tree(nullptr),
	m_discretization_state(0), m_constraining(MeshData::CONSTRAIN_NONE)
{
	m_points = new DataContainer<MeshPoint2d>(part_size);
	m_elements = new DataContainer<MeshElement>(part_size);
}

/////////////////////////////////////////////////////////////////////////////
// Destruktor
MeshContainer2d::~MeshContainer2d()
{
	deleteAll();
	delete m_elements;
	delete m_points;
	if(m_quad_tree) delete m_quad_tree;
}

void MeshContainer2d::removeAllTags(TagExtended::TagType tag_type)
{
	for(int i = 0; i < m_points->countInt(); i++) // clear all tags for points, edges and elements
		getPointAt(i)->removeTag(tag_type);
	for(int i = 0; i < m_elements->countInt(); i++) 
		getElementAt(i)->removeTag(tag_type);
	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge())
		it.getEdge()->removeTag(tag_type);
}

/////////////////////////////////////////////////////////////////////////////
// PrzyPIsuje ka¿demu punktowi numer wynikaj¹cy z jego indeksu w tablicy
void MeshContainer2d::renumeratePoints()
{
	int count = m_points->countInt();
	for(int i = 0; i < count; i++){
		m_points->getDataAt(i)->setIntTag(TagExtended::TAG_ID, i);
	}
}

/////////////////////////////////////////////////////////////////////////////
// Nadaje numery punktom (oraz odpowiednio je przemieszcza) zgodnie z tablicami
//	incydencji punktów
void MeshContainer2d::renumeratePointsByNeighbours(bool and_reverse, bool distance_sort)
{
	int count = m_points->countInt();
	if(count < 1) return;
	int i, k = 1, start = 0;

	for(i = 0; i < count; i++){
		MeshPoint2d *pt2, *pt1 = m_points->getDataAt(i);
		pt1->setIntTag(TagExtended::TAG_ID, i);
		int rank = pt1->getRank();
		if(k <= i){	// Wa¿ne dla niespójnych siatek
			k = i+1;
			if(distance_sort) sortByDistance(start+1, i, start);
			start = i+1;
		}
		// Dla wszystkich punktów z s¹siedztwa
		if(rank > 1){
			for(int n = 0; n < rank; n++){
				MeshEdge2d* edge = pt1->getEdge(n);
/*
				// Punkty wewnêtrzne krawêdzi
				int inner_count = edge->getInnerPointsCount();
				for(int m = 0; m < inner_count; m++){
					pt2 = edge->getInnerMeshPoint(m);
					int id2 = pt2->getIndex();
					if(id2 > k){
						// Zamieñ punkty - tzn. przydziel kolejny numer s¹siadowi
						m_points->switchDataItems(id2, k++);
					}else if(id2 == k){
						k++;
					}
				}
*/
				pt2 = edge->getOtherPoint(pt1);
				int id2 = pt2->getIndex();
				if(id2 > k){
					// Zamieñ punkty - tzn. przydziel kolejny numer s¹siadowi
					m_points->switchDataItems(id2, k++);
				}else if(id2 == k){
					k++;
				}
			}
		}
	}
	if(distance_sort && start < count-1) sortByDistance(start+1, count-1, start);

	if(!and_reverse) return;
	// Druga faza - wstecz
	start = count - 1;
	k = count - 2;
	for(i = count - 1; i >= 0; i--){
		MeshPoint2d *pt2, *pt1 = m_points->getDataAt(i);
		pt1->setIntTag(TagExtended::TAG_ID, i);
		int rank = pt1->getRank();
		if(k >= i){	// Wa¿ne dla niespójnych siatek
			k = i - 1;
			if(distance_sort) sortByDistance(i, start-1, start);
			start = i - 1;
		}
		// Dla wszystkich punktów z s¹siedztwa
		if(rank > 1){
			for(int n = 0; n < rank; n++){
				MeshEdge2d* edge = pt1->getEdge(n);
/*
				// Punkty wewnêtrzne krawêdzi
				int inner_count = edge->getInnerPointsCount();
				for(int m = 0; m < inner_count; m++){
					pt2 = edge->getInnerMeshPoint(m);
					int id2 = pt2->getIndex();
					if(id2 < k){
						// Zamieñ punkty - tzn. przydziel s¹siedni numer s¹siadowi
						m_points->switchDataItems(id2, k--);
					}else if(id2 == k){
						k--;
					}
				}
*/
				pt2 = edge->getOtherPoint(pt1);
				int id2 = pt2->getIndex();
				if(id2 < k){
					// Zamieñ punkty - tzn. przydziel s¹siedni numer s¹siadowi
					m_points->switchDataItems(id2, k--);
				}else if(id2 == k){
					k--;
				}
			}
		}
	}
	if(distance_sort && start > 0) sortByDistance(0, start-1, start);
}


/////////////////////////////////////////////////////////////////////////////
// Usuwa dane dotycz¹ce przechowywanej geometrii obszaru
void MeshContainer2d::deleteAll()
{
	m_elements->deleteAll();
	m_points->deleteAll();
}

//////////////////////////////////////////////////////////////////////
// Zwraca wspó³rzêdne prostok¹ta obejmuj¹cego wszystkie punkty 
DRect MeshContainer2d::getBoundingRect() const
{
	DRect rect;

	int count = m_points->countInt();
	for(int i = 0; i < count; i++){
		MeshPoint2d* point = m_points->getDataAt(i);
		if(point->getRank() == 0) // no incident edge
			rect.addPoint(point->getCoordinates());
	}

	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge())
		it.getEdge()->addToBoundingRect(rect);

	return rect;
}

/*
//////////////////////////////////////////////////////////////////////
// Insert given number of inner nodes for all edges within this mesh
void MeshContainer2d::setEdgeInnerPoints(Metric2dContext& mc, int count, bool use_shapes)
{
	// W trakcie dodawania punktów wewnêtrznych ta wartoœæ bêdzie
	//	zwiêkszana, ale przegl¹danie nowych punktów by³oby bezu¿yteczne
	int pct = m_points->countInt();	
									
	int id = pct;
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = m_points->getDataAt(i);
		int ect = point->getRank();
		if(ect < 2) continue;
		for(int k = 0; k < ect; k++){	// Dla wszystkich krawêdzi incydentnych do kolejnego punktu
			MeshEdge2d* edge = point->getEdge(k);
			if(edge->getPointIndex(point) == 0){	// Ka¿da krawêdŸ tylko raz
				MeshPoint2d **points, **old_points;
				if(count > 0){
					points = new MeshPoint2d*[count];	
				}else{
					points = nullptr;
				}
				// Tworzenie tablicy nowych punktów
				int n;
				for(n = 0; n < count; n++){
					double t;
					DPoint2d pt = edge->getAspectPoint(mc, ((n+1) / (double)(count+1)), t, use_shapes);
					m_points->addDataItem(points[n] = new MeshPoint2d(++id, pt.x, pt.y));
					points[n]->addEdge(edge);
					points[n]->setAsInnerNode();
				}
				// Zamiana tablic
				int old_ct = edge->addInnerPoints(count, points, &old_points);
				// Usuwanie starych punktów (o ile istniej¹)
				for(n = 0; n < old_ct; n++){
					delete removeMeshPoint(old_points[n]->getIndex());
				}
				if(old_ct > 0) delete[] old_points;
				if(old_ct > count) pct -= (old_ct - count);
			}
		}
	}

	m_inner_edges_ct = count;

	renumeratePoints();
}
*/
struct SortData{ 
	double dist; 
	MeshPoint2d* point; 
};

int compareSortData( const void *arg1, const void *arg2 )
{
	return ((*(const SortData*)arg1).dist > (*(const SortData*)arg2).dist) ? 1 : -1;
}

/////////////////////////////////////////////////////////////////////////////
//	Porz¹dkowanie punktów wed³ug odleg³oœci od PIerwszego
void MeshContainer2d::sortByDistance(int start, int end, int target)
{
	int sort_size = end - start + 1;
	if(sort_size < 2) return;

	int i, direction;
	if(target > end){
		direction = -1;	// malej¹ce odleg³oœci w kierunku target
	}else{
		direction = 1;	// rosn¹ce odleg³oœci od target
	}
	SortData *points = new SortData[sort_size];

	DPoint2d target_point = m_points->getDataAt(target)->getCoordinates();
	int index = 0;
	for(i = start; i <= end; i++){
		points[index].point = m_points->getDataAt(i);
		points[index].dist = direction * points[index].point->getCoordinates().distance2(target_point);
		index++;
	}

	// Sortuj rosn¹co tablicê
	qsort(points, (size_t)sort_size, sizeof(SortData), compareSortData);

	//Przemieœæ punkty
	for(i = 0; i < sort_size; i++){
		MeshPoint2d* pt = points[i].point;
		int k = i + start;
		int id = pt->getIndex();
		if(k != id){
			// Zamieñ
			m_points->switchDataItems(k, id);
		}
		pt->setIntTag(TagExtended::TAG_ID, k);
	}

	delete[] points;
}

void MeshContainer2d::setControlSpace(CS2dPtr space)
{
	m_control = space;
}

int MeshContainer2d::getElementsCount(int edge_count) const
{
	int ect = getElementsCount();
	int count = 0;
	for(int i = 0; i < ect; i++) 
		if(getElementAt(i)->getEdgeCount() == edge_count) ++count;
	return count;
}

bool MeshContainer2d::anyEdgeWithShape() const
{
	int pct = getPointsCount();
	for(int i = 0; i < pct; i++){
		const MeshPoint2d* point = getPointAt(i);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			if(point->getEdge(j)->getType() == EDGE_CURVE) return true;
		}
	}
	return false;
}

bool MeshContainer2d::createControlSpaceFromTriangles(double factor)
{
	if(!m_surface) return false;

	START_CLOCK("MC2d::createControlSpaceFromTriangles");

	DRect rect = getBoundingRect();
	rect.inflate(ControlSpace2d::param_inflate_box_factor);

	switch(ControlSpace2d::param_control_type){
	case MeshData::CONTROL_UNIFORM:
		m_control = std::make_shared<ControlSpace2dMatrixUniform>(m_surface, rect, 
					ControlSpace2dMatrixUniform::param_uniform_nx, 
					ControlSpace2dMatrixUniform::param_uniform_nx);
		break;
	case MeshData::CONTROL_MESH:
		m_control = std::make_shared<ControlSpace2dMesh>(m_surface, rect,
					ControlSpace2dMesh::param_interpolation_method);
		break;
	case MeshData::CONTROL_QUADTREE:
		m_control = std::make_shared<ControlSpace2dQuadTree>(m_surface, rect,
					ControlSpace2dAdaptive::param_control_nxy);
		break;
	default:
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Unknown CONTROL_TYPE!");
		break;
	}

	int tct = getElementsCount();
	auto space = m_control->getAsAdaptive();
	if (space) {
		if (false) {
			for (int i = 0; i < tct; i++) {
				MeshTriangle2d* triangle = (MeshTriangle2d*)getElementAt(i);
				if (triangle->getType() != ELEMENT_MESH_TRIANGLE) continue;
				ControlDataMatrix2d cdm_simplex = ControlDataMatrix2d::countMetric(
					triangle->getPoint(0)->getCoordinates(),
					triangle->getPoint(1)->getCoordinates(),
					triangle->getPoint(2)->getCoordinates());
				space->addControlNode(ControlNode2d(triangle->getMiddlePoint(), cdm_simplex * factor));
			}
		}
		else {
			int pct = getPointsCount();
			ControlDataMatrix2d *cdms = new ControlDataMatrix2d[pct];
			int *cdms_ct = new int[pct];
			for (int i = 0; i < pct; i++) cdms_ct[i] = 0;
			for (int i = 0; i < tct; i++) {
				MeshTriangle2d* triangle = (MeshTriangle2d*)getElementAt(i);
				if (triangle->getType() != ELEMENT_MESH_TRIANGLE) continue;
				ControlDataMatrix2d cdm = ControlDataMatrix2d::countMetric(
					triangle->getPoint(0)->getCoordinates(),
					triangle->getPoint(1)->getCoordinates(),
					triangle->getPoint(2)->getCoordinates());
				for (int j = 0; j < 3; j++) {
					int id = triangle->getPoint(j)->getIndex();
					cdms[id] += cdm;
					cdms_ct[id]++;
				}
			}
			for (int i = 0; i < pct; i++)
				if (cdms_ct[i] > 0)
					space->addControlNode(ControlNode2d(getPointAt(i)->getCoordinates(),
						cdms[i] * (factor / cdms_ct[i])));
			delete[] cdms_ct;
			delete[] cdms;
		}
		space->interpolate();
		space->adaptToParameterization();

		if (ControlSpace2dAdaptive::param_gradation_ratio > 0.0) space->smoothen();
	}
	STOP_CLOCK("MC2d::createControlSpaceFromTriangles");
//	space->logDescription();
//	space->storeEPS("control-space");
	return true;
}

bool MeshContainer2d::createControlSpace(CS2dPtr user_space,
	CS3dPtr user_space_3d_1, CS3dPtr user_space_3d_2)
{
	if(!m_surface) return false;
	m_control.reset();

	START_CLOCK("MC2d::createControlSpace");

	DRect rect = getBoundingRect();
	rect.inflate(ControlSpace2d::param_inflate_box_factor);

	bool will_use_surface_curvature = 
		ControlSpace2dAdaptive::param_use_surface_curvature &&
		m_surface->getType() != SURFACE_PLANE;
	bool will_use_contour_curvature = 
		ControlSpace2dAdaptive::param_use_contour_curvature &&
		(m_surface->getType() != SURFACE_PLANE || anyEdgeWithShape());
	bool will_use_adaptive_user_space = user_space && 
		user_space->isAdaptive();

	if(will_use_surface_curvature ||
		(!will_use_adaptive_user_space && will_use_contour_curvature) ||
		!user_space)
	{
		m_control = MeshGenerator2d::createNewControlSpace(m_surface, rect);
	}

#ifdef _DEBUG
//#define STORE_CONTROL_EPS
#endif

	auto space = m_control->getAsAdaptive();

	if(will_use_surface_curvature){
		assert(space);
		// Curvature of surface - base control data
		space->setSurfaceCurvatureControlData();
	#ifdef STORE_CONTROL_EPS
		space->storeEPS("control-with-surface");
	#endif

	}

	if(will_use_adaptive_user_space){	// if exists and is adaptive
		if(space)
			space->applyAsMinimum(user_space);
		else if(!m_control){
			m_control = user_space;
		}
	}

	//-----------------------------------
	// Curvature of boundaries - additional control data

	if(will_use_contour_curvature){
		assert(space);
		if((space->initializationState() == 0) && !space->interpolate())
			space->setMaxMetric();
		space->setContourCurvatureControlData(this);
	#ifdef STORE_CONTROL_EPS
		space->storeEPS("control-with-contours");
	#endif

	}

	if(user_space && !user_space->isAdaptive()){	// if exists and is not adaptive
		if(space)
			space->applyAsMinimum(user_space);
		else if(!m_control){
			LOG4CPLUS_WARN(MeshLog::logger_console, "Using non-adaptive user control space!");
			m_control = user_space;
		}
	}

	assert(m_control);
	if(!m_control) m_control = std::make_shared<ControlSpace2dQuadTree>(m_surface, rect, 10);

	if(!space->isAdaptive()) return true;

	if((space->initializationState() == 0) && !space->interpolate())
		space->setMaxMetric();
	space->adaptToParameterization();

#ifdef STORE_CONTROL_EPS
	space->storeEPS("control-with-adapt-parameterization");
#endif

	if(ControlSpace2dAdaptive::param_gradation_ratio > 0.0)
		space->smoothen();

	if(user_space_3d_1)
		space->applyAsMinimum(user_space_3d_1);
	if(user_space_3d_2)
		space->applyAsMinimum(user_space_3d_2);

	STOP_CLOCK("MC2d::createControlSpace");

/*
	if(true){
		DPoint2d left_pt(rect.left, 0.5*(rect.top + rect.bottom));
		DPoint2d right_pt(rect.right, left_pt.y);
		const DVector2d dvx(1.0, 0.0);
		ofstream ftest("ctest.txt");
		for(double t = 0.0; t <= 1.0; t += 0.01){
			DPoint2d pt = left_pt * (1.0-t) + right_pt * t;
			ControlDataMatrix2d cdm = space->getMetricAtPoint(pt);
			double len = (cdm * dvx).length();
			ftest << pt.x << "\t" << len << endl;
		}
	}
*/

//	space->logDescription();
	return true;
}

StatData MeshContainer2d::getElementsIDSpan() const
{
	StatData stat;
	int count = m_elements->countInt();
	if(count < 1) return stat;

	int q = m_elements->getDataAt(0)->getIdSpan();
	stat.minimum = stat.maximum = stat.average = q;

	for(int i = 1; i < count; i++){
		q = m_elements->getDataAt(i)->getIdSpan();
		stat.average += q;
		if(q > stat.maximum) stat.maximum = q;
		if(q < stat.minimum) stat.minimum = q;
	}
	stat.average /= count;

	return stat;

}

bool MeshContainer2d::isValid()
{
	// check for proper element orientation;
	int ect = getElementsCount();
	int pct = getPointsCount();
	for(int i = 0; i < ect; i++){
		MeshElement* element = getElementAt(i);
		if(!element){
			LOG4CPLUS_INFO(MeshLog::logger_console, "invalid> empty element?");
			assert(false);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> empty element?");
			return false;
		}
		if(element->isInverted()){
			LOG4CPLUS_INFO(MeshLog::logger_console, "invalid> inverted element?");
//			MeshView::showDebugMesh("inverted element", this, element);
			assert(false);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> inverted element?");
			return false;
		}
		int edge_ct = element->getEdgeCount();
		for(int j = 0; j < edge_ct; j++){
			MeshPoint2d* point = element->getPoint(j);
			if(!point || point->getIndex() < 0 || point->getIndex() >= pct){
				LOG4CPLUS_INFO(MeshLog::logger_console, "invalid> erroneous point?");
				assert(false);
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> erroneous point?");
				return false;
			}
			MeshEdge2d* edge = element->getEdge(j);
			if(!edge) return false;
			if(!edge->isBorder() && (edge->getMeshElement(0) == nullptr || edge->getMeshElement(1) == nullptr)){
				LOG4CPLUS_INFO(MeshLog::logger_console, "invalid> inconsistent edge?");
				assert(false);
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> inconsistent edge?");
				return false;
			}
		}
	}
	for(int i = 0; i < pct; i++){
		assert(getPointAt(i)->getRank() > 0);
	}

/*
	// check for inner-edges adjacency
	int pct = getPointsCount();
	for(i = 0; i < pct; i++){
		MeshPoint2d* point = getPointAt(i);
		if(!point->isBorder()){
			int rank = point->getRank();
			if(rank > 1){
				for(int j = 0; j < rank; j++){
					MeshEdge2d* edge = point->getEdge(j);
					if(edge->getMeshElement(0) == nullptr ||
						edge->getMeshElement(1) == nullptr)
						return false;
				}
			}
		}
	}
*/
	return true;
}

void MeshContainer2d::clearSearchTree()
{
	if(m_quad_tree) { 
		delete m_quad_tree; 
		m_quad_tree = nullptr; 
	} 
}

int MeshContainer2d::getMaxSearchTreeLevel() const
{
	return m_quad_tree ? m_quad_tree->getMaxLevel() : 0;
}

void MeshContainer2d::setSearchTree(QuadTree *tree)
{
	if(m_quad_tree) delete m_quad_tree;
	m_quad_tree = tree;
}

int MeshContainer2d::addMeshTriangle(MeshTriangle2d *triangle, bool no_heap)
{
	m_last_triangle = triangle;
	if(m_quad_tree) m_quad_tree->insertTriangleLink(triangle);
	return m_elements->addDataItem(triangle, no_heap);
}

MeshTriangle2d* MeshContainer2d::removeMeshTriangle(MeshTriangle2d* triangle)
{ 
	if(m_quad_tree) m_quad_tree->removeTriangleLink(triangle);
	return (MeshTriangle2d*)m_elements->removeDataItem(triangle->getIndex());
}

MeshTriangle2d* MeshContainer2d::getNearTriangle(const DPoint2d& pt)
{ 
	if(m_quad_tree)
		return m_quad_tree->getNearestTriangle(pt, m_last_triangle);
	else return m_last_triangle; 
}

bool MeshContainer2d::storeTxt(const string& fname, int id) const
{
	int pct = getPointsCount();
	int ect = getElementsCount();
	if(pct < 1 || ect < 1) return false;

	ostringstream fname1,fname2,fname3;
	fname1 << fname << '-' << id << "-p.txt";
	fname2 << fname << '-' << id << "-i.txt";
	fname3 << fname << '-' << id << "-e.txt";

	ofstream fpoints(fname1.str().c_str());
	ofstream felements(fname2.str().c_str());
	ofstream fedges(fname3.str().c_str());

	fpoints << pct << endl;
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = getPointAt(i);
		DPoint2d pt = point->getCoordinates();
		fpoints << i << "\t" << pt.x << "\t" << pt.y << endl;
	}

	felements << ect << endl;
	for(int i = 0; i < ect; i++){
		MeshElement* element = getElementAt(i);
		felements << i;
		int edges_ct = element->getEdgeCount();
		for(int j = 0; j < edges_ct; j++)
			felements << '\t' << element->getPoint(j)->getIndex();
		felements << endl;
	}

	DataVector<MeshEdge2d*> all_edges(4*pct);
	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge()){
		all_edges.add(it.getEdge());
	}
	fedges << all_edges.countInt() << endl;
	for(size_t i = 0; i < all_edges.countInt(); i++){
		fedges << i;
		fedges << '\t' << all_edges[i]->getMeshPoint(0)->getIndex();
		fedges << '\t' << all_edges[i]->getMeshPoint(1)->getIndex();
		fedges << endl;
	}

	return true;
}

bool MeshContainer2d::storeSurf(const string& fname, bool proper_orientation, int id) const
{
	int pct = getPointsCount();
	int ect = getElementsCount();
	if(pct < 1 || ect < 1 || !m_surface) return false;

	ostringstream fname1,fname2;
	fname1 << fname << '-' << id << "-p.txt";
	fname2 << fname << '-' << id << "-i.txt";

	ofstream fpoints(fname1.str().c_str());
	ofstream felements(fname2.str().c_str());

	fpoints << pct << endl;
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = getPointAt(i);
		const DPoint3d pt = m_surface->getPoint(point->getCoordinates());
		fpoints << i << "\t" << pt.x << "\t" << pt.y  << "\t" << pt.z << endl;
	}

	felements << ect << endl;
	for(int i = 0; i < ect; i++){
		MeshElement* element = getElementAt(i);
		int edges_ct = element->getEdgeCount();
		felements << i << '\t' << edges_ct;
		for(int j = 0; j < edges_ct; j++)
			if(proper_orientation)
				felements << '\t' << element->getPoint(j)->getIndex();
			else
				felements << '\t' << element->getPoint(edges_ct-j-1)->getIndex();
		felements << '\t' << element->getAreaID() << endl;
	}

	return true;
}

MeshData::StatData MeshContainer2d::getAlphaQuality(bool on_surface) const
{
	MeshData::StatData data(1.0, 0.0, 1.0);
	int ect = getElementsCount();
	double n_power = 1.0 / (double)ect;
	double x = 1.0;

	for(int i = 0; i < ect; i++){
		MeshElement* element = getElementAt(i);
		double q;
		if(on_surface){
			const DPoint3d& pt0 = m_surface->getPoint(element->getPoint(0)->getCoordinates());
			const DPoint3d& pt1 = m_surface->getPoint(element->getPoint(1)->getCoordinates());
			const DPoint3d& pt2 = m_surface->getPoint(element->getPoint(2)->getCoordinates());
			if(element->getEdgeCount() == 3){
				q = DTriangle3d::alphaQuality(pt0, pt1, pt2);
			}else{
				const DPoint3d& pt3 = m_surface->getPoint(element->getPoint(3)->getCoordinates());
				q = DQuad3d::alphaQuality(pt0, pt1, pt2, pt3);
			}
		}else{
			const DPoint2d& pt0 = element->getPoint(0)->getCoordinates();
			const DPoint2d& pt1 = element->getPoint(1)->getCoordinates();
			const DPoint2d& pt2 = element->getPoint(2)->getCoordinates();
			if(element->getEdgeCount() == 3){
				q = DTriangle2d::alphaQuality(pt0, pt1, pt2);
			}else{
				const DPoint2d& pt3 = element->getPoint(3)->getCoordinates();
				q = DQuad2d::alphaQuality(pt0, pt1, pt2, pt3);
			}
		}
		if(q < data.minimum) data.minimum = q;
		if(q > data.maximum) data.maximum = q;
		if(data.minimum > 0.0){
			x *= q;
			if(x < 0.1){
				data.average *= pow(x, n_power);
				x = 1.0;
			}
		}
	}

	if(data.minimum > 0.0)
		data.average *= pow(x, n_power);
	else
		data.average = 0.0;

	return data;
}

MeshData::StatData MeshContainer2d::getAngleQuality(bool on_surface) const
{
	MeshData::StatData data(PI, 0.0, 0.0);
	int ect = getElementsCount();
	int count = 0;

	for(int i = 0; i < ect; i++){
		MeshElement* element = getElementAt(i);
		int edge_ct = element->getEdgeCount();
		if(on_surface){
			DPoint3d pt[4];
			for(int j = 0; j < edge_ct; j++)
				pt[j] = m_surface->getPoint(element->getPoint(j)->getCoordinates());
			for(int j = 0; j < edge_ct; j++){
				double angle = (pt[(j+1)%edge_ct] - pt[j]).getAngle(pt[(j+2)%edge_ct] - pt[j]);
				if(angle < data.minimum) data.minimum = angle;
				if(angle > data.maximum) data.maximum = angle;
				data.average += angle;
				++count;
			}
		}else{
			DPoint2d pt[4];
			for(int j = 0; j < edge_ct; j++)
				pt[j] = element->getPoint(j)->getCoordinates();
			for(int j = 0; j < edge_ct; j++){
				double angle = (pt[(j+1)%edge_ct] - pt[j]).getAngle(pt[(j+2)%edge_ct] - pt[j]);
				if(angle < data.minimum) data.minimum = angle;
				if(angle > data.maximum) data.maximum = angle;
				data.average += angle;
				++count;
			}
		}
	}

	data.average /= count;
	return data;
}

bool MeshContainer2d::statMetricDifference(Metric2dContext& mc, DataStatistics& stats) const
{
	int ect = getElementsCount();
	for(int i = 0; i < ect; i++){
		const MeshTriangle2d* triangle = (MeshTriangle2d*)getElementAt(i);
		if(triangle->getType() != ELEMENT_MESH_TRIANGLE) continue; // triangles only
			//  - 2D metric from simplex (three points)
		double diff = triangle->countMetricDiff(mc);
		stats.add(diff);
	}
	return true;
}

bool MeshContainer2d::statMeanRatio(Metric2dContext& mc, DataStatistics& stats, bool ext_metric) const
{
	int ect = getElementsCount();
	for(int i = 0; i < ect; i++){
		const MeshTriangle2d* triangle = (MeshTriangle2d*)getElementAt(i);
		if(triangle->getType() != ELEMENT_MESH_TRIANGLE) continue; // triangles only
			//  - 2D metric from simplex (three points)
		double diff = triangle->getMeanRatio(mc, ext_metric);
		stats.add(diff);
	}
	return true;
}

bool MeshContainer2d::statMetricEdgeLength(Metric2dContext& mc, DataStatistics& stats, bool ext_metric) const
{
	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge())
		stats.add(it.getEdge()->getLengthQuality(mc, ext_metric));
	return true;
}

bool MeshContainer2d::statMetricAlphaQuality(Metric2dContext& mc, DataStatistics& stats, bool /*ext_metric*/) const
{
	int ect = getElementsCount();
	for(int i = 0; i < ect; i++)
		stats.add(getElementAt(i)->getAlphaQuality(mc, true));
	return true;
}

bool MeshContainer2d::statRealAngles(DataStatistics& stats) const
{
	int ect = getElementsCount();
	for(int i = 0; i < ect; i++){
		const MeshElement* element = getElementAt(i);
		DataVector<double> angles = element->getRealAngles(m_surface);
		for(size_t j = 0; j < angles.countInt(); j++)
			stats.add(angles[j]);
	}
	return true;
}

void MeshContainer2d::storeEPS(const char* name, int id) const
{
	ostringstream fname;
	fname << name << '-' << id << ".eps";

	DRect box = getBoundingRect();
	EPSFile eps(fname.str(), box.x0, box.x1, box.y0, box.y1);

	int pct = getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = getPointAt(i);
		const DPoint2d& pt1 = point->getCoordinates();
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

MeshViewSet* MeshContainer2d::getViewSet(MeshViewSet* view_set, bool with_surface, bool proper_orientation, bool with_elements) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;

	int ect = getElementsCount();
	int pct = getPointsCount();
	if(ect + pct == 0) return view_set;

	if(view_set){
		view_set->prepareFreePlace(pct, with_elements?2*ect:0, ect, 0);
	}else{
		view_set = new MeshViewSet(pct, with_elements?2*ect:0, ect, 0);
	}
	view_set->setMesh(this);

	// elements
	if(with_elements)
		for(int i = 0; i < ect; i++)
			view_set->addElement(getElementAt(i), with_surface?m_surface:nullptr, proper_orientation);

	// others
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = getPointAt(i);
		view_set->addPoint(point, with_surface?m_surface:nullptr);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0){
				view_set->addEdge(edge, with_surface?m_surface:nullptr);
			}
		}
	}

	return view_set;
}

MeshViewSet* MeshContainer2d::getDebugViewSet(const MeshElement* el1, const MeshElement* el2, double radius) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;
	if(!el1) return nullptr;

	ControlSpace2dIdentity cs;
	Metric2dContext mc_identity(&cs);
	double max_len = 0.0;
	int edge_count = el1->getEdgeCount();
	for(int i = 0; i < edge_count; i++)
		max_len = std::max(max_len, el1->getEdge(i)->getLength(mc_identity));
	if(el2){
		edge_count = el2->getEdgeCount();
		for(int i = 0; i < edge_count; i++)
			max_len = std::max(max_len, el2->getEdge(i)->getLength(mc_identity));
	}
	const DPoint2d middle = el1->getMiddlePoint();
	double r2 = sqr(radius*max_len);

	int ect = getElementsCount();
	int pct = getPointsCount();

	ofstream of("debug_view.log");

	MeshViewSet *view_set = new MeshViewSet(pct, 2*ect, ect, 0);
	view_set->setMesh(this);

	// elements
	for(int i = 0; i < ect; i++){
		MeshElement* el = getElementAt(i);
		if(middle.distance2(el->getMiddlePoint()) > r2) continue;
		if(el == el1)
			view_set->addElementWithEdges(el, nullptr, true, 1);
		else if(el == el2)
			view_set->addElementWithEdges(el, nullptr, true, 2);
		else
			view_set->addElementWithEdges(el, nullptr, true, 0);
		// log
		int edct = el->getEdgeCount();
		of << "Element index=" << el->getIndex() << " pts=" << edct << " [";
		for(int j = 0; j < edct; j++) of << el->getPoint(j)->getIntTag(TagExtended::TAG_ID) << " ";
		of << "]" << endl;
		of << "\tinverted=" << (el->isInverted()?"true":"false") << endl;
		of << "\tarea-012=" << DTriangle2d::area(
			el->getPoint(0)->getCoordinates(),
			el->getPoint(1)->getCoordinates(),
			el->getPoint(2)->getCoordinates()) << endl;
		if(edct == 4){
			of << "\tarea-123=" << DTriangle2d::area(
				el->getPoint(1)->getCoordinates(),
				el->getPoint(2)->getCoordinates(),
				el->getPoint(3)->getCoordinates()) << endl;
			of << "\tarea-230=" << DTriangle2d::area(
				el->getPoint(2)->getCoordinates(),
				el->getPoint(3)->getCoordinates(),
				el->getPoint(0)->getCoordinates()) << endl;
			of << "\tarea-301=" << DTriangle2d::area(
				el->getPoint(3)->getCoordinates(),
				el->getPoint(0)->getCoordinates(),
				el->getPoint(1)->getCoordinates()) << endl;
		}
	}

	// others
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = getPointAt(i);
		if(middle.distance2(point->getCoordinates()) > r2) continue;
		view_set->addPoint(point, nullptr);
	}

	return view_set;
}

MeshViewSet* MeshContainer2d::getDebugViewSet(const MeshPoint2d* pt1, const MeshPoint2d* pt2, double radius) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;
	if(!pt1) return nullptr;

	ControlSpace2dIdentity cs;
	Metric2dContext mc_identity(&cs);
	const DPoint2d middle = pt1->getCoordinates();
	double r2 = 0.0;
	if(pt2) r2 = middle.distance2(pt2->getCoordinates());
	else{
		int rank = pt1->getRank();
		for(int i = 0; i < rank; i++)
			r2 = std::max(r2, sqr(pt1->getEdge(i)->getLength(mc_identity)));
	}
	r2 *= sqr(radius);

	int ect = getElementsCount();
	int pct = getPointsCount();

	ofstream of("debug_view.log");

	MeshViewSet *view_set = new MeshViewSet(pct, 2*ect, ect, 0);
	view_set->setMesh(this);

	// elements
	for(int i = 0; i < ect; i++){
		MeshElement* el = getElementAt(i);
		if(el->getType() == ELEMENT_MESH_AREA) continue;
		if(middle.distance2(el->getMiddlePoint()) > r2) continue;
		view_set->addElementWithEdges(el, nullptr, true, 0);
		// log
		int edct = el->getEdgeCount();
		of << "Element index=" << el->getIndex() << " pts=" << edct << " [";
		for(int j = 0; j < edct; j++) of << el->getPoint(j)->getIntTag(TagExtended::TAG_ID) << " ";
		of << "]" << endl;
		of << "\tinverted=" << (el->isInverted()?"true":"false") << endl;
		of << "\tarea-012=" << DTriangle2d::area(
			el->getPoint(0)->getCoordinates(),
			el->getPoint(1)->getCoordinates(),
			el->getPoint(2)->getCoordinates()) << endl;
		if(edct == 4){
			of << "\tarea-123=" << DTriangle2d::area(
				el->getPoint(1)->getCoordinates(),
				el->getPoint(2)->getCoordinates(),
				el->getPoint(3)->getCoordinates()) << endl;
			of << "\tarea-230=" << DTriangle2d::area(
				el->getPoint(2)->getCoordinates(),
				el->getPoint(3)->getCoordinates(),
				el->getPoint(0)->getCoordinates()) << endl;
			of << "\tarea-301=" << DTriangle2d::area(
				el->getPoint(3)->getCoordinates(),
				el->getPoint(0)->getCoordinates(),
				el->getPoint(1)->getCoordinates()) << endl;
		}
	}

	// others
	for(int i = 0; i < pct; i++){
		MeshPoint2d* point = getPointAt(i);
		if(middle.distance2(point->getCoordinates()) > r2) continue;
		if(point == pt1)
			view_set->addPoint(point, nullptr, 1);
		else if(point == pt2)
			view_set->addPoint(point, nullptr, 2);
		else
			view_set->addPoint(point, nullptr, 0);
	}

	return view_set;
}

MeshEdge2d* MeshContainer2d::getCoupledEdge(const MeshEdge3d* edge3d) const
{
	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge()){
		if(it.getEdge()->getPtrTag(TagExtended::TAG_ME_2D_3D) == edge3d)
			return it.getEdge();
	}
	return nullptr;
}

void MeshContainer2d::countQuality(Metric2dContext& mc, int criterion){
	int ect = getElementsCount();
	for(int i = 0; i < ect; i++)
		getElementAt(i)->countQuality(mc, false, criterion);
}

void MeshContainer2d::countMetricDifferenceQuality(){
	int ect = getElementsCount();
	Metric2dContext mc(m_control);
	for(int i = 0; i < ect; i++)
		getElementAt(i)->countMetricDiffQuality(mc);
}

bool MeshContainer2d::storePJM(const string& fname, int grid_type, int id) const
{
	int pct = getPointsCount();
	int ect = getElementsCount();
	if(pct < 1 || ect < 1){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Empty mesh - nothing to write!");
		return false;
	}
	int tct = getElementsCount(3);
	int qct = getElementsCount(4);
	assert(tct + qct == ect);
	if(tct > 0 && qct > 0){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Mixed mesh (tri and quad) - don't know how to write!");
		return false;
	}

//	MeshElement* first_element = getElementAt(0);
	int edge_inner_nodes_count = 0; //first_element->getEdge(0)->getInnerPointsCount(); // TODO for inner nodes
	if(edge_inner_nodes_count > 1){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Too many inner nodes for edges - don't know how to write!");
		return false;
	}

	if(grid_type == 1){
		LOG4CPLUS_INFO(MeshLog::logger_console, "Storing grid GR_2D_Q4");
		if(tct > 0){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Mesh contains some triangles! Leaving."); return false; }
		if(edge_inner_nodes_count > 0){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Non-zero number of inner nodes for edges! Leaving."); return false; }
	}else if(grid_type == 3){
		LOG4CPLUS_INFO(MeshLog::logger_console, "Storing grid GR_2D_Q8");
		if(tct > 0){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Mesh contains some triangles! Leaving."); return false; }
		if(edge_inner_nodes_count != 1){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Number of inner nodes for edges different from 1! Leaving."); return false; }
	}else if(grid_type == 5){
		LOG4CPLUS_INFO(MeshLog::logger_console, "Storing grid GR_2D_Q9");
		if(tct > 0){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Mesh contains some triangles! Leaving."); return false; }
		if(edge_inner_nodes_count != 1){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Number of inner nodes for edges different from 1! Leaving."); return false; }
	}else if(grid_type == 7){
		LOG4CPLUS_INFO(MeshLog::logger_console, "Storing grid GR_2D_T3");
		if(qct > 0){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Mesh contains some quadrangles! Leaving."); return false; }
		if(edge_inner_nodes_count > 0){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Non-zero number of inner nodes for edges! Leaving."); return false; }
	}else if(grid_type == 8){
		LOG4CPLUS_INFO(MeshLog::logger_console, "Storing grid GR_2D_T6");
		if(qct > 0){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Mesh contains some quadrangles! Leaving."); return false; }
		if(edge_inner_nodes_count != 1){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Number of inner nodes for edges different from 1! Leaving."); return false; }
	}else{
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Unknown grid type! Leaving."); 
		return false;
	}

	ostringstream fname_desc, fname_data;
	fname_data << fname;
	if(id > 0) fname_data << '-' << id;
	fname_desc << fname;
	if(id > 0) fname_desc << '-' << id;
	fname_desc << ".desc";

	ofstream fdata(fname_data.str().c_str());
	ofstream fdesc(fname_desc.str().c_str());

	fdata << "GRID_TYPE " << grid_type << endl;
	DRect rect = getBoundingRect();
	fdata << "DIMENSIONS\t" << rect.getDX() << '\t' << rect.getDY() << endl;
//	fdata << "GRID_LEVELS\t0\t2" << endl << endl;

	fdata << "POINTS " << pct << endl << endl;
	for(int i = 0; i < pct; i++){
		DPoint2d pt = getPointAt(i)->getCoordinates();
		fdata << (abs(pt.x) < mesh_data.relative_small_number ? 0.0 : pt.x) << "\t" 
			  << (abs(pt.y) < mesh_data.relative_small_number ? 0.0 : pt.y) << endl;
	}

	if(grid_type == 5){ // GR_2D_Q9 -> add additional points within quads
		for(int i = 0; i < qct; i++){
			DPoint2d pt = getElementAt(i)->getMiddlePoint();
			fdata << pt.x << "\t" << pt.y << endl;
		}
	}


	fdata << endl << endl << "ELEMENTS \t" << ect << '\t' << pct << endl << endl;
	for(int i = 0; i < ect; i++){
		MeshElement* element = getElementAt(i);
		int edges_ct = element->getEdgeCount();
		for(int j = 0; j < edges_ct; j++){
			fdata << element->getPoint(j)->getIndex() << '\t';
// TODO for inner points
//			MeshEdge2d* edge = element->getEdge(j);
//			int ipct = edge->getInnerPointsCount();
//			if(ipct != edge_inner_nodes_count){
//				LOG4CPLUS_WARN(MeshLog::logger_console, "Different number of inner nodes for some edges!!!");
//			}
//			for(int k = 0; k < ipct; k++)
//				fdata << edge->getInnerMeshPoint(k)->getIndex() << '\t';
		}
		if(grid_type == 5) fdata << (i+pct) << '\t'; // GR_2D_Q9 -> add additional points within quads
		for(int j = 0; j < edges_ct; j++){
			MeshEdge2d* edge = element->getEdge(j);
			MeshElement* other_element = edge->getOtherElement(element);
			if(other_element){
				fdata << '\t' << other_element->getIndex();
			}else{
				MeshEdge3d* edge3d = (MeshEdge3d*)edge->getPtrTag(TagExtended::TAG_ME_2D_3D);
				if(edge3d){
					MeshBoundaryCondition* condition = (MeshBoundaryCondition*)edge3d->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
					if(condition) fdata << "\t" << condition->getCondition();
					else fdata << "\t?";
				}else fdata << "\t?";
			}
		}
		fdata << "\t" << element->getAreaID() << endl;
	}

	fdata << endl << "COLOR 1" << endl << ect << '\t';
	for(int i = 0; i < ect; i++) fdata << i << ' ';
	fdata << endl;

	fdesc << pct << " nodes" << endl;
	fdesc << tct << " triangles" << endl;
	fdesc << qct << " quads" << endl << endl;
	MeshData::StatData alpha = getAlphaQuality();
	fdesc << "Alpha quality (ave): " << alpha.average << endl;
	fdesc << "Alpha quality (min): " << alpha.minimum << endl;
	MeshData::StatData angles = getAngleQuality();
	fdesc << "Inner angle (min)  : " << angles.minimum*180.0/PI << "o" << endl;

	return true;
}

int MeshContainer2d::logMetricQuality(const string& caption)
{
	Metric2dContext mc(m_control);
	DataStatistics stats;

	DataVector<double> ranges(10);
	ranges.add(0.2);
	ranges.add(0.5);
	ranges.add(0.8);
	ranges.add(0.95);
	ranges.add(1.05);
	ranges.add(1.2);
	ranges.add(1.5);
	ranges.add(2.0);

	MeshGenerator2d::statMetricQuality(mc, this, stats, MeshData::QUALITY_SPACE);
	if(stats.calculate()) stats.logStats(caption+"MetricQuality - Triangle Area", "MQ-TS", ranges);
	stats.clear();
	MeshGenerator2d::statMetricQuality(mc, this, stats, MeshData::QUALITY_CIRCLE_SPACE);
	if(stats.calculate()) stats.logStats(caption+"MetricQuality - Circumcircle Area", "MQ-CS", ranges);
	stats.clear();
	MeshGenerator2d::statMetricQuality(mc, this, stats, MeshData::QUALITY_ALPHA);
	if(stats.calculate()) stats.logStats(caption+"MetricQuality - Alpha", "MQ-AL", ranges);
	stats.clear();
	statMetricEdgeLength(mc, stats, false);
	if(stats.calculate()) stats.logStats(caption+"MetricQuality - Edge Length", "MQ-EL", ranges);

	return 1;
}

bool MeshContainer2d::storeMatlabFile(const string& fname, const DPoint3d& clip, double quality_clip) const
{
	ofstream os(fname.c_str());

	os << "clear all; figure(1); clf; hold on; axis equal; axis off; view(3);" << endl;
	os << "c(1,:) = [0.9,0.9,0.9];" << endl;
	os << "c(2,:) = [0.7,0.7,0.7];" << endl;
	os << "colormap(c);" << endl;
	os << "load 'mesh_vertices.txt' -ascii" << endl;
	os << "load 'mesh_faces.txt' -ascii" << endl;
	os << "load 'mesh_ci.txt' -ascii" << endl;
	os << "patch('Vertices',mesh_vertices,'faces',mesh_faces,'FaceVertexCData',mesh_ci(:,1),'FaceColor','flat');" << endl;

	// store vertices
	ofstream ofsv("mesh_vertices.txt");
	int pct = getPointsCount();
	for(int j = 0; j < pct; j++){
		const DPoint3d& pt = m_surface->getPoint(getPointAt(j)->getCoordinates());
		ofsv << pt.x << ' ' << pt.y << ' ' << pt.z << endl;
	}

	// check elements
	int ect = getElementsCount();
	DataVector<bool> hidden(ect, true);
	for(int i = 0; i < ect; i++){
		MeshElement* element = getElementAt(i);
		if(quality_clip < 1.05 && element->getQuality() > quality_clip) continue;
		bool clipped = false;
		int bpct = element->getEdgeCount();
		for(int j = 0; j < bpct; j++){
			const DPoint3d& bpt = m_surface->getPoint(element->getPoint(j)->getCoordinates());
			if(bpt.x > clip.x || bpt.y > clip.y || bpt.z > clip.z){
				clipped = true; break;
			}
		}
		if(!clipped) hidden[i] = false;
	}

	ofstream ofsf("mesh_faces.txt");
	ofstream ofsc("mesh_ci.txt");

	// store visible (not-clipped) elements
	int total_fct = 0;
	for(int i = 0; i < ect; i++){
		if(hidden[i]) continue;
		MeshElement* element = getElementAt(i);
		// store faces
		if(total_fct++ > 0) os << "; ";
		int epct = element->getEdgeCount();
		for(int k = 0; k < epct; k++)
			ofsf << (1+ element->getPoint(k)->getIndex()) << ' ';
		ofsf << endl;
		ofsc << 1 << endl;
	}

	return true;
}

/// Partition mesh (remove marked elements and transfer them to other mesh)
MeshContainer2d* MeshContainer2d::splitByElements(int split_area_id)
{
	int element_count = getElementsCount();
	MeshContainer2d* other_mesh = new MeshContainer2d(element_count/2);
	bool heap_order = isHeapOrder();
	if(heap_order) setHeapOrder(false);
	// split elements
	for(int i = 0; i < element_count; ){
		MeshElement* element = getElementAt(i);
		if(element->getAreaID() == split_area_id){
			other_mesh->addMeshElement(removeMeshElement(i));
			--element_count;
		}else ++i;
	}
	int pct = getPointsCount();
	DataHashTable<MeshEdge2d*> cut_edges(3*pct, 0);
	DataHashTable<MeshPoint2d*> near_cut_points(pct, 0);

	// identify split contour
	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge()){
		MeshEdge2d* edge = it.getEdge();
		MeshElement* element0 = edge->getMeshElement(0);
		MeshElement* element1 = edge->getMeshElement(1);
		MeshPoint2d* points[2] = { edge->getMeshPoint(0), edge->getMeshPoint(1) };
		if(element0 && element1 && (element0->getAreaID() != element1->getAreaID())){
			cut_edges.insert(edge);
			MeshPoint2d* mpoints[2] = { nullptr, nullptr };
			for(int j = 0; j < 2; j++){
				mpoints[j] = (MeshPoint2d*) points[j]->getPtrTag(TagExtended::TAG_CUT_MP_2D);
				if(!mpoints[j]){ // clone point
					mpoints[j] = new MeshPoint2d(points[j]);
					other_mesh->addMeshPoint(mpoints[j]);
					points[j]->setPtrTag(TagExtended::TAG_CUT_MP_2D, mpoints[j]);
					mpoints[j]->setPtrTag(TagExtended::TAG_CUT_MP_2D, points[j]);
				}
			}
			MeshEdge2d* medge = new MeshEdge2d(mpoints[0], mpoints[1]);
			// TODO copy edge parameters or create some new type of edge for inter-partition connections
			if(edge->isBorder())
				medge->copyBorderFlagsFrom(edge);
			else{
				edge->setBorder(TagBorder::OUTER | TagBorder::FIXED);
				medge->setBorder(TagBorder::OUTER | TagBorder::FIXED);
			}
			points[0]->setBorder(TagBorder::OUTER | TagBorder::FIXED);
			points[1]->setBorder(TagBorder::OUTER | TagBorder::FIXED);
			mpoints[0]->setBorder(TagBorder::OUTER | TagBorder::FIXED);
			mpoints[1]->setBorder(TagBorder::OUTER | TagBorder::FIXED);
		}
	}
	// move inner points (of split area)
	DataSimpleList<EdgeInfo> border_edges_near_split;
	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge()){
		MeshEdge2d* edge = it.getEdge();
		MeshElement* element0 = edge->getMeshElement(0);
		MeshElement* element1 = edge->getMeshElement(1);
		MeshPoint2d* points[2] = { edge->getMeshPoint(0), edge->getMeshPoint(1) };
		if(edge->isBorder() && 
			((points[0]->availableTag(TagExtended::TAG_CUT_MP_2D) && !points[1]->availableTag(TagExtended::TAG_CUT_MP_2D)) ||
			(!points[0]->availableTag(TagExtended::TAG_CUT_MP_2D) && points[1]->availableTag(TagExtended::TAG_CUT_MP_2D))))
		{
			border_edges_near_split.insert(EdgeInfo(points[0], points[1], edge->getBorderFlags()));
			edge->clearBorder();
		}
		if(!cut_edges.contains(edge)){
			int id = element0 ? element0->getAreaID() : element1->getAreaID();
			if(id == split_area_id){
				for(int j = 0; j < 2; j++){
					if(!points[j]->availableTag(TagExtended::TAG_CUT_MP_2D)){
						other_mesh->addMeshPoint(removeMeshPoint(points[j]));
						near_cut_points.insert(points[j]);
					}
				}
			}
		}
	}

	// adjust incidency for split contour
	int split_element_count = other_mesh->getElementsCount();
	for(int i = 0; i < split_element_count; i++){
		MeshElement* element = other_mesh->getElementAt(i);
		int ect = element->getEdgeCount();
		for(int j = 0; j < ect; j++){
			MeshPoint2d* point = element->getPoint(j);
			if(point->availableTag(TagExtended::TAG_CUT_MP_2D)){ // i.e. split contour (but original points)
				element->switchPointsWithEdges(point, (MeshPoint2d*)(point->getPtrTag(TagExtended::TAG_CUT_MP_2D)));
			}
		}
	}

	// adjust border edges near split contour
	while(border_edges_near_split.notEmpty()){
		EdgeInfo info = border_edges_near_split.removeFirst();
		if(info.pt0->availableTag(TagExtended::TAG_CUT_MP_2D)){ // first point is on contour
			if(near_cut_points.contains(info.pt1)) // other_mesh
				info.pt0 = (MeshPoint2d*)(info.pt0->getPtrTag(TagExtended::TAG_CUT_MP_2D));
			info.pt0->getEdgeToPoint(info.pt1)->setBorder(info.btype);
		}else{ // second point is on the contour
			if(near_cut_points.contains(info.pt0)) // other_mesh
				info.pt1 = (MeshPoint2d*)(info.pt1->getPtrTag(TagExtended::TAG_CUT_MP_2D));
			info.pt0->getEdgeToPoint(info.pt1)->setBorder(info.btype);
		}
	}

	other_mesh->setConstrainingPhase(m_constraining);
	other_mesh->setControlSpace(m_control);
	other_mesh->setDiscretizationState(m_discretization_state);
	other_mesh->setSurface(m_surface);

	if(heap_order){
		setHeapOrder(true);
		other_mesh->setHeapOrder(true);
	}

	return other_mesh;
}

	/// Partition mesh (remove marked elements and transfer them to other mesh)
MeshContainer2d* MeshContainer2d::splitByElements(TagExtended::TagType tag_type, int tag_value)
{
	int element_count = getElementsCount();
	MeshContainer2d* other_mesh = new MeshContainer2d(element_count/2);
	bool heap_order = isHeapOrder();
	if(heap_order) setHeapOrder(false);
	// split elements
	for(int i = 0; i < element_count; ){
		MeshElement* element = getElementAt(i);
		if(element->getIntTag(tag_type, tag_value-1) == tag_value){
			other_mesh->addMeshElement(removeMeshElement(i));
			--element_count;
		}else ++i;
	}
	int pct = getPointsCount();
	DataHashTable<MeshEdge2d*> cut_edges(3*pct, 0);
	DataHashTable<MeshPoint2d*> near_cut_points(pct, 0);

	// identify split contour
	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge()){
		MeshEdge2d* edge = it.getEdge();
		MeshElement* element0 = edge->getMeshElement(0);
		MeshElement* element1 = edge->getMeshElement(1);
		if(!element0 || !element1) continue;
		bool cut_element0 = (element0->getIntTag(tag_type, tag_value-1) == tag_value);
		bool cut_element1 = (element1->getIntTag(tag_type, tag_value-1) == tag_value);
		MeshPoint2d* points[2] = { edge->getMeshPoint(0), edge->getMeshPoint(1) };
		if(cut_element0 != cut_element1){
			cut_edges.insert(edge);
			MeshPoint2d* mpoints[2] = { nullptr, nullptr };
			for(int j = 0; j < 2; j++){
				mpoints[j] = (MeshPoint2d*) points[j]->getPtrTag(TagExtended::TAG_CUT_MP_2D);
				if(!mpoints[j]){ // clone point
					mpoints[j] = new MeshPoint2d(points[j]);
					other_mesh->addMeshPoint(mpoints[j]);
					points[j]->setPtrTag(TagExtended::TAG_CUT_MP_2D, mpoints[j]);
					mpoints[j]->setPtrTag(TagExtended::TAG_CUT_MP_2D, points[j]);
				}
			}
			MeshEdge2d* medge = new MeshEdge2d(mpoints[0], mpoints[1]);
			medge->setPtrTag(TagExtended::TAG_CUT_MP_2D, edge);
			edge->setPtrTag(TagExtended::TAG_CUT_MP_2D, medge);
			// TODO copy edge parameters or create some new type of edge for inter-partition connections
			edge->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, edge->getBorderFlags());
			medge->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, edge->getBorderFlags());
			if(edge->isBorder())
				medge->copyBorderFlagsFrom(edge);
			else{
				edge->setBorder();
				medge->setBorder();
			}
			// save
			points[0]->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, points[0]->getBorderFlags());
			points[1]->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, points[1]->getBorderFlags());
			mpoints[0]->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, points[0]->getBorderFlags());
			mpoints[1]->setIntTag(TagExtended::TAG_CUT_ORYG_DATA, points[1]->getBorderFlags());
			// change
			points[0]->setBorder(TagBorder::OUTER | TagBorder::FIXED);
			points[1]->setBorder(TagBorder::OUTER | TagBorder::FIXED);
			mpoints[0]->setBorder(TagBorder::OUTER | TagBorder::FIXED);
			mpoints[1]->setBorder(TagBorder::OUTER | TagBorder::FIXED);
		}
	}
	// move inner points (of split area)
	DataSimpleList<EdgeInfo> border_edges_near_split;
	DataSimpleList<MeshPoint2d*> points_to_move;
	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge()){
		MeshEdge2d* edge = it.getEdge();
		MeshElement* element0 = edge->getMeshElement(0);
		MeshElement* element1 = edge->getMeshElement(1);
		MeshPoint2d* points[2] = { edge->getMeshPoint(0), edge->getMeshPoint(1) };
		if(edge->isBorder() && 
			((points[0]->availableTag(TagExtended::TAG_CUT_MP_2D) && !points[1]->availableTag(TagExtended::TAG_CUT_MP_2D)) ||
			(!points[0]->availableTag(TagExtended::TAG_CUT_MP_2D) && points[1]->availableTag(TagExtended::TAG_CUT_MP_2D))))
		{
			border_edges_near_split.insert(EdgeInfo(points[0], points[1], edge->getBorderFlags()));
			edge->clearBorder();
		}
		if(!cut_edges.contains(edge)){
			MeshElement* valid_element = element0 ? element0 : element1;
			if(valid_element->getIntTag(tag_type, tag_value-1) == tag_value){
				for(int j = 0; j < 2; j++){
					if(!points[j]->availableTag(TagExtended::TAG_CUT_MP_2D) &&
						!near_cut_points.contains(points[j])){
						near_cut_points.insert(points[j]);
						points_to_move.append(points[j]);
					}
				}
			}
		}
	}
	while(!points_to_move.empty()){
		MeshPoint2d* point = points_to_move.removeFirst();
		other_mesh->addMeshPoint(removeMeshPoint(point));
	}

	// adjust incidency for split contour
	int split_element_count = other_mesh->getElementsCount();
	for(int i = 0; i < split_element_count; i++){
		MeshElement* element = other_mesh->getElementAt(i);
		int ect = element->getEdgeCount();
		for(int j = 0; j < ect; j++){
			MeshPoint2d* point = element->getPoint(j);
			if(point->availableTag(TagExtended::TAG_CUT_MP_2D)){ // i.e. split contour (but original points)
				element->switchPointsWithEdges(point, (MeshPoint2d*)(point->getPtrTag(TagExtended::TAG_CUT_MP_2D)));
			}
		}
	}

	// adjust border edges near split contour
	while(border_edges_near_split.notEmpty()){
		EdgeInfo info = border_edges_near_split.removeFirst();
		if(info.pt0->availableTag(TagExtended::TAG_CUT_MP_2D)){ // first point is on contour
			if(near_cut_points.contains(info.pt1)) // other_mesh
				info.pt0 = (MeshPoint2d*)(info.pt0->getPtrTag(TagExtended::TAG_CUT_MP_2D));
			info.pt0->getEdgeToPoint(info.pt1)->setBorder(info.btype);
		}else{ // second point is on the contour
			if(near_cut_points.contains(info.pt0)) // other_mesh
				info.pt1 = (MeshPoint2d*)(info.pt1->getPtrTag(TagExtended::TAG_CUT_MP_2D));
			info.pt0->getEdgeToPoint(info.pt1)->setBorder(info.btype);
		}
	}

	other_mesh->setConstrainingPhase(m_constraining);
	other_mesh->setControlSpace(m_control);
	other_mesh->setDiscretizationState(m_discretization_state);
	other_mesh->setSurface(m_surface);

	//if(heap_order){
	//	setHeapOrder(true);
	//	other_mesh->setHeapOrder(true);
	//}

	return other_mesh;
}

/// converts the mesh into a boundary-mesh only (removes all inner nodes and elements)
void MeshContainer2d::convertToBoundaryMesh()
{
	DataHashTableKeyValue<int, MeshArea*> areas(100, -1);
	DataVector<MeshArea*> new_areas(100);
	// remove obsolete elements
	while(getElementsCount() > 0){
		MeshElement* element = removeMeshElement(0);
		int area_id = element->getAreaID();
		MeshArea* area = areas.getValue(area_id, nullptr);
		if(!area){
			area = new MeshArea(); // dummy
			area->setAreaID(area_id);
			areas.insert(area_id, area);
			new_areas.add(area);
		}
		for(int i = 0; i < element->getEdgeCount(); i++){
			MeshEdge2d* edge = element->getEdge(i);
			if(edge->isBorder()){
				edge->switchElementLink(element, area);
			}
		}
		delete element;
	}
	// insert new (dummy) areas
	for(size_t i = 0; i < new_areas.countInt(); i++)
		addMeshElement(new_areas[i]);

	// remove obsolete (inner) points
	for(int i = 0; i < getPointsCount(); ){
		MeshPoint2d* point = getPointAt(i);
		if(point->isBorder()){
			++i;
		}else{
			removeMeshPoint(i);
		}
	}
}

/// merge two meshes into one (the mesh from argument will be erased!)
bool MeshContainer2d::merge(MeshContainer2d* &mesh)
{
	// move elements
	while(mesh->getElementsCount() > 0){
		MeshElement* element = mesh->removeMeshElement(0);
		addMeshElement(element);
		int ect = element->getEdgeCount();
		// check edges
		for(int i = 0; i < ect; i++){
			MeshEdge2d* edge = element->getEdge(i);
			if(edge->isBorder()){
				MeshEdge2d* bedge = (MeshEdge2d*)edge->getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
				if(!bedge) continue;
				MeshEdge2d* cut_edge = (MeshEdge2d*)bedge->getPtrTag(TagExtended::TAG_CUT_MP_2D);
				if(cut_edge){
					edge->clearBorder(); // marked as non-border in order to be removed automatically later
					cut_edge->setBorder(cut_edge->getIntTag(TagExtended::TAG_CUT_ORYG_DATA, TagBorder::NONE));
				}
			}
		}
		// check points
		for(int i = 0; i < ect; i++){
			MeshPoint2d* point = element->getPoint(i);
			if(!point->isBorder()) continue;
			MeshPoint2d* bpoint = (MeshPoint2d*)point->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
			if(!bpoint) continue;
			MeshPoint2d* cut_point = (MeshPoint2d*)bpoint->getPtrTag(TagExtended::TAG_CUT_MP_2D);
			if(cut_point){
				cut_point->setBorder((char)cut_point->getIntTag(TagExtended::TAG_CUT_ORYG_DATA, TagBorder::NONE));
				element->switchPointsWithEdges(point, cut_point);
				if(point->getRank() == 0)
					delete mesh->removeMeshPoint(point);
			}
		}
	}
	// move inner and other non-cut points
	while(mesh->getPointsCount() > 0){
		addMeshPoint(mesh->removeMeshPoint(0));
	}
	// clean
	delete mesh;
	mesh = nullptr;
	return true;
}
