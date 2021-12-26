/////////////////////////////////////////////////////////////////////////////
// MeshEdge2d.cpp
// Reprezentuje krawêdŸ ³¹cz¹ca dwa wierzcho³ki w siatce
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////

#include "DPoint.h"
#include "MeshEdge2d.h"
#include "MeshPoint3d.h"
#include "MeshPoint2d.h"
#include "MeshElement.h"
#include "MeshTriangle2d.h"
#include "SurfaceParametric.h"
#include "Curve2dSegment.h"
#include "ControlSpace2d.h"
#include "MeshViewSet.h"
#include "MeshEdge3d.h"

//////////////////////////////////////////////////////////////////////
// Standardowy konstruktor
MeshEdge2d::MeshEdge2d(MeshPoint2d *p1, MeshPoint2d *p2, MeshElement* e1, MeshElement* e2)
{
	assert(p1);
	assert(p2);
	assert(p1 != p2);
	points[0] = p1;
	points[1] = p2;
	elements[0] = e1;
	elements[1] = e2;
	p1->addEdge(this);
	p2->addEdge(this);
}

//////////////////////////////////////////////////////////////////////
// Destruktor
MeshEdge2d::~MeshEdge2d()
{
	if(points[0]){ // no delete-all mode
		points[0]->removeEdge(this);
		points[1]->removeEdge(this);
	}
}

/// Removes link to this point (fo delete-all phase) - returns true if last one
bool MeshEdge2d::removePointLink(MeshPoint2d* point)
{
	if(points[0] == point) points[0] = 0;
	else if(points[1] == point) points[1] = 0;

	return (points[0] == 0 && points[1] == 0);
}

/// Create new edge with similar geometry
MeshEdge2d* MeshEdge2d::cloneGeometric(MeshPoint2d *p0, MeshPoint2d *p1, double /* ksi0 */, double /* ksi1 */, MeshElement* e1, MeshElement* e2) const
{
	return new MeshEdge2d(p0, p1, e1, e2);
}

//////////////////////////////////////////////////////////////////////
// Usuwa dowi¹zanie do podanego elementu, zwraca true jeœli po tej
//	operacji krawêdŸ mo¿na usun¹æ
bool MeshEdge2d::removeElementLink(const MeshElement *element)
{
	int i = (elements[0] == element)?0:1;
	if(elements[i] != element) return false;

	elements[i] = nullptr;

	// Jeœli krawêdŸ nie jest brzegowa i nie jest zwi¹zana ju¿ z ¿adnym obiektem, 
	//		mo¿na j¹ usun¹æ (true)
	return (!elements[1-i] && !isBorder());
}

//////////////////////////////////////////////////////////////////////
// Dodaje dowi¹zanie do danego elementu po lewej stronie tej krawêdzi, 
//	przyjmuj¹c ¿e krawêdŸ jest skierowana od punktu _p1_.
void MeshEdge2d::addElementLink(MeshElement *element, const MeshPoint2d *p1)
{
	if(points[0] == p1){
		elements[0] = element;
	}else{
		elements[1] = element;
	}
}

//////////////////////////////////////////////////////////////////////
// Ustawia kierunek krawêdzi wzglêdem punktu _point_ 
//	1 - skierowana dodatnio (od punktu)
//	0 - nieskierowana
// -1 - skierowana ujemnie (do punktu)
// -2 - nieokreœlona (nie brzegowa)
void MeshEdge2d::setDirection(int dir, const MeshPoint2d *point)
{
	if(dir == -2){
		clearBorder();
		return;
	}
	if(dir == 0){
		setBorder(TagBorder::INNER);
	}else{
		setBorder(TagBorder::OUTER);
		int i = (dir == 1)?0:1;
		if(point != points[i]){
			assert(point == points[1-i]);
			// Zamiana stron
			switchSide();
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Ustawia kierunek krawêdzi z uwzglêdnieniem jej poprzedniego
//	 kierunku (krawêdŸ skierowana w obie strony staje siê nieskierowana)
void MeshEdge2d::addDirection(int dir, const MeshPoint2d *point)
{
	if(dir != 0){
		int i = (dir == 1)?0:1;
		if(isBorder(TagBorder::OUTER)){
			if(point == points[1-i]) setBorder(TagBorder::INNER);
		}else{
			setBorder(TagBorder::OUTER);
			if(point != points[i]){
				assert(point == points[1-i]);
				switchSide();
			}
		}
	}
}

double MeshEdge2d::getLengthQuality(Metric2dContext& mc, bool ext_metric) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->getLengthQuality(mc, ext_metric);
	}
	int qms = MeshTriangle2d::param_quality_metric_selection;
	MeshTriangle2d::param_quality_metric_selection = 
		ext_metric ? MeshData::QM_VERTICES_AVE : MeshData::QM_MIDDLE;
	mc.countMetricAtPoints(points[0], points[1]);
	MeshTriangle2d::param_quality_metric_selection = qms;

	return points[0]->getMetricCoordinates(mc).distance(
			points[1]->getMetricCoordinates(mc));
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca d³ugoœæ ca³kowit¹ krawêdzi 
double MeshEdge2d::getLength(Metric2dContext& mc, bool local_metric) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->getLength(mc, local_metric);
	}
	if(local_metric) mc.countMetricAtPoint(getPoint(0.5));
	return points[0]->getMetricCoordinates(mc).distance(
		points[1]->getMetricCoordinates(mc));
}

double MeshEdge2d::getLengthMetricAdapted(Metric2dContext& mc) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if (boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) {
			return boundary_edge->getLengthMetricAdapted(mc);
		}
	}
	return getLengthMetricAdapted(mc, 0.0, 1.0);
}

double MeshEdge2d::getLengthMetricAdapted(Metric2dContext& mc, double ksi0, double ksi1, int lev) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->getLengthMetricAdapted(mc, ksi0, ksi1, lev);
	}
	const DPoint2d pt0 = getPoint(ksi0);
	const DPoint2d pt1 = getPoint(ksi1);
	double ksi_m = 0.5 * (ksi0 + ksi1);
	const DPoint2d pt_m = getPoint(ksi_m);
	mc.countMetricAtPoint(DPoint2d::average(pt0, pt_m));
	double len_0 = mc.transformPStoMS(pt_m-pt0).length();
	mc.countMetricAtPoint(DPoint2d::average(pt1, pt_m));
	len_0 += mc.transformPStoMS(pt_m-pt1).length();
	if(lev > 5) return len_0;
	mc.countMetricAtPoint(pt_m);
	double len_1 = mc.transformPStoMS(pt1-pt0).length();
	if((len_0 - len_1) < 1e-3*len_1) return len_0;
	else return getLengthMetricAdapted(mc, ksi0, ksi_m, lev+1) +
		getLengthMetricAdapted(mc, ksi_m, ksi1, lev+1);
}

/// Returns the length of this edge (for variable metric)
double MeshEdge2d::getRequiredLength(Metric2dContext& mc, double ksi0, double &ksi1, double req_length) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->getRequiredLength(mc, ksi0, ksi1, req_length);
	}
	assert(ksi1 > ksi0);
	const int nct = 10;
	const double dksi = (ksi1-ksi0) / nct;
	double len = 0.0;
	ksi1 = ksi0;
	for(int i = 0; i < nct; i++){
		double len_i = getLengthMetricAdapted(mc, ksi1, ksi1+dksi);
		if(len+len_i > req_length){
			if(i > 0) return len;
			else{
				ksi1 = ksi0+dksi;
				return getRequiredLength(mc, ksi0, ksi1, req_length);
			}
		}else{
			ksi1 += dksi;
			len += len_i;
		}
	}
	return len;
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca d³ugoœæ czêœci krawêdzi ograniczonej przez parametry t0 i t1.
//	Parametry nale¿¹ do przedzia³u [0, 1]
double MeshEdge2d::getLength(Metric2dContext& mc, double ksi0, double ksi1, bool local_metric) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->getLength(mc, ksi0, ksi1, local_metric);
	}

	if(local_metric) 
		mc.countMetricAtPoint(getPoint((ksi0+ksi1)*0.5));
	return mc.transformPStoMS(getPoint(ksi0)-getPoint(ksi1)).length();
}

double MeshEdge2d::getLengthMax(Metric2dContext& mc, double ksi0, double ksi1) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->getLengthMax(mc, ksi0, ksi1);
	}

	const DPoint2d pt1 = getPoint(ksi0);
	const DPoint2d pt2 = getPoint(ksi1);
	mc.countMetricAtPoint(pt1);
	double len1 = mc.transformPStoMS(pt1-pt2).length();
	mc.countMetricAtPoint(pt2);
	double len2 = mc.transformPStoMS(pt1-pt2).length();
	return std::max(len1, len2);
}

double MeshEdge2d::checkAndGetLength(Metric2dContext& mc, double ksi0, double& ksi1, double max_len, bool local_metric) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->checkAndGetLength(mc, ksi0, ksi1, max_len, local_metric);
	}

	if(local_metric) mc.countMetricAtPoint(getPoint((ksi0+ksi1)*0.5));
	const DMPoint2d pt0 = mc.transformPStoMS(getPoint(ksi0));
	const DMPoint2d pt1 = mc.transformPStoMS(getPoint(ksi1));
	double len = pt0.distance(pt1);
	if(len > max_len){
		ksi1 = ksi0 + max_len / len * (ksi1-ksi0);
		len = max_len;
		double diff = pt0.distance(mc.transformPStoMS(getPoint(ksi1))) - max_len;
		assert(abs(diff*diff) < 1e-5);
		if(abs(diff) < mesh_data.relative_small_number)
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Warning. Edge discretize -> diff = " << diff);
	}
	return len;
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca punkt nale¿¹cy do krawêdzi na podstawie wartoœci parametru t [0, 1]
DPoint2d MeshEdge2d::getPoint(double ksi) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->getPoint(ksi);
	}

	return DPoint2d(points[0]->getCoordinates(), points[1]->getCoordinates(), ksi);
}

//////////////////////////////////////////////////////////////////////
// Zwraca identyfikator obszaru po zadanej stronie krawêdzi
int MeshEdge2d::getIncidentAreaID(int side) const
{
	int i = (side == AREA_LEFT)?0:1;
	if(!elements[i]){
		return -1;
	}else{
		return elements[i]->getAreaID();
	}
}

//////////////////////////////////////////////////////////////////////
// Zwraca kierunek krawêdzi wzglêdem zadanego punktu
int MeshEdge2d::getDirection(const MeshPoint2d *point) const
{
	if(!isBorder()){
		return -2;
	}else if(isBorder(TagBorder::OUTER)){
		return (point == points[0])?1:-1;
	}else{
		return 0;
	}
}

//////////////////////////////////////////////////////////////////////
// Ustawia identyfikator obszaru dla elementu po zadanej stronie
//	krawêdzi
void MeshEdge2d::setIncidentAreaID(int side, int id)
{
	int i = (side == AREA_LEFT)?0:1;
	if(elements[i]){
		elements[i]->setAreaID(id);
	}
}

//////////////////////////////////////////////////////////////////////
// Zwraca pozycjê podanego punktu w tej krawêdzi
int MeshEdge2d::getPointIndex(const MeshPoint2d * point) const
{
	if(points[0] == point){
		return 0;
	}else if(points[1] == point){
		return 1;
	}else
		return -1;
}

//////////////////////////////////////////////////////////////////////
// Zmienia kierunek krawêdzi na przeciwny (zamienia elementy zwi¹zane
//	z obywdwoma stronami)
void MeshEdge2d::switchSide()
{
	MeshPoint2d* point = points[0];
	points[0] = points[1];
	points[1] = point;
	MeshElement* e = elements[0];
	elements[0] = elements[1];
	elements[1] = e;
	//clearSideTag();	// Wyczyszczenie znaczników
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca tablicê punktów definiuj¹cych krzyw¹ ³aman¹ oPIsuj¹c¹ kszta³t tej 
//	krawêdzi
void MeshEdge2d::getPolyLine(DataVector<DPoint2d> & polyline) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE){
			boundary_edge->getPolyLine(polyline);
			return;
		}
	}

	polyline.add(points[0]->getCoordinates());
	polyline.add(points[1]->getCoordinates());
}

/// Returns the approximation of this edge via polyline (array of points)
void MeshEdge2d::getPolyLine(DataVector<DPoint3d> & polyline, SurfaceConstPtr surface) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE){
			boundary_edge->getPolyLine(polyline, surface);
			return;
		}
		if(surface->getType() == SURFACE_PLANE){
			MeshEdge3d* edge3d = (MeshEdge3d*)getPtrTag(TagExtended::TAG_ME_2D_3D);
			if(edge3d){
				polyline.add(edge3d->getMeshPoint(0)->getCoordinates());
				polyline.add(edge3d->getMeshPoint(1)->getCoordinates());
				return;
			}
		}
	}

	surface->getPolyLine(polyline, points[0]->getCoordinates(), points[1]->getCoordinates());
}

void MeshEdge2d::switchPoints(const MeshPoint2d *point1, MeshPoint2d *point2)
{
	if(points[0] == point1)
		points[0] = point2;
	else if(points[1] == point1)
		points[1] = point2;
}

MeshElement* MeshEdge2d::getMeshElement(const MeshPoint2d *point) const
{
	if(points[0] == point){
		return elements[0];
	}else if(points[1] == point){
		return elements[1];
	}

	return nullptr;
}

double MeshEdge2d::getElementsWeight() const
{
	double w = 0.0;
	int ct = 0;

	int weights[] = {0, 0, 0, 1, 2};
	if(elements[0]){ 
		w += weights[elements[0]->getEdgeCount()];
		++ct;
	}
	if(elements[1]){
		w += weights[elements[1]->getEdgeCount()];
		++ct;
	}
	if(ct > 0) w /= ct;

	return w;
}

double MeshEdge2d::getLength(SurfaceConstPtr surface, double ksi0, double ksi1) const
{
	double len = 0.0;
	
	double factor = 1.0 - MEASURE_PRECISION;
	double heap[20];
	DPoint3d pt0 = surface->getPoint(getPoint(ksi0));
	DPoint3d pt1 = surface->getPoint(getPoint(ksi1));
	heap[0] = ksi1;
	int top = 0;
	while(true){
		double ksi_middle = (ksi0+ksi1)*0.5;
		DPoint3d pt_middle = surface->getPoint(getPoint(ksi_middle));
		double distance = pt0.distance(pt1);
		double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
		if((distance >= factor * real_distance) || (top == 19)){
			len += real_distance;
			if(top == 0) break;
			ksi0 = ksi1;
			pt0 = pt1;
			pt1 = surface->getPoint(getPoint(ksi1 = heap[--top]));
		}else{						
			heap[++top] = ksi1 = ksi_middle;
			pt1 = pt_middle;
		}
	}	
	return len;
}

std::shared_ptr<MeshViewEdgeData> MeshEdge2d::getViewData(SurfaceConstPtr surface) const
{
	DPoint3d dpts[2];
	for(int i = 0; i < 2; i++){
		if(!surface)
			dpts[i] = DPoint3d(points[i]->getCoordinates(), 0.0);
		else{
			MeshPoint3d* mpt3 = (MeshPoint3d*)points[i]->getPtrTag(TAG_MP_2D_3D);
			if(mpt3) dpts[i] = mpt3->getCoordinates();
			else dpts[i] = surface->getPoint(points[i]->getCoordinates());
		}
	}
	return std::make_shared<MeshViewEdgeData>(dpts[0], dpts[1], isBorder());
}

void MeshEdge2d::addToBoundingRect(DRect& rect) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE){
			boundary_edge->addToBoundingRect(rect);
			return;
		}
	}

	rect.addPoint(points[0]->getCoordinates());
	rect.addPoint(points[1]->getCoordinates());
}

/// numerical approximation -> calculate parameter for point on surface+curve
DPoint2d MeshEdge2d::surfaceParameters(SurfaceConstPtr surface, const DPoint3d& pt, double & ksi, bool same_direction) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->surfaceParameters(surface, pt, ksi, same_direction);
	}

	const DPoint2d & pt_0 = getMeshPoint(same_direction ? 0 : 1)->getCoordinates();
	const DPoint2d & pt_1 = getMeshPoint(same_direction ? 1 : 0)->getCoordinates();
	ksi = surface->getSegmentParameters(pt, pt_0, pt_1, ksi, 0.0, 1.0);
	return DPoint2d(pt_0, pt_1, ksi);
}

/// Returns the parameter ksi for a given point (numerical approximation)
double MeshEdge2d::getParameter(const DPoint2d& pt) const
{
	return Curve2dSegment( points[0]->getCoordinates(), points[1]->getCoordinates()).getParameter(pt, 0.5);
}

double MeshEdge2d::getNonPlanarCurvature(SurfaceConstPtr surface, double ksi, double* cdt_len) const
{
	if(isBorder()){
		MeshEdge2d* boundary_edge = (MeshEdge2d*)getPtrTag(TagExtended::TAG_BOUNDARY_EDGE);
		if(boundary_edge && boundary_edge->getType() != EDGE_SIMPLE) 
			return boundary_edge->getNonPlanarCurvature(surface, ksi, cdt_len);
	}

	const DPoint2d pt  = getPoint(ksi);
	const DVector3d fs  = surface->getDerivative(DEquation::deriv_ds, pt);
	const DVector3d ft  = surface->getDerivative(DEquation::deriv_dt, pt);

	// check
	double l1 = fs.length2();
	double l2 = ft.length2();
	double g_ratio = std::min(l1,l2) / std::max(l1,l2);
	if(g_ratio < MIN_PARAM_GRATIO) return 0.0;

	const DVector3d fss = surface->getDerivative(DEquation::deriv_dss, pt);
	const DVector3d fst = surface->getDerivative(DEquation::deriv_dst, pt);
	const DVector3d ftt = surface->getDerivative(DEquation::deriv_dtt, pt);
	const DVector2d dt  = points[1]->getCoordinates() - points[0]->getCoordinates();

	const DVector3d cdt = fs * dt.x + ft * dt.y;
	const DVector3d cdtt = fss * (dt.x*dt.x) + fst * (2*dt.x*dt.y) + ftt * (dt.y*dt.y);

	// curvature
	double ca = abs(cdt.crossProduct(cdtt).length());
	double cb = cdt.length2();
	assert(cb > 0.0);
	if(cdt_len){
		cb *= (*cdt_len = sqrt(cb));
	}else cb *= sqrt(cb);
	return /* cr= */ ca / cb;
}

/// Replaces one element link with another
bool MeshEdge2d::switchElementLink(const MeshElement* old_element, MeshElement* new_element)
{
	if(elements[0] == old_element){
		elements[0] = new_element; 
		return true;
	}else if(elements[1] == old_element){
		elements[1] = new_element;
		return true;
	}else
		return false;
}
