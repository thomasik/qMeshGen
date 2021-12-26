/////////////////////////////////////////////////////////////////////////////
// MeshTriangle2d.cpp
// Derived from MeshElement - implements triangular mesh element
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "MeshTriangle2d.h"

#include "common.h"
#include "MeshEdge2d.h"
#include "MeshData.h"
#include "MeshGenerator2d.h"
#include "ControlSpace2d.h"
#include "DMetric2d.h"
#include "MeshViewSet.h"
#include "DTriangle.h"

int MeshTriangle2d::param_quality_criterion = MeshData::QUALITY_CIRCLE_SPACE_AND_EDGES;
int MeshTriangle2d::param_quality_metric_selection = MeshData::QM_MIDDLE;
//unsigned int MeshTriangle2d::m_counters[3] = {0, 0, 0};

//////////////////////////////////////////////////////////////////////
// Konstruktor
MeshTriangle2d::MeshTriangle2d(MeshPoint2d *p1, MeshPoint2d *p2, MeshPoint2d *p3) :
	MeshElement(3, triangle_edges, triangle_points)
{
	points[0] = p1;
	points[1] = p2;
	points[2] = p3;

	assert(this->isOK());

	edges[0] = p1->getEdgeToPoint(p2);
	if(!edges[0]) edges[0] = new MeshEdge2d(p1, p2, this);
	else edges[0]->addElementLink(this, p1);

	edges[1] = p2->getEdgeToPoint(p3);
	if(!edges[1]) edges[1] = new MeshEdge2d(p2, p3, this);
	else edges[1]->addElementLink(this, p2);

	edges[2] = p3->getEdgeToPoint(p1);
	if(!edges[2]) edges[2] = new MeshEdge2d(p3, p1, this);
	else edges[2]->addElementLink(this, p3);
}

//////////////////////////////////////////////////////////////////////
// Destructor
MeshTriangle2d::~MeshTriangle2d()
{
	if(edges[0]){ // i.e. no delete-all
		// If an edge was connected with this element only, it should be removed
		if(edges[0]->removeElementLink(this)) delete edges[0];
		if(edges[1]->removeElementLink(this)) delete edges[1];
		if(edges[2]->removeElementLink(this)) delete edges[2];
	}
}

/// clear some data for faster all-delete process
void MeshTriangle2d::preDeleteAll()
{
	edges[0] = 0; // block normal adjacency update
}

//////////////////////////////////////////////////////////////////////
// Zwraca parametr "alpha" oPIsuj¹cy "równobocznoœæ" trójk¹ta
//	Parametr ten przyjmuje wartoœæ od 0(linia) do 1 (trójk¹t równoboczny)
double MeshTriangle2d::getAlphaQuality(Metric2dContext& mc, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint2d pt1 = getPoint(0)->getMetricCoordinates(mc);
	const DMPoint2d pt2 = getPoint(1)->getMetricCoordinates(mc);
	const DMPoint2d pt3 = getPoint(2)->getMetricCoordinates(mc);

	double denominator = pt1.distance2(pt2) + pt2.distance2(pt3) + pt3.distance2(pt1);
	return 2*SQRT3 * (pt2-pt1).crossProduct(pt3-pt1) / denominator;
}

double MeshTriangle2d::getAngle(Metric2dContext& mc, int i, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint2d pt1 = getPoint(i)->getMetricCoordinates(mc);
	const DMPoint2d pt2 = getPoint((i+1)%3)->getMetricCoordinates(mc);
	const DMPoint2d pt3 = getPoint((i+2)%3)->getMetricCoordinates(mc);

	return acos( (pt2-pt1).scalarProductNormalized(pt3-pt1));
}

//////////////////////////////////////////////////////////////////////
//
bool MeshTriangle2d::isPointInOuterCircle(Metric2dContext& mc, const DPoint2d& point, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint2d pt = mc.transformPStoMS(point);
	const DMPoint2d pt1 = getPoint(0)->getMetricCoordinates(mc);
	const DMPoint2d pt2 = getPoint(1)->getMetricCoordinates(mc);
	const DMPoint2d pt3 = getPoint(2)->getMetricCoordinates(mc);

	return DTriangle2d::inSphereCheck(pt, pt1, pt2, pt3) >= 0.0;
}

//////////////////////////////////////////////////////////////////////
// Zwraca d³ugoœæ promienia ko³a wPIsanego
double MeshTriangle2d::getInnerCircleRadius(Metric2dContext& mc, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint2d p1 = getPoint(0)->getMetricCoordinates(mc);
	const DMPoint2d p2 = getPoint(1)->getMetricCoordinates(mc);
	const DMPoint2d p3 = getPoint(2)->getMetricCoordinates(mc);

	double A = p1.distance(p2);
	double B = p2.distance(p3);
	double C = p3.distance(p1);
	double G = A + B + C;
	double D = p1.x*(p2.y - p3.y) - p2.x*(p1.y - p3.y) + p3.x*(p1.y - p2.y);

	return D / G;
}

double MeshTriangle2d::countMetricDiff(Metric2dContext& mc) const
{
	ControlDataMatrix2d cdmm_simplex = ControlDataMatrix2d::countMetricTensor(
			points[0]->getCoordinates(), points[1]->getCoordinates(), points[2]->getCoordinates());
	ControlDataMatrix2d cdmm_space;
	if(mc.getControlSpace()->getBaseSurface()){
		ControlDataMatrix2d p_cdm;
		ControlDataMatrix2d s_cdm = mc.getControlSpace()->getMetricAndParameterizationAtPoint(
				getMiddlePoint(), p_cdm);
		const DMatrix2d sp = s_cdm * p_cdm.inverse();
		const DMatrix2d sp2 = sp * sp;
		cdmm_space = ControlDataMatrix2d(sp2.m[0][0], sp2.m[0][1], sp2.m[1][1]);
	}else{
		cdmm_space = mc.getControlSpace()->getMetricAtPoint(getMiddlePoint()).transformationToTensor();
	}
	return cdmm_simplex.countDifferenceRR(cdmm_space);
}


double MeshTriangle2d::countMetricDiffQuality(Metric2dContext& mc)
{
	return quality = (1.0 - 0.25 * countMetricDiff(mc));
}

//////////////////////////////////////////////////////////////////////
// Counts quality of triangle according to given quality coefficient
double MeshTriangle2d::countQuality(Metric2dContext& mc, bool mind_area, int criterion)
{
	mc.countMetricAtPoints(getPoint(0), getPoint(1), getPoint(2));

	const DMPoint2d p1 = getPoint(0)->getMetricCoordinates(mc);
	const DMPoint2d p2 = getPoint(1)->getMetricCoordinates(mc);
	const DMPoint2d p3 = getPoint(2)->getMetricCoordinates(mc);

	double area = MeshElement::last_area = DTriangle2d::area(p1, p2, p3);

	int quality_criterion = criterion;
	if(criterion == -1) quality_criterion = param_quality_criterion;
	if(mind_area && area < 1.2*getMinArea()) 
		return quality = mesh_data.relative_infinity;

	switch(quality_criterion){
	case MeshData::QUALITY_RADIA:
		return quality = 2.0 * DTriangle2d::innerCircleRadius(p1, p2, p3) / 
			DTriangle2d::outerCircleRadius(p1, p2, p3);
	case MeshData::QUALITY_ALPHA:
		return quality = DTriangle2d::alphaQuality(p1, p2, p3);
	case MeshData::QUALITY_QUOTIENT_r_p:
		{
			// p - po³owa obwodu
			double p =  p1.distance(p2) + p2.distance(p3) + p3.distance(p1);
			return quality = (18.0 / SQRT3) * DTriangle2d::innerCircleRadius(p1, p2, p3) / p;
		}
	case MeshData::QUALITY_QUOTIENT_AREAS:
		{
			double r = DTriangle2d::innerCircleRadius(p1, p2, p3);
			return quality = (9.0 / SQRT3) * r*r / area;
		}
	case MeshData::QUALITY_QUOTIENT_ppp:
		{
			double a = p1.distance(p2);
			double b = p2.distance(p3);
			double c = p3.distance(p1);
			double p = 0.5*(a+b+c);	// half of circumference
			return quality = 27.0 * (p-a)*(p-b)*(p-c) / (p*p*p);
		}
	case MeshData::QUALITY_QUOTIENT_r_h:
		{
			double h = std::max(p1.distance(p2), p2.distance(p3));
			h = std::max(h, p3.distance(p1)); // longest edge of the triangle
			return quality = (6.0 / SQRT3) * DTriangle2d::innerCircleRadius(p1, p2, p3) / h;
		}
	case MeshData::QUALITY_SPACE:
		{
			double wanted_area = getOptimumArea();
			return quality = wanted_area / area;
		}
	case MeshData::QUALITY_CIRCLE_SPACE:
		{
			double wanted_radius = getOptimumCircumradius2();
			double circle_radius = DTriangle2d::outerCircleRadius2(p1, p2, p3);
			return quality = wanted_radius / circle_radius;
		}
	case MeshData::QUALITY_CIRCLE_SPACE_AND_EDGES:
		{
			double wanted_radius = getOptimumCircumradius2();
			double circle_radius = DTriangle2d::outerCircleRadius2(p1, p2, p3);
			quality = wanted_radius / circle_radius;
			if(quality > MeshGenerator2d::param_quality_threshold){
				// additionally check metric length of edges
				if( (!edges[0]->isBorder() && edges[0]->getLength(mc, true) > 1.5) ||
					(!edges[1]->isBorder() && edges[1]->getLength(mc, true) > 1.5) || 
					(!edges[2]->isBorder() && edges[2]->getLength(mc, true) > 1.5))
				{	// needs refining
					quality = 0.95 * MeshGenerator2d::param_quality_threshold;
				}
			}
			return quality;
		}
	case MeshData::QUALITY_SMALL_ANGLES:
		{
			double smallest_angle = std::min((p2-p1).getAngle(p3-p1), (p3-p2).getAngle(p1-p2));
			smallest_angle = std::min(smallest_angle, (p1-p3).getAngle(p2-p3));
			return quality = smallest_angle / (PI/3);
		}
	}

	LOG4CPLUS_WARN(MeshLog::logger_console, "Unknown quality_criterion in MeshTriangle2d::countQuality()");
	return 0.0;
}

/// Returns "mean ratio" coefficient quality
double MeshTriangle2d::getMeanRatio(Metric2dContext& mc, bool ext_metric) const
{
	int qms = MeshTriangle2d::param_quality_metric_selection;
	MeshTriangle2d::param_quality_metric_selection = 
		ext_metric ? MeshData::QM_VERTICES_AVE : MeshData::QM_MIDDLE;
	mc.countMetricAtPoints(getPoint(0), getPoint(1), getPoint(2));
	MeshTriangle2d::param_quality_metric_selection = qms;

	const DMPoint2d p1 = getPoint(0)->getMetricCoordinates(mc);
	const DMPoint2d p2 = getPoint(1)->getMetricCoordinates(mc);
	const DMPoint2d p3 = getPoint(2)->getMetricCoordinates(mc);

	double area = DTriangle2d::area(p1, p2, p3);

	return sgn(area) * 4.0 * SQRT3 * area / (p1.distance2(p2) + p1.distance2(p3) + p2.distance2(p3));
}

//////////////////////////////////////////////////////////////////////
// ZaPIsuje dane o trójk¹cie do pliku tekstowego
void MeshTriangle2d::exportToGRD(FILE *file) const
{
	int id = 0; //edges[0]->getInnerPointsCount() == 0 ? 0 : 1;
	fprintf(file, "%d %2d\t%3d %3d %3d", id, area_id, getPoint(0)->getIndex(), 
		getPoint(1)->getIndex(), getPoint(2)->getIndex());
//	for(int i = 0; i < 3; i++) edges[i]->exportInnerToGRD(file);
	fprintf(file, "\t");
	//for(int i = 0; i < 3; i++) fprintf(file, "%2d ", edges[i]->getBorderType());
	//fprintf(file, "\t");
	for(int i = 0; i < 3; i++){
		MeshElement* element = edges[i]->getOtherElement(this);
		fprintf(file, "%3d ", (element)?element->getIndex():-1);
	}
	fprintf(file, "\n");
}

//////////////////////////////////////////////////////////////////////
// Okreœla czy zadany punkt znajduje siê wewn¹trz trójk¹ta
bool MeshTriangle2d::isPointInside(const DPoint2d& pt0) const
{
	const DPoint2d& pt1 = getPoint(0)->getCoordinates();
	const DPoint2d& pt2 = getPoint(1)->getCoordinates();
	const DPoint2d& pt3 = getPoint(2)->getCoordinates();

	return  DTriangle2d::det(pt1, pt2, pt0) >= -mesh_data.relative_small_number &&
			DTriangle2d::det(pt2, pt3, pt0) >= -mesh_data.relative_small_number &&
			DTriangle2d::det(pt3, pt1, pt0) >= -mesh_data.relative_small_number;
}

//////////////////////////////////////////////////////////////////////
// Zwraca trójk¹t s¹siaduj¹cy z danym trójk¹tem w kierunku okreœlonym
//	przez zadany punkt
MeshTriangle2d* MeshTriangle2d::getNeighbourInDirection(const DPoint2d& dpt, bool mind_border) const
{
	double max = -1.0;
	int best = -1;
	for(int j = 0; j < 3; j++){
		// Wyznacz cosinus k¹ta pomiêdzy wektorem normalnym do j-tej 
		//	krawêdzi a wektorem kierunkowym
		const DPoint2d& dpt0 = getPoint(j)->getCoordinates();
		const DPoint2d& dpt1 = getPoint((j+1)%3)->getCoordinates();
		const DPoint2d mpt = DPoint2d::average(dpt0, dpt1); // Middle of edge
		DVector2d dv = dpt - mpt; // direction vector
		double len = dv.length();
		if(len > mesh_data.relative_small_number){
			// Count normalized vectors
			dv /= len;
			DVector2d nv(dpt0.y - dpt1.y, dpt1.x - dpt0.x); // Normal to edge
			const DPoint2d& dpt2 = getPoint((j+2)%3)->getCoordinates();
			double d = (DTriangle2d::det(dpt0, dpt1, dpt2) * 
				DTriangle2d::det(dpt0, dpt1, mpt + nv) > 0) ? -1.0 : 1.0;
			len = nv.length();
			if(len > mesh_data.relative_small_number){
				nv /= (d*len);
				double cos_alpha = dv.scalarProduct(nv); // for cosinus of angle 
				if(cos_alpha > max){
					best = j;
					max = cos_alpha;
				}
			}else continue;
		}else{
			best = j;
			break;
		}
	}

	if (best < 0) return nullptr;
	if(mind_border && edges[best]->isBorder()) return nullptr;
	return (MeshTriangle2d*)(edges[best]->getOtherElement(this));
}

//////////////////////////////////////////////////////////////////////
// Zamienia przek¹tn¹ w czworok¹cie zbudowanym z dwóch podanych 
//		(s¹siednich) trójk¹tów 
//	(jeœli ustawiony jest parametr improve_only - zmiana jest wykonywana
//		tylko jeœli polepszy to jakoœæ siatki)	
bool MeshTriangle2d::swapWithNeighbour(Metric2dContext& mc, int ind, bool improve_only, 
									 bool local_metric, TagExtended::TagType side_tag)
{
	if(edges[ind]->isBorder() || edges[ind]->getIntTag(TagExtended::TAG_SIDE_EDGE) > 0) return false;
	MeshTriangle2d* triangle = (MeshTriangle2d*)edges[ind]->getOtherElement(this);
	if(!triangle || triangle->getEdgeCount() != 3) return false;

	// Punkty "tworz¹ce przek¹tn¹"
	MeshPoint2d * pt1 = getPoint(ind);
	MeshPoint2d * pt2 = getPoint((ind+1)%3);
	// Punkty przeciwne
	MeshPoint2d * pt3 = getPoint((ind+2)%3);
	int j = (triangle->getPointIndex(pt1)+1)%3;
	MeshPoint2d * pt4 = triangle->getPoint(j);	

	if(local_metric){
		const DPoint2d middle = DPoint2d::average(pt1->getCoordinates(), pt2->getCoordinates(), 
			pt3->getCoordinates(), pt4->getCoordinates());
		mc.countMetricAtPoint(middle);
	}

	double area1 = getArea(mc, false);
	double area2 = triangle->getArea(mc, false);

	if(area1 <= 0.0){
		for(int i = 0; i < 5; i++){
			MeshGenerator2d::movePointByLaplace(getPoint(0));
			MeshGenerator2d::movePointByLaplace(getPoint(1));
			MeshGenerator2d::movePointByLaplace(getPoint(2));
			area1 = getArea(mc, false);
//			if(area1 > -mesh_data.relative_small_number) break;
			if(area1 > 0.0) break;
		}
	}
	if(area2 <= 0.0){
		for(int i = 0; i < 5; i++){
			MeshGenerator2d::movePointByLaplace(triangle->getPoint(0));
			MeshGenerator2d::movePointByLaplace(triangle->getPoint(1));
			MeshGenerator2d::movePointByLaplace(triangle->getPoint(2));
			area2 = triangle->getArea(mc, false);
//			if(area2 > -mesh_data.relative_small_number) break;
			if(area2 > 0.0) break;
		}
	}

	if(area1 <= 0.0 || area2 <= 0.0) return false;

	// Check if swap is possible (if the quad is convex)
	if(DTriangle2d::det(
		pt1->getMetricCoordinates(mc), 
		pt4->getMetricCoordinates(mc), 
		pt3->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER) return false;
	if(DTriangle2d::det(
		pt2->getMetricCoordinates(mc), 
		pt3->getMetricCoordinates(mc), 
		pt4->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER) return false;

	if(improve_only){
		const DMPoint2d mdpt1 = pt1->getMetricCoordinates(mc);
		const DMPoint2d mdpt2 = pt2->getMetricCoordinates(mc);
		const DMPoint2d mdpt3 = pt3->getMetricCoordinates(mc);
		const DMPoint2d mdpt4 = pt4->getMetricCoordinates(mc);

		if(MeshGenerator2d::param_swap_criterion == MeshData::SWAP_ANGLES){
			double cos_1 = (mdpt1-mdpt3).scalarProductNormalized(mdpt2-mdpt3);
			double cos_2 = (mdpt2-mdpt4).scalarProductNormalized(mdpt1-mdpt4);
			if((cos_1+cos_2) > -METRIC_SMALL_NUMBER) return false;
		}else{
			if(!checkSwapCriterion(mdpt1, mdpt2, mdpt3, mdpt4)) return false;
		}
	}

	MeshEdge2d* edge13 = edges[(ind+2)%3];
	MeshEdge2d* edge24 = triangle->edges[j];

	if(side_tag != TagExtended::TAG_NONE){
		edges[(ind+1)%3]->setIntTag(side_tag);
		edge13->setIntTag(side_tag);
		edge24->setIntTag(side_tag);
		triangle->edges[(j+2)%3]->setIntTag(side_tag);
	}

	// Swap diagonals ...
	delete edges[ind];
	MeshEdge2d *edge34 = new MeshEdge2d(pt3, pt4, this, triangle);
	// Switch points
	points[ind] = pt4;
	triangle->points[(j+1)%3] = pt3;
	// Update edges
	edges[ind] = edge24;
	edges[(ind+2)%3] = edge34;
	triangle->edges[j] = edge34;
	triangle->edges[(j+1)%3] = edge13;
	edge13->removeElementLink(this);
	edge13->addElementLink(triangle, pt3);
	edge24->removeElementLink(triangle);
	edge24->addElementLink(this, pt4);

	return true;
}

//////////////////////////////////////////////////////////////////////
// Swaps the diagonal in the quad formed by two given (adjacent) triangles
bool MeshTriangle2d::swapWithNeighbourNoCheck(int i)
{
	assert(!edges[i]->isBorder() && edges[i]->zeroIntTag(TagExtended::TAG_SIDE_EDGE));
	MeshElement* element = edges[i]->getOtherElement(this);
	assert(element && element->getEdgeCount() == 3);
	MeshTriangle2d* triangle = (MeshTriangle2d*)element;

	// Punkty "tworz¹ce przek¹tn¹"
	MeshPoint2d * pt1 = getPoint(i);
	//MeshPoint2d * pt2 = points[(i+1)%3];
	// Punkty przeciwne
	MeshPoint2d * pt3 = getPoint((i+2)%3);
	int j = (triangle->getPointIndex(pt1)+1)%3;
	MeshPoint2d * pt4 = triangle->getPoint(j);	
	
	// Zamieñ przek¹tne ...
	delete edges[i];
	MeshEdge2d *edge34 = new MeshEdge2d(pt3, pt4, this, triangle);
	// Zamieñ punkty
	points[i] = pt4;
	triangle->points[(j+1)%3] = pt3;
	// Uaktualnij krawêdzie
	MeshEdge2d* edge13 = edges[(i+2)%3];
	MeshEdge2d* edge24 = triangle->edges[j];
	edges[i] = edge24;
	edges[(i+2)%3] = edge34;
	triangle->edges[j] = edge34;
	triangle->edges[(j+1)%3] = edge13;
	edge13->removeElementLink(this);
	edge13->addElementLink(triangle, pt3);
	edge24->removeElementLink(triangle);
	edge24->addElementLink(this, pt4);

	assert(this->isOK());
	assert(triangle->isOK());

	return true;
}

/////////////////////////////////////////////////////////////////////////////
// Sprawdza czy zamiana krawêdzi _i_ jest mo¿liwa i korzystna
bool MeshTriangle2d::shouldSwap(Metric2dContext& mc, int i, bool local_metric) const
{
	if(edges[i]->isBorder()) return false;
	MeshElement* element = edges[i]->getOtherElement(this);
	if(!element || element->getEdgeCount() != 3) return false;
	MeshTriangle2d* triangle = (MeshTriangle2d*)element;
	// Punkty "tworz¹ce przek¹tn¹"
	MeshPoint2d * pt1 = getPoint(i);
	MeshPoint2d * pt2 = getPoint((i+1)%3);
	// Punkty przeciwne
	MeshPoint2d * pt3 = getPoint((i+2)%3);
	MeshPoint2d * pt4 = triangle->getOtherPoint(pt1, pt2);

	if(local_metric)
		mc.countMetricAtPoint(
			DPoint2d::average(pt1->getCoordinates(), pt2->getCoordinates(),
								pt3->getCoordinates(), pt4->getCoordinates()));

	// Check if swap is possible (if the quad is convex)
	if(DTriangle2d::det(
		pt1->getMetricCoordinates(mc), 
		pt4->getMetricCoordinates(mc), 
		pt3->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER) return false;
	if(DTriangle2d::det(
		pt2->getMetricCoordinates(mc), 
		pt3->getMetricCoordinates(mc), 
		pt4->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER) return false;

	const DMPoint2d dpt1 = pt1->getMetricCoordinates(mc);
	const DMPoint2d dpt2 = pt2->getMetricCoordinates(mc);
	const DMPoint2d dpt3 = pt3->getMetricCoordinates(mc);
	const DMPoint2d dpt4 = pt4->getMetricCoordinates(mc);

	if(MeshGenerator2d::param_swap_criterion == MeshData::SWAP_ANGLES){
		// Czy kryterium k¹tów jest naruszone ?
		double cos_1 = (dpt1-dpt3).scalarProductNormalized(dpt2-dpt3);
		double cos_2 = (dpt2-dpt4).scalarProductNormalized(dpt1-dpt4);
		return ((cos_1+cos_2) < -METRIC_SMALL_NUMBER);
	}else{
		return checkSwapCriterion(dpt1, dpt2, dpt3, dpt4);
	}
}

bool MeshTriangle2d::swapWithNeighbour(Metric2dContext& mc, const MeshEdge2d *edge, bool improve_only, bool local_metric)
{
	if(edges[0] == edge) 
		return swapWithNeighbour(mc, 0, improve_only, local_metric);
	else if(edges[1] == edge) 
		return swapWithNeighbour(mc, 1, improve_only, local_metric);
	else if(edges[2] == edge) 
		return swapWithNeighbour(mc, 2, improve_only, local_metric);
	else
		return false;
}

// Zwraca prawdê, jeœli jest spe³niony warunek zamiany (tzn. trzeba zamieniæ krawêdŸ)
bool MeshTriangle2d::checkSwapCriterion(const DMPoint2d &dpt1, const DMPoint2d &dpt2, 
									  const DMPoint2d &dpt3, const DMPoint2d &dpt4) const
{
	switch(MeshGenerator2d::param_swap_criterion){
	case MeshData::SWAP_ALPHA:
		{
			double sum1 = DTriangle2d::alphaQuality(dpt1, dpt2, dpt3) + DTriangle2d::alphaQuality(dpt2, dpt1, dpt4);
			double sum2 = DTriangle2d::alphaQuality(dpt1, dpt4, dpt3) + DTriangle2d::alphaQuality(dpt2, dpt3, dpt4);
			return (sum2 - sum1 > METRIC_SMALL_NUMBER);
		}
	case MeshData::SWAP_RR:
		{
			DPoint2d opt;
			double sum1 = DTriangle2d::innerCircleRadius(dpt1, dpt2, dpt3) / 
					DTriangle2d::outerCircleRadius(dpt1, dpt2, dpt3) +
				DTriangle2d::innerCircleRadius(dpt2, dpt1, dpt4) / 
					DTriangle2d::outerCircleRadius(dpt2, dpt1, dpt4);
			double sum2 = DTriangle2d::innerCircleRadius(dpt1, dpt4, dpt3) / 
					DTriangle2d::outerCircleRadius(dpt1, dpt4, dpt3) +
				DTriangle2d::innerCircleRadius(dpt2, dpt3, dpt4) / 
					DTriangle2d::outerCircleRadius(dpt2, dpt3, dpt4);
			return ((sum2 - sum1) > METRIC_SMALL_NUMBER);
		}
	case MeshData::SWAP_RP:
		{
			double p12 = dpt1.distance(dpt2);
			double p13 = dpt1.distance(dpt3);
			double p14 = dpt1.distance(dpt4);
			double p23 = dpt2.distance(dpt3);
			double p24 = dpt2.distance(dpt4);
			double p34 = dpt3.distance(dpt4);
			double sum1 = DTriangle2d::innerCircleRadius(dpt1, dpt2, dpt3) / (p12+p23+p13) +
				DTriangle2d::innerCircleRadius(dpt2, dpt1, dpt4) / (p12+p14+p24);
			double sum2 = DTriangle2d::innerCircleRadius(dpt1, dpt4, dpt3) / (p14+p34+p13) +
				DTriangle2d::innerCircleRadius(dpt2, dpt3, dpt4) / (p23+p34+p24);
			return (sum2 - sum1 > METRIC_SMALL_NUMBER);
		}
	case MeshData::SWAP_SCST:
		{
			double r123 = DTriangle2d::innerCircleRadius(dpt1, dpt2, dpt3);
			double r214 = DTriangle2d::innerCircleRadius(dpt2, dpt1, dpt4);
			double r143 = DTriangle2d::innerCircleRadius(dpt1, dpt4, dpt3);
			double r234 = DTriangle2d::innerCircleRadius(dpt2, dpt3, dpt4);

			double sum1 =  r123 * r123 / DTriangle2d::area(dpt1, dpt2, dpt3) +
				r214 * r214 / DTriangle2d::area(dpt2, dpt1, dpt4);
			double sum2 = r143 * r143 / DTriangle2d::area(dpt1, dpt4, dpt3) +
				r234 * r234 / DTriangle2d::area(dpt2, dpt3, dpt4);
			return (sum2 - sum1 > METRIC_SMALL_NUMBER);
		}
	case MeshData::SWAP_PPP:
		{
			double p12 = dpt1.distance(dpt2);
			double p13 = dpt1.distance(dpt3);
			double p14 = dpt1.distance(dpt4);
			double p23 = dpt2.distance(dpt3);
			double p24 = dpt2.distance(dpt4);
			double p34 = dpt3.distance(dpt4);

			double p123 = 0.5 * (p12 + p23 + p13);
			double p214 = 0.5 * (p12 + p14 + p24);
			double p143 = 0.5 * (p14 + p34 + p13);
			double p234 = 0.5 * (p23 + p34 + p24);

			double q1_before = (p123-p12)*(p123-p23)*(p123-p13) / (p123*p123*p123);
			double q2_before = (p214-p12)*(p214-p14)*(p214-p24) / (p214*p214*p214);
			double q1_after = (p143-p14)*(p143-p34)*(p143-p13) / (p143*p143*p143);
			double q2_after = (p234-p23)*(p234-p34)*(p234-p24) / (p234*p234*p234);
			return (std::min(q1_after, q2_after) - std::min(q1_before, q2_before) > METRIC_SMALL_NUMBER);
		}
	case MeshData::SWAP_RH:
		{
			double p12 = dpt1.distance(dpt2);
			double p13 = dpt1.distance(dpt3);
			double p14 = dpt1.distance(dpt4);
			double p23 = dpt2.distance(dpt3);
			double p24 = dpt2.distance(dpt4);
			double p34 = dpt3.distance(dpt4);

			double r123 = DTriangle2d::innerCircleRadius(dpt1, dpt2, dpt3);
			double r214 = DTriangle2d::innerCircleRadius(dpt2, dpt1, dpt4);
			double r143 = DTriangle2d::innerCircleRadius(dpt1, dpt4, dpt3);
			double r234 = DTriangle2d::innerCircleRadius(dpt2, dpt3, dpt4);

			double sum1 = r123 / std::max(p12, std::max(p23, p13)) + 
				r214 / std::max(p12, std::max(p14, p24));
			double sum2 = r143 / std::max(p14, std::max(p34, p13)) + 
				r234 / std::max(p23, std::max(p34, p24));
			return (sum2 - sum1 > METRIC_SMALL_NUMBER);
		}
	}

	return false;
}

MeshTriangle2d* MeshTriangle2d::findTriangleByNeighbours(const DPoint2d& mpt, bool mind_border)
{
	MeshTriangle2d* triangle = this;
	MeshTriangle2d* prev_triangle = nullptr;
	unsigned int counter = 0;

	while(true){
		// Czy punkt jest wewn¹trz trójk¹ta?
		if(!triangle){
//			LOG4CPLUS_INFO(MeshLog::logger_console, "MeshTriangle2d: error findTriangleByNeighbours - null triangle encountered.");
			return nullptr;
		}
		if(++counter > 10000){
			LOG4CPLUS_INFO(MeshLog::logger_console, "MeshTriangle2d: error findTriangleByNeighbours - 10000 steps exceeded.");
			return nullptr;
		}
		if(triangle->isPointInside(mpt)) break;
		// Jeœli to jeszcze nie ten trójk¹t ...
		// Wybór jednego z trójk¹tów s¹siednich - w kierunku szukanego punktu
		MeshTriangle2d* t = triangle->getNeighbourInDirection(mpt, mind_border);
		if(!t) return nullptr;

		if(prev_triangle == t){
			LOG4CPLUS_INFO(MeshLog::logger_console, "MeshTriangle2d: error findTriangleByNeighbours - loop encountered.");
			return nullptr;
		}
		prev_triangle = triangle;
		triangle = t;
	}

	return triangle;
}

double MeshTriangle2d::getArea(Metric2dContext& mc, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	return DTriangle2d::area(
		getPoint(0)->getMetricCoordinates(mc), 
		getPoint(1)->getMetricCoordinates(mc), 
		getPoint(2)->getMetricCoordinates(mc));
}

double MeshTriangle2d::getAreaNoMetric() const
{
	return DTriangle2d::area(
		getPoint(0)->getCoordinates(), 
		getPoint(1)->getCoordinates(), 
		getPoint(2)->getCoordinates());
}

DMPoint2d MeshTriangle2d::getOuterCircleCenter(Metric2dContext& mc, bool local_metric)
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	return DTriangle2d::outerCircleCenter(
		getPoint(0)->getMetricCoordinates(mc), 
		getPoint(1)->getMetricCoordinates(mc), 
		getPoint(2)->getMetricCoordinates(mc));
}

bool MeshTriangle2d::isInverted() const
{
	return DTriangle2d::det(
		getPoint(0)->getCoordinates(), 
		getPoint(1)->getCoordinates(), 
		getPoint(2)->getCoordinates()) <= 0.0;	
}

MeshPoint2d* MeshTriangle2d::getOtherPoint(const MeshPoint2d* p0, const MeshPoint2d* p1) const
{
	if(getPoint(0) != p0 && getPoint(0) != p1) return getPoint(0);
	if(getPoint(1) != p0 && getPoint(1) != p1) return getPoint(1);
	assert(getPoint(2) != p0 && getPoint(2) != p1);
	return getPoint(2);
}

/// add this element as a lump of mass (respecting control space)
void MeshTriangle2d::addForInertialCenter(Metric2dContext& mc, DPoint2d& inertial_center, double& total_mass) const
{
	addForInertialCenterAdaptiveSplitLongestEdge(mc, 
		getPoint(0)->getCoordinates(), 
		getPoint(1)->getCoordinates(),
		getPoint(2)->getCoordinates(),
		inertial_center, total_mass);
}

/// add this element as a lump of mass (respecting control space)
void MeshTriangle2d::addForInertialCenterAdaptiveSplitLongestEdge(Metric2dContext& mc, 
	const DPoint2d& pt0, const DPoint2d& pt1, const DPoint2d& pt2, 
	DPoint2d& inertial_center, double& total_mass, int level)
{
	const DPoint2d middle = DPoint2d::average(pt0, pt1, pt2);
	mc.countMetricAtPoint(middle);
	const DMPoint2d mpt0 = mc.transformPStoMS(pt0);
	const DMPoint2d mpt1 = mc.transformPStoMS(pt1);
	const DMPoint2d mpt2 = mc.transformPStoMS(pt2);
	double len01 = mpt0.distance2(mpt1);
	double len02 = mpt0.distance2(mpt2);
	double len12 = mpt1.distance2(mpt2);
	if(level > 10 || std::max(len01, std::max(len02, len12)) < 4.0){ // len^2 !!!
		double area = DTriangle2d::area(mpt0, mpt1, mpt2);
		total_mass += area;
		inertial_center.add(middle, area);
	}else{
		if(len01 > len02 && len01 > len12){
			const DPoint2d pt = DPoint2d::average(pt0, pt1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt0, pt, pt2, 
				inertial_center, total_mass, level+1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt, pt1, pt2, 
				inertial_center, total_mass, level+1);
		}else if (len02 > len12){
			const DPoint2d pt = DPoint2d::average(pt0, pt2);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt0, pt1, pt, 
				inertial_center, total_mass, level+1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt1, pt2, pt, 
				inertial_center, total_mass, level+1);
		}else{
			const DPoint2d pt = DPoint2d::average(pt1, pt2);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt0, pt1, pt, 
				inertial_center, total_mass, level+1);
			addForInertialCenterAdaptiveSplitLongestEdge(mc, pt0, pt, pt2, 
				inertial_center, total_mass, level+1);
		}
	}
}

/// add this element as a lump of mass (respecting control space)
void MeshTriangle2d::addForInertialMoments(Metric2dContext& mc, const DPoint2d& inertial_center, DMatrix2d& inertial_moments) const
{
	addForInertialMomentsAdaptiveSplitLongestEdge(mc, 
		getPoint(0)->getCoordinates(), 
		getPoint(1)->getCoordinates(),
		getPoint(2)->getCoordinates(),
		inertial_center, inertial_moments);
}

/// add this element as a lump of mass (respecting control space)
void MeshTriangle2d::addForInertialMomentsAdaptiveSplitLongestEdge(Metric2dContext& mc, 
	const DPoint2d& pt0, const DPoint2d& pt1, const DPoint2d& pt2, 
	const DPoint2d& inertial_center, DMatrix2d& inertial_moments, int level)
{
	const DPoint2d middle = DPoint2d::average(pt0, pt1, pt2);
	mc.countMetricAtPoint(middle);
	const DMPoint2d mpt0 = mc.transformPStoMS(pt0);
	const DMPoint2d mpt1 = mc.transformPStoMS(pt1);
	const DMPoint2d mpt2 = mc.transformPStoMS(pt2);
	double len01 = mpt0.distance2(mpt1);
	double len02 = mpt0.distance2(mpt2);
	double len12 = mpt1.distance2(mpt2);
	if(level > 10 || std::max(len01, std::max(len02, len12)) < 4.0){ // len^2 < 4.0 (len < 2.0)
		double area = DTriangle2d::area(mpt0, mpt1, mpt2);
		const DVector2d dv = middle - inertial_center;
		inertial_moments.m[0][0] +=  area * (dv.y*dv.y); // Ixx
		inertial_moments.m[1][1] +=  area * (dv.x*dv.x); // Iyy
		inertial_moments.m[0][1] += -area * (dv.x*dv.y); // Ixy
		inertial_moments.m[1][0] += -area * (dv.x*dv.y); // Iyx = Ixy
	}else{
		if(len01 > len02 && len01 > len12){
			const DPoint2d pt = DPoint2d::average(pt0, pt1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt0, pt, pt2, 
				inertial_center, inertial_moments, level+1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt, pt1, pt2, 
				inertial_center, inertial_moments, level+1);
		}else if (len02 > len12){
			const DPoint2d pt = DPoint2d::average(pt0, pt2);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt0, pt1, pt, 
				inertial_center, inertial_moments, level+1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt1, pt2, pt, 
				inertial_center, inertial_moments, level+1);
		}else{
			const DPoint2d pt = DPoint2d::average(pt1, pt2);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt0, pt1, pt, 
				inertial_center, inertial_moments, level+1);
			addForInertialMomentsAdaptiveSplitLongestEdge(mc, pt0, pt, pt2, 
				inertial_center, inertial_moments, level+1);
		}
	}
}

