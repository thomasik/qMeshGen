// MeshQuad2d.cpp: implementation of the MeshQuad2d class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshQuad2d.h"
#include "ControlSpace2d.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "Metric2dContext.h"
#include "DTriangle.h"
#include "DQuad.h"

MeshQuad2d::MeshQuad2d(MeshPoint2d *p1, MeshPoint2d *p2, MeshPoint2d *p3, MeshPoint2d* p4) :
	MeshElement(4, quad_edges, quad_points)
{
	points[0] = p1;
	points[1] = p2;
	points[2] = p3;
	points[3] = p4;

	edges[0] = p1->getEdgeToPoint(p2);
	if(!edges[0]) edges[0] = new MeshEdge2d(p1, p2, this);
	else edges[0]->addElementLink(this, p1);

	edges[1] = p2->getEdgeToPoint(p3);
	if(!edges[1]) edges[1] = new MeshEdge2d(p2, p3, this);
	else edges[1]->addElementLink(this, p2);

	edges[2] = p3->getEdgeToPoint(p4);
	if(!edges[2]) edges[2] = new MeshEdge2d(p3, p4, this);
	else edges[2]->addElementLink(this, p3);

	edges[3] = p4->getEdgeToPoint(p1);
	if(!edges[3]) edges[3] = new MeshEdge2d(p4, p1, this);
	else edges[3]->addElementLink(this, p4);
}

MeshQuad2d::~MeshQuad2d()
{
	if(edges[0]){ // i.e. no delete-all
		// If an edge was connected with this element only, it should be removed
		if(edges[0]->removeElementLink(this)) delete edges[0];
		if(edges[1]->removeElementLink(this)) delete edges[1];
		if(edges[2]->removeElementLink(this)) delete edges[2];
		if(edges[3]->removeElementLink(this)) delete edges[3];
	}
}

/// clear some data for faster all-delete process
void MeshQuad2d::preDeleteAll()
{
	edges[0] = 0; // block normal adjacency update
}

//////////////////////////////////////////////////////////////////////
// Zwraca parametr "alpha" oPIsuj¹cy "kwadratowoœæ" czworok¹ta
//	Parametr ten przyjmuje wartoœæ od 0(linia) do 1 (kwadrat)
double MeshQuad2d::getAlphaQuality(Metric2dContext& mc, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint2d pt1 = getPoint(0)->getMetricCoordinates(mc);
	const DMPoint2d pt2 = getPoint(1)->getMetricCoordinates(mc);
	const DMPoint2d pt3 = getPoint(2)->getMetricCoordinates(mc);
	const DMPoint2d pt4 = getPoint(3)->getMetricCoordinates(mc);

	return DQuad2d::alphaQuality(pt1, pt2, pt3, pt4);
}

double MeshQuad2d::getAngle(Metric2dContext& mc, int i, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint2d pt1 = getPoint(i)->getMetricCoordinates(mc);
	const DMPoint2d pt2 = getPoint((i+1)%4)->getMetricCoordinates(mc);
	const DMPoint2d pt3 = getPoint((i+3)%4)->getMetricCoordinates(mc);

	return DMVector2d::angle(pt1, pt2, pt3);
}

//////////////////////////////////////////////////////////////////////
// Wylicza jakoœæ czworok¹ta wg. aktualnie przyjêtego kryterium ???
double MeshQuad2d::countQuality(Metric2dContext& mc, bool /* mind_area */, int /* criterion */)
{
	mc.countMetricAtPoint(getMiddlePoint());
	return DQuad2d::alphaQuality(
		getPoint(0)->getMetricCoordinates(mc),
		getPoint(1)->getMetricCoordinates(mc),
		getPoint(2)->getMetricCoordinates(mc), 
		getPoint(3)->getMetricCoordinates(mc));
}

//////////////////////////////////////////////////////////////////////
// ZaPIsuje dane o czworok¹cie do pliku tekstowego
void MeshQuad2d::exportToGRD(FILE *file) const
{
	int id = 2; //edges[0]->getInnerPointsCount() == 0 ? 2 : 3;
	fprintf(file, "%d %2d\t%3d %3d %3d %3d", id, area_id, getPoint(0)->getIndex(), 
		getPoint(1)->getIndex(), getPoint(2)->getIndex(), getPoint(3)->getIndex());
	int i;
	//for(i = 0; i < 4; i++) edges[i]->exportInnerToGRD(file);
	fprintf(file, "\t");
	//for(i = 0; i < 4; i++){
	//	fprintf(file, "%2d ", edges[i]->getBorderType());
	//}
	//fprintf(file, "\t");
	for(i = 0; i < 4; i++){
		MeshElement* element = edges[i]->getOtherElement(this);
		fprintf(file, "%3d ", (element)?element->getIndex():-1);
	}
	fprintf(file, "\n");
}

//////////////////////////////////////////////////////////////////////
// Okreœla czy zadany punkt znajduje siê wewn¹trz czworok¹ta
bool MeshQuad2d::isPointInside(const DPoint2d& pt0) const
{
	const DPoint2d& pt1 = getPoint(0)->getCoordinates();
	const DPoint2d& pt2 = getPoint(1)->getCoordinates();
	const DPoint2d& pt3 = getPoint(2)->getCoordinates();
	const DPoint2d& pt4 = getPoint(3)->getCoordinates();

	return  (DTriangle2d::det(pt1, pt2, pt0) >= 0.0 &&
		DTriangle2d::det(pt2, pt3, pt0) >= 0.0 && 
		DTriangle2d::det(pt3, pt4, pt0) >= 0.0 && 
		DTriangle2d::det(pt4, pt1, pt0) >= 0.0);
}

double MeshQuad2d::getArea(Metric2dContext& mc, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getMiddlePoint());

	const DMPoint2d pt1 = getPoint(0)->getMetricCoordinates(mc);
	const DMPoint2d pt2 = getPoint(1)->getMetricCoordinates(mc);
	const DMPoint2d pt3 = getPoint(2)->getMetricCoordinates(mc);
	const DMPoint2d pt4 = getPoint(3)->getMetricCoordinates(mc);

	double area1 = DTriangle2d::area(pt1, pt2, pt3);
	double area2 = DTriangle2d::area(pt3, pt4, pt1);

	double area = area1 + area2;
	return (area >= 0.0)?area:-area;
}

double MeshQuad2d::getAreaNoMetric() const
{
	const DPoint2d& pt1 = getPoint(0)->getCoordinates();
	const DPoint2d& pt2 = getPoint(1)->getCoordinates();
	const DPoint2d& pt3 = getPoint(2)->getCoordinates();
	const DPoint2d& pt4 = getPoint(3)->getCoordinates();

	double area1 = DTriangle2d::area(pt1, pt2, pt3);
	double area2 = DTriangle2d::area(pt3, pt4, pt1);

	double area = area1 + area2;
	return (area >= 0.0)?area:-area;
}

bool MeshQuad2d::isInverted() const
{
	const DPoint2d& pt1 = getPoint(0)->getCoordinates();
	const DPoint2d& pt2 = getPoint(1)->getCoordinates();
	const DPoint2d& pt3 = getPoint(2)->getCoordinates();
	const DPoint2d& pt4 = getPoint(3)->getCoordinates();

	double area1 = DTriangle2d::det(pt1, pt2, pt3);
	double area2 = DTriangle2d::det(pt3, pt4, pt1);

	if(area1 * area2 < 0.0){
		area1 = DTriangle2d::det(pt2, pt3, pt4);
		area2 = DTriangle2d::det(pt4, pt1, pt2);
	}

	if(area1 * area2 < 0.0) return true;
	if(area1 < 0.0 && area2 < 0.0) return true;

	return false;
}

