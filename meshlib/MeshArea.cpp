/////////////////////////////////////////////////////////////////////////////
// MeshArea.cpp
// Element pomocniczy reprezentuj¹cy obszar geometryczny
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////

#include "DPoint.h"
#include "MeshPoint2d.h"
#include "MeshPoint3d.h"
#include "MeshEdge2d.h"
#include "MeshArea.h"
#include "Curve2dParametric.h"
#include "SurfaceParametric.h"
#include "MeshData.h"
#include "MeshViewSet.h"

//////////////////////////////////////////////////////////////////////
// Standardowy konstruktor
MeshArea::MeshArea(int ct, MeshPoint2d** new_points, bool _filled, double k) : 
	MeshElement(ct, nullptr, new_points), filled(_filled), k_material(k)
{
	if(ct > 0){
		edges = new MeshEdge2d*[ct];
		for(int i = 0; i < ct; i++){
			// The edges should already be there
			MeshPoint2d* next_point = new_points[(i+1)%count];
			edges[i] = new_points[i]->getEdgeToPoint(next_point);
			assert(edges[i]);
			if(filled){
				edges[i]->addElementLink(this, new_points[i]);
				edges[i]->addDirection(1, new_points[i]);
			}else{
				edges[i]->addElementLink(this, next_point);
				edges[i]->addDirection(1, next_point);
			}
		}
	}
}

MeshArea::MeshArea(DataVector<MeshPoint2d*> &new_points, bool _filled, double k) :
	MeshElement((int)new_points.countInt(), nullptr, nullptr), filled(_filled), k_material(k)
{
	if(count > 0){
		points = new MeshPoint2d*[count];
		for(int i = 0; i < count; i++) points[i] = new_points[i];
		edges = new MeshEdge2d*[count];
		for(int i = 0; i < count; i++){
			// The edges should already be there
			MeshPoint2d* next_point = new_points[(i+1)%count];
			edges[i] = new_points[i]->getEdgeToPoint(next_point);
			assert(edges[i]);
			if(filled){
				edges[i]->addElementLink(this, new_points[i]);
				edges[i]->addDirection(1, new_points[i]);
			}else{
				edges[i]->addElementLink(this, next_point);
				edges[i]->addDirection(1, next_point);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Destruktor
MeshArea::~MeshArea()
{
	for(int i = 0; i < count; i++){	
		// Jeœli dana krawêdŸ by³a po³¹czona tylko z tym elementem, nale¿y j¹ usun¹æ
		if(edges[i]->removeElementLink(this)) delete edges[i];
	}
	if(count > 0){
		delete[] points;
		delete[] edges;
	}
}

/// clear some data for faster all-delete process
void MeshArea::preDeleteAll()
{
	if(count > 0){
		delete[] points;
		delete[] edges;
		count = 0; // block normal adjacency update
	}
}

void MeshArea::drawGL() const
{
	// TODO ?
	assert(false);
}

DRect MeshArea::getBoundingRect() const
{
	DRect rect;
	for(int i = 0; i < count; i++)
		edges[i]->addToBoundingRect(rect);
	return rect;
}

double MeshArea::getShortestEdge(std::shared_ptr<const SurfaceParametric> surface) const
{
	double shortest = mesh_data.relative_infinity;
	for(int i = 0; i < count; i++){
		//Curve2dParametric* figure = edges[i]->getShape();
		double factor = 1.0 - MEASURE_PRECISION;
		double heap[100];
		double t0 = 0.0;
		double t1 = 1.0;
		DPoint2d pt2d = edges[i]->getPoint(t0);
		DPoint3d pt0 = surface->getPoint(pt2d);
		pt2d = edges[i]->getPoint(t1);
		DPoint3d pt1 = surface->getPoint(pt2d);
		heap[0] = t1;
		int top = 0;
		while(true){
			double t_middle = (t0+t1)*0.5;
			pt2d = edges[i]->getPoint(t_middle);
			DPoint3d pt_middle = surface->getPoint(pt2d);
			double distance = pt0.distance(pt1);
			double real_distance = pt0.distance(pt_middle) + pt_middle.distance(pt1);
			if((distance >= factor * real_distance) || (top == 99)){
				double dist = (edges[i]->getPoint(t0)).distance(edges[i]->getPoint(t1));
				if(dist > 0.0 && dist < shortest) shortest = dist;
				if(top == 0) break;
				t0 = t1;
				pt0 = pt1;
				pt1 = surface->getPoint(edges[i]->getPoint(t1 = heap[--top]));
			}else{						
				heap[++top] = t1 = t_middle;
				pt1 = pt_middle;
			}
		}
	}

	return shortest;
}

//////////////////////////////////////////////////////////////////////
// Ustala orientacjê wielok¹ta (zadanego tablic¹ wierzcho³ków)
//	false - zorientowany zgodnie z ruchem wskazówek (brak materia³u)
//	true - zorientowany przeciwnie do kierunku wskazówek (obszar wype³niony)
bool MeshArea::getOrientation(const DataVector<MeshPoint2d*> &points)
{
	int pct = (int)points.countInt();
	if(pct < 3) return false;

	double d = 0.0;

	for(int i = 0; i < pct; i++){
		const DPoint2d& p0 = points[i]->getCoordinates();
		const DPoint2d& p1 = points[(i + 1) % pct]->getCoordinates();

		d += (p1.x-p0.x)*(p1.y+p0.y);
	}
	return d < 0.0;
}

MeshViewSet* MeshArea::getViewSet(MeshViewSet* view_set, std::shared_ptr<const SurfaceParametric> surface) const
{
	if(count == 0) return view_set;

	if(view_set){
		view_set->prepareFreePlace(count, count, 0, 0);
	}else{
		view_set = new MeshViewSet(count, count, 0, 0);
	}

	// others
	for(int i = 0; i < count; i++){
		view_set->addPoint(getPoint(i), surface);
		view_set->addEdge(edges[i], surface);
	}

	return view_set;
}
