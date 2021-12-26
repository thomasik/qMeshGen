/////////////////////////////////////////////////////////////////////////////
// MeshElement.cpp
// Klasa abstrakcyjna - reprezentuje element siatki
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////

#include <algorithm>

#include "MeshElement.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "SurfaceParametric.h"
#include "MeshViewSet.h"
#include "MeshPoint3d.h"
#include "SurfacePlane.h"
#include "DTriangle.h"

int MeshElement::param_point_distance = MeshData::RP_VERTICES;

double MeshElement::last_area = 0.0;

//////////////////////////////////////////////////////////////////////
// Standardowy konstruktor
MeshElement::MeshElement(int _count, MeshEdge2d**_edges, MeshPoint2d** _points)
	: quality(-1.0), area_id(-1), count(_count), edges(_edges), points(_points)
{ }

//////////////////////////////////////////////////////////////////////
// Funkcja zwraca kolejn¹ krawêdŸ wzglêdem podanej (wszystkie elementy
//	skierowane s¹ w lewo - jeœli zawieraj¹ materia³)
MeshEdge2d* MeshElement::getNextEdge(const MeshEdge2d *edge) const
{
	for(int i = 0; i < count; i++){
		if(edges[i] == edge) return edges[(i+1) % count];
	}
	assert(false);
	return nullptr;
}

MeshEdge2d* MeshElement::getOtherEdge(const MeshEdge2d* edge, const MeshPoint2d* point) const
{
	for(int i = 0; i < count; i++)
		if(edges[i] != edge && edges[i]->incidentTo(point)) 
			return edges[i];
	assert(false);
	return nullptr;
}

//////////////////////////////////////////////////////////////////////
// Porównuje jakoœæ dwóch elementów
short MeshElement::compareTo(MeshElement *item)
{
	double q1 = getQuality();
	double q2 = item->getQuality();
	if(q1 > q2){
		return 1;
	}else if(q1 < q2){
		return -1;
	}else return 0;
}

//////////////////////////////////////////////////////////////////////
// Funkcja zwraca poprzedni¹ krawêdŸ wzglêdem podanej
MeshEdge2d* MeshElement::getPrevEdge(const MeshEdge2d *edge) const
{
	for(int i = 0; i < count; i++){
		if(edges[i] == edge) return edges[(i + count - 1) % count];
	}
	assert(false);
	return nullptr;
}

//////////////////////////////////////////////////////////////////////
// Funkcja zamienia podan¹ krawêdŸ na now¹
void MeshElement::switchEdges(const MeshEdge2d *edge1, MeshEdge2d *edge2)
{
	for(int i = 0; i < count; i++){
		if(edges[i] == edge1){
			edges[i] = edge2;
			return;
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Funkcja zamienia podany stary punkt na nowy
void MeshElement::switchPoints(const MeshPoint2d *point1, MeshPoint2d *point2)
{
	for(int i = 0; i < count; i++){
		if(points[i] == point1){
			points[i] = point2;
			return;
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Funkcja zwraca indeks zadanego punktu w tablicy wierzcho³ków
int MeshElement::getPointIndex(const MeshPoint2d *point) const
{
	for(int i = 0; i < count; i++)
		if(getPoint(i) == point) return i;

	return -1;
}

//////////////////////////////////////////////////////////////////////
// Funkcja zamienia stary punkt na nowy, z uaktualnieniem s¹siednich 
//	krawêdzi
void MeshElement::switchPointsWithEdges(const MeshPoint2d *point1, MeshPoint2d *point2)
{
	for(int i = 0; i < count; i++){
		if(getPoint(i) == point1){
			points[i] = point2;
			int j = (i + count - 1) % count;
			// KrawêdŸ przed punktem
			if(edges[j]->removeElementLink(this)) delete edges[j];
			MeshPoint2d* mp_j = getPoint(j);
			edges[j] = mp_j->getEdgeToPoint(point2);
			if(!edges[j]){
				// Nowa krawêdŸ
				edges[j] = new MeshEdge2d(mp_j, point2, this);
			}else{
				edges[j]->addElementLink(this, mp_j);
			}
			// KrawêdŸ po punkcie
			if(edges[i]->removeElementLink(this)) delete edges[i];
			j = (i + 1) % count;
			mp_j = getPoint(j);
			edges[i] = mp_j->getEdgeToPoint(point2);
			if(!edges[i]){
				// Nowa krawêdŸ
				edges[i] = new MeshEdge2d(point2, mp_j, this);
			}else{
				edges[i]->addElementLink(this, point2);
			}

//			assert(!isInverted());

			return;
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Funkcja zwraca kolejny punkt wzglêdem podanego
MeshPoint2d* MeshElement::getNextPoint(const MeshPoint2d *point) const
{
	for(int i = 0; i < count; i++){
		if(getPoint(i) == point) return getPoint((i+1) % count);
	}
	return nullptr;
}

//////////////////////////////////////////////////////////////////////
// Funkcja zwraca poprzedni punkt wzglêdem podanego
MeshPoint2d* MeshElement::getPrevPoint(const MeshPoint2d *point) const
{
	for(int i = 0; i < count; i++){
		if(getPoint(i) == point) return getPoint((i + count - 1) % count);
	}
	return nullptr;
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca œrodek trójk¹ta (jako œrodek masy wierzcho³ków)
DPoint2d MeshElement::getMiddlePoint() const
{
	DPoint2d middle = getPoint(0)->getCoordinates();
	for(int i = 1; i < count; i++){
		middle.add(getPoint(i)->getCoordinates());
	}
	middle /= (double)count;
	return middle;
}

//////////////////////////////////////////////////////////////////////
// Zwraca maksymaln¹ ró¿nicê wartoœci identyfikatorów wierzcho³ków
int MeshElement::getIdSpan() const
{
	if(count < 1) return 0;
	int min_id = getPoint(0)->getIntTag(TagExtended::TAG_ID);
	int max_id = min_id;
	for(int i = 0; i < count; i++){
		// Wierzcho³ek
		int id = getPoint(i)->getIntTag(TagExtended::TAG_ID);
		if(id < min_id) min_id = id;
		if(id > max_id) max_id = id;
/*
		// Punkty wewnêtrzne krawêdzi
		int ict = edges[i]->getInnerPointsCount();
		for(int j = 0; j < ict; j++){
			id = edges[i]->getInnerMeshPoint(j)->getID();
			if(id < min_id) min_id = id;
			if(id > max_id) max_id = id;
		}
*/
	}
	return max_id - min_id + 1;
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca wierzcho³ek (w elemencie) o indeksie wiêkszym ni¿ k, przy
//	czym zwracany jest ten z wierzcho³ków, który ma najni¿sz¹ krotnoœæ
MeshPoint2d* MeshElement::getMinVertexBiggerThen(int k) const
{
	MeshPoint2d *point = nullptr;
	for(int i = 0; i < count; i++){
		MeshPoint2d* mp = getPoint(i);
		if(mp->getIndex() >= k && ((point == nullptr) || (mp->getRank() < point->getRank()))){
			point = mp;
		}
	}

	return point;
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca wierzcho³ek le¿¹cy na krawêdzi przylegaj¹cej do zadanego punktu, 
//	o indeksie wiêkszym ni¿ k, (o ile taki istnieje)
/*
MeshPoint2d* MeshElement::getInnerVertexBiggerThen(const MeshPoint2d* point, int k) const
{
	for(int i = 0; i < count; i++){	// Dla ka¿dej krawêdzi
		if(edges[i]->getPointIndex(point) > -1){	// Jeœli zawiera ten punkt 
			int ict = edges[i]->getInnerPointsCount();
			for(int j = 0; j < ict; j++){
				MeshPoint2d* pt = edges[i]->getInnerMeshPoint(j);
				if(pt->getIndex() >= k) return pt;
			}
		}
	}

	return nullptr;
}
*/

std::shared_ptr<MeshViewFaceData> MeshElement::getViewData(double shrink_ratio,
	SurfaceConstPtr surface, bool proper_orientation) const
{
	auto data = std::make_shared<MeshViewFaceData>(count);
	data->area_id = area_id;
	data->quality = (float)quality;
	FPoint3d middle = FPoint3d::zero;

	int shift = 0;
	if(count == 4){
		// check for convex quads
		if((DTriangle2d::det(getPoint(0)->getCoordinates(), 
				getPoint(1)->getCoordinates(), getPoint(2)->getCoordinates()) < 0.0) ||
			DTriangle2d::det(getPoint(2)->getCoordinates(), 
				getPoint(3)->getCoordinates(), getPoint(0)->getCoordinates()) < 0.0)
			shift = 1;
	}

	float f = 1.0f / count;
	for(int i = 0; i < count; i++){
		const MeshPoint2d* mpt = getPoint((i+shift)%count);
		if(!surface){
			const DPoint2d& ptd = mpt->getCoordinates();
			data->pts.add( FPoint3d((float)ptd.x, (float)ptd.y) );
		}else{
			MeshPoint3d* mpt3 = (MeshPoint3d*)mpt->getPtrTag(TagExtended::TAG_MP_2D_3D);
			if(mpt3) data->pts.add( mpt3->getCoordinates() );
			else data->pts.add( surface->getPoint(mpt->getCoordinates()) );
		}
		if(mpt->availableTag(TagExtended::TAG_VISUALIZATION)){
			data->wpts.add( (float) mpt->getDoubleTag(TagExtended::TAG_VISUALIZATION) );
		}else{
			data->wpts.add( 1.0f );
		}
		middle.add(data->pts[i], f);
	}

	for(int i = 2; i < count; i++){
		data->indices.add( 0 );
		data->indices.add( i-1 );
		data->indices.add( i );
	}

	if(!proper_orientation){
		std::reverse( data->pts.begin(), data->pts.end() );
		std::reverse( data->wpts.begin(), data->wpts.end() );
		std::reverse( data->indices.begin(), data->indices.end() );
	}

	if(shrink_ratio < 1.0f){
		for(int i = 0; i < count; i++)
			data->pts[i]  = middle + (data->pts[i] - middle) * (float)shrink_ratio;
	}
	data->countNormal();
	return data;
}

bool MeshElement::properOrientation(const MeshPoint2d *mp1, const MeshPoint2d *mp2) const
{	
	for(int i = 0; i < count; i++){
		if(getPoint(i) == mp1)
			return (getPoint((i+1)%count) == mp2);
	}
	assert(false);
	return false;
}

MeshPoint2d* MeshElement::closeToVertex(Metric2dContext& mc, MeshPoint2d* point) const
{
	const DMPoint2d dpt = point->getMetricCoordinates(mc);
	for(int i = 0; i < count; i++){
		MeshPoint2d* mpt = getPoint(i);
		if(dpt.distance2(mpt->getMetricCoordinates(mc)) < METRIC_SMALL_NUMBER)
			return mpt;
	}
	return nullptr;
}

double MeshElement::distance2(const DPoint2d& pt) const
{
	double min_dist2 = mesh_data.relative_infinity;
	if((param_point_distance & MeshData::RP_MIDDLE) != 0){
		double dist2 = pt.distance2(getMiddlePoint());
		if(dist2 < min_dist2) min_dist2 = dist2;
	}
	if((param_point_distance & MeshData::RP_VERTICES) != 0)
		for(int i = 0; i < count; i++){
			double dist2 = pt.distance2(getPoint(i)->getCoordinates());
			if(dist2 < min_dist2) min_dist2 = dist2;
		}
	if((param_point_distance & MeshData::RP_MIDEDGES) != 0)
		for(int i = 0; i < count; i++){
			double dist2 = pt.distance2(getEdge(i)->getPoint(0.5));
			if(dist2 < min_dist2) min_dist2 = dist2;
		}
	return min_dist2;
}

void MeshElement::getRepresentativePoints(DataVector<DPoint2d> & pts) const
{
	if((param_point_distance & MeshData::RP_MIDDLE) != 0)
		pts.add(getMiddlePoint());
	if((param_point_distance & MeshData::RP_VERTICES) != 0)
		for(int i = 0; i < count; i++)
			pts.add(getPoint(i)->getCoordinates());
	if((param_point_distance & MeshData::RP_MIDEDGES) != 0)
		for(int i = 0; i < count; i++)
			pts.add(getEdge(i)->getPoint(0.5));
}

/// Get real angles (projected onto plane)
DataVector<double> MeshElement::getRealAngles(SurfaceConstPtr surface) const
{
	DataVector<DPoint3d> pts3d(count);
	for(int i = 0; i < count; i++)
		pts3d.add(surface->getPoint(points[i]->getCoordinates()));

	const DVector3d e0 = (pts3d[1]-pts3d[0]).normalized();
	const DVector3d e1 = e0.crossProduct(pts3d[2]-pts3d[0]).crossProduct(e0).normalized();
	SurfacePlane plane(pts3d[0], e0, e1);

	DataVector<DPoint2d> pts2d(count);
	for(int i = 0; i < count; i++)
		pts2d.add(plane.getParameters(pts3d[i]));

	DataVector<double> angles(count);
	for(int i = 0; i < count; i++)
		angles.add((pts2d[(i+1)%count] - pts2d[i]).getAngle(pts2d[(i+count-1)%count] - pts2d[i]));

	return angles;
}

bool MeshElement::isPointInside(MeshPoint2d *point) const 
{ 
	return isPointInside(point->getCoordinates()); 
}

/// add this element as a lump of mass (respecting control space)
void MeshElement::addForInertialCenter(Metric2dContext& mc, DPoint2d& inertial_center, double& total_mass) const
{
	double area = getArea(mc);
	total_mass += area;
	inertial_center.add(getMiddlePoint(), area);
}

/// add this element as a lump of mass (respecting control space)
void MeshElement::addForInertialMoments(Metric2dContext& mc, const DPoint2d& inertial_center, DMatrix2d& inertial_moments) const
{
	double area = getArea(mc);
	const DVector2d dv = getMiddlePoint() - inertial_center;
	inertial_moments.m[0][0] += area * (dv.y*dv.y); // Ixx
	inertial_moments.m[1][1] += area * (dv.x*dv.x); // Iyy
	inertial_moments.m[0][1] += -area * (dv.x*dv.y); // Ixy
	inertial_moments.m[1][0] += -area * (dv.x*dv.y); // Iyx = Ixy
}

void MeshElement::setTagForEdges(TagExtended::TagType type, int t) const 
{ 
	for(int i = 0; i < count; i++) 
		edges[i]->setIntTag(type, t); 
}

/// Returns the reference to the element incident through the i-th edge (can be nullptr)
MeshElement* MeshElement::getNeighbour(int i) const 
{ 
	return edges[i]->getOtherElement(this); 
}
