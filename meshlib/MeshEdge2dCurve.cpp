/////////////////////////////////////////////////////////////////////////////
// MeshEdge2dCurve.cpp
// 2D MeshEdge2d with curvilinear shape
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "MeshEdge2dCurve.h"
#include "Curve2dParametric.h"
#include "SurfaceParametric.h"
#include "ControlSpace2dAdaptive.h"
#include "MeshPoint2d.h"
#include "MeshViewSet.h"

//////////////////////////////////////////////////////////////////////
// Standardowy konstruktor
MeshEdge2dCurve::MeshEdge2dCurve(MeshPoint2d *p1, MeshPoint2d *p2, 
		Curve2dConstPtr new_shape, double t0, double t1, 
		MeshElement* e1, MeshElement* e2)
	: MeshEdge2d(p1, p2, e1, e2), m_shape(new_shape), m_t0(t0), m_t1(t1)
{
	assert(new_shape != nullptr);
	assert(m_t0 != m_t1);
}

/// Create new edge with similar geometry (but possibly adjusted points topology and location)
MeshEdge2d* MeshEdge2dCurve::cloneGeometric(MeshPoint2d *p0, MeshPoint2d *p1, double ksi0, double ksi1, MeshElement* e1, MeshElement* e2) const
{
	double dt = m_t1 - m_t0;
	return new MeshEdge2dCurve(p0, p1, m_shape, m_t0 + dt * ksi0, m_t0 + dt * ksi1, e1, e2);
}
/////////////////////////////////////////////////////////////////////////////
// Zwraca d³ugoœæ ca³kowit¹ krawêdzi 
double MeshEdge2dCurve::getLength(Metric2dContext& mc, bool local_metric) const
{
	return m_shape->getLength(mc, m_t0, m_t1, local_metric);
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca d³ugoœæ czêœci krawêdzi ograniczonej przez parametry t0 i t1.
//	Parametry nale¿¹ do przedzia³u [0, 1]
double MeshEdge2dCurve::getLength(Metric2dContext& mc, double ksi0, double ksi1, bool local_metric) const
{
	double dt = m_t1-m_t0;
	return m_shape->getLength(mc, m_t0 + dt * ksi0, m_t0 + dt * ksi1, local_metric);
}

double MeshEdge2dCurve::getLengthMax(Metric2dContext& mc, double ksi0, double ksi1) const
{
	double dt = m_t1-m_t0;
	return m_shape->getLength(mc, m_t0 + dt * ksi0, m_t0 + dt * ksi1, true);
}

double MeshEdge2dCurve::checkAndGetLength(Metric2dContext& mc, double ksi0, double& ksi1, double max_len, bool local_metric) const
{
	double dt = m_t1-m_t0;
	double t = m_t0 + dt * ksi1;
	double len = m_shape->checkAndGetLength(mc, m_t0 + dt * ksi0, t, max_len, local_metric);
	ksi1 = (t - m_t0) / dt;
	return len;
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca punkt nale¿¹cy do krawêdzi na podstawie wartoœci parametru t [0, 1]
DPoint2d MeshEdge2dCurve::getPoint(double ksi) const
{
	return m_shape->getPoint(m_t0 + (m_t1 - m_t0) * ksi);
}

//////////////////////////////////////////////////////////////////////
// Zmienia kierunek krawêdzi na przeciwny (zamienia elementy zwi¹zane
//	z obywdwoma stronami)
void MeshEdge2dCurve::switchSide()
{
	MeshEdge2d::switchSide();

	double t = m_t0;
	m_t0 = m_t1;
	m_t1 = t;
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca tablicê punktów definiuj¹cych krzyw¹ ³aman¹ oPIsuj¹c¹ kszta³t tej 
//	krawêdzi
void MeshEdge2dCurve::getPolyLine(DataVector<DPoint2d> & polyline) const
{
	DataVector<double> poly_double;
	m_shape->getPolyLineInRange(m_t0, m_t1, poly_double);
	size_t count = poly_double.countInt();
	polyline.prepare(count);
	for(size_t i = 0; i < count; i++)
		polyline.add(m_shape->getPoint(poly_double[i]));
}

/// Returns the approximation of this edge via polyline (array of points)
void MeshEdge2dCurve::getPolyLine(DataVector<DPoint3d> & polyline, SurfaceConstPtr surface) const
{
	surface->getPolyLine(polyline, m_shape, m_t0, m_t1);
}


/// Includes this edge into bounding rectangle
void MeshEdge2dCurve::addToBoundingRect(DRect& rect) const
{
	rect.addRect(m_shape->getBoundingRect(m_t0, m_t1));
}

/// numerical approximation -> calculate parameter for point on surface+curve
DPoint2d MeshEdge2dCurve::surfaceParameters(SurfaceConstPtr surface, 
			const DPoint3d& pt, double & ksi, bool /* same_direction */) const
{
	double t = surface->getShapeParameters(pt, m_shape, m_t0 + (m_t1-m_t0) * ksi, 
		std::min(m_t0, m_t1), std::max(m_t0, m_t1));
	ksi = (t - m_t0) / (m_t1 - m_t0);

	const DPoint2d result = m_shape->getPoint(t);

	double diff = pt.distance(surface->getPoint(result));
	if(diff > 10*mesh_data.relative_small_number){
		MeshViewSet* set = new MeshViewSet;
		DataVector<DPoint3d> polyline(100);
		getPolyLine(polyline, surface);
		for(size_t i = 1; i < polyline.countInt(); i++)
			set->addEdge(polyline[i-1], polyline[i]);
		set->addPoint(surface->getPoint(m_shape->getPoint(m_t0)), 1, 0);
		set->addPoint(surface->getPoint(m_shape->getPoint(m_t1)), 1, 1);
		set->addPoint(pt, 2, 2);
		set->addPoint(surface->getPoint(result), 3, 3);
		ostringstream oss;
		oss << "MeshEdge2dCurve::surfaceParameters (" << diff << ")";
		SHOW_MESH(oss.str(), set);
	}

	return result;
}

/// Returns the parameter ksi for a given point (numerical approximation)
double MeshEdge2dCurve::getParameter(const DPoint2d& pt) const
{
	double t = m_shape->getParameter(pt, 0.5*(m_t0+m_t1));
	return (t - m_t0) / (m_t1 - m_t0);
}


/// Update given ACS for curvature of edge contour
bool MeshEdge2dCurve::updateACSwithCurvature(CS2dAPtr space, SurfaceConstPtr surface, double d2_threshold) const
{
	const DPoint3d mpt0 = surface->getPoint(points[0]->getCoordinates());
	const DPoint3d mpt1 = surface->getPoint(points[1]->getCoordinates());
	const DPoint3d p0 = surface->getPoint(MeshEdge2d::getPoint(0.5));	// without shape
	const DPoint3d p1 = surface->getPoint(getPoint(0.5));				// with shape
	double d2 = p0.distance2(p1) / mpt0.distance2(mpt1);
	if(d2 > d2_threshold)
		return space->updateForBoundaryShape(this);
	else
		return false;
}

double MeshEdge2dCurve::getNonPlanarCurvature(SurfaceConstPtr surface, double ksi, double* cdt_len) const
{
	return  m_shape->getNonPlanarCurvature(surface, m_t0 + (m_t1-m_t0) * ksi, cdt_len);
}

/// Returns curvature of edge shape on plane
double MeshEdge2dCurve::getPlanarCurvature(double ksi) const
{
	return m_shape->getPlanarCurvature(m_t0 + (m_t1-m_t0) * ksi);
}
