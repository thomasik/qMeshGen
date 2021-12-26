// SurfaceCylinder.cpp: implementation of the SurfaceCylinder class.
// Tomasz Jurczyk, 2006-
// Generation of unstructured meshes
//
//////////////////////////////////////////////////////////////////////

#include "SurfaceCylinder.h"
#include "DLeastSquaresFitting.h"
#include "MeshFace.h"
#include "MeshPoint3d.h"

#define CYLINDER_PARAM_MARGIN_X 0.05

/// Standard constructor

SurfaceCylinder::SurfaceCylinder(const DPoint3d& pnt0, const DPoint3d& pnt1, double radius) 
	: SurfaceParametric(), m_cylinder(pnt0, (pnt1-pnt0).normalized(), radius), m_non_periodic(false)
{
	double l1 = 1.0; //m_cylinder.m_axis_vt.length2();
	double l2 = radius*radius;
	m_dratio = std::min(l1, l2) / std::max(l1, l2);

	m_valid = true;
}

SurfaceCylinder::SurfaceCylinder(const DPoint3d& pnt0, const DVector3d& vt, double radius) 
	: SurfaceParametric(), m_cylinder(pnt0, vt, radius), m_non_periodic(false)
{
	double l1 = 1.0; // vt.length2();
	double l2 = radius*radius;
	m_dratio = std::min(l1, l2) / std::max(l1, l2);

	m_valid = true;
}

SurfaceCylinder::SurfaceCylinder(const DCylinder & cylinder)
	: SurfaceParametric(), m_cylinder(cylinder), m_non_periodic(false)
{
	double l1 = 1.0; //m_cylinder.m_axis_vt.length2();
	double l2 = sqr(m_cylinder.m_radius);
	m_dratio = std::min(l1, l2) / std::max(l1, l2);

	m_valid = true;
}

DPoint2d SurfaceCylinder::shiftedParams(const DPoint2d& param) const
{
	assert( m_non_periodic );
	double x = param.x + m_param_shift_x;
	while( x < 0.0 ) x += 1.0;
	while( x > 1.0 ) x -= 1.0;
	return DPoint2d( x, param.y );
}

DPoint2d SurfaceCylinder::unshiftedParams(const DPoint2d& param) const
{
	assert( m_non_periodic );
	double x = param.x - m_param_shift_x;
	while( x < 0.0 ) x += 1.0;
	while( x > 1.0 ) x -= 1.0;
	return DPoint2d( x, param.y );
}

const DPoint3d SurfaceCylinder::getPoint(const DPoint2d& param) const
{
	return m_cylinder.getPoint( m_non_periodic ? shiftedParams(param) : param );
}

/// Returns the principle curvature for this plane and parameters [s,t]
const SurfaceCurvature SurfaceCylinder::getCurvature(const DPoint2d& param, double & g_ratio) const
{ 
	g_ratio = m_dratio; 
	double tmp;
	return m_cylinder.getCurvature( m_non_periodic ? shiftedParams(param) : param, tmp );
}

const DPoint2d SurfaceCylinder::getParameters(const DPoint3d& point) const
{
	return m_non_periodic ? 
		unshiftedParams( m_cylinder.getParam( point ) ) :
		m_cylinder.getParam( point );
}

/// Returns the parameters of the surface for the given point (numerical approximation) and signed distance (according to normal vecor)
bool SurfaceCylinder::getParametersAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z ) const 
{
	bool ok = m_cylinder.getParamAndSignedDistance( point, param, z );
	if(m_non_periodic) param = unshiftedParams( param );
	return ok;
}

const DVector3d SurfaceCylinder::getNormalVector(const DPoint2d& param) const
{
	return m_cylinder.getNormalVector( m_non_periodic ? shiftedParams(param) : param );
}

const DVector3d SurfaceCylinder::getNormalVectorDerivative(int deriv, const DPoint2d& param) const
{
	return m_cylinder.getNormalVectorDerivative(deriv, m_non_periodic ? shiftedParams(param) : param);
}

/// Returns the normal vector to surface for the given parameters
const DVector3d SurfaceCylinder::getNormalVectorForPoint3d(const DPoint3d& pt) const
{
	return m_cylinder.getNormalVectorForPoint3d( pt );
}

const DVector3d SurfaceCylinder::getDerivative(int deriv, const DPoint2d& param) const
{
	return m_cylinder.getDerivative( deriv, m_non_periodic ? shiftedParams(param) : param );
}

/// invert the orientation of the surface (change diretion of normal vector) if possible
bool SurfaceCylinder::invertOrientation()
{
	m_non_periodic = false; // reset
	return m_cylinder.invertOrientation();
}

/// Store XML description to stream
ostream& SurfaceCylinder::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<cylinder>\n";
	os << prefix << "\t<pt0> "		<< m_cylinder.m_center		<< " </pt0>\n";
	os << prefix << "\t<axis> "		<< m_cylinder.m_axis_vt		<< " </axis>\n";
	os << prefix << "\t<radius> "	<< m_cylinder.m_radius		<< " </radius>\n";
	return os << prefix << "</cylinder>\n";
}

std::shared_ptr<SurfaceParametric> SurfaceCylinder::adjustedForFaces( const DataVector< MeshFace* > & mfaces,
		const DataVector< MeshPoint3d* > & mpoints, MeshFace* central_face  ) const
{
	assert( mfaces.notEmpty() );
	assert( mpoints.notEmpty() );
	if( mfaces.empty() || mpoints.empty() ) return nullptr;

	size_t mfct = mfaces.countInt();
	DataVector< DVector3d > cvectors( mfct );
	for(size_t i = 0; i < mfct; i++)
		cvectors.add( mfaces[i]->getBaseNormal() );
	DVector3d caxis;

	double res = DLeastSquaresFitting::fitVectorOrthonormal( cvectors, caxis );
	assert( res != LS_FIT_ERROR );

	size_t mpct = mpoints.countInt();
	DataVector< DPoint3d > cpoints( mpct );
	for(size_t i = 0; i < mpct; i++)
		cpoints.add( mpoints[i]->getCoordinates() );

	DCylinder cylinder;

	res = DLeastSquaresFitting::fitCylinder( cpoints, caxis, cylinder ); 
	assert( res != LS_FIT_ERROR );

	// orientation...

	auto surf = std::make_shared<SurfaceCylinder>( cylinder );
	if( central_face != nullptr ){
		surf->setOrientationLikeFace( central_face );
		surf->fixCenterAndRangeForPeriodic( central_face->getMiddlePoint() );
	}

	return surf;
}

// fix center and range (for parameters) to fix periodic surfaces
bool SurfaceCylinder::fixCenterAndRangeForPeriodic(const DPoint3d& center_point ) 
{ 
	DPoint2d center_param = m_cylinder.getParam( center_point );

//	DPoint3d check_point = m_cylinder.getPoint( center_param );
//	DPoint3d cc_point = getPoint( center_param );
//	double check_len = center_point.distance( check_point );

	// calculate param_shift and param_range
	m_param_shift_x = center_param.x - 0.5;
	m_non_periodic = true;

//	DPoint3d check_point_after = getPoint( DPoint2d( 0.5, center_param.y ) );
//	double check_len_after = center_point.distance( check_point_after );

	return m_non_periodic; 
}

bool SurfaceCylinder::withinParamRange( const DPoint2d& param ) const 
{ 
	return !m_non_periodic || (
		(param.x >= CYLINDER_PARAM_MARGIN_X) && (param.x <= (1.0-CYLINDER_PARAM_MARGIN_X)) );	
}
