// SurfaceSphere.cpp: implementation of the SurfaceSphere class.
// Tomasz Jurczyk, 2015-
// Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include "SurfaceSphere.h"
#include "DLeastSquaresFitting.h"
#include "DataVector.h"
#include "MeshPoint3d.h"

#define SPHERE_PARAM_MARGIN_X 0.05
#define SPHERE_PARAM_MARGIN_Y 0.05

/// Standard constructor

SurfaceSphere::SurfaceSphere(const DPoint3d& pnt, double radius) 
	: SurfaceParametric(), m_sphere(pnt, radius), m_param_shift_x(0.0), m_non_periodic(false)
{
	m_valid = (radius != 0.0);
}

SurfaceSphere::SurfaceSphere(const DSphere & sphere) 
	: SurfaceParametric(), m_sphere( sphere ), m_param_shift_x(0.0), m_non_periodic(false)
{
	m_valid = (sphere.m_radius != 0.0);
}

DPoint2d SurfaceSphere::shiftedParams(const DPoint2d& param) const
{
	assert( m_non_periodic );
	double x = param.x + m_param_shift_x;
	while( x < 0.0 ) x += 1.0;
	while( x > 1.0 ) x -= 1.0;
	return DPoint2d( x, param.y );
}

DPoint2d SurfaceSphere::unshiftedParams(const DPoint2d& param) const
{
	assert( m_non_periodic );
	double x = param.x - m_param_shift_x;
	while( x < 0.0 ) x += 1.0;
	while( x > 1.0 ) x -= 1.0;
	return DPoint2d( x, param.y );
}

const DPoint3d SurfaceSphere::getPoint(const DPoint2d& param) const
{
	return m_sphere.getPoint( m_non_periodic ? shiftedParams(param) : param );
}

/// Returns the principle curvature for this sphere and parameters [s,t]
const SurfaceCurvature SurfaceSphere::getCurvature(const DPoint2d& param, double & g_ratio) const
{ 
	return m_sphere.getCurvature( m_non_periodic ? shiftedParams(param) : param, g_ratio );
}

const DPoint2d SurfaceSphere::getParameters(const DPoint3d& point) const
{
	return m_non_periodic ? unshiftedParams( m_sphere.getParam( point ) )
		: m_sphere.getParam( point );
}

/// Returns the parameters of the surface for the given point (numerical approximation) and signed distance (according to normal vecor)
bool SurfaceSphere::getParametersAndSignedDistance(const DPoint3d& point, DPoint2d& param, double & z ) const 
{
	bool ok = m_sphere.getParamAndSignedDistance( point, param, z );
	if(m_non_periodic) param = unshiftedParams( param );
	return ok;
}

const DVector3d SurfaceSphere::getNormalVector(const DPoint2d& param) const
{
	return m_sphere.getNormalVector( m_non_periodic ? shiftedParams(param) : param );
}

const DVector3d SurfaceSphere::getNormalVectorForPoint3d(const DPoint3d& pt) const
{
	return m_sphere.getNormalVectorForPoint3d( pt );
}

const DVector3d SurfaceSphere::getDerivative(int deriv, const DPoint2d& param) const
{
	return m_sphere.getDerivative( deriv, m_non_periodic ? shiftedParams(param) : param );
}

/// invert the orientation of the surface (change diretion of normal vector) if possible
bool SurfaceSphere::invertOrientation()
{
	m_non_periodic = false; // reset
	return m_sphere.invertOrientation();
}

/// Store XML description to stream
ostream& SurfaceSphere::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<sphere>\n";
	os << prefix << "\t<midpt> "	<< m_sphere.m_center	<< " </midpt>\n";
	os << prefix << "\t<radius> "	<< m_sphere.m_radius	<< " </radius>\n";
	return os << prefix << "</sphere>\n";
}

std::shared_ptr<SurfaceParametric>  SurfaceSphere::adjustedForFaces(
	const DataVector< MeshFace* > & mfaces, 
	const DataVector< MeshPoint3d* > & mpoints,
	MeshFace* central_face ) const
{
	assert( mfaces.notEmpty() );
	assert( mpoints.notEmpty() );
	if( mfaces.empty() || mpoints.empty() ) return nullptr;

	DSphere dsphere = m_sphere;
	size_t mpct = mpoints.countInt();
	DataVector< DPoint3d > dpoints( mpct );
	for(size_t i = 0; i < mpct; i++)
		dpoints.add( mpoints[i]->getCoordinates() );

	DLeastSquaresFitting::fitSphereIter( dpoints, dsphere, false, false );

	auto surf = std::make_shared<SurfaceSphere>( dsphere );
	if( central_face != nullptr ) surf->setOrientationLikeFace( central_face );

	return surf;
}

// fix center and range (for parameters) to fix periodic surfaces
bool SurfaceSphere::fixCenterAndRangeForPeriodic(const DPoint3d& center_point ) 
{ 
	// -> reorientate sphere poles ??? from OZ to OX or OY, if better...
	DVector3d cnormal = getNormalVectorForPoint3d( center_point );

	double adx = abs(cnormal.x);
	double ady = abs(cnormal.y);
	double adz = abs(cnormal.z);
	if( adx < ady && adx < adz ) // set OX
		m_sphere.setPoleAxis( DSphere::POLE_OX );
	else if( ady < adz ) // set OY
		m_sphere.setPoleAxis( DSphere::POLE_OY );
	// else leave OZ...
	else m_sphere.setPoleAxis( DSphere::POLE_OZ );

	// calculate param_shift
	DPoint2d center_param = m_sphere.getParam( center_point );
	m_param_shift_x = center_param.x - 0.5;

	return (m_non_periodic = true); 
}

bool SurfaceSphere::withinParamRange( const DPoint2d& param ) const 
{ 
	return !m_non_periodic || (
		(param.x >= SPHERE_PARAM_MARGIN_X) && 
		(param.x <= (1.0-SPHERE_PARAM_MARGIN_X)) &&
		(param.y >= SPHERE_PARAM_MARGIN_Y) && 
		(param.y <= (1.0-SPHERE_PARAM_MARGIN_Y)) );	
}
