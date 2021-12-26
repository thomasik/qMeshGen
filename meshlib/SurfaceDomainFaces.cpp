// SurfaceDomainFaces.cpp: implementation of the SurfaceDomainFaces class.
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include "SurfaceDomainFaces.h"
#include "GeometricPredicates.h"

#include "MeshViewSet.h"
#include "SurfaceParametric.h"
#include "DTriangle.h"
#include "DLine.h"
#include "Metric3dContext.h"
#include "DEquation.h"

SurfaceDomainFaces::SurfaceDomainFaces( 
	const DataVector<DPoint3d> & local_points, 	const DataVector<ITriangle>& local_faces, 
	const DataVector<double>& f_dist2_tol, const DataVector<double>& p_approx_err )
{
	size_t lpct = local_points.countInt();
	size_t lfct = local_faces.countInt();
	// vertices
	m_points.reserve( lpct );
	m_papprox_err.reserve( lpct );
	assert( p_approx_err.countInt() >= lpct );
	for(size_t i = 0; i < lpct; i++) {
		m_points.push_back( local_points[i] );
		m_papprox_err.push_back( p_approx_err[i] );
	}
	// faces
	m_tri_indices.reserve( lfct );
	m_fdist2_tol.reserve( lfct );
	assert( f_dist2_tol.countInt() >= lfct );
	for (size_t i = 0; i < lfct; i++){
		m_tri_indices.push_back( local_faces[i] );
		m_fdist2_tol.push_back( f_dist2_tol[i] );
	}
}

double SurfaceDomainFaces::getInsideQuality(const DPoint3d& point, const DPoint2d& /* param */) const
{
	bool is_within = false;
	auto it = m_tri_indices.begin();
	auto it_end = m_tri_indices.end();
	auto it_dist2 = m_fdist2_tol.begin();
	double quality = -1.0;
	while( it != it_end ){
		const ITriangle tri = *it++;
		double dist2_tol = *it_dist2++;

		if( DTriangle3d::containsProjectedPoint( point, m_points[tri.a], m_points[tri.b], m_points[tri.c], 1e-2 ) ||
			DTriangle3d::distance2ToPoint( point, m_points[tri.a], m_points[tri.b], m_points[tri.c] ) <= dist2_tol ) 
		{
			quality = 0.0;
			if( w_inside_ratio > 0.0 ) {
				quality += w_inside_ratio; // * 1.0;
			}
			if( w_approx_error > 0.0 ) {
				double appr = ( m_papprox_err[tri.a] + m_papprox_err[tri.b] + m_papprox_err[tri.c] ) / 3.0;
				quality += w_approx_error * appr;
			}
			break;
		}
	}

//	if( !is_within ) {
	if( false ) {
		MeshViewSet* set = new MeshViewSet;
		draw(set);
		set->addPoint( point, 1 );
		set->addInfo("is_within", is_within? "true" : "false");
		set->addInfo("quality", quality);
		SHOW_MESH( "SDF::isInside", set );
	}

	return quality;
}

void SurfaceDomainFaces::draw(MeshViewSet* set, const SurfaceParametric* surface, int id) const
{
	if(set == nullptr) return;

	auto it = m_tri_indices.begin();
	auto it_end = m_tri_indices.end();
	while( it != it_end ){
		const ITriangle tri = *it++;
		set->addFace( m_points[tri.a], m_points[tri.b], m_points[tri.c], id );
	}
	set->addInfo("domain type", "SD-Faces" );
}

/// Store XML description to stream
ostream& SurfaceDomainFaces::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<domain-faces>\n";
	//for(auto itr = m_hull.begin(), end = m_hull.end(); itr != end; itr++){
	//	os << prefix << "\t<pt>" << *itr << "</pt>\n";
	//}
	return os << prefix << "</domain-faces>\n";
}
