// SurfaceDomainFacesPlanar.cpp: implementation of the SurfaceDomainFacesPlanar class.
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include "SurfaceDomainFacesPlanar.h"
#include "GeometricPredicates.h"

#include "MeshViewSet.h"
#include "SurfaceParametric.h"
#include "DTriangle.h"
#include "DSegment.h"
#include "DLine.h"
#include "Metric3dContext.h"
#include "DEquation.h"

SurfaceDomainFacesPlanar::SurfaceDomainFacesPlanar( 
	const DataVector<DPoint2d> & local_points, 	
	const DataVector<ITriangle>& local_faces, 
	const DataVector<double>& p_approx_quality )
{
	size_t lpct = local_points.countInt();
	size_t lfct = local_faces.countInt();
	// vertices
	m_points.reserve( lpct );
	m_papprox_quality.reserve( lpct );
	DRect box;
	assert( p_approx_quality.countInt() >= lpct );
	for(size_t i = 0; i < lpct; i++) {
		box.addPoint( local_points[i] );
		m_points.push_back( local_points[i] );
		m_papprox_quality.push_back( p_approx_quality[i] );
	}
	m_eps = 1e-3 * box.getDiameter();

	// faces
	m_tri_indices.reserve( lfct );
	inv_count = 0;
	double inv_min = 0.0;
	for (size_t i = 0; i < lfct; i++){
		const ITriangle& tri = local_faces[i];
		const DPoint2d& a = m_points[tri.a];
		const DPoint2d& b = m_points[tri.b];
		const DPoint2d& c = m_points[tri.c];
		double det = GeometricPredicates::orient2d( a, b , c );
		if( det <= 0.0 ) {
			// stat
			if( inv_count == 0 ) inv_min = det;
			else if( inv_min > det ) inv_min = det;
			inv_count++;
			m_tri_indices.push_back( ITriangle(tri.a, tri.c, tri.b) ); // invert inverted...
		}else
			m_tri_indices.push_back( tri );
	}

	if( inv_count == lfct ) {
		inv_count = 0; // all were inverted, so it should be ok...
	}

	if( inv_count > 0 ){
//	if( true ){
		MeshViewSet * set = new MeshViewSet;
		draw(set);
		set->addInfo("inv_count", inv_count );
		set->addInfo("inv_min", inv_min );
		SHOW_MESH("SDFP::create", set );
	}
}

DPoint2d SurfaceDomainFacesPlanar::getMiddleParam() const
{
	if( m_points.empty() ) return DPoint2d::zero;
	double f = 1.0 / m_points.size();
	DPoint2d mid = DPoint2d::zero;
	for(auto itr = m_points.begin(), end = m_points.end(); itr != end; itr++)
		mid.add(*itr, f);
	return mid;
}

double SurfaceDomainFacesPlanar::getInsideQuality(const DPoint3d& /* point */, const DPoint2d& param) const
{
	bool is_within = false;
	auto it = m_tri_indices.begin();
	auto it_end = m_tri_indices.end();

	//static int local_counter = 0;
	//local_counter++;

	double quality = AQ_INVALID;
	while( it != it_end ){
		const ITriangle& tri = *it++;

		const DPoint2d& a = m_points[tri.a];
		const DPoint2d& b = m_points[tri.b];
		const DPoint2d& c = m_points[tri.c];

		assert( GeometricPredicates::orient2d( a, b, c ) >=  0.0 );

		double ab = GeometricPredicates::orient2d( a, b, param);
		double bc = GeometricPredicates::orient2d( b, c, param);
		double ca = GeometricPredicates::orient2d( c, a, param);

		is_within = (ab > -m_eps) && (bc > -m_eps) && (ca > -m_eps);

		if( false ){
			MeshViewSet* set = new MeshViewSet;
			for(unsigned int i = 0; i < m_tri_indices.size(); i++){
				ITriangle itri = m_tri_indices[i];
				set->addEdge( DPoint3d( m_points[ itri.a], 0.0), DPoint3d( m_points[ itri.b], 0.0), -1 );
				set->addEdge( DPoint3d( m_points[ itri.b], 0.0), DPoint3d( m_points[ itri.c], 0.0), -1 );
				set->addEdge( DPoint3d( m_points[ itri.c], 0.0), DPoint3d( m_points[ itri.a], 0.0), -1 );
			}
			set->addEdge( DPoint3d( a, 0.0), DPoint3d( b, 0.0), 1 );
			set->addEdge( DPoint3d( b, 0.0), DPoint3d( c, 0.0), 1 );
			set->addEdge( DPoint3d( c, 0.0), DPoint3d( a, 0.0), 1 );
			set->addPoint( DPoint3d( param, 0.1 ), 1 );
			set->addInfo("ab", ab);
			set->addInfo("bc", bc);
			set->addInfo("ca", ca);
			set->addInfo("eps", m_eps);
			set->addInfo("is_within", is_within);
			SHOW_MESH("SDFP::isInside", set);
		}

		//is_within = DTriangle2d::containsPoint(a, b, c, param, 1e-2 );

		if( is_within )
		{
			assert( w_approx_error == 1.0 );
			if( ab < 0.0 ) ab = 0.0;
			if( bc < 0.0 ) bc = 0.0;
			if( ca < 0.0 ) ca = 0.0;
			double abc = ab+bc+ca;
			if( abc > 0.0 ) {
				quality = ( bc * m_papprox_quality[tri.a] + ca * m_papprox_quality[tri.b] + 
							ab * m_papprox_quality[tri.c] ) / abc;
			}else{
				double dist2_a = param.distance2(a);
				double dist2_b = param.distance2(b);
				double dist2_c = param.distance2(c);
				// select nearest...
				if( dist2_a < dist2_b )
					quality = m_papprox_quality[(dist2_c < dist2_a ) ? tri.c : tri.a];
				else
					quality = m_papprox_quality[(dist2_c < dist2_b ) ? tri.c : tri.b];
			}
			//assert( abc > 0.0 );
			//quality += w_approx_error * appr;
			AQ_ASSERT( quality );
			break;
		}
	}

//	if( !is_within ) {
	if( false ) {
		MeshViewSet* set = new MeshViewSet;
		draw(set);
		set->addPoint( DPoint3d(param, 0.01), 1 );
		set->addInfo("is_within", is_within? "true" : "false");
		set->addInfo("quality", quality);
		SHOW_MESH( "SDFP::isInside", set );
	}

	return quality;
}

void SurfaceDomainFacesPlanar::draw(MeshViewSet* set, const SurfaceParametric* surface, int id) const
{
	if(set == nullptr) return;

	auto it = m_tri_indices.begin();
	auto it_end = m_tri_indices.end();
	while( it != it_end ){
		const ITriangle tri = *it++;
		if( surface ){
			set->addFace( 
				surface->getPoint(m_points[tri.a]), 
				surface->getPoint(m_points[tri.b]), 
				surface->getPoint(m_points[tri.c]), id );
		}else{
			set->addFace( 
				DPoint3d(m_points[tri.a], 0.0), 
				DPoint3d(m_points[tri.b], 0.0), 
				DPoint3d(m_points[tri.c], 0.0), id );
		}
	}

	if(surface){
		DataVector< DPoint2d > params( m_points.size() );
		auto itp = m_points.begin();
		auto itp_end = m_points.end();
		DRect box;
		while( itp != itp_end ){
			const DPoint2d& p = *itp++;
			params.add( p );
			box.addPoint( p );
		}
		surface->createViewSetForPoints( set, params );
	}

	set->addInfo("domain type", "SD-FacesPlanar" );
}

/// Store XML description to stream
ostream& SurfaceDomainFacesPlanar::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<domain-faces-planar>\n";
	os << prefix << "\t<points>\n";
	size_t count = m_points.size();
	for(size_t i = 0; i < count; i++){
		os << prefix << "\t<pt id=\"" << i << "\">" << m_points[i] 
			<< "</pt> <approx-quality>" << m_papprox_quality[i] << "</approx-quality>\n";
	}
	os << prefix << "\t</points>\n";
	os << prefix << "\t<triangles>\n";
	count = m_tri_indices.size();
	for(size_t i = 0; i < count; i++){
		const ITriangle& t = m_tri_indices[i];
		os << prefix << "\t<triangle>" << t.a << " " << t.b << " " << t.c << "</triangle>\n";
	}
	os << prefix << "\t</triangles>\n";
	return os << prefix << "</domain-faces-planar>\n";
}
