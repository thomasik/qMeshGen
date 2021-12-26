// SurfaceDomainHull.cpp: implementation of the SurfaceDomainHull class.
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2013-
//	Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include<algorithm>

#include "SurfaceDomainHull.h"
#include "GeometricPredicates.h"

#include "MeshViewSet.h"
#include "SurfaceParametric.h"
#include "DTriangle.h"
#include "DLine.h"
#include "Metric3dContext.h"
#include "DEquation.h"

SurfaceDomainHull::SurfaceDomainHull(const DataVector<DPoint2d> & hull)
{
	for(size_t i = 0; i < hull.countInt(); i++)
		m_hull.push_back( hull[i] );
}

DPoint2d SurfaceDomainHull::getMiddleParam() const { 
	if( m_hull.empty() ) return DPoint2d::zero; 
	else if(!m_hdist.empty() ) return m_middle;
	else{
		double f = 1.0 / m_hull.size();
		DPoint2d mid = DPoint2d::zero;
		for(auto itr = m_hull.begin(), end = m_hull.end(); itr != end; itr++)
			mid.add(*itr, f);
		return mid;
	}
}

double SurfaceDomainHull::getInsideQuality(const DPoint3d& /* point */, const DPoint2d& param) const
{
	if( m_hdist.empty() ){
		assert( !m_hull.empty() );
		for(unsigned int i = 1; i < m_hull.size(); i++){
			double det = GeometricPredicates::orient2d( m_hull[i-1], m_hull[i], param );
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "det_" << i << "\t" << det);
			if(det < -SMALL_NUMBER) return AQ_INVALID;
		}
		return true;
	}else{
		// ... calculate "angle" for the given param
		const DVector2d vp = param - m_middle;
		double vp_len = vp.length();
		if( vp_len < SMALL_NUMBER ) {
			double quality = 0.0;
			if( w_inside_ratio > 0.0 ) {
				quality = w_inside_ratio; // *1.0;
			}
			if( w_approx_error > 0.0 ) {
				assert( false ); // not available?
			}
			return quality;
		}
		const DVector2d vpn = vp.normalized();
		double sp = (vpn.x > 0) ? (1 + vpn.y) : (3 - vpn.y); // 0..4
		// ... find proper pair of values in m_hdist
		int l = 0, r = m_hdist.countInt()-1;
		while( r - l > 1 ) {
			int c = (l+r)/2;
			if( m_hdist[c] > sp ) r = c; else l = c;
		}
		// ... calculate distances
		const int N1 = m_hdist.countInt();
		const int N2 = (int)m_hull.size();
		assert( N1 == N2+2 );
		int hl = ( ((l > 0) ? l : N2) + m_hoffset - 1) % N2;
		//int hr = ( ((r > N2) ? 1 : r) + m_hoffset) % N2;
		int hr = (hl+1) % N2;
		DPoint2d cpt = DLine2d::crossPoint( m_hull[hl], m_hull[hr], m_middle, param );
		if(true){
			DLine2d line_s(m_hull[hl], m_hull[hr]);
			double s = line_s.paramForPoint( cpt );
			double max_s = line_s.paramForPoint( m_hull[hr] );
			DLine2d line_t(m_middle, param);
			double t = line_t.paramForPoint( cpt );
			double max_t = line_t.paramForPoint( param );
//			if( s < 0.0 || s > max_s || t < 0.0 ) {
			if( false ) {
				MeshViewSet* set = new MeshViewSet;
				DPoint3d prev( m_hull[ N2 - 1 ], 0.0 );
				for(int i = 0; i < N2; i++){
					DPoint3d next( m_hull[ i ], 0.0 );
					set->addEdge( prev, next );
					set->addPoint( next, 0, i );
					prev = next;
				}
				set->addInfo("hl-hr", to_string(hl) + " - " + to_string(hr) );
				set->addInfo("s / max_s", s / max_s);
				set->addInfo("t / max_t", t / max_t);
				set->addLabel( DPoint3d( m_middle, 0.0), "mid" );
				set->addLabel( DPoint3d( param, 0.0), "pt" );
				set->addLabel( DPoint3d( cpt, 0.0), "c" );
				SHOW_MESH(" SDH::isInside / crossPoint ", set );
			}
		}
		double cr = vp_len / m_middle.distance( cpt );
		// check
		bool within = (cr <= (1.0 + SMALL_NUMBER) );
		if( within ) {
			double quality = 0.0;
			if( w_inside_ratio > 0.0 ) {
				quality = w_inside_ratio * (1.0 - cr);
			}
			if( w_approx_error > 0.0 ) {
				assert( false ); // not available?
			}
			AQ_ASSERT( quality );
			return quality;
		}else
			return AQ_INVALID;
	}
}

void SurfaceDomainHull::addPoint(const DPoint2d& pt)
{
	m_points.push_back(pt);
}

bool SurfaceDomainHull::createHull()
{
	if(m_points.size() < 3) return false;

	// Implementation of Andrew's monotone chain 2D convex hull algorithm.
	// Asymptotic complexity: O(n log n).

	int n = (int)m_points.size();
	int k = 0;

	m_hull.clear();
	m_hull.resize(2*n);

	// Sort points lexicographically
	sort( m_points.begin(), m_points.end() );

	// Build lower hull
	for (int i = 0; i < n; i++) {
		while (k >= 2 && GeometricPredicates::orient2d(m_hull[k-2], m_hull[k-1], m_points[i]) <= SMALL_NUMBER) k--;
		m_hull[k++] = m_points[i];
	}

	// Build upper hull
	for (int i = n-2, t = k+1; i >= 0; i--) {
		while (k >= t && GeometricPredicates::orient2d(m_hull[k-2], m_hull[k-1], m_points[i]) <= SMALL_NUMBER) k--;
		m_hull[k++] = m_points[i];
	}
	
	//if(true){
	//	MeshViewSet* set = new MeshViewSet();
	//	for(int i = 0; i < n; i++)
	//		set->addPoint( DPoint3d(m_points[i], 0.0) );
	//	for(int i = 1; i < k; i++)
	//		set->addPoint( DPoint3d(m_hull[i], 0.1), 1, i);
	//	for(int i = 1; i < k; i++)
	//		set->addEdge( DPoint3d(m_hull[i-1], 0.1), DPoint3d(m_hull[i], 0.1), 2);
	//	SHOW_MESH("2d-hull of k-points", set);
	//}

	m_hull.resize(k-1);

	/*
	int inside_count = 0;
	for(int i = 0; i < n; i++) {
		bool inside = isInside( m_points[i] );
		if(inside) ++inside_count;
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "DomainHull - inside test: " << inside_count << "/" << n);
	*/

	return true;
}

DPoint2d SurfaceDomainHull::calculateHullMiddle( Metric3dContext & mc, const SurfaceParametric* surface ) const
{
	if( m_hull.empty() ) return DPoint2d::zero;

	DPoint2d mass_mid = DPoint2d::zero;

	int hct = (int)m_hull.size();
	assert( hct >= 3 );
	auto itr = m_hull.begin(), end = m_hull.end();
	DPoint2d first = *itr++;
	DPoint2d prev = *itr++;
	double total_area = 0.0;
	while( itr != end ) {
		DPoint2d next = *itr++;
		DPoint2d tmid( first, prev, next ); // average of three points
		mc.countMetricAtPoint( surface->getPoint( tmid ) );
		double f = DMTriangle3d::area(
			mc.transformRStoMS( surface->getPoint( first ) ),
			mc.transformRStoMS( surface->getPoint( prev ) ),
			mc.transformRStoMS( surface->getPoint( next ) ) );

		total_area += f;
		mass_mid.add( tmid, f );
		prev = next;
	}

	mass_mid /= total_area;

	return mass_mid;
}

bool SurfaceDomainHull::createHullDist( Metric3dContext & mc, const SurfaceParametric* surface )
{
	if( m_hull.empty() ) return false;

	m_middle = calculateHullMiddle( mc, surface );

	int hct = (int)m_hull.size();
	DataVector<double> temp_hdist(hct);
	m_hoffset = 0;

	const DVector3d ds = surface->getDerivative( DEquation::deriv_ds, m_middle );
	const DVector3d dt = surface->getDerivative( DEquation::deriv_dt, m_middle );
	mc.countMetricAtPoint( surface->getPoint( m_middle ) );

	m_slen2 = mc.transformRStoMS( ds ).length2();
	m_tlen2 = mc.transformRStoMS( dt ).length2();

	int i = -1;
	for(auto itr = m_hull.begin(), end = m_hull.end(); itr != end; itr++){
		const DVector2d vpn = (*itr - m_middle).normalized();
		double sp = (vpn.x > 0) ? (1 + vpn.y) : (3 - vpn.y); // 0..4
		int index = temp_hdist.add( sp );
		if( (++i > 0) && (sp < temp_hdist[m_hoffset]) ) m_hoffset = i;
	}

	// sort and wrap angles
	m_hdist.clear();
	m_hdist.prepare( hct + 2 );
	m_hdist.add( temp_hdist[ (m_hoffset+hct-1)%hct] - 4.0 );
	for(int i = 0; i < hct; i++) 
		m_hdist.add( temp_hdist[ (m_hoffset+i)%hct] );
	m_hdist.add( temp_hdist[ m_hoffset ] + 4.0 );

	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "===== hdist [angle ] ===========");
	//for(int i = 0; i < m_hdist.countInt(); i++) {
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, i << "\t" << m_hdist[i]);
	//}
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "==============================================");

	return true;
}

void SurfaceDomainHull::draw(MeshViewSet* set, const SurfaceParametric* surface, int id) const
{
	if(set == nullptr) return;
	if(surface == nullptr) return;
	if(m_hull.empty()) return;

	DPoint3d last_pt = surface->getPoint( *m_hull.rbegin() );
	DPoint3d pt;
	DPoint2d middle;
	if( m_hdist.empty() ){
		double f = 1.0 / m_hull.size();
		for(auto itr = m_hull.begin(), end = m_hull.end(); itr != end; itr++)
			middle.add(*itr, f);
	} else
		middle = m_middle;

	DPoint3d smiddle = surface->getPoint(middle);
	set->addPoint( smiddle, id, id );
	for(auto itr = m_hull.begin(), end = m_hull.end(); itr != end; itr++){
		pt = surface->getPoint(*itr);
		set->addEdge( last_pt, pt, id );
		DPoint3d hpt = surface->getPoint( DPoint2d( *itr, middle, 0.5) );
		set->addEdge( smiddle, hpt, id );
		set->addEdge( hpt, pt, id );
		last_pt = pt;
	}
	set->addInfo("domain type", "SD-Hull" );
}

/// Store XML description to stream
ostream& SurfaceDomainHull::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<domain-hull>\n";
	for(auto itr = m_hull.begin(), end = m_hull.end(); itr != end; itr++){
		os << prefix << "\t<pt>" << *itr << "</pt>\n";
	}
	return os << prefix << "</domain-hull>\n";
}
