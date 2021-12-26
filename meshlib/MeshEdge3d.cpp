// MeshEdge3d.cpp: implementation of the MeshEdge3d class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshEdge3d.h"

#include "MeshPoint3d.h"
#include "MeshEdge2d.h"
#include "MeshFace.h"
#include "MeshDomainSurface.h"
#include "MeshContainer2d.h"
#include "MeshArea.h"
#include "Curve2dParametric.h"
#include "Curve3dParametric.h"
#include "SurfaceParametric.h"
#include "MeshBoundaryCondition.h"
#include "MeshViewSet.h"
#include "MeshTetrahedron.h"
#include "Metric3dContext.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MeshEdge3d::MeshEdge3d(MeshPoint3d* point1, MeshPoint3d* point2, MeshFace* face)
	: local_curve(nullptr)
{
	assert(point1);
	assert(point2);
	assert( point1 != point2 );
	points[0] = point1;
	points[1] = point2;
	point1->addEdgeLink(this);
	point2->addEdgeLink(this);
	if(face)
		faces.add(face);
}

MeshEdge3d::~MeshEdge3d()
{
	assert( !availableTag(TagExtended::TAG_COLLAPSE_3D) );
	assert( !availableTag(TagExtended::TAG_ACTIVE) );
	if(points[0]){ // normal update
		points[0]->removeEdgeLink(this);
		points[1]->removeEdgeLink(this);
		assert(faces.countInt() == 0);
	}else{
		for(size_t i = 0; i < faces.countInt(); i++){
			MeshFace* face = faces[i];
			if(face->removeEdgeLink(this)) delete face;
		}
	}
}

/// Removes link to this point (fo delete-all phase) - returns true if last one
bool MeshEdge3d::removePointLink(MeshPoint3d* point)
{
	if(points[0] == point) points[0] = 0;
	else if(points[1] == point) points[1] = 0;

	return (points[0] == 0 && points[1] == 0);
}

void MeshEdge3d::addFaceLink(MeshFace * const face)
{
	assert(!faces.contains(face));
	faces.add(face);
}

bool MeshEdge3d::removeFaceLink(MeshFace * const face)
{
	faces.remove(face);
	return (faces.empty() && !isBorder());
}

bool MeshEdge3d::removeFaceLinkDisregardBoundary(MeshFace * const face)
{
	faces.remove(face);
	return faces.empty();
}

int MeshEdge3d::getPointIndex(const MeshPoint3d *point) const
{
	if(points[0] == point){
		return 0;
	}else if(points[1] == point){
		return 1;
	}else{
		return -1;
	}
}

double MeshEdge3d::getLength(Metric3dContext& mc, bool local_metric) const
{
	if(local_metric) mc.countMetricAtPoint(getPoint(0.5));
	return points[0]->getMetricCoordinates(mc).distance(
			points[1]->getMetricCoordinates(mc));
}

double MeshEdge3d::getLengthNoMetric() const
{
	return points[0]->getCoordinates().distance(points[1]->getCoordinates());
}

DPoint3d MeshEdge3d::getPoint(double t) const
{
	const DPoint3d& p0 = points[0]->getCoordinates();
	const DPoint3d& p1 = points[1]->getCoordinates();
	return DPoint3d(p0, p1, t);
}

int MeshEdge3d::getFaceIndex( MeshFace* const face) const 
{
	return (int)faces.find(face);
}

std::shared_ptr<MeshViewEdgeData> MeshEdge3d::getViewData() const
{
	char b = 0;
	if(isBorder(TagBorder::FIXED)) b = 3;
	else if(isBorder(TagBorder::RIDGE)) b = 2;
	else if(isBorder()) b = 1;
	return std::make_shared<MeshViewEdgeData>(points[0]->getCoordinates(), points[1]->getCoordinates(), b);
}

/////////////////////////////////////////////////////////////////////////////
// Zwraca punkt nale¿¹cy do krawêdzi na podstawie wartoœci parametru t, 
//	okreœlaj¹cego stosunek d³ugoœci wydzielonego odcinka
DPoint3d MeshEdge3d::getAspectPoint(double ksi) const
{
	return DPoint3d(points[0]->getCoordinates(), points[1]->getCoordinates(), ksi);
}

MeshFace* MeshEdge3d::getFaceToPoint(const MeshPoint3d* mpt) const
{
	for(size_t i = 0; i < faces.countInt(); i++)
		if(faces[i]->getOtherPoint(points[0], points[1]) == mpt) 
			return faces[i];

	assert(false);
	return nullptr;
}

double MeshEdge3d::getLengthQuality(Metric3dContext& mc, bool ext_metric) const
{
	int qms = MeshTetrahedron::param_quality_metric_selection;
	MeshTetrahedron::param_quality_metric_selection = 
		ext_metric ? MeshData::QM_VERTICES_AVE : MeshData::QM_MIDDLE;
	mc.countMetricAtPoints(points[0], points[1]);
	MeshTetrahedron::param_quality_metric_selection = qms;

	auto ret = points[0]->getMetricCoordinates(mc).distance(
			points[1]->getMetricCoordinates(mc));
	return ret;
}

double MeshEdge3d::getLengthMetricAdapted(Metric3dContext& mc) const
{
	return getLengthMetricAdapted(mc, 0.0, 1.0);
}

double MeshEdge3d::getLengthMetricAdapted(Metric3dContext& mc, double t0, double t1, int lev) const
{
	const DPoint3d pt0 = getPoint(t0);
	const DPoint3d pt1 = getPoint(t1);
	double tm = 0.5*(t0+t1);
	const DPoint3d ptm = getPoint(tm);
	mc.countMetricAtPoint(DPoint3d::average(pt0, ptm));
	double len_0 = mc.transformRStoMS(ptm-pt0).length();
	mc.countMetricAtPoint(DPoint3d::average(pt1, ptm));
	len_0 += mc.transformRStoMS(ptm-pt1).length();
	if(lev > 5) return len_0;
	mc.countMetricAtPoint(ptm);
	double len_1 = mc.transformRStoMS(pt1-pt0).length();
	if((len_0 - len_1) < 1e-3*len_1) return len_0;
	else return getLengthMetricAdapted(mc, t0, tm, lev+1) +
		getLengthMetricAdapted(mc, tm, t1, lev+1);
}

double MeshEdge3d::getLength() const
{
	return points[0]->getCoordinates().distance(points[1]->getCoordinates());
}

double MeshEdge3d::getLength2() const
{
	return points[0]->getCoordinates().distance2(points[1]->getCoordinates());
}

double MeshEdge3d::getLength2(Metric3dContext& mc) const
{
	return points[0]->getMetricCoordinates(mc).distance2(points[1]->getMetricCoordinates(mc));
}

int MeshEdge3d::approximatePolyline(DataVector<DPoint3d> & polyline) const
{
	polyline.clear();
	polyline.add(points[0]->getCoordinates());
	polyline.add(points[1]->getCoordinates());
	return 2;
}

void MeshEdge3d::switchPoints(const MeshPoint3d *point1, MeshPoint3d *point2)
{
	if(points[0] == point1) points[0] = point2;
	else if(points[1] == point1) points[1] = point2;
}

bool MeshEdge3d::adjacentBlocks(DataVector<MeshBlock*> & blocks, bool boundary_check) const
{
	int fct = faces.countInt();
	blocks.prepare(2*fct);
	for(int i = 0; i < fct; i++){
		MeshFace* face = getFaceAt(i);
		for(int m = 0; m < 2; m++){
			MeshBlock* block = face->getBlock(m);
			if(boundary_check){
				assert(isBorder() || block);
				if(!isBorder() && !block) return false;
			}
			if(block) blocks.addIfNew(block);
		}
	}
	return true;
}

DataVector<const MeshBlock*> MeshEdge3d::adjacentBlocks(bool boundary_check) const
{
	int fct = faces.countInt();
	DataVector<const MeshBlock*> blocks(2 * fct);

	for (int i = 0; i < fct; i++) {
		MeshFace* face = getFaceAt(i);
		for (int m = 0; m < 2; m++) {
			MeshBlock* block = face->getBlock(m);
			if (boundary_check) {
				assert(isBorder() || block);
				if (!isBorder() && !block) return false;
			}
			if (block) blocks.addIfNew(block);
		}
	}

	return blocks;
}

DataVector<MeshBlock*> MeshEdge3d::adjacentBlocks(bool boundary_check)
{
	int fct = faces.countInt();
	DataVector<MeshBlock*> blocks(2 * fct);

	for (int i = 0; i < fct; i++) {
		MeshFace* face = getFaceAt(i);
		for (int m = 0; m < 2; m++) {
			MeshBlock* block = face->getBlock(m);
			if (boundary_check) {
				assert(isBorder() || block);
				if (!isBorder() && !block) return false;
			}
			if (block) blocks.addIfNew(block);
		}
	}

	return blocks;
}

MeshPoint3d* MeshEdge3d::commonVertex(const MeshEdge3d* edge) const
{
	if(points[0] == edge->points[0]) return points[0];
	else if(points[0] == edge->points[1]) return points[0];
	else if(points[1] == edge->points[0]) return points[1];
	else if(points[1] == edge->points[1]) return points[1];
	else return nullptr;
}

/// Returns the other boundary face of the edge (supposing there are only two!)
MeshFace* MeshEdge3d::getOtherBorderFace(const MeshFace* face) const
{
	for(size_t i = 0; i < faces.countInt(); i++){
		if(faces[i] != face && faces[i]->isBorder()) return faces[i];
	}
	return nullptr;
}

/*
//#define COUNT_GET_POINT_ON_SURFACE_CASES
#ifdef COUNT_GET_POINT_ON_SURFACE_CASES
#define INC_COUNTER( x ) x++
#else
#define INC_COUNTER( x )
#endif

DPoint3d MeshEdge3d::getPointOnSurface( Metric3dContext& mc, double t) const
{
#ifdef COUNT_GET_POINT_ON_SURFACE_CASES
	static int counter_curve = 0;
	static int counter_same_surface = 0;
	static int counter_surface_applicable = 0;
	static int counter_one_candidate_ok = 0;
	static int counter_one_candidate_fail = 0;
	static int counter_best_candidate_ok = 0;
	static int counter_best_candidate_fail = 0;
	static int counter_no_candidate = 0;
#endif

	DataVector<DPoint3d> polyline(100);
	if( isBorder() ) {
		Curve3dConstPtr* curve = getLocalCurve();
		assert( curve != nullptr );
		double t0 = curve->getParameter( points[0]->getCoordinates() );
		double t1 = curve->getParameter( points[1]->getCoordinates() );
		curve->getPolyLine( t0, t1, polyline );
		INC_COUNTER(counter_curve);
	}else{
		// -> get common surface
		SurfaceParametric* surface0 = points[0]->getLocalSurface();
		SurfaceParametric* surface1 = points[1]->getLocalSurface();
		const DPoint3d& dpt0 = points[0]->getCoordinates();
		const DPoint3d& dpt1 = points[1]->getCoordinates();
		DPoint2d pt_2d_0, pt_2d_1;
		SurfaceParametric* surface = nullptr;
		if( surface0 != nullptr ){
			if( surface0 == surface1 ) {
				surface = surface0;
				pt_2d_0 = surface->getParameters( dpt0 );
				pt_2d_1 = surface->getParameters( dpt1 );
				INC_COUNTER(counter_same_surface);
			}else{
				mc.countMetricAtPoint( dpt1 );
				if( surface0->withinDomain(mc, dpt1, points[1]->getBaseNormal(), pt_2d_1 ) ) {
					surface = surface0;
					pt_2d_0 = surface->getParameters( dpt0 );
					INC_COUNTER(counter_surface_applicable);
				}
			}
		} else if (surface1 != nullptr ){
			mc.countMetricAtPoint( dpt0 );
			if( surface1->withinDomain(mc, dpt0, points[0]->getBaseNormal(), pt_2d_0 ) ) {
				surface = surface1;
				pt_2d_1 = surface->getParameters( dpt1 );
				INC_COUNTER(counter_surface_applicable);
			}
		}
		if( ! surface ) {
			// gather candidates
			assert( getFaceCount() > 0 );
			MeshFace* e_face = getFaceAt(0);
			assert( e_face != nullptr );
			assert( e_face->availableTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN ) );
			int domain_tag = e_face->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );
			DataVector<SurfaceParametric*> surf_candidates( 100 );
			points[0]->gatherLocalSurfaceCandidates( mc, domain_tag, surf_candidates );
			points[1]->gatherLocalSurfaceCandidates( mc, domain_tag, surf_candidates );

			// check
			int sc_ct = surf_candidates.countInt();
			if( sc_ct == 0 ){
				INC_COUNTER(counter_no_candidate);
				return getPoint(t);
			}

			mc.countMetricAtPoints( points[0], points[1] );
			const DVector3d& dvn0 = points[0]->getBaseNormal();
			const DVector3d& dvn1 = points[1]->getBaseNormal();

			if( sc_ct == 1 ) { // if only one candidate, make it simpler
				surface = surf_candidates[0];
				if( !surface->withinDomain(mc, dpt0, dvn0, pt_2d_0) ||
					!surface->withinDomain(mc, dpt1, dvn1, pt_2d_1)) 
				{
					INC_COUNTER(counter_one_candidate_fail);
					surface = nullptr;
				}else
					INC_COUNTER(counter_one_candidate_ok);
			}else{
				DPoint2d best_pt_2d_0, best_pt_2d_1;
				double best_value, q0, q1;
				for(int i = 0; i < sc_ct; i++){
					if( surf_candidates[i]->withinDomain(mc, dpt0, dvn0, pt_2d_0, &q0 ) &&  
						surf_candidates[i]->withinDomain(mc, dpt1, dvn1, pt_2d_1, &q1 ) )
					{
						double q_min = std::min( q0, q1 );
						if( !surface || q_min > best_value ) {
							surface = surf_candidates[i];
							best_pt_2d_0 = pt_2d_0;
							best_pt_2d_1 = pt_2d_1;
							best_value = q_min;
						}						
					}
				}
				if( surface ) {
					pt_2d_0 = best_pt_2d_0;
					pt_2d_1 = best_pt_2d_1;
					INC_COUNTER(counter_best_candidate_ok);
				}else
					INC_COUNTER(counter_best_candidate_fail);
			}
		}

		if( ! surface ) return getPoint(t);
		// -> get point params
		surface->getPolyLine( polyline, pt_2d_0, pt_2d_1 );
	}
	// from polyline -> count point
	double total_len = 0.0;
	int poly_count = polyline.countInt();
	DataVector<double> lengths( poly_count );
	for(int i = 1; i < poly_count; i++) {
		double len = polyline[i].distance( polyline[i-1] );
		lengths.add(len);
		total_len += len;
	}
	double req_len = t * total_len;
	int lcount = poly_count-1;
	int k = 0;
	while( k < lcount && req_len > lengths[k] ) req_len -= lengths[k++];
	if( k < lcount ) {
		DPoint3d result( polyline[k], polyline[k+1], req_len / lengths[k] );
		if( false ){
			MeshViewSet* set = new MeshViewSet;
			for(int i = 1; i < poly_count; i++)
				set->addEdge( polyline[i-1], polyline[i], 0 );
			set->addLabel( result, "mid" );
			set->addPoint( points[0] );
			set->addPoint( points[1] );

#ifdef COUNT_GET_POINT_ON_SURFACE_CASES
			set->addInfo("counter_curve", counter_curve );
			set->addInfo("counter_same_surface", counter_same_surface );
			set->addInfo("counter_surface_applicable", counter_surface_applicable );
			set->addInfo("counter_one_candidate_ok", counter_one_candidate_ok );
			set->addInfo("counter_one_candidate_fail", counter_one_candidate_fail );
			set->addInfo("counter_best_candidate_ok", counter_best_candidate_ok );
			set->addInfo("counter_best_candidate_fail", counter_best_candidate_fail );
			set->addInfo("counter_no_candidate", counter_no_candidate );
#endif
			set->addInfo("poly_count", poly_count );

			SHOW_MESH( "ME3d::getPointOnSurface", set );
		}
		return result;
	}else
		return polyline.last();
}
*/