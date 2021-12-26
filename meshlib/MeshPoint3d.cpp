// MeshPoint3d.cpp: implementation of the MeshPoint3d class.
//
//////////////////////////////////////////////////////////////////////


#include "MeshPoint2d.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshFace.h"
#include "MeshViewSet.h"
#include "DataHashTable.h"
#include "Metric3dContext.h"
#include "SurfaceParametric.h"
#include "Curve3dParametric.h"
#include "ControlSpace3dAdaptive.h"
#include "DataStatistics.h"
#include "MeshLog.h"

//////////////////////////////////////////////////////////////////////
MeshPoint3d::MeshPoint3d(double x, double y, double z) :
	coord(x,y,z), active(0), local_surface(nullptr), local_curve(nullptr) /* base_coord(x,y,z) */ { }

//////////////////////////////////////////////////////////////////////
MeshPoint3d::MeshPoint3d(const DPoint3d& pt) :
	coord(pt), active(0), local_surface(nullptr), local_curve(nullptr) /* base_coord(pt) */ { }

//////////////////////////////////////////////////////////////////////
MeshPoint3d::MeshPoint3d(const MeshPoint3d& point) : 
	coord(point.coord), active(0), local_surface(nullptr), local_curve(nullptr) /* base_coord(point.coord) */ { }

/// Surface-point constructor
MeshPoint3d::MeshPoint3d(SurfaceConstPtr surface, const DPoint2d& param, double aq, const DPoint3d * point3d)
	 : coord( point3d?*point3d:surface->getPoint(param)), active(0), 
		local_surface( new SurfaceData(surface, param, aq, coord) ), local_curve(nullptr) 
{ 
	AQ_ASSERT( aq );
}

/// Curve-point constructor
MeshPoint3d::MeshPoint3d(Curve3dConstPtr curve, const double& t, const DPoint3d * point3d)
	 : coord( point3d?*point3d:curve->getPoint(t)), active(0), 
		local_surface( new SurfaceData() ), local_curve( new CurveData( curve, t, coord) ) { }


MeshPoint3d::~MeshPoint3d()
{
//	edge_count = 0;	//aby unikn¹æ zbêdnego zagnie¿d¿onego usuwania odnoœników
	while(edges.notEmpty())
		delete edges[0];	// Destructor MeshEdge3d will take care of removing this association
	
	if(local_curve != nullptr) delete local_curve;
	if(local_surface != nullptr) delete local_surface;
}

void MeshPoint3d::preDeleteAll()
{
	for(int i = 0; i < edges.countInt(); i++){
		MeshEdge3d* edge = edges[i];
		if(edge->removePointLink(this)) delete edge;
	}
	edges.clear();
}

void MeshPoint3d::copyDataFrom(const MeshPoint3d* point)
{
	coord = point->coord;
	border = point->border;
}

MeshEdge3d* MeshPoint3d::getEdgeToPoint(const MeshPoint3d *point) const
{
	for(int i = 0; i < edges.countInt(); i++)
		if(edges[i]->getOtherPoint(this) == point)
			return edges[i];

	return nullptr;
}

void MeshPoint3d::addEdgeLink(MeshEdge3d *edge)
{
	assert(!edges.contains(edge));
	edges.add(edge);
}

MeshFace* MeshPoint3d::getFaceToPoints(const MeshPoint3d *point1, const MeshPoint3d *point2) const
{
	MeshEdge3d* edge = getEdgeToPoint(point1);
	if(!edge) return nullptr;
	int fct = edge->getFaceCount();
	for(int i = 0; i < fct; i++){
		MeshFace* face = edge->getFaceAt(i);
		if(face->incidentToPoint(point2)) return face;
	}
	return nullptr;
}

std::shared_ptr<MeshViewPointData> MeshPoint3d::getViewData() const
{
	char b = 0;
	if(isBorder(TagBorder::FIXED | TagBorder::CORNER)) b = 3;
	else if(isBorder(TagBorder::RIDGE)) b = 2;
	else if(isBorder()) b = 1;
	return std::make_shared<MeshViewPointData>(coord, b, 
		getIntTag(TagExtended::TAG_ID, index), true);
}

const DMPoint3d MeshPoint3d::getMetricCoordinates(const Metric3dContext& mc)
{
	return mc.transformRStoMS(coord);
}

bool MeshPoint3d::adjacentBlocks(DataVector<MeshBlock*> & blocks, DataHashTable<MeshBlock*> & visited_blocks, bool boundary_check) const
{
	blocks.prepare(2*edges.countInt());
	for(int j = 0; j < edges.countInt(); j++){
		int fct = edges[j]->getFaceCount();
		for(int i = 0; i < fct; i++){
			MeshFace* face = edges[j]->getFaceAt(i);
			for(int m = 0; m < 2; m++){
				MeshBlock* block = face->getBlock(m);
				if(boundary_check){
					assert(isBorder() || block);
					if(!isBorder() && !block) return false;
				}
				if(block && !visited_blocks.insert(block)) blocks.add(block);
			}
		}
	}
	return true;
} 

bool MeshPoint3d::adjacentBlocks(DataVector<MeshBlock*> & blocks, bool boundary_check) const
{
	blocks.prepare(2*edges.countInt());
	for(int j = 0; j < edges.countInt(); j++){
		int fct = edges[j]->getFaceCount();
		for(int i = 0; i < fct; i++){
			MeshFace* face = edges[j]->getFaceAt(i);
			for(int m = 0; m < 2; m++){
				MeshBlock* block = face->getBlock(m);
				if(boundary_check){
					assert(isBorder() || block);
					if(!isBorder() && !block) return false;
				}
				if(block) blocks.addIfNew(block);
			}
		}
	}
	return true;
}

bool MeshPoint3d::adjacentFaces(DataVector<MeshFace*> & faces) const
{
	faces.prepare(2*edges.countInt());
	for(int j = 0; j < edges.countInt(); j++){
		int fct = edges[j]->getFaceCount();
		for(int i = 0; i < fct; i++){
			faces.addIfNew(edges[j]->getFaceAt(i));
		}
	}
	return true;
}

bool MeshPoint3d::adjacentBlocksWithHash(DataVector<MeshBlock*> & blocks) const
{
	blocks.prepare(2*edges.countInt());
	DataHashTable<MeshBlock*> hash((unsigned int)(2*edges.countInt()), nullptr);

	for(int j = 0; j < edges.countInt(); j++){
		int fct = edges[j]->getFaceCount();
		for(int i = 0; i < fct; i++){
			MeshFace* face = edges[j]->getFaceAt(i);
			for(int m = 0; m < 2; m++){
				MeshBlock* block = face->getBlock(m);
				assert(isBorder() || block);
				if(!isBorder() && !block) return false;
				if(block && hash.insert(block)) blocks.add(block);
			}
		}
	}
	return true;
}

/// Returns the number of adjacent boundary edges
int MeshPoint3d::getBorderEdgesCount() const
{
	int count = 0;
	for(int i = 0; i < edges.countInt(); i++)
		if(edges[i]->isBorder()) count++;
	return count;
}

/*
/// Sets coordinates of the point after fitting to its ascribed local surface/curve
bool MeshPoint3d::setToLocalShape(Metric3dContext & mc)
{
	return recalculateForLocalShape( mc );
}

/// Returns coordinates of the point after fitting to its ascribed local surface/curve
DPoint3d MeshPoint3d::getCoordinatesWithLocalShape(Metric3dContext & mc, const DPoint3d& new_pt)
{
	DPoint3d pt = new_pt;
	recalculateForLocalShape(mc, pt);
	return pt;
}
*/

// Tries to move point to fit its ascribed local surface/curve
bool MeshPoint3d::moveToLocalShape(Metric3dContext& mc)
{
	assert(false); // is it really used?
	if( hasLocalCurve() )
		return tryMovingPoint( mc, getLocalCurve()->getPoint( getLocalCurveParam() ) );
	else if( hasLocalSurface() )
		return tryMovingPoint( mc, getLocalSurface()->getPoint(	getLocalSurfaceParam() ) );
	else 
		return false;
}

/// Tries to move point to new coordinates (gradually, checking for inverted faces)
bool MeshPoint3d::tryMovingPoint(Metric3dContext& mc, const DPoint3d& new_pt)
{
	//assert(false); // is it really used?
	//if(point->isBorder()) return false;

	int rank = getRank();
	if(rank < 3) return false;

	const DPoint3d old_pt = getCoordinates();
	const DVector3d vector = new_pt - old_pt;

	// gather adjacent faces
	DataVector<MeshFace*> faces;
	if(!adjacentFaces(faces)) return false;

	size_t fct = faces.countInt();
	double factor = 1.0;
	for(int k = 0; k < 10; k++){	// 10 tries
		setCoordinates(old_pt + vector * factor);
		mc.countMetricAtPoint( coord, true );
		bool valid = true;
		// Check normal vectors for faces
		for(size_t i = 0; valid && i < fct; i++){
			valid = faces[i]->validDirect(mc);
		}
		if(valid){
//			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Normal check ok, min sc = " << min_sc);
			return true;
		}else factor *= 0.5;
	}
	// 10x fault -> cancel movement
	setCoordinates(old_pt);
	return false;
}

/// Sets the coordinates for this point
void MeshPoint3d::setCoordinates( SurfaceConstPtr surface, const DPoint2d& param, 
								 double aq, const DPoint3d * point3d)
{
	AQ_ASSERT( aq );
	coord = point3d ? *point3d : surface->getPoint( param );
	local_surface->setSurfaceParam( surface, param, aq, coord );
}

/// Sets the coordinates for this point
void MeshPoint3d::setCoordinates( Curve3dConstPtr curve, const double& t, const DPoint3d * point3d)
{
	coord = point3d ? *point3d : curve->getPoint( t );
	local_curve->setCurveParam( curve, t, coord );
}

/// Recheck best surface
//bool MeshPoint3d::SurfaceData::recheckBestSurface( const DPoint3d& pt )
//{
//	if( best_surface != nullptr) {
//		if( best_surface->withinDomain( pt, best_surface_param, &best_surface_param ) )
//			return true;
//		else {
//			if(surface_params){
//				unsigned int i;
//				if( surface_params->contains( best_surface, i ) ) 
//					surface_params->slotValue( i ).counter = -param_counter;
//			}
//			best_surface = nullptr; // apparently, after change of coordinates, it is no longer valid
//		}
//	}
//
//	if( surface_params == nullptr ) return false;
//
//	unsigned int hcount = surface_params->countInt();
//	for(unsigned int i = 0; i < hcount; i++) {
//		SurfaceConstPtr surface = surface_params->slotKey( i );
//		if( surface != nullptr){
//			auto& pwt = surface_params->slotValue( i );
//			if( (pwt.counter != -param_counter) && surface->withinDomain( pt, pwt.param, &pwt.param ) ) {
//				pwt.counter = param_counter;
//				best_surface = surface;
//				best_surface_param = pwt.param;
//				return true;
//			}else
//				pwt.counter = -param_counter; // checked and invalid
//		}
//	}
//
//	return false;
//}

/// Tries to move point to new coordinates (gradually, checking for inverted faces)
bool MeshPoint3d::tryMovingPoint( /*Metric3dContext& mc, */ SurfaceConstPtr surface, 
		const DPoint2d& new_pt, const DataVector< MeshFace* > & faces )
{
	assert( faces.notEmpty() );
	int rank = getRank();
	if(rank < 3) return false;

	assert( local_surface != nullptr );
	assert( surface != nullptr );
	DPoint3d old_pt = coord;
	ParamAndQuality old_pt_2d = local_surface->getSurfaceParamQuality( surface, coord );
	assert( old_pt_2d.isValid() );
	const DVector2d vt = new_pt - old_pt_2d.param;

	//static int counter = 0;
	//counter++;

	// cancel if vt.length() < eps ?
	if( vt.length2() < VERY_SMALL_NUMBER ) return true;

	//mc.countMetricAtPoint( coord, true );

	size_t fct = faces.countInt();
	double factor = 1.0;
	for(int k = 0; k < 10; k++){	// 10 attempts
		DPoint2d new_pt_2d = old_pt_2d.param + vt * factor;
		setCoordinates( surface, new_pt_2d, AQ_UNKNOWN );

		bool valid = true;
		for(size_t i = 0; valid && i < fct; i++){
			valid = faces[i]->valid( surface );
		}
		if(valid) return true;
		else factor *= 0.5;
	}
	// 10x fault -> cancel movement
	setCoordinates( surface, old_pt_2d.param, old_pt_2d.quality, &old_pt );
	return false;
}

/// Tries to move point to new coordinates (gradually, checking for inverted faces)
bool MeshPoint3d::tryMovingPoint( /*Metric3dContext& mc, */ Curve3dConstPtr curve, const double& t)
{
	int rank = getRank();
	if(rank < 3) return false;

	// gather adjacent faces
	DataVector<MeshFace*> faces;
	if(!adjacentFaces(faces)) return false;

	assert( local_curve != nullptr );
	assert( curve != nullptr );
	DPoint3d old_pt = coord;
	double old_t = getLocalCurveParam( curve );
	const double dt = t - old_t;

	// cancel if dt < eps ?
	if( dt < VERY_SMALL_NUMBER ) return true;

	size_t fct = faces.countInt();
	double factor = 1.0;
	for(int k = 0; k < 10; k++){	// 10 attempts
		double new_t = old_t + dt * factor;
		setCoordinates( curve, new_t );

		bool valid = true;
		for(size_t i = 0; valid && i < fct; i++){
			valid = faces[i]->valid( faces[i]->getOptLocalSurface() );
		}
		if(valid) return true;
		else factor *= 0.5;
	}
	// 10x fault -> cancel movement
	setCoordinates( curve, old_t, &old_pt );
	return false;
}

/*

//#define SHOW_SURFACE_SELECTION
#define SHOW_NO_VALID_LOCAL_SURFACE

/// Recalculates coordinates of this point according to local shape(s) - with possible change of active shape!
bool MeshPoint3d::recalculateForLocalShape( Metric3dContext & mc )
{
	if(isBorder()){
		if( isBorder(TagBorder::CORNER) ) return false;

		if(local_curve && local_curve->withinDomain( local_curve->getParameter(coord) ) )
			return true;

		double t;
		for(int i = 0; i < getRank(); i++){
			MeshEdge3d* other_edge = getEdge(i);
			if(!other_edge->isBorder()) continue;
			Curve3dConstPtr * edge_curve = other_edge->getLocalCurve();
			assert(edge_curve);
			if( (edge_curve != local_curve) && edge_curve->withinDomain( t = edge_curve->getParameter(coord) ) ) {
				local_curve = edge_curve;
				coord = edge_curve->getPoint(t);
				return true;
			}
			MeshPoint3d* other_point = other_edge->getOtherPoint(this);
			if( other_point->isBorder(TagBorder::CORNER) ) continue;
			Curve3dConstPtr * other_curve = other_point->local_curve;
			if( (other_curve != local_curve) && other_curve->withinDomain( t = other_curve->getParameter(coord) ) ) {
				local_curve = other_curve;
				local_surface = other_point->local_surface;
				coord = other_curve->getPoint(t);
				return true;
			}
		}
		// found nothing adequate...
		return false;
	}else{
		DPoint2d pt_2d;
		static const bool SELECT_FIRST_VALID = false;
		static const bool KEEP_INITIAL_VALID = true;
		double q;
		const DVector3d pnorm = getBaseNormal();

		if(KEEP_INITIAL_VALID && local_surface && local_surface->withinDomain( mc, coord, pnorm, pt_2d, &q) && (q > 0.3)){
			moveLocal(mc, local_surface->getPoint( pt_2d ) );
			return true;
		}

		// else choose and update

		assert( getRank() > 0 );
		MeshFace* pt_face = getEdge(0)->getFaceAt(0);
		assert( pt_face != nullptr );
		assert( pt_face->availableTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN ) );
		int domain_tag = pt_face->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );

		DataVector<SurfaceParametric*> surf_candidates( 100 );

		gatherLocalSurfaceCandidates( mc, domain_tag, surf_candidates );

		int sc_ct = surf_candidates.countInt();

		if( sc_ct == 1 ) { // if only one candidate, make it simpler
			if( surf_candidates[0]->withinDomain(mc, coord, pnorm, pt_2d ) ) {
				local_surface = surf_candidates[0];
				moveLocal( mc, local_surface->getPoint( pt_2d ) );
				return true;
			}else{
#ifdef SHOW_NO_VALID_LOCAL_SURFACE
				MeshViewSet* set = new MeshViewSet;
				surf_candidates[0]->drawDomain(set);
				set->addPoint( this , 1 );
				set->addInfo("local surf. count", 1 );
				SHOW_MESH(" MP3d::recalculateForLocalShape(mc) - single fail", set );
#endif				
				return false;
			}
		}

#ifdef SHOW_SURFACE_SELECTION
		DataVector<double> domain_q( sc_ct, 0.0 );
		int ct_in_domain = 0;
#endif

		SurfaceParametric* best_surface = nullptr;
		DPoint2d best_pt_2d;
		double best_value;
		mc.countMetricAtPoint( coord );
		for(int i = 0; i < sc_ct; i++){
			if( surf_candidates[i]->withinDomain(mc, coord, pnorm, pt_2d, &q ) ){
				if( SELECT_FIRST_VALID ) {
					local_surface = surf_candidates[i];
					coord = local_surface->getPoint( pt_2d );
					return true;
				}else{
					if( !best_surface || q > best_value ) {
						best_surface = surf_candidates[i];
						best_pt_2d = pt_2d;
						best_value = q;
					}

#ifdef SHOW_SURFACE_SELECTION
					domain_q[i] = q;
					ct_in_domain++;
#endif
				}
			}
		}

#ifdef SHOW_SURFACE_SELECTION
//		if(local_surface == nullptr) {
		if( false ) {
			//if( ct_in_domain < sc_ct ){ // show outside-domain
			//	MeshViewSet* set = new MeshViewSet;
			//	set->addPoint( this, 2 );
			//	DataVector<MeshFace*> faces;
			//	if( this->adjacentFaces( faces ) ){
			//		for(int i = 0; i < faces.countInt(); i++)
			//			set->addFaceWithEdges( faces[i] );
			//	}
			//	int k = 0;
			//	for(int i = 0; i < sc_ct; i++){
			//		if( surf_dist[i] < 0.0 ){
			//			surf_candidates[i]->drawDomain( set, k++ );
			//		}
			//	}
			//	set->addInfo("surfaces inside", ct_in_domain );
			//	set->addInfo("surfaces outside", sc_ct - ct_in_domain );
			//	SHOW_MESH( "surfaces outside domain", set );
			//}
			if( ct_in_domain > 0 ){ // show inside-domain
				MeshViewSet* set = new MeshViewSet;
				set->addPoint( this, 2 );
				DataVector<MeshFace*> faces;
				if( this->adjacentFaces( faces ) ){
					for(int i = 0; i < faces.countInt(); i++)
						set->addFaceWithEdges( faces[i] );
				}
				int k = 0;
				for(int i = 0; i < sc_ct; i++){
					if( surf_dist[i] >= 0.0 ){
						surf_candidates[i]->drawDomain( set, (surf_candidates[i] == best_surface) ? 1 : 0 );
						string s = "quality";
						if( surf_candidates[i] == best_surface ) s+= "*";
						set->addInfo(s, to_string(domain_q[i]) );
					}
				}
				if( best_surface ) {
					set->addPoint( best_surface->getPoint( best_pt_2d ), 3 );
				}
				set->addInfo("surfaces inside", ct_in_domain );
				SHOW_MESH( "surfaces inside domain", set );
			}
		}
#endif

		if(best_surface){
			local_surface = best_surface;
			moveLocal( mc, best_surface->getPoint( best_pt_2d ) );

			//assert( best_surface && best_surface->withinDomain( coord, pt_2d ) );

			return true;
		}else{
#ifdef SHOW_NO_VALID_LOCAL_SURFACE
				MeshViewSet* set = new MeshViewSet;
				for(int i = 0; i < surf_candidates.countInt(); i++)
					surf_candidates[i]->drawDomain( set, i );
				set->addPoint( this , 1 );
				set->addInfo("local surf. count", surf_candidates.countInt() );
				SHOW_MESH(" MP3d::recalculateForLocalShape(mc) - multi fail", set );
				//---
				//for(int i = 0; i < surf_candidates.countInt(); i++) {
				//	set = new MeshViewSet;
				//	SurfaceParametric* i_surf = surf_candidates[i];
				//	i_surf->drawDomain( set, i );
				//	set->addPoint( this , 1 );
				//	set->addLabel( i_surf->getPoint( i_surf->getParameters( coord ) ), "surf-pt" );
				//	set->addInfo("local surf. #", i );
				//	SHOW_MESH(" MP3d::multi-fail i-th surf-candidate", set );
				//	bool res = i_surf->withinDomain(mc, coord, pnorm, pt_2d, &q );
				//	LOG4CPLUS_INFO(MeshLog::logger_console, "local-surf-q", q);
				//}
#endif	
			return false;
		}
	}

	return true;
}

/// Recalculates coordinates of the given point according to local shape(s) - with possible change of active shape!
bool MeshPoint3d::recalculateForLocalShape(Metric3dContext & mc, DPoint3d & pt)
{
	if( recalculateForLocalShape( mc ) ){
		// local shape is available and correct (within domain for "this" point ...
		// ... now try to check if it is correct for the other point "pt" as well...
		if( isBorder() ) {
			assert( local_curve != nullptr );
			double t;
			if(local_curve->withinDomain(t = local_curve->getParameter(pt) ) ) {
				pt = local_curve->getPoint(t);
				return true;
			}
			// try to find other curve ?

			// ... gather candidate curves
			DataVector<Curve3dConstPtr*> curv_candidates( 2 * getRank() );
			for(int i = 0; i < getRank(); i++){
				MeshEdge3d* other_edge = getEdge(i);
				if(!other_edge->isBorder()) continue;
				Curve3dConstPtr * edge_curve = other_edge->getLocalCurve();
				assert(edge_curve);
				if(edge_curve != local_curve) curv_candidates.addIfNew( edge_curve );

				MeshPoint3d* other_point = other_edge->getOtherPoint(this);
				if( other_point->isBorder(TagBorder::CORNER) ) continue;
				Curve3dConstPtr * other_curve = other_point->local_curve;
				if( (other_curve != local_curve) ) curv_candidates.addIfNew( other_curve );
			}

			// ... check if any contains both this and "pt"
			Curve3dConstPtr* valid_pt_curve = nullptr;
			double valid_pt_curve_t;
			for(int i = 0; i < curv_candidates.countInt(); i++){
				if( curv_candidates[i]->withinDomain( t = curv_candidates[i]->getParameter( pt ) ) ) {
					if( !valid_pt_curve ) {
						valid_pt_curve = curv_candidates[i];
						valid_pt_curve_t = t;
					}
					if( curv_candidates[i]->withinDomain( curv_candidates[i]->getParameter( coord ) ) ) {
						pt = curv_candidates[i]->getPoint( t );
						return true; // found curve valid fot both ...
					}
				}
			}

			// ... else if any contains "pt"
			if( valid_pt_curve ){
				pt = valid_pt_curve->getPoint( valid_pt_curve_t );
				return true;
			}

			// ... else failed
			return false;
		}else {
			assert( local_surface != nullptr );
			DPoint2d pt_2d, pt_2d_c;
			const DVector3d pnorm = getBaseNormal();
			if(local_surface->withinDomain(mc, pt, pnorm, pt_2d ) ) {
				pt = local_surface->getPoint( pt_2d );
				return true;
			}
			// try to find other surface ?

			// ... gather candidate surfaces
			assert( getRank() > 0 );
			MeshFace* pt_face = getEdge(0)->getFaceAt(0);
			assert( pt_face != nullptr );
			assert( pt_face->availableTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN ) );
			int domain_tag = pt_face->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );
			DataVector<SurfaceParametric*> surf_candidates( 100 );

			gatherLocalSurfaceCandidates( mc, domain_tag, surf_candidates, &pt );
		
			// ... check if any contains both this and "pt"
			SurfaceParametric* valid_pt_surface = nullptr;
			DPoint2d valid_pt_surface_param;
			for(int i = 0; i < surf_candidates.countInt(); i++){
				if( surf_candidates[i]->withinDomain(mc, pt, pnorm, pt_2d ) ) {
					if( !valid_pt_surface ) {
						valid_pt_surface = surf_candidates[i];
						valid_pt_surface_param = pt_2d;
					}
					if( surf_candidates[i]->withinDomain( mc, coord, pnorm, pt_2d_c ) ) {
						pt = surf_candidates[i]->getPoint( pt_2d );
						return true; // found surface valid fot both ...
					}
				}
			}

			// ... else if any contains "pt"
			if( valid_pt_surface ){
				pt = valid_pt_surface->getPoint( valid_pt_surface_param );
				return true;
			}

		// ... else failed
			return false;
		}
	} else
		return false;
}
*/

/// Gather local surface candidates for this point
//bool MeshPoint3d::gatherLocalSurfaceCandidates( Metric3dContext & mc, int domain_tag, 
//			DataVector< SurfaceConstPtr > & scandidates, DPoint3d * dpt ) const
//{
//	SurfaceParametricSet* cs_set = nullptr;
//	CS3dPtr cs = (CS3dPtr)mc.getControlSpace();
//	if(cs && cs->isAdaptive()) cs_set = cs->getLocalSurfaceSetAtPoint( dpt ? *dpt : coord );
//
//	int max_surf_ct = 1 + getRank() + ( cs_set ? cs_set->countInt() : 0 );
//	scandidates.prepare( max_surf_ct );
//
//	bool any_new = false;
//
//	// - acs-surfaces
//	if( cs_set ){ 
//		for(int i = 0; i < cs_set->countInt(); i++){
//			SurfaceConstPtr surf = cs_set->getSurface(i);
//			if( (surf != local_surface->best_surface) && surf->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, domain_tag ) )
//				any_new |= scandidates.addIfNew( surf );
//		}
//	}
//
//	// - local-surface
//	if( local_surface ) any_new |= scandidates.addIfNew( local_surface->best_surface );
//
//	// TODO -> ??? plus surface_params hash-table from vicinity points?
//	// - vicinity-surfaces
//	for(int i = 0; i < getRank(); i++){ 
//		SurfaceConstPtr other_surface = getEdge(i)->getOtherPoint(this)->local_surface->best_surface;
//		if(other_surface && (other_surface != local_surface->best_surface) && 
//				other_surface->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, domain_tag ) ) 
//			any_new |= scandidates.addIfNew( other_surface );
//	}
//
//	return any_new;
//}

/// Check if point is planar with respect to adjacent faces (according to normals)
bool MeshPoint3d::planar(const Metric3dContext & mc) const
{
	double skip_len2 = mc.getMinSkipLen2();
	if(skip_len2 > 0.0){
		int ect = getRank();
		for(int i = 0; i < ect; i++)
			if( getEdge(i)->getLength2() < skip_len2) return true; // no point in checking, will be removed later anyway...
	}

	assert( !base_normal.isZero() );

	DataVector<MeshFace*> pfaces(edges.countInt());
	bool result = adjacentFaces(pfaces);
	assert( result );

	static double min_sp = 0.9;

	for(int i = 0; i < pfaces.countInt(); i++){
		if(pfaces[i]->hasBaseNormal()){
			double sp = base_normal.scalarProduct( pfaces[i]->getBaseNormal() );
			if(sp < min_sp){
				LOG4CPLUS_INFO(MeshLog::logger_console, "New min sp for MP3d::planar: " << sp);
				//if(true){
				//	MeshViewSet* set = new MeshViewSet();
				//	set->addPoint(this);
				//	for(int j = 0; j < pfaces.countInt(); j++){
				//		set->addFace(pfaces[j], (j==i) ? 1 : 0);
				//	}
				//	SHOW_MESH("MeshPoint3d::planar", set);
				//}
				min_sp = sp;
			}
			if(sp < 0.75) return false;
		}
	}

	return true;
}

/*
/// optimize selection of local surface 
bool MeshPoint3d::selectBestLocalSurface( Metric3dContext& mc, int local_surface_tag)
{
	if( isBorder() ) return false;
	if( local_surface_tag < 0 ) {
		assert( getRank() > 0 );
		MeshFace* pt_face = getEdge(0)->getFaceAt(0);
		assert( pt_face != nullptr );
		assert( pt_face->availableTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN ) );
		local_surface_tag = pt_face->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );
	}
	if( local_surface->best_surface_quality == 1.0 ) return false; // not really possible to improve

	DataVector<SurfaceConstPtr> surf_candidates;
	if( !gatherLocalSurfaceCandidates( mc, local_surface_tag, surf_candidates ) ) return false;

	DPoint2d pt_2d;
	double quality;
	SurfaceConstPtr old_surf = local_surface->best_surface;
	for(int i = 0; i < surf_candidates.countInt(); i++) {
		SurfaceConstPtr csurf = surf_candidates[i];
		if( csurf == old_surf ) continue;
		if( csurf->withinDomain( mc, this, &pt_2d, &quality ) && (quality > local_surface->best_surface_quality) ){
			setLocalSurface( csurf, pt_2d, quality );
		}
	}

	return old_surf != local_surface->best_surface;
}
*/

/// Move point to new coordinates, with constraint of the maximum metric distance from base coordinates
/*
bool MeshPoint3d::moveLocal(Metric3dContext& mc, const DPoint3d& new_coord, double max_dist, bool max_cut)
{
	mc.countMetricAtPoint( base_coord, true );
	double dist = mc.transformRStoMS( new_coord - base_coord ).length();
	if(dist > max_dist) {

#ifdef _DEBUG
		if( false ){
//		if( dist > 1.5 * max_dist ){
			MeshViewSet* set = new MeshViewSet;
			set->addPoint( this, 2 );
			DataVector<MeshFace*> faces;
			if( this->adjacentFaces( faces ) ){
				for(int i = 0; i < faces.countInt(); i++)
					set->addFaceWithEdges( faces[i] );
			}
			if( local_surface ) {
				local_surface->drawDomain( set, 1 );
				
				DPoint2d pt_2d;
				double q;
				if( local_surface->withinDomain( mc, new_coord, getBaseNormal(), pt_2d, &q ) ) {
					set->addInfo("local surface quality", q );
				}else{
					set->addInfo("local surface quality", "invalid");
				}
			}

			set->addLabel( base_coord, "base" );
			set->addLabel( new_coord, "new" );
			set->addLabel( DPoint3d( base_coord, new_coord, max_dist / dist ), "max" );

			set->addInfo("metric dist new/max", to_string(dist) + " / " + to_string(max_dist) );
			SHOW_MESH( "moveLocal too far ...", set );
		}
#endif
		if( max_cut ) {
			coord = DPoint3d( base_coord, new_coord, max_dist / dist );
		}else return false;
	}else
		coord = new_coord;
	return true;
}
*/

/// Returns arbitrary face adjacent to this point
MeshFace* MeshPoint3d::anyAdjacentFace() const
{
	for(int i = 0; i < edges.countInt(); i++) {
		MeshEdge3d* edge = edges[i];
		if( edge->getFaceCount() == 0 ) continue;
		return edge->getFaceAt(0);
	}

	assert(false);
	return nullptr;
}

/// Get local vicinity set of mesh points (with this point)
DataVector<MeshPoint3d*> MeshPoint3d::localVicinityPoints()
{
	size_t rank = edges.countInt();
	DataVector<MeshPoint3d*> vicinity(rank+1);
	vicinity.add(this);
	for(size_t i = 0; i < rank; i++)
		vicinity.add( edges[i]->getOtherPoint( this ) );
	return vicinity;
}

/// Get local vicinity set of mesh points (with this point)
DataVector<MeshPoint3d*> MeshPoint3d::localVicinityPointsBorderEdges( bool border )
{
	assert( isBorder() == border );
	size_t rank = edges.countInt();
	DataVector<MeshPoint3d*> vicinity(rank+1);
	vicinity.add(this);
	for(size_t i = 0; i < rank; i++) {
		if( edges[i]->isBorder() == border) 
			vicinity.add( edges[i]->getOtherPoint( this ) );
	}
	return vicinity;
}

/// Standard constructor
MeshPoint3d::SurfaceData::SurfaceData( SurfaceConstPtr surface, const DPoint2d& param, 
									  double aq, const DPoint3d& pt )
{ 
	assert( surface != nullptr );
	assert( surface->isFixed() );
	AQ_ASSERT( aq );
	//DPoint2d check_param = param;
	//bool check_valid = best_surface->withinDomain( pt, check_param );
	//assert( check_valid == true );
	// ParamData* ref = &
	surface_params.append( ParamData(surface, pt, param, aq ) ); 
	//hsparams.insert( surface, ref );
}

MeshPoint3d::SurfaceData::ParamData::ParamData(SurfaceConstPtr _surface, 
				const DPoint3d& _ref_pt,	const DPoint2d& _param, double aq )
	: surface(_surface), ref_pt(_ref_pt), param(_param), quality(aq) 
{
	assert( surface != nullptr );
	assert( surface->isFixed() );
	AQ_ASSERT( aq );

	if( quality == AQ_UNKNOWN ) {
		param = surface->getDomainMiddleParam();
		quality = surface->withinDomainQuality( ref_pt, param );
	}
}

double MeshPoint3d::SurfaceData::ParamData::updateGetQuality( const DPoint3d& pt )
{
	return quality = surface->withinDomainQuality( ref_pt = pt, param ); 
}

double MeshPoint3d::SurfaceData::ParamData::updateIfNeededGetQuality( const DPoint3d& pt )
{
	if( (ref_pt != pt) || (quality == AQ_UNKNOWN) ) {
		quality = surface->withinDomainQuality( ref_pt = pt, param ); 
	}else{
		//DPoint2d check_param = param;
		//bool check_valid = surface->withinDomain( pt, check_param );
		//if( valid != check_valid ) {
		//	MeshViewSet* set = new MeshViewSet;
		//	set->addPoint(pt, 1, 0);
		//	surface->drawDomain( set );
		//	SHOW_MESH("MP3d::SD::PD::updateIfNeeded inconsistent...", set);
		//}
		//assert( valid == check_valid );
	}
	return quality;
}


MeshPoint3d::CurveData::ParamData::ParamData(Curve3dConstPtr _curve, const DPoint3d& _ref_pt, double _param )
		: curve(_curve), ref_pt(_ref_pt), param(_param), valid(true) 
{ 
	assert( curve->isFixed() ); 
}

MeshPoint3d::CurveData::ParamData::ParamData(Curve3dConstPtr _curve, const DPoint3d& _ref_pt )
	: curve(_curve), ref_pt(_ref_pt), param( curve->getMiddleParam() )
{
	assert( curve->isFixed() );
	valid = curve->withinDomain( _ref_pt, param );
}


void MeshPoint3d::CurveData::ParamData::set( const double& _param, const DPoint3d& pt) 
{ 
	param = _param; 
	ref_pt = pt; 
	valid = true; 
}

bool MeshPoint3d::CurveData::ParamData::update( const DPoint3d& pt )
{
	return valid = curve->withinDomain( ref_pt = pt, param ); 
}

bool MeshPoint3d::CurveData::ParamData::updateIfNeeded( const DPoint3d& pt )
{
	if( ref_pt != pt ) {
		valid = curve->withinDomain( ref_pt = pt, param ); 
	}
	return valid;
}

MeshPoint3d::CurveData::CurveData( Curve3dConstPtr local_curve, const double& t, const DPoint3d& pt)
{ 
	assert( local_curve != nullptr ); 
}

/// Get parameters for this curve
double MeshPoint3d::CurveData::getCurveParam( Curve3dConstPtr curve, const DPoint3d& pt )
{
	assert( curve != nullptr );
	ParamData* pd = findCurveParamData( curve );
	if( pd != nullptr ) {
		pd->updateIfNeeded( pt );
		return pd->param;
	}else if( curve->isFixed() ) {
		pd = &curve_params.append( ParamData( curve, pt ) );
		return pd->param;
	}else{
		return curve->getParameter( pt, curve->getMiddleParam() );
	}
}

MeshPoint3d::CurveData::ParamData* MeshPoint3d::CurveData::findCurveParamData( Curve3dConstPtr curve ) const
{
	return curve_params.findFirst( [curve] ( const ParamData& pd ) { return pd.curve == curve; } );
}

MeshPoint3d::SurfaceData::ParamData* MeshPoint3d::SurfaceData::findSurfaceParamData( SurfaceConstPtr surface ) const
{
	return surface_params.findFirst( [surface] ( const ParamData& pd ) { return pd.surface == surface; } );
}

/// Get parameters for this surface
DPoint2d MeshPoint3d::SurfaceData::getSurfaceParam( SurfaceConstPtr surface, const DPoint3d& pt )
{
	assert( surface != nullptr );
	ParamData* pd = findSurfaceParamData( surface );
	if( pd != nullptr ) {
		pd->updateIfNeededGetQuality( pt );
		return pd->param;
	}else if( surface->isFixed() ) {
		pd = &surface_params.append( ParamData( surface, pt ) );
		return pd->param;
	}else{
		return surface->getParameters( pt );
	}
}

/// Get parameters for this surface
ParamAndQuality MeshPoint3d::SurfaceData::getSurfaceParamQuality( SurfaceConstPtr surface, const DPoint3d& pt )
{
	assert( surface != nullptr );
	ParamData* pd = findSurfaceParamData( surface );
	if( pd != nullptr ) {
		pd->updateIfNeededGetQuality( pt );
		return pd->getParamAndQuality();
	}else if( surface->isFixed() ) {
		pd = &surface_params.append( ParamData( surface, pt ) );
		return pd->getParamAndQuality();
	}else{
		return ParamAndQuality( surface->getParameters( pt ), AQ_VALID_MAX );
	}
}

ParamAndQuality::ParamAndQuality( SurfaceConstPtr surface, const DPoint3d& pt )
{
	assert( surface != nullptr );
	quality = surface->withinDomainQuality( pt, param );
}

/// Get parameters for this surface
ParamAndQuality MeshPoint3d::SurfaceData::checkAndGetSurfaceParamQuality( SurfaceConstPtr surface, const DPoint3d& pt )
{
	assert( surface != nullptr );
	ParamData* pd = findSurfaceParamData( surface );
	if( pd != nullptr ) {
		pd->updateIfNeeded( pt );
		return pd->getParamAndQuality();
	}else if( surface->isFixed() ) {
		pd = &surface_params.append( ParamData( surface, pt ) );
		return pd->getParamAndQuality();
	}else{
		return ParamAndQuality( surface, pt );
	}
}

bool MeshPoint3d::SurfaceData::checkAndGetSurfaceParam( SurfaceConstPtr surface, DPoint2d& param, const DPoint3d& pt )
{
	assert( surface != nullptr );
	ParamData* pd = findSurfaceParamData( surface );
	if( pd != nullptr ) {
		pd->updateIfNeeded( pt );
		param = pd->param;
		return pd->isValid();
	}else if( surface->isFixed() ) {
		pd = &surface_params.append( ParamData( surface, pt ) );
		param = pd->param;
		return pd->isValid();
	}else{
		return surface->withinDomain( pt, param );
	}
}

/// Get parameters for this surface
bool MeshPoint3d::CurveData::checkAndGetCurveParam( Curve3dConstPtr curve, const DPoint3d& pt, double& t )
{
	assert( curve != nullptr );
	ParamData* pd = findCurveParamData( curve );
	if( pd != nullptr ) {
		pd->updateIfNeeded( pt );
		t = pd->param;
		return pd->valid;
	}else if( curve->isFixed() ) {
		pd = &curve_params.append( ParamData( curve, pt ) );
		t = pd->param;
		return pd->valid;
	}else{
		return curve->withinDomain( pt, t );
	}
}

/// Changes local curve
void MeshPoint3d::setLocalCurve(Curve3dConstPtr curve, const double & t) { 
	assert( isBorder() );
	assert( curve->isFixed() );
	if( local_curve != nullptr )
		local_curve->setCurveParam(curve, t, coord);
	else {
		local_curve = new CurveData(curve, t, coord);
		if( local_surface == nullptr )
			local_surface = new SurfaceData();
	}
}

/// Changes local surface
void MeshPoint3d::setLocalSurface(SurfaceConstPtr surface, const DPoint2d& param, double aq ) 
{ 
	assert( surface->isFixed() );
	AQ_ASSERT( aq );
	if( local_surface != nullptr )
		local_surface->setSurfaceParam(surface, param, aq, coord);
	else
		local_surface = new SurfaceData(surface, param, aq, coord); 
}

/// Returns local param for surface
ParamAndQuality MeshPoint3d::getLocalSurfaceParamQuality( SurfaceConstPtr surface ) const { 
	if( local_surface != nullptr)
		return local_surface->getSurfaceParamQuality( surface, coord ); 
	else
		return ParamAndQuality( surface->getParameters( coord ), AQ_VALID_MAX );
}

DPoint2d MeshPoint3d::getLocalSurfaceParam( SurfaceConstPtr surface ) const { 
	if( local_surface != nullptr)
		return local_surface->getSurfaceParam( surface, coord ); 
	else
		return surface->getParameters( coord );
}

/// Returns local param for surface
DPoint2d MeshPoint3d::getLocalSurfaceParam( ) const { 
	assert( local_surface != nullptr); 
	assert( local_surface->surface_params.notEmpty() );
	return local_surface->getSurfaceParam( local_surface->getValidSurface(coord), coord ); 
}

/// Get parameters for this surface
bool MeshPoint3d::checkAndGetSurfaceParam( SurfaceConstPtr surface, DPoint2d & param )
{
	if( local_surface != nullptr)
		return local_surface->checkAndGetSurfaceParam( surface, param, coord );
	else
		return surface->withinDomainQuality( coord, param ) >= AQ_VALID_MIN;
}

/// Get parameters for this curve
bool MeshPoint3d::checkAndGetCurveParam(  Curve3dConstPtr curve, double& t )
{
	if( local_curve != nullptr)
		return local_curve->checkAndGetCurveParam( curve, coord, t );
	else 
		return curve->withinDomain( coord, t );
}

/// Returns local param for surface
double MeshPoint3d::getLocalCurveParam( Curve3dConstPtr curve ) const { 
	assert( local_curve != nullptr); 
	return local_curve->getCurveParam( curve, coord ); 
}

/// Returns local param for surface
double MeshPoint3d::getLocalCurveParam( ) const { 
	assert( local_curve != nullptr); 
	return local_curve->getCurveParam( local_curve->getValidCurve(coord), coord ); 
}

void MeshPoint3d::clearLocalShapes()
{
	if( local_surface != nullptr ) { delete local_surface; local_surface = nullptr; }
	if( local_curve != nullptr )   { delete local_curve;   local_curve = nullptr; }
}

/// Checks validity of local surface params
bool MeshPoint3d::checkLocalSurfaceParams() const
{
	return (local_surface == nullptr) || local_surface->checkSurfaceParams();
}

/// Checks validity of local surface params
bool MeshPoint3d::SurfaceData::checkSurfaceParams() const
{
	int invalid = 0;
	DPoint2d check_param;
	double max_quality_diff = 0.0;

	surface_params.forEach( [&] ( const ParamData& pd ) {
		assert( pd.surface->isFixed() );
		check_param = pd.param;
		double q = pd.surface->withinDomainQuality( pd.ref_pt, check_param );
		if( (q >= AQ_VALID_MIN) != pd.isValid() ){
			if( true ) {
				MeshViewSet* set = new MeshViewSet;
				pd.surface->drawDomain( set );
				set->addInfo("q", q);
				set->addInfo("pd.isValid", pd.isValid() ? "true":"false" );
				set->addPoint( pd.ref_pt, 1, 0 );
				SHOW_MESH( "domain check failed", set );
			}
			++invalid;
		}
		if( q >= AQ_VALID_MIN && pd.quality >= AQ_VALID_MIN ) {
			double dq = std::abs( q - pd.quality );
			if( dq > max_quality_diff ) { 
				max_quality_diff = dq;
			}
		}

	} );

	return invalid == 0;
}

/// Add local surface
void MeshPoint3d::SurfaceData::setSurfaceParam(SurfaceConstPtr surface, const DPoint2d& param, 
											   double aq, const DPoint3d& pt ) 
{ 
	assert( surface != nullptr );
	assert( surface->isFixed() );
	AQ_ASSERT( aq );

	ParamData* pd = findSurfaceParamData( surface );
	if( pd != nullptr ) {
		pd->set( param, aq, pt );
	}else{
		surface_params.append( ParamData( surface, pt, param, aq ) );
	}

	//DPoint2d check_param = param;
	//bool check_valid = surface->withinDomain( pt, check_param );
	//if( pd.valid != check_valid ) {
	//	MeshViewSet* set = new MeshViewSet;
	//	set->addPoint(pt, 1, 0);
	//	surface->drawDomain( set );
	//	SHOW_MESH("MP3d::SD::setSurfaceParam inconsistent...", set);
	//}
	//assert( pd.valid == check_valid );
}

/// Set/change local curve param
void MeshPoint3d::CurveData::setCurveParam( Curve3dConstPtr curve, const double& t, const DPoint3d& pt ) 
{ 
	assert( curve != nullptr );
	assert( curve->isFixed() );

	ParamData* pd = findCurveParamData( curve );
	if( pd != nullptr ) {
		pd->set( t, pt );
	}else{
		curve_params.append( ParamData( curve, pt, t ) );
	}
}

/// Get neighboring point
MeshPoint3d* MeshPoint3d::incidentPoint(int i ) const 
{ 
	return getEdge(i)->getOtherPoint(this); 
}

//double MeshPoint3d::checkMaxLocalSurfaceParamDist2( DataStatistics& stats ) const 
//{
//	return (local_surface != nullptr) ? local_surface->checkMaxLocalSurfaceParamDist2(stats, coord) : 0.0;
//}
//
//double MeshPoint3d::SurfaceData::checkMaxLocalSurfaceParamDist2( DataStatistics& stats, const DPoint3d& pt ) const 
//{
//	int counter = 0;
//	double max_dist2 = -1.0;
//
//	assert( best_surface != nullptr );
//
//	if( best_surface_param.isValid( pt) ) {
//		counter++;
//		max_dist2 = best_surface_param.param.distance2( best_surface->getParameters( pt ) );
//	}
//
//	if( surface_params != nullptr ){
//		DataVector< SurfaceConstPtr > hsurfaces;
//		surface_params->getKeys( hsurfaces );
//		for(int i = 0; i < hsurfaces.countInt(); i++ ) {
//			if( hsurfaces[i] != best_surface ) { 
//				counter++;
//				const auto& pd = surface_params->getValue( hsurfaces[i] );
//				if( pd.isValid( pt ) ) {
//					double dist2 = pd.param.distance2( hsurfaces[i]->getParameters( pt ) );
//					if( dist2 > max_dist2 ) max_dist2 = dist2;
//				}
//			}
//		}
//	}
//
//	stats.add(counter);
//	return max_dist2;
//}

/// Checks the specific local surface
bool MeshPoint3d::hasLocalSurface(SurfaceConstPtr surface ) const
{
	return (local_surface != nullptr ) && local_surface->hasLocalSurface( surface );
}

/// Checks the specific local surface
bool MeshPoint3d::SurfaceData::hasLocalSurface(SurfaceConstPtr surface ) const
{
	return findSurfaceParamData( surface ) != nullptr;
}

/// Returns the local curve
Curve3dConstPtr MeshPoint3d::getLocalCurve() const {
	return (local_curve == nullptr ) ? nullptr : local_curve->getLocalCurve();
}

SurfaceConstPtr MeshPoint3d::getLocalSurface() const 
{ 
	return (local_surface == nullptr ) ? nullptr : local_surface->getLocalSurface();
}


SurfaceConstPtr MeshPoint3d::getLocalValidSurface() const 
{ 
	return (local_surface == nullptr ) ? nullptr : local_surface->getValidSurface( coord );
}

/// Returns the local surface, valid and common to all points
SurfaceConstPtr MeshPoint3d::getOptLocalCommonSurface( int ls_tag, 
			const DataVector< MeshPoint3d* > & mpoints,
			DataHashTable< SurfaceConstPtr > & hchecked_surfaces,
			SurfaceConstPtr & opt_surface, double& opt_quality) const
{
	return (local_surface == nullptr ) ? opt_surface : local_surface->getOptLocalCommonSurface( ls_tag,
		this, mpoints, hchecked_surfaces, opt_surface, opt_quality );
}

Curve3dConstPtr MeshPoint3d::CurveData::getLocalCurve(  ) 
{
	return curve_params.notEmpty() ? curve_params.getFirst().curve : nullptr;
}

SurfaceConstPtr MeshPoint3d::SurfaceData::getLocalSurface(  ) 
{
	return surface_params.notEmpty() ? surface_params.getFirst().surface : nullptr;
}

SurfaceConstPtr MeshPoint3d::SurfaceData::getValidSurface( const DPoint3d& pt ) 
{
	ParamData* pd = surface_params.findFirst( [&pt] ( ParamData& pd ) {
		return pd.updateIfNeeded(pt); } );

	return pd ? pd->surface : nullptr;
} 

Curve3dConstPtr MeshPoint3d::CurveData::getValidCurve( const DPoint3d& pt ) 
{
	ParamData* pd = curve_params.findFirst( [&pt] ( ParamData& pd ) {
		return pd.updateIfNeeded(pt); } );

	return pd ? pd->curve : nullptr;
} 

/// Returns the local surface, valid and common to all points
SurfaceConstPtr MeshPoint3d::SurfaceData::getOptLocalCommonSurface( int ls_tag, 
	const MeshPoint3d* mpoint, const DataVector< MeshPoint3d* > & mpoints,
	DataHashTable< SurfaceConstPtr > & hchecked_surfaces,
	SurfaceConstPtr & opt_surface, double & opt_quality)
{
	if( opt_quality == AQ_VALID_MAX ) return opt_surface;

	const DPoint3d& pt = mpoint->getCoordinates();
	size_t mct = mpoints.countInt();

	surface_params.forEach( [&] ( ParamData& pd ) {
		if( hchecked_surfaces.insert( pd.surface) &&
			pd.surface->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, ls_tag) ) 
		{
			double worst_quality = pd.updateIfNeededGetQuality( pt );
			for(size_t i = 0; (worst_quality >=AQ_VALID_MIN) && (i < mct); i++) {
				if( mpoints[i] == mpoint ) continue;
				double q = mpoints[i]->getLocalSurfaceQuality( pd.surface );
				if( q < worst_quality ) worst_quality = q;
			}
			LOG_GETLOCALSURFACE( "MPOINT-SURFACE-QUALITY " << worst_quality );
			if( worst_quality > opt_quality ) {
				opt_quality = worst_quality;
				opt_surface = pd.surface;
			}
		}
	} );

	return opt_surface;
}

/// Checks the specific local surface
double MeshPoint3d::getLocalSurfaceQuality(SurfaceConstPtr surface ) const
{
	assert( local_surface != nullptr );
	if(local_surface != nullptr ) 
		return local_surface->getLocalSurfaceQuality( surface, coord );
	else
		return AQ_INVALID;
}

/// Checks the specific local surface
double MeshPoint3d::SurfaceData::getLocalSurfaceQuality(SurfaceConstPtr surface, const DPoint3d& pt )
{
	assert( surface != nullptr );

	ParamData* pd = findSurfaceParamData( surface );
	if( pd != nullptr ) {
		return pd->updateIfNeededGetQuality( pt );
	}else{
		ParamData pdata( surface, pt );
		if( surface->isFixed() && pdata.quality >= AQ_VALID_MIN )
			surface_params.append( pdata );
		return pdata.quality;
	}
}

