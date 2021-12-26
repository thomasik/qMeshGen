// MeshPoly3d.cpp: implementation of the MeshPoly3d class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshPoly3d.h"
#include "MeshEdge3d.h"
#include "MeshTriangle3d.h"
#include "DTriangle.h"
#include "DQuad.h"
#include "DPlane.h"
#include "DSegment.h"
#include "DLeastSquaresFitting.h"
#include "SurfaceParametric.h"
#include "SurfacePlane.h"
#include "ControlSpace2d.h"
#include "ControlSpace2dIdentity.h"
#include "MeshPoint2d.h"
#include "MeshContainer2d.h"
#include "MeshGenerator2d.h"
#include "MeshArea.h"
#include "MeshEdge2d.h"
#include "MeshTriangle2d.h"
#include "MeshQuad3d.h"
#include "MeshViewSet.h"
#include "Metric3dContext.h"

MeshPoly3d::MeshPoly3d(const DataVector<MeshPoint3d*> & poly_points,
	MeshBlock* bounded_block) : MeshFace((int)poly_points.countInt(), 
									new MeshEdge3d*[poly_points.countInt()], 
									new MeshPoint3d*[poly_points.countInt()])
{
	blocks[0] = bounded_block;
	for(int i = 0; i < count; i++)
		points[i] = poly_points[i];

	attachToEdges();
}

MeshPoly3d::MeshPoly3d(int poly_count, MeshPoint3d** poly_points, 
	MeshBlock* bounded_block) : MeshFace(poly_count, new MeshEdge3d*[poly_count], new MeshPoint3d*[poly_count])
{
	blocks[0] = bounded_block;
	for(int i = 0; i < count; i++)
		points[i] = poly_points[i];

	attachToEdges();
}

/// Creates and returns a copy of this face (+ whole connectivity)
MeshFace* MeshPoly3d::clone() const
{
	MeshPoly3d* face = new MeshPoly3d(count, points, blocks[0]);
	face->blocks[1] = blocks[1];
	face->copyAllExtraData(this);
	return face;
}

MeshPoly3d::~MeshPoly3d()
{
	detachFromEdges();
	delete[] points;
	delete[] edges;
}

#define FIT_ERR 2e-2
/// Returns the vector normal to this face if possible
//bool MeshPoly3d::checkAndGetNormalVector(DVector3d& vn) const
//{
//	assert( count > 2);
//
//	DPoint3d middle = getMiddlePoint();
//	DVector3d n;
//
//	for(int i = 0; i < count; i++){
//		n += (points[i]->getCoordinates() - middle).crossProduct( points[(i+1)%count]->getCoordinates() - middle);
//	}
//	if(n.isZero()) return false;
//	vn = n.normalized();
//	return true;
///*
//	DataMatrix<DPoint3d> poly(count);
//	DBox box;
//	for(int i = 0; i < count; i++){
//		const DPoint3d& pt = points[i]->getCoordinates();
//		poly.add( pt );
//		box.addPoint(pt);
//	}
//	// fit plane
//	DPlane plane;
//	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(poly, plane, true);
//
//	if(plane_max_dist < 0.0){
//		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MeshPoly3d - error fitting plane");
//		return false;
//	}else if(plane_max_dist > FIT_ERR * box.getDiameter()){
//		LOG4CPLUS_WARN(MeshLog::logger_console, "MeshPoly3d - problem fitting plane, rel-dist", plane_max_dist/box.getDiameter());
//	}
//	vn = plane.vn;
//	return true;
//*/
//}
//
///// Returns the vector normal to this face
//DVector3d MeshPoly3d::getNormalVector() const
//{
//	assert( count > 2);
//
//	DPoint3d middle = getMiddlePoint();
//	DVector3d n;
//
//	for(int i = 0; i < count; i++){
//		n += (points[i]->getCoordinates() - middle).crossProduct( points[(i+1)%count]->getCoordinates() - middle);
//	}
//	if(n.isZero()) return DVector3d::zero;
//	return n.normalized();
///*
//	DataMatrix<DPoint3d> poly(count);
//	DBox box;
//	for(int i = 0; i < count; i++){
//		const DPoint3d& pt = points[i]->getCoordinates();
//		poly.add( pt );
//		box.addPoint(pt);
//	}
//	// fit plane
//	DPlane plane;
//	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(poly, plane, true);
//
//	if(plane_max_dist < 0.0){
//		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MeshPoly3d - error fitting plane");
//		assert(false);
//		return DVector3d::zero;
//	}else if(plane_max_dist > FIT_ERR * box.getDiameter()){
//		LOG4CPLUS_WARN(MeshLog::logger_console, "MeshPoly3d - problem fitting plane, rel-dist", plane_max_dist/box.getDiameter());
//	}
//	return plane.vn;
//*/
//}

/// remove point, return resulting face (the same, or a new one)
MeshFace* MeshPoly3d::removePoint(const MeshPoint3d* point)
{
	if(count < 4) return nullptr;
	else if(count == 4) {
		DataVector<MeshPoint3d*> new_pts(3);
		for(int i = 0; i < count; i++)
			if(points[i] != point) new_pts.add(points[i]);
		MeshFace* new_face = new MeshTriangle3d( new_pts[0], new_pts[1], new_pts[2], blocks[0] );
		new_face->copyAllExtraData(this);
		return new_face;
	}else if(count == 5) {
		DataVector<MeshPoint3d*> new_pts(4);
		for(int i = 0; i < count; i++)
			if(points[i] != point) new_pts.add(points[i]);
		MeshFace* new_face = new MeshQuad3d( new_pts[0], new_pts[1], new_pts[2], new_pts[3], blocks[0] );
		new_face->copyAllExtraData(this);
		return new_face;
	}else {
		// remove two edges, add one edge ...
		if(edges[0]){ // unless detached already ...
			MeshPoint3d* ept0 = nullptr;
			MeshPoint3d* ept1 = nullptr;
			for(int i = 0, j = 0; i < count; i++){
				if(edges[i]->incidentTo(point)) {
					// mark points
					if(ept0){ // second edge
						ept1 = edges[i]->getOtherPoint(point);
					}else{ // first edge
						ept0 = edges[i]->getOtherPoint(point);
					}
					// detach edges
					if(edges[i]->removeFaceLink(this))
						delete edges[i];
					// attach new edge
					if(ept1){
						edges[j] = ept0->getEdgeToPoint(ept1);
						if(!edges[j]){
							// New edge
							edges[j] = new MeshEdge3d(ept0, ept1, this);
						}else{
							edges[j]->addFaceLink(this);
						}
						j++;
					}
				}else{
					edges[j++] = edges[i];
				}
			}
		}
		// remove point
		int i = 0;
		while( (i < count) && (points[i] != point)) ++i;
		assert(i < count); // if points is part of this face...
		while( ++i < count) points[i-1] = points[i];
		--count;

		if(edges[0] == nullptr) // ie. if detached
			attachToEdges();

		// return the same face, reduced
		return this;
	}
}

std::shared_ptr<MeshViewFaceData> MeshPoly3d::getViewData(double shrink_ratio, bool proper_orientation) const
{
	if(count == 3 || count == 4)
		return MeshFace::getViewData(shrink_ratio, proper_orientation);

	DataVector<DPoint3d> poly(count);
	DBox box;
	DPlane plane;
	double f = 1.0 / count;
	for(int i = 0; i < count; i++){
		const DPoint3d& pt = points[i]->getCoordinates();
		poly.add( pt );
		plane.p0.add(pt, f);
		box.addPoint(pt);
	}
	plane.vn = getNormalVector();
	if(plane.vn.isZero()){
		//LOG4CPLUS_WARN(MeshLog::logger_console, "MeshPoly3d::view - error fitting plane");
//		assert(false);
		return nullptr;
	}
	plane.vn.orthonormalVectors( plane.e0, plane.e1 );

	// opposite orientation - for GL visualisation, normals should be outside ...

//	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(poly, plane, true);
//	if(plane_max_dist < 0.0){
//		LOG4CPLUS_WARN(MeshLog::logger_console, "MeshPoly3d::view - error fitting plane");
////		assert(false);
//		return false;
//	}else if(plane_max_dist > FIT_ERR * box.getDiameter()){
//		LOG4CPLUS_WARN(MeshLog::logger_console, "MeshPoly3d - problem fitting plane, rel-dist", plane_max_dist/box.getDiameter());
//	}

	DataVector<DPoint2d> points2d(count);
	if(proper_orientation){
		for(int i = 0; i < count; i++)
			points2d.add( plane.projectToPlane(plane.p0 + (points[i]->getCoordinates() - plane.p0) * shrink_ratio) );
	}else{
		for(int i = 0; i < count/2; i++){
			DPoint3d tpoint = poly[i];
			poly[i] = poly[count-i+1];
			poly[count-i+1] = tpoint;
		}
		for(int i = 0; i < count; i++)
			points2d.add( plane.projectToPlane(plane.p0 + (points[count-i-1]->getCoordinates() - plane.p0) * shrink_ratio) );
	}

	DataVector<int> trans_tab;
	int pct = DPoint2d::clearZeroEdges(points2d, trans_tab);

	if(pct < 3) return nullptr; // no triangle-face, really

	auto data = std::make_shared<MeshViewFaceData>(pct);
	for(int i = 0; i < pct; i++) data->pts.add(poly[ trans_tab[i] ]);
	data->quality = 1.0;
	data->normal = plane.vn;
	data->area_id = 0;

	// ... split into triangles
	DataVector<int> l(pct), r(pct);
	for(int i = 0; i < pct; i++) {
		l.add( (i+pct-1) % pct );
		r.add( (i+1) % pct);
	}
	int k = pct-1;
	int loop_counter = 0;
	while(pct > 3){
		k = r[k];
		if(loop_counter++ > pct){ // shouldn't happen ...
			if(true){
				MeshViewSet* set = new MeshViewSet();
				int ki = k;
				for(int i = 0; i < pct; i++){
					set->addPoint( plane.projectToSpace(points2d[ki]), 0, i);
					set->addEdge( plane.projectToSpace(points2d[ki]), 
									plane.projectToSpace(points2d[r[ki]]));
					ki = r[ki];
				}
				SHOW_MESH("error ear-cutting for poly", set);
			}
			int k0 = k;
			for(int i = 2; i < pct; i++){
				k = r[k];
				if( DTriangle2d::det(points2d[k0], points2d[k], points2d[r[k]]) > VERY_SMALL_NUMBER) {
					data->indices.add(k0);
					data->indices.add(k);
					data->indices.add(r[k]);
				}
			}
			return data;
		}
		// -- is triangle ok?
		DPoint2d tpts[3] = { points2d[l[k]], points2d[k], points2d[r[k]] };
		double tdet = DTriangle2d::det(tpts[0], tpts[1], tpts[2]);
		if(tdet < 0.0) continue;
		// -- does it contain any other points
		bool valid = true;
		int jk =  r[r[k]]; // first to the right, not part of the tested triangle
		for(int j = 3; valid && (j < pct); j++, jk = r[jk]) { // skipPIng the three selected vertices
			valid = !DTriangle2d::containsPoint(tpts[0], tpts[1], tpts[2], points2d[jk]);
		}
		if(valid) {
			--pct;

			if(tdet > VERY_SMALL_NUMBER) {
				data->indices.add(l[k]);
				data->indices.add(k);
				data->indices.add(r[k]);
			}
			l[ r[k] ] = l[k];
			r[ l[k] ] = r[k];

			loop_counter = 0;
		}
	}
	k = r[k];
	// ... and the last three vertices
	if( DTriangle2d::det(points2d[l[k]], points2d[k], points2d[r[k]]) > VERY_SMALL_NUMBER) {
		data->indices.add(l[k]);
		data->indices.add(k);
		data->indices.add(r[k]);
	}
	return data;
}

/// Split face into triangles
bool MeshPoly3d::splitToTriangles( DataVector<MeshFace*> & split_faces ) const
{
	if(count == 3){
		MeshFace* f = new MeshTriangle3d( points[0], points[1], points[2], blocks[0] );
		f->copyAllExtraData(this);
		split_faces.add( f );
		return true;
	}else if( count < 3) return true;
	
	int sfct_start = split_faces.countInt();

	// ear-cutting
	DataVector<DPoint3d> poly(count);
	DBox box;
	for(int i = 0; i < count; i++){
		const DPoint3d& pt = points[i]->getCoordinates();
		poly.add( pt );
		box.addPoint(pt);
	}
	DPlane plane;
	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(poly, plane, true);
	if(plane_max_dist < 0.0){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MeshPoly3d - error fitting plane");
		assert(false);
		return false;
	}else if(plane_max_dist > FIT_ERR * box.getDiameter()){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"MeshPoly3d - problem fitting plane, rel-dist: "
			<< plane_max_dist/box.getDiameter());
	}

	DataVector<DPoint2d> points2d(count);
	SurfaceConstPtr surface = getCurrentLocalSurface();
	if( surface != nullptr ) {
		DPoint2d param;
		for(int i = 0; (surface != nullptr) && (i < count); i++ ){
			if( points[i]->checkAndGetSurfaceParam( surface, param ) )
				points2d.add( param );
			else {
				points2d.clear();
				surface = nullptr;
			}
		}
	}

	if( surface == nullptr ) {
		for(int i = 0; i < count; i++)
			points2d.add( plane.projectToPlane( points[i]->getCoordinates() ) );
	}

	// ... split into triangles
	int pct = count;
	DataVector<int> l(pct), r(pct);
	for(int i = 0; i < pct; i++) {
		l.add( (i+pct-1) % pct );
		r.add( (i+1) % pct);
	}
	int k = pct-1;
	int loop_counter = 0;
	DataVector<MeshEdge3d*> inner_edges( pct-3 );
	while(pct > 3){
		k = r[k];
		if(loop_counter++ > pct){ // shouldn't happen ...
			LOG4CPLUS_WARN(MeshLog::logger_console, "Error ear-cutting for poly");
			split_faces.leaveOnly(sfct_start);
			return MeshFace::splitToTriangles(split_faces);
		}
		// -- is triangle ok?
		DPoint2d tpts[3] = { points2d[l[k]], points2d[k], points2d[r[k]] };
		double tdet = DTriangle2d::det(tpts[0], tpts[1], tpts[2]);
		if(tdet < 0.0) continue;
		// -- does it contain any other points
		bool valid = true;
		int jk =  r[r[k]]; // first to the right, not part of the tested triangle
		for(int j = 3; valid && (j < pct); j++, jk = r[jk]) { // skipPIng the three selected vertices
			valid = !DTriangle2d::containsPoint(tpts[0], tpts[1], tpts[2], points2d[jk]);
		}
		if(valid) {
			--pct;
			split_faces.add( new MeshTriangle3d( points[l[k]], points[k], points[r[k]], blocks[0] ) );
			inner_edges.add( points[l[k]]->getEdgeToPoint( points[r[k]] ) );

			l[ r[k] ] = l[k];
			r[ l[k] ] = r[k];

			loop_counter = 0;
		}
	}
	k = r[k];
	// ... and the last three vertices
	split_faces.add( new MeshTriangle3d( points[l[k]], points[k], points[r[k]], blocks[0] ) );

	for(int i = sfct_start; i < split_faces.countInt(); i++)
		split_faces[i]->copyAllExtraData(this);

	if( surface != nullptr ) {
		MeshViewSet* set1 = nullptr;
		MeshViewSet* set2 = nullptr;
		if(false){
			set1 = new MeshViewSet();
			for(int i = sfct_start; i < split_faces.countInt(); i++)
				set1->addFaceWithEdgesAndPoints(split_faces[i]);
		}
	
		bool any_swap;
		int swap_counter = 0;
		ControlSpace3dIdentity space( box.getDiameter() );
		Metric3dContext mc(&space);
		do{
			any_swap = false;
			for(int i = 0; i < inner_edges.countInt(); i++){
				MeshEdge3d* edge = inner_edges[i];
				assert( edge->getFaceCount() == 2 );
				MeshTriangle3d* t0 = (MeshTriangle3d*)edge->getFaceAt(0);
				MeshTriangle3d* t1 = (MeshTriangle3d*)edge->getFaceAt(1);
				assert( t0->getType() == FACE_TRIANGLE );
				assert( t1->getType() == FACE_TRIANGLE );
				MeshEdge3d* swap_edge = t0->swapWithNeighbour(mc, surface, t0->getNeighbourIndex(t1), true);
				if(swap_edge){
					inner_edges[i] = swap_edge;
					any_swap = true;
					++swap_counter;
				}
			}
		}while(any_swap);

		if(set1){
			if(swap_counter > 0){
				set2 = new MeshViewSet();
				for(int i = sfct_start; i < split_faces.countInt(); i++)
					set2->addFaceWithEdgesAndPoints(split_faces[i]);
				LOG4CPLUS_INFO(MeshLog::logger_console, "Split MeshPoly3d, swaps " << swap_counter);
				SHOW_MESH("poly before", set1);
				SHOW_MESH_NORESET("poly after", set2);
			}else{
				delete set1;
			}
		}
	}

	return true;
}

bool MeshPoly3d::splitToTriangles( DataVector< DTriangle3d > & triangles ) const
{
	if(count == 3){
		triangles.add( DTriangle3d( points[0]->getCoordinates(), 
			points[1]->getCoordinates(), points[2]->getCoordinates() ) );
		return true;
	}else if( count < 3) return true;
	
	int sfct_start = triangles.countInt();

	// ear-cutting
	DataVector<DPoint3d> poly(count);
	DBox box;
	for(int i = 0; i < count; i++){
		const DPoint3d& pt = points[i]->getCoordinates();
		poly.add( pt );
		box.addPoint(pt);
	}
	DPlane plane;
	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(poly, plane, true);
	if(plane_max_dist < 0.0){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MeshPoly3d - error fitting plane");
		assert(false);
		return false;
	}else if(plane_max_dist > FIT_ERR * box.getDiameter()){
		LOG4CPLUS_WARN(MeshLog::logger_console, "MeshPoly3d - problem fitting plane, rel-dist: " << plane_max_dist/box.getDiameter());
	}

	DataVector<DPoint2d> points2d(count);
	for(int i = 0; i < count; i++)
		points2d.add( plane.projectToPlane( points[i]->getCoordinates() ) );

	// ... split into triangles
	int pct = count;
	DataVector<int> l(pct), r(pct);
	for(int i = 0; i < pct; i++) {
		l.add( (i+pct-1) % pct );
		r.add( (i+1) % pct);
	}
	int k = pct-1;
	int loop_counter = 0;
	//DataVector<MeshEdge3d*> inner_edges( pct-3 );
	while(pct > 3){
		k = r[k];
		if(loop_counter++ > pct){ // shouldn't happen ...
			LOG4CPLUS_WARN(MeshLog::logger_console, "Error ear-cutting for poly");
			triangles.leaveOnly(sfct_start);
			return MeshFace::splitToTriangles(triangles);
		}
		// -- is triangle ok?
		DPoint2d tpts[3] = { points2d[l[k]], points2d[k], points2d[r[k]] };
		double tdet = DTriangle2d::det(tpts[0], tpts[1], tpts[2]);
		if(tdet < 0.0) continue;
		// -- does it contain any other points
		bool valid = true;
		int jk =  r[r[k]]; // first to the right, not part of the tested triangle
		for(int j = 3; valid && (j < pct); j++, jk = r[jk]) { // skipPIng the three selected vertices
			valid = !DTriangle2d::containsPoint(tpts[0], tpts[1], tpts[2], points2d[jk]);
		}
		if(valid) {
			--pct;
			triangles.add( DTriangle3d( points[l[k]]->getCoordinates(), 
				points[k]->getCoordinates(), points[r[k]]->getCoordinates() ) );
			//inner_edges.add( points[l[k]]->getEdgeToPoint( points[r[k]] ) );

			l[ r[k] ] = l[k];
			r[ l[k] ] = r[k];

			loop_counter = 0;
		}
	}
	k = r[k];
	// ... and the last three vertices
	triangles.add( DTriangle3d( points[l[k]]->getCoordinates(), 
		points[k]->getCoordinates(), points[r[k]]->getCoordinates() ) );

	//bool any_swap;
	//int swap_counter = 0;
	//ControlSpace3dIdentity space( box.getDiameter() );
	//Metric3dContext mc(&space);
	//do{
	//	any_swap = false;
	//	for(int i = 0; i < inner_edges.countInt(); i++){
	//		MeshEdge3d* edge = inner_edges[i];
	//		assert( edge->getFaceCount() == 2 );
	//		MeshTriangle3d* t0 = (MeshTriangle3d*)edge->getFaceAt(0);
	//		MeshTriangle3d* t1 = (MeshTriangle3d*)edge->getFaceAt(1);
	//		assert( t0->getType() == FACE_TRIANGLE );
	//		assert( t1->getType() == FACE_TRIANGLE );
	//		MeshEdge3d* swap_edge = t0->swapWithNeighbour(mc, t0->getNeighbourIndex(t1), true);
	//		if(swap_edge){
	//			inner_edges[i] = swap_edge;
	//			any_swap = true;
	//			++swap_counter;
	//		}
	//	}
	//}while(any_swap);

	return true;
}

//		-1.0 -> if self-intersections
//		0.0  -> if equilinear edges, or same vertices...
//		//0.5  -> if concave
//		1.0	 -< if convex
double MeshPoly3d::shapeQuality( const DataVector<DPoint3d> & poly )
{
	DBox box;
	int count = poly.countInt();
	for(int i = 0; i < count; i++)
		box.addPoint(poly[i]);

	DPlane plane;
	double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(poly, plane);
	if(plane_max_dist < 0.0){
//		LOG4CPLUS_ERROR(MeshLog::logger_console,   "MeshPoly3d - error fitting plane");
//		assert(false);
		return 0.0;
	}else if(plane_max_dist > FIT_ERR * box.getDiameter()){
		LOG4CPLUS_WARN(MeshLog::logger_console, "MeshPoly3d - problem fitting plane, rel-dist: " 
			<< plane_max_dist/box.getDiameter());
		// return 0.0; ???
	}

	DataVector<DPoint2d> points2d(count);
	for(int i = 0; i < count; i++)
		points2d.add( plane.projectToPlane( poly[i] ) );

	DataVector<std::shared_ptr<DSegment2d>> segments(count);
	for(int i = 0; i < count; i++){
		auto psegment = std::make_shared<DSegment2d>(points2d[i], points2d[(i+1)%count]);
		if(psegment->length2() > 0.0) 
			segments.add( psegment );
	}
	int sct = segments.countInt();
	for(int i = 0; i < sct-2; i++){
		int sct_i = (i==0) ? sct-1 : sct;
		for(int j = i+2; j < sct_i; j++){
			if( segments[i]->crossingSegment( *segments[j] ) ) return -1.0;
		}
	}

	return (sct < count) ? 0.0 : 1.0;
}

/// returns shape quality for a face
double MeshPoly3d::getShapeQuality() const
{
	// ... return:
	//		-1.0 -> if self-intersections
	//		0.0  -> if equilinear edges, or same vertices...
	//		//0.5  -> if concave
	//		1.0	 -< if convex

	if(count < 3) return 0.0;
	else if(count == 3)
		return DTriangle3d::alphaQuality(points[0]->getCoordinates(),
				points[1]->getCoordinates(), points[2]->getCoordinates());
	else if(count == 4){
		return DQuad3d::shapeQuality(
				points[0]->getCoordinates(), points[1]->getCoordinates(),
				points[2]->getCoordinates(), points[3]->getCoordinates());
	}

	DataVector<DPoint3d> poly(count);
	for(int i = 0; i < count; i++){
		poly.add( points[i]->getCoordinates() );
	}

	if(false){
		MeshViewSet* set = new MeshViewSet;
		set->addFaceWithEdgesAndPoints(this);
		SHOW_MESH("MeshPoly3d::getShapeQuality", set);
	}

	return MeshPoly3d::shapeQuality( poly );
}

/// returns shape quality for a face
double MeshPoly3d::getShapeQuality(Metric3dContext& mc) const
{
	if(count < 3) return 0.0;
	else if(count == 3)
		return DTriangle3d::alphaQuality(
				points[0]->getMetricCoordinates(mc).toRealSpace(),
				points[1]->getMetricCoordinates(mc).toRealSpace(), 
				points[2]->getMetricCoordinates(mc).toRealSpace());
	else if(count == 4){
		return DQuad3d::shapeQuality(
				points[0]->getMetricCoordinates(mc).toRealSpace(), 
				points[1]->getMetricCoordinates(mc).toRealSpace(),
				points[2]->getMetricCoordinates(mc).toRealSpace(), 
				points[3]->getMetricCoordinates(mc).toRealSpace());
	}

	DataVector<DPoint3d> poly(count);
	for(int i = 0; i < count; i++){
		poly.add( points[i]->getMetricCoordinates(mc).toRealSpace() );
	}

	if(false){
		MeshViewSet* set = new MeshViewSet;
		set->addFaceWithEdgesAndPoints(this);
		SHOW_MESH("MeshPoly3d::getShapeQuality", set);
	}

	return MeshPoly3d::shapeQuality( poly );
}

// check, whether this face is valid with respect to the given surface
bool MeshPoly3d::valid( SurfaceConstPtr surface ) const
{
	// ... return:
	//		false	-> if self-intersections
	//		false	-> if equilinear edges, or same vertices...
	//		//true	-> if concave
	//		true	-> if convex

	if(count < 3) return false;
	else if(count == 3)
		return DTriangle2d::valid(
			getPoint(0)->getLocalSurfaceParam( surface ),
			getPoint(1)->getLocalSurfaceParam( surface ), 
			getPoint(2)->getLocalSurfaceParam( surface ) );
	else if(count == 4){
		return DQuad2d::valid(
			getPoint(0)->getLocalSurfaceParam( surface ),
			getPoint(1)->getLocalSurfaceParam( surface ), 
			getPoint(2)->getLocalSurfaceParam( surface ), 
			getPoint(3)->getLocalSurfaceParam( surface )
		);
	}else{
		DataVector<std::shared_ptr<DSegment2d>> segments(count);
		for(int i = 0; i < count; i++){
			int j = (i+1)%count;
			auto psegment = std::make_shared<DSegment2d>(
				points[i]->getLocalSurfaceParam( surface ), 
				points[j]->getLocalSurfaceParam( surface ) );
			if(psegment->length2() > 0.0) 
				segments.add( psegment );
			else return false;
		}
		assert( segments.countInt() == count );
		for(int i = 0; i < count-2; i++){
			int sct_i = (i==0) ? count-1 : count;
			for(int j = i+2; j < sct_i; j++){
				if( segments[i]->crossingSegment( *segments[j] ) ) return false;
			}
		}

		return true;
	}
}

/*

// Copyright 2001 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
 

// Assume that classes are already given for the objects:
//    Point with 2D coordinates {float x, y;}
//    Polygon with n vertices {int n; Point *V;} with V[n]=V[0]
//    Tnode is a node element structure for a BBT
//    BBT is a class for a Balanced Binary Tree
//        such as an AVL, a 2-3, or a  red-black tree
//        with methods given by the  placeholder code:

typedef struct _BBTnode Tnode;
struct _BBTnode {
        void* val;
        // plus node mgmt info ...
};

class BBT {
        Tnode   *root;
public:
                 BBT() {root = (Tnode*)0;}   // constructor
                 ~BBT() {freetree();}        // destructor

        Tnode*  insert( void* ){};       // insert data into the tree
        Tnode*  find( void* ){};         // find data from the tree
        Tnode*  next( Tnode* ){};        // get next tree node
        Tnode*  prev( Tnode* ){};        // get previous tree node
        void    remove( Tnode*  ){};     // remove node from the tree
        void    freetree(){};            // free all tree data structs
};
// NOTE:
// Code for these methods must be provided for the algorithm to work.
// We have not provided it since binary tree algorithms are well-known
// and code is widely available. Further, we want to reduce the clutter
// accompanying the essential sweep line algorithm.
//===================================================================
 

#define FALSE   0
#define TRUE    1
#define LEFT    0
#define RIGHT   1

extern void
qsort(void*, unsigned, unsigned, int(*)(const void*,const void*));

// xyorder(): determines the xy lexicographical order of two points
//      returns: (+1) if p1 > p2; (-1) if p1 < p2; and  0 if equal
int xyorder( Point* p1, Point* p2 )
{
    // test the x-coord first
    if (p1->x > p2->x) return 1;
    if (p1->x < p2->x) return (-1);
    // and test the y-coord second
    if (p1->y > p2->y) return 1;
    if (p1->y < p2->y) return (-1);
    // when you exclude all other possibilities, what remains  is...
    return 0;  // they are the same point
}

// isLeft(): tests if point P2 is Left|On|Right of the line P0 to P1.
//      returns: >0 for left, 0 for on, and <0 for  right of the line.
//      (see Algorithm 1 on Area of Triangles)
inline float
isLeft( Point P0, Point P1, Point P2 )
{
    return (P1.x - P0.x)*(P2.y - P0.y) - (P2.x - P0.x)*(P1.y -  P0.y);
}
//===================================================================
 


// EventQueue Class

// Event element data struct
typedef struct _event Event;
struct _event {
    int      edge;          // polygon edge i is V[i] to V[i+1]
    int      type;          // event type: LEFT or RIGHT vertex
    Point*   eV;            // event vertex
};

int E_compare( const void* v1, const void* v2 ) // qsort compare two events
{
    Event**    pe1 = (Event**)v1;
    Event**    pe2 = (Event**)v2;

    return xyorder( (*pe1)->eV, (*pe2)->eV );
}

// the EventQueue is a presorted array (no insertions needed)
class EventQueue {
    int      ne;                // total number of events in array
    int      ix;                // index of next event on queue
    Event*   Edata;             // array of all events
    Event**  Eq;                // sorted list of event pointers
public:
              EventQueue(Polygon P);     // constructor
             ~EventQueue(void)           // destructor
                  { delete Eq; delete Edata;}

    Event*   next();                     // next event on queue
};

// EventQueue Routines
EventQueue::EventQueue( Polygon P )
{
    ix = 0;
    ne = 2 * P.n;           // 2 vertex events for each edge
    Edata = (Event*)new Event[ne];
    Eq = (Event**)new (Event*)[ne];
    for (int i=0; i < ne; i++)           // init Eq array pointers
        Eq[i] = &Edata[i];

    // Initialize event queue with edge segment endpoints
    for (int i=0; i < P.n; i++) {        // init data for edge i
        Eq[2*i]->edge = i;
        Eq[2*i+1]->edge = i;
        Eq[2*i]->eV   = &(P.V[i]);
        Eq[2*i+1]->eV = &(P.V[i+1]);
        if (xyorder( &P.V[i], &P.V[i+1]) < 0)  { // determine type
            Eq[2*i]->type    = LEFT;
             Eq[2*i+1]->type = RIGHT;
        }
        else {
            Eq[2*i]->type    = RIGHT;
             Eq[2*i+1]->type = LEFT;
        }
    }
    // Sort Eq[] by increasing x and y
    qsort( Eq, ne, sizeof(Event*), E_compare );
}

Event* EventQueue::next()
{
    if (ix >= ne)
        return (Event*)0;
    else
        return Eq[ix++];
}
//===================================================================
 


// SweepLine Class

// SweepLine segment data struct
typedef struct _SL_segment SLseg;
struct _SL_segment {
    int      edge;          // polygon edge i is V[i] to V[i+1]
    Point    lP;            // leftmost vertex point
    Point    rP;            // rightmost vertex point
    SLseg*   above;         // segment above this one
    SLseg*   below;         // segment below this one
};

// the Sweep Line itself
class SweepLine {
    int      nv;            // number of vertices in polygon
    Polygon* Pn;            // initial Polygon
    BBT      Tree;          // balanced binary tree
public:
              SweepLine(Polygon P)            // constructor
                  { nv = P.n; Pn = &P; }
             ~SweepLine(void)                 // destructor
                  { Tree.freetree();}

    SLseg*   add( Event* );
    SLseg*   find( Event* );
    int      intersect( SLseg*, SLseg*  );
    void     remove( SLseg* );
};

SLseg* SweepLine::add( Event* E )
{
    // fill in SLseg element data
    SLseg* s = new SLseg;
    s->edge  = E->edge;

    // if it is being added, then it must be a LEFT edge event
    // but need to determine which endpoint is the left one
    Point* v1 = &(Pn->V[s->edge]);
    Point* v2 = &(Pn->V[s->edge+1]);
    if (xyorder( v1, v2) < 0) { // determine which is leftmost
        s->lP = *v1;
        s->rP = *v2;
    }
    else {
        s->rP = *v1;
        s->lP = *v2;
    }
    s->above = (SLseg*)0;
    s->below = (SLseg*)0;

    // add a node to the balanced binary tree
    Tnode* nd = Tree.insert(s);
    Tnode* nx = Tree.next(nd);
    Tnode* np = Tree.prev(nd);
    if (nx != (Tnode*)0) {
        s->above = (SLseg*)nx->val;
        s->above->below = s;
    }
    if (np != (Tnode*)0) {
        s->below = (SLseg*)np->val;
        s->below->above = s;
    }
    return s;
}

SLseg* SweepLine::find( Event* E )
{
    // need a segment to find it in the tree
    SLseg* s = new SLseg;
    s->edge  = E->edge;
    s->above = (SLseg*)0;
    s->below = (SLseg*)0;

    Tnode* nd = Tree.find(s);
    delete s;
    if (nd == (Tnode*)0)
        return (SLseg*)0;

    return (SLseg*)nd->val;
}

void SweepLine::remove( SLseg* s )
{
    // remove the node from the balanced binary tree
    Tnode* nd = Tree.find(s);
    if (nd == (Tnode*)0)
        return;       // not there

    // get the above and below segments pointing to each other
    Tnode* nx = Tree.next(nd);
    if (nx != (Tnode*)0) {
        SLseg* sx = (SLseg*)(nx->val);
        sx->below = s->below;
    }
    Tnode* np = Tree.prev(nd);
    if (np != (Tnode*)0) {
        SLseg* sp = (SLseg*)(np->val);
        sp->above = s->above;
    }
    Tree.remove(nd);       // now  can safely remove it
    delete s;
}

// test intersect of 2 segments and return: 0=none, 1=intersect
int SweepLine::intersect( SLseg* s1, SLseg* s2)
{
    if (s1 == (SLseg*)0 || s2 == (SLseg*)0)
        return FALSE;       // no intersect if either segment doesn't exist

    // check for consecutive edges in polygon
    int e1 = s1->edge;
    int e2 = s2->edge;
    if (((e1+1)%nv == e2) || (e1 == (e2+1)%nv))
        return FALSE;       // no non-simple intersect since consecutive

    // test for existence of an intersect point
    float lsign, rsign;
    lsign = isLeft(s1->lP, s1->rP, s2->lP);    //  s2 left point sign
    rsign = isLeft(s1->lP, s1->rP, s2->rP);    //  s2 right point sign
    if (lsign * rsign > 0) // s2 endpoints have same sign  relative to s1
        return FALSE;       // => on same side => no intersect is possible
    lsign = isLeft(s2->lP, s2->rP, s1->lP);    //  s1 left point sign
    rsign = isLeft(s2->lP, s2->rP, s1->rP);    //  s1 right point sign
    if (lsign * rsign > 0) // s1 endpoints have same sign  relative to s2
        return FALSE;       // => on same side => no intersect is possible
    // the segments s1 and s2 straddle each other
    return TRUE;            // => an intersect exists
}
//===================================================================
 


// simple_Polygon(): test if a Polygon is simple or not
//     Input:  Pn = a polygon with n vertices V[]
//     Return: FALSE(0) = is NOT simple
//             TRUE(1)  = IS simple
int
simple_Polygon( Polygon Pn )
{
    EventQueue  Eq(Pn);
    SweepLine   SL(Pn);
    Event*      e;                  // the current event
    SLseg*      s;                  // the current SL segment

    // This loop processes all events in the sorted queue
    // Events are only left or right vertices since
    // No new events will be added (an intersect => Done)
    while (e = Eq.next()) {         // while there are events
        if (e->type == LEFT) {      // process a left vertex
            s = SL.add(e);          // add it to the sweep line
            if (SL.intersect(  s, s->above))
                 return FALSE;      // Pn is NOT simple
            if (SL.intersect(  s, s->below))
                 return FALSE;      // Pn is NOT simple
        }
        else {                      // processs a right vertex
            s = SL.find(e);
            if (SL.intersect(  s->above, s->below))
                 return FALSE;      // Pn is NOT simple
            SL.remove(s);           // remove it from the sweep line
        }
    }
    return TRUE;      // Pn IS simple
}
//===================================================================
*/