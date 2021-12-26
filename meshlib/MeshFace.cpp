// MeshFace.cpp: implementation of the MeshFace class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshData.h"
#include "MeshFace.h"
#include "MeshPoint3d.h"
#include "MeshViewSet.h"
#include "MeshTetrahedron.h"
#include "MeshBoundaryCondition.h"
#include "MeshEdge3d.h"
#include "Metric3dContext.h"
#include "DTriangle.h"
#include "MeshGenerator3dSurface.h"
#include "SurfaceParametric.h"
#include "MeshTriangle3d.h"
#include "MeshGenerator2d.h"
#include "MeshContainer2d.h"
#include "MeshTriangle2d.h"
#include "DPlane.h"
#include "GeometricPredicates.h"

MeshFace::MeshFace(int _count, MeshEdge3d** _edges, MeshPoint3d** _points) 
	: count(_count), edges(_edges), points(_points), sdata(nullptr)
{
	blocks[0] = blocks[1] = nullptr;
}

/// copy additional data for cloning
void MeshFace::copyAllExtraData(const MeshFace* source)
{
	copyBorderFlagsFrom(source);
	copyAllTags(source);
	copySurfaceData(source);
}

MeshEdge3d* MeshFace::createEdge(MeshPoint3d* p1, MeshPoint3d* p2)
{
	return new MeshEdge3d(p1, p2, this);
}


/// Creates and returns a copy of this face (+ whole connectivity)
MeshFace* MeshFace::clone() const
{
	MeshFace* face = new MeshFace(count, edges, points);
	face->copyAllExtraData(this);
	return face;
}

bool MeshFace::clearBlockLink(MeshBlock *block)
{
	if(blocks[0] == block){
		blocks[0] = nullptr;
	}else if(blocks[1] == block){
		blocks[1] = nullptr;
	}else assert(false);

	return (!blocks[0] && !blocks[1] && !isBorder());
}

void MeshFace::setBlockLink(MeshBlock *block, MeshFace::BlockOrientation side)
{
	int i = 1-side;
	if(blocks[i] == nullptr){
		blocks[i] = block;
	}else{
		assert(!blocks[1-i]);
		blocks[1-i] = block;
	}
}

DBox MeshFace::getBoundingBox() const
{
	DBox box;
	for(int i = 0; i < count; i++){
		box.addPoint(points[i]->getCoordinates());
	}
	return box;
}

DPoint3d MeshFace::getMiddlePoint() const
{
	DPoint3d middle = points[0]->getCoordinates();
	for(int i = 1; i < count; i++){
		middle.add(points[i]->getCoordinates());
	}
	middle /= (double)count;
	return middle;
}

bool MeshFace::incidentToPoint(const MeshPoint3d *point) const
{
	for(int i = 0; i < count; i++) 
		if(points[i] == point) return true;
	return false;
}

bool MeshFace::incidentToPoints(const MeshPoint3d *point1, const MeshPoint3d *point2) const
{
	bool one_ok = false;
	for(int i = 0; i < count; i++) 
		if(points[i] == point1 || points[i] == point2) 
			if(one_ok) return true; // i.e. found second one
			else one_ok = true;
	return false;
}

int MeshFace::getBlockIndex(const MeshBlock *block) const
{
	if(blocks[0] == block) return 0;
	else if(blocks[1] == block) return 1;
	else return -1;
}

int MeshFace::getOrientation(const DPoint3d &pt) const
{
	const DVector3d normal = 
		(points[1]->getCoordinates() - points[0]->getCoordinates()).crossProduct(
		 points[2]->getCoordinates() - points[0]->getCoordinates());
	if(normal.scalarProduct(pt - points[0]->getCoordinates()) < 0.0)
		return BLOCK_DOWN;
	else return BLOCK_UP;
}

void MeshFace::setBlockLink(MeshBlock *block, const DPoint3d &top_point)
{
	setBlockLink(block, getOrientation(top_point)?MeshFace::BLOCK_DOWN:MeshFace::BLOCK_UP);
}

/// Returns the vector normal to this face if possible
bool MeshFace::checkAndGetNormalVector(DVector3d& vn) const
{
	assert( count > 2);

	const DPoint3d& ptA = points[0]->getCoordinates();
	DVector3d n;

	for(int i = 2; i < count; i++){
		n += (points[i-1]->getCoordinates() - ptA).crossProduct( points[i]->getCoordinates() - ptA);
	}

	if(n.isZero()) return false;
	vn = n.normalized();
	return true;

	//if(count < 3) return false;
	//vn = (points[1]->getCoordinates() - points[0]->getCoordinates()).crossProduct(
	//	 points[2]->getCoordinates() - points[0]->getCoordinates());
	//double len = vn.length();
	//if(len < VERY_SMALL_NUMBER) return false;
	//vn /= len;
	//return true;
}

DVector3d MeshFace::getNormalVector() const
{
	assert( count > 2);

	const DPoint3d& ptA = points[0]->getCoordinates();
	DVector3d n;

	for(int i = 2; i < count; i++){
		n += (points[i-1]->getCoordinates() - ptA).crossProduct( points[i]->getCoordinates() - ptA);
	}

	if(n.isZero()) return DVector3d::zero;
	return n.normalized();

	//assert(count > 2);
	//return (points[1]->getCoordinates() - points[0]->getCoordinates()).crossProduct(
	//	 points[2]->getCoordinates() - points[0]->getCoordinates()).normalized();
}

int MeshFace::getEdgeIndex(const MeshEdge3d *edge) const
{
	for(int i = 0; i < count; i++) 
		if(edges[i] == edge) return i;
	return -1;
}

void MeshFace::setBlockLink(MeshBlock *block, const MeshPoint3d* mpt1, const MeshPoint3d* mpt2)
{
	int i = 0;
	for(i = 0; i < count; i++){
		if(points[i] == mpt1){
			if(points[(i+1)%count] == mpt2){
				setBlockLink(block, BLOCK_UP);
				return;
			}else{
				assert(points[(i+count-1)%count] == mpt2);
				setBlockLink(block, BLOCK_DOWN);
				return;
			}
		}
	}
	assert(false);
}

bool MeshFace::properOrientation(const MeshPoint3d *mp1, const MeshPoint3d *mp2) const
{	
	for(int i = 0; i < count; i++){
		if(points[i] == mp1)
			return (points[(i+1)%count] == mp2);
	}
	assert(false);
	return false;
}

/// Returns true if both faces are oriented the same way (checking vertices)
bool MeshFace::sameOrientation(const MeshFace* face) const
{
	for(int i = 0; i < count; i++){
		if(face->incidentToPoints(points[i], points[(i+1)%count]))
			return face->properOrientation(points[(i+1)%count], points[i]);
	}
	assert(false);
	return false;
}


//////////////////////////////////////////////////////////////////////
// Funkcja zwraca kolejny punkt wzglêdem podanego
MeshPoint3d* MeshFace::getNextPoint(const MeshPoint3d *point) const
{
	for(int i = 0; i < count; i++){
		if(points[i] == point) return points[(i+1) % count];
	}
	return nullptr;
}

//////////////////////////////////////////////////////////////////////
// Funkcja zwraca poprzedni punkt wzglêdem podanego
MeshPoint3d* MeshFace::getPrevPoint(const MeshPoint3d *point) const
{
	for(int i = 0; i < count; i++){
		if(points[i] == point) 
			return points[(i + count - 1) % count];
	}
	return nullptr;
}

bool MeshFace::incidentToEdge(const MeshEdge3d *edge) const
{
	for(int i = 0; i < count; i++) if(edges[i] == edge) 
		return true;
	return false;
}

MeshPoint3d* MeshFace::getOtherPoint(const MeshPoint3d *pt0, const MeshPoint3d *pt1)
{
	for(int i = 0; i < count; i++) 
		if(points[i] != pt0 && points[i] != pt1) 
			return points[i];
	return nullptr;
}

// Checks, whether both points are "visible" from the same side of this face
bool MeshFace::sameSide(const DPoint3d& ptA, const DPoint3d& ptB) const
{
	const DPoint3d& pt0 = points[0]->getCoordinates();
	const DPoint3d& pt1 = points[1]->getCoordinates();
	const DPoint3d& pt2 = points[2]->getCoordinates();

	return (DTriangle3d::orient3d(pt0, pt1, pt2, ptA) *
		DTriangle3d::orient3d(pt0, pt1, pt2, ptB)) > 0.0;
}

bool MeshFace::outsideOrientation(MeshBlock * block) const
{
	MeshPoint3d* opp_point = block->getOppositePoint(this);
	assert(opp_point != nullptr);
	assert(!this->incidentToPoint(opp_point));

	double det = GeometricPredicates::orient3d_fast(
		points[0]->getCoordinates(),
		points[1]->getCoordinates(),
		points[2]->getCoordinates(),
		opp_point->getCoordinates());

	return det < 0.0;
}

std::shared_ptr<MeshViewFaceData> MeshFace::getViewData(double shrink_ratio, bool proper_orientation) const
{
	
	auto data = std::make_shared<MeshViewFaceData>(count);
	data->area_id = 0;
	data->quality = (float)getDoubleTag( TagExtended::TAG_QUALITY, 1.0 );
	DPoint3d middle = getMiddlePoint();

	if(proper_orientation){
		for(int i = 0; i < count; i++)
			data->pts.add( middle + (points[i]->getCoordinates() - middle) * shrink_ratio );
	}else{
		for(int i = 0; i < count; i++)
			data->pts.add( middle + (points[count-i-1]->getCoordinates() - middle) * shrink_ratio );
	}

	for(int i = 2; i < count; i++){
		data->indices.add( 0 );
		data->indices.add( i-1 );
		data->indices.add( i );
	}

	data->normal = getNormalVector();
	if (!proper_orientation)
		data->normal.invert();
	if(data->normal.isZero()){
		return nullptr;
	}else{
		return data;
	}
}

void MeshFace::setBorderWithEdges(char border_flags)
{
	setBorder(border_flags);
	for(int i = 0; i < count; i++)
		if(edges[i]) edges[i]->setBorder(border_flags);
}

void MeshFace::setBorderFlagsWithEdges(char border_flags)
{
	setBorderFlags(border_flags);
	for(int i = 0; i < count; i++)
		if(edges[i]) edges[i]->setBorderFlags(border_flags);
}

void MeshFace::setBorderWithEdgesAndPoints(char border_flags)
{
	setBorder(border_flags);
	for(int i = 0; i < count; i++){
		if(edges[i]) edges[i]->setBorder(border_flags);
		if(points[i]) points[i]->setBorder(border_flags);
	}
}

bool MeshFace::removeBlockLink(const MeshBlock *block)
{
	int i = (blocks[0] == block)?0:1;
	if(blocks[i] != block) return false;

	blocks[i] = nullptr;

	// If face is not boundary and is not incident to any blocks, it can be removed (return true)
	return (!blocks[1-i] && !isBorder());
}

bool MeshFace::removeBlockLinkDisregardBoundary(const MeshBlock *block)
{
	int i = (blocks[0] == block)?0:1;
	if(blocks[i] != block) return false;

	blocks[i] = nullptr;

	// If face is not incident to any blocks, it can be removed (return true)
	return !blocks[1-i];
}

/// Returns edge from this face adjacent to the given edge through point
MeshEdge3d* MeshFace::getOtherEdge(const MeshEdge3d* edge, const MeshPoint3d* point)
{
	for(int i = 0; i < count; i++){
		if(points[i] == point){
			if(edges[i] == edge)
				return edges[(i+count-1)%count];
			else
				return edges[i];
		}
	}
	return nullptr;
}

//////////////////////////////////////////////////////////////////////
// Switch old point into new one, with update of adjacent edges
void MeshFace::switchPointsWithEdges(const MeshPoint3d *point1, MeshPoint3d *point2)
{
	for(int i = 0; i < count; i++){
		if(points[i] == point1){
			points[i] = point2;
			break;
		}
	}
	for(int i = 0; i < count; i++){
		MeshEdge3d* edge = edges[i];
		if(edge->getMeshPoint(0) == point1 || edge->getMeshPoint(1) == point1){
			MeshPoint3d* other_point = edge->getOtherPoint(point1);
			if(edge->removeFaceLink(this)) delete edge;
			edges[i] = other_point->getEdgeToPoint(point2);
			if(!edges[i]){
				// New edge
				edges[i] = createEdge(other_point, point2);
			}else{
				edges[i]->addFaceLink(this);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Switch old point into new one, with update of adjacent edges
void MeshFace::switchPointsWithEdgesBoundary(const MeshPoint3d *point1, MeshPoint3d *point2)
{
	for(int i = 0; i < count; i++){
		if(points[i] == point1){
			points[i] = point2;
			break;
		}
	}
	for(int i = 0; i < count; i++){
		MeshEdge3d* edge = edges[i];
		if(edge->getMeshPoint(0) == point1 || edge->getMeshPoint(1) == point1){
			MeshPoint3d* other_point = edge->getOtherPoint(point1);
			char is_border = edge->getBorderFlags();
			if(edge->removeFaceLinkDisregardBoundary(this)) delete edge;
			edges[i] = other_point->getEdgeToPoint(point2);
			if(!edges[i]){
				// New edge
				edges[i] = createEdge(other_point, point2);
			}else{
				edges[i]->addFaceLink(this);
			}
			if(is_border) edges[i]->setBorderFlags(is_border);
		}
	}
}

/// Create face-edges connections
void MeshFace::attachToEdges()
{
	for(int i = 0; i < count; i++){
		edges[i] = points[i]->getEdgeToPoint(points[(i+1)%count]);
		if(!edges[i]){
			// New edge
			edges[i] = createEdge(points[i], points[(i+1)%count]);
		}else{
			edges[i]->addFaceLink(this);
		}
	}
}

/// Removes face-edges connections
void MeshFace::detachFromEdges()
{
	if(count > 0 && edges[0]){
		for(int i = 0; i < count; i++){
			if(edges[i]){
				if(edges[i]->removeFaceLink(this))
					delete edges[i];
				edges[i] = nullptr;
			}
		}
	}
}

/// Removes face-edges connections
void MeshFace::detachFromEdgesDisregardBoundary()
{
	if(count > 0 && edges[0]){
		for(int i = 0; i < count; i++){
			if(edges[i]){
				if(edges[i]->removeFaceLinkDisregardBoundary(this))
					delete edges[i];
				edges[i] = nullptr;
			}
		}
	}
}


DVector3d MeshFace::getLocalSurfaceNormal( SurfaceConstPtr surface) const
{
	int reversed = 1;
	if( surface == nullptr ) {
		if( sdata != nullptr ) {
			surface = sdata->local_surface;
			reversed = sdata->local_surface_orientation;
		}else {
			assert( false );
			return DVector3d::zero;
		}
	}

	double w = 1.0 / count;
	DPoint2d mid_param;
	for(int i = 0; i < count; i++)
		mid_param.add( points[i]->getLocalSurfaceParam( surface ), w );

	return surface->getNormalVector( mid_param ) * reversed;
}

// check, whether this triangular face is valid (second version)
bool MeshFace::validDirect(Metric3dContext& mc, double last_q) const
{
	double skip_len2 = mc.getMinSkipLen2();
	if(skip_len2 > 0.0){
		int ect = getEdgeCount();
		for(int i = 0; i < ect; i++)
			if( getEdge(i)->getLength2() < skip_len2) return true; // no point in checking, will be removed later anyway...
	}

	static const double MIN_ALPHA_QUALITY = 0.01;

	// 1. Quality
	if(last_q >= 0.0){
		double q = getShapeQuality(mc);
		if( ( (q <= 0.0) && (last_q > 0.0) ) ||
			( (q < MIN_ALPHA_QUALITY) && (last_q > MIN_ALPHA_QUALITY) ) ) {
	//	if((q < MIN_ALPHA_QUALITY) && (q < last_q)){
//			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Invalid quality (face #" << getIndex() << "): q= " << q << ", last_q= " << last_q);
			return false;
		}
	}

	// 2. Normal
	static const double MAX_ANGLE_PRODUCT = SQRT2*0.5;

	if(sdata != nullptr) {
		const DVector3d vn = getNormalVector();

		if(sdata->local_surface){
			// check using normal for face and normal for surface
			double sp = vn.scalarProduct( getLocalSurfaceNormal() );
			if(sp <= 0.0){
//				LOG4CPLUS_INFO(MeshLog::logger_mesh, "Invalid surf normal (face #" << getIndex() << "): sp= " << sp);
				return false;
			}
		}else if(! sdata->base_normal.isZero()){
			// compare with base value of normal
			double sp = vn.scalarProduct(sdata->base_normal);
			if(sp < MAX_ANGLE_PRODUCT){
//				LOG4CPLUS_INFO(MeshLog::logger_mesh, "Invalid last normal (face #" << getIndex() << "): sp= " << sp);
				return false;
			}

			// ? compare with average normals of adjacent faces ?
		}
	}

	return true;
}

void MeshFace::setTagForEdges(TagExtended::TagType type, int t) const 
{ 
	for(int i = 0; i < count; i++) 
		edges[i]->setIntTag(type, t); 
}

/// Set tag for all adjacent edges
void MeshFace::setIntTagForPoints(TagExtended::TagType type, int t) const
{
	for(int i = 0; i < count; i++)
		points[i]->setIntTag(type, t);
}

/// Returns the edge between the two faces, or nullptr if not adjacent
MeshEdge3d* MeshFace::getCommonEdge(MeshFace* face) const
{
	for(int i = 0; i < count; i++){
		for(int j = 0; j < count; j++){
			if(edges[i] == face->edges[j]) return edges[i];
		}
	}
	return nullptr;
}

/// Removes link to this edge (fo delete-all phase) - returns true if last one
bool MeshFace::removeEdgeLink(MeshEdge3d* edge)
{
	bool no_more_edges = true;
	for(int i = 0; i < count; i++){
		if(edges[i] == edge) edges[i] = 0;
		else if(edges[i]) no_more_edges = false;
	}

	return no_more_edges;
}

/// clears surface data
void MeshFace::clearSurfaceData()
{
	if(sdata){
		delete sdata;
		sdata = nullptr;
	}
}

/// Copy all surface data
void MeshFace::copySurfaceData(const MeshFace* source)
{
	if(source->sdata){
		if(!sdata) sdata = new SurfaceData();
		sdata->base_normal = source->sdata->base_normal;
		sdata->local_surface = source->sdata->local_surface;
		sdata->local_surface_orientation = source->sdata->local_surface_orientation;
		sdata->ce0 = source->sdata->ce0;
		sdata->ce1 = source->sdata->ce1;
		sdata->cv0 = source->sdata->cv0;
		sdata->cv1 = source->sdata->cv1;
	}else{
		if(sdata) clearSurfaceData();
	}
}

/// Returns local surface
SurfaceConstPtr MeshFace::getOptLocalSurface() {
	SurfaceConstPtr opt_surface = nullptr;
	if( sdata!= nullptr ) opt_surface = sdata->local_surface;

	if( opt_surface != nullptr ) {
		double worst_quality = AQ_VALID_MAX;
		for(int i = 0; (worst_quality >= 0.0) && (i < count); i++ ) {
			double q = points[i]->getLocalSurfaceQuality( opt_surface );
			if( q < worst_quality ) worst_quality = q;
		}
		if( worst_quality >= AQ_VALID_MIN ) return opt_surface; // the current one is ok.
	}

	int local_surface_tag = getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );

	const unsigned int HCT = 100;
	DataHashTable< SurfaceConstPtr > hchecked_surfaces( HCT, nullptr );

	// All points have the same valid surface ?
	DataVector< MeshPoint3d* > mpoints( count ) ;
	for(int i = 0 ; i < count; i++) mpoints.add( points[i] );
	double opt_quality = AQ_INVALID;

	for(int i = 0; (i < count); i++ ) {
		mpoints[i]->getOptLocalCommonSurface( local_surface_tag, mpoints, hchecked_surfaces, opt_surface, opt_quality );
	}
	assert( opt_surface != nullptr );
	if( opt_surface != nullptr ) setLocalSurface( opt_surface );
	
	return opt_surface;
}

/// Split face into triangles
bool MeshFace::splitToTriangles( DataVector<MeshFace*> & split_faces ) const
{
	assert(count > 2);
	if(count < 3) return true; // ie. wasn't able to fit a plane...

	size_t sfct_start = split_faces.countInt();

	if(count == 3){
		split_faces.add( new MeshTriangle3d( points[0], points[1], points[2], blocks[0] ) );
	}else if( count == 4) {
		// like quad...
		const DPoint3d& dpt0 = points[0]->getCoordinates();
		const DPoint3d& dpt1 = points[1]->getCoordinates();
		const DPoint3d& dpt2 = points[2]->getCoordinates();
		const DPoint3d& dpt3 = points[3]->getCoordinates();

		const DVector3d vn1 = DVector3d::crossProduct(dpt0, dpt1, dpt2);
		const DVector3d vn2 = DVector3d::crossProduct(dpt0, dpt2, dpt3);
		const DVector3d vn3 = DVector3d::crossProduct(dpt0, dpt1, dpt3);
		const DVector3d vn4 = DVector3d::crossProduct(dpt1, dpt2, dpt3);

		double s12 = (vn1.isZero() || vn2.isZero()) ? 0.0 : vn1.scalarProductNormalized(vn2);
		double s34 = (vn3.isZero() || vn4.isZero()) ? 0.0 : vn3.scalarProductNormalized(vn4);

		if( s12 > 0.0){
			split_faces.add( new MeshTriangle3d( points[0], points[1], points[2], blocks[0] ) );
			split_faces.add( new MeshTriangle3d( points[0], points[2], points[3], blocks[0] ) );
		}else{
			assert( s34 >= 0.0 );
			split_faces.add( new MeshTriangle3d( points[0], points[1], points[3], blocks[0] ) );
			split_faces.add( new MeshTriangle3d( points[1], points[2], points[3], blocks[0] ) );
		}
	}else{
		DataVector<DPoint3d> poly(count);
		for(int i = 0; i < count; i++)
			poly.add( points[i]->getCoordinates() );

		//if(true){
		//	MeshViewSet* set = new MeshViewSet();
		//	for(int i = 0; i < count; i++){
		//		set->addPoint(poly[i], 0, i);
		//		set->addEdge(poly[i], poly[(i+1)%count], 0);
		//	}
		//	SHOW_MESH("poly to split", set);
		//}

		MeshContainer2d* mesh = MeshGenerator2d::triangulatePoly(poly);
		assert(mesh != nullptr);

		int tct = mesh->getElementsCount();
		assert(tct == count-2);
		for(int i = 0; i < tct; i++){
			MeshTriangle2d* triangle = (MeshTriangle2d*) mesh->getElementAt(i);
			int ids[] = {
				triangle->getPoint(0)->getIntTag(TagExtended::TAG_ID, -1),
				triangle->getPoint(1)->getIntTag(TagExtended::TAG_ID, -1),
				triangle->getPoint(2)->getIntTag(TagExtended::TAG_ID, -1) 
			};
			assert( ids[0] >= 0 );
			assert( ids[1] >= 0 );
			assert( ids[2] >= 0 );
			split_faces.add( new MeshTriangle3d(points[ids[0]], points[ids[1]], points[ids[2]], blocks[0]) );
		}

		//if(true){
		//	MeshViewSet* set = new MeshViewSet();
		//	for(int i = 0; i < count; i++)
		//		set->addPoint(points[i]->getCoordinates(), 0, i);
		//	for(int i = 0; i < tct; i++){
		//		set->addFaceWithEdges( split_faces[i]);
		//	}
		//	SHOW_MESH("triangulated poly", set);
		//}

		delete mesh;
	}

	for(size_t i = sfct_start; i < split_faces.countInt(); i++)
		split_faces[i]->copyAllExtraData(this);

	return true;
}

/// Split face into triangles
bool MeshFace::splitToTriangles( DataVector< DTriangle3d> & triangles ) const
{
	assert(count > 2);
	if(count < 3) return true; // ie. wasn't able to fit a plane...

	if(count == 3){
		triangles.add( DTriangle3d( points[0]->getCoordinates(), 
			points[1]->getCoordinates(), points[2]->getCoordinates() ) );
	}else if( count == 4) {
		// like quad...
		const DPoint3d& dpt0 = points[0]->getCoordinates();
		const DPoint3d& dpt1 = points[1]->getCoordinates();
		const DPoint3d& dpt2 = points[2]->getCoordinates();
		const DPoint3d& dpt3 = points[3]->getCoordinates();

		const DVector3d vn1 = DVector3d::crossProduct(dpt0, dpt1, dpt2);
		const DVector3d vn2 = DVector3d::crossProduct(dpt0, dpt2, dpt3);
		const DVector3d vn3 = DVector3d::crossProduct(dpt0, dpt1, dpt3);
		const DVector3d vn4 = DVector3d::crossProduct(dpt1, dpt2, dpt3);

		double s12 = (vn1.isZero() || vn2.isZero()) ? 0.0 : vn1.scalarProductNormalized(vn2);
		double s34 = (vn3.isZero() || vn4.isZero()) ? 0.0 : vn3.scalarProductNormalized(vn4);

		if( s12 > 0.0){
			triangles.add( DTriangle3d( points[0]->getCoordinates(), 
				points[1]->getCoordinates(), points[2]->getCoordinates() ) );
			triangles.add( DTriangle3d(	points[0]->getCoordinates(), 
				points[2]->getCoordinates(), points[3]->getCoordinates() ) );
		}else{
			assert( s34 >= 0.0 );
			triangles.add( DTriangle3d(	points[0]->getCoordinates(), 
				points[1]->getCoordinates(), points[3]->getCoordinates() ) );
			triangles.add( DTriangle3d(	points[1]->getCoordinates(), 
				points[2]->getCoordinates(), points[3]->getCoordinates() ) );
		}
	}else{
		DataVector<DPoint3d> poly(count);
		for(int i = 0; i < count; i++)
			poly.add( points[i]->getCoordinates() );

		MeshContainer2d* mesh = MeshGenerator2d::triangulatePoly(poly);
		assert(mesh != nullptr);

		int tct = mesh->getElementsCount();
		assert(tct == count-2);
		for(int i = 0; i < tct; i++){
			MeshTriangle2d* triangle = (MeshTriangle2d*) mesh->getElementAt(i);
			int ids[] = {
				triangle->getPoint(0)->getIntTag(TagExtended::TAG_ID, -1),
				triangle->getPoint(1)->getIntTag(TagExtended::TAG_ID, -1),
				triangle->getPoint(2)->getIntTag(TagExtended::TAG_ID, -1) 
			};
			assert( ids[0] >= 0 );
			assert( ids[1] >= 0 );
			assert( ids[2] >= 0 );
			triangles.add( DTriangle3d(points[ids[0]]->getCoordinates(), 
				points[ids[1]]->getCoordinates(), points[ids[2]]->getCoordinates() ) );
		}

		delete mesh;
	}

	return true;
}

/// switch orientation of the face (direction of vertices)
void MeshFace::switchOrientation()
{
	// blocks
	//MeshBlock* btemp = blocks[0];
	//blocks[0] = blocks[1];
	//blocks[1] = btemp;

	// points
	int pct = getPointCount();
	for(int i = 0; i < pct/2; i++){
		MeshPoint3d* ptemp = points[i];
		points[i] = points[pct-i-1];
		points[pct-i-1] = ptemp;
	}

	// edges
	for(int i = 0; i < (pct-1)/2; i++){
		MeshEdge3d* etemp = edges[i];
		edges[i] = edges[pct-i-2];
		edges[pct-i-2] = etemp;
	}
}

/// index of neighbouring face
int MeshFace::getNeighbourIndex(const MeshFace* other_face) const
{
	for(int i = 0; i < count; i++){
		if( edges[i]->getFaceCount() == 2 && 
			edges[i]->getOtherFace(this) == other_face) return i;
	}
	return -1;
}

/// count base plane and 2d bounding box
bool MeshFace::countPlaneBBox( DPlane& plane, DRect& box2d ) const
{
	assert( count > 2 );
	// plane...
	plane.vn = getNormalVector();
	if( plane.vn.isZero() ) return false;
	plane.p0 = points[0]->getCoordinates();
	plane.vn.orthonormalVectors( plane.e0, plane.e1 );
	// bbox
	box2d.clear();
	for(int i = 0; i < count; i++)
		box2d.addPoint( plane.projectToPlane( points[i]->getCoordinates() ) );

	return true;
}

/// get all points for the given set of faces
void MeshFace::getMeshPointsForFaces( 
	const DataVector<MeshFace*> & mfaces, 
	DataVector< MeshPoint3d* > & mpoints )
{
	size_t mct = mfaces.countInt();
	assert( mct > 0 );

	int pct = 0;
	for(size_t i = 0; i < mct; i++) pct += mfaces[i]->getPointCount();

	mpoints.clear();
	mpoints.prepare( pct );

	for(size_t i = 0; i < mct; i++){
		const MeshFace* mf = mfaces[i];
		int fpct = mf->getPointCount();
		for(int j = 0; j < fpct; j++)
			mpoints.addIfNew( mf->getPoint(j) );
	}
}

/// count inverted-rank for this face
int MeshFace::getInvertedRank(int * inv_index) const
{
	DVector3d dn0, dn1;
	if( !checkAndGetNormalVector(dn0) ) return -1;
	
	int inverted_rank = 0;
	for(int j = 0; j < count; j++ ) {
		if( edges[j]->isBorder() ) continue;
		MeshFace* other_face = edges[j]->getOtherFace(this);
		if( other_face->checkAndGetNormalVector(dn1) ) {
			double sp = dn0.scalarProduct(dn1);
			if( sp < -0.9 ) { 
				inverted_rank += SP_BAD; 
				if(inv_index != nullptr) *inv_index = j; 
			}else if( sp > 0.9 ) inverted_rank += SP_GOOD;
			else inverted_rank += SP_OK;
		}
	}

	return inverted_rank;
}


DVector3d MeshFace::getOuterNormalForEdge(int i, double nvalue)
{
	DVector3d n = getNormalVector();
	DVector3d v_edge = ( edges[i]->getFaceIndex(this) == 0 ) ?
		(points[1]->getCoordinates() - points[0]->getCoordinates()) :
		(points[0]->getCoordinates() - points[1]->getCoordinates());
	assert( ! v_edge.isZero() );
	return v_edge.crossProduct(n).normalized( nvalue );
}
