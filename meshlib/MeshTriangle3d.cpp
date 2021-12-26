// MeshTriangle3d.cpp: implementation of the MeshTriangle3d class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshTriangle3d.h"
#include "MeshEdge3d.h"
#include "DTriangle.h"
#include "Metric3dContext.h"
#include "SurfaceParametric.h"
#include "MeshViewSet.h"

MeshTriangle3d::MeshTriangle3d(MeshPoint3d *p1, MeshPoint3d *p2, MeshPoint3d *p3, MeshBlock* bounded_block)
	: MeshFace(3, triangle_edges, triangle_points)
{
	blocks[0] = bounded_block;
	points[0] = p1;
	points[1] = p2;
	points[2] = p3;

//	double area = pt1.triangleArea(pt2, pt3);
//	if(area < 0.0){	// Zmiana orientacji
//		points[1] = p3;	
//		points[2] = p2;
//	}
	attachToEdges();
}

/// Creates and returns a copy of this face (+ whole connectivity)
MeshFace* MeshTriangle3d::clone() const
{
	MeshTriangle3d* face = new MeshTriangle3d(points[0], points[1], points[2], blocks[0]);
	face->blocks[1] = blocks[1];
	
	face->copyBorderFlagsFrom(this);
	face->copyAllTags(this);
	face->copySurfaceData(this);
	return face;
}

MeshTriangle3d::~MeshTriangle3d()
{
	detachFromEdges();
}

/// Area of a triangular face
double MeshTriangle3d::area() const
{
	return DTriangle3d::area(
		points[0]->getCoordinates(),
		points[1]->getCoordinates(),
		points[2]->getCoordinates());
}

/// Area of a triangular face
double MeshTriangle3d::area(Metric3dContext& mc) const
{
	return DMTriangle3d::area(
		points[0]->getMetricCoordinates(mc),
		points[1]->getMetricCoordinates(mc),
		points[2]->getMetricCoordinates(mc));
}

/// Alpha quality coefficient of a triangular face
double MeshTriangle3d::alphaQuality() const
{
	return DTriangle3d::alphaQuality(
		points[0]->getCoordinates(),
		points[1]->getCoordinates(),
		points[2]->getCoordinates());
}

MeshEdge3d* MeshTriangle3d::swapWithNeighbour(Metric3dContext& mc, 
		SurfaceConstPtr surface, int ind, bool improve_only, 
		bool local_metric, TagExtended::TagType side_tag)
{
	if(edges[ind]->isBorder() || edges[ind]->getIntTag(TagExtended::TAG_SIDE_EDGE) > 0) return nullptr;
	MeshTriangle3d* triangle = (MeshTriangle3d*)edges[ind]->getOtherFace(this);
	if(!triangle || triangle->getEdgeCount() != 3) return nullptr;

	MeshPoint3d* p0 = edges[ind]->getMeshPoint(0);
	MeshPoint3d* p1 = edges[ind]->getMeshPoint(1);
	MeshPoint3d* p2 = this->getOtherPoint(p0, p1);
	MeshPoint3d* p3 = triangle->getOtherPoint(p0, p1);

	if(p2->getEdgeToPoint(p3)) return nullptr; // such edge already exists...

	if(local_metric){
		const DPoint3d middle = DPoint3d::average(p0->getCoordinates(), p1->getCoordinates(), 
			p2->getCoordinates(), p3->getCoordinates());
		mc.countMetricAtPoint(middle);
	}

	DMPoint3d dmp0 = p0->getMetricCoordinates(mc);
	DMPoint3d dmp1 = p1->getMetricCoordinates(mc);
	DMPoint3d dmp2 = p2->getMetricCoordinates(mc);
	DMPoint3d dmp3 = p3->getMetricCoordinates(mc);

	double q1_pre = abs(DMTriangle3d::alphaQuality(dmp0, dmp1, dmp2));
	double q2_pre = abs(DMTriangle3d::alphaQuality(dmp0, dmp1, dmp3));

	double q1_after = abs(DMTriangle3d::alphaQuality(dmp1, dmp2, dmp3));
	double q2_after = abs(DMTriangle3d::alphaQuality(dmp0, dmp2, dmp3));


	if(improve_only){
		if(std::min(q1_after, q2_after) <= std::min(q1_pre, q2_pre)) return nullptr; // current min quality is better than after swap
	}

	// ... make swap
	this->switchPointsWithEdges(p0, p3);
	triangle->switchPointsWithEdges(p1, p2);

	// ... check
	bool acceptable = false;
	if( surface != nullptr )
		acceptable = this->valid( surface ) && triangle->valid( surface );
	else
		acceptable = this->validDirect(mc) && triangle->validDirect(mc);

	if(!acceptable){
		// reverse operation
		this->switchPointsWithEdges(p3, p0);
		triangle->switchPointsWithEdges(p2, p1);
	}

	return acceptable ? p2->getEdgeToPoint(p3) : nullptr;
}

/// returns shape quality for a face
double MeshTriangle3d::getShapeQuality(Metric3dContext& mc) const
{
	return DMTriangle3d::alphaQuality(
			points[0]->getMetricCoordinates(mc),
			points[1]->getMetricCoordinates(mc),
			points[2]->getMetricCoordinates(mc));
}

// check, whether this face is valid with respect to the given surface
bool MeshTriangle3d::valid( SurfaceConstPtr surface) const
{
	//static int counter = 0;
	//counter++;
	DPoint2d lspt0, lspt1, lspt2;
	bool d0 = points[0]->checkAndGetSurfaceParam( surface, lspt0 );
	bool d1 = points[1]->checkAndGetSurfaceParam( surface, lspt1 );
	bool d2 = points[2]->checkAndGetSurfaceParam( surface, lspt2 );

	//DPoint2d spt0 = lspt0, spt1 = lspt1, spt2 = lspt2;
	//bool b0 = surface->withinDomain( points[0]->getCoordinates(), spt0 );
	//bool b1 = surface->withinDomain( points[1]->getCoordinates(), spt1 );
	//bool b2 = surface->withinDomain( points[2]->getCoordinates(), spt2 );

	//if( d0 != b0 || d1 != b1 || d2 != b2 ) {
	//	MeshViewSet* set = new MeshViewSet;
	//	set->addFaceWithEdgesAndPoints(this);
	//	surface->drawDomain( set );
	//	set->addInfo( to_string( points[0]->getIndex() ),
	//		string( d0 ? "d-tak" : "d-nie") + " | " + string( b0 ? "b-tak" : "b-nie") );
	//	set->addInfo( to_string( points[1]->getIndex() ),
	//		string( d1 ? "d-tak" : "d-nie") + " | " + string( b1 ? "b-tak" : "b-nie") );
	//	set->addInfo( to_string( points[2]->getIndex() ),
	//		string( d2 ? "d-tak" : "d-nie") + " | " + string( b2 ? "b-tak" : "b-nie") );
	//	SHOW_MESH("MeshTriangle3d::valid inconsistent...", set);
	//}

	return d0 && d1 && d2 && DTriangle2d::valid( lspt0, lspt1, lspt2 );
}

/// returns shape quality for a face
double MeshTriangle3d::getShapeQuality() const
{
	return DTriangle3d::alphaQuality(
			points[0]->getCoordinates(),
			points[1]->getCoordinates(),
			points[2]->getCoordinates());
}

/// Split face into (mesh) triangles
bool MeshTriangle3d::splitToTriangles( DataVector< MeshFace* > & split_faces ) const
{
	assert(false); // shouldn't be needed anyway...
	MeshFace* f = new MeshTriangle3d( points[0], points[1], points[2], blocks[0] );
	f->copyAllExtraData(this);
	split_faces.add( f );
	return true;
}

/// Split face into triangles
bool MeshTriangle3d::splitToTriangles( DataVector< DTriangle3d > & triangles ) const
{
	triangles.add( DTriangle3d( 
		points[0]->getCoordinates(), 
		points[1]->getCoordinates(), 
		points[2]->getCoordinates() ) );
	return true;
}
