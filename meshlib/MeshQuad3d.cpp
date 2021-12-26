// MeshQuad3d.cpp: implementation of the MeshQuad3d class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshQuad3d.h"
#include "MeshEdge3d.h"
#include "MeshTriangle3d.h"
#include "DTriangle.h"
#include "MeshViewSet.h"
#include "DQuad.h"

MeshQuad3d::MeshQuad3d(MeshPoint3d *p1, MeshPoint3d *p2, MeshPoint3d *p3, MeshPoint3d *p4, 
	MeshBlock* bounded_block) : MeshFace(4, quad_edges, quad_points)
{
	blocks[0] = bounded_block;
	points[0] = p1;
	points[1] = p2;
	points[2] = p3;
	points[3] = p4;

//	double area1 = pt1.triangleArea(pt2, pt3);
//	double area2 = pt3.triangleArea(pt4, pt1);
//	if(area1 * area2 < 0){
		// punkty nie tworz¹ wypuk³ego czworok¹ta ...
//	}
//	if(area1 < 0.0 && area2 < 0.0){	// Zmiana orientacji
//		points[1] = p4;
//		points[2] = p3;
//		points[3] = p2;
//	}
	attachToEdges();
}

/// Creates and returns a copy of this face (+ whole connectivity)
MeshFace* MeshQuad3d::clone() const
{
	MeshQuad3d* face = new MeshQuad3d(points[0], points[1], points[2], points[3], blocks[0]);
	face->blocks[1] = blocks[1];
	face->copyAllExtraData(this);
	return face;
}

MeshQuad3d::~MeshQuad3d()
{
	detachFromEdges();
}

/// remove point, return resulting face (the same, or a new one)
MeshFace* MeshQuad3d::removePoint(const MeshPoint3d* point) { 
	DataVector<MeshPoint3d*> tpoints(4);
	for(int i = 0; i < count; i++)
		if( points[i] != point) tpoints.add(points[i]);
	assert(tpoints.countInt() == 3);
	
	MeshFace* new_face = new MeshTriangle3d( tpoints[0], tpoints[1], tpoints[2], blocks[0] ); 
	new_face->copyAllExtraData(this);
	return new_face;
}

// check, whether this face is valid with respect to the given surface
bool MeshQuad3d::valid( SurfaceConstPtr surface ) const
{
	return DQuad2d::valid( 
		getPoint(0)->getLocalSurfaceParam( surface ),
		getPoint(1)->getLocalSurfaceParam( surface ), 
		getPoint(2)->getLocalSurfaceParam( surface ), 
		getPoint(3)->getLocalSurfaceParam( surface ) );
}

/// returns shape quality for a face
double MeshQuad3d::getShapeQuality(Metric3dContext& mc) const
{
	return DQuad3d::shapeQuality(
						points[0]->getMetricCoordinates(mc).toRealSpace(),
						points[1]->getMetricCoordinates(mc).toRealSpace(),
						points[2]->getMetricCoordinates(mc).toRealSpace(),
						points[3]->getMetricCoordinates(mc).toRealSpace());
}

/// returns shape quality for a face
double MeshQuad3d::getShapeQuality() const
{
	return DQuad3d::shapeQuality(
				points[0]->getCoordinates(), points[1]->getCoordinates(),
				points[2]->getCoordinates(), points[3]->getCoordinates());

	//double alpha[4];
	//int i, j;
	//for(i = 0; i < 4; i++){
	//	double aq = DMTriangle3d::alphaQuality(
	//		points[i]->getMetricCoordinates(mc), 
	//		points[(i+1)%4]->getMetricCoordinates(mc), 
	//		points[(i+2)%4]->getMetricCoordinates(mc));

	//	// Sorting of alpha-values of 4 triangles
	//	for(j = i-1; j >= 0; j--){
	//		if(alpha[j] > aq) alpha[j+1] = alpha[j];
	//		else break;
	//	}
	//	alpha[j+1] = aq;
	//}

	//double denominator = alpha[2] * alpha[3];
	//if(abs(denominator) > SMALL_NUMBER){
	//	return alpha[0] * alpha[1] / denominator;
	//}else{
	//	return 0.0;
	//}
}

/// Split face into (mesh) triangles
bool MeshQuad3d::splitToTriangles( DataVector< MeshFace* > & split_faces ) const
{
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

	MeshFace *f0 = nullptr, *f1 = nullptr;

	if( s12 > 0.0){
		f0 = new MeshTriangle3d( points[0], points[1], points[2], blocks[0] );
		f1 = new MeshTriangle3d( points[0], points[2], points[3], blocks[0] );
	}else{
		assert( s34 >= 0.0 );
		f0 = new MeshTriangle3d( points[0], points[1], points[3], blocks[0] );
		f1 = new MeshTriangle3d( points[1], points[2], points[3], blocks[0] );
	}

	f0->copyAllExtraData(this);
	f1->copyAllExtraData(this);

	split_faces.add( f0 );
	split_faces.add( f1 );

	return true;
}

/// Split face into triangles
bool MeshQuad3d::splitToTriangles( DataVector< DTriangle3d > & triangles ) const
{
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

	MeshFace *f0 = nullptr, *f1 = nullptr;

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

	return true;
}
