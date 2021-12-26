/////////////////////////////////////////////////////////////////////////////
// SurfaceQuadricQTree.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "SurfaceQuadricQTree.h"
#include "DPlanarQuadric.h"
#include "MeshFace.h"
#include "DLeastSquaresFitting.h"
#include "MeshContainer3dSurface.h"
#include "GeometricPredicates.h"
#include "Metric3dContext.h"
#include "MeshViewSet.h"
#include "Metric2dContext.h"
#include "ControlSpace2dIdentity.h"
#include "ControlSpace2dProjected.h"

QuadricLeaf::QuadVertexWhich	QuadricLeaf::mid_to_nodes[4][2] = 
			{{VX2D_SW, VX2D_SE}, {VX2D_SE, VX2D_NE},
			{VX2D_NW, VX2D_NE}, {VX2D_SW, VX2D_NW}};
QuadricLeaf::QuadMidVertexWhich	QuadricLeaf::nodes_to_mid[4][2] = 
			{{MVX2D_DOWN, MVX2D_LEFT}, {MVX2D_DOWN, MVX2D_RIGHT},
			{MVX2D_UP, MVX2D_LEFT}, {MVX2D_UP, MVX2D_RIGHT}};
QuadricLeaf QuadricLeaf::empty;


QuadricLeaf::QuadricLeaf() :m_level(0), m_valid(false)
{ 
	m_leaves[0] = nullptr; 
	m_neighbours[0] = m_neighbours[1] = m_neighbours[2] = m_neighbours[3] = nullptr;
	m_mid_vertices[0] = m_mid_vertices[1] = m_mid_vertices[2] = m_mid_vertices[3] = nullptr;
}

QuadricLeaf::QuadricLeaf(double mid_x, double mid_y, double dx, double dy, int level) 
	: m_middle(mid_x, mid_y), m_dl(dx, dy), m_level(level), m_valid(false)
{
	m_leaves[0] = nullptr; 
	m_neighbours[0] = m_neighbours[1] = m_neighbours[2] = m_neighbours[3] = nullptr;
	m_mid_vertices[0] = m_mid_vertices[1] = m_mid_vertices[2] = m_mid_vertices[3] = nullptr;

}

QuadricLeaf::~QuadricLeaf()
{ 
	if(m_leaves[0]) 
		for(int i = 0; i < 4; i++) delete m_leaves[i]; 
}

const QuadricLeaf* QuadricLeaf::findLastLeaf(const DPoint2d& pt) const
{
	if(m_leaves[LF_FIRST] == nullptr) return this;
	if(pt.x < m_middle.x) 
		if(pt.y < m_middle.y) return m_leaves[LF_SW]->findLastLeaf(pt);
		else return m_leaves[LF_NW]->findLastLeaf(pt);
	else if(pt.y < m_middle.y) return m_leaves[LF_SE]->findLastLeaf(pt);
		else return m_leaves[LF_NE]->findLastLeaf(pt);
}

const DPoint3d QuadricLeaf::getPoint(const DPoint2d& pt, SurfaceConstPtr surface) const
{
	assert(!isSplit());
	assert(m_valid);

	double tx = (pt.x - m_middle.x + 2*m_dl.x) / (4*m_dl.x);
	double ty = (pt.y - m_middle.y + 2*m_dl.y) / (4*m_dl.y);

	double coeff[4] = {(1-tx)*(1-ty), // VX2D_SW
		tx *(1-ty), // VX2D_SE
		(1-tx)* ty, // VX2D_NW
		tx * ty};	// VX2D_NE
	DVectorN<6> ave_vq = m_vertices[0]->vq * coeff[0];
	for(int i = 1; i < 4; i++)
		ave_vq += m_vertices[i]->vq * coeff[i];
	return DQuadricOnSurface::getPoint(pt, surface, ave_vq);
}

bool QuadricLeaf::withinParamRange( const DPoint2d& pt ) const
{
	assert(!isSplit());

	return m_valid;

	//return	
	//	m_vertices[0]->valid() || m_vertices[1]->valid() || 
	//	m_vertices[2]->valid() || m_vertices[3]->valid();
}

const DVector3d QuadricLeaf::getDerivative(int deriv, const DPoint2d& pt, SurfaceConstPtr surface) const
{
	assert(!isSplit());
	assert(m_valid);

	double tx = (pt.x - m_middle.x + 2*m_dl.x) / (4*m_dl.x);
	double ty = (pt.y - m_middle.y + 2*m_dl.y) / (4*m_dl.y);

	double coeff[4] = {(1-tx)*(1-ty), // VX2D_SW
		tx *(1-ty), // VX2D_SE
		(1-tx)* ty, // VX2D_NW
		tx * ty};	// VX2D_NE
	DVectorN<6> ave_vq = m_vertices[0]->vq * coeff[0];
	for(int i = 1; i < 4; i++)
		ave_vq += m_vertices[i]->vq * coeff[i];
	return DQuadricOnSurface::getDerivative(deriv, pt, surface, ave_vq);
}

SurfaceQuadricQTree::SurfaceQuadricQTree(
	SurfacePtr base_surface,
	const DRect& box, int nx, int ny ) 
		: m_base_surface(base_surface), m_box(box), m_nx(nx), m_ny(ny),
			m_grid( nx, ny, QuadricLeaf::empty )
{
	double dx = box.getDX() / nx;
	double dy = box.getDY() / ny;

	// add regular grid vertices
	QuadricNode qv(box.x0, box.y0);
	for(int j = 0; j <= ny; j++, qv.coord.y += dy){
		qv.coord.x = box.x0;
		for(int i = 0; i <= nx; i++, qv.coord.x += dx){
			m_grid_vertices.add(qv);
		}
	}

	// init regular grid boxes with coordinates and vertices indices
	DVector2d dpt(dx/4, dy/4);
	DPoint2d pt(m_box.x0 + dx / 2, m_box.y0 + dy / 2);
	for(int j = 0, ix = 0; j < m_ny; j++, ix++){
		for(int i = 0; i < m_nx; i++, ix++){
			QuadricLeaf& qleaf = m_grid.get(j,i);
			qleaf.setCoordinates(pt, dpt);
			pt.x += dx;
			qleaf.setNeighbours(
				((j>0)		? (&m_grid.get(j-1,i)) : nullptr),  // down
				((i<m_nx-1) ? (&m_grid.get(j,i+1)) : nullptr),  // right
				((j<m_ny-1) ? (&m_grid.get(j+1,i)) : nullptr),  // up
				((i>0)		? (&m_grid.get(j,i-1)) : nullptr)); // left
			qleaf.setVertices(
				&(m_grid_vertices[ix]), 
				&(m_grid_vertices[ix+1]), 
				&(m_grid_vertices[ix+m_nx+1]), 
				&(m_grid_vertices[ix+m_nx+2]));
		}
		pt.y += dy;
		pt.x = m_box.x0 + dx / 2;
	}
}

#define CONTAINING_EPS SMALL_NUMBER

bool QuadricNode::checkFaceDist( const DMPoint2d &qpt, 
	const DMPoint2d &fa, const DMPoint2d &fb, const DMPoint2d &fc,
	MeshFace* face)
{
	double dist2 = DTriangle2d::distance2ToPoint( qpt, fa, fb, fc);

	if( dist2 < min_dist2 ) {
		min_dist2 = dist2;
		fnearest = face;
		if( min_dist2 < CONTAINING_EPS ) {
			setQuality( AQ_VALID_MAX );// already found
			return true;
		}
	}

	return false;
}

bool SurfaceQuadricQTree::initializeWithFaces( 
	Metric3dContext &mc, 
	MeshContainer3dSurface* mesh, 
	const DataVector<MeshFace*> & sfaces,
	double tolerance )
{
	int sfct = sfaces.countInt();
	for(int i = 0; i < sfct; i++)
		sfaces[i]->setIntTag( TagExtended::TAG_SQQTREE, 1 );

	auto cs3d = mesh->getControlSpace();
	CS2dConstPtr csp(new ControlSpace2dProjected( cs3d, m_base_surface ));
	Metric2dContext mc2d( csp );

	struct FaceParams {
		FaceParams(bool v = true) : valid(v) {}
		bool valid;
		DPoint2d params[3];
	};

	DataVector< FaceParams > fparams( sfct, FaceParams(true) );

	for(int i = 0; i < sfct; i++) {
		MeshFace* face = sfaces[i];

		DPoint2d fmid = DPoint2d::zero;
		DRect fbox;
		FaceParams& fp = fparams[i];
		for(int k = 0; fp.valid && k < 3; k++) {
			fp.valid = face->getPoint(k)->checkAndGetSurfaceParam( 
				m_base_surface, fp.params[k] );
			fbox.addPoint( fp.params[k] );
			fmid.add( fp.params[k], 1.0/3.0 );
		}
		if( !fp.valid) continue;

		mc2d.countMetricAtPoint( fmid );

		const DMPoint2d a = mc2d.transformPStoMS( fp.params[0] );
		const DMPoint2d b = mc2d.transformPStoMS( fp.params[1] );
		const DMPoint2d c = mc2d.transformPStoMS( fp.params[2] );

		int ix0 = (int) ceil(m_nx * (fbox.x0 - m_box.x0) / (m_box.x1 - m_box.x0));
		int iy0 = (int) ceil(m_ny * (fbox.y0 - m_box.y0) / (m_box.y1 - m_box.y0));
		int ix1 = (int) (m_nx * (fbox.x1 - m_box.x0) / (m_box.x1 - m_box.x0));
		int iy1 = (int) (m_ny * (fbox.y1 - m_box.y0) / (m_box.y1 - m_box.y0));

		for(int iy = iy0; iy <= iy1; iy++){
			int k = (m_nx+1) * iy;
			for(int ix = ix0; ix <= ix1; ix++){
				int kq = k + ix;
				QuadricNode& qnode = m_grid_vertices[kq];
				if( !qnode.valid() ) {
					qnode.checkFaceDist( 
						mc2d.transformPStoMS( qnode.coord ),
						a, b, c, face);
				}
			}
		}
	}

	// testowo: dla ka¿dego wêz³a -> znaleŸæ najbli¿szy trójk¹t?

	for(size_t i = 0; i <  m_grid_vertices.countInt(); i++){
		QuadricNode& qnode = m_grid_vertices[i];
		if( qnode.valid() ) continue;
		mc2d.countMetricAtPoint( qnode.coord );
		DMPoint2d qpt = mc2d.transformPStoMS( qnode.coord );

		for(int j = 0; j < sfct; j++){
			FaceParams& fp = fparams[j];
			if( !fp.valid) continue;

			if( qnode.checkFaceDist( qpt, 
					mc2d.transformPStoMS( fp.params[0] ),
					mc2d.transformPStoMS( fp.params[1] ),
					mc2d.transformPStoMS( fp.params[2] ),
					sfaces[j]) ) break;
		}
	}

	m_grid.forEach( [ this, &mc2d, &sfaces ] (QuadricLeaf& ql ) {
		ql.adaptForFaces( m_grid_vertices, m_base_surface, mc2d, sfaces );
	} );

	if (true) {
		MeshViewSet *set = new MeshViewSet();
		// 3d version (or use parametric one?)
		sfaces.forEach([&set](MeshFace* f) { set->addEdges(f); });
		m_grid.forEach([&set,this](QuadricLeaf& ql) { ql.addToViewSet(set, m_base_surface); });
		m_grid_vertices.forEach([&set, this](QuadricNode& qn) {
			DPoint3d pt3d = m_base_surface->getPoint(qn.getCooord());
			set->addPoint(pt3d, (qn.valid(AQ_VALID_MAX)?2:(qn.valid()?1:0)));
			//DVector3d vn = m_base_surface.ptr()->getNormalVector(qn.getCooord());
			//set->addEdge(pt3d, pt3d + vn, 1);
		});
		SHOW_MESH("Adapted QuadricQTree", set);
	}

	//sfaces.forEach( [this,mesh] (MeshFace* face ) {
	//
	//	int fpct = face->getPointCount();
	//
	//	// gather (some) layer of vertices (not crossing border edges/points)
	//	DataVector<DPoint3d> points;
	//	if(mesh->gatherLayeredVerticesTopological(face, points, 3, 
	//		false, 10, nullptr, TagExtended::TAG_SQQTREE, 1)) 
	//	{
	//		// approximate surface  -> z(1,x,y,xy,xx,yy)
	//		DPlanarQuadric pquadric;
	//		double max_dist = DLeastSquaresFitting::fitPlanarQuadric(points, m_base_plane, pquadric);
	//	}
	//
	//} );

//DataVector< DPoint2d > param_points( 2*sfct );
//	DataHashTableKeyValue< MeshPoint3d*, int > hppoints( 4*sfct, nullptr );
//	DRect sbox;
	//for(int j = 0; j < sfct; j++ ) {
	//	MeshFace* face = sfaces[j];
	//	assert( face->getPointCount() == 3 );
	//	for(int k = 0; k < 3; k++ ) {
	//		MeshPoint3d* mpoint = face->getPoint(k);
	//		int pPI = hppoints.getValue( mpoint, -1 );
	//		if( pPI == -1 ) {
	//			pPI = param_points.add(
	//				m_base_surface.ptr()->getParameters( mpoint->getCoordinates() ) );
	//			hppoints.insert( mpoint,  pPI );
	//			sbox.addPoint( param_points[pPI] );
	//		}
	//	}
	//}
	//double eps = 1e-3 * sbox.getDiameter();
	
	DataVector< DQuadricOnSurface > face_quadrics( sfct );
	DataHashTableKeyValue< MeshFace*,int > hface_quadrics( (unsigned int)(2*sfct), nullptr );

	for(size_t i = 0; i < m_grid_vertices.countInt(); i++ ) {
		// face for qnode.coord
		QuadricNode& qnode = m_grid_vertices[i];
		if (!qnode.valid()) continue;

		// metric for qnode.coord
		// DPoint3d sdpt = m_base_surface.ptr()->getPoint(qnode.coord);
		//mc.countMetricAtPoint( sdpt );
		//ControlDataMatrix2d scdm = m_base_surface.ptr()->projectTransformationTensor( qnode.coord, 
		//	cs3d->getMetricAtPoint( sdpt ) );
		// ... or create ControlSpace2dIdentity(scdm) 
		//ControlSpace2dIdentity csi( scdm, m_base_surface.ptr() ); 
		//Metric2dContext mc2d( &csi );
		//DMetric2d dm2d(scdm, m_base_surface, qnode.coord);

		// face for qnode.coord
		// TODO - triangle walking from seed ...
		//int cfi = -1;
		//bool found_containing = false;
		//double ab, bc, ca;
		//for(int j = 0; (cfi < 0) && (j < sfct); j++ ) {
		//	MeshFace* face = sfaces[j];
		//	assert( face->getPointCount() == 3 );
		//	
		//	const DPoint2d& a = param_points[ hppoints.getValue( face->getPoint(0), -1 ) ];
		//	const DPoint2d& b = param_points[ hppoints.getValue( face->getPoint(1), -1 ) ];
		//	const DPoint2d& c = param_points[ hppoints.getValue( face->getPoint(2), -1 ) ];
		//
		//	ab = GeometricPredicates::orient2d( a, b, qnode.coord);
		//	bc = GeometricPredicates::orient2d( b, c, qnode.coord);
		//	ca = GeometricPredicates::orient2d( c, a, qnode.coord);
		//
		//	if( (ab > -tolerance) && (bc > -tolerance) && (ca > -tolerance) )
		//		cfi = j; // index of containing face
		//}

		//if( cfi < 0 ) { // search for nearest triangle
		//double min_dist2 = LARGE_NUMBER;
		//if( true ) { // search for nearest triangle
		//	//mc.countMetricAtPoint( m_base_surface->getPoint( qnode.coord) );
		//	double eps_dist2 = 1.0;
		//	mc2d.countMetricAtPoint( qnode.coord );
		//	DMPoint2d qpt = mc2d.transformPStoMS( qnode.coord );
		//	for(int j = 0; j < sfct; j++ ) {
		//		MeshFace* face = sfaces[j];
		//		assert( face->getPointCount() == 3 );
		//		DPoint2d params[3];
		//		bool params_valid = true;
		//		for(int k = 0; params_valid && k < 3; k++)
		//			params_valid = face->getPoint(k)->checkAndGetSurfaceParam( m_base_surface.ptr(), params[k] );
		//		if( !params_valid) continue;
		//		const DMPoint2d a = mc2d.transformPStoMS( params[0] );
		//		const DMPoint2d b = mc2d.transformPStoMS( params[1] );
		//		const DMPoint2d c = mc2d.transformPStoMS( params[2] );
		//		double mdist2 = DTriangle2d::distance2ToPoint( qpt, a, b, c);
		//		if( mdist2 < min_dist2 ) {
		//			cfi = j;
		//			min_dist2 = mdist2;
		//			found_containing = (mdist2 <= 0.0);
		//			if(found_containing) break; // containing triangle found...
		//		}
		//	}
		//}

		MeshFace* cface = qnode.fnearest;
		assert(cface != nullptr);
		int cfi = hface_quadrics.getValue(cface, -1);

		LOG4CPLUS_INFO(MeshLog::logger_mesh, "QNode #" << i << " \tcfi=" << cfi << " \tmindist=" << std::sqrt(qnode.min_dist2));
		if( false ){
			MeshViewSet* set = new MeshViewSet;
			for(size_t ii = 0; ii < m_grid_vertices.countInt(); ii++){
				QuadricNode& qqnode = m_grid_vertices[ii];
				DPoint3d dpt = m_base_surface->getPoint(qqnode.coord);
				set->addPoint( dpt, (ii == i ) ? 3 : 1);
				set->addLabel( dpt, to_string( qqnode.approx_quality ) );
			}
			for(int jj = 0; jj < sfct; jj++){
				set->addFaceWithEdgesAndPoints( sfaces[jj],  (sfaces[jj] == cface) ? 3 : 1 );
			}
			set->addInfo("min_mdist", to_string(std::sqrt(qnode.min_dist2)) );
			set->addInfo("cfi", to_string(cfi));
			SHOW_MESH("SQQT::initializeWithFaces", set);
		}

		if( cfi == -1) {
			// if not already approximated, approximate pquadric
			// gather (some) layer of vertices (not crossing border edges/points)
			DataVector<DPoint3d> points;
			if(mesh->gatherLayeredVerticesTopological(cface, points, 3,
				false, 10, nullptr, TagExtended::TAG_SQQTREE, 1)) 
			{
				DQuadricOnSurface dquadric;
				// approximate surface  -> z(1,x,y,xy,xx,yy)
				double max_dist = DLeastSquaresFitting::fitQuadricOnSurface(points, 
							m_base_surface, dquadric);
				cfi = (int)face_quadrics.add( dquadric );
				hface_quadrics.insert(cface, cfi);

				if(true){
					MeshViewSet* set = new MeshViewSet;
					// draw m_base_surface
					createViewSetForBaseGrid(set, -1);
					// draw face ...
					set->addFace(cface);
					sfaces.forEach( [set] (MeshFace * face) { set->addEdges( face ); } );
					// draw local quadric-on-surface
					DQuadricOnSurface::createViewSetForPoints(set, points, 
						dquadric.base_surface, dquadric.vq, 3 );

					for (size_t i = 0; i < points.countInt(); i++) {
						//set->addLabel(points[i], "points-0");
						DPoint2d param = dquadric.base_surface->getParameters(points[i]);
						//DPoint3d test_pt = dquadric.getPoint(param);
						DVector3d vn = dquadric.getDerivative(DEquation::deriv_ds, param).crossProduct(
							dquadric.getDerivative(DEquation::deriv_dt, param)).normalized(0.25);
						set->addEdge(points[i], points[i] + vn, 1);
					}
						
					SHOW_MESH("local quadric-on-surface for quadric-qtree initialization", set);
				}
			}
		}

		// set pquadric, and approx-error and inside/outside
		if( cfi >= 0 ) {			
			qnode.vq = face_quadrics[cfi].vq;
			DPoint3d pq_point = face_quadrics[cfi].getPoint( qnode.coord );
			DPoint3d f_point = DPoint3d::zero;
			if( qnode.approx_quality == AQ_VALID_MAX ) { // i.e. is inside
				DPoint2d params[3];
				bool params_valid = true;
				for(int k = 0; params_valid && k < 3; k++)
					params_valid = sfaces[cfi]->getPoint(k)->checkAndGetSurfaceParam( 
																m_base_surface, params[k] );
				assert(params_valid);

				double ab = GeometricPredicates::orient2d( params[0], params[1], qnode.coord);
				double bc = GeometricPredicates::orient2d( params[1], params[2], qnode.coord);
				double ca = GeometricPredicates::orient2d( params[2], params[0], qnode.coord);
				double abc = ab+bc+ca;
				f_point.add( sfaces[cfi]->getPoint(0)->getCoordinates(),  bc/abc ); 
				f_point.add( sfaces[cfi]->getPoint(1)->getCoordinates(),  ca/abc ); 
				f_point.add( sfaces[cfi]->getPoint(2)->getCoordinates(),  ab/abc ); 

				mc.countMetricAtPoint( m_base_surface->getPoint( qnode.coord ) );
				double approx_dist = mc.transformRStoMS( pq_point - f_point ).length();
				qnode.approx_quality = (approx_dist > tolerance) ? 0.0 : (1.0 - approx_dist/tolerance);
			}
		}
		else {
			// leave vq as zeros -> revert to project surface...
			qnode.approx_quality = AQ_VALID_MIN;
		}
	}

	// extrapolate one layer ...
	// extrapolateGridLayer();

	for(int i = 0; i < sfct; i++)
		sfaces[i]->removeTag( TagExtended::TAG_SQQTREE );

	if( true ) {
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < sfct; i++)
			set->addFaceWithEdgesAndPoints( sfaces[i] );
		createViewSet( set );

		SHOW_MESH("quadric-on-surface qtree", set );
	}

	m_valid = true;
	return true;
}

bool SurfaceQuadricQTree::invertOrientation() 
{
	assert(false); // not valid, grid-nodes coordinates would have to be recalculated as well ...

	if( ! m_base_surface->invertOrientation() )
		return false;

	m_grid_vertices.forEach( [] ( QuadricNode& qn ) {
		DPlanarQuadric::invertSurfaceQuadricOrientation( qn.vq );
	} );

	return true;
}

void SurfaceQuadricQTree::extrapolateGridLayer()
{
	m_grid.forEach( [ ] ( QuadricLeaf & ql ) {
		int counter = 0;
		for(int i = QuadricLeaf::VX2D_FIRST; i <= QuadricLeaf::VX2D_LAST; i++ )
			if( ql.getVertex(i)->approx_quality < AQ_VALID_MIN ) counter++;
		if( counter > 0 && counter < 4 ) {
			DVectorN<6> pq(0.0);
			// gather from valid
			for(int i = QuadricLeaf::VX2D_FIRST; i <= QuadricLeaf::VX2D_LAST; i++ )
				if( ql.getVertex(i)->approx_quality >= AQ_VALID_MIN ) 
					pq += ql.getVertex(i)->vq;
			// average if more than one valid
			if( counter < 3 ) pq *= 1.0/(4-counter);
			// set to all invalid
			for(int i = QuadricLeaf::VX2D_FIRST; i <= QuadricLeaf::VX2D_LAST; i++ )
				if( ql.getVertex(i)->approx_quality < AQ_VALID_MIN ) 
					ql.getVertex(i)->vq = pq;
		}
	} );
}

bool SurfaceQuadricQTree::showNode( 
	MeshViewSet* set, int k, 
	DataVector< int > & shown_nodes,
	DataVector< DPoint3d > & pt3d_nodes ) const
{
	if( shown_nodes[k] != -1 ) return (shown_nodes[k] == 1); // 1 - valid, 0 - invalid, -1 - not known yet
	const QuadricNode& qnode = m_grid_vertices[k];
	if (withinParamRange(qnode.coord)) {
		shown_nodes[k] = 1;
		pt3d_nodes[k] = getPoint(qnode.coord);
		set->addLabel(pt3d_nodes[k], to_string(qnode.approx_quality));
		return true;
	}else{
		shown_nodes[k] = 0;
		return false;
	}
}

MeshViewSet * SurfaceQuadricQTree::createViewSetForBaseGrid( MeshViewSet* set, int id ) const
{
	if( set == nullptr ) set = new MeshViewSet;

	DPoint3d p1, p2;

	// horizontal
	for(int iy = 0; iy <= m_ny; iy++){
		int k1 = iy * (m_nx+1);
		p2 = m_base_surface->getPoint( m_grid_vertices[k1].coord );
		for(int ix = 0; ix < m_nx; ix++, k1++) {
			p1 = p2;
			p2 = m_base_surface->getPoint( m_grid_vertices[k1+1].coord );
			set->addEdge( p1, p2, id );
		}
	}
	// vertical
	for(int ix = 0; ix <= m_nx; ix++){
		int k1 = ix;
		p2 = m_base_surface->getPoint( m_grid_vertices[k1].coord );
		for(int iy = 0; iy < m_ny; iy++, k1 += (m_nx+1)) {
			p1 = p2;
			p2 = m_base_surface->getPoint( m_grid_vertices[k1+m_nx+1].coord );
			set->addEdge( p1, p2, id );
		}
	}

	return set;
}

/// create sketchy representation of this surface for the area with the given points
MeshViewSet * SurfaceQuadricQTree::createViewSet( MeshViewSet* set ) const
{
	if( set == nullptr ) set = new MeshViewSet;

	DataVector< int > shown_nodes( m_grid_vertices.countInt(), -1 );
	DataVector< DPoint3d > pt3d_nodes ( m_grid_vertices.countInt(), DPoint3d::zero );

	// horizontal
	for(int iy = 0; iy <= m_ny; iy++){
		for(int ix = 0; ix < m_nx; ix++) {
			int k1 = iy * (m_nx+1) + ix;
			int k2 = k1+1;

			bool ok1 = showNode( set, k1, shown_nodes, pt3d_nodes );
			bool ok2 = showNode( set, k2, shown_nodes, pt3d_nodes );

			const QuadricNode& qnode1 = m_grid_vertices[k1];
			const QuadricNode& qnode2 = m_grid_vertices[k2];

			if( ok1 || ok2 ) {
				DPoint2d mid_param(qnode1.coord, qnode2.coord, 0.5);
				if (withinParamRange(mid_param)) {
					DPoint3d mid = getPoint(mid_param);
					if(ok1)	set->addEdge(pt3d_nodes[k1], mid, 3);
					if(ok2) set->addEdge(pt3d_nodes[k2], mid, 3);
				}
			}
		}
	}
	// vertical
	for(int ix = 0; ix <= m_nx; ix++){
		for(int iy = 0; iy < m_ny; iy++) {
			int k1 = iy * (m_nx+1) + ix;
			int k2 = k1 + m_nx + 1;

			bool ok1 = showNode( set, k1, shown_nodes, pt3d_nodes );
			bool ok2 = showNode( set, k2, shown_nodes, pt3d_nodes );

			const QuadricNode& qnode1 = m_grid_vertices[k1];
			const QuadricNode& qnode2 = m_grid_vertices[k2];

			if( ok1 || ok2 ) {
				DPoint2d mid_param(qnode1.coord, qnode2.coord, 0.5);
				if (withinParamRange(mid_param)) {
					DPoint3d mid = getPoint(mid_param);
					if(ok1)	set->addEdge(pt3d_nodes[k1], mid, 3);
					if(ok2) set->addEdge(pt3d_nodes[k2], mid, 3);
				}
			}
		}
	}

	return set;
}

SurfacePtr SurfaceQuadricQTree::fitQuadricQTreeSurface(
		Metric3dContext &mc,
		SurfacePtr base_surface, 
		const DRect& box, 
		double tolerance,
		MeshContainer3dSurface* mesh, 
		const DataVector<MeshFace*> &sfaces,
		int nxy)
{
	assert( nxy >= 10 );

	double real_dx = base_surface->segmentLength(box.getX0Y0(), box.getX1Y0() );
	double real_dy = base_surface->segmentLength(box.getX0Y0(), box.getX0Y1() );
	double ratio = real_dx / real_dy;	// adjust m_nx:m_ny ratio to real lengths of domain
	int nx = std::max(4,(int)(sqrt(ratio*nxy)));
	int ny = std::max(4, (int)(nx/ratio));
	nx = std::min(nx, nxy/2);
	ny = std::min(ny, nxy/2);

	auto surface = std::make_shared<SurfaceQuadricQTree>( base_surface, box, nx, ny );

	// initialize with sfaces here...
	// 1. for each face -> list of nodes ?
	// 2. for each node -> find face
	surface->initializeWithFaces( mc, mesh, sfaces, tolerance );

	return surface;
}

void SurfaceQuadricQTree::countLocalCoordinates(const DPoint2d &pt, int &ix, int &iy) const
{
	ix = (int) (m_nx * (pt.x - m_box.x0) / (m_box.x1 - m_box.x0));
	iy = (int) (m_ny * (pt.y - m_box.y0) / (m_box.y1 - m_box.y0));
	if(ix > (m_nx - 1)) ix = m_nx - 1;
	if(iy > (m_ny - 1)) iy = m_ny - 1;
	if(ix < 0) ix = 0;
	if(iy < 0) iy = 0;
}

const QuadricLeaf* SurfaceQuadricQTree::findLastLeaf(const DPoint2d& pt) const
{
	int ix,iy;
	countLocalCoordinates(pt, ix, iy);
	return m_grid.get(iy,ix).findLastLeaf(pt);
}

const DPoint3d SurfaceQuadricQTree::getPoint(const DPoint2d& param) const
{
	//assert(m_initialized>0);
	const DPoint2d pt_fit = m_box.fitInPoint(param);
	return findLastLeaf(pt_fit)->getPoint(pt_fit, m_base_surface);
}

const DVector3d SurfaceQuadricQTree::getDerivative(int deriv, const DPoint2d& param) const
{
	//assert(m_initialized>0);
	const DPoint2d pt_fit = m_box.fitInPoint(param);
	return findLastLeaf(pt_fit)->getDerivative(deriv, pt_fit, m_base_surface);
}

bool SurfaceQuadricQTree::withinParamRange( const DPoint2d& param ) const
{
	if( m_box.contains( param ) ) 
		return findLastLeaf(param)->withinParamRange(param);
	else
		return false;
}

void QuadricLeaf::setNeighbour(int side, QuadricLeaf* nb, bool skip_first)
{
	if(!skip_first) m_neighbours[side] = nb;
	if(isSplit()){
		m_leaves[mid_to_nodes[side][0]]->setNeighbour(side, nb);
		m_leaves[mid_to_nodes[side][1]]->setNeighbour(side, nb); 
	}
}

#define MAX_SPLIT_LEVEL 10

bool QuadricLeaf::split(QuadricVertices& grid_vertices, 
		SurfaceConstPtr surface, QuadricVerticesPtr* split_vertices)
{
	if(m_level >= MAX_SPLIT_LEVEL)
		return false;

	assert(m_leaves[LF_FIRST] == nullptr);

	// Balance before split...
	//balance(grid_vertices, surface, split_vertices);
	
	m_leaves[LF_SW] = new QuadricLeaf(m_middle.x-m_dl.x, m_middle.y-m_dl.y, m_dl.x/2, m_dl.y/2, m_level+1);
	m_leaves[LF_SE] = new QuadricLeaf(m_middle.x+m_dl.x, m_middle.y-m_dl.y, m_dl.x/2, m_dl.y/2, m_level+1);
	m_leaves[LF_NW] = new QuadricLeaf(m_middle.x-m_dl.x, m_middle.y+m_dl.y, m_dl.x/2, m_dl.y/2, m_level+1);
	m_leaves[LF_NE] = new QuadricLeaf(m_middle.x+m_dl.x, m_middle.y+m_dl.y, m_dl.x/2, m_dl.y/2, m_level+1);

	QuadricNode* qv_middle;

	// middle
	size_t index = grid_vertices.add(QuadricNode(m_middle));
	qv_middle = &(grid_vertices[index]);
	if(split_vertices) split_vertices->add(qv_middle);

	// check middle nodes at edges
	QuadMidVertexWhich opposite[4] = {MVX2D_UP, MVX2D_LEFT, MVX2D_DOWN, MVX2D_RIGHT};
	for(int i = MVX2D_FIRST; i <= MVX2D_LAST; i++){
		if(!m_mid_vertices[i]){	// vertex missing
			QuadricNode *q0 = m_vertices[mid_to_nodes[i][0]];
			QuadricNode *q1 = m_vertices[mid_to_nodes[i][1]];
			index = grid_vertices.add(QuadricNode(DPoint2d::average(q0->coord, q1->coord)));
			m_mid_vertices[i] = &(grid_vertices[index]);
			if(split_vertices) split_vertices->add(m_mid_vertices[i]);
			if(m_neighbours[i] && (m_neighbours[i]->m_mid_vertices[opposite[i]] == nullptr))
				m_neighbours[i]->m_mid_vertices[opposite[i]] = m_mid_vertices[i];
		}
	}

	// neighbours
	// - between leaves
	m_leaves[LF_SW]->m_neighbours[NB_UP] = m_leaves[LF_NW];
	m_leaves[LF_NW]->m_neighbours[NB_DOWN] = m_leaves[LF_SW];
	m_leaves[LF_SE]->m_neighbours[NB_UP] = m_leaves[LF_NE];
	m_leaves[LF_NE]->m_neighbours[NB_DOWN] = m_leaves[LF_SE];
	m_leaves[LF_SW]->m_neighbours[NB_RIGHT] = m_leaves[LF_SE];
	m_leaves[LF_SE]->m_neighbours[NB_LEFT] = m_leaves[LF_SW];
	m_leaves[LF_NW]->m_neighbours[NB_RIGHT] = m_leaves[LF_NE];
	m_leaves[LF_NE]->m_neighbours[NB_LEFT] = m_leaves[LF_NW];
	// - outside
	for(int i = NB_FIRST; i <= NB_LAST; i++){
		if(!m_neighbours[i]) continue;
		int l_i0 = mid_to_nodes[i][0];
		int l_i1 = mid_to_nodes[i][1];
		if(m_neighbours[i]->isSplit()){
			int l_j0 = mid_to_nodes[opposite[i]][0];
			int l_j1 = mid_to_nodes[opposite[i]][1];
			m_neighbours[i]->m_leaves[l_j0]->setNeighbour(opposite[i], m_leaves[l_i0]);
			m_neighbours[i]->m_leaves[l_j1]->setNeighbour(opposite[i], m_leaves[l_i1]);
			m_leaves[l_i0]->m_neighbours[i] = m_neighbours[i]->m_leaves[l_j0];
			m_leaves[l_i1]->m_neighbours[i] = m_neighbours[i]->m_leaves[l_j1];
		}else{
			m_leaves[l_i0]->m_neighbours[i] = m_neighbours[i];
			m_leaves[l_i1]->m_neighbours[i] = m_neighbours[i];
		}
	}
	// vertices
	m_leaves[LF_SW]->setVertices(m_vertices[VX2D_SW], m_mid_vertices[MVX2D_DOWN], 
		m_mid_vertices[MVX2D_LEFT], qv_middle);
	m_leaves[LF_SE]->setVertices(m_mid_vertices[MVX2D_DOWN], m_vertices[VX2D_SE], 
		qv_middle, m_mid_vertices[MVX2D_RIGHT]);
	m_leaves[LF_NW]->setVertices(m_mid_vertices[MVX2D_LEFT], qv_middle,
		m_vertices[VX2D_NW], m_mid_vertices[MVX2D_UP]);
	m_leaves[LF_NE]->setVertices(qv_middle, m_mid_vertices[MVX2D_RIGHT],
		m_mid_vertices[MVX2D_UP], m_vertices[VX2D_NE]);
	// mid_vertices
	for(int i = NB_FIRST; i <= NB_LAST; i++){
		if(!m_neighbours[i] || !m_neighbours[i]->isSplit()) continue;
		for(int j = 0; j < 2; j++){
			int l_j = mid_to_nodes[opposite[i]][j];
			if(m_neighbours[i]->m_leaves[l_j]->isSplit()){
				int l_i = mid_to_nodes[i][j];
				m_leaves[l_i]->m_mid_vertices[i] = 
					m_neighbours[i]->m_leaves[l_j]->m_mid_vertices[opposite[i]];
			}
		}
	}
	// done
	return true;
}

#define BALANCE_LEVEL 1

void QuadricLeaf::balance(QuadricVertices& grid_vertices,
						  SurfaceConstPtr surface, 
						  QuadricVerticesPtr* split_vertices)
{
	for(int i = NB_FIRST; i <= NB_LAST; i++){
		if(!m_neighbours[i]) continue;	// border of qtree
		while((m_level - m_neighbours[i]->m_level) >= BALANCE_LEVEL)
			if(!m_neighbours[i]->split(grid_vertices, surface, split_vertices)) break;
	}
}

#define NEAR_EPS2 4.0

bool QuadricLeaf::adaptForFaces( QuadricVertices & qv, SurfaceConstPtr surface,
						Metric2dContext& mc2d, const DataVector<MeshFace*> & sfaces )
{
	int containing_ct = 0;

	for(int i = 0; i < 4; i++)
		if( m_vertices[i]->valid(AQ_VALID_MAX) )
			containing_ct++;

	if (containing_ct == 4) { // fully valid
		m_valid = true;
		return true; 
	}

	if (containing_ct == 0) {

		double min_min_dist2 = std::min(
			std::min(m_vertices[0]->min_dist2, m_vertices[1]->min_dist2),
			std::min(m_vertices[2]->min_dist2, m_vertices[3]->min_dist2));

		mc2d.countMetricAtPoint(m_middle);
		double min_cell_dist2 = 0.25 * std::min(
			mc2d.transformPStoMS( DVector2d(getWidth(),  0.0) ).length2(),
			mc2d.transformPStoMS( DVector2d(0.0, getHeight()) ).length2());

		if (min_min_dist2 > min_cell_dist2) { // fully invalid, too far from any 
			m_valid = false;
			return true;
		}
	}

	double max_min_dist2 = std::max(
		std::max(m_vertices[0]->min_dist2, m_vertices[1]->min_dist2),
		std::max(m_vertices[2]->min_dist2, m_vertices[3]->min_dist2));

	if( max_min_dist2 < NEAR_EPS2 ) {
		for(int i = 0; i < 4; i++)
			if( !m_vertices[i]->valid() )
				m_vertices[i]->setQuality( AQ_VALID_MIN );			

		m_valid = true;
		return true; // partially valid, sufficiently close to face-domain
	}
	
	// split
	QuadricVerticesPtr split_vertices; // split vertics to initialize 
	if (!split(qv, surface, &split_vertices)) return false;

	// initialize new vertices
	// ... from all faces (or only from the faces of valid vertices of this leaf) ...
	split_vertices.forEach( [ &mc2d, sfaces, surface ] (QuadricNode* qn) {
		mc2d.countMetricAtPoint( qn->coord );
		DMPoint2d qpt = mc2d.transformPStoMS( qn->coord );

		size_t sfct = sfaces.countInt();
		DPoint2d a,b,c;
		for(size_t j = 0; j < sfct; j++){
			bool valid = 
				sfaces[j]->getPoint(0)->checkAndGetSurfaceParam( surface, a ) &&
				sfaces[j]->getPoint(1)->checkAndGetSurfaceParam( surface, b ) &&
				sfaces[j]->getPoint(2)->checkAndGetSurfaceParam( surface, c );
			if( valid && qn->checkFaceDist( qpt, mc2d.transformPStoMS( a ),
					mc2d.transformPStoMS( b ), mc2d.transformPStoMS( c ), sfaces[j]) ) 
				break;
		}
	} );

	// check recursively
	int ok = 0;
	for (int i = 0; i < 4; i++)
		if (m_leaves[i]->adaptForFaces(qv, surface, mc2d, sfaces)) ok++;

	// m_valid is irrelevant, since not a leaf anymore...
	return (ok == 4);
}

void QuadricLeaf::addToViewSet(MeshViewSet* set, SurfaceConstPtr surface,
	TagExtended::TagType tag) const
{
	if (isSplit()) {
		for (int i = LF_FIRST; i <= LF_LAST; i++)
			m_leaves[i]->addToViewSet(set, surface, tag);
	} else if (m_valid) {
		auto data = std::make_shared<MeshViewFaceData>(4);
		data->area_id = m_level;
		data->quality = 1.0;

		DPoint3d dpts[4];
		int order[4] = { 0, 1, 3, 2 };
		for (int i = 0; i < 4; i++) {
			if (surface) {
				dpts[i] = surface->getPoint(m_vertices[order[i]]->coord);
			}else {
				dpts[i] = DPoint3d(m_vertices[order[i]]->coord, 0.0);
			}
		}

		DPoint3d middle(dpts[0], dpts[2], 0.5);

		for (int i = 0; i < 4; i++)
			data->pts.add( middle + (dpts[i] - middle) * 0.95 );

		data->normal = (data->pts[1] - data->pts[0]).crossProduct(
						data->pts[2] - data->pts[0]).normalized();

		data->indices.add(0);	data->indices.add(1);	data->indices.add(2);
		data->indices.add(0);	data->indices.add(2);	data->indices.add(3);

		set->m_faces.add(data);
	}

}
