/////////////////////////////////////////////////////////////////////////////
// MeshContainer3dSurface.cpp
// Class responsible for storing points and faces of a 3D surface mesh
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2010-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "MeshContainer3dSurface.h"
#include "MeshViewSet.h"
#include "MeshEdge3d.h"
#include "MeshFace.h"
#include "MeshDomainVolume.h"
#include "Metric3dContext.h"
#include "DLeastSquaresFitting.h"
#include "MeshGenerator3dSurface.h"
#include "DataList.h"
#include "DataStatistics.h"
#include "SurfaceParametric.h"
#include "SurfacePlane.h"
#include "SurfacePlanarQuadric.h"
#include "SurfaceCylinder.h"
#include "SurfaceBSplinePlanar.h"
#include "SurfaceBSplineCylindrical.h"
#include "Curve3dParametric.h"
#include "Curve3dBSpline.h"
#include "Curve2dSegment.h"
#include "Curve3dSurfaceParametric.h"
#include "DPlane.h"
#include "DLine.h"
#include "DTriangle.h"
#include "DPlanarQuadric.h"
#include "SurfaceDomainHull.h"
#include "SurfaceDomainFaces.h"
#include "SurfaceDomainFacesPlanar.h"
#include "ControlSpace3dAdaptive.h"
#include "ControlSpace2dIdentity.h"
#include "DOrientedBox.h"
#include "DLinearQuadric.h"
#include "Curve3dLinearQuadric.h"
#include "SurfaceQuadricQTree.h"

#define SHOW_GATHER_LAYERED
//#define SHOW_COMMON_SURFACE

int MeshContainer3dSurface::param_surface_domain_type = MeshContainer3dSurface::SDOMAIN_FACES_PLANAR;

/////////////////////////////////////////////////////////////////////////////
// Standard constructor, part_size defines step of increasing the capacity of arrays
MeshContainer3dSurface::MeshContainer3dSurface(int part_size) 
	: m_part_size(part_size), m_border_stage(BORDER_KNOWN), m_box_diameter(-1.0)
{
	m_points = new DataContainer<MeshPoint3d>(part_size);
	m_faces = new DataContainer<MeshFace>(part_size);
}

/////////////////////////////////////////////////////////////////////////////
// Destructor
MeshContainer3dSurface::~MeshContainer3dSurface()
{
	deleteAll();
	delete m_faces; m_faces = nullptr;
	delete m_points; m_points = nullptr;
}

void MeshContainer3dSurface::removeAllTags(TagExtended::TagType tag_type)
{
	for(int i = 0; i < m_points->countInt(); i++) // clear all tags for points, edges and elements
		getPointAt(i)->removeTag(tag_type);
	for(int i = 0; i < m_faces->countInt(); i++) 
		getFaceAt(i)->removeTag(tag_type);
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge())
		it.getEdge()->removeTag(tag_type);
}

/////////////////////////////////////////////////////////////////////////////
// Removes all data
void MeshContainer3dSurface::deleteAll()
{
	m_faces->deleteAll();
	m_points->deleteAll();
}

//////////////////////////////////////////////////////////////////////
// Bounding box containing all points in this container
const DBox& MeshContainer3dSurface::getBoundingBox(bool force_update)
{
	if( force_update || !m_box.valid ) {
		m_box.clear();

		int count = m_points->countInt();
		for(int i = 0; i < count; i++)
			m_box.addPoint(m_points->getDataAt(i)->getCoordinates());

//	for(IteratorEdge2d it = getFirstEdge2d(); it.isValid(); it.nextEdge())
//		it.getEdge()->addToBoundingRect(rect);
	}

	return m_box;
}

double  MeshContainer3dSurface::getBoundingBoxDiameter(bool force_update)
{
	if( force_update || !m_box.valid) {
		m_box_diameter = getBoundingBox( force_update ).getDiameter();
	}else if( m_box_diameter <= 0.0 ) {
		m_box_diameter = m_box.getDiameter();
	}

	return m_box_diameter;
}

void MeshContainer3dSurface::setControlSpace(CS3dPtr space)
{
	m_control = space;
}

int MeshContainer3dSurface::getFacesCount(int edge_count) const
{
	int ect = getFacesCount();
	int count = 0;
	for(int i = 0; i < ect; i++) 
		if(getFaceAt(i)->getEdgeCount() == edge_count) ++count;
	return count;
}

MeshViewSet* MeshContainer3dSurface::getViewSet(MeshViewSet* view_set, TagExtended::TagType color_tag_type) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;

	int fct = getFacesCount();
	int pct = getPointsCount();
	if(fct + pct == 0) return view_set;

	if(view_set){
		view_set->prepareFreePlace(pct, 2*fct, fct, 0);
	}else{
		view_set = new MeshViewSet(pct, 2*fct, fct, 0);
	}
	//view_set->setMesh(this);

	// faces
	for(int i = 0; i < fct; i++){
		MeshFace* face = getFaceAt(i);
		view_set->addFace(face, 
			(color_tag_type != TagExtended::TAG_NONE) ? face->getIntTag(color_tag_type, -2) : -2, 
			MeshViewSet::param_shrink, face->getBlock(0) != nullptr);
	}

	// others
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if (color_tag_type != TagExtended::TAG_NONE)
			view_set->addPoint( point->getCoordinates(), point->getIntTag(color_tag_type, 0), point->getIntTag(color_tag_type, 0) );
		else
			view_set->addPoint(point);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0){
				view_set->addEdge(edge);
			}
		}
	}

	if(true){ // additional quality/stat info
		view_set->addInfo("mesh type", "3d surface");
		view_set->addInfo("valid", isValid(true) ? "true" : "false");

		int f3 = 0, f4 = 0, fn = 0, fmax = 0;
		for(int i = 0; i < fct; i++){
			MeshFace* face = getFaceAt(i);
			int fpct = face->getPointCount();
			if(fpct == 3) f3++;
			else if(fpct == 4) f4 ++;
			else fn++;
			if(fpct > fmax) fmax = fpct;
		}
		if(f3 < fct) view_set->addInfo("all face count", fct);
		if(f3 > 0) view_set->addInfo("tri-face count", f3);
		if(f4 > 0) view_set->addInfo("quad-face count", f4);
		if(fn > 0) view_set->addInfo("poly-face count", fn);
		if(fmax > 4) view_set->addInfo("max poly-face", fmax);

		double min_elen, max_elen;
		int ect = 0;
		int e1 = 0;
		for(auto it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
			MeshEdge3d* edge = it.getEdge();
			if(edge->getFaceCount() == 1) e1++;
			double len = edge->getLength();
			if(ect == 0) min_elen = max_elen = len;
			else{
				if(len > max_elen) max_elen = len;
				if(len < min_elen) min_elen = len;
			}
			ect++;
		}
		view_set->addInfo("closed surface", (e1 == 0) ? "yes" : "no");
		view_set->addInfo("edge count", ect);
		view_set->addInfo("min edge len", min_elen);
		view_set->addInfo("max edge len", max_elen);
		view_set->addInfo("node count", pct);
	}

	return view_set;
}

MeshViewSet* MeshContainer3dSurface::getDebugViewSet(const MeshPoint3d* pt1, const MeshPoint3d* pt2, double radius, TagExtended::TagType color_tag_type) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;
	if(!pt1) return nullptr;

	const DPoint3d middle = pt1->getCoordinates();
	double r2 = 0.0;
	if(pt2) r2 = middle.distance2(pt2->getCoordinates());
	else{
		int rank = pt1->getRank();
		for(int i = 0; i < rank; i++)
			r2 = std::max(r2, sqr(pt1->getEdge(i)->getLength()));
	}
	r2 *= sqr(radius);

	int fct = getFacesCount();
	int pct = getPointsCount();

//	ofstream of("debug_view.log");

	MeshViewSet *view_set = new MeshViewSet(pct, 2*fct, fct, 0);
//	view_set->setMesh(this);

	// faces
	for(int i = 0; i < fct; i++){
		MeshFace* face = getFaceAt(i);
		//if(el->getType() == ELEMENT_MESH_AREA) continue;
		if(middle.distance2(face->getMiddlePoint()) > r2) continue;
		view_set->addFaceWithEdges(face, 
			(color_tag_type != TagExtended::TAG_NONE) ? face->getIntTag(color_tag_type, -2) : -2, 
			MeshViewSet::param_shrink, face->getBlock(0) != nullptr);
	}

	// others
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(middle.distance2(point->getCoordinates()) > r2) continue;
		if(point == pt1)
			view_set->addPoint(point, 1);
		else if(point == pt2)
			view_set->addPoint(point, 2);
		else
			view_set->addPoint(point, 0);
	}

	return view_set;
}

MeshViewSet* MeshContainer3dSurface::getDebugViewSetTopological(const MeshPoint3d* start_point, int radius, TagExtended::TagType color_tag_type) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;
	if(!start_point) return nullptr;

	int pct = getPointsCount();
	int fct = getFacesCount();

	DataVector<int> p_layers(pct, -1);
	DataVector<int> f_layers(fct, -1);
	DataVector<int> f_list(fct);
	DataVector<int> p_list(pct);
	DataVector<int> f_layers_count;

	gatherFirstLayerTopologicalForPoint(start_point, f_layers_count, p_list, p_layers, f_list, f_layers);
	if(radius > 1)
		gatherLayeredVerticesTopological(f_layers_count, p_list, p_layers, f_list, f_layers, radius);
	int l_outside = 1+f_layers[ f_list[ f_layers_count[ f_layers_count.countInt()-2] ] ];

//	ofstream of("debug_view.log");

	int vpct = p_list.countInt();
	int vfct = f_list.countInt();
	MeshViewSet *view_set = new MeshViewSet(vpct, 2*vfct, vfct, 0);
//	view_set->setMesh(this);

	// faces
	for(int i = 0; i < vfct; i++){
		int fid = f_list[i];
		if(f_layers[fid] >= l_outside) continue;
		MeshFace* face = getFaceAt(fid);
		view_set->addFaceWithEdges(face, 
			(color_tag_type != TagExtended::TAG_NONE) ? face->getIntTag(color_tag_type, -2) : -2, 
			MeshViewSet::param_shrink, face->getBlock(0) != nullptr);
	}

	// others
	for(int i = 0; i < vpct; i++){
		int pid = p_list[i];
		if(p_layers[pid] >= l_outside) continue;
		MeshPoint3d* point = getPointAt(pid);
		if(point == start_point)
			view_set->addPoint(point, 1);
		else
			view_set->addPoint(point, 0);
	}

	return view_set;
}

MeshViewSet* MeshContainer3dSurface::getDebugViewSetTopological(MeshFace* start_face, int radius, TagExtended::TagType color_tag_type) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;
	if(!start_face) return nullptr;

	int pct = getPointsCount();
	int fct = getFacesCount();

	DataVector<int> p_layers(pct, -1);
	DataVector<int> f_layers(fct, -1);
	DataVector<int> f_list(fct);
	DataVector<int> p_list(pct);
	DataVector<int> f_layers_count;

	DataVector<MeshFace*> start_faces;
	start_faces.add(start_face);

	gatherFirstLayerTopological(start_faces, f_layers_count, p_list, p_layers, f_list, f_layers);
	if(radius > 1)
		gatherLayeredVerticesTopological(f_layers_count, p_list, p_layers, f_list, f_layers, radius);
	int l_outside = 1;
	if(f_layers_count.countInt() >= 2)
		l_outside = 1+f_layers[ f_list[ f_layers_count[ f_layers_count.countInt()-2] ] ];

//	ofstream of("debug_view.log");

	int vpct = p_list.countInt();
	int vfct = f_list.countInt();
	MeshViewSet *view_set = new MeshViewSet(vpct, 2*vfct, vfct, 0);
//	view_set->setMesh(this);

	// faces
	for(int i = 0; i < vfct; i++){
		int fid = f_list[i];
		if(f_layers[fid] >= l_outside) continue;
		MeshFace* face = getFaceAt(fid);
		view_set->addFaceWithEdgesAndPoints(face, 
			(color_tag_type != TagExtended::TAG_NONE) ? face->getIntTag(color_tag_type, -2) : ((face==start_face)?1:-2), 
			MeshViewSet::param_shrink, face->getBlock(0) != nullptr);
	}

	return view_set;
}

/// Clears all boundary flags from surface entities
void MeshContainer3dSurface::clearBoundaryFlags()
{
	int fct = getFacesCount();
	for(int i = 0; i < fct; i++) getFaceAt(i)->clearBorder();
	int pct = getPointsCount();
	for(int i = 0; i < pct; i++) getPointAt(i)->clearBorder();
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge())
		it.getEdge()->clearBorder();
	setBorderStage(BORDER_UNKNOWN);
}

/// Sets boundary flags (edges and points) according to the given tag_type
void MeshContainer3dSurface::setBoundaryTagEdges(TagExtended::TagType tag_type)
{
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(edge->isBorder()) continue; // already set
		int rank = edge->getFaceCount();
		if(rank < 2) continue;

		MeshFace* face = edge->getFaceAt(0);
		if(face->availableTag(tag_type)){
			int tag_value = face->getIntTag(tag_type);
			for(int i = 1; i < rank; i++){
				if(edge->getFaceAt(i)->getIntTag(tag_type, tag_value-1) != tag_value){
					edge->setBorder( TagBorder::OUTER | TagBorder::RIDGE ); // ok, found at least one different, mark as boundary
					edge->getMeshPoint(0)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
					edge->getMeshPoint(1)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
					break;
				}
			}
		}else{
			for(int i = 1; i < rank; i++){
				if(edge->getFaceAt(i)->availableTag(tag_type)){
					edge->setBorder( TagBorder::OUTER | TagBorder::RIDGE ); // ok, found at least one different, mark as boundary
					edge->getMeshPoint(0)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
					edge->getMeshPoint(1)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
					break;
				}
			}
		}
	}
	// additionally search for single non-boundary points of contact
	int pct = getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(!point->availableTag(tag_type)) continue;
		int tag = point->getIntTag(tag_type);
		int counter = 0;
		while(tag && counter < 2){
			if((tag & 1) != 0) counter++;
			tag >>= 1;
		}
		if(counter > 1) point->setBorder(TagBorder::OUTER | TagBorder::CORNER);  // if more than one tag-flag
	}

/*	const char IS_TAG = 1;
	const char VISITED = 2;
	DataVector<char> points_istag(pct, 2); // 0 - no tag, 1 - is tag, 2 - not visited yet
	DataVector<int> points_tag(pct, 0); // tag value, if any
	int fct = getFacesCount();
	for(int i = 0; i < fct; i++){
		MeshFace* face = getFaceAt(i);
		bool face_istag = face->availableTag(tag_type);
		int face_vtag = face_istag ? face->getIntTag(tag_type) : 0;
		// ... check points of this face
		int fpct = face->getPointCount();
		for(int j = 0; j < fpct; j++){
			MeshPoint3d* point = face->getPoint(j);
			if(point->isBorder()) continue;
			int id = point->getIndex();
			if(points_istag[id] & VISITED != 0){
				// already initialized, so check if ok
				bool point_istag = points_istag[id] & IS_TAG != 0;
				if( point_istag != face_istag){
					point->setBorder();
				}else if( point_istag && (points_tag[id] != face_vtag)){
					point->setBorder();
				}
			}else{ 
				// first time
				points_istag[id] |= VISITED;
				if(face_istag){
					points_istag[id] |= IS_TAG;
					points_tag[id] = face_vtag;
				}
			}
		}
	}
*/
	setBorderStage(BORDER_TOPOLOGY);
}

/// Sets boundary flags (edges and points) for "feature" edges
void MeshContainer3dSurface::setBoundaryFeatureEdges()
{
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(edge->isBorder()) continue; // already set
		if(edge->getFaceCount() != 2){
			edge->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
			edge->getMeshPoint(0)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
			edge->getMeshPoint(1)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
		}
	}
	int pct = getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(point->isBorder( TagBorder::CORNER) ) continue;
		int bect = point->getBorderEdgesCount();
		if(bect == 0) continue;
		if(!point->isBorder()) point->setBorderFlags( TagBorder::OUTER | TagBorder::RIDGE );
		if(bect != 2) point->setBorderFlags( TagBorder::CORNER );
	}
	setBorderStage(BORDER_TOPOLOGY);
}

/// Sets boundary flags (edges and points) for sharp edges (fmax - max value for scalar product of normal vectors)
void MeshContainer3dSurface::setBoundarySharpEdges(double fmax, TagExtended::TagType tag_type, int tag_value)
{
	int pct = getPointsCount();
	DataVector<MeshEdge3d*> sharp_edges(100);
	DataVector<int> sharp_vertices(pct, 0);

	DataHashTableKeyValue<MeshEdge3d*, double> hash_bedges(6*pct, nullptr);

	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(edge->isBorder()){
			hash_bedges.insert(edge, -3.0);
			continue; // already set
		}
		if(edge->getFaceCount() != 2){
			edge->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
			edge->getMeshPoint(0)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
			edge->getMeshPoint(1)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
			hash_bedges.insert(edge, -2.0);
			continue; // not applicable
		}

		MeshFace* f0 = edge->getFaceAt(0);
		MeshFace* f1 = edge->getFaceAt(1);
		if(tag_type != TagExtended::TAG_NONE){ // if tag_type should be used, both faces have to be tagged properly
			if(f0->getIntTag(tag_type, tag_value-1) != tag_value) continue;
			if(f1->getIntTag(tag_type, tag_value-1) != tag_value) continue;
		}

		DVector3d dn0, dn1;
		if( ! f0->checkAndGetNormalVector(dn0) ) continue;
		if( ! f1->checkAndGetNormalVector(dn1) ) continue;

			double sp = dn0.scalarProduct(dn1);
		edge->setDoubleTag(TagExtended::TAG_NORMAL_SP, sp);
//		if(!f0->sameOrientation(f1)) sp = -sp;
		if(sp < fmax){
			sharp_edges.add(edge);
			sharp_vertices[ edge->getMeshPoint(0)->getIndex() ]++;
			sharp_vertices[ edge->getMeshPoint(1)->getIndex() ]++;
			hash_bedges.insert(edge, sp);
//			edge->setBorder( TagBorder::OUTER | TagBorder::RIDGE ); // ok, sharp enough, mark as boundary
//			edge->getMeshPoint(0)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
//			edge->getMeshPoint(1)->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
		}
	}

#ifdef _DEBUG
	if(sharp_edges.notEmpty()){
		MeshViewSet* set = new MeshViewSet(2*sharp_edges.countInt(), 10*sharp_edges.countInt(), 2*sharp_edges.countInt(), 0);
		for(int i = 0; i < sharp_edges.countInt(); i++){
			MeshEdge3d* edge = sharp_edges[i];
			set->addEdge(edge);
			set->addPoint(edge->getMeshPoint(0));
			set->addPoint(edge->getMeshPoint(1));
			//set->addFaceWithEdges( edge->getFaceAt(0) );
			//set->addFaceWithEdges( edge->getFaceAt(1) );
		}
		SHOW_MESH("Sharp edges (to check)", set);
	}
#endif

	for(int i = 0; i < sharp_edges.countInt(); i++){
		MeshEdge3d* edge = sharp_edges[i];
		for(int j = 0; j < 2; j++){
			MeshFace* face = edge->getFaceAt(j);
			bool all_sharp = true;
			int ect = face->getEdgeCount();
			for(int k = 0; all_sharp && (k < ect); k++)
				all_sharp = hash_bedges.contains( face->getEdge(k) );
			if(all_sharp){
				SHOW_MESH("face with all sharp edges", getDebugViewSetTopological(face, 2));
			}
		}
	}

	// extend open border edges - if fairly linear and with softened "sharpness" of normals for faces
	const double MIN_V_SP = 0.866;
//	const double MAX_N_SP = 0.866; // -> cos(30o)  //(1.0 + fmax) * 0.5;
	for(int i = 0; i < sharp_edges.countInt(); i++){
		MeshEdge3d* edge = sharp_edges[i];
		MeshPoint3d* ep[2] = { edge->getMeshPoint(0), edge->getMeshPoint(1) };
		int ep0ct = sharp_vertices[ ep[0]->getIndex() ];
		int ep1ct = sharp_vertices[ ep[1]->getIndex() ];
		int idx = -1;
		if( ep0ct == 1 && ep1ct > 1) idx = 0;
		else if( ep0ct > 1 && ep1ct == 1) idx = 1;

		if( idx >= 0 ){ // open boundary edge
			MeshPoint3d* epx = ep[idx];
			// - find best edge
			const DVector3d v0 = ( ep[idx]->getCoordinates() - ep[1-idx]->getCoordinates() ).normalized() ;
			int rank = epx->getRank();
			MeshEdge3d* best_edge = nullptr;
			double best_v_sp = -2.0;
			for(int j = 0; j < rank; j++){
				MeshEdge3d* other_edge = epx->getEdge(j);
				if( other_edge == edge ) continue;
				MeshPoint3d* other_point = other_edge->getOtherPoint(epx);
				double v_sp = v0.scalarProduct( ( other_point->getCoordinates() - ep[idx]->getCoordinates() ).normalized() );
				if( !best_edge || v_sp > best_v_sp ){
					best_edge = other_edge;
					best_v_sp = v_sp;
				}
			}

			const double MAX_N_SP = best_v_sp + 0.5;

			//if( best_edge ){
			//	MeshViewSet* set = new MeshViewSet;
			//	set->addEdge(edge);
			//	set->addEdge(best_edge, 1);
			//	for(int k = 0; k < edge->getFaceCount(); k++)
			//		set->addFace( edge->getFaceAt(k) );
			//	set->addInfo( "v_sp / min", (best_v_sp > MIN_V_SP ? "OK " : "NO ") + to_string( best_v_sp ) 
			//			+ " | " + to_string( MIN_V_SP ) );
			//	double n_sp = best_edge->getDoubleTag(TagExtended::TAG_NORMAL_SP, -2.0);
			//	set->addInfo( "n_sp / max", (n_sp < MAX_N_SP ? "OK " : "NO ") + to_string( n_sp ) 
			//			+ " | " + to_string( MAX_N_SP ) );
			//	SHOW_MESH("sharp edge straight-extension", set);
			//}

			if( best_edge && best_v_sp > MIN_V_SP ) {
				assert( !best_edge->isBorder() );
				assert( best_edge->getFaceCount() == 2 );
				
				assert( best_edge->availableTag( TagExtended::TAG_NORMAL_SP ) );
				double n_sp = best_edge->getDoubleTag(TagExtended::TAG_NORMAL_SP, -1.0);

				if( n_sp < MAX_N_SP ){
					sharp_edges.add(best_edge);
					sharp_vertices[ best_edge->getMeshPoint(0)->getIndex() ]++;
					sharp_vertices[ best_edge->getMeshPoint(1)->getIndex() ]++;
				}
			}
		}
	}

	// join interrupted chains of sharp edges
	for(int i = 0; i < sharp_edges.countInt(); i++){
		MeshEdge3d* edge = sharp_edges[i];
		MeshPoint3d* ep[2] = { edge->getMeshPoint(0), edge->getMeshPoint(1) };
		int ep0ct = sharp_vertices[ ep[0]->getIndex() ];
		int ep1ct = sharp_vertices[ ep[1]->getIndex() ];
		int idx = -1;
		if( ep0ct == 1 && ep1ct > 1) idx = 0;
		else if( ep0ct > 1 && ep1ct == 1) idx = 1;

		if(idx >= 0){
			MeshPoint3d* epx = ep[idx];
			int rank = epx->getRank();
			for(int j = 0; j < rank; j++){
				MeshEdge3d* other_edge = epx->getEdge(j);
				MeshPoint3d* other_point = other_edge->getOtherPoint(epx);
				if( other_point == ep[1-idx] ) continue;
				if( sharp_vertices[ other_point->getIndex() ] == 1 ){
					// check direction ?

					//MeshViewSet* set = new MeshViewSet;
					//set->addEdge(edge);
					//set->addEdge(other_edge, 1);
					//for(int k = 0; k < edge->getFaceCount(); k++){
					//	set->addFace( edge->getFaceAt(k) );
					//}
					//SHOW_MESH("sharp edge extension", set);
		
					sharp_edges.add(other_edge);
					sharp_vertices[ epx->getIndex() ]++;
					sharp_vertices[ other_point->getIndex() ]++;
					break;
				}
			}
		}
	}

	// for now - simple, remove all single sharp edges, if there is no chain...
	for(int i = 0; i < sharp_edges.countInt(); ){
		MeshEdge3d* edge = sharp_edges[i];
		MeshPoint3d* ep0 = edge->getMeshPoint(0);
		MeshPoint3d* ep1 = edge->getMeshPoint(1);
		int ep0ct = sharp_vertices[ ep0->getIndex() ];
		int ep1ct = sharp_vertices[ ep1->getIndex() ];

		if( (ep0ct == 1) &&	(ep1ct == 1)) sharp_edges.removeAt(i);
		else{
			edge->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
			ep0->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
			ep1->setBorder( TagBorder::OUTER | TagBorder::RIDGE );
			if(ep0ct != 2) ep0->setBorderFlags( TagBorder::CORNER );
			if(ep1ct != 2) ep1->setBorderFlags( TagBorder::CORNER );
			i++;
		}
	}

	for(int i = 0; i < pct; i++){
		if( sharp_vertices[i] != 2)	continue;
		MeshPoint3d* point = getPointAt(i);
		MeshEdge3d* edge0 = nullptr;
		MeshEdge3d* edge1 = nullptr;
		int pect = point->getRank();
		for(int j = 0; j < pect; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if (edge->isBorder()) { if (edge0) edge1 = edge; else edge0 = edge; }
		}
		assert(edge0 && edge1);
		DVector3d v0 = (point->getCoordinates() - edge0->getOtherPoint(point)->getCoordinates()).normalized();
		DVector3d v1 = (edge1->getOtherPoint(point)->getCoordinates() - point->getCoordinates()).normalized();
		if(v0.scalarProduct(v1) < fmax){
			point->setBorderFlags( TagBorder::CORNER ); // ok, sharp enough, mark as corner-boundary
		}
	}

	int corner_points = 0;
	for(int i = 0; i < pct; i++)
		if(getPointAt(i)->isBorder(TagBorder::CORNER)) corner_points++;
	if(corner_points > 0){
		LOG4CPLUS_INFO(MeshLog::logger_console, "Corner points: " << corner_points);
	}

	setBorderStage(BORDER_KNOWN);
}

/// Removes all local surface definitions and counters
void MeshContainer3dSurface::clearLocalShapes()
{
	if(m_local_surfaces.countInt() > 0 || m_local_curves.countInt() > 0){
		m_local_surfaces.clear();
		m_local_curves.clear();
		// clear info in points' tags
		int pct = getPointsCount();
		for(int i = 0; i < pct; i++){
			getPointAt(i)->clearLocalShapes();
		}
		// and in faces' tags
		int fct = getFacesCount();
		for(int i = 0; i < fct; i++){
			getFaceAt(i)->clearSurfaceData();
		}
		// and in edges' tags
		for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge())
			it.getEdge()->clearLocalCurve();
	}
}

/// Search for local reparameterization curves (for boundary contours)
//int MeshContainer3dSurface::identifyLocalCurves(Metric3dContext& mc, double tolerance, 
//	TagExtended::TagType tag_type, int tag_value1, int tag_value2)
//{
//	int pct = getPointsCount();
//	for(int i = 0; i < pct; i++){
//		MeshPoint3d* point = getPointAt(i);
//		if(!point->isBorder()) continue;
//		if((tag_type != TagExtended::TAG_NONE) && !point->hasAllIntFlags(tag_type, tag_value1 | tag_value2)) continue;
//		for(int j = 0; j < point->getRank(); j++){
//			MeshEdge3d* edge = point->getEdge(j);
//			if( edge->isBorder() && !edge->hasLocalCurve() ) {
//				if(! approximateLocalCurve(mc, edge, tolerance)){
//					LOG4CPLUS_WARN(MeshLog::logger_console, "Failed to approximate local curve for edge at point #", point->getIndex());
//				}
//			}
//		}	
//	}
//	return m_local_curves.countInt();
//}

// Search for local reparameterization surfaces (for whole mesh, or only faces with the given tag)
//int MeshContainer3dSurface::identifyLocalSurfaces(Metric3dContext& mc, double tolerance, TagExtended::TagType tag_type, int tag_value)
//{
//	static const double TOL_F[] = { 0.5, 0.8, 1.0, 1.5, 2.0 }; // tolerance coefficient
//	static const int    LAY_F[] = {   3,   3,   3,   2,   1 }; // number of topological layers
//	static const double MTR_F[] = { 2.0, 1.5, 1.0, 0.5, 0.0 }; // metric radius for geometrical layer
//	static const int TOL_F_CT = sizeof(TOL_F) / sizeof(double);
//
//	int fct = getFacesCount();
//	int fct_left = fct;
//
//	for(int it = 0; (fct_left > 0) && (it < TOL_F_CT); it++){
//
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "===> Identifying local surfaces for TOL_F= " << TOL_F[it] << " LAY_F= " << LAY_F[it] << " ====");
//
//#ifdef T_DEBUG_
//		{
//			ostringstream ostr;
//			ostr << "Identifying local surfaces for TOL_F= " << TOL_F[it] << " LAY_F= " << LAY_F[it] << " ...";
//			LOG4CPLUS_INFO(MeshLog::logger_console,ostr.str());
//		}
//#endif // T_DEBUG_
//
//		// First, try match a surface for ALL points/faces (with the given tag)
//		fct_left -= approximateLocalSurface(mc, nullptr, TOL_F[it] * tolerance, LAY_F[it], MTR_F[it], tag_type, tag_value);
//		if(fct_left == 0) break;
//
//		// TODO? 
//		//		- For each patch of faces (surrounded by border edges)
//		//		- Start from the middle of each patch? Or other approach???
//
//		// Browse faces - if any left
//		for(int i = 0; i < fct; i++){
//			MeshFace* face = getFaceAt(i);
//			if(face->hasLocalSurface()) continue;
//			if(tag_type != TagExtended::TAG_NONE && !face->hasAnyIntFlags(tag_type, tag_value))
//				continue;
//			int ct = approximateLocalSurface(mc, face, TOL_F[it] * tolerance, LAY_F[it], MTR_F[it], tag_type, tag_value);
//			if(ct > 0) {
//				fct_left -= ct;
//#ifdef T_DEBUG_
//				{
//					ostringstream ostr;
//					ostr << "Local surfaces set (" << m_local_surfaces.countInt() << ") for " << (fct-fct_left) << " out of " << fct << " faces ...";
//					LOG4CPLUS_INFO(MeshLog::logger_console,ostr.str());
//				}
//#endif // T_DEBUG_
//				if(fct_left == 0) break;
//			}
//		}
//	}
//
//	// TODO? Finally, if any point left without any local surface, try match for this point all available surfaces...
//
//	if(fct_left > 0) {
//		LOG4CPLUS_WARN(MeshLog::logger_console, "Failed to approximate surface for some faces, fct_left", fct_left);
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Failed to approximate surface for some faces, fct_left= " << fct_left);
//	}
//
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Created local surfaces set: " << m_local_surfaces.countInt() << " sufaces.");
//
//	return m_local_surfaces.countInt();
//}

Curve3dConstPtr MeshContainer3dSurface::approximateLocalCurveForChain( Metric3dContext& mc, 
		SurfaceConstPtr surface, const DataVector<MeshEdge3d*> chain_edges, 
		const DataVector<MeshPoint3d*> chain_points, double tolerance ) 
{
	int cect = chain_edges.countInt();
	int cpct = chain_points.countInt();
	assert( cpct == cect+1 );
	DataVector<double> cparams(cpct);
	Curve3dPtr curve;

	if( surface != nullptr ) {
		DataVector<DPoint2d> points2d(cpct);

		for(int j = 0; j < cpct; j++)
			points2d.add( chain_points[j]->getLocalSurfaceParam(surface) );

		curve = Curve3dParametric::fitCurveOnSurface( mc, surface, points2d, cparams, tolerance );
	}else{
		// use 3D bspline
		DataVector<DPoint3d> dmpoints(cpct);
		for(int i = 0; i < cpct; i++)
			dmpoints.add(chain_points[i]->getCoordinates());

		assert(  chain_points[0]->isBorder(TagBorder::CORNER) && chain_points.last()->isBorder(TagBorder::CORNER) );

		DataVector< Curve3dPtr> cand_curves;

		cand_curves.add( Curve3dBSpline::fitToPointSequence( dmpoints, mc, tolerance ) );
		DLinearQuadric3d lquadric3d;
		if( DLeastSquaresFitting::fitLinearQuadric( dmpoints, lquadric3d, false ) >= 0.0 ) {
			cand_curves.add( std::make_shared<Curve3dLinearQuadric>( lquadric3d ) );
		}
		// + fit quadric-curve-3d (especially if pct < 5)?

		int n = cand_curves.countInt();
		assert( n > 0 );
		DataVector< DataVector<double> > cand_params(n, DataVector<double>());
		DataVector< double > cand_max_dist2(n, 0.0); 
		DataVector< double > cand_min_param(n, -LARGE_NUMBER );
		DataVector< double > cand_max_param(n,  LARGE_NUMBER );

		for(int k = 0; k < cpct; k++){
			const DPoint3d& pt_3d = chain_points[k]->getCoordinates();
			mc.countMetricAtPoint( pt_3d );

			for(int j = 0; j < n; j++ ){
				double last_t = cand_params[j].empty() ? 0.0 : cand_params[j].last();
				double t = cand_curves[j]->getParameter( chain_points[k]->getCoordinates(), last_t );
				if( k == 0 ) {
					cand_min_param[j] = t; 
					cand_max_param[j] = t;
				}else {
					if( t < cand_min_param[j] ) cand_min_param[j] = t;
					if( t > cand_max_param[j] ) cand_max_param[j] = t;
				}
				cand_params[j].add(t); 
				double d = mc.transformRStoMS( cand_curves[j]->getPoint( t ) - pt_3d ).length2();
				if(d > cand_max_dist2[j]) cand_max_dist2[j] = d;
			}
		}

		int best_curve = 0;
		// -> select best ?
		double tol2 = tolerance*tolerance;
		for(int i = 1; i < n; i++)
			if( cand_max_dist2[i] < cand_max_dist2[best_curve] ) best_curve = i;
//		if( cand_max_dist2[best_curve] > tol2 ) best_curve = -1; // none is good enough

		if( best_curve >= 0 ) {
			cparams.clear();
			for(int k = 0; k < cpct; k++){
				cparams.add( cand_params[best_curve][k] );
			}
		}

		if( best_curve >= 0 ){
			cand_curves[best_curve]->setMinMaxParam( cand_min_param[best_curve], cand_max_param[best_curve] );
			LOG4CPLUS_INFO(MeshLog::logger_mesh, 
				" * Local curve best=#" << best_curve << " for " << cpct << "pts (" 
				<< cand_curves[best_curve]->getSimpleDescription() 
				<< "), max_dist=" << sqrt(cand_max_dist2[best_curve]));
			curve = cand_curves[best_curve];
		}
	}

	if(curve){
		curve->setFixed();
		m_local_curves.add(curve);

		for(int i = 0; i < cpct; i++){
			chain_points[i]->setLocalCurve( curve, cparams[i] );
		}
		for(int i = 0; i < cect; i++){
			MeshEdge3d* cedge = chain_edges[i];
			assert(!cedge->hasLocalCurve());
			cedge->setLocalCurve(curve);
		}
	}

	return curve;
}

/// Search for local reparameterization curve for the given (boundary) edge, within the local surface patch
Curve3dConstPtr MeshContainer3dSurface::approximateLocalCurve(Metric3dContext& mc, 
			MeshEdge3d* edge, SurfaceConstPtr surface, 
			DataHashTable<MeshEdge3d*> * visited_edges, double tolerance)
{
	assert( edge->isBorder() );
	assert( (surface == nullptr) || edge->getMeshPoint(0)->hasLocalSurface( surface ) );
	assert( (surface == nullptr) || edge->getMeshPoint(1)->hasLocalSurface( surface ) );

	// try with parametric curve in parametric surface
	DataVector<MeshEdge3d*> chain_edges;
	MeshPoint3d* ccpoint = gatherBorderContourChainForSurface(edge, surface, chain_edges);

	assert ( ! chain_edges.empty() );

	int cect = chain_edges.countInt();

	DataVector<MeshPoint3d*> chain_points(cect+1);
	chain_points.add(ccpoint);
	for(int j = 0; j < cect; j++) {
		MeshEdge3d* cedge = chain_edges[j];
		chain_points.add(ccpoint = cedge->getOtherPoint(ccpoint));
		if(visited_edges != nullptr) 
			visited_edges->insert( cedge );
	}

	MeshPoint3d* cp_first = chain_points[0];
	MeshPoint3d* cp_last  = chain_points.last();
	assert( cp_first != cp_last );

	bool full_chain = ( cp_first->isBorder(TagBorder::CORNER) && cp_last->isBorder(TagBorder::CORNER) );
	if( full_chain )
		return approximateLocalCurveForChain( mc, surface, chain_edges, chain_points, tolerance );
	else
		return nullptr;
}

// Search for local reparameterization surface for the given face
//   and the mesh in the vicinity (without crossing boundaries!)
SurfaceConstPtr MeshContainer3dSurface::approximateLocalSurface(Metric3dContext& mc, 
		const DataVector<MeshFace*> &faces, int & ascribed_faces,
		double tolerance, int top_layers, double metric_radius,
		TagExtended::TagType tag_type, int tag_value)
{
//	START_CLOCK("lsurf-begin");

	ascribed_faces = 0;
	// 1. Gather local neigborhood (without boundary and crossing boundary)
	int pct = getPointsCount();
	int fct = getFacesCount();

	DataVector<int> p_layers(pct, -1);
	DataVector<int> f_layers(fct, -1);

	DataVector<int> f_list(fct);
	DataVector<int> p_list(pct);

	DataVector<int> f_layers_count;

	//static const int LAYERS_COUNT = 3;
	//static const double GEOM_METRIC_RADIUS = 2.0;

	gatherFirstLayerTopological(faces, f_layers_count, p_list, p_layers, f_list, f_layers, tag_type, tag_value);

	bool other_layers_available = false;

	if( faces.notEmpty() ){
		// a) topological (neighbours of neighbors...)
		other_layers_available = gatherLayeredVerticesTopological(f_layers_count, 
			p_list, p_layers, f_list, f_layers, top_layers, tag_type, tag_value);
		// b) geometrical (by neighbours, but using distance in metric space)
		if((metric_radius > 0.0) && other_layers_available) {
			DPoint3d middle = DPoint3d::zero;
			double w = 1.0 / faces.countInt();
			for(int i = 0; i < faces.countInt(); i++)
				middle.add( faces[i]->getMiddlePoint(), w );

			other_layers_available = gatherLayeredVerticesGeometrical(mc, middle, 
				f_layers_count, p_list, p_layers, f_list, f_layers, metric_radius, tag_type, tag_value);
		}
	}

	assert(f_layers_count.countInt() >= 2);
	int l_outside = 1+f_layers[ f_list[ f_layers_count[ f_layers_count.countInt()-2] ] ];

//	STOP_CLOCK("lsurf-begin");
	// 2. Try fitting surface  (tolerance in local metric)

//	START_CLOCK("lsurf-layers-loop");
	SurfacePtr temp_surface;

	int old_outside = l_outside;
	int surf_iterations = 0;
	int lpfct = p_list.countInt() + f_list.countInt();
	DataVector< double > approx_quality_1(lpfct), approx_quality_2(lpfct);
	DataVector< double > * approx_quality[] = { &approx_quality_1, &approx_quality_2 };
	DataVector<DPoint3d> local_points_1(lpfct), local_points_2(lpfct);
	DataVector<DPoint3d> * local_points[] = { &local_points_1, &local_points_2 };
	DataVector<DPoint2d> local_params_1(lpfct), local_params_2(lpfct);
	DataVector<DPoint2d> * local_params[] = { &local_params_1, &local_params_2 };
	DataVector<DVector3d> local_normals_1(lpfct), local_normals_2(lpfct);
	DataVector<DVector3d> * local_normals[] = { &local_normals_1, &local_normals_2 };

	int fit_index = 1;
	while(true){

		//showDebugLocalFaces( " MC3dS:: surf_iter #"+to_string(surf_iterations), p_list, p_layers, f_list, f_layers, l_outside );

		fit_index = 1 - fit_index;

		local_points[fit_index]->clear();
		for(int i = 0; i < p_list.countInt(); i++) {
			int pid = p_list[i];
			if(p_layers[pid] < l_outside) {
				MeshPoint3d* point = getPointAt(pid);
				local_points[fit_index]->add( point->getCoordinates() );
				local_normals[fit_index]->add( point->getBaseNormal() );
			}
		}

		// add extra points -> middle of faces
		for(int i = 0; i < f_list.countInt(); i++) {
			int fid = f_list[i];
			if(f_layers[fid] < l_outside) {
				MeshFace* face = getFaceAt(fid);
				local_points[fit_index]->add( face->getMiddlePoint() );
				local_normals[fit_index]->add( face->getBaseNormal() );
			}
		}
		assert(local_points[fit_index]->countInt() >= 3); // since starting from face

		if( local_points[fit_index]->countInt() < 10 ) { // add even more extra points -> middle of edges
			DataHashTable<MeshEdge3d*> hedges(100, nullptr);
			for(int i = 0; i < f_list.countInt(); i++) {
				int fid = f_list[i];
				if(f_layers[fid] < l_outside) {
					MeshFace* face = getFaceAt(fid);
					int fect = face->getEdgeCount();
					for(int j = 0; j < fect; j++){
						MeshEdge3d* edge = face->getEdge(j);
						if( hedges.insert( edge ) ){
							local_points[fit_index]->add( edge->getPoint(0.5) );
							local_normals[fit_index]->add(
								(edge->getMeshPoint(0)->getBaseNormal() + 
								 edge->getMeshPoint(1)->getBaseNormal() ).normalized() );
						}
					}
				}
			}
		}
//		local_count[fit_index] = local_points[fit_index]->countInt();

		auto fit_surface = SurfaceParametric::fitSurface(mc, 
			*local_points[fit_index], *local_normals[fit_index], *local_params[fit_index],
			tolerance, faces.empty() ? getFaceAt(0) : faces[0],
			0.1, approx_quality[fit_index]);

		if(fit_surface){
			temp_surface = fit_surface;
		}else if(temp_surface) {
			// revert
			l_outside = old_outside;
			fit_index = 1 - fit_index;
			break;
		}
		if(!temp_surface) return nullptr;

#ifdef SHOW_GATHER_LAYERED
		if(false){
			MeshViewSet* set = new MeshViewSet;
			temp_surface->createViewSetForPoints( set, *local_params[fit_index] );
			ostringstream ostr;
			ostr << "surface fit for tolerance " << tolerance;
			SHOW_MESH(ostr.str(), set);
		}

		if(false){
//		if( face->getIndex() == 153 ){
			MeshViewSet* set = new MeshViewSet;
			for(int i = 0; i < f_list.countInt(); i++) {
				int fid = f_list[i];
				MeshFace* sface = getFaceAt( fid );
				set->addFaceWithEdges( sface, f_layers[fid] );
			}
			for(int i = 0; i < p_list.countInt(); i++)
				set->addPoint( getPointAt(p_list[i]), p_layers[p_list[i]]);
			SHOW_MESH("Layers before MC3dS::gatherLayeredVerticesForSurface", set);
		}
#endif

		++surf_iterations;

		if( faces.notEmpty() ) {
			if(other_layers_available){
				other_layers_available = gatherLayeredVerticesForSurface(mc, temp_surface, 
											f_layers_count, p_list, p_layers, 
											f_list, f_layers, tolerance, tag_type, tag_value);
				int new_l_outside = 1+f_layers[ f_list[ f_layers_count[ f_layers_count.countInt()-2] ] ];
				if(new_l_outside > l_outside){
					old_outside = l_outside;
					l_outside = new_l_outside;
				}else
					break;  // no additional points found -> stop
			}else
				break; // no additonal points for sure -> stop
		}else
			break; // no starting point -> all points already included -> stop
	}

	//if(true) {
	//	int p_list_count = 0;
	//	for( int i = 0 ; i < p_list.countInt(); i++ )
	//		if( p_layers[ p_list[i] ] < l_outside ) p_list_count++;
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh, " p_list_count=" << p_list_count << ", local_points=" << local_points.countInt() << ", local_count=" << local_count);
	//}

	DataVector<DPoint3d> & lpoints = *local_points[fit_index];
	DataVector<DPoint2d> & lparams = *local_params[fit_index];
	DataVector<double> & laquality = *approx_quality[fit_index];


	cleanAndPackPFLists(l_outside, p_list, p_layers, f_list, f_layers,
		lpoints, lparams, laquality);

	DataVector<int> pref( pct, -1 );
	int lpct = p_list.countInt();
	for(int i = 0; i < lpct; i++) {
		int pid = p_list[i];
		assert( p_layers[ pid ] < l_outside );
		pref[ pid ] = i;
	}

	fixLocalSurfaceOrientation(temp_surface, f_list, pref, lpoints, lparams); // ?

	DBox sbox;
	for(int i = 0; i < lpct; i++) {
		sbox.addPoint( lpoints[i] );
		sbox.addPoint( temp_surface->getPoint( lparams[i] ) );
	}
	mc.countMetricAtPoint( sbox.getMiddlePoint() );
	sbox.growDM( mc, 0.5 );

	if( faces.notEmpty() && other_layers_available ) {
		createLocalSurfaceDomain(mc, temp_surface, sbox, tolerance, f_list, pref, 
			lpoints, lparams, laquality );
	}else{
		for(int i = 0; i < lpct; i++) 
			laquality[i] = 1.0;
	}

	temp_surface->setIntTag( tag_type, tag_value ); // mark surface with local-sub-domain id...
	temp_surface->setFixed( true );

	SurfaceConstPtr surface = temp_surface;
	m_local_surfaces.add(surface);	


//	START_CLOCK("lsurf-set");
	ascribed_faces += ascribeLocalSurface( mc, surface, sbox, f_list, p_list, lparams, laquality ); 

	int cct = createContoursForLocalSurface(mc, surface, tolerance, f_list);

	ostringstream log_info;
	log_info << "NP= " << (lpct) << " NF= " << ascribed_faces << " ITER= " << surf_iterations;
	if( cct > 0 ) log_info << " CONTOURS= " << cct;
	LOG4CPLUS_INFO(MeshLog::logger_mesh, log_info.str());

// ORIENTATION ???
	fixLocalSurfaceOrientationForInvertedFaces(surface, f_list, pref, lparams );

//	STOP_CLOCK("lsurf-set");

	return surface;
}

// Search for local reparameterization surface for the given face
//   and the mesh in the vicinity (without crossing boundaries!)
SurfaceConstPtr MeshContainer3dSurface::approximateLocalSurfaceViaNormals(
		Metric3dContext& mc, const DataVector<MeshFace*> &faces, int & ascribed_faces,
		double sp_min, double tolerance, TagExtended::TagType tag_type, int tag_value)
{
//	START_CLOCK("lsurf-begin");

	//- face

	assert( faces.notEmpty() );
	ascribed_faces = 0;

	//- base_plane (from base_normal)

	DVector3d ave_normal = faces[0]->getBaseNormal();
	if( faces.countInt() > 1 ) {
		for(int i = 1; i < faces.countInt(); i++) 
			ave_normal += faces[i]->getBaseNormal();
		ave_normal.normalize();
	}
	SurfaceConstPtr base_plane( new SurfacePlane (faces[0]->getMiddlePoint(), ave_normal ));

	//- vicinity_via_normals: gather local neigborhood (without boundary and crossing boundary)

	int pct = getPointsCount();
	int fct = getFacesCount();

	DataVector<int> p_layers(pct, -1);
	DataVector<int> f_layers(fct, -1);

	DataVector<int> f_list(fct);
	DataVector<int> p_list(pct);

	DataVector<int> f_layers_count;

	gatherFirstLayerTopological(faces, f_layers_count, p_list, p_layers, f_list, f_layers, tag_type, tag_value);

	double actual_sp_min = gatherLayeredVerticesViaNormals( 
		base_plane, sp_min, 
		f_layers_count, p_list, p_layers, f_list, f_layers, 
		tag_type, tag_value);
	assert(f_layers_count.countInt() >= 2);
	int l_outside = 1+f_layers[ f_list[ f_layers_count[ f_layers_count.countInt()-2] ] ];

	DataVector<MeshFace*> plane_faces( f_list.countInt() );
	for(int i = 0; i < f_list.countInt(); i++) {
		int fid = f_list[i];
		if(f_layers[fid] < l_outside)
			plane_faces.add( getFaceAt( fid ) );
	}

	//- base_surface (for a wider set of faces, gathered for plane-projection)

	double model_diameter = getBoundingBoxDiameter();
	double cmin_eps = ControlSpace2dAdaptive::param_curvature_ratio / ( 0.5 * model_diameter );
	auto base_surface = SurfaceParametric::calculateBestBaseSurface( plane_faces, cmin_eps );

	//- check_fit (if good enough - STOP)
	 //- vicinity_via_normals
	 //- quadric_on_surface
	 //- check_fit
	 //- DONE

	if( base_surface->getType() != ElementType::SURFACE_PLANE ) {

		// check normal-fit of plane_faces for base_surface ...
		if(false){

			DataVector<bool> check_ok[3];

			for(int i = 0; i < plane_faces.countInt(); i++) {
				MeshFace* face = plane_faces[i];

				// check normal
				double face_sp = base_surface->getNormalVectorForPoint3d( face->getMiddlePoint() 
								).scalarProduct( face->getBaseNormal() );
				check_ok[0].add( face_sp >= sp_min );

				// check param and point-normal
				bool all_params_in_range = true;
				double min_point_sp = 1.0;
				int ofpct = face->getPointCount();
				for(int k = 0; all_params_in_range && (k < ofpct); k++){
					MeshPoint3d* ep = face->getPoint(k);
					all_params_in_range = base_surface->withinParamRange( 
						base_surface->getParameters( 
							ep->getCoordinates() ) );
					double point_sp = base_surface->getNormalVectorForPoint3d( ep->getCoordinates() 
								).scalarProduct( face->getBaseNormal() );
					if( point_sp < min_point_sp ) min_point_sp = point_sp;
				}
				
				check_ok[1].add( min_point_sp >= sp_min );
				check_ok[2].add( all_params_in_range );
			}

			string labels[3] = { "face_sp", "point_sp", "param" };
			for(int l = 0; l < 3; l++){
				MeshViewSet* set = new MeshViewSet;
				DataVector<DPoint2d> params( plane_faces.countInt() * 3 );
				int fail_count = 0;
				for(int i = 0; i < plane_faces.countInt(); i++) {
					MeshFace* face = plane_faces[i];
					bool ok = check_ok[l][i];
					if(!ok) fail_count++;
					set->addFaceWithEdgesAndPoints( face, ok ? 0 : 1 );
					int fpct = face->getPointCount();
					for(int j = 0; j < fpct; j++)
						params.add( base_surface->getParameters( face->getPoint(j)->getCoordinates()) );
				}
				set->addInfo("check_type", labels[l]);
				set->addInfo("fail count", fail_count);
				set->addInfo("base surface", base_surface->getSimpleDescription());
				base_surface->createViewSetForPoints(set, params);
				SHOW_MESH("plane-faces back-check for base-surface", set);
			}
		}

		// ... for setting center value for param.u and range checking...

		p_layers.setAll( -1 );
		f_layers.setAll( -1 );

		f_list.clear();
		p_list.clear();
		f_layers_count.countInt();

		// faces: oryginal or updated after first gatherLayered... ?                                        
		//gatherFirstLayerTopological(faces, f_layers_count, p_list, p_layers, f_list, f_layers, tag_type, tag_value);
		gatherFirstLayerTopological(plane_faces, f_layers_count, p_list, p_layers, f_list, f_layers, tag_type, tag_value);

		actual_sp_min = gatherLayeredVerticesViaNormals( 
			base_surface, sp_min, 
			f_layers_count, p_list, p_layers, f_list, f_layers, 
			tag_type, tag_value);

		assert(f_layers_count.countInt() >= 2);
		l_outside = 1+f_layers[ f_list[ f_layers_count[ f_layers_count.countInt()-2] ] ];
	}

	int lpfct = p_list.countInt() + f_list.countInt();
	DataVector<DPoint3d> lpoints(lpfct);
	DataVector<DPoint2d> lparams(lpfct);
	DataVector<DVector3d> lnormals(lpfct);
	DataVector< double > laquality(lpfct);

	cleanAndPackPFLists(l_outside, p_list, p_layers, f_list, f_layers,
		lpoints, lparams, laquality);

	assert( lpoints.empty() );
	assert( lparams.empty() );

	int lpct = p_list.countInt();
	int lfct = f_list.countInt();

	for(int i = 0; i < lpct; i++) {
		int pid = p_list[i];
		assert( p_layers[pid] < l_outside);
		MeshPoint3d* point = getPointAt(pid);
		lpoints.add( point->getCoordinates() );
		lnormals.add( point->getBaseNormal() );
	}

	DataVector< MeshFace* > sfaces( lfct );
	for(int i = 0; i < lfct; i++) sfaces.add( getFaceAt( f_list[i] ) );

	SurfacePtr temp_surface;
	if( actual_sp_min > 0.9 ) {
		DataVector< MeshPoint3d* > spoints( lpct );
		for(int i = 0; i < lpct; i++) spoints.add( getPointAt( p_list[i] ) );

		temp_surface  = base_surface->adjustedForFaces( sfaces, spoints, faces[0] );
		assert( lparams.empty() ); // since they were calculated for other surface
		SurfaceParametric::checkSurfaceFit( mc, temp_surface, lpoints, lnormals, lparams, tolerance, &laquality, true);
	}

	if( !temp_surface ) {
		// 2. Try fitting quadric-on-surface  (tolerance in local metric)

		DRect brect;
		for(int i = 0; i < p_list.countInt(); i++) {
			const DPoint3d& dpt = getPointAt( p_list[i] )->getCoordinates();
			DPoint2d dparam = base_surface->getParameters( dpt );
			lparams.add( dparam );
			brect.addPoint( dparam );
		}
		brect.inflate(0.05);

		temp_surface = SurfaceQuadricQTree::fitQuadricQTreeSurface(mc,
			base_surface, brect, tolerance, this, sfaces);

		if (true) {
			MeshViewSet* set = new MeshViewSet;
			temp_surface->createViewSetForPoints(set, lparams);
			SHOW_MESH("set 1", set);
		}

		SurfaceParametric::checkSurfaceFit( mc, temp_surface, lpoints, lnormals, lparams, tolerance, &laquality, false);
		// ... check? but no tolerance and no surface-removing ...
		if (true) {
			MeshViewSet* set = new MeshViewSet;
			temp_surface->createViewSetForPoints(set, lparams);
			SHOW_MESH("set 3", set);
		}

	}

	if(temp_surface == nullptr) return nullptr;

	if(true){
		MeshViewSet* set = new MeshViewSet;
		temp_surface->createViewSetForPoints( set, lparams );
		set->addInfo("tolerance", tolerance);
		set->addInfo("surface-type", temp_surface->getSimpleDescription() );
		SHOW_MESH("surface fit - via normals", set);
	}

	DataVector<int> pref( pct, -1 );
	for(int i = 0; i < lpct; i++) {
		int pid = p_list[i];
		assert( p_layers[ pid ] < l_outside );
		pref[ pid ] = i;
	}

	fixLocalSurfaceOrientation(temp_surface, f_list, pref, lpoints, lparams); // ?

	DBox sbox;
	for(int i = 0; i < lpct; i++) {
		sbox.addPoint( lpoints[i] );
		sbox.addPoint( temp_surface->getPoint( lparams[i] ) );
	}
	mc.countMetricAtPoint( sbox.getMiddlePoint() );
	sbox.growDM( mc, 0.5 );

	// no domain necessary - it is built-in, but add OBB !!

	//DOrientedBox sobox = temp_surface->getOrientedBox( lpoints );
	//DOrientedBox sobox_rot = temp_surface->getOrientedBoxOpt( lpoints );

	//if( sobox.getVolume() > sobox_rot.getVolume() )
	//	sobox = sobox_rot;
	//sobox.growDM( mc, 0.5 );

	{ // set marks for local-surface-inner-edges
		int lfct = f_list.countInt();
		// triangle->edges
		DataHashTableKeyValue< MeshEdge3d*, int > hedges( lfct*3, nullptr );
		for(int i = 0; i < lfct; i++) {
			MeshFace* f = getFaceAt( f_list[i] );
			int ect = f->getEdgeCount();
			for(int j = 0; j < ect; j++)
				hedges.incValue( f->getEdge(j), 1 );
		}

		DataVector<MeshEdge3d*> hkeys;
		hedges.getKeys( hkeys );
		for(int i = 0; i < hkeys.countInt(); i++){			
			if( hedges.getValue( hkeys[i], 0 ) == 2 ){ // edge inner for local-surface-domain
				hkeys[i]->setIntTag( TagExtended::TAG_LOCAL_SURFACE_INNER, 1 );
			}
		}
	}

	temp_surface->setIntTag( tag_type, tag_value ); // mark surface with local-sub-domain id...
	temp_surface->setFixed( true );

	SurfaceConstPtr surface = temp_surface;
	m_local_surfaces.add(surface);	

//	START_CLOCK("lsurf-set");
	ascribed_faces += ascribeLocalSurface( mc, surface, sbox, f_list, p_list, lparams, laquality ); 

	int cct = createContoursForLocalSurface(mc, surface, tolerance, f_list);
	{
		ostringstream log_info;
		log_info << "NP= " << (lpct) << " NF= " << ascribed_faces;
		if (cct > 0) log_info << " CONTOURS= " << cct;
		LOG4CPLUS_INFO(MeshLog::logger_mesh, log_info.str());
	}
// ORIENTATION ???
	fixLocalSurfaceOrientationForInvertedFaces(surface, f_list, pref, lparams );

	//	STOP_CLOCK("lsurf-set");

	return surface;
}

void MeshContainer3dSurface::createLocalSurfaceDomain(
		Metric3dContext& mc,
		SurfacePtr surface,
		const DBox& sbox,
		double tolerance,
		const DataVector<int> & f_list,
		const DataVector<int> & pref,
		const DataVector<DPoint3d>& lpoints, 
		const DataVector<DPoint2d>& lparams,
		const DataVector<double>& laquality)
{
//	START_CLOCK("lsurf-domain");

	int lfct = f_list.countInt();
	// triangle->edges
	DataHashTableKeyValue< MeshEdge3d*, int > hedges( lfct*3, nullptr );
	for(int i = 0; i < lfct; i++) {
		MeshFace* f = getFaceAt( f_list[i] );
		int ect = f->getEdgeCount();
		for(int j = 0; j < ect; j++)
			hedges.incValue( f->getEdge(j), 1 );
	}

	DataVector<MeshEdge3d*> hkeys;
	hedges.getKeys( hkeys );
	for(int i = 0; i < hkeys.countInt(); i++){			
		if( hedges.getValue( hkeys[i], 0 ) == 2 ){ // edge inner for local-surface-domain
			hkeys[i]->setIntTag( TagExtended::TAG_LOCAL_SURFACE_INNER, 1 );
		}
	}

	int lpct = lparams.countInt();
	// Create and attach domain info ...

	DataVector< DVector3d > poffs(lpct, DVector3d::zero );

	SurfaceDomain* sdomain = nullptr;

	if( param_surface_domain_type == MeshContainer3dSurface::SDOMAIN_HULL ) { // if hull
		SurfaceDomainHull* shd = new SurfaceDomainHull;

		for(int i = 0; i < lpct; i++) shd->addPoint( lparams[i] );

		shd->createHull();
		shd->createHullDist( mc, surface.get() );
		sdomain = shd;
	}else if( param_surface_domain_type == MeshContainer3dSurface::SDOMAIN_FACES) { // if faces
		// local_points -> coord

		// local_faces -> 3xindices of vertices
		DataVector< ITriangle > local_faces( lfct );
		DataVector<double> local_dist2( lfct );
		double tol2 = tolerance * tolerance;
		//showDebugLocalFaces( " MC3dS:: -> domain faces", p_list, p_layers, f_list, f_layers, l_outside, &pref);

		for(int i = 0; i < lfct; i++) {
			int fid = f_list[i];
			MeshFace* f = getFaceAt( fid );
			mc.countMetricAtPoint( f->getMiddlePoint(), true);
			//const DVector3d fn = f->getNormalVector();
			//double mlen2 = tol2 * fn.length2() / mc.transformRStoMS( fn ).length2();
			double mlen2 = tol2 * sqr( mc.minLength() );

			int fpct = f->getPointCount();
			if( fpct == 3 ){
				local_dist2.add( mlen2 );
				ITriangle itri(
					pref[ f->getPoint(0)->getIndex() ],
					pref[ f->getPoint(1)->getIndex() ],
					pref[ f->getPoint(2)->getIndex() ] );
				local_faces.add( itri );
			}else{
				DataVector<MeshFace*> split_faces;
				if( f->splitToTriangles(split_faces) ) {
					for(int j = 0; j < split_faces.countInt(); j++) {
						MeshFace* sf = split_faces[j];
						assert( sf->getPointCount() == 3 );
						local_dist2.add( mlen2 );
						local_faces.add( ITriangle(
							pref[ sf->getPoint(0)->getIndex() ],
							pref[ sf->getPoint(1)->getIndex() ],
							pref[ sf->getPoint(2)->getIndex() ] ) );
						delete sf;
					}
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "Failed splitting poly-face to triangles...");
				}
			}
		}
		sdomain = new SurfaceDomainFaces( lpoints, local_faces, local_dist2, laquality );
	}else if( param_surface_domain_type == MeshContainer3dSurface::SDOMAIN_FACES_PLANAR) { // if faces
		// local_points -> coord

		// identify border-inner faces
		// edge->points
		DataVector<bool> bpoints(getPointsCount(), false);
		for(int i = 0; i < hkeys.countInt(); i++){			
			if( hedges.checkValue( hkeys[i], 1 ) ){
				bpoints[ hkeys[i]->getMeshPoint(0)->getIndex() ] = true;
				bpoints[ hkeys[i]->getMeshPoint(1)->getIndex() ] = true;
			}
		}

		// for each face -> for each border edge -> calculate extension vector [n x e], and add to vertices
		for(int i = 0; i < lfct; i++) {
			int fid = f_list[i];
			MeshFace* f = getFaceAt( fid );
			DVector3d n = f->getNormalVector();
			int fpct = f->getPointCount();
			for(int j = 0; j < fpct; j++) {
				MeshPoint3d* mp0 = f->getPoint(j);
				MeshPoint3d* mp1 = f->getNextPoint( mp0 );
				MeshEdge3d* edge = f->getEdge(j);
				assert( edge == mp0->getEdgeToPoint(mp1) );
				if( hedges.checkValue( edge, 1 ) ) {
					DVector3d v_edge = mp1->getCoordinates() - mp0->getCoordinates();
					assert( ! v_edge.isZero() );
					DVector3d v_off = v_edge.crossProduct(n).normalized();
					//if( edge->getFaceIndex(f) == 1 ) 
					//	v_off.reverse(); 
					poffs[ pref[mp0->getIndex()] ] += v_off;
					poffs[ pref[mp1->getIndex()] ] += v_off;
				}
			}
		}
		// for each border vertex -> normalize vertices with tolerance length
		for(int i = 0; i < lpct; i++) {
			DVector3d& voffs = poffs[i];
			if( !voffs.isZero() ) {
				double d = voffs.length();
				double dx = 2.0 / d;
				double fx = (dx < 2.0 ) ? 2.0 : (4.0/dx); // avoid too long offsets, for sharp edges...
				mc.countMetricAtPoint( lpoints[i] );
				double mlen = mc.transformRStoMS( voffs ).length();
				double req_len = fx * tolerance / mlen;
				voffs.normalize( req_len );
			}
		}
		if( false ) {
			MeshViewSet* set = new MeshViewSet;
			for(int i = 0; i < lfct; i++) {
				int fid = f_list[i];
				set->addFaceWithEdges( getFaceAt( fid ) );
			}
			for(int i = 0; i < lpct; i++) {
				const DVector3d& voffs = poffs[i];
				if( !voffs.isZero() ) {
					set->addEdge( lpoints[i], lpoints[i]+voffs, 2 );
				}
			}
			SHOW_MESH("border offset vectors", set );
		}

		// local_faces -> 3xindices of vertices
		DataVector< ITriangle > local_faces( lfct );
		DataVector< bool > local_faces_blayer( lfct );

		for(int i = 0; i < lfct; i++) {
			int fid = f_list[i];
			MeshFace* f = getFaceAt( fid );

			int fpct = f->getPointCount();
			if( fpct == 3 ){
				int p0i = f->getPoint(0)->getIndex();
				int p1i = f->getPoint(1)->getIndex();
				int p2i = f->getPoint(2)->getIndex();
				local_faces.add( ITriangle(	pref[ p0i ], pref[ p1i ], pref[ p2i ] ) );
				local_faces_blayer.add( bpoints[p0i] || bpoints[p1i] || bpoints[p2i] );
			}else{
				DataVector<MeshFace*> split_faces;
				if( f->splitToTriangles(split_faces) ) {
					for(int j = 0; j < split_faces.countInt(); j++) {
						MeshFace* sf = split_faces[j];
						assert( sf->getPointCount() == 3 );
						int p0i = sf->getPoint(0)->getIndex();
						int p1i = sf->getPoint(1)->getIndex();
						int p2i = sf->getPoint(2)->getIndex();
						local_faces.add( ITriangle( pref[ p0i ], pref[ p1i ], pref[ p2i ] ) );
						local_faces_blayer.add( bpoints[p0i] || bpoints[p1i] || bpoints[p2i] );
						delete sf;
					}
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "Failed splitting poly-face to triangles...");
				}
			}
		}
		DataVector< DPoint2d > sd_params(lpct);
		for(int i = 0; i < lpct; i++) {
			const DVector3d& voffs = poffs[i];
			sd_params.add( voffs.isZero() ?
				lparams[i] : surface->getParameters( lpoints[i] + voffs ) );
		}
		SurfaceDomainFacesPlanar * sdfp = 
			new SurfaceDomainFacesPlanar( sd_params, local_faces, laquality );
		if( sdfp->invCount() > 0 ) {
			MeshViewSet* set = new MeshViewSet;
			surface->createViewSetForPoints( set, sd_params );
			for(int i = 0; i < lfct; i++) {
				int fid = f_list[i];
				set->addFaceWithEdgesAndPoints( getFaceAt( fid ) );
			}
			SHOW_MESH("inv faces for surf-domain", set);
		}
		sdomain = sdfp;
	}

	assert( sdomain != nullptr); 

	DOrientedBox sobox = surface->getOrientedBox( lpoints );
	DOrientedBox sobox_rot = surface->getOrientedBoxOpt( lpoints );

	//int sobox_not_within = 0;
	//int sobox_rot_not_within = 0;
	//for(int i = 0; i < lpoints.countInt(); i++){
	//	const DPoint3d& dpt = lpoints[i];
	//	if( !sobox.contains( dpt ) ) ++sobox_not_within;
	//	if( !sobox_rot.contains( dpt ) ) ++sobox_rot_not_within;
	//}

	if( sobox.getVolume() > sobox_rot.getVolume() )
		sobox = sobox_rot;

	sobox.growDM( mc, 0.5 );

	//int sobox_grow_not_within = 0;
	//for(int i = 0; i < lpoints.countInt(); i++){
	//	const DPoint3d& dpt = lpoints[i];
	//	if( !sobox.contains( dpt ) ) ++sobox_grow_not_within;
	//}

	//if( sobox_not_within > 0 || sobox_rot_not_within > 0 || sobox_grow_not_within > 0 ) {
	//	LOG4CPLUS_WARN(MeshLog::logger_console, "obb-box not within", to_string(sobox_not_within) + "/" +
	//		to_string(sobox_rot_not_within) + "/" + to_string(sobox_grow_not_within));
	//}

	sdomain->setOBBox( sobox );
	surface->setDomain( sdomain );

//	STOP_CLOCK("lsurf-domain");
}

int MeshContainer3dSurface::ascribeLocalSurface( 
	Metric3dContext& mc, 
	SurfaceConstPtr surface, 
	const DBox& sbox,
	const DataVector<int> & f_list,
	const DataVector<int> & p_list,
	const DataVector<DPoint2d>& lparams,
	const DataVector<double>& laquality)
{
	// 4a. Include surface into faces
	int lfct = f_list.countInt();
	int ascribed_faces = 0;
	for(int i = 0; i < lfct; i++){
		int fid = f_list[i];
		MeshFace* face = getFaceAt( fid );
		if( ! face->hasLocalSurface()){

			ascribed_faces++;
			face->setLocalSurface(surface);
		}
	}

	// TODO? Check if this surface include any/some of other, less vast surfaces

	// 4. Include sid in tag-info for all "points_id"

	//if(acs && !acs->isAdaptive()) acs = nullptr;
	//DataHashTableKeyValue<SurfaceParametricSet*, SurfaceSetPtrHolder> local_surface_hash( 2*p_list.countInt(), (SurfaceParametricSet*)1);

	if(m_control && m_control->isAdaptive()) {
		// "zero"-value different from nullptr, since nullptr is also possible as a key
		DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> local_surface_hash( 
			(unsigned int)(128*p_list.countInt()), std::make_shared<SurfaceParametricSet>() );
		m_control->getAsAdaptive()->addLocalSurfaceAtBBox( sbox, surface, local_surface_hash );
	}

	int lpct = p_list.countInt();
	for(int i = 0; i < lpct; i++){
		int pid = p_list[i];
		MeshPoint3d* point = getPointAt(pid);

		//bool is_ok1 = point->checkLocalSurfaceParams();

		point->setLocalSurface( surface, lparams[i], laquality[i] );

		//if( true ) {
		//	ParamAndQuality pq1 = point->getLocalSurfaceParamQuality( surface );
		//	DPoint2d param = lparams[p_index];
		//	double q2 = surface->withinDomainQuality( point->getCoordinates(), param );
		//	double dq = std::abs( pq1.quality - q2 );
		//	if( dq > 0.1 ) {
		//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "xxx " << dq);
		//	}
		//}

		//bool is_ok2 = point->checkLocalSurfaceParams();
		//if( ! is_ok2 ) {
		//	MeshViewSet * set = new MeshViewSet;
		//	for(int i = 0; i < lpct; i++ )
		//		set->addPoint( getPointAt( p_list[i] ) );
		//	for(int i = 0; i < lpct; i++ )
		//		set->addPoint( lpoints[i] );
		//	for(int i = 0; i < lfct; i++ )
		//		set->addFaceWithEdges( getFaceAt( f_list[i] ) );
		//	for(int i = 0; i < lpct; i++) {
		//		const DVector3d& voffs = poffs[i];
		//		if( !voffs.isZero() ) {
		//			set->addEdge( lpoints[i], lpoints[i]+voffs, 2 );
		//		}
		//	}
		//	SHOW_MESH("p_list/f_list", set);
		//}
		//if( !is_ok1 || !is_ok2 ) { 
		//	DPoint2d check_param = lparams[p_index];
		//	bool check_valid = surface->withinDomain( point->getCoordinates(), check_param );
		//	LOG4CPLUS_INFO(MeshLog::logger_mesh, " check_valid -> " << check_valid);
		//}
		//assert(is_ok1 && is_ok2);
	}
	
	return ascribed_faces;
}

int MeshContainer3dSurface::createContoursForLocalSurface(
		Metric3dContext& mc, 
		SurfaceConstPtr surface, 
		double tolerance, 
		const DataVector<int> & f_list)
{
	int lfct = f_list.countInt();
	int cct = 0;
	// create contours for all border-edges in this local surface
	DataHashTable< MeshEdge3d* > visited_edges( 2*lfct, nullptr );
	for(int i = 0; i < lfct; i++) {
		int fid = f_list[i];
		MeshFace* face = getFaceAt( fid );
		int fect = face->getEdgeCount();
		for(int j = 0; j < fect; j++){
			MeshEdge3d* edge = face->getEdge(j);
			if( !edge->isBorder() || edge->hasLocalCurve() || visited_edges.contains(edge) ) continue;
			auto curve = approximateLocalCurve(mc, edge, surface, &visited_edges, tolerance);
			if( curve != nullptr ){
				cct++;
				//if(false){
				//	MeshViewSet* set = new MeshViewSet;
				//	surface->createViewSetForPoints( set, lparams );
				//	DataVector< DPoint3d > polyline;
				//	curve->getPolyLine(polyline);
				//	for(int i = 1; i < polyline.countInt(); i++)
				//		set->addEdge( polyline[i-1], polyline[i], 1 );
				//	set->addInfo("curve-type", curve->getSimpleDescription() );
				//	SHOW_MESH("Curve contour (full-chain) for surface", set);
				//}
			}
		}
	}
	return cct;
}

void MeshContainer3dSurface::fixLocalSurfaceOrientationForInvertedFaces(
		SurfaceConstPtr surface,
		DataVector<int> & f_list,
		const DataVector<int> & pref,
		const DataVector<DPoint2d>& lparams )
{
	// ... remove faces with other surface already set
	for(int i = f_list.countInt()-1; i >= 0; i--)
		if(getFaceAt(f_list[i])->getCurrentLocalSurface() != surface)
			f_list.removeAt(i);

	DataVector<int> face_info(getFacesCount(), 0); 
	// -1 -> face orientation inv, 1 -> face orientation ok, 
	// 0 -> face orientation unknown

	// set surface orientation for faces
	while(! f_list.empty() ){
		// choose best start face (closest to surface with respect to normal vector)
		int best_face_i = -1;
		double best_face_sc = 0.0;
		for(int i = 0; i < f_list.countInt(); i++){
			MeshFace* fface = getFaceAt(f_list[i]);
			const DVector3d fvn = fface->getNormalVector();
			DPoint2d fpt_2d;
			int fpct = fface->getPointCount();
			double fratio = 1.0 / fpct;
			for(int j = 0; j < fpct; j++)
				fpt_2d.add( lparams[ pref[ fface->getPoint(j)->getIndex() ] ], fratio );
			const DVector3d svn = surface->getNormalVector( fpt_2d );
			double sc = fvn.scalarProduct(svn);
			if((best_face_i < 0) || (std::abs(sc) > std::abs(best_face_sc))){
				best_face_i = i;
				best_face_sc = sc;
			}
		}
		// set orientation for best face
		face_info[f_list[best_face_i]] = ((best_face_sc >= 0.0) ? 1 : -1);
		// ... and propagate through vicinity
		DataVector<int> local_faces(f_list.countInt());
		local_faces.add(f_list[best_face_i]);
		int inverted_count = 0;
		for(int i = 0; i < local_faces.countInt(); i++){
			int fid = local_faces[i];
			MeshFace* fface = getFaceAt(fid);
			fface->setLocalSurfaceOrientation(face_info[fid]);
			// ... check if inverted
			const DVector3d fvn = fface->getNormalVector();
			DPoint2d fpt_2d;
			int fpct = fface->getPointCount();
			double fratio = 1.0 / fpct;
			for(int i = 0; i < fpct; i++)
				fpt_2d.add( lparams[ pref[ fface->getPoint(i)->getIndex() ] ], fratio );
			const DVector3d svn = surface->getNormalVector( fpt_2d );
			double sc = face_info[fid] * fvn.scalarProduct(svn);
			if(sc <= 0.0) inverted_count++;
			// ... check neighbours
			for(int j = 0; j < fface->getEdgeCount(); j++){
				MeshEdge3d* edge = fface->getEdge(j);
				if(edge->isBorder()) continue;
				MeshFace* other_face = edge->getOtherFace(fface);
				if( other_face->getCurrentLocalSurface() != surface ) 
					continue;
				int other_fid = other_face->getIndex();
				if(face_info[other_fid] == 0){ // i.e. if neighbour is unsure, use this face
					assert( ! other_face->hasLocalSurfaceOrientation() );
					face_info[other_fid] = (fface->sameOrientation(other_face) ? face_info[fid] : -face_info[fid]);
					local_faces.add(other_fid);
				}
			}
		}
		if(inverted_count > 0)
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Inverted faces: " << inverted_count << " out of local " << local_faces.countInt());
		// ... if most faces are inverted, reverse the orientation for all
		if( 2*inverted_count > local_faces.countInt() ){
			for(int i = 0; i < local_faces.countInt(); i++){
				getFaceAt(local_faces[i])->setLocalSurfaceOrientation(-face_info[local_faces[i] ]);
			}
		}
		// ... remove ready faces from surface_faces
		for(int i = f_list.countInt()-1; i >= 0; i--)
			if(std::abs(face_info[f_list[i]]) == 1)
				f_list.removeAt(i);

		//if(true){
		//	MeshViewSet* dset = new MeshViewSet(0, 10*fct, fct);
		//	for(int i = 0; i < fct; i++){
		//		MeshFace* face = getFaceAt(i);
		//		dset->addFaceWithEdges(face, std::abs(face_info[i]), MeshViewSet::param_shrink, face->getBlock(0) != nullptr);
		//	}
		//	SHOW_MESH("Orientation of faces", dset);
		//}
	}
}

void MeshContainer3dSurface::fixLocalSurfaceOrientation(
		SurfacePtr surface,
		const DataVector<int> & f_list,
		const DataVector<int> & pref,
		const DataVector<DPoint3d>& lpoints, 
		DataVector<DPoint2d>& lparams)
{
	MeshFace* fface = getFaceAt(f_list[0]);
	//const DVector3d fnormal = fface->getBaseNormal();
	const DVector3d fnormal = fface->getNormalVector();
	DPoint2d fpt_2d;
	int fpct = fface->getPointCount();
	double fratio = 1.0 / fpct;
	for(int i = 0; i < fpct; i++)
		fpt_2d.add( lparams[ pref[ fface->getPoint(i)->getIndex() ] ], fratio );
	const DVector3d snormal = surface->getNormalVector( fpt_2d );
	double sp = fnormal.scalarProduct(snormal);
	if(sp < 0.0 ){
		LOG4CPLUS_INFO(MeshLog::logger_console, "Fixing orientation in MC3d...");
		surface->invertOrientation();
		// have to recalculate params
		for(int i = 0; i < lparams.countInt(); i++) 
			lparams[i] = surface->getParameters( lpoints[i] );
	}
}

// clean and pack p_list, f_list ...
void MeshContainer3dSurface::cleanAndPackPFLists( int l_outside,
		DataVector<int> & p_list, DataVector<int> & p_layers, 
		DataVector<int> & f_list, DataVector<int> & f_layers,
		DataVector<DPoint3d>& lpoints, DataVector<DPoint2d>& lparams, 
		DataVector<double>& laquality)
{
	int j = -1;
	for(int i = 0; i < p_list.countInt(); i++) {
		if( p_layers[ p_list[i] ] < l_outside ){
			j++;
			if( j != i )
				p_list[j] = p_list[i];
		}
	}
	j++;
	p_list.leaveOnly(j);
	lpoints.leaveOnly(j);
	lparams.leaveOnly(j);
	laquality.leaveOnly(j);

	j = -1;
	for(int i = 0; i < f_list.countInt(); i++) {
		if( f_layers[ f_list[i] ] < l_outside ){
			j++;
			if( j != i ) f_list[j] = f_list[i];
		}
	}
	j++;
	if( j < f_list.countInt() )
		f_list.leaveOnly(j);
}

/// Tries to move points to fit theirs ascribed local surfaces
//int MeshContainer3dSurface::moveInnerPointsToLocalShape(Metric3dContext& mc, 
//			TagExtended::TagType tag_type, int tag_value, int forbid_tag_value)
//{
//	if(m_local_surfaces.empty() && m_local_curves.empty()) return 0;
//
//	int moved_count = 0;
//	int pct = getPointsCount();
//	for(int i = 0; i < pct; i++){
//		MeshPoint3d* point = getPointAt(i);
//		if(point->isBorder()) continue;
//		if(tag_type != TagExtended::TAG_NONE && !point->hasAnyIntFlags(tag_type, tag_value)) 
//			continue;
//		if(forbid_tag_value > 0 && point->hasAnyIntFlags(tag_type, forbid_tag_value)) 
//			continue;
//		if(point->moveToLocalShape(mc)) moved_count++;
//	}
//	return moved_count;
//}

/// Tries to move points to fit theirs ascribed local surfaces
//int MeshContainer3dSurface::moveBorderPointsToLocalShape(Metric3dContext& mc, 
//			TagExtended::TagType tag_type, int tag_value, int forbid_tag_value)
//{
//	if(m_local_surfaces.empty() && m_local_curves.empty()) return 0;
//
//	int moved_count = 0;
//	int pct = getPointsCount();
//	for(int i = 0; i < pct; i++){
//		MeshPoint3d* point = getPointAt(i);
//		if(!point->isBorder()) continue;
//		if(tag_type != TagExtended::TAG_NONE && !point->checkIntTag(tag_type, tag_value)) 
//			continue;
//		if(forbid_tag_value > 0 && point->hasAnyIntFlags(tag_type, forbid_tag_value)) 
//			continue;
//		if(point->moveToLocalShape(mc)) moved_count++;
//	}
//	return moved_count;
//}

/// Checks validity of local surface params
bool MeshContainer3dSurface::checkLocalSurfaceParams() const
{
	int pct = getPointsCount();
	int invalid = 0;
	for(int i = 0; i < pct; i++ ) {
		bool valid = getPointAt(i)->checkLocalSurfaceParams();
		if(!valid) invalid++;
	}
	return invalid == 0;
}

#define COUNT_ALL_INVALID

#ifdef COUNT_ALL_INVALID
#define IS_INVALID(x) x++
#else
#define IS_INVALID(x) return false
#endif

bool MeshContainer3dSurface::isValid(bool allow_degenerate_faces) const
{

#ifdef COUNT_ALL_INVALID
	unsigned int counter_multi_edges = 0;
	unsigned int counter_low_sp = 0;
	unsigned int counter_bad_shape = 0;
#endif // COUNT_ALL_INVALID

	// check adjacent faces for non-boundary edges
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(edge->isBorder()) continue;
		if(edge->getFaceCount() != 2){
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "MC3dS::isValid() => false, edge->getFaceCount() = " << edge->getFaceCount());
			IS_INVALID(counter_multi_edges);
		}

		MeshFace* f0 = edge->getFaceAt(0);
		MeshFace* f1 = edge->getFaceAt(1);

		DVector3d dn0, dn1;
		if( ! f0->checkAndGetNormalVector(dn0) ) continue;
		if( ! f1->checkAndGetNormalVector(dn1) ) continue;

		double sp = dn0.scalarProduct(dn1);
		if(sp < -0.9) {
			//MeshViewSet* set = new MeshViewSet;
			//set->addFaceWithEdgesAndPoints(f0);
			//set->addFaceWithEdgesAndPoints(f1);
			//SHOW_MESH("MC3dS::isValid() => false, sp < -0.9", set);
			//SHOW_MESH("MC3dS::isValid() => false, sp < -0.9", getDebugViewSetTopological( edge->getMeshPoint(0), 2));
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "MC3dS::isValid() => false, sp = " << sp);
			IS_INVALID(counter_low_sp);
		}
	}

	// check shape of faces
	int fct = getFacesCount();
	for(int i = 0; i < fct; i++){
		double q = getFaceAt(i)->getShapeQuality();
		bool invalid = allow_degenerate_faces ? (q < 0.0) : (q <= 0.0);
		if(invalid){
//			SHOW_MESH("MC3dS::isValid() => false, shape quality < 0", getDebugViewSetTopological( getFaceAt(i)->getPoint(0), 2));
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "MC3dS::isValid() => false, shape quality = " << q);
			IS_INVALID(counter_bad_shape);
		}
	}

#ifdef COUNT_ALL_INVALID
	if( counter_multi_edges + counter_low_sp + counter_bad_shape > 0 ) {
		LOG4CPLUS_INFO(MeshLog::logger_mesh, 
			"==== MESH INVALID ===="
			<< "=> MULTI_EDGES: " << counter_multi_edges << endl
			<< "=> LOW_SP:      " << counter_low_sp << endl
			<< "=> BAD_SHAPE:   " << counter_bad_shape << endl
			<< "======================");
		return false;
	}else{
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "==== MESH VALID ! ====");
	}
#endif // COUNT_ALL_INVALID

	// no errors found, then return ok
	return true;
}

/// Check average length of inner edges (with metric)
double MeshContainer3dSurface::checkInnerEdgesLength(Metric3dContext& mc, TagExtended::TagType tag_type, int tag_value)
{
	DataStatistics stat;
	DataVector<std::shared_ptr<DataStatistics>> sub_stats;
	DataHashTableKeyValue<int, int> sub_hash(1000, -1);
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(edge->isBorder() || edge->getFaceCount() != 2) continue;

		MeshFace* face0 = edge->getFaceAt(0);
		MeshFace* face1 = edge->getFaceAt(1);
		assert(face0 && face1); // since count == 2

		if(tag_type != TagExtended::TAG_NONE){
			if(!face0->checkIntTag(tag_type, tag_value)) continue;
			if(!face1->checkIntTag(tag_type, tag_value)) continue;
		}

		MeshPoint3d* p0 = edge->getMeshPoint(0);
		MeshPoint3d* p1 = edge->getMeshPoint(1);
		mc.countMetricAtPoints(p0, p1);
		double mlen = edge->getLength(mc);
		stat.add(mlen);

		int sub_domain_id = face0->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );
		if( sub_domain_id >= 0 ){
			int sub_index = sub_hash.getValue( sub_domain_id, -1 );
			if( sub_index < 0 ){
				sub_index = sub_stats.add( 
					std::make_shared<DataStatistics>() );
				sub_hash.insert( sub_domain_id, sub_index );
			}
			sub_stats[sub_index]->add(mlen);
		}
	}

	DataVector<int> keys;
	sub_hash.getKeys( keys );
	for(int i = 0; i < keys.countInt(); i++){
		int sub_index = sub_hash.getValue( keys[i], -1 );
		assert(sub_index >= 0);
		ostringstream log_info;
		log_info << " sub-domain #" << keys[i] << " - ";
		if(sub_stats[sub_index]->calculate()){
			log_info << sub_stats[sub_index]->average() 
					<< " [ " << sub_stats[sub_index]->minimum()
					<< " / " << sub_stats[sub_index]->maximum();
		}
		else
			log_info << " xxx";
		LOG4CPLUS_INFO(MeshLog::logger_mesh, log_info.str());
	}

	if(stat.calculate()){
		return stat.average();
	}else
		return -1.0;
}

/// Check average length of border edges (with metric)
double MeshContainer3dSurface::checkBorderEdgesLength(Metric3dContext& mc, 
			TagExtended::TagType tag_type, int tag_value1, int tag_value2)
{
	DataStatistics stat;
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(!edge->isBorder()) continue;
		
		if( (tag_type != TagExtended::TAG_NONE) && (edge->getFaceCount() == 2)) {
			MeshFace* face0 = edge->getFaceAt(0);
			MeshFace* face1 = edge->getFaceAt(1);
			assert(face0 && face1); // since count == 2

			if( ( !(face0->checkIntTag(tag_type, tag_value1) &&
					face1->checkIntTag(tag_type, tag_value2)) )
					&&
				( !(face1->checkIntTag(tag_type, tag_value1) &&
					face0->checkIntTag(tag_type, tag_value2)) ) ) continue;
		}

		MeshPoint3d* p0 = edge->getMeshPoint(0);
		MeshPoint3d* p1 = edge->getMeshPoint(1);
		mc.countMetricAtPoints(p0, p1);
		stat.add(edge->getLength(mc));
	}

	if(stat.calculate()){
		return stat.average();
	}else
		return -1.0;
}

/// Gahter first layer of faces around the given point
void MeshContainer3dSurface::gatherFirstLayerTopologicalForPoint(
			const MeshPoint3d* point, DataVector<int> & f_layers_count,
			DataVector<int> & p_list, DataVector<int> & p_layers, 
			DataVector<int> & f_list, DataVector<int> & f_layers,
			TagExtended::TagType tag_type, int tag_value) const
{
	if(point){
		DataVector<MeshFace*> faces;
		point->adjacentFaces(faces); 

		// faces for 0-layer
		if(point->isBorder()){
			int si = 0;
			if(tag_type != TagExtended::TAG_NONE){
				while(si < faces.countInt() && !faces[si]->checkIntTag(tag_type, tag_value)) si++;

			}
			if( si >= faces.countInt()) return; // no available faces ...
			DataVector<MeshFace*> xfaces;
			xfaces.add(faces[si]);
			for(int i = 0; i < xfaces.countInt(); i++){
				const MeshFace* face = xfaces[i];
				int efct = face->getEdgeCount();
				for(int j = 0; j < efct; j++){
					MeshEdge3d* edge = face->getEdge(j);
					if( edge->isBorder() || !edge->incidentTo(point) ) continue;
					xfaces.addIfNew(edge->getOtherFace(face));
				}
			}
			gatherFirstLayerTopological(xfaces, f_layers_count, p_list, p_layers, f_list, f_layers, tag_type, tag_value);
		}else
			gatherFirstLayerTopological(faces, f_layers_count, p_list, p_layers, f_list, f_layers, tag_type, tag_value);
	}else{ // no point -> gather all 
		gatherFirstLayerTopological( DataVector<MeshFace*>() , f_layers_count, p_list, p_layers, f_list, f_layers, tag_type, tag_value);
	}
}

/// Gahter first layer of faces for the given face (actually it's only this face and its vertices)
void MeshContainer3dSurface::gatherFirstLayerTopological(
			const DataVector<MeshFace*> &faces, 
			DataVector<int> & f_layers_count,
			DataVector<int> & p_list, DataVector<int> & p_layers, 
			DataVector<int> & f_list, DataVector<int> & f_layers,
			TagExtended::TagType tag_type, int tag_value) const
{
	f_layers_count.clear();
	f_layers_count.add(0);

	int fct = faces.countInt();
	if(fct == 0) fct = getFacesCount();	// given faces or all

	for(int i = 0; i < fct; i++){
		const MeshFace* f = faces.notEmpty() ? faces[i] : getFaceAt(i);
		if( tag_type == TagExtended::TAG_NONE || f->checkIntTag(tag_type, tag_value) ) {
			int fi = f->getIndex();
			f_list.add(fi);
			f_layers[fi] = 0;
			int fpct = f->getPointCount();
			for(int j = 0; j < fpct; j++){
				MeshPoint3d* point = f->getPoint(j);
				int pi = point->getIndex();
				if(p_layers[pi] < 0){
					p_layers[pi] = 0;
					p_list.add(pi);
				}
			}
		}
	}

	if( faces.empty() ) f_layers_count.add( f_list.countInt() );
}

/// Gather set of nodes through topological neighbourhood, with face-face adjacency, no border-crossing
bool MeshContainer3dSurface::gatherLayeredVerticesTopological(DataVector<int> & f_layers_count, 
				DataVector<int> & p_list, DataVector<int> & p_layers, 
				DataVector<int> & f_list, DataVector<int> & f_layers,
				int layers, TagExtended::TagType tag_type, int tag_value) const
{
	assert( !f_list.empty() ); // has to start from something ...
	int l_last = f_layers[f_list.last()]; // layer number for the last face in the list
	int fk = f_layers_count[l_last];

//	bool show_case = (m_local_surfaces.countInt() == 1323);
	bool show_case = false;

	for(int l = l_last; (l < layers) && (fk < f_list.countInt()); l++) {

		if(show_case){
			MeshViewSet* set = new MeshViewSet();
			for(int i = 0; i < f_list.countInt(); i++){
				set->addFaceWithEdges( getFaceAt( f_list[i]), f_layers[ f_list[i] ] );
			}
			for(int i = 0; i < p_list.countInt(); i++){
				set->addPoint( getPointAt( p_list[i] )->getCoordinates(), p_layers[ p_list[i] ], p_layers[ p_list[i]] );
			}
			SHOW_MESH("Next layer: gatherLayeredVerticesTopological", set);
		}

		for(int k = fk; k < f_list.countInt(); k++){
			int fi = f_list[k];
			if(f_layers[fi] != l) continue;
			else if(k != fk) f_list.switchData(k, fk);
			fk++;

			MeshFace* face = getFaceAt(fi); // face to check neighbours
			int ect = face->getEdgeCount();

			for(int i = 0; i < ect; i++){
				MeshEdge3d* edge = face->getEdge(i);
				if(edge->isBorder()) continue;
				int efct = edge->getFaceCount();
				for(int j = 0; j < efct; j++){
					MeshFace* other_face = edge->getFaceAt(j);
					int ofi = other_face->getIndex();
					if(f_layers[ofi] >= 0) continue; // already checked ...
					if( tag_type != TagExtended::TAG_NONE && ! other_face->checkIntTag( tag_type, tag_value ))
						continue; // wrong tag_type

					// calculate layer for face ( min_point_layer + 1)
					int opct = other_face->getPointCount();
					int min_point_layer = l+2;
					int min_cpoint_layer = l+2;
					for(int pi = 0; pi < opct; pi++){
						MeshPoint3d* point = other_face->getPoint(pi);
						int pl = p_layers[point->getIndex()];
						if(pl >= 0){
							if(point->isBorder(TagBorder::CORNER)){ // ... avoid gathering layers around single corner points
								if(pl < min_cpoint_layer) min_cpoint_layer = pl;
							}else{
								if(pl < min_point_layer) min_point_layer = pl;
							}
						}
					}
					if(min_point_layer >= l+2 && min_cpoint_layer < l+2)
					{
						min_point_layer = l;
						//min_point_layer = min_cpoint_layer;
					}

					assert(min_point_layer < l+2);

					f_layers[ofi] = min_point_layer+1;
					f_list.add(ofi);

					// update not-assigned points (or too-high assigned - for polymesh)
					for(int pi = 0; pi < opct; pi++){
						int opi = other_face->getPoint(pi)->getIndex();
						if( p_layers[opi] < 0) {
							p_layers[opi] = min_point_layer+1;
							p_list.add(opi);
						}else if(p_layers[opi] > min_point_layer+1){
							p_layers[opi] = min_point_layer+1;
						}
					}
				}
			}
		}
		assert( f_layers_count.last() != fk );
		f_layers_count.add(fk);
	}

	if(show_case){
		MeshViewSet* set = new MeshViewSet();
		for(int i = 0; i < f_list.countInt(); i++){
			set->addFaceWithEdges( getFaceAt( f_list[i]), f_layers[ f_list[i] ] );
		}
		for(int i = 0; i < p_list.countInt(); i++){
			set->addPoint( getPointAt( p_list[i] )->getCoordinates(), p_layers[ p_list[i] ], p_layers[ p_list[i]] );
		}
		SHOW_MESH("Next layer: gatherLayeredVerticesTopological", set);
	}


#ifdef SHOW_GATHER_LAYERED
	if(false){
		MeshViewSet* set = new MeshViewSet(p_list.countInt(), 0, f_list.countInt(), 0);
		for(int i = 0; i < f_list.countInt(); i++)
			if(f_layers[ f_list[i]] < layers)
				set->addFaceWithEdges(getFaceAt( f_list[i] ), f_layers[ f_list[i] ]);
			else
				set->addEdges(getFaceAt( f_list[i] ));
		for(int i = 0; i < p_list.countInt(); i++)
			if(p_layers[ p_list[i] ] < layers)
				set->addPoint(getPointAt( p_list[i] ), p_layers[ p_list[i] ]);
		SHOW_MESH("Layered vicinity - topological", set);
	}
#endif

	return fk < f_list.countInt();
/*
	for(int k = 0; k < layers; k++){
		int last_PId = points_id.countInt();
		for(int i = first_PId; i < last_PId; i++){
			MeshPoint3d* mp = getPointAt(points_id[i]);
			if(mp->isBorder()) continue; // don't cross border
			int rank = mp->getRank();
			for(int j = 0; j < rank; j++){
				MeshEdge3d* edge = mp->getEdge(j);
				MeshPoint3d* other_mp = edge->getOtherPoint(mp);
				int pid = other_mp->getIndex();
				if(included_points[pid]) continue;
//				if(other_mp->hasLocalSurface()) continue;
				if(tag_type != TagExtended::TAG_NONE && !other_mp->hasAnyIntFlags(tag_type, tag_value)) 
					continue;
				// ok, include new point in the neighborhood
				points_id.add(pid);
				included_points[pid] = true;
			}
		}
		first_PId = last_PId;
	}
	return first_PId;
*/
}

/// Gather set of nodes through neighbourhood (point-point adjacency) with surface fitting
bool MeshContainer3dSurface::gatherLayeredVerticesForSurface(Metric3dContext& mc, 
		SurfaceConstPtr surface, DataVector<int> & f_layers_count, 
		DataVector<int> & p_list, DataVector<int> & p_layers, 
		DataVector<int> & f_list, DataVector<int> & f_layers,
		double tolerance,  TagExtended::TagType tag_type, int tag_value) const
{
	assert( ! f_list.empty() );
	const int l_last = f_layers[f_list.last()]; // layer number for the last face in the list
	int fk = f_layers_count[l_last];
	const int l_outside = l_last+1;

	double tol2 = tolerance*tolerance;

	int pct = getPointsCount();

	// ... recheck last (extra) layer points
	int p_found = 0;
	for(int i = 0; i < p_list.countInt(); i++){
		int pid = p_list[i];
		if( p_layers[pid] != l_last ) continue;
		MeshPoint3d* mp = getPointAt(pid);
		const DPoint3d& other_pt = mp->getCoordinates();
		const DPoint2d& other_param = mp->getLocalSurfaceParam(surface);
		const DPoint3d surf_pt = surface->getPoint(other_param);
		mc.countMetricAtPoint(other_pt);
		double dist2 = mc.transformRStoMS(surf_pt - other_pt).length2();
		if(dist2 > tol2) p_layers[pid] = l_outside;
		else p_found++;
	}
	if(p_found == 0){ // nothing new, revert the changes, and quit
		for(int i = 0; i < p_list.countInt(); i++){
			int pid = p_list[i];
			if( p_layers[pid] == l_outside )
				p_layers[pid] = l_last;
		}
		return fk < f_list.countInt();
	}

	// ... recheck last (extra) layer triangles
	int faces_in = f_list.countInt() - fk;
	for(int k = fk; k < f_list.countInt(); k++){
		int fi = f_list[k];
		assert( f_layers[fi] == l_last );
		MeshFace* face = getFaceAt(fi);
		int fpct = face->getPointCount();
		for(int i = 0; i < fpct; i++){
			if(p_layers[face->getPoint(i)->getIndex()] == l_outside){
				f_layers[fi] = l_outside;
				faces_in--;
				break;
			}
		}
	}

	if(faces_in == 0) { // actually, no new faces (desPIte some new points) -> revert and quit
		for(int i = 0; i < p_list.countInt(); i++){
			int pid = p_list[i];
			if( p_layers[pid] == l_outside )
				p_layers[pid] = l_last;
		}
		for(int k = fk; k < f_list.countInt(); k++){
			assert( f_layers[f_list[k]] == l_outside );
			f_layers[f_list[k]] = l_last;
		}
		return fk < f_list.countInt();
	}
			
	for(int k = fk; k < f_list.countInt(); k++){
		int fi = f_list[k];
		assert( f_layers[fi] == l_last || f_layers[fi] == l_outside);
		if(f_layers[fi] == l_outside) continue;
		else if(k != fk) f_list.switchData(k, fk);
		fk++;

		MeshFace* face = getFaceAt(fi); // face to check neighbours
		int ect = face->getEdgeCount();

		for(int i = 0; i < ect; i++){
			MeshEdge3d* edge = face->getEdge(i);
			if(edge->isBorder()) continue;
			int efct = edge->getFaceCount();
			for(int j = 0; j < efct; j++){
				MeshFace* other_face = edge->getFaceAt(j);
				int ofi = other_face->getIndex();
				if(f_layers[ofi] >= 0) continue; // already checked ...
				if( tag_type != TagExtended::TAG_NONE && ! other_face->checkIntTag( tag_type, tag_value ))
					continue; // wrong tag_type

//				if(other_face->hasLocalSurface())
//					continue;	// already set for another surface ...

				int opct = other_face->getPointCount();
				bool all_other_points_in = true;
				for(int pi = 0; pi < opct; pi++){
					MeshPoint3d* other_point = other_face->getPoint(pi);
					int pid = other_point->getIndex();
					if(p_layers[pid] < 0){
						const DPoint3d& other_pt = other_point->getCoordinates();
						const DPoint2d& other_param = other_point->getLocalSurfaceParam(surface);
						const DPoint3d surf_pt = surface->getPoint(other_param);
						mc.countMetricAtPoint(other_pt);
						double dist2 = mc.transformRStoMS(surf_pt - other_pt).length2();
						if(dist2 > tol2){
							p_layers[pid] = l_outside;
							all_other_points_in = false;
						} else 
							p_layers[pid] = l_last;
						p_list.add(pid);
					}else if(p_layers[pid] == l_outside) {
						all_other_points_in = false;
					}
				}

				f_layers[ofi] = all_other_points_in ? l_last : l_outside;
				f_list.add(ofi);

/* Previous version -> for triangle only

				MeshPoint3d* other_point = other_face->getOtherPoint( 
					edge->getMeshPoint(0), edge->getMeshPoint(1));
				int pid = other_point->getIndex();
				if(p_layers[pid] < 0) { 
					const DPoint3d& other_pt = other_point->getCoordinates();
					const DPoint2d other_param = surface->getParameters(other_pt);
					const DPoint3d surf_pt = surface->getPoint(other_param);
					mc.countMetricAtPoint(other_pt);
					double dist2 = mc.transformRStoMS(surf_pt - other_pt).length2();
					if(dist2 > tol2) p_layers[pid] = l_outside;
					else p_layers[pid] = l_last;
					p_list.add(pid);
				}

				f_layers[ofi] = p_layers[pid];
				f_list.add(ofi);
*/
			}
		}
	}

	assert( f_layers_count.last() != fk );
	f_layers_count.add(fk);

#ifdef SHOW_GATHER_LAYERED
	if(false){
		MeshViewSet* set = new MeshViewSet(p_list.countInt(), 0, f_list.countInt(), 0);
		for(int i = 0; i < f_list.countInt(); i++)
			if(f_layers[ f_list[i]] < l_outside)
				set->addFaceWithEdges(getFaceAt( f_list[i] ), f_layers[ f_list[i] ]);
			else
				set->addEdges(getFaceAt( f_list[i] ));
		for(int i = 0; i < p_list.countInt(); i++)
			if(p_layers[ p_list[i] ] < l_outside)
				set->addPoint(getPointAt( p_list[i] ), p_layers[ p_list[i] ]);
		SHOW_MESH("Layered vicinity - surface", set);
	}
#endif

	return fk < f_list.countInt();
}

/// Gather set of nodes through topological neighbourhood, with point-point adjacency
bool MeshContainer3dSurface::gatherLayeredVerticesGeometrical(Metric3dContext& mc, 
		const DPoint3d& dpt, DataVector<int> & f_layers_count, 
		DataVector<int> & p_list, DataVector<int> & p_layers, 
		DataVector<int> & f_list, DataVector<int> & f_layers,
		double radius, TagExtended::TagType tag_type, int tag_value) const
{
	assert( ! f_list.empty() );
	const int l_last = f_layers[f_list.last()]; // layer number for the last face in the list
	int fk = f_layers_count[l_last];
	const int l_outside = l_last+1;

//	bool show_case = (m_local_surfaces.countInt() == 1323);
	bool show_case = false;

	if(show_case){
		MeshViewSet* set = new MeshViewSet();
		for(int i = 0; i < f_list.countInt(); i++){
			set->addFaceWithEdges( getFaceAt( f_list[i]), f_layers[ f_list[i] ] );
		}
		for(int i = 0; i < p_list.countInt(); i++){
			set->addPoint( getPointAt( p_list[i] )->getCoordinates(), p_layers[ p_list[i] ], p_layers[ p_list[i]] );
//			set->addPoint( getPointAt( p_list[i] ), p_layers[ p_list[i] ] );
		}
		SHOW_MESH("Start: gatherLayeredVerticesGeometrical", set);
	}

	// ... prepare local metric
	mc.countMetricAtPoint(dpt);
	const DMPoint3d dmpoint = mc.transformRStoMS(dpt);
	const double R2 = radius * radius;

	// ... recheck last (extra) layer points
	int p_found = 0;
	for(int i = 0; i < p_list.countInt(); i++){
		int pid = p_list[i];
		if( p_layers[pid] != l_last ) continue;
		MeshPoint3d* mp = getPointAt(pid);
		double dist2 = (mp->getMetricCoordinates(mc) - dmpoint).length2();
		if(dist2 > R2) p_layers[pid] = l_outside;
		else p_found++;
	}
	if(p_found == 0){ // nothing new, revert the changes, and quit
		for(int i = 0; i < p_list.countInt(); i++){
			int pid = p_list[i];
			if( p_layers[pid] == l_outside )
				p_layers[pid] = l_last;
		}
		return fk < f_list.countInt();
	}

	// ... recheck last (extra) layer triangles
	int faces_in = f_list.countInt() - fk;
	for(int k = fk; k < f_list.countInt(); k++){
		int fi = f_list[k];
		assert( f_layers[fi] == l_last );
		MeshFace* face = getFaceAt(fi);
		if( tag_type != TagExtended::TAG_NONE && ! face->checkIntTag( tag_type, tag_value ))
			continue; // wrong tag_type
		int fpct = face->getPointCount();
		for(int i = 0; i < fpct; i++){
			if(p_layers[face->getPoint(i)->getIndex()] == l_outside){
				f_layers[fi] = l_outside;
				faces_in--;
				break;
			}
		}
	}

	if(faces_in == 0) { // actually, no new faces (desPIte some new points) -> revert and quit
		for(int i = 0; i < p_list.countInt(); i++){
			int pid = p_list[i];
			if( p_layers[pid] == l_outside )
				p_layers[pid] = l_last;
		}
		for(int k = fk; k < f_list.countInt(); k++){
			assert( f_layers[f_list[k]] == l_outside );
			f_layers[f_list[k]] = l_last;
		}
		return fk < f_list.countInt();
	}

	if(show_case){
		MeshViewSet* set = new MeshViewSet();
		for(int i = 0; i < f_list.countInt(); i++){
			set->addFaceWithEdges( getFaceAt( f_list[i]), f_layers[ f_list[i] ] );
		}
		for(int i = 0; i < p_list.countInt(); i++){
			set->addPoint( getPointAt( p_list[i] )->getCoordinates(), p_layers[ p_list[i] ], p_layers[ p_list[i]] );
//			set->addPoint( getPointAt( p_list[i] ), p_layers[ p_list[i] ] );
		}
		SHOW_MESH("After recheck: gatherLayeredVerticesGeometrical", set);
	}

	for(int k = fk; k < f_list.countInt(); k++){
		int fi = f_list[k];
		assert( f_layers[fi] == l_last || f_layers[fi] == l_outside);
		if(f_layers[fi] == l_outside) continue;
		else if(k != fk) f_list.switchData(k, fk);
		fk++;

		MeshFace* face = getFaceAt(fi); // face to check neighbours
		int ect = face->getEdgeCount();

		for(int i = 0; i < ect; i++){
			MeshEdge3d* edge = face->getEdge(i);
			if(edge->isBorder()) continue;
			int efct = edge->getFaceCount();
			for(int j = 0; j < efct; j++){
				MeshFace* other_face = edge->getFaceAt(j);
				int ofi = other_face->getIndex();
				if(f_layers[ofi] >= 0) continue; // already checked ...
				if( tag_type != TagExtended::TAG_NONE && ! other_face->checkIntTag( tag_type, tag_value ))
					continue; // wrong tag_type

				int opct = other_face->getPointCount();
				bool all_other_points_in = true;
				for(int pi = 0; pi < opct; pi++){
					MeshPoint3d* other_point = other_face->getPoint(pi);
					int pid = other_point->getIndex();
					if(p_layers[pid] < 0){
						double dist2 = (other_point->getMetricCoordinates(mc) - dmpoint).length2();
						if(dist2 > R2){
							p_layers[pid] = l_outside;
							all_other_points_in = false;
						} else 
							p_layers[pid] = l_last;
						p_list.add(pid);
					}else if(p_layers[pid] == l_outside) {
						all_other_points_in = false;
					}
				}

				f_layers[ofi] = all_other_points_in ? l_last : l_outside;
				f_list.add(ofi);
			}
		}
	}

	assert(f_layers_count.last() != fk );
	f_layers_count.add(fk);

	if(show_case){
		MeshViewSet* set = new MeshViewSet();
		for(int i = 0; i < f_list.countInt(); i++){
			set->addFaceWithEdges( getFaceAt( f_list[i]), f_layers[ f_list[i] ] );
		}
		for(int i = 0; i < p_list.countInt(); i++){
			set->addPoint( getPointAt( p_list[i] )->getCoordinates(), p_layers[ p_list[i] ], p_layers[ p_list[i]] );
		}
		SHOW_MESH("Start: gatherLayeredVerticesGeometrical", set);
	}

#ifdef SHOW_GATHER_LAYERED
	if(false){
		MeshViewSet* set = new MeshViewSet(p_list.countInt(), 0, f_list.countInt(), 0);
		for(int i = 0; i < f_list.countInt(); i++)
			if(f_layers[ f_list[i]] < l_outside)
				set->addFaceWithEdges(getFaceAt( f_list[i] ), f_layers[ f_list[i] ]);
			else
				set->addEdges(getFaceAt( f_list[i] ));
		for(int i = 0; i < p_list.countInt(); i++)
			if(p_layers[ p_list[i] ] < l_outside)
				set->addPoint(getPointAt( p_list[i] ), p_layers[ p_list[i] ]);
		SHOW_MESH("Layered vicinity - geometrical", set);
	}
#endif

	return fk < f_list.countInt();
/*
	// ... search
	for(int i = first_PId; i < points_id.countInt(); i++){
		MeshPoint3d* mp = getPointAt(points_id[i]);
		if(mp->isBorder()) continue; // don't cross border
		int rank = mp->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = mp->getEdge(j);
			MeshPoint3d* other_mp = edge->getOtherPoint(mp);
			int pid = other_mp->getIndex();
			if(included_points[pid]) continue;
//			if(other_mp->hasLocalSurface()) continue;
			if(tag_type != TagExtended::TAG_NONE && !other_mp->hasAnyIntFlags(tag_type, tag_value)) 
				continue;
			double dist2 = (other_mp->getMetricCoordinates(mc) - dmpoint).length2();
			if(dist2 > R2) continue;
			// ok, include new point in the neighborhood
			points_id.add(pid);
			included_points[pid] = true;
		}
	}
*/
}

/// Gather set of nodes through neighbourhood, via normals, with face-face adjacency
double MeshContainer3dSurface::gatherLayeredVerticesViaNormals(
		SurfaceConstPtr base_surface, double sp_min,
		DataVector<int> & f_layers_count, 
		DataVector<int> & p_list, DataVector<int> & p_layers, 
		DataVector<int> & f_list, DataVector<int> & f_layers,
		TagExtended::TagType tag_type, int tag_value) const
{
	assert( !f_list.empty() ); // has to start from something ...
	int l_last = f_layers[f_list.last()]; // layer number for the last face in the list
	int fk = f_layers_count[l_last];

//	bool show_case = (m_local_surfaces.countInt() == 1323);
	bool show_case = false;

	double show_normals_len = 0.0;
	if( show_case ) {
		DBox bbox;
		for(int i = 0; i < p_list.countInt(); i++){
			bbox.addPoint( getPointAt( p_list[i] )->getCoordinates() );
		}
		show_normals_len = 0.2 * bbox.getDiameter();
	}

	double actual_min_sp = 1.0;

	for(int l = l_last; fk < f_list.countInt(); l++) {

		if(show_case){
			MeshViewSet* set = new MeshViewSet();
			for(int i = 0; i < f_list.countInt(); i++){
				MeshFace* face = getFaceAt( f_list[i]);
				set->addFaceWithEdges( face, f_layers[ f_list[i] ] );
				DPoint3d fpmid = face->getMiddlePoint();
				set->addEdge( fpmid, 
					fpmid + face->getBaseNormal().normalized( show_normals_len ), 2 );
			}
			for(int i = 0; i < p_list.countInt(); i++){
				set->addPoint( getPointAt( p_list[i] )->getCoordinates(), 
					p_layers[ p_list[i] ], p_layers[ p_list[i]] );
			}
			SHOW_MESH("Next layer: gatherLayeredVerticesViaNormals", set);
		}

		for(int k = fk; k < f_list.countInt(); k++){
			int fi = f_list[k];
			if(f_layers[fi] != l) continue;
			else if(k != fk) f_list.switchData(k, fk);
			fk++;

			MeshFace* face = getFaceAt(fi); // face to check neighbours
			int ect = face->getEdgeCount();

			for(int i = 0; i < ect; i++){
				MeshEdge3d* edge = face->getEdge(i);
				if(edge->isBorder()) continue;
				int efct = edge->getFaceCount();
				for(int j = 0; j < efct; j++){
					MeshFace* other_face = edge->getFaceAt(j);
					int ofi = other_face->getIndex();
					if(f_layers[ofi] >= 0) continue; // already checked ...
					if( tag_type != TagExtended::TAG_NONE && ! other_face->checkIntTag( tag_type, tag_value ))
						continue; // wrong tag_type

					// check normal
					double sp = base_surface->getNormalVectorForPoint3d( other_face->getMiddlePoint() 
								).scalarProduct( other_face->getBaseNormal() );
					if( sp < sp_min ) continue;

					// check param
					MeshPoint3d* ep0 = edge->getMeshPoint(0);
					MeshPoint3d* ep1 = edge->getMeshPoint(1);
					bool all_params_in_range = true;
					int ofpct = other_face->getPointCount();
					for(int k = 0; all_params_in_range && (k < ofpct); k++){
						MeshPoint3d* ep = other_face->getPoint(k);
						if( ep == ep0 || ep == ep1 ) continue;
						all_params_in_range = base_surface->withinParamRange( 
							base_surface->getParameters( 
								ep->getCoordinates() ) );
					}
					if( !all_params_in_range ) continue;

					if( sp < actual_min_sp ) actual_min_sp = sp;

					f_layers[ofi] = l+1;
					f_list.add(ofi);

					// update not-assigned points (or too-high assigned - for polymesh)
					int opct = other_face->getPointCount();
					for(int pi = 0; pi < opct; pi++){
						int opi = other_face->getPoint(pi)->getIndex();
						if( p_layers[opi] < 0) {
							p_layers[opi] = l+1;
							p_list.add(opi);
						}else if(p_layers[opi] > l+1){
							p_layers[opi] = l+1;
						}
					}
				}
			}
		}
		assert( f_layers_count.last() != fk );
		f_layers_count.add(fk);
	}

	//show_case = true;
	show_case = ( base_surface->getType() != ElementType::SURFACE_PLANE );

	if(show_case){
		MeshViewSet* set = new MeshViewSet();
		for(int i = 0; i < f_list.countInt(); i++){
			MeshFace* face = getFaceAt( f_list[i]);
			set->addFaceWithEdges( face, f_layers[ f_list[i] ] );
			DPoint3d fpmid = face->getMiddlePoint();
			set->addEdge( fpmid, 
				fpmid + face->getBaseNormal().normalized( show_normals_len ), 2 );
		}
		if(true){ // show normals for base surface and arbitrary face
			MeshFace* face = getFaceAt( f_list[0]);
			DPoint3d fpmid = face->getMiddlePoint();
			if( show_normals_len == 0.0 ){
				DBox bbox;
				for(int i = 0; i < p_list.countInt(); i++){
					bbox.addPoint( getPointAt( p_list[i] )->getCoordinates() );
				}
				show_normals_len = 0.2 * bbox.getDiameter();
			}
			set->addEdge( fpmid, 
				fpmid + face->getBaseNormal().normalized( show_normals_len ), 2 );
			set->addEdge( fpmid,
				fpmid + base_surface->getNormalVectorForPoint3d(fpmid).normalized( show_normals_len), 1 );

			// center
			DRect param_rect;
			DataVector<DPoint2d> params(p_list.countInt());
			for(int i = 0; i < p_list.countInt(); i++){
				DPoint2d param = base_surface->getParameters( getPointAt( p_list[i] )->getCoordinates() );
				param_rect.addPoint( param );
				params.add( param );
			}
			base_surface->createViewSetForPoints( set, params );
			set->addInfo("x range", to_string(param_rect.x0) + " | " + to_string(param_rect.x1) );
			set->addInfo("y range", to_string(param_rect.y0) + " | " + to_string(param_rect.y1) );
			set->addLabel( base_surface->getPoint( DPoint2d( 0.5, (param_rect.y0+param_rect.y1)*0.5) ), "midx-0.5" );
		}
		for(int i = 0; i < p_list.countInt(); i++){
			set->addPoint( getPointAt( p_list[i] )->getCoordinates(), p_layers[ p_list[i] ], p_layers[ p_list[i]] );
		}
		SHOW_MESH("All layers: gatherLayeredVerticesViaNormals", set);
	}

	return actual_min_sp;
}

/// Gather set of nodes through topological neighbourhood, with face-face adjacency
bool MeshContainer3dSurface::gatherLayeredVerticesTopological(
		MeshFace* start_face,	DataVector<DPoint3d> & points, 
		int layers, bool cross_borders, int min_points, 
		DataVector<MeshFace*> * layer_faces, 
		TagExtended::TagType tag_type, int tag_value) const
{
	assert(start_face);
	assert( tag_type == TagExtended::TAG_NONE || start_face->checkIntTag( tag_type, tag_value ) );

	int fct = getFacesCount();
	int pct = getPointsCount();

	DataVector<int> p_layers(pct, -1);
	DataVector<int> f_layers(fct, -1);

	DataVector<int> f_list(fct);
	DataVector<int> p_list(pct);

	int sfi = start_face->getIndex();
	f_list.add(sfi);
	f_layers[sfi] = 0;

	int fpct = start_face->getPointCount();
	for(int i = 0; i < fpct; i++){
		int pi = start_face->getPoint(i)->getIndex();
		p_layers[pi] = 0;
		p_list.add(pi);
	}

	int fk = 0;
	for(int l = 0; (l < layers || p_list.countInt() < min_points) && (fk < f_list.countInt()); l++) {

		for(int k = fk; k < f_list.countInt(); k++){
			int fi = f_list[k];
			if(f_layers[fi] != l) continue;
			else if(k != fk) f_list.switchData(k, fk);
			fk++;

			MeshFace* face = getFaceAt(fi); // face to check neighbours
			int ect = face->getEdgeCount();

			for(int i = 0; i < ect; i++){
				MeshEdge3d* edge = face->getEdge(i);
				if(!cross_borders && edge->isBorder()) continue;
				int efct = edge->getFaceCount();
				for(int j = 0; j < efct; j++){
					MeshFace* other_face = edge->getFaceAt(j);
					int ofi = other_face->getIndex();
					if(f_layers[ofi] >= 0) continue; // already checked ...
					if( tag_type != TagExtended::TAG_NONE && ! other_face->checkIntTag( tag_type, tag_value ))
						continue; // wrong tag_type

					// calculate layer for face ( min_point_layer + 1)
					int opct = other_face->getPointCount();
					int min_point_layer = l+2;
					for(int pi = 0; pi < opct; pi++){
						int pl = p_layers[other_face->getPoint(pi)->getIndex()];
						if( pl >= 0 && pl < min_point_layer) min_point_layer = pl;
					}
					assert(min_point_layer < l+2);

					f_layers[ofi] = min_point_layer+1;
					f_list.add(ofi);

					// update not-assigned points
					for(int pi = 0; pi < opct; pi++){
						int opi = other_face->getPoint(pi)->getIndex();
						if( p_layers[opi] >= 0) continue;
						p_layers[opi] = min_point_layer+1;
						p_list.add(opi);
					}
				}
			}
		}
	}

	if(layer_faces){
		for(int i = 0; i < f_list.countInt(); i++){
			int fid = f_list[i];
			MeshFace* face = getFaceAt(fid);
			face->setIntTag(TagExtended::TAG_LAYER_ID, f_layers[fid]);
			layer_faces->add(face);
		}
	}

	if(false){
		MeshViewSet* set = new MeshViewSet(p_list.countInt(), 0, f_list.countInt(), 0);
		for(int i = 0; i < f_list.countInt(); i++)
			if(f_layers[ f_list[i]] < layers)
				set->addFaceWithEdges(getFaceAt( f_list[i] ), f_layers[ f_list[i] ]);
		for(int i = 0; i < p_list.countInt(); i++)
			if(p_layers[ p_list[i] ] < layers)
				set->addPoint(getPointAt( p_list[i] ), p_layers[ p_list[i] ]);
		SHOW_MESH("Layered vicinity", set);
	}

	// + sort p_list

	for(int i = 0; i < p_list.countInt(); i++) {
		int pid = p_list[i];
		if( p_layers[ pid ] < layers )
			points.add( getPointAt( pid )->getCoordinates() );
	}
	for(int i = 0; i < f_list.countInt(); i++) {
		int fid = f_list[i];
		if(f_layers[ fid ] < layers)
			points.add( getFaceAt( fid )->getMiddlePoint() );
	}

	return p_list.countInt() >= min_points;
}

/// Gather chain of boundary edges forming border contour, not crossing "corner"-boundary-points
MeshPoint3d* MeshContainer3dSurface::gatherBorderContourChain(
		MeshEdge3d* start_edge, MeshPoint3d* start_point, 
		DataVector<MeshEdge3d*> & edges, 
		int layers, bool only_without_curve)
{
	assert(start_edge && start_edge->isBorder());
	MeshPoint3d* point = start_point;
	while(! point->isBorder(TagBorder::CORNER) && (layers-- != 0)) {
		MeshEdge3d* other_edge = nullptr;
		int i = 0;
		assert( point->getBorderEdgesCount() == 2 ); // since not corner-point
		do{
			other_edge = point->getEdge(i++);
		}while((other_edge == start_edge) || (!other_edge->isBorder()));
		if( only_without_curve && other_edge->hasLocalCurve() )
			return point;
		edges.add(other_edge);
		start_edge = other_edge;
		point = other_edge->getOtherPoint(point);
		if(point == start_point) return point;
	}
	return point;
}

/// Gather chain of boundary edges forming border contour, not crossing "corner"-boundary-points
MeshPoint3d* MeshContainer3dSurface::gatherBorderContourChainForSurface(
		MeshEdge3d* start_edge, MeshPoint3d* start_point, 
		SurfaceConstPtr surface, DataVector<MeshEdge3d*> & edges)
{
	assert(start_edge && start_edge->isBorder());
	MeshPoint3d* point = start_point;
	assert( (surface == nullptr) || (start_point->hasLocalSurface( surface ) ) );
	while(! point->isBorder(TagBorder::CORNER) ) {
		MeshEdge3d* other_edge = nullptr;
		int i = 0;
		assert( point->getBorderEdgesCount() == 2 ); // since not corner-point
		do{
			other_edge = point->getEdge(i++);
		}while((other_edge == start_edge) || (!other_edge->isBorder()));
		MeshPoint3d* other_point = other_edge->getOtherPoint(point);
		if( (surface != nullptr) && !other_point->hasLocalSurface( surface ) )
			return point;
		edges.add(other_edge);
		start_edge = other_edge;
		point = other_point;
		if(point == start_point) return point;
	}
	return point;
}

/// Gather chain of boundary edges forming border contour, not crossing "corner"-boundary-points
MeshPoint3d* MeshContainer3dSurface::gatherBorderContourChain(
		MeshEdge3d* start_edge, 
		DataVector<MeshEdge3d*> & chain_edges,
		int layers, bool only_without_curve)
{
	MeshPoint3d* ept0 = start_edge->getMeshPoint(0);
	MeshPoint3d* ept1 = start_edge->getMeshPoint(1);

	DataVector<MeshEdge3d*> edges0, edges1;
	MeshPoint3d* end_pt1 = gatherBorderContourChain(start_edge, ept1, edges1, layers, only_without_curve);
	if((end_pt1 == ept1) && edges1.notEmpty()){ // closed chain of edges...
		int ct1 = edges1.countInt();

		// prepare full and ordered chain of edges
		chain_edges.clear();
		chain_edges.prepare(ct1);

		for(int i = 0; i < ct1; i++)
			chain_edges.add(edges1[i]);
		assert(chain_edges.notEmpty());
		return end_pt1;
	}else{
		MeshPoint3d* end_pt0 = gatherBorderContourChain(start_edge, ept0, edges0, layers, only_without_curve);
		int ct0 = edges0.countInt();
		int ct1 = edges1.countInt();

		// prepare full and ordered chain of edges
		chain_edges.clear();
		chain_edges.prepare(ct0+ct1+1);

		for(int i = ct0-1; i >= 0; i--)
			chain_edges.add(edges0[i]);
		chain_edges.add(start_edge);
		for(int i = 0; i < ct1; i++)
			chain_edges.add(edges1[i]);

		assert(chain_edges.notEmpty());
		return end_pt0;
	}
}

/// Gather chain of boundary edges forming border contour, not crossing "corner"-boundary-points
MeshPoint3d* MeshContainer3dSurface::gatherBorderContourChainForSurface(
		MeshEdge3d* start_edge, SurfaceConstPtr surface,
		DataVector<MeshEdge3d*> & chain_edges )
{
	MeshPoint3d* ept0 = start_edge->getMeshPoint(0);
	MeshPoint3d* ept1 = start_edge->getMeshPoint(1);

	DataVector<MeshEdge3d*> edges0, edges1;
	MeshPoint3d* end_pt1 = gatherBorderContourChainForSurface(start_edge, ept1, surface, edges1);
	if((end_pt1 == ept1) && edges1.notEmpty()){ // closed chain of edges...
		int ct1 = edges1.countInt();

		// split chain into two ...
		int ct1_half = ct1/2;

		// prepare full and ordered chain of edges
		chain_edges.clear();
		chain_edges.prepare(ct1_half);

		chain_edges.add( start_edge );
		for(int i = 0; i < ct1_half; i++){
			end_pt1 = edges1[i]->getOtherPoint( end_pt1 );
			chain_edges.add(edges1[i]);
		}

		ept0->setBorderFlags( TagBorder::CORNER );
		end_pt1->setBorderFlags( TagBorder::CORNER );

		return ept0;
	}else{
		MeshPoint3d* end_pt0 = gatherBorderContourChainForSurface(start_edge, ept0, surface, edges0);
		int ct0 = edges0.countInt();
		int ct1 = edges1.countInt();

		// prepare full and ordered chain of edges
		chain_edges.clear();
		chain_edges.prepare(ct0+ct1+1);

		for(int i = ct0-1; i >= 0; i--)
			chain_edges.add(edges0[i]);
		chain_edges.add(start_edge);
		for(int i = 0; i < ct1; i++)
			chain_edges.add(edges1[i]);

		assert(chain_edges.notEmpty());
		return end_pt0;
	}
}

/// Convert all non-triangular faces to triangles
void MeshContainer3dSurface::convertPolysToTriangles()
{
	START_CLOCK("MC3dS:convertPolysToTriangles");

	DataVector<MeshFace*> split_faces(100);

	int fct = getFacesCount();
	DataVector<MeshFace*> polys(fct);

	for(int i = 0; i < fct; i++){
		MeshFace* face = getFaceAt(i);
		if(face->getType() == FACE_QUAD || face->getType() == FACE_POLY){
			polys.add(face);
		}
	}

//	assert(isValid());
	for(int i = 0; i < polys.countInt(); i++){

#ifdef T_DEBUG_
		if( i % 100 == 0 ){
			ostringstream ostr;
			ostr << "Converting polys to triangles, face " << (i+1) << " out of " << polys.countInt();
			LOG4CPLUS_INFO(MeshLog::logger_console,ostr.str());
		}
#endif // T_DEBUG_

		MeshFace* face = polys[i];
		//if(i==3562){
		//	MeshViewSet* set = new MeshViewSet();
		//	set->addFaceWithEdgesAndPoints(face);
		//	SHOW_MESH("face to split", set);
		//}
		if( face->splitToTriangles(split_faces) ) {
			//if(i==3562){
			//	MeshViewSet* set = new MeshViewSet();
			//	for(int j = 0; j < split_faces.countInt(); j++)
			//		set->addFaceWithEdgesAndPoints(split_faces[j]);
			//	SHOW_MESH("face after split", set);
			//}
			delete removeMeshFace(face->getIndex());
			for(int j = 0; j < split_faces.countInt(); j++)
				addMeshFace( split_faces[j] );
//			assert(isValid());
			split_faces.clear();
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Failed splitting poly-face to triangles...");
		}
	}

/*
	assert(!m_faces->isHeapOrder());

	// find faces to split (and remove)
	int i = 0;
	while( i < getFacesCount() ){
		MeshFace* face = getFaceAt(i);
		switch(face->getType()) {
		case FACE_TRIANGLE:
			i++; 
			break;
		case FACE_QUAD:
		case FACE_POLY:
			if( face->splitToTriangles(split_faces) ) {
				delete removeMeshFace(i);
			}else{
				i++;
				LOG4CPLUS_WARN(MeshLog::logger_console, "Failed splitting poly-face to triangles...");
			}
			break;
		default:
			assert(false);
		}
	}

	// add new faces
	for(int i = 0; i < split_faces.countInt(); i++)
		addMeshFace( split_faces[i] );
*/
	STOP_CLOCK("MC3dS:convertPolysToTriangles");
}

/// Store surface mesh to .XML file
bool MeshContainer3dSurface::storeXML(const string& fname, const string& description) const
{
	ofstream fxml(fname);
	if(!fxml){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening file for write: " << fname);
		return false;
	}
	fxml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
	fxml << "<meshdoc xmlns=\"http://www.icsr.agh.edu.pl\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">" << endl;
	fxml << "\t<header>" << endl;
	fxml << "\t\t<creator>Tomasz Jurczyk</creator>" << endl;
	fxml << "\t\t<version>" << mesh_data.version() << "</version>" << endl;
	if(description != "")
		fxml << "\t\t<description>" << description << "</description>" << endl;
	fxml << "\t</header>" << endl;
	fxml << "\t<surface-mesh>" << endl;

	// points
	int pct = getPointsCount();
	fxml << "\t\t<pointarray count=\"" << pct << "\">" << endl;
	for(int i = 0; i < pct; i++){
		const DPoint3d& pt = getPointAt(i)->getCoordinates();
		fxml << "\t\t\t" << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
	}
	fxml << "\t\t</pointarray>" << endl;

	// faces
	int fct = getFacesCount();
	fxml << "\t\t<facearray count=\"" << fct<< "\">" << endl;
	for(int i = 0; i < fct; i++){
		MeshFace* face = getFaceAt(i);
		int fpct = face->getPointCount();
		fxml << "\t\t\t" << fpct;
		for(int j = 0; j < fpct; j++){
			fxml << " " << face->getPoint(j)->getIndex();
		}
		MeshBlock* bl = face->getBlock(0);
		fxml << "\t" << (bl ? bl->getAreaID() : -1);
		bl = face->getBlock(1);
		fxml << "\t" << (bl ? bl->getAreaID() : -1);
		fxml << endl;
	}
	fxml << "\t\t</facearray>" << endl;

	// local surfaces
	int lsct = m_local_surfaces.countInt();
	if(lsct > 0){
		fxml << "\t\t<local-surfaces>" << endl;
		for(int i = 0; i < lsct; i++){
			SurfaceConstPtr surface = m_local_surfaces[i];
			DataVector<int> PIds(pct);
			for(int j = 0; j < pct; j++)
				if(getPointAt(j)->getLocalValidSurface() == surface) PIds.add(j);
			if(PIds.notEmpty()){
				fxml << "\t\t\t<local-surface>" << endl;
				surface->storeXML(fxml, "\t\t\t\t");
				surface->storeDomainXML(fxml, "\t\t\t\t");
				fxml << "\t\t\t\t<local-points>" << endl;
				fxml << "\t\t\t\t";
				for(int j = 0; j < PIds.countInt(); j++)
					fxml << " " << PIds[j];
				fxml << endl;
				fxml << "\t\t\t\t</local-points>" << endl;
				fxml << "\t\t\t</local-surface>" << endl;
			}
		}
		fxml << "\t\t</local-surfaces>" << endl;
	}

	// local curves
	int lcct = m_local_curves.countInt();
	if(lcct > 0){
		fxml << "\t\t<local-curves>" << endl;
		for(int i = 0; i < lcct; i++){
			Curve3dConstPtr curve = m_local_curves[i];
			DataVector<int> PIds(pct);
			for(int j = 0; j < pct; j++)
				if(getPointAt(j)->getLocalCurve() == curve) PIds.add(j);
			if(PIds.notEmpty()){
				fxml << "\t\t\t<local-curve>" << endl;
				curve->storeXML(fxml, "\t\t\t\t");
				fxml << "\t\t\t\t<local-points>" << endl;
				fxml << "\t\t\t\t";
				for(int j = 0; j < PIds.countInt(); j++)
					fxml << " " << PIds[j];
				fxml << endl;
				fxml << "\t\t\t\t</local-points>" << endl;
				fxml << "\t\t\t</local-curve>" << endl;
			}
		}
		fxml << "\t\t</local-curves>" << endl;
	}

	// border edges/vertices
	fxml << "\t\t<border-points>" << endl;
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(point->isBorder())
			fxml << "\t\t\t" << i << "\t" << (int)point->getBorderFlags() << endl;
	}
	fxml << "\t\t</border-points>" << endl;

	fxml << "\t\t<border-edges>" << endl;
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(edge->isBorder())
			fxml << "\t\t\t" << edge->getMeshPoint(0)->getIndex() << "\t" 
					<< edge->getMeshPoint(1)->getIndex() << "\t" 
					<< (int)edge->getBorderFlags() << endl;
	}
	fxml << "\t\t</border-edges>" << endl;

	fxml << "\t</surface-mesh>" << endl;
	fxml << "</meshdoc>" << endl;
	return true;
}

/// Store surface mesh to .OFF file
bool MeshContainer3dSurface::storeOFF(const string& fname) const
{
	LOG4CPLUS_INFO(MeshLog::logger_console,
		"Storing surface mesh to OFF file: " << fname);
	ofstream off(fname);
	if(!off){
		LOG4CPLUS_ERROR(MeshLog::logger_console, 
			"Error opening file for write: " << fname);
		return false;
	}
	off << "OFF" << endl;
	int pct = getPointsCount();
	int fct = getFacesCount();
	off << pct << ' ' << fct << " 0" << endl;
	const int PREC = 8;
	off.precision(PREC);
	off << fixed;
	for(int i = 0; i < pct; i++){
		const DPoint3d& pt = getPointAt(i)->getCoordinates();
		off << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
	}
	for(int i = 0; i < fct; i++){
		MeshFace* face = getFaceAt(i);
		int fct = face->getPointCount();
		off << fct;
		for(int j = 0; j < fct; j++){
			off << "\t" << face->getPoint(j)->getIndex();
		}
		off << endl;
	}
	return true;
}

/// check local surface approximation quality
void MeshContainer3dSurface::checkLocalSurfaces(Metric3dContext& mc) const
{

//	int pct = getPointsCount();
//	double max_param_dist2 = 0.0;
//	DataStatistics surf_count;
//	for(int i = 0; i < pct; i++){
//		double param_dist2 = getPointAt(i)->checkMaxLocalSurfaceParamDist2( surf_count );
//	}

/*
	CS3dPtr acs = (CS3dPtr)mc.getControlSpace();
	assert(acs->isAdaptive());

	int pct = getPointsCount();

	int surf_inconsistent_count = 0;
	int surf_not_best_count = 0;
	int surf_not_in_domain_count = 0;
	int empty_set_count = 0;
	int ipct = 0;

	DataStatistics stat_dist;
	DataStatistics lstat_norm;
	DataStatistics pstat_norm;
	DataStatistics lstat_size;

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if( point->isBorder() ) continue;
		ipct++;

		const DPoint3d& dpt = point->getCoordinates();
		SurfaceParametricSet* lset = acs->getLocalSurfaceSetAtPoint( dpt );
		SurfaceParametric* psurf = point->getLocalSurface();
		const DVector3d& pnormal = point->getBaseNormal();


		SurfaceParametric* best_surf = nullptr;

		if( lset != nullptr ){
			lstat_size.add( lset->countInt() );

			const DPoint3d spt = lset->fitToNearest( mc, dpt, pnormal, &best_surf);

			if( ! lset->containsSurface( psurf ) ) surf_inconsistent_count++;
			else if( psurf != best_surf) surf_not_best_count++;
		
			if( best_surf != nullptr ){
				double dist = mc.transformRStoMS(spt - dpt).length();
				stat_dist.add( dist );
			}
		}else
			++empty_set_count;

		DPoint2d param;
		if( psurf && psurf->withinDomain( mc, dpt, pnormal, param ) ){
			double sp = pnormal.scalarProduct( psurf->getNormalVector( param ) );
			pstat_norm.add(sp);
		}else
			++surf_not_in_domain_count;

		if( lset != nullptr ) {
			for(int j = 0; j < lset->countInt(); j++){
				SurfaceParametric* lsurf = lset->getSurface(j);
				if( lsurf->withinDomain( mc, dpt, pnormal, param ) ){
					double sp = pnormal.scalarProduct( lsurf->getNormalVector( param ) );
					lstat_norm.add(sp);
				}
			}
		}
	}

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "============ checkLocalSurfaces for " << pct << " points (" << ipct  << " inner) =========");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, " inconsistent  : " << surf_inconsistent_count);
	LOG4CPLUS_INFO(MeshLog::logger_mesh, " not in domain : " << surf_not_in_domain_count);
	LOG4CPLUS_INFO(MeshLog::logger_mesh, " empty set     : " << empty_set_count);
	LOG4CPLUS_INFO(MeshLog::logger_mesh, " not_best      : " << surf_not_best_count);
	if(stat_dist.calculate()){
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " dist min    : " << stat_dist.minimum());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " dist ave    : " << stat_dist.average());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " dist max    : " << stat_dist.maximum());
	}
	if(pstat_norm.calculate()){
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " p-sp-norm min : " << pstat_norm.minimum());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " p-sp-norm ave : " << pstat_norm.average());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " p-sp-norm max : " << pstat_norm.maximum());
	}
	if(lstat_norm.calculate()){
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " l-sp-norm min : " << lstat_norm.minimum());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " l-sp-norm ave : " << lstat_norm.average());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " l-sp-norm max : " << lstat_norm.maximum());
	}
	if(lstat_size.calculate()){
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " l-count min : " << lstat_size.minimum());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " l-count ave : " << lstat_size.average());
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " l-count max : " << lstat_size.maximum());
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "===============================================================");
*/
}

/// check surface/surface-set counters for mesh points
void MeshContainer3dSurface::checkSurfaceSetCounters(Metric3dContext& mc, bool log_data) const
{
	assert(m_control->isAdaptive());
	int pct = getPointsCount();

	static DataHashTableKeyValue< SurfaceSetConstPtr, int> set_indices(1000, nullptr);
	static DataHashTableKeyValue< SurfaceConstPtr, int> surf_indices(1000, nullptr);
	static DataVector< std::shared_ptr<DataVector<int>> > set_tables(10);
	static DataVector< std::shared_ptr<DataVector<int>> > surf_tables(10);

	bool first_time = set_tables.empty();

	auto set_counters = first_time ? std::make_shared<DataVector<int>>(1000) : 
		std::make_shared<DataVector<int>>(set_tables[set_tables.countInt()-1]->countInt(), 0);
	auto surf_counters = first_time ? std::make_shared<DataVector<int>>(1000) : 
		std::make_shared<DataVector<int>>(surf_tables[surf_tables.countInt()-1]->countInt(), 0);

	set_tables.add(set_counters);
	surf_tables.add(surf_counters);

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if( point->isBorder() ) continue;

		const DPoint3d& dpt = point->getCoordinates();
		const DVector3d& pnormal = point->getBaseNormal();
		auto lset = m_control->getAsAdaptive()->getLocalSurfaceSetAtPoint( dpt );
		auto psurf = point->getLocalValidSurface();

		int surf_index = -1;
		for(int j = 0; j < lset->countInt(); j++) {
			if(lset->getSurface(j) == psurf) {
				surf_index = j;
				break;
			}
		}
		//assert( surf_index > -1 );

		if( surf_index > -1 ) {
			int set_offset = set_indices.getValue(lset, -1);
			if( set_offset == -1 ){ // first time
				set_indices.insert( lset, set_offset = set_counters->countInt() ); 
				set_counters->addItems(lset->countInt(), 0);
			}
			(*set_counters)[ set_offset + surf_index ] ++;
		}else{
			DPoint2d param;
			if(true){
				MeshViewSet* set = new MeshViewSet;
				psurf->drawDomain( set );
				set->addPoint( point );
				SHOW_MESH( "problem", set );
			}
			assert( psurf && psurf->withinDomain( dpt, param ) );
			LOG4CPLUS_WARN(MeshLog::logger_console, "p-surf not in surf-set");
		}

		int surf_offset = surf_indices.getValue( psurf, -1 );
		if( surf_offset == -1 ){ // first time
			surf_indices.insert( psurf, surf_offset = surf_counters->countInt() );
			surf_counters->add(1);
		}else{
			(*surf_counters)[ surf_offset ] ++;
		}
	}

	if(log_data){
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
			" ===[*" << set_tables.countInt() 
			<< "]===== count point-call for surface in sets ======== ");
		DataVector<SurfaceSetConstPtr> all_sets;
		set_indices.getKeys( all_sets );
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, " === total " << all_sets.countInt() << " sets ");
		for(int i = 0; i < all_sets.countInt(); i++){
			auto set = all_sets[i];
			int set_offset = set_indices.getValue( set, -1 );
			assert( set_offset > -1 );
			for(int j = 0; j < set->countInt(); j++){
//				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, ((unsigned int)set) << "\t" 
//					<< ((unsigned int) (set->getSurface(j)));
				ostringstream log_info;
				for(int k = 0; k < set_tables.countInt(); k++)
					log_info << "\t" << ( ((set_offset+j) < set_tables[k]->countInt()) ? set_tables[k]->get(set_offset+j) : -1 );
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, log_info.str());
			}
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "--------------------");
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, " =========== count point-call for surfaces ======== ");
		DataVector<SurfaceConstPtr> all_surfaces;
		surf_indices.getKeys( all_surfaces );
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, " === total " << all_surfaces.countInt() << " surfaces ");
		for(int i = 0; i < all_surfaces.countInt(); i++){
			SurfaceConstPtr surf = all_surfaces[i];
			int surf_offset = surf_indices.getValue( surf, -1 );
			assert( surf_offset > -1 );
//			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, ((unsigned int)surf);
			ostringstream log_info;
			for(int k = 0; k < surf_tables.countInt(); k++)
				log_info << "\t" << ( (surf_offset < surf_tables[k]->countInt()) ? surf_tables[k]->get(surf_offset) : -1 );
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, log_info.str());
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "===============================================================");
	}
}

/// show surface mesh with marking local surfaces
void MeshContainer3dSurface::showLocalSurfacesForMesh(const string& label) const
{
	LOG4CPLUS_INFO(MeshLog::logger_console,
		"Preparing visualization of the mesh...");
	int fct = getFacesCount();
	DataHashTableKeyValue<SurfaceConstPtr, int> surf_ids(100, nullptr);
	int surf_id_counter = 1;
	for(int i = 0; i < fct; i++){
		MeshFace* face = getFaceAt(i);
		assert(face->hasLocalSurface());
		SurfaceConstPtr surf0 = nullptr;
		bool all_alike = true;
		int fpct = face->getPointCount();
		for(int j = 0; all_alike && (j < fpct); j++){
			MeshPoint3d* fpoint = face->getPoint(j);
			SurfaceConstPtr surfj = fpoint->getLocalSurface();
			if(fpoint->isBorder()) continue;
			if(surf0 == nullptr) surf0 = surfj;
			else all_alike = (surfj == surf0);
		}
		int id = -1;
		if(!surf0) surf0 = face->getCurrentLocalSurface();
		if(surf0) {
			if(all_alike) {
				id = surf_ids.getValue(surf0, 0);
				if(id == 0) surf_ids.insert(surf0, id = surf_id_counter++);
			}else{
				id = -1;
			}
		}
		face->setIntTag(TagExtended::TAG_VISUALIZATION, id);
	}

	for(int i = 0; i < getPointsCount(); i++) {
		MeshPoint3d* point = getPointAt(i);
		if( point->isBorder() )
			point->setIntTag(TagExtended::TAG_VISUALIZATION, 0);
		else{
			SurfaceConstPtr surf = point->getLocalSurface();
			point->setIntTag(TagExtended::TAG_VISUALIZATION, surf_ids.getValue(surf, 0) );
		}
	}
	SHOW_MESH(label, getViewSet(nullptr, TagExtended::TAG_VISUALIZATION));
}

/*
/// show surface mesh with marking local surfaces quality for vertices
void MeshContainer3dSurface::showLocalSurfacesQuality(Metric3dContext& mc, const string& label) const
{
	LOG4CPLUS_INFO(MeshLog::logger_console,"Preparing visualization of the mesh...");

	int pct = getPointsCount();
	DataVector<int> qtable( pct );
	double quality;
	DPoint2d pt_2d;
	double worst_quality = 1.0;
	for(int i = 0; i < pct; i++) {
		int q100 = -1;
		MeshPoint3d* point = getPointAt(i);
		if( point->isBorder() ) q100 = -1;
		else{
			SurfaceConstPtr surf = point->getLocalSurface();
			if( surf && surf->withinDomain( mc, point->getCoordinates(), point->getBaseNormal(), pt_2d, &quality ) ) {
				q100 = (int)(100*quality);
				if(q100 < worst_quality) worst_quality = q100;
			}else
				LOG4CPLUS_WARN(MeshLog::logger_console, "MC3S::showLocalSIR - point not within local-surf. domain" );
		}

		point->setIntTag(TagExtended::TAG_VISUALIZATION, q100 );
		qtable.add(q100);
	}

	LOG4CPLUS_INFO(MeshLog::logger_console, "MC3S::showLocalSQuality - worst quality", worst_quality );

	int fct = getFacesCount();
	for(int i = 0; i < fct; i++){
		MeshFace* face = getFaceAt(i);
		int min_q = 100;
		int fpct = face->getPointCount();
		for(int j = 0; j < fpct; j++){
			MeshPoint3d* fpoint = face->getPoint(j);
			int qp = qtable[ fpoint->getIndex() ];
			if( qp > min_q ) min_q = qp;
		}
		face->setDoubleTag(TagExtended::TAG_QUALITY, 0.01*min_q);
	}

	MeshViewSet* set = getViewSet( nullptr, TagExtended::TAG_VISUALIZATION );
	set->setPolygonFillMode( MeshViewSet::FILL_QUALITY );
	SHOW_MESH( label, set );
}
*/

/// show surface mesh with marking faces with sharp edges
void MeshContainer3dSurface::showFacesWithSharpEdges(const string& label) const
{
	LOG4CPLUS_INFO(MeshLog::logger_console,"Preparing visualization of the mesh...");
	int pct = getPointsCount();
	int fct = getFacesCount();
	for(int i = 0; i < fct; i++)
		getFaceAt(i)->setIntTag( TagExtended::TAG_VISUALIZATION, 0 );

	int sharp_counter = 0;
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(edge->isBorder()) continue; 
		assert( edge->getFaceCount() == 2); // since it is not marked as border

		MeshFace* f0 = edge->getFaceAt(0);
		MeshFace* f1 = edge->getFaceAt(1);

		DVector3d dn0, dn1;
		if( ! f0->checkAndGetNormalVector(dn0) ) continue;
		if( ! f1->checkAndGetNormalVector(dn1) ) continue;

		double sp = dn0.scalarProduct(dn1);
		if(sp < MeshGenerator3dSurface::param_sharp_edge_threshold){
			f0->setIntTag( TagExtended::TAG_VISUALIZATION, 3 );
			f1->setIntTag( TagExtended::TAG_VISUALIZATION, 3 );
			sharp_counter++;
		}
	}

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(point->isBorder()) continue;
		SurfaceConstPtr psurf = point->getLocalSurface();
		if(psurf == nullptr) continue;

		const DVector3d & pbase_normal = point->getBaseNormal();
		const DVector3d   psurf_normal = psurf->getNormalVector( point->getLocalSurfaceParam(psurf) );

		DVector3d vbase_normal;
		DVector3d vsurf_normal;
		int vbct = 0, vsct = 0;
		for(int j = 0; j < point->getRank(); j++){
			MeshPoint3d* other_point = point->getEdge(j)->getOtherPoint( point );
			if( other_point->isBorder() ) continue;
			vbase_normal += other_point->getBaseNormal();
			vbct++;
			SurfaceConstPtr osurf = other_point->getLocalSurface();
			if(osurf == nullptr) continue;
			vsurf_normal += osurf->getNormalVector( other_point->getLocalSurfaceParam(osurf) );
			vsct++;
		}
		if(vbct > 0) vbase_normal /= vbct;
		if(vsct > 0) vsurf_normal /= vsct;

		double sp_v = vbase_normal.scalarProduct( vsurf_normal );
		double sp_p = pbase_normal.scalarProduct( psurf_normal );

		double sp_b_s = pbase_normal.scalarProduct( vsurf_normal );
		double sp_s_s = psurf_normal.scalarProduct( vsurf_normal );

		const double MIN_SP = 0.9;
		if( sp_v < MIN_SP || sp_p < MIN_SP || sp_b_s < MIN_SP || sp_s_s < MIN_SP ){
			MeshViewSet* set = getDebugViewSetTopological(point, 3, TagExtended::TAG_VISUALIZATION);
			set->addInfo("sp_v | sp_p", to_string(sp_v) + " | " + to_string(sp_p) );
			set->addInfo("sp_b_s | sp_s_s", to_string(sp_b_s) + " | " + to_string(sp_s_s) );
			SHOW_MESH( label + " (sp-normals)", set );
		}
	}

	MeshViewSet* set = getViewSet( nullptr, TagExtended::TAG_VISUALIZATION );
	set->addInfo("sharp edges count", sharp_counter);
	SHOW_MESH( label, set );
}

/// Show faces/points in the local sub-domain - just for test
void MeshContainer3dSurface::showDebugLocalFaces(const string& caption, 
								DataVector<int> & p_list, DataVector<int> & p_layers,
								DataVector<int> & f_list, DataVector<int> & f_layers,
								int l_outside, DataVector<int> * pref ) const
{
	DataVector<int> * local_pref = nullptr;
	if( pref == nullptr ){
		local_pref = new DataVector<int>( getPointsCount(), -1 );
		int lpct = p_list.countInt();
		for(int i = 0; i < lpct; i++) 
			if( p_layers[ p_list[i] ] < l_outside )
				local_pref->set( p_list[i], i );
		pref = local_pref;
	}

	MeshViewSet* dset = new MeshViewSet();
	int lfct = f_list.countInt();
	for(int i = 0; i < lfct; i++) {
		int fid = f_list[i];
		if( f_layers[fid] >= l_outside ) continue;
		MeshFace* f = getFaceAt( fid );
		if( f->getPointCount() != 3 ) continue;
		int ta = pref->get( f->getPoint(0)->getIndex() );
		int tb = pref->get( f->getPoint(1)->getIndex() );
		int tc = pref->get( f->getPoint(2)->getIndex() );
		bool valid = ta >= 0 && tb >= 0 && tc >= 0;
		dset->addFaceWithEdges( f, valid ? 0 : 1 );
	}
	SHOW_MESH( caption, dset);

	if(local_pref) delete local_pref;
}

/// for-each face
void MeshContainer3dSurface::forEachFace( std::function <void(MeshFace* face)> f )
{
	int fct = m_faces->countInt();
	for(int i = 0; i < fct; i++) f( m_faces->getDataAt(i) );
}

/// for-each point
void MeshContainer3dSurface::forEachPoint( std::function <void(MeshPoint3d* point)> f )
{
	int pct = m_points->countInt();
	for(int i = 0; i < pct; i++) f( m_points->getDataAt(i) );
}
