// MeshGenerator3dSurface.cpp: implementation of the MeshGenerator3dSurface class.
//
//////////////////////////////////////////////////////////////////////

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "MeshGenerator3dSurface.h"

#include "MeshContainer3d.h"
#include "MeshContainer3dSurface.h"
#include "MeshFace.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "DataVector.h"
#include "DataHashTable.h"
#include "MeshDomainVolume.h"
#include "DataMatrix.h"
#include "DLeastSquaresFitting.h"
#include "SurfaceParametric.h"
#include "SurfacePlane.h"
#include "SurfacePlanarQuadric.h"
#include "Curve3dParametric.h"
#include "DTriangle.h"
#include "Metric3dContext.h"
#include "MeshTriangle3d.h"
#include "MeshQuad3d.h"
#include "MeshPoly3d.h"
#include "ControlSpace3dOctree.h"
#include "ControlSpace2dIdentity.h"
#include "DPlanarQuadric.h"
#include "DLine.h"
#include "DPlane.h"
#include "DTriangle.h"
#include "MeshGenerator2d.h"
#include "MeshContainer2d.h"
#include "MeshTriangle2d.h"
#include "SurfaceAnalytic.h"
#include "DataHeapVector.h"
#include "DataPtrMatrix.h"

#include <vector>
#include <algorithm>

#include "MeshViewSet.h"

/// Parameter for sharp edges identification (maximum scalar product value for normals)
double MeshGenerator3dSurface::param_sharp_edge_threshold = 0.5;
/// Parameter for approximation of local surfaces and curves
double MeshGenerator3dSurface::param_local_shape_tolerance = 0.15;
/// Method of local surface identification
int MeshGenerator3dSurface::param_local_surface_identification_method;
// -> initialized in MeshData::setDefaultParameters ...
//	= MeshGenerator3dSurface::LSI_HIGH_CURVATURE_FIRST;
//	= MeshGenerator3dSurface::LSI_VIA_NORMALS;

const double MeshGenerator3dSurface::MIN_EDGE_LEN = 0.68;
const double MeshGenerator3dSurface::MIN_BEDGE_LEN = 0.65;
const double MeshGenerator3dSurface::MAX_BEDGE_LEN_JOINED = 1.40;
const double MeshGenerator3dSurface::MAX_EDGE_LEN = 1.5;
const double MeshGenerator3dSurface::MAX_BEDGE_LEN = 1.80;

const double MeshGenerator3dSurface::LOCAL_SURFACE_TOL_F[] = { 0.5, 0.8, 1.0, 1.5, 2.0 }; // tolerance coefficient
const int    MeshGenerator3dSurface::LOCAL_SURFACE_LAY_F[] = {   3,   3,   3,   2,   1 }; // number of topological layers
const double MeshGenerator3dSurface::LOCAL_SURFACE_MTR_F[] = { 2.0, 1.5, 1.0, 0.5, 0.0 }; // metric radius for geometrical layer
const int	 MeshGenerator3dSurface::LOCAL_SURFACE_TOL_F_CT = sizeof(MeshGenerator3dSurface::LOCAL_SURFACE_TOL_F) / sizeof(double);
const int	 MeshGenerator3dSurface::LOCAL_SURFACE_TOL_F_CT_WHOLE = 3;	// for whole-surface approximation


#ifdef _DEBUG
	#define REMESH_SHOW
	#define REMESH_STATS
#endif

/// Modify crack surface mesh according to the given CS3d, with local operators (edge collapsing and smoothing)
bool MeshGenerator3dSurface::remeshSurfaceMeshWithLocalTransformations(Metric3dContext& mc, MeshContainer3dSurface* mesh)
{
#ifdef REMESH_SHOW
	bool any_curve = false;
	if(true){
		mesh->showLocalSurfacesForMesh( "Surface mesh - WLT start" );
	//	mesh->showLocalSurfacesQuality(mc, "Surface mesh local-surf-quality - WLT start" );

		MeshViewSet* set = new MeshViewSet();
		DataHashTableKeyValue<Curve3dConstPtr, int> curve_ids(100, nullptr);
		int curve_id_counter = 1;
		int pct = mesh->getPointsCount();
		for(int i = 0; i < pct; i++){
			MeshPoint3d* point = mesh->getPointAt(i);
			if(! point->isBorder() ) continue;
			auto curve0 = point->getLocalCurve();
			int id = 0;
			if(curve0){
				id = curve_ids.getValue(curve0, 0);
				if(id == 0) curve_ids.insert(curve0, id = curve_id_counter++);
			}
			set->addPoint( point->getCoordinates(), id, id);
			for(int j = 0; j < point->getRank(); j++){
				MeshEdge3d* edge = point->getEdge(j);
				if(! edge->isBorder() ) continue;
				MeshPoint3d* other_point = edge->getOtherPoint(point);
				auto curve1 = other_point->getLocalCurve();
				set->addEdge(point->getCoordinates(), other_point->getCoordinates(), (curve0==curve1) ? id : -2);
			}
		}
		any_curve = (curve_id_counter > 1);
		if(any_curve)
			SHOW_MESH("Curve edges", set);
		else
			delete set; // no border edge, actually...
	}
#endif

#ifdef REMESH_STATS
	if(true){
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges: " <<  mesh->checkInnerEdgesLength(mc));
		double be_len = mesh->checkBorderEdgesLength(mc);
		if(be_len >= 0.0)  // if any...
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges: " << be_len);
	}
#endif

#ifdef REMESH_STATS
	mesh->checkLocalSurfaces(mc);
#endif

	START_CLOCK("MG3dS:remeshSurfaceWLT");
	auto sm_laplace = MeshData::SM_LAPLACE_MIXED;

	for(int ii = 0; ii < 2; ii++){
		int iter = 0;
		bool refill_edges = false;
		static const int MIN_RE = 5;
		int re = MIN_RE+1; // removed edges during remeshing
		while(re > MIN_RE){

			iter++;
	#ifdef REMESH_STATS
			ostringstream log_info;
			log_info << "Iteration " << iter;
	#endif

			LOG4CPLUS_DEBUG(MeshLog::logger_console,
				"Checking and fixing inverted faces for local surfaces...");
			int inv_count = MeshGenerator3dSurface::checkAndFixInvertedFacesForLocalSurfaces(mc, mesh);
	#ifdef REMESH_STATS
			log_info << ", inv=" << inv_count;
	#endif

	//		mesh->showLocalSurfacesForMesh( " TEST after fix inverted" );

//			sm_laplace = (iter < 3) ? MeshData::SM_LAPLACE : MeshData::SM_LAPLACE_METRIC;

			MeshGenerator3dSurface::smoothen(mc, mesh, 1, TagExtended::TAG_NONE, 0, 
				sm_laplace | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE); 
	//		assert(mesh->isValid());
	//		SHOW_MESH_NORESET("after smoothen-1", mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_DOMAIN));
			MeshGenerator3dSurface::smoothenLaplaceBoundary(mc, mesh); 
	//		SHOW_MESH_NORESET("after smoothen-b", mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_DOMAIN));
	//		assert(mesh->isValid());
			MeshGenerator3dSurface::smoothen(mc, mesh, 1, TagExtended::TAG_NONE, 0, 
				sm_laplace | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE); 

	//		mesh->showLocalSurfacesForMesh( " TEST after smoothen" );

			SHOW_MESH_NORESET("after smoothen-ibi", 
				mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_DOMAIN));

			int inner_re = MeshGenerator3dSurface::collapseInnerEdges(mc, mesh, false, refill_edges);
	#ifdef REMESH_STATS
			log_info << ", inner_re=" << inner_re;
	#endif

	//		mesh->showLocalSurfacesForMesh( " TEST after collapse inner edges" );

			SHOW_MESH_NORESET("after collapse inner", 
				mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));

			int border_re = MeshGenerator3dSurface::collapseBoundaryEdges(mc, mesh, false, refill_edges);
	#ifdef REMESH_STATS
			log_info << ", border_re=" << border_re;
			LOG4CPLUS_INFO(MeshLog::logger_mesh, log_info.str());

			//mesh->showLocalSurfacesForMesh( " TEST after collapse boundary edges" );

			if(true){
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges: " <<  mesh->checkInnerEdgesLength(mc));
				double be_len = mesh->checkBorderEdgesLength(mc);
				if(be_len >= 0.0)  // if any...
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges: " << be_len);
			}
	#endif
			re = inner_re + border_re;

		#ifdef REMESH_STATS
			mesh->checkLocalSurfaces(mc);
		#endif

			if(refill_edges) break;
			if(re < 50) refill_edges = true;

	#ifdef REMESH_SHOW
			if(false){
				LOG4CPLUS_INFO(MeshLog::logger_console,"Preparing visualization of the mesh...");
				int fct = mesh->getFacesCount();
				DataHashTableKeyValue<SurfaceConstPtr, int> surf_ids(100, nullptr);
				int surf_id_counter = 1;
				for(int i = 0; i < fct; i++){
					MeshFace* face = mesh->getFaceAt(i);
					SurfaceConstPtr surf0 = face->getPoint(0)->getLocalSurface();
					SurfaceConstPtr surf1 = face->getPoint(1)->getLocalSurface();
					SurfaceConstPtr surf2 = face->getPoint(2)->getLocalSurface();
					int id = 0;
					if(surf0) {
						if(surf0 == surf1 && surf0 == surf2) {
							id = surf_ids.getValue(surf0, 0);
							if(id == 0) surf_ids.insert(surf0, id = surf_id_counter++);
						}else{
							id = -1;
						}
					}
					face->setIntTag(TagExtended::TAG_VISUALIZATION, id);
				}
				ostringstream ostr;
				ostr << "Surface mesh after iteration #" << iter << ", removed edges: " << re;
				SHOW_MESH_NORESET(ostr.str(), mesh->getViewSet(nullptr, TagExtended::TAG_VISUALIZATION));

				if(any_curve){
					MeshViewSet* set = new MeshViewSet();
					DataHashTableKeyValue<Curve3dConstPtr, int> curve_ids(100, nullptr);
					int curve_id_counter = 1;
					int pct = mesh->getPointsCount();
					for(int i = 0; i < pct; i++){
						MeshPoint3d* point = mesh->getPointAt(i);
						if(! point->isBorder() ) continue;
						auto curve0 = point->getLocalCurve();
						int id = 0;
						if(curve0){
							id = curve_ids.getValue(curve0, 0);
							if(id == 0) curve_ids.insert(curve0, id = curve_id_counter++);
						}
						set->addPoint( point->getCoordinates(), id, id);
						for(int j = 0; j < point->getRank(); j++){
							MeshEdge3d* edge = point->getEdge(j);
							if(! edge->isBorder() ) continue;
							MeshPoint3d* other_point = edge->getOtherPoint(point);
							auto curve1 = other_point->getLocalCurve();
							set->addEdge(point->getCoordinates(), other_point->getCoordinates(), (curve0==curve1) ? id : -2);
						}
					}
					SHOW_MESH("Curve edges", set);
				}
			}
	#endif
		}

		//mesh->showLocalSurfacesForMesh( " TEST before split + smoothen" );

		MeshGenerator3dSurface::smoothen(mc, mesh, 1, TagExtended::TAG_NONE, 0, 
			sm_laplace | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE); 
		int inserted = MeshGenerator3dSurface::splitEdges(mc, mesh);
		MeshGenerator3dSurface::smoothen(mc, mesh, 1, TagExtended::TAG_NONE, 0, 
			sm_laplace | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE); 

		mesh->showLocalSurfacesForMesh( " TEST after split + smoothen" );

	#ifdef REMESH_STATS
		if(true){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Inserted nodes: " << inserted);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges: " <<  mesh->checkInnerEdgesLength(mc));
			double be_len = mesh->checkBorderEdgesLength(mc);
			if(be_len >= 0.0)  // if any...
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges: " << be_len);
		}
	#endif

		//optimizeLocalSurfaces( mc, mesh );

	#ifdef REMESH_STATS
		mesh->checkLocalSurfaces(mc);
	#endif

		MeshGenerator3dSurface::optimizeBoundaryEdges(mc, mesh);
		MeshGenerator3dSurface::smoothen(mc, mesh, 1, TagExtended::TAG_NONE, 0, 
			sm_laplace | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE); 

		mesh->showLocalSurfacesForMesh( " TEST after optimize + smoothen" );


	#ifdef REMESH_STATS
		if(true){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges: " <<  mesh->checkInnerEdgesLength(mc));
			double be_len = mesh->checkBorderEdgesLength(mc);
			if(be_len >= 0.0)  // if any...
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges: " << be_len);
		}
	#endif

	}

#ifdef _DEBUG
	//mesh->showLocalSurfacesQuality(mc, "Surface mesh local-surface-quality - WLT finish" );
	//mesh->showFacesWithSharpEdges( "Surface mesh sharp edges - WLT finish" );
#endif

	#ifdef REMESH_SHOW
	if(false){
		mesh->showLocalSurfacesForMesh( "Surface mesh, final" );

		if(any_curve){
			MeshViewSet* set = new MeshViewSet();
			DataHashTableKeyValue<Curve3dConstPtr, int> curve_ids(100, nullptr);
			int curve_id_counter = 1;
			int pct = mesh->getPointsCount();
			for(int i = 0; i < pct; i++){
				MeshPoint3d* point = mesh->getPointAt(i);
				if(! point->isBorder() ) continue;
				auto curve0 = point->getLocalCurve();
				int id = 0;
				if(curve0){
					id = curve_ids.getValue(curve0, 0);
					if(id == 0) curve_ids.insert(curve0, id = curve_id_counter++);
				}
				set->addPoint( point->getCoordinates(), id, id);
				for(int j = 0; j < point->getRank(); j++){
					MeshEdge3d* edge = point->getEdge(j);
					if(! edge->isBorder() ) continue;
					MeshPoint3d* other_point = edge->getOtherPoint(point);
					auto curve1 = other_point->getLocalCurve();
					set->addEdge(point->getCoordinates(), other_point->getCoordinates(), (curve0==curve1) ? id : -2);
				}
			}
			SHOW_MESH("Curve edges, final", set);
		}
	}
	#endif

#ifdef _DEBUG
	if(false){
		const int PIds[] = { 1, 3, 5 };
		for(int pid : PIds)
			checkPointForLocalSurfaces( mc, mesh, mesh->getPointAt( pid ) );
	}
#endif

#ifdef REMESH_SHOW
	int bad_counter = 0;
	for(int i = 0; i < mesh->getFacesCount(); i++){
		MeshFace* face = mesh->getFaceAt(i);
		mc.countMetricAtPoint(face->getMiddlePoint());
		double q = DMTriangle3d::alphaQuality(
			face->getPoint(0)->getMetricCoordinates(mc),
			face->getPoint(1)->getMetricCoordinates(mc),
			face->getPoint(2)->getMetricCoordinates(mc));
		if(q < 0.1){
			LOG4CPLUS_DEBUG(MeshLog::logger_console, "alpha quality: " << q);
			++bad_counter;
			//SHOW_MESH("Bad surface element", mesh->getDebugViewSet(face->getPoint(0), face->getPoint(1), 2.0));
		}
	}
	if(bad_counter > 0) LOG4CPLUS_DEBUG(MeshLog::logger_console, "faces with bad alpha quality: " << bad_counter);
#endif

#ifdef _DEBUG
	if(true){ // show edges and bad face-edges
		LOG4CPLUS_INFO(MeshLog::logger_console,"Preparing visualization of the mesh...");
		MeshViewSet* set = new MeshViewSet;
		int bad_ct = 0;
		for(IteratorEdge3d it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
			MeshEdge3d* edge = it.getEdge();
			set->addEdge(edge);
			if(edge->isBorder() && edge->getFaceCount() == 2){
				MeshFace* face0 = edge->getFaceAt(0);
				MeshFace* face1 = edge->getFaceAt(1);
				double sc = face0->getNormalVector().scalarProduct(face1->getNormalVector());
				if(!face0->sameOrientation(face1))
					sc = -sc;
				if(sc <= -0.9){ // 6o
					bad_ct++;
					set->addFace(face0, -2, MeshViewSet::param_shrink, face0->getBlock(0) != nullptr);
					set->addFace(face1, -2, MeshViewSet::param_shrink, face0->getBlock(0) != nullptr);
				}
			}
		}
		if(bad_ct > 0){
			SHOW_MESH("Bad face-angles", set);
		}else
			delete set;
	}
#endif

	STOP_CLOCK("MG3dS:remeshSurfaceWLT");

	return true;
}

/// Modify crack surface mesh according to the given CS3d, with local operators (edge collapsing and smoothing)
bool MeshGenerator3dSurface::remeshCrackSurfaceMeshWithLocalTransformations(Metric3dContext& mc, MeshContainer3dSurface* mesh)
{
#ifdef REMESH_SHOW
	SHOW_MESH("Surface mesh", mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));
#endif

#ifdef REMESH_STATS
	if(true){
		double ave_len = mesh->checkInnerEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 1);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges (1): " << ave_len);
		ave_len = mesh->checkInnerEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 2);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges (2): " << ave_len);
		ave_len = mesh->checkBorderEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 1, 2);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges (1-2): " << ave_len);
		ave_len = mesh->checkBorderEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 2, 2);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges (2-2): " << ave_len);
	}
#endif

	START_CLOCK("MG3dS:remeshSurfaceWLT");

	int iter = 0;
	bool refill_edges = false;
	static const int MIN_RE = 5;
	int re = MIN_RE+1; // removed edges during remeshing
	while(re > MIN_RE){

		iter++;
#ifdef REMESH_STATS
		ostringstream log_info;
		log_info << "Iteration " << iter;
#endif

		int inv_count = MeshGenerator3dSurface::checkAndFixInvertedFacesForLocalSurfaces(mc, mesh, TagExtended::TAG_ADAPT_SURF, 1);
		inv_count += MeshGenerator3dSurface::checkAndFixInvertedFacesForLocalSurfaces(mc, mesh, TagExtended::TAG_ADAPT_SURF, 2);
#ifdef REMESH_STATS
		log_info << ", inv=" << inv_count;
#endif

//		mesh->moveInnerPointsToLocalSurfaces(mc, TagExtended::TAG_ADAPT_SURF, 1, 4);
//		mesh->moveBorderPointsToLocalSurfaces(mc, TagExtended::TAG_ADAPT_SURF, 1+2, 4);

//		SHOW_MESH_NORESET("after move", mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));

//		MeshGenerator3dSurface::improveNearBorder(mc, mesh, TagExtended::TAG_ADAPT_SURF, 1);
//		MeshGenerator3dSurface::improveNearBorder(mc, mesh, TagExtended::TAG_ADAPT_SURF, 2);

//		SHOW_MESH_NORESET("after improve-near-border", mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));

		MeshGenerator3dSurface::smoothen(mc, mesh, 1, TagExtended::TAG_ADAPT_SURF, 1+2); // for points, tags are treated as flags
//		assert(mesh->isValid());
		MeshGenerator3dSurface::smoothenLaplaceBoundary(mc, mesh, TagExtended::TAG_ADAPT_SURF, 1+2, 4); // for points, tags are treated as flags
//		assert(mesh->isValid());
		MeshGenerator3dSurface::smoothen(mc, mesh, 1, TagExtended::TAG_ADAPT_SURF, 1+2); // for points, tags are treated as flags

//		SHOW_MESH_NORESET("after smoothen", mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));

		int inner_re = 0;
		inner_re += MeshGenerator3dSurface::collapseInnerEdges(mc, mesh, false, refill_edges, TagExtended::TAG_ADAPT_SURF, 1);
		inner_re += MeshGenerator3dSurface::collapseInnerEdges(mc, mesh, false, refill_edges, TagExtended::TAG_ADAPT_SURF, 2);
#ifdef REMESH_STATS
		log_info << ", inner_re=" << inner_re;
#endif

//		SHOW_MESH_NORESET("after collapse inner", mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));

		int border_re = 0;
		border_re += MeshGenerator3dSurface::collapseBoundaryEdges(mc, mesh, false, refill_edges, TagExtended::TAG_ADAPT_SURF, 1, 2);
		border_re += MeshGenerator3dSurface::collapseBoundaryEdges(mc, mesh, false, refill_edges, TagExtended::TAG_ADAPT_SURF, 2, 2);
#ifdef REMESH_STATS
		log_info << ", border_re=" << border_re;
		LOG4CPLUS_INFO(MeshLog::logger_mesh, log_info.str());

		if(true){
			double ave_len = mesh->checkInnerEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 1);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges (1): " << ave_len);
			ave_len = mesh->checkInnerEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 2);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges (2): " << ave_len);
			ave_len = mesh->checkBorderEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 1, 2);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges (1-2): " << ave_len);
			ave_len = mesh->checkBorderEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 2, 2);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges (2-2): " << ave_len);
		}
#endif
		re = inner_re + border_re;

		if(refill_edges) break;
		if(re < 50) refill_edges = true;

#ifdef REMESH_SHOW
		ostringstream ostr;
		ostr << "Surface mesh after iteration #" << iter << ", removed edges: " << re;
		SHOW_MESH_NORESET(ostr.str(), mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));
#endif
	}

	//mesh->moveInnerPointsToLocalSurfaces(mc, TagExtended::TAG_ADAPT_SURF, 1);
	//mesh->moveBorderPointsToLocalSurfaces(mc, TagExtended::TAG_ADAPT_SURF, 1+2);

	//MeshGenerator3dSurface::improveNearBorder(mc, mesh, TagExtended::TAG_ADAPT_SURF, 1);
	//MeshGenerator3dSurface::improveNearBorder(mc, mesh, TagExtended::TAG_ADAPT_SURF, 2);

	MeshGenerator3dSurface::smoothen(mc, mesh, 2, TagExtended::TAG_ADAPT_SURF, 1+2); // for points, tags are treated as flags

#ifdef REMESH_STATS
	if(true){
		double ave_len = mesh->checkInnerEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 1);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges (1): " << ave_len);
		ave_len = mesh->checkInnerEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 2);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for inner edges (2): " << ave_len);
		ave_len = mesh->checkBorderEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 1, 2);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges (1-2): " << ave_len);
		ave_len = mesh->checkBorderEdgesLength(mc, TagExtended::TAG_ADAPT_SURF, 2, 2);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Average length for border edges (2-2): " << ave_len);
	}
#endif

#ifdef REMESH_SHOW
	SHOW_MESH_NORESET("Surface mesh, final", mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));
#endif

	//if(true){
	//	MeshPoint3d* tpt1 = mesh->getPointAt(1927);
	//	MeshPoint3d* tpt2 = mesh->getPointAt(3043);
	//	//MeshEdge3d* edge = tpt1->getEdgeToPoint(tpt2);
	//	//mc.countMetricAtPoints(tpt1, tpt2);
	//	//double len = edge->getLength(mc);
	//	//LOG4CPLUS_INFO(MeshLog::logger_console, "44-1676 metric len", len);
	//	SHOW_MESH("test points", mesh->getDebugViewSet(tpt1, tpt2, 
	//			2.0, TagExtended::TAG_ADAPT_SURF));
	//}

#ifdef REMESH_SHOW
	for(int i = 0; i < mesh->getFacesCount(); i++){
		MeshFace* face = mesh->getFaceAt(i);
		mc.countMetricAtPoint(face->getMiddlePoint());
		double q = DMTriangle3d::alphaQuality(
			face->getPoint(0)->getMetricCoordinates(mc),
			face->getPoint(1)->getMetricCoordinates(mc),
			face->getPoint(2)->getMetricCoordinates(mc));
		if(q < 0.1){
			LOG4CPLUS_DEBUG(MeshLog::logger_console, "alpha quality: " << q);
			SHOW_MESH("Bad surface element", mesh->getDebugViewSet(face->getPoint(0), face->getPoint(1), 
				2.0, TagExtended::TAG_ADAPT_SURF));
		}
	}
#endif

#ifdef _DEBUG
	if(true){ // show edges and bad face-edges
		MeshViewSet* set = new MeshViewSet;
		int bad_ct = 0;
		for(IteratorEdge3d it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
			MeshEdge3d* edge = it.getEdge();
			set->addEdge(edge);
			if(edge->isBorder() && edge->getFaceCount() == 2){
				MeshFace* face0 = edge->getFaceAt(0);
				MeshFace* face1 = edge->getFaceAt(1);
				double sc = face0->getNormalVector().scalarProduct(face1->getNormalVector());
				if(!face0->sameOrientation(face1))
					sc = -sc;
				if(sc <= -0.9){ // 6o
					bad_ct++;
					int face_tag = face0->getIntTag(TagExtended::TAG_ADAPT_SURF, 0);
					set->addFace(face0, face_tag, MeshViewSet::param_shrink, face0->getBlock(0) != nullptr);
					face_tag = face1->getIntTag(TagExtended::TAG_ADAPT_SURF, 0);
					set->addFace(face1, face_tag, MeshViewSet::param_shrink, face0->getBlock(0) != nullptr);
				}
			}
		}
		if(bad_ct > 0){
			SHOW_MESH("Bad face-angles", set);
		}else
			delete set;
	}
#endif

	STOP_CLOCK("MG3dS:remeshSurfaceWLT");
	return true;
}

/// Collapse (all) vertices rank two
int MeshGenerator3dSurface::collapseAllInnerVerticesRankLow(Metric3dContext& mc, MeshContainer3dSurface* mesh)
{
	int counter = 0;
	int i = 0;
	int pct = mesh->getPointsCount();
	bool modified;
	do{
		modified = false;
		while(i < pct){
			MeshPoint3d* point = mesh->getPointAt(i);
			if((!point->isBorder()) && 
				( ((point->getRank() == 2) && collapseInnerVertexRankTwo(mc, mesh, point)) ||
				  ((point->getRank() == 3) && 
				   ( (mesh->getBorderStage() == BORDER_KNOWN) || point->planar(mc) ) &&
				   removeTriangulationPointSimple(mesh, point))) ){
				counter++;
				pct--;
				modified = true;
			}else
				i++;
		}
	}while(modified);

	return counter;
}

/// Collapse too short edge, special case for poly-mesh and vertex rank two
bool MeshGenerator3dSurface::collapseInnerVertexRankTwo(Metric3dContext& mc, MeshContainer3dSurface* mesh,
	MeshPoint3d* removed_point, bool no_border_points, bool refill_edges, 
	DataContainer<MeshEdge3d::ActiveEdge> * act_edges, TagExtended::TagType act_tag)
{
	assert(removed_point != nullptr);
	if(removed_point->isBorder()) return false;
	if(removed_point->getRank() != 2) return false;

//	static int special_counter = 0;

	MeshEdge3d* edges[2] = { removed_point->getEdge(0), removed_point->getEdge(1) };

	MeshPoint3d* p0 = edges[0]->getOtherPoint(removed_point);
	MeshPoint3d* p1 = edges[1]->getOtherPoint(removed_point);

	MeshFace* f[2] = { edges[0]->getFaceAt(0), edges[0]->getFaceAt(1) };

	bool join_faces = false;

//	join_faces = (f0->getEdgeCount() == 4) || (f1->getEdgeCount() == 4);

	if(!join_faces){
		// check, if collapsing would produce inner edge with four incident faces ...
		MeshEdge3d* temp_edge = p0->getEdgeToPoint(p1);
		join_faces = (temp_edge && !f[0]->incidentToEdge(temp_edge) && !f[1]->incidentToEdge(temp_edge));
	}

	if(!join_faces){

		double sp = f[0]->getBaseNormal().scalarProduct( f[1]->getBaseNormal() );
		if(sp >= 0.9) {
			// check if some other vertex of one of the faces lies within ...
			DTriangle3d tri(p0->getCoordinates(), removed_point->getCoordinates(), p1->getCoordinates());
			double tdet = tri.det();
			if(std::abs(tdet) > VERY_SMALL_NUMBER) {
				DPlane plane(tri.pt_b, tri.pt_a - tri.pt_b, tri.pt_c - tri.pt_b);
				DPoint2d pt_a = plane.projectToPlane(tri.pt_a);
				DPoint2d pt_b = plane.projectToPlane(tri.pt_b);
				DPoint2d pt_c = plane.projectToPlane(tri.pt_c);

				if( DTriangle2d::det(pt_a, pt_b, pt_c) < 0.0 ) {
					DPoint2d pt_x = pt_c;
					pt_c = pt_b; pt_b = pt_x;
				}

				for(int j = 0; (!join_faces) && (j < 2); j++){
					int fpct = f[j]->getPointCount();
					for(int i = 0; (!join_faces) && (i < fpct); i++){
						MeshPoint3d* point = f[j]->getPoint(i);
						if(point == removed_point) continue;
						if(point == p0 || point == p1) continue;
						join_faces = DTriangle2d::containsPoint( pt_a, pt_b, pt_c, 
							plane.projectToPlane(point->getCoordinates()) );
					}
				}
			}
		}
	}

//	++special_counter;
	bool show_case = false;

	if(show_case){
		ostringstream ostr;
		ostr << "Selected edge to collapse - special";
		if(join_faces) ostr << " (join)";
//		ostr << " #" << special_counter;
		ostr << ", real-len=" << p0->getCoordinates().distance(p1->getCoordinates()) << endl;
		ostr << ": p" << p0->getIndex() << " - " << (p0->isBorder() ? "border" : "inner");
		ostr << ", p" << p1->getIndex() << " - " << (p1->isBorder() ? "border" : "inner");
		SHOW_MESH(ostr.str(), mesh->getDebugViewSetTopological(removed_point, 2));
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, ostr.str());
	}

	// remove active-tag from edges (if available)
	if(act_edges){
		for(int j = 0; j < 2; j++){
			MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)edges[j]->getPtrTag(act_tag);
			if(aact){
				edges[j]->removeTag(act_tag);
				delete act_edges->removeDataItem(aact->getIndex());
			}
		}
	}

	if(join_faces){
		int fpct[2] = { f[0]->getPointCount(), f[1]->getPointCount() };
		int poly_ct = fpct[0]+fpct[1]-4;
		assert(poly_ct > 2);
		// create chain of points for a new (poly)face
		DataVector<MeshPoint3d*> poly_points( poly_ct );
		for(int j = 0; j < 2; j++){
			int ir = 0;
			while( f[j]->getPoint(ir) != removed_point) ir++;
			for(int i = 2; i < fpct[j]; i++) {
				ir = (ir+1)%fpct[j];
				poly_points.add(f[j]->getPoint(ir));
			}
		}
		// create new face
		MeshFace* fn = nullptr;
		switch(poly_ct){
		case 3:
			fn = new MeshTriangle3d(poly_points[0], poly_points[1], poly_points[2], f[0]->getBlock(0));
			break;
		case 4:
			fn = new MeshQuad3d(poly_points[0], poly_points[1], poly_points[2], poly_points[3], f[0]->getBlock(0));
			break;
		default:
			fn = new MeshPoly3d(poly_points, f[0]->getBlock(0));
			break;
		}
		fn->copySurfaceData(f[0]);
		fn->copyAllTags(f[0]);
		// replace faces and point with the new face
		delete mesh->removeMeshFace(f[0]);
		delete mesh->removeMeshFace(f[1]);
		delete mesh->removeMeshPoint(removed_point);
		mesh->addMeshFace(fn);

		//if(true){
		//	MeshViewSet* set = new MeshViewSet();
		//	set->addFaceWithEdges(fn);
		//	for(int i = 0; i < fn->getPointCount(); i++){
		//		set->addPoint( fn->getPoint(i) );
		//	}
		//	SHOW_MESH("Polygon after join-edge", set);
		//}

	}else{
		// transform mesh - remove this point
		for(int j = 0; j < 2; j++){
			MeshFace* reduced_face = f[j]->removePoint(removed_point);
		// triangle->nullptr, quad->triangle, poly(5)->quad, poly(n)->poly(n-1)
			if(reduced_face != f[j]){
				delete mesh->removeMeshFace(f[j]);
				if(reduced_face != nullptr) mesh->addMeshFace(reduced_face);
				f[j] = reduced_face;
			}
		}

		delete mesh->removeMeshPoint(removed_point);

		// check created edge
		MeshEdge3d* temp_edge = p0->getEdgeToPoint(p1);
		assert(temp_edge);
		assert(temp_edge->getFaceCount() == 2);
		if(act_edges){
			MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(act_tag);
			if(aact){
				mc.countMetricAtPoints(p0, p1);
				double len = temp_edge->getLength(mc);
				if(len != aact->len){
					if(len >= MIN_EDGE_LEN){
						temp_edge->removeTag(act_tag);
						delete act_edges->removeDataItem(aact->index);
					}else{
						aact->len = len;
						act_edges->updateDataItemPosition(aact);
					}
				}
			}else if(refill_edges && !temp_edge->isBorder() && (temp_edge->getFaceCount() == 2)){
				if(! (no_border_points && (p0->isBorder() || p1->isBorder()) )) {
					mc.countMetricAtPoints(p0, p1);
					double len = temp_edge->getLength(mc);
					if(len < MIN_EDGE_LEN){
						aact = new MeshEdge3d::ActiveEdge(temp_edge, len);
						act_edges->addDataItem(aact);
						temp_edge->setPtrTag(act_tag, aact);
					}
				}
			}
		}
	}

	if(show_case){
		ostringstream ostr;
		ostr << "After collapse-special";
//		ostr << "#" << special_counter;
		if(!mesh->isValid()) ostr << endl << "MESH INVALID!";
		SHOW_MESH_NORESET(ostr.str(), mesh->getDebugViewSetTopological(p0, 2));
	}else{
//		assert(mesh->isValid());
	}

	return true;
}

/// Collapse too short edges (inner)
int MeshGenerator3dSurface::collapseInnerEdges(
		Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		bool no_border_points, bool refill_edges,
		TagExtended::TagType tag_type, int tag_value)
{
//	START_CLOCK("MeshGenerator3dSurface::collapseInnerEdges");

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Mesh transformation: collapsing inner edges ...");
	// check ranks
	int removed_count = collapseAllInnerVerticesRankLow(mc, mesh);

	// gather edges
	DataContainer<MeshEdge3d::ActiveEdge> active_edges(mesh->getPointsCount(), true);
	for(IteratorEdge3d it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(edge->isBorder() || edge->getFaceCount() != 2) continue;
		MeshPoint3d* p0 = edge->getMeshPoint(0);
		MeshPoint3d* p1 = edge->getMeshPoint(1);
		if(no_border_points){
			if(p0->isBorder() || p1->isBorder()) continue;
		}
		MeshFace* face0 = edge->getFaceAt(0);
		MeshFace* face1 = edge->getFaceAt(1);
		assert(face0 && face1); // since count == 2

		if(tag_type != TagExtended::TAG_NONE){
			if(!face0->checkIntTag(tag_type, tag_value)) continue;
			if(!face1->checkIntTag(tag_type, tag_value)) continue;
		}
		
		{  // check if artificial extra "between-domain-edge" ?
			int sid0 = face0->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN );
			int sid1 = face1->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN );
			if( sid0 != sid1 ) {
				LOG4CPLUS_WARN(MeshLog::logger_console, "MG3d::collapseInnerEdges - extra-border-collapse-cancel");
				continue;
			}
		}

		mc.countMetricAtPoints(p0, p1);
		double len = edge->getLength(mc);
		if(len < MIN_EDGE_LEN){
			MeshEdge3d::ActiveEdge* act = new MeshEdge3d::ActiveEdge(edge, len);
			active_edges.addDataItem(act);
			edge->setPtrTag(TagExtended::TAG_COLLAPSE_3D, act);
		}
	}

	bool check_edges_len = false;

	if(check_edges_len){
		double min_len = LARGE_NUMBER;
		double max_len = 0.0;
		int min_counter = 0;
		int missing_counter = 0;
		for(auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
			MeshEdge3d* edge = it.getEdge();
			if(edge->isBorder()) continue;
			mc.countMetricAtPoints(edge->getMeshPoint(0), edge->getMeshPoint(1));
			double len = edge->getLength(mc);
			if(min_counter == 0 || len < min_len){ min_len = len; min_counter = 1; }
			else if (len == min_len) min_counter++;
			if( len > max_len ) max_len = len;
			if(len < MIN_EDGE_LEN && !edge->availableTag(TagExtended::TAG_COLLAPSE_3D))
				missing_counter++;
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
				"Min. metric inner edge length: " << min_len << "  (#" << min_counter 
				<< "), missing edges = " << missing_counter);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Max. metric inner edge length: " << max_len);
	}

//	mesh->setHeapOrder(false); // should be false anyway...

	//if( mesh->hasAnyLocalShape() )
	//	assert(mesh->isValid());
	//int counter = 0;

	bool show_case = false;

#ifdef T_DEBUG_
	double last_act_len_ratio = -1.0;
#endif // T_DEBUG_

	while(active_edges.countInt() > 0){

		//counter++;
		//if( mesh->hasAnyLocalShape() && !mesh->isValid() ) {
		//	MeshViewSet* set = mesh->getViewSet();
		//	set->addInfo("counter", to_string(counter) );
		//	SHOW_MESH("mesh invalid!", set);
		//}

		//if( counter >= 39 ) 
		//	show_case = true;

		MeshEdge3d::ActiveEdge* act = active_edges.removeDataItem(0);
		act->edge->removeTag(TagExtended::TAG_COLLAPSE_3D);

#ifdef T_DEBUG_
		double len_ratio = act->len / MIN_EDGE_LEN;
		if( len_ratio - last_act_len_ratio > 0.01 ){
			last_act_len_ratio = len_ratio;
			ostringstream ostr;
			ostr.precision(2);
			ostr << "Collapsing inner edges with ratio " << fixed << len_ratio;
			ostr << " [#" << removed_count << "]";
			LOG4CPLUS_INFO(MeshLog::logger_console,ostr.str());
		}
#endif // T_DEBUG_

		if(act->len > MIN_EDGE_LEN || act->edge->isBorder() || act->edge->getFaceCount() != 2){
			delete act;
			continue;
		}

//		show_case = (act->len > 0.0);
//		show_case = (removed_count >= 1528);

		MeshPoint3d* p[2] = { act->edge->getMeshPoint(0), act->edge->getMeshPoint(1) };

		bool p0_border = p[0]->isBorder();
		bool p1_border = p[1]->isBorder();

		if(no_border_points){
			if(p0_border || p1_border){ // none at all
				delete act;
				continue;
			}
		}else{
			if(p0_border && p1_border){ // one at most
				delete act;
				continue;
			}
		}

		MeshFace* f[2] = { act->edge->getFaceAt(0), act->edge->getFaceAt(1) };
		MeshPoint3d* removed_point = nullptr;

		if(!p0_border && p[0]->getRank() == 2) { // special case
			removed_point = p[0];
		}else if(!p1_border && p[1]->getRank() == 2) { // special case, symmetrical
			removed_point = p[1];
		}

		if(removed_point && 
			collapseInnerVertexRankTwo(mc, mesh, removed_point, no_border_points, refill_edges, 
										&active_edges, TagExtended::TAG_COLLAPSE_3D))
		{
			removed_count++;
			delete act;
			continue;
		}

		// check length of adjacent edges
		//const double VICINITY_LEN_THRESHOLD = 0.8;
		//double max_vlen = 0.0;
		//if(!p0->isBorder()){
		//	for(int i = 0; max_vlen < VICINITY_LEN_THRESHOLD && i < p0->getRank(); i++){
		//		MeshEdge3d* edge = p0->getEdge(i);
		//		if(edge == act->edge) continue;
		//		mc.countMetricAtPoints(p0, edge->getOtherPoint(p0));
		//		double len = edge->getLength(mc);
		//		if(len > max_vlen) max_vlen = len;
		//	}
		//}
		//if(!p1->isBorder()){
		//	for(int i = 0; max_vlen < VICINITY_LEN_THRESHOLD && i < p1->getRank(); i++){
		//		MeshEdge3d* edge = p1->getEdge(i);
		//		if(edge == act->edge) continue;
		//		mc.countMetricAtPoints(p1, edge->getOtherPoint(p1));
		//		double len = edge->getLength(mc);
		//		if(len > max_vlen) max_vlen = len;
		//	}
		//}
		//if(max_vlen >= VICINITY_LEN_THRESHOLD){
		//	delete act;
		//	continue;
		//}

		// for non-triangle, arbitrary from the rest ...
		MeshPoint3d* px[2] = { f[0]->getOtherPoint(p[0], p[1]), f[1]->getOtherPoint(p[0], p[1]) };

		// check extra conditions

		for(int j = 0; j < 2; j++) {
			while(!px[j]->isBorder() && (px[j]->getRank() == 3)){
				if( removeTriangulationPointSimple(mesh, px[j], &active_edges, TagExtended::TAG_COLLAPSE_3D) ){
					f[j] = act->edge->getOtherFace(f[1-j]);
					px[j] = f[j]->getOtherPoint(p[0], p[1]);
				}else break;
			}
		}

		bool valid = true;
		for(int j = 0; valid && (j < 2); j++){
			for(int i = 0; i < p[j]->getRank(); i++){
				MeshPoint3d* other_point = p[j]->getEdge(i)->getOtherPoint(p[j]);
				if(other_point != p[1-j] && other_point != px[0] && other_point != px[1] &&
					other_point->getEdgeToPoint(p[1-j]) != nullptr)
				{
					valid = false;
					break;
				}
			}
		}
		if(!valid){
			if(show_case){
				ostringstream ostr;
				ostr << "Selected edge to collapse, metric-len=" << act->len;
				ostr << ", real-len=" << p[0]->getCoordinates().distance(p[1]->getCoordinates());
				ostr << endl << "#" << removed_count;
				ostr << ": p" << p[0]->getIndex() << " - " << (p[0]->isBorder() ? "border" : "inner");
				ostr << ", p" << p[1]->getIndex() << " - " << (p[1]->isBorder() ? "border" : "inner");
				ostr << endl << "EXTRA CONDITIONS: NOT VALID";
				SHOW_MESH(ostr.str(), mesh->getDebugViewSetTopological(p[0], 2));
			}
			delete act;
			continue;
		}

		SurfaceConstPtr surface = nullptr;

		if( mesh->hasAnyLocalShape() ){
			DataVector< MeshFace* > mfaces(2);
			mfaces.add( f[0] ); mfaces.add( f[1] );
			DataVector< MeshPoint3d*> mpoints(4);
			MeshFace::getMeshPointsForFaces( mfaces, mpoints );
			surface = getLocalSurfaceForPoints( mc, mesh, mfaces, mpoints );
			if( surface == nullptr ) {
				if(true){
					ostringstream ostr;
					ostr << "Selected edge to collapse, metric-len=" << act->len;
					ostr << ", real-len=" << p[0]->getCoordinates().distance(p[1]->getCoordinates());
					ostr << endl << "#" << removed_count;
					ostr << ": p" << p[0]->getIndex() << " - " << (p[0]->isBorder() ? "border" : "inner");
					ostr << ", p" << p[1]->getIndex() << " - " << (p[1]->isBorder() ? "border" : "inner");
					ostr << endl << "NO COMMON LOCAL SURFACE";
					SHOW_MESH(ostr.str(), mesh->getDebugViewSetTopological(p[0], 2));
				}
				delete act;
				continue;
			}
		}


		removed_count++;

		//show_case = f[0]->availableTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN )
		//	&& (f[0]->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, 0 ) ||
		//		f[0]->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, 9 ) );
//		show_case = (removed_count >= 17565);
//		show_case = false;

		if(show_case){
			ostringstream ostr;
			ostr << "Selected edge to collapse, metric-len=" << act->len;
			ostr << ", real-len=" << p[0]->getCoordinates().distance(p[1]->getCoordinates());
			ostr << endl << "#" << removed_count;
			ostr << ": p" << p[0]->getIndex() << " - " << (p[0]->isBorder() ? "border" : "inner");
			ostr << ", p" << p[1]->getIndex() << " - " << (p[1]->isBorder() ? "border" : "inner");
			if(!mesh->isValid( mc.getMinSkipLen2() >= 0.0 )) ostr << endl << "MESH INVALID!";

			if(check_edges_len){
				double min_len = LARGE_NUMBER;
				int min_counter = 0;
				int missing_counter = -1;
				for(auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
					MeshEdge3d* edge = it.getEdge();
					if(edge->isBorder()) continue;
					mc.countMetricAtPoints(edge->getMeshPoint(0), edge->getMeshPoint(1));
					double len = edge->getLength(mc);
					if(min_counter == 0 || len < min_len){ min_len = len; min_counter = 1; }
					else if (len == min_len) min_counter++;
					if(len < MIN_EDGE_LEN && !edge->availableTag(TagExtended::TAG_COLLAPSE_3D))
						missing_counter++;
				}
				ostr << endl << "Min. metric inner edge length: " << min_len << "  (#" 
						<< min_counter << "), missing edges = " << missing_counter;
			}
			
//			SHOW_MESH(ostr.str(), mesh->getDebugViewSetTopological(p0, 2, TagExtended::TAG_ADAPT_SURF));
			MeshViewSet* set = mesh->getDebugViewSetTopological(f[0], 2);
			set->addInfo("sub-domain id", f[0]->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1) );
			SHOW_MESH(ostr.str(), set);
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, ostr.str());
		}

//		assert(mesh->isValid());

		bool zero_len = (act->len == 0.0);
		delete act;

		// In case of poly-faces, there may be other edges between f0 and f1 - need to remove active-flag
		bool other_common_act_edge = false;
		int f0_ect = f[0]->getEdgeCount();
		for(int i = 0; i < f0_ect; i++){
			MeshEdge3d* temp_edge = f[0]->getEdge(i);
			if(temp_edge->isBorder()) continue;
			assert(temp_edge->getFaceCount() == 2);
			if(temp_edge->getOtherFace(f[0]) == f[1]){
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					other_common_act_edge = true;
					temp_edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
					delete active_edges.removeDataItem(aact->getIndex());
				}
			}
		}

		f[0]->detachFromEdges();
		f[1]->detachFromEdges();

		DPoint3d old_pt;
		ParamAndQuality old_pt_2d( DPoint2d::zero, AQ_UNKNOWN );

		int sfct = std::max( p[0]->getRank(), p[1]->getRank() );
		DataVector<MeshFace*> switched_faces(sfct);
		DataVector<double>	  switched_quality(sfct);

		double fq[2] = { 1.0, 1.0 };
		DVector3d fn[2];
		for(int j = 0; j < 2; j++) {
			if(f[j]->getPointCount() > 3) {
				fq[j] = f[j]->getShapeQuality(mc);
			}
		}

		bool both_inner = false;

		if(p[0]->isBorder()){ // --> move towards p0
			// remove tag from edges which will be removed
			for(int i = 0; i < p[1]->getRank(); i++){
				MeshEdge3d* temp_edge = p[1]->getEdge(i);
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					temp_edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
					delete active_edges.removeDataItem(aact->getIndex());
				}
			}
			// switch all elements incident to point
			while(p[1]->getRank() > 0){ // edges will be removed automatically with no adjacent faces...
				MeshFace* face = p[1]->getEdge(0)->getFaceAt(0);
				switched_faces.add(face);
				switched_quality.add(face->getShapeQuality(mc));
				face->switchPointsWithEdges(p[1], p[0]);
			}
			// set point to remove later
			removed_point = p[1];
			if( surface != nullptr ) old_pt_2d = p[0]->getLocalSurfaceParamQuality( surface );
			else old_pt = p[0]->getCoordinates();
		}else if(p[1]->isBorder()){ // --> move towards p1
			// remove tag from edges which will be removed
			for(int i = 0; i < p[0]->getRank(); i++){
				MeshEdge3d* temp_edge = p[0]->getEdge(i);
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					temp_edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
					delete active_edges.removeDataItem(aact->getIndex());
				}
			}
			// switch all elements incident to point
			while(p[0]->getRank() > 0){ // edges will be removed automatically with no adjacent faces...
				MeshFace* face = p[0]->getEdge(0)->getFaceAt(0);
				switched_faces.add(face);
				switched_quality.add(face->getShapeQuality(mc));
				face->switchPointsWithEdges(p[0], p[1]);
			}
			removed_point = p[0];
			p[0] = p[1]; p[1] = removed_point;
			if( surface != nullptr ) old_pt_2d = p[0]->getLocalSurfaceParamQuality( surface );
			else old_pt = p[0]->getCoordinates();
		}else{ // --> move using Laplace
			both_inner = true;
			// remove tag from edges which will be removed
			for(int i = 0; i < p[1]->getRank(); i++){
				MeshEdge3d* temp_edge = p[1]->getEdge(i);
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					temp_edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
					delete active_edges.removeDataItem(aact->getIndex());
				}
			}
			// switch all elements incident to point
			while(p[1]->getRank() > 0){ // edges will be removed automatically with no adjacent faces...
				MeshFace* face = p[1]->getEdge(0)->getFaceAt(0);
				switched_faces.add(face);
				switched_quality.add(face->getShapeQuality(mc));
				face->switchPointsWithEdges(p[1], p[0]);
			}
			removed_point = p[1];
			if( surface != nullptr ) old_pt_2d = p[0]->getLocalSurfaceParamQuality( surface );
			else old_pt = p[0]->getCoordinates();

			if(p[0]->hasLocalSurface())
				movePointByLaplace(mesh, mc, p[0] );
		}

		// Check if the mesh is OK
		// ... check normal vectors for switched faces
		valid = true;
		if(!zero_len){
			sfct = switched_faces.countInt();
			for(int i = 0; valid && i < sfct; i++){
				//valid = switched_faces[i]->validForLocalVicinitySurface(mc);
				valid = ( surface != nullptr ) ? switched_faces[i]->valid(surface)
					: switched_faces[i]->validDirect( mc, switched_quality[i] );
			}

			if(!valid && both_inner){
				const double tp[] = { 2.0, -1.0, 1.0, 0.0, 0.5 };
				int tPI = sizeof(tp) / sizeof(double);
				if( surface != nullptr ) {
					const ParamAndQuality param0 = p[0]->getLocalSurfaceParamQuality( surface );
					const DPoint2d param1 = p[1]->getLocalSurfaceParam( surface );
					do{
						valid = true;
						tPI--;
						const DPoint2d temp_pt(old_pt_2d.param, param1, tp[tPI]);
						p[0]->setCoordinates( surface, DPoint2d(old_pt_2d.param, param1, tp[tPI]), false );

						for(int i = 0; valid && i < sfct; i++){
							//valid = switched_faces[i]->validForLocalVicinitySurface(mc);
							valid = switched_faces[i]->valid(surface);
						}			
					}while(!valid && (tPI > 0));

					if( valid ) p[0]->setCoordinates( surface, DPoint2d(old_pt_2d.param, param1, tp[tPI]), AQ_UNKNOWN );
					else p[0]->setCoordinates( surface, param0.param, param0.quality );
				}else{
					const DPoint3d dpt0 = old_pt;
					do{
						valid = true;
						tPI--;
						const DPoint3d temp_pt(old_pt, p[1]->getCoordinates(), tp[tPI]);
						p[0]->setCoordinates( temp_pt );

						for(int i = 0; valid && i < sfct; i++){
							//valid = switched_faces[i]->validForLocalVicinitySurface(mc);
							valid = switched_faces[i]->validDirect(mc, switched_quality[i]);
						}			
					}while( !valid && (tPI > 0) );
					if( !valid ) p[0]->setCoordinates( old_pt );
				}
			}
		}

		if(valid){
			for(int j = 0; j < 2; j++){
				MeshFace* reduced_face = f[j]->removePoint(removed_point);
				// triangle->nullptr, quad->triangle, poly(5)->quad, poly(n)->poly(n-1)
				if(reduced_face != f[j]){
					delete mesh->removeMeshFace(f[j]);
					if(reduced_face != nullptr) mesh->addMeshFace(reduced_face);
					f[j] = reduced_face;
				}
			}

			delete mesh->removeMeshPoint(removed_point);

			DataVector<MeshEdge3d*> edges_to_check(50);
			int p0_rank = p[0]->getRank();
			for(int i = 0; i < p0_rank; i++)
				edges_to_check.add(p[0]->getEdge(i));
			if(other_common_act_edge && f[0] && f[1]){
				// both f0 and f1 are poly-faces, modified but not replaced
				f0_ect = f[0]->getEdgeCount();
				for(int i = 0; i < f0_ect; i++){
					MeshEdge3d* temp_edge = f[0]->getEdge(i);
					if(temp_edge->isBorder()) continue;
					assert(temp_edge->getFaceCount() == 2);
					if(temp_edge->getOtherFace(f[0]) == f[1])
						edges_to_check.addIfNew( temp_edge );
				}
			}

			bool fmodified[2] = { false, false };
			for(int j = 0; j < 2; j++){
				if(f[j]){ // i.e. if it was not a triangle
					if( surface ? f[j]->valid(surface) : f[j]->validDirect(mc, fq[j]) ) continue;
					//// 1. Laplace - if all neighbours with similar normal vectors
					//improveFace(mc, f[j]);
					//if( (fmodified[j] = f[j]->valid(mc, fn[j], fq[j]) ) ) continue;
					//// maybe try again ?
					//improveFace(mc, f[j]);
					//if( (fmodified[j] = f[j]->valid(mc, fn[j], fq[j]) ) ) continue;
					//LOG4CPLUS_INFO(MeshLog::logger_console, "Failed improving f[0/1] after collapse", j);
					// 2. try join?

					// 3. try join vertices -> prepare for collapse ? ...
					int fect = f[j]->getEdgeCount();
					MeshEdge3d* shortest_edge = f[j]->getEdge(0);
					double min_len = shortest_edge->getLength(mc);
					for(int k = 1; k < fect; k++){
						MeshEdge3d* edge = f[j]->getEdge(k);
						double len = edge->getLength(mc);
						if(len < min_len){
							min_len = len;
							shortest_edge = edge;
						}
					}
					DPoint3d middle = shortest_edge->getPoint(0.5);
					MeshPoint3d* se_p0 = shortest_edge->getMeshPoint(0);
					MeshPoint3d* se_p1 = shortest_edge->getMeshPoint(1);
					//LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Step #" << removed_count << " Fixing f" << j << " -> moving together points " 
					//	<< se_p0->getIndex() << " and " << se_p1->getIndex() << endl;
					if( surface ) {
						DPoint2d param;
						double pquality = surface->withinDomainQuality( middle, param );
						assert( pquality >= 0.0 );
						middle = surface->getPoint( param );
						se_p0->setCoordinates( surface, param, pquality, &middle );
						se_p1->setCoordinates( surface, param, pquality, &middle );
					}else{
						se_p0->setCoordinates( middle );
						se_p1->setCoordinates( middle );
					}
					fmodified[j] = true;
				}
			}

//			show_case = (fmodified[0] || fmodified[1]);

			for(int j = 0; j < 2; j++){
				if(fmodified[j]) {
					int fpct = f[j]->getPointCount();
					for(int k = 0; k < fpct; k++){
						MeshPoint3d* point = f[j]->getPoint(k);
						int rank = point->getRank();
						for(int l = 0; l < rank; l++){
							edges_to_check.addIfNew( point->getEdge(l) );
						}
					}
				}
			}

			int echeck_count = edges_to_check.countInt();
			// update active edges incident to p0
			for(int i = 0; i < echeck_count; i++){
				MeshEdge3d* edge = edges_to_check[i];
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)edge->getPtrTag(TagExtended::TAG_COLLAPSE_3D);
				if(aact){
					mc.countMetricAtPoints(p[0], edge->getOtherPoint(p[0]));
					double len = edge->getLength(mc);
					if(len != aact->len){
						if(len >= MIN_EDGE_LEN){
							edge->removeTag(TagExtended::TAG_COLLAPSE_3D);
							delete active_edges.removeDataItem(aact->index);
						}else{
							aact->len = len;
							active_edges.updateDataItemPosition(aact);
						}
					}
				}else{
					if(!refill_edges) continue;
					if(edge->isBorder() || edge->getFaceCount() != 2)
						continue;
					MeshPoint3d* ep0 = edge->getOtherPoint(p[0]);
					if(no_border_points && (p[0]->isBorder() || ep0->isBorder()))
						continue;

					MeshFace* face0 = edge->getFaceAt(0);
					MeshFace* face1 = edge->getFaceAt(1);
					assert(face0 && face1); // since count == 2

					if(tag_type != TagExtended::TAG_NONE){
						if(!face0->checkIntTag(tag_type, tag_value)) continue;
						if(!face1->checkIntTag(tag_type, tag_value)) continue;
					}

					mc.countMetricAtPoints(p[0], ep0);
					double len = edge->getLength(mc);
					if(len < MIN_EDGE_LEN){
						aact = new MeshEdge3d::ActiveEdge(edge, len);
						active_edges.addDataItem(aact);
						edge->setPtrTag(TagExtended::TAG_COLLAPSE_3D, aact);
					}
				}
			}

			if(show_case){
				ostringstream ostr;
				ostr << "After collapse #" << removed_count;
				if(!mesh->isValid()) ostr << endl << "MESH INVALID!";
				SHOW_MESH_NORESET(ostr.str(), mesh->getDebugViewSetTopological(p[0], 2, TagExtended::TAG_ADAPT_SURF));
			}else{
//				assert(mesh->isValid());
			}

			DataVector<MeshPoint3d*> vlow; // list of vertices with rank two
			if((p0_rank == 2) || (p0_rank == 3)) vlow.add(p[0]);
			for(int i = 0; i < p0_rank; i++){
				MeshPoint3d* point = p[0]->getEdge(i)->getOtherPoint(p[0]);
				if(!point->isBorder() && ( (point->getRank() == 2) || (point->getRank() == 3)))
					vlow.add(point);
			}

			if(vlow.notEmpty()){
				for(int i = 0; i < vlow.countInt(); i++){
					MeshPoint3d* rpoint = vlow[i];
					int rank = rpoint->getRank();
					DataVector<MeshPoint3d*> rp(3);
					for(int j = 0; j < rank; j++)
						rp.add(rpoint->getEdge(j)->getOtherPoint(rpoint));
					bool result = false;
					if(rank == 2)
						result = collapseInnerVertexRankTwo(mc, mesh, rpoint, 
									no_border_points, refill_edges, &active_edges, TagExtended::TAG_COLLAPSE_3D);
					else if(rank == 3)
						result = removeTriangulationPointSimple(mesh, rpoint, &active_edges, TagExtended::TAG_COLLAPSE_3D);
					if(result){
						removed_count++;
						for(int j = 0; j < rank; j++)
							if( !rp[j]->isBorder() && ( (rp[j]->getRank() == 2) || rp[j]->getRank() == 3)) 
								vlow.addIfNew(rp[j]);
					}
				}
			}
			//if(true){
			//	// check ...
			//	int pct = mesh->getPointsCount();
			//	vtwo.clear();
			//	for(int i = 0; i < pct; i++){
			//		MeshPoint3d* point = mesh->getPointAt(i);
			//		if(!point->isBorder() && (point->getRank() == 2))
			//			vtwo.add(point);
			//	}
			//	assert(vtwo.empty());
			//}

		}else{
			if(show_case){
//			if(true){
				MeshViewSet* set = new MeshViewSet();
				DataVector<MeshPoint3d*> spoints(100);
				spoints.add(p[0]);
				for(int i = 0; i < sfct; i++){
					MeshFace* sf = switched_faces[i];
					set->addFaceWithEdges(sf);
					int sfpct = sf->getPointCount();
					for(int j = 0; j < sfpct; j++)
						spoints.addIfNew(sf->getPoint(j));
				}
				set->addPoint(p[0], 1);
				for(int i = 0; i < spoints.countInt(); i++)
					set->addPoint(spoints[i]);
				ostringstream ostr;
				ostr << "Invalid switched faces #" << removed_count;
				SHOW_MESH(ostr.str(), set);
			}
			//SHOW_MESH("Error after collapse?", mesh->getDebugViewSet(p0, nullptr, 2.0, TagExtended::TAG_ADAPT_SURF));
			// ... restore 
			if(!p[0]->isBorder()) {
				if( surface != nullptr) p[0]->setCoordinates( surface, old_pt_2d.param, old_pt_2d.quality );
				else p[0]->setCoordinates( old_pt );
			}
			for(int i = 0; i < sfct; i++){
				switched_faces[i]->switchPointsWithEdges(p[0], removed_point);
			}
			f[0]->attachToEdges();
			f[1]->attachToEdges();
			//SHOW_MESH("Restored after collapse?", mesh->getDebugViewSet(p0, nullptr, 2.0, TagExtended::TAG_ADAPT_SURF));
//			if(true){
			if(show_case){
				ostringstream ostr;
				ostr << "Restored after collapse #" << removed_count;
				if(!mesh->isValid()) ostr << endl << "MESH INVALID!";
				SHOW_MESH_NORESET(ostr.str(), mesh->getDebugViewSetTopological(p[0], 2, TagExtended::TAG_ADAPT_SURF));
			}else{
//				assert(mesh->isValid());
			}
		}

	}

	if(check_edges_len){
		double min_len = LARGE_NUMBER;
		double max_len = 0.0;
		int min_counter = 0;
		int missing_counter = 0;
		for(auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
			MeshEdge3d* edge = it.getEdge();
			if(edge->isBorder()) continue;
			mc.countMetricAtPoints(edge->getMeshPoint(0), edge->getMeshPoint(1));
			double len = edge->getLength(mc);
			if(min_counter == 0 || len < min_len){ min_len = len; min_counter = 1; }
			else if (len == min_len) min_counter++;
			if( len > max_len ) max_len = len;
			if(len < MIN_EDGE_LEN && !edge->availableTag(TagExtended::TAG_COLLAPSE_3D))
				missing_counter++;
		}
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
			"*Min. metric inner edge length: " << min_len << "  (#" 
				<< min_counter << "), missing edges = " << missing_counter);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "*Max. metric inner edge length: " << max_len);
	}
		
	//	STOP_CLOCK("MeshGenerator3dSurface::collapseInnerEdges");
	return removed_count;
}

/// Collapse too short edges (boundary)
int MeshGenerator3dSurface::collapseBoundaryEdges(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		bool no_corner_points, bool refill_edges,
		TagExtended::TagType tag_type, int tag_value1, int tag_value2)
{
//	START_CLOCK("MeshGenerator3dSurface::collapseBoundaryEdges");

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Mesh transformation: collapsing boundary edges ...");
	// gather edges
	DataContainer<MeshEdge3d::ActiveEdge> active_edges(mesh->getPointsCount(), true);
	for(IteratorEdge3d it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if( !edge->isBorder() ) continue;
		if(tag_type != TagExtended::TAG_NONE){
			bool all_with_tags = true;
			for(int i = 0; all_with_tags && (i < edge->getFaceCount()); i++){
				MeshFace* face = edge->getFaceAt(i);
				all_with_tags = face->checkIntTag(tag_type, tag_value1) ||
					face->checkIntTag(tag_type, tag_value2);
			}
			if(!all_with_tags) continue;
		}

		mc.countMetricAtPoints(edge->getMeshPoint(0), edge->getMeshPoint(1));
		double len = edge->getLength(mc);
		if(len < MIN_BEDGE_LEN){
			MeshEdge3d::ActiveEdge* act = new MeshEdge3d::ActiveEdge(edge, len);
			active_edges.addDataItem(act);
			edge->setPtrTag(TagExtended::TAG_COLLAPSE_3D, act);
		}
	}

//	assert(mesh->isValid());

	int removed_count = 0;
//	int counter = 0;

	while(active_edges.countInt() > 0){
//		counter++;
		MeshEdge3d::ActiveEdge* act = active_edges.removeDataItem(0);
		act->edge->removeTag(TagExtended::TAG_COLLAPSE_3D);

		MeshEdge3d* act_edge = act->edge;
		double act_len = act->len;
		delete act;

		if(act_len < MIN_BEDGE_LEN)
			if(collapseBoundaryEdge(mc, mesh, act_edge, true, no_corner_points, refill_edges, &active_edges, TagExtended::TAG_COLLAPSE_3D))
				removed_count++;

		// debug check
//		assert(mesh->isValid());
		//for(int i = 0; i < active_edges.countInt(); i++){
		//	MeshEdge3d::ActiveEdge* act = active_edges.getDataAt(i);
		//	bool ok = act->edge->availableTag(TagExtended::TAG_COLLAPSE_3D);
		//}
	}
//	STOP_CLOCK("MeshGenerator3dSurface::collapseBoundaryEdges");
	return removed_count;
}

bool MeshGenerator3dSurface::collapseBoundaryEdge(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
					MeshEdge3d* col_edge, bool check_length, 
					bool no_corner_points, bool refill_edges,
					DataContainer<MeshEdge3d::ActiveEdge> * act_edges, TagExtended::TagType act_tag,
					TagExtended::TagType tag_type, int tag_value1, int tag_value2)
{
	if( ! col_edge->isBorder() ) return false;

	MeshPoint3d* p0 = col_edge->getMeshPoint(0);
	MeshPoint3d* p1 = col_edge->getMeshPoint(1);
	MeshEdge3d* p0_edge = nullptr;
	MeshEdge3d* p1_edge = nullptr;
	MeshPoint3d* p00 = nullptr;
	MeshPoint3d* p11 = nullptr;

	bool p0_corner = p0->isBorder(TagBorder::CORNER) || (p0->getBorderEdgesCount() != 2);
	bool p1_corner = p1->isBorder(TagBorder::CORNER) || (p1->getBorderEdgesCount() != 2);

	if(no_corner_points){
		if(p0_corner || p1_corner) // none at all
			return false;
	}else{
		if(p0_corner && p1_corner) // one at most
			return false;
	}

	int efct = col_edge->getFaceCount();
	DataVector<MeshFace*> c_faces(efct);
	DataVector<MeshPoint3d*> cf_points(efct);
	DataHashTableKeyValue< int, int > htags(10, -1);
	DataVector< DataVector<MeshFace*> > cf_tag_faces(10);

	for(int i = 0; i < efct; i++){
		MeshFace* face = col_edge->getFaceAt(i);
		c_faces.add(face);
		cf_points.add(face->getOtherPoint(p0, p1));
		int ftag = face->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );
		unsigned int hi;
		int cfi;
		if(htags.contains( ftag, hi) ) {
			cfi = htags.slotValue(hi);
		}else{
			cfi = (int)cf_tag_faces.add( DataVector< MeshFace*>() );
			htags.insert( ftag, cfi );
		}
		cf_tag_faces[cfi].add( face );
	}

	// check extra conditions
	for(int i = 0; i < p0->getRank(); i++){
		MeshPoint3d* other_point = p0->getEdge(i)->getOtherPoint(p0);
		if(other_point == p1) continue;
		if(other_point->getEdgeToPoint(p1) == nullptr) continue;
		bool no_cf_point = true;
		for(int j = 0; no_cf_point && (j < efct); j++)
			no_cf_point = (other_point != cf_points[j]);
		if( no_cf_point ) return false;
	}
	for(int i = 0; i < p1->getRank(); i++){
		MeshPoint3d* other_point = p1->getEdge(i)->getOtherPoint(p1);
		if(other_point == p0) continue;
		if(other_point->getEdgeToPoint(p0) == nullptr) continue;
		bool no_cf_point = true;
		for(int j = 0; no_cf_point && (j < efct); j++)
			no_cf_point = (other_point != cf_points[j]);
		if( no_cf_point ) return false;
	}

	//for(int i = 0; i < p1->getRank(); i++){
	//	MeshPoint3d* other_point = p1->getEdge(i)->getOtherPoint(p1);
	//	if(other_point != p0 && other_point != p2 && other_point != p3 &&
	//		other_point->getEdgeToPoint(p0) != nullptr)
	//	{
	//		valid = false;
	//		break;
	//	}
	//}
	//if(!valid) continue;

	if(check_length){
		// *** check lengths
		mc.countMetricAtPoints( p0, p1 );
		double col_len = col_edge->getLength(mc);
		if( !p0_corner ){
			for(int i = 0; i < p0->getRank(); i++){
				MeshEdge3d* e = p0->getEdge(i);
				if(e == col_edge || !e->isBorder()) continue;
				p0_edge = e;
				mc.countMetricAtPoints(p0, e->getOtherPoint(p0));
				double len = e->getLength(mc);
				if( p1_corner && (col_len + len >= MAX_BEDGE_LEN_JOINED) ) return false;
			}
		}
		if( !p1_corner ){
			for(int i = 0; i < p1->getRank(); i++){
				MeshEdge3d* e = p1->getEdge(i);
				if(e == col_edge || !e->isBorder()) continue;
				p1_edge = e;
				mc.countMetricAtPoints(p1, e->getOtherPoint(p1));
				double len = e->getLength(mc);
				if( p0_corner && (col_len + len >= MAX_BEDGE_LEN_JOINED) ) return false;
			}
		}
	}

	char bf = col_edge->getBorderFlags();

	if( !p0_corner ) {
		for(int i = 0; (p0_edge == nullptr) && (i < p0->getRank()); i++){
			MeshEdge3d* e = p0->getEdge(i);
			if( e != col_edge && e->isBorder() ) p0_edge = e;
		}
		p00 = p0_edge->getOtherPoint( p0 );
	}
	if( !p1_corner ) {
		for(int i = 0; (p1_edge == nullptr) && (i < p1->getRank()); i++){
			MeshEdge3d* e = p1->getEdge(i);
			if( e != col_edge && e->isBorder() ) p1_edge = e;
		}
		p11 = p1_edge->getOtherPoint( p1 );
	}

	auto col_shape = col_edge->getLocalCurve();
	double t0 = p0->getLocalCurveParam( col_shape );
	double t1 = p1->getLocalCurveParam( col_shape );

	bool show_case = false;
//	bool show_case = (p1->getIndex() == 226);

	if(show_case){
		MeshViewSet* set = new MeshViewSet;
		set->addInfo("p0-type", p0_corner ? "corner" : "ridge");
		set->addInfo("p1-type", p1_corner ? "corner" : "ridge");
		set->addInfo("p0-p1", to_string(p0->getIndex())+"-"+to_string(p1->getIndex()));
		set->addInfo("t0|t1", to_string(t0)+" | "+to_string(t1));
		DataVector< MeshFace* > pfaces;
		p0->adjacentFaces( pfaces );
		p1->adjacentFaces( pfaces );
		set->addEdge( col_edge, 2 );
		set->addPoint( p0, 2 );
		set->addPoint( p1, 2 );
		for(int i = 0; i < pfaces.countInt(); i++)
			set->addFaceWithEdgesAndPoints( pfaces[i], 
				pfaces[i]->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 ) );
		SHOW_MESH("Selected bedge to collapse", set);
	}

	col_edge->clearBorder();
	for(int i = 0; i < efct; i++)
		c_faces[i]->detachFromEdges();

	MeshPoint3d* removed_point = nullptr;
	DataVector<MeshFace*> switched_faces(10);
	MeshPoint3d* other_bp = nullptr;

	if(p0_corner){ // --> move towards p0
		assert( !p1_corner);
		other_bp = p11;
		// remove tag from edges which will be removed
		for(int i = 0; i < p1->getRank(); i++){
			MeshEdge3d* temp_edge = p1->getEdge(i);
			if(temp_edge->isBorder())
				temp_edge->clearBorder();
			if(act_edges){
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(act_tag);
				if(aact){
					temp_edge->removeTag(act_tag);
					delete act_edges->removeDataItem(aact->getIndex());
				}
			}
		}
		// switch all elements incident to point
		p1->clearBorder();
		while(p1->getRank() > 0){ // edges will be removed automatically with no adjacent faces...
			MeshFace* f = p1->getEdge(0)->getFaceAt(0);
			switched_faces.add(f);
			f->switchPointsWithEdges(p1, p0);
		}
		removed_point = p1;
		MeshEdge3d* temp_edge = p0->getEdgeToPoint(other_bp);
		assert(temp_edge);
		temp_edge->setBorder(bf);
		temp_edge->setLocalCurve(col_shape);
	}else if(p1_corner){ // --> move towards p1
		assert( !p0_corner);
		other_bp = p00;
		// remove tag from edges which will be removed
		for(int i = 0; i < p0->getRank(); i++){
			MeshEdge3d* temp_edge = p0->getEdge(i);
			if(temp_edge->isBorder())
				temp_edge->clearBorder();
			if(act_edges){
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(act_tag);
				if(aact){
					temp_edge->removeTag(act_tag);
					delete act_edges->removeDataItem(aact->getIndex());
				}
			}
		}
		// switch all elements incident to point
		p0->clearBorder();
		while(p0->getRank() > 0){ // edges will be removed automatically with no adjacent faces...
			MeshFace* f = p0->getEdge(0)->getFaceAt(0);
			switched_faces.add(f);
			f->switchPointsWithEdges(p0, p1);
		}
		removed_point = p0;
		MeshEdge3d* temp_edge = p1->getEdgeToPoint(other_bp);
		assert(temp_edge);
		temp_edge->setBorder(bf);
		temp_edge->setLocalCurve(col_shape);
		p0 = p1;
		t0 = t1;
	}else{ // --> move using Laplace
		other_bp = p11;
		// remove tag from edges which will be removed
		for(int i = 0; i < p1->getRank(); i++){
			MeshEdge3d* temp_edge = p1->getEdge(i);
			if(temp_edge->isBorder())
				temp_edge->clearBorder();
			if(act_edges){
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)temp_edge->getPtrTag(act_tag);
				if(aact){
					temp_edge->removeTag(act_tag);
					delete act_edges->removeDataItem(aact->getIndex());
				}
			}
		}
		// switch all elements incident to point
		p1->clearBorder();
		while(p1->getRank() > 0){ // edges will be removed automatically with no adjacent faces...
			MeshFace* f = p1->getEdge(0)->getFaceAt(0);
			switched_faces.add(f);
			f->switchPointsWithEdges(p1, p0);
		}
		removed_point = p1;
		MeshEdge3d* temp_edge = p0->getEdgeToPoint(other_bp);
		assert(temp_edge);
		temp_edge->setBorder(bf);
		temp_edge->setLocalCurve(col_shape);
		if(show_case){
			MeshViewSet* set = new MeshViewSet;
			DataVector< MeshFace* > pfaces;
			p0->adjacentFaces( pfaces );
			set->addPoint( p0, 2 );
			for(int i = 0; i < pfaces.countInt(); i++) {
				set->addFaceWithEdgesAndPoints( pfaces[i], 
					pfaces[i]->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 ) );
				DPoint3d pmid = pfaces[i]->getMiddlePoint();
				SurfaceConstPtr psurface = pfaces[i]->getOptLocalSurface();
				DVector3d fn[4] = {
					pfaces[i]->getLocalSurfaceNormal(),
					psurface->getNormalVector( psurface->getParameters( pmid ) ),
					pfaces[i]->getNormalVector(),
					pfaces[i]->getBaseNormal() };
				double ave_len = (pfaces[i]->getEdge(0)->getLength() + pfaces[i]->getEdge(1)->getLength()
					+ pfaces[i]->getEdge(2)->getLength() ) / 3.0;
				for(int j = 0; j < 4; j++)
					set->addEdge( pmid, pmid + fn[j]*ave_len, j );
			}
			SHOW_MESH_NORESET("after collapse, before boundary-laplace", set);
		}
		moveBoundaryPointByLaplace(mesh, mc, p0);
	}

	// Check if the mesh is OK
	bool valid = true;
	// ... check if the boundary edges after collapse are not too large
	if(check_length && other_bp){
		MeshEdge3d* e = p0->getEdgeToPoint(other_bp);
		if(e){
			mc.countMetricAtPoints(p0, other_bp);
			double len = e->getLength(mc);
			if(len > MAX_BEDGE_LEN_JOINED) valid = false;
		}
	}

	if( valid ){ 
		for(int i = 0; valid && (i < cf_tag_faces.countInt()); i++){
			const DataVector<MeshFace*> & dmfaces = cf_tag_faces[i];
			DataVector<MeshPoint3d*> dmpoints;
			MeshFace::getMeshPointsForFaces( dmfaces, dmpoints );
			SurfaceConstPtr surface = 
				getLocalSurfaceForPoints( mc, mesh, dmfaces, dmpoints);
			for(int j = 0; valid && (j < dmfaces.countInt()); j++ ) {
				valid = dmfaces[j]->valid( surface );
			}
		}
	}
	// ... check normal vectors for switched faces
	//for(int i = 0; valid && i < sfct; i++){
	//	valid = switched_faces[i]->validDirect(mc);
	//}
	if(valid){
		for(int i = 0; i < efct; i++)
			delete mesh->removeMeshFace(c_faces[i]);
		delete mesh->removeMeshPoint(removed_point);
		if(other_bp) moveBoundaryPointByLaplace(mesh, mc, other_bp);

		if(act_edges){
			// update active edges incident to p0
			for(int i = 0; i < p0->getRank(); i++){
				MeshEdge3d* edge = p0->getEdge(i);
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)edge->getPtrTag(act_tag);
				if(aact){
					mc.countMetricAtPoints(p0, edge->getOtherPoint(p0));
					double len = edge->getLength(mc);
					if(len != aact->len){
						if(len >= MIN_BEDGE_LEN){
							edge->removeTag(act_tag);
							delete act_edges->removeDataItem(aact->index);
						}else{
							aact->len = len;
							act_edges->updateDataItemPosition(aact);
						}
					}
				}else{
					if(!refill_edges) continue;
					if(!edge->isBorder()) continue;
					MeshPoint3d* ep0 = edge->getMeshPoint(0);
					MeshPoint3d* ep1 = edge->getMeshPoint(1);

					bool ep0_corner = ep0->isBorder(TagBorder::CORNER) || (ep0->getBorderEdgesCount() != 2);
					bool ep1_corner = ep1->isBorder(TagBorder::CORNER) || (ep1->getBorderEdgesCount() != 2);

					if(no_corner_points){
						if( ep0_corner || ep1_corner) continue;
					}else{
						if( ep0_corner && ep1_corner) continue;
					}

					if(tag_type != TagExtended::TAG_NONE){
						bool all_with_tags = true;
						for(int i = 0; all_with_tags && (i < edge->getFaceCount()); i++){
							MeshFace* face = edge->getFaceAt(i);
							all_with_tags = face->checkIntTag(tag_type, tag_value1) ||
								face->checkIntTag(tag_type, tag_value2);
						}
						if(!all_with_tags) continue;
					}

					mc.countMetricAtPoints(ep0, ep1);
					double len = edge->getLength(mc);
					if(len < MIN_BEDGE_LEN){
						aact = new MeshEdge3d::ActiveEdge(edge, len);
						act_edges->addDataItem(aact);
						edge->setPtrTag(act_tag, aact);
					}
				}
			}
		}

		if(show_case){
			SHOW_MESH_NORESET("After collapse", 
				mesh->getDebugViewSetTopological(p0, 2, TagExtended::TAG_ADAPT_SURF));
		}

		return true;
	}else{
		if(!p0_corner && !p1_corner) {
			assert( col_shape != nullptr );
			p0->setCoordinates( col_shape, t0 );
		}
		MeshEdge3d* e = p0->getEdgeToPoint(other_bp);
		if(e) e->clearBorder();
		for(int i = 0; i < switched_faces.countInt(); i++){
			switched_faces[i]->switchPointsWithEdges(p0, removed_point);
		}
		removed_point->setBorder();
		for(int i = 0; i < efct; i++)
			c_faces[i]->attachToEdges();
		e = removed_point->getEdgeToPoint(other_bp);
		if(e){
			e->setBorder(bf);
			e->setLocalCurve(col_shape);
		}
		e = p0->getEdgeToPoint(removed_point);
		if(e){
			e->setBorder(bf);
			e->setLocalCurve(col_shape);
		}
		if(show_case){
			SHOW_MESH_NORESET("Restored after collapse", 
				mesh->getDebugViewSetTopological(p0, 2, TagExtended::TAG_ADAPT_SURF));
		}
	}

	return false;
}

/// Test Laplace methods
bool MeshGenerator3dSurface::testLaplace(Metric3dContext& mc,  MeshContainer3dSurface* mesh, MeshPoint3d* point)
{
	// calculate coordinates for
	// - simple Laplace
	// - variable Laplace
	// - mass center of triangles
	// - ...

	const DPoint3d oryg_coord = point->getCoordinates();

	bool simple_res = movePointByLaplace( mesh, mc, point );
	const DPoint3d simple_coord = point->getCoordinates();

	point->setCoordinates( oryg_coord );
	bool var_res = movePointByLaplaceForVariableMetric(mesh, mc, point );
	const DPoint3d var_coord = point->getCoordinates();

	point->setCoordinates( oryg_coord );
	bool mass_res = movePointByMassCenter(mesh, mc, point );
	const DPoint3d mass_coord = point->getCoordinates();

	point->setCoordinates( oryg_coord );
	double mg = mc.getMetricGradation( oryg_coord );

	static int test_counter = 0;
	test_counter++;

	// generate info
	if( test_counter >= 397 ){		// 199/397/515/605/669

		MeshViewSet* set = new MeshViewSet;
		set->addPoint( point );
		DataVector<MeshFace*> faces(point->getRank());
		if( point->adjacentFaces( faces )) {
			for(int i = 0; i < faces.countInt(); i++)
				set->addFaceWithEdgesAndPoints( faces[i] );
		}
		if(simple_res) set->addPoint( simple_coord, 1, 1 );
		if(var_res) set->addPoint( var_coord, 2, 2 );
		if(mass_res) set->addPoint( mass_coord, 3, 3 );

		DataStatistics stats[4];

		for(int i = 0; i < point->getRank(); i++) {
			MeshPoint3d* other_point = point->getEdge(i)->getOtherPoint(point);
			mc.countMetricAtPoints( point, other_point );
			const DPoint3d& other_coord = other_point->getCoordinates();
			stats[0].add( mc.transformRStoMS( oryg_coord - other_coord ).length() );
			if(simple_res) stats[1].add( mc.transformRStoMS( simple_coord - other_coord ).length() );
			if(var_res) stats[2].add( mc.transformRStoMS( var_coord - other_coord ).length() );
			if(mass_res) stats[3].add( mc.transformRStoMS( mass_coord - other_coord ).length() );
		}

		if( stats[0].calculate() && stats[1].calculate() && stats[2].calculate() && stats[3].calculate() ){
			set->addInfo("metric edge lengths", "oryg[" + to_string(stats[0].countInt()) 
									+ "] | simple[" + to_string(stats[1].countInt())
									+ "] | var[" + to_string(stats[2].countInt())
									+ "] | mass[" + to_string(stats[3].countInt()) );
			double vmin[4] = { stats[0].minimum(), stats[1].minimum(), stats[2].minimum(), stats[3].minimum()  };
			double vmin_max = max( max( vmin[0], vmin[1] ), max( vmin[2], vmin[3] ) );
			set->addInfo("minimum",			  to_string(vmin[0]) + ( (vmin[0] == vmin_max) ? "*":" ")
									+ " | " + to_string(vmin[1]) + ( (vmin[1] == vmin_max) ? "*":" ")
									+ " | " + to_string(vmin[2]) + ( (vmin[2] == vmin_max) ? "*":" ")
									+ " | " + to_string(vmin[3]) + ( (vmin[3] == vmin_max) ? "*":" ") );
			double vmax[4] = { stats[0].maximum(), stats[1].maximum(), stats[2].maximum(), stats[3].maximum() };
			double vmax_min = min( min( vmax[0], vmax[1] ), min( vmax[2], vmax[3] ) );
			set->addInfo("maximum",			  to_string(vmax[0]) + ( (vmax[0] == vmax_min) ? "*":" ")
									+ " | " + to_string(vmax[1]) + ( (vmax[1] == vmax_min) ? "*":" ")
									+ " | " + to_string(vmax[2]) + ( (vmax[2] == vmax_min) ? "*":" ")
									+ " | " + to_string(vmax[3]) + ( (vmax[3] == vmax_min) ? "*":" ") );
			set->addInfo("average", to_string(stats[0].average()) + " | " + to_string(stats[1].average())
									+ " | " + to_string(stats[2].average()) + " | " + to_string(stats[3].average()) );
		}
		set->addInfo("metric grad. ", mg);

		// show on one image
		SHOW_MESH(" test Laplace", set);
	}

	return true;
}

bool MeshGenerator3dSurface::movePointByLaplace(MeshContainer3dSurface* mesh, 
		Metric3dContext& mc, MeshPoint3d *point, SurfaceConstPtr surface )
{
	assert( mesh != nullptr );
	if(point->isBorder()) return false;

	int rank = point->getRank();
	DataVector< MeshFace* > faces(rank);
	if( ! point->adjacentFaces( faces ) ) return false;

	if( surface == nullptr ) {
		surface = getLocalSurfaceForPoints( mc, mesh, faces, point->localVicinityPoints() );
		if( surface == nullptr ) return false;
	}

//	mc.countMetricAtPoint(point->getCoordinates());

	//static int counter = 0;
	//counter++;

	double w = 1.0 / rank;
	DPoint2d new_pt;
	for(int j = 0; j < rank; j++){
		//MeshPoint3d* other_point = point->incidentPoint(j);
		//const DPoint2d other_param = other_point->getLocalSurfaceParam( surface );
		//const DPoint3d pt0 = other_point->getCoordinates();
		//const DPoint2d test_param = surface->getParameters( pt0 ); // TODO -> has x inverted... compared to other_param...
		//new_pt.add( other_param, w );
		new_pt.add( point->incidentPoint(j)->getLocalSurfaceParam( surface ), w );
	}

	//if( counter == 5779 ) {
	//	MeshViewSet* set = new MeshViewSet;
	//	for(int i = 0; i < rank; i++)
	//		set->addFaceWithEdgesAndPoints( faces[i] );
	//	set->addLabel( surface->getPoint( new_pt ), "newpt" );
	//	SHOW_MESH("movePointByLaplace", set);
	//}

	return point->tryMovingPoint( surface, new_pt, faces );
}

bool MeshGenerator3dSurface::movePointByMassCenter(MeshContainer3dSurface* mesh, 
		Metric3dContext& mc, MeshPoint3d *point)
{
	if(point->isBorder()) return false;

	DPoint2d new_pt;
	int rank = point->getRank();
	DataVector< MeshFace* > faces(rank);
	if( ! point->adjacentFaces( faces ) ) return false;

	SurfaceConstPtr surface = getLocalSurfaceForPoints( mc, mesh, faces, point->localVicinityPoints() );
	if( surface == nullptr ) return false;

	DataHashTableKeyValue<MeshFace*, double> fareas( rank*2, nullptr );
	DPoint3d fmid;
	for(int i = 0; i < faces.countInt(); i++) {
		fmid = faces[i]->getMiddlePoint();
		mc.countMetricAtPoint( fmid );
		if( faces[i]->getType() != FACE_TRIANGLE ) return false;
		fareas.setValue( faces[i], ((MeshTriangle3d*)faces[i])->area( mc ) );
	}

	double total_mass = 0.0;
	for(int j = 0; j < rank; j++){
		MeshEdge3d* edge = point->getEdge(j);
		int efct = edge->getFaceCount();
		double emass = 0.0;
		for(int k = 0; k < efct; k++ )
			emass += fareas.getValue( edge->getFaceAt(k), 0.0 );
		emass /= efct;

		new_pt.add( edge->getOtherPoint(point)->getLocalSurfaceParam( surface ), emass);
		total_mass += emass;
	}

	new_pt /= total_mass;
	return point->tryMovingPoint( surface, new_pt, faces );
}

bool MeshGenerator3dSurface::movePointByLaplaceForVariableMetric(MeshContainer3dSurface* mesh, 
		Metric3dContext& mc, MeshPoint3d *point)
{
	//LOG4CPLUS_WARN(MeshLog::logger_console, "MG3dS::movePointByLaplaceFVM - to be rewritten");
	//return MeshGenerator3dSurface::movePointByLaplace( mesh, mc, point );

	assert( mesh != nullptr );
	if(point->isBorder()) return false;

	int rank = point->getRank();
	DataVector< MeshFace* > faces(rank);
	if( ! point->adjacentFaces( faces ) ) return false;

	SurfaceConstPtr surface = getLocalSurfaceForPoints( mc, mesh, 
		faces, point->localVicinityPoints() );
	if( surface == nullptr ) return false;

	const DPoint3d& sdpt_3d = point->getCoordinates();

	DataVector<double> mlens( rank );
	double len_ave = 0.0;

	for(int j = 0; j < rank; j++){
		MeshPoint3d* other_point = point->incidentPoint(j);
		mc.countMetricAtPoints(point, other_point);
		
		double len = mc.transformRStoMS( other_point->getCoordinates() - sdpt_3d ).length();
		mlens.add( len );
		len_ave += len;

		if( len < VERY_SMALL_NUMBER )
			return MeshGenerator3dSurface::movePointByLaplace( mesh, mc, point, surface );
	}

	len_ave /= rank;

	if( len_ave < 0.5 || len_ave > 2.0 ) {
		return MeshGenerator3dSurface::movePointByLaplace( mesh, mc, point, surface );
	}

	DPoint2d new_pt = DPoint2d::zero;
	//DPoint2d new_pt_normal = DPoint2d::zero;
	DPoint2d sparam = point->getLocalSurfaceParam( surface );
	double inv_rank = 1.0 / rank;

	//DataVector< DPoint3d > pref_points( rank );

	for(int j = 0; j < rank; j++) {
		MeshPoint3d* other_point = point->incidentPoint(j);
		DPoint2d other_param = other_point->getLocalSurfaceParam( surface );
		DPoint2d pref_param( other_param, sparam, len_ave / mlens[j] );
		new_pt.add( pref_param, inv_rank );

		//new_pt_normal.add( other_param, inv_rank );
		//pref_points.add( surface->getPoint( pref_param ) );
	}

	//if( surface->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, 0 ) ) {
	//	MeshViewSet* set = new MeshViewSet;
	//	set->addPoint( point, 1 );
	//	for(int i = 0; i < faces.countInt(); i++)
	//		set->addFaceWithEdgesAndPoints( faces[i] );
	//	for(int i = 0; i < rank; i++){
	//		set->addLabel( point->getEdge(i)->getPoint(0.5), to_string( mlens[i] ) );
	//		set->addPoint( pref_points[i], 2, i );
	//	}
	//	set->addLabel( surface->getPoint( new_pt ), "new_pt" );
	//	set->addLabel( surface->getPoint( new_pt_normal ), "norm_pt" );
	//	set->addInfo( "surface-tag-domain", surface->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 ) );
	//	set->addInfo( "ave_mlen", len_ave );
	//	SHOW_MESH("MG3d::movePointBLFVM", set);
	//}

	return point->tryMovingPoint( surface, new_pt, faces );
}

bool MeshGenerator3dSurface::moveBoundaryPointByLaplace( MeshContainer3dSurface* mesh, 
		Metric3dContext& mc, MeshPoint3d *point )
{
	if(!point->isBorder() || point->isBorder(TagBorder::FIXED | TagBorder::CORNER)) return false;

	DataVector<MeshPoint3d*> mpoints = point->localVicinityPointsBorderEdges( true );
	assert( mpoints.countInt() == 3 ); // ???

	const MeshEdge3d* medge01 = mpoints[0]->getEdgeToPoint( mpoints[1] );

	auto curve = medge01->getLocalCurve();
	assert( curve );

	double t = 0.0;


	for(int j = 1; j < mpoints.countInt(); j++){
		t += mpoints[j]->getLocalCurveParam( curve );
	}

	t /= (mpoints.countInt()-1);
 
	return point->tryMovingPoint( curve, t );

}

bool MeshGenerator3dSurface::moveBoundaryPointByLaplaceForVariableMetric(MeshContainer3dSurface* mesh, 
		Metric3dContext& mc, MeshPoint3d *point)
{
	if(!point->isBorder() || point->isBorder(TagBorder::FIXED | TagBorder::CORNER)) return false;
	if(!point->hasLocalCurve()) return false;

	DataVector<MeshPoint3d*> mpoints = point->localVicinityPointsBorderEdges( true );
	assert( mpoints.countInt() == 3 ); // ???

	auto medge01 = mpoints[0]->getEdgeToPoint( mpoints[1] );
	auto curve = medge01->getLocalCurve();
	assert( curve != nullptr );

	DPoint3d sdpt = point->getCoordinates();
	int brank = mpoints.countInt()-1;

	DataVector<double> mlens( brank );
	double len_ave = 0.0;

	for(int j = 1; j <= brank; j++){
		mc.countMetricAtPoints( point, mpoints[j] );
		double len = mc.transformRStoMS( mpoints[j]->getCoordinates() - sdpt ).length();

		mlens.add( len );
		len_ave += len;

		if( len < VERY_SMALL_NUMBER )
			return MeshGenerator3dSurface::moveBoundaryPointByLaplace( mesh, mc, point );
	}

	len_ave *= 0.5;

	if( len_ave < 0.5 || len_ave > 2.0 ) {
		return MeshGenerator3dSurface::moveBoundaryPointByLaplace( mesh, mc, point );
	}

	double new_t = 0.0;
	//DPoint2d new_pt_normal = DPoint2d::zero;
	double st = point->getLocalCurveParam( curve );

	//DataVector< DPoint3d > pref_points( rank );

	for(int j = 0; j < brank; j++) {
		double ot = mpoints[j+1]->getLocalCurveParam( curve );
		double k = len_ave / mlens[j];
		double pref_t =  ot * (1-k) + st * k;;
		new_t += pref_t * 0.5;
	}

	return point->tryMovingPoint( curve, new_t );
}

bool MeshGenerator3dSurface::smoothen(Metric3dContext& mc, MeshContainer3dSurface* mesh, int steps, 
		TagExtended::TagType tag_type, int tag_value, int method)
{
	LOG4CPLUS_INFO(MeshLog::logger_console,"Smoothing mesh (mixture of methods, inner) ...");

	bool ok = mesh->checkLocalSurfaceParams();

	for(int i = 0; i < steps; i++){

		if((method & MeshData::SM_TOP_SWAP) != 0){
			MeshGenerator3dSurface::smoothenTopologicalSwap(mc, mesh, tag_type, tag_value);
			bool ok = mesh->checkLocalSurfaceParams();
//			assert(mesh->isValid());
//			SHOW_MESH_NORESET("after smoothenTopologicalSwap", mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_DOMAIN));
		}
		if((method & MeshData::SM_LAPLACE) != 0){
			MeshGenerator3dSurface::smoothenLaplace(mc, mesh, false, tag_type, tag_value); 
			bool ok = mesh->checkLocalSurfaceParams();
//			assert(mesh->isValid());
//			SHOW_MESH_NORESET("after smoothenLaplace", mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_DOMAIN));
		}
		if((method & MeshData::SM_LAPLACE_MIXED) != 0){
			MeshGenerator3dSurface::smoothenLaplaceMixed(mc, mesh, tag_type, tag_value); 
			bool ok = mesh->checkLocalSurfaceParams();
//			assert(mesh->isValid());
//			SHOW_MESH_NORESET("after smoothenLaplaceMixed", mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));
		}
		if((method & MeshData::SM_LAPLACE_METRIC) != 0){
			MeshGenerator3dSurface::smoothenLaplace(mc, mesh, true, tag_type, tag_value); 
			bool ok = mesh->checkLocalSurfaceParams();
//			assert(mesh->isValid());
//			SHOW_MESH_NORESET("after smoothenLaplace+metric", mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_DOMAIN));
		}
		//if((method & MeshData::SM_METRIC) != 0){
		//	MeshGenerator3dSurface::smoothenMetric(mc, mesh, tag_type, tag_value); 
		//}
		if((method & MeshData::SM_DEL_SWAP) != 0){
			MeshGenerator3dSurface::smoothenQualitySwap(mc, mesh, tag_type, tag_value); 
			bool ok = mesh->checkLocalSurfaceParams();
//			assert(mesh->isValid());
			//MeshGenerator3dSurface::smoothenEdgeLengthSwap(mc, mesh, tag_type, tag_value); 
//			SHOW_MESH_NORESET("after smoothenQualitySwap", mesh->getViewSet(nullptr, TagExtended::TAG_ADAPT_SURF));
		}
		if((method & MeshData::SM_DEL_SWAP_COMPLETE) != 0){
			MeshGenerator3dSurface::smoothenQualitySwapComplete(mc, mesh, tag_type, tag_value);
			bool ok = mesh->checkLocalSurfaceParams();
//			assert(mesh->isValid());
//			SHOW_MESH_NORESET("after smoothenQualitySwapComplete", mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_DOMAIN));
		}
	}
	//MeshGenerator3dSurface::smoothenPostCheck(mc, mesh);

	return true;
}

bool MeshGenerator3dSurface::smoothenTopologicalSwap(Metric3dContext& mc, MeshContainer3dSurface *mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return false;

	for(int i = 0; i < count; i++){
		MeshPoint3d* p0 = mesh->getPointAt(i);
		if(p0->isBorder()) continue;
		if((tag_type != TagExtended::TAG_NONE) && !p0->hasAnyIntFlags(tag_type, tag_value)) 
			continue;

		for(int j = 0; j < p0->getRank(); j++){
			MeshEdge3d* edge = p0->getEdge(j);
			if(edge->getPointIndex(p0) == 0 && !edge->isBorder() && (edge->getFaceCount() == 2)){
				MeshTriangle3d* triangle1 = (MeshTriangle3d*)edge->getFaceAt(0);
				MeshTriangle3d* triangle2 = (MeshTriangle3d*)edge->getFaceAt(1);
				if(triangle1->getEdgeCount() != 3 || triangle2->getEdgeCount() != 3){
					continue;
				}
				MeshPoint3d* p1 = edge->getOtherPoint(p0);
				if(p1->isBorder()) continue;
				MeshPoint3d* p2 = triangle1->getOtherPoint(p0, p1);
				if(p2->isBorder()) continue;
				MeshPoint3d* p3 = triangle2->getOtherPoint(p0, p1);
				if(p3->isBorder()) continue;

				if(p2->getEdgeToPoint(p3)) continue; // such edge already exists...

				int diff = p0->getRank() + p1->getRank() - p2->getRank() - p3->getRank();
				if(diff > 2){

					DataVector<MeshFace*> mfaces(2);
					mfaces.add(triangle1); mfaces.add(triangle2);
					DataVector<MeshPoint3d*> mpoints(4);
					mpoints.add(p0); mpoints.add(p1);
					mpoints.add(p2); mpoints.add(p3);
					SurfaceConstPtr surface = getLocalSurfaceForPoints( mc, mesh, mfaces, mpoints );

					//if( surface != nullptr ) {
					//	bool valid1 = triangle1->valid( surface );
					//	bool valid2 = triangle2->valid( surface );
					//	if( !valid1  || !valid2 ) {
					//		MeshViewSet* set = new MeshViewSet;
					//		set->addFaceWithEdgesAndPoints( triangle1, 1 );
					//		set->addFaceWithEdgesAndPoints( triangle2, 2 );
					//		set->addInfo("valid-1", valid1 ? "yes" : "no" );
					//		set->addInfo("valid-2", valid2 ? "yes" : "no" );
					//		double nlen = (triangle1->getBoundingBox().getDiameter()*0.5);
					//		DPoint3d fmid = triangle1->getMiddlePoint();
					//		set->addEdge( fmid, fmid+triangle1->getNormalVector() * nlen, 1 );
					//		fmid = triangle2->getMiddlePoint();
					//		set->addEdge( fmid, fmid+triangle2->getNormalVector() * nlen, 1 );
					//		for(int j = 0; j < mpoints.countInt(); j++){
					//			const DPoint2d& param = mpoints[j]->getLocalSurfaceParam( surface );
					//			const DPoint3d& dpt = mpoints[j]->getCoordinates();
					//			//const DPoint2d param_test = surface->getParameters( dpt );
					//			set->addEdge( dpt, dpt + surface->getNormalVector( param ) * nlen, 0 );
					//		}
					//		SHOW_MESH("non-valid triangles?", set);
					//	}
					//}

					bool show_case = false;

					if(show_case)
						SHOW_MESH("Edge top-swap", mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));

					triangle1->switchPointsWithEdges(p0, p3);
					triangle2->switchPointsWithEdges(p1, p2);

					//mc.countMetricAtPoints(p0, p1);
					bool acceptable = (surface != nullptr) ? (triangle1->valid( surface ) && triangle2->valid( surface ) )
						: ( triangle1->validDirect(mc) && triangle2->validDirect(mc) );

					if(show_case)
						SHOW_MESH_NORESET(
							acceptable ? "Edge swap OK" : "Edge swap failed", 
							mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));

					if(!acceptable){
						// reverse operation
						triangle1->switchPointsWithEdges(p3, p0);
						triangle2->switchPointsWithEdges(p2, p1);
						if(show_case)
							SHOW_MESH_NORESET("Edge swap reversed", mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));
					}
				}
			}
		}
	}

	for(int i = 0; i < count; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if(point->isBorder()) continue;
		if((tag_type != TagExtended::TAG_NONE) && !point->hasAnyIntFlags(tag_type, tag_value)) 
			continue;

		int rank = point->getRank();
		if(rank < 5){
			if(MeshGenerator3dSurface::removeTriangulationPointSimple(mesh, point)){
				count = mesh->getPointsCount();
				i = 0;
			}
		}
	}
	
	return true;
}

bool MeshGenerator3dSurface::smoothenLaplaceMixed(Metric3dContext& mc, MeshContainer3dSurface *mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return true;

	for(int i = 0; i < count; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if((tag_type != TagExtended::TAG_NONE) && !point->hasAnyIntFlags(tag_type, tag_value)) 
			continue;

		double mg = mc.getMetricGradation(point->getCoordinates());
		if(mg > 1.5) movePointByLaplaceForVariableMetric(mesh, mc, point);
		else movePointByLaplace(mesh, mc, point);
	}

	return true;
}

bool MeshGenerator3dSurface::smoothenLaplace(Metric3dContext& mc, MeshContainer3dSurface* mesh, bool variable_metric,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return true;

	bool any_moved = false;
	for(int i = 0; i < count; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if((tag_type != TagExtended::TAG_NONE) && !point->hasAnyIntFlags(tag_type, tag_value)) 
			continue;

		//if( ! point->isBorder() && point->getEdge(0)->getFaceAt(0)->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, 9 ) ){
		//	// test
		//	testLaplace( mc, point );
		//}

		if(variable_metric)
			any_moved |= movePointByLaplaceForVariableMetric(mesh, mc, point);
		else
			any_moved |= movePointByLaplace(mesh, mc, point);
	}

	return any_moved;
}

bool MeshGenerator3dSurface::smoothenLaplaceBoundary(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		TagExtended::TagType tag_type, int tag_value, int forbid_tag_value)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return true;

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Smoothing mesh (Laplace for boundary) ...");

	bool any_moved = false;
	for(int i = 0; i < count; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if( !point->isBorder()) continue;
		if( point->isBorder(TagBorder::FIXED | TagBorder::CORNER) ) continue;

		if((tag_type != TagExtended::TAG_NONE) && !point->hasAnyIntFlags(tag_type, tag_value)) 
			continue;
		if((tag_type != TagExtended::TAG_NONE) && point->hasAnyIntFlags(tag_type, forbid_tag_value)) 
			continue;

		any_moved |= moveBoundaryPointByLaplace(mesh, mc, point);
	}

	return any_moved;
}

//////////////////////////////////////////////////////////////////////
// Removes point (and adjacent triangles) from the mesh, updates mesh
// by inserting new triangles within thus created cavity
// Elements adjacent to point have to be triangular.
// Delaunay criterion is not guaranteed.
bool MeshGenerator3dSurface::removeTriangulationPointSimple(
			MeshContainer3dSurface* mesh, MeshPoint3d *point,
			DataContainer<MeshEdge3d::ActiveEdge> * act_edges, TagExtended::TagType act_tag) 
{
//	assert(mesh->isValid());

	assert(point);
	if(point->isBorder()) return false;
	int rank = point->getRank();
	if(rank > 4 || rank < 3) return false;	// for now, only these simple cases are implemented

	if(rank == 3){
		// remove two triangles and modify the third one
		DataVector<MeshFace*> faces;
		if(!point->adjacentFaces(faces)) return false;

		for(int i = 0; i < rank; i++){ // triangles required ...
			if( faces[i]->getType() != FACE_TRIANGLE ) 
				return false;
		}

		if(act_edges){
			for(int i = 0; i < rank; i++){
				MeshEdge3d* edge = point->getEdge(i);
				MeshEdge3d::ActiveEdge* aact = (MeshEdge3d::ActiveEdge*)edge->getPtrTag(act_tag);
				if(aact){
					edge->removeTag(act_tag);
					delete act_edges->removeDataItem(aact->getIndex());
				}		
			}
		}

		MeshPoint3d* other_mp = nullptr;
		for(int i = 0; i < 3; i++){
			other_mp = point->getEdge(i)->getOtherPoint(point);
			if(!faces[0]->incidentToPoint(other_mp)) break;
		}
		delete mesh->removeMeshFace(faces.removeLast());
		delete mesh->removeMeshFace(faces.removeLast());
		faces[0]->switchPointsWithEdges(point, other_mp);
		delete mesh->removeMeshPoint(point);
	}else{
		return false;
		//DataVector<MeshFace*> faces;
		//if(!point->adjacentFaces(faces)) return false;
		//// choose the shortest edge
		//const DPoint3d& dpt = point->getCoordinates();
		//MeshEdge3d* min_edge = point->getEdge(0);
		//assert(!min_edge->isBorder());
		//double MIN_EDGE_LEN = dpt.distance2(min_edge->getOtherPoint(point)->getCoordinates());
		//for(int i = 1; i < 4; i++){
		//	MeshEdge3d* edge = point->getEdge(i);
		//	assert(!edge->isBorder());
		//	double len = dpt.distance2(edge->getOtherPoint(point)->getCoordinates());
		//	if(len < MIN_EDGE_LEN){
		//		min_edge = edge;
		//		MIN_EDGE_LEN = len;
		//	}
		//}
		//MeshPoint3d* other_mp = min_edge->getOtherPoint(point);
		//MeshFace* face1 = nullptr;
		//MeshFace* face2 = nullptr;
		//for(int i = 0; i < 4; i++){
		//	if(faces[i]->incidentToEdge(min_edge)) 
		//		delete mesh->removeMeshFace(faces[i]);
		//	else if(face1) face2 = faces[i];
		//	else face1 = faces[i];
		//}
		//face1->switchPointsWithEdges(point, other_mp);
		//face2->switchPointsWithEdges(point, other_mp);
		//delete mesh->removeMeshPoint(point);

		//// just in case..., check normals and swap edge if necessary
		//DVector3d vn1 = face1->getNormalVector();
		//DVector3d vn2 = face2->getNormalVector();
		//if(face1->getBlock(0) != face2->getBlock(0)) vn2.reverse();
		//if(vn1.scalarProduct(vn2) <= param_sharp_edge_threshold){
		//	//SHOW_MESH("ERROR After remove point x4", mesh->getDebugViewSet(other_mp, nullptr, 2.0, TagExtended::TAG_ADAPT_SURF));
		//	MeshEdge3d* edge = nullptr;
		//	for(int i = 0; i < face1->getEdgeCount(); i++){
		//		edge = face1->getEdge(i);
		//		if(edge->getOtherFace(face1) == face2) break;
		//	}
		//	MeshPoint3d* p0 = other_mp;
		//	MeshPoint3d* p1 = edge->getOtherPoint(p0);
		//	MeshPoint3d* p2 = face1->getOtherPoint(p0, p1);
		//	MeshPoint3d* p3 = face2->getOtherPoint(p0, p1);

		//	// ... make swap
		//	face1->switchPointsWithEdges(p0, p3);
		//	face2->switchPointsWithEdges(p1, p2);

		//	//SHOW_MESH("SWAP After remove point x4", mesh->getDebugViewSet(other_mp, nullptr, 2.0, TagExtended::TAG_ADAPT_SURF));
		//}
	}

//	assert(mesh->isValid());

	return true;
}

bool MeshGenerator3dSurface::smoothenEdgeLengthSwap(Metric3dContext& mc, MeshContainer3dSurface *mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return false;

	int total_count = 0;

	for(int i = 0; i < count; i++){
		MeshPoint3d* p0 = mesh->getPointAt(i);
		if((tag_type != TagExtended::TAG_NONE) && !p0->hasAnyIntFlags(tag_type, tag_value)) 
			continue;

		for(int j = 0; j < p0->getRank(); j++){
			MeshEdge3d* edge = p0->getEdge(j);
			if(edge->isBorder() || edge->getPointIndex(p0) != 0) continue;
			if(edge->getFaceCount() != 2) continue;
			MeshTriangle3d* triangle1 = (MeshTriangle3d*)edge->getFaceAt(0);
			MeshTriangle3d* triangle2 = (MeshTriangle3d*)edge->getFaceAt(1);
			if(triangle1->getEdgeCount() != 3 || triangle2->getEdgeCount() != 3)
				continue;

			MeshPoint3d* p1 = edge->getOtherPoint(p0);
			MeshPoint3d* p2 = triangle1->getOtherPoint(p0, p1);
			MeshPoint3d* p3 = triangle2->getOtherPoint(p0, p1);

			if(p2->getEdgeToPoint(p3)) continue; // such edge already exists...

			// check quality condition
			DPoint3d middle;
			middle.add(p0->getCoordinates(), 0.25);
			middle.add(p1->getCoordinates(), 0.25);
			middle.add(p2->getCoordinates(), 0.25);
			middle.add(p3->getCoordinates(), 0.25);
			mc.countMetricAtPoint(middle);
			double dlen01 = abs(1.0 - mc.transformRStoMS( p1->getCoordinates() - p0->getCoordinates()).length2());
			double dlen23 = abs(1.0 - mc.transformRStoMS( p3->getCoordinates() - p2->getCoordinates()).length2());

			if(dlen01 <= dlen23) continue; // current edge is closer to 1.0

			DataVector<MeshFace*> mfaces(2);
			mfaces.add(triangle1); mfaces.add(triangle2);
			DataVector<MeshPoint3d*> mpoints(4);
			mpoints.add(p0); mpoints.add(p1);
			mpoints.add(p2); mpoints.add(p3);

			SurfaceConstPtr surface = getLocalSurfaceForPoints(mc, mesh, mfaces, mpoints );
			bool show_case = false;

			// else swap
			if(show_case)
				SHOW_MESH("Edge len-swap", mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));

			// ... make swap
			triangle1->switchPointsWithEdges(p0, p3);
			triangle2->switchPointsWithEdges(p1, p2);

			// ... check
			bool acceptable = (surface != nullptr) ? (triangle1->valid( surface ) && triangle2->valid( surface ) )
				: ( triangle1->validDirect(mc) && triangle2->validDirect(mc) );

			if(show_case)
				SHOW_MESH_NORESET(
					acceptable ? "Edge swap OK" : "Edge swap failed", 
					mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));

			if(!acceptable){
				// reverse operation
				triangle1->switchPointsWithEdges(p3, p0);
				triangle2->switchPointsWithEdges(p2, p1);
				if(show_case)
					SHOW_MESH_NORESET("Edge swap reversed", mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));
			}
		}
	}
	return true;
}

bool MeshGenerator3dSurface::smoothenQualitySwapComplete(Metric3dContext& mc, MeshContainer3dSurface *mesh,
		TagExtended::TagType tag_type, int tag_value)
{

	if(!mesh) return false;
	int tct = mesh->getFacesCount();
	if(tct < 1) return false;
	int total_count = 0;

//	START_CLOCK("MG3dS::smoothenQualitySwapComplete");

	DataVector<MeshTriangle3d*> active_triangles(tct);
	for(int i = 0; i < tct; i++){
		MeshTriangle3d* triangle = (MeshTriangle3d*)mesh->getFaceAt(i);
		if(triangle->getEdgeCount() != 3) continue;
		if((tag_type != TagExtended::TAG_NONE) &&
			! triangle->hasAnyIntFlags(tag_type, tag_value)) continue;

		triangle->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
		triangle->setTagForEdges(TagExtended::TAG_MG2D_SM_SWAP);
		active_triangles.add(triangle);
	}

//	int counter = 0;
	int active_count[3] = {-1, -1, -1};
	DataVector<MeshFace*> mfaces( 2, nullptr );
	DataVector<MeshPoint3d*> mpoints(4);
	while(active_triangles.countInt() > 0){
		//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Del-swap-complete, step " << ++counter << ", active " << active_triangles.countInt());
		for(int i = 0; i < active_triangles.countInt(); ){
			MeshTriangle3d * triangle = active_triangles[i];
			mfaces[0] = triangle;
			triangle->setZeroTag(TagExtended::TAG_MG2D_SM_SWAP);
			for(int j = 0; j < 3; j++){
				MeshEdge3d* edge = triangle->getEdge(j);
				if(edge->zeroIntTag(TagExtended::TAG_MG2D_SM_SWAP) || edge->isBorder()) continue;
				MeshTriangle3d* other_triangle = (MeshTriangle3d*)edge->getOtherFace(triangle);
				mfaces[1] = other_triangle;
				MeshFace::getMeshPointsForFaces(mfaces, mpoints);
//				bool ok1 = mesh->checkLocalSurfaceParams();
				SurfaceConstPtr surface = getLocalSurfaceForPoints(mc, mesh, mfaces, mpoints );
				if(triangle->swapWithNeighbour(mc, surface, j, true, true, TagExtended::TAG_MG2D_SM_SWAP)){
//					bool ok3 = mesh->checkLocalSurfaceParams();
					++total_count;
					triangle->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
					if(other_triangle->zeroIntTag(TagExtended::TAG_MG2D_SM_SWAP)){
						if((tag_type == TagExtended::TAG_NONE) ||
							other_triangle->checkIntTag(tag_type, tag_value))
						{
							other_triangle->setIntTag(TagExtended::TAG_MG2D_SM_SWAP);
							active_triangles.add(other_triangle);
						}
					}
					break;
				}else{
//					bool ok4 = mesh->checkLocalSurfaceParams();
					triangle->getEdge(j)->setZeroTag(TagExtended::TAG_MG2D_SM_SWAP);
				}
			}
			if(triangle->zeroIntTag(TagExtended::TAG_MG2D_SM_SWAP)) 
				active_triangles.removeAt(i);
			else ++i;
		}
		active_count[0] = active_count[1];
		active_count[1] = active_count[2];
		active_count[2] = active_triangles.countInt();
		if(active_count[2] == active_count[1] || active_count[2] == active_count[0]) break;
	}

	mesh->removeAllTags(TagExtended::TAG_MG2D_SM_SWAP);
//	STOP_CLOCK("MG3dS::smoothenQualitySwapComplete");
	LOG4CPLUS_DEBUG(MeshLog::logger_console, "Quality-swapped edges count= " << total_count);

	return true;
}

bool MeshGenerator3dSurface::smoothenQualitySwap(Metric3dContext& mc, MeshContainer3dSurface *mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int count = mesh->getPointsCount();
	if(count < 1) return false;

	for(int i = 0; i < count; i++){
		MeshPoint3d* p0 = mesh->getPointAt(i);
		if((tag_type != TagExtended::TAG_NONE) && !p0->hasAnyIntFlags(tag_type, tag_value)) 
			continue;

		for(int j = 0; j < p0->getRank(); j++){
			MeshEdge3d* edge = p0->getEdge(j);
			if(edge->isBorder() || edge->getPointIndex(p0) != 0) continue;
			if(edge->getFaceCount() != 2) continue;
			MeshTriangle3d* triangle1 = (MeshTriangle3d*)edge->getFaceAt(0);
			MeshTriangle3d* triangle2 = (MeshTriangle3d*)edge->getFaceAt(1);
			if(triangle1->getEdgeCount() != 3 || triangle2->getEdgeCount() != 3)
				continue;

			MeshPoint3d* p1 = edge->getOtherPoint(p0);
			MeshPoint3d* p2 = triangle1->getOtherPoint(p0, p1);
			MeshPoint3d* p3 = triangle2->getOtherPoint(p0, p1);

			if(p2->getEdgeToPoint(p3)) continue; // such edge already exists...

			// check quality condition
			DPoint3d middle;
			middle.add(p0->getCoordinates(), 0.25);
			middle.add(p1->getCoordinates(), 0.25);
			middle.add(p2->getCoordinates(), 0.25);
			middle.add(p3->getCoordinates(), 0.25);
			mc.countMetricAtPoint(middle);

			DMPoint3d dmp0 = p0->getMetricCoordinates(mc);
			DMPoint3d dmp1 = p1->getMetricCoordinates(mc);
			DMPoint3d dmp2 = p2->getMetricCoordinates(mc);
			DMPoint3d dmp3 = p3->getMetricCoordinates(mc);

			double q1_pre = abs(DMTriangle3d::alphaQuality(dmp0, dmp1, dmp2));
			double q2_pre = abs(DMTriangle3d::alphaQuality(dmp0, dmp1, dmp3));

			double q1_after = abs(DMTriangle3d::alphaQuality(dmp1, dmp2, dmp3));
			double q2_after = abs(DMTriangle3d::alphaQuality(dmp0, dmp2, dmp3));

			if(std::min(q1_after, q2_after) <= std::min(q1_pre, q2_pre)) continue; // current min quality is better than after swap

			DataVector<MeshFace*> mfaces(2);
			mfaces.add(triangle1); mfaces.add(triangle2);
			DataVector<MeshPoint3d*> mpoints(4);
			mpoints.add(p0); mpoints.add(p1);
			mpoints.add(p2); mpoints.add(p3);

			SurfaceConstPtr surface = getLocalSurfaceForPoints(mc, mesh, mfaces, mpoints );
			bool show_case = false;

			// else swap
			if(show_case)
				SHOW_MESH("Edge len-swap", mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));

			// ... make swap
			triangle1->switchPointsWithEdges(p0, p3);
			triangle2->switchPointsWithEdges(p1, p2);

			// ... check
			bool acceptable = 
				triangle1->valid(surface) && 
				triangle2->valid(surface);

			if(show_case)
				SHOW_MESH_NORESET(
					acceptable ? "Edge swap OK" : "Edge swap failed", 
					mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));

			if(!acceptable){
				// reverse operation
				triangle1->switchPointsWithEdges(p3, p0);
				triangle2->switchPointsWithEdges(p2, p1);
				if(show_case)
					SHOW_MESH_NORESET("Edge swap reversed", mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));
			}
		}
	}
	return true;
}

bool MeshGenerator3dSurface::improveNearBorder(Metric3dContext& mc, MeshContainer3dSurface *mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	if(!mesh) return false;
	int fct = mesh->getFacesCount();
	if(fct < 1) return false;

	for(int i = 0; i < fct; i++){
		MeshFace* face = mesh->getFaceAt(i);
		if((tag_type != TagExtended::TAG_NONE) && !face->checkIntTag(tag_type, tag_value)) 
			continue;
		if(face->getEdgeCount() != 3) continue;

		DataVector<MeshEdge3d*> bedges(3);
		MeshEdge3d *edge;
		MeshEdge3d *iedge = nullptr;

		for(int j = 0; j < 3; j++)
			if((edge=face->getEdge(j))->isBorder()) 
				bedges.add(edge);
			else iedge = edge;
		if(bedges.countInt() != 2) continue;
		if(!iedge || iedge->getFaceCount() != 2) continue;

		MeshPoint3d* cpt = bedges[0]->commonVertex(bedges[1]);
		MeshPoint3d* p0 = bedges[0]->getOtherPoint(cpt);
		MeshPoint3d* p1 = bedges[1]->getOtherPoint(cpt);

		const DVector3d vt0 = (p0->getCoordinates() - cpt->getCoordinates()).normalized();
		const DVector3d vt1 = (p1->getCoordinates() - cpt->getCoordinates()).normalized();

		if(vt0.scalarProduct(vt1) > -0.5) continue;

		//SHOW_MESH("Bad boundary alignment", mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));

		MeshTriangle3d* triangle1 = (MeshTriangle3d*) iedge->getFaceAt(0);
		MeshTriangle3d* triangle2 = (MeshTriangle3d*) iedge->getFaceAt(1);

		MeshPoint3d* p2 = triangle1->getOtherPoint(p0, p1);
		MeshPoint3d* p3 = triangle2->getOtherPoint(p0, p1);

		DataVector<MeshFace*> mfaces(2);
		mfaces.add(triangle1); mfaces.add(triangle2);
		DataVector<MeshPoint3d*> mpoints(4);
		mpoints.add(p0); mpoints.add(p1);
		mpoints.add(p2); mpoints.add(p3);

		SurfaceConstPtr surface = getLocalSurfaceForPoints(mc, mesh, mfaces, mpoints );

		// ... make swap
		triangle1->switchPointsWithEdges(p0, p3);
		triangle2->switchPointsWithEdges(p1, p2);

		// ... check
		bool acceptable = 
			(surface != nullptr ) ? (triangle1->valid( surface ) && triangle2->valid( surface ))
			: ( triangle1->validDirect(mc) && triangle2->validDirect(mc) );

		//SHOW_MESH_NORESET(
		//	acceptable ? "Edge swap OK" : "Edge swap failed", 
		//	mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));

		if(!acceptable){
			// reverse operation
			triangle1->switchPointsWithEdges(p3, p0);
			triangle2->switchPointsWithEdges(p2, p1);
			//SHOW_MESH_NORESET("Edge swap reversed", mesh->getDebugViewSet(p0, p1, 2.0, TagExtended::TAG_ADAPT_SURF));
		}
	}
	return true;
}

/// Check and fix faces with normals inverted with respect to local surfaces
int MeshGenerator3dSurface::checkAndFixInvertedFacesForLocalSurfaces(
		Metric3dContext& mc, MeshContainer3dSurface* mesh,
		TagExtended::TagType tag_type, int tag_value)
{
	int fct = mesh->getFacesCount();
	int inverted_count = 0;

	for(int i = 0; i < fct; i++){
		MeshFace* face = mesh->getFaceAt(i);
		if((tag_type != TagExtended::TAG_NONE) && !face->checkIntTag(tag_type, tag_value)) 
			continue;
		
		if(! face->hasLocalSurface() ) continue;
		int orient = face->getLocalSurfaceOrientation();
		DVector3d svn = face->getLocalSurfaceNormal( );
		DVector3d fvn = face->getNormalVector();
		double sc = svn.scalarProduct(fvn);
		if(sc <= 0.0){
			inverted_count++;
			bool all_points_boundary = true;
			for(int j = 0; all_points_boundary && (j < face->getPointCount()); j++)
				all_points_boundary = face->getPoint(j)->isBorder();
			if(all_points_boundary){ // can't be helped, so it's better to live it alone...
				face->setLocalSurfaceOrientation(-orient);
				continue;
			}
//			LOG4CPLUS_INFO(MeshLog::logger_console, "Inverted face, n.n", sc);
//			SHOW_MESH("Inverted face", mesh->getDebugViewSet(face->getPoint(0), face->getPoint(1), 2.0,
//				TagExtended::TAG_ADAPT_SURF));
			if(sc > -1.0){
				SurfacePlane sp(face->getMiddlePoint(), svn, svn.crossProduct(fvn));
				for(int j = 0; j < face->getPointCount(); j++){
					MeshPoint3d* point = face->getPoint(j);
					if(point->isBorder()) continue;
					const DPoint3d& dpt = point->getCoordinates();
					DPoint3d spt = sp.getPoint( sp.getParameters( dpt ) );
					DVector3d dv = spt - dpt;
					mc.countMetricAtPoint(dpt);
					double fs = 2 * std::max(1.0, 0.2 / std::abs(sc));
					if( face->hasLocalSurface() ){
						SurfaceConstPtr face_surface = face->getOptLocalSurface();
						DPoint2d param = face_surface->getParametersNear( dpt + dv*fs, point->getLocalSurfaceParam( face_surface ) );
						DataVector< MeshFace* > faces;
						point->adjacentFaces( faces );
						point->tryMovingPoint( face_surface, param, faces );
					}else
						point->tryMovingPoint(mc, dpt + dv*fs);
				}
			fvn = face->getNormalVector();
			sc = svn.scalarProduct(fvn);
//			LOG4CPLUS_INFO(MeshLog::logger_console, "After fix, n.n", sc);
//			SHOW_MESH_NORESET("Inverted face after fix", mesh->getDebugViewSet(face->getPoint(0), face->getPoint(1), 2.0,
//				TagExtended::TAG_ADAPT_SURF));
			}
		}
	}

	return inverted_count;
}

/// Check and fix faces with normals inverted with respect to base normals
int MeshGenerator3dSurface::checkAndFixInvertedFacesForBaseNormal( MeshContainer3dSurface* mesh )
{
	int fct = mesh->getFacesCount();
	int inverted_count = 0;

	for(int i = 0; i < fct; i++){
		MeshFace* face = mesh->getFaceAt(i);
		if(!face->hasBaseNormal()) continue;

		DVector3d dn;
		if( ! face->checkAndGetNormalVector(dn) ) continue;
		
		double sp = dn.scalarProduct( face->getBaseNormal() );
		if(sp > 0.0) continue;
		
		inverted_count++;
		SHOW_MESH("Inverted face", mesh->getDebugViewSetTopological( face->getPoint(0), 2));

		int fpct = face->getPointCount();
		bool is_ok = false;
		for(int j = 0; (!is_ok) && (j < fpct); j++){
			MeshPoint3d* point = face->getPoint(j);
//			LOG4CPLUS_INFO(MeshLog::logger_console, "Pushing point", point->getIndex());
			DataVector<MeshFace*> faces;
			if(!point->adjacentFaces(faces)) continue;
			int pfct = faces.countInt();
//			DataVector<double> faces_sc(pfct);
			DVector3d fk_dn;
			//for(int k = 0; k < pfct; k++){
			//	if(faces[k]->hasBaseNormal() && faces[k]->checkAndGetNormalVector(fk_dn)){
			//		faces_sc.add( fk_dn.scalarProduct( faces[k]->getBaseNormal() ) );
			//	}else{
			//		faces_sc.add( -1.5 );
			//	}
			//}
			// move point
			const DPoint3d old_pt = point->getCoordinates();
			const DPoint3d prev_pt = face->getPoint((j-1+fpct)%fpct)->getCoordinates();
			const DPoint3d next_pt = face->getPoint((j+1)%fpct)->getCoordinates();
			const DPoint3d mid_pt(prev_pt, next_pt, 0.5);
			const DPoint3d new_pt(mid_pt, old_pt, -1.0);
			point->setCoordinates(new_pt);
			// check again sc for adjacent faces ...
			is_ok = true;
			for(int k = 0; is_ok && (k < pfct); k++){
				is_ok = (faces[k]->getShapeQuality() >= 0.0);
				if(!is_ok) break;
				if(!faces[k]->hasBaseNormal()) continue;
				double sc = -1.5;
				if(faces[k]->checkAndGetNormalVector(fk_dn))
					sc = fk_dn.scalarProduct( faces[k]->getBaseNormal() );
				is_ok = (sc > 0.0); // || (faces_sc[k] <= 0.0);
			}
			if(!is_ok){
//				SHOW_MESH("Inverted face - better now? NO", mesh->getDebugViewSetTopological( face->getPoint(0), 2));
				point->setCoordinates(old_pt);
			}else{
//				SHOW_MESH("Inverted face - better now? YES", mesh->getDebugViewSetTopological( face->getPoint(0), 2));
				//for(int k = 0; k < pfct; k++){
				//	MeshViewSet* set = new MeshViewSet;
				//	set->addFaceWithEdgesAndPoints(faces[k]);
				//	SHOW_MESH("Polygon after pushing", set);
				//}
				inverted_count--;
			}
		}
	}

	return inverted_count;
}

/// Creates surface mesh from volume mesh (copy boundary faces and entities)
MeshContainer3dSurface* MeshGenerator3dSurface::copySurfaceMeshFromVolumeMesh(MeshContainer3d* mesh3d, MeshDomainVolume* mdv)
{
	if(!mesh3d) return nullptr;

	int pct = mesh3d->getPointsCount();
	int pct_border = 0;
	for(int i = 0; i < pct; i++)
		if(mesh3d->getPointAt(i)->isBorder()) ++pct_border;
	int fct_border = 0;
	for(IteratorFace it = mesh3d->getFirstFace(); it.isValid(); it.nextFace())
		if(it.getFace()->isBorder()) ++fct_border;

	if(pct_border < 3 || fct_border < 1) return nullptr;

	MeshContainer3dSurface* surface_mesh = new MeshContainer3dSurface(pct_border);

	DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> hash_points(2*pct_border, nullptr);

	for(int i = 0; i < pct; i++){
		MeshPoint3d* mp = mesh3d->getPointAt(i);
		if(!mp->isBorder()) continue;
		MeshPoint3d* mp_copy = new MeshPoint3d(*mp);
		hash_points.insert(mp, mp_copy);
		surface_mesh->addMeshPoint(mp_copy);
	}

	if(!mdv){
		// create additional "dummy" volume
		auto mdvs = std::make_shared<MeshDomainVolume>();
		surface_mesh->addDomainVolume(mdvs);
		mdv = mdvs.get();
		mdv->setAreaID(1);
	}

	// discrete faces
	for(IteratorFace it = mesh3d->getFirstFace(); it.isValid(); it.nextFace()){
		MeshFace* face = it.getFace();
		if(!face->isBorder() || face->getType() != FACE_TRIANGLE) continue;
		assert(face->getType() == FACE_TRIANGLE);
		MeshFace* face_copy;
		if(face->getBlock(0)){
			face_copy = new MeshTriangle3d(
				hash_points.getValue(face->getPoint(2), nullptr),
				hash_points.getValue(face->getPoint(1), nullptr),
				hash_points.getValue(face->getPoint(0), nullptr),
				mdv);
		}else{
			face_copy = new MeshTriangle3d(
				hash_points.getValue(face->getPoint(0), nullptr),
				hash_points.getValue(face->getPoint(1), nullptr),
				hash_points.getValue(face->getPoint(2), nullptr),
				mdv);
		}
		surface_mesh->addMeshFace(face_copy);
	}

	surface_mesh->setControlSpace(mesh3d->getControlSpace());

	return surface_mesh;
}

/// Create new adaptive control space with sizing taken from mesh faces
CS3dAPtr MeshGenerator3dSurface::createACSfromMeshFaces(MeshContainer3dSurface* mesh)
{
	if(!mesh) return nullptr;
	int pct = mesh->getPointsCount();
	if(pct == 0) return nullptr;

	DBox box;
	for(int i = 0; i < pct; i++)
		box.addPoint(mesh->getPointAt(i)->getCoordinates());

	box.inflate(ControlSpace2d::param_inflate_box_factor);

	auto space = std::make_shared<ControlSpace3dOctree>(box);
	space->setMaxMetric();
	if(space->addSimpleBoundaryControlData(mesh))
		space->smoothen();

	return space;
}

/// Update sizing data in the acs with curvature information aproximated from local parameterization surfaces (if available)
bool MeshGenerator3dSurface::updateACSwithLocalCurvature(CS3dPtr cs, MeshContainer3dSurface* mesh)
{
	if(!cs || !cs->isAdaptive()) return false;
	if(!mesh) return false;

	int pct = mesh->getPointsCount();
	if(pct == 0) return false;

	auto acs = cs->getAsAdaptive();

	double model_diameter = acs->getBoundingBoxDiameter();
	bool any_change = false;
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		// ... local surface
		auto surface = point->getLocalValidSurface();
		if(surface) 
			any_change |= acs->setMinControlFromSurfaceCurvature(point->getCoordinates(), surface, model_diameter);

		// ... local curves
		auto curve = point->getLocalCurve();
		if(curve)
			any_change |= acs->setMinControlFromCurveCurvature(point->getCoordinates(), curve, model_diameter);
	}

	return any_change;
}


/// Create surface mesh from point/face data, with small fixes
MeshContainer3dSurface* MeshGenerator3dSurface::createMeshFromFaces(
	DataVector<DPoint3d> & points, DataSimpleList< DataVector<int> > & faces)
{
	int pct = points.countInt();
	auto surface_mesh = new MeshContainer3dSurface(pct);

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = new MeshPoint3d(points[i]);
		surface_mesh->addMeshPoint(point);
	}

	for(auto it = faces.iterator(); it.valid(); it.moveNext()){
		const DataVector<int> & face = it.item();
		int fpct = face.countInt();
		assert(fpct > 2);
		MeshFace* mface = nullptr;
		if(fpct == 3){
			mface = new MeshTriangle3d(
				surface_mesh->getPointAt(face[0]), 
				surface_mesh->getPointAt(face[1]), 
				surface_mesh->getPointAt(face[2]));
			surface_mesh->addMeshFace(mface);
		}else if(fpct == 4){
			mface = new MeshQuad3d(
				surface_mesh->getPointAt(face[0]), 
				surface_mesh->getPointAt(face[1]), 
				surface_mesh->getPointAt(face[2]), 
				surface_mesh->getPointAt(face[3]));
			surface_mesh->addMeshFace(mface);
		}else{
			DataVector<MeshPoint3d*> spoints(fpct);
			for(int i = 0; i < fpct; i++)
				spoints.add( surface_mesh->getPointAt(face[i]) );
			mface = new MeshPoly3d(spoints);
			surface_mesh->addMeshFace(mface);
		}
	}

	surface_mesh->setControlSpace( std::make_shared<ControlSpace3dIdentity>() );
	auto mdv = std::make_shared<MeshDomainVolume>(surface_mesh);
	mdv->setAreaID(0);
	surface_mesh->addDomainVolume(mdv);

	return surface_mesh;
}

/// Check and fix all surface meshes for the domain mesh
int MeshGenerator3dSurface::validateSurfaceMesh(MeshContainer3dSurface* surface_mesh, double min_dist)
{
	assert(surface_mesh);
	if(!surface_mesh) return 0;
	int result = 0;

//	START_CLOCK("MG3dS:validateSurfaceMesh");

	// remove too short edges
	//SHOW_MESH("validateSurfaceMesh-start", surface_mesh->getViewSet() );
	if(min_dist > 0.0){
		CS3dPtr space(new ControlSpace3dIdentity(min_dist / MIN_EDGE_LEN));
		Metric3dContext mc(space);
		mc.setMinSkipLen2( 0.25 * min_dist*min_dist );
		int loop_result = 0;
		int loop_max = 3;
		do{
			loop_result = 0;
			int local_result = MeshGenerator3dSurface::collapseInnerEdges(mc, surface_mesh);
			local_result += MeshGenerator3dSurface::collapseBoundaryEdges(mc, surface_mesh);
			if(local_result > 0) {
				loop_result += local_result;
				//SHOW_MESH("validateSurfaceMesh-before fix-thin", surface_mesh->getViewSet() );
			}
			local_result = MeshGenerator3dSurface::fixThinTriangles(mc, surface_mesh);
			if(local_result > 0) {
				loop_result += local_result;
				//SHOW_MESH("validateSurfaceMesh-before inverted", surface_mesh->getViewSet() );
			}
			local_result = MeshGenerator3dSurface::fixInvertedTriangles(mc, surface_mesh);
			if(local_result > 0) {
				loop_result += local_result;
				//SHOW_MESH("validateSurfaceMesh-after inverted", surface_mesh->getViewSet() );
			}

			result += loop_result;
		}while( (loop_result > 0) && (--loop_max >= 0) );
	}
	//SHOW_MESH("validateSurfaceMesh-before check-fix-inverted", surface_mesh->getViewSet() );
	MeshGenerator3dSurface::checkAndFixInvertedFacesForBaseNormal(surface_mesh);
	//SHOW_MESH("validateSurfaceMesh-after check-fix-inverted", surface_mesh->getViewSet() );

//	STOP_CLOCK("MG3dS:validateSurfaceMesh");

	//double min_len2 = LARGE_NUMBER;
	//for(auto it = surface_mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
	//	double len2 = it.getEdge()->getLength2();
	//	if(len2 < min_len2) min_len2 = len2;
	//}
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Min. edge length after validation: " << sqrt(min_len2) << "  (" << ( sqrt(min_len2) / min_dist ) << " of min_dist)");

	return result;
}

/// Collapse inverted faces, return number of fixed faces
int MeshGenerator3dSurface::fixInvertedTriangles(Metric3dContext& mc, MeshContainer3dSurface* mesh)
{
	int fct = mesh->getFacesCount();
	vector< MeshFace* > inv_faces;
	for(int i = 0; i < fct; i++){
		MeshTriangle3d* tri = (MeshTriangle3d*)mesh->getFaceAt(i);
		if(tri->getType() != FACE_TRIANGLE) continue; // inverted triangles only, for polys - join?

		MeshEdge3d* edges[3] = { tri->getEdge(0), tri->getEdge(1), tri->getEdge(2) };
		if( edges[0]->isBorder() && edges[1]->isBorder() && edges[2]->isBorder() ) continue; // what can be done

		int inverted_rank = tri->getInvertedRank();
		if( inverted_rank >= MeshFace::SP_BAD ) { // at least one such edge
			tri->setIntTag( TagExtended::TAG_INV_RANK, inverted_rank );
			inv_faces.push_back( tri );
		}
	}

	if( inv_faces.empty() ) return 0;


	//if(true) {
	//	DVector3d vn;
	//	double min_sp = 1.0;
	//	for(int i = 0; i < fct; i++ ) {
	//		MeshFace* face = mesh->getFaceAt(i);
	//		int iq = 0;
	//		if( face->checkAndGetNormalVector( vn ) ){
	//			DVector3d vbn = face->getBaseNormal();
	//			if( !vbn.isZero() ) {
	//				double sp = vbn.scalarProduct( vn );
	//				if(sp < min_sp) min_sp = sp;
	//				iq = (int) ((1.0+sp)*50);
	//			}
	//		}
	//		face->setIntTag( TagExtended::TAG_VISUALIZATION, iq );
	//	}
	//	MeshViewSet* set = mesh->getViewSet(nullptr, TagExtended::TAG_VISUALIZATION);
	//	set->addInfo("min_sp", min_sp);
	//	SHOW_MESH("fixInvertedTriangles::sp base_normal x normal", set);
	//}

	sort( inv_faces.begin(), inv_faces.end(), 
		[] ( MeshFace *f1, MeshFace *f2 ) { 
			return f1->getIntTag( TagExtended::TAG_INV_RANK, 0 ) > f2->getIntTag( TagExtended::TAG_INV_RANK, 0 );
		} 
	);

	int fixed_count = 0;
	bool show_case = false;

	for( auto face : inv_faces ) {
		int inv_index;
		int inverted_rank = face->getInvertedRank( &inv_index );
		if( inverted_rank < MeshFace::SP_BAD ) continue; // was corrected in the meantime, apparently ...

		//show_case = (fixed_count >= 46);

		MeshEdge3d* inv_edge = face->getEdge(inv_index);
		if(inv_edge->getFaceCount() != 2){
			LOG4CPLUS_WARN(MeshLog::logger_console, "MG3dS::fixThinTriangles, multi-face edge, todo...");
			continue;
		}
		MeshFace* other_face = inv_edge->getOtherFace(face);
		if( other_face->getType() != FACE_TRIANGLE ) continue;

		DVector3d tri_vn = face->getBaseNormal();
		DVector3d otri_vn = other_face->getBaseNormal();

		//if( true ) {
		//	double min_sp = 1.0;
		//	DVector3d vn;
		//	if( face->checkAndGetNormalVector( vn ) ){
		//		double sp = vn.scalarProduct( tri_vn );
		//		if(sp < min_sp) min_sp = sp;
		//	}
		//	if( other_face->checkAndGetNormalVector( vn ) ){
		//		double sp = vn.scalarProduct( otri_vn );
		//		if(sp < min_sp) min_sp = sp;
		//	}
		//	if( min_sp < 0.5 ) {
		//		LOG4CPLUS_WARN(MeshLog::logger_console, "min_sp", min_sp);
		//	}
		//}

		if( show_case ){
			if(true) { 
				MeshViewSet* set = new MeshViewSet;
				set->addFaceWithEdgesAndPoints(face, 1);
				DPoint3d fmid = face->getMiddlePoint();
				set->addEdge( fmid, fmid+face->getNormalVector() * (face->getBoundingBoxDiameter()*0.5), 1 );
				for( int j = 0; j < 3; j++){
					MeshFace* face_j = face->getEdge(j)->getOtherFace( face );
					set->addFaceWithEdgesAndPoints(face_j, (face_j==other_face) ? 2 : 0);
					fmid = face_j->getMiddlePoint();
					set->addEdge( fmid, fmid+face_j->getNormalVector() * (face_j->getBoundingBoxDiameter()*0.5), 0 );
					set->addLabel( fmid, to_string(j) );
				}
				set->addInfo("inverted_rank", inverted_rank);
				set->addInfo("inverted_rank - other", other_face->getInvertedRank() );
				SHOW_MESH("inverted? - getInvertedRank", set);
			}
			if(false) {
				MeshViewSet* set = mesh->getDebugViewSetTopological( face, 2);
				set->addInfo("inverted_rank", inverted_rank);
				set->addInfo("inverted_rank - other", other_face->getInvertedRank() );
				set->addInfo("inv_edge", to_string(inv_edge->getMeshPoint(0)->getIndex()) + "-" +
					to_string(inv_edge->getMeshPoint(1)->getIndex()) );
				set->addInfo( "points", to_string( face->getPoint(0)->getIndex() ) + "-" +
									to_string( face->getPoint(1)->getIndex() ) + "-" +
									to_string( face->getPoint(2)->getIndex() ) );
				SHOW_MESH("inverted?", set);
			}
		}

		// method first - as in fixThinTriangles ... 

		mc.countMetricAtPoint( face->getMiddlePoint() );

		MeshPoint3d* point_max = face->getPoint((inv_index+2)%3);
		
		const DPoint3d& ept0 = inv_edge->getMeshPoint(0)->getCoordinates();
		const DPoint3d& ept1 = inv_edge->getMeshPoint(1)->getCoordinates();
		DVector3d evt = ept1-ept0;
		double t = evt.scalarProduct(point_max->getCoordinates() - ept0) / evt.length2();

		if(t > 0 && t < 1.0) {
			DPoint3d opp_pt(ept0, ept1, t);

			double area_1 = DTriangle3d::area( face->getPoint(0)->getCoordinates(),
					face->getPoint(1)->getCoordinates(), face->getPoint(2)->getCoordinates() );
			double area_2 = DTriangle3d::area( other_face->getPoint(0)->getCoordinates(),
					other_face->getPoint(1)->getCoordinates(), other_face->getPoint(2)->getCoordinates() );
								
			if( area_1 < area_2) face->setBaseNormal( otri_vn );
			else other_face->setBaseNormal( tri_vn );

			SurfaceConstPtr surface = nullptr;
			if( face->hasLocalSurface() ) {
				DataVector<MeshFace*> mfaces(2);
				mfaces.add(face); mfaces.add(other_face);
				DataVector<MeshPoint3d*> mpoints(4);
				MeshFace::getMeshPointsForFaces( mfaces, mpoints );

				surface = getLocalSurfaceForPoints(mc, mesh, mfaces, mpoints );
			}

			if( face->swapWithNeighbour(mc, surface, inv_index, false) ) {
				//point_max->setCoordinates( opp_pt );
				++fixed_count;
				if( surface == nullptr ){
					face->setBaseNormal( face->getNormalVector() );
					other_face->setBaseNormal( other_face->getNormalVector() );
				}

				if( show_case ){
					MeshViewSet* set = mesh->getDebugViewSetTopological( face, 2);
					set->addInfo("inverted_rank", face->getInvertedRank() );
					set->addInfo("inverted_rank - other", other_face->getInvertedRank() );
					set->addInfo( "points", to_string( face->getPoint(0)->getIndex() ) + "-" +
										to_string( face->getPoint(1)->getIndex() ) + "-" +
										to_string( face->getPoint(2)->getIndex() ) );
					set->addInfo("fixed count", fixed_count);
					SHOW_MESH_NORESET("inverted triangle after fix-1", set);
				}

				continue;
			}else{
				face->setBaseNormal(tri_vn); // restore ...
				other_face->setBaseNormal(otri_vn);

				if(show_case) { 
					MeshViewSet* set = new MeshViewSet;
					set->addFaceWithEdgesAndPoints(face, 1);
					DPoint3d fmid = face->getMiddlePoint();
					set->addEdge( fmid, fmid+face->getNormalVector() * (face->getBoundingBox().getDiameter()*0.5), 1 );
					for( int j = 0; j < 3; j++){
						MeshFace* face_j = face->getEdge(j)->getOtherFace( face );
						set->addFaceWithEdgesAndPoints(face_j, (face_j==other_face) ? 2 : 0);
						fmid = face_j->getMiddlePoint();
						set->addEdge( fmid, fmid+face_j->getNormalVector() * (face_j->getBoundingBox().getDiameter()*0.5), 0 );
						set->addLabel( fmid, to_string(j) );
					}
					set->addInfo("inverted_rank", inverted_rank);
					set->addInfo("inverted_rank - other", other_face->getInvertedRank() );
					set->addInfo("face index", face->getIndex() );
					SHOW_MESH("inverted - fix-1 FAILED", set);
				}
			}
		}else{
			if(show_case) { 
				MeshViewSet* set = new MeshViewSet;
				set->addFaceWithEdgesAndPoints(face, 1);
				DPoint3d fmid = face->getMiddlePoint();
				set->addEdge( fmid, fmid+face->getNormalVector() * (face->getBoundingBox().getDiameter()*0.5), 1 );
				for( int j = 0; j < 3; j++){
					MeshFace* face_j = face->getEdge(j)->getOtherFace( face );
					set->addFaceWithEdgesAndPoints(face_j, (face_j==other_face) ? 2 : 0);
					fmid = face_j->getMiddlePoint();
					set->addEdge( fmid, fmid+face_j->getNormalVector() * (face_j->getBoundingBox().getDiameter()*0.5), 0 );
					set->addLabel( fmid, to_string(j) );
				}
				set->addInfo("inverted_rank", inverted_rank);
				set->addInfo("inverted_rank - other", other_face->getInvertedRank() );
				set->addInfo("face index", face->getIndex() );
				set->addInfo("t", t );
				SHOW_MESH("inverted - fix-1 FAILED (t out of range)", set);
			}
		}

		// method second - else, try harder with both faces...
		DataVector<MeshFace*> swap_tri(4);
		DataVector<int> swap_k(4);
		DataVector<double> swap_sp(4);
		for(int k = 0; k < 3; k++){
			MeshEdge3d* nedge = face->getEdge(k);
			if(nedge->getFaceCount() != 2) continue;
			MeshFace* nface = nedge->getOtherFace(face);
			if(nface == other_face) continue;
			DVector3d ndn;
			if( ! nface->checkAndGetNormalVector(ndn) ) continue;
			//double nsp = 
			tri_vn.scalarProduct(ndn);
			double onsp = otri_vn.scalarProduct(ndn);
//			LOG4CPLUS_INFO(MeshLog::logger_mesh, "f0-nf" << k << " - sp= " << nsp << ", other_sp= " << onsp);
			swap_tri.add(face);
			swap_k.add(k);
			swap_sp.add(onsp);
		}
		for(int k = 0; k < 3; k++){
			MeshEdge3d* nedge = other_face->getEdge(k);
			if(nedge->getFaceCount() != 2) continue;
			MeshFace* nface = nedge->getOtherFace(other_face);
			if(nface == face) continue;
			DVector3d ndn;
			if( ! nface->checkAndGetNormalVector(ndn) ) continue;
			double nsp = tri_vn.scalarProduct(ndn);
			//double onsp = 
			otri_vn.scalarProduct(ndn);
//						LOG4CPLUS_INFO(MeshLog::logger_mesh, "f1-nf" << k << " - sp= " << onsp << ", other_sp= " << nsp);
			swap_tri.add(other_face);
			swap_k.add(k);
			swap_sp.add(nsp);
		}

		bool swap_ok = false;
		while(!swap_ok && swap_tri.notEmpty()){
			int best_i = 0;
			for(int l = 1; l < swap_tri.countInt(); l++){
				if( swap_sp[l] > swap_sp[best_i] )
					best_i = l;
			}
			MeshFace* swap_t = swap_tri[best_i];
			MeshFace* swap_other = swap_t->getEdge(swap_k[best_i])->getOtherFace(swap_t);

			swap_t->setBaseNormal( swap_other->getBaseNormal() );
			swap_ok = ( swap_t->swapWithNeighbour(mc, nullptr, swap_k[best_i], false) != nullptr );
			swap_tri.removeAt(best_i);
			swap_k.removeAt(best_i);
			swap_sp.removeAt(best_i);
		}

		if(swap_ok) {
			face->setBaseNormal( face->getNormalVector() );
			other_face->setBaseNormal( other_face->getNormalVector() );
//			SHOW_MESH_NORESET("inverted? after extra swap", mesh->getDebugViewSetTopological( face, 2));
			++fixed_count;

			if( show_case ){
				MeshViewSet* set = mesh->getDebugViewSetTopological( face, 2);
				set->addInfo("inverted_rank", face->getInvertedRank() );
				set->addInfo("inverted_rank - other", other_face->getInvertedRank() );
				set->addInfo( "points", to_string( face->getPoint(0)->getIndex() ) + "-" +
									to_string( face->getPoint(1)->getIndex() ) + "-" +
									to_string( face->getPoint(2)->getIndex() ) );
				set->addInfo("fixed count", fixed_count);
				SHOW_MESH_NORESET("inverted triangle after fix-2", set);
			}

			continue;
		}

		// failed, restore...
		face->setBaseNormal( otri_vn );
		other_face->setBaseNormal( otri_vn );

		if( show_case ){
			MeshViewSet* set = mesh->getDebugViewSetTopological( face, 2);
			set->addInfo("inverted_rank", face->getInvertedRank() );
			set->addInfo("inverted_rank - other", other_face->getInvertedRank() );
			set->addInfo( "points", to_string( face->getPoint(0)->getIndex() ) + "-" +
								to_string( face->getPoint(1)->getIndex() ) + "-" +
								to_string( face->getPoint(2)->getIndex() ) );
			set->addInfo("fixed count", fixed_count);
			SHOW_MESH_NORESET("inverted triangle fix - failed", set);
		}

		// SHOW_MESH_NORESET("inverted? swap failed", mesh->getDebugViewSetTopological( face, 2));	
	}
			
	return fixed_count;
}

/// Collapse too thin (or inverted) faces, return number of fixed faces
int MeshGenerator3dSurface::fixThinTriangles(Metric3dContext& mc, MeshContainer3dSurface* mesh)
{
	static const double MIN_H = 0.1;
	static const double SF2 = 4.0 / (MIN_H * MIN_H);
	int fixed_count = 0;

	for(int i = 0; i < mesh->getFacesCount(); i++){
		MeshTriangle3d* tri = (MeshTriangle3d*)mesh->getFaceAt(i);
		if(tri->getType() != FACE_TRIANGLE) continue;
		MeshEdge3d* edges[3] = { tri->getEdge(0), tri->getEdge(1), tri->getEdge(2) };
		if( edges[0]->isBorder() && edges[1]->isBorder() && edges[2]->isBorder() ) continue;

		mc.countMetricAtPoint( tri->getMiddlePoint() );

		DMPoint3d mpts[3] = { 
			tri->getPoint(0)->getMetricCoordinates(mc),
			tri->getPoint(1)->getMetricCoordinates(mc),
			tri->getPoint(2)->getMetricCoordinates(mc) };

		double s = DMTriangle3d::area( mpts[0], mpts[1], mpts[2] );
		double a2_threshold = SF2 * s * s;
		int i_max = 0;
		double a2_max = mpts[0].distance2(mpts[1]);
		double a2 = mpts[1].distance2(mpts[2]);
		if(a2 > a2_max) { a2_max = a2; i_max = 1; }
		a2 = mpts[2].distance2(mpts[0]);
		if(a2 > a2_max) { a2_max = a2; i_max = 2; }
		if( a2_max < a2_threshold) continue;

		MeshEdge3d* opposite_edge = edges[i_max];
		if(opposite_edge->getFaceCount() != 2){
			LOG4CPLUS_WARN(MeshLog::logger_console, "MG3dS::fixThinTriangles, multi-face edge, todo...");
			continue;
		}
		MeshPoint3d* point_max = tri->getPoint((i_max+2)%3);
		
		const DPoint3d& pt0 = opposite_edge->getMeshPoint(0)->getCoordinates();
		const DPoint3d& pt1 = opposite_edge->getMeshPoint(1)->getCoordinates();
		DVector3d evt = pt1-pt0;
		double t = evt.scalarProduct(point_max->getCoordinates() - pt0) / evt.length2();

		bool show_case = false; 
		if( show_case ){
			MeshViewSet* set = new MeshViewSet;
			DataVector<MeshFace*> faces;
			tri->getPoint(0)->adjacentFaces(faces);
			tri->getPoint(1)->adjacentFaces(faces);
			tri->getPoint(2)->adjacentFaces(faces);
			for(int k = 0; k < faces.countInt(); k++)
				set->addFaceWithEdgesAndPoints(faces[k]);
			set->addInfo("pt1", tri->getPoint(0)->getIndex());
			set->addInfo("pt2", tri->getPoint(1)->getIndex());
			set->addInfo("pt3", tri->getPoint(2)->getIndex());
			set->addInfo("triangle mc-quality", DMTriangle3d::alphaQuality(mpts[0], mpts[1], mpts[2]));
			set->addInfo("h len", 2.0 * s / sqrt(a2_max));
			set->addInfo("pt max", point_max->getIndex());
			set->addInfo("t opposite", t);
			SHOW_MESH("thin triangle to fix?", set);
		}
		if(t <= 0 || t >= 1.0) {
			//LOG4CPLUS_WARN(MeshLog::logger_console, "MG3dS::fixThinTriangles, t outside [0,1], why?");
			continue;
		}

		DPoint3d opp_pt(pt0, pt1, t);
		MeshFace* other_face = opposite_edge->getOtherFace(tri);
		if(other_face->getType() == FACE_TRIANGLE){
			DVector3d tri_vn = tri->getBaseNormal();
			DVector3d otri_vn = other_face->getBaseNormal();
			tri->setBaseNormal( otri_vn );

			SurfaceConstPtr surface = nullptr;
			if( tri->hasLocalSurface() ) {
				DataVector<MeshFace*> mfaces(2);
				mfaces.add(tri); mfaces.add(other_face);
				DataVector<MeshPoint3d*> mpoints(4);
				MeshFace::getMeshPointsForFaces( mfaces, mpoints );

				surface = getLocalSurfaceForPoints(mc, mesh, mfaces, mpoints );
			}

			if( tri->swapWithNeighbour(mc, surface, i_max, false) ) {
				point_max->setCoordinates( opp_pt );
				++fixed_count;
				if( surface == nullptr ){
					tri->setBaseNormal( tri->getNormalVector() );
					other_face->setBaseNormal( other_face->getNormalVector() );
				}
				if( show_case ){
					MeshViewSet* set = new MeshViewSet;
					DataVector<MeshFace*> faces;
					tri->getPoint(0)->adjacentFaces(faces);
					tri->getPoint(1)->adjacentFaces(faces);
					tri->getPoint(2)->adjacentFaces(faces);
					for(int k = 0; k < faces.countInt(); k++)
						set->addFaceWithEdgesAndPoints(faces[k]);
					set->addInfo("pt_max", point_max->getIndex());
					set->addInfo("fixed count", fixed_count);
					SHOW_MESH_NORESET("thin triangle after fix", set);
				}
			}else{
				tri->setBaseNormal(tri_vn); // restore ...
				if( show_case ){
					MeshViewSet* set = new MeshViewSet;
					DataVector<MeshFace*> faces;
					tri->getPoint(0)->adjacentFaces(faces);
					tri->getPoint(1)->adjacentFaces(faces);
					tri->getPoint(2)->adjacentFaces(faces);
					for(int k = 0; k < faces.countInt(); k++)
						set->addFaceWithEdgesAndPoints(faces[k]);
					set->addInfo("pt1", tri->getPoint(0)->getIndex());
					set->addInfo("pt2", tri->getPoint(1)->getIndex());
					set->addInfo("pt3", tri->getPoint(2)->getIndex());
					set->addInfo("pt max", point_max->getIndex());
					SHOW_MESH("FAILED?", set);
				}
			}
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "MG3dS::fixThinTriangles, combine with quad/poly?");			
		}
	}

	return fixed_count;
}


CS3dPtr MeshGenerator3dSurface::createACSFromApproxCurvatureDirect(MeshContainer3dSurface* mesh)
{
	if(!mesh) return nullptr;
	int pct = mesh->getPointsCount();
	if(pct == 0) return nullptr;

//	DataVector<double> edge_ave_lens(pct, 0.0);
	START_CLOCK("MG3dS:createACSFromApproxCurvatureDirect");

	DBox box;
	for(int i = 0; i < pct; i++){
		MeshPoint3d* mp = mesh->getPointAt(i);
		box.addPoint(mp->getCoordinates());
		//int rank = mp->getRank();
		//assert(rank > 0);
		//double sum = 0.0;
		//for(int j = 0; j < rank; j++){
		//	sum += mp->getEdge(j)->getLength();
		//}
		//edge_ave_lens[i] = sum/rank;
	}

	box.inflate(ControlSpace2d::param_inflate_box_factor);

	auto space = std::make_shared<ControlSpace3dOctree>(box);
	space->setMaxMetric();


	double p = box.getDiameter();
	int fct = mesh->getFacesCount();

	// ======== test ACS ========
	//CS3dPtr space_f = new ControlSpace3dOctree(box);
	//space_f->setMaxMetric();
	//CS3dPtr space_01 = new ControlSpace3dOctree(box);
	//space_01->setMaxMetric();
	//CS3dPtr space_012 = new ControlSpace3dOctree(box);
	//space_012->setMaxMetric();

	//CS3dPtr space_ave_f = new ControlSpace3dOctree(box);
	//CS3dPtr space_ave_01 = new ControlSpace3dOctree(box);
	//CS3dPtr space_ave_012 = new ControlSpace3dOctree(box);

	//DataVector<int> used_faces_01(fct, -1);
	//DataVector<int> used_faces_012(fct, -1);

	//SurfaceAnalytic surface_torus("(5cm+2cm*cos(v))*cos(u)", "(7cm+1cm*cos(v))*sin(u)", "3cm*sin(v)");

//	DataVector<int> trans_table(fct);
//	for(int i = 0; i < fct; i++) trans_table.add(i);
//	trans_table.shuffle();
	// ======== test ACS ========

	bool any_changes = false;
	int approx_error_ct = 0;

	const int LAYERS_COUNT = 3;

//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "CURV_X" << "\t" << "CURV_Y" << "\t"
//		<< "DATA_LX" << "\t" << "DATA_LY" << "\t" << "DATA_A" << endl;

	// control-check
	//DataVector<ControlDataMatrix3d> cdm_f0(fct, ControlDataMatrix3d::identity);
	//DataVector< DataVector<ControlDataMatrix3d> > cdm_f1(fct, DataVector<ControlDataMatrix3d>());
	//DataVector< DataVector<ControlDataMatrix3d> > cdm_f2(fct, DataVector<ControlDataMatrix3d>());

	for(int i = 0; i < fct; i++){

// ======== test ACS ========
//	for(int ri = 0; ri < fct; ri++){
//		int i = trans_table[ri];

//		if( used_faces_01[i] >= 0) continue; // curvature already set for this face
// ======== test ACS ========

#ifdef T_DEBUG_
		if( i % 1000 == 0 ){
			ostringstream ostr;
			ostr << "Calculating approx. curvature for ACS, face " << (i+1) << " out of " << fct;
			LOG4CPLUS_INFO(MeshLog::logger_console,ostr.str());
		}
#endif // T_DEBUG_

		MeshFace* face = mesh->getFaceAt(i);
		int fpct = face->getPointCount();

		// gather (some) layer of vertices (not crossing border edges/points)
		DataVector<DPoint3d> points;
		// control-check
//		DataVector<MeshFace*> lfaces;
//		if(mesh->gatherLayeredVerticesTopological(face, points, LAYERS_COUNT, false, 6, &lfaces)) {
		if(mesh->gatherLayeredVerticesTopological(face, points, LAYERS_COUNT, false, 6)) {
			// approximate surface  -> z(1,x,y,xy,xx,yy) or z(1,x,y,xy,xx,yy, xyy, xxy, xxx, yyy)
			DPlanarQuadric pquadric;
			double max_dist = DLeastSquaresFitting::fitPlanarQuadric(points, pquadric);
			
			DBox vicinity_box;
			for(int j = 0; j < points.countInt(); j++)
				vicinity_box.addPoint(points[j]);
			double diff_threshold = vicinity_box.getDiameter() / (2*LAYERS_COUNT - 1);
			diff_threshold *= MeshGenerator3dSurface::param_local_shape_tolerance;

			double local_max_dist = pquadric.distance( face->getPoint(0)->getCoordinates() );
			for(int j = 1; j < fpct; j++){
				double d = pquadric.distance( face->getPoint(j)->getCoordinates() );
				if(d < local_max_dist) local_max_dist = d;
			}

			bool approx_error = ( local_max_dist > diff_threshold );

//			if(approx_error){ // if invalid
			if(false){ 
				MeshViewSet* set = new MeshViewSet(points.countInt(), 2*DPlanarQuadric::SKETCH_LINES, 1, 0);
				set->addFace(face);
				pquadric.createViewSetForPoints(set, points);
				set->addInfo("diff threshold", diff_threshold);
				set->addInfo("local max dist", local_max_dist);
				set->addInfo("max dist", max_dist);
				set->addInfo("diff ratio", (max_dist / diff_threshold) );
				set->addInfo("local ratio", (local_max_dist / diff_threshold) );
				SHOW_MESH("Approximation error for face during ACS creation", set);
			}

			if(approx_error){
				approx_error_ct++;
				continue;
			}
			//if(true){
			//	MeshViewSet* set = new MeshViewSet(points.countInt(), 2*DPlanarQuadric::SKETCH_LINES, 1, 0);
			//	set->addFace(face);
			//	pquadric.createViewSetForPoints(set, points);
			//	SHOW_MESH("Approx. surface for ACS creation", set);
			//}

			// create parametric surface
			SurfaceConstPtr surface(new SurfacePlanarQuadric(pquadric));

			// control-check
			ControlDataMatrix3d cdm;

			// update control space with surface curvature
			const DPoint3d fpmiddle = face->getMiddlePoint();
			SurfaceCurvature curvature2d;

			any_changes |= space->setMinControlFromSurfaceCurvature(face, surface, p, true, &cdm, &curvature2d);
			if( !cdm.isZero() ) {
				double max_curvature = std::max( curvature2d.c1, curvature2d.c2 );
				if( max_curvature > 0.0 ) face->setDoubleTag( TagExtended::TAG_MAX_CURVATURE, max_curvature );
				DVector2d v2dx = DVector2d::v_ox.turned( sin(curvature2d.angle), cos(curvature2d.angle) );
				DVector2d v2dy( v2dx.y, -v2dx.x );
				face->setCurvatureData(
					surface->getPoint(DPoint2d::zero + v2dx ) - surface->getPoint(DPoint2d::zero),
					surface->getPoint(DPoint2d::zero + v2dy ) - surface->getPoint(DPoint2d::zero),
					curvature2d.c1, curvature2d.c2 );
			}
			// ======== test ACS ========
			//ControlDataMatrix3d cdmf;

			//bool free_f_01  = (used_faces_01[i] < 0);
			//assert(free_f_01);
			//bool free_f_012 = (used_faces_012[i] < 0);

			//DataStatistics stat_max_dist;
			//double stat_fdiff_0 = -1.0;
			//DataStatistics stat_fdiff_1;
			//DataStatistics stat_fdiff_2;

			//for(int j = 0; j < lfaces.countInt(); j++){
			//	MeshFace* lface = lfaces[j];
			//	int lid = lface->getIntTag(TagExtended::TAG_LAYER_ID, -1);
			//	assert(lid >= 0);
			//	if(lid > LAYERS_COUNT) continue;
			//	const DPoint3d lfpmiddle = lface->getMiddlePoint();

			//	if(ControlSpace3dAdaptive::calculateControlFromSurfaceCurvature(lfpmiddle, &surface, p, cdm)){

			//		if(free_f_01 && (used_faces_01[j] < 0) && (lid <= 1)){
			//			used_faces_01[j] = lid;
			//			//space_ave_01->addControlPoint( lfpmiddle, cdm );
			//			any_changes |= space_01->setMinControl( lfpmiddle, cdm );
			//		}
			//		//if(free_f_012 && (used_faces_012[j] < 0) && (lid <= 2)){
			//		//	used_faces_012[j] = lid;
			//		//	space_ave_012->addControlPoint( lfpmiddle, cdm );
			//		//	space_012->setMinControl( lfpmiddle, cdm );
			//		//}
			//	}

			//	if(true){ // + directly from f(u,v)
			//		double sin_v = lfpmiddle.z / 0.03;
			//		double v1 = (sin_v <= -1.0) ? -PI/2 : (sin_v >= 1.0) ? PI/2 : asin( sin_v );
			//		double v2 = PI - v1;

			//		double cos_v1 = cos(v1);
			//		double cos_v2 = cos(v2);

			//		double cos_u1 = lfpmiddle.x / (0.05 + 0.02 * cos_v1);
			//		double cos_u2 = lfpmiddle.x / (0.05 + 0.02 * cos_v2);
			//		double sin_u1 = lfpmiddle.y / (0.07 + 0.01 * cos_v1);
			//		double sin_u2 = lfpmiddle.y / (0.07 + 0.01 * cos_v2);

			//		double u1_c1 = (cos_u1 <= -1.0) ? PI : (cos_u1 >= 1.0) ? 0 : acos( cos_u1 );
			//		double u1_c2 = -u1_c1;
			//		double u2_c1 = (cos_u2 <= -1.0) ? PI : (cos_u2 >= 1.0) ? 0 : acos( cos_u2 );
			//		double u2_c2 = -u2_c1;

			//		double u1_s1 = (sin_u1 <= -1.0) ? -PI/2 : (sin_u1 >= 1.0) ? PI/2 : asin( sin_u1 );
			//		double u1_s2 = PI - u1_s1;
			//		double u2_s1 = (sin_u2 <= -1.0) ? -PI/2 : (sin_u2 >= 1.0) ? PI/2 : asin( sin_u2 );
			//		double u2_s2 = PI - u2_s1;

			//		double v[] = { v1, v2 };
			//		double u[] = { u1_c1, u1_c2, u2_c1, u2_c2, u1_s1, u1_s2, u2_s1, u2_s2 };

			//		DPoint2d best_pt2d;
			//		double min_dist2 = LARGE_NUMBER;
			//		for(int iv = 0; iv < 2; iv++){
			//			for(int iu = 0; iu < 8; iu++){
			//				DPoint2d pt2d(u[iu], v[iv]);
			//				DPoint3d pt_uv = surface_torus.getPoint( pt2d ); 
			//				double dist2 = pt_uv.distance2( lfpmiddle );
			//				if(dist2 < min_dist2){
			//					min_dist2 = dist2;
			//					best_pt2d = pt2d;
			//				}
			//			}
			//		}
			//		stat_max_dist.add( sqrt(min_dist2) );

			//		if(ControlSpace3dAdaptive::calculateControlFromSurfaceCurvature(lfpmiddle, &surface_torus, p, cdmf, &best_pt2d)){
			//			space_f->setMinControl( lfpmiddle, cdmf );
			//			space_ave_f->addControlPoint( lfpmiddle, cdmf );

			//			double fdiff = cdm.countDifferenceRR(cdmf);
			//			if(lid == 0) stat_fdiff_0 = fdiff;
			//			else if(lid == 1) stat_fdiff_1.add(fdiff);
			//			else stat_fdiff_2.add(fdiff);
			//		}
			//	}

			//}
			//if(stat_max_dist.calculate() && stat_fdiff_1.calculate() && stat_fdiff_2.calculate()) {
			//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, i << "\t" << stat_max_dist.maximum() 
			//		<< "\t" << stat_fdiff_0 
			//		<< "\t" << stat_fdiff_1.minimum() 
			//		<< "\t" << stat_fdiff_1.average() 
			//		<< "\t" << stat_fdiff_1.maximum() 
			//		<< "\t" << stat_fdiff_2.minimum() 
			//		<< "\t" << stat_fdiff_2.average() 
			//		<< "\t" << stat_fdiff_2.maximum() 
			//		<< endl;
			//}
			// ======== test ACS ========

			// control-check
			//cdm_f0[i] = cdm;
			//for(int j = 0; j < lfaces.countInt(); j++){
			//	MeshFace* lface = lfaces[j];
			//	int lid = lface->getIntTag(TagExtended::TAG_LAYER_ID, -1);
			//	assert(lid >= 0);
			//	if(lid == 0 || lid > 2) continue;
			//	if(ControlSpace3dAdaptive::calculateControlFromSurfaceCurvature(lface->getMiddlePoint(), &surface, p, cdm)){
			//		if(lid == 1) cdm_f1[ lfaces[j]->getIndex() ].add( cdm );
			//		else cdm_f2[ lfaces[j]->getIndex() ].add( cdm );
			//	}
			//}

			// check border edges for contour curvature
			int ect = face->getEdgeCount();
			for(int j = 0; j < ect; j++){
				MeshEdge3d* edge = face->getEdge(j);
				if(!edge->isBorder()) continue;
				MeshPoint3d* ept0 = edge->getMeshPoint(0);
				MeshPoint3d* ept1 = edge->getMeshPoint(1);
				if(ept0->isBorder(TagBorder::CORNER) && ept1->isBorder(TagBorder::CORNER)) continue;

				DataVector<MeshEdge3d*> cedges;
				MeshPoint3d* cpoint = MeshContainer3dSurface::gatherBorderContourChain(edge, cedges, 2);
				int cpct = cedges.countInt()+1;
				DataVector<DPoint2d> points2d(cpct);
				//DataVector<MeshPoint3d*> points3d(cpct);
				//points3d.add(cpoint);
				points2d.add( surface->getParameters( cpoint->getCoordinates() ) );
				for(int k = 0; k < cedges.countInt(); k++){
					cpoint = cedges[k]->getOtherPoint(cpoint);
					//points3d.add(cpoint);
					points2d.add( surface->getParameters( cpoint->getCoordinates() ) );
				}
				// ... straight line -> nothing
				DLine2d line;;
				DLeastSquaresFitting::fitLine(points2d, line, false);
				int lpct = points2d.countInt();
				DPoint2d pt0 = line.getPoint( line.paramForPoint(points2d[0]) );
				DPoint2d pt1 = line.getPoint( line.paramForPoint(points2d[lpct-1]));

				double real_line_max_dist2 = 0.0;
				for(int k = 0; k < lpct; k++){
					double d = surface->getPoint( line.getPoint( line.paramForPoint(points2d[k]) ) ).distance2( 
									surface->getPoint( points2d[k] ) );
					if(d > real_line_max_dist2) real_line_max_dist2 = d;
				}

				if(real_line_max_dist2 < diff_threshold) continue;	//  close enough to line, no curvature to calculate

				if(true){
					MeshViewSet * set = new MeshViewSet;
					DPoint3d last_pt;
					for(int k = 0; k < lpct; k++){
						DPoint3d curr_pt = surface->getPoint( points2d[k] );
						set->addPoint( curr_pt );
						if(k > 0) set->addEdge( last_pt, curr_pt, 1);
						last_pt = curr_pt;
					}
					set->addInfo("chain points", lpct);
					set->addInfo("real diff", sqrt(real_line_max_dist2));
					SHOW_MESH("nonlinear border contour chain for ACS" , set);
				}

				// ... else, contour curvature for both vertices
				// TODO_CURVE
				LOG4CPLUS_WARN(MeshLog::logger_console,
					"Missing curve approximation for contour cuvrature sizing... - TODO");
			}
		}
	}

	if( approx_error_ct > 0 ){
		LOG4CPLUS_DEBUG(MeshLog::logger_console, 
			"Curvature approx. success" << ( fct - approx_error_ct )
			<< "/" << fct);
	}

	if(any_changes) space->smoothen();


	// ======== test ACS ========
	//if(any_changes) space_01->smoothen();
	//space_012->smoothen();
	//space_f->smoothen();
	//space_ave->interpolate();
	//space_ave_01->interpolate();
	//space_ave_012->interpolate();
	//space_ave_f->interpolate();

	//const int NS = 8;
	//DataStatistics stats[NS][NS];
	//string labels[NS] = { "min-f", "min", "min-01", "min-012", "ave-f", "ave", "ave-01", "ave-012" };
	//for(int i = 0; i < pct; i++){
	//	const DPoint3d& dpt = mesh->getPointAt(i)->getCoordinates();
	//	ControlDataMatrix3d cdm[NS] = {
	//		space_f->getMetricAtPoint(dpt),
	//		space->getMetricAtPoint(dpt),
	//		space_01->getMetricAtPoint(dpt),
	//		space_012->getMetricAtPoint(dpt),
	//		space_ave_f->getMetricAtPoint(dpt),
	//		space_ave->getMetricAtPoint(dpt),
	//		space_ave_01->getMetricAtPoint(dpt),
	//		space_ave_012->getMetricAtPoint(dpt) 
	//	};
	//	for(int j = 0; j < NS; j++){
	//		for(int k = 0; k < j; k++){
	//			double diffRR = cdm[j].countDifferenceRR(cdm[k]);
	//			stats[j][k].add(diffRR);
	//		}
	//	}
	//}
	//LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "=================================================================\n";
	//LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " ACS1 \t ACS2 \t MIN \t MAX \t AVE \t STDDEV\t 0.5 \t 1.0 \t 2.0 \t 3.0 \t 4.0 \t 6.0 \t + \n";
	//for(int j = 0; j < NS; j++){
	//	for(int k = 0; k < j; k++){
	//		if(stats[j][k].calculate()){
	//			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, labels[j] << "\t" << labels[k] << "\t" 
	//				<< stats[j][k].minimum() << "\t"
	//				<< stats[j][k].maximum() << "\t"
	//				<< stats[j][k].average() << "\t"
	//				<< stats[j][k].stdDev() << "\t"
	//				<< stats[j][k].getDataCountBottom(0.5) << "\t"
	//				<< stats[j][k].getDataCountInRange(0.5, 1.0) << "\t"
	//				<< stats[j][k].getDataCountInRange(1.0, 2.0) << "\t"
	//				<< stats[j][k].getDataCountInRange(2.0, 3.0) << "\t"
	//				<< stats[j][k].getDataCountInRange(3.0, 4.0) << "\t"
	//				<< stats[j][k].getDataCountInRange(4.0, 6.0) << "\t"
	//				<< stats[j][k].getDataCountTop(6.0) << "\n";
	//		}
	//	}
	//}
	//LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "=================================================================\n";
	// ======== test ACS ========

	// control-check
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "==================== CONTROL-CHECK-START =================");
	//for(int i = 0; i < fct; i++){
	//	int f1_ct = cdm_f1[i].countInt();
	//	int f2_ct = cdm_f2[i].countInt();
	//	if(cdm_f0[i].isIdentity()){
	//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "** (" << f1_ct << ") (" << f2_ct << ")");
	//		continue;
	//	}else
	//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "ok (" << f1_ct << ") " ;
	//	if(f1_ct > 0){
	//		double min_rr = cdm_f0[i].countDifferenceRR(cdm_f1[i][0]);
	//		double max_rr = min_rr;
	//		for(int j = 1; j < f1_ct; j++){
	//			double rr = cdm_f0[i].countDifferenceRR(cdm_f1[i][j]);
	//			if(rr > max_rr) max_rr = rr;
	//			if(rr < min_rr) min_rr = rr;
	//		}
	//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "\t" << min_rr << " \t" << max_rr << " \t";
	//	}
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "(" << f2_ct << ") " ;
	//	if(f2_ct > 0){
	//		double min_rr = cdm_f0[i].countDifferenceRR(cdm_f2[i][0]);
	//		double max_rr = min_rr;
	//		for(int j = 1; j < f2_ct; j++){
	//			double rr = cdm_f0[i].countDifferenceRR(cdm_f2[i][j]);
	//			if(rr > max_rr) max_rr = rr;
	//			if(rr < min_rr) min_rr = rr;
	//		}
	//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "\t" << min_rr << " \t" << max_rr << " \t";
	//	}
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;
	//}
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "==================== CONTROL-CHECK-END ===================");

	STOP_CLOCK("MG3dS:createACSFromApproxCurvatureDirect");

//	SHOW_MESH("ACS", space->getViewSet());

	return space;
}

CS3dPtr MeshGenerator3dSurface::createACSFromApproxCurvatureAverage(MeshContainer3dSurface* mesh)
{
	if(!mesh) return nullptr;
	int pct = mesh->getPointsCount();
	if(pct == 0) return nullptr;

	START_CLOCK("MG3dS:createACSFromApproxCurvatureAverage");

	DBox box;
	for(int i = 0; i < pct; i++){
		MeshPoint3d* mp = mesh->getPointAt(i);
		box.addPoint(mp->getCoordinates());
	}

	box.inflate(ControlSpace2d::param_inflate_box_factor);

	auto space = std::make_shared<ControlSpace3dOctree>(box);
	space->setMaxMetric();

	double p = box.getDiameter();
	int fct = mesh->getFacesCount();

	DataVector< ControlDataMatrix3d> cdms( fct, ControlDataMatrix3d() );

	bool any_changes = false;
	int approx_error_ct = 0;

	const int LAYERS_COUNT = 3;


	for(int i = 0; i < fct; i++){

#ifdef T_DEBUG_
		if( i % 1000 == 0 ){
			ostringstream ostr;
			ostr << "Calculating approx. curvature for ACS, face " << (i+1) << " out of " << fct;
			LOG4CPLUS_INFO(MeshLog::logger_console,ostr.str());
		}
#endif // T_DEBUG_

		MeshFace* face = mesh->getFaceAt(i);
		int fpct = face->getPointCount();

		// gather (some) layer of vertices (not crossing border edges/points)
		DataVector<DPoint3d> points;
		if(mesh->gatherLayeredVerticesTopological(face, points, LAYERS_COUNT, false, 10)) {
			// approximate surface  -> z(1,x,y,xy,xx,yy) or z(1,x,y,xy,xx,yy, xyy, xxy, xxx, yyy)
			DPlanarQuadric pquadric;
			double max_dist = DLeastSquaresFitting::fitPlanarQuadric(points, pquadric);
			
			DBox vicinity_box;
			for(int j = 0; j < points.countInt(); j++)
				vicinity_box.addPoint(points[j]);
			double diff_threshold = vicinity_box.getDiameter() / (2*LAYERS_COUNT - 1);
			diff_threshold *= MeshGenerator3dSurface::param_local_shape_tolerance;

			double local_max_dist = pquadric.distance( face->getPoint(0)->getCoordinates() );
			for(int j = 1; j < fpct; j++){
				double d = pquadric.distance( face->getPoint(j)->getCoordinates() );
				if(d < local_max_dist) local_max_dist = d;
			}

			bool approx_error = ( local_max_dist > diff_threshold );

//			if(approx_error){ // if invalid
			if(false){ // if invalid
				MeshViewSet* set = new MeshViewSet(points.countInt(), 2*DPlanarQuadric::SKETCH_LINES, 1, 0);
				set->addFace(face);
				pquadric.createViewSetForPoints(set, points);
				set->addInfo("diff threshold", diff_threshold);
				set->addInfo("local max dist", local_max_dist);
				set->addInfo("max dist", max_dist);
				set->addInfo("diff ratio", (max_dist / diff_threshold) );
				set->addInfo("local ratio", (local_max_dist / diff_threshold) );
				SHOW_MESH("Approximation error for face during ACS creation", set);
			}

			if(approx_error){
				approx_error_ct++;
				continue;
			}

			// create parametric surface
			SurfaceConstPtr surface(new SurfacePlanarQuadric(pquadric));

			// control-check
			ControlDataMatrix3d cdm;

			// update control space with surface curvature
			const DPoint3d fpmiddle = face->getMiddlePoint();
			SurfaceCurvature scurvature2d;

			//any_changes |= space->setMinControlFromSurfaceCurvature(face, &surface, p, true, &cdm, &max_curvature);

			if(ControlSpace3dAdaptive::calculateControlFromSurfaceCurvature(fpmiddle, surface, 
				p, cdm, nullptr, &scurvature2d)) 
			{
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, i << "\t" << scurvature2d.c1 << "\t" << scurvature2d.c2);
				cdms[i] = cdm;
				DVector2d v2dx = DVector2d::v_ox.turned( sin(scurvature2d.angle), cos(scurvature2d.angle) );
				DVector2d v2dy( v2dx.y, -v2dx.x );
				face->setCurvatureData(
					surface->getPoint(DPoint2d::zero + v2dx ) - surface->getPoint(DPoint2d::zero),
					surface->getPoint(DPoint2d::zero + v2dy ) - surface->getPoint(DPoint2d::zero),
					scurvature2d.c1 * scurvature2d.s1, 
					scurvature2d.c2 * scurvature2d.s2 );
				//face->setApproxPQuadric( pquadric );
				face->setBaseNormal( 
					surface->getNormalVector( 
						surface->getParameters( fpmiddle ) ) );
			}

			// check border edges for contour curvature
			int ect = face->getEdgeCount();
			for(int j = 0; j < ect; j++){
				MeshEdge3d* edge = face->getEdge(j);
				if(!edge->isBorder()) continue;
				MeshPoint3d* ept0 = edge->getMeshPoint(0);
				MeshPoint3d* ept1 = edge->getMeshPoint(1);
				if(ept0->isBorder(TagBorder::CORNER) && ept1->isBorder(TagBorder::CORNER)) continue;

				DataVector<MeshEdge3d*> cedges;
				MeshPoint3d* cpoint = MeshContainer3dSurface::gatherBorderContourChain(edge, cedges, 2);
				int cpct = cedges.countInt()+1;
				DataVector<DPoint2d> points2d(cpct);
				//DataVector<MeshPoint3d*> points3d(cpct);
				//points3d.add(cpoint);
				points2d.add( surface->getParameters( cpoint->getCoordinates() ) );
				for(int k = 0; k < cedges.countInt(); k++){
					cpoint = cedges[k]->getOtherPoint(cpoint);
					//points3d.add(cpoint);
					points2d.add( surface->getParameters( cpoint->getCoordinates() ) );
				}
				// ... straight line -> nothing
				DLine2d line;
				DLeastSquaresFitting::fitLine(points2d, line);
				int lpct = points2d.countInt();
				DPoint2d pt0 = line.getPoint( line.paramForPoint(points2d[0]) );
				DPoint2d pt1 = line.getPoint( line.paramForPoint(points2d[lpct-1]));

				double real_line_max_dist2 = 0.0;
				for(int k = 0; k < lpct; k++){
					double d = surface->getPoint( line.getPoint( line.paramForPoint(points2d[k]) ) ).distance2( 
									surface->getPoint( points2d[k] ) );
					if(d > real_line_max_dist2) real_line_max_dist2 = d;
				}

				if(real_line_max_dist2 < diff_threshold) continue;	//  close enough to line, no curvature to calculate

				if(true){
					MeshViewSet * set = new MeshViewSet;
					DPoint3d last_pt;
					for(int k = 0; k < lpct; k++){
						DPoint3d curr_pt = surface->getPoint( points2d[k] );
						set->addPoint( curr_pt );
						if(k > 0) set->addEdge( last_pt, curr_pt, 1);
						last_pt = curr_pt;
					}
					set->addInfo("chain points", lpct);
					set->addInfo("real diff", sqrt(real_line_max_dist2));
					SHOW_MESH("nonlinear border contour chain for ACS" , set);
				}

				// ... else, contour curvature for both vertices
				// TODO_CURVE
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"Missing curve approximation for contour cuvrature sizing... - TODO");
			}
		}
	}

	if( approx_error_ct > 0 ){
		LOG4CPLUS_DEBUG(MeshLog::logger_console, 
			"Curvature approx. success" << ( fct - approx_error_ct ) << "/" << fct);
	}

	if( false ) {
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < fct; i++) {
			if(! cdms[i].isZero() )
				set->addFaceWithEdges( mesh->getFaceAt(i) );
		}
		SHOW_MESH("Faces with available curvature", set );
	}
	if( false ) {
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < fct; i++) {
			MeshFace* face = mesh->getFaceAt(i);
			set->addEdges( face );
			ControlDataMatrix3d cdm = cdms[i];
			if( cdms[i].isZero() ) continue;
			const DPoint3d fmid = face->getMiddlePoint();
			set->addLabel( fmid, to_string(i) );
			double flen = 0.3 * face->getBoundingBox().getDiameter();
			set->addEdge( fmid, fmid + face->getCurvatureDirection0().normalized( flen ), 1 );
			set->addEdge( fmid, fmid + face->getCurvatureDirection1().normalized( flen ), 2 );
		}
		SHOW_MESH("Faces with main directions", set );
	}

//	ofstream f_cdm("cdm-sources-off.txt");

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Averaging and introducing metric data...");
	for(int i = 0; i < fct; i++ ){
		MeshFace* face = mesh->getFaceAt(i);
		int ect = face->getEdgeCount();
		ControlDataMatrix3d cdm = cdms[i];
		int w = cdm.isZero() ? 0 : 1;
		for(int j = 0; j < ect; j++){
			MeshEdge3d* edge = face->getEdge(j);
			if( edge->isBorder() ) continue;
			assert( edge->getFaceCount() == 2 );
			int ofid = edge->getOtherFace(face)->getIndex();
			if( !cdms[ofid].isZero() ){
				cdm += cdms[ofid];
				w++;
			}
		}
		if( w > 1 ) cdm /= w;
		if( w > 0 ) {

//			face->getMiddlePoint().storeSimple(f_cdm, '\t');
//			cdm.storeSimple(f_cdm);

			if( ect == 3 ){
				const DPoint3d& pt_a = face->getPoint(0)->getCoordinates();
				any_changes |= space->getAsAdaptive()->setMinControl( 
					ControlDataExtMatrix3dTriangle( pt_a, 
						face->getPoint(1)->getCoordinates() - pt_a, 
						face->getPoint(2)->getCoordinates() - pt_a, 0.0, cdm ) );
			}else{
				DataVector< DTriangle3d > triangles;
				if( face->splitToTriangles( triangles ) ) {
					for(int fi = 0; fi < triangles.countInt(); fi++) {
						const DTriangle3d& dtri = triangles[fi];
						any_changes |= space->getAsAdaptive()->setMinControl( 
							ControlDataExtMatrix3dTriangle( dtri.pt_a, 
								dtri.pt_b - dtri.pt_a, dtri.pt_c - dtri.pt_a, 0.0, cdm ) );
					}
				}
			}
		}
	}

	if(any_changes) space->smoothen();

	STOP_CLOCK("MG3dS:createACSFromApproxCurvatureAverage");

	return space;
}

/// Split the given edge (normal or boundary) by inserting a new point in the middle
MeshPoint3d* MeshGenerator3dSurface::splitEdgeWLT(Metric3dContext& mc, 
					MeshContainer3dSurface* mesh, MeshEdge3d* act_edge,
					DataContainer<MeshEdge3d::ActiveEdge> * act_edges, TagExtended::TagType act_tag)
{
	MeshPoint3d* p0 = act_edge->getMeshPoint(0);
	MeshPoint3d* p1 = act_edge->getMeshPoint(1);
	MeshPoint3d* point = nullptr;
	Curve3dConstPtr local_curve;
	SurfaceConstPtr local_surface;

	// rearange the mesh (i.e. "split" faces)
	int fct = act_edge->getFaceCount();
	DataVector<MeshFace*> act_faces(fct);
	for(int i = 0; i < fct; i++)
		act_faces.add(act_edge->getFaceAt(i));

	DPoint2d local_surface_param;

	if( act_edge->isBorder() ){ // new point following local curve
		local_curve = act_edge->getLocalCurve();
		assert( local_curve != nullptr );
		double t0 = p0->getLocalCurveParam( local_curve );
		double t1 = p1->getLocalCurveParam( local_curve );
		double t = (t0+t1)*0.5;
		point = new MeshPoint3d( local_curve, t );
		DVector3d pn = p0->getBaseNormal() + p1->getBaseNormal();
		point->setBaseNormal( pn.normalized() );
		//local_surface = act_faces[0]->getValidLocalSurface();
		//assert( local_surface != nullptr );
		//DPoint2d local_surface_param;
		//local_surface_param.add( p0->getLocalSurfaceParam( local_surface ), 0.5);
		//local_surface_param.add( p1->getLocalSurfaceParam( local_surface ), 0.5);
		//local_surface_param = local_surface->getParametersNear( point->getCoordinates(), local_surface_param );
		//point->setLocalSurface( local_surface, local_surface_param );
	}else{
		assert( fct == 2 );
		DataVector<MeshPoint3d*> mpoints(4);
		MeshFace::getMeshPointsForFaces( act_faces, mpoints ); 
		local_surface = getLocalSurfaceForPoints( mc, mesh, act_faces, mpoints);
		assert( local_surface != nullptr );
		local_surface_param.add( p0->getLocalSurfaceParam( local_surface ), 0.5);
		local_surface_param.add( p1->getLocalSurfaceParam( local_surface ), 0.5);
		point = new MeshPoint3d( local_surface, local_surface_param );
		point->setBaseNormal( local_surface->getNormalVector( local_surface_param ) );
	}

	mesh->addMeshPoint(point);
	// set border info
	char bflags = act_edge->getBorderFlags();
	point->setBorder( bflags ); // ie. ridge or normal boundary

	for(int i = 0; i < fct; i++){
		MeshFace* face = act_faces[i];
		MeshFace* face_clone = face->clone();
		face->switchPointsWithEdges(p1, point);
		face_clone->switchPointsWithEdges(p0, point);
		mesh->addMeshFace(face_clone);
	}

	if(bflags != 0){
		assert( local_curve != nullptr );
		delete act_edge;
		act_edge = p0->getEdgeToPoint(point);
		assert(act_edge);
		act_edge->setBorderFlags(bflags);
		act_edge->setLocalCurve( local_curve );
		act_edge = p1->getEdgeToPoint(point);
		assert(act_edge);
		act_edge->setBorderFlags(bflags);
		act_edge->setLocalCurve( local_curve );
	}

	// no more act_edge

	bool show_case = false;

	if(show_case){
		SHOW_MESH("After split", mesh->getDebugViewSet(point));
	}

	// project new point onto associated surface + set proper surface !
	if(point->isBorder())
		MeshGenerator3dSurface::moveBoundaryPointByLaplace( mesh, mc, point );
	else
		MeshGenerator3dSurface::movePointByLaplace( mesh, mc, point);

	bool any_swap = false;
	act_faces.clear();
	if(point->adjacentFaces(act_faces)){
		for(int i = 0; i < act_faces.countInt(); i++){
			MeshTriangle3d* triangle = (MeshTriangle3d*) act_faces[i];
			assert(triangle->getType() == FACE_TRIANGLE);
			int j = 0;
			do{
				act_edge = triangle->getEdge(j++);
			}while( act_edge->incidentTo(point) );
			MeshEdge3d::ActiveEdge* act = nullptr;
			MeshPoint3d *act_pt0, *act_pt1;
			if(act_edges){
				act = (MeshEdge3d::ActiveEdge*) act_edge->getPtrTag(act_tag);
				if(act){
					act_edge->removeTag(act_tag);
					act_pt0 = act_edge->getMeshPoint(0);
					act_pt1 = act_edge->getMeshPoint(1);
				}
			}
			if(MeshGenerator3dSurface::swap22(mc, mesh, act_edge) != nullptr){
				if(act) act_edges->removeDataItemValue( act );
				any_swap = true;
			}else{
				if(act){
					act->edge = act_edge = act_pt0->getEdgeToPoint(act_pt1); // update ...
					act_edge->setPtrTag(act_tag, act);
				}
			}
		}
	}

	if(any_swap){
//			if(show_case){
		//if(true){
		//	SHOW_MESH_NORESET("After split & swap", mesh->getDebugViewSet(point));
		//}
		if(point->isBorder())
			MeshGenerator3dSurface::moveBoundaryPointByLaplace( mesh, mc, point );
		else
			MeshGenerator3dSurface::movePointByLaplace( mesh, mc, point );
	}

	if(act_edges){
		// add new edges to active_edges
		for(int i = 0; i < point->getRank(); i++){
			MeshEdge3d* edge = point->getEdge(i);

			mc.countMetricAtPoints(point, edge->getOtherPoint(point));
			double len = edge->getLength(mc);
			bool too_long = len > (edge->isBorder() ? MAX_BEDGE_LEN : MAX_EDGE_LEN);
			if(too_long){
				MeshEdge3d::ActiveEdge* act = new MeshEdge3d::ActiveEdge(edge, 1.0/len);
				act_edges->addDataItem(act);
				edge->setPtrTag(act_tag, act);
			}
		}
	}

	if(show_case){
		SHOW_MESH_NORESET("After split & smoothing", mesh->getDebugViewSet(point));
	}

	return point;
}

/// Split the given edge (normal or boundary) by inserting a new point in the middle
MeshPoint3d* MeshGenerator3dSurface::splitEdgeSimple( MeshContainer3dSurface* mesh, MeshEdge3d* act_edge )
{
	assert( ! act_edge->hasLocalCurve() );

	MeshPoint3d* p0 = act_edge->getMeshPoint(0);
	MeshPoint3d* p1 = act_edge->getMeshPoint(1);
	MeshPoint3d* point = nullptr;

	assert( ! p0->hasLocalSurface() );
	assert( ! p1->hasLocalSurface() );

	DPoint3d new_pt_coord(p0->getCoordinates(), p1->getCoordinates(), 0.5);
	point = new MeshPoint3d( new_pt_coord );
	point->setBaseNormal( (( p0->getBaseNormal() + p1->getBaseNormal() ) * 0.5).normalized() );

	mesh->addMeshPoint(point);
	// set border info
	char bflags = act_edge->getBorderFlags();
	point->setBorder( bflags ); // ie. ridge or normal boundary

	// rearange the mesh (i.e. "split" faces)
	int fct = act_edge->getFaceCount();
	DataVector<MeshFace*> act_faces(fct);
	for(int i = 0; i < fct; i++)
		act_faces.add(act_edge->getFaceAt(i));

	for(int i = 0; i < fct; i++){
		MeshFace* face = act_faces[i];
		MeshFace* face_clone = face->clone();
		face->switchPointsWithEdges(p1, point);
		face_clone->switchPointsWithEdges(p0, point);
		mesh->addMeshFace(face_clone);
	}

	if(bflags != 0){
		delete act_edge;
		act_edge = p0->getEdgeToPoint(point);
		assert(act_edge);
		act_edge->setBorderFlags(bflags);
		act_edge = p1->getEdgeToPoint(point);
		assert(act_edge);
		act_edge->setBorderFlags(bflags);
	}

	if(false){
		SHOW_MESH("After split", mesh->getDebugViewSet(point));
	}

	return point;
}
//#define SHOW_ACTIVE_EDGES

/// Insert nodes at boundary, for too long edges, return number of inserted points
int MeshGenerator3dSurface::splitEdges(Metric3dContext& mc, MeshContainer3dSurface* mesh)
{

	LOG4CPLUS_INFO(MeshLog::logger_console,"Mesh transformation: splitting edges ...");

	START_CLOCK("MeshGenerator3dSurface::splitEdges");
	// gather edges

#ifdef SHOW_ACTIVE_EDGES
	MeshViewSet* set = new MeshViewSet;
#endif

	DataContainer<MeshEdge3d::ActiveEdge> active_edges(mesh->getPointsCount(), true);
	for(IteratorEdge3d it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();

		MeshPoint3d* p0 = edge->getMeshPoint(0);
		MeshPoint3d* p1 = edge->getMeshPoint(1);

		mc.countMetricAtPoints(p0, p1);
		double len = edge->getLength(mc);
		bool too_long = len > (edge->isBorder() ? MAX_BEDGE_LEN : MAX_EDGE_LEN);
		if(too_long){
			MeshEdge3d::ActiveEdge* act = new MeshEdge3d::ActiveEdge(edge, 1.0/len);
			active_edges.addDataItem(act);
			edge->setPtrTag(TagExtended::TAG_ACTIVE, act);
#ifdef SHOW_ACTIVE_EDGES
			set->addEdge(edge);
#endif
		}else
			edge->removeTag(TagExtended::TAG_ACTIVE);
	}

#ifdef SHOW_ACTIVE_EDGES
	if(active_edges.countInt() > 0)
		SHOW_MESH("Edges to split", set);
	else delete set;
#endif


	int inserted_count = 0;

	while(active_edges.countInt() > 0){
		MeshEdge3d::ActiveEdge* act = active_edges.removeDataItem(0);
		MeshEdge3d* act_edge = act->edge;

		inserted_count++;

		bool show_case = false;
		if(show_case){
			MeshPoint3d* p0 = act_edge->getMeshPoint(0);
			MeshPoint3d* p1 = act_edge->getMeshPoint(1);
			ostringstream ostr;
			MeshViewSet* set = mesh->getDebugViewSet(p0, p1, 4.0);
			set->addInfo("inserted count", inserted_count);
			set->addInfo("edge length", (1.0/act->len) );
			set->addInfo("p0", p0->getIndex() );
			set->addInfo("p1", p1->getIndex() );
			set->addInfo("edge type", (act_edge->isBorder(TagBorder::RIDGE) ? "ridge" : "normal"));
			SHOW_MESH("Selected edge to split", set);
		}

		act_edge->removeTag(TagExtended::TAG_ACTIVE);
		delete act;

		splitEdgeWLT(mc, mesh, act_edge, &active_edges, TagExtended::TAG_ACTIVE);
	}

	STOP_CLOCK("MeshGenerator3dSurface::splitEdges");

	return inserted_count;
}


MeshEdge3d* MeshGenerator3dSurface::swap22(Metric3dContext& mc, MeshContainer3dSurface* mesh, MeshEdge3d *edge)
{
	if(edge->isBorder() || edge->getFaceCount() != 2) return nullptr;

	DataVector<MeshFace*> mfaces(2);
	mfaces.add( edge->getFaceAt(0) );
	mfaces.add( edge->getFaceAt(1) );

	DataVector<MeshPoint3d*> mpoints(4);
	mpoints.add( edge->getMeshPoint(0) );
	mpoints.add( edge->getMeshPoint(1) );
	mpoints.add( mfaces[0]->getOtherPoint( mpoints[0], mpoints[1] ) );
	mpoints.add( mfaces[1]->getOtherPoint( mpoints[0], mpoints[1] ) );

	if(mpoints[2]->getEdgeToPoint(mpoints[3])){
		//SHOW_MESH("crossing edge", mesh->getDebugViewSet(epts[0], epts[1]));
		return nullptr;
	}

	SurfaceConstPtr surface = getLocalSurfaceForPoints( mc, mesh, mfaces, mpoints );
	if(((MeshTriangle3d*)mfaces[0])->swapWithNeighbour(mc, surface, mfaces[0]->getEdgeIndex(edge), true, true))
		return mpoints[2]->getEdgeToPoint( mpoints[3] );
	else return nullptr;
}

/// Rearrange boundary nodes for optimum length of border edges
void MeshGenerator3dSurface::optimizeBoundaryEdges(Metric3dContext& mc, MeshContainer3dSurface* mesh)
{
	assert(mesh);

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Mesh transformation: optimizing boundary edges ...");
	// find all border chains

	int pct = mesh->getPointsCount();

	struct ChainData{
		DataVector<MeshEdge3d*> cedges;
		MeshPoint3d* cpoint;
	};

	DataSimpleList< ChainData > chains;
	DataHashTable< MeshEdge3d* > hash_bedges(pct, nullptr);

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if(!point->isBorder()) continue;
		for(int j = 0; j < point->getRank(); j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(!edge->isBorder() || hash_bedges.contains(edge) ) continue;
			ChainData chain;
			chain.cpoint = MeshContainer3dSurface::gatherBorderContourChain(edge, chain.cedges, -1, false);

			// split chain into two - starting from both ends ?

			chains.append( chain );
			for(int k = 0; k < chain.cedges.countInt(); k++){
				hash_bedges.insert( chain.cedges[k] );
			}
		}
	}

	bool show_case = false;

	// foreach chain -> check & rearrange
	while( chains.notEmpty() ){
		ChainData cp = chains.removeFirst();

		double total_mlen = 0.0;

		int ect = cp.cedges.countInt();
		for(int i = 0; i < ect; i++)
			total_mlen += cp.cedges[i]->getLength(mc, true);

		int opt_ect = (int)( total_mlen + 0.5 );
		if(opt_ect < 2) continue;

		double ave_mlen = total_mlen / opt_ect;

		// prepare queue of edges
		DataSimpleList<MeshEdge3d*> edge_list;
		for(int i = 0; i < ect; i++)
			edge_list.append(cp.cedges[i]);

		// process
		auto chain_curve = edge_list.getFirst()->getLocalCurve();
		MeshPoint3d* epoint = cp.cpoint;
		double t0 = epoint->getLocalCurveParam( chain_curve );
		while( edge_list.notEmpty() ){
			MeshEdge3d* edge = edge_list.removeFirst();
			assert( edge->getLocalCurve() == chain_curve );
			MeshPoint3d* next_point = edge->getOtherPoint(epoint);
			double t1 = next_point->getLocalCurveParam( chain_curve );
			double mlen = edge->getLength(mc, true);

			if(show_case){
				ostringstream ostr;
				ostr << "Boundary edge to optimize, ave_mlen= " << ave_mlen << ", mlen= " << mlen;
				SHOW_MESH(ostr.str(), mesh->getDebugViewSet(epoint, next_point, 4.0));
			}

			//==================================================================================
			//for(int i = 0; i < mesh->getPointsCount(); i++){
			//	MeshPoint3d* point = mesh->getPointAt(i);
			//	if(! point->isBorder() ) continue;
			//	for(int j = 0; j < point->getRank(); j++){
			//		MeshEdge3d* edge = point->getEdge(j);
			//		if( edge->getPointIndex(point) != 0 || !edge->isBorder() ) continue;

			//		assert( edge->hasLocalCurve() );
			//	}
			//}
			//==================================================================================

			if(mlen > 1.8*ave_mlen) { // split
				MeshPoint3d* new_point = MeshGenerator3dSurface::splitEdgeWLT(mc, mesh, edge);
				// insert two new edges into queue
				edge = new_point->getEdgeToPoint(next_point);
				assert(edge);
				edge_list.insert(edge);
				edge = epoint->getEdgeToPoint(new_point);
				assert(edge);
				edge_list.insert(edge);

				//==================================================================================
				//for(int i = 0; i < mesh->getPointsCount(); i++){
				//	MeshPoint3d* point = mesh->getPointAt(i);
				//	if(! point->isBorder() ) continue;
				//	for(int j = 0; j < point->getRank(); j++){
				//		MeshEdge3d* edge = point->getEdge(j);
				//		if( edge->getPointIndex(point) != 0 || !edge->isBorder() ) continue;

				//		assert( edge->hasLocalCurve() );
				//	}
				//}
				//==================================================================================

				continue;
			}
			
			MeshEdge3d* next_edge = (edge_list.empty() ? nullptr : edge_list.getFirst());
			double next_mlen = (next_edge ? next_edge->getLength(mc, true) : 0.0);

			if( next_edge && ((mlen + next_mlen) < 1.2*ave_mlen) ){ // join
				edge_list.removeFirst(); // remove next_edge;
				// storing info about next_next edge, since collapsing may replace that edge...
				MeshPoint3d* nn_point = next_edge->getOtherPoint(next_point);
				char pbf = epoint->getBorderFlags();
				epoint->setBorderFlags( TagBorder::CORNER );
				if(MeshGenerator3dSurface::collapseBoundaryEdge(mc, mesh, edge, false, false, false)){
					edge = epoint->getEdgeToPoint(nn_point);
					assert(edge);
					edge_list.insert(edge); // once more, another collapsing might be required...
					epoint->setBorder(pbf);
					//==================================================================================
					//for(int i = 0; i < mesh->getPointsCount(); i++){
					//	MeshPoint3d* point = mesh->getPointAt(i);
					//	if(! point->isBorder() ) continue;
					//	for(int j = 0; j < point->getRank(); j++){
					//		MeshEdge3d* edge = point->getEdge(j);
					//		if( edge->getPointIndex(point) != 0 || !edge->isBorder() ) continue;

					//		assert( edge->hasLocalCurve() );
					//	}
					//}
					//==================================================================================

					continue;
				}else{
					edge_list.insert(next_edge);
					epoint->setBorder(pbf);
				}
			}

			if( ! next_point->isBorder( TagBorder::CORNER ) ) {
				// move next_point accordingly
				double k = ave_mlen/mlen;
				double t = t0 * (1.0-k) + t1 * k;
				//mc.countMetricAtPoint( epoint->getCoordinates() );
				//DPoint3d sdpt = next_point->getCoordinatesWithLocalShape(mc, 
				//	DPoint3d(epoint->getCoordinates(), next_point->getCoordinates(), ) );
				next_point->tryMovingPoint( chain_curve, t );			
			}

			// ...
			if(show_case){
				ostringstream ostr;
				ostr << "Boundary edge - optimized, ave_mlen= " << ave_mlen 
						<< ", mlen= " << edge->getLength(mc, true);
				SHOW_MESH_NORESET(ostr.str(), mesh->getDebugViewSet(epoint, next_point, 4.0));
			}

			if(edge_list.notEmpty()) {
				epoint = next_point;
				t0 = t1;
			}
		}
		if(!epoint->isBorder( TagBorder::CORNER ))
			moveBoundaryPointByLaplaceForVariableMetric(mesh, mc, epoint);
	}
}

/// Try improving quality of the given face
void MeshGenerator3dSurface::improveFace(Metric3dContext & mc, MeshContainer3dSurface* mesh, MeshFace* face)
{
	assert(face);

	int fpct = face->getPointCount();
	SurfaceConstPtr surface = face->getOptLocalSurface();
	assert( surface != nullptr );
	if(surface){
		for(int i = 0; i < fpct; i++){
			MeshPoint3d* point = face->getPoint(i);
			if(point->isBorder())
				moveBoundaryPointByLaplace( mesh, mc, point );
			else movePointByLaplace(mesh, mc, point);
		}
	}else{

		if(true){
			MeshViewSet* set = new MeshViewSet();
			DataVector<MeshFace*> pfaces(100);
			for(int i = 0; i < fpct; i++){
				MeshPoint3d* point = face->getPoint(i);
				point->adjacentFaces(pfaces);
				set->addPoint(point);
			}
			for(int i = 0; i < pfaces.countInt(); i++){
				set->addFaceWithEdges(pfaces[i], (pfaces[i] == face) ? 2 : 1);
			}
			SHOW_MESH("MG3dS::improveFace before", set);
		}

		for(int i = 0; i < fpct; i++){
			MeshPoint3d* point = face->getPoint(i);
			if(point->isBorder()) continue;
			DVector3d pnormal;
			int prank = point->getRank();
			DataVector<MeshFace*> pfaces(prank);
			DataVector<DVector3d> pnormals(prank);
			if(!point->adjacentFaces(pfaces)) continue;
			int pfct = pfaces.countInt();
			DVector3d n;
			for(int j = 0; j < pfct; j++){
				if(pfaces[j]->checkAndGetNormalVector(n)){
					pnormals.add(n);
					pnormal += n;
				}
			}
			if(pnormals.empty()) continue;
			n /= pnormals.countInt();
			double min_sp = 1.0;
			for(int j = 0; j < pnormals.countInt(); j++){
				double sp = pnormal.scalarProduct( pnormals[j] );
				if(sp < min_sp) min_sp = sp;
			}
			if( min_sp < MeshGenerator3dSurface::param_sharp_edge_threshold ) 
				continue; // too large normal difference from average ...
			DPoint3d new_pt;
			double pd = 1.0 / prank;
			for(int j = 0; j < prank; j++)
				new_pt.add(point->getEdge(j)->getOtherPoint(point)->getCoordinates(), pd);
			point->tryMovingPoint(mc, new_pt);
		}

		if(true){
			MeshViewSet* set = new MeshViewSet();
			DataVector<MeshFace*> pfaces(100);
			for(int i = 0; i < fpct; i++){
				MeshPoint3d* point = face->getPoint(i);
				point->adjacentFaces(pfaces);
				set->addPoint(point);
			}
			for(int i = 0; i < pfaces.countInt(); i++){
				set->addFaceWithEdges(pfaces[i], (pfaces[i] == face) ? 2 : 1);
			}
			SHOW_MESH("MG3dS::improveFace - after", set);
		}
	}
}

/// Calculate base normals for surface faces in the given subdomain
bool MeshGenerator3dSurface::updateNormals(MeshContainer3dSurface* mesh, 
	std::shared_ptr<const DataVector<MeshFace*>> sub_faces)
{
	DVector3d vn;
	DataVector< MeshFace* > missing_fnormals;
	int sfct = sub_faces->countInt();
	for(int i = 0; i < sfct; i++) {
		MeshFace* face = sub_faces->get(i);
		if( face->checkAndGetNormalVector( vn ) ) {
			face->setBaseNormal(vn);
		}else{
			missing_fnormals.add( face );
		}
	}
	if( missing_fnormals.notEmpty() ) // missing normals for faces or points -> leaving previous values of base normals
		LOG4CPLUS_WARN(MeshLog::logger_console,  
			"MG3dS::updateNormals: missing f-normals count = " << missing_fnormals.countInt() );

	int pct = mesh->getPointsCount();
	DataVector<int> pw(pct, 0);
	DataVector<DVector3d> pn(pct, DVector3d::zero);

	// faces -> vertices for average normals
	for(int i = 0; i < sfct; i++){
		MeshFace* face = sub_faces->get(i);
		if( face->hasBaseNormal() ) {
			const DVector3d& vni = face->getBaseNormal();
			assert(!vni.isZero());
			int fpct = face->getPointCount();
			for(int j = 0; j < fpct; j++){
				MeshPoint3d* point = face->getPoint(j);
				int idx = point->getIndex();
				pw[idx]++;
				pn[idx] += vni;
			}
		}
	}

	for(int i = 0; i < pct; i++){
		if(pw[i] > 0){
			if(pw[i] > 1) pn[i].normalize();
			mesh->getPointAt(i)->setBaseNormal(pn[i]);
		}
	}

	// additional averaging of base normals for faces ?
	if( true ) {
		for(int i = 0; i < sfct; i++) {
			MeshFace* face = sub_faces->get(i);
			vn = DVector3d::zero;
			int fpct = face->getPointCount();
			for(int j = 0; j < fpct; j++)
				vn += face->getPoint(j)->getBaseNormal();
			if(!vn.isZero())
				face->setBaseNormal( vn.normalized() );
		}
	}

	return true;
}

/// Calculate base normals for surface faces
bool MeshGenerator3dSurface::calculateNormals(MeshContainer3dSurface* mesh)
{
	assert(mesh);
	int fct = mesh->getFacesCount();
	DVector3d vn;
	DataVector<int> missing_fnormals(fct);
	DataVector<bool> available_fnormal(fct, true);

	// ... normals for faces
	for(int i = 0; i < fct; i++){
		MeshFace* face = mesh->getFaceAt(i);
		if(face->checkAndGetNormalVector(vn)){
			face->setBaseNormal(vn);
		}else{
			available_fnormal[i] = false;
			missing_fnormals.add(i);
		}
	}

	// ... (average) normals for points
	int pct = mesh->getPointsCount();
	DataVector<int> pw(pct, 0);
	DataVector<DVector3d> pn(pct, DVector3d::zero);
	// faces -> vertices for average normals
	for(int i = 0; i < fct; i++){
		if(available_fnormal[i]){
			MeshFace* face = mesh->getFaceAt(i);
			const DVector3d& vni = face->getBaseNormal();
			assert(!vni.isZero());
			int fpct = face->getPointCount();
			for(int j = 0; j < fpct; j++){
				MeshPoint3d* point = face->getPoint(j);
				int idx = point->getIndex();
				pw[idx]++;
				pn[idx] += vni;
			}
		}
	}
	DataVector<int> missing_pnormals(pct);
	for(int i = 0; i < pct; i++){
		if( !pn[i].isZero() ){
			if( pw[i] > 1 ) pn[i].normalize();
			mesh->getPointAt(i)->setBaseNormal(pn[i]);
		}else
			missing_pnormals.add(i);
	}

	// ... fill missing normals for faces/points
	if(missing_fnormals.notEmpty()){
		bool any_change = true;
		while( missing_fnormals.notEmpty() && any_change){
			any_change = false;
			int i = 0;
			while( i < missing_fnormals.countInt() ){
				MeshFace* face = mesh->getFaceAt(missing_fnormals[i]);
				vn = DVector3d::zero;
				int count = 0;
				int fpct = face->getPointCount();
				for(int j = 0; j < fpct; j++){
					int idx = face->getPoint(j)->getIndex();
					if(pw[idx] > 0) {
						vn += pn[idx];
						count++;
					}
				}
				if(count > 0 && !vn.isZero()){
					vn.normalize();
					face->setBaseNormal(vn);
					missing_fnormals.removeAt(i);
					any_change = true;
				}else i++;
			}
			i = 0;
			while( i < missing_pnormals.countInt() ){
				MeshPoint3d* point = mesh->getPointAt( missing_pnormals[i] );
				vn = DVector3d::zero;
				int count = 0;
				DataVector<MeshFace*> pfaces;
				if(point->adjacentFaces(pfaces)){
					for(int j = 0; j < pfaces.countInt(); j++){
						if(available_fnormal[pfaces[j]->getIndex()]){
							vn += pfaces[j]->getBaseNormal();
							count++;
						}
					}
				}
				if(count > 0 && !vn.isZero() ){
					vn.normalize();
					point->setBaseNormal(vn);
					missing_pnormals.removeAt(i);
					any_change = true;
				}else i++;
			}
		}
	}

	return missing_fnormals.empty();
}

/// Check and fix faces using topological orientation info
bool MeshGenerator3dSurface::checkAndFixInvertedFacesTopological(MeshContainer3dSurface* surface_mesh)
{
	assert(surface_mesh);
	int fct = surface_mesh->getFacesCount();

	DataVector<int> other_count(fct, 0);
	DataVector<int> bad_orient_count(fct, 0);
	DataVector<int> bad_faces(fct);

	// check
	for(int i = 0; i < fct; i++){
		MeshFace* face = surface_mesh->getFaceAt(i);
		int ect = face->getEdgeCount();
		for(int j = 0; j < ect; j++){
			MeshEdge3d* edge = face->getEdge(j);
			if(edge->getFaceCount() != 2) continue;	// not interested in 1 or 3+ cases
			other_count[i]++;
			MeshFace* other_face = edge->getOtherFace(face);
			int other_i = other_face->getIndex();
			if(other_i < i) continue; // already checked
			if(! face->sameOrientation(other_face)){
				bad_orient_count[i]++;
				bad_orient_count[other_i]++;
			}
		}
		if(bad_orient_count[i] > 0) bad_faces.add(i);
	}

	if( bad_faces.empty() ) return true;

	// fix 
	int iteration_limit = 2 * bad_faces.countInt();
	int switch_counter = 0;
	while( bad_faces.notEmpty() && (--iteration_limit > 0)){

		int best_i = -1;
		int max_gain = -1000;
		for(int i = 0; i < bad_faces.countInt(); i++){
			int fi = bad_faces[i];
			int boc = bad_orient_count[fi];
			if(boc == 0) {
				bad_faces.removeAt(i);
				i--;
				continue;
			}
			int gain = 2*boc - other_count[fi];
			if(gain > max_gain){
				max_gain = gain;
				best_i = i;
			}
		}
		if(bad_faces.empty()) break;

		if(max_gain < 0){
			SHOW_MESH("OFF: max_gain < 0", surface_mesh->getViewSet());
		}
		assert(max_gain >= 0);

		int fi = bad_faces[best_i];
		MeshFace* face = surface_mesh->getFaceAt(fi);
		face->switchOrientation();
		switch_counter++;
		bad_orient_count[fi] = other_count[fi] - bad_orient_count[fi];
		int ect = face->getEdgeCount();
		for(int i = 0; i < ect; i++){
			MeshEdge3d* edge = face->getEdge(i);
			if(edge->getFaceCount() != 2) continue;	// not interested in 1 or 3+ cases
			MeshFace* other_face = edge->getOtherFace(face);
			int other_i = other_face->getIndex();
			if(face->sameOrientation(other_face)){
				bad_orient_count[other_i]--;
				assert( bad_orient_count[other_i] >= 0 );
			}else{
				bad_orient_count[other_i]++;
				if(bad_orient_count[other_i] == 1){
					//assert(!bad_faces.contains(other_i));
					bad_faces.addIfNew(other_i);
				}
			}
		}
	}

	if(switch_counter > 0){
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Switched orientation for faces count = " << switch_counter);
	}

	return bad_faces.empty();
}

/// Search for local reparameterization surfaces (for whole mesh, or only faces with the given tag)
void MeshGenerator3dSurface::identifyLocalSurfaces(Metric3dContext& mc, MeshContainer3dSurface* mesh, 
		TagExtended::TagType tag_type, int tag_value)
{
	DataVector< std::shared_ptr<DataVector<MeshFace*>> > sub_domains(100);
	int sub_ct = splitToSubdomains(mesh, sub_domains, 0, tag_type, tag_value);

	sub_ct = fixOpenBorders(mesh, sub_domains, sub_ct);

	//if(true){
	//	int fct = mesh->getFacesCount();
	//	int stfct = 0;
	//	for(int i = 0; i < sub_domains.countInt(); i++)
	//		stfct += sub_domains[i]->countInt();
	//	assert( fct == stfct );
	//	for(int i = 0; i < sub_domains.countInt(); i++){
	//		auto sd = sub_domains[i];
	//		for(int j = 0; j < sd->countInt(); j++){
	//			int fid = sd->get(j)->getIntTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1);
	//			assert( fid == i );
	//		}
	//	}
	//}

	int surf_count = 0;
	for(int i = 0; i < sub_ct; i++){
		//identifyLocalSurfacesSubdivide(mc, mesh, sub_domains[i], TagExtended::TAG_LOCAL_SURFACE_DOMAIN, tolerance);
		updateNormals( mesh, sub_domains[i] );

		switch(param_local_surface_identification_method){
		case LSI_GREEDY:
			surf_count += identifyLocalSurfacesGreedy(mc, mesh, sub_domains[i] );
			break;
		case LSI_HIGH_CURVATURE_FIRST:
			surf_count += identifyLocalSurfacesHighCurvatureFirst(mc, mesh, sub_domains[i] );
			break;
		case LSI_VIA_NORMALS:
			surf_count += identifyLocalSurfacesViaNormals(mc, mesh, sub_domains[i] );
			break;
		default:
			assert(false);
			break;
		}
	}

	string label = "Identified local surfaces";
	switch(param_local_surface_identification_method){
	case LSI_GREEDY: label += " (GREEDY)"; break;
	case LSI_HIGH_CURVATURE_FIRST: label += " (HIGH-CURV)"; break;
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_console, label << " -> " << surf_count);

//	mesh->checkSurfaceSetCounters( mc, false );
	fixLocalSurfaceSeams( mc, mesh );
//	optimizeLocalSurfaces( mc, mesh );
//	mesh->checkSurfaceSetCounters( mc, false );

//	int extra_surf_count = complementLocalSurfaces( mc, mesh, TagExtended::TAG_LOCAL_SURFACE_DOMAIN, tolerance );
//	LOG4CPLUS_INFO(MeshLog::logger_console, "Extra local surfaces", extra_surf_count);

	identifyLocalCurves( mc, mesh, param_local_shape_tolerance );

//#ifdef REMESH_STATS
//	mesh->checkLocalSurfaces( mc );
//#endif
}

/// try to split domains with open borders
int MeshGenerator3dSurface::fixOpenBorders(MeshContainer3dSurface* mesh, 
		DataVector< std::shared_ptr<DataVector<MeshFace*>> > & sub_domains, int sub_tag_init)
{
	const bool visualisation_on = false;

	int base_sub_ct = sub_domains.countInt();
	for(int si = 0; si < base_sub_ct; si++){
		auto sub_faces = sub_domains[si];
		int sfct = sub_faces->countInt();
		// does this subdomain has any open border edge?
		DataHashTable<int> hpoints(2*sfct, -1);
		DataVector<MeshPoint3d*> open_border_points(100);
		for(int i = 0; i < sfct; i++){
			MeshFace* face = sub_faces->get(i);
			int fpct = face->getPointCount();
			for(int j = 0; j < fpct; j++){
				MeshPoint3d* point = face->getPoint(j);
				if( hpoints.insert(point->getIndex()) && point->isBorder(TagBorder::CORNER) && (point->getBorderEdgesCount() == 1))
					open_border_points.add(point);
			}
		}

		int spct = hpoints.countInt();
		DataHashTable<MeshPoint3d*> hsb_points(2*spct, nullptr);
		DataHashTable<MeshEdge3d*>  hsb_edges(2*spct, nullptr);

		if( open_border_points.notEmpty() ){
			if(visualisation_on){
				MeshViewSet* set = new MeshViewSet();
				for(int i = 0; i < sfct; i++)
					set->addFaceWithEdgesAndPoints(sub_faces->get(i));
				SHOW_MESH("sub-domain with open border", set);
			}

			//if(true) {
			//	MeshViewSet* set = new MeshViewSet();
			//	for(int i = 0; i < sfct; i++) {
			//		MeshFace* face = sub_faces->get(i);
			//		int ect = face->getEdgeCount();
			//		for(int j = 0; j < ect; j++){
			//			MeshEdge3d* edge = face->getEdge(j);
			//			if(edge->isBorder()) 
			//				set->addEdge( edge, 1 );
			//			else if( hsb_edges.contains( edge ) )
			//				set->addEdge( edge, 2 );
			//		}
			//	}
			//	SHOW_MESH("borders / sub-borders for local sub-domain", set);
			//}

			// join single border points ...
			while( open_border_points.notEmpty() ) {
				//  - identify start points, forbiden points (border, too close), goal points (border)
				MeshPoint3d* start_point = open_border_points.removeLast();
				DataHashTableKeyValue<MeshPoint3d*, double> hpw(2*spct, nullptr);
				DataVector<MeshPoint3d*> forbidden_points(spct);
				DataVector<int> forbidden_layers(spct);
				const int MAX_DIST = 10;
				forbidden_points.add(start_point);
				forbidden_layers.add(0);
				hpw.insert(start_point, 0.0);  // start_point -> weight == 0.0

				for(int i = 0; i < forbidden_points.countInt(); i++ ){
					MeshPoint3d* fpoint = forbidden_points[i];
					int layer = forbidden_layers[i];
					if(layer >= MAX_DIST) break;
					int ect = fpoint->getRank();
					for(int j = 0; j < ect; j++){
						MeshEdge3d* edge = fpoint->getEdge(j);
						if( !edge->isBorder() && !hsb_edges.contains(edge) ) continue; // edge border or marked as surf-border
						MeshPoint3d* other_point = edge->getOtherPoint(fpoint);
						if( ! hpoints.contains(other_point->getIndex()) ) continue;	// not in this sub-domain
						if( hpw.contains(other_point) ) continue;  // already checked
						// else -> should be forbidden
						forbidden_points.add(other_point);
						forbidden_layers.add(i+1);
						hpw.insert(other_point, -1.0); // forbidden_point
					}
				}

				DataVector<MeshPoint3d*> flow_processed(spct);
				DataVector<MeshPoint3d*> flow_active(spct);
				DataHashTableKeyValue<MeshPoint3d*, MeshPoint3d*> hprev(2*spct, nullptr);
				flow_active.add( start_point );

				// calculate flow (min-heap for calculating flow...)
				MeshPoint3d* flow_end_point = nullptr;
				while( flow_active.notEmpty() ){
					// find min ...
					int min_i = 0;
					double min_w = hpw.getValue( flow_active[0], -2.0);
					assert(min_w >= 0.0);
					for(int i = 1; i < flow_active.countInt(); i++){
						double w = hpw.getValue( flow_active[i], -2.0);
						assert(w > 0.0);
						if( w < min_w){ min_i = i; min_w = w; }
					}
					MeshPoint3d* act_point = flow_active.removeAt(min_i);
					flow_processed.add(act_point);

					//if(true){
					//	MeshViewSet* set = mesh->getDebugViewSetTopological(act_point, 4);
					//	set->addLabel( start_point->getCoordinates(), "SP");
					//	set->addInfo("min_w", to_string(min_w));
					//	SHOW_MESH("min-flow, act_point", set); 
					//}

					if(act_point != start_point && act_point->isBorder()){
						// found ending point!
						flow_end_point = act_point;
						break;
					}else{
						int ect = act_point->getRank();
						for(int i = 0; i < ect; i++){
							MeshEdge3d* edge = act_point->getEdge(i);
							if(edge->isBorder() || hsb_edges.contains(edge) ) continue; 
							MeshPoint3d* other_point = edge->getOtherPoint( act_point );
							double other_w = hpw.getValue( other_point, -2.0 );
							if(other_w == -1.0) continue; // forbidden							
							if(other_w != -2.0 && other_w < min_w){
								assert( flow_processed.contains( other_point ) );
								continue; // lower flow-weight, must have been checked already
							}
							// calculate weight for edge
							double sp = edge->getDoubleTag(TagExtended::TAG_NORMAL_SP, 0.5);
							double w_sp = std::max( sp, 0.0); // 0.0 .. 1.0
							double w_e = 0.2 + 0.8 * w_sp;

							double new_w = min_w + w_e;
							if( other_w == -2.0 ){ // first time reached
								hpw.insert( other_point, new_w );
								hprev.insert( other_point, act_point );
								flow_active.add( other_point );
							}else if( new_w < other_w ) { // already reached, but this route is better
								hpw.setValue( other_point, new_w );
								hprev.setValue( other_point, act_point );
							}
						}
					}
				}

				if( flow_end_point ){
					// is flow_end_point in open_border_points?
					if( open_border_points.contains( flow_end_point ) ){
						open_border_points.remove( flow_end_point );
					}
					// mark edges and points, backtracking...
					while( flow_end_point != start_point ){
						MeshPoint3d* prev_point = hprev.getValue( flow_end_point, nullptr );
						assert( prev_point != nullptr );
						MeshEdge3d* edge = flow_end_point->getEdgeToPoint( prev_point );
						assert( edge != nullptr );
						hsb_edges.insert( edge );
						hsb_points.insert( flow_end_point );
						flow_end_point = prev_point;
					}
				}else{
					LOG4CPLUS_ERROR(MeshLog::logger_console,   
						"MG3dS::fixOpenBorders - couldn't join open border");
					break;
				}
			}

			if(visualisation_on) {
				MeshViewSet* set = new MeshViewSet();
				for(int i = 0; i < sfct; i++) {
					MeshFace* face = sub_faces->get(i);
					int ect = face->getEdgeCount();
					for(int j = 0; j < ect; j++){
						MeshEdge3d* edge = face->getEdge(j);
						if(edge->isBorder()) 
							set->addEdge( edge, 1 );
						else if( hsb_edges.contains( edge ) )
							set->addEdge( edge, 2 );
					}
				}
				SHOW_MESH("borders / sub-borders for local sub-domain", set);
			}

			// reassign faces, create additional sub_domains
			DataHashTableKeyValue< MeshFace*, int > face_ids( 2*sfct, nullptr );
			for(int i = 0; i < sfct; i++)
				face_ids.insert( sub_faces->get(i), 1 );

			while( true ) {
				// ... start with first face
				assert( sub_faces->notEmpty() );
				MeshFace* face = sub_faces->get(0);
				assert( face_ids.getValue( face, 0 ) == 1 );
				DataVector<MeshFace*> afaces( sub_faces->countInt() );
				afaces.add( face );
				face_ids.setValue( face, 2 );
				int k = 0;
				// ... mark all adjacent faces, by adjacency, not crossing borders and sborders
				while( k < afaces.countInt() ) {
					face = afaces[k++];
					int ect = face->getEdgeCount();
					for(int i = 0; i < ect; i++){
						MeshEdge3d* edge = face->getEdge(i);
						if(edge->isBorder() || hsb_edges.contains(edge) ) continue;
						MeshFace* other_face = edge->getOtherFace( face );
						if( face_ids.getValue( other_face, 0 ) == 1 ) {
							afaces.add( other_face );
							face_ids.setValue( other_face, 2 );
						}
					}
				}
				// ... if number == total, no sub_sub_domain, no splitting -> end
				if( k == sub_faces->countInt() ) break;
				// ... else create new sub_domain, remove from this, repeat
				auto sub_sub_faces = std::make_shared<DataVector<MeshFace*>> ( k );
				int i = 0;
				while( i < sub_faces->countInt() ){
					face = sub_faces->get(i);
					if( face_ids.getValue( face, 0 ) == 2 ){ // move
						sub_faces->removeAt(i);
						sub_sub_faces->add( face );
						face->setIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_init );
					}else i++;
				}
				assert( sub_sub_faces->notEmpty() );
				sub_domains.add( sub_sub_faces );
				sub_tag_init++;
			}
		}
	}

	return sub_domains.countInt();
}

/// Identify local surfaces using greedy approach
int MeshGenerator3dSurface::identifyLocalSurfacesGreedy(Metric3dContext& mc, 
	MeshContainer3dSurface* mesh, 
	std::shared_ptr<const DataVector<MeshFace*>> sub_faces )
{
	int sub_tag_value = sub_faces->get(0)->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );
	assert( sub_tag_value >= 0 );

	int sfct_left = sub_faces->countInt(); 

	int ls_ct = 0;

	DataVector<MeshFace*> start_faces(1, nullptr);
	DataVector<MeshFace*> start_faces_empty(1);

	for(int it = 0; (sfct_left > 0) && (it < LOCAL_SURFACE_TOL_F_CT); it++){

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
			"===> Identifying local surfaces for TOL_F= " << LOCAL_SURFACE_TOL_F[it]
				<< " LAY_F= " << LOCAL_SURFACE_LAY_F[it] << " ====");

		// First, try match a surface for ALL points/faces (for sub-domain)
		int sfct = 0;
		mesh->approximateLocalSurface(mc, start_faces_empty, sfct, 
			LOCAL_SURFACE_TOL_F[it] * param_local_shape_tolerance, 
			LOCAL_SURFACE_LAY_F[it], LOCAL_SURFACE_MTR_F[it], 
			TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value);
		if(sfct > 0){
			++ls_ct;
			sfct_left -= sfct;
		}

		if(sfct_left == 0) break;

		// Browse faces - if any left
		for(int i = 0; i < sfct; i++){
			start_faces[0] = sub_faces->get(i);
			if(start_faces[0]->hasLocalSurface()) continue;
			assert(start_faces[0]->checkIntTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value));
			mesh->approximateLocalSurface(mc, start_faces, sfct,  
				LOCAL_SURFACE_TOL_F[it] * param_local_shape_tolerance, LOCAL_SURFACE_LAY_F[it], 
				LOCAL_SURFACE_MTR_F[it], TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value);
			if(sfct > 0) {
				++ls_ct;
				sfct_left -= sfct;
				if(sfct_left == 0) break;
			}
		}
	}


	if(sfct_left > 0) {
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"Failed to approximate surface for some faces, fct_left= " << sfct_left);
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
		"Created local surfaces set for subdomain: " << ls_ct << " surfaces.");

	return ls_ct;
}

/// Identify local surfaces using higher-curvature first
int MeshGenerator3dSurface::identifyLocalSurfacesHighCurvatureFirst( Metric3dContext& mc, 
	MeshContainer3dSurface* mesh, 
	std::shared_ptr<const DataVector<MeshFace*>> sub_faces )
{
	int sub_tag_value = sub_faces->get(0)->getIntTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1);
	assert( sub_tag_value >= 0 );

	int sfct = 0;

	DataVector<MeshFace*> start_faces_empty(1);

	for(int it = 0; it < LOCAL_SURFACE_TOL_F_CT_WHOLE; it++){

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
			"===> Identifying local whole-surface for TOL_F= " << LOCAL_SURFACE_TOL_F[it]
				<< " LAY_F= " << LOCAL_SURFACE_LAY_F[it] << " ====");

		// with empty "start_faces" - try fit for all faces in this sub-domain
		SurfaceConstPtr surface = mesh->approximateLocalSurface(mc, 
				start_faces_empty, sfct, 
				LOCAL_SURFACE_TOL_F[it] * param_local_shape_tolerance, 
				LOCAL_SURFACE_LAY_F[it], LOCAL_SURFACE_MTR_F[it], 
				TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value);

		if( surface != nullptr )
		{
			assert( sfct == sub_faces->countInt() );
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Created local surface for whole subdomain.");
			return 1;
		}
	}

	DataVector<MeshFace*> start_faces(1, nullptr);


	vector< MeshFace* > face_table;
	int sfct_left = sub_faces->countInt();
	face_table.reserve( sfct_left );

	for(int i = 0; i < sfct_left; i++) {
		face_table.push_back( sub_faces->get(i) );
	}

	// ... sort
	sort( face_table.begin(), face_table.end(), 
		[]( MeshFace* a, MeshFace* b) {	return a->getMaxCurvature() > b->getMaxCurvature();	} );

	int ls_ct = 0;
	//int pct = mesh->getPointsCount();

	for(int it = 0; (sfct_left > 0) && (it < LOCAL_SURFACE_TOL_F_CT); it++){

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
			"===> Identifying local surfaces for TOL_F= " << LOCAL_SURFACE_TOL_F[it]
				<< " LAY_F= " << LOCAL_SURFACE_LAY_F[it] << " ====");

		// Browse faces - if any left
		for(auto face : face_table) {
			start_faces[0] = face;
			if(face->hasLocalSurface()) continue;
			assert(face->checkIntTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value));
			mesh->approximateLocalSurface(mc, start_faces, sfct, 
				LOCAL_SURFACE_TOL_F[it] * param_local_shape_tolerance, LOCAL_SURFACE_LAY_F[it], 
				LOCAL_SURFACE_MTR_F[it], TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value);
			if(sfct > 0) {
				++ls_ct;
				sfct_left -= sfct;
				if(sfct_left == 0) break;
			}
		}
	}

	if(sfct_left > 0) {
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"Failed to approximate surface for some faces, fct_left= " << sfct_left);
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Created local surfaces set for subdomain: " << ls_ct << " surfaces.");

	return ls_ct;
}

/// Identify local surfaces using normal-clustering, then quadtree of quadrics
int MeshGenerator3dSurface::identifyLocalSurfacesViaNormals( Metric3dContext& mc, 
	MeshContainer3dSurface* mesh, 
	std::shared_ptr<DataVector<MeshFace*>> sub_faces )
{
	assert( sub_faces->notEmpty() );
	int sub_tag_value = sub_faces->get(0)->getIntTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1);
	assert( sub_tag_value >= 0 );

	int sfct = 0;
	int sfct_left = sub_faces->countInt();

	// ... sort
	sort( sub_faces->begin(), sub_faces->end(), 
		[]( MeshFace* a, MeshFace* b) {	return a->getMaxCurvature() > b->getMaxCurvature();	} );

	int ls_ct = 0;
	static const double SP_MIN = 0.5; // alfa: 0-60o

	sub_faces->forEach( 
		[&] ( MeshFace* face ) { 
			if(!face->hasLocalSurface()) {
				DataVector<MeshFace*> start_faces(1, face);
				assert(face->checkIntTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value));
				//mesh->gatherFirstLayerTopological( start_faces, 
				//	layers_count, p_list, p_layers,	f_list, f_layers,
				//	TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value);
				//mesh->gatherLayeredVerticesViaNormal( min_sp,
				//	layers_count, p_list, p_layers, f_list, f_layers, 
				//	TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value );
				mesh->approximateLocalSurfaceViaNormals(mc, start_faces, sfct, 
					SP_MIN, param_local_shape_tolerance,
					TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_value);
				if(sfct > 0) {
					++ls_ct;
					sfct_left -= sfct;
				}
			}
		} 
	);

	if(sfct_left > 0) {
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"Failed to approximate surface for some faces, fct_left= " << sfct_left);
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Created local surfaces set for subdomain: " << ls_ct << " surfaces.");

	return ls_ct;
}

/// Split surface mesh into closed sub-domains (using boundary and tag info
int MeshGenerator3dSurface::splitToSubdomains(MeshContainer3dSurface* mesh, 
	DataVector< std::shared_ptr<DataVector<MeshFace*>> > & sub_domains, 
	int sub_tag_init, TagExtended::TagType tag_type, int tag_value)
{
	int fct = mesh->getFacesCount();
	DataVector<MeshFace*> f_list(fct);

	for(int i = 0; i < fct; i++){
		MeshFace* face = mesh->getFaceAt(i);
		if(tag_type != TagExtended::TAG_NONE && face->hasIntFlag(tag_type, tag_value)) continue; // not applicable
		if(face->availableTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN)) continue; // alredy in some subdomain

		// gather full neighbourhood (not crossing borders)
		f_list.clear();
		f_list.add(face);
		face->setIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_init );
		int k = -1;
		while(++k < f_list.countInt()){
			face = f_list[k];
			int ect = face->getEdgeCount();
			for(int j = 0; j < ect; j++){
				MeshEdge3d* edge = face->getEdge(j);
				if(edge->isBorder()) continue;
				MeshFace* other_face = edge->getOtherFace(face);
				if(tag_type != TagExtended::TAG_NONE && other_face->checkIntTag(tag_type, tag_value)) continue; // not applicable
				if(other_face->availableTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN)){
					assert(other_face->checkIntTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_init));
					continue; // alredy included
				}
				f_list.add(other_face);
				other_face->setIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, sub_tag_init );
			}
		}

		//mesh->gatherFirstLayerTopological(face, f_layers_count, p_list, p_layers, f_list, f_layers, tag_type, tag_value);
		//bool other_layers_available = mesh->gatherLayeredVerticesTopological(f_layers_count, 
		//	p_list, p_layers, f_list, f_layers, -1, tag_type, tag_value);

		// create and mark subdomain
		auto sub_faces = std::make_shared<DataVector<MeshFace*>> ( f_list.countInt() );
		for(int j = 0; j < f_list.countInt(); j++){
			sub_faces->add( f_list[j] );
		}

		sub_domains.add( sub_faces );
		++sub_tag_init;
	}

	return sub_tag_init;
}

/*
/// Create additional local surfaces for points with low locals-surface-inside-ratio
int MeshGenerator3dSurface::complementLocalSurfaces( Metric3dContext& mc, MeshContainer3dSurface* mesh )
{
	// create heap of nodes with worst inside-ratio criterion
	// ... they will not be worse -> no heap really necessary ...

	// ... gather points with inside-ratio lower than threshold
	const double QUALITY_THRESHOLD = 0.3;
	int pct = mesh->getPointsCount();
	vector< pair<double, MeshPoint3d*> > qtable;

	for(int i = 0; i < pct; i++) {
		MeshPoint3d* point = mesh->getPointAt(i);
		if( point->isBorder() ) continue;

		const double & quality = point->getLocalSurfaceQuality();
		if( quality < QUALITY_THRESHOLD ){
			qtable.push_back( pair<double, MeshPoint3d*>( quality, point ) );
		}
	}
	// ... sort
	sort( qtable.begin(), qtable.end(), 
		[]( const pair<double, MeshPoint3d*>& a, const pair<double, MeshPoint3d*>& b) {
			return a.first < b.first;
		} );

	int surf_count = 0;
	DataVector<MeshFace*> start_faces(1);
	start_faces.add(nullptr);

	for(auto a : qtable) {
		SurfaceParametric* surf = a.second->getLocalSurface();
		const double& quality = a.second->getLocalSurfaceQuality();
		if( !surf || quality < QUALITY_THRESHOLD ) {
			start_faces[0] = a.second->getEdge(0)->getFaceAt(0);
			int sub_tag_value = start_faces[0]->getIntTag( sub_tag_type );
			for(int it = 0; it < LOCAL_SURFACE_TOL_F_CT; it++){
				int ct = mesh->approximateLocalSurface(mc, start_faces, 
					LOCAL_SURFACE_TOL_F[it] * tolerance, LOCAL_SURFACE_LAY_F[it], 
					LOCAL_SURFACE_MTR_F[it], sub_tag_type, sub_tag_value);
				if(ct >= 0) {
					++surf_count;
					break;
				}
			}
		}
	}

	return surf_count;
}
*/

/// improve local surface approximation quality
bool MeshGenerator3dSurface::identifyLocalCurves( Metric3dContext& mc, 
		MeshContainer3dSurface* mesh, double tolerance )
{
	// bsplines3d, always full-chain
	for( auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge() ) {
		MeshEdge3d* edge = it.getEdge();
		if(!edge->isBorder() || edge->hasLocalCurve()) continue;

		mesh->approximateLocalCurve(mc, edge, nullptr, nullptr, tolerance );
		assert( edge->hasLocalCurve() );
		if(false){
			auto curve = edge->getLocalCurve();
			MeshViewSet* set = new MeshViewSet;
			DataVector< DPoint3d > polyline;
			curve->getPolyLine(polyline);
			for(int j = 1; j < polyline.countInt(); j++)
				set->addEdge( polyline[j-1], polyline[j], 1 );
			set->addInfo("curve-type", curve->getSimpleDescription() );
			SHOW_MESH("Curve contour (full-chain) direct-3d", set);
		}
	}

	return true;
}

/// improve local surface approximation quality
int MeshGenerator3dSurface::fixLocalSurfaceSeams(Metric3dContext& mc, 
		MeshContainer3dSurface* mesh)
{
	if(false){
		MeshViewSet* set = new MeshViewSet;
		int to_fix = 0;
		for( auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge() ){
			MeshEdge3d* edge = it.getEdge();
			if( !edge->isBorder() && !edge->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_INNER, 1 ) ) {
				assert( edge->getFaceCount() == 2 );
				MeshFace* f0 = edge->getFaceAt(0);
				MeshFace* f1 = edge->getFaceAt(1);
				int stv0 = f0->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1);
				int stv1 = f1->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1);
				if( stv0 != stv1 ){
					//set->addEdge( edge, 3 );
					continue;
				}
				++to_fix;
				set->addEdge( edge, 2 );
			}
			else set->addEdge( edge, -1 );
		}
		if( to_fix > 0 ){
			set->addInfo( "edges to_fix", to_fix );
			SHOW_MESH( "MG3dS::fixLSSeams - check edges", set );
		}else
			delete set;
	}

	// find and fix non-border edge with no overlayed local surface
	DataVector<MeshFace*> faces(2, nullptr);
	int extra_surfaces = 0;
	int extra_failed = 0;
	for( auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge() ){
		MeshEdge3d* edge = it.getEdge();
		if( !edge->isBorder() && !edge->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_INNER, 1 ) ) {
			assert( edge->getFaceCount() == 2 );
			faces[0] = edge->getFaceAt(0);
			faces[1] = edge->getFaceAt(1);
			int stv0 = faces[0]->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1);
			int stv1 = faces[1]->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1);
			if( stv0 == stv1 ){
				SurfaceConstPtr surface = nullptr;
				int sfct = 0;
				for(int iti = 0; ( surface == nullptr) && (iti < LOCAL_SURFACE_TOL_F_CT); iti++){
					surface = mesh->approximateLocalSurface( mc, faces, sfct, 
						LOCAL_SURFACE_TOL_F[iti] * param_local_shape_tolerance, 
						LOCAL_SURFACE_LAY_F[iti], LOCAL_SURFACE_MTR_F[iti], 
						TagExtended::TAG_LOCAL_SURFACE_DOMAIN, stv0);
				}
				if( !edge->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_INNER, 1 ) ){
					LOG4CPLUS_WARN(MeshLog::logger_console, 
						"MG3dS::fixLSSeams - failed approx extra local surface");
					++extra_failed;
				}
				if( surface != nullptr ) ++extra_surfaces;
			}
		}
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "MG3dS::fixLSSeams - extra surfaces: " << extra_surfaces);
	if( extra_failed > 0) 
		LOG4CPLUS_WARN(MeshLog::logger_console, "MG3dS::fixLSSeams - extra failed: " << extra_failed);

	return extra_surfaces;
}

/*
/// improve local surface approximation quality
int MeshGenerator3dSurface::optimizeLocalSurfaces(Metric3dContext& mc, MeshContainer3dSurface* mesh)
{
	int pct = mesh->getPointsCount();
	int change_counter = 0;

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = mesh->getPointAt(i);
		if( point->isBorder() ) continue;

		if( point->selectBestLocalSurface( mc ) ) 
			change_counter++;
	}

	string label = "Reassigned surfaces";
	LOG4CPLUS_INFO(MeshLog::logger_console, label, change_counter);
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "*** " << label << " => " << change_counter);
	return change_counter;
}
*/

/// Check local surface information for the given point
bool MeshGenerator3dSurface::checkPointForLocalSurfaces( Metric3dContext& mc, 
		MeshContainer3dSurface* mesh, MeshPoint3d* point )
{
	SurfaceSetConstPtr cs_set;
	CS3dPtr cs = mesh->getControlSpace();
	if(cs && cs->isAdaptive())
		cs_set = cs->getAsAdaptive()->getLocalSurfaceSetAtPoint( point->getCoordinates() );

	MeshFace* pt_face = point->getEdge(0)->getFaceAt(0);
	assert( pt_face != nullptr );
	assert( pt_face->availableTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN ) );
	int domain_tag = pt_face->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );

	for(int i = 0; i < cs_set->countInt(); i++){
		SurfaceConstPtr surf = cs_set->getSurface(i);
		if( ! surf->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, domain_tag ) ) continue;

		MeshViewSet* set = mesh->getDebugViewSetTopological( point, 2, TagExtended::TAG_LOCAL_SURFACE_DOMAIN );
		set->addInfo( "l-surf #", to_string(i+1) + "/" + to_string(cs_set->countInt()) );

		SHOW_MESH("checkPoint-surface", set);
	}

	return true;
}

/// get local surface valid for the given set of (local) mesh points
SurfaceConstPtr MeshGenerator3dSurface::getLocalSurfaceForPoints( Metric3dContext& mc, 
		MeshContainer3dSurface* mesh, const DataVector<MeshFace*> & mfaces, 
		const DataVector< MeshPoint3d* > & mpoints )
{
	assert( mfaces.notEmpty() );
	assert( mpoints.notEmpty() );

	int local_surface_tag = mfaces[0]->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 );
	int mfct = mfaces.countInt();
	for(int i = 1; i < mfct; i++){
		if( ! mfaces[i]->checkIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, local_surface_tag ) )
			return nullptr; // all faces should be in the same local-subdomain
	}

	int mpct = mpoints.countInt();

	static int run_counter = 0;
	run_counter++;

	static unsigned int counters[3] = { 0 };

	const unsigned int HCT = 100;
	DataHashTable< SurfaceConstPtr > hchecked_surfaces( HCT, nullptr );

	SurfaceConstPtr opt_surface = nullptr;
	double opt_quality = AQ_INVALID;

	LOG_GETLOCALSURFACE( "START " << mpct );

	// 1. All points have the same valid surface ?
	for(int i = 0; i < mpct; i++ ) {
		mpoints[i]->getOptLocalCommonSurface( local_surface_tag, 
			mpoints, hchecked_surfaces, opt_surface, opt_quality );
	}

	if( opt_quality < AQ_VALID_MAX ) {
		// 2. Check acs-surfaces?
		CS3dPtr cs = mesh->getControlSpace();
		if (cs) {
			auto acs = cs->getAsAdaptive();
			if (acs != nullptr) {
				for (int i = 0; (opt_quality < AQ_VALID_MAX) && (i < mpct); i++) {
					auto lset = acs->getLocalSurfaceSetAtPoint(mpoints[i]->getCoordinates());
					for (int j = 0; (opt_quality < AQ_VALID_MAX) && (j < lset->countInt()); j++) {
						SurfaceConstPtr surface = lset->getSurface(j);
						if (surface->checkIntTag(TagExtended::TAG_LOCAL_SURFACE_DOMAIN, local_surface_tag) &&
							hchecked_surfaces.insert(surface))
						{
							double worst_q = 1.0;
							for (int k = 0; k < mpct; k++) {
								double q = mpoints[k]->getLocalSurfaceQuality(surface);
								if (q < worst_q) worst_q = q;
							}
							LOG_GETLOCALSURFACE("ACS-SURFACE-QUALITY " << worst_q);
							if (worst_q > opt_quality) {
								opt_surface = surface;
								opt_quality = worst_q;
							}
						}
					}
				}
			}
		}
	}

	if(opt_surface == nullptr) {
		LOG_GETLOCALSURFACE( "CREATE-NEW " << mpct );
		// 3. Create new local_surface ...
		//bool ok3 = mesh->checkLocalSurfaceParams();
		int sfct = 0;
		for(int it = 0; ( opt_surface == nullptr ) && (it < LOCAL_SURFACE_TOL_F_CT); it++){
			opt_surface = mesh->approximateLocalSurface( mc, mfaces, sfct, 
				LOCAL_SURFACE_TOL_F[it] * param_local_shape_tolerance, 
				LOCAL_SURFACE_LAY_F[it], LOCAL_SURFACE_MTR_F[it], 
				TagExtended::TAG_LOCAL_SURFACE_DOMAIN, local_surface_tag);
			if(opt_surface != nullptr) counters[2]++;
		}

		//bool ok4 = mesh->checkLocalSurfaceParams();
	}

	LOG_GETLOCALSURFACE( "END " << mpct );

	if( opt_surface != nullptr ) {
		for(int i = 0; i < mfct; i++)
			mfaces[i]->setLocalSurface( opt_surface );
		//for(int i = 0; i < mpct; i++)
		//	mpoints[i]->selectBestLocalSurface( opt_surface );
	}

	return opt_surface;
}

/// get local curve valid for the given set of (local) mesh points
//const Curve3dConstPtr* MeshGenerator3dSurface::getLocalCurveForPoints( Metric3dContext& mc, 
//		MeshContainer3dSurface* mesh, const DataVector<MeshEdge3d*> & medges, 
//		const DataVector< MeshPoint3d* > & mpoints )
//{
//	assert( mpoints.notEmpty() );
//	assert( medges.notEmpty() );
//
//	int mpct = mpoints.countInt();
//
//	static unsigned int counters[5] = { 0 };
//
//	// 1. All points have the same first curve ?
//	const Curve3dConstPtr* curve = mpoints[0]->getLocalCurve();
//	if( curve != nullptr ){
//		for(int i = 1; (curve != nullptr) && (i < mpct); i++)
//			if( mpoints[i]->getLocalCurve() != curve )
//				curve = nullptr;
//	
//		if( curve != nullptr ) {
//			for(int i = 0; i < medges.countInt(); i++)
//				medges[i]->setLocalCurve( curve );
//			counters[0]++;
//			return curve;
//		}
//	}
//
//	// 2. All points have some common valid curve ? (no domain-check)
//	const unsigned int HCT = 100;
//	DataHashTableKeyValue< const Curve3dConstPtr*, int > hcurve_counter(HCT, nullptr);
//	for( int i = 0; i < mct; i++ )
//		mpoints[i]->incForValidCurves( hcurve_counter );
//	for( unsigned int i = 0; i < hcurve_counter.countInt(); i++ ) {
//		curve = hcurve_counter.slotKey(i);
//		if( curve != nullptr && hcurve_counter.slotValue(i) == mct) {
//			for(int i = 0; i < medges.countInt(); i++)
//				medges[i]->setLocalCurve( curve );
//			counters[1]++;
//			return curve;
//		}
//	}
//
//	// 3. Check all available curves (with domain-check)
//	DataHashTable< const Curve3dConstPtr* > hcurve_check( HCT, nullptr );
//	DataVector< const Curve3dConstPtr* > curves( HCT );
//	double t;
//	for( int i = 0; i < mct; i++ ){
//		mpoints[i]->getAllLocalCurves( curves );
//		for( int j = 0; j < curves.countInt(); j++ ){
//			curve = curves[j];
//			if( hcurve_check.insert( curve ) ){
//				for(int k = 0; (curve != nullptr) && (k < mct); k++){
//					if( i != k ) {
//						if( ! mpoints[k]->checkAndGetCurveParam( curve, t ) )
//							curve = nullptr;
//					}
//				}
//				if( curve != nullptr ){
//					for(int i = 0; i < medges.countInt(); i++)
//						medges[i]->setLocalCurve( curve );
//					counters[2]++;
//					return curve;
//				}
//			}
//		}
//	}
//
//
//	// 5. Create new local_curve ...
//	//curve = nullptr;
//	//int sfct = 0;
//	//for(int it = 0; ( surface == nullptr ) && (it < LOCAL_SURFACE_TOL_F_CT); it++){
//	//	surface = mesh->approximateLocalSurface( mc, mfaces, sfct, 
//	//		LOCAL_SURFACE_TOL_F[it] * param_local_shape_tolerance, 
//	//		LOCAL_SURFACE_LAY_F[it], LOCAL_SURFACE_MTR_F[it], 
//	//		TagExtended::TAG_LOCAL_SURFACE_DOMAIN, local_surface_tag);
//	//	if(surface != nullptr){
//	//	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "getLocalSurf-#5: " << local_surface_tag << "/" 
//	//	//		<<  surface->getIntTag( TagExtended::TAG_LOCAL_SURFACE_DOMAIN, -1 ) << endl;
//	//		counters[4]++;
//	//	}
//	//}
//
//	//if( mfaces.countInt() == 2 ){
//	//	MeshViewSet* set = new MeshViewSet;
//	//	for(int i = 0; i < mfaces.countInt(); i++)
//	//		set->addFaceWithEdgesAndPoints( mfaces[i] );
//	//	if( surface != nullptr ) {
//	//		DataMatrix<DPoint2d> params( mpoints.countInt() );
//	//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "=================");
//	//		for(int j = 0; j < mpoints.countInt(); j++) {
//	//			const DPoint2d& mparam = mpoints[j]->getLocalSurfaceParam( surface );
//	//			params.add( mparam );
//	//			DPoint2d test_param = surface->getParameters( mpoints[j]->getCoordinates() );
//	//			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, j << " -> dist2d= " << mparam.distance(test_param) << "   /    "
//	//				<< "dist3d= " << surface->getPoint( mparam ).distance( mpoints[j]->getCoordinates() ) << endl;
//	//		}
//	//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "=================");
//	//		set = surface->createViewSetForPoints( set, params );			
//	//	}
//
//	//	SHOW_MESH( (surface == nullptr) ? "MG3dS::getLocalSurfaceForPoints - Couldn't find a surface?"
//	//		: "MG3dS::getLocalSurfaceForPoints - new local surface", set );
//	//}
//
//	//assert( surface != nullptr );
//
//	return curve;
//}

void MeshGenerator3dSurface::remeshSurfaceMesh( MeshContainer3dSurface* surface_mesh, double min_dist )
{
	MeshGenerator3dSurface::checkAndFixInvertedFacesTopological(surface_mesh);

	if( min_dist < 0.0 ) {
		DBox bbox = surface_mesh->getBoundingBox();
		double diameter = bbox.getDiameter();
		min_dist = diameter * 1e-3;
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Setting boundary by topology ...");
	//surface_mesh->clearBoundaryFlags();
	//surface_mesh->setBoundaryFeatureEdges();

	if(true){
		SHOW_MESH("remesh: top-boundary", surface_mesh->getViewSet());
	}

	MeshGenerator3dSurface::calculateNormals(surface_mesh);

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Mesh validation & removing tiny edges ...");
	MeshGenerator3dSurface::validateSurfaceMesh(surface_mesh, min_dist);

	if( !surface_mesh->isValid() ) LOG4CPLUS_WARN(MeshLog::logger_console, "Invalid mesh after validation");

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Setting boundary by sharp edges ...");
	//surface_mesh->setBoundarySharpEdges(MeshGenerator3dSurface::param_sharp_edge_threshold);

#ifdef _DEBUG
	if(true){
		int bp_count = 0, bpc_count = 0;
		int be_count = 0;
		MeshViewSet* bset = new MeshViewSet();
		int pct = surface_mesh->getPointsCount();
		for(int i = 0; i < pct; i++){
			MeshPoint3d* point = surface_mesh->getPointAt(i);
			if(point->isBorder()) { bset->addPoint(point); bp_count++; }
			if(point->isBorder(TagBorder::CORNER)) bpc_count++;
			int rank = point->getRank();
			for(int j = 0; j < rank; j++){
				MeshEdge3d* edge = point->getEdge(j);
				if(edge->getPointIndex(point) != 0) continue;
				if(edge->isBorder()) { bset->addEdge(edge); be_count++; }
			}
		}
		SHOW_MESH("OFF: create + validate", surface_mesh->getViewSet());
		if(bp_count > 0) {
			bset->addInfo("boundary points", bp_count);
			bset->addInfo("boundary corners", bpc_count);
			bset->addInfo("boundary edges", be_count);
			SHOW_MESH("OFF: boundary", bset);
		}else
			delete bset;
	}
#endif

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Creating ACS from approximate curvature ...");
	CS3dPtr cs = MeshGenerator3dSurface::createACSFromApproxCurvatureAverage(surface_mesh);
	//CS3dPtr cs = MeshGenerator3dSurface::createACSFromApproxCurvatureDirect(surface_mesh);
	surface_mesh->setControlSpace(cs);
	Metric3dContext mc(cs);

	surface_mesh->clearLocalShapes();
	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Identifying local surfaces & curves ...");
	MeshGenerator3dSurface::identifyLocalSurfaces( mc, surface_mesh );
	//LOG4CPLUS_INFO(MeshLog::logger_console,"Identifying local curves ...");
	//surface_mesh->identifyLocalCurves(mc,MeshGenerator3dSurface::param_local_shape_tolerance);

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Converting polys to triangles ...");
	surface_mesh->convertPolysToTriangles();

	MeshGenerator3dSurface::smoothen(mc, surface_mesh, 2, TagExtended::TAG_NONE, 0, 
		MeshData::SM_LAPLACE | MeshData::SM_TOP_SWAP | MeshData::SM_DEL_SWAP_COMPLETE); 

	//surface_mesh->checkSurfaceSetCounters( mc, false );

	if( !surface_mesh->isValid() ) LOG4CPLUS_WARN(MeshLog::logger_console, "Invalid mesh after poly-to-triangle");

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Surface remeshing with local transformations ...");
	MeshGenerator3dSurface::remeshSurfaceMeshWithLocalTransformations(mc, surface_mesh);

	//surface_mesh->checkSurfaceSetCounters( mc, true );
	//int removed_points = MeshGenerator3dSurface::validateSurfaceMesh(surface_mesh, min_dist);

	SHOW_MESH("OFF: create + remesh", surface_mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_DOMAIN));
	// check & process

	if( !surface_mesh->isValid() ) LOG4CPLUS_WARN(MeshLog::logger_console, "Invalid mesh after remeshing");

	LOG4CPLUS_DEBUG(MeshLog::logger_console,"Storing final mesh ...");
	// store data - OFF
	surface_mesh->storeOFF( "surface-mesh.remeshed.off" );
}
