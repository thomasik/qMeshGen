#include "MarchingCubes.h"
#include "common.h"
#include "MeshContainer3dSurface.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshTriangle3d.h"
#include "MeshViewSet.h"
#include "ControlSpace2dIdentity.h"
#include "MeshDomainVolume.h"
#include "MeshGenerator3dSurface.h"
#include "Metric3dContext.h"

MeshContainer3dSurface* MarchingCubes::createSurfaceMesh(
	MarchingModelDiscreteLayered& model)
{
	int nvx, nvy;
	model.getCellCount(nvx, nvy); // # cell-count -> number of vertices for marching-cells
	assert(nvx > 0);
	assert(nvy > 0);
	int nex = nvx - 1;
	int ney = nvy - 1;
	int nvxy = nvx * nvy;

	DRect brect = model.getLayerBoundingRect();
	double dx = brect.getDX() / nvx;
	double dy = brect.getDY() / nvy;
	assert(dx > 0.0 && dy > 0.0);

	double z0 = model.getZ0();
	double last_dz = 0.0;
	double dz = 0.0;

	DataMatrix< VertexData > gv0(nvy, nvx, VertexData()), gv1(nvy, nvx, VertexData());
	auto p_gv0 = &gv0;
	auto p_gv1 = &gv1;
	DataMatrix< EdgeData > gex0(nvy, nex, EdgeData()), gex1(nvy, nex, EdgeData());
	auto p_gex0 = &gex0;
	auto p_gex1 = &gex1;
	DataMatrix< EdgeData > gey0(ney, nvx, EdgeData()), gey1(ney, nvx, EdgeData());
	auto p_gey0 = &gey0;
	auto p_gey1 = &gey1;
	DataMatrix< EdgeData > gez(nvy, nvx, EdgeData());
	auto p_gez = &gez;

	// create coordinates for grid vertices
	DVector3d gdv(0.0, 0.0, dz / 2);
	for (int iy = 0; iy < nvy; iy++) {
		DPoint3d gpt(brect.x0+dx/2, brect.y0+dy/2+dy*iy, z0);
		for (int ix = 0; ix < nvx; ix++, gpt.x += dx) {
			gv0(iy, ix).coord = gpt;
			gv1(iy, ix).coord = gpt;
		}
	}

	MeshContainer3dSurface* mesh = new MeshContainer3dSurface(nvxy);

	int lct = 0;
	while (model.prepareNextLayer(dz)) {
		// set vertices
		double dz2 = last_dz + dz;
		for (int iy = 0; iy < nvy; iy++) {
			for (int ix = 0; ix < nvx; ix++) {
				auto& vd = p_gv1->get(iy, ix);
				vd.isoValue = model.isoValue(iy, ix, &vd.surfId);
				vd.coord.z += dz2;
				vd.mpoint = nullptr;
			}
		}
		// clear mpoints for edges
		p_gex1->forEach([](EdgeData& ed) { ed.mpoint = nullptr; });
		p_gey1->forEach([](EdgeData& ed) { ed.mpoint = nullptr; });
		p_gez->forEach([](EdgeData& ed) { ed.mpoint = nullptr; });

		// create triangles
		VertexData* cell_vertices[8];
		EdgeData* cell_edges[12];

		static const int CVI[8] = { 3, 2, 0, 1, 7, 6, 4, 5 };
		static const int CEI[12] = {/* ex */ 2, 0, 6, 4, 
									/* ey */ 3, 1, 7, 5,  
									/* ez */ 11, 10, 8, 9 };

		for (int iy = 0; iy < ney; iy++) {
			for (int ix = 0; ix < nex; ix++) {
				// prepare cell vertices/edges
				for (int j = 0; j < 4; j++) {
					cell_vertices[CVI[j]]   = &p_gv0->get(iy + (j/2), ix + (j%2));
					cell_vertices[CVI[j+4]] = &p_gv1->get(iy + (j/2), ix + (j%2));
					cell_edges[CEI[j + 8]]  = &gez.get(iy + (j / 2), ix + (j % 2));
				}
				for (int j = 0; j < 2; j++) {
					cell_edges[CEI[j]]     = &p_gex0->get(iy + j, ix);
					cell_edges[CEI[j + 2]] = &p_gex1->get(iy + j, ix);
					cell_edges[CEI[j + 4]] = &p_gey0->get(iy, ix + j);
					cell_edges[CEI[j + 6]] = &p_gey1->get(iy, ix + j);
				}
				// make triangles
				polygonise(cell_vertices, cell_edges, mesh, 128);
			}
		}
		// switch
		std::swap(p_gv0, p_gv1);
		std::swap(p_gex0, p_gex1);
		std::swap(p_gey0, p_gey1);
		last_dz = dz;

		++lct;
		//mesh->storeOFF("dicom-mesh-" + to_string(lct) + ".off");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Layer " << lct << " done.");
	}

	mesh->setControlSpace(std::make_shared<ControlSpace3dIdentity>());
	auto mdv = std::make_shared<MeshDomainVolume>(mesh);
	mdv->setAreaID(0);
	mesh->addDomainVolume(mdv);

	//SHOW_MESH(" marching cubes - surface mesh (sub_id)", mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_ID) );

	//markBorders(mesh, model);

	//int pct = mesh->getPointsCount();
	//for (int i = 0; i < pct; i++) {
	//	MeshPoint3d* mp = mesh->getPointAt(i);
	//	if (mp->isBorder()) continue;
	//	int sid = mp->getIntTag(TagExtended::TAG_LOCAL_SURFACE_ID, -1);
	//	model.moveToSurface(mp, sid);
	//}
	//
	//markFaces(mesh);

	return mesh;
}

MeshContainer3dSurface* MarchingCubes::createSurfaceMesh( 
	const MarchingModel& model, double isoLevel)
{
	DBox box = model.getBoundingBox();
	assert( box.valid );
	double dx = box.getDX();
	double dy = box.getDY();
	double dz = box.getDZ();
	assert( dx > 0.0 && dy > 0.0 && dz > 0.0 );

	double res = model.getResolution();
	assert( res > 0.0 );

	int nex = 1 + (int)(dx / res); // number of edges in dimension x
	int ney = 1 + (int)(dy / res); // number of edges in dimension y
	int nez = 1 + (int)(dz / res); // number of edges in dimension z
	int nvx = nex+1; // number of vertices in dimension x
	int nvy = ney+1; // number of vertices in dimension y
	int nvz = nez+1; // number of vertices in dimension z

	dx /= nex;
	dy /= ney;
	dz /= nez;

	MeshContainer3dSurface* mesh = new MeshContainer3dSurface( nex * ney );

	// prepare grid-vertices (and edges)
	int vct = nvx * nvy * nvz;
	DataVector< VertexData > grid_vertices( vct );
	VertexData vd;
	for(int iz = 0; iz < nvz; iz++) {
		double tz = (double) iz / nez; // [0.0, 1.0]
		vd.coord.z = box.z0 * (1.0 - tz) + box.z1 * tz;
		for(int iy = 0; iy < nvy; iy++) {
			double ty = (double) iy / ney; // [0.0, 1.0]
			vd.coord.y = box.y0 * (1.0 - ty) + box.y1 * ty;
			for(int ix = 0; ix < nvx; ix++) {
				double tx = (double) ix / nex; // [0.0, 1.0]
				vd.coord.x = box.x0 * (1.0 - tx) + box.x1 * tx;
				grid_vertices.add( vd );
			}
		}
	}

	int ectx = nex*nvy*nvz;
	int ecty = ney*nvx*nvz;
	int ectz = nez*nvx*nvy;
	int ect = ectx + ecty + ectz;
	DataVector< EdgeData > grid_edges( ect, EdgeData() );

	// calculate isovalues
	assert( vct == grid_vertices.countInt() );
	for(int i = 0; i < vct; i++) {
		VertexData & vd = grid_vertices[i];
		vd.isoValue = model.isoValue( vd.coord, &vd.surfId );
	}

	// create triangles
	int cell_ct = nex*ney*nez;
	VertexData* cell_vertices[8];
	EdgeData* cell_edges[12];
	int ncxy = nex * ney;
	int nvxy = nvx*nvy;
	int nexy = nvx*ney + nex*nvy;
	static const int CVI[8] = { 3, 2, 0, 1, 7, 6, 4, 5 }; 
	static const int CVD[8] = { 0, 1, nvx, nvx+1, nvxy, nvxy+1, nvxy+nvx, nvxy+nvx+1 }; 
	static const int CEI[12] = { /* ex */ 2, 0, 6, 4, /* ey */ 3, 1, 7, 5,  /* ez */ 11, 10, 8, 9 };
	static const int CED[12] = { 0, nex, nex*nvy, nex*(nvy+1), 	 0, 1, ney*nvx, ney*nvx+1,	 0, 1, nvx, nvx+1 };
	for(int i = 0; i < cell_ct; i++) {
		int iz = i / ncxy;
		int ixy = i % ncxy;
		int iy = ixy / nex;
		int ix = ixy % nex;
		// cell_vertices <- ...
		int icv = ix + iy*nvx + iz*nvxy;
		for(int j = 0; j < 8; j++)
			cell_vertices[ CVI[j] ] = &grid_vertices[icv+CVD[j]];
		// cell_edges <- ....
		int icex = ix + iy * nex + iz * nex*nvy;
		int icey = ectx + ix + iy * nvx + iz * ney*nvx;
		int icez = ectx + ecty + ix + iy * nvx + iz * nvx*nvy;
		for(int j = 0; j < 4; j++) {
			cell_edges[ CEI[j] ]   = &grid_edges[icex+CED[j]];
			cell_edges[ CEI[j+4] ] = &grid_edges[icey+CED[j+4]];
			cell_edges[ CEI[j+8] ] = &grid_edges[icez+CED[j+8]];
		}

		if(false){
			MeshViewSet* set = new MeshViewSet;
			for(int j = 0; j < 8; j++) {
				set->addLabel( cell_vertices[ CVI[j] ]->coord, to_string( icv+CVD[j] ) );
			}
			set->addInfo("nex ney nez", to_string(nex)+"x"+to_string(ney)+"x"+to_string(nez) );
			set->addInfo("cell nr", i);
			int enr[12];
			for(int j = 0; j < 4; j ++){
				enr[ CEI[j] ]   = icex+CED[j];
				enr[ CEI[j+4] ] = icey+CED[j+4];
				enr[ CEI[j+8] ] = icez+CED[j+8];
			}
			string estr;
			for(int j = 0; j < 12; j++) estr += to_string(enr[j])+" ";
			set->addInfo("edges", estr);
			SHOW_MESH( "marching cube cell", set );
		}

		// make triangles
		polygonise( cell_vertices, cell_edges, mesh, isoLevel );
	}

	mesh->setControlSpace(std::make_shared<ControlSpace3dIdentity>());
	auto mdv = std::make_shared<MeshDomainVolume>(mesh);
	mdv->setAreaID(0);
	mesh->addDomainVolume(mdv);

	//SHOW_MESH(" marching cubes - surface mesh (sub_id)", mesh->getViewSet(nullptr, TagExtended::TAG_LOCAL_SURFACE_ID) );

	markBorders( mesh, model );

	int pct = mesh->getPointsCount();
	for(int i = 0; i < pct; i++) {
		MeshPoint3d* mp = mesh->getPointAt(i);
		if( mp->isBorder() ) continue;
		int sid = mp->getIntTag( TagExtended::TAG_LOCAL_SURFACE_ID, -1 );
		model.moveToSurface( mp, sid );
	}

	markFaces( mesh );

	MeshGenerator3dSurface::remeshSurfaceMesh( mesh, 0.1 * res );

	return mesh;
}

void MarchingCubes::markBorders( MeshContainer3dSurface* mesh, const MarchingModel& model )
{
	// check edges ...
	for( auto it = mesh->getFirstEdge3d(); it.isValid(); it.nextEdge() ) {
		MeshEdge3d* edge = it.getEdge();
		if( edge->isBorder() ) continue;
		MeshPoint3d* mp0 = edge->getMeshPoint(0);
		MeshPoint3d* mp1 = edge->getMeshPoint(1);
		if( mp0->isBorder() || mp1->isBorder() ) continue;
		int sid0 = mp0->getIntTag( TagExtended::TAG_LOCAL_SURFACE_ID, -1 );
		int sid1 = mp1->getIntTag( TagExtended::TAG_LOCAL_SURFACE_ID, -1 );
		if( sid0 != sid1 ) {
			// split
			MeshPoint3d* split_point = MeshGenerator3dSurface::splitEdgeSimple( mesh, edge );
			//SHOW_MESH("markBorders::after split", mesh->getDebugViewSetTopological( split_point, 3 ) );
			// move iteratively
			int counter = 0;
			DPoint3d old_coord = split_point->getCoordinates();
			double dist2 = 0.0;
			do {
				model.moveToSurface( split_point, sid0 );
				model.moveToSurface( split_point, sid1 );
				const DPoint3d& coord = split_point->getCoordinates();
				dist2 = old_coord.distance2( coord );
				old_coord = coord;
				counter++;
				//LOG4CPLUS_INFO(MeshLog::logger_mesh, "dist2 == " << dist2);
			}while( (dist2 > SMALL_NUMBER) && (counter < 100) );
			// mark
			split_point->setBorder( TagBorder::INNER | TagBorder::RIDGE );
			//SHOW_MESH("markBorders::after move", mesh->getDebugViewSetTopological( split_point, 3 ) );
			// mark adjacent edges - in created triangles...
			int rank = split_point->getRank();
			for(int i = 0; i < rank; i++){
				MeshEdge3d* medge = split_point->getEdge(i);
				MeshPoint3d* mpoint = medge->getOtherPoint( split_point );
				if( mpoint == mp0 || mpoint == mp1 ) continue;
				if( mpoint->isBorder() )
					medge->setBorder( TagBorder::INNER | TagBorder::RIDGE );
			}
		}
	}
}

void MarchingCubes::markFaces( MeshContainer3dSurface* mesh)
{
	int fct = mesh->getFacesCount();
	for(int i = 0; i < fct; i++) {
		MeshFace* mf = mesh->getFaceAt(i);
		int fpct = mf->getPointCount();
		int surf_id = -2;
		for(int j = 0; (j < fpct) && (surf_id != -1); j++){
			MeshPoint3d* mp = mf->getPoint(j);
			if(mp->isBorder()) continue;
			int sid = mp->getIntTag( TagExtended::TAG_LOCAL_SURFACE_ID, -1 );
			assert( sid > -1 );
			if(sid == -1) surf_id = -1; // shouldn't happen
			if( surf_id == -2 ) surf_id = sid;
			else if( surf_id != sid) surf_id = -1; // different
		}
		if( surf_id > -2 )
			mf->setIntTag( TagExtended::TAG_LOCAL_SURFACE_ID, surf_id );
	}
}

/*
   http://paulbourke.net/geometry/polygonise/
   Given a grid cell and an isolevel, calculate the triangular facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles" will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above of totally below the isolevel.
*/
int MarchingCubes::polygonise(
		MarchingCubes::VertexData* cell_vertices[8], 
		MarchingCubes::EdgeData* cell_edges[12],
		MeshContainer3dSurface* mesh, 
		const double & isoLevel )
{

	// Determine the index into the edge table which
	// tells us which vertices are inside of the surface
	int cubeIndex = 0;
	for(int i = 0; i < 8; i++)
		if (cell_vertices[i]->isoValue < isoLevel) cubeIndex |= TWO_POWER[i];

	/* Cube is entirely in/out of the surface */
	if (EDGE_TABLE[cubeIndex] == 0)
		return 0;

	if(false){
		MeshViewSet* set = new MeshViewSet;
		for(int i = 0; i < 8; i++) {
			set->addPoint( cell_vertices[i]->coord, (cell_vertices[i]->isoValue < isoLevel) ? 1 : 0, i);
		}
		SHOW_MESH( "marching cube partial cell", set );
	}

	/* Find the vertices where the surface intersects the cube */
	for(int i = 0; i < 12; i++)
		if (EDGE_TABLE[cubeIndex] & TWO_POWER[i]) {
			VertexData * vd[2] = { cell_vertices[ EDGE_VERTICES[i][0] ], cell_vertices[ EDGE_VERTICES[i][1] ] };
			if( cell_edges[i]->mpoint == nullptr && vd[0]->mpoint == nullptr && vd[1]->mpoint == nullptr) {
				DPoint3d ept;
				int surfId;
				int res = vertexInterpolate( *vd[0], *vd[1], ept, isoLevel, &surfId );
				MeshPoint3d* mp = new MeshPoint3d( ept );
				mp->setIntTag( TagExtended::TAG_LOCAL_SURFACE_ID, surfId );
				mesh->addMeshPoint( mp );
				if( res < 2 ) vd[res]->mpoint = mp;
				else cell_edges[i]->mpoint = mp;
			}
		}
	/* Create the triangles */
	int ntriang = 0;
	for (int i=0; TRI_TABLE[cubeIndex][i]!=-1; i+=3) {
		MeshPoint3d* mpts[3] = { getRefMeshPoint(cubeIndex, i, cell_vertices, cell_edges ),
			getRefMeshPoint(cubeIndex, i+1, cell_vertices, cell_edges ),
			getRefMeshPoint(cubeIndex, i+2, cell_vertices, cell_edges ) };

		if( mpts[0] != mpts[1] && mpts[0] != mpts[2] && mpts[1] != mpts[2] ) { // all different
			mesh->addMeshFace( new MeshTriangle3d( mpts[0], mpts[1], mpts[2] ) );
			ntriang++;
		}
	}

	return ntriang;
}

MeshPoint3d* MarchingCubes::getRefMeshPoint(int ci, int i, 
				MarchingCubes::VertexData* cell_vertices[8], 
				MarchingCubes::EdgeData* cell_edges[12] )
{
	int ei = TRI_TABLE[ci][i];
	if( cell_edges[ei]->mpoint != nullptr) return cell_edges[ei]->mpoint;
	else if( cell_vertices[ EDGE_VERTICES[ei][0] ]->mpoint != nullptr) return cell_vertices[ EDGE_VERTICES[ei][0] ]->mpoint;
	else {
		assert( cell_vertices[ EDGE_VERTICES[ei][1] ]->mpoint != nullptr );
		return cell_vertices[ EDGE_VERTICES[ei][1] ]->mpoint;
	}
}

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
int MarchingCubes::vertexInterpolate(const MarchingCubes::VertexData& vd0, 
						const MarchingCubes::VertexData& vd1, DPoint3d& pt, 
						const double& isoLevel, int * surfId )
{
	//if( std::abs(isoLevel - vd0.isoValue) < 1e-7 ) {
	//	pt = vd0.coord;
	//	return 0;
	//}
	//if( std::abs(isoLevel - vd1.isoValue) < 1e-7 ) {
	//	pt = vd1.coord;
	//	return 1;
	//}
	//if( std::abs(vd0.isoValue - vd1.isoValue) < 1e-7 ) {
	//	assert(false);
	//	pt = vd0.coord;
	//	return 0;
	//}

	if( std::abs(vd0.isoValue - vd1.isoValue) < VERY_SMALL_NUMBER ) {
		pt = DPoint3d( vd0.coord, vd1.coord, 0.5 );
		if( surfId != nullptr) *surfId = vd0.surfId;
		return 2;
	}

	double mu = (isoLevel - vd0.isoValue) / (vd1.isoValue - vd0.isoValue);

	static const double THRESHOLD = 0.01;
	if( mu < THRESHOLD ) {	
		pt = vd0.coord; 		
		if( surfId != nullptr) *surfId = vd0.surfId;
		return 0; 
	}else if ( mu > (1-THRESHOLD) ) { 
		pt = vd1.coord; 
		if( surfId != nullptr) *surfId = vd1.surfId;
		return 1; 
	}

	pt = DPoint3d( vd0.coord, vd1.coord, mu );
	if( surfId != nullptr) *surfId = (mu > 0.5) ? vd1.surfId : vd0.surfId;
	return 2;
}

int MarchingCubes::EDGE_TABLE[256]={
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   
};

int MarchingCubes::TRI_TABLE[256][16] = 
{
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
	{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
	{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
	{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
	{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
	{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
	{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
	{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
	{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
	{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
	{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
	{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
	{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
	{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
	{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
	{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
	{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
	{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
	{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
	{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
	{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
	{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
	{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
	{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
	{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
	{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
	{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
	{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
	{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
	{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
	{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
	{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
	{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
	{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
	{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
	{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
	{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
	{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
	{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
	{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
	{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
	{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
	{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
	{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
	{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
	{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
	{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
	{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
	{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
	{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
	{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
	{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
	{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
	{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
	{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
	{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
	{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
	{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
	{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
	{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
	{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
	{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
	{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
	{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
	{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
	{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
	{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
	{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
	{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
	{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
	{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
	{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
	{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
	{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
	{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
	{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
	{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
	{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
	{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
	{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
	{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
	{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
	{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
	{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
	{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
	{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
	{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
	{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
	{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
	{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
	{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
	{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
	{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
	{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
	{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
	{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
	{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};

int MarchingCubes::TWO_POWER[12] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 };
int MarchingCubes::EDGE_VERTICES[12][2] = 
			{ {0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7} };
