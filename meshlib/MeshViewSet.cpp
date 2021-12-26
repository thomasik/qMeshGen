/////////////////////////////////////////////////////////////////////////////
// MeshViewSet.cpp
// Class for storing data for visualisation
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "MeshViewSet.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshElement.h"
#include "MeshFace.h"
#include "MeshBlock.h"
#include "SurfaceParametric.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshTetrahedron.h"

#include "DLeastSquaresFitting.h"
#include "SurfacePlane.h"
#include "DTriangle.h"
#include "DPlane.h"

#define USE_BOARD

#ifdef USE_BOARD
#include "Board.h"
using namespace LibBoard;
#endif

/// face/block shrink for visualization
float MeshViewSet::param_shrink = 0.95f;
/// whether any visualization should be prepared/shown
int MeshViewSet::param_show_visualization = 1; // 1 -> yes, 0 -> no
/// reference to the visualization module
MeshViewExt* MeshViewSet::m_view = nullptr;

void MeshViewSet::addTetra(const FPoint3d& pt, float dx, float dy, float dz, double q)
{
	auto data = std::make_shared<MeshViewBlockData>(4, 4);
	data->area_id = -2;
	data->quality = q;

	data->pts.add(pt);
	data->pts.add(pt + FVector3d(0.0, dy, 0.0));
	data->pts.add(pt + FVector3d(dx, 0.0, 0.0));
	data->pts.add(pt + FVector3d(0.0, 0.0, dz));

	//------------
	for (int i = 0; i < 12; i++)
		data->indices.add(MeshTetrahedron::PCONF[i / 3][i % 3]);

	for (int i = 0; i < 4; i++)
		data->normals.add(FVector3d::crossProduct(
			data->pts[MeshTetrahedron::PCONF[i][0]],
			data->pts[MeshTetrahedron::PCONF[i][1]], 
			data->pts[MeshTetrahedron::PCONF[i][2]]).normalized());

	m_blocks.add(data);
}

void MeshViewSet::addBlock(const MeshBlock* block, int id)
{
	if(m_hash_blocks.contains(block)) return;
	auto data = block->getViewData(param_shrink);
	if(id > -2) data->area_id = id;
	if(block->getVolumeNoMetric() < 0.0) data->area_id = -10;
	m_blocks.add(data);
	m_hash_blocks.insert(block, data);
	// adjacency info
	int fct = block->getFaceCount();
	for(int i = 0; i < fct; i++){
		MeshBlock* other_block = block->getNeighbour(i);
		if (other_block == nullptr) continue;
		auto other_data = m_hash_blocks.getValue(other_block, nullptr);
		if(other_data){ // already inserted
			other_data->adjacent.add(data);
			data->adjacent.add(other_data);
		}
	}
}

void MeshViewSet::addBlockWithEdges(const MeshBlock* block, int id)
{
	addBlock(block, id);
	int ect = block->getEdgeCount();
	for(int i = 0; i < ect; i++){
		addEdge(block->getEdge(i));
	}
}

void MeshViewSet::addEmptyBlockWithEdges(const MeshBlock* block, int id)
{
	int area_id = block->getAreaID();
	if(id > -2) area_id = id;
	if(block->getVolumeNoMetric() < 0.0) area_id = -10;

	int ect = block->getEdgeCount();
	for(int i = 0; i < ect; i++){
		addEdge(block->getEdge(i), area_id);
	}
}

void MeshViewSet::addEdges(const MeshBlock* block)
{
	int ect = block->getEdgeCount();
	for(int i = 0; i < ect; i++){
		addEdge(block->getEdge(i));
	}
}

void MeshViewSet::addEdges(const MeshFace* face)
{
	int ect = face->getEdgeCount();
	for(int i = 0; i < ect; i++){
		addEdge(face->getEdge(i));
	}
}

void MeshViewSet::addElement(const MeshElement* element, SurfaceConstPtr surface, bool proper_orientation, int id)
{
	if(element->getEdgeCount() == 0) return;
	if(m_hash_faces.contains(element)) return;
	auto data = element->getViewData(param_shrink, surface, proper_orientation);
	if(id > -2) data->area_id = id;
	//if(element->isTagged()) data->part = 1;
	m_faces.add(data);
	m_hash_faces.insert(element, data);
}

void MeshViewSet::addElementWithEdges(const MeshElement* element, SurfaceConstPtr surface, bool proper_orientation, int id)
{
	addElement(element, surface, proper_orientation, id);
	int ect = element->getEdgeCount();
	for(int i = 0; i < ect; i++){
		addEdge(element->getEdge(i), surface);
	}
}

void MeshViewSet::addPoint(const MeshPoint3d* point, int part)
{
	if(m_hash_points.contains(point)) return;
	auto data = point->getViewData();
	data->part = part;
	m_points.add(data);
	m_hash_points.insert(point, data);
}

void MeshViewSet::addPoint(const FPoint3d& point, int part)
{
	auto data = std::make_shared<MeshViewPointData>(point);
	data->part = part;
	m_points.add(data);
}

void MeshViewSet::addPoint(const FPoint3d& point, int part, int id)
{
	auto data = std::make_shared<MeshViewPointData>(point, false, id, true);
	data->part = part;
	m_points.add(data);
}

void MeshViewSet::addPoint(const MeshPoint2d* point, SurfaceConstPtr surface, int part)
{
	if(m_hash_points.contains(point)) return;
	auto data = point->getViewData(surface);
	data->part = part;
	m_points.add(data);
	m_hash_points.insert(point, data);
}

void MeshViewSet::addEdge(const MeshEdge3d* edge, int id)
{
	if(m_hash_edges.contains(edge)) return;
	auto data = edge->getViewData();
	if(id > -2) data->part = id;
	m_edges.add(data);
	m_hash_edges.insert(edge, data);
	// adjacency info
	DataVector<MeshBlock*> blocks;
	if(edge->adjacentBlocks(blocks, false)){
		for(int i = 0; i < blocks.countInt(); i++){
			auto block_data = m_hash_blocks.getValue(blocks[i], nullptr);
			if(block_data){ // already inserted
				data->adjacent_blocks.add(block_data);
			}
		}
	}
}

void MeshViewSet::addEdge(const MeshEdge2d* edge, SurfaceConstPtr surface)
{
	if(m_hash_edges.contains(edge)) return;
	auto data = edge->getViewData(surface);
	m_edges.add(data);
	m_hash_edges.insert(edge, data);
}

void MeshViewSet::addEdge(const FPoint3d& p0, const FPoint3d& p1, int id)
{
	auto data = std::make_shared<MeshViewEdgeData>(p0, p1);
	if(id > -2) data->part = id;
	m_edges.add(data);
}

void MeshViewSet::addFace(const MeshFace* face, int id, double shrink, bool proper_orientation)
{
	if(m_hash_faces.contains(face)) return;
	auto data = face->getViewData(shrink, proper_orientation);
	if(data){
		if(id > -2) data->area_id = id;
		m_faces.add(data);
		m_hash_faces.insert(face, data);
	}
}

void MeshViewSet::addFaceWithEdges(const MeshFace* face, int id, double shrink, bool proper_orientation)
{
	addFace(face, id, shrink, proper_orientation);
	int ect = face->getEdgeCount();
	for(int i = 0; i < ect; i++){
		addEdge(face->getEdge(i));
	}
}

void MeshViewSet::addFaceWithEdgesAndPoints(const MeshFace* face, int id, double shrink, bool proper_orientation)
{
	addFace(face, id, shrink, proper_orientation);
	int ect = face->getEdgeCount();
	for(int i = 0; i < ect; i++){
		addEdge(face->getEdge(i));
	}
	int pct = face->getPointCount();
	for(int i = 0; i < pct; i++){
		addPoint(face->getPoint(i));
	}
}

void MeshViewSet::addFace(const FPoint3d& p0, const FPoint3d& p1, const FPoint3d& p2, int id)
{
	auto data = std::make_shared<MeshViewFaceData>(3);
	data->quality = 1.0;
	data->pts.add(p0);
	data->pts.add(p1);
	data->pts.add(p2);
	data->indices.add(0); data->indices.add(1); data->indices.add(2);
	data->normal = (p1-p0).crossProduct(p2-p0).normalized();
	data->area_id = (id > -2) ? id : 0;
	m_faces.add(data);
}

void MeshViewSet::addFace(const MeshFace* face, const FVector3d& shift, int id)
{
	if(m_hash_faces.contains(face)) return;
	auto data = face->getViewData();
	if(data){
		for(int i = 0; i < data->pts.countInt(); i++)
			data->pts[i] += shift;
		if(id > -2) data->area_id = id;
		m_faces.add(data);
		m_hash_faces.insert(face, data);
	}
}

void MeshViewSet::addPolygonConvex(const DataVector<DPoint3d> & poly, int id)
{
	int pct = poly.countInt();

	auto data = std::make_shared<MeshViewFaceData>(pct);
	for(int i = 0; i < pct; i++) data->pts.add(poly[i]);
	data->quality = 1.0;
	data->area_id = (id > -2) ? id : 0;
	for(int i = 2; i < pct; i++){
		data->indices.add( 0 );
		data->indices.add( i-1 );
		data->indices.add( i );
	}
	data->countNormal();

	m_faces.add(data);
}

void MeshViewSet::addPolygon(const DataVector<DPoint3d> & poly, int id)
{
	int pct = poly.countInt();
	if(pct == 3){
		// simple triangular face
		addFace(poly[0], poly[1], poly[2], id);
	}else if(pct > 3){
		// polygon, assuming planar - to triangulate (ear-clipPIng...)
		DPlane plane;
		double res = DLeastSquaresFitting::fitHyperplaneOrthogonal(poly, plane, true);

		auto data = std::make_shared<MeshViewFaceData>(pct);
		for(int i = 0; i < pct; i++) data->pts.add(poly[i]);
		data->quality = 1.0;
		data->normal = plane.vn;
		data->area_id = (id > -2) ? id : 0;

		DataVector<DPoint2d> points2d(pct);
		for(int i = 0; i < pct; i++) points2d.add( plane.projectToPlane(poly[i]) );
		// ... split into triangles
		DataVector<int> l(pct), r(pct);
		for(int i = 0; i < pct; i++) {
			l.add( (i+pct-1) % pct );
			r.add( (i+1) % pct);
		}
		int k = pct-1;
		while(pct > 3){
			k = r[k];
			// -- is triangle ok?
			DPoint2d tpts[3] = { points2d[l[k]], points2d[k], points2d[r[k]] };
			if( DTriangle2d::det(tpts[0], tpts[1], tpts[2]) <= 0.0) continue;
			// -- does it containe any other points
			bool valid = true;
			int jk =  r[r[k]]; // first to the right, not part of the tested triangle
			for(int j = 3; valid && (j < pct); j++, jk = r[jk]) { // skipPIng the three selected vertices
				valid = !DTriangle2d::containsPoint(tpts[0], tpts[1], tpts[2], points2d[jk]);
			}
			if(valid) {
				--pct;
				data->indices.add( l[k] );
				data->indices.add( k );
				data->indices.add( r[k] );
				l[ r[k] ] = l[k];
				r[ l[k] ] = r[k];
			}
		}
		k = r[k];
		// ... and the last three vertices
		data->indices.add( l[k] );
		data->indices.add( k );
		data->indices.add( r[k] );

		m_faces.add(data);
	}
	// else - nothing
}


DBox MeshViewSet::getBoundingBox() const
{
	DBox box;
	// points
	int count = m_points.countInt();
	for(int i = 0; i < count; i++)
		box.addPoint(m_points[i]->pt);
	// edges
	count = m_edges.countInt();
	for(int i = 0; i < count; i++){
		const auto& data = m_edges[i];
		box.addPoint(data->pt1);
		box.addPoint(data->pt2);
	}
	// faces
	count = m_faces.countInt();
	for(int i = 0; i < count; i++){
		const auto& data = m_faces[i];
		for(int j = 0; j < data->pts.countInt(); j++)
			box.addPoint(data->pts[j]);
	}
	// blocks
	count = m_blocks.countInt();
	for(int i = 0; i < count; i++){
		const auto& data = m_blocks[i];
		for(int j = 0; j < data->pts.countInt(); j++)
			box.addPoint(data->pts[j]);
	}

	return box;
}

void MeshViewSet::log()
{
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "MeshViewSet, bb=" << getBoundingBox());
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Blocks:");
	int count = m_blocks.countInt();
	for(int i = 0; i < count; i++){
		const auto& data = m_blocks[i];
		ostringstream log_line;
		log_line <<  data->pts.countInt() << ") -> ";
		for(int j = 0; j < data->pts.countInt(); j++)
			log_line << data->pts[j] << " ";
		LOG4CPLUS_INFO(MeshLog::logger_mesh, log_line.str());
	}

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Faces:");
	count = m_faces.countInt();
	for(int i = 0; i < count; i++){
		const auto& data = m_faces[i];
		ostringstream log_line;
		log_line << data->pts.countInt() << ") -> ";
		for(int j = 0; j < data->pts.countInt(); j++)
			log_line << data->pts[j] << " ";
		LOG4CPLUS_INFO(MeshLog::logger_mesh, log_line.str());
	}

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Edges:");
	count = m_edges.countInt();
	for(int i = 0; i < count; i++){
		const auto& data = m_edges[i];
		LOG4CPLUS_INFO(MeshLog::logger_mesh, data->pt1 << " " << data->pt2);
		}

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Vertices:");
	count = m_points.countInt();
	for(int i = 0; i < count; i++){
		const auto& data = m_points[i];
		LOG4CPLUS_INFO(MeshLog::logger_mesh, data->pt);
	}
}

bool MeshViewSet::storeMatlabFile(const string& fname, bool use_mesh, const DPoint3d& clip, double quality_clip) const
{
	if(use_mesh){
		if(m_mesh2d) return m_mesh2d->storeMatlabFile(fname, clip, quality_clip);
		else if(m_mesh3d) return m_mesh3d->storeMatlabFile(fname, clip, quality_clip);
		else LOG4CPLUS_WARN(MeshLog::logger_console, "No mesh available for storeMatlabFile");
	}

	ofstream os(fname.c_str());

	os << "clear all; figure(1); clf; hold on; axis equal; axis off; view(3);" << endl;
	os << "c(1,:) = [0.9,0.9,0.9];" << endl;
	os << "c(2,:) = [0.7,0.7,0.7];" << endl;
	os << "colormap(c);" << endl;
	os << "load 'mesh_vertices.txt' -ascii" << endl;
	os << "load 'mesh_faces.txt' -ascii" << endl;
	os << "load 'mesh_ci.txt' -ascii" << endl;
	os << "patch('Vertices',mesh_vertices,'faces',mesh_faces,'FaceVertexCData',mesh_ci(:,1),'FaceColor','flat');" << endl;

	ofstream ofsv("mesh_vertices.txt");
	ofstream ofsf("mesh_faces.txt");
	ofstream ofsc("mesh_ci.txt");

	// blocks
	int count = m_blocks.countInt();
	// check visibility
	int clip_count = 0;
	int hidden_count = 0;
	// 1. check clip
	for(int i = 0; i < count; i++){
		const auto& data = m_blocks[i];
		data->hidden = 0;
		if(quality_clip < 1.05 && data->quality > quality_clip)
			continue;
		for(int j = 0; j < data->pts.countInt(); j++)
			if(data->pts[j].x > clip.x || data->pts[j].y > clip.y || data->pts[j].z > clip.z){
				data->hidden = 1; 
				++clip_count;
				break;
			}
	}
	// 2. check sight-line
	for(int i = 0; i < count; i++){
		const auto& data = m_blocks[i];
		if(data->hidden) continue;
		data->hidden = 2;
		for(int j = 0; j < data->adjacent.countInt(); j++){
			// if boundary or near clip-plane
			if((data->adjacent[j] == nullptr) || (data->adjacent[j]->hidden == 1)){
				data->hidden = 0; // make visible
				break;
			}
		}
		if(data->hidden == 2) ++hidden_count;
	}

	if(true){
		ostringstream osstr;
		osstr << "Store blocks total/visible/clipped/hidden = " << count << '/' 
			<< (count - clip_count - hidden_count) << '/'
			<< clip_count << '/' << hidden_count;
		LOG4CPLUS_INFO(MeshLog::logger_console, osstr.str());
	}

	int vert_index = 1;
	for(int i = 0; i < count; i++){
		const auto& data = m_blocks[i];
		if(data->normals.countInt() != 4) continue; // only tetrahedra so far...
		if(data->hidden) continue;
		// store vertices
		for(int j = 0; j < data->pts.countInt(); j++){
			ofsv << data->pts[j].x << ' ' << data->pts[j].y << ' ' << data->pts[j].z << endl;
		}
		// store faces
		for(int j = 0; j < data->normals.countInt(); j++){
			for(int k = 0; k < 3; k++)
				ofsf << (vert_index + data->indices[j*3 + k]) << ' ';
			ofsf << endl;
			ofsc << (data->part >= 0 ? 2 : 1) << endl;
		}
		vert_index += data->pts.countInt();
	}

	// faces
	count = m_faces.countInt();
	bool extra_inner_faces = m_blocks.countInt() == 0;
	if(extra_inner_faces) LOG4CPLUS_INFO(MeshLog::logger_console, "Storing extra inner faces");
	const float INNER_FACTOR = (float) (-1e-5 * getBoundingBox().getDiameter() );
	for(int i = 0; i < count; i++){
		const auto& data = m_faces[i];
		if(quality_clip < 1.05 && data->quality > quality_clip)
			continue;
		bool clipped = false;
		for(int j = 0; j < data->pts.countInt(); j++)
			if(data->pts[j].x > clip.x || data->pts[j].y > clip.y || data->pts[j].z > clip.z){
				clipped = true; break;
			}
		if(clipped) continue;
		// store vertices
		for(int j = 0; j < data->pts.countInt(); j++){
			ofsv << data->pts[j].x << ' ' << data->pts[j].y << ' ' << data->pts[j].z << endl;
		}
		// store face
		for(int j = 0; j < data->pts.countInt(); j++){
			ofsf << (vert_index+j) << ' ';
		}
		ofsf << endl;
		ofsc << (data->part >= 0 ? 2 : 1) << endl;
		vert_index += data->pts.countInt();
		if(extra_inner_faces){
			const FVector3d nv = data->normal * INNER_FACTOR;
			// store vertices
			for(int j = 0; j < data->pts.countInt(); j++){
				ofsv << (data->pts[j].x + nv.x) << ' ' 
					 << (data->pts[j].y + nv.y) << ' ' 
					 << (data->pts[j].z + nv.z) << endl;
			}
			// store face
			for(int j = 0; j < data->pts.countInt(); j++){
				ofsf << (vert_index+j) << ' ';
			}
			ofsf << endl;
			ofsc << 0 << endl;
			vert_index += data->pts.countInt();
		}
	}

	return true;
}

bool MeshViewSet::storeMatlabFile(const string& fpath, const string& fname) const
{
	ofstream os((fpath+"/"+fname).c_str());

	os << "clear all; figure(1); clf; hold on; axis equal; axis off; view(3);" << endl;
	os << "c(1,:) = [0.9,0.9,0.9];" << endl;
	os << "c(2,:) = [0.7,0.7,0.7];" << endl;
	os << "colormap(c);" << endl;
	os << "load 'mesh_vertices.txt' -ascii" << endl;
	os << "load 'mesh_faces.txt' -ascii" << endl;
	os << "load 'mesh_ci.txt' -ascii" << endl;
	os << "patch('Vertices',mesh_vertices,'faces',mesh_faces,'FaceVertexCData',mesh_ci(:,1),'FaceColor','flat');" << endl;

	ofstream ofsv((fpath+"/mesh_vertices.txt").c_str());
	ofstream ofsf((fpath+"/mesh_faces.txt").c_str());
	ofstream ofsc((fpath+"/mesh_ci.txt").c_str());

	// blocks
	int count = m_blocks.countInt();
	int clipped = 0;
	int hidden = 0;

	int vert_index = 1;
	for(int i = 0; i < count; i++){
		const auto& data = m_blocks[i];
		if(data->normals.countInt() != 4) continue; // only tetrahedra so far...

		if(data->hidden == 1) ++clipped;
		if(data->hidden == 2) ++hidden;
		if(data->hidden) continue;

		// store vertices
		for(int j = 0; j < data->pts.countInt(); j++){
			ofsv << data->pts[j].x << ' ' << data->pts[j].y << ' ' << data->pts[j].z << endl;
		}
		// store faces
		for(int j = 0; j < data->normals.countInt(); j++){
			for(int k = 0; k < 3; k++)
				ofsf << (vert_index + data->indices[j*3 + k]) << ' ';
			ofsf << endl;
			//ofsc << (data->part >= 0 ? 2 : 1) << endl;
			ofsc << (data->adjacent[j] ? 2 : 1) << endl;
		}
		vert_index += data->pts.countInt();
	}

	if(true){
		ostringstream osstr;
		osstr << "Store blocks total/visible/clipped/hidden = " << count << '/' 
			<< (count - clipped - hidden) << '/' << clipped << '/' << hidden;
		LOG4CPLUS_INFO(MeshLog::logger_console, osstr.str());
	}

	// faces
	count = m_faces.countInt();
	bool extra_inner_faces = m_blocks.countInt() == 0;
	if(extra_inner_faces) LOG4CPLUS_INFO(MeshLog::logger_console, "Storing extra inner faces");
	const float INNER_FACTOR = (float) (-1e-5 * getBoundingBox().getDiameter() );
	for(int i = 0; i < count; i++){
		const auto& data = m_faces[i];
		if(data->hidden) continue;

		// store vertices
		for(int j = 0; j < data->pts.countInt(); j++){
			ofsv << data->pts[j].x << ' ' << data->pts[j].y << ' ' << data->pts[j].z << endl;
		}
		// store face
		for(int j = 0; j < data->pts.countInt(); j++){
			ofsf << (vert_index+j) << ' ';
		}
		ofsf << endl;
		ofsc << (data->part >= 0 ? 2 : 1) << endl;
		vert_index += data->pts.countInt();
		if(extra_inner_faces){
			const FVector3d nv = data->normal * INNER_FACTOR;
			// store vertices
			for(int j = 0; j < data->pts.countInt(); j++){
				ofsv << (data->pts[j].x + nv.x) << ' ' 
					 << (data->pts[j].y + nv.y) << ' ' 
					 << (data->pts[j].z + nv.z) << endl;
			}
			// store face
			for(int j = 0; j < data->pts.countInt(); j++){
				ofsf << (vert_index+j) << ' ';
			}
			ofsf << endl;
			ofsc << 0 << endl;
			vert_index += data->pts.countInt();
		}
	}

	return true;
}

void MeshViewSet::showViewSet(const string& desc, MeshViewSet* set, int t)
{
	if(param_show_visualization == 0) return;
	if(m_view) m_view->showViewSet(desc, set, t);
}

void MeshViewSet::showViewSetNoReset(const string& desc, MeshViewSet* set, int t)
{
	if(param_show_visualization == 0) return;
	if(m_view) m_view->showViewSetNoReset(desc, set, t);
}

void MeshViewSet::showDebugMesh(const string& desc, const MeshContainer3d* mesh, 
		const MeshPoint3d* pt1, const MeshPoint3d* pt2, double radius, int t)
{
	if(param_show_visualization == 0) return;
	showViewSet(desc, mesh->getDebugViewSet(pt1, pt2, radius), t);
}

void MeshViewSet::showDebugMesh(const string& desc, const MeshContainer3d* mesh, 
		const MeshBlock* el1, const MeshBlock* el2, double radius, int t)
{
	if(param_show_visualization == 0) return;
	showViewSet(desc, mesh->getDebugViewSet(el1, el2, radius), t);
}

void MeshViewSet::showDebugMesh(const string& desc, 
		const MeshContainer2d* mesh, const MeshPoint2d* pt1, 
		const MeshPoint2d* pt2, double radius, int t)
{
	if(param_show_visualization == 0) return;
	showViewSet(desc, mesh->getDebugViewSet(pt1, pt2, radius), t);
}

void MeshViewSet::showDebugMesh(const string& desc, const MeshContainer2d* mesh, 
		const MeshElement* el1, const MeshElement* el2, double radius, int t)
{
	if(param_show_visualization == 0) return;
	showViewSet(desc, mesh->getDebugViewSet(el1, el2, radius), t);
}

void MeshViewSet::showDebugMesh(const MeshContainer3d* boundary, 
								const MeshContainer3d* mesh, int phase, int t)
{
	if(param_show_visualization == 0) return;
	stringstream sstr;
	sstr << "Boundary recovery (phase " << phase << ")";
	showViewSet(sstr.str(), mesh->getDebugViewSet(boundary), t);
}

//MeshViewSet::MeshViewSet(const MeshViewSet& set)
//	: m_points(set.m_points.countInt()), m_edges(set.m_edges.countInt()), 
//	  m_faces(set.m_faces.countInt()), m_blocks(set.m_blocks.countInt()), 
//	  m_labels(set.m_labels), 
//	  m_hash( std::max(1000, (set.m_points.countInt() + set.m_edges.countInt() + set.m_faces.countInt() + set.m_blocks.countInt())), nullptr),
//	  m_mesh2d(set.m_mesh2d), m_mesh3d(set.m_mesh3d), 
//	  polygon_fill(set.polygon_fill)
//{
//	for(int i = 0; i < set.m_points.countInt(); i++){
//		m_points.add(new MeshViewPointData(*set.m_points.getconst(i)));
//	}
//	for(int i = 0; i < set.m_edges.countInt(); i++){
//		m_edges.add(new MeshViewEdgeData(*set.m_edges.getconst(i)));
//	}
//	for(int i = 0; i < set.m_faces.countInt(); i++){
//		m_faces.add(new MeshViewFaceData(*set.m_faces.getconst(i)));
//	}
//	for(int i = 0; i < set.m_blocks.countInt(); i++){
//		m_blocks.add(new MeshViewBlockData(*set.m_blocks.getconst(i)));
//	}
//}
//
void MeshViewSet::clearLabels()
{
	m_labels.clear();
}

void MeshViewSet::addLabel(const FPoint3d& pt, const string& label)
{
	m_labels.add(LabelInfo(pt, label));
}

FPoint3d MeshViewSet::transProject(float trans_matrix[], const FPoint3d& pt)
{
	float trans_pt[4];
    //Modelview transform
    trans_pt[0]=trans_matrix[0]*pt.x+trans_matrix[4]*pt.y+trans_matrix[8]*pt.z+trans_matrix[12];  //w is always 1
    trans_pt[1]=trans_matrix[1]*pt.x+trans_matrix[5]*pt.y+trans_matrix[9]*pt.z+trans_matrix[13];
    trans_pt[2]=trans_matrix[2]*pt.x+trans_matrix[6]*pt.y+trans_matrix[10]*pt.z+trans_matrix[14];
    trans_pt[3]=trans_matrix[3]*pt.x+trans_matrix[7]*pt.y+trans_matrix[11]*pt.z+trans_matrix[15];
    //The result normalizes between -1 and 1
    if(trans_pt[3] == 0.0)        //The w value
        return FPoint3d::zero;
    trans_pt[3]=1.0f/trans_pt[3];
    //Perspective division
    //trans_pt[0]*=trans_pt[3];
    //trans_pt[1]*=trans_pt[3];
    //trans_pt[2]*=trans_pt[3];
    //Window coordinates
    //Map x, y to range 0-1
    //windowCoordinate[0]=(trans_pt[0]*0.5+0.5)*viewport[2]+viewport[0];
    //windowCoordinate[1]=(trans_pt[1]*0.5+0.5)*viewport[3]+viewport[1];
    //This is only correct when glDepthRange(0.0, 1.0)
    //windowCoordinate[2]=(1.0+trans_pt[2])*0.5;  //Between 0 and 1
    //return 1;
	return FPoint3d( trans_pt[0] * trans_pt[3], trans_pt[1] * trans_pt[3], trans_pt[2] * trans_pt[3] );
}

bool MeshViewSet::storeEPSFile(float trans_matrix[], float model_matrix[], int view_mode, const string& fname) const
{
#ifdef USE_BOARD
	LOG4CPLUS_INFO(MeshLog::logger_console,"Storing EPS file ...");

	Board board;
	board.setClippingRectangle(-1, 1, 2, 2);
	board.setLineWidth(0.1);
	//board.setFontSize(10.0);

	enum {VIEW_BLOCKS = 1, VIEW_FACES = 2, VIEW_EDGES = 4, VIEW_NODES = 8, VIEW_WHITE = 16, STORE_BACK_FACES = 32, VIEW_LABELS = 64 };
	const Color colors_faces[]  = { Color::Blue, Color::Green, Color::Red, Color::Purple, Color::Yellow, Color::Cyan };
    const Color colors_points[] = { Color::Black, Color::Yellow, Color::Cyan, Color::Red, Color::Green, Color::Blue };
	const Color colors_edges[]  = { Color::Yellow, Color::Yellow, Color::Cyan, Color::Purple };
	const Color color_wedge  = Color::Gray;
	const Color color_light_gray(200, 200, 200);
	const int MDEPTH = 10000;
	const double DT = 0.1;
	// light (diffuse) ...
	const FVector3d lv = FVector3d(-2.0, -2.0, 5.0).normalized();

	// blocks
	if( (view_mode & VIEW_BLOCKS) != 0){
		int count = m_blocks.countInt();

		auto mode = getPolygonFillMode();

		for(int i = 0; i < count; i++){
			const auto& data = m_blocks[i];
			if(data->hidden) continue;

			DataVector<FPoint3d> block_v(data->pts.countInt());
			for(int j = 0; j < data->pts.countInt(); j++)
				block_v.add( MeshViewSet::transProject(trans_matrix, data->pts[j] ) );

			// -> check if the whole block outside clip
			int outside[4] = {0};
			for(int j = 0; j < data->pts.countInt(); j++){
				const FPoint3d& pt = block_v[j];
				if(pt.x < -1) outside[0]++;
				else if(pt.x >  1) outside[1]++;
				if(pt.y < -1) outside[2]++;
				else if(pt.y >  1) outside[3]++;
			}
			if( outside[0] == data->pts.countInt() ||
				outside[1] == data->pts.countInt() ||
				outside[2] == data->pts.countInt() || 
				outside[3] == data->pts.countInt()) continue;

			// -> set color
			Color block_color = Color::Gray;
			if(mode == MeshViewSet::FILL_QUALITY){
				float q = (float)data->quality;
				if(q > 0.6){ // blue - green
					block_color.setRGBf(0.0f, (1.0f - q)*5.0f, (q - 0.8f)*5.0f );
				}else if(q > 0.4){ /// green - yellow
					block_color.setRGBf( (0.8f - q)*(10.0f/3.0f), 1.0f, 0.0f );
				}else if(q > 0.2){ // yellow - orange
					block_color.setRGBf( 1.0f, 0.5f + (q - 0.2f)*(5.0f/3.0f), 0.0f );
				}else if(q > -1e-10){ // orange - red
					block_color.setRGBf( 1.0f,  q * 2.5f, 0.0f );
				}else{ // red
					block_color = Color::Red;
				}
			}else if(mode == MeshViewSet::FILL_NODES){
				assert(false);
			}else{
				assert(mode == MeshViewSet::FILL_AREA);
				if(data->part > 0){
					block_color =  Color::Lime;
				}else if(data->area_id < 0){
					block_color =  Color::Gray;
				}else{
					int cid = data->area_id;//+data->count-3;
					if(cid < 6){
						block_color = colors_faces[ cid % 6 ];
					}else{
						RandomGen rg(cid+2);
						block_color.setRGBf( (float)rg.doub(), (float)rg.doub(), (float)rg.doub() );
					}
				}
			}

			for(int k = 0; k < data->normals.countInt(); k++){
				// normal
				FVector3d fv = data->normals[k];
				fv = FVector3d(
						model_matrix[0]*fv.x+model_matrix[4]*fv.y+model_matrix[8]*fv.z,
						model_matrix[1]*fv.x+model_matrix[5]*fv.y+model_matrix[9]*fv.z,
						model_matrix[2]*fv.x+model_matrix[6]*fv.y+model_matrix[10]*fv.z ).normalized();
				if(fv.z < 0.0)	continue; // no need to draw...

				// face
				float sc = lv.scalarProduct(fv);
				if(sc < 0.0f) sc = -sc;
				sc = 0.5f + 0.5f * sc;

				Color face_color(
					(unsigned char) (sc * block_color.red()),
					(unsigned char) (sc * block_color.green()),
					(unsigned char) (sc * block_color.blue()));

				board.setPenColor( face_color );

				DataVector<FPoint3d> poly(3);
				std::vector<Point> polyline;
				for(int j = 0; j < 3; j++){
					const FPoint3d pt = MeshViewSet::transProject(trans_matrix, 
						data->pts[data->indices[k*3 + j]] );
					poly.add( pt );
					polyline.push_back( Point(pt.x, pt.y) );
				}

				double dmax = poly[0].z;
				double dmin = dmax;
				for(int j = 1; j < poly.countInt(); j++){
					if(poly[j].z > dmax) dmax = poly[j].z;
					if(poly[j].z < dmin) dmin = poly[j].z;
				}

				//double d = dmax;
				double d = (1-DT) * dmax + DT * dmin;

				// -> draw filled polyline
				board.fillPolyline( polyline, (int) (MDEPTH*(d+1.0)));			
			}
		}
	}

	// faces
	bool store_back_faces = ( ( view_mode & STORE_BACK_FACES) != 0 );
	if( (view_mode & VIEW_FACES) != 0 ){
		int count = m_faces.countInt();
		auto mode = getPolygonFillMode();
		for(int i = 0; i < count; i++){
			const auto& data = m_faces[i];
			if(data->hidden) continue;

			DataVector<FPoint3d> poly(data->pts.countInt());
			std::vector<Point> polyline;
			for(int j = 0; j < data->pts.countInt(); j++){
				FPoint3d pt = MeshViewSet::transProject(trans_matrix, data->pts[j] );
				poly.add( pt );
				polyline.push_back( Point(pt.x, pt.y) );
			}

			// -> check if the whole face outside clip
			int outside[4] = {0};
			for(int j = 0; j < data->pts.countInt(); j++){
				const FPoint3d& pt = poly[j];
				if(pt.x < -1) outside[0]++;
				else if(pt.x >  1) outside[1]++;
				if(pt.y < -1) outside[2]++;
				else if(pt.y >  1) outside[3]++;
			}
			if( outside[0] == data->pts.countInt() || 
				outside[1] == data->pts.countInt() ||
				outside[2] == data->pts.countInt() || 
				outside[3] == data->pts.countInt()) continue;

			// -> set color
			Color face_color = Color::Gray;
			if(mode == MeshViewSet::FILL_QUALITY){
				float q = (float)data->quality;
				if(q > 0.6){ // blue - green
					face_color.setRGBf(0.0f, (1.0f - q)*5.0f, (q - 0.8f)*5.0f );
				}else if(q > 0.4){ /// green - yellow
					face_color.setRGBf( (0.8f - q)*(10.0f/3.0f), 1.0f, 0.0f );
				}else if(q > 0.2){ // yellow - orange
					face_color.setRGBf( 1.0f, 0.5f + (q - 0.2f)*(5.0f/3.0f), 0.0f );
				}else if(q > -1e-10){ // orange - red
					face_color.setRGBf( 1.0f,  q * 2.5f, 0.0f );
				}else{ // red
					face_color = Color::Red;
				}
			}else if(mode == MeshViewSet::FILL_NODES){
				assert(false);
			}else if (mode == MeshViewSet::FILL_LGRAY){
				face_color = color_light_gray;
			}else{
				assert(mode == MeshViewSet::FILL_AREA);
				if (data->part > 0) {
					face_color = Color::Lime;
				}
				else if (data->area_id < 0) {
					face_color = Color::Gray;
				}
				else {
					int cid = data->area_id;//+data->count-3;
					if (cid < 6) {
						face_color = colors_faces[cid % 6];
					}
					else {
						RandomGen rg(cid + 2);
						face_color.setRGBf((float)rg.doub(), (float)rg.doub(), (float)rg.doub());
					}
				}
			}

			FVector3d fv = data->normal;
			fv = FVector3d(
					model_matrix[0]*fv.x+model_matrix[4]*fv.y+model_matrix[8]*fv.z,
					model_matrix[1]*fv.x+model_matrix[5]*fv.y+model_matrix[9]*fv.z,
					model_matrix[2]*fv.x+model_matrix[6]*fv.y+model_matrix[10]*fv.z ).normalized();

			if(fv.z < 0.0) {
				if( store_back_faces ) face_color = Color::Gray;
				else continue;
			}

			float sc = lv.scalarProduct(fv);

//			if(sc > 0.5) face_color = Color::Lime;

			if(sc < 0.0f) sc = -sc;
			sc = 0.5f + 0.5f * sc;

			face_color.setRGBi(
				(unsigned char) (sc * face_color.red()),
				(unsigned char) (sc * face_color.green()),
				(unsigned char) (sc * face_color.blue()));

			board.setPenColor( face_color );

			//float d = poly[0].z;
			//for(int j = 1; j < data->pts.countInt(); j++){
			//	if(poly[j].z > d) d = poly[j].z;
			//}

			double dmax = poly[0].z;
			double dmin = dmax;
			for(int j = 1; j < poly.countInt(); j++){
				if(poly[j].z > dmax) dmax = poly[j].z;
				if(poly[j].z < dmin) dmin = poly[j].z;
			}

			//double d = dmax;
			double d = (1-DT) * dmax + DT * dmin;

			// -> draw filled polyline
			board.fillPolyline( polyline, (int) (MDEPTH*(d+1.0)));
		}
	}

	if((view_mode & VIEW_EDGES) != 0){
		int count = m_edges.countInt();
		for(int i = 0; i < count; i++){
			const auto& data = m_edges[i];
			if(data->hidden) continue;

			if (getPolygonFillMode() == MeshViewSet::FILL_LGRAY) {
				board.setPenColor(Color::Black);
			}else if(data->part > 0){
				board.setPenColor( colors_points[ data->part % 6 ] );
			}else if(data->part < 0){
				board.setPenColor( color_wedge );
			}else{
				assert(data->border < 4);
				board.setPenColor( colors_edges[ data->border ] );
			}
			FPoint3d pt1 = MeshViewSet::transProject(trans_matrix, data->pt1 );
			FPoint3d pt2 = MeshViewSet::transProject(trans_matrix, data->pt2 );
			if( pt1.x < -1 && pt2.x < -1) continue;
			if( pt1.x >  1 && pt2.x >  1) continue;
			if( pt1.y < -1 && pt2.y < -1) continue;
			if( pt1.y >  1 && pt2.y >  1) continue;

			//double d = (pt1.z + pt2.z)/2;
			double d = (pt1.z > pt2.z) ? ((1-DT) * pt1.z + DT * pt2.z) : ((1-DT) * pt2.z + DT * pt1.z);
			board.drawLine( pt1.x, pt1.y, pt2.x, pt2.y, (int)(MDEPTH * (d+1.0)) );
		}
	}

	// vertices
	if((view_mode & VIEW_NODES) != 0){
		board.setLineWidth(1);
		int count = m_points.countInt();
		for(int i = 0; i < count; i++){
			const auto& data = m_points[i];
			if(data->hidden) continue;

			if(data->part > 0){
				board.setPenColor( colors_points[ data->part % 6 ] );
			}else if(data->part < 0){
				board.setPenColor( color_wedge );
			}else{
				assert(data->border < 4);
				board.setPenColor( colors_edges[ data->border ] );
			}

			FPoint3d pt = MeshViewSet::transProject(trans_matrix, data->pt );
			if( pt.x < -1 || pt.x > 1 || pt.y < -1 || pt.y > 1) continue;
			board.drawDot( pt.x, pt.y, (int) (MDEPTH * (pt.z + 1.0)) );
			if(data->numbered && (view_mode & VIEW_LABELS) != 0 ){
				board.drawText( pt.x, pt.y, to_string(data->id), (int) (MDEPTH * (pt.z + 1.0)) );
			}
		}
	}
	if( (view_mode & VIEW_LABELS) != 0 ) {
		if(m_labels.notEmpty()){
			board.setPenColor( Color::Red );
			for(int i = 0; i < m_labels.countInt(); i++){
				FPoint3d pt = MeshViewSet::transProject(trans_matrix, m_labels[i].pt);
				if( pt.x < -1 || pt.x > 1 || pt.y < -1 || pt.y > 1) continue;
				board.drawDot( pt.x, pt.y );
				board.drawText( pt.x, pt.y, m_labels[i].label, (int) (MDEPTH * (pt.z + 1.0)) );
			}
		}
	}

	board.scale(100);
	board.saveEPS( fname.c_str() );

	LOG4CPLUS_INFO(MeshLog::logger_console,"Ready");

#else
	LOG4CPLUS_ERROR(MeshLog::logger_console,   "Board library not included in this release...");
#endif // USE_BOARD
	return true;
}

void MeshViewSet::setHiddenByClipPlane(const MeshViewSet::ClipPlane& cp)
{
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Blocks:");
	int count = m_blocks.countInt();
	for (int i = 0; i < count; i++) {
		const auto& data = m_blocks[i];
		if (data->hidden) continue;
		for (int j = 0; j < data->pts.countInt(); j++)
			if (cp.clipped(data->pts[j])) {
				data->hidden = true;
				break;
			}
	}

	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Faces:");
	count = m_faces.countInt();
	for (int i = 0; i < count; i++) {
		const auto& data = m_faces[i];
		if (data->hidden) continue;
		for (int j = 0; j < data->pts.countInt(); j++)
			if (cp.clipped(data->pts[j])) {
				data->hidden = true;
				break;
			}
	}

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Edges:");
	count = m_edges.countInt();
	for (int i = 0; i < count; i++) {
		const auto& data = m_edges[i];
		if (data->hidden) continue;
		data->hidden = cp.clipped(data->pt1) ||	cp.clipped(data->pt2);
	}

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Vertices:");
	count = m_points.countInt();
	for (int i = 0; i < count; i++) {
		const auto& data = m_points[i];
		if (data->hidden) continue;
		data->hidden = cp.clipped(data->pt);
	}
}
