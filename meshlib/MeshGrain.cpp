/////////////////////////////////////////////////////////////////////////////
// MeshGrain.cpp: implementation of the MeshGrain class
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include <algorithm>

#include "MeshGrain.h"
#include "MeshViewSet.h"
#include "DMetric3d.h"
#include "DataMatrix.h"
#include "SurfacePlane.h"
#include "SurfaceCylinder.h"
#include "SurfaceBSplinePlanar.h"
#include "SurfaceBSplineCylindrical.h"
#include "SurfaceMulti.h"
#include "SurfaceCorrected.h"
#include "DataHashTable.h"
#include "DLeastSquaresFitting.h"
#include "DLine.h"
#include "Curve2dSegment.h"
#include "Curve2dBSpline.h"
#include "DSegment.h"
#include "DataList.h"
#include "MeshPoint3d.h"
#include "MeshGenerator3d.h"
#include "MeshGenerator3dDelaunayBoundary.h"
#include "MeshContainer3d.h"
#include "ControlSpace2dIdentity.h"
#include "ControlSpace3d.h"
#include "Metric3dContext.h"
#include "DPlane.h"

// tolerance for geometric approximation (plane, quadric, etc.)
double MeshGrain::param_grain_tolerance = 2.0;
// tolerance for how close can nodes be to be treated as different vertices
double MeshGrain::param_grain_node_identity_tolerance = 1.0;
// what to show
int MeshGrain::param_grain_view = 0;
// preset bounding box
DBox MeshGrain::param_bounding_box;

/// Transforms set of text files into .msh format
bool MeshGrain::parseGrainFiles(const string& dir, const string& /* msh_file */)
{
	DataVector<DPoint3d> nodes;
	DataVector<int> trans_nodes;
	DataVector<LineInfo> lines;
	DataVector<FaceInfo> faces;
	DataVector<BlockInfo> blocks;

	if(parseGrainNodeFile(dir, nodes, trans_nodes) == 0)
		return false;

	if(parseGrainLineFile(dir, nodes, trans_nodes, lines, faces) == 0)
		return false;
	if(parseGrainFaceFile(dir, nodes, trans_nodes, lines, faces, blocks) == 0)
		return false;

	checkTopology(nodes, lines, faces, blocks);

//	createSurfaceMeshes(faces, lines, nodes);

	classifyGrainSurfaces(faces, lines, nodes);

	//DataVector<MeshViewSet*> before(faces.countInt(), nullptr);
	//if(true){
	//	for(int i = 0; i < faces.countInt(); i++){
	//		const FaceInfo& face = faces[i];
	//		if(!face.valid) continue;
	//		before[i] = getFaceViewSet(face, nodes, lines);
	//	}
	//}

	//if(true){
	//	MeshViewSet* set = nullptr;
	//	BlockInfo& bl = blocks[31];
	//	for(int i = 0; i < bl.faces.countInt(); i++){
	//		const FaceInfo& face = faces[bl.faces[i]];
	//		if(face.valid){
	//			set = getFaceViewSet(face, nodes, lines, set);
	//		}
	//	}
	//	SHOW_MESH("grain #31", set);
	//}

	while(true){
		optimizeNodePlacement(faces, nodes); //-> optimize nodes
		if(recheckCloseNodes(blocks, faces, lines, nodes, trans_nodes) == 0)
			break;
	}
//	recheckPlanarFaceVertices(blocks, faces, lines, nodes);
	adjustSurfacesForNodes(faces, nodes);

	if((param_grain_view & 1) == 1){
	//if(true){
		for(int i = 0; i < faces.countInt(); i++){

	//		//if(i < before.countInt() && before[i]){
	//		//	ostringstream oss;
	//		//	oss << "face #" << i << " before optimized nodes";
	//		//	SHOW_MESH(oss.str(), before[i]);
	//		//}

			const FaceInfo& face = faces[i];
			if(!face.valid) continue;

			if(true){
				ostringstream oss;
				oss << "face #" << i << " with optimized nodes";
				for(int j = 0; j < face.lines.countInt(); j++){
					LineInfo& linfo = lines[face.lines[j]];
					if(linfo.adjacent_nonplanar > 0)
						oss << " [" << linfo.node0 << "-" << linfo.node1 
							<< "]x" << linfo.adjacent_nonplanar;
				}
				SHOW_MESH(oss.str(),
					getFaceViewSet(face, nodes, lines));
			}
		}
	}

	classifyGrainLines(faces, lines, nodes);
	classifyGrainBlocks(blocks, faces, nodes);

	storeGrainXml(dir, nodes, trans_nodes, lines, faces, blocks);
	//storeGrainXmlWithPlanes(dir, nodes, trans_nodes, lines, faces, blocks);
	return false;
}

/// Load list of nodes
int MeshGrain::parseGrainNodeFile(const string& dir, DataVector<DPoint3d> & nodes, DataVector<int> & trans_nodes)
{
	string fname = dir+"Nodes.txt";
	ifstream ifs(fname.c_str());
	string line;
	getline(ifs, line);
	if(!ifs){
		//LOG4CPLUS_WARN(MeshLog::logger_console, "Error opening file", fname);
		return parseGrainNodesShortFile(dir, nodes, trans_nodes);
	}

	int node_id = -1, node_ct;
	nodes.clear();
	DataVector<bool> available_nodes;
	while(ifs){
		string text;
		if(line.find("Number of Nodes:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> text >> node_ct;
			LOG4CPLUS_DEBUG(MeshLog::logger_console, "grain nodes, count= " << node_ct);
			nodes.addItems(node_ct+1, DPoint3d::zero);
			available_nodes.addItems(node_ct+1, false);
			trans_nodes.addItems(node_ct+1, 0);
		}
		if(line.find("Node number:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> node_id;
			if(node_id < 1 || node_id > node_ct)
				LOG4CPLUS_ERROR(MeshLog::logger_console, "Invalid node id: " << node_id);
		}
		if(line.find("Node coordinates:") != string::npos){
			if(node_id > 0){
				istringstream iss(line);
				DPoint3d pt;
				iss >> text >> text >> pt.x >> pt.y >> pt.z;
				nodes[node_id] = pt;
				available_nodes[node_id] = true;
				param_bounding_box.addPoint(pt);
			}else
				LOG4CPLUS_ERROR(MeshLog::logger_console,   "Node coordinates without proper id");
			node_id = -1;
		}
		getline(ifs, line);
	}

	for(int i = 1; i <= node_ct; i++)
		if(!available_nodes[i]) 
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Missing node coordinates for id: " << i);

	double DIST_TOL2 = param_grain_node_identity_tolerance * param_grain_node_identity_tolerance;
	trans_nodes[1] = 1;
	int close_node_count = 0;
	for(int i = 2; i <= node_ct; i++){
		trans_nodes[i] = i;
		// check if close to other nodes ?
		for(int j = 1; j < i; j++){
			double dist2 = nodes[i].distance2(nodes[j]);
			if(dist2 < DIST_TOL2){
				++close_node_count;
				if(trans_nodes[i] > trans_nodes[j])
					trans_nodes[i] = trans_nodes[j];
				else
					trans_nodes[j] = trans_nodes[i];
			}
		}
	}

	DataVector<int> pts_counts(node_ct+1, 0);
	DataVector<DPoint3d> pts_ave(node_ct+1, DPoint3d(0.0,0.0,0.0));

	for(int i = 1; i <= node_ct; i++){
		pts_counts[trans_nodes[i]]++;
		pts_ave[trans_nodes[i]].add(nodes[i]);
		if(trans_nodes[i] != i){
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Fixing close nodes " << i << " -> " << trans_nodes[i]);
		}
	}

	if(close_node_count){
		LOG4CPLUS_WARN(MeshLog::logger_console, "Total number of joined nodes: " << close_node_count);
	}

	for(int i = 1; i <= node_ct; i++){
		if(pts_counts[i] > 1){
			nodes[i] = pts_ave[i];
			nodes[i] /= pts_counts[i];
		}
	}

	return node_ct;
}

/// Load list of nodes
int MeshGrain::parseGrainNodesShortFile(const string& dir, DataVector<DPoint3d> & nodes, DataVector<int> & trans_nodes)
{
	string fname = dir+"NodesShort.txt";
	ifstream ifs(fname.c_str());
	string line;
	getline(ifs, line);
	if(!ifs){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening node file: " << fname);
		return 0;
	}

	int node_id = -1, node_ct;
	nodes.clear();
	DataVector<bool> available_nodes;
	while(ifs){
		string text;
		if(line.find("Number of Nodes:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> text >> node_ct;
			LOG4CPLUS_DEBUG(MeshLog::logger_console, "grain nodes count: " << node_ct);
			nodes.addItems(node_ct+1, DPoint3d());
			available_nodes.addItems(node_ct+1, false);
			trans_nodes.addItems(node_ct+1, 0);
		}else{
			// id + x + y + z
			istringstream iss(line);
			DPoint3d pt;
			iss >> node_id >> pt.x >> pt.y >> pt.z;
			if(iss){
				nodes[node_id] = pt;
				available_nodes[node_id] = true;
				param_bounding_box.addPoint(pt);
			}
		}
		getline(ifs, line);
	}

	for(int i = 1; i <= node_ct; i++)
		if(!available_nodes[i]) 
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Missing node coordinates for id: " << i);

	double DIST_TOL2 = param_grain_node_identity_tolerance * param_grain_node_identity_tolerance;
	trans_nodes[1] = 1;
	int close_node_count = 0;

	for(int i = 2; i <= node_ct; i++){
		trans_nodes[i] = i;
		// check if close to other nodes ?
		for(int j = 1; j < i; j++){
			double dist2 = nodes[i].distance2(nodes[j]);
			if(dist2 < DIST_TOL2){
				++close_node_count;
				if(trans_nodes[i] > trans_nodes[j])
					trans_nodes[i] = trans_nodes[j];
				else
					trans_nodes[j] = trans_nodes[i];
			}
		}
	}

	DataVector<int> pts_counts(node_ct+1, 0);
	DataVector<DPoint3d> pts_ave(node_ct+1, DPoint3d(0.0,0.0,0.0));

	for(int i = 1; i <= node_ct; i++){
		pts_counts[trans_nodes[i]]++;
		pts_ave[trans_nodes[i]].add(nodes[i]);
		if(trans_nodes[i] != i){
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Fixing close nodes " << i << " -> " << trans_nodes[i]);
		}
	}

	if(close_node_count){
		LOG4CPLUS_WARN(MeshLog::logger_console, "Total number of joined nodes: " << close_node_count);
	}

	for(int i = 1; i <= node_ct; i++){
		if(pts_counts[i] > 1){
			nodes[i] = pts_ave[i];
			nodes[i] /= pts_counts[i];
		}
	}

	return node_ct;
}

/// Load list of lines
int MeshGrain::parseGrainLineFile(const string& dir, 		
			DataVector<DPoint3d> & nodes, const DataVector<int> & trans_nodes, 
			DataVector<LineInfo> & lines, DataVector<FaceInfo> & faces)
{
	string fname = dir+"Lines.txt";
	ifstream ifs(fname.c_str());
	string line;
	getline(ifs, line);
	if(!ifs){
		//LOG4CPLUS_WARN(MeshLog::logger_console, "Error opening file", fname);
		return parseGrainLinesShortFile(dir, nodes, trans_nodes, lines, faces);
	}

	int line_id = -1, line_ct = -1, edge_ct = -1;
	LineInfo linfo;
	lines.clear();
	while(ifs){
		string text;
		if(line_ct < 1 && line.find("Number of Lines:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> text >> line_ct;
			LOG4CPLUS_DEBUG(MeshLog::logger_console, "grain lines count: " << line_ct);
			lines.addItems(line_ct+1);
		}
		if(line.find("Line number:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> line_id;
			linfo = LineInfo();
			if(line_id < 1){
				LOG4CPLUS_ERROR(MeshLog::logger_console, "Invalid line id: " << line_id);
				line_id = -1;
				linfo.valid = false;
			}else{
				if(line_id > line_ct){
					//LOG4CPLUS_WARN(MeshLog::logger_console, "Line number outside the declared count", line_id);
					lines.addItems(line_id - line_ct);
					line_ct = line_id;
				}
				linfo.valid = true;
			}
		}
		if(line_id > 0 && line.find("Nodes:") != string::npos){
			istringstream iss(line);
			iss >> text >> linfo.node0 >> linfo.node1;
			if(!iss || linfo.node0 < 1 || linfo.node1 < 1) 
				linfo.valid = false;
			else{
				int n0 = trans_nodes[linfo.node0];
				int n1 = trans_nodes[linfo.node1];
				if(n0 == n1) linfo.valid = false;
				else{
					linfo.node0 = std::min(n0, n1);
					linfo.node1 = std::max(n0, n1);
				}
			}
		}
		if(line_id > 0 && linfo.valid && line.find("Number of Planes:") != string::npos){
			istringstream iss(line);
			int plane_ct = 0;
			iss >> text >> text >> text >> plane_ct;
			getline(ifs, line);
			getline(ifs, line);
			istringstream iss_p(line);
			int p_ct = 0;
			int p_id = 0;
			iss_p >> p_id;
			while(iss_p){
				linfo.incident_faces.add(p_id);
				++p_ct;
				iss_p >> p_id;
			}
			if(p_ct < plane_ct){
				LOG4CPLUS_WARN(MeshLog::logger_console, "Incomplete list of planes for line, id: " << line_id);
			}else if(p_ct > plane_ct){
				LOG4CPLUS_WARN(MeshLog::logger_console, "Overfull list of planes for line, id: " << line_id);
			}
		}
		if(linfo.node0 > 0 && linfo.valid && line.find("Number of Edges:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> text >> edge_ct;
			linfo.inner_points.clear();
			linfo.inner_points.prepare(edge_ct);
			getline(ifs, line);
			getline(ifs, line);
			DPoint3d pt;
			for(int i = 0; i < edge_ct; i++){
				getline(ifs, line);
				istringstream iss_p(line);
				iss_p >> pt.x >> pt.y >> pt.z;
				if(!iss_p){
					LOG4CPLUS_WARN(MeshLog::logger_console, "Incomplete edge-data for line, id: " << line_id);
					//linfo.valid = false;
					break;
				}
				linfo.inner_points.add(pt);
			}
			// update
			lines[line_id] = linfo;
			// next line-data
			line_id = -1;
			linfo.node0 = -1;
			linfo.valid = false;
		}
		getline(ifs, line);
	}

	// check for duplicate lines (due to node-collapsing - but not necessary!)
	checkForDuplicateLines(lines, nodes);

	for(int i = 0; i < lines.countInt(); i++){
		if(lines[i].valid || lines[i].trans_line >= 0){
			// insert incidency data for faces
			for(int j = 0; j < lines[i].incident_faces.countInt(); j++){
				int p_id = lines[i].incident_faces[j];
				if(faces.countInt() <= p_id){
					faces.addItems(p_id - faces.countInt() + 1);
				}
				int tline_id = (lines[i].trans_line < 0) ? i : lines[i].trans_line; 
				faces[p_id].lines.add(tline_id);
			}
		}
	}

	return line_ct;
}

/// Load list of lines
int MeshGrain::parseGrainLinesShortFile(const string& dir, 
		DataVector<DPoint3d> & nodes, const DataVector<int> & trans_nodes, 
		DataVector<LineInfo> & lines, DataVector<FaceInfo> & faces)
{
	string fname = dir+"LinesShort.txt";
	ifstream ifs(fname.c_str());
	string line;
	getline(ifs, line);
	if(!ifs){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening edge file: " << fname);
		return 0;
	}

	int line_id = -1, line_ct = -1, edge_ct = -1;
	LineInfo linfo;
	lines.clear();
	while(ifs){
		string text;
		if(line_ct < 1 && line.find("Number of Lines:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> text >> line_ct;
			LOG4CPLUS_DEBUG(MeshLog::logger_console, "grain lines count: " << line_ct);
			lines.addItems(line_ct+1);
		}
		if(line.find("Line number:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> line_id;
			linfo = LineInfo();
			if(line_id < 1){
				LOG4CPLUS_ERROR(MeshLog::logger_console, "Invalid line id: " << line_id);
				line_id = -1;
				linfo.valid = false;
			}else{
				if(line_id > line_ct){
					//LOG4CPLUS_WARN(MeshLog::logger_console, "Line number outside the declared count", line_id);
					lines.addItems(line_id - line_ct);
					line_ct = line_id;
				}
				linfo.valid = true;
			}
		}
		if(line_id > 0 && line.find("Nodes:") != string::npos){
			istringstream iss(line);
			iss >> text >> linfo.node0 >> linfo.node1;
			if(!iss || linfo.node0 < 1 || linfo.node1 < 1) 
				linfo.valid = false;
			else{
				int n0 = trans_nodes[linfo.node0];
				int n1 = trans_nodes[linfo.node1];
				if(n0 == n1) linfo.valid = false;
				else{
					linfo.node0 = std::min(n0, n1);
					linfo.node1 = std::max(n0, n1);
				}
			}
		}
		if(line_id > 0 && linfo.valid && line.find("Number of Planes:") != string::npos){
			istringstream iss(line);
			int plane_ct = 0;
			iss >> text >> text >> text >> plane_ct;
			getline(ifs, line);
			istringstream iss_p(line);
			int p_ct = 0;
			int p_id = 0;
			iss_p >> text >> p_id;
			while(iss_p){
				linfo.incident_faces.add(p_id);
				++p_ct;
				iss_p >> p_id;
			}
			if(p_ct < plane_ct){
				LOG4CPLUS_WARN(MeshLog::logger_console, "Incomplete list of planes for line, id: " << line_id);
			}else if(p_ct > plane_ct){
				LOG4CPLUS_WARN(MeshLog::logger_console, "Overfull list of planes for line, id: " << line_id);
			}
		}
		if(linfo.node0 > 0 && linfo.valid && line.find("Number of Edges:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> text >> edge_ct;
			linfo.inner_points.clear();
			linfo.inner_points.prepare(edge_ct);
			getline(ifs, line);
			getline(ifs, line);
			DPoint3d pt;
			for(int i = 0; i < edge_ct; i++){
				getline(ifs, line);
				istringstream iss_p(line);
				iss_p >> pt.x >> pt.y >> pt.z;
				if(!iss_p){
					LOG4CPLUS_WARN(MeshLog::logger_console, "Incomplete edge-data for line, id: " << line_id);
					//linfo.valid = false;
					break;
				}
				linfo.inner_points.add(pt);
			}
			// update
			lines[line_id] = linfo;
			// next line-data
			line_id = -1;
			linfo.node0 = -1;
			linfo.valid = false;
		}
		getline(ifs, line);
	}

	// check for duplicate lines (due to node-collapsing - but not necessary!)
	checkForDuplicateLines(lines, nodes);

	int valid_lines = 0;
	for(int i = 0; i < lines.countInt(); i++){
		if(lines[i].valid || lines[i].trans_line >= 0){
			++valid_lines;
			// insert incidency data for faces
			for(int j = 0; j < lines[i].incident_faces.countInt(); j++){
				int p_id = lines[i].incident_faces[j];
				if(faces.countInt() <= p_id){
					faces.addItems(p_id - faces.countInt() + 1);
				}
				int tline_id = (lines[i].trans_line < 0) ? i : lines[i].trans_line; 
				faces[p_id].lines.add(tline_id);
			}
		}
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "valid grain edges: " << valid_lines);

	return line_ct;
}

// check for duplicate lines (due to node-collapsing - but not necessary!)
bool MeshGrain::checkForDuplicateLines(
		DataVector<LineInfo> & lines, DataVector<DPoint3d> & nodes)
{
	DataHashTableKeyValue<int, int> node_pairs(2*lines.countInt(), -1);

	int lct = lines.countInt();
	for(int i = 0; i < lct; i++){
		LineInfo& line1 = lines[i];
		if(!line1.valid) continue;
		int pair_id = line1.node0 * lct + line1.node1;
		int i2 = node_pairs.getValue(pair_id, -1);
		if(i2 < 0){
			node_pairs.insert(pair_id, i);
			continue;
		}
		// found double
		LineInfo& line2 = lines[i2];
		// ... select: collapse/split
		// ... check max inner-node distance from vertex-segment for both edges
		const DLine3d segment(nodes[line1.node0], nodes[line1.node1]);
		double max_dist2_1 = 0.0;
		int max_i1 = 0;
		for(int j = 0; j < line1.inner_points.countInt(); j++){
			double dist2 = segment.distanceToPoint2(line1.inner_points[j]);
			if(dist2 > max_dist2_1){
				max_dist2_1 = dist2;
				max_i1 = j;
			}
		}
		double max_dist2_2 = 0.0;
		int max_i2 = 0;
		for(int j = 0; j < line2.inner_points.countInt(); j++){
			double dist2 = segment.distanceToPoint2(line2.inner_points[j]);
			if(dist2 > max_dist2_2){
				max_dist2_2 = dist2;
				max_i2 = j;
			}
		}
		bool linear_1 = (line1.inner_points.countInt() < 5) || (max_dist2_1 < param_grain_tolerance);
		bool linear_2 = (line2.inner_points.countInt() < 5) || (max_dist2_2 < param_grain_tolerance);
		// ... show
		if((param_grain_view & 4) > 0){
			MeshViewSet* set = new MeshViewSet;
			set->addEdge(nodes[line1.node0], nodes[line1.node1], 0);
			set->addPoint(nodes[line1.node0], 0, line1.node0);
			set->addPoint(nodes[line1.node1], 0, line1.node1);
			// -- inner nodes
			for(int j = 0; j < line1.inner_points.countInt(); j++)
				set->addPoint(line1.inner_points[j], 1);
			for(int j = 0; j < line2.inner_points.countInt(); j++)
				set->addPoint(line2.inner_points[j], 2);
			// that's all
			string caption = "edge collision ";
			caption += linear_1 ? "(1 linear) " : "(1 nonlinear) ";
			caption += linear_2 ? "(2 linear) " : "(2 nonlinear) ";
			SHOW_MESH(caption, set);
		}
		// ... decide
		if( linear_1 && linear_2){
			// both can be treated as linear -> join
			// ... mark translation-line
			line1.trans_line = i2;
			// ... join inner data
			for(int j = 0; j < line1.inner_points.countInt(); j++){
				line2.inner_points.add(line1.inner_points[j]);
			}
			// ... invalidate this one
			line1.valid = false;
		}else{
			// split the non-linear ones
			if(!linear_2){
				LineInfo& line = lines[i2];
				node_pairs.setValue(pair_id, -2); // non-empty, but also non-valid
				// split line2
				// ... select new node (most distant one)
				int node_id = nodes.add(line.inner_points.removeOrderedAt(max_i2));
				// ... create new line
				LineInfo new_line(line.node1, node_id, true);
				new_line.incident_faces = line.incident_faces; // the same incident faces
				line.node1 = node_id;
				// ... split inner points
				double tm = segment.paramForPoint(nodes[node_id]);
				for(int j = 0; j < line.inner_points.countInt(); ){
					if(segment.paramForPoint(line.inner_points[j]) > tm){
						new_line.inner_points.add(line.inner_points.removeOrderedAt(j));
					}else ++j;
				}
				lines.add(new_line);
			}
			if(linear_1){
				node_pairs.insert(pair_id, i);
			}else{
				// split line1
				LineInfo& line = lines[i]; // refresh, since might have changes due to earlier addition to lines
				// ... select new node (most distant one)
				int node_id = nodes.add(line.inner_points.removeOrderedAt(max_i1));
				// ... create new line
				LineInfo new_line(line.node1, node_id, true);
				new_line.incident_faces = line.incident_faces; // the same incident faces
				line.node1 = node_id;
				// ... split inner points
				double tm = segment.paramForPoint(nodes[node_id]);
				for(int j = 0; j < line.inner_points.countInt(); ){
					if(segment.paramForPoint(line.inner_points[j]) > tm){
						new_line.inner_points.add(line.inner_points.removeOrderedAt(j));
					}else ++j;
				}
				lines.add(new_line);
			}
		}
	}

	return true;
}

/// Load list of faces
int MeshGrain::parseGrainFaceFile(const string& dir, 
		const DataVector<DPoint3d> & nodes, const DataVector<int> & trans_nodes, 
		const DataVector<LineInfo> & lines, DataVector<FaceInfo> & faces, 
		DataVector<BlockInfo> & blocks, const string& filename)
{
	string fname = dir + filename;
	ifstream ifs(fname.c_str());
	string line;
	getline(ifs, line);
	if(!ifs){
		if(filename == "Planes.txt")
			return parseGrainFaceFile(dir, nodes, trans_nodes, lines, faces, blocks, "PlanesShort.txt");
		else
			LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening face file: " << fname);
	}

	int face_id = -1, face_ct = -1, points_ct = -1;
	MeshViewSet *set = nullptr;
	int empty_faces = 0;
	int valid_faces = 0;
	while(ifs){
		string text;
		if(face_ct < 1 && line.find("Number of planes:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> text >> face_ct;
			//if(face_ct != faces.countInt()-1){
			//	LOG4CPLUS_WARN(MeshLog::logger_console, "Number of faces in Plane.txt different from that in Lines.txt");
			//}
			LOG4CPLUS_DEBUG(MeshLog::logger_console, "grain faces count= " << face_ct);
			//faces.add(face_ct+1, FaceInfo());
		}
		if(line.find("Plane number:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> face_id;
			if(face_id < 1){
				LOG4CPLUS_ERROR(MeshLog::logger_console, "Invalid face id: " << face_id);
				face_id = -1;
			}else{
				if(face_id >= faces.countInt()){
					faces.addItems(face_id - faces.countInt() + 1);
				}
				FaceInfo& finfo = faces[face_id];
				finfo.valid = true;
				// fix lines 
				// ... (apply translation info and remove duplicates)
				DataVector<bool> used_line(lines.countInt(), false);
				for(int i = 0; i < finfo.lines.countInt(); ){
					if(lines[finfo.lines[i]].trans_line >= 0)
						finfo.lines[i] = lines[finfo.lines[i]].trans_line;
					if(used_line[finfo.lines[i]])
						finfo.lines.removeAt(i);
					else 
						used_line[finfo.lines[i++]] = true;
				}
				// fix nodes
				int line_count = finfo.lines.countInt();
				// check cycle
				DataVector<int> nodes_check(nodes.countInt(), 0);
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "*** FACE " << face_id << " ***");
				for(int i = 0; i < line_count; i++){
					int id = finfo.lines[i];
					nodes_check[lines[id].node0]++;
					nodes_check[lines[id].node1]++;
					LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "i = " << i << " -> " << lines[id].node0 << "," << lines[id].node1);
				}
				if((param_grain_view & 4) > 0){
					// show face
					set = new MeshViewSet;
					// -- face edges
					for(int j = 0; j < finfo.lines.countInt(); j++){
						const LineInfo& tline = lines[finfo.lines[j]];
						if(!tline.valid) continue;
						set->addEdge(nodes[tline.node0], nodes[tline.node1], 0);
						set->addPoint(nodes[tline.node0], 0, tline.node0);
						set->addPoint(nodes[tline.node1], 0, tline.node1);
						for(int k = 0; k < tline.inner_points.countInt(); k++)
							set->addPoint(tline.inner_points[k], 1);
					}
					//--------------------
				}
				bool any_change = true;
				while(any_change){ // find vertices with single edge
					any_change = false;
					for(int i = 0; i < nodes.countInt(); i++){
						if(nodes_check[i] == 1){
							LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
								" non-cyclic boundary edges for face id" <<  face_id);
							any_change = true;
							// remove edge with this node
							for(int j = 0; j < line_count; j++){
								int id = finfo.lines[j];
								if(lines[id].node0 == i || lines[id].node1 == i){
									nodes_check[lines[id].node0]--;
									nodes_check[lines[id].node1]--;
									line_count--;
									finfo.lines.removeAt(j);
									break;
								}
							}
						}
					}
				}
				any_change = true;
				while(any_change){ // find edges with odd vertex-ranks
					any_change = false;
					for(int j = 0; j < line_count; j++){
						int id = finfo.lines[j];
						if( (nodes_check[lines[id].node0] % 2 == 1) &&
							(nodes_check[lines[id].node1] % 2 == 1))
						{
							LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
								" non-cyclic boundary edges for face id" <<  face_id);
							nodes_check[lines[id].node0]--;
							nodes_check[lines[id].node1]--;
							line_count--;
							finfo.lines.removeAt(j);
							any_change = true;
							break;
						}
					}
				}
				for(int j = 0; j < line_count; ){
					int id = finfo.lines[j];
					if(lines[id].node0 == lines[id].node1){
						nodes_check[lines[id].node0] -= 2;
						line_count--;
						finfo.lines.removeAt(j);
					}else
						j++;
				}
				bool any_invalid_node = false;
				for(int i = 0; i < nodes.countInt(); i++){
					if(nodes_check[i] != 2 && nodes_check[i] != 0){
						LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
							"non-cyclic boundary edges, face " << 
							face_id << ", node " << i << ", count " << nodes_check[i]);
						any_invalid_node = true;
					}
				}
				if(any_invalid_node){
					finfo.valid = false;
				}else{
					// sort
					if(finfo.lines.countInt() > 2){
						int line_id = finfo.lines[0];
						int last_node = lines[line_id].node0;
						finfo.nodes.add(last_node);
						for(int i = 1; i < line_count; i++){
							for(int j = i; j < line_count; j++){
								line_id = finfo.lines[j];
								if(lines[line_id].node0 == last_node){
									// switch lines
									if(j > i) {
										finfo.lines[j] = finfo.lines[i];
										finfo.lines[i] = line_id;
									}
									last_node = lines[line_id].node1;
									finfo.nodes.add(last_node);
									break;
								}else if(lines[line_id].node1 == last_node){
									// switch lines
									if(j > i) {
										finfo.lines[j] = finfo.lines[i];
										finfo.lines[i] = line_id;
									}
									last_node = lines[line_id].node0;
									finfo.nodes.add(last_node);
									break;
								}
							}
							assert(finfo.nodes.countInt() == i+1);
						}
					}else{
						++empty_faces;
						finfo.valid = false;
					}
				}
			}
		}
		if(face_id > 0 && (line.find("Grain Number") != string::npos ||
			line.find("Grain number") != string::npos)){
			FaceInfo& finfo = faces[face_id];
			if(finfo.valid){ // insert information into blocks
				istringstream iss(line);
				iss >> text >> text >> finfo.block0 >> finfo.block1;
				if(finfo.block0 > 0) {
					if(blocks.countInt() <= finfo.block0) 
						blocks.addItems(finfo.block0 - blocks.countInt() + 1);
					blocks[finfo.block0].faces.add(face_id);
				}
				if(finfo.block1 > 0) {
					if(blocks.countInt() <= finfo.block1) 
						blocks.addItems(finfo.block1 - blocks.countInt() + 1);
					blocks[finfo.block1].faces.add(face_id);
				}
			}
		}
		//if(face_id > 0 && line.find("Number of lines:") != string::npos){
		//	istringstream iss(line);
		//	iss >> text >> text >> text >> lines_ct;
		//	finfo.lines.clear();
		//	finfo.lines.prepare(lines_ct);
		//	getline(ifs, line);
		//	getline(ifs, line);
		//	int id;
		//	istringstream iss_l(line);
		//	for(int i = 0; i < lines_ct; i++){
		//		iss_l >> id;
		//		if(!iss_l){
		//			finfo.valid = false;
		//			break;
		//		}
		//		finfo.lines.add(id);
		//	}
		//}
		//if(face_id > 0 && line.find("Number of Nodes:") != string::npos){
		//	istringstream iss(line);
		//	iss >> text >> text >> text >> nodes_ct;
		//	finfo.nodes.clear();
		//	finfo.nodes.prepare(nodes_ct);
		//	getline(ifs, line);
		//	getline(ifs, line);
		//	int id;
		//	istringstream iss_n(line);
		//	for(int i = 0; i < nodes_ct; i++){
		//		iss_n >> id;
		//		if(!iss_n){
		//			finfo.valid = false;
		//			break;
		//		}
		//		finfo.nodes.add(id);
		//	}
		//}
		if(face_id > 0 && line.find("Number of Faces:") != string::npos){
			istringstream iss(line);
			iss >> text >> text >> text >> points_ct;
			FaceInfo& finfo = faces[face_id];
			finfo.inner_points.clear();
			finfo.inner_points.prepare(points_ct);
			getline(ifs, line);
			getline(ifs, line);
			DPoint3d pt;
			for(int i = 0; i < points_ct; i++){
				getline(ifs, line);
				istringstream iss_p(line);
				iss_p >> pt.x >> pt.y >> pt.z;
				if(!iss_p){
					LOG4CPLUS_WARN(MeshLog::logger_console, "Incomplete face-data for plane, id:" << face_id);
					//finfo.valid = false;
					break;
				}
				finfo.inner_points.add(pt);
			}

			if((param_grain_view & 4) > 0){
				// -- face inner nodes
				for(int j = 0; j < finfo.inner_points.countInt(); j++)
					set->addPoint(finfo.inner_points[j], 3);

				ostringstream oss;
				oss << (finfo.valid ? "" : "invalid ") << "grain face #" << face_id << " (initial)";
				SHOW_MESH(oss.str(), set);
				//--------------------
			}

			if(finfo.valid) ++valid_faces;
			faces[face_id] = finfo;
			face_id = -1;

		}
		getline(ifs, line);
	}

	if(empty_faces > 0)
		LOG4CPLUS_WARN(MeshLog::logger_console, "Empty faces (after pruning), total: " << empty_faces);

	//for(int i = 1; i <= face_ct; i++){
	//	if(!faces[i].valid){
	//		LOG4CPLUS_ERROR(MeshLog::logger_console, "Invalid or missing data for face, id", i);
	//		continue;
	//	}
	//	//if(faces[i].nodes.countInt() != faces[i].lines.countInt()){
	//	//	LOG4CPLUS_ERROR(MeshLog::logger_console, "Inconsistent node/edge count for plane", i);
	//	//	LOG4CPLUS_WARN(MeshLog::logger_console, "Nodes", faces[i].nodes.countInt());
	//	//	LOG4CPLUS_WARN(MeshLog::logger_console, "Lines", faces[i].lines.countInt());
	//	//}
	//}

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "valid grain faces: " << valid_faces);

	//if(true){
	//	MeshViewSet* set = nullptr;
	//	for(int i = 1; i <= face_ct; i++){
	//		if(faces[i].valid){
	//			set = getFaceViewSet(faces[i], nodes, lines, set);
	//		}
	//	}
	//	SHOW_MESH("model faces", set);
	//}

	return valid_faces;
}

// adjust surfaces to better fit face-vertices
void MeshGrain::adjustSurfacesForNodes(
		DataVector<FaceInfo> & faces, 
		DataVector<DPoint3d> & nodes)
{
	// check boundary
	for(int i = 1; i < faces.countInt(); i++){
		FaceInfo& finfo = faces[i];
		if(!finfo.valid) continue;
		int nct = finfo.nodes.countInt();
		// check boundary ...
		if(finfo.block0 < 0 || finfo.block1 < 0){
			assert(finfo.surface->getType() == SURFACE_PLANE);
			const DVector3d vn = finfo.surface->getNormalVector(DPoint2d::zero);
			const DPoint3d ptzero = finfo.surface->getPoint(DPoint2d::zero);
			if(abs(vn.x) > abs(vn.y) && abs(vn.x) > abs(vn.z)){
				double xb = (abs(ptzero.x - param_bounding_box.x0) < abs(ptzero.x - param_bounding_box.x1)) ?
					param_bounding_box.x0 : param_bounding_box.x1;
				for(int j = 0; j < nct; j++)
					nodes[finfo.nodes[j]].x = xb;
			}else if(abs(vn.y) > abs(vn.z)){
				double yb = (abs(ptzero.y - param_bounding_box.y0) < abs(ptzero.y - param_bounding_box.y1)) ?
					param_bounding_box.y0 : param_bounding_box.y1;
				for(int j = 0; j < nct; j++)
					nodes[finfo.nodes[j]].y = yb;
			}else{
				double zb = (abs(ptzero.z - param_bounding_box.z0) < abs(ptzero.z - param_bounding_box.z1)) ?
					param_bounding_box.z0 : param_bounding_box.z1;
				for(int j = 0; j < nct; j++)
					nodes[finfo.nodes[j]].z = zb;
			}
		}
	}
	// adjust surfaces
	for(int i = 1; i < faces.countInt(); i++){
		FaceInfo& finfo = faces[i];
		if(!finfo.valid) continue;
		int nct = finfo.nodes.countInt();
		DataVector<DPoint3d> all_nodes(nct);
		for(int j = 0; j < nct; j++){
			all_nodes.add(nodes[finfo.nodes[j]]);
		}
		auto surface = finfo.surface;
		if(surface->getType() == SURFACE_PLANE){
			// check max distance
			double max_dist2 = 0.0;
			for(int j = 0; j < nct; j++){
				double dist2 = nodes[finfo.nodes[j]].distance2(
					surface->getPoint(surface->getParameters(nodes[finfo.nodes[j]])));
				if(dist2 > max_dist2) max_dist2 = dist2;
			}
			if(max_dist2 > 0.25){ // re-fit plane
				DPlane plane;
				DLeastSquaresFitting::fitHyperplaneOrthogonal(all_nodes, plane);
				finfo.surface = std::make_shared<SurfacePlane>(plane.p0, plane.e0, plane.e1);
			}
			if(nct > 3){
				max_dist2 = 0.0;
				for(int j = 0; j < nct; j++){
					double dist2 = nodes[finfo.nodes[j]].distance2(
						surface->getPoint(surface->getParameters(nodes[finfo.nodes[j]])));
					if(dist2 > max_dist2) max_dist2 = dist2;
				}
				if(max_dist2 > 0.25){
					// TODO -> transform to bspline-surface and fit
				}
			}
		}else if(surface->getType() == SURFACE_BSPLINE){
			((SurfaceBSplinePlanar*)(surface.get()))->fitAdditionally(all_nodes);
		}
	}
	// add corrections
	const double MIN_CORRECTION_DISTANCE = 0.2;
	for(int i = 1; i < faces.countInt(); i++){
		FaceInfo& finfo = faces[i];
		if(!finfo.valid) continue;
		std::shared_ptr<SurfaceCorrected> surf;
		// - vertices
		for(int j = 0; j < finfo.nodes.countInt(); j++){
			const DPoint3d& pt_real = nodes[finfo.nodes[j]];
			DPoint2d pt_param = finfo.surface->getParameters(pt_real);
			DPoint3d pt_surface = finfo.surface->getPoint(pt_param);
			double dist = pt_real.distance(pt_surface);
			if(dist > MIN_CORRECTION_DISTANCE){
				if(!surf) surf = std::make_shared<SurfaceCorrected>(finfo.surface);
				surf->insertCorrectionVector(pt_param, 20*dist, pt_real - pt_surface);
			}
		}
		if(surf){
			finfo.surface = surf;
		}
	}
}

/// Classify lines (linear, arcs, other - approximation)
void MeshGrain::classifyGrainLines(
		DataVector<FaceInfo> & faces, 
		DataVector<LineInfo> & lines, 
		const DataVector<DPoint3d> & nodes)
{
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "======================================");

	const double DIST_TOL2 = 0.15 * param_grain_node_identity_tolerance * param_grain_node_identity_tolerance;

	// -> adjust inner nodes in edges with regard to surfaces
	for(int i = 1; i < lines.countInt(); i++){
		LineInfo& linfo = lines[i];

		// ... gather adjacent non-planar surfaces
		DataVector<SurfacePtr> surfaces;
		DataVector<SurfacePtr> surfaces_planar;
		DataVector<SurfacePtr> surfaces_nonplanar;
		int surf_ct = linfo.incident_faces.countInt();
		for(int j = 0; j < surf_ct; j++){
			FaceInfo& finfo = faces[linfo.incident_faces[j]];
			if(finfo.valid){
				surfaces.add(finfo.surface);
				if(finfo.planar())
					surfaces_planar.add(finfo.surface);
				else
					surfaces_nonplanar.add(finfo.surface);
			}
		}
		linfo.adjacent_nonplanar = surfaces_nonplanar.countInt(); // refresh the counter ...
		if(linfo.adjacent_nonplanar == 0){
			linfo.inner_points.clear(); // all adjacent faces are planar, so edge must be linear
			continue;
		}
		// else ->

		surf_ct = surfaces.countInt();
		assert(surf_ct > 1);
		int ipt_ct = linfo.inner_points.countInt();

		// ! correction
		DataVector<DPoint3d> ipts(ipt_ct);
		double max_dist = 0.0;
		double max_err = 0.0;
		for(int j = 0; j < ipt_ct; j++){
			const DPoint3d& ipt_oryg = linfo.inner_points[j];
			DPoint3d ipt_new = ipt_oryg;
			DataVector<DVector3d> plane_normals(surf_ct+1);
			DataVector<double> plane_ds(surf_ct+1);
			double err = 0.0;
			double dist = 0.0;
			for(int is = 0; is < 5; is++){
				// ... gather surface normals
				for(int k = 0; k < surf_ct; k++){
					DPoint2d pt_param = surfaces[k]->getParameters(ipt_oryg);
					DPoint3d pt_surf = surfaces[k]->getPoint(pt_param);
					DVector3d vn_surf = surfaces[k]->getNormalVector(pt_param).normalized();
					plane_normals.add(vn_surf);
					double d = vn_surf.x * pt_surf.x + vn_surf.y * pt_surf.y + vn_surf.z * pt_surf.z;
					plane_ds.add(d);
				}
				// ... and one more - orthogonal
				int min_k = 1;
				double min_dot = abs(plane_normals[0].scalarProduct(plane_normals[1]));
				for(int k = 2; k < surf_ct; k++){
					double dot = abs(plane_normals[0].scalarProduct(plane_normals[k]));
					if(dot < min_dot){
						min_k = k;
						min_dot = dot;
					}
				}
				if(min_dot < 0.95){
					DVector3d vn_ext = plane_normals[0].crossProduct(plane_normals[min_k]).normalized();
					plane_normals.add(vn_ext);
					double d = vn_ext.x * ipt_oryg.x + vn_ext.y * ipt_oryg.y + vn_ext.z * ipt_oryg.z;
					plane_ds.add(d);
					err = DLeastSquaresFitting::fitPointToPlanes(ipt_new, plane_normals, plane_ds);
					dist = ipt_new.distance(ipt_oryg);
					if(err < 0.01){
						//if(dist > 0.5){
						//	ipt_new = ipt_oryg + (ipt_new-ipt_oryg) * (0.5 / dist);
						//}
						break;
					}
				}else break;
			}
			ipts.add(ipt_new);
			if(err > max_err) max_err = err;
			if(dist > max_dist) max_dist = dist;
		}

		// --- show 
		if(false){
			MeshViewSet* set = new MeshViewSet;
			// ... common edge vertices
			set->addPoint(nodes[linfo.node0], 1, 0);
			set->addPoint(nodes[linfo.node1], 1, 1);
			// ... common edge inner points
			for(int j = 0; j < ipt_ct; j++){
				set->addPoint(linfo.inner_points[j], 2);
				set->addPoint(ipts[j], 3);
				if(j > 0){
					set->addEdge(linfo.inner_points[j-1], linfo.inner_points[j], 2);
					set->addEdge(ipts[j-1], ipts[j], 3);
				}
			}
			// ... surfaces
			for(int j = 0; j < surf_ct; j++){
				FaceInfo& finfo = faces[linfo.incident_faces[j]];
				if(!finfo.valid) continue;
				auto surface = finfo.surface;
				if(surface->getType() == SURFACE_BSPLINE){
					set = ((SurfaceBSplinePlanar*)(surface.get()))->getViewSet(set);
				}else{
					DRect br;
					for(int k = 0; k < finfo.nodes.countInt(); k++){
						br.addPoint(surface->getParameters(nodes[finfo.nodes[k]]));
					}
					double RES = 20;
					double dx = br.getDX() / RES;
					double dy = br.getDY() / RES;
					for(int ix = 0; ix <= RES; ix++){
						for(int iy = 1; iy <= RES; iy++){
							DPoint3d ptA = surface->getPoint(DPoint2d(br.x0 + ix * dx, br.y0 + iy * dy));
							DPoint3d ptB = surface->getPoint(DPoint2d(br.x0 + ix * dx, br.y0 + (iy-1) * dy));
							DPoint3d ptC = surface->getPoint(DPoint2d(br.x0 + iy * dx, br.y0 + ix * dy));
							DPoint3d ptD = surface->getPoint(DPoint2d(br.x0 + (iy-1) * dx, br.y0 + ix * dy));
							set->addEdge(ptA, ptB, 0);
							set->addEdge(ptC, ptD, 0);
						}
					}
				}
			}
			// ... show
			ostringstream oss;
			oss << "common edge, line #" << i 
				<< ", with " << surf_ct << " valid surfaces"
				<< ", max_err = " << max_err
				<< ", max dist = " << max_dist;
			SHOW_MESH(oss.str(), set);
		}
		// ---

		const DPoint3d& node_pt0 = nodes[linfo.node0];
		const DPoint3d& node_pt1 = nodes[linfo.node1];
		const DPoint3d middle(node_pt0, node_pt1, 0.5);
		const DVector3d vn = (node_pt1 - node_pt0).normalized();
		double d = vn.scalarProduct(middle - DPoint3d::zero); // plane-equation
		double node_dist = vn.scalarProduct(node_pt1 - DPoint3d::zero) - d;
		linfo.inner_points.clear();
		for(int j = 0; j < ipts.countInt(); j++){
			if( ipts[j].distance2(node_pt0) < DIST_TOL2 || // if  close to one vertex
				ipts[j].distance2(node_pt1) < DIST_TOL2 || // ... or close to the other one
				abs(vn.scalarProduct(ipts[j] - DPoint3d::zero ) - d) > node_dist ) // ... or outside
			{
				// remove inner points too close or outside the line-vertices
				// --- or actually, do nothing...
				//ipts.removeOrderedAt(j);
				//j--;
			}else{
				linfo.inner_points.add(ipts[j]);
			}
		}

	}

	int bspline_edges = 0;

	for(int i = 1; i < faces.countInt(); i++){
		FaceInfo& finfo = faces[i];
		if(!finfo.valid) continue;
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "=== Classify lines for face -> " << i);
		auto surface = finfo.surface;
		// check orientation
		finfo.local_node_params.clear();
		finfo.local_node_params.prepare(finfo.nodes.countInt());
		int pct = finfo.nodes.countInt();
		for(int j = 0; j < pct; j++){
			finfo.local_node_params.add(surface->getParameters(nodes[finfo.nodes[j]]));
		}
		if(!DPoint2d::properOrientation(finfo.local_node_params)){ // reverse
			// .. nodes
			std::reverse( finfo.nodes.begin(), finfo.nodes.end() );
			std::reverse( finfo.local_node_params.begin(), finfo.local_node_params.end() );
			//finfo.nodes.reverse();
			//finfo.local_node_params.reverse();
			// .. lines
			int last_id = finfo.lines.removeLast();
			std::reverse( finfo.lines.begin(), finfo.lines.end() );
			//finfo.lines.reverse();
			finfo.lines.add(last_id);
		}
		// classify edges
		finfo.linear.clear();
		finfo.linear.addItems(finfo.lines.countInt(), false);
		finfo.curves.clear();
		finfo.curves.addItems(finfo.lines.countInt());
		DataVector<double> max_distances(finfo.lines.countInt(), 0.0);
		for(int j = 0; j < finfo.lines.countInt(); j++){
			LineInfo& linfo = lines[finfo.lines[j]];
			const DPoint2d param_point0 = surface->getParameters(nodes[linfo.node0]);
			const DPoint2d param_point1 = surface->getParameters(nodes[linfo.node1]);
			// -> if to few inner nodes or all adjacent surfaces are planar, then nothing to do ...
			if(linfo.inner_points.countInt() < 3 || linfo.adjacent_nonplanar == 0){
				finfo.linear[j] = true;
				finfo.curves[j] = std::make_shared<Curve2dSegment>(param_point0, param_point1);
				continue;
			}
			// project points
			DataVector<DPoint2d> param_points(linfo.inner_points.countInt());
			Curve2dSegment line(param_point0, param_point1);
			double line_max_dist2 = 0.0;
			double last_t = 0.0;
			for(int k = 0; k < linfo.inner_points.countInt(); k++){
				DPoint2d param = surface->getParameters(linfo.inner_points[k]);
				param_points.add(param);
				DPoint2d line_param = line.getPoint(last_t = line.getParameter(param, last_t));
				double dist2 = surface->getPoint(param).distance2(surface->getPoint(line_param));
				if(dist2 > line_max_dist2) line_max_dist2 = dist2;
			}
			double max_dist = sqrt(line_max_dist2);
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, " * max point-to-line distance -> " << max_dist);
			//if(true){
			//	//--------------------
			//	// show line
			//	MeshViewSet *set = new MeshViewSet;
			//	ostringstream oss;
			//	// -- nodes
			//	set->addPoint(nodes[linfo.node0], 2, linfo.node0);
			//	set->addPoint(nodes[linfo.node1], 2, linfo.node1);
			//	// -- points 3d and projected to surface
			//	for(int k = 0; k < linfo.inner_points.countInt(); k++){
			//		set->addPoint(surface->getPoint(param_points[k]), 2);
			//	}
			//	// -- fit line
			//	set->addEdge(
			//		surface->getPoint(param_point0), 
			//		surface->getPoint(param_point1), 0);
			//	oss << "fit line for surface, max_dist=" << max_dist ;			
			//	SHOW_MESH(oss.str(), set);
			//}

			if(max_dist < 0.1*param_grain_tolerance){
				finfo.linear[j] = true;
				finfo.curves[j] = std::make_shared<Curve2dSegment>(param_point0, param_point1);
			}else{
				// fit bspline
				auto bspline = std::make_shared<Curve2dBSpline>();
				double ml = param_grain_tolerance;
				DMetric2d dm(ControlDataMatrix2d(ml, 0.0, ml), surface, param_point0);
				if(finfo.planar()){
					max_dist = bspline->fitToPoints(param_points, 
						param_point0, param_point1, dm);
				}else{
					bspline->fitThroughPoints(param_points, 
						param_point0, param_point1, dm);
				}
				finfo.curves[j] = bspline;
				bspline_edges++;
			}
			max_distances[j] = max_dist * param_grain_tolerance;
		}
		if(true){
			//--------------------
			// show line
			ostringstream oss;
			for(int j = 0; j < finfo.nodes.countInt(); j++){
				oss << "[" << finfo.nodes[j] << "-" << finfo.nodes[(j+1)%finfo.nodes.countInt()]
					<< " " << (finfo.linear[j] ? "linear" : "bspline")
					<< ", max_dist=" << max_distances[j] << "] ";
			}
			SHOW_MESH(oss.str(), getFaceViewSet(finfo, nodes, lines));
		}
	}
	if(bspline_edges > 0)
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "bspline edges, total= " << bspline_edges);
}

/// Create surface meshes from points, for each face separately
void MeshGrain::createSurfaceMeshes(DataVector<FaceInfo> & faces, 
	DataVector<LineInfo> & lines, 
	const DataVector<DPoint3d> & nodes)
{
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "======================================");
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "=== Create surface meshes, total = " << faces.countInt());

	for(int i = 1; i < faces.countInt(); i++){
		FaceInfo & face = faces[i];
		if(!face.valid){
			if((param_grain_view & 4) > 0){
				//--------------------
				// show surface
				ostringstream oss;
				oss << "(invalid) grain face #" << i;
				SHOW_MESH(oss.str(), 
					getFaceViewSet(face, nodes, lines));
				//--------------------
			}
			continue;
		}

		int pct = face.nodes.countInt() + face.inner_points.countInt();
		for(int j = 0; j < face.lines.countInt(); j++){
			if(lines[face.lines[j]].valid)
				pct += lines[face.lines[j]].inner_points.countInt();
		}

		DataVector<MeshPoint3d*> face_mesh_points(pct);
		DBox fbox;
		// gather all points
		// -- face nodes
		for(int j = 0; j < face.nodes.countInt(); j++){
			MeshPoint3d* point = new MeshPoint3d(nodes[face.nodes[j]]);
			point->setIntTag(TagExtended::TAG_GRAIN_NODE_ID, j);
			point->setIntTag(TagExtended::TAG_GRAIN_NODE_TYPE, 1); // 1 - vertex
 			fbox.addPoint(point->getCoordinates());
			face_mesh_points.add(point);
		}
		// -- face edges
		for(int j = 0; j < face.lines.countInt(); j++){
			LineInfo& line = lines[face.lines[j]];
			if(!line.valid) continue;
			for(int k = 0; k < line.inner_points.countInt(); k++){
				MeshPoint3d* point = new MeshPoint3d(line.inner_points[k]);
				point->setIntTag(TagExtended::TAG_GRAIN_NODE_ID, j);
				point->setIntTag(TagExtended::TAG_GRAIN_NODE_ID2, k);
				point->setIntTag(TagExtended::TAG_GRAIN_NODE_TYPE, 2); // 2 - inner node for line
 				fbox.addPoint(point->getCoordinates());
				face_mesh_points.add(point);
			}
		}
		// -- face inner nodes
		for(int j = 0; j < face.inner_points.countInt(); j++){
			MeshPoint3d* point = new MeshPoint3d(face.inner_points[j]);
			point->setIntTag(TagExtended::TAG_GRAIN_NODE_ID, j);
			point->setIntTag(TagExtended::TAG_GRAIN_NODE_TYPE, 4); // 4 - inner node for face
 			fbox.addPoint(point->getCoordinates());
			face_mesh_points.add(point);
		}
		assert(pct == face_mesh_points.countInt());

		auto csi = std::make_shared<ControlSpace3dIdentity>(2.0*param_grain_tolerance);
		fbox.inflate(2.0);
		MeshContainer3d* vmesh = MeshGenerator3dDelaunayBoundary::createInitialMesh(pct, fbox, csi);

		DataVector<int> sequence(pct);
		for(int i = 0; i < pct; i++) sequence.add(i);
		std::random_shuffle( sequence.begin(), sequence.end() );

		// * insert new points into the mesh
		Metric3dContext mc(csi);
		for(int i = 0; i < pct; i++){
			// Select one of points in random
			mc.countMetricAtPoint(face_mesh_points[i]->getCoordinates());
			bool success = MeshGenerator3d::addPointToTriangulation(mc, vmesh, face_mesh_points[i]);
			assert(success);
		}
	}
}

/// Classify surfaces (planar, quadratic, other - approximation)
void MeshGrain::classifyGrainSurfaces(DataVector<FaceInfo> & faces, 
									  DataVector<LineInfo> & lines, 
									  const DataVector<DPoint3d> & nodes)
{
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "======================================");
	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "=== Classify surfaces, total = " << faces.countInt());

	int planar_count    = 0;
	int cylinder_count  = 0;
	int bspline_count   = 0;
	int multi_count		= 0;
	int unknown_count   = 0;
	int invalid_count   = 0;


	for(int i = 1; i < faces.countInt(); i++){
		FaceInfo & face = faces[i];
		if(!face.valid){
			if((param_grain_view & 4) > 0){
				//--------------------
				// show surface
				ostringstream oss;
				oss << "(invalid) grain face #" << i;
				SHOW_MESH(oss.str(), 
					getFaceViewSet(face, nodes, lines));
				//--------------------
			}
			++invalid_count;
			continue;
		}

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Plane => " << i);
		// gather all points
		DataVector<DPoint3d> all_points(face.nodes.countInt());
		// -- face nodes
		for(int j = 0; j < face.nodes.countInt(); j++)
			all_points.add(nodes[face.nodes[j]]);
//		if(face.inner_points.countInt() < 10){
			// -- face edges
			for(int j = 0; j < face.lines.countInt(); j++){
				LineInfo& line = lines[face.lines[j]];
				if(!line.valid) continue;
				for(int k = 0; k < line.inner_points.countInt(); k++)
					all_points.add(line.inner_points[k]);
			}
//		}
		// -- face inner nodes
		for(int j = 0; j < face.inner_points.countInt(); j++)
			all_points.add(face.inner_points[j]);

		int all_pct = all_points.countInt();
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "Plane fitting, number of points => " << all_pct);
		// Least Squares Fitting
		//	* hyperplanar fitting of points using orthogonal regression
		// -- A -> average of the sample points
		DPlane plane;
		double plane_max_dist = DLeastSquaresFitting::fitHyperplaneOrthogonal(
			all_points, plane);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "max point-to-plane distance = " <<  plane_max_dist);
		bool found_surface = false;

		// check plane
		if(plane_max_dist < 0.5*param_grain_tolerance || all_pct < SurfaceBSplinePlanar::FIT_PROBE_MIN){
			face.surface = std::make_shared<SurfacePlane>(plane);
			++planar_count;
			found_surface = true;
		}

		DCylinder cylinder;
		double cyl_max_dist;

		if(!found_surface){ // -> fit cylinder
			cyl_max_dist = DLeastSquaresFitting::fitCylinder(all_points, plane, cylinder);
			if(cyl_max_dist < 0.5*param_grain_tolerance){
				face.surface = std::make_shared<SurfaceCylinder>(cylinder);
				++cylinder_count;
				// -> mark as an additional nonplanar surface for all its edges
				for(int j = 0; j < face.lines.countInt(); j++)
					lines[face.lines[j]].adjacent_nonplanar++;
				found_surface = true;
			}
		}

		if(!found_surface){	// -> fit bspline-surface
			// planar
			SurfacePlane splane(plane);
			DRect brect;
			for(int j = 0; j < all_pct; j++)
				brect.addPoint(splane.getParameters(all_points[j]));
			brect.inflate(0.05);
			auto bsurf2 = std::make_shared<SurfaceBSplinePlanar>(splane, brect);
			double bsurf_max_dist2 = bsurf2->fitToPoints(all_points, 2);
			auto bsurf3 = std::make_shared<SurfaceBSplinePlanar>(splane, brect);
			double bsurf_max_dist3 = bsurf3->fitToPoints(all_points, 3);
			double bsurf_max_dist_min = std::min(bsurf_max_dist2, bsurf_max_dist3);
			// cylindrical
			SurfaceCylinder cylinder(cylinder);
			DRect brect_cyl;
			for(int j = 0; j < all_pct; j++)
				brect_cyl.addPoint(cylinder.getParameters(all_points[j]));
			brect_cyl.inflate(0.05);
			auto bcsurf2 = std::make_shared<SurfaceBSplineCylindrical>(cylinder, brect_cyl);
			double bcsurf_max_dist2 = bcsurf2->fitToPoints(all_points, 2);
			auto bcsurf3 = std::make_shared<SurfaceBSplineCylindrical>(cylinder, brect_cyl);
			double bcsurf_max_dist3 = bcsurf3->fitToPoints(all_points, 3);
			double bcsurf_max_dist_min = std::min(bcsurf_max_dist2, bcsurf_max_dist3);
			// check
			if(plane_max_dist <= bsurf_max_dist_min && plane_max_dist < param_grain_tolerance){
				// stay with planar fit
				face.surface = std::make_shared<SurfacePlane>(plane);
				++planar_count;
				found_surface = true;
			}else if(cyl_max_dist <= bcsurf_max_dist_min && cyl_max_dist < param_grain_tolerance){
				// stay with cylinder fit
				face.surface = std::make_shared<SurfaceCylinder>(cylinder);
				++cylinder_count;
				// -> mark as an additional nonplanar surface for all its edges
				for(int j = 0; j < face.lines.countInt(); j++)
					lines[face.lines[j]].adjacent_nonplanar++;
				found_surface = true;
			}else{
				// -> mark as an additional nonplanar surface for all its edges
				for(int j = 0; j < face.lines.countInt(); j++)
					lines[face.lines[j]].adjacent_nonplanar++;
				// find best approximation
				double max_dist = 0.0;
				if(bsurf_max_dist_min < param_grain_tolerance){
					// use better of the two bsurfaces
					if(bsurf_max_dist2 <= bsurf_max_dist3){
						face.surface = bsurf2;
						max_dist = bsurf_max_dist2;
					}else{
						face.surface = bsurf3;
						max_dist = bsurf_max_dist3;
					}
					LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
						"max point-to-bspline (direct-fit) distance = " <<  max_dist);
					++bspline_count;
				}else if(bcsurf_max_dist_min < param_grain_tolerance){
					// use better of the two bcsurfaces
					if(bcsurf_max_dist2 <= bcsurf_max_dist3){
						face.surface = bcsurf2;
						max_dist = bcsurf_max_dist2;
					}else{
						face.surface = bcsurf3;
						max_dist = bcsurf_max_dist3;
					}
					LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
						"max point-to-bcspline (direct-fit) distance = " <<  max_dist);
					++bspline_count;
				}else{
					// try higher level

					double max_dist = bsurf2->fitToPoints(all_points, 4);
					SurfacePtr bsurf_best;
					if(max_dist > param_grain_tolerance){
						max_dist = bsurf2->fitToPoints(all_points, 5);
						double max_dist_m = 0.0;
						auto msurface = SurfaceBSplinePlanar::multiFitToPointsRegular(
							bsurf2, all_points, max_dist_m);
						if(max_dist_m > max_dist){
							bsurf_best = bsurf2;
						}else{
							max_dist = max_dist_m;
							bsurf_best = msurface;
						}
					}
					// ... and ?
					if(max_dist < 2*param_grain_tolerance){
						face.surface = bsurf_best;
						if(bsurf_best->getType() == SURFACE_BSPLINE){
							LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
								"max point-to-bspline (direct-fit) distance = " <<  max_dist);
							++bspline_count;
						}else{
							LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
								"max point-to-multisurface (direct-fit) distance = " <<  max_dist);
							++multi_count;
						}
					}else{
//						face.surface.set(SurfaceBSplinePlanar::multiFitToPoints(
//							bsurf2, all_points, 2, 3, 
//							param_grain_tolerance, max_dist));
//						if(max_dist < param_grain_tolerance){
//							LOG4CPLUS_INFO(MeshLog::logger_mesh, "max point-to-multisurface (direct-fit) distance = " <<  max_dist);
//							++multi_count;
//						}else{
							++unknown_count;
//						}
					}
				}
			}
		}

		if((param_grain_view & (face.planar() ? 2 : 1)) > 0){
//		if(!face.planar()){
			//--------------------
			// show surface
			MeshViewSet *set = getFaceViewSet(face, nodes, lines);

			// -- fitting plane vectors
			set->addPoint(plane.p0, 4);
			set->addEdge(plane.p0, plane.p0 + plane.e0, 3);
			set->addEdge(plane.p0, plane.p0 + plane.e1, 3);
			set->addEdge(plane.p0, plane.p0 + plane.e0.crossProduct(plane.e1) * 0.2, 4);

			ostringstream oss;
			if(!face.surface) oss << "complex";
			else if(face.surface->getType() == SURFACE_PLANE) oss << "planar";	
			else if(face.surface->getType() == SURFACE_CYLINDER) oss << "cylinder";	
			else if(face.surface->getType() == SURFACE_BSPLINE) oss << "bspline";
			else if(face.surface->getType() == SURFACE_CORRECTED) oss << "corrected";
			else if(face.surface->getType() == SURFACE_MULTI) oss << "multi-bspline";
			else oss << "unknown";

			oss << " grain face #" << i;
			SHOW_MESH(oss.str(), set);
			//--------------------
		}
	}
	if(planar_count > 0)   LOG4CPLUS_DEBUG(MeshLog::logger_console, "Planar faces: " << planar_count);
	if(cylinder_count > 0) LOG4CPLUS_DEBUG(MeshLog::logger_console, "Cylinder faces: " << cylinder_count);
	if(bspline_count > 0)  LOG4CPLUS_DEBUG(MeshLog::logger_console, "BSPline faces: " << bspline_count);
	if(unknown_count > 0)  LOG4CPLUS_DEBUG(MeshLog::logger_console, "Complex faces: " << unknown_count);
}

bool MeshGrain::storeGrainXml(const string& dir, const DataVector<DPoint3d> & nodes,
		const DataVector<int> & trans_nodes, const DataVector<LineInfo> & lines,
		const DataVector<FaceInfo> & faces, const DataVector<BlockInfo> & blocks)
{
	string fname = dir+"grain.xml";
	ofstream ofs(fname.c_str());
	if(!ofs){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Error opening file for write: " << fname);
		return false;
	}

	const char* header = 
		"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		"<meshdoc xmlns=\"http://www.icsr.agh.edu.pl\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n"
		"\t<header>\n"
		"\t\t<creator>Tomasz Jurczyk</creator>\n"
		"\t\t<version>1.1</version>\n"
		"\t\t<description>auto-generated grain brep-model</description>\n"
		"\t</header>\n\n";

	ofs << header;

	//const char* hconst = 
	//	"\t<constants>\n"
	//	"\t\t<const>\n"
	//	"\t\t\t<name>CR</name> <value>2.0</value>\n"
	//	"\t\t</const>\n"
	//	"\t</constants>\n\n";

	//ofs << hconst;

	ofs << "\t<model>" << endl;
	// write nodes
	ofs << "\t\t<vertices>" << endl;
	int vct = 0;
	for(int i = 1; i < nodes.countInt(); i++){
		if(i >= trans_nodes.countInt() || trans_nodes[i] == i){
			++vct;
			ofs << "\t\t\t<vertex vid=\"" << i << "\"> " 
				<< nodes[i]	<< " </vertex>" << endl;
		}
	}
	ofs << "\t\t</vertices>" << endl << endl;
	LOG4CPLUS_INFO(MeshLog::logger_console, "Stored grain nodes: " << vct);

	int sct = 0;
	DataVector<bool> stored_surfaces(faces.countInt(), false);
	// write surfaces
	for(int i = 1; i < faces.countInt(); i++){
		const FaceInfo& finfo = faces[i];
		if(!finfo.valid || !finfo.surface) continue;
		bool linear_edges = true;
		for(int j = 0; j < finfo.linear.countInt(); j++)
			if(!finfo.linear[j]){
				linear_edges = false;
				break;
			}
		if(finfo.planar() && linear_edges){
			// for planes with linear edges, not necessary really ...
		}else{
			if(++sct == 1) ofs << "\t\t<surfaces>\n";
			// ... surface info here ...
			stored_surfaces[i] = true;
			ofs << "\t\t\t<surface sid=\"" << i << "\">\n";
			finfo.surface->storeXML(ofs, "\t\t\t\t");
			ofs << "\t\t\t</surface>\n";
		}
	}
	if(sct > 0){
		ofs << "\t\t</surfaces>\n\n";
		LOG4CPLUS_INFO(MeshLog::logger_console, "Stored grain surfaces: " << sct);
	}

	int cct = 0;
	if(sct > 0){
		// write curves
		ofs << "\t\t<curves>" << endl;
		for(int i = 1; i < faces.countInt(); i++){
			const FaceInfo& finfo = faces[i];
			if(!finfo.valid) continue;
			for(int j = 0; j < finfo.curves.countInt(); j++){
				if(!finfo.linear[j] && finfo.curves[j]){
					auto curve = finfo.curves[j];
					ofs << "\t\t\t<curve cid=\"" << cct << "\">\n";
					curve->setIntTag(TagExtended::TAG_ID, cct);
					curve->storeXML(ofs, "\t\t\t\t");
					ofs << "\t\t\t</curve>\n";
					cct++;
				}
			}
		}
		ofs << "\t\t</curves>" << endl << endl;
		LOG4CPLUS_INFO(MeshLog::logger_console, "Stored grain curves: " << cct);
	}
	if(cct > 0){
		// write edges
		ofs << "\t\t<edges>" << endl;
		for(int i = 1; i < faces.countInt(); i++){
			const FaceInfo& finfo = faces[i];
			if(!finfo.valid) continue;
			for(int j = 0; j < finfo.curves.countInt(); j++){
				if(!finfo.linear[j] && finfo.curves[j]){
					ofs << "\t\t\t<edge sid=\"" << i << "\""
						<< " vid0=\"" << lines[finfo.lines[j]].node0 << "\""
						<< " vid1=\"" << lines[finfo.lines[j]].node1 << "\">"
						<< endl;
					ofs << "\t\t\t\t<cid>" << finfo.curves[j]->getIntTag(TagExtended::TAG_ID)
						<< "</cid> <t0>0</t0> <t1>1</t1>\n";
					ofs << "\t\t\t</edge>\n";
				}
			}
		}
		ofs << "\t\t</edges>" << endl << endl;
	}

	// write faces
	ofs << "\t\t<faces>" << endl;
	int fct = 0;
	for(int i = 1; i < faces.countInt(); i++){
		const FaceInfo& finfo = faces[i];
		if(!finfo.valid) continue;
		++fct;
		if(stored_surfaces[i]){
			ofs << "\t\t\t<face fid=\"" << i << "\" sid=\"" 
				<< i << "\">"; 
		}else{
			ofs << "\t\t\t<face fid=\"" << i << "\">\n\t\t\t\t"; 
		}
		for(int j = 0; j < finfo.nodes.countInt(); j++){
			if(stored_surfaces[i]){
				ofs << "\n\t\t\t\t<vertex vid=\"" << finfo.nodes[j] << "\"> "
					<< finfo.local_node_params[j] << " </vertex>";
			}else{
				ofs << "<vertex vid=\"" << finfo.nodes[j] << "\" />\t";
			}
		}
		ofs << "\n\t\t\t</face>" << endl << endl;
	}
	ofs << "\t\t</faces>" << endl << endl;
	LOG4CPLUS_INFO(MeshLog::logger_console, "Stored grain faces: " << fct);

	// write blocks
	ofs << "\t\t<blocks>" << endl;
	int bct = 0;
	for(int i = 1; i < blocks.countInt(); i++){
		const BlockInfo& binfo = blocks[i];
		if(binfo.faces.empty()) continue;
		++bct;
		ofs << "\t\t\t<block bid=\"" << i << "\">\n";
		for(int j = 0; j < binfo.faces.countInt(); j++){
			if(faces[binfo.faces[j]].valid){
				ofs << "\t\t\t\t<face fid=\"" <<  binfo.faces[j] << "\" ";
				if(binfo.inverted[j]) ofs << "inverted=\"true\" ";
				ofs << "/>\n";
			}
		}
		ofs << "\t\t\t</block>" << endl << endl;
	}
	ofs << "\t\t</blocks>" << endl;
	LOG4CPLUS_INFO(MeshLog::logger_console, "Stored grain blocks: " << bct);

	ofs << "\t</model>" << endl;
	ofs << "\t<sizing>" << endl;
	ofs << "\t\t<params>" << endl;
	ofs << "\t\t\t<param>" << endl;
	ofs << "\t\t\t\t<name>ACS_MIN_LENGTH</name>" << endl;
	ofs << "\t\t\t\t<value>1.0</value>" << endl;
	ofs << "\t\t\t</param>" << endl;
	ofs << "\t\t</params>" << endl;
	ofs << "\t</sizing>" << endl;
	ofs << "</meshdoc>" << endl;
	return true;
}

void MeshGrain::classifyGrainBlocks(DataVector<BlockInfo> & blocks,
		const DataVector<FaceInfo> & faces, 
		const DataVector<DPoint3d> & nodes)
{
	int valid_blocks = 0;
	for(int i = 1; i < blocks.countInt(); i++){
		BlockInfo& binfo = blocks[i];
		if(binfo.faces.empty()) continue;
		++valid_blocks;
		binfo.inverted.addItems(binfo.faces.countInt(), false);
		// count middle point for block
		DataVector<bool> used(nodes.countInt()+1, false);
		DPoint3d middle;
		int middle_ct = 0;
		for(int j = 0; j < binfo.faces.countInt(); j++){
			const FaceInfo& finfo = faces[binfo.faces[j]];
			for(int k = 0; k < finfo.nodes.countInt(); k++){
				int nid = finfo.nodes[k];
				if(!used[nid]){
					middle.add(nodes[nid]);
					used[nid] = true;
					++middle_ct;
				}
			}
		}
		middle /= middle_ct;
		// check orientation
		for(int j = 0; j < binfo.faces.countInt(); j++){
			const FaceInfo& finfo = faces[binfo.faces[j]];
			const DVector3d face_vn = finfo.surface->getNormalVector(DPoint2d(0.5, 0.5));
			const DVector3d face_vt = middle - finfo.surface->getPoint(DPoint2d(0.5, 0.5));
			if(face_vn.scalarProduct(face_vt) > 0.0)
				binfo.inverted[j] = true;
		}
	}
	LOG4CPLUS_DEBUG(MeshLog::logger_console, "valid grain blocks: " << valid_blocks);
}

/// Check topology of model entities
bool MeshGrain::checkTopology(const DataVector<DPoint3d> & nodes, 
	DataVector<LineInfo> & lines, 
	DataVector<FaceInfo> & faces, 
	DataVector<BlockInfo> & blocks)
{

	// check face-vertex incidency
	DataVector<int> rank(nodes.countInt(), 0);
	for(int i = 1; i < faces.countInt(); i++){
		for(int j = 0; j < faces[i].nodes.countInt(); j++)
			rank[faces[i].nodes[j]]++;
	}
	for(int i = 0; i < nodes.countInt(); i++){
		if(rank[i] > 0 && rank[i] < 3){
			LOG4CPLUS_WARN(MeshLog::logger_console, "Sparse vertex id: " << i);
		}
	}

	// check face-line incidency for block
	for(int i = 1; i < blocks.countInt(); i++){
		DataVector< DataVector<int> > line_faces(lines.countInt(), DataVector<int>());

		for(int j = 0; j < blocks[i].faces.countInt(); j++){
			int fid = blocks[i].faces[j];
			const FaceInfo& face = faces[fid];
			assert(face.lines.countInt() == face.nodes.countInt());
			for(int k = 0; k < face.lines.countInt(); k++){
				line_faces[face.lines[k]].add(fid);
			}
		}
		int count = 0;
		for(int j = 0; j < lines.countInt(); j++){
			int trank = line_faces[j].countInt();
			if(trank != 0 && trank != 2){
				if(count++ == 0){
					LOG4CPLUS_WARN(MeshLog::logger_console, "===> face-continuity, grain id" << i);
				}
				LOG4CPLUS_WARN(MeshLog::logger_console, 
					"--> line (" << j << ": " << 
					lines[j].node0 << "-" << lines[j].node1 <<
					"), rank = " << trank);
				for(int k = 0; k < trank; k++){
					LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "----> face id = " << line_faces[j][k]);
				}
			}
		}
	}

	// find border faces, mark "fixed coordinates" of theirs vertices

	return true;
}

/// Optimize placement of nodes according to face-surfaces
void MeshGrain::optimizeNodePlacement(DataVector<FaceInfo> & faces, 
	DataVector<DPoint3d> & nodes)
{
	// gather list of faces for each node
	DataVector< DataVector<int> > node_faces(nodes.countInt(), DataVector<int>());
	for(int i = 0; i < faces.countInt(); i++){
		if(!faces[i].valid) continue;
		for(int j = 0; j < faces[i].nodes.countInt(); j++)
			node_faces[faces[i].nodes[j]].add(i);
	}
	// calculate new coordinates
	for(int i = 0; i < nodes.countInt(); i++){
		if(node_faces[i].empty()) continue;
		// create matrix
		DMatrix3d A;
		DVector3d b;
		int nfct = node_faces[i].countInt();
		DataVector<SurfacePtr> surfaces(nfct);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
			"=============> Calculating better coordinates for node #" << i 
			<< " with " << nfct << " incident faces");
		MeshViewSet* set = nullptr;
		for(int j = 0; j < nfct; j++){
			FaceInfo& f = faces[node_faces[i][j]];
			//set = getFaceViewSet(f, nodes, lines, set);
			if(!f.planar())
				surfaces.add(f.surface);
			//  -> project point3d
			DPoint2d param_point = f.surface->getParameters(nodes[i]);
			//	-> for non-planar use tangent plane as an approximation
			DVector3d vn = f.surface->getNormalVector(param_point);
			A.m[0][0] += vn.x * vn.x;
			A.m[1][1] += vn.y * vn.y;
			A.m[2][2] += vn.z * vn.z;
			A.m[0][1] += vn.x * vn.y;
			A.m[0][2] += vn.x * vn.z;
			A.m[1][2] += vn.y * vn.z;
			DPoint3d plane_pt = f.surface->getPoint(param_point);
			double d = vn.x * plane_pt.x + vn.y * plane_pt.y + vn.z * plane_pt.z;
			b.x += d * vn.x;
			b.y += d * vn.y;
			b.z += d * vn.z;
			// debug info
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
				" --> face #" << node_faces[i][j] << ", oryginal distance = "
				<< (vn.scalarProduct(nodes[i] - DPoint3d::zero) - d));
		}
		A.m[1][0] = A.m[0][1];
		A.m[2][0] = A.m[0][2];
		A.m[2][1] = A.m[1][2];

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, fixed << b.x << " | " << A.m[0][0]);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, fixed << b.y << " | " << A.m[1][0] << "\t" << A.m[1][1]);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, fixed << b.z << " | " << A.m[2][0] << "\t" << A.m[2][1] << "\t" << A.m[2][2]);
		// solve
		const DPoint3d& old_pt = nodes[i];
		DPoint3d new_pt = countOptimizedPoint(old_pt, A, b);
		// stat
		double dist = new_pt.distance(old_pt);
		if(dist > 0.5){
			new_pt = old_pt + (new_pt-old_pt) * (0.5 / dist);
		}
		// -> fix for non-planar surfaces
		if(!surfaces.empty()){
			double w = 1.0 / surfaces.countInt();
			DPoint3d new_spt = DPoint3d::zero;
			for(int k = 0; k < surfaces.countInt(); k++){
				new_spt.add(surfaces[k]->getPoint(surfaces[k]->getParameters(new_pt)), w);
			}
			new_pt = new_spt;
		}
		// adjust with param_bounding_box
		if(param_bounding_box.valid){
			if(new_pt.x - param_bounding_box.x0 < 0.5*param_grain_node_identity_tolerance)
				new_pt.x = param_bounding_box.x0;
			if(param_bounding_box.x1 - new_pt.x < 0.5*param_grain_node_identity_tolerance)
				new_pt.x = param_bounding_box.x1;

			if(new_pt.y - param_bounding_box.y0 < 0.5*param_grain_node_identity_tolerance)
				new_pt.y = param_bounding_box.y0;
			if(param_bounding_box.y1 - new_pt.y < 0.5*param_grain_node_identity_tolerance)
				new_pt.y = param_bounding_box.y1;

			if(new_pt.z - param_bounding_box.z0 < 0.5*param_grain_node_identity_tolerance)
				new_pt.z = param_bounding_box.z0;
			if(param_bounding_box.z1 - new_pt.z < 0.5*param_grain_node_identity_tolerance)
				new_pt.z = param_bounding_box.z1;
		}
		// stat
		dist = new_pt.distance(old_pt);
		// show
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "optimal point relocation distance (" << i << ") => " << dist);
		if(set){
			if(dist > 0.4){
				set->addPoint(new_pt, 5, i);
				ostringstream oss;
				oss << "node #" << i << " placement optimization, dist=" << dist;
				SHOW_MESH(oss.str(), set);
			}else
				delete set;
		}
		// apply
		nodes[i] = new_pt;
	}
}

/// Check topology of model entities
MeshViewSet* MeshGrain::getFaceViewSet(const FaceInfo& face,
	const DataVector<DPoint3d> & nodes, 
	const DataVector<LineInfo> & lines,
	MeshViewSet* set)
{
	// show surface
	if(!set) set = new MeshViewSet;
	// -- face nodes
	for(int j = 0; j < face.nodes.countInt(); j++){
		set->addPoint(nodes[face.nodes[j]], 0, face.nodes[j]);
	}
	// -- face edges
	for(int j = 0; j < face.lines.countInt(); j++){
		const LineInfo& line = lines[face.lines[j]];
		if(!line.valid) continue;
		bool linear_j = (j < face.linear.countInt()) && face.linear[j];
		bool cvalid_j = (j < face.curves.countInt()) && face.curves[j];
		if(linear_j || !face.surface || !cvalid_j){
			set->addEdge(nodes[line.node0], nodes[line.node1], 0);
		}else{
			DataVector<DPoint3d> polyline;
			face.surface->getPolyLine(polyline, face.curves[j], 0.0, 1.0);
			for(int k = 0; k < polyline.countInt()-1; k++){
				set->addEdge(polyline[k], polyline[k+1], 2);
			}
		}
		set->addPoint(nodes[line.node0], 2, line.node0);
		set->addPoint(nodes[line.node1], 2, line.node1);
		for(int k = 0; k < line.inner_points.countInt(); k++)
			set->addPoint(line.inner_points[k], 1);
	}
	// -- face inner nodes
	for(int j = 0; j < face.inner_points.countInt(); j++)
		set->addPoint(face.inner_points[j], face.planar() ? 3 : 4);
	// that's all
	return set;
}

/// Solve optimization problem for node placement
DPoint3d MeshGrain::countOptimizedPoint(const DPoint3d& old_pt, const DMatrix3d& A, const DVector3d& b)
{
	bool deg_x = abs(b.x) < SMALL_NUMBER && abs(A.m[0][0]) < SMALL_NUMBER && 
		abs(A.m[0][1]) < SMALL_NUMBER && abs(A.m[0][2]) < SMALL_NUMBER;
	bool deg_y = abs(b.y) < SMALL_NUMBER && abs(A.m[1][0]) < SMALL_NUMBER && 
		abs(A.m[1][1]) < SMALL_NUMBER && abs(A.m[1][2]) < SMALL_NUMBER;
	bool deg_z = abs(b.z) < SMALL_NUMBER && abs(A.m[2][0]) < SMALL_NUMBER && 
		abs(A.m[2][1]) < SMALL_NUMBER && abs(A.m[2][2]) < SMALL_NUMBER;

	if(deg_x && deg_y && deg_z){
		return old_pt;
	}else if(deg_x && deg_y){
		// count z
		return DPoint3d(old_pt.x, old_pt.y, b.z / A.m[2][2]);
	}else if(deg_x && deg_z){
		// count y
		return DPoint3d(old_pt.x, b.y / A.m[1][1], old_pt.z);
	}else if(deg_y && deg_z){
		// count x
		return DPoint3d(b.x / A.m[0][0], old_pt.y, old_pt.z);
	}else if(deg_x){
		// count y and z
		DMatrix2d A2(A.m[1][1], A.m[1][2], A.m[2][1], A.m[2][2]);
		DVector2d b2(b.y, b.z), res;
		if(A2.solve(b2, res)){
			return DPoint3d(old_pt.x, res.x, res.y);
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Error solving (reduced for x) optimization system for node");
			return old_pt;
		}
	}else if(deg_y){
		// count x and z
		DMatrix2d A2(A.m[0][0], A.m[0][2], A.m[2][0], A.m[2][2]);
		DVector2d b2(b.x, b.z), res;
		if(A2.solve(b2, res)){
			return DPoint3d(res.x, old_pt.y, res.y);
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Error solving (reduced for y) optimization system for node");
			return old_pt;
		}
	}else if(deg_z){
		// count x and y
		DMatrix2d A2(A.m[0][0], A.m[0][1], A.m[1][0], A.m[1][1]);
		DVector2d b2(b.x, b.y), res;
		if(A2.solve(b2, res)){
			return DPoint3d(res.x, res.y, old_pt.z);
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Error solving (reduced for z) optimization system for node");
			return old_pt;
		}
	}
	// else - count all
	DVector3d res;
	if(A.solve(b, res)){
		return DPoint3d::zero + res;
	}else{
		LOG4CPLUS_WARN(MeshLog::logger_console, "Error solving optimization system for node");
		return old_pt;
	}
}

/// Close line with too close vertices
bool MeshGrain::closeLine(int i, DataVector<BlockInfo> & blocks,
		DataVector<FaceInfo> & faces, DataVector<LineInfo> & lines, 
		DataVector<int> & trans_nodes)
{
	LineInfo& line = lines[i];
	if(!line.valid) return false;
	// mark line
	line.valid = false;
	// fix nodes (remove node1)
	trans_nodes[line.node1] = line.node0;
	for(int j = 0; j < trans_nodes.countInt(); j++)
		if(trans_nodes[j] == line.node1)
			trans_nodes[j] = line.node0;
	// fix lines
	for(int j = 0; j < lines.countInt(); j++){
		if(lines[j].valid){
			if(lines[j].node0 == line.node1) lines[j].node0 = line.node0;
			if(lines[j].node1 == line.node1) lines[j].node1 = line.node0;
		}
	}
	// fix faces
	for(int j = 0; j < faces.countInt(); j++){
		FaceInfo& face = faces[j];
		if(!face.valid) continue;
		// ... lines
		bool incident_line = false;
		for(int k = 0; k < face.lines.countInt(); k++){
			if(face.lines[k] == i){
				if(k==0){
					face.lines[0] = face.lines[face.lines.countInt()-1];
					face.lines.removeLast();
				}else{
					face.lines.removeOrderedAt(k);
				}
				incident_line = true;
				break;
			}
		}
		// ... nodes if incident 
		for(int k = 0; k < face.nodes.countInt(); k++){
			if(face.nodes[k] == line.node1){
				if(incident_line)
					face.nodes.removeOrderedAt(k);
				else
					face.nodes[k] = line.node0;
				break;
			}
		}
		if(face.nodes.countInt() <= 2){
			// ... invalidate face and remove from all blocks
			face.valid = false;
			for(int k = 0; k < blocks.countInt(); k++){
				BlockInfo& block = blocks[k];
				for(int l = 0; l < block.faces.countInt(); l++){
					if(block.faces[l] == j){
						block.faces.removeAt(l);
						break;
					}
				}
			}
			// ... for degenerated cases
			if(face.nodes.countInt() == 2){
				// .. remove and replace one of the edges
				assert(face.lines.countInt() == 2);
				int line0 = face.lines[0];
				int line1 = face.lines[1];
				lines[line1].valid = false;
				lines[line1].trans_line = line0;
				for(int k = 0; k < faces.countInt(); k++){
					FaceInfo& kface = faces[k];
					if(!kface.valid) continue;
					for(int l = 0; l < kface.lines.countInt(); l++){
						if(kface.lines[l] == line1)
							kface.lines[l] = line0;
					}
				}
			}
		}
	}
	return true;
}

/// Check close nodes after optimization
int MeshGrain::recheckCloseNodes(DataVector<BlockInfo> & blocks,
		DataVector<FaceInfo> & faces, DataVector<LineInfo> & lines, 
		DataVector<DPoint3d> & nodes, DataVector<int> & trans_nodes)
{
	int close_node_count = 0;

	bool any_change = true;
	double DIST_TOL2 = 0.25 * param_grain_node_identity_tolerance * param_grain_node_identity_tolerance;
	while(any_change){
		any_change = false;
		// check lines directly
		for(int i = 0; i < lines.countInt(); i++){
			LineInfo& line = lines[i];
			if(!line.valid) continue;
			// check if vertices are too close to each other ?
			double dist2 = nodes[line.node0].distance2(nodes[line.node1]);
			if(dist2 >= DIST_TOL2) continue;
			// else, join ...
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, " rechecking, joining nodes: " << line.node0 << " and " << line.node1);
			LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "      (squared) distance = " << dist2);
			if(closeLine(i, blocks, faces, lines, trans_nodes)){
				++close_node_count;
				any_change = true;
			}
		}
		// check lines projected to (planar?) faces
		for(int i = 0; i < faces.countInt(); i++){
			FaceInfo& face = faces[i];
			if(!face.valid || !face.surface) continue;
			for(int j = 0; j < face.lines.countInt(); j++){
				LineInfo& line = lines[face.lines[j]];
				assert(line.valid);
				const DPoint3d pt0 = face.surface->getPoint(
					face.surface->getParameters( nodes[line.node0]));
				const DPoint3d pt1 = face.surface->getPoint(
					face.surface->getParameters( nodes[line.node1]));
				double dist2 = pt0.distance2(pt1);
				if(dist2 >= DIST_TOL2) continue;
				// else, join ...
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, " rechecking, joining nodes: " << line.node0 << " and " << line.node1);
				LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "      (squared projected) distance = " << dist2);
				if(closeLine(face.lines[j], blocks, faces, lines, trans_nodes)){
					++close_node_count;
					any_change = true;
					break;
				}
			}
		}
		// check lines projected to (planar?) faces
		for(int i = 0; i < faces.countInt(); i++){
			FaceInfo& face = faces[i];
			if(!face.valid || !face.surface) continue;
			int flct = face.lines.countInt();
			for(int j = 0; j < flct; j++){
				LineInfo& line = lines[face.lines[j]];
				LineInfo& prev_line = lines[face.lines[(j+flct-1)%flct]];
				LineInfo& next_line = lines[face.lines[(j+1)%flct]];
				assert(line.valid && prev_line.valid && next_line.valid);
				const DPoint3d pt0 = face.surface->getPoint(
					face.surface->getParameters( nodes[line.node0]));
				const DPoint3d pt1 = face.surface->getPoint(
					face.surface->getParameters( nodes[line.node1]));
				DPoint3d pt0_prev;
				DPoint3d pt1_next;
				if(prev_line.node0 == line.node0)
					pt0_prev = face.surface->getPoint(face.surface->getParameters( nodes[prev_line.node1]));
				else if(prev_line.node1 == line.node0)
					pt0_prev = face.surface->getPoint(face.surface->getParameters( nodes[prev_line.node0]));
				else if(next_line.node0 == line.node0)
					pt0_prev = face.surface->getPoint(face.surface->getParameters( nodes[next_line.node1]));
				else if(next_line.node1 == line.node0)
					pt0_prev = face.surface->getPoint(face.surface->getParameters( nodes[next_line.node0]));
				else assert(false);

				if(prev_line.node0 == line.node1)
					pt1_next = face.surface->getPoint(face.surface->getParameters( nodes[prev_line.node1]));
				else if(prev_line.node1 == line.node1)
					pt1_next = face.surface->getPoint(face.surface->getParameters( nodes[prev_line.node0]));
				else if(next_line.node0 == line.node1)
					pt1_next = face.surface->getPoint(face.surface->getParameters( nodes[next_line.node1]));
				else if(next_line.node1 == line.node1)
					pt1_next = face.surface->getPoint(face.surface->getParameters( nodes[next_line.node0]));
				else assert(false);

				const DVector3d vp = (pt0_prev - pt0).normalized();
				const DVector3d vt = (pt1 - pt0).normalized();
				const DVector3d vn = (pt1 - pt1_next).normalized();

				double s_pt = vp.scalarProduct(vt);
				double s_tn = vt.scalarProduct(vn);

				bool prev_ok = (abs(1.0 - s_pt) > SMALL_NUMBER);
				bool next_ok = (abs(1.0 - s_tn) > SMALL_NUMBER);

				if(prev_ok && next_ok) continue;

				double dist2 = pt0.distance2(pt1);
				
				if(!prev_ok && pt0_prev.distance2(pt0) < dist2) continue;
				if(!next_ok && pt1_next.distance2(pt1) < dist2) continue;

				if((param_grain_view & 1) == 1){
//				if(true){
					ostringstream oss;
					oss << "face #" << i << " line " << line.node0 << "-" << line.node1 
						<< ", dot_prev = " << s_pt
						<< ", dot_next = " << s_tn << " "
						<< ", dist2 = " << dist2;
					SHOW_MESH(oss.str(),
						getFaceViewSet(face, nodes, lines));
				}

				if(dist2 < 2*DIST_TOL2){
					// join ...
					LOG4CPLUS_DEBUG(MeshLog::logger_mesh, " rechecking, joining nodes: " << line.node0 << " and " << line.node1);
					LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "      (squared projected) distance = " << dist2);
					if(closeLine(face.lines[j], blocks, faces, lines, trans_nodes)){
						++close_node_count;
						any_change = true;
						break;
					}
				}else{
					// move away ...
					if(!prev_ok){
						// ... line.node1
						LOG4CPLUS_DEBUG(MeshLog::logger_mesh, " rechecking, moving node: " << line.node1);
						DVector3d dv = vp.crossProduct(face.surface->getNormalVector(DPoint2d::zero));
						if(dv.scalarProduct(vn) < 0.0) dv *= -1.0;
						dv *= param_grain_node_identity_tolerance;
						nodes[line.node1] += dv;
					}else{
						// ... line.node0
						LOG4CPLUS_DEBUG(MeshLog::logger_mesh, " rechecking, moving node: " << line.node0);
						DVector3d dv = vp.crossProduct(face.surface->getNormalVector(DPoint2d::zero));
						if(dv.scalarProduct(vp) < 0.0) dv *= -1.0;
						dv *= param_grain_node_identity_tolerance;
						nodes[line.node0] += dv;
					}
					if((param_grain_view & 1) == 1){
//					if(true){
						ostringstream oss;
						oss << "face #" << i << " line " << line.node0 << "-" << line.node1 
							<< ", after slight relocation";
						SHOW_MESH(oss.str(),
							getFaceViewSet(face, nodes, lines));
					}
				}
			}
		}
	}

	if(close_node_count > 0){
		LOG4CPLUS_WARN(MeshLog::logger_console, "Total number of extra joined nodes: " << close_node_count);
	}

	return close_node_count;
}

/// Check vertices for planar faces
int MeshGrain::recheckPlanarFaceVertices(DataVector<BlockInfo> & blocks,
		DataVector<FaceInfo> & faces, DataVector<LineInfo> & lines, 
		DataVector<DPoint3d> & nodes)
{
	int modifications = 0;

	for(int i = 0; i < faces.countInt(); i++){
		FaceInfo& face = faces[i];
		if(!face.valid || !face.planar() || face.nodes.countInt() < 4) continue;

		// ... create plane from vertices only and measure max distance
		int nct = face.nodes.countInt();
		DataVector<DPoint3d> all_points(nct);
		for(int j = 0; j < nct; j++)
			all_points.add(nodes[face.nodes[j]]);
		DPlane plane;
		double max_dist_v = DLeastSquaresFitting::fitHyperplaneOrthogonal(all_points, plane, true);

//		if(i == 469) SHOW_MESH("face 469", getFaceViewSet(face, nodes, lines));
		
		if(max_dist_v < 0.25) continue;
		// ... calculate max distance from pre-calculated plane
		DVector3d vn = face.surface->getNormalVector(DPoint2d::zero);
		double d = (face.surface->getPoint(DPoint2d::zero) - DPoint3d::zero).scalarProduct(vn);
		double max_dist_f = 0.0;
		int max_f = 0;
		for(int j = 0; j < nct; j++){
			double dist = abs((nodes[face.nodes[j]] - DPoint3d::zero).scalarProduct(vn) - d);
			if(dist > max_dist_f){
				max_dist_f = dist;
				max_f = j;
			}
		}
		// ... create plane from vertices only and without the most distant
		DataVector<DPoint3d> almost_all_points(nct-1);
		for(int j = 0; j < nct; j++){
			if(j != max_f)
				almost_all_points.add(nodes[face.nodes[j]]);
		}
		double max_dist_without = DLeastSquaresFitting::fitHyperplaneOrthogonal(
			almost_all_points, plane);

		vn = plane.vn;
		d = (plane.p0 - DPoint3d::zero).scalarProduct(vn);
		double node_dist = abs((nodes[face.nodes[max_f]] - DPoint3d::zero).scalarProduct(vn) - d);

		if(node_dist - max_dist_without < 0.5) continue;

		// add line joining the previous and next vertex
		int max_f_p = (max_f+nct-1)%nct;
		int max_f_n = (max_f+1)%nct;
		LineInfo nline(face.nodes[max_f_p], face.nodes[max_f_n], true);
		int line_id = lines.add(nline);

		if((param_grain_view & 1) == 1){
			ostringstream oss;
			oss << "face #" << i << " with optimized nodes, max_dist_v = " << max_dist_v
				<< ", max_dist_f = " << max_dist_f
				<< ", max_dist_v_without [node #" << face.nodes[max_f] << " "
				<< nline.node0 << "-" << nline.node1 << "] = (" 
				<< max_dist_without << " + " << node_dist << ")";
			SHOW_MESH(oss.str(),
				getFaceViewSet(face, nodes, lines));
		}

		// create additional face
		FaceInfo nface(face.block0, face.block1, true);
		// ... nodes
		nface.nodes.add(face.nodes[max_f_p]); // three nodes for new face
		nface.nodes.add(face.nodes[max_f]);
		nface.nodes.add(face.nodes[max_f_n]);
		face.nodes.removeOrderedAt(max_f); // remove one node from old face
		// ... lines
		nface.lines.add(line_id); // three lines for the new face
		nface.lines.add(face.lines[max_f]);
		nface.lines.add(face.lines[max_f_n]);
		face.lines.removeOrderedAt(max_f_n); // replace two edges in old face with one new
		face.lines[max_f_n > max_f ? max_f : max_f-1] = line_id;
		// ... calculate plane parameters for the new face
		//nface.plane_pt = nodes[nface.nodes[0]];
		//nface.plane_e0 = (nodes[nface.nodes[1]] - nodes[nface.nodes[0]]).normalized();
		//nface.plane_e1 = (nodes[nface.nodes[2]] - nodes[nface.nodes[0]]).normalized();
		//vn = nface.plane_e0.crossProduct(nface.plane_e1);
		//nface.plane_e1 = vn.crossProduct(nface.plane_e0).normalized();
		nface.surface = std::make_shared<SurfacePlane>(
			face.surface->getPoint(DPoint2d::zero),
			face.surface->getDerivative(DEquation::deriv_ds, DPoint2d::zero),
			face.surface->getDerivative(DEquation::deriv_dt, DPoint2d::zero));
		// ... add
		int face_id = faces.add(nface);
		// ... add to blocks
		if(nface.block0 >= 0){
			blocks[nface.block0].faces.add(face_id);
		}
		if(nface.block1 >= 0){
			blocks[nface.block1].faces.add(face_id);
		}

		//SHOW_MESH("old face", getFaceViewSet(faces[i], nodes, lines));
		//SHOW_MESH("new face", getFaceViewSet(nface, nodes, lines));
	}

	return modifications;
}

bool MeshGrain::FaceInfo::planar() const 
{ 
	return surface && 
		((surface->getType() == SURFACE_PLANE) ||
		((surface->getType() == SURFACE_CORRECTED) && 
		 ((SurfaceCorrected*)(surface.get()))->getBaseSurface()->getType() == SURFACE_PLANE)); 
}
