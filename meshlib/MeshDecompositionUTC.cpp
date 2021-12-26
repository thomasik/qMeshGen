/////////////////////////////////////////////////////////////////////////////
// MeshDecompositionUTC.cpp: implementation of the MeshDecompositionUTC class.
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	21-24 September 2005
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "MeshDecompositionUTC.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshPoint3d.h"
#include "MeshPoint2d.h"
#include "MeshArea.h"
#include "MeshFace.h"
#include "MeshBlock.h"
#include "MeshDomainVolume.h"
#include "DataVector.h"
#include "MeshViewSet.h"
#include "MeshTriangle3d.h"
#include "MeshQuad3d.h"
#include "MeshEdge2d.h"
#include "MeshEdge3d.h"
#include "SurfacePlane.h"
#include "ControlSpace2dMatrixUniform.h"
#include "ControlSpace3dOctree.h"
#include "MeshGenerator2d.h"
#include "DataHashTable.h"

MeshDecompositionUTC::MeshDecompositionUTC(const string& fname)
{
	m_mesh = 0;
	m_epsilon = 0.1;
	m_steep_edge_threshold = 0.1;
	m_visualization = 1;
	m_orientation_switch = false;

	ifstream file(fname.c_str());
	if(!file){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Can't open file: " << fname);
		return;
	}
	m_cut_mesh_name = "cut-mesh";
	while(file){
		string t;
		file >> t;
		if(!file) break;
		if(t.compare("MESH_NAME") == 0){
			file >> t;
			m_mesh = readIncidenceMeshDescription(t);
		}else if(t.compare("PLANE_NODE") == 0){
			file >> t;
			if(t.compare("inertial") == 0){
				m_plane_center = m_mesh->getInertialCenter();
			}else{
				assert(t.compare("explicit") == 0);
				file >> m_plane_center;
			}
		}else if(t.compare("PLANE_NORMAL") == 0){
			file >> m_plane_normal;
			m_plane_normal.normalize();
		}else if(t.compare("CUT_MESH_NAME") == 0){
			file >> m_cut_mesh_name;
		}else if(t.compare("VISUALIZATION") == 0){
			file >> m_visualization;
		}else if(t.compare("STEEP_EDGE_THRESHOLD") == 0){
			file >> t;
			DEquation eq;
			if(eq.parse(t.c_str()))	m_steep_edge_threshold = eq.getValue(0.0);
		}else if(t.compare("PENALTY_NODE_DIST") == 0){
			file >> t;
			DEquation eq;
			if(eq.parse(t.c_str()))	m_penalty_node_dist = eq.getValue(0.0);
		}else if(t.compare("PENALTY_NODE_CLASS") == 0){
			file >> t;
			DEquation eq;
			if(eq.parse(t.c_str()))	m_penalty_node_class = eq.getValue(0.0);
		}else if(t.compare("PENALTY_EDGE_LENGTH") == 0){
			file >> t;
			DEquation eq;
			if(eq.parse(t.c_str()))	m_penalty_edge_length = eq.getValue(0.0);
		}else if(t.compare("PENALTY_EDGE_ORIENT") == 0){
			file >> t;
			DEquation eq;
			if(eq.parse(t.c_str()))	m_penalty_edge_orient = eq.getValue(0.0);
		}else if(t.compare("PENALTY_EDGE_BAD_ORIENT") == 0){
			file >> t;
			DEquation eq;
			if(eq.parse(t.c_str()))	m_penalty_edge_bad_orient = eq.getValue(0.0);
		}else if(t.compare("PENALTY_EDGE_CROSS") == 0){
			file >> t;
			DEquation eq;
			if(eq.parse(t.c_str()))	m_penalty_edge_cross = eq.getValue(0.0);
		}else if(t.compare("NODE_DIST_EPSILON") == 0){
			file >> t;
			DEquation eq;
			if(eq.parse(t.c_str()))	m_epsilon = eq.getValue(0.0);
		}else{
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "CUT_MESH -> unknown property");
			LOG4CPLUS_ERROR(MeshLog::logger_console,   t.c_str());
		}
	}
	m_log_file.open((m_cut_mesh_name + ".log").c_str());
}


MeshDecompositionUTC::~MeshDecompositionUTC()
{ 
	delete m_mesh; 
}

double MeshDecompositionUTC::countNodePenalty(CutNode& cn){
	cn.penalty = 
		m_penalty_node_dist * cn.dist +
		m_penalty_node_class * std::abs(cn.mpoint->getIntTag(TagExtended::TAG_UTC_CUT_PENALTY));

	return cn.penalty;
}

double MeshDecompositionUTC::countEdgePenalty(CutEdge& ce){
	ce.penalty = 
			m_penalty_edge_length * ce.length +
			m_penalty_edge_orient * (1.0 - ce.orientation) + 
			m_penalty_edge_cross  * ce.crossing;

	if(ce.orientation <= 0.0)
		ce.penalty += m_penalty_edge_bad_orient;

	return ce.penalty;
}

/**
 * Finds loops of edges for conforming cutting contour
 * \param mesh - 3d surface triangular mesh (closed)
 * \return true if successful (and marks nodes and edges in mesh)
 */
DecompEdgeArrayList MeshDecompositionUTC::findContour3D(DecompNodeList* node_list, DecompFaceList* cut_faces)
{
	int node_list_count = node_list->countInt();
	CutNode* cut_nodes = new CutNode[node_list_count];
	int cut_faces_count = cut_faces->countInt();
	CutEdge* cut_edges = new CutEdge[3*cut_faces_count];

	DataHashTableKeyValue<MeshPoint3d*, int> kv_points(node_list_count, 0);

	double max_node_dist = 0.0;
	for(int i = 0; i < node_list_count; i++){
		cut_nodes[i].mpoint = node_list->get(i);
		cut_nodes[i].dist = abs(cut_nodes[i].mpoint->getDoubleTag(TagExtended::TAG_UTC_CUT_DIST));
		kv_points.insert(cut_nodes[i].mpoint, i); // for back-translation
		if(cut_nodes[i].dist > max_node_dist)
			max_node_dist = cut_nodes[i].dist;
	}

	// Scaling of distance
	for(int i = 0; i < node_list_count; i++){
		cut_nodes[i].dist /= max_node_dist;
	}

	// Gather edges
	for(int i = 0; i < cut_faces_count; i++){
		MeshFace* face = cut_faces->get(i);
		int ect = face->getEdgeCount();
		for(int j = 0; j < ect; j++) 
			face->getEdge(j)->setIntTag(TagExtended::TAG_UTC_CUT);
	}

	// both directed/undirected graph 
	double max_edge_len = 0.0;
	int cut_edges_count = 0;
	for(int i = 0; i < cut_faces_count; i++){
		MeshFace* face = cut_faces->get(i);
		const DVector3d direction = m_plane_normal.crossProduct(face->getNormalVector());
		int ect = face->getEdgeCount();
		MeshPoint3d* last_point = face->getPoint(0);
		for(int j = 0; j < ect; j++){
			MeshEdge3d* edge = face->getEdge(j);
			MeshPoint3d* point = face->getPoint((j+1)%ect);
			if(edge->availableTag(TagExtended::TAG_UTC_CUT)){
				const DVector3d edge_vec = point->getCoordinates() - last_point->getCoordinates();
				cut_edges[cut_edges_count].medge = edge;
				cut_edges[cut_edges_count].length = edge_vec.length();
				max_edge_len = std::max(max_edge_len, cut_edges[cut_edges_count].length);
				cut_edges[cut_edges_count].orientation = 
					direction.scalarProduct(edge_vec / cut_edges[cut_edges_count].length);
				cut_edges[cut_edges_count].crossing = 
					((last_point->getIntTag(TagExtended::TAG_UTC_CUT_PENALTY) * 
						   point->getIntTag(TagExtended::TAG_UTC_CUT_PENALTY) < 0) ? 1 : 0);
				if(cut_edges[cut_edges_count].orientation >= 0.0){
					cut_edges[cut_edges_count].node0 = kv_points.getValue(last_point, 0);
					cut_edges[cut_edges_count].node1 = kv_points.getValue(point, 0);
				}else{
					cut_edges[cut_edges_count].orientation *= -1.0;
					cut_edges[cut_edges_count].node0 = kv_points.getValue(point, 0);
					cut_edges[cut_edges_count].node1 = kv_points.getValue(last_point, 0);
				}
				cut_nodes[cut_edges[cut_edges_count].node0].edges_out.add(cut_edges_count);
				cut_nodes[cut_edges[cut_edges_count].node1].edges_in.add(cut_edges_count);
				// mark edge
				edge->removeTag(TagExtended::TAG_UTC_CUT);
				++cut_edges_count;
			}
			last_point = point;
		}
	}

	// Scaling of edge lengths
	for(int i = 0; i < cut_edges_count; i++){
		cut_edges[i].length /= max_edge_len;
	}


	// graph improvement: each node at least one edge in and one out
	for(int i = 0; i < node_list_count; i++){
		int rank_out = cut_nodes[i].edges_out.countInt();
		int rank_in  = cut_nodes[i].edges_in.countInt();
		assert(rank_out + rank_in > 1);
		if(rank_out == 0){
			// find min orientation
			//int min_orient = 0;
			int min_ce_i = cut_nodes[i].edges_in.get(0);
			for(int j = 1; j < rank_in; j++){
				int ce_i = cut_nodes[i].edges_in.get(j);
				if(cut_edges[ce_i].orientation < cut_edges[min_ce_i].orientation){
					min_ce_i = ce_i;
				}
			}
			// switch edge
			assert(i == cut_edges[min_ce_i].node1);
			cut_edges[min_ce_i].orientation *= -1.0;
			int i_other = cut_edges[min_ce_i].node0;
			cut_edges[min_ce_i].node0 = i;
			cut_edges[min_ce_i].node1 = i_other;
			cut_nodes[i].edges_in.remove(min_ce_i);
			cut_nodes[i].edges_out.add(min_ce_i);
			cut_nodes[i_other].edges_out.remove(min_ce_i);
			cut_nodes[i_other].edges_in.add(min_ce_i);
		}else if(rank_in == 0){
			// find min orientation
			// int min_orient = 0;
			int min_ce_i = cut_nodes[i].edges_out.get(0);
			for(int j = 1; j < rank_out; j++){
				int ce_i = cut_nodes[i].edges_out.get(j);
				if(cut_edges[ce_i].orientation < cut_edges[min_ce_i].orientation){
					min_ce_i = ce_i;
				}
			}
			// switch edge
			assert(i == cut_edges[min_ce_i].node0);
			cut_edges[min_ce_i].orientation *= -1.0;
			int i_other = cut_edges[min_ce_i].node1;
			cut_edges[min_ce_i].node0 = i_other;
			cut_edges[min_ce_i].node1 = i;
			cut_nodes[i].edges_out.remove(min_ce_i);
			cut_nodes[i].edges_in.add(min_ce_i);
			cut_nodes[i_other].edges_in.remove(min_ce_i);
			cut_nodes[i_other].edges_out.add(min_ce_i);
		}
	}

	// find separate cycles
	int colored_count = 0;
	int color_index = -1;
	DataVector<int> color_counters;
	while(colored_count < node_list_count){
		int i = 0;
		while(cut_nodes[i].color >= 0) i++;
		DataVector<int> not_visited(node_list_count);
		cut_nodes[i].color = ++color_index;
		not_visited.add(i);
		int current_color_count = 1;
		while(not_visited.countInt() > 0){
			int k = not_visited.removeLast();
			int rank_out = cut_nodes[k].edges_out.countInt();
			for(int j = 0; j < rank_out; j++){
				int k1 = cut_edges[cut_nodes[k].edges_out[j]].node1;
				if(cut_nodes[k1].color < 0){
					++current_color_count;
					cut_nodes[k1].color = color_index;
					not_visited.add(k1);
				}
			}
			int rank_in = cut_nodes[k].edges_in.countInt();
			for(int j = 0; j < rank_in; j++){
				int k1 = cut_edges[cut_nodes[k].edges_in[j]].node0;
				if(cut_nodes[k1].color < 0){
					++current_color_count;
					cut_nodes[k1].color = color_index;
					not_visited.add(k1);
				}
			}
		}
		colored_count += current_color_count;
		color_counters.add(current_color_count);
	}

	int cut_count = color_counters.countInt();
	m_log_file << cut_count << " connected subgraphs in cut-contour.\n";

	DecompEdgeArrayList cut_contours(cut_count);

	for(int i = 0; i < node_list_count; i++){
		countNodePenalty(cut_nodes[i]);
		cut_nodes[i].min_flow_penalty = 1e20;
	}
	for(int i = 0; i < cut_edges_count; i++){
		countEdgePenalty(cut_edges[i]);
	}

	// find min edge for each cycle (+ make it a bridge-edge) -> for directed graph
	for(int k = 0; k < cut_count; k++){
		int min_edge_i = -1;
		double min_penalty = 0.0;

		for(int i = 0; i < node_list_count; i++){
			if(cut_nodes[i].color == k){
				double node_penalty = cut_nodes[i].penalty;
				int rank = cut_nodes[i].edges_out.countInt();
				for(int j = 0; j < rank; j++){
					int ce_i = cut_nodes[i].edges_out[j];
					double penalty = node_penalty + 
						cut_edges[ce_i].penalty + 
						cut_nodes[cut_edges[ce_i].node1].penalty;
					if(min_edge_i < 0 || min_penalty > penalty){
						min_edge_i = ce_i;
						min_penalty = penalty;
					}
				}
			}
		}

		cut_edges[min_edge_i].medge->setIntTag(TagExtended::TAG_UTC_CUT);
		cut_edges[min_edge_i].medge->setBorder();
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Minimum edge with penalty " << min_penalty);
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, " - node0: " << cut_nodes[cut_edges[min_edge_i].node0].penalty);
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, " - node1: " << cut_nodes[cut_edges[min_edge_i].node1].penalty);
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, " - edge : " << cut_edges[min_edge_i].penalty);

		cut_edges[min_edge_i].tag = 1; // part of cut-contour

//		int rank = cut_graph[min_edge_i].directed_edges.countInt();
//		for(int j = 0; j < rank; j++){
//			if(j != min_edge_j) 
//				cut_graph[min_edge_i].directed_edges[j].tag = 2; // edge forbidden
//		}

		// Mark all nodes with the minimum penalty flow from the begin point (directed graph)
		int end_point_i = cut_edges[min_edge_i].node0;
		int begin_point_i = cut_edges[min_edge_i].node1;
		cut_nodes[begin_point_i].min_flow_penalty = cut_nodes[begin_point_i].penalty;
		int current_color_count = color_counters[k];
		DataVector<int> to_visit(current_color_count);
		to_visit.add(begin_point_i);
		int visited_count = 0;
		while(to_visit.countInt() > 0){
			int min_i = 0;
			double lowest_penalty = cut_nodes[to_visit[min_i]].min_flow_penalty;
			for(int i = 1; i < to_visit.countInt(); i++){
				if(cut_nodes[to_visit[i]].min_flow_penalty < lowest_penalty){
					min_i = i; 
					lowest_penalty = cut_nodes[to_visit[i]].min_flow_penalty;
				}
			}
			int visited_i = to_visit.removeAt(min_i);
			++visited_count;
			cut_nodes[visited_i].visited = true;
			int rank = cut_nodes[visited_i].edges_out.countInt();
			for(int i = 0; i < rank; i++){
				int ce_i = cut_nodes[visited_i].edges_out[i];
				//cut_edges[ce_i].medge->setTag(); //--------
				int n1 = cut_edges[ce_i].node1;
				if(!cut_nodes[n1].visited){
					double penalty = lowest_penalty + cut_edges[ce_i].penalty + cut_nodes[n1].penalty;
					if(penalty < cut_nodes[n1].min_flow_penalty){
						cut_nodes[n1].min_flow_penalty = penalty;
					}
					to_visit.addIfNew(n1);
				}
			}
		}
		if(visited_count != current_color_count){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   "Couldn't visit all nodes for this color!");
		}

		int contour_index = color_counters[k];
		MeshEdge3d** contour = new MeshEdge3d*[contour_index+1];

		// Backtrack the directed graph and find (one of) the shortest flow
		int visit_i = end_point_i;
		do{
			assert(cut_nodes[visit_i].visited);
			// select correct back-route
			int rank_in = cut_nodes[visit_i].edges_in.countInt();
			assert(rank_in > 0);
			int best_ce = cut_nodes[visit_i].edges_in.get(0);
			double best_diff = cut_nodes[visit_i].min_flow_penalty
				- cut_nodes[visit_i].penalty
				- cut_edges[best_ce].penalty 
				- cut_nodes[cut_edges[best_ce].node0].min_flow_penalty;
			for(int j = 1; j < rank_in; j++){
				int ce = cut_nodes[visit_i].edges_in.get(j);
				double diff = cut_nodes[visit_i].min_flow_penalty
					- cut_nodes[visit_i].penalty
					- cut_edges[ce].penalty 
					- cut_nodes[cut_edges[ce].node0].min_flow_penalty;
				if(diff > best_diff){
					best_diff = diff;
					best_ce = ce;
				}
			}
			if(cut_edges[best_ce].orientation < 0)
                cut_edges[best_ce].medge->setIntTag(TagExtended::TAG_UTC_CUT, 2); //2 for back-oriented
			else if(abs(cut_edges[best_ce].orientation) < m_steep_edge_threshold)
				cut_edges[best_ce].medge->setIntTag(TagExtended::TAG_UTC_CUT, 3); //3 for steep
			else cut_edges[best_ce].medge->setIntTag(TagExtended::TAG_UTC_CUT, 1); //1 for normal
			contour[--contour_index] = cut_edges[best_ce].medge;
			visit_i = cut_edges[best_ce].node0;
		}while(visit_i != begin_point_i);
		int rank_in = cut_nodes[visit_i].edges_in.countInt();
		for(int j = 0; j < rank_in; j++){
			int ce = cut_nodes[visit_i].edges_in.get(j);
			if(cut_edges[ce].node0 == end_point_i){
				cut_edges[ce].medge->setIntTag(TagExtended::TAG_UTC_CUT, 
					(cut_edges[ce].orientation >= 0) ? 1 : 2); //2 for back-oriented
				contour[--contour_index] = cut_edges[ce].medge;
				break;
			}
		}

		int last_i = color_counters[k];
		if(contour_index > 0){
			int ind = 0;
			for(int i = contour_index; i < last_i; contour[ind++] = contour[i++]);
			contour[ind] = nullptr;
		}else contour[last_i] = nullptr;

		m_log_file << "Contour " << k << ": " << (last_i-contour_index);
		m_log_file << " nodes selected out of " << last_i << " total in this cut.\n";

		cut_contours.add(contour);
	}
	// check for connected cycles

	// log graph
	/*
	for(int i = 0; i < node_list_count; i++){
		int rank_out = cut_nodes[i].edges_out.countInt();
		int rank_in  = cut_nodes[i].edges_in.countInt();
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, i << "[" << cut_nodes[i].mpoint->getIndex() << "] (" << rank_in << "/" << rank_out << ") -> ";
		if(rank_out == 0) node_list->get(i)->setTag(3);
		for(int j = 0; j < rank_out; j++){
			int ce_i = cut_nodes[i].edges_out[j];
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "[" << cut_edges[ce_i].node1;
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, ", orient=" << cut_edges[ce_i].orientation;
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, ", len=" << cut_edges[ce_i].length << "] ";
		}
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;
	}
	*/

	delete[] cut_nodes;
	delete[] cut_edges;

	return cut_contours;
}

/**
 * Marks points according to the cutting plane
 */
void MeshDecompositionUTC::markPointsOrientation()
{
	int p_count = m_mesh->getPointsCount();
	for(int i = 0; i < p_count; i++){
		MeshPoint3d* point = m_mesh->getPointAt(i);
		// min length of incident edges
		double min_length = point->getEdge(0)->getLengthNoMetric();
		int rank = point->getRank();
		for(int j = 1; j < rank; j++)
		    min_length = std::min(min_length, point->getEdge(j)->getLengthNoMetric());
		// distance of point[i] from plane
		double dist = m_plane_normal.scalarProduct(point->getCoordinates()-m_plane_center);
		point->setDoubleTag(TagExtended::TAG_UTC_CUT_DIST, dist);
		if(abs(dist) < m_epsilon*min_length){
			point->setIntTag(TagExtended::TAG_UTC_CUT_PENALTY, 0);
		}else{
			point->setIntTag(TagExtended::TAG_UTC_CUT_PENALTY, (dist>0) ? 1 : -1);
		}
	}
}

DecompFaceList* MeshDecompositionUTC::selectCutTriangles()
{
	MeshBlock* block = m_mesh->getBlockAt(0);
	assert(block);
	int fct = block->getFaceCount();
	DecompFaceList* face_list = new DecompFaceList(10+(int)sqrt((double)fct));
	double min_angle = PI;
	double max_angle = 0.0;
	int ranges[4] = {0, 0, 0, 0};
	for(int i = 0; i < fct; i++){
		MeshFace* face = block->getFace(i);
		int types[] = {0, 0, 0};
		int ect = face->getEdgeCount();
		for(int j = 0; j < ect; j++){
			types[face->getPoint(j)->getIntTag(TagExtended::TAG_UTC_CUT_PENALTY)+1]++;
		}
		if((types[0] == 2) && (types[2] == 1))		face->setIntTag(TagExtended::TAG_UTC_CUT, 0);
		else if((types[0] == 1) && (types[2] == 2)) face->setIntTag(TagExtended::TAG_UTC_CUT, 1);
		else if((types[1] == 2) && (types[0] == 1)) face->setIntTag(TagExtended::TAG_UTC_CUT, 2);
		else if((types[1] == 2) && (types[2] == 1)) face->setIntTag(TagExtended::TAG_UTC_CUT, 3);
		else if((types[0] == 1) && (types[2] == 1)) face->setIntTag(TagExtended::TAG_UTC_CUT, 4);
		else if(types[2] > 0) face->setIntTag(TagExtended::TAG_UTC_CUT, 5);  // not interesting
		else face->setIntTag(TagExtended::TAG_UTC_CUT, 6); // not interesting

		if(face->getIntTag(TagExtended::TAG_UTC_CUT) < 5){
	 	    face_list->add(face);
			double angle = m_plane_normal.getAngle(face->getNormalVector());
			if(angle > PI/2) angle = PI - angle;
			if(angle > max_angle) max_angle = angle;
			if(angle < min_angle) min_angle = angle;
			if(angle < PI/6) ranges[0]++;
			else if(angle < PI/4) ranges[1]++;
			else if(angle < PI/3) ranges[2]++;
			else ranges[3]++;
		}
	}
	m_log_file << face_list->countInt() << " faces cut by plane.\n";
	m_log_file << " - min angle = " << (min_angle * 180/PI) << endl;
	m_log_file << " - max angle = " << (max_angle * 180/PI) << endl;
	m_log_file << " -- number of faces with angle [ 0,30): " << ranges[0] << endl;
	m_log_file << " -- number of faces with angle [30,45): " << ranges[1] << endl;
	m_log_file << " -- number of faces with angle [45,60): " << ranges[2] << endl;
	m_log_file << " -- number of faces with angle [60,90): " << ranges[3] << endl;

	storeCutFaces(face_list);

	return face_list;
}

DecompNodeList* MeshDecompositionUTC::selectPoints(DecompFaceList* cut_faces)
{
	int pct = m_mesh->getPointsCount();
	DecompNodeList* pts_list = new DecompNodeList(10+(int)sqrt((double)pct));
	DataVector<bool> pts_near_plane(pct, false);

	int fct = cut_faces->countInt();
	for(int i = 0; i < fct; i++){
		MeshFace* face = cut_faces->get(i);
		int ect = face->getEdgeCount();
		for(int j = 0; j < ect; j++)
		    pts_near_plane[face->getPoint(j)->getIndex()] = true;
	}
	for(int i = 0; i < pct; i++){
		if(pts_near_plane[i])
	 	    pts_list->add(m_mesh->getPointAt(i));
		else
 		    m_mesh->getPointAt(i)->setIntTag(TagExtended::TAG_UTC_CUT_PENALTY, -10);
	}

	m_log_file << pts_list->countInt() << " total nodes in faces cut by plane.\n";
	return pts_list;
}

MeshContainer3d* MeshDecompositionUTC::readIncidenceMeshDescription(const string& fname)
{
	ostringstream fname1, fname2, fname_cs;
	fname1 << fname << "-p.txt";
	fname2 << fname << "-i.txt";
	fname_cs << fname << "-cs.txt";

	ifstream fpoints(fname1.str().c_str());
	ifstream felements(fname2.str().c_str());

	int p_count;
	fpoints >> p_count;
	MeshContainer3d* mesh = new MeshContainer3d(p_count);

	DataVector<MeshPoint3d*> points(p_count);
	DataHashTableKeyValue<int, MeshPoint3d*> kv_PIds(p_count, -1);

	for(int i = 0; i < p_count; i++){
		int id;
		double x,y,z;
		fpoints >> id >> x >> y >> z;
		MeshPoint3d* point = new MeshPoint3d(x, y, z);
		points.add(point);
		mesh->addMeshPoint(point);
		kv_PIds.insert(id, point);
	}

	string buffer_line;
	int f_count;
	getline(felements, buffer_line);
	{
		istringstream iss(buffer_line);
		iss >> f_count;
	}
	DataVector<MeshFace*> faces(f_count);
	for(int i = 0; i < f_count; i++){
		getline(felements, buffer_line);
		istringstream iss(buffer_line);
		int id, ect;
		int PIds[4];
		iss >> id >> ect;
		assert((ect == 3) || (ect == 4));
		for(int j = 0; j < ect; j++)
			iss >> PIds[j];
		MeshFace* face = nullptr;
		if(ect == 3){
			face = new MeshTriangle3d(
				kv_PIds.getValue(PIds[0], 0),
				kv_PIds.getValue(PIds[1], 0),
				kv_PIds.getValue(PIds[2], 0));
		}else{
			face = new MeshQuad3d(
				kv_PIds.getValue(PIds[0], 0),
				kv_PIds.getValue(PIds[1], 0),
				kv_PIds.getValue(PIds[2], 0),
				kv_PIds.getValue(PIds[3], 0));
		}
		faces.add(face);
		//face->setTag(id);
	}

	MeshDomainVolume* block = new MeshDomainVolume(faces, points);
	mesh->addMeshBlock(block);

	auto cs_space = std::make_shared<ControlSpace3dOctree>(mesh->getBoundingBox());
	if(cs_space->loadTXT(fname_cs.str().c_str())){
		mesh->setControlSpace(cs_space);
		LOG4CPLUS_INFO(MeshLog::logger_console, "Control space read successfully.");
	}

	return mesh;
}

MeshViewSet* MeshDecompositionUTC::getViewSet(MeshViewSet* view_set)
{
	if(!m_mesh) return view_set;
	int pct = m_mesh->getPointsCount();
	MeshBlock* block = m_mesh->getBlockAt(0);
	assert(block);
	int fct = block->getFaceCount();

	if(view_set){
		view_set->prepareFreePlace(pct, 2*fct, fct, 0);
	}else{
		view_set = new MeshViewSet(pct, 2*fct, fct, 0);
	}

	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = m_mesh->getPointAt(i);
		if(point->getIntTag(TagExtended::TAG_UTC_CUT_PENALTY) > -10){
		    view_set->addPoint(point, point->getIntTag(TagExtended::TAG_UTC_CUT_PENALTY)+1);

		}
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(edge->availableTag(TagExtended::TAG_UTC_CUT) && (edge->getPointIndex(point) == 0)){
				view_set->addEdge(edge);
			}
		}
	}

	const DVector3d shift = m_plane_normal * (0.05 * m_mesh->getBoundingBox().getDiameter());
	for(int i = 0; i < fct; i++){
		MeshFace* face = block->getFace(i);
		int tag = face->getIntTag(TagExtended::TAG_UTC_CUT);
		if(tag == 5)
		    view_set->addFace(face, shift, tag);
		else if(tag == 6)
		    view_set->addFace(face, shift * -1.0, tag);
		else
		    view_set->addFace(face, tag);
	}

	return view_set;
}

MeshContainer2d* MeshDecompositionUTC::createCutMesh(DecompEdgeArrayList& contours)
{
	// for each contour create seam-mesh
	int contour_count = contours.countInt();
	int total_node_count = 0;
	DataVector<int> contour_node_count(contour_count);
	for(int i = 0; i < contour_count; i++){
		int count = 0;
		for(MeshEdge3d** edges = contours[i]; *edges; edges++) 
			++count;
		contour_node_count.add(count);
		total_node_count += count;
	}

	MeshContainer2d* boundary = new MeshContainer2d(total_node_count);
	DVector3d vs;
	if(m_plane_normal.x == 0) vs = DVector3d(1.0, 0.0, 0.0);
	else if(m_plane_normal.y == 0) vs = DVector3d(0.0, 1.0, 0.0);
	else if(m_plane_normal.z == 0) vs = DVector3d(0.0, 0.0, 1.0);
	else if(abs(m_plane_normal.y) > abs(m_plane_normal.z)) 
		vs = DVector3d(1.0, -m_plane_normal.x/m_plane_normal.y, 0.0);
	else vs = DVector3d(1.0, 0.0, -m_plane_normal.x/m_plane_normal.z);
	DVector3d vt = m_plane_normal.crossProduct(vs);

	SurfacePtr plane(new SurfacePlane(m_plane_center, vs, vt));
	boundary->setSurface(plane);

	struct ControlInfo{
		DPoint2d pt;
		ControlDataStretch2d stretch;
	} * control_data = new ControlInfo[total_node_count];

	DRect bounding_rect;

	int control_i = 0;
	for(int i = 0; i < contour_count; i++){
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "===> Contour << " << i);
		DataVector<MeshPoint2d*> plane_points(contour_node_count[i]);
		MeshEdge3d** edges = contours[i];

		MeshPoint3d* mp = edges[0]->getMeshPoint(0);
		if(edges[1]->getPointIndex(mp) >= 0){
            mp = edges[0]->getMeshPoint(1);
			assert(edges[1]->getPointIndex(mp) < 0);
		}
		int cct = contour_node_count[i];
		DataVector<MeshPoint3d*> mesh_points3d(cct);
		for(int j = 0; j < cct; j++){
			mesh_points3d.add(mp);
//			LOG4CPLUS_INFO(MeshLog::logger_mesh, " * node " << mp->getID());
			DPoint2d pt;
			int prev_j = (j+cct-1)%cct;
/*
			DPoint3d pt3 = mp->getCoordinates();
			DPoint3d pt3_back = edges[prev_j]->getOtherPoint(mp)->getCoordinates();
			DPoint3d pt3_forward = edges[j]->getOtherPoint(mp)->getCoordinates();
			DPoint3d vt_forward = pt3_forward - pt3;
			DPoint3d vt_back = pt3_back - pt3;
			DPoint3d pt3_end = pt3 + vt_forward + vt_back;
			// DPoint3d pt3_end = pt3 + vt_forward.normalized() + vt_back.normalized();
			bool result = plane->getSegmentCrossingPointParam(pt3, pt3_end, pt);
			assert(result == true);
*/

			if(edges[prev_j]->getIntTag(TagExtended::TAG_UTC_CUT) == 1){ // normal edge
				pt = plane->getParameters(mp->getCoordinates());
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, " -- backward edge fixing");
				m_log_file << " --> node " << j << " moved (backward edge)" << endl;

				//assert(edges[prev_j]->getTag() == 2); // backward edge
				DPoint3d back_pt = edges[prev_j]->getOtherPoint(mp)->getCoordinates();
				DPoint3d forward_pt = edges[j]->getOtherPoint(mp)->getCoordinates();
				pt = plane->getParameters(DPoint3d::average(back_pt, forward_pt));
			}

			bounding_rect.addPoint(pt);
			MeshPoint2d* mptx = new MeshPoint2d(pt);
			plane_points.add(mptx);
			mptx->setBorder();
			mptx->setPtrTag(TagExtended::TAG_MP_2D_3D, mp);
			boundary->addMeshPoint(mptx);
			mp = edges[j]->getOtherPoint(mp);
			control_data[control_i+j].stretch.ly = edges[j]->getLengthNoMetric();
		}

		// show surface
		MeshViewSet *view_set = new MeshViewSet(cct*2, cct*2, 0, 0);
		MeshPoint2d*   last_mpt  = plane_points[cct-1];
		MeshPoint3d* last_mpt3 = mesh_points3d[cct-1];
		for(int j = 0; j < cct; j++){
			MeshPoint2d*   mpt  = plane_points[j];
			MeshPoint3d* mpt3 = mesh_points3d[j];
			view_set->addPoint(plane->getPoint(mpt->getCoordinates()), 1);
			view_set->addPoint(mpt3->getCoordinates(), 2);
			view_set->addEdge(plane->getPoint(last_mpt->getCoordinates()),
				plane->getPoint(mpt->getCoordinates()), 0);
			view_set->addEdge(last_mpt3->getCoordinates(), 
				mpt3->getCoordinates(), 1);
			last_mpt = mpt;
			last_mpt3 = mpt3;
		}
		SHOW_MESH("contour projection", view_set);

/*
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "=========================================");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Testing boundary nodes for cycle " << i);
		double min_dist = 1e20;
		for(int j = 0; j < cct; j++){
			const DPoint2d pt = plane_points[j]->getCoordinates();
			const DPoint2d prev_pt = plane_points[(j+cct-1)%cct]->getCoordinates();
			const DPoint2d next_pt = plane_points[(j+1)%cct]->getCoordinates();
			double dist = pt.distance(prev_pt);
			if(dist < min_dist) min_dist = dist;
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "* node " << j << " back-dist = " << dist;
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "  angle = " << pt.getAngle(prev_pt, next_pt));
		}
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "*** Min distance = " << min_dist);
*/
		for(int j = 0; j < cct; j++){
			MeshPoint2d* mp0 = plane_points[j];
			MeshPoint2d* mp1 = plane_points[(j+1)%cct];
			// set edge
			MeshEdge2d* edge = mp0->getEdgeToPoint(mp1);
			if(edge){
				LOG4CPLUS_WARN(MeshLog::logger_console, "Edge incident to two contours ?!");
			}else{
				edge = new MeshEdge2d(mp0, mp1);
				edge->setBorder(TagBorder::OUTER); // i.e. normal boundary
			}
			// count control measure
			const DPoint2d& pt0 = mp0->getCoordinates();
			const DPoint2d& pt1 = mp1->getCoordinates();
			const DVector2d vdt = pt1 - pt0;
			control_data[control_i+j].stretch.lx = vdt.length();
			control_data[control_i+j].stretch.angle = 
				(abs(vdt.y) > mesh_data.relative_small_number) ? atan(vdt.y/vdt.x) : PI/2;
			control_data[control_i+j].pt = DPoint2d::average(pt0, pt1);
		}

		MeshArea* area = new MeshArea(plane_points);
		area->setFilled(MeshArea::getOrientation(plane_points));
		area->setAreaID(0);
		boundary->addMeshElement(area);
		control_i += cct;

		if(m_visualization > 0){
			SHOW_MESH("Contour area (2d)", area->getViewSet(0, 0));
		}

		if(m_visualization > 0){
			SHOW_MESH("Contour area (3d)", area->getViewSet(0, plane));
		}

	}

	auto space = MeshGenerator2d::createNewControlSpace(plane, bounding_rect);
	auto acs = space->getAsAdaptive();
	if (acs) {
		acs->setMaxMetric();
		acs->adaptToParameterization();
		if (m_mesh->getControlSpace())
			acs->applyAsMinimum(m_mesh->getControlSpace());
		for (int i = 0; i < total_node_count; i++) {
			acs->setMinControl(control_data[i].pt, DMetric2d::stretchToMatrix(control_data[i].stretch));
		}
		acs->smoothen();
	}
	boundary->setControlSpace(space);
	
	// triangulate
	MeshContainer2d *cut_mesh = MeshGenerator2d::createInitialMesh(boundary, &bounding_rect);
	if(cut_mesh == 0){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Error crating the initial cut-mesh");
		return 0;
	}
	cut_mesh->setControlSpace(space);
	Metric2dContext mc(space);
	int ct = MeshGenerator2d::addBoundaryNodes(mc, cut_mesh, boundary);
	if(ct < total_node_count){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Couldn't add all boundary nodes for cut-mesh");
		return 0;
	}

	bool success = MeshGenerator2d::constrainToBorder(mc, boundary, cut_mesh);
	if(!success){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Failed constraining to boundary for cut-mesh");
		return 0;
	}
	delete boundary;

	if(m_visualization > 1){
		SHOW_MESH("After constraining (2d)", cut_mesh->getViewSet(0, false));
	}

	if(m_visualization > 1){
		SHOW_MESH("After constraining (3d)", cut_mesh->getViewSet(0, true));
	}

	cut_mesh->clearSearchTree();
	MeshGenerator2d::addInnerNodes(mc, cut_mesh);
	MeshGenerator2d::smoothen(mc, cut_mesh, 2);

	// TODO check for "boundary triangles already existing"
	int ect = cut_mesh->getElementsCount();
	for(int i = 0; i < ect; i++){
		MeshElement* element = cut_mesh->getElementAt(i);
		MeshPoint3d* mpt0 = (MeshPoint3d*)element->getPoint(0)->getPtrTag(TagExtended::TAG_MP_2D_3D);
		if(!mpt0) continue;
		MeshPoint3d* mpt1 = (MeshPoint3d*)element->getPoint(1)->getPtrTag(TagExtended::TAG_MP_2D_3D);
		if(!mpt1) continue;
		MeshPoint3d* mpt2 = (MeshPoint3d*)element->getPoint(2)->getPtrTag(TagExtended::TAG_MP_2D_3D);
		if(!mpt2) continue;
		m_log_file << " ! cut-element with all boundary-nodes: [";
		m_log_file << mpt0->getIndex() << ',';
		m_log_file << mpt1->getIndex() << ',';
		m_log_file << mpt2->getIndex() << "] -> ";
		if(mpt0->getFaceToPoints(mpt1, mpt2)){
			m_log_file << "face already in 3D mesh" << endl;
		}else{
			m_log_file << "no such face in 3D mesh" << endl;
		}
	}


	if(m_visualization > 0){
		SHOW_MESH("Final mesh", cut_mesh->getViewSet());
	}

	// clear obsolete data
	for(int i = 0; i < contour_count; i++){
		delete[] contours[i];
	}

	delete[] control_data;

	return cut_mesh;
}

void MeshDecompositionUTC::storeOpenCutMesh(int part_id, 
		int part_pct, MeshPoint3d** part_points, 
		int part_fct, MeshFace** part_faces, bool with_cut_faces)
{
	ostringstream fname1,fname2;
	fname1 << m_cut_mesh_name << (with_cut_faces?"-open-part":"-fullopen-part") << part_id << "-p.txt";
	fname2 << m_cut_mesh_name << (with_cut_faces?"-open-part":"-fullopen-part") << part_id << "-i.txt";

	ofstream fpoints(fname1.str().c_str());
	ofstream felements(fname2.str().c_str());

	DataHashTableKeyValue<void*, int> kv_index_ref(part_pct, 0);

	fpoints << part_pct << endl;
	for(int i = 0; i < part_pct; i++){
		DPoint3d pt = part_points[i]->getCoordinates();
		kv_index_ref.insert(part_points[i], i);
		fpoints << i << "\t";
		fpoints << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
	}

	if(with_cut_faces){
		felements << part_fct << endl;
	}else{
		int ct = 0;
		for(int i = 0; i < part_fct; i++){
			if(part_faces[i]->getIntTag(TagExtended::TAG_UTC_CUT) >= 5) ct++;
		}
		felements << ct << endl;
	}
	for(int i = 0; i < part_fct; i++){
		if(!with_cut_faces && (part_faces[i]->getIntTag(TagExtended::TAG_UTC_CUT) < 5)) continue;
		felements << i;
		int edges_ct = part_faces[i]->getEdgeCount();
		felements << '\t' << edges_ct;
		for(int j = 0; j < edges_ct; j++){
			int ind = m_orientation_switch?(edges_ct-1-j):j;
			felements << '\t' << kv_index_ref.getValue(part_faces[i]->getPoint(ind), 0);
		}
		felements << '\t' << ((part_faces[i]->getIntTag(TagExtended::TAG_UTC_CUT) < 5)?2:1);
		felements << endl;
	}
}

void MeshDecompositionUTC::storeCutContours(DecompEdgeArrayList& contours)
{
	ostringstream fname1,fname2;
	fname1 << m_cut_mesh_name << "-contours-p.txt";
	fname2 << m_cut_mesh_name << "-contours-i.txt";

	ofstream fpoints(fname1.str().c_str());
	ofstream felements(fname2.str().c_str());

	int contour_count = contours.countInt();
	int pct = m_mesh->getPointsCount();
	DataVector<int> cut_points_tags(pct, -1);
	int cut_fct = 0;
	for(int i = 0; i < contour_count; i++){
		for(MeshEdge3d** edges = contours[i]; *edges; edges++){
			MeshEdge3d* edge = *edges;
			cut_points_tags[edge->getMeshPoint(0)->getIndex()] = 1;
			cut_points_tags[edge->getMeshPoint(1)->getIndex()] = 1;
			++cut_fct;
		}
	}
	int cut_pct = 0;
	for(int i = 0; i < pct; i++){
		if(cut_points_tags[i] > 0)
			cut_points_tags[i] = cut_pct++;
	}

	fpoints << cut_pct << endl;
	for(int i = 0; i < pct; i++){
		if(cut_points_tags[i] < 0) continue;
		DPoint3d pt = m_mesh->getPointAt(i)->getCoordinates();
		fpoints << cut_points_tags[i] << "\t";
		fpoints << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
	}

	felements << cut_fct << endl;
	int i_fct = 0;
	for(int i = 0; i < contour_count; i++){
		for(MeshEdge3d** edges = contours[i]; *edges; edges++){
			MeshEdge3d* edge = *edges;
			MeshPoint3d* mpt0 = edge->getMeshPoint(0);
			MeshPoint3d* mpt1 = edge->getMeshPoint(1);
			felements << i_fct++;
			felements << "\t3\t" << cut_points_tags[mpt0->getIndex()];
			felements << '\t' << cut_points_tags[mpt1->getIndex()];
			felements << '\t' << cut_points_tags[mpt1->getIndex()];
			felements << "\t1" << endl;
		}
	}
}

void MeshDecompositionUTC::storeCutFaces(DecompFaceList* cut_faces)
{
	ostringstream fname1,fname2;
	fname1 << m_cut_mesh_name << "-cutfaces-p.txt";
	fname2 << m_cut_mesh_name << "-cutfaces-i.txt";

	ofstream fpoints(fname1.str().c_str());
	ofstream felements(fname2.str().c_str());

	int pct = m_mesh->getPointsCount();
	DataVector<int> cut_points_tags(pct, -1);
	int cut_fct = cut_faces->countInt();
	for(int i = 0; i < cut_fct; i++){
		MeshFace* face = cut_faces->get(i);
		int ect = face->getEdgeCount();
		for(int j = 0; j < ect; j++){
			cut_points_tags[face->getPoint(j)->getIndex()] = 1;
		}
	}
	int cut_pct = 0;
	for(int i = 0; i < pct; i++){
		if(cut_points_tags[i] > 0)
			cut_points_tags[i] = cut_pct++;
	}

	fpoints << cut_pct << endl;
	for(int i = 0; i < pct; i++){
		if(cut_points_tags[i] < 0) continue;
		const DPoint3d& pt = m_mesh->getPointAt(i)->getCoordinates();
		fpoints << cut_points_tags[i] << "\t";
		fpoints << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
	}

	felements << cut_fct << endl;
	for(int i = 0; i < cut_fct; i++){
		felements << i;
		int edges_ct = cut_faces->get(i)->getEdgeCount();
		felements << '\t' << edges_ct;
		for(int j = 0; j < edges_ct; j++){
			int ind = m_orientation_switch?(edges_ct-1-j):j;
			felements << '\t' << cut_points_tags[cut_faces->get(i)->getPoint(ind)->getIndex()];
		}
		felements << '\t' << (1+cut_faces->get(i)->getIntTag(TagExtended::TAG_UTC_CUT));
		felements << endl;
	}
}

void MeshDecompositionUTC::storeClosedCutMesh(int part_id, 
		int part_pct, MeshPoint3d** part_points, int part_fct, MeshFace** part_faces, 
		int cut_pct, MeshPoint2d** cut_points, int cut_fct, MeshElement** cut_faces,
		SurfaceParametric* surface, bool orientation_valid)
{
	if(m_mesh->getControlSpace() && m_mesh->getControlSpace()->getType() == MeshData::CONTROL_OCTREE_3D){
		ostringstream fname_cs;
		fname_cs << m_cut_mesh_name << "-closed-part" << part_id << "-cs.txt";
		((ControlSpace3dOctree*)(m_mesh->getControlSpace().get()))->storeTXT(fname_cs.str().c_str());
	}

	ostringstream fname1,fname2;
	fname1 << m_cut_mesh_name << "-closed-part" << part_id << "-p.txt";
	fname2 << m_cut_mesh_name << "-closed-part" << part_id << "-i.txt";

	ofstream fpoints(fname1.str().c_str());
	ofstream felements(fname2.str().c_str());

	int boundary_count = 0;
	for(int i = 0; i < cut_pct; i++){
		if(cut_points[i]->availableTag(TagExtended::TAG_MP_2D_3D)){
			++boundary_count;
		}
	}

	DataHashTableKeyValue<void*, int> kv_index(part_pct+cut_pct, 0);

	fpoints << (part_pct+cut_pct-boundary_count) << endl;
	for(int i = 0; i < part_pct; i++){
		const DPoint3d& pt = part_points[i]->getCoordinates();
		kv_index.insert(part_points[i], i);
		fpoints << i << "\t";
		fpoints << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
	}
	for(int i = 0, j = part_pct; i < cut_pct; i++){
		if(cut_points[i]->availableTag(TagExtended::TAG_MP_2D_3D)) continue;
		const DPoint3d& pt = surface->getPoint(cut_points[i]->getCoordinates());
		kv_index.insert(cut_points[i], j);
		fpoints << j++ << "\t";
		fpoints << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
	}

	felements << (part_fct+cut_fct) << endl;
	for(int i = 0; i < part_fct; i++){
		felements << i;
		int edges_ct = part_faces[i]->getEdgeCount();
		felements << '\t' << edges_ct;
		for(int j = 0; j < edges_ct; j++){
			int ind = m_orientation_switch?(edges_ct-1-j):j;
			felements << '\t' << kv_index.getValue(part_faces[i]->getPoint(ind), 0);
		}
		felements << '\t' << ((part_faces[i]->getIntTag(TagExtended::TAG_UTC_CUT) < 5)?2:1);
		felements << endl;
	}
	for(int i = 0; i < cut_fct; i++){
		felements << (part_fct+i);
		int edges_ct = cut_faces[i]->getEdgeCount();
		felements << '\t' << edges_ct;
		bool any_mpt3d = false;
		for(int j = 0; j < edges_ct; j++){
			int ind = orientation_valid?j:(edges_ct-j-1);
			if(m_orientation_switch) ind = edges_ct-1-ind;
			MeshPoint2d* mpt = cut_faces[i]->getPoint(ind);
			MeshPoint3d* mpt3d = (MeshPoint3d*)mpt->getPtrTag(TagExtended::TAG_MP_2D_3D);
			if(mpt3d){
				felements << '\t' << kv_index.getValue(mpt3d, 0);
				any_mpt3d = true;
			}else{
				felements << '\t' << kv_index.getValue(mpt, 0);
			}
		}
		felements << '\t' << (any_mpt3d ? 3 : 4);
		felements << endl;
	}
}

/// !!! Stores additional edges !!! 
void MeshDecompositionUTC::storeInterfaceCutContours(int part_id, int cut_pct, MeshPoint2d** cut_points)
{
	ostringstream fname1,fname2;
	fname1 << m_cut_mesh_name << "-interface-contour-part" << part_id << "-p.txt";
	fname2 << m_cut_mesh_name << "-interface-contour-part" << part_id << "-i.txt";

	ofstream fpoints(fname1.str().c_str());
	ofstream felements(fname2.str().c_str());

	int store_pct = 0;
	for(int i = 0; i < cut_pct; i++){
		if(cut_points[i]->availableTag(TagExtended::TAG_MP_2D_3D))
			++store_pct;
	}
	DataHashTableKeyValue<void*, int> kv_index(cut_pct, 0);

	fpoints << store_pct << endl;
	for(int i = 0, j = 0; i < cut_pct; i++){
		if(!cut_points[i]->availableTag(TagExtended::TAG_MP_2D_3D)) continue;
		const DPoint2d& pt = cut_points[i]->getCoordinates();
		kv_index.insert(cut_points[i], j);
		fpoints << j++ << "\t";
		fpoints << pt.x << "\t" << pt.y << "\t" << '0' << endl;
	}

	int store_fct = 0;
	for(int i = 0; i < cut_pct; i++){
		if(!cut_points[i]->availableTag(TagExtended::TAG_MP_2D_3D)) continue;
		int rank = cut_points[i]->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = cut_points[i]->getEdge(j);
			if(edge->getPointIndex(cut_points[i]) != 0) continue;
			MeshPoint2d* mpt = edge->getOtherPoint(cut_points[i]);
			if(mpt->availableTag(TagExtended::TAG_MP_2D_3D))
				++store_fct;
		}
	}

	felements << store_fct << endl; 
	int i_fct = 0;
	for(int i = 0; i < cut_pct; i++){
		if(!cut_points[i]->availableTag(TagExtended::TAG_MP_2D_3D)) continue;
		int rank = cut_points[i]->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge2d* edge = cut_points[i]->getEdge(j);
			if(edge->getPointIndex(cut_points[i]) != 0) continue;
			MeshPoint2d* mpt = edge->getOtherPoint(cut_points[i]);
			if(!mpt->availableTag(TagExtended::TAG_MP_2D_3D)) continue;
			felements << i_fct++;
			felements << '\t' << 3;
			felements << '\t' << kv_index.getValue(cut_points[i], 0);
			felements << '\t' << kv_index.getValue(mpt, 0);
			felements << '\t' << kv_index.getValue(mpt, 0);
			felements << '\t' << 1 << endl;
		}
	}
}

void MeshDecompositionUTC::storeInterfaceCutMesh(int part_id, int cut_pct, MeshPoint2d** cut_points, 
		int cut_fct, MeshElement** cut_faces, SurfaceParametric* surface)
{
	ostringstream fname1,fname2;
	fname1 << m_cut_mesh_name << "-interface" << (surface?"3D":"2D") << "-part" << part_id << "-p.txt";
	fname2 << m_cut_mesh_name << "-interface" << (surface?"3D":"2D") << "-part" << part_id << "-i.txt";

	ofstream fpoints(fname1.str().c_str());
	ofstream felements(fname2.str().c_str());

	DataHashTableKeyValue<void*, int> kv_index(cut_pct, 0);

	fpoints << cut_pct << endl;
	for(int i = 0; i < cut_pct; i++){
		if(surface){
			MeshPoint3d* mpt3d = (MeshPoint3d*)cut_points[i]->getPtrTag(TagExtended::TAG_MP_2D_3D);
			const DPoint3d pt = mpt3d ? mpt3d->getCoordinates() : surface->getPoint(cut_points[i]->getCoordinates());
			kv_index.insert(mpt3d ? (void*)mpt3d : (void*)cut_points[i], i);
			fpoints << i << "\t";
			fpoints << pt.x << "\t" << pt.y << "\t" << pt.z << endl;
		}else{
			const DPoint2d& pt = cut_points[i]->getCoordinates();
			kv_index.insert(cut_points[i], i);
			fpoints << i << "\t";
			fpoints << pt.x << "\t" << pt.y << "\t" << '0' << endl;
		}
	}

	felements << cut_fct << endl;
	for(int i = 0; i < cut_fct; i++){
		felements << i;
		int edges_ct = cut_faces[i]->getEdgeCount();
		felements << '\t' << edges_ct;
		bool any_mpt3d = false;
		for(int j = 0; j < edges_ct; j++){
			int ind = m_orientation_switch?(edges_ct-j-1):j;
			MeshPoint2d* mpt = cut_faces[i]->getPoint(ind);
			MeshPoint3d* mpt3d = (MeshPoint3d*)mpt->getPtrTag(TagExtended::TAG_MP_2D_3D);
			if(surface && mpt3d){
				felements << '\t' << kv_index.getValue(mpt3d, 0);
				any_mpt3d = true;
			}else{
				felements << '\t' << kv_index.getValue(mpt, 0);
			}
		}
		felements << '\t' << (any_mpt3d ? 3 : 4);
		felements << endl;
	}
}

/*
bool MeshDecompositionUTC::storeCutSurfaceMeshes(MeshContainer2d* cut_mesh)
{
	MeshBlock* block = m_mesh->getBlockAt(0);
	assert(block);
	int fct = block->getFaceCount();
	for(int i = 0; i < fct; i++)
		block->getFace(i)->setBorder();
	int pct = block->getPointCount();
	for(int i = 0; i < pct; i++)
		block->getPoint(i)->removeTag(TagExtended::TAG_UTC_CUT);
	int cut_pct = cut_mesh->getPointsCount();
	for(int i = 0; i < cut_pct; i++){
		MeshPoint2d* mpt = cut_mesh->getPointAt(i);
		mpt->removeTag(TagExtended::TAG_UTC_CUT);
	}
	int cut_fct = cut_mesh->getElementsCount();
	for(int i = 0; i < cut_fct; i++)
		cut_mesh->getElementAt(i)->setAreaID(0);

	m_log_file << "=== CUT MESHES ===" << endl;
	
	int part_id = 1;
	int visited_fct = 0;
	MeshFace** part_faces = new MeshFace*[fct];
	MeshPoint3d** part_points = new MeshPoint3d*[pct];
	MeshElement** cut_faces = new MeshElement*[cut_fct];
	MeshPoint2d** cut_points = new MeshPoint2d*[cut_pct];
	while(visited_fct < fct){
		for(int i = 0; i < fct; i++){
			part_faces[0] = block->getFace(i);
			if(part_faces[0]->isBorder())	break;
		}
		assert(part_faces[0]);
		part_faces[0]->setBorder(); // differene way of marking cut-interface is necassary
		int part_fct = 1;
		int part_pct = 0;
		for(int i = 0; i < part_fct; i++){
			int rank = part_faces[i]->getEdgeCount();
			for(int j = 0; j < rank; j++){
				MeshEdge3d* edge = part_faces[i]->getEdge(j);
				MeshPoint3d* mpt = part_faces[i]->getPoint(j);
				if(mpt->getIntTag(TagExtended::TAG_UTC_CUT) < part_id){
					part_points[part_pct++] = mpt;
					mpt->setIntTag(TagExtended::TAG_UTC_CUT, part_id);
				}
				if(edge->availableTag(TagExtended::TAG_UTC_CUT)) continue; // don't cross the cut-line
				int e_rank = edge->getFaceCount();
				for(int k = 0; k < e_rank; k++){
					MeshFace* face = edge->getFaceAt(k);
					if(face->getBorderType() == 0){
						face->setBorderType(part_id);
						part_faces[part_fct++] = face;
					}
				}
			}
		}
		m_log_file << "PART " << part_id << ": " << part_pct << " oryginal nodes, ";
		m_log_file << part_fct << " oryginal faces." << endl;

		storeOpenCutMesh(part_id, part_pct, part_points, part_fct, part_faces, true);
		storeOpenCutMesh(part_id, part_pct, part_points, part_fct, part_faces, false);

		// check 2d cut mesh
		int cut_part_fct = 0;
		int cut_part_pct = 0;
		bool orientation_valid = true;
		for(int i = 0; i < cut_pct; i++){
			MeshPoint2d* cut_mp = cut_mesh->getPointAt(i);
			if(cut_mp->getIntTag(TagExtended::TAG_UTC_CUT) == part_id) continue;
			MeshPoint3d* cut_mp3d = (MeshPoint3d*)cut_mp->getPtrTag(TagExtended::TAG_MP_2D_3D);
			if(!cut_mp3d || (cut_mp3d->getIntTag(TagExtended::TAG_UTC_CUT) < part_id)) continue;
			MeshEdge2d* first_edge = cut_mp->getEdge(0);
			cut_faces[cut_part_fct] = first_edge->getMeshElement(0);
			if(!cut_faces[cut_part_fct]) cut_faces[cut_part_fct] = first_edge->getMeshElement(1);
			assert(cut_faces[cut_part_fct]);
			cut_faces[cut_part_fct]->setAreaID(part_id);
			++cut_part_fct;
			// set orientation (triangle -> points3d -> edge -> other triangle with current part_id)
			if(cut_part_fct == 1){ // only for first triangle
				int first_rank = cut_mp->getRank();
				for(int k = 0; k < first_rank; k++){
					first_edge = cut_mp->getEdge(k);
					MeshPoint2d* other_mp = first_edge->getOtherPoint(cut_mp);
					MeshPoint3d* other_mp3d = (MeshPoint3d*)other_mp->getPtrTag(TagExtended::TAG_MP_2D_3D);
					if(!other_mp3d) continue;
					MeshEdge3d* edge3d = cut_mp3d->getEdgeToPoint(other_mp3d);
					if(!edge3d) continue;
					MeshFace* face = 0;
					for(int l = 0; l < edge3d->getFaceCount(); l++)
						if((face = edge3d->getFaceAt(l))->getBorderType() == part_id) break;
					bool face_orientation = face->properOrientation(cut_mp3d, other_mp3d);
					bool element_orientation = cut_faces[0]->properOrientation(cut_mp, other_mp);
					if(face_orientation == element_orientation) 
						orientation_valid = false; // for neighbouring elements direction of two nodes should differ
					break;
				}
			}
			// follow triangles incident to this node
			int local_fct = 1;
			int local_pct = 0;
			for(int k = cut_part_fct-1; k < cut_part_fct; k++){
				int rank = cut_faces[k]->getEdgeCount();
				for(int j = 0; j < rank; j++){
					MeshEdge2d* edge = cut_faces[k]->getEdge(j);
					MeshPoint2d* mpt = cut_faces[k]->getPoint(j);
					if(mpt->getIntTag(TagExtended::TAG_UTC_CUT) < part_id){
						++local_pct;
						mpt->setIntTag(TagExtended::TAG_UTC_CUT, part_id);
//						if(!mpt->getPoint3dLink()) 
						cut_points[cut_part_pct++] = mpt;
					}
					for(int l = 0; l < 2; l++){
						MeshElement* face = edge->getMeshElement(l);
						if(!face || face->getAreaID() == part_id) continue;
						face->setAreaID(part_id);
						cut_faces[cut_part_fct++] = face;
						++local_fct;
					}
				}
			}
			m_log_file << " * CUT with: " << local_pct << "  nodes (including boundary), ";
			m_log_file << local_fct << " faces." << endl;
		}

		storeClosedCutMesh(part_id, part_pct, part_points, part_fct, part_faces,
			cut_part_pct, cut_points, cut_part_fct, cut_faces, 
			cut_mesh->getSurface(), orientation_valid);

		if(orientation_valid){	// store interface only once
			storeInterfaceCutMesh(part_id, 	cut_part_pct, cut_points,	// in 3D
				cut_part_fct, cut_faces, cut_mesh->getSurface());
			storeInterfaceCutMesh(part_id, 	cut_part_pct, cut_points,	// in 2D
				cut_part_fct, cut_faces, 0);
			storeInterfaceCutContours(part_id, cut_part_pct, cut_points);
		}

		visited_fct += part_fct;
		++part_id;
	}

	delete[] cut_points;
	delete[] cut_faces;
	delete[] part_points;
	delete[] part_faces;
	return true;
}
*/

bool MeshDecompositionUTC::run()
{
	if(!m_mesh) return false;
	START_CLOCK("cut mesh");
	markPointsOrientation();
	DecompFaceList* cut_faces = selectCutTriangles();
	DecompNodeList* node_list = selectPoints(cut_faces);
	DecompEdgeArrayList contours = findContour3D(node_list, cut_faces);

	// ********** show
	if(m_visualization > 0){
		SHOW_MESH("cut mesh", getViewSet());
	}

	storeCutContours(contours);
	MeshContainer2d* cut_mesh = createCutMesh(contours);
	STOP_CLOCK("cut mesh");
//	storeCutSurfaceMeshes(cut_mesh);


	delete cut_faces;
	delete node_list;

	return true;
}
