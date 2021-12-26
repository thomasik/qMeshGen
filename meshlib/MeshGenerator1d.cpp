// MeshGenerator1d.cpp: implementation of the MeshGenerator1d class.
//
//////////////////////////////////////////////////////////////////////

#include "MeshGenerator1d.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
#include "MeshPoint2d.h"
#include "MeshPoint3d.h"
#include "MeshEdge2d.h"
#include "MeshDomainEdge3d.h"
#include "SurfaceParametric.h"
#include "MeshDomainSurface.h"
#include "MeshDomainVolume.h"
#include "MeshArea.h"
#include "ControlSpace2d.h"
#include "ControlSpace2dAdaptive.h"
#include "ControlSpace3d.h"
#include "Metric3dContext.h"

int MeshGenerator1d::param_even_node_count = 0;

/// Whether to smoothen discretization of edge (for last segment, which may be smaller than required)
int MeshGenerator1d::param_smoothen_last_node = 5;

/// Maximum difference for metric-length at edge vertices during discretization
double MeshGenerator1d::param_boundary_metric_conform_rato = 1.8;

// Distributes vertexes along boundaries, conforming to the control space
int MeshGenerator1d::discretizeEdgesMin(MeshContainer3d *boundary)
{
	int pcount = boundary->getPointsCount();
	if(pcount < 1) return 0;
	int total_points = 0;


	START_CLOCK("MG1d::discretizeEdges");

	boundary->clearDiscretization();

	for(int i = 0; i < pcount; i++){
		MeshPoint3d* point3d = boundary->getPointAt(i);
		int rank = point3d->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge3d = point3d->getEdge(j);
			if(edge3d->getPointIndex(point3d) == 0){
				if(edge3d->getType() != EDGE_DOMAIN_3D)
					LOG4CPLUS_ERROR(MeshLog::logger_console,   "Gen1d: discretization of non-domain edge?");
				else
					total_points += discretizeEdgeMin((MeshDomainEdge3d*)edge3d);
			}
		}
	}

	STOP_CLOCK("MG1d::discretizeEdges");

	return total_points;
}

// Distributes vertexes along boundaries, conforming to the control space
int MeshGenerator1d::discretizeEdgeMin(MeshDomainEdge3d* edge3d)
{
	double max_metric_confrom_ratio = param_boundary_metric_conform_rato;
	double max_metric_confrom_ratio2 = max_metric_confrom_ratio*max_metric_confrom_ratio;

	auto& disc_points = edge3d->getDiscretization();
	disc_points.clear();
	edge3d->setValidDiscretization(true); // actually, it will be valid soon

	// create 
	int face_count = edge3d->getFaceCount();
	assert(face_count > 0);
	int selected_face = 0;
	double selected_length = -1.0;
	DataVector<CS2dPtr> control_spaces(face_count);
	DataVector<MeshEdge2d*> mesh_edges(face_count);
	DataVector<SurfaceConstPtr> surfaces(face_count);
	// select one
	for(int k = 0; k < face_count; k++){
		MeshDomainSurface* face = (MeshDomainSurface*)edge3d->getFaceAt(k);
		assert(face->getType() == FACE_DOMAIN);
		SurfaceConstPtr surface = face->getBaseSurface();
		surfaces.add(surface);
		MeshContainer2d* mesh = face->getBoundary();
		CS2dPtr control = mesh->getControlSpace();
		if(!control){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "No control, creating new one...");
			// Prepare control space (introductory factors)
			MeshDomainVolume* dvolume0 = (MeshDomainVolume*)face->getBlock(0);
			MeshDomainVolume* dvolume1 = (MeshDomainVolume*)face->getBlock(1);
			CS3dPtr cs_0, cs_1;
			if(dvolume0){
				cs_0 = dvolume0->getControlSpace();
				if(!cs_0) cs_0 = dvolume0->getUserControlSpace();
			}
			if(dvolume1 && (dvolume0 != dvolume1)){
				cs_1 = dvolume1->getControlSpace();
				if(!cs_1) cs_1 = dvolume1->getUserControlSpace();
			}
			if(face->isBoundedBothSides()){ // for inner boundary adjust sizing
				ControlSpace2dAdaptive::param_max_diameter_ratio *= ControlSpace2dAdaptive::param_inner_boundary_ratio;
			}
			mesh->createControlSpace(face->getUserControlSpace(), cs_0, cs_1);
			if(face->isBoundedBothSides()){ // for inner boundary revert sizing adjustment
				ControlSpace2dAdaptive::param_max_diameter_ratio /= ControlSpace2dAdaptive::param_inner_boundary_ratio;
			}
			control = mesh->getControlSpace();
		}
		assert(control);
		control_spaces.add(control);
		MeshEdge2d* mesh_edge = mesh->getCoupledEdge(edge3d);
		mesh_edges.add(mesh_edge);
		assert(mesh_edge);
		Metric2dContext mc(control);
		double edge_length = mesh_edge->getLengthMetricAdapted(mc);
//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "xxx edge_length = " << edge_length);
		if(edge_length > selected_length){
			selected_length = edge_length;
			selected_face = k;
		}
	}

	LOG4CPLUS_TRACE(MeshLog::logger_mesh, 
		face_count << " incident faces, selected: " << selected_face
		 << " with metric length: " << selected_length);
    MeshDomainSurface* face = (MeshDomainSurface*)edge3d->getFaceAt(selected_face);
	assert(face->getType() == FACE_DOMAIN);
	SurfaceConstPtr surface = face->getBaseSurface();
	MeshContainer2d* mesh = face->getBoundary();
	Metric2dContext mc(mesh->getControlSpace());
	MeshEdge2d* mesh_edge = mesh_edges[selected_face];

	auto freepoints = edge3d->getFreePoints();
	int fcount = freepoints ? (int)freepoints->countInt() : 0;
	if(selected_length < 1.5 && (fcount == 0)){ // small enough already
		if(param_even_node_count){
			DPoint3d pt = surface->getPoint(mesh_edge->getPoint(0.5));
			disc_points.add(std::make_shared<MeshPoint3d>(pt));
			disc_points[0]->setDoubleTag(TagExtended::TAG_DISCRETIZATION_PARAM, 0.5);
			return 1;
		}else
			return 0;
	}
	// else, insert some points ...
	DataVector<double> fnodes(fcount+1);

	// first, insert the "free" points
	for(int i = 0; i < fcount; i++){
		DPoint2d param = surface->getParameters(freepoints->get(i)->getCoordinates());
		double ksi = mesh_edge->getParameter(param);
		fnodes.add(ksi);
	}
	for(int i = 1; i < fcount; i++){ // just in case... -> sort
		int i_start = i-1;
		int i_min = i_start;
		for(int j = i; j < fcount; j++)
			if(fnodes[j] < fnodes[i_min]) i_min = j;
		if(i_min != i_start){
			fnodes.switchData(i_start, i_min);
			freepoints->switchData(i_start, i_min);
		}
	}

	const double PRECISION = 1e-10;
	for(int i = 1; i < fcount; ){ // ... and remove duplicates
		if(abs(fnodes[i] - fnodes[i-1]) < PRECISION){
			fnodes.removeOrderedAt(i);
			freepoints->removeOrderedAt(i);
			--fcount;
		}else
			++i;
	}
	// + 1.0
	fnodes.add(1.0);

	DataVector<double> nodes((int)(2*selected_length) + fcount); // length ~ number of discretized edges

	double t_begin = 0.0;
	double t_end = 0.0;

	for(int ti = 0; ti < (int)fnodes.countInt(); ti++){
		t_begin = t_end;
		t_end = fnodes[ti];
		int i_begin = (int)nodes.countInt();

		// (possibly) discretize segment t_begin - t_end
		double last_t = t_begin;
		double next_t = t_begin;
		double edge_length = (fcount == 0) ? selected_length : mesh_edge->getLengthMetricAdapted(mc, t_begin, t_end);
		double dt = (t_end - t_begin) / edge_length;
		double max_len = mesh_edge->getLengthMax(mc, t_begin, t_begin+dt);
		if(max_len > 1.0) dt /= max_len;
		double required_length = 1.0;
		// count all segments
		while(true){
			// count next segment
			double local_last_t = last_t;
			while(required_length > PRECISION){ 
				next_t = std::min(local_last_t + dt, t_end);
				double real_length = mesh_edge->getRequiredLength(mc, local_last_t, next_t, required_length);
				//double real_length_1 = mesh_edge->getLengthForVariableMetric(local_last_t, next_t, 1e-3);
				assert(required_length >= real_length);
				required_length -= real_length;
				if(t_end - next_t < PRECISION) break; // reached end of the segment
				if(required_length > PRECISION){ // next bit of segment required
					dt = (next_t-local_last_t) * required_length / real_length;
					local_last_t = next_t;
				}
			}
			if(t_end - next_t > PRECISION){ // it's not the end of whole edge
				nodes.add(next_t);
				if (nodes.countInt() >= MeshPoint2d::param_max_node_count) {
					LOG4CPLUS_ERROR(MeshLog::logger_console, "Maximum number of mesh nodes exceeded.");
					return 0;
				}
				dt = next_t - last_t;
				last_t = next_t;
				required_length = 1.0;
			}else break;
		}
		// Adjust ending length proportionally
		double real_length = 1.0 - required_length;
		int nodes_ct = (int)nodes.countInt() - i_begin;
		if(real_length >= 1.0 || (param_smoothen_last_node == 0) || (nodes_ct == 0)){
			// nothing
		}else if(real_length >= 0.5){
			double t0 = t_begin;
			int first_node = 0;
			
			if(param_smoothen_last_node > 0){
				first_node = std::max(0, nodes_ct-param_smoothen_last_node);
				if(first_node > 0) t0 = nodes[i_begin+first_node-1];
			}

			DPoint3d last_pt = surface->getPoint(mesh_edge->getPoint(t0));
			double total_len = 0.0;
			for(int i = first_node; i < nodes_ct; i++){
				DPoint3d pt = surface->getPoint(mesh_edge->getPoint(nodes[i_begin+i]));
				total_len += last_pt.distance(pt);
				last_pt = pt;
			}
			double rest_len = last_pt.distance(surface->getPoint(mesh_edge->getPoint(t_end)));
			double rest_len_x = rest_len * (1.0-real_length) / real_length;
			double ratio = (total_len + rest_len) / (total_len + rest_len + rest_len_x);
			/* Old method
			double dt = (1.0 - real_length) / (nodes_ct+1);
			double ratio = 1.0 / (1.0 + dt);
			*/
			double *dts = new double[nodes_ct];
			double temp_t = t0;
			for(int i = first_node; i < nodes_ct; i++){
				dts[i] = ratio*(nodes[i_begin+i]-temp_t);
				temp_t = nodes[i_begin+i];
			}
			temp_t = t0;
			for(int i = first_node; i < nodes_ct; i++)
				temp_t = nodes[i_begin+i] = temp_t + dts[i];
			delete[] dts;
		}else{
			int nodes_ct = (int)nodes.countInt() - i_begin;
			double t0 = t_begin;
			int first_node = 0;
			
			if(param_smoothen_last_node > 0){
				if (nodes_ct > param_smoothen_last_node)
					first_node = nodes_ct - param_smoothen_last_node;
				if(first_node > 0) t0 = nodes[i_begin+first_node-1];
			}

			DPoint3d last_pt = surface->getPoint(mesh_edge->getPoint(t0));
			double total_len = 0.0;
			for(int i = first_node; i < nodes_ct; i++){
				DPoint3d pt = surface->getPoint(mesh_edge->getPoint(nodes[i_begin+i]));
				total_len += last_pt.distance(pt);
				last_pt = pt;
			}
			double rest_len = last_pt.distance(surface->getPoint(mesh_edge->getPoint(t_end)));
			double ratio = (total_len + rest_len) / total_len;

			/* Old method
			double ratio = 1.0 + real_length / nodes_ct;
			*/
			double *dts = new double[nodes_ct];
			double temp_t = t0;
			for(int i = first_node; i < nodes_ct-1; i++){
				dts[i] = ratio*(nodes[i_begin+i]-temp_t);
				temp_t = nodes[i_begin+i];
			}
			nodes.removeLast();
			temp_t = t0;
			for(int i = first_node; i < nodes_ct-1; i++)
				temp_t = nodes[i_begin+i] = temp_t + dts[i];
			if(nodes_ct > 1)
				while(nodes[nodes.countInt()-1] >= t_end)
					nodes.removeLast();
		}
		// add "t_end" temporarily
		nodes.add(t_end);
	}

	int node_count = nodes.countInt();
	if(node_count > 1){
		assert(nodes[0] > 0.0);
		assert(nodes[node_count-2] < 1.0);
	}

	// Check conformity of control space for incident faces
	if(face_count > 1){
		for(int k = 0; k < face_count; k++){
			if(k == selected_face) continue;
			MeshEdge2d* mesh_edge_k = mesh_edges[k];
			Metric2dContext incident_mc(control_spaces[k]);
			SurfaceConstPtr surface_k = surfaces[k];
			bool same_dir = (
				mesh_edge_k->getMeshPoint(0)->getPtrTag(TagExtended::TAG_MP_2D_3D) ==
				mesh_edge->getMeshPoint(0)->getPtrTag(TagExtended::TAG_MP_2D_3D));
			DPoint2d pt_0 = mesh_edge_k->getMeshPoint(same_dir?0:1)->getCoordinates();
			double ksi = same_dir ? 0.0 : 1.0;
			double min_ksi = 0.0;
			for(int i = 0; i < node_count; ){
				const DPoint3d pt3 = surface->getPoint(mesh_edge->getPoint(nodes[i]));
				const DPoint2d pt_1 = mesh_edge_k->surfaceParameters(surface_k, pt3, ksi, same_dir);
				if(nodes[i] - ((i>0)?nodes[i-1]:min_ksi) > 1e-3){///???????????????????
					incident_mc.countMetricAtPoint(DPoint2d::average(pt_0, pt_1));
					const DVector2d dpt = pt_1-pt_0;
					double dist = incident_mc.transformPStoMS(dpt).length2();
					if(dist > max_metric_confrom_ratio2){
						incident_mc.countMetricAtPoint(pt_0);
						double dist_0 = incident_mc.transformPStoMS(dpt).length2();
						incident_mc.countMetricAtPoint(pt_1);
						double dist_1 = incident_mc.transformPStoMS(dpt).length2();
						double dist_min = std::min(dist, std::min(dist_0, dist_1));
						double dist_max = std::max(dist, std::max(dist_0, dist_1));
						if(dist_max / dist_min < 2.0){
							LOG4CPLUS_WARN(MeshLog::logger_console, 
								"too large incident control-space nonconformity.");
							double new_node = 0.5*(nodes[i] + ((i>0)?nodes[i-1]:min_ksi));
							nodes.insertAt(i, new_node);
							++node_count;
						}else{
							LOG4CPLUS_WARN(MeshLog::logger_console, 
								"too large metric fluctuations.");
							pt_0 = pt_1;
							++i;
						}
					}else{
						pt_0 = pt_1;
						++i;
					}
				}else{
					pt_0 = pt_1;
					++i;
				}
			}
		}
	}
	if(node_count > 1){
		assert(nodes[0] > 0.0);
		assert(nodes[node_count-2] < 1.0);
	}

	node_count = (int)nodes.countInt();
	// Check parity (if required) -- for now, not used with "free-points" !!!
	if(param_even_node_count && (fcount == 0)){
		if(node_count == 1){ //(with the additional "1.0")
			nodes.insertAt(0, 0.5);
			++node_count;
		}else if(node_count % 2 == 1){ //(with the additional "1.0")
			// Insert an additional node (or delete?)
			int mid_node = node_count/2;
			int first_node = 0;
			int last_node = node_count - 1;
			int nct = node_count;
			if(param_smoothen_last_node > 0){
				nct = std::min(param_smoothen_last_node, nct);
				first_node = std::max(first_node, mid_node-nct);
				last_node  = std::min(last_node,  mid_node+nct);
			}
			double *dts = new double[last_node-first_node];
			double first_t = (first_node > 0)?nodes[0]:0.0;
			double temp_t = first_t;
			double short_len = 0.0;
			for(int i = first_node, j=0; i <= last_node; i++){
				if(i != mid_node) short_len += (dts[j++] = nodes[i]-temp_t);
				temp_t = nodes[i];
			}
			double ratio = (nodes[last_node] - first_t)/ short_len;
			nodes.removeOrderedAt(mid_node);
			temp_t = first_t;
			for(int i = first_node; i < last_node; i++){
				temp_t = nodes[i] = temp_t + ratio * dts[i-first_node];
			}
			--node_count;
			delete[] dts;
		}
	}
	if(node_count > 1){
		assert(nodes[0] > 0.0);
		assert(nodes[node_count-2] < 1.0);
	}

	nodes.removeLast(); // remove the obsolete "1.0"
	disc_points.prepare(node_count);
	bool same_dir = (mesh_edge->getMeshPoint(0)->getPtrTag(TagExtended::TAG_MP_2D_3D) ==
		edge3d->getMeshPoint(0));
	MeshPoint3d* last_point = edge3d->getMeshPoint(0);
	--node_count;
//		int ict = 0;
	for(int i = 0; i < node_count; i++){
		int inode = same_dir?i:(node_count-i-1);
		double ksi = nodes[inode];
		DPoint3d pt = surface->getPoint(mesh_edge->getPoint(ksi));
		double dist = last_point->getCoordinates().distance2(pt);
		if(dist > mesh_data.relative_small_number){
			auto dp = std::make_shared<MeshPoint3d>(pt);
			disc_points.add(dp);
			last_point = dp.get();
			last_point->setDoubleTag(TagExtended::TAG_DISCRETIZATION_PARAM, ksi);
			for(int j = 0; j < fcount; j++){
				if(abs(ksi - fnodes[j]) < PRECISION){
					void* mbc = freepoints->get(j)->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
					if(mbc) last_point->setPtrTag(TagExtended::TAG_BOUNDARY_COND, mbc);
					break;
				}
			}
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Boundary nodes too close");
		}
	}

	LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "MG1d::disretizeEdgeMin (" << disc_points.countInt() << ") done.");

	return disc_points.countInt();
}

// --- old version ---
//// Distributes vertexes along boundaries, conforming to the control space
//int MeshGenerator1d::discretizeEdgeMin(MeshDomainEdge3d* edge3d)
//{
//	double max_metric_confrom_ratio = param_boundary_metric_conform_rato;
//	double max_metric_confrom_ratio2 = max_metric_confrom_ratio*max_metric_confrom_ratio;
//
//	DataPtrVector<MeshPoint3d> & disc_points = edge3d->getDiscretization();
//	disc_points.clear();
//	edge3d->setValidDiscretization(true); // actually, it will be valid soon
//
//	// create 
//	int face_count = edge3d->getFaceCount();
//	assert(face_count > 0);
//	int selected_face = 0;
//	double selected_length = -1.0;
//	DataVector<CS2dPtr> control_spaces(face_count);
//	DataVector<MeshEdge2d*> mesh_edges(face_count);
//	DataVector<SurfaceParametric*> surfaces(face_count);
//	// select one
//	for(int k = 0; k < face_count; k++){
//		MeshDomainSurface* face = (MeshDomainSurface*)edge3d->getFaceAt(k);
//		assert(face->getType() == FACE_DOMAIN);
//		SurfaceParametric* surface = face->getBaseSurface();
//		surfaces.add(surface);
//		MeshContainer2d* mesh = face->getBoundary();
//		CS2dPtr control = mesh->getControlSpace();
//		if(!control){
//			// Prepare control space (introductory factors)
//			MeshDomainVolume* dvolume0 = (MeshDomainVolume*)face->getBlock(0);
//			MeshDomainVolume* dvolume1 = (MeshDomainVolume*)face->getBlock(1);
//			CS3dPtr ucs_0 = 
//				dvolume0 ? dvolume0->getUserControlSpace() : nullptr;
//			CS3dPtr ucs_1 = 
//				(dvolume1 && (dvolume0 != dvolume1)) ? dvolume1->getUserControlSpace() : nullptr;
//			mesh->createControlSpace(face->getUserControlSpace(), ucs_0, ucs_1);
//			control = mesh->getControlSpace();
//		}
//		assert(control);
//		control_spaces.add(control);
//		MeshEdge2d* mesh_edge = mesh->getCoupledEdge(edge3d);
//		mesh_edges.add(mesh_edge);
//		assert(mesh_edge);
//		Metric2dContext mc(control);
//		double edge_length = mesh_edge->getLengthMetricAdapted(mc);
//		if(edge_length > selected_length){
//			selected_length = edge_length;
//			selected_face = k;
//		}
//	}
////	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, face_count << " incident faces, selected: " << selected_face;
////	LOG4CPLUS_INFO(MeshLog::logger_mesh, " with metric length: " << selected_length);
//    MeshDomainSurface* face = (MeshDomainSurface*)edge3d->getFaceAt(selected_face);
//	assert(face->getType() == FACE_DOMAIN);
//	SurfaceParametric* surface = face->getBaseSurface();
//	MeshContainer2d* mesh = face->getBoundary();
//	Metric2dContext mc(mesh->getControlSpace());
//	MeshEdge2d* mesh_edge = mesh_edges[selected_face];
//	if(selected_length > 1.5){ // this length ~ number of discretized edges
//		DataVector<double> nodes((int)(2*selected_length));
//		double last_t = 0.0;
//		double next_t = 0.0;
//		double dt = 1.0 / selected_length;
//		double max_len = mesh_edge->getLengthMax(mc, 0.0, dt);
//		if(max_len > 1.0) dt /= max_len;
//		double required_length = 1.0;
//		// count all segments
//		const double PRECISION = 1e-10;
//		while(true){
//			// count next segment
//			double local_last_t = last_t;
//			while(required_length > PRECISION){ 
//				next_t = std::min(local_last_t + dt, 1.0);
//				double real_length = mesh_edge->getRequiredLength(mc, local_last_t, next_t, required_length);
//				//double real_length_1 = mesh_edge->getLengthForVariableMetric(local_last_t, next_t, 1e-3);
//				assert(required_length >= real_length);
//				required_length -= real_length;
//				if(1.0 - next_t < PRECISION) break; // reached end of the segment
//				if(required_length > PRECISION){ // next bit of segment required
//					dt = (next_t-local_last_t) * required_length / real_length;
//					local_last_t = next_t;
//				}
//			}
//			if(1.0 - next_t > PRECISION){ // it's not the end of whole edge
//				nodes.add(next_t);
//				LOG4CPLUS_ASSERT(MeshLog::logger_mesh, (nodes.countInt() < MeshPoint2d::param_max_node_count),
//					MeshingException("Maximum number of mesh nodes exceeded."), 0);
//				dt = next_t - last_t;
//				last_t = next_t;
//				required_length = 1.0;
//			}else break;
//		}
//		// Adjust ending length proportionally
//		double real_length = 1.0 - required_length;
//		if(real_length >= 1.0 || (param_smoothen_last_node == 0)){
//			// nothing
//		}else if(real_length >= 0.5){
//			int nodes_ct = nodes.countInt();
//			double t0 = 0.0;
//			int first_node = 0;
//			
//			if(param_smoothen_last_node > 0){
//				first_node = std::max(0, nodes_ct-param_smoothen_last_node);
//				if(first_node > 0) t0 = nodes[first_node-1];
//			}
//
//			DPoint3d last_pt = surface->getPoint(mesh_edge->getPoint(t0));
//			double total_len = 0.0;
//			for(int i = first_node; i < nodes_ct; i++){
//				DPoint3d pt = surface->getPoint(mesh_edge->getPoint(nodes[i]));
//				total_len += last_pt.distance(pt);
//				last_pt = pt;
//			}
//			double rest_len = last_pt.distance(surface->getPoint(mesh_edge->getPoint(1.0)));
//			double rest_len_x = rest_len * (1.0-real_length) / real_length;
//			double ratio = (total_len + rest_len) / (total_len + rest_len + rest_len_x);
//			/* Old method
//			double dt = (1.0 - real_length) / (nodes_ct+1);
//			double ratio = 1.0 / (1.0 + dt);
//			*/
//			double *dts = new double[nodes_ct];
//			double temp_t = t0;
//			for(int i = first_node; i < nodes_ct; i++){
//				dts[i] = ratio*(nodes[i]-temp_t);
//				temp_t = nodes[i];
//			}
//			temp_t = t0;
//			for(int i = first_node; i < nodes_ct; i++)
//				temp_t = nodes[i] = temp_t + dts[i];
//			delete[] dts;
//		}else{
//			int nodes_ct = nodes.countInt();
//			double t0 = 0.0;
//			int first_node = 0;
//			
//			if(param_smoothen_last_node > 0){
//				first_node = std::max(0, nodes_ct-param_smoothen_last_node);
//				if(first_node > 0) t0 = nodes[first_node-1];
//			}
//
//			DPoint3d last_pt = surface->getPoint(mesh_edge->getPoint(t0));
//			double total_len = 0.0;
//			for(int i = first_node; i < nodes_ct; i++){
//				DPoint3d pt = surface->getPoint(mesh_edge->getPoint(nodes[i]));
//				total_len += last_pt.distance(pt);
//				last_pt = pt;
//			}
//			double rest_len = last_pt.distance(surface->getPoint(mesh_edge->getPoint(1.0)));
//			double ratio = (total_len + rest_len) / total_len;
//
//			/* Old method
//			double ratio = 1.0 + real_length / nodes_ct;
//			*/
//			double *dts = new double[nodes_ct];
//			double temp_t = t0;
//			for(int i = first_node; i < nodes_ct-1; i++){
//				dts[i] = ratio*(nodes[i]-temp_t);
//				temp_t = nodes[i];
//			}
//			nodes.removeLast();
//			temp_t = t0;
//			for(int i = first_node; i < nodes_ct-1; i++)
//				temp_t = nodes[i] = temp_t + dts[i];
//			if(nodes_ct > 1)
//				while(nodes[nodes.countInt()-1] >= 1.0)
//					nodes.removeLast();
//		}
//		// add "1.0" temporarily
//		nodes.add(1.0);
//		int node_count = nodes.countInt();
//		if(node_count > 1){
//			assert(nodes[0] > 0.0);
//			assert(nodes[node_count-2] < 1.0);
//		}
//
///*
//		// Check max ratio of lenghts for neighbouring edges
//		DPoint3d last_pt = mesh_edge->getMeshPoint(0)->getPoint3dLink()->getCoordinates();
//		DPoint3d mid_pt = surface->getPoint(mesh_edge->getPoint(nodes[0]));					
//		for(int i = 1; i < node_count; ){
//			DPoint3d next_pt = surface->getPoint(mesh_edge->getPoint(nodes[i]));
//			double len1 = mid_pt.distance(last_pt);
//			double len2 = mid_pt.distance(next_pt);
//			// length ratio
//			if(len2 > len1){
//				if(len2 / len1 > ControlSpace2dAdaptive::param_stretch_max_ratio){
//					double new_node = nodes[i-1]*0.7 + nodes[i]*0.3;
//					nodes.insertAt(i, new_node);
//					++node_count;
//					LOG4CPLUS_WARN(MeshLog::logger_console, "too large incident edge ratio.");
//				}else{
//					last_pt = mid_pt;
//					mid_pt = next_pt;
//					++i;
//				}
//			}else{
//				if(len1 / len2 > ControlSpace2dAdaptive::param_stretch_max_ratio){
//					double n2 = (i>1)?nodes[i-2]:0.0;
//					double new_node = nodes[i-1]*0.7 + n2*0.3;
//					nodes.insertAt(i-1, new_node);
//					++node_count;
//					last_pt = surface->getPoint(mesh_edge->getPoint(new_node));
//					++i;
//					LOG4CPLUS_WARN(MeshLog::logger_console, "too large incident edge ratio.");
//				}else{
//					last_pt = mid_pt;
//					mid_pt = next_pt;
//					++i;
//				}
//			}
//		}
//		if(node_count > 1){
//			assert(nodes[0] > 0.0);
//			assert(nodes[node_count-2] < 1.0);
//		}
//*/
//		// Check conformity of control space for incident faces
//		if(face_count > 1){
//			for(int k = 0; k < face_count; k++){
//				if(k == selected_face) continue;
//				MeshEdge2d* mesh_edge_k = mesh_edges[k];
//				Metric2dContext incident_mc(control_spaces[k]);
//				SurfaceParametric* surface_k = surfaces[k];
//				bool same_dir = (
//					mesh_edge_k->getMeshPoint(0)->getPtrTag(TagExtended::TAG_MP_2D_3D) ==
//					mesh_edge->getMeshPoint(0)->getPtrTag(TagExtended::TAG_MP_2D_3D));
//				DPoint2d pt_0 = mesh_edge_k->getMeshPoint(same_dir?0:1)->getCoordinates();
//				double ksi = same_dir ? 0.0 : 1.0;
//				double min_ksi = 0.0;
//				for(int i = 0; i < node_count; ){
//					const DPoint3d pt3 = surface->getPoint(mesh_edge->getPoint(nodes[i]));
//					const DPoint2d pt_1 = mesh_edge_k->surfaceParameters(surface_k, pt3, ksi, same_dir);
//					if(nodes[i] - ((i>0)?nodes[i-1]:min_ksi) > 1e-3){///???????????????????
//						incident_mc.countMetricAtPoint(DPoint2d::average(pt_0, pt_1));
//						const DVector2d dpt = pt_1-pt_0;
//						double dist = incident_mc.transformPStoMS(dpt).length2();
//						if(dist > max_metric_confrom_ratio2){
//							incident_mc.countMetricAtPoint(pt_0);
//							double dist_0 = incident_mc.transformPStoMS(dpt).length2();
//							incident_mc.countMetricAtPoint(pt_1);
//							double dist_1 = incident_mc.transformPStoMS(dpt).length2();
//							double dist_min = std::min(dist, std::min(dist_0, dist_1));
//							double dist_max = std::max(dist, std::max(dist_0, dist_1));
//							if(dist_max / dist_min < 2.0){
//								LOG4CPLUS_WARN(MeshLog::logger_console, "too large incident control-space nonconformity.");
//								double new_node = 0.5*(nodes[i] + ((i>0)?nodes[i-1]:min_ksi));
//								nodes.insertAt(i, new_node);
//								++node_count;
//							}else{
//								LOG4CPLUS_WARN(MeshLog::logger_console, "too large metric fluctuations.");
//								pt_0 = pt_1;
//								++i;
//							}
//						}else{
//							pt_0 = pt_1;
//							++i;
//						}
//					}else{
//						pt_0 = pt_1;
//						++i;
//					}
//				}
//			}
//		}
//		if(node_count > 1){
//			assert(nodes[0] > 0.0);
//			assert(nodes[node_count-2] < 1.0);
//		}
//
//		node_count = nodes.countInt();
//		// Check parity (if required)
//		if(param_even_node_count){
//			if(node_count == 1){ //(with the additional "1.0")
//				nodes.insertAt(0, 0.5);
//				++node_count;
//			}else if(node_count % 2 == 1){ //(with the additional "1.0")
//				// Insert an additional node (or delete?)
//				int mid_node = node_count/2;
//				int first_node = 0;
//				int last_node = node_count - 1;
//				int nct = node_count;
//				if(param_smoothen_last_node > 0){
//					nct = std::min(param_smoothen_last_node, nct);
//					first_node = std::max(first_node, mid_node-nct);
//					last_node  = std::min(last_node,  mid_node+nct);
//				}
//				double *dts = new double[last_node-first_node];
//				double first_t = (first_node > 0)?nodes[0]:0.0;
//				double temp_t = first_t;
//				double short_len = 0.0;
//				for(int i = first_node, j=0; i <= last_node; i++){
//					if(i != mid_node) short_len += (dts[j++] = nodes[i]-temp_t);
//					temp_t = nodes[i];
//				}
//				double ratio = (nodes[last_node] - first_t)/ short_len;
//				nodes.removeOrderedAt(mid_node);
//				temp_t = first_t;
//				for(int i = first_node; i < last_node; i++){
//					temp_t = nodes[i] = temp_t + ratio * dts[i-first_node];
//				}
//				--node_count;
//				delete[] dts;
//			}
//		}
//		if(node_count > 1){
//			assert(nodes[0] > 0.0);
//			assert(nodes[node_count-2] < 1.0);
//		}
//
//		nodes.removeLast(); // remove the obsolete "1.0"
//		disc_points.prepare(node_count);
//		bool same_dir = (mesh_edge->getMeshPoint(0)->getPtrTag(TagExtended::TAG_MP_2D_3D) ==
//			edge3d->getMeshPoint(0));
//		MeshPoint3d* last_point = edge3d->getMeshPoint(0);
//		--node_count;
////		int ict = 0;
//		for(int i = 0; i < node_count; i++){
//			int inode = same_dir?i:(node_count-i-1);
//			DPoint3d pt = surface->getPoint(mesh_edge->getPoint(nodes[inode]));
//			double dist = last_point->getCoordinates().distance2(pt);
//			if(dist > mesh_data.relative_small_number){
//				disc_points.add(last_point = new MeshPoint3d(pt));
//				last_point->setDoubleTag(TagExtended::TAG_DISCRETIZATION_PARAM, nodes[inode]);
//			}else{
//				LOG4CPLUS_WARN(MeshLog::logger_console, "Boundary nodes too close");
//			}
//		}
//		return disc_points.countInt();
//	}else if(param_even_node_count){
//		DPoint3d pt = surface->getPoint(mesh_edge->getPoint(0.5));
//		disc_points.add(new MeshPoint3d(pt));
//		disc_points[0]->setDoubleTag(TagExtended::TAG_DISCRETIZATION_PARAM, 0.5);
//		return 1;
//	}
//	return 0;
//}
