///////////////////////////////////////////////////////////////////////////
// ControlSpace3dOctree.cpp
// This class implements the control space (sizing and stretching info)
//  in the form of an Octree
/////////////////////////////////////////////////////////////////////////////
//  Tomasz Jurczyk, 2006-
//  Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include "ControlSpace3dOctree.h"
#include "MeshPoint3d.h"
#include "MeshEdge3d.h"
#include "MeshFace.h"
#include "MeshContainer3d.h"
#include "ControlSpace2dAdaptive.h"
#include "common.h"
#include "DTriangle.h"
#include "DSegment.h"
#include "DTetrahedron.h"
#include "MeshViewSet.h"
#include "MeshContainer3dSurface.h"

#ifdef _DEBUG
//#define STORE_CONTROL_EPS
#endif

int ControlSpace3dOctree::param_max_depth = 20;

double ControlSpace3dOctree::OctLeaf::param_octree_reduce_threshold = 0.2;
/// balance level
int ControlSpace3dOctree::OctLeaf::param_balance_level = 1;

const int ControlSpace3dOctree::OctLeaf::FEDATA[12][4] = {
	{ VX3D_LSW, VX3D_LSE, FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL }, // lower edges OXY
	{ VX3D_LNW, VX3D_LNE, FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL },
	{ VX3D_LSW, VX3D_LNW, FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL },
	{ VX3D_LSE, VX3D_LNE, FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL },
	{ VX3D_HSW, VX3D_HSE, FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL }, // higher edges OXY
	{ VX3D_HNW, VX3D_HNE, FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL },
	{ VX3D_HSW, VX3D_HNW, FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL },
	{ VX3D_HSE, VX3D_HNE, FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL },
	{ VX3D_LSW, VX3D_HSW, FORBIDDEN_EDGE_ZH, FORBIDDEN_EDGE_ZL }, // vertical edges
	{ VX3D_LSE, VX3D_HSE, FORBIDDEN_EDGE_ZH, FORBIDDEN_EDGE_ZL },
	{ VX3D_LNW, VX3D_HNW, FORBIDDEN_EDGE_ZH, FORBIDDEN_EDGE_ZL },
	{ VX3D_LNE, VX3D_HNE, FORBIDDEN_EDGE_ZH, FORBIDDEN_EDGE_ZL }
};

ControlSpace3dOctree::ControlSpace3dOctree(const DBox& box, int nxyz) 
	: ControlSpace3dKdTree(box), m_grid_vertices(3*nxyz)
{
	LOG4CPLUS_INFO(MeshLog::logger_console, "Creating octree control space.");
	assert(nxyz >= 8);
	// length
	double len[3] = { box.getDX(), box.getDY(), box.getDZ() };
	double total_len = len[0]+len[1]+len[2];
	double fl[3] = { len[0] / total_len, len[1] / total_len, len[2] / total_len };

	double ratio = pow(nxyz / (fl[0]*fl[1]*fl[2]), 1.0/3.0);
	for(int i = 0; i < 3; i++)
		m_n[i] = std::min(std::max(2,(int)(fl[i]*ratio)), nxyz/3);

	double d[3] = { len[0]/m_n[0], len[1]/m_n[1], len[2]/m_n[2] };

	m_grid = new OctLeaf[m_n[0]*m_n[1]*m_n[2]];
	
	// add regular grid vertices
	ControlNode3d qv(m_box.x0, m_box.y0, m_box.z0);
	for(int iz = 0; iz <= m_n[2]; iz++, qv.coord.z += d[2]){
		qv.coord.y = m_box.y0;
		for(int iy = 0; iy <= m_n[1]; iy++, qv.coord.y += d[1]){
			qv.coord.x = m_box.x0;
			for(int ix = 0; ix <= m_n[0]; ix++, qv.coord.x += d[0])
				m_grid_vertices.add( std::make_shared<ControlNode3d>(qv.coord, qv.control_data) );
		}
	}

	// init regular grid boxes with coordinates and vertices indices
	const DVector3d dpt(d[0]/4, d[1]/4, d[2]/4);
	const int nxy = m_n[0]*m_n[1];
	const int nx1 = m_n[0]+1;
	const int nxy1 = (m_n[0]+1)*(m_n[1]+1);
	DPoint3d pt(m_box.x0 + d[0]/2, m_box.y0 + d[1]/2, m_box.z0 + d[2]/2);
	for(int iz = 0, k = 0, in = 0; iz < m_n[2]; iz++, in+=nx1){
		for(int iy = 0; iy < m_n[1]; iy++, in++){
			for(int ix = 0; ix < m_n[0]; ix++, in++, k++){
				m_grid[k].init(pt, dpt);
				pt.x += d[0];
				// neighbours
				m_grid[k].setSingleNeighbour(FC3D_LOW,   (iz > 0) ? &m_grid[k-nxy] : nullptr);
				m_grid[k].setSingleNeighbour(FC3D_HIGH,  (iz < m_n[2]-1) ? &m_grid[k+nxy] : nullptr);
				m_grid[k].setSingleNeighbour(FC3D_SOUTH, (iy > 0) ? &m_grid[k-m_n[0]] : nullptr);
				m_grid[k].setSingleNeighbour(FC3D_NORTH, (iy < m_n[1]-1) ? &m_grid[k+m_n[0]] : nullptr);
				m_grid[k].setSingleNeighbour(FC3D_WEST,  (ix > 0) ? &m_grid[k-1] : nullptr);
				m_grid[k].setSingleNeighbour(FC3D_EAST,  (ix < m_n[0]-1) ? &m_grid[k+1] : nullptr);
				// vertices
				m_grid[k].setSingleVertex(VX3D_LSW, m_grid_vertices[in]);
				m_grid[k].setSingleVertex(VX3D_LSE, m_grid_vertices[in+1]);
				m_grid[k].setSingleVertex(VX3D_LNW, m_grid_vertices[in+nx1]);
				m_grid[k].setSingleVertex(VX3D_LNE, m_grid_vertices[in+nx1+1]);
				m_grid[k].setSingleVertex(VX3D_HSW, m_grid_vertices[in+nxy1]);
				m_grid[k].setSingleVertex(VX3D_HSE, m_grid_vertices[in+nxy1+1]);
				m_grid[k].setSingleVertex(VX3D_HNW, m_grid_vertices[in+nxy1+nx1]);
				m_grid[k].setSingleVertex(VX3D_HNE, m_grid_vertices[in+nxy1+nx1+1]);
			}
			pt.y += d[1];
			pt.x = m_box.x0 + d[0]/2;
		}
		pt.z += d[2];
		pt.y = m_box.y0 + d[1]/2;
		pt.x = m_box.x0 + d[0]/2;
	}

	assert(valid());
}

ControlSpace3dOctree::~ControlSpace3dOctree()
{
	if(m_grid) delete[] m_grid;
}

/// Returns number of control nodes in adaptive control structure
int ControlSpace3dOctree::getControlNodesCount() const
{
	return m_grid_vertices.countInt();
}

/// Invoke function for all control nodes of this space (read-only)
void ControlSpace3dOctree::forEachControlNode(const std::function<void(const ControlNode3d&node)>& fg) const {
	for (int i = m_grid_vertices.countInt() - 1; i >= 0; i--)
		fg( *m_grid_vertices[i] );
}

/// Invoke function for all control nodes of this space
void ControlSpace3dOctree::forEachControlNode(const std::function<void(ControlNode3d&node)>& fg) {
	for (int i = m_grid_vertices.countInt() - 1; i >= 0; i--)
		fg(*m_grid_vertices[i]);
}

/////////////////////////////////////////////////////////////////
// Interpolate metric data for not-initialized Oct tree nodes	
bool ControlSpace3dOctree::interpolate()
{
	// calculate data from source points
	int leaf_count = m_n[0] * m_n[1] * m_n[2];
	for(int i = 0; i < leaf_count; i++)
		m_grid[i].gatherDataFromSourcePoints();

	int gv_count = m_grid_vertices.countInt();
	int ready_count = 0;
	for(int i = 0; i < gv_count; i++){
		auto qv = m_grid_vertices[i];
		if(qv->w > 0.0){
			qv->control_data /= qv->w;
			qv->w = -1.0; // proper initialization
		}
		if(qv->w < 0.0) ++ready_count;
	}
	// + extended metric sources
//	int ct = m_ext_source_points.countInt();
//	for(int i = 0; i < ct; i++)
//		ready_count += setMinControlValue(*(m_ext_source_points.get(i)));

	// propagate control values from within leaves
	for(int i = 0; (ready_count < gv_count) && (i < leaf_count); i++)
		ready_count += m_grid[i].propagateUp();
	// extrapolate to other top-grid nodes
	if(ready_count < gv_count)
		ready_count += extrapolateMainControlNodes();
	// propagate control values into leaves
	for(int i = 0; (ready_count < gv_count) && (i < leaf_count); i++)
		ready_count += m_grid[i].propagateDown();

	if(ready_count == gv_count) m_initialized = 1;

	return (ready_count == gv_count);
}

int ControlSpace3dOctree::extrapolateMainControlNodes()
{
	// Start from the initialized nodes
	int v[3] = { m_n[0]+1, m_n[1]+1, m_n[2]+1 };
	int vx = v[0];
	int vxy = v[0]*v[1];
	int main_count = vxy*v[2];
	DataVector<int> ready_nodes(main_count);
	DataVector<int> levels(main_count, 0);
	for(int i = 0; i < main_count; i++){
		auto qv = m_grid_vertices[i];
		if(qv->w < 0.0){
			ready_nodes.add(i);
		}else
			levels[i] = main_count;
	}
	int new_nodes[6];
	int nearest_nodes[6];
	int start_ready_nodes = ready_nodes.countInt();
	if(start_ready_nodes == 0) return false;
	for(int i = 0; ready_nodes.countInt() < main_count; i++){
		int k = ready_nodes[i];
		int cur_level = levels[k];
		int kxy = k % vxy;
		int iz = k / vxy;
		int iy = kxy / vx;
		int ix = kxy % vx;
		int ct = 0;
		// check neighbours of 'k'
		if(ix > 0 && levels[k-1] == main_count)			new_nodes[ct++] = k-1;
		if(ix < m_n[0] && levels[k+1] == main_count)	new_nodes[ct++] = k+1;
		if(iy > 0 && levels[k-vx] == main_count)		new_nodes[ct++] = k-vx;
		if(iy < m_n[1] && levels[k+vx] == main_count)	new_nodes[ct++] = k+vx;
		if(iz > 0 && levels[k-vxy] == main_count)		new_nodes[ct++] = k-vxy;
		if(iz < m_n[2] && levels[k+vxy] == main_count)	new_nodes[ct++] = k+vxy;
		// if any new, calculate new value and add to ready_nodes
		for(int l = 0; l < ct; l++){
			k = new_nodes[l];
			kxy = k % vxy;
			iz = k / vxy;
			iy = kxy / vx;
			ix = kxy % vx;
			int j, nct = 0;
			for(j = 1; ix-j >= 0 && levels[k-j] > cur_level; j++);	// left
			if(ix-j >= 0) nearest_nodes[nct++] = k-j;
			for(j = 1; ix+j < v[0] && levels[k+j] > cur_level; j++);	// right
			if(ix+j < v[0]) nearest_nodes[nct++] = k+j;
			for(j = 1; iy-j >= 0 && levels[k-vx*j] > cur_level; j++);	// south
			if(iy-j >= 0) nearest_nodes[nct++] = k-vx*j;
			for(j = 1; iy+j < v[1] && levels[k+vx*j] > cur_level; j++);	// north
			if(iy+j < v[1]) nearest_nodes[nct++] = k+vx*j;
			for(j = 1; iz-j >= 0 && levels[k-vxy*j] > cur_level; j++);	// low
			if(iz-j >= 0) nearest_nodes[nct++] = k-vxy*j;
			for(j = 1; iz+j < v[2] && levels[k+vxy*j] > cur_level; j++);	// high
			if(iz+j < v[2]) nearest_nodes[nct++] = k+vxy*j;
			assert(nct > 0);
			auto qv = m_grid_vertices[k];
			assert(qv->w == 0.0);
			for(int m = 0; m < nct; m++){
				auto nqv = m_grid_vertices[nearest_nodes[m]];
				double dist = qv->coord.distance2(nqv->coord);
				double nw = (dist < mesh_data.relative_small_number)?mesh_data.relative_infinity:(1.0/dist);
				qv->control_data += nqv->control_data * nw;
				qv->w += nw;
			}
			qv->control_data /= qv->w;
			qv->w = -3.0;
			levels[k] = cur_level+1;
			ready_nodes.add(k);
		}
	}
	return (main_count - start_ready_nodes);
}

/////////////////////////////////////////////////////////////////////
// Add discrete metric information
void ControlSpace3dOctree::addControlNode(const ControlNode3d& node)
{
	assert(!m_initialized);
	findLastLeaf(node.coord)->adaptAndInsertControlPoint(m_grid_vertices, 
		std::make_shared<ControlNode3d>(node.coord, node.control_data));
}

ControlDataMatrix3d ControlSpace3dOctree::getMetricAtPoint(const DPoint3d& pt) const
{
	assert(m_initialized>0);
	return findLastLeaf(pt)->getMetricAtPoint(pt);
}

bool ControlSpace3dOctree::setMinControl(const DPoint3d& pt, const ControlDataMatrix3d& cdm, bool min_value_set)
{
	assert(m_initialized);
	assert(valid());
	return findLastLeaf(pt)->adaptAndSetControlPoint(m_grid_vertices, 
		std::make_shared<ControlNode3d>(pt, cdm), min_value_set);
}

/**
 * Calculates the int coordinates of left-bottom box (of matrix) containg the given point
 */
void ControlSpace3dOctree::countLocalCoordinates(const DPoint3d &pt, int in[]) const
{
	in[0] = (int) (m_n[0] * (pt.x - m_box.x0) / m_box.getDX());
	in[1] = (int) (m_n[1] * (pt.y - m_box.y0) / m_box.getDY());
	in[2] = (int) (m_n[2] * (pt.z - m_box.z0) / m_box.getDZ());
	for(int i = 0; i < 3; i++){
		if(in[i] > (m_n[i] - 1)) in[i] = m_n[i] - 1;
		if(in[i] < 0) in[i] = 0;
	}
}

/**
 * Calculates the int coordinates of left-bottom box (of matrix) containg the given point
 */
int ControlSpace3dOctree::countLocalCoordinates(const DPoint3d &pt) const
{
	int in[3];
	countLocalCoordinates(pt, in);
	return in[2]*m_n[1]*m_n[0] + in[1] * m_n[0] + in[0];
}

ControlDataMatrix3d ControlSpace3dOctree::OctLeaf::getMetricAtPoint(const DPoint3d& pt) const
{
	assert(!isSplit());

	double tx = (pt.x - m_middle.x + 2*m_dl.x) / (4*m_dl.x);
	double ty = (pt.y - m_middle.y + 2*m_dl.y) / (4*m_dl.y);
	double tz = (pt.z - m_middle.z + 2*m_dl.z) / (4*m_dl.z);

	//enum OctVertexWhich {VX3D_LSW = 0, VX3D_LSE = 1, VX3D_LNW = 2, VX3D_LNE = 3,
	//		VX3D_HSW = 4, VX3D_HSE = 5, VX3D_HNW = 6, VX3D_HNE = 7,

	double coeff[8] = {
		(1-tx)* (1-ty)* (1-tz),	// VX3D_LSW
		   tx * (1-ty)* (1-tz),	// VX3D_LSE
		(1-tx)*    ty * (1-tz),	// VX3D_LNW
		   tx *    ty * (1-tz),	// VX3D_LNE
		(1-tx)* (1-ty)*    tz ,	// VX3D_HSW
		   tx * (1-ty)*    tz ,	// VX3D_HSE
		(1-tx)*    ty *    tz ,	// VX3D_HNW
		   tx *    ty *    tz 	// VX3D_HNE
	};	
	ControlDataMatrix3d ave = m_vertices[0]->control_data * coeff[0];
	for(int i = 1; i < 8; i++)
		ave += m_vertices[i]->control_data * coeff[i];
	return ave;
}

double ControlSpace3dOctree::OctLeaf::interpolateDoubleTag(const DPoint3d& pt, TagExtended::TagType type) const
{
	assert(!isSplit());

	double tx = (pt.x - m_middle.x + 2*m_dl.x) / (4*m_dl.x);
	double ty = (pt.y - m_middle.y + 2*m_dl.y) / (4*m_dl.y);
	double tz = (pt.z - m_middle.z + 2*m_dl.z) / (4*m_dl.z);

	double coeff[8] = {
		(1-tx)* (1-ty)* (1-tz),	// VX3D_LSW
		   tx * (1-ty)* (1-tz),	// VX3D_LSE
		(1-tx)*    ty * (1-tz),	// VX3D_LNW
		   tx *    ty * (1-tz),	// VX3D_LNE
		(1-tx)* (1-ty)*    tz ,	// VX3D_HSW
		   tx * (1-ty)*    tz ,	// VX3D_HSE
		(1-tx)*    ty *    tz ,	// VX3D_HNW
		   tx *    ty *    tz 	// VX3D_HNE
	};	

	double ave = m_vertices[0]->getDoubleTag(type) * coeff[0];
	for(int i = 1; i < 8; i++)
		ave += m_vertices[i]->getDoubleTag(type) * coeff[i];
	return ave;
}

double ControlSpace3dOctree::OctLeaf::getMetricGradationRatio(const DPoint3d& pt) const
{
	assert(!isSplit());

	double tx = (pt.x - m_middle.x + 2*m_dl.x) / (4*m_dl.x);
	double ty = (pt.y - m_middle.y + 2*m_dl.y) / (4*m_dl.y);
	double tz = (pt.z - m_middle.z + 2*m_dl.z) / (4*m_dl.z);

	double coeff[8] = {
		(1-tx)* (1-ty)* (1-tz),	// VX3D_LSW
		   tx * (1-ty)* (1-tz),	// VX3D_LSE
		(1-tx)*    ty * (1-tz),	// VX3D_LNW
		   tx *    ty * (1-tz),	// VX3D_LNE
		(1-tx)* (1-ty)*    tz ,	// VX3D_HSW
		   tx * (1-ty)*    tz ,	// VX3D_HSE
		(1-tx)*    ty *    tz ,	// VX3D_HNW
		   tx *    ty *    tz 	// VX3D_HNE
	};	
	assert(!m_vertices[0]->gradationUnknown());
	double ave = m_vertices[0]->max_gradation_ratio * coeff[0];
	for (int i = 1; i < 8; i++) {
		assert(!m_vertices[i]->gradationUnknown());
		ave += m_vertices[i]->max_gradation_ratio * coeff[i];
	}
	return ave;
}

ControlSpace3dOctree::OctLeafExtra::OctLeafExtra(int level) 
	: m_level(level), m_source_points(nullptr) 
{ 
	for(int i = FC3D_FIRST; i <= FC3D_LAST; i++){
		m_neighbours[i] = nullptr; 
		m_midface_vertices[i] = nullptr;
	}
	for(int i = ED3D_FIRST; i <= ED3D_LAST; i++) m_midedge_vertices[i] = nullptr;
}

ControlSpace3dOctree::OctLeaf::OctLeaf() : m_extra(nullptr) 
{ 
	m_leaves[0] = nullptr;
}

ControlSpace3dOctree::OctLeaf::~OctLeaf() { 
	if(m_leaves[0]) 
		for(int i = 0; i < 8; i++) 
			delete m_leaves[i]; 
	if(m_extra) 
		delete m_extra;
}

ControlSpace3dOctree::OctLeaf::OctLeaf(const DPoint3d& middle, const DVector3d& dl, int level) 
	: m_middle(middle), m_dl(dl)
{
	m_leaves[0] = nullptr;
	m_extra = new OctLeafExtra(level);
}

void ControlSpace3dOctree::OctLeaf::init(const DPoint3d& middle, const DVector3d& dl, int level)
{
	m_middle = middle;
	m_dl = dl;
	assert(!m_extra);
	m_extra = new OctLeafExtra(level);
}

void ControlSpace3dOctree::OctLeaf::adaptAndInsertControlPoint(
	OctVertexList& grid_vertices, std::shared_ptr<ControlNode3d> mqv)
{
	if(isSplit()){
		findLastLeaf(mqv->coord)->adaptAndInsertControlPoint(grid_vertices, mqv);
		return;
	}
	assert(mqv->control_data.det() > 0.0);
	DMetric3d dmp(mqv->control_data);
	double len_x = dmp.transformRStoMS(DVector3d(4*m_dl.x, 0.0, 0.0)).length2();
	double len_y = dmp.transformRStoMS(DVector3d(0.0, 4*m_dl.y, 0.0)).length2();
	double len_z = dmp.transformRStoMS(DVector3d(0.0, 0.0, 4*m_dl.z)).length2();
	bool split_needed = (len_x > METRIC_LENGTH_RATIO2 || 
		len_y > METRIC_LENGTH_RATIO2 || len_z > METRIC_LENGTH_RATIO2);

	auto qv = std::make_shared<ControlNode3d>(mqv->coord, mqv->control_data);

	if(!split_needed && m_extra->m_source_points){	// further check
		// compare new source node with already existing
		// check with:
		//	 all points (used now)
		//	(?) average, minimum, some kind of approximation
		double max_diff = 0.0;
		for(int i = 0; i < m_extra->m_source_points->countInt(); i++){
			auto sqv = m_extra->m_source_points->get(i);
			if(dmp.transformRStoMS(qv->coord - sqv->coord).length() < METRIC_SMALL_NUMBER){
				// combine
				qv->control_data.setMinimum(sqv->control_data);
				m_extra->m_source_points->removeAt(i);
				max_diff = 0.0;
				i = -1; // start again
				continue;
			}
			double diff = qv->control_data.countDifferenceRR(sqv->control_data);
			if(diff > max_diff){
				max_diff = diff;
				if(max_diff > ControlSpace2dAdaptive::param_threshold_diff) break;
			}
		}
		split_needed = (max_diff > ControlSpace2dAdaptive::param_threshold_diff);
	}

	if(!m_extra->m_source_points) m_extra->m_source_points = new DataVector<std::shared_ptr<ControlNode3d>>;
	m_extra->m_source_points->add(qv);

	if(split_needed){
		OctVertexList split_vertices(100);
		split(grid_vertices, split_vertices); // split and redistribute source points
	}
}

bool ControlSpace3dOctree::OctLeaf::adaptAndSetControlPoint(
	OctVertexList& grid_vertices, std::shared_ptr<ControlNode3d> qv, bool min_value_set)
{
	if(isSplit()){
		return findLastLeaf(qv->coord)->adaptAndSetControlPoint(grid_vertices, qv, min_value_set);
	}
	// compare new source node with already existing
	// check with calculated value:
	const ControlDataMatrix3d cdm = getMetricAtPoint(qv->coord);
	double diff = qv->control_data.countDifferenceRR(cdm);
	bool any_changes = false;
	if(diff > ControlSpace2dAdaptive::param_threshold_diff){
		// split may be unnecessary, if qv.control_data is greater than already set for this leaf...
		ControlDataMatrix3d min_cdm = cdm;
		min_cdm.setMinimum(qv->control_data);
		if(cdm.countDifferenceRR(min_cdm) < ControlSpace2dAdaptive::param_threshold_diff){
			// qv different from current metric, but introduces no changes - so skip
			return false;
		}
		DMetric3d dmp(qv->control_data);
		double len_x = dmp.transformRStoMS(DVector3d(4*m_dl.x, 0.0, 0.0)).length2();
		double len_y = dmp.transformRStoMS(DVector3d(0.0, 4*m_dl.y, 0.0)).length2();
		double len_z = dmp.transformRStoMS(DVector3d(0.0, 0.0, 4*m_dl.z)).length2();
		bool split_needed = (len_x > METRIC_LENGTH_RATIO2 || 
			len_y > METRIC_LENGTH_RATIO2 || len_z > METRIC_LENGTH_RATIO2 );
		if(split_needed){
			OctVertexList split_vertices(100);
			if(split(grid_vertices, split_vertices, true)){ // split (+ set values for new nodes)
				findLastLeaf(qv->coord)->adaptAndSetControlPoint(grid_vertices, qv, min_value_set);
				return true;
			}
		}
		if(min_value_set){
			for(int i = VX3D_FIRST; i <= VX3D_LAST; i++)
				any_changes |= m_vertices[i]->control_data.setMinimum(qv->control_data);
		}
	}

	return any_changes;
}

void ControlSpace3dOctree::OctLeaf::setNeighbour(int side, OctLeaf* nb, bool skip_first)
{
	if(!skip_first) m_extra->m_neighbours[side] = nb;
	if(isSplit())
		for(int i = 0; i < 4; i++)
			m_leaves[DBox::face_to_vertex[side][i]]->setNeighbour(side, nb);
}

/// Set mid-edge vertex for selected edge
bool ControlSpace3dOctree::OctLeaf::setMidEdgeVertex(CN3dPtr q0, CN3dPtr q1, CN3dPtr node)
{
	for(int i = ED3D_FIRST; i <= ED3D_LAST; i++){
		auto cn0 = m_vertices[DBox::edge_to_vertex[i][0]];
		auto cn1 = m_vertices[DBox::edge_to_vertex[i][1]];
		if((q0 == cn0 && q1 == cn1) || (q0 == cn1 && q1 == cn0)){
			assert(m_extra->m_midedge_vertices[i] == nullptr);
			m_extra->m_midedge_vertices[i] = node;
			return true;
		}
	}
	assert(false);
	return false;
}


bool ControlSpace3dOctree::OctLeaf::split(OctVertexList& grid_vertices, 
		OctVertexList& split_vertices, bool init_values, bool with_balancing)
{
	if(m_extra->m_level >= ControlSpace3dOctree::param_max_depth)
		return false;

	assert(m_leaves[VX3D_FIRST] == nullptr); // not split yet

	assert(valid());

	// Balance before split...
	if(with_balancing)
		balance(grid_vertices, split_vertices, init_values);

	assert(valid());

	const DVector3d mid_dl = m_dl/2;
	m_leaves[VX3D_LSW] = new OctLeaf(DPoint3d(m_middle.x-m_dl.x, m_middle.y-m_dl.y, 
		m_middle.z-m_dl.z), mid_dl, m_extra->m_level+1);
	m_leaves[VX3D_LSE] = new OctLeaf(DPoint3d(m_middle.x+m_dl.x, m_middle.y-m_dl.y, 
		m_middle.z-m_dl.z), mid_dl, m_extra->m_level+1);
	m_leaves[VX3D_LNW] = new OctLeaf(DPoint3d(m_middle.x-m_dl.x, m_middle.y+m_dl.y, 
		m_middle.z-m_dl.z), mid_dl, m_extra->m_level+1);
	m_leaves[VX3D_LNE] = new OctLeaf(DPoint3d(m_middle.x+m_dl.x, m_middle.y+m_dl.y, 
		m_middle.z-m_dl.z), mid_dl, m_extra->m_level+1);
	m_leaves[VX3D_HSW] = new OctLeaf(DPoint3d(m_middle.x-m_dl.x, m_middle.y-m_dl.y, 
		m_middle.z+m_dl.z), mid_dl, m_extra->m_level+1);
	m_leaves[VX3D_HSE] = new OctLeaf(DPoint3d(m_middle.x+m_dl.x, m_middle.y-m_dl.y, 
		m_middle.z+m_dl.z), mid_dl, m_extra->m_level+1);
	m_leaves[VX3D_HNW] = new OctLeaf(DPoint3d(m_middle.x-m_dl.x, m_middle.y+m_dl.y, 
		m_middle.z+m_dl.z), mid_dl, m_extra->m_level+1);
	m_leaves[VX3D_HNE] = new OctLeaf(DPoint3d(m_middle.x+m_dl.x, m_middle.y+m_dl.y, 
		m_middle.z+m_dl.z), mid_dl, m_extra->m_level+1);

	// middle
	auto qv_middle = std::make_shared<ControlNode3d>(m_middle);
	grid_vertices.add(qv_middle);
	split_vertices.add(qv_middle);
	if(init_values){
		for(int i = VX3D_FIRST; i <= VX3D_LAST; i++)
			qv_middle->control_data += m_vertices[i]->control_data;
		qv_middle->control_data *= 0.125;
		qv_middle->w = -5.0;
	}

	// check middle nodes at faces
	OctFaceWhich opposite_face[6] = {FC3D_NORTH, FC3D_SOUTH, FC3D_EAST, FC3D_WEST, FC3D_HIGH, FC3D_LOW};
	for(int i = FC3D_FIRST; i <= FC3D_LAST; i++){
		if(!m_extra->m_midface_vertices[i]){	// vertex missing
			CN3dPtr qf[4];
			for(int j = 0; j < 4; j++)
				qf[j] = m_vertices[DBox::face_to_vertex[i][j]];
			m_extra->m_midface_vertices[i] = std::make_shared<ControlNode3d>(
				DPoint3d::average(qf[0]->coord, qf[1]->coord, qf[2]->coord, qf[3]->coord));
			grid_vertices.add(m_extra->m_midface_vertices[i]);
			split_vertices.add(m_extra->m_midface_vertices[i]);
			if(m_extra->m_neighbours[i]){
				assert(m_extra->m_neighbours[i]->m_extra->m_midface_vertices[opposite_face[i]] == nullptr);
				m_extra->m_neighbours[i]->m_extra->m_midface_vertices[opposite_face[i]] = m_extra->m_midface_vertices[i];
			}
			if(init_values){
				m_extra->m_midface_vertices[i]->control_data = (qf[0]->control_data + qf[1]->control_data + 
					qf[2]->control_data + qf[3]->control_data)*0.25;
				m_extra->m_midface_vertices[i]->w = -5.0;
			}
		}
	}

	// check middle nodes at edges
	for(int i = ED3D_FIRST; i <= ED3D_LAST; i++){
		if(!m_extra->m_midedge_vertices[i]){	// vertex missing
			auto q0 = m_vertices[DBox::edge_to_vertex[i][0]];
			auto q1 = m_vertices[DBox::edge_to_vertex[i][1]];
			m_extra->m_midedge_vertices[i] = std::make_shared<ControlNode3d>(DPoint3d::average(q0->coord, q1->coord));
			grid_vertices.add(m_extra->m_midedge_vertices[i]);
			split_vertices.add(m_extra->m_midedge_vertices[i]);
			// 3 neighbours of edge
			int side[2] = { DBox::edge_to_face[i][0], DBox::edge_to_face[i][1] };
			OctLeaf* nbs[3] = { m_extra->m_neighbours[side[0]], m_extra->m_neighbours[side[1]], 
				m_extra->m_neighbours[side[0]] ? m_extra->m_neighbours[side[0]]->m_extra->m_neighbours[side[1]] : nullptr};
			for(int j = 0; j < 3; j++)
				if(nbs[j]) nbs[j]->setMidEdgeVertex(q0, q1, m_extra->m_midedge_vertices[i]);
			if(init_values){
				m_extra->m_midedge_vertices[i]->control_data = (q0->control_data + q1->control_data)*0.5;
				m_extra->m_midedge_vertices[i]->w = -5.0;
			}
		}
	}

	// neighbours
	// - between leaves
	m_leaves[VX3D_LSW]->m_extra->m_neighbours[FC3D_EAST] = m_leaves[VX3D_LSE];
	m_leaves[VX3D_LSW]->m_extra->m_neighbours[FC3D_NORTH] = m_leaves[VX3D_LNW];
	m_leaves[VX3D_LSW]->m_extra->m_neighbours[FC3D_HIGH] = m_leaves[VX3D_HSW];

	m_leaves[VX3D_LSE]->m_extra->m_neighbours[FC3D_WEST] = m_leaves[VX3D_LSW];
	m_leaves[VX3D_LSE]->m_extra->m_neighbours[FC3D_NORTH] = m_leaves[VX3D_LNE];
	m_leaves[VX3D_LSE]->m_extra->m_neighbours[FC3D_HIGH] = m_leaves[VX3D_HSE];

	m_leaves[VX3D_LNW]->m_extra->m_neighbours[FC3D_EAST] = m_leaves[VX3D_LNE];
	m_leaves[VX3D_LNW]->m_extra->m_neighbours[FC3D_SOUTH] = m_leaves[VX3D_LSW];
	m_leaves[VX3D_LNW]->m_extra->m_neighbours[FC3D_HIGH] = m_leaves[VX3D_HNW];

	m_leaves[VX3D_LNE]->m_extra->m_neighbours[FC3D_WEST] = m_leaves[VX3D_LNW];
	m_leaves[VX3D_LNE]->m_extra->m_neighbours[FC3D_SOUTH] = m_leaves[VX3D_LSE];
	m_leaves[VX3D_LNE]->m_extra->m_neighbours[FC3D_HIGH] = m_leaves[VX3D_HNE];

	m_leaves[VX3D_HSW]->m_extra->m_neighbours[FC3D_EAST] = m_leaves[VX3D_HSE];
	m_leaves[VX3D_HSW]->m_extra->m_neighbours[FC3D_NORTH] = m_leaves[VX3D_HNW];
	m_leaves[VX3D_HSW]->m_extra->m_neighbours[FC3D_LOW] = m_leaves[VX3D_LSW];

	m_leaves[VX3D_HSE]->m_extra->m_neighbours[FC3D_WEST] = m_leaves[VX3D_HSW];
	m_leaves[VX3D_HSE]->m_extra->m_neighbours[FC3D_NORTH] = m_leaves[VX3D_HNE];
	m_leaves[VX3D_HSE]->m_extra->m_neighbours[FC3D_LOW] = m_leaves[VX3D_LSE];

	m_leaves[VX3D_HNW]->m_extra->m_neighbours[FC3D_EAST] = m_leaves[VX3D_HNE];
	m_leaves[VX3D_HNW]->m_extra->m_neighbours[FC3D_SOUTH] = m_leaves[VX3D_HSW];
	m_leaves[VX3D_HNW]->m_extra->m_neighbours[FC3D_LOW] = m_leaves[VX3D_LNW];

	m_leaves[VX3D_HNE]->m_extra->m_neighbours[FC3D_WEST] = m_leaves[VX3D_HNW];
	m_leaves[VX3D_HNE]->m_extra->m_neighbours[FC3D_SOUTH] = m_leaves[VX3D_HSE];
	m_leaves[VX3D_HNE]->m_extra->m_neighbours[FC3D_LOW] = m_leaves[VX3D_LNE];

	// - outside
	for(int i = FC3D_FIRST; i <= FC3D_LAST; i++){
		if(!m_extra->m_neighbours[i]) continue;
		int li[4] = { 
			DBox::face_to_vertex[i][0], DBox::face_to_vertex[i][1],
			DBox::face_to_vertex[i][2], DBox::face_to_vertex[i][3] };
		if(m_extra->m_neighbours[i]->isSplit()){
			int lj[4] = { 
				DBox::face_to_vertex[opposite_face[i]][0], 
				DBox::face_to_vertex[opposite_face[i]][1],
				DBox::face_to_vertex[opposite_face[i]][2], 
				DBox::face_to_vertex[opposite_face[i]][3] };
			for(int j = 0; j < 4; j++){
				m_extra->m_neighbours[i]->m_leaves[lj[j]]->setNeighbour(opposite_face[i], m_leaves[li[j]]);
				m_leaves[li[j]]->m_extra->m_neighbours[i] = m_extra->m_neighbours[i]->m_leaves[lj[j]];
			}
		}else{
			for(int j = 0; j < 4; j++)
				m_leaves[li[j]]->m_extra->m_neighbours[i] = m_extra->m_neighbours[i];
		}
	}
	// vertices 
	m_leaves[VX3D_LSW]->setSingleVertex(VX3D_LSW, m_vertices[VX3D_LSW]);
	m_leaves[VX3D_LSW]->setSingleVertex(VX3D_LSE, m_extra->m_midedge_vertices[ED3D_LS]);
	m_leaves[VX3D_LSW]->setSingleVertex(VX3D_LNW, m_extra->m_midedge_vertices[ED3D_LW]);
	m_leaves[VX3D_LSW]->setSingleVertex(VX3D_LNE, m_extra->m_midface_vertices[FC3D_LOW]);
	m_leaves[VX3D_LSW]->setSingleVertex(VX3D_HSW, m_extra->m_midedge_vertices[ED3D_SW]);
	m_leaves[VX3D_LSW]->setSingleVertex(VX3D_HSE, m_extra->m_midface_vertices[FC3D_SOUTH]);
	m_leaves[VX3D_LSW]->setSingleVertex(VX3D_HNW, m_extra->m_midface_vertices[FC3D_WEST]);
	m_leaves[VX3D_LSW]->setSingleVertex(VX3D_HNE, qv_middle);

	m_leaves[VX3D_LSE]->setSingleVertex(VX3D_LSW, m_extra->m_midedge_vertices[ED3D_LS]);
	m_leaves[VX3D_LSE]->setSingleVertex(VX3D_LSE, m_vertices[VX3D_LSE]);
	m_leaves[VX3D_LSE]->setSingleVertex(VX3D_LNW, m_extra->m_midface_vertices[FC3D_LOW]);
	m_leaves[VX3D_LSE]->setSingleVertex(VX3D_LNE, m_extra->m_midedge_vertices[ED3D_LE]);
	m_leaves[VX3D_LSE]->setSingleVertex(VX3D_HSW, m_extra->m_midface_vertices[FC3D_SOUTH]);
	m_leaves[VX3D_LSE]->setSingleVertex(VX3D_HSE, m_extra->m_midedge_vertices[ED3D_SE]);
	m_leaves[VX3D_LSE]->setSingleVertex(VX3D_HNW, qv_middle);
	m_leaves[VX3D_LSE]->setSingleVertex(VX3D_HNE, m_extra->m_midface_vertices[FC3D_EAST]);

	m_leaves[VX3D_LNW]->setSingleVertex(VX3D_LSW, m_extra->m_midedge_vertices[ED3D_LW]);
	m_leaves[VX3D_LNW]->setSingleVertex(VX3D_LSE, m_extra->m_midface_vertices[FC3D_LOW]);
	m_leaves[VX3D_LNW]->setSingleVertex(VX3D_LNW, m_vertices[VX3D_LNW]);
	m_leaves[VX3D_LNW]->setSingleVertex(VX3D_LNE, m_extra->m_midedge_vertices[ED3D_LN]);
	m_leaves[VX3D_LNW]->setSingleVertex(VX3D_HSW, m_extra->m_midface_vertices[FC3D_WEST]);
	m_leaves[VX3D_LNW]->setSingleVertex(VX3D_HSE, qv_middle);
	m_leaves[VX3D_LNW]->setSingleVertex(VX3D_HNW, m_extra->m_midedge_vertices[ED3D_NW]);
	m_leaves[VX3D_LNW]->setSingleVertex(VX3D_HNE, m_extra->m_midface_vertices[FC3D_NORTH]);

	m_leaves[VX3D_LNE]->setSingleVertex(VX3D_LSW, m_extra->m_midface_vertices[FC3D_LOW]);
	m_leaves[VX3D_LNE]->setSingleVertex(VX3D_LSE, m_extra->m_midedge_vertices[ED3D_LE]);
	m_leaves[VX3D_LNE]->setSingleVertex(VX3D_LNW, m_extra->m_midedge_vertices[ED3D_LN]);
	m_leaves[VX3D_LNE]->setSingleVertex(VX3D_LNE, m_vertices[VX3D_LNE]);
	m_leaves[VX3D_LNE]->setSingleVertex(VX3D_HSW, qv_middle);
	m_leaves[VX3D_LNE]->setSingleVertex(VX3D_HSE, m_extra->m_midface_vertices[FC3D_EAST]);
	m_leaves[VX3D_LNE]->setSingleVertex(VX3D_HNW, m_extra->m_midface_vertices[FC3D_NORTH]);
	m_leaves[VX3D_LNE]->setSingleVertex(VX3D_HNE, m_extra->m_midedge_vertices[ED3D_NE]);

	m_leaves[VX3D_HSW]->setSingleVertex(VX3D_LSW, m_extra->m_midedge_vertices[ED3D_SW]);
	m_leaves[VX3D_HSW]->setSingleVertex(VX3D_LSE, m_extra->m_midface_vertices[FC3D_SOUTH]);
	m_leaves[VX3D_HSW]->setSingleVertex(VX3D_LNW, m_extra->m_midface_vertices[FC3D_WEST]);
	m_leaves[VX3D_HSW]->setSingleVertex(VX3D_LNE, qv_middle);
	m_leaves[VX3D_HSW]->setSingleVertex(VX3D_HSW, m_vertices[VX3D_HSW]);
	m_leaves[VX3D_HSW]->setSingleVertex(VX3D_HSE, m_extra->m_midedge_vertices[ED3D_HS]);
	m_leaves[VX3D_HSW]->setSingleVertex(VX3D_HNW, m_extra->m_midedge_vertices[ED3D_HW]);
	m_leaves[VX3D_HSW]->setSingleVertex(VX3D_HNE, m_extra->m_midface_vertices[FC3D_HIGH]);

	m_leaves[VX3D_HSE]->setSingleVertex(VX3D_LSW, m_extra->m_midface_vertices[FC3D_SOUTH]);
	m_leaves[VX3D_HSE]->setSingleVertex(VX3D_LSE, m_extra->m_midedge_vertices[ED3D_SE]);
	m_leaves[VX3D_HSE]->setSingleVertex(VX3D_LNW, qv_middle);
	m_leaves[VX3D_HSE]->setSingleVertex(VX3D_LNE, m_extra->m_midface_vertices[FC3D_EAST]);
	m_leaves[VX3D_HSE]->setSingleVertex(VX3D_HSW, m_extra->m_midedge_vertices[ED3D_HS]);
	m_leaves[VX3D_HSE]->setSingleVertex(VX3D_HSE, m_vertices[VX3D_HSE]);
	m_leaves[VX3D_HSE]->setSingleVertex(VX3D_HNW, m_extra->m_midface_vertices[FC3D_HIGH]);
	m_leaves[VX3D_HSE]->setSingleVertex(VX3D_HNE, m_extra->m_midedge_vertices[ED3D_HE]);

	m_leaves[VX3D_HNW]->setSingleVertex(VX3D_LSW, m_extra->m_midface_vertices[FC3D_WEST]);
	m_leaves[VX3D_HNW]->setSingleVertex(VX3D_LSE, qv_middle);
	m_leaves[VX3D_HNW]->setSingleVertex(VX3D_LNW, m_extra->m_midedge_vertices[ED3D_NW]);
	m_leaves[VX3D_HNW]->setSingleVertex(VX3D_LNE, m_extra->m_midface_vertices[FC3D_NORTH]);
	m_leaves[VX3D_HNW]->setSingleVertex(VX3D_HSW, m_extra->m_midedge_vertices[ED3D_HW]);
	m_leaves[VX3D_HNW]->setSingleVertex(VX3D_HSE, m_extra->m_midface_vertices[FC3D_HIGH]);
	m_leaves[VX3D_HNW]->setSingleVertex(VX3D_HNW, m_vertices[VX3D_HNW]);
	m_leaves[VX3D_HNW]->setSingleVertex(VX3D_HNE, m_extra->m_midedge_vertices[ED3D_HN]);

	m_leaves[VX3D_HNE]->setSingleVertex(VX3D_LSW, qv_middle);
	m_leaves[VX3D_HNE]->setSingleVertex(VX3D_LSE, m_extra->m_midface_vertices[FC3D_EAST]);
	m_leaves[VX3D_HNE]->setSingleVertex(VX3D_LNW, m_extra->m_midface_vertices[FC3D_NORTH]);
	m_leaves[VX3D_HNE]->setSingleVertex(VX3D_LNE, m_extra->m_midedge_vertices[ED3D_NE]);
	m_leaves[VX3D_HNE]->setSingleVertex(VX3D_HSW, m_extra->m_midface_vertices[FC3D_HIGH]);
	m_leaves[VX3D_HNE]->setSingleVertex(VX3D_HSE, m_extra->m_midedge_vertices[ED3D_HE]);
	m_leaves[VX3D_HNE]->setSingleVertex(VX3D_HNW, m_extra->m_midedge_vertices[ED3D_HN]);
	m_leaves[VX3D_HNE]->setSingleVertex(VX3D_HNE, m_vertices[VX3D_HNE]);

	// mid_vertices
	for(int i = FC3D_FIRST; i <= FC3D_LAST; i++){
		if(!m_extra->m_neighbours[i] || !m_extra->m_neighbours[i]->isSplit()) continue;
		// i.e. leaves from ith side have theirs counterparts
		for(int j = 0; j < 4; j++){
			int vi = DBox::face_to_vertex[i][j];
			OctLeaf* nb = m_leaves[vi]->m_extra->m_neighbours[i];
			if(nb->m_extra->m_midface_vertices[opposite_face[i]])
				m_leaves[vi]->m_extra->m_midface_vertices[i] = nb->m_extra->m_midface_vertices[opposite_face[i]];
			for(int fi = 0; fi < 4; fi++){
				int ei = DBox::face_to_edge[opposite_face[i]][fi];
				if(nb->m_extra->m_midedge_vertices[ei])
					m_leaves[vi]->setMidEdgeVertex(nb->m_vertices[DBox::edge_to_vertex[ei][0]],
						nb->m_vertices[DBox::edge_to_vertex[ei][1]], nb->m_extra->m_midedge_vertices[ei]);
			}
		}
	}
	// done
	if(m_extra->m_source_points){	// split too
		int ct = m_extra->m_source_points->countInt();
		for(int i = VX3D_FIRST; i <= VX3D_LAST; i++)
			m_leaves[i]->m_extra->m_source_points = 
				new DataVector<std::shared_ptr<ControlNode3d>>(std::min(10,ct/4));
		for(int i = 0; i < ct; i++){
			auto qv = m_extra->m_source_points->get(i);
			findLastLeaf(qv->coord)->adaptAndInsertControlPoint(grid_vertices, qv);
		}
		delete m_extra->m_source_points;
		m_extra->m_source_points = nullptr;
	}

	assert(valid(2));
	for(int i = VX3D_FIRST; i <= VX3D_LAST; i++){
		assert(m_leaves[i]->valid());
	}

	return true;
}

void ControlSpace3dOctree::OctLeaf::balance(OctVertexList& grid_vertices,
	OctVertexList& split_vertices, bool init_values)
{
	// through faces
	for(int i = FC3D_FIRST; i <= FC3D_LAST; i++){
		if(!m_extra->m_neighbours[i]) continue;	// border of qtree
		while((m_extra->m_level - m_extra->m_neighbours[i]->m_extra->m_level) >= 1)
			if(!m_extra->m_neighbours[i]->split(grid_vertices, split_vertices, init_values)) break;
	}
	// through edges
	for(int i = ED3D_FIRST; i <= ED3D_LAST; i++){
		int side[2] = { DBox::edge_to_face[i][0], DBox::edge_to_face[i][1] };
		OctLeaf* ql = m_extra->m_neighbours[side[0]] 
			? m_extra->m_neighbours[side[0]]->m_extra->m_neighbours[side[1]] : nullptr;
		if(!ql) continue;
		if((m_extra->m_level - ql->m_extra->m_level) >= 1)
			ql->split(grid_vertices, split_vertices, init_values);
	}
}

ControlSpace3dOctree::OctLeaf* ControlSpace3dOctree::OctLeaf::findLastLeaf(const DPoint3d& pt)
{
	if(m_leaves[VX3D_FIRST] == nullptr) return this;

	if(pt.x < m_middle.x) 
		if(pt.y < m_middle.y) 
			if(pt.z < m_middle.z) return m_leaves[VX3D_LSW]->findLastLeaf(pt);
			else return m_leaves[VX3D_HSW]->findLastLeaf(pt);
		else 
			if(pt.z < m_middle.z) return m_leaves[VX3D_LNW]->findLastLeaf(pt);
			else return m_leaves[VX3D_HNW]->findLastLeaf(pt);
	else if(pt.y < m_middle.y) 
			if(pt.z < m_middle.z) return m_leaves[VX3D_LSE]->findLastLeaf(pt);
			else return m_leaves[VX3D_HSE]->findLastLeaf(pt);
		else 
			if(pt.z < m_middle.z) return m_leaves[VX3D_LNE]->findLastLeaf(pt);
			else return m_leaves[VX3D_HNE]->findLastLeaf(pt);
}

ControlSpace3dOctree::OctLeaf* ControlSpace3dOctree::findLastLeaf(const DPoint3d& pt) const
{
	int k = countLocalCoordinates(pt);
	return m_grid[k].findLastLeaf(pt);
}

bool ControlSpace3dOctree::valid() const 
{
	// temporary, since validation is not necessary until some problems arise
	return true;

	int count = m_n[0]*m_n[1]*m_n[2];
	for(int i = 0; i < count; i++)
		if(!m_grid[i].valid()) return false;
	return true;
}

bool ControlSpace3dOctree::OctLeaf::valid(int depth) const
{
	// check for repeated pairs of nodes per edge
	for(int i = ED3D_FIRST; i < ED3D_LAST; i++){
		auto cni0 = m_vertices[DBox::edge_to_vertex[i][0]];
		auto cni1 = m_vertices[DBox::edge_to_vertex[i][1]];
		for(int j = i+1; j <= ED3D_LAST; j++){
			auto cnj0 = m_vertices[DBox::edge_to_vertex[j][0]];
			auto cnj1 = m_vertices[DBox::edge_to_vertex[j][1]];

			bool ok = !((cni0 == cnj0) && (cni1 == cnj1));
			assert(ok);
			ok &= !((cni0 == cnj1) && (cni1 == cnj0));
			assert(ok);
			if(!ok) return false;
		}
	}
	if(!m_extra) return true; // if readonly, finish here
	// check edge-nodes neighbor consistency
	for(int i = ED3D_FIRST; i <= ED3D_LAST; i++){
		auto q0 = m_vertices[DBox::edge_to_vertex[i][0]];
		auto q1 = m_vertices[DBox::edge_to_vertex[i][1]];
		// 3 neighbours of edge
		int side[2] = { DBox::edge_to_face[i][0], DBox::edge_to_face[i][1] };
		OctLeaf* nbs[3] = { m_extra->m_neighbours[side[0]], m_extra->m_neighbours[side[1]], 
			m_extra->m_neighbours[side[0]] ? m_extra->m_neighbours[side[0]]->m_extra->m_neighbours[side[1]] : nullptr};
		for(int j = 0; j < 3; j++){
			if(!nbs[j] || nbs[j]->m_extra->m_level != m_extra->m_level) continue;
			int found = 0;
			for(int k = ED3D_FIRST; k <= ED3D_LAST; k++){
				auto nq0 = nbs[j]->m_vertices[DBox::edge_to_vertex[k][0]];
				auto nq1 = nbs[j]->m_vertices[DBox::edge_to_vertex[k][1]];
				if((nq0 == q0) && (nq1 == q1)){
					found++;
				}else if((nq0 == q1) && (nq1 == q0)){
					found++;
				}
			}
			assert(found == 1);
			if(found != 1) return false;
		}
	}

	if(depth > 0)
		for(int i = FC3D_FIRST; i <= FC3D_LAST; i++)
			if(m_extra->m_neighbours[i] && !m_extra->m_neighbours[i]->valid(depth-1))
				return false;

	return true;
}

void ControlSpace3dOctree::compact()
{
	// remove obsolete leaves (and extra data)
	int leaves_total = 0;
	int leaves_obsolete = 0;

	int count = m_n[0]*m_n[1]*m_n[2];
	for(int i = 0; i < count; i++)
		m_grid[i].optimize(leaves_total, leaves_obsolete);

	// scan through list of vertices and delete vertices with ref=1

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "ACS::Octree [total] " << leaves_total);
	LOG4CPLUS_DEBUG(MeshLog::logger_console, "ACS::Octree [obs  ] " << leaves_obsolete);

	count = m_grid_vertices.countInt();
	LOG4CPLUS_DEBUG(MeshLog::logger_console, "ACS::Octree [v total] " << count);
	for(int i = 0; i < count; )
		if(m_grid_vertices[i].use_count() == 1){
			m_grid_vertices.removeAt(i);
			--count;
		}else ++i;

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "ACS::Octree [v left ] " << count);

	// done
	ControlSpace3dAdaptive::compact(); // set readonly mode
}

bool ControlSpace3dOctree::OctLeaf::optimize(int & leaves_total, int & leaves_obsolete)
{
	bool no_subleaves = true;
	if(isSplit()){
		for(int i = 0; i < 8; i++)
			no_subleaves &= m_leaves[i]->optimize(leaves_total, leaves_obsolete);
		if(no_subleaves && param_octree_reduce_threshold > 0.0){
			// check middle point only
			// ... interpolated
			ControlDataMatrix3d cdm_middle = m_vertices[0]->control_data;
			for(int i = 1; i < 8; i++)
				cdm_middle += m_vertices[i]->control_data;
			cdm_middle *= (1.0/8.0);
			// ... from middle point
			double diff = cdm_middle.countDifferenceRR(
				m_leaves[0]->m_vertices[VX3D_HNE]->control_data);
			if(diff < param_octree_reduce_threshold * ControlSpace2dAdaptive::param_threshold_diff){ 
				leaves_obsolete += 8;
				for(int i = 0; i < 8; i++) delete m_leaves[i];
				m_leaves[0] = nullptr;
			}else{
				no_subleaves = false;
			}
			// check middle point + midedges
//			/*
//			*/
		}
	}	
	++leaves_total;
	delete m_extra;
	m_extra = nullptr;

	return no_subleaves;
}

void ControlSpace3dOctree::testOctree()
{
/*
	storeEPS("control", 0);

	DPoint2d test_point(0.001, 0.001);

	for(int j = 1; j < 10; j++){
		int k = countLocalCoordinates(test_point);
		OctLeaf* leaf = m_grid[k].findLastLeaf(test_point);
		OctVertexList split_vertices(100);
		leaf->split(m_grid_vertices, split_vertices, base_surface);
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Control OctTree split -> " << split_vertices.countInt() << 
			" new vertices." << endl;
		storeEPS("control", j);
	}
*/
}

void ControlSpace3dOctree::OctLeaf::gatherDataFromSourcePoints()
{
	if(isSplit()){
		assert(m_extra->m_source_points == nullptr);
		for(int i = VX3D_FIRST; i <= VX3D_LAST; i++)
			m_leaves[i]->gatherDataFromSourcePoints();
	}else{
		if(m_extra->m_source_points){
			int sct = m_extra->m_source_points->countInt();
			for(int i = 0; i < sct; i++){
				auto qv = m_extra->m_source_points->get(i);
				for(int j = VX3D_FIRST; j <= VX3D_LAST; j++){
					double dist = qv->coord.distance2(m_vertices[j]->coord);
					double nw = (dist < mesh_data.relative_small_number)?mesh_data.relative_infinity:(1.0/dist);
					m_vertices[j]->control_data += qv->control_data * nw;
					m_vertices[j]->w += nw;
				}
			}
			delete m_extra->m_source_points;
			m_extra->m_source_points = nullptr;
		}
	}
}

/// Propagate control values from within leaves
int ControlSpace3dOctree::OctLeaf::propagateUp()
{
	int ct = 0;
	if(isSplit()){
		for(int i = VX3D_FIRST; i <= VX3D_LAST; i++)
			ct += m_leaves[i]->propagateUp();
	}
	for(int i = VX3D_FIRST; i <= VX3D_LAST; i++){
		if(m_vertices[i]->w == 0.0){
			// TODO improve?
			for(int j = 0; j < 3; j++){
				auto cn = m_extra->m_midedge_vertices[DBox::vertex_to_edge[i][j]];
				if(cn && cn->w < 0.0){ // i.e. already initialized
					m_vertices[i]->control_data += cn->control_data;
					m_vertices[i]->w += 1.0;
				}
			}
			for(int j = 0; j < 3; j++){
				auto cn = m_extra->m_midface_vertices[DBox::vertex_to_face[i][j]];
				if(cn && cn->w < 0.0){ // i.e. already initialized
					m_vertices[i]->control_data += cn->control_data;
					m_vertices[i]->w += 1.0;
				}
			}
			if(m_vertices[i]->w > 1.0)
				m_vertices[i]->control_data /= m_vertices[i]->w;
			if(m_vertices[i]->w > 0.0){
				m_vertices[i]->w = -2.0;
				++ct;
			}
		}
	}
	return ct;
}

/// Propagate control values into leaves
int ControlSpace3dOctree::OctLeaf::propagateDown()
{
	int ct = 0;
	for(int i = ED3D_FIRST; i <= ED3D_LAST; i++){	// mid-vertices of edges
		if(m_extra->m_midedge_vertices[i] && (m_extra->m_midedge_vertices[i]->w == 0.0)){
			m_extra->m_midedge_vertices[i]->control_data = 
				(m_vertices[DBox::edge_to_vertex[i][0]]->control_data +
				 m_vertices[DBox::edge_to_vertex[i][1]]->control_data) * 0.5;
			m_extra->m_midedge_vertices[i]->w = -4.0;
			++ct;
		}
	}
	for(int i = FC3D_FIRST; i <= FC3D_LAST; i++){	// mid-vertices of faces
		if(m_extra->m_midface_vertices[i] && (m_extra->m_midface_vertices[i]->w == 0.0)){
			m_extra->m_midface_vertices[i]->control_data = 
				(m_vertices[DBox::face_to_vertex[i][0]]->control_data +
				 m_vertices[DBox::face_to_vertex[i][1]]->control_data +
				 m_vertices[DBox::face_to_vertex[i][2]]->control_data +
				 m_vertices[DBox::face_to_vertex[i][3]]->control_data) * 0.25;
			m_extra->m_midface_vertices[i]->w = -4.0;
			++ct;
		}
	}
	if(isSplit()){
		// middle-vertex
		auto node = m_leaves[VX3D_LSW]->m_vertices[VX3D_HNE];
		if(node->w == 0.0){
			for(int i = VX3D_FIRST; i <= VX3D_LAST; i++)
				node->control_data += m_vertices[i]->control_data;
			node->control_data *= (1.0/(VX3D_LAST+1));
			node->w = -4.0;
			++ct;
		}
		// and check lower
		for(int i = VX3D_FIRST; i <= VX3D_LAST; i++)
			ct += m_leaves[i]->propagateDown();
	}
	return ct;
}

void ControlSpace3dOctree::OctLeaf::statOctree(int & max_level, int & leaf_counter, int level) const
{
	if(isSplit()){
		for(int i = 0; i < 8; i++)
			m_leaves[i]->statOctree(max_level, leaf_counter, level+1);
	}else{
		++leaf_counter;
		if(level > max_level) max_level = level;
	}
}

/// Log basic information about this control space
void ControlSpace3dOctree::logDescription() const
{
	int max_level = 0;
	int leaf_counter = 0;
	const int main_leaf_count = m_n[0] * m_n[1] * m_n[2];
	for(int i = 0; i < main_leaf_count; i++)
		m_grid[i].statOctree(max_level, leaf_counter);

	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"Octree [" << m_n[0] << "," << m_n[1] << "," << m_n[2] << "], nodes= " <<
		m_grid_vertices.countInt() << ", max depth= " << 
		max_level << ", leaves= " << leaf_counter);
/*
	ofstream file("qtree-test.txt");
	const DPoint2d pt0 = m_box.getLeftBottom();
	const DPoint2d pt1 = m_box.getLeftTop();
	const DPoint2d up(0.0, 1.0);
	for(double t = 0.0; t <= 0.5; t += 0.002){
		const DPoint2d pt = pt0*(1.0-t) + pt1*t;
		const ControlDataMatrix2d cdm = findLastLeaf(pt)->getMetricAtPoint(pt);
		double len = (cdm * up).length();
		file << pt.y << '\t' << len << endl;
	}
*/
}

/// Smoothen variance of metric within the control space
bool ControlSpace3dOctree::smoothen()
{
	assert(m_initialized == 1); // with parameterization matrix
	// uptree
	int leaf_count = m_n[0] * m_n[1] * m_n[2];
	int up_count = 0;
	for(int i = 0; i < leaf_count; i++){
		up_count += m_grid[i].smoothenMetricAtNodes(true);
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"ControlSpace3dOctree smoothing (up) - " << up_count << " modifications.");
	// at the top level ??
	/*
	int k = 0;
	for(int iy = 0; iy <= m_ny; iy++){
		for(int ix = 0; ix <= m_nx; ix++){
			ControlNode2d& node = m_grid_vertices.getDataAt(k);
		}
	}
	*/
	// downtree
	int down_count = 0;
	for(int i = 0; i < leaf_count; i++){
		down_count += m_grid[i].smoothenMetricAtNodes(false);
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"ControlSpace3dOctree smoothing (down) - " << down_count << " modifications.");

	return up_count+down_count > 0;
}

int ControlSpace3dOctree::OctLeaf::smoothenMetricAtNodes(bool uptree)
{
	int count = 0;
	if(uptree && isSplit()){	// first leaves then this node
		for(int i = VX3D_FIRST; i <= VX3D_LAST; i++)
			count += m_leaves[i]->smoothenMetricAtNodes(uptree);
	}
	// check along leaf-edges
	const DVector3d dv[3] = { DVector3d(2*m_dl.x, 0.0, 0.0),
		DVector3d(0.0, 2*m_dl.y, 0.0), DVector3d(0.0, 0.0, 2*m_dl.z) };

	for(int i = ED3D_FIRST; i <= ED3D_LAST; i++){
		int v0 = DBox::edge_to_vertex[i][0];
		int v1 = DBox::edge_to_vertex[i][1];

		int dir = 0;
		switch(i){
			case ED3D_LS:
			case ED3D_LN:
			case ED3D_HS:
			case ED3D_HN:
				dir = 0; break; // along x
			case ED3D_LW:
			case ED3D_LE:
			case ED3D_HW:
			case ED3D_HE:
				dir = 1; break; // along y
			case ED3D_SW:
			case ED3D_SE:
			case ED3D_NW:
			case ED3D_NE:
				dir = 2; break; // along z
			default:
				assert(false);
		}
		int result = ControlSpace3dAdaptive::smoothenMetricForNodes(
			m_vertices[v0].get(), m_vertices[v1].get(), dv[dir]);

		if((result & 1) != 0) ++count;
		if((result & 2) != 0) ++count;
	}

	if(!uptree && isSplit()){	// first this node then leaves
		for(int i = VX3D_FIRST; i <= VX3D_LAST; i++)
			count += m_leaves[i]->smoothenMetricAtNodes(uptree);
	}

	return count;
}

void ControlSpace3dOctree::OctLeaf::setSingleVertex(int side, CN3dPtr vert)
{
#ifdef _DEBUG
	for (int i = VX3D_FIRST; i <= VX3D_LAST; i++)
		if (i != side)
			assert(m_vertices[i] != vert);
#endif // _DEBUG

	m_vertices[side] = vert;
}

double ControlSpace3dOctree::getMetricGradationRatio(const DPoint3d& pt) const
{
	assert(m_initialized>0);
	return findLastLeaf(pt)->getMetricGradationRatio(pt);
}

/// Store octree grid to text file
bool ControlSpace3dOctree::storeTXT(const char* fname)
{
	ofstream ofs(fname);
	if(!ofs){
		LOG4CPLUS_ERROR(MeshLog::logger_console, "Couldn't open [" << fname << "] file for write!");
		return false;
	}

	ofs << m_box << endl;
	ofs << "ControlSpace3dOctree " << m_n[0] << " " << m_n[1] << " " << m_n[2];
	ofs << " " << (m_initialized < 10 ? 0 : 1) << endl; // 1 -> readonly (i.e. no balancing)
	int cn_count = m_grid_vertices.countInt();
	ofs << "ControlNodes " << cn_count << endl;
	int main_leaves_count = m_n[0]*m_n[1]*m_n[2];

	for(int i = 0; i < cn_count; i++){
		auto node = m_grid_vertices[i];
		ofs << node->coord << "\t";
		ofs << node->control_data << "\t";
		ofs << node->max_gradation_ratio << endl;
		node->w = i;
	}
	ofs << "ControlLeaves" << endl;
	for(int i = 0; i < main_leaves_count; i++){
		m_grid[i].storeTxtSplit(ofs);
		ofs << endl;
	}

	for(int i = 0; i < cn_count; i++)
		m_grid_vertices[i]->w = -1.0;

	return true;
}

void ControlSpace3dOctree::OctLeaf::storeTxtSplit(ostream& os) const
{
	os << 'v';
	for(int i = 0; i < 8; i++) os << (int)m_vertices[i]->w << ',';
	if(isSplit()){
		os << '1';
		for(int i = 0; i < 8; i++) m_leaves[i]->storeTxtSplit(os);
	}else{
		os << '0';
	}
}

/// Load octree grid from text file
bool ControlSpace3dOctree::loadTXT(const char* fname)
{
	ifstream ifs(fname);
	if(!ifs){
		LOG4CPLUS_ERROR(MeshLog::logger_console, 
			"Couldn't open [" << fname << "] file for read!");
		return false;
	}

	DBox box;
	string buffer, text;
	getline(ifs, buffer);
	if(ifs){
		istringstream str_file(buffer);
		str_file >> box;
		if(!str_file){
			LOG4CPLUS_ERROR(MeshLog::logger_console,   
				"Couldn't read box dimensions");
			return false;
		}
	}
	if(box.x0 > m_box.x0 + SMALL_NUMBER ||
		box.x1 < m_box.x1 - SMALL_NUMBER ||
		box.y0 > m_box.y0 + SMALL_NUMBER ||
		box.y1 < m_box.y1 - SMALL_NUMBER ||
		box.z0 > m_box.z0 + SMALL_NUMBER ||
		box.z1 < m_box.z1 - SMALL_NUMBER)
	{
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"Loaded box dimensions don't quite contain oryginal box!");
		//return false;
	}
	m_box = box;

	if(m_grid) delete[] m_grid;
	getline(ifs, buffer);
	int readonly;
	if(ifs){
		istringstream str_file(buffer);
		str_file >> text >> m_n[0] >> m_n[1] >> m_n[2];
		if(!str_file || text != "ControlSpace3dOctree"){
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"CS3dOctree loading: syntax error -> " << buffer);
			return false;
		}
		str_file >> readonly;
		if(!str_file) readonly = 0;
	}

	const int nxyz = m_n[0] * m_n[1] * m_n[2]; // number of leaves for initial grid

	// load control nodes
	m_grid_vertices.clear();
	int cn_count = 0;
	getline(ifs, buffer);
	if(ifs){
		istringstream str_file(buffer);
		str_file >> text >> cn_count;
		if(!str_file || text != "ControlNodes"){
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"CS3dOctree loading: syntax error -> " << buffer);
			return false;
		}
	}
	ControlNode3d cn;
	for(int i = 0; i < cn_count; i++){
		getline(ifs, buffer);
		istringstream str_file(buffer);
		str_file >> cn.coord >> cn.control_data >> cn.max_gradation_ratio;
		if(ifs && str_file){
			auto scn = std::make_shared<ControlNode3d>(cn.coord, cn.control_data);
			scn->max_gradation_ratio = cn.max_gradation_ratio;
			m_grid_vertices.add(scn);
		}else{
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"CS3dOctree loading: syntax error -> " << buffer);
			return false;
		}
	}

	// like in constructor
	double len[3] = { box.getDX(), box.getDY(), box.getDZ() };
//	double total_len = len[0]+len[1]+len[2];
//	double fl[3] = { len[0] / total_len, len[1] / total_len, len[2] / total_len };

//	double ratio = pow(nxyz / (fl[0]*fl[1]*fl[2]), 1.0/3.0);
//	for(int i = 0; i < 3; i++)
//		m_n[i] = std::min(std::max(2,(int)(fl[i]*ratio)), nxyz/3);

	double d[3] = { len[0]/m_n[0], len[1]/m_n[1], len[2]/m_n[2] };

	m_grid = new OctLeaf[nxyz];
	
	// init regular grid boxes with coordinates and vertices indices
	const DVector3d dpt(d[0]/4, d[1]/4, d[2]/4);
//	const int nxy = m_n[0]*m_n[1];
	const int nx1 = m_n[0]+1;
//	const int nxy1 = (m_n[0]+1)*(m_n[1]+1);
	DPoint3d pt(m_box.x0 + d[0]/2, m_box.y0 + d[1]/2, m_box.z0 + d[2]/2);
	for(int iz = 0, k = 0, in = 0; iz < m_n[2]; iz++, in+=nx1){
		for(int iy = 0; iy < m_n[1]; iy++, in++){
			for(int ix = 0; ix < m_n[0]; ix++, in++, k++){
				m_grid[k].init(pt, dpt);
				pt.x += d[0];
			}
			pt.y += d[1];
			pt.x = m_box.x0 + d[0]/2;
		}
		pt.z += d[2];
		pt.y = m_box.y0 + d[1]/2;
		pt.x = m_box.x0 + d[0]/2;
	}


	// split leaves
	getline(ifs, buffer);
	if(ifs){
		istringstream str_file(buffer);
		str_file >> text;
		if(!str_file || text != "ControlLeaves"){
			LOG4CPLUS_ERROR(MeshLog::logger_console, 
				"CS3dOctree loading: syntax error -> " << buffer);
			return false;
		}
	}
	for(int i = 0; i < nxyz; i++){
		getline(ifs, buffer);
		if(ifs){
			istringstream str_file(buffer);
			m_grid[i].loadTxtSplit(str_file, m_grid_vertices);

		}
	}

	m_initialized = 10; // always set as readonly

	LOG4CPLUS_DEBUG(MeshLog::logger_console, 
		"CS3dOctree loaded as readonly, valid: " << valid());

	return true;
}

void ControlSpace3dOctree::OctLeaf::loadTxtSplit(istream& is, OctVertexList& grid_vertices)
{
	char z;
	int index;
	// load vertices
	is >> z; assert(z == 'v');
	for(int i = 0; i < 8; i++){
		is >> index >> z; assert(z == ',');
		m_vertices[i] = grid_vertices[index];
	}

	// readonly ?
	delete m_extra;
	m_extra = nullptr;

	// load split
	is >> z;
	if(z == '1'){ // should be split
		if(!isSplit()){
			const DVector3d mid_dl = m_dl/2;
			m_leaves[VX3D_LSW] = new OctLeaf(DPoint3d(m_middle.x-m_dl.x, m_middle.y-m_dl.y, 
				m_middle.z-m_dl.z), mid_dl);
			m_leaves[VX3D_LSE] = new OctLeaf(DPoint3d(m_middle.x+m_dl.x, m_middle.y-m_dl.y, 
				m_middle.z-m_dl.z), mid_dl);
			m_leaves[VX3D_LNW] = new OctLeaf(DPoint3d(m_middle.x-m_dl.x, m_middle.y+m_dl.y, 
				m_middle.z-m_dl.z), mid_dl);
			m_leaves[VX3D_LNE] = new OctLeaf(DPoint3d(m_middle.x+m_dl.x, m_middle.y+m_dl.y, 
				m_middle.z-m_dl.z), mid_dl);
			m_leaves[VX3D_HSW] = new OctLeaf(DPoint3d(m_middle.x-m_dl.x, m_middle.y-m_dl.y, 
				m_middle.z+m_dl.z), mid_dl);
			m_leaves[VX3D_HSE] = new OctLeaf(DPoint3d(m_middle.x+m_dl.x, m_middle.y-m_dl.y, 
				m_middle.z+m_dl.z), mid_dl);
			m_leaves[VX3D_HNW] = new OctLeaf(DPoint3d(m_middle.x-m_dl.x, m_middle.y+m_dl.y, 
				m_middle.z+m_dl.z), mid_dl);
			m_leaves[VX3D_HNE] = new OctLeaf(DPoint3d(m_middle.x+m_dl.x, m_middle.y+m_dl.y, 
				m_middle.z+m_dl.z), mid_dl);
		}
		for(int i = 0; i < 8; i++)
			m_leaves[i]->loadTxtSplit(is, grid_vertices);
	}else if(z == '0'){
		assert(!isSplit());
	}else{
		assert(false);
	}
}

/// Calculate inside information for nodes/elements
void ControlSpace3dOctree::markInsideNodes(const MeshContainer3dSurface* surface_mesh)
{
	START_CLOCK("CS3dOT::markInsideNodes");
	//MeshViewSet* set = new MeshViewSet;
	START_CLOCK("CS3dOT::markInsideNodes-1");

	double dx = m_box.getDX() / m_n[0];
	double dy = m_box.getDY() / m_n[1];
	double dz = m_box.getDZ() / m_n[2];

	const double EPSILON = 0.01;
	int nxy = m_n[0] * m_n[1];

	//int counter = 0;

	int fct = surface_mesh->getFacesCount();
	for(int i = 0; i < fct; i++){
		MeshFace* face = surface_mesh->getFaceAt(i);

		const DPoint3d pt0 = m_box.fitInPoint(face->getPoint(0)->getCoordinates());
		const DPoint3d pt1 = m_box.fitInPoint(face->getPoint(1)->getCoordinates());
		const DPoint3d pt2 = m_box.fitInPoint(face->getPoint(2)->getCoordinates());

		double d01 = pt0.distance(pt1);
		double d12 = pt1.distance(pt2);
		double d20 = pt2.distance(pt0);

		double diameter = std::max(std::max(d01, d12), d20);
		const double MIN_DL = EPSILON * diameter;

		DBox face_box;
		face_box.addPoint(pt0);
		face_box.addPoint(pt1);
		face_box.addPoint(pt2);

		double xmin = (face_box.x0 - m_box.x0) / dx - EPSILON; // x
		double xmax = (face_box.x1 - m_box.x0) / dx + EPSILON;
		double ymin = (face_box.y0 - m_box.y0) / dy - EPSILON; // y
		double ymax = (face_box.y1 - m_box.y0) / dy + EPSILON;
		double zmin = (face_box.z0 - m_box.z0) / dz - EPSILON; // z
		double zmax = (face_box.z1 - m_box.z0) / dz + EPSILON;

		int x0_grid = (int)(xmin+0.5);
		int y0_grid = (int)(ymin+0.5);
		int z0_grid = (int)(zmin+0.5);

		int ixmin = std::max((int)xmin, 0);
		int ixmax = std::min((int)xmax, m_n[0]-1);
		int iymin = std::max((int)ymin, 0);
		int iymax = std::min((int)ymax, m_n[1]-1);
		int izmin = std::max((int)zmin, 0);
		int izmax = std::min((int)zmax, m_n[2]-1);

		// debug - only
		// ... temporary clean all nodes
		//for(int i = 0; i < m_grid_vertices.countInt(); i++){
		//	ControlNode3d* node = m_grid_vertices[i];
		//	node->removeTag(TagExtended::TAG_INSIDE_AREA);
		//	node->removeTag(TagExtended::TAG_INSIDE_AREA_DIST);
		//}

		// debug only find non-OXY faces close to control nodes (any node)
		//DataVector<ControlNode3d*> close_nodes;
		//if(face_box.getDZ() > MIN_DL){
		//	DTriangle3d face_triangle(pt0, pt1, pt2);
		//	const double MIN_DL2 = MIN_DL * MIN_DL;
		//	for(int i = 0; i < m_grid_vertices.countInt(); i++){
		//		if(face_triangle.distance2ToPoint(m_grid_vertices[i]->coord) < MIN_DL2)
		//			close_nodes.add(m_grid_vertices[i]);
		//	}
		//	if(close_nodes.countInt() > 0){
		//		MeshViewSet* set = new MeshViewSet();
		//		set->addFace(it.getFace());
		//		for(int i = 0; i < close_nodes.countInt(); i++){
		//			int v = close_nodes[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -2);
		//			if(v > -1) set->addPoint(close_nodes[i]->coord, v, v);
		//			else set->addPoint(close_nodes[i]->coord);
		//		}
		//		SHOW_MESH("near-face", set);
		//	}
		//}

		if(face_box.getDX() < MIN_DL){

			//MeshViewSet* set = new MeshViewSet();
			//set->addFace(it.getFace());

			//---------
			if(abs(0.5*(xmin+xmax) - x0_grid) < EPSILON){ // along YZ grid plane
				for(int iy = iymin; iy <= iymax; iy++){
					for(int iz = izmin; iz <= izmax; iz++){
						int k = iz * nxy + iy * m_n[0];
						if(x0_grid > 0){
							m_grid[k+x0_grid-1].markInsideNodesYZ(pt0, pt1, pt2, face, true, MIN_DL); // plane_x1
							//for(int i = 0; i < 8; i++) set->addPoint(m_grid[k+x0_grid-1].m_vertices[i]->coord);
						}if(x0_grid < m_n[0]){
							m_grid[k+x0_grid].markInsideNodesYZ(pt0, pt1, pt2, face, false, MIN_DL); // plane_x0
							//for(int i = 0; i < 8; i++) set->addPoint(m_grid[k+x0_grid].m_vertices[i]->coord);
						}
					}
				}
				//SHOW_MESH("leaf along CS grid-plane YZ", set);
			}else{ // along one YZ plane-stripe
				for(int iy = iymin; iy <= iymax; iy++){
					for(int iz = izmin; iz <= izmax; iz++){
						int k = iz * nxy + iy * m_n[0] + ixmin;
						m_grid[k].markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
						//for(int i = 0; i < 8; i++) set->addPoint(m_grid[k].m_vertices[i]->coord);
					}
				}
				//SHOW_MESH("leaf between CS grid-planes YZ", set);
			}
		}else if(face_box.getDY() < MIN_DL){
			//---------
			if(abs(0.5*(ymin+ymax) - y0_grid) < EPSILON){ // along XZ grid plane
				for(int ix = ixmin; ix <= ixmax; ix++){
					for(int iz = izmin; iz <= izmax; iz++){
						int k = iz * nxy + ix;
						if(y0_grid > 0)
							m_grid[k+(y0_grid-1)*m_n[0]].markInsideNodesXZ(pt0, pt1, pt2, face, true, MIN_DL); // plane_y1
						if(y0_grid < m_n[1])
							m_grid[k+y0_grid*m_n[0]].markInsideNodesXZ(pt0, pt1, pt2, face, false, MIN_DL); // plane_y0
					}
				}
			}else{ // along one XZ plane-stripe
				for(int ix = ixmin; ix <= ixmax; ix++){
					for(int iz = izmin; iz <= izmax; iz++){
						int k = iz * nxy + iymin * m_n[0] + ix;
						m_grid[k].markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
					}
				}
			}
		}else if(face_box.getDZ() < MIN_DL){
			//---------
			if(abs(0.5*(zmin+zmax) - z0_grid) < EPSILON){ // along XY grid plane
				for(int ix = ixmin; ix <= ixmax; ix++){
					for(int iy = iymin; iy <= iymax; iy++){
						int k = iy * m_n[0] + ix;
						if(z0_grid > 0)
							m_grid[k+(z0_grid-1)*nxy].markInsideNodesXY(pt0, pt1, pt2, face, true, MIN_DL); // plane_z1
						if(z0_grid < m_n[2])
							m_grid[k+z0_grid*nxy].markInsideNodesXY(pt0, pt1, pt2, face, false, MIN_DL); // plane_z0
					}
				}
			}else{ // along one XY plane-stripe
				for(int ix = ixmin; ix <= ixmax; ix++){
					for(int iy = iymin; iy <= iymax; iy++){
						int k = izmin * nxy + iy * m_n[0] + ix;
						m_grid[k].markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
					}
				}
			}
		}else{ // several X strips
			//--
			//if(xmax-0.1 > 0.0) continue;
			//--
			//it.getFace()->setIntTag(TagExtended::TAG_INSIDE_AREA, 1); // just for debug-view

			//MeshViewSet* set = new MeshViewSet();
			//set->addFace(it.getFace());
			//---------
			for(int iz = izmin; iz <= izmax; iz++){
				for(int iy = iymin; iy <= iymax; iy++){
					int k = iz * nxy + iy * m_n[0];
					for(int ix = ixmin; ix <= ixmax; ix++){
						m_grid[k + ix].markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
						//for(int i = 0; i < 8; i++) set->addPoint(m_grid[k+ix].m_vertices[i]->coord);
					}
				}
			}
			//SHOW_MESH("CS cross-leaf", set);
		}

		//if(close_nodes.countInt() > 0){
		//	MeshViewSet* set = new MeshViewSet();
		//	set->addFace(it.getFace());
		//	for(int i = 0; i < close_nodes.countInt(); i++){
		//		int v = close_nodes[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -2);
		//		if(v > -1) set->addPoint(close_nodes[i]->coord, v, v);
		//		else set->addPoint(close_nodes[i]->coord);
		//	}
		//	SHOW_MESH("near-face after", set);
		//}

		//LOG4CPLUS_INFO(MeshLog::logger_console, "Inside-nodes after face number", counter);
		//showMarkInside(boundary_mesh);
		//++counter;
	}

	STOP_CLOCK("CS3dOT::markInsideNodes-1");

//	showMarkInside(bfaces);

	// mark boundary CS-nodes as 0 unless already set to "inside"
	START_CLOCK("CS3dOT::markInsideNodes-1a");

	// switch "-1" into "no tag"
	int cct = m_grid_vertices.countInt();
	for(int i = 0; i < cct; i++){
		if(m_grid_vertices[i]->getIntTag(TagExtended::TAG_INSIDE_AREA) == -1)
			m_grid_vertices[i]->removeTag(TagExtended::TAG_INSIDE_AREA);
	}
	// mark at boundary
	markInsideNodesAtBoundary();

	STOP_CLOCK("CS3dOT::markInsideNodes-1a");

//	showMarkInside(bfaces);

	START_CLOCK("CS3dOT::markInsideNodes-2");

	int nxyz = nxy * m_n[2];

	// propagate control values from within leaves
	DataVector<OctLeaf*> not_ready(nxyz);
	for(int i = 0; i < nxyz; i++){
		if(!m_grid[i].propagateInsideMark())
			not_ready.add(m_grid+i);
	}

//	showMarkInside(boundary_mesh);

	if(!not_ready.empty()){
		// extrapolate to other top-grid nodes
		propagateMainGridInsideMark();
//		showMarkInside(bfaces);

		int last_not_ready_count = (int)not_ready.countInt();
		while(!not_ready.empty()){
			// propagate control values into leaves			
			for(int i = 0; i < not_ready.countInt(); ){
				if(not_ready[i]->propagateInsideMark())
					not_ready.removeAt(i);
				else ++i;
			}
//			showMarkInside(bfaces);
			if(last_not_ready_count == not_ready.countInt()){
				// nothing done -> something wrong
				//	... set rest to 1
				for(int i = 0; i < not_ready.countInt(); i++){
					not_ready[i]->setInsideMark(1);
				}
				not_ready.clear();
			}else{
				// else update
				last_not_ready_count = (int)not_ready.countInt();
			}
		}
	}

	// clear forbidden-tag
	for(int i = 0; i < m_grid_vertices.countInt(); i++){
		m_grid_vertices[i]->removeTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN);
	}

	STOP_CLOCK("CS3dOT::markInsideNodes-2");

//	showMarkInside(bfaces);

	START_CLOCK("CS3dOT::markInsideNodes-3");
	// calculate prediction
	ControlPrediction cp = calculateControlPrediction();

	STOP_CLOCK("CS3dOT::markInsideNodes-3");
	STOP_CLOCK("CS3dOT::markInsideNodes");

	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"========== Control Prediction ==============");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		" * Leaves [0-1-2-3-4-5-6-7-8]: [" 
		<< cp.counter[0] << ","
		<< cp.counter[1] << "," << cp.counter[2] << ","
		<< cp.counter[3] << "," << cp.counter[4] << ","
		<< cp.counter[5] << "," << cp.counter[6] << ","
		<< cp.counter[7] << "," << cp.counter[8] << "]");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"CS-prediction: " << cp.value[2] << " (" << cp.value[0] << " - " << cp.value[1] << ")");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"============================================");

	LOG4CPLUS_INFO(MeshLog::logger_console, 
		"Predicted (CS) number of tetrahedra: " << (int)cp.value[2]);

//	SHOW_MESH("control space (markInsideNodes + extend)", set);
}

ControlSpace3dOctree::ControlPrediction ControlSpace3dOctree::calculateControlPrediction() const
{
	ControlPrediction cp;
	START_CLOCK("CS3dOT::calculateCSPrediction");
	int nxyz = m_n[0] * m_n[1] * m_n[2];
	for(int i = 0; i < nxyz; i++)
		m_grid[i].calculateControlPrediction(cp);
	STOP_CLOCK("CS3dOT::calculateCSPrediction");

	return cp;
}

/// Calculate mesh prediction using CS
void ControlSpace3dOctree::OctLeaf::calculateControlPrediction(ControlSpace3dOctree::ControlPrediction & cp) const
{
	if(isSplit()){
		for(int i = 0; i < 8; i++)
			m_leaves[i]->calculateControlPrediction(cp);
	}else{
		int inside_count = 0;
		int all_nodes[] = {0, 1, 2, 3, 4, 5, 6, 7};
		int inside_nodes[8];
		for(int i = 0; i < 8; i++){
			if(m_vertices[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1) > 0)
				inside_nodes[inside_count++] = i;
		}
		cp.counter[inside_count]++;
		if(inside_count == 8) { // for min-approx -> all nodes inside
			calculateControlPrediction(cp, 0, 8, all_nodes);
		}
		if(inside_count > 0){ // for max-approx and ave-approx -> at least one node inside
			calculateControlPrediction(cp, 1, 8, all_nodes);
			calculateControlPrediction(cp, 2, inside_count, inside_nodes);
		}
	}
}

/// Calculate mesh prediction using CS
void ControlSpace3dOctree::OctLeaf::calculateControlPrediction(
	ControlSpace3dOctree::ControlPrediction & cp, int ind, int inside_count, int inside[]) const
{
	if(inside_count < 1) return;

	ControlDataMatrix3d cdm_ave_arytm = m_vertices[inside[0]]->control_data;
	assert(!m_vertices[inside[0]]->gradationUnknown());
	double gradation_ave = m_vertices[inside[0]]->max_gradation_ratio;
	for(int i = 1; i < inside_count; i++){
		cdm_ave_arytm += m_vertices[inside[i]]->control_data;
		assert(!m_vertices[inside[i]]->gradationUnknown());
		gradation_ave += m_vertices[inside[i]]->max_gradation_ratio;
	}
	cdm_ave_arytm *= 1.0 / inside_count;
	gradation_ave *= 1.0 / inside_count;

	gradation_ave = 1.0 + 0.25 * (gradation_ave - 1.0);

	double f = 0.125 * inside_count * gradation_ave;

	DMetric3d dm_ave_arytm(cdm_ave_arytm);

	static const double TETRA_VOLUME = SQRT2 / 12.0;
/*
	const DVector3d vx(4*m_dl.x, 0.0, 0.0);
	const DVector3d vy(0.0, 4*m_dl.y, 0.0);
	const DVector3d vz(0.0, 0.0, 4*m_dl.z);

	double volume = 
		dm_ave_arytm.transformRStoMS(vx).length() * 
		dm_ave_arytm.transformRStoMS(vy).length() * 
		dm_ave_arytm.transformRStoMS(vz).length();
*/

	const DMPoint3d dmpoints[8] = {
		DMPoint3d(0.0, 0.0, 0.0),
		dm_ave_arytm.transformRStoMS(DPoint3d(4*m_dl.x, 0.0,	  0.0)),
		dm_ave_arytm.transformRStoMS(DPoint3d(4*m_dl.x, 0.0,	  4*m_dl.z)),
		dm_ave_arytm.transformRStoMS(DPoint3d(0.0,		0.0,	  4*m_dl.z)),
		dm_ave_arytm.transformRStoMS(DPoint3d(0.0,		4*m_dl.y, 0.0)),
		dm_ave_arytm.transformRStoMS(DPoint3d(4*m_dl.x, 4*m_dl.y, 0.0)),
		dm_ave_arytm.transformRStoMS(DPoint3d(4*m_dl.x, 4*m_dl.y, 4*m_dl.z)),
		dm_ave_arytm.transformRStoMS(DPoint3d(0.0,		4*m_dl.y, 4*m_dl.z))
	};
	// Five tetrahedra for cube
	const int conf[5][4] = {{0, 1, 2, 5}, {0, 2, 3, 7}, {0, 5, 7, 4}, {2, 5, 6, 7}, {0, 5, 2, 7}};
	double volume = 0.0;
	for(int i = 0; i < 5; i++)
		volume += DTetrahedron::volume(
				dmpoints[conf[i][0]], dmpoints[conf[i][1]],
				dmpoints[conf[i][2]], dmpoints[conf[i][3]]);

	cp.value[ind] += f * volume / TETRA_VOLUME;
}

int ControlSpace3dOctree::OctLeaf::markInsideNodes(
	const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, 
	MeshFace* face, const double& MIN_DL)
{
	const double EPSILON = 0.01;
	const int FEDGES3D[12][5] = {
		//  0       1         2               3                  4
		{ VX3D_LSW, VX3D_LSE, EDGE_LSWE, FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL },
		{ VX3D_LNW, VX3D_LNE, EDGE_LNWE, FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL },
		{ VX3D_HSW, VX3D_HSE, EDGE_HSWE, FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL },
		{ VX3D_HNW, VX3D_HNE, EDGE_HNWE, FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL },
		{ VX3D_LSW, VX3D_LNW, EDGE_LSNW, FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL },
		{ VX3D_LSE, VX3D_LNE, EDGE_LSNE, FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL },
		{ VX3D_HSW, VX3D_HNW, EDGE_HSNW, FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL },
		{ VX3D_HSE, VX3D_HNE, EDGE_HSNE, FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL },
		{ VX3D_LSW, VX3D_HSW, EDGE_LHSW, FORBIDDEN_EDGE_ZH, FORBIDDEN_EDGE_ZL },
		{ VX3D_LSE, VX3D_HSE, EDGE_LHSE, FORBIDDEN_EDGE_ZH, FORBIDDEN_EDGE_ZL },
		{ VX3D_LNW, VX3D_HNW, EDGE_LHNW, FORBIDDEN_EDGE_ZH, FORBIDDEN_EDGE_ZL },
		{ VX3D_LNE, VX3D_HNE, EDGE_LHNE, FORBIDDEN_EDGE_ZH, FORBIDDEN_EDGE_ZL },
	};

	int res = 0;

	if(isSplit()){
		// recursive
		double ix0 = (pt0.x - m_middle.x) / m_dl.x;
		double ix1 = (pt1.x - m_middle.x) / m_dl.x;
		double ix2 = (pt2.x - m_middle.x) / m_dl.x;
		double iy0 = (pt0.y - m_middle.y) / m_dl.y;
		double iy1 = (pt1.y - m_middle.y) / m_dl.y;
		double iy2 = (pt2.y - m_middle.y) / m_dl.y;
		double iz0 = (pt0.z - m_middle.z) / m_dl.z;
		double iz1 = (pt1.z - m_middle.z) / m_dl.z;
		double iz2 = (pt2.z - m_middle.z) / m_dl.z;

		bool xlow  = 
			(ix0 < EPSILON || ix1 < EPSILON || ix2 < EPSILON) && 
			(ix0 > -2.0-EPSILON || ix1 > -2.0-EPSILON || ix2 > -2.0-EPSILON);
		bool xhigh = 
			(ix0 < 2.0+EPSILON || ix1 < 2.0+EPSILON || ix2 < 2.0+EPSILON) && 
			(ix0 > -EPSILON || ix1 > -EPSILON || ix2 > -EPSILON);
		bool ylow  = 
			(iy0 < EPSILON || iy1 <  EPSILON || iy2 <  EPSILON) && 
			(iy0 > -2.0-EPSILON || iy1 >  -2.0-EPSILON || iy2 >  -2.0-EPSILON);
		bool yhigh = 
			(iy0 < 2.0+EPSILON || iy1 < 2.0+EPSILON || iy2 < 2.0+EPSILON) && 
			(iy0 > -EPSILON || iy1 > -EPSILON || iy2 > -EPSILON);
		bool zlow  = 
			(iz0 < EPSILON || iz1 <  EPSILON || iz2 <  EPSILON) && 
			(iz0 > -2.0-EPSILON || iz1 >  -2.0-EPSILON || iz2 >  -2.0-EPSILON);
		bool zhigh = 
			(iz0 < 2.0+EPSILON || iz1 < 2.0+EPSILON || iz2 < 2.0+EPSILON) && 
			(iz0 > -EPSILON || iz1 > -EPSILON || iz2 > -EPSILON);

		int fres[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
		if(xlow  && ylow  && zlow)  fres[VX3D_LSW] = m_leaves[VX3D_LSW]->markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
		if(xlow  && yhigh && zlow)  fres[VX3D_LNW] = m_leaves[VX3D_LNW]->markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
		if(xhigh && ylow  && zlow)  fres[VX3D_LSE] = m_leaves[VX3D_LSE]->markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
		if(xhigh && yhigh && zlow)  fres[VX3D_LNE] = m_leaves[VX3D_LNE]->markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
		if(xlow  && ylow  && zhigh) fres[VX3D_HSW] = m_leaves[VX3D_HSW]->markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
		if(xlow  && yhigh && zhigh) fres[VX3D_HNW] = m_leaves[VX3D_HNW]->markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
		if(xhigh && ylow  && zhigh) fres[VX3D_HSE] = m_leaves[VX3D_HSE]->markInsideNodes(pt0, pt1, pt2, face, MIN_DL);
		if(xhigh && yhigh && zhigh) fres[VX3D_HNE] = m_leaves[VX3D_HNE]->markInsideNodes(pt0, pt1, pt2, face, MIN_DL);

		//if(true){
		//	MeshViewSet* set = new MeshViewSet();
		//	set->addFace(face);
		//	drawToViewSet(set, 1);
		//	if(xlow  && ylow  && zlow)  m_leaves[VX3D_LSW]->drawToViewSet(set, 0, 0.95);
		//	if(xlow  && yhigh && zlow)  m_leaves[VX3D_LNW]->drawToViewSet(set, 0, 0.95);
		//	if(xhigh && ylow  && zlow)  m_leaves[VX3D_LSE]->drawToViewSet(set, 0, 0.95);
		//	if(xhigh && yhigh && zlow)  m_leaves[VX3D_LNE]->drawToViewSet(set, 0, 0.95);
		//	if(xlow  && ylow  && zhigh) m_leaves[VX3D_HSW]->drawToViewSet(set, 0, 0.95);
		//	if(xlow  && yhigh && zhigh) m_leaves[VX3D_HNW]->drawToViewSet(set, 0, 0.95);
		//	if(xhigh && ylow  && zhigh) m_leaves[VX3D_HSE]->drawToViewSet(set, 0, 0.95);
		//	if(xhigh && yhigh && zhigh) m_leaves[VX3D_HNE]->drawToViewSet(set, 0, 0.95);

		//	SHOW_MESH("sub-leaves for general face crossing", set);
		//}

		// propagate-up "forbidden-edges" info
		for(int i = 0; i < 12; i++){
			if(((fres[FEDGES3D[i][0]] & FEDGES3D[i][2]) != 0) || ((fres[FEDGES3D[i][1]] & FEDGES3D[i][2]) != 0)) {
				m_vertices[FEDGES3D[i][0]]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FEDGES3D[i][3]);
				m_vertices[FEDGES3D[i][1]]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FEDGES3D[i][4]);
				res |= FEDGES3D[i][2];
			}
		}
	}else{
		MeshBlock* bl0 = face->getBlock(0);
		MeshBlock* bl1 = face->getBlock(1);
		int inside0 = -1;
		if(bl0) inside0 = (bl0->getAreaID() > -1) ? 1 : 0;
		int inside1 = -1;
		if(bl1) inside1 = (bl1->getAreaID() > -1) ? 1 : 0;

		int vtags[8], oryg_vtags[8];
		for(int i = 0; i < 8; i++)
			oryg_vtags[i] = vtags[i] = m_vertices[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -2);

		//MeshViewSet* set = new MeshViewSet();
		//set->addFace(face, 2);
		//const DVector3d nv = face->getNormalVector() * pt0.distance(pt1);
		//const DPoint3d mpt = face->getMiddlePoint();
		//if(bl0){
		//	set->addEdge(mpt, mpt-nv, inside0);
		//	set->addPoint(mpt-nv, inside0, inside0);
		//}
		//if(bl1){
		//	set->addEdge(mpt, mpt+nv, inside1);
		//	set->addPoint(mpt+nv, inside1, inside1);
		//}
		//int cross_edges = 0;

		// for each edge of the oct-block 
		//   check face-segment crossing
		//	... cross     -> mark nodes + mark as forbidden
		//  ... not-cross -> nothing
		//  ... undecided -> mark as forbidden
		const double MIN_DL2 = MIN_DL * MIN_DL;
		for(int i = 0; i < 12; i++){
			int i0 = FEDGES3D[i][0];
			int i1 = FEDGES3D[i][1];
			const DPoint3d& ept0 = m_vertices[i0]->coord;
			const DPoint3d& ept1 = m_vertices[i1]->coord;
			const DTriangle3d face_triangle(pt0, pt1, pt2);

			// control vertices too close to the face plane?
			if(vtags[i0] < 2){
				if(face_triangle.distance2ToPoint(ept0) < MIN_DL2){
					m_vertices[i0]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
					m_vertices[i0]->setIntTag(TagExtended::TAG_INSIDE_AREA, vtags[i0] = 2);
				}
			}
			if(vtags[i1] < 2){
				if(face_triangle.distance2ToPoint(ept1) < MIN_DL2){
					m_vertices[i1]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
					m_vertices[i1]->setIntTag(TagExtended::TAG_INSIDE_AREA, vtags[i1] = 2);
				}
			}
			if(vtags[i0] > 1  || vtags[i1] > 1) continue;

			// control segment too close to the face edges?
			if( DSegment3d::distance2ToSegment(ept0, ept1, pt0, pt1) < MIN_DL2 ||
				DSegment3d::distance2ToSegment(ept0, ept1, pt1, pt2) < MIN_DL2 ||
				DSegment3d::distance2ToSegment(ept0, ept1, pt2, pt0) < MIN_DL2)
			{
				// if close to edges -> too ambiguous, mark edge as forbidden and skip
				m_vertices[i0]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FEDGES3D[i][3]);
				m_vertices[i1]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FEDGES3D[i][4]);
				res |= FEDGES3D[i][2];
				continue;
			}

			// else -> crossing or not ?
			if(face_triangle.crossSegment(ept0, ept1)){
				// view
				//set->addEdge(ept0, ept1);
				//cross_edges++;
				// log
				//LOG4CPLUS_INFO(MeshLog::logger_console, "=======> MIN_DL2", MIN_DL2);
				//LOG4CPLUS_INFO(MeshLog::logger_console, "======== 0-1", DSegment3d::distance2ToSegment(ept0, ept1, pt0, pt1));
				//LOG4CPLUS_INFO(MeshLog::logger_console, "======== 1-2", DSegment3d::distance2ToSegment(ept0, ept1, pt1, pt2));
				//LOG4CPLUS_INFO(MeshLog::logger_console, "======== 2-0", DSegment3d::distance2ToSegment(ept0, ept1, pt2, pt0));
				// which side ?
				double v0 = face_triangle.orient3d(ept0);
				int ind0 = (v0 > 0.0) ? 0 : 1;
				double fv0 = abs(v0);
				double fv1 = abs(face_triangle.orient3d(ept1));
				// distance
				double dist = ept0.distance(ept1);
				double edist[2] = { dist * fv0 / (fv0+fv1), 0.0 };
				edist[1] = dist - edist[0];
				// mark
				assert(edist[0] >= MIN_DL);
				assert(edist[1] >= MIN_DL);
				// ... first side
				if(m_vertices[FEDGES3D[i][ind0]]->getDoubleTag(
					TagExtended::TAG_INSIDE_AREA_DIST, 2*edist[ind0]) > edist[ind0]) // if not available will also be greater
				{
					m_vertices[FEDGES3D[i][ind0]]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, edist[ind0]);
					m_vertices[FEDGES3D[i][ind0]]->setIntTag(TagExtended::TAG_INSIDE_AREA, inside0);
					// if(!bl0) inside0 == -1 which means it is unknown and will be removed later
					vtags[FEDGES3D[i][ind0]] = inside0;
				}
				// ... second side
				ind0 = 1 - ind0; // switch 0<->1
				if(m_vertices[FEDGES3D[i][ind0]]->getDoubleTag(
					TagExtended::TAG_INSIDE_AREA_DIST, 2*edist[ind0]) > edist[ind0])
				{
					m_vertices[FEDGES3D[i][ind0]]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, edist[ind0]);
					m_vertices[FEDGES3D[i][ind0]]->setIntTag(TagExtended::TAG_INSIDE_AREA, inside1);
					vtags[FEDGES3D[i][ind0]] = inside1;
				}

				// mark as forbidden
				m_vertices[i0]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FEDGES3D[i][3]);
				m_vertices[i1]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FEDGES3D[i][4]);
				res |= FEDGES3D[i][2];
			}
		}

		//if(cross_edges > 0){
		//	for(int i = 0; i < 8; i++){
		//		if(vtags[i] > -2){
		//			set->addPoint(m_vertices[i]->coord, (vtags[i] != oryg_vtags[i]) ? 1 : -2, vtags[i]);
		//		}else{
		//			set->addPoint(m_vertices[i]->coord);
		//		}
		//	}
		//	SHOW_MESH("leaf", set);
		//}else{
		//	delete set;
		//}
	}

	return res;
}

void ControlSpace3dOctree::OctLeaf::markInsideNodesXY(
	const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, 
	MeshFace* face, bool plane_z1, const double& MIN_DL)
{
	static const int v0[4] = {VX3D_LSE, VX3D_LSW, VX3D_LNE, VX3D_LNW};
	static const int v1[4] = {VX3D_HSE, VX3D_HSW, VX3D_HNE, VX3D_HNW};

	if(isSplit()){
		// recursive
		for(int i = 0; i < 4; i++)
			m_leaves[plane_z1 ? v1[i] : v0[i]]->markInsideNodesXY(pt0, pt1, pt2, face, plane_z1, MIN_DL);
	}else{
		const DPoint2d pt2d0(pt0.x, pt0.y);
		const DPoint2d pt2d1(pt1.x, pt1.y);
		const DPoint2d pt2d2(pt2.x, pt2.y);

		MeshBlock* b0 = face->getBlock(0);
		MeshBlock* b1 = face->getBlock(1);
		int inside0 = -1;
		if(b0) inside0 = (b0->getAreaID() > -1) ? 1 : 0;
		int inside1 = -1;
		if(b1) inside1 = (b1->getAreaID() > -1) ? 1 : 0;

		//MeshViewSet* set = new MeshViewSet();
		//set->addFace(face, 2);
		//const DVector3d nv = face->getNormalVector() * pt0.distance(pt1);
		//const DPoint3d mpt = face->getMiddlePoint();
		//if(b0){
		//	set->addEdge(mpt, mpt-nv, inside0);
		//	set->addPoint(mpt-nv, inside0, inside0);
		//}
		//if(b1){
		//	set->addEdge(mpt, mpt+nv, inside1);
		//	set->addPoint(mpt+nv, inside1, inside1);
		//}
		//int setview_points[8] = { -2, -2, -2, -2, -2, -2, -2, -2 };

		const double MIN_DL2 = MIN_DL * MIN_DL;

		for(int i = 0; i < 4; i++){
			auto node = m_vertices[plane_z1 ? v1[i] : v0[i]];
			const DPoint2d pt2d(node->coord.x, node->coord.y);
			if( DSegment2d::distanceToPoint2(pt2d0, pt2d1, pt2d) < MIN_DL2 ||
				DSegment2d::distanceToPoint2(pt2d1, pt2d2, pt2d) < MIN_DL2 ||
				DSegment2d::distanceToPoint2(pt2d2, pt2d0, pt2d) < MIN_DL2)
			{ // if node lies too near to one of the triangle edges
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
				node->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
			}else if(DTriangle2d::containsPoint(pt2d0, pt2d1, pt2d2, pt2d, MIN_DL2)){
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
				node->setIntTag(TagExtended::TAG_INSIDE_AREA, 4);
				//setview_points[plane_z1 ? v1[i] : v0[i]] = 4;
				// other node
				node = m_vertices[plane_z1 ? v0[i] : v1[i]];
				double dist = 4*m_dl.z;
				if(node->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 2*dist) > dist)
				{
					node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
					int v = (DTriangle3d::orient3d(pt0, pt1, pt2, node->coord) > 0.0) ? inside0 : inside1;
					if(v >= 0){
						node->setIntTag(TagExtended::TAG_INSIDE_AREA, v);
						//setview_points[plane_z1 ? v0[i] : v1[i]] = v;
					}else{
						node->removeTag(TagExtended::TAG_INSIDE_AREA);
						//setview_points[plane_z1 ? v0[i] : v1[i]] = -1;
					}
				}else{
					//setview_points[plane_z1 ? v0[i] : v1[i]] = 
					//	node->getIntTag(TagExtended::TAG_INSIDE_AREA, -2);
				}
			}
		}

		//for(int i = 0; i < 8; i++){
		//	if(setview_points[i] > -2)
		//		set->addPoint(m_vertices[i]->coord, setview_points[i], setview_points[i]);
		//	else
		//		set->addPoint(m_vertices[i]->coord);
		//}
		//SHOW_MESH("leaf XY", set);

	}
}

void ControlSpace3dOctree::OctLeaf::markInsideNodesXZ(
	const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, 
	MeshFace* face, bool plane_y1, const double& MIN_DL)
{
	static const int v0[4] = {VX3D_LSE, VX3D_LSW, VX3D_HSE, VX3D_HSW};
	static const int v1[4] = {VX3D_LNE, VX3D_LNW, VX3D_HNE, VX3D_HNW};

	if(isSplit()){
		// recursive
		for(int i = 0; i < 4; i++)
			m_leaves[plane_y1 ? v1[i] : v0[i]]->markInsideNodesXZ(pt0, pt1, pt2, face, plane_y1, MIN_DL);
	}else{
		const DPoint2d pt2d0(pt0.x, pt0.z);
		const DPoint2d pt2d1(pt1.x, pt1.z);
		const DPoint2d pt2d2(pt2.x, pt2.z);

		const double MIN_DL2 = MIN_DL * MIN_DL;

		MeshBlock* b0 = face->getBlock(0);
		MeshBlock* b1 = face->getBlock(1);
		int inside0 = -1;
		if(b0) inside0 = (b0->getAreaID() > -1) ? 1 : 0;
		int inside1 = -1;
		if(b1) inside1 = (b1->getAreaID() > -1) ? 1 : 0;

		for(int i = 0; i < 4; i++){
			auto node = m_vertices[plane_y1 ? v1[i] : v0[i]];
			const DPoint2d pt2d(node->coord.x, node->coord.z);
			if( DSegment2d::distanceToPoint2(pt2d0, pt2d1, pt2d) < MIN_DL2 ||
				DSegment2d::distanceToPoint2(pt2d1, pt2d2, pt2d) < MIN_DL2 ||
				DSegment2d::distanceToPoint2(pt2d2, pt2d0, pt2d) < MIN_DL2)
			{ // if node lies too near to one of the triangle edges
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
				node->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
			}else if(DTriangle2d::containsPoint(pt2d0, pt2d1, pt2d2, pt2d, MIN_DL2)){
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
				node->setIntTag(TagExtended::TAG_INSIDE_AREA, 4);
				// other node
				node = m_vertices[plane_y1 ? v0[i] : v1[i]];
				double dist = 4*m_dl.y;
				if(node->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 2*dist) > dist)
				{
					node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
					int v = (DTriangle3d::orient3d(pt0, pt1, pt2, node->coord) > 0.0) ? inside0 : inside1;
					if(v >= 0) node->setIntTag(TagExtended::TAG_INSIDE_AREA, v);
					else node->removeTag(TagExtended::TAG_INSIDE_AREA);
				}
			}
		}
	}
}

void ControlSpace3dOctree::OctLeaf::markInsideNodesYZ(
	const DPoint3d& pt0, const DPoint3d& pt1, const DPoint3d& pt2, 
	MeshFace* face, bool plane_x1, const double& MIN_DL)
{
	static const int v0[4] = {VX3D_LSW, VX3D_HSW, VX3D_LNW, VX3D_HNW};
	static const int v1[4] = {VX3D_LSE, VX3D_HSE, VX3D_LNE, VX3D_HNE};

	if(isSplit()){
		// recursive
		for(int i = 0; i < 4; i++)
			m_leaves[plane_x1 ? v1[i] : v0[i]]->markInsideNodesYZ(pt0, pt1, pt2, face, plane_x1, MIN_DL);
	}else{
		const DPoint2d pt2d0(pt0.y, pt0.z);
		const DPoint2d pt2d1(pt1.y, pt1.z);
		const DPoint2d pt2d2(pt2.y, pt2.z);

		const double MIN_DL2 = MIN_DL * MIN_DL;

		MeshBlock* b0 = face->getBlock(0);
		MeshBlock* b1 = face->getBlock(1);
		int inside0 = -1;
		if(b0) inside0 = (b0->getAreaID() > -1) ? 1 : 0;
		int inside1 = -1;
		if(b1) inside1 = (b1->getAreaID() > -1) ? 1 : 0;

		for(int i = 0; i < 4; i++){
			auto node = m_vertices[plane_x1 ? v1[i] : v0[i]];
			const DPoint2d pt2d(node->coord.y, node->coord.z);
			if( DSegment2d::distanceToPoint2(pt2d0, pt2d1, pt2d) < MIN_DL2 ||
				DSegment2d::distanceToPoint2(pt2d1, pt2d2, pt2d) < MIN_DL2 ||
				DSegment2d::distanceToPoint2(pt2d2, pt2d0, pt2d) < MIN_DL2)
			{ // if node lies too near to one of the triangle edges
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
				node->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
			}else if(DTriangle2d::containsPoint(pt2d0, pt2d1, pt2d2, pt2d, MIN_DL2)){
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
				node->setIntTag(TagExtended::TAG_INSIDE_AREA, 4);
				// other node
				node = m_vertices[plane_x1 ? v0[i] : v1[i]];
				double dist = 4*m_dl.x;
				if(node->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 2*dist) > dist)
				{
					node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
					int v = (DTriangle3d::orient3d(pt0, pt1, pt2, node->coord) > 0.0) ? inside0 : inside1;
					if(v >= 0) node->setIntTag(TagExtended::TAG_INSIDE_AREA, v);
					else node->removeTag(TagExtended::TAG_INSIDE_AREA);
				}
			}
		}
	}
}

void ControlSpace3dOctree::markInsideNodesAtBoundary()
{
	int nxy = m_n[0] * m_n[1];
	// XY
	int k0 = 0;
	int k1 = nxy * (m_n[2] - 1);
	for(int i = 0; i < nxy; i++){
		m_grid[k0++].markInsideNodesAtBoundary(FC3D_LOW);
		m_grid[k1++].markInsideNodesAtBoundary(FC3D_HIGH);
	}
	// XZ
	for(int i = 0; i < m_n[2]; i++){
		k0 = nxy * i;
		k1 = k0 + nxy - m_n[0];
		for(int j = 0; j < m_n[0]; j++){
			m_grid[k0++].markInsideNodesAtBoundary(FC3D_SOUTH);
			m_grid[k1++].markInsideNodesAtBoundary(FC3D_NORTH);
		}
	}
	// YZ
	for(int i = 0; i < m_n[2]; i++){
		k0 = nxy * i;
		k1 = k0 + m_n[0] - 1;
		for(int j = 0; j < m_n[1]; j++){
			m_grid[k0].markInsideNodesAtBoundary(FC3D_WEST);
			m_grid[k1].markInsideNodesAtBoundary(FC3D_EAST);
			k0 += m_n[0];
			k1 += m_n[0];
		}
	}
}

/// mark inside-information for boundary CS-leaves
void ControlSpace3dOctree::OctLeaf::markInsideNodesAtBoundary(OctFaceWhich side)
{
	if(isSplit()){
		for(int i = 0; i < 4; i++)
			m_leaves[DBox::face_to_vertex[side][i]]->markInsideNodesAtBoundary(side);
	}else{
		for(int i = 0; i < 4; i++){
			int k = DBox::face_to_vertex[side][i];
//			m_vertices[k]->setIntTag(TagExtended::TAG_ID, side);
			if(m_vertices[k]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1) < 0){
				m_vertices[k]->setIntTag(TagExtended::TAG_INSIDE_AREA, 0); // i.e. "not-inside"
			}
		}
	}
}

bool ControlSpace3dOctree::OctLeaf::propagateInsideMark()
{
	if(isSplit()){
		for(int i = 0; i < 8; i++)
			m_leaves[i]->propagateInsideMark();
	}

	int inside_mark[] = {
		m_vertices[0]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[1]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[2]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[3]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[4]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[5]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[6]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[7]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1) };

	int empty_count = 0;
	for(int i = 0; i < 8; i++)
		if(inside_mark[i] < 0) ++empty_count;
	if(empty_count == 0 || empty_count == 8) return (empty_count == 0);

	int forbidden[] = {
		m_vertices[0]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[1]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[2]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[3]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[4]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[5]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[6]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[7]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN) 
	};

//	static int counter = 10;

	// show degug if: at least one "0" (+ at least one forbidden ...)
	//if(inside_mark[0] == 0 || inside_mark[1] == 0 || inside_mark[2] == 0 || inside_mark[3] == 0){
	//	MeshViewSet* set = new MeshViewSet;
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX3D_SW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX3D_SE]->coord, 0.0), 
	//		(m_vertices[VX3D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) ||
	//		 m_vertices[VX3D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX3D_NW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX3D_NE]->coord, 0.0),  
	//		(m_vertices[VX3D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) ||
	//		 m_vertices[VX3D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX3D_SW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX3D_NW]->coord, 0.0),  
	//		(m_vertices[VX3D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) ||
	//		 m_vertices[VX3D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX3D_SE]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX3D_NE]->coord, 0.0),  
	//		(m_vertices[VX3D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) ||
	//		 m_vertices[VX3D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL)) ? 1 : 0);
	//	for(int i = 0; i < 4; i++){
	//		int v = m_vertices[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1);
	//		if(v >= 0) set->addPoint(DPoint3d(m_vertices[i]->coord, 0.0), v, v);
	//	}

	//	++counter;
	//	ostringstream oss;
	//	oss << "Before, counter = " << counter;
	//	SHOW_MESH(oss.str().c_str(), set);
	//}

	/// TODO -> 4 edges to 12 edges

	while(empty_count > 0){
		int changes = 0;
		for(int i = 0; (i < 12) && (empty_count > 0); i++){
			// SW - SE
			if( ((inside_mark[FEDATA[i][0]] < 0) != (inside_mark[FEDATA[i][1]] < 0)) &&
		 		 (inside_mark[FEDATA[i][0]] + inside_mark[FEDATA[i][1]] < 1) &&
		 		 (((forbidden[FEDATA[i][0]] & FEDATA[i][2]) == 0) ||
				  ((forbidden[FEDATA[i][1]] & FEDATA[i][3]) == 0)))
			{
				if(inside_mark[FEDATA[i][0]] == -1){
					m_vertices[FEDATA[i][0]]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
						inside_mark[FEDATA[i][0]] = inside_mark[FEDATA[i][1]]);
	//				m_vertices[FEDATA[i][0]]->setIntTag(TagExtended::TAG_ID, ++counter);
				}else{
					m_vertices[FEDATA[i][1]]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
						inside_mark[FEDATA[i][1]] = inside_mark[FEDATA[i][0]]);
	//				m_vertices[FEDATA[i][1]]->setIntTag(TagExtended::TAG_ID, ++counter);
				}
				--empty_count;
				++changes;
			}
		}
		// no applicable modification ? skip the loop
		if(changes == 0) break;
	}

	// show degug if: at least one "0" (+ at least one forbidden ...)
	//if(inside_mark[0] == 0 || inside_mark[1] == 0 || inside_mark[2] == 0 || inside_mark[3] == 0){
	//	MeshViewSet* set = new MeshViewSet;
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX3D_SW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX3D_SE]->coord, 0.0), 
	//		(m_vertices[VX3D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) &&
	//		 m_vertices[VX3D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX3D_NW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX3D_NE]->coord, 0.0),  
	//		(m_vertices[VX3D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) &&
	//		 m_vertices[VX3D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX3D_SW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX3D_NW]->coord, 0.0),  
	//		(m_vertices[VX3D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) &&
	//		 m_vertices[VX3D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX3D_SE]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX3D_NE]->coord, 0.0),  
	//		(m_vertices[VX3D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) &&
	//		 m_vertices[VX3D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL)) ? 1 : 0);
	//	for(int i = 0; i < 4; i++){
	//		int v = m_vertices[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1);
	//		if(v >= 0) set->addPoint(DPoint3d(m_vertices[i]->coord, 0.0), v, v);
	//	}

	//	ostringstream oss;
	//	oss << "After, counter = " << counter;
	//	SHOW_MESH(oss.str().c_str(), set);
	//}

	return (empty_count == 0);
}

bool ControlSpace3dOctree::OctLeaf::setInsideMark(int value)
{
	if(isSplit()){
		for(int i = 0; i < 8; i++)
			m_leaves[i]->setInsideMark(value);
	}

	for(int i = 0; i < 8; i++)
		if(m_vertices[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1) < 0)
			m_vertices[i]->setIntTag(TagExtended::TAG_INSIDE_AREA, value);

	return true;
}

/// Propagate inside mark for the main grid
void ControlSpace3dOctree::propagateMainGridInsideMark()
{
	// Start from the initialized nodes
	int vx = m_n[0]+1;
	int vy = m_n[1]+1;
	int vz = m_n[2]+1;
	int vxy = vx * vy;
	int main_count = vx*vy*vz;
	DataVector<int> ready_nodes(main_count);
	for(int i = 0; i < main_count; i++){
		if(m_grid_vertices[i]->availableTag(TagExtended::TAG_INSIDE_AREA))
			ready_nodes.add(i);
	}

	int start_ready_nodes = (int)ready_nodes.countInt();
	if(start_ready_nodes == 0) return;
	for(int i = 0; i < ready_nodes.countInt() && ready_nodes.countInt() < main_count; i++){
		int k = ready_nodes[i];
		auto qvk = m_grid_vertices[k];
		int inside_tag = qvk->getIntTag(TagExtended::TAG_INSIDE_AREA);
		if(inside_tag > 1) continue; // only 0 and 1 are used to extend info
		int iz = k / vxy;
		int kxy = k % vxy;
		int ix = kxy % vx;
		int iy = kxy / vx;
		// check neighbours of 'k'
		int nid[]   = { (k-1),    (k+1),       (k-vx),   (k+vx),      (k-vxy),  (k+vxy)     };
		bool nbrd[] = { (ix > 0), (ix < vx-1), (iy > 0), (iy < vy-1), (iz > 0), (iz < vz-1) };
		int nflags[]  = { 
			FORBIDDEN_EDGE_XL, FORBIDDEN_EDGE_XH, 
			FORBIDDEN_EDGE_YL, FORBIDDEN_EDGE_YH,
			FORBIDDEN_EDGE_ZL, FORBIDDEN_EDGE_ZH 
		};
		int nrflags[] = { 
			FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL, 
			FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL, 
			FORBIDDEN_EDGE_ZH, FORBIDDEN_EDGE_ZL 
		};

		for(int j = 0; j < 6; j++){
			if(nbrd[j]){
				auto qv = m_grid_vertices[nid[j]];
				if( !qv->availableTag(TagExtended::TAG_INSIDE_AREA) &&
					(!qv->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, nflags[j]) ||
					 !qvk->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, nrflags[j])))
				{
					qv->setIntTag(TagExtended::TAG_INSIDE_AREA, inside_tag);
					ready_nodes.add(nid[j]);
				}
			}
		}
	}
}

/// Show CS-grid with inside info for debug
void ControlSpace3dOctree::showMarkInside(const DataVector<MeshFace*>& bfaces) const
{
	// start with CS-grid edges
	// MeshViewSet* set = getViewSet();
	MeshViewSet* set = new MeshViewSet();

	// border edges
	for(int i = 0; i < bfaces.countInt(); i++){
		MeshFace* face = bfaces[i];
		//if(face->availableTag(TagExtended::TAG_INSIDE_AREA))
		//	set->addFace(face);
		//else
			for(int j = 0; j < 3; j++)
				set->addEdge(face->getEdge(j), 1);
	}

	// add inside-marks for vertices 
	int gv_count = m_grid_vertices.countInt();
	for(int i = 0; i < gv_count; i++){
		auto gv = m_grid_vertices[i];
		if(gv->availableTag(TagExtended::TAG_INSIDE_AREA)){
			int v = gv->getIntTag(TagExtended::TAG_INSIDE_AREA);
			set->addPoint(gv->coord, v, gv->getIntTag(TagExtended::TAG_ID, v));
		}else 
			set->addPoint(gv->coord);
	}


	// + forbidden edges
	DataSimpleList<ControlSpace3dOctree::OctLeaf*> leaves;
	int nxyz = m_n[0] * m_n[1] * m_n[2];
	for(int i = 0; i < nxyz; i++)
		leaves.append(m_grid+i);

	while(!leaves.empty()){
		OctLeaf* leaf = leaves.removeFirst();

		if(leaf->isSplit()){
			for(int i = 0; i < 8; i++)
				leaves.append(leaf->m_leaves[i]);
		}

		for(int i = 0; i < 12; i++){
			if( leaf->m_vertices[OctLeaf::FEDATA[i][0]]->hasIntFlag(
					TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, OctLeaf::FEDATA[i][2]) &&
				leaf->m_vertices[OctLeaf::FEDATA[i][1]]->hasIntFlag(
					TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, OctLeaf::FEDATA[i][3]))
			{
				DPoint3d pt0(leaf->m_middle, leaf->m_vertices[OctLeaf::FEDATA[i][0]]->coord, 0.95);
				DPoint3d pt1(leaf->m_middle, leaf->m_vertices[OctLeaf::FEDATA[i][1]]->coord, 0.95);

				set->addEdge(pt0, pt1, 3);
			}
		}
	}

	// show
	SHOW_MESH("CS-grid-3d with inside marks", set);
}

/// draw octree element
void ControlSpace3dOctree::OctLeaf::drawToViewSet(MeshViewSet* set, int part, double ratio) const
{
	DPoint3d pts[8];
	for(int i = 0; i < 8; i++){
		pts[i] = DPoint3d(m_middle, m_vertices[i]->coord, ratio);
	}

	set->addEdge(pts[VX3D_LSW], pts[VX3D_LSE], part);
	set->addEdge(pts[VX3D_LNW], pts[VX3D_LNE], part);
	set->addEdge(pts[VX3D_LSW], pts[VX3D_LNW], part);
	set->addEdge(pts[VX3D_LSE], pts[VX3D_LNE], part);

	set->addEdge(pts[VX3D_HSW], pts[VX3D_HSE], part);
	set->addEdge(pts[VX3D_HNW], pts[VX3D_HNE], part);
	set->addEdge(pts[VX3D_HSW], pts[VX3D_HNW], part);
	set->addEdge(pts[VX3D_HSE], pts[VX3D_HNE], part);

	set->addEdge(pts[VX3D_LSW], pts[VX3D_HSW], part);
	set->addEdge(pts[VX3D_LNW], pts[VX3D_HNW], part);
	set->addEdge(pts[VX3D_LSE], pts[VX3D_HSE], part);
	set->addEdge(pts[VX3D_LNE], pts[VX3D_HNE], part);
}

double ControlSpace3dOctree::getLocalResolution(const DPoint3d& pt) const
{
	OctLeaf* leaf = findLastLeaf(pt);
	return std::min(std::min(leaf->m_dl.x, leaf->m_dl.y), leaf->m_dl.z);
}

double ControlSpace3dOctree::interpolateDoubleTag(const DPoint3d& pt, TagExtended::TagType type) const
{
	const DPoint3d pt_fit = m_box.fitInPoint(pt);
	return findLastLeaf(pt_fit)->interpolateDoubleTag(pt_fit, type);
}

void ControlSpace3dOctree::setMinDoubleTag(const DPoint3d& pt, TagExtended::TagType type, double q)
{
	const DPoint3d pt_fit = m_box.fitInPoint(pt);
	OctLeaf* leaf = findLastLeaf(pt_fit);
	if(!leaf) return;
	for(int i = 0; i < 8; i++){
		auto node = leaf->getVertex(i);
		double v = node->getDoubleTag(type);
		if(q > v) node->setDoubleTag(type, q);
	}
}

void ControlSpace3dOctree::setMaxDoubleTag(const DPoint3d& pt, TagExtended::TagType type, double q)
{
	const DPoint3d pt_fit = m_box.fitInPoint(pt);
	OctLeaf* leaf = findLastLeaf(pt_fit);
	if(!leaf) return;
	for(int i = 0; i < 8; i++){
		auto node = leaf->getVertex(i);
		double v = node->getDoubleTag(type);
		if(q < v) node->setDoubleTag(type, q);
	}
}

// stat	

int ControlSpace3dOctree::getElementCount() const 
{ 
	int grid_count = m_n[0] * m_n[1] * m_n[2];
	int element_count = 0;
	for (int i = 0; i < grid_count; i++)
		element_count += m_grid[i].getElementCount();
	return element_count;
}

int ControlSpace3dOctree::getMaxDepth() const 
{ 
	int grid_count = m_n[0] * m_n[1] * m_n[2];
	int max_depth = 0;
	for (int i = 0; i < grid_count; i++) {
		int depth = m_grid[i].getMaxDepth();
		if (depth > max_depth) max_depth = depth;
	}
	return max_depth;
}

int ControlSpace3dOctree::getMaxBalance() const 
{ 
	int max_balance = 0;
	int grid_count = m_n[0] * m_n[1] * m_n[2];
	for (int i = 0; i < grid_count; i++) {
		m_grid[i].getMaxDepthAndBalance(max_balance);
	}
	return max_balance;
}

int ControlSpace3dOctree::getTotalBytes() const 
{ 
	return sizeof(*this) // main struct
		+ sizeof(OctLeaf) * getElementCount() // tree nodes
		+ sizeof(ControlNode3d) * getValueCount(); // values
}

/// Returns the "screenshot" of this control space for visualization
MeshViewSet* ControlSpace3dOctree::getViewSet(MeshViewSet* set, bool with_points) const
{
	if(set){
		set->prepareFreePlace(m_grid_vertices.countInt(), 0, 4*m_n[0]*m_n[1]*m_n[2], 0);
	}else{
		set = new MeshViewSet(m_grid_vertices.countInt(), 0, 4*m_n[0]*m_n[1]*m_n[2], 0);
	}

	double dx = m_box.getDX() / m_n[0];
	double dy = m_box.getDY() / m_n[1];
	double dz = m_box.getDZ() / m_n[2];

	// main grid
	double z = m_box.z0;
	for(int i = 0; i <= m_n[2]; i++){
		DPoint3d ptx0(m_box.x0, m_box.y0, z);
		DPoint3d ptx1(m_box.x0, m_box.y1, z);
		for(int j = 0; j <= m_n[0]; j++){
			set->addEdge(ptx0, ptx1);
			ptx0.x += dx; ptx1.x += dx;
		}
		DPoint3d pty0(m_box.x0, m_box.y0, z);
		DPoint3d pty1(m_box.x1, m_box.y0, z);
		for(int j = 0; j <= m_n[1]; j++){
			set->addEdge(pty0, pty1);
			pty0.y += dy; pty1.y += dy;
		}
		z += dz;
	}
	DPoint3d ptz0(m_box.x0, m_box.y0, m_box.z0);
	DPoint3d ptz1(m_box.x0, m_box.y0, m_box.z1);
	for(int i = 0; i <= m_n[0]; i++){
		ptz0.y = m_box.y0;
		ptz1.y = m_box.y0;
		for(int j = 0; j <= m_n[1]; j++){
			set->addEdge(ptz0, ptz1);
			ptz0.y += dy; ptz1.y += dy;
		}
		ptz0.x += dx; ptz1.x += dx;
	}

	// inside
	int count = m_n[0]*m_n[1]*m_n[2];
	for(int k = 0; k < count; k++){
		m_grid[k].addToViewSet(set);
	}

	if(with_points){
		int gv_count = m_grid_vertices.countInt();
		for(int i = 0; i < gv_count; i++){
			auto qv = m_grid_vertices[i];
			if(qv->w < 0) 
				set->addPoint(qv->coord, (int)(-qv->w));
		}
	}

	return set;
}

void ControlSpace3dOctree::OctLeaf::addToViewSet(MeshViewSet* set) const
{
	if(isSplit()){
		DVector3d v[3] = {
			DVector3d(2*m_dl.x, 0.0, 0.0),
			DVector3d(0.0, 2*m_dl.y, 0.0),
			DVector3d(0.0, 0.0, 2*m_dl.z) };

		for(int i = 0; i < 3; i++){
			set->addEdge( m_middle-v[i], m_middle+v[i] );
			for(int j = 0; j < 3; j++){
				if(i == j) continue;
				set->addEdge( m_middle-v[i]-v[j], m_middle+v[i]-v[j] );
				set->addEdge( m_middle-v[i]+v[j], m_middle+v[i]+v[j] );
			}
		}
		for(int i = 0; i < 8; i++)
			m_leaves[i]->addToViewSet(set);
	}else{
		//MeshViewFaceData* data = new MeshViewFaceData(4);
		//data->area_id = 0;
		//data->quality = 1.0;

		//DPoint3d dpts[4];
		//int order[4] = {0, 1, 3, 2};
		//for(int i = 0; i < 4; i++){
		//	if(surface){
		//		dpts[i] = surface->getPoint(m_vertices[order[i]]->coord);
		//	}else{
		//		dpts[i] = DPoint3d(m_vertices[order[i]]->coord, 0.0);
		//	}
		//	if(tag == TagExtended::TAG_NONE){
		//		data->wpts[i] = 1.0 - 0.1 * m_vertices[order[i]]->max_gradation_ratio;
		//	}else{
		//		if(m_vertices[order[i]]->availableTag(tag)) 
		//			data->wpts[i] = m_vertices[order[i]]->getDoubleTag(tag);
		//		else 
		//			data->wpts[i] = 1.0;
		//	}
		//}

		//DPoint3d middle(dpts[0], dpts[2], 0.5);

		//for(int i = 0; i < 4; i++)
		//	data->pts[i] = middle + (dpts[i] - middle) * 0.95;

		//data->normal = (data->pts[1] - data->pts[0]).crossProduct(data->pts[2] - data->pts[0]).normalized();

		//set->m_faces.add(data);
	}

	if(m_extra && m_extra->m_source_points){
		int ct = (int)m_extra->m_source_points->countInt();
		for(int i = 0; i < ct; i++){
			auto node = m_extra->m_source_points->get(i);
			set->addPoint(node->coord, -1);
		}
	}
}

/// Insert additional info -> about local surface available for the given point
bool ControlSpace3dOctree::addLocalSurfaceAtPoint(
		const DPoint3d& pt, SurfaceConstPtr local_surface, 
		DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & local_surface_hash)
{
	OctLeaf* leaf = findLastLeaf(pt);
	assert(leaf != nullptr);
	leaf->addToLocalSurfaceSet( local_surface, local_surface_hash );
	return true;
}

/// Insert the given local surface into surface set (false -> already there)
bool ControlSpace3dOctree::OctLeaf::addToLocalSurfaceSet(
		SurfaceConstPtr local_surface, 
		DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & local_surface_hash)
{
	auto local_set = getLocalSurfaceSet();
	if(local_set && local_set->containsSurface( local_surface )) return false; // already set

	// nullptr or valid
	auto ext_local_set = local_surface_hash.getValue( local_set, nullptr );
	if( ! ext_local_set ){
		ext_local_set = std::make_shared<SurfaceParametricSet>( local_set, local_surface );
		local_surface_hash.setValue( local_set, ext_local_set );
	}

	setLocalSurfaceSet( ext_local_set );
	return true;
}

int ControlSpace3dOctree::OctLeaf::getElementCount() const {
	if (!isSplit()) return 1;
	int result = 1;
	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
		result += m_leaves[vi]->getElementCount();
	}
	return result;
}

int ControlSpace3dOctree::OctLeaf::getMaxDepth() const {
	if (!isSplit()) return 1;
	int dmax = 0;
	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
		int d = m_leaves[vi]->getMaxDepth();
		if (d > dmax) dmax = d;
	}
	return 1 + dmax;
}

int ControlSpace3dOctree::OctLeaf::getMaxDepthAndBalance(int & max_balance) const {
	if (!isSplit()) return 1;
	int dmax, dmin;
	dmax = dmin = m_leaves[VX3D_LAST]->getMaxDepthAndBalance(max_balance);
	for (int vi = VX3D_FIRST; vi < VX3D_LAST; vi++) {
		int d = m_leaves[vi]->getMaxDepthAndBalance(max_balance);
		if (d > dmax) dmax = d;
		if (d < dmin) dmin = d;
	}
	int balance = dmax - dmin;
	if (balance > max_balance) max_balance = balance;
	return 1 + dmax;
}

/// Insert additional info -> about local surface available for the given bounding box
bool ControlSpace3dOctree::addLocalSurfaceAtBBox(const DBox& box, 
	SurfaceConstPtr local_surface, 
	DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & local_surface_hash)
{
	double dx = m_box.getDX() / m_n[0];
	double dy = m_box.getDY() / m_n[1];
	double dz = m_box.getDZ() / m_n[2];
	int nxy = m_n[0] * m_n[1];

	double xmin = (box.x0 - m_box.x0) / dx; // x
	double xmax = (box.x1 - m_box.x0) / dx;
	double ymin = (box.y0 - m_box.y0) / dy; // y
	double ymax = (box.y1 - m_box.y0) / dy;
	double zmin = (box.z0 - m_box.z0) / dz; // z
	double zmax = (box.z1 - m_box.z0) / dz;

	int ixmin = std::max((int)xmin, 0);
	int ixmax = std::min((int)xmax, m_n[0]-1);
	int iymin = std::max((int)ymin, 0);
	int iymax = std::min((int)ymax, m_n[1]-1);
	int izmin = std::max((int)zmin, 0);
	int izmax = std::min((int)zmax, m_n[2]-1);

	// for leaves and sub-leaves -> recursively check collision box-leaf
	for(int iz = izmin; iz <= izmax; iz++){
		for(int iy = iymin; iy <= iymax; iy++){
			int k = iz * nxy + iy * m_n[0];
			for(int ix = ixmin; ix <= ixmax; ix++){
				m_grid[k + ix].addLocalSurfaceAtBBox(box, local_surface, local_surface_hash);
			}
		}
	}

	return true;
}

/// Insert additional info -> about local surface available for the given bounding box
void ControlSpace3dOctree::OctLeaf::addLocalSurfaceAtBBox(const DBox& box, 
		SurfaceConstPtr local_surface, 
		DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & local_surface_hash)
{
	if( isSplit() ){
		// check where to recurse
		bool is_x0 = (box.x0 <= m_middle.x);
		bool is_x1 = (box.x1 >= m_middle.x);
		bool is_y0 = (box.y0 <= m_middle.y);
		bool is_y1 = (box.y1 >= m_middle.y);
		bool is_z0 = (box.z0 <= m_middle.z);
		bool is_z1 = (box.z1 >= m_middle.z);

		if( is_x0 && is_y0 && is_z0) m_leaves[VX3D_LSW]->addLocalSurfaceAtBBox(box, local_surface, local_surface_hash);
		if( is_x1 && is_y0 && is_z0) m_leaves[VX3D_LSE]->addLocalSurfaceAtBBox(box, local_surface, local_surface_hash);
		if( is_x0 && is_y1 && is_z0) m_leaves[VX3D_LNW]->addLocalSurfaceAtBBox(box, local_surface, local_surface_hash);
		if( is_x1 && is_y1 && is_z0) m_leaves[VX3D_LNE]->addLocalSurfaceAtBBox(box, local_surface, local_surface_hash);
		if( is_x0 && is_y0 && is_z1) m_leaves[VX3D_HSW]->addLocalSurfaceAtBBox(box, local_surface, local_surface_hash);
		if( is_x1 && is_y0 && is_z1) m_leaves[VX3D_HSE]->addLocalSurfaceAtBBox(box, local_surface, local_surface_hash);
		if( is_x0 && is_y1 && is_z1) m_leaves[VX3D_HNW]->addLocalSurfaceAtBBox(box, local_surface, local_surface_hash);
		if( is_x1 && is_y1 && is_z1) m_leaves[VX3D_HNE]->addLocalSurfaceAtBBox(box, local_surface, local_surface_hash);
	}else{
		// check if within? should be anyhow..
		// store
		addToLocalSurfaceSet( local_surface, local_surface_hash );
	}
}

/// Returns additional info -> about local surfaces (set) available for the given point
SurfaceSetConstPtr ControlSpace3dOctree::getLocalSurfaceSetAtPoint(const DPoint3d& pt) 
{
	return findLastLeaf(pt)->getLocalSurfaceSet();
}
