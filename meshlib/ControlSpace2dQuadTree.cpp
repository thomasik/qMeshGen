///////////////////////////////////////////////////////////////////////////
// ControlSpace2dQuadTree.cpp
// This class implements the control space (sizing and stretching info)
//  in the form of the QuadTree
/////////////////////////////////////////////////////////////////////////////
//  Tomasz Jurczyk, 2005-
//  Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "ControlSpace2dQuadTree.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "SurfaceParametric.h"
#include "MeshContainer2d.h"
#include "Curve2dSegment.h"
#include "EPSFile.h"
#include "MeshViewSet.h"
#include "ReparameterizationData.h"
#include "SurfacePlane.h"
#include "DTriangle.h"
#include "DSegment.h"
#include "ControlSpace3d.h"
#include "ControlSpace3dAdaptive.h"
#include "MeshGenerator2d.h"

#ifdef _DEBUG
//#define STORE_CONTROL_EPS
#endif

int ControlSpace2dQuadTree::param_max_depth = 20;
int ControlSpace2dQuadTree::param_max_depth_surface_adaptation = 10;

/// Whether control qtree mid-nodes should be used for metric interpolation
int ControlSpace2dQuadTree::QuadLeaf::param_use_midnodes = 0;
/// balance level
int ControlSpace2dQuadTree::QuadLeaf::param_balance_level = 1;

ControlSpace2dQuadTree::QuadLeaf::QuadVertexWhich 
	ControlSpace2dQuadTree::QuadLeaf::mid_to_nodes[4][2] = 
			{{VX2D_SW, VX2D_SE}, {VX2D_SE, VX2D_NE},
			{VX2D_NW, VX2D_NE}, {VX2D_SW, VX2D_NW}};
ControlSpace2dQuadTree::QuadLeaf::QuadMidVertexWhich 
	ControlSpace2dQuadTree::QuadLeaf::nodes_to_mid[4][2] = 
			{{MVX2D_DOWN, MVX2D_LEFT}, {MVX2D_DOWN, MVX2D_RIGHT},
			{MVX2D_UP, MVX2D_LEFT}, {MVX2D_UP, MVX2D_RIGHT}};

ControlSpace2dQuadTree::ControlSpace2dQuadTree(SurfaceConstPtr surface, 
			const DRect& box, int nxy) 
	: ControlSpace2dAdaptive(surface, box), m_grid_vertices(2*nxy)
{
	LOG4CPLUS_INFO(MeshLog::logger_console, "Creating quadtree control space.");
	assert(nxy >= 4);
	// length x
	Curve2dSegment middle_xline(DPoint2d(box.x0, (box.y0+box.y1)*0.5), 
		DPoint2d(box.x1, (box.y0+box.y1)*0.5));
	double real_dx = middle_xline.getLengthOnSurface(0.0, 1.0, surface);
	// length y
	Curve2dSegment middle_yline(DPoint2d((box.x0+box.x1)*0.5, box.y0), 
		DPoint2d((box.x0+box.x1)*0.5, box.y1));
	double real_dy = middle_yline.getLengthOnSurface(0.0, 1.0, surface);

	double ratio = real_dx / real_dy;	// adjust m_nx:m_ny ratio to real lengths of domain
	m_nx = std::max(2,(int)(sqrt(ratio*nxy)));
	m_ny = std::max(2, (int)(m_nx/ratio));
	m_nx = std::min(m_nx, nxy/2);
	m_ny = std::min(m_ny, nxy/2);

	double dx = box.getDX() / m_nx;
	double dy = box.getDY() / m_ny;

	m_grid = new QuadLeaf[m_nx*m_ny];
	
	// add regular grid vertices
	ControlNode2d qv(m_box.x0, m_box.y0);
	for(int j = 0; j <= m_ny; j++, qv.coord.y += dy){
		qv.coord.x = m_box.x0;
		for(int i = 0; i <= m_nx; i++, qv.coord.x += dx){
			m_grid_vertices.add(qv);
		}
	}

	// init regular grid boxes with coordinates and vertices indices
	DVector2d dpt(dx/4, dy/4);
	DPoint2d pt(m_box.x0 + dx / 2, m_box.y0 + dy / 2);
	for(int j = 0, k = 0, ix = 0; j < m_ny; j++, ix++){
		for(int i = 0; i < m_nx; i++, ix++){
			m_grid[k].setCoordinates(pt, dpt);
			pt.x += dx;
			m_grid[k].setNeighbours(((j>0)?(&m_grid[k-m_nx]):nullptr), ((i<m_nx-1)?(&m_grid[k+1]):nullptr),
				((j<m_ny-1)?(&m_grid[k+m_nx]):nullptr), ((i>0)?(&m_grid[k-1]):nullptr));
			m_grid[k++].setVertices(&(m_grid_vertices[ix]), &(m_grid_vertices[ix+1]), 
					&(m_grid_vertices[ix+m_nx+1]), &(m_grid_vertices[ix+m_nx+2]));
		}
		pt.y += dy;
		pt.x = m_box.x0 + dx / 2;
	}
}

ControlSpace2dQuadTree::~ControlSpace2dQuadTree()
{
	if(m_grid) delete[] m_grid;
}


/// Returns number of control nodes in adaptive control structure
int ControlSpace2dQuadTree::getControlNodesCount() const
{
	return m_grid_vertices.countInt();
}

/// Invoke function for all control nodes of this space (read-only)
void ControlSpace2dQuadTree::forEachControlNode(const std::function<void(const ControlNode2d & node)>& fg) const {
	for (int i = m_grid_vertices.countInt() - 1; i >= 0; i--)
		fg(m_grid_vertices[i]);
}

/// Invoke function for all control nodes of this space
void ControlSpace2dQuadTree::forEachControlNode(const std::function<void(ControlNode2d & node)>& fg) {
	for (int i = m_grid_vertices.countInt() - 1; i >= 0; i--)
		fg(m_grid_vertices[i]);
}

/////////////////////////////////////////////////////////////////
// Interpolate metric data for not-initialized quad tree nodes	
bool ControlSpace2dQuadTree::interpolate()
{
#ifdef STORE_CONTROL_EPS
	storeEPS("control-interp-before-gather");
//	SHOW_MESH("control space before gather", getViewSet());
#endif
	// calculate data from source points
	int leaf_count = m_nx * m_ny;
	for(int i = 0; i < leaf_count; i++)
		m_grid[i].gatherDataFromSourcePoints(base_surface);

	int gv_count = (int)m_grid_vertices.countInt();
	int ready_count = 0;
	for(int i = 0; i < gv_count; i++){
		ControlNode2d& qv = m_grid_vertices.get(i);
		if(qv.w > 0.0){
			qv.control_data /= qv.w;
			qv.w = -1.0; // proper initialization
		}
		if(qv.w < 0.0) ++ready_count;
	}
#ifdef STORE_CONTROL_EPS
	storeEPS("control-interp-gather");
#endif
	// + extended metric sources
	size_t ct = m_ext_source_points.countInt();
	for(size_t i = 0; i < ct; i++)
		ready_count += setMinControlValue(*(m_ext_source_points.get(i)));
#ifdef STORE_CONTROL_EPS
	storeEPS("control-interp-extended");
#endif
	// propagate control values from within leaves
	for(int i = 0; (ready_count < gv_count) && (i < leaf_count); i++)
		ready_count += m_grid[i].propagateUp();
#ifdef STORE_CONTROL_EPS
	storeEPS("control-interp-propagate-up");
#endif
	// extrapolate to other y1-grid nodes
	if(ready_count < gv_count)
		ready_count += extrapolateMainControlNodes();
#ifdef STORE_CONTROL_EPS
	storeEPS("control-interp-extrapolate");
#endif
	// propagate control values into leaves
	for(int i = 0; (ready_count < gv_count) && (i < leaf_count); i++)
		ready_count += m_grid[i].propagateDown();
#ifdef STORE_CONTROL_EPS
	storeEPS("control-interp-propagate-down");
#endif

	if(ready_count == gv_count) m_initialized = 1;

	return (ready_count == gv_count);
}

int ControlSpace2dQuadTree::extrapolateMainControlNodes()
{
	// Start from the initialized nodes
	int vx = m_nx+1;
	int vy = m_ny+1;
	int main_count = vx*vy;
	DataVector<int> ready_nodes(main_count);
	int *levels = new int[main_count];
	for(int i = 0; i < main_count; i++){
		ControlNode2d& qv = m_grid_vertices.get(i);
		if(qv.w < 0.0){
			ready_nodes.add(i);
			levels[i] = 0;
		}else
			levels[i] = main_count;
	}
	int new_nodes[4];
	int nearest_nodes[4];
	int start_ready_nodes = ready_nodes.countInt();
	if(start_ready_nodes == 0) return false;
	for(int i = 0; ready_nodes.countInt() < main_count; i++){
		int j, k = ready_nodes[i];
		int cur_level = levels[k];
		int ix = k % vx;
		int iy = k / vx;
		int ct = 0;
		// check neighbours of 'k'
		if(ix > 0 && levels[k-1] == main_count)		new_nodes[ct++] = k-1;
		if(ix < vx-1 && levels[k+1] == main_count)	new_nodes[ct++] = k+1;
		if(iy > 0 && levels[k-vx] == main_count)	new_nodes[ct++] = k-vx;
		if(iy < vy-1 && levels[k+vx] == main_count)	new_nodes[ct++] = k+vx;
		// if any new, calculate new value and add to ready_nodes
		for(int l = 0; l < ct; l++){
			k = new_nodes[l];
			ix = k % vx;
			iy = k / vx;
			int nct = 0;
			for(j = 1; ix-j >= 0 && levels[k-j] > cur_level; j++);	// x0
			if(ix-j >= 0) nearest_nodes[nct++] = k-j;
			for(j = 1; ix+j < vx && levels[k+j] > cur_level; j++);	// x1
			if(ix+j < vx) nearest_nodes[nct++] = k+j;
			for(j = 1; iy-j >= 0 && levels[k-vx*j] > cur_level; j++);	// down
			if(iy-j >= 0) nearest_nodes[nct++] = k-vx*j;
			for(j = 1; iy+j < vy && levels[k+vx*j] > cur_level; j++);	// x1
			if(iy+j < vy) nearest_nodes[nct++] = k+vx*j;
			assert(nct > 0);
			ControlNode2d& qv = m_grid_vertices.get(k);
			assert(qv.w == 0.0);
			DMetric2d dmp(base_surface, qv.coord);
			const DPoint2d middle = dmp.transformPStoRS(qv.coord);
			for(int m = 0; m < nct; m++){
				ControlNode2d& nqv = m_grid_vertices.get(nearest_nodes[m]);
				double dist = middle.distance2(dmp.transformPStoRS(nqv.coord));
				double nw = (dist < mesh_data.relative_small_number)?mesh_data.relative_infinity:(1.0/dist);
				qv.control_data += nqv.control_data * nw;
				qv.w += nw;
			}
			qv.control_data /= qv.w;
			qv.w = -3.0;
			levels[k] = cur_level+1;
			ready_nodes.add(k);
		}
	}
	delete[] levels;
	return (main_count - start_ready_nodes);
}

/////////////////////////////////////////////////////////////////////
// Add discrete metric information
void ControlSpace2dQuadTree::addControlNode(const ControlNode2d& node)
{
	assert(!m_initialized);
	findLastLeaf(node.coord)->adaptAndInsertControlPoint(
		m_grid_vertices, node, base_surface);
}

double ControlSpace2dQuadTree::getMetricGradationRatio(const DPoint2d& pt) const
{
	assert(m_initialized>0);
	const DPoint2d pt_fit = m_box.fitInPoint(pt);
	return findLastLeaf(pt_fit)->getMetricGradationRatio(pt_fit);
}

ReparameterizationData* ControlSpace2dQuadTree::getReparameterization(const DPoint2d& pt) const
{
	assert(m_initialized>0);
	const DPoint2d pt_fit = m_box.fitInPoint(pt);
	return findLastLeaf(pt_fit)->getReparameterization(pt_fit);
}

ControlDataMatrix2d ControlSpace2dQuadTree::getMetricAtPoint(const DPoint2d& pt) const
{
	assert(m_initialized>0);
	const DPoint2d pt_fit = m_box.fitInPoint(pt);
	return findLastLeaf(pt_fit)->getMetricAtPoint(pt_fit);
}

double ControlSpace2dQuadTree::interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const
{
	const DPoint2d pt_fit = m_box.fitInPoint(pt);
	return findLastLeaf(pt_fit)->interpolateDoubleTag(pt_fit, type);
}

void ControlSpace2dQuadTree::setMinDoubleTag(const DPoint2d& pt, TagExtended::TagType type, double q)
{
	const DPoint2d pt_fit = m_box.fitInPoint(pt);
	QuadLeaf* leaf = findLastLeaf(pt_fit);
	if(!leaf) return;
	for(int i = 0; i < 4; i++){
		ControlNode2d* node = leaf->getVertex(i);
		double v = node->getDoubleTag(type);
		if(q > v) node->setDoubleTag(type, q);
	}
}

void ControlSpace2dQuadTree::setMaxDoubleTag(const DPoint2d& pt, TagExtended::TagType type, double q)
{
	const DPoint2d pt_fit = m_box.fitInPoint(pt);
	QuadLeaf* leaf = findLastLeaf(pt_fit);
	if(!leaf) return;
	for(int i = 0; i < 4; i++){
		ControlNode2d* node = leaf->getVertex(i);
		double v = node->getDoubleTag(type);
		if(q < v) node->setDoubleTag(type, q);
	}
}

ControlDataMatrix2d ControlSpace2dQuadTree::getMetricAndParameterizationAtPoint(
	const DPoint2d& pt, ControlDataMatrix2d& p_cdm) const 
{
	assert(m_initialized>1);
	const DPoint2d pt_fit = m_box.fitInPoint(pt);
	return findLastLeaf(pt_fit)->getMetricAndParameterizationAtPoint(pt_fit, p_cdm, base_surface);
}

void ControlSpace2dQuadTree::setSurfaceCurvatureControlData()
{
	if(!base_surface) return;

	assert(!m_initialized);

	double model_diameter = m_box.getDiameter();
	double max_len = model_diameter * std::min(param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
	double min_len = std::max(model_diameter * param_min_diameter_ratio, param_min_length);

	// count control data for base grid vertices
	int vertex_count = (int)m_grid_vertices.countInt();
	for(int i = 0; i < vertex_count; i++){
		double g_ratio;
		const SurfaceCurvature curvature = base_surface->getCurvature(m_grid_vertices[i].coord, g_ratio);
		if(g_ratio > MIN_PARAM_GRATIO){ 
			ControlDataStretch2d data = adjustCurvatureData(curvature, param_curvature_ratio,
				model_diameter, min_len, max_len, param_stretch_max_ratio);
			// Convert (lx, ly, alfa) -> (dxx, dxy, dyy)
			m_grid_vertices[i].control_data = DMetric2d::stretchToMatrix(data);
			m_grid_vertices[i].w = -1.0;
		}else{
			LOG4CPLUS_WARN(MeshLog::logger_console, "Improper parameterization: " << g_ratio);
		}
	}

	int leaf_count = m_nx * m_ny;
	for(int i = 0; i < leaf_count; i++)
		m_grid[i].adaptToSurfaceCurvature(base_surface, m_grid_vertices, 
				param_curvature_ratio, param_stretch_max_ratio, 
				model_diameter, min_len, max_len);

#ifdef STORE_CONTROL_EPS
	storeEPS("control-surface-curvature");
#endif

	interpolate();
}

bool ControlSpace2dQuadTree::setMinControl(const DPoint2d& pt, const ControlDataMatrix2d& cdm, bool min_value_set)
{
	assert(m_initialized);
	if(!m_box.contains(pt)) return false;
	return findLastLeaf(pt)->adaptAndSetControlPoint(m_grid_vertices, 
				ControlNode2d(pt, cdm), base_surface, min_value_set);
}

/**
 * Calculates the int coordinates of x0-y0 box (of matrix) containg the given point
 */
void ControlSpace2dQuadTree::countLocalCoordinates(const DPoint2d &pt, int &ix, int &iy) const
{
	ix = (int) (m_nx * (pt.x - m_box.x0) / (m_box.x1 - m_box.x0));
	iy = (int) (m_ny * (pt.y - m_box.y0) / (m_box.y1 - m_box.y0));
	if(ix > (m_nx - 1)) ix = m_nx - 1;
	if(iy > (m_ny - 1)) iy = m_ny - 1;
	if(ix < 0) ix = 0;
	if(iy < 0) iy = 0;
}

/**
 * Calculates the int coordinates of x0-y0 box (of matrix) containg the given point
 */
int ControlSpace2dQuadTree::countLocalCoordinates(const DPoint2d &pt) const
{
	int ix = (int) (m_nx * (pt.x - m_box.x0) / (m_box.x1 - m_box.x0));
	int iy = (int) (m_ny * (pt.y - m_box.y0) / (m_box.y1 - m_box.y0));
	if(ix > (m_nx - 1)) ix = m_nx - 1;
	if(iy > (m_ny - 1)) iy = m_ny - 1;
	if(ix < 0) ix = 0;
	if(iy < 0) iy = 0;
	return iy * m_nx + ix;
}

/// Draws the structure of the control space
void ControlSpace2dQuadTree::storeEPS(const char* name, int id)
{
	ostringstream fname;
	fname << name << '-' << id << ".eps";

	DRect rect = m_box;
	rect.inflate(0.05);
	EPSFile eps(fname.str(), rect.x0, rect.x1, rect.y0, rect.y1);

	double dx = m_box.getDX() / m_nx;
	double dy = m_box.getDY() / m_ny;

	// main grid
	DPoint2d ptx0(m_box.x0, m_box.y0);
	DPoint2d ptx1(m_box.x0, m_box.y1);
	for(int j = 0; j <= m_nx; j++){
		eps.drawLine(ptx0, ptx1);
		ptx0.x += dx;
		ptx1.x += dx;
	}
	DPoint2d pty0(m_box.x0, m_box.y0);
	DPoint2d pty1(m_box.x1, m_box.y0);
	for(int j = 0; j <= m_ny; j++){
		eps.drawLine(pty0, pty1);
		pty0.y += dy;
		pty1.y += dy;
	}
	// inside
	int count = m_nx * m_ny;
	for(int k = 0; k < count; k++){
		m_grid[k].storeEPS(&eps);
	}
	// nodes
	int gv_count = m_grid_vertices.countInt();
	int main_count = (m_nx+1)*(m_ny+1);
	for(int i = 0; i < gv_count; i++){
		const ControlNode2d& qv = m_grid_vertices.get(i);
		double gray_color = 1.0;
		if(qv.w < 0) gray_color = -qv.w*0.2-0.2;
		eps.drawPoint(qv.coord, (i<main_count)?0.4:0.1, gray_color);
	}
}

/// Returns the "screenshot" of this control space for visualization
MeshViewSet* ControlSpace2dQuadTree::getViewSet(MeshViewSet* set, bool with_surface, TagExtended::TagType tag) const
{
	if(set){
		set->prepareFreePlace((int)m_grid_vertices.countInt(), 0, 4*m_nx*m_ny, 0);
	}else{
		set = new MeshViewSet((int)m_grid_vertices.countInt(), 0, 4*m_nx*m_ny, 0);
	}

	double dx = m_box.getDX() / m_nx;
	double dy = m_box.getDY() / m_ny;

	// main grid
	DPoint2d ptx0(m_box.x0, m_box.y0);
	DPoint2d ptx1(m_box.x0, m_box.y1);
	for(int j = 0; j <= m_nx; j++){
		if(with_surface){
			set->addEdge(
				base_surface->getPoint(ptx0),
				base_surface->getPoint(ptx1));
		}else{
			set->addEdge(DPoint3d(ptx0, 0.0), DPoint3d(ptx1, 0.0));
		}
		ptx0.x += dx;
		ptx1.x += dx;
	}
	DPoint2d pty0(m_box.x0, m_box.y0);
	DPoint2d pty1(m_box.x1, m_box.y0);
	for(int j = 0; j <= m_ny; j++){
		if(with_surface){
			set->addEdge(
				base_surface->getPoint(pty0),
				base_surface->getPoint(pty1));
		}else{
			set->addEdge(DPoint3d(pty0, 0.0), DPoint3d(pty1, 0.0));
		}
		pty0.y += dy;
		pty1.y += dy;
	}
	// inside
	int count = m_nx * m_ny;
	for(int k = 0; k < count; k++){
		m_grid[k].addToViewSet(set, with_surface ? base_surface : nullptr, tag);
	}

	m_grid_vertices.forEach([set](const ControlNode2d& qv) {
		if (qv.w < 0)
			set->addPoint(DPoint3d(qv.coord.x, qv.coord.y, 0.0), (int)(-qv.w));
	});

	return set;
}

ControlDataMatrix2d ControlSpace2dQuadTree::QuadLeaf::getMetricAtPoint(const DPoint2d& pt) const
{
	assert(!isSplit());

	double tx = (pt.x - m_middle.x + 2*m_dl.x) / (4*m_dl.x);
	double ty = (pt.y - m_middle.y + 2*m_dl.y) / (4*m_dl.y);
	if(param_use_midnodes){
		double coeff[8] = {	// Shape functions for eight-node quadrilateral
			(1-tx)*(1-ty)*(-2*tx-2*ty+1),	// VX2D_SW
			tx*(1-ty)*(2*tx-2*ty-1),		// VX2D_SE
			(1-tx)*ty*(2*ty-2*tx-1),		// VX2D_NW
			tx*ty*(2*tx+2*ty-3),			// VX2D_NE
			4*(1-tx)*(1-ty)*tx,			// MVX2D_DOWN
			4*(1-ty)*tx*ty,				// MVX2D_RIGHT 
			4*(1-tx)*tx*ty,				// MVX2D_UP 
			4*(1-ty)*(1-tx)*ty			// MVX2D_LEFT
		};
		ControlDataMatrix2d ave;
		// Add contribution from middle nodes
		//	- if any missing, redistribute its weight to neighbours
		for(int i = MVX2D_FIRST; i <= MVX2D_LAST; i++){
			if(m_mid_vertices[i])
				ave += m_mid_vertices[i]->control_data * coeff[4+i];
			else{
				coeff[mid_to_nodes[i][0]] += 0.5*coeff[4+i];
				coeff[mid_to_nodes[i][1]] += 0.5*coeff[4+i];
			}
		}
		for(int i = VX2D_FIRST; i <= VX2D_LAST; i++)
			ave += m_vertices[i]->control_data * coeff[i];
		return ave;
	}else{
		double coeff[4] = {(1-tx)*(1-ty), // VX2D_SW
			tx *(1-ty), // VX2D_SE
			(1-tx)* ty, // VX2D_NW
			tx * ty};	// VX2D_NE
		ControlDataMatrix2d ave = m_vertices[0]->control_data * coeff[0];
		for(int i = 1; i < 4; i++)
			ave += m_vertices[i]->control_data * coeff[i];
		return ave;
	}
}

double ControlSpace2dQuadTree::QuadLeaf::interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const
{
	assert(!isSplit());

	double tx = (pt.x - m_middle.x + 2*m_dl.x) / (4*m_dl.x);
	double ty = (pt.y - m_middle.y + 2*m_dl.y) / (4*m_dl.y);
	if(param_use_midnodes){
		double coeff[8] = {	// Shape functions for eight-node quadrilateral
			(1-tx)*(1-ty)*(-2*tx-2*ty+1),	// VX2D_SW
			tx*(1-ty)*(2*tx-2*ty-1),		// VX2D_SE
			(1-tx)*ty*(2*ty-2*tx-1),		// VX2D_NW
			tx*ty*(2*tx+2*ty-3),			// VX2D_NE
			4*(1-tx)*(1-ty)*tx,			// MVX2D_DOWN
			4*(1-ty)*tx*ty,				// MVX2D_RIGHT 
			4*(1-tx)*tx*ty,				// MVX2D_UP 
			4*(1-ty)*(1-tx)*ty			// MVX2D_LEFT
		};
		double ave = 0.0;
		// Add contribution from middle nodes
		//	- if any missing, redistribute its weight to neighbours
		for(int i = MVX2D_FIRST; i <= MVX2D_LAST; i++){
			if(m_mid_vertices[i])
				ave += m_mid_vertices[i]->getDoubleTag(type) * coeff[4+i];
			else{
				coeff[mid_to_nodes[i][0]] += 0.5*coeff[4+i];
				coeff[mid_to_nodes[i][1]] += 0.5*coeff[4+i];
			}
		}
		for(int i = VX2D_FIRST; i <= VX2D_LAST; i++)
			ave += m_vertices[i]->getDoubleTag(type) * coeff[i];
		return ave;
	}else{
		double coeff[4] = {(1-tx)*(1-ty), // VX2D_SW
			tx *(1-ty), // VX2D_SE
			(1-tx)* ty, // VX2D_NW
			tx * ty};	// VX2D_NE
		double ave = m_vertices[0]->getDoubleTag(type) * coeff[0];
		for(int i = 1; i < 4; i++)
			ave += m_vertices[i]->getDoubleTag(type) * coeff[i];
		return ave;
	}
}

double ControlSpace2dQuadTree::QuadLeaf::getMetricGradationRatio(const DPoint2d& pt) const
{
	assert(!isSplit());

	double tx = (pt.x - m_middle.x + 2*m_dl.x) / (4*m_dl.x);
	double ty = (pt.y - m_middle.y + 2*m_dl.y) / (4*m_dl.y);

	double coeff[4] = {(1-tx)*(1-ty), // VX2D_SW
		tx *(1-ty), // VX2D_SE
		(1-tx)* ty, // VX2D_NW
		tx * ty};	// VX2D_NE
	double ave = m_vertices[0]->max_gradation_ratio * coeff[0];
	for(int i = 1; i < 4; i++)
		ave += m_vertices[i]->max_gradation_ratio * coeff[i];
	return ave;
}

ReparameterizationData* ControlSpace2dQuadTree::QuadLeaf::getReparameterization(const DPoint2d& /* pt */) const
{
	assert(!isSplit());
	return m_reparameterization;
}

ControlSpace2dQuadTree::QuadLeaf::QuadLeaf() 
	: m_reparameterization(nullptr), m_level(0), m_source_points(nullptr) 
{ 
	m_leaves[0] = nullptr; 
	m_neighbours[0] = m_neighbours[1] = m_neighbours[2] = m_neighbours[3] = nullptr;
	m_mid_vertices[0] = m_mid_vertices[1] = m_mid_vertices[2] = m_mid_vertices[3] = nullptr;
}

ControlSpace2dQuadTree::QuadLeaf::QuadLeaf(double mid_x, double mid_y, double dx, double dy, int level) 
	: m_middle(mid_x, mid_y), m_reparameterization(nullptr), 
	m_level(level), m_dl(dx, dy), m_source_points(nullptr) 
{
	m_leaves[0] = nullptr; 
	m_neighbours[0] = m_neighbours[1] = m_neighbours[2] = m_neighbours[3] = nullptr;
	m_mid_vertices[0] = m_mid_vertices[1] = m_mid_vertices[2] = m_mid_vertices[3] = nullptr;

}

ControlSpace2dQuadTree::QuadLeaf::~QuadLeaf()
{ 
	if(m_leaves[0]) 
		for(int i = 0; i < 4; i++) delete m_leaves[i]; 
	if(m_reparameterization) delete m_reparameterization;
}

ControlDataMatrix2d ControlSpace2dQuadTree::QuadLeaf::getMetricAndParameterizationAtPoint(
	const DPoint2d& pt, ControlDataMatrix2d& p_cdm, SurfaceConstPtr /* surface */) const
{
	if(ControlSpace2dAdaptive::param_cached_parameterization_matrix == 1 ||
		m_reparameterization != nullptr){	// or don't assume as constant ? (i.e. if not planar reparameterization surface)
		assert(m_middle_cdm.det() > 0.0);
		p_cdm = m_middle_cdm;
		return getMetricAtPoint(pt);
	}
	assert(ControlSpace2dAdaptive::param_cached_parameterization_matrix == 2); // count param_matrix from nodes
	assert(!isSplit());

	double tx = (pt.x - m_middle.x + 2*m_dl.x) / (4*m_dl.x);
	double ty = (pt.y - m_middle.y + 2*m_dl.y) / (4*m_dl.y);
	if(param_use_midnodes){
		double coeff[8] = {	// Shape functions for eight-node quadrilateral
			(1-tx)*(1-ty)*(-2*tx-2*ty+1),	// VX2D_SW
			tx*(1-ty)*(2*tx-2*ty-1),		// VX2D_SE
			(1-tx)*ty*(2*ty-2*tx-1),		// VX2D_NW
			tx*ty*(2*tx+2*ty-3),			// VX2D_NE
			4*(1-tx)*(1-ty)*tx,			// MVX2D_DOWN
			4*(1-ty)*tx*ty,				// MVX2D_RIGHT 
			4*(1-tx)*tx*ty,				// MVX2D_UP 
			4*(1-ty)*(1-tx)*ty			// MVX2D_LEFT
		};
		ControlDataMatrix2d ave;
		p_cdm.reset();
		// Add contribution from middle nodes
		//	- if any missing, redistribute its weight to neighbours
		for(int i = MVX2D_FIRST; i <= MVX2D_LAST; i++){
			if(m_mid_vertices[i]){
				ave += m_mid_vertices[i]->control_data * coeff[4+i];
				p_cdm += m_mid_vertices[i]->param_data * coeff[4+i];
			}else{
				coeff[mid_to_nodes[i][0]] += 0.5*coeff[4+i];
				coeff[mid_to_nodes[i][1]] += 0.5*coeff[4+i];
			}
		}
		for(int i = VX2D_FIRST; i <= VX2D_LAST; i++){
			ave += m_vertices[i]->control_data * coeff[i];
			p_cdm += m_vertices[i]->param_data * coeff[i];
		}
		return ave;
	}else{
		double coeff[4] = {(1-tx)*(1-ty), // VX2D_SW
			tx *(1-ty), // VX2D_SE
			(1-tx)* ty, // VX2D_NW
			tx * ty};	// VX2D_NE
		ControlDataMatrix2d ave = m_vertices[0]->control_data * coeff[0];
		p_cdm = m_vertices[0]->param_data * coeff[0];
		//assert(m_vertices[0]->param_data.det() > 0.0);
		double total_p_coeff = (m_vertices[0]->param_data.det() > 0.0) ? coeff[0] : 0.0;
		for(int i = 1; i < 4; i++){
			ave += m_vertices[i]->control_data * coeff[i];
			p_cdm += m_vertices[i]->param_data * coeff[i];
			//assert(m_vertices[i]->param_data.det() > 0.0);
			if(m_vertices[i]->param_data.det() > 0.0)
				total_p_coeff += coeff[i];
		}
		if(total_p_coeff > 0.0) p_cdm /= total_p_coeff;
		else p_cdm = m_middle_cdm;
		return ave;
	}
}

bool ControlSpace2dQuadTree::QuadLeaf::checkAndSetReparameterization(SurfaceConstPtr surface)
{
	if(m_reparameterization) return true; // already there

	return false;
/*
	// calculate plane
	const DPoint3d pts[4] = {
		surface->getPoint(m_vertices[0]->coord),
		surface->getPoint(m_vertices[1]->coord),
		surface->getPoint(m_vertices[2]->coord),
		surface->getPoint(m_vertices[3]->coord)
	};
	const DVector3d dvs_s[2] = { pts[VX2D_SE] - pts[VX2D_SW], pts[VX2D_NE] - pts[VX2D_NW] };
	const DVector3d dvs_t[2] = { pts[VX2D_NE] - pts[VX2D_SE], pts[VX2D_NW] - pts[VX2D_SW] };

	const DVector3d dv_s = ((dvs_s[0].length2() > dvs_s[1].length2()) ? dvs_s[0] : dvs_s[1]).normalized();
	const DVector3d dv_t = ((dvs_t[0].length2() > dvs_t[1].length2()) ? dvs_t[0] : dvs_t[1]);
	const DVector3d dv_td = dv_s.crossProduct(dv_t).crossProduct(dv_s).normalized();
	SurfacePlane* plane = new SurfacePlane(surface->getPoint(m_middle), dv_s, dv_td);

	// check if OK
	double total_min_h = mesh_data.relative_infinity;
	double cf = sqr(ControlSpace2dAdaptive::param_curvature_ratio);
	for(int i = VX2D_FIRST; i <= VX2D_LAST; i++){
		double min_h = sqr(m_vertices[i]->control_data.minEigenvalue());
		const DPoint2d pt_planar = plane->getParameters(pts[i]);
		const DPoint3d pt_approximated = plane->getPoint(pt_planar);
		double dist = pts[i].distance2(pt_approximated);
		if(cf * dist > min_h){
			delete plane;
			return false;
		}
		total_min_h = std::min(min_h, total_min_h);
	}

	m_reparameterization = new ReparameterizationData(
		plane, m_middle, DPoint2d(0.0, 0.0), total_min_h);
	double g_ratio;
	m_middle_cdm = plane->countParameterizationMatrix(m_middle, g_ratio);
	return true;
*/
}

void ControlSpace2dQuadTree::QuadLeaf::adaptToSurfaceCurvature(
	SurfaceConstPtr surface, QuadVertexArray& grid_vertices, double c_ratio, 
	double stretch_max_ratio, double model_diameter, double min_len, double max_len)
{
	for(int i = LF_FIRST; i <= LF_LAST; i++)
		if(m_vertices[i]->w > -1.0) return;	// 
	for(int i = MVX2D_FIRST; i <= MVX2D_LAST; i++)
		if(m_mid_vertices[i] && m_mid_vertices[i]->w > -1.0) return;

	if(!isSplit()){
		double g_ratio;
		const SurfaceCurvature curvature = surface->getCurvature(m_middle, g_ratio);
		if(g_ratio < MIN_PARAM_GRATIO){
			LOG4CPLUS_WARN(MeshLog::logger_console, "Improper parameterization: " << g_ratio);
			return;
		}
		ControlDataStretch2d data = ControlSpace2dAdaptive::adjustCurvatureData(curvature, c_ratio,
			model_diameter, min_len, max_len, stretch_max_ratio);

		// Przeliczenie (lx, ly, alfa) -> (dxx, dxy, dyy)
		ControlDataMatrix2d data_calculated = DMetric2d::stretchToMatrix(data);
		ControlDataMatrix2d data_assessed = getMetricAtPoint(m_middle);

		double diff = data_calculated.countDifferenceRR(data_assessed);
		if(diff > param_threshold_diff && 
			m_level < ControlSpace2dQuadTree::param_max_depth_surface_adaptation)
		{
			// split
			QuadVertexList split_vertices(100);
			split(grid_vertices, surface, &split_vertices);
			size_t split_vertices_count = split_vertices.countInt();
			for(size_t i = 0; i < split_vertices_count; i++){ 
				const SurfaceCurvature curv = surface->getCurvature(split_vertices[i]->coord, g_ratio);
				if(g_ratio > MIN_PARAM_GRATIO){
					ControlDataStretch2d d = ControlSpace2dAdaptive::adjustCurvatureData(curv, c_ratio,
						model_diameter, min_len, max_len, stretch_max_ratio);
					split_vertices[i]->control_data = DMetric2d::stretchToMatrix(d);
					split_vertices[i]->w = -1.0;
				}else{
					LOG4CPLUS_WARN(MeshLog::logger_console, "Improper parameterization: " << g_ratio);
				}
			}
		}
	}
	if(isSplit()){ // which means also all m_leaves[i] != nullptr 
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			m_leaves[i]->adaptToSurfaceCurvature(surface, 
				grid_vertices, c_ratio, stretch_max_ratio,
				model_diameter, min_len, max_len);
	}
}

void ControlSpace2dQuadTree::QuadLeaf::adaptAndInsertControlPoint(
	QuadVertexArray& grid_vertices, const ControlNode2d& mqv, 
	SurfaceConstPtr surface)
{
	if(isSplit()){
		findLastLeaf(mqv.coord)->adaptAndInsertControlPoint(grid_vertices, mqv, surface);
		return;
	}
	DMetric2d dmp(mqv.control_data, surface, mqv.coord);
	double len_x = dmp.transformPStoMS(DVector2d(4*m_dl.x, 0.0)).length2();
	double len_y = dmp.transformPStoMS(DVector2d(0.0, 4*m_dl.y)).length2();
	bool split_needed = (len_x > METRIC_LENGTH_RATIO2 || len_y > METRIC_LENGTH_RATIO2 );

	ControlNode2d qv = mqv;
	if(!split_needed && m_source_points){	// further check
		// compare new source node with already existing
		// check with:
		//	 all points (used now) 
		//	(?) average, minimum, some kind of approximation
		double max_diff = 0.0;
		for(size_t i = 0; i < m_source_points->countInt(); i++){
			const ControlNode2d& sqv = m_source_points->get(i);
			if(dmp.transformPStoMS(qv.coord - sqv.coord).length() < METRIC_SMALL_NUMBER){
				// combine
				qv.control_data.setMinimum(sqv.control_data);
				m_source_points->removeAt(i);
				max_diff = 0.0;
				i = -1; // start again
				continue;
			}
			double diff = qv.control_data.countDifferenceRR(sqv.control_data);
			if(diff > max_diff){
				max_diff = diff;
				if(max_diff > param_threshold_diff) break;
			}
		}
		split_needed = (max_diff > param_threshold_diff);
	}

	if(!m_source_points) m_source_points = new DataVector<ControlNode2d>;
	m_source_points->add(qv);

	if(split_needed){
		split(grid_vertices, surface); // split and redistribute source points
	}
}

bool ControlSpace2dQuadTree::QuadLeaf::adaptAndSetControlPoint(
	QuadVertexArray& grid_vertices, 
	const ControlNode2d& qv, SurfaceConstPtr surface, bool min_value_set)
{
	if(isSplit()){
		return findLastLeaf(qv.coord)->adaptAndSetControlPoint(grid_vertices, qv, surface, min_value_set);
	}
	// compare new source node with already existing
	// check with calculated value:
	const ControlDataMatrix2d cdm = getMetricAtPoint(qv.coord);
	double diff = qv.control_data.countDifferenceRR(cdm);
	bool any_changes = false;
	if(diff > param_threshold_diff){
		// split may be unnecessary, if qv.control_data is greater than already set for this leaf...
		ControlDataMatrix2d min_cdm = cdm;
		min_cdm.setMinimum(qv.control_data);
		if(cdm.countDifferenceRR(min_cdm) < param_threshold_diff){
			// qv different from current metric, but introduces no changes - so skip
			return false;
		}
		DMetric2d dmp(qv.control_data, surface, qv.coord);
		double len_x = dmp.transformPStoMS(DVector2d(4*m_dl.x, 0.0)).length2();
		double len_y = dmp.transformPStoMS(DVector2d(0.0, 4*m_dl.y)).length2();
		bool split_needed = (len_x > METRIC_LENGTH_RATIO2 || len_y > METRIC_LENGTH_RATIO2 );
		if(split_needed){
			if(split(grid_vertices, surface, nullptr, true)){ // split (+ set values for new nodes)
				findLastLeaf(qv.coord)->adaptAndSetControlPoint(grid_vertices, qv, surface, min_value_set);
				return true;
			}
		}
		if(min_value_set){
			for(int i = VX2D_FIRST; i <= VX2D_LAST; i++)
				any_changes |= m_vertices[i]->control_data.setMinimum(qv.control_data);
		}
	}

	return any_changes;
}

void ControlSpace2dQuadTree::QuadLeaf::setNeighbour(int side, QuadLeaf* nb, bool skip_first)
{
	if(!skip_first) m_neighbours[side] = nb;
	if(isSplit()){
		m_leaves[mid_to_nodes[side][0]]->setNeighbour(side, nb);
		m_leaves[mid_to_nodes[side][1]]->setNeighbour(side, nb); 
	}
}

bool ControlSpace2dQuadTree::QuadLeaf::split(QuadVertexArray& grid_vertices, 
		SurfaceConstPtr surface, QuadVertexList* split_vertices, bool init_values)
{
	if(m_level >= ControlSpace2dQuadTree::param_max_depth)
		return false;

	assert(m_leaves[LF_FIRST] == nullptr);

	// Balance before split...
	balance(grid_vertices, surface, split_vertices, init_values);
	
	m_leaves[LF_SW] = new QuadLeaf(m_middle.x-m_dl.x, m_middle.y-m_dl.y, m_dl.x/2, m_dl.y/2, m_level+1);
	m_leaves[LF_SE] = new QuadLeaf(m_middle.x+m_dl.x, m_middle.y-m_dl.y, m_dl.x/2, m_dl.y/2, m_level+1);
	m_leaves[LF_NW] = new QuadLeaf(m_middle.x-m_dl.x, m_middle.y+m_dl.y, m_dl.x/2, m_dl.y/2, m_level+1);
	m_leaves[LF_NE] = new QuadLeaf(m_middle.x+m_dl.x, m_middle.y+m_dl.y, m_dl.x/2, m_dl.y/2, m_level+1);

	ControlNode2d* qv_middle;

	// middle
	size_t index = grid_vertices.add(ControlNode2d(m_middle));
	qv_middle = &(grid_vertices[index]);
	if(split_vertices) split_vertices->add(qv_middle);
	if(init_values){
		double g_ratio;
		qv_middle->control_data = (m_vertices[VX2D_SW]->control_data + m_vertices[VX2D_SE]->control_data +
			m_vertices[VX2D_NW]->control_data + m_vertices[VX2D_NE]->control_data)*0.25;
		qv_middle->param_data = surface->countParameterizationMatrix(m_middle, g_ratio);
		qv_middle->w = -5.0;
		for(int i = LF_FIRST; i <= LF_LAST; i++){
			m_leaves[i]->m_middle_cdm = surface->countParameterizationMatrix(m_middle, g_ratio);
			if(g_ratio < MIN_PARAM_GRATIO) m_leaves[i]->m_middle_cdm = m_middle_cdm;
		}
	}

	// check middle nodes at edges
	QuadMidVertexWhich opposite[4] = {MVX2D_UP, MVX2D_LEFT, MVX2D_DOWN, MVX2D_RIGHT};
	for(int i = MVX2D_FIRST; i <= MVX2D_LAST; i++){
		if(!m_mid_vertices[i]){	// vertex missing
			ControlNode2d *q0 = m_vertices[mid_to_nodes[i][0]];
			ControlNode2d *q1 = m_vertices[mid_to_nodes[i][1]];
			index = grid_vertices.add(ControlNode2d(DPoint2d::average(q0->coord, q1->coord)));
			m_mid_vertices[i] = &(grid_vertices[index]);
			if(split_vertices) split_vertices->add(m_mid_vertices[i]);
			if(m_neighbours[i] && (m_neighbours[i]->m_mid_vertices[opposite[i]] == nullptr))
				m_neighbours[i]->m_mid_vertices[opposite[i]] = m_mid_vertices[i];
			if(init_values){
				double g_ratio;
				m_mid_vertices[i]->control_data = (q0->control_data + q1->control_data)*0.5;
				m_mid_vertices[i]->param_data = 
					surface->countParameterizationMatrix(m_mid_vertices[i]->coord, g_ratio);
				m_mid_vertices[i]->w = -5.0;
			}
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
	if(m_source_points){	// split too
		size_t ct = m_source_points->countInt();
		size_t ct_min = std::min((size_t)10, ct / 4);
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			m_leaves[i]->m_source_points = new DataVector<ControlNode2d>(ct_min);
		for(size_t i = 0; i < ct; i++){
			const ControlNode2d& qv = m_source_points->get(i);
			findLastLeaf(qv.coord)->adaptAndInsertControlPoint(grid_vertices, qv, surface);
		}
		delete m_source_points;
		m_source_points = nullptr;
	}
	return true;
}

void ControlSpace2dQuadTree::QuadLeaf::balance(QuadVertexArray& grid_vertices,
		SurfaceConstPtr surface, QuadVertexList* split_vertices, bool init_values)
{
	for(int i = NB_FIRST; i <= NB_LAST; i++){
		if(!m_neighbours[i]) continue;	// border of qtree
		while((m_level - m_neighbours[i]->m_level) >= param_balance_level)
			if(!m_neighbours[i]->split(grid_vertices, surface, split_vertices, init_values)) break;
	}
}

ControlSpace2dQuadTree::QuadLeaf* ControlSpace2dQuadTree::QuadLeaf::findLastLeaf(const DPoint2d& pt)
{
	if(m_leaves[LF_FIRST] == nullptr) return this;
	if(pt.x < m_middle.x) 
		if(pt.y < m_middle.y) return m_leaves[LF_SW]->findLastLeaf(pt);
		else return m_leaves[LF_NW]->findLastLeaf(pt);
	else if(pt.y < m_middle.y) return m_leaves[LF_SE]->findLastLeaf(pt);
		else return m_leaves[LF_NE]->findLastLeaf(pt);
}

void ControlSpace2dQuadTree::QuadLeaf::storeEPS(EPSFile* eps_file) const
{
	if(isSplit()){
		eps_file->drawLine(DPoint2d(m_middle.x, m_middle.y - 2*m_dl.y), 
			DPoint2d(m_middle.x, m_middle.y + 2*m_dl.y));
		eps_file->drawLine(DPoint2d(m_middle.x - 2*m_dl.x, m_middle.y), 
			DPoint2d(m_middle.x + 2*m_dl.x, m_middle.y));
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			m_leaves[i]->storeEPS(eps_file);
	}

	if(m_source_points){
		m_source_points->forEach([eps_file](const ControlNode2d& cn) {
			eps_file->drawPointCross(cn.coord, 0.4);
		});
	}
}

void ControlSpace2dQuadTree::QuadLeaf::addToViewSet(MeshViewSet* set, SurfaceConstPtr surface, 
				TagExtended::TagType tag) const
{
	if(isSplit()){
		if(surface){
			set->addEdge(
				surface->getPoint(DPoint2d(m_middle.x, m_middle.y - 2*m_dl.y)), 
				surface->getPoint(DPoint2d(m_middle.x, m_middle.y + 2*m_dl.y)));
			set->addEdge(
				surface->getPoint(DPoint2d(m_middle.x - 2*m_dl.x, m_middle.y)), 
				surface->getPoint(DPoint2d(m_middle.x + 2*m_dl.x, m_middle.y)));
		}else{
			set->addEdge(DPoint3d(m_middle.x, m_middle.y - 2*m_dl.y, 0.0), 
					DPoint3d(m_middle.x, m_middle.y + 2*m_dl.y, 0.0));
			set->addEdge(DPoint3d(m_middle.x - 2*m_dl.x, m_middle.y, 0.0), 
					DPoint3d(m_middle.x + 2*m_dl.x, m_middle.y, 0.0));
		}
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			m_leaves[i]->addToViewSet(set, surface, tag);
	}else{
		auto data = std::make_shared<MeshViewFaceData>(4);
		data->area_id = 0;
		data->quality = 1.0;

		DPoint3d dpts[4];
		int order[4] = {0, 1, 3, 2};
		for(int i = 0; i < 4; i++){
			if(surface){
				dpts[i] = surface->getPoint(m_vertices[order[i]]->coord);
			}else{
				dpts[i] = DPoint3d(m_vertices[order[i]]->coord, 0.0);
			}
			if(tag == TagExtended::TAG_NONE){
				data->wpts[i] = 1.0f - 0.1f * (float)m_vertices[order[i]]->max_gradation_ratio;
			}else{
				if(m_vertices[order[i]]->availableTag(tag)) 
					data->wpts[i] = (float)m_vertices[order[i]]->getDoubleTag(tag);
				else 
					data->wpts[i] = 1.0f;
			}
		}

		DPoint3d middle(dpts[0], dpts[2], 0.5);

		for(int i = 0; i < 4; i++)
			data->pts[i] = middle + (dpts[i] - middle) * 0.95;

		data->normal = (data->pts[1] - data->pts[0]).crossProduct(data->pts[2] - data->pts[0]).normalized();

		set->m_faces.add(data);
	}

	if(m_source_points){
		m_source_points->forEach([set](const ControlNode2d& cn) {
			set->addPoint(DPoint3d(cn.coord.x, cn.coord.y, 0.0), -1);
		});
	}
}

ControlSpace2dQuadTree::QuadLeaf* ControlSpace2dQuadTree::findLastLeaf(const DPoint2d& pt) const
{
	int k = countLocalCoordinates(pt);
	return m_grid[k].findLastLeaf(pt);
}

void ControlSpace2dQuadTree::testQuadTree()
{
/*
	storeEPS("control", 0);

	DPoint2d test_point(0.001, 0.001);

	for(int j = 1; j < 10; j++){
		int k = countLocalCoordinates(test_point);
		QuadLeaf* leaf = m_grid[k].findLastLeaf(test_point);
		QuadVertexList split_vertices(100);
		leaf->split(m_grid_vertices, split_vertices, base_surface.ptr());
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "Control QuadTree split -> " << split_vertices.countInt() << 
			" new vertices." << endl;
		storeEPS("control", j);
	}
*/
}

void ControlSpace2dQuadTree::adaptToParameterization()
{
	assert(m_initialized);
	///
	size_t gv_count = m_grid_vertices.countInt();
	for(size_t i = 0; i < gv_count; i++){
		ControlNode2d& node = m_grid_vertices.get(i);
		double g_ratio;
		node.param_data = base_surface->countParameterizationMatrix(node.coord, g_ratio);
		if(g_ratio < MIN_PARAM_GRATIO) 
			node.param_data.reset();
	}
	// check and recursive adaptation
	int leaf_count = m_nx * m_ny;
	for(int i = 0; i < leaf_count; i++)
		m_grid[i].adaptToParameterization(base_surface, m_grid_vertices);

#ifdef STORE_CONTROL_EPS
	storeEPS("control-adapted-to-parameterization");
#endif
	m_initialized = 2;
}

void ControlSpace2dQuadTree::QuadLeaf::adaptToParameterization(
	SurfaceConstPtr surface, QuadVertexArray& grid_vertices)
{
	if(!isSplit()){
		double g_ratio;
		m_middle_cdm = surface->countParameterizationMatrix(m_middle, g_ratio);
		if(g_ratio < MIN_PARAM_GRATIO){
			m_middle_cdm = m_vertices[VX2D_FIRST]->param_data;
			checkAndSetReparameterization(surface);
			return;
		}
		double max_diff = 0.0;
		for(int i = VX2D_FIRST; i <= VX2D_LAST; i++){
			if(m_vertices[i]->param_data.det() > 0.0){
				double diff = m_middle_cdm.countDifferenceRR(m_vertices[i]->param_data);
				if(diff > max_diff) max_diff = diff;
			}else{
				checkAndSetReparameterization(surface);
				return;
			}
		}
		if(max_diff > param_threshold_diff && 
			m_level < ControlSpace2dQuadTree::param_max_depth_surface_adaptation)
		{// split
			split(grid_vertices, surface, nullptr, true);
		}
	}
	if(isSplit()){
		m_leaves[LF_SW]->adaptToParameterization(surface, grid_vertices);
		m_leaves[LF_SE]->adaptToParameterization(surface, grid_vertices);
		m_leaves[LF_NW]->adaptToParameterization(surface, grid_vertices);
		m_leaves[LF_NE]->adaptToParameterization(surface, grid_vertices);
	}
}

/// Refines control space with respect to reference 3d control space
bool ControlSpace2dQuadTree::adaptToControlSpace3d(CS3dConstPtr space)
{
	// check and recursive adaptation
	int leaf_count = m_nx * m_ny;
	bool any_changes = false;
	for(int i = 0; i < leaf_count; i++)
		any_changes |= m_grid[i].adaptToControlSpace3d(base_surface, m_grid_vertices, space);
	return any_changes;
}

/// Adapt quad leaf with respect to reference 3d control space (split if required)
bool ControlSpace2dQuadTree::QuadLeaf::adaptToControlSpace3d(SurfaceConstPtr surface, 
	QuadVertexArray& grid_vertices, CS3dConstPtr space)
{
	bool any_changes = false;

	if(!isSplit()){
		// calcuate reference metric from 3d space in the middle of this leaf
		const DPoint3d spt3d = surface->getPoint(m_middle);
		if(!space->containsPoint(spt3d)) return false;
		const ControlDataMatrix3d cdm3d = space->getMetricAtPoint(spt3d);
		const ControlDataMatrix2d cdm_ref = surface->projectTransformationTensor(m_middle, cdm3d);

//#ifdef _DEBUG
//		double d3_min = cdm3d.minEigenvalue();
//		double d2_min = cdm_ref.minEigenvalue();
//#endif

		// compare with already existing
		const ControlDataMatrix2d cdm = getMetricAtPoint(m_middle);
		double diff = cdm_ref.countDifferenceRR(cdm);
		if(diff > param_threshold_diff){
			// split may be unnecessary, if cdm_ref is greater than already set for this leaf...
			ControlDataMatrix2d min_cdm = cdm;
			min_cdm.setMinimum(cdm_ref);
			if(cdm.countDifferenceRR(min_cdm) < param_threshold_diff){
				// different from current metric, but introduces no changes - so skip
				return false;
			}
			DMetric2d dmp(cdm_ref, surface, m_middle);
			double len_x = dmp.transformPStoMS(DVector2d(4*m_dl.x, 0.0)).length2();
			double len_y = dmp.transformPStoMS(DVector2d(0.0, 4*m_dl.y)).length2();
			bool split_needed = (len_x > METRIC_LENGTH_RATIO2 || len_y > METRIC_LENGTH_RATIO2 );
			if(split_needed){
				QuadVertexList split_vertices(100);
				split(grid_vertices, surface, &split_vertices, true);
				size_t split_vertices_count = split_vertices.countInt();
				for(size_t i = 0; i < split_vertices_count; i++){
					const DPoint3d local_spt3d = surface->getPoint(split_vertices[i]->coord);
					if(space->containsPoint(local_spt3d)){
						const ControlDataMatrix3d local_cdm3d = space->getMetricAtPoint(local_spt3d);
						ControlDataMatrix2d local_cdm_ref = surface->projectTransformationTensor(split_vertices[i]->coord, local_cdm3d);
//#ifdef _DEBUG
//						double local_d3_min = local_cdm3d.minEigenvalue();
//						double local_d2_min = local_cdm_ref.minEigenvalue();
//#endif
						split_vertices[i]->control_data.setMinimum(local_cdm_ref);
					}else{
						split_vertices[i]->control_data = cdm_ref;
					}
				}
				any_changes = true;
			}
		}
	}
	
	if(isSplit()){
		any_changes |= m_leaves[LF_SW]->adaptToControlSpace3d(surface, grid_vertices, space);
		any_changes |= m_leaves[LF_SE]->adaptToControlSpace3d(surface, grid_vertices, space);
		any_changes |= m_leaves[LF_NW]->adaptToControlSpace3d(surface, grid_vertices, space);
		any_changes |= m_leaves[LF_NE]->adaptToControlSpace3d(surface, grid_vertices, space);
	}

	return any_changes;
}

void ControlSpace2dQuadTree::QuadLeaf::gatherDataFromSourcePoints(
	SurfaceConstPtr surface)
{
	if(isSplit()){
		assert(m_source_points == nullptr);
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			m_leaves[i]->gatherDataFromSourcePoints(surface);
	}else{
		if(m_source_points){
			size_t sct = m_source_points->countInt();
			for(size_t i = 0; i < sct; i++){
				const ControlNode2d& qv = m_source_points->get(i);
				DMetric2d dmp(surface, qv.coord);
				const DPoint2d middle = dmp.transformPStoRS(qv.coord);
				for(int j = VX2D_FIRST; j <= VX2D_LAST; j++){
					double dist = middle.distance2(dmp.transformPStoRS(m_vertices[j]->coord));
					double nw = (dist < mesh_data.relative_small_number)?mesh_data.relative_infinity:(1.0/dist);
					m_vertices[j]->control_data += qv.control_data * nw;
					m_vertices[j]->w += nw;
				}
			}
			delete m_source_points;
			m_source_points = nullptr;
		}
	}
}

/// Propagate control values from within leaves
int ControlSpace2dQuadTree::QuadLeaf::propagateUp()
{
	int ct = 0;
	if(isSplit()){
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			ct += m_leaves[i]->propagateUp();
	}
	for(int i = VX2D_FIRST; i <= VX2D_LAST; i++){
		if(m_vertices[i]->w == 0.0){
			// TODO improve?
			for(int j = 0; j < 2; j++){
				if(m_mid_vertices[nodes_to_mid[i][j]] && 
					m_mid_vertices[nodes_to_mid[i][j]]->w < 0.0){ // i.e. already initialized
					m_vertices[i]->control_data += m_mid_vertices[nodes_to_mid[i][j]]->control_data;
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
int ControlSpace2dQuadTree::QuadLeaf::propagateDown()
{
	int ct = 0;
	for(int i = MVX2D_FIRST; i <= MVX2D_LAST; i++){	// mid-vertices of edges
		if(m_mid_vertices[i] && (m_mid_vertices[i]->w == 0.0)){
			m_mid_vertices[i]->control_data = 
				(m_vertices[mid_to_nodes[i][0]]->control_data + 
				m_vertices[mid_to_nodes[i][1]]->control_data) * 0.5;
			m_mid_vertices[i]->w = -4.0;
			++ct;
		}
	}
	if(isSplit()){
		// middle-vertex
		ControlNode2d* node = m_leaves[LF_NW]->m_vertices[VX2D_SE];
		if(node->w == 0.0){
			for(int i = MVX2D_FIRST; i <= MVX2D_LAST; i++)
				node->control_data += m_mid_vertices[i]->control_data;
			node->control_data *= 0.25;
			node->w = -4.0;
			++ct;
		}
		// and check lower
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			ct += m_leaves[i]->propagateDown();
	}
	return ct;
}

void ControlSpace2dQuadTree::QuadLeaf::statQuadTree(int & max_level, int & leaf_counter, int level) const
{
	if(isSplit()){
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			m_leaves[i]->statQuadTree(max_level, leaf_counter, level+1);
	}else{
		++leaf_counter;
		if(level > max_level) max_level = level;
	}
}

/// Log basic information about this control space
void ControlSpace2dQuadTree::logDescription() const
{
	int max_level = 0;
	int leaf_counter = 0;
	const int main_leaf_count = m_nx * m_ny;
	for(int i = 0; i < main_leaf_count; i++)
		m_grid[i].statQuadTree(max_level, leaf_counter);

	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"QTree [" << m_nx << "," << m_ny << "], nodes= " <<
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
bool ControlSpace2dQuadTree::smoothen()
{
	assert(m_initialized == 2); // with parameterization matrix
	// reset gradation ratio in nodes
	for(size_t i = 0; i < m_grid_vertices.countInt(); i++)
		m_grid_vertices[i].max_gradation_ratio = 1.0;
	// uptree
	int leaf_count = m_nx * m_ny;
	int up_count = 0;
	for(int i = 0; i < leaf_count; i++){
		up_count += m_grid[i].smoothenMetricAtNodes(true);
	}
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "ControlSpace2dQuadTree smoothing (up) - " << up_count << " modifications.");

	// downtree
	int down_count = 0;
	for(int i = 0; i < leaf_count; i++){
		down_count += m_grid[i].smoothenMetricAtNodes(false);
	}
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "ControlSpace2dQuadTree smoothing (down) - " << down_count << " modifications.");
	return up_count+down_count > 0;
}

int ControlSpace2dQuadTree::QuadLeaf::smoothenMetricAtNodes(bool uptree)
{
	int count = 0;
	if(uptree){
		if(isSplit()){	// first leaves then this node
			for(int i = LF_FIRST; i <= LF_LAST; i++)
				count += m_leaves[i]->smoothenMetricAtNodes(uptree);
		}
		if(m_real_dl.x == 0){ // not initialized
			if(isSplit()){
				for(int i = LF_FIRST; i <= LF_LAST; i++)
					m_real_dl += m_leaves[i]->m_real_dl;
				m_real_dl *= 0.5;
			}else{
				assert(m_middle_cdm.det() > 0.0);
				m_real_dl.x = 2*(m_middle_cdm * DVector2d(m_dl.x, 0.0)).length();
				m_real_dl.y = 2*(m_middle_cdm * DVector2d(0.0, m_dl.y)).length();
			}
		}
	}
	const DVector2d dvx(m_real_dl.x, 0.0);
	const DVector2d dvy(0.0, m_real_dl.y);
	// check along leaf-edges
	for(int i = MVX2D_FIRST; i <= MVX2D_LAST; i++){
		int v0 = mid_to_nodes[i][0];
		int v1 = mid_to_nodes[i][1];

		bool along_x = (i == MVX2D_DOWN) || (i == MVX2D_UP);
		int result = ControlSpace2dAdaptive::smoothenMetricForNodes(
			m_vertices[v0], m_vertices[v1], along_x ? dvx : dvy);

		if((result & 1) != 0) ++count;
		if((result & 2) != 0) ++count;
	}

	if(!uptree && isSplit()){	// first this node then leaves
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			count += m_leaves[i]->smoothenMetricAtNodes(uptree);
	}

	return count;
}

/// Calculate inside information for nodes/elements
void ControlSpace2dQuadTree::markInsideNodes(const MeshContainer2d* boundary_mesh)
{
	//MeshViewSet* set = new MeshViewSet;

	START_CLOCK("CSQT::markInsideNodes-1");

	double dx = (m_box.x1 - m_box.x0) / m_nx;
	double dy = (m_box.y1 - m_box.y0) / m_ny;

	const double EPSILON = 0.01;
	double diameter = std::max(dx, dy);
	double z0 = EPSILON * diameter;


	m_grid_vertices.forEach([](ControlNode2d& cn) {
		cn.removeTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN);
	});

	for(IteratorEdge2d it = boundary_mesh->getFirstEdge2d(); it.isValid(); it.nextEdge()){
		const DPoint2d pt0 = m_box.fitInPoint(it.getEdge()->getMeshPoint(0)->getCoordinates());
		const DPoint2d pt1 = m_box.fitInPoint(it.getEdge()->getMeshPoint(1)->getCoordinates());
		// for visualization only
		//set->addEdge( DPoint3d(pt0, z0), DPoint3d(pt1, z0), 1);
		//set->addPoint( DPoint3d(pt0, z0), 1);
		//set->addPoint( DPoint3d(pt1, z0), 1);

		double x0 = (pt0.x - m_box.x0) / dx;
		double x1 = (pt1.x - m_box.x0) / dx;
		double y0 = (pt0.y - m_box.y0) / dy;
		double y1 = (pt1.y - m_box.y0) / dy;
		double xmin = std::min(x0, x1) - EPSILON; // x
		double xmax = std::max(x0, x1) + EPSILON;
		double ymin = std::min(y0, y1) - EPSILON; // y
		double ymax = std::max(y0, y1) + EPSILON;
		DVector2d dpt = pt1 - pt0;

		int x0_grid = (int)(x0+0.5);
		int y0_grid = (int)(y0+0.5);

		int ixmin = std::max((int)xmin, 0);
		int ixmax = std::min((int)xmax, m_nx-1);
		int iymin = std::max((int)ymin, 0);
		int iymax = std::min((int)ymax, m_ny-1);

		if(abs(dpt.x) < z0){
			//---------
			if(abs(x0 - x0_grid) < EPSILON){ // along Y grid line
				for(int iy = iymin; iy <= iymax; iy++){
					int k = iy * m_nx;
					if(x0_grid > 0)
						m_grid[k+x0_grid-1].markInsideNodesY(pt0, pt1, it.getEdge(), true); // x1-line
					if(x0_grid < m_nx)
						m_grid[k+x0_grid].markInsideNodesY(pt0, pt1, it.getEdge(), false); // x0-line
				}
			}else{ // along one Y stripe
				for(int iy = iymin; iy <= iymax; iy++){
					int k = iy * m_nx + ixmin;
					m_grid[k].markInsideNodes(pt0, pt1, it.getEdge());
				}
			}
		}else if(abs(dpt.y) < z0){
			//---------
			if(abs(y0 - y0_grid) < EPSILON){ // along X grid line
				for(int ix = ixmin; ix <= ixmax; ix++){
					if(y0_grid > 0)
						m_grid[(y0_grid-1)*m_nx+ix].markInsideNodesX(pt0, pt1, it.getEdge(), true); // y1-line
					if(y0_grid < m_ny)
						m_grid[y0_grid*m_nx+ix].markInsideNodesX(pt0, pt1, it.getEdge(), false); // y0-line
				}
			}else{ // along one X stripe
				for(int ix = ixmin; ix <= ixmax; ix++){
					int k = iymin * m_nx + ix;
					m_grid[k].markInsideNodes(pt0, pt1, it.getEdge());
				}
			}
		}else{ // several Y strips
			double fdpt = (dpt.y / dpt.x);
			//---------

			//MeshViewSet* set = new MeshViewSet;
			//set->addEdge(DPoint3d(pt0, 0.0), DPoint3d(pt1, 0.0), 1);

			for(int ix = ixmin; ix <= ixmax; ix++){
				x0 = m_box.x0 + ix * dx;
				x1 = x0+dx;

				double pty0 = pt0.y + (x0 - pt0.x) * fdpt;
				double pty1 = pt0.y + (x1 - pt0.x) * fdpt;

				double ptymin = std::min(pt0.y, pt1.y);
				double ptymax = std::max(pt0.y, pt1.y);
				if(pty0 < ptymin) pty0 = ptymin;
				else if (pty0 > ptymax) pty0 = ptymax;
				if(pty1 < ptymin) pty1 = ptymin;
				else if (pty1 > ptymax) pty1 = ptymax;

				double yy0 = (pty0 - m_box.y0) / dy;
				double yy1 = (pty1 - m_box.y0) / dy;
				double yymin = std::min(yy0, yy1) - EPSILON; // y
				double yymax = std::max(yy0, yy1) + EPSILON;

				int iymin = std::max((int)yymin, 0);
				int iymax = std::min((int)yymax, m_ny-1);
				for(int iy = iymin; iy <= iymax; iy++){
					int k = iy * m_nx + ix;
					m_grid[k].markInsideNodes(pt0, pt1, it.getEdge());

					//int pairs[4][2] = { 
					//	{QuadLeaf::VX2D_SW, QuadLeaf::VX2D_NW}, 
					//	{QuadLeaf::VX2D_SE, QuadLeaf::VX2D_NE}, 
					//	{QuadLeaf::VX2D_SW, QuadLeaf::VX2D_SE}, 
					//	{QuadLeaf::VX2D_NW, QuadLeaf::VX2D_NE} };  
					//for(int i = 0; i < 4; i++){
					//	set->addEdge(
					//		DPoint3d(m_grid[k].m_vertices[pairs[i][0]]->coord, 0.0),
					//		DPoint3d(m_grid[k].m_vertices[pairs[i][1]]->coord, 0.0));
					//}
				}
			}
			//SHOW_MESH("scan edge", set);
		}
	}

	STOP_CLOCK("CSQT::markInsideNodes-1");
	START_CLOCK("CSQT::markInsideNodes-2");

	int nxy = m_nx * m_ny;

	// propagate control values from within leaves
	DataVector<QuadLeaf*> not_ready(nxy);
	for(int i = 0; i < nxy; i++){
		if(!m_grid[i].propagateInsideMark())
			not_ready.add(m_grid+i);
	}

//	showMarkInside(boundary_mesh);

	if(!not_ready.empty()){
		// extrapolate to other y1-grid nodes
		propagateMainGridInsideMark();
//		showMarkInside(boundary_mesh);

		size_t last_count = not_ready.countInt();
		while(!not_ready.empty()){
			// propagate control values into leaves			
			for(size_t i = 0; i < not_ready.countInt(); ){
				if(not_ready[i]->propagateInsideMark())
					not_ready.removeAt(i);
				else ++i;
			}

			if(not_ready.countInt() == last_count){
				// TODO: check why?
				while(!not_ready.empty()){
					QuadLeaf *ql = not_ready.removeLast();
					ql->propagateDown();
					ql->propagateUp();
				}
				break;
			}else
				last_count = not_ready.countInt();

//			showMarkInside(boundary_mesh);
		}
	}

	STOP_CLOCK("CSQT::markInsideNodes-2");

	if(MeshGenerator2d::show_prediction){
		START_CLOCK("CSQT::markInsideNodes-3");
		// calculate prediction
		ControlPrediction cp = calculateControlPrediction();

		STOP_CLOCK("CSQT::markInsideNodes-3");

		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, 
			"========== Control Prediction ==============");
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh,
			" * Leaves [0-1-2-3-4]: [" << cp.counter[0] << ","
			<< cp.counter[1] << "," << cp.counter[2] << ",");
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, cp.counter[3] << "," << cp.counter[4] << "]");
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "LEAVES [MIN/AVE/MAX] - " 
			<< cp.value[0] << " / " << cp.value[2] << " / " << cp.value[1]);
		LOG4CPLUS_DEBUG(MeshLog::logger_mesh, "============================================");

		LOG4CPLUS_INFO(MeshLog::logger_console, "Predicted (CS) number of triangles: " << (int)cp.value[2]);
	}
//	SHOW_MESH("control space (markInsideNodes + extend)", set);
}

ControlSpace2dQuadTree::ControlPrediction ControlSpace2dQuadTree::calculateControlPrediction() const
{
	ControlPrediction cp;
	int nxy = m_nx * m_ny;
	for(int i = 0; i < nxy; i++)
		m_grid[i].calculateControlPrediction(cp);

	return cp;
}

/// Calculate mesh prediction using CS
void ControlSpace2dQuadTree::QuadLeaf::calculateControlPrediction(ControlSpace2dQuadTree::ControlPrediction & cp) const
{
	if(isSplit()){
		for(int i = 0; i < 4; i++)
			m_leaves[i]->calculateControlPrediction(cp);
	}else{
		int inside[4] = {
			m_vertices[VX2D_SW]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
			m_vertices[VX2D_SE]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
			m_vertices[VX2D_NW]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
			m_vertices[VX2D_NE]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1) 
		};
		int inside_count = 0;
		int all_nodes[] = {0, 1, 2, 3};
		int inside_nodes[4];
		for(int i = 0; i < 4; i++)
			if(inside[i] > 0)
				inside_nodes[inside_count++] = i;
		cp.counter[inside_count]++;
		if(inside_count == 4) { // for min-approx -> all nodes inside
			calculateControlPrediction(cp, 0, 4, all_nodes);
		}
		if(inside_count > 0){ // for max-approx and ave-approx -> at least one node inside
			calculateControlPrediction(cp, 1, 4, all_nodes);
			calculateControlPrediction(cp, 2, inside_count, inside_nodes);
		}
	}
}

/// Calculate mesh prediction using CS
void ControlSpace2dQuadTree::QuadLeaf::calculateControlPrediction(
	ControlSpace2dQuadTree::ControlPrediction & cp, int ind, int inside_count, int inside[]) const
{
	static const double TRI_AREA = 0.92 * 0.25 * SQRT3;
//	const DVector2d vx(4*m_dl.x, 0.0);
//	const DVector2d vy(0.0, 4*m_dl.y);
	if(inside_count < 1) return;

	ControlDataMatrix2d cdm_ave_arytm = m_vertices[inside[0]]->control_data;
	ControlDataMatrix2d p_cdm = m_vertices[inside[0]]->param_data;
	for(int i = 1; i < inside_count; i++){
		cdm_ave_arytm += m_vertices[inside[i]]->control_data;
		p_cdm += m_vertices[inside[i]]->param_data;
	}
	double w = 1.0 / inside_count;
	cdm_ave_arytm *= w;
	p_cdm *= w;

	double f = 0.25 * inside_count;

	DMetric2d dm_ave_arytm(cdm_ave_arytm, p_cdm);
	//cp.value[ind] += f * (dm_ave_arytm.transformPStoMS(vx).length() * dm_ave_arytm.transformPStoMS(vy).length()) / TRI_AREA;
	cp.value[ind] += f * (dm_ave_arytm.detPStoMS() * (16.0 * m_dl.x * m_dl.y)) / TRI_AREA;
}

int ControlSpace2dQuadTree::QuadLeaf::markInsideNodes(const DPoint2d& pt0, const DPoint2d& pt1, MeshEdge2d* edge)
{
	const double EPSILON = 1e-2;
	const double MIN_DL = EPSILON * std::max(m_dl.x, m_dl.y);

	int res = 0;

	if(isSplit()){
		// recursive
		double ix0 = (pt0.x - m_middle.x) / m_dl.x;
		double ix1 = (pt1.x - m_middle.x) / m_dl.x;
		double iy0 = (pt0.y - m_middle.y) / m_dl.y;
		double iy1 = (pt1.y - m_middle.y) / m_dl.y;

		bool xlow  = 
			(ix0 < EPSILON || ix1 < EPSILON) && 
			(ix0 > -2.0-EPSILON || ix1 > -2.0-EPSILON);
		bool xhigh = 
			(ix0 < 2.0+EPSILON || ix1 < 2.0+EPSILON) && 
			(ix0 > -EPSILON || ix1 > -EPSILON);
		bool ylow  = 
			(iy0 < EPSILON || iy1 <  EPSILON) && 
			(iy0 > -2.0-EPSILON || iy1 >  -2.0-EPSILON);
		bool yhigh = 
			(iy0 < 2.0+EPSILON || iy1 < 2.0+EPSILON) && 
			(iy0 > -EPSILON || iy1 > -EPSILON);

		//MeshViewSet* set = new MeshViewSet;
		//set->addEdge(DPoint3d(pt0, 0.0), DPoint3d(pt1, 0.0), 1);
		//int pairs[4][2] = { {VX2D_SW, VX2D_NW}, {VX2D_SE, VX2D_NE}, {VX2D_SW, VX2D_SE}, {VX2D_NW, VX2D_NE} };  
		//for(int i = 0; i < 4; i++){
		//		set->addEdge(	DPoint3d(m_leaves[VX2D_SW]->m_vertices[pairs[i][0]]->coord, 0.0),
		//						DPoint3d(m_leaves[VX2D_SW]->m_vertices[pairs[i][1]]->coord, 0.0),
		//						(xlow  && ylow) ? 1 : -2);
		//}
		//for(int i = 0; i < 4; i++){
		//		set->addEdge(	DPoint3d(m_leaves[VX2D_NW]->m_vertices[pairs[i][0]]->coord, 0.0),
		//						DPoint3d(m_leaves[VX2D_NW]->m_vertices[pairs[i][1]]->coord, 0.0),
		//						(xlow  && yhigh) ? 1 : -2);
		//}
		//for(int i = 0; i < 4; i++){
		//		set->addEdge(	DPoint3d(m_leaves[VX2D_SE]->m_vertices[pairs[i][0]]->coord, 0.0),
		//						DPoint3d(m_leaves[VX2D_SE]->m_vertices[pairs[i][1]]->coord, 0.0),
		//						(xhigh && ylow) ? 1 : -2);
		//}
		//for(int i = 0; i < 4; i++){
		//		set->addEdge(	DPoint3d(m_leaves[VX2D_NE]->m_vertices[pairs[i][0]]->coord, 0.0),
		//						DPoint3d(m_leaves[VX2D_NE]->m_vertices[pairs[i][1]]->coord, 0.0),
		//						(xhigh && yhigh) ? 1 : -2);
		//}
		//SHOW_MESH("leaf recursion", set);

		int fres[4] = { 0, 0, 0, 0 };
		if(xlow  && ylow)  fres[VX2D_SW] = m_leaves[VX2D_SW]->markInsideNodes(pt0, pt1, edge);
		if(xlow  && yhigh) fres[VX2D_NW] = m_leaves[VX2D_NW]->markInsideNodes(pt0, pt1, edge);
		if(xhigh && ylow)  fres[VX2D_SE] = m_leaves[VX2D_SE]->markInsideNodes(pt0, pt1, edge);
		if(xhigh && yhigh) fres[VX2D_NE] = m_leaves[VX2D_NE]->markInsideNodes(pt0, pt1, edge);

		// propagate-up "forbidden-edges" info
		if(((fres[VX2D_SW] & FORBIDDEN_EDGE_XL) != 0) || ((fres[VX2D_NW] & FORBIDDEN_EDGE_XL) != 0)) {
			m_vertices[VX2D_SW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH);
			m_vertices[VX2D_NW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL);
			res |= FORBIDDEN_EDGE_XL;
		}
		if(((fres[VX2D_SE] & FORBIDDEN_EDGE_XH) != 0) || ((fres[VX2D_NE] & FORBIDDEN_EDGE_XH) != 0)) {
			m_vertices[VX2D_SE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH);
			m_vertices[VX2D_NE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL);
			res |= FORBIDDEN_EDGE_XH;
		}
		if(((fres[VX2D_SW] & FORBIDDEN_EDGE_YL) != 0) || ((fres[VX2D_SE] & FORBIDDEN_EDGE_YL) != 0)) {
			m_vertices[VX2D_SW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH);
			m_vertices[VX2D_SE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL);
			res |= FORBIDDEN_EDGE_YL;
		}
		if(((fres[VX2D_NW] & FORBIDDEN_EDGE_YH) != 0) || ((fres[VX2D_NE] & FORBIDDEN_EDGE_YH) != 0)) {
			m_vertices[VX2D_NW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH);
			m_vertices[VX2D_NE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL);
			res |= FORBIDDEN_EDGE_YH;
		}
	}else{
		// 
		MeshElement* el0 = edge->getMeshElement(0);
		MeshElement* el1 = edge->getMeshElement(1);
		int inside0 = -1;
		if(el0) inside0 = (el0->getAreaID() > -1) ? 1 : 0;
		int inside1 = -1;
		if(el1) inside1 = (el1->getAreaID() > -1) ? 1 : 0;

		double x0 = m_middle.x - 2*m_dl.x;
		double x1 = m_middle.x + 2*m_dl.x;
		double y0 = m_middle.y - 2*m_dl.y;
		double y1 = m_middle.y + 2*m_dl.y;

		const double MIN_DL2 = MIN_DL * MIN_DL;

		if(DSegment2d::distanceToPoint2(pt0, pt1, DPoint2d(x0, y0)) < MIN_DL2){
			m_vertices[VX2D_SW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
			m_vertices[VX2D_SW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 2);
		}
		if(DSegment2d::distanceToPoint2(pt0, pt1, DPoint2d(x1, y0)) < MIN_DL2){
			m_vertices[VX2D_SE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
			m_vertices[VX2D_SE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 2);
		}
		if(DSegment2d::distanceToPoint2(pt0, pt1, DPoint2d(x0, y1)) < MIN_DL2){
			m_vertices[VX2D_NW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
			m_vertices[VX2D_NW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 2);
		}
		if(DSegment2d::distanceToPoint2(pt0, pt1, DPoint2d(x1, y1)) < MIN_DL2){
			m_vertices[VX2D_NE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
			m_vertices[VX2D_NE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 2);
		}

		DVector2d dpt = pt1 - pt0;

		if(abs(dpt.x) > MIN_DL){ // check x0 and x1
			double s = (x0 - pt0.x) / dpt.x; // x0
			if(s > -EPSILON && s < 1.0 + EPSILON){
				double y = pt0.y + s * dpt.y;
				if(y > y0 - MIN_DL && y < y1 + MIN_DL){ 
					if(abs(s) < EPSILON || abs(1.0-s) < EPSILON){ // edge vertex on the x0 control edge
						m_vertices[VX2D_SW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH);
						m_vertices[VX2D_NW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL);
						res |= FORBIDDEN_EDGE_XL;
					}else{
						// if y within [y0,y1] > mark nodes
						double dist = y - y0;
						int v0 = -1, v1 = -1;
						if(dist < MIN_DL){
							m_vertices[VX2D_SW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							m_vertices[VX2D_SW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
						}else if(!m_vertices[VX2D_SW]->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
							(m_vertices[VX2D_SW]->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist))
						{
							m_vertices[VX2D_SW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							v0 = (DTriangle2d::det(pt0, pt1, DPoint2d(x0, y0)) > 0.0) ? inside0 : inside1;
							if(v0 >= 0) m_vertices[VX2D_SW]->setIntTag(TagExtended::TAG_INSIDE_AREA, v0);
						}
						dist = y1 - y;
						if(dist < MIN_DL){
							m_vertices[VX2D_NW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							m_vertices[VX2D_NW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
						}else if(!m_vertices[VX2D_NW]->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
							(m_vertices[VX2D_NW]->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist))
						{
							m_vertices[VX2D_NW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							v1 = (DTriangle2d::det(pt0, pt1, DPoint2d(x0, y1)) > 0.0) ? inside0 : inside1;
							if(v1 >= 0) m_vertices[VX2D_NW]->setIntTag(TagExtended::TAG_INSIDE_AREA, v1);
						}
						if(v0 < 0 || v1 < 0){
							m_vertices[VX2D_SW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH);
							m_vertices[VX2D_NW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL);
							res |= FORBIDDEN_EDGE_XL;
						}
					}
				}
			}
			s = (x1 - pt0.x) / dpt.x; // x1
			if(s > -EPSILON && s < 1.0 + EPSILON){
				double y = pt0.y + s * dpt.y;
				if(y > y0 - MIN_DL && y < y1 + MIN_DL){ 
					if(abs(s) < EPSILON || abs(1.0-s) < EPSILON){ // edge vertex on the x1 control edge
						m_vertices[VX2D_SE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH);
						m_vertices[VX2D_NE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL);
						res |= FORBIDDEN_EDGE_XH;
					}else{
						// if y within [y0,y1] > mark nodes
						double dist = y - y0;
						int v0 = -1, v1 = -1;
						if(dist < MIN_DL){
							m_vertices[VX2D_SE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							m_vertices[VX2D_SE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
						}else if(!m_vertices[VX2D_SE]->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
							(m_vertices[VX2D_SE]->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist))
						{
							m_vertices[VX2D_SE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							v0 = (DTriangle2d::det(pt0, pt1, DPoint2d(x1, y0)) > 0.0) ? inside0 : inside1;
							if(v0 >= 0) m_vertices[VX2D_SE]->setIntTag(TagExtended::TAG_INSIDE_AREA, v0);
						}
						dist = y1 - y;
						if(dist < MIN_DL){
							m_vertices[VX2D_NE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							m_vertices[VX2D_NE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
						}else if(!m_vertices[VX2D_NE]->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
							(m_vertices[VX2D_NE]->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist))
						{
							m_vertices[VX2D_NE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							v1 = (DTriangle2d::det(pt0, pt1, DPoint2d(x1, y1)) > 0.0) ? inside0 : inside1;
							if(v1 >= 0) m_vertices[VX2D_NE]->setIntTag(TagExtended::TAG_INSIDE_AREA, v1);
						}
						if(v0 < 0 || v1 < 0){
							m_vertices[VX2D_SE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH);
							m_vertices[VX2D_NE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL);
							res |= FORBIDDEN_EDGE_XH;
						}
					}
				}
			}
		}

		if(abs(dpt.y) > MIN_DL){ // check y0 and y1
			double s = (y0 - pt0.y) / dpt.y; // y0
			if(s > -EPSILON && s < 1.0 + EPSILON){
				double x = pt0.x + s * dpt.x;
				if(x > x0 - MIN_DL && x < x1 + MIN_DL){ 
					if(abs(s) < EPSILON || abs(1.0-s) < EPSILON){ // edge vertex on the y0 control edge
						m_vertices[VX2D_SW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH);
						m_vertices[VX2D_SE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL);
						res |= FORBIDDEN_EDGE_YL;
					}else{
						// if x within [x0,x1] > mark nodes
						double dist = x - x0;
						int v0 = -1, v1 = -1;
						if(dist < MIN_DL){
							m_vertices[VX2D_SW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							m_vertices[VX2D_SW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
						}else if(!m_vertices[VX2D_SW]->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
							(m_vertices[VX2D_SW]->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist))
						{
							m_vertices[VX2D_SW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							v0 = (DTriangle2d::det(pt0, pt1, DPoint2d(x0, y0)) > 0.0) ? inside0 : inside1;
							if(v0 >= 0) m_vertices[VX2D_SW]->setIntTag(TagExtended::TAG_INSIDE_AREA, v0);
						}
						dist = x1 - x;
						if(dist < MIN_DL){
							m_vertices[VX2D_SE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							m_vertices[VX2D_SE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
						}else if(!m_vertices[VX2D_SE]->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
							(m_vertices[VX2D_SE]->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist))
						{
							m_vertices[VX2D_SE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							v1 = (DTriangle2d::det(pt0, pt1, DPoint2d(x1, y0)) > 0.0) ? inside0 : inside1;
							if(v1 >= 0) m_vertices[VX2D_SE]->setIntTag(TagExtended::TAG_INSIDE_AREA, v1);
						}
						if(v0 < 0 || v1 < 0){
							m_vertices[VX2D_SW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH);
							m_vertices[VX2D_SE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL);
							res |= FORBIDDEN_EDGE_YL;
						}
					}
				}
			}
			s = (y1 - pt0.y) / dpt.y; // y1
			if(s > -EPSILON && s < 1.0 + EPSILON){
				double x = pt0.x + s * dpt.x;
				if(x > x0 - MIN_DL && x < x1 + MIN_DL){ 
					if(abs(s) < EPSILON || abs(1.0-s) < EPSILON){ // edge vertex on the y1 control edge
						m_vertices[VX2D_NW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH);
						m_vertices[VX2D_NE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL);
						res |= FORBIDDEN_EDGE_YH;
					}else{
						// if y within [y0,y1] > mark nodes
						double dist = x - x0;
						int v0 = -1, v1 = -1;
						if(dist < MIN_DL){
							m_vertices[VX2D_NW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							m_vertices[VX2D_NW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
						}else if(!m_vertices[VX2D_NW]->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
							(m_vertices[VX2D_NW]->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist))
						{
							m_vertices[VX2D_NW]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							v0 = (DTriangle2d::det(pt0, pt1, DPoint2d(x0, y1)) > 0.0) ? inside0 : inside1;
							if(v0 >= 0) m_vertices[VX2D_NW]->setIntTag(TagExtended::TAG_INSIDE_AREA, v0);
						}
						dist = x1 - x;
						if(dist < MIN_DL){
							m_vertices[VX2D_NE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							m_vertices[VX2D_NE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 3);
						}else if(!m_vertices[VX2D_NE]->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
							(m_vertices[VX2D_NE]->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist))
						{
							m_vertices[VX2D_NE]->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
							v1 = (DTriangle2d::det(pt0, pt1, DPoint2d(x1, y1)) > 0.0) ? inside0 : inside1;
							if(v1 >= 0) m_vertices[VX2D_NE]->setIntTag(TagExtended::TAG_INSIDE_AREA, v1);
						}
						if(v0 < 0 || v1 < 0){
							m_vertices[VX2D_NW]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH);
							m_vertices[VX2D_NE]->setIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL);
							res |= FORBIDDEN_EDGE_YH;
						}
					}
				}
			}
		}

		//MeshViewSet* set = new MeshViewSet;
		//set->addEdge(DPoint3d(x0, y0, 0.0), DPoint3d(x1, y0, 0.0));
		//set->addEdge(DPoint3d(x1, y0, 0.0), DPoint3d(x1, y1, 0.0));
		//set->addEdge(DPoint3d(x1, y1, 0.0), DPoint3d(x0, y1, 0.0));
		//set->addEdge(DPoint3d(x0, y1, 0.0), DPoint3d(x0, y0, 0.0));
		//set->addEdge(DPoint3d(pt0, 0.0), DPoint3d(pt1, 0.0), 1);
		//int pairs[4][2] = { {VX2D_SW, VX2D_NW}, {VX2D_SE, VX2D_NE}, {VX2D_SW, VX2D_SE}, {VX2D_NW, VX2D_NE} };  
		//for(int i = 0; i < 4; i++){
		//	int v = m_vertices[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1);
		//	if(v >= 0) set->addPoint(DPoint3d(m_vertices[i]->coord, 0.0), v, v);
		//	if(forbidden_edges.contains(ControlEdgeSimple(m_vertices[pairs[i][0]], m_vertices[pairs[i][1]]))) {
		//		set->addEdge(
		//			DPoint3d(m_vertices[pairs[i][0]]->coord, 0.0),
		//			DPoint3d(m_vertices[pairs[i][1]]->coord, 0.0), 2);
		//	}
		//}
		//SHOW_MESH("leaf", set);
	}

	return res;
}

void ControlSpace2dQuadTree::QuadLeaf::markInsideNodesY(const DPoint2d& pt0, const DPoint2d& pt1, MeshEdge2d* edge, bool right_line)
{
	if(isSplit()){
		// recursive
		if(right_line){
			m_leaves[VX2D_SE]->markInsideNodesY(pt0, pt1, edge, right_line);
			m_leaves[VX2D_NE]->markInsideNodesY(pt0, pt1, edge, right_line);
		}else{
			m_leaves[VX2D_SW]->markInsideNodesY(pt0, pt1, edge, right_line);
			m_leaves[VX2D_NW]->markInsideNodesY(pt0, pt1, edge, right_line);
		}
	}else{
		const double EPSILON = 1e-2;
		const double MIN_DL = EPSILON * std::max(m_dl.x, m_dl.y);

		MeshElement* el0 = edge->getMeshElement(0);
		MeshElement* el1 = edge->getMeshElement(1);
		int inside0 = -1;
		if(el0) inside0 = (el0->getAreaID() > -1) ? 1 : 0;
		int inside1 = -1;
		if(el1) inside1 = (el1->getAreaID() > -1) ? 1 : 0;

		double y0 = m_middle.y - 2*m_dl.y;
		double y1 = m_middle.y + 2*m_dl.y;
		double ymin = std::min(pt0.y, pt1.y) - MIN_DL;
		double ymax = std::max(pt0.y, pt1.y) + MIN_DL;

		if(ymin < y0 && ymax > y0){ // y0
			ControlNode2d* node = m_vertices[right_line ? VX2D_SE : VX2D_SW];
			node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
			node->setIntTag(TagExtended::TAG_INSIDE_AREA, 4);
			// other node
			node = m_vertices[right_line ? VX2D_SW : VX2D_SE];
			double dist = 2*m_dl.x;
			if(!node->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
				node->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist)
			{
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
				int v = (DTriangle2d::det(pt0, pt1, node->coord) > 0.0) ? inside0 : inside1;
				if(v >= 0) node->setIntTag(TagExtended::TAG_INSIDE_AREA, v);
			}
		}

		if(ymin < y1 && ymax > y1){ // y1
			ControlNode2d* node = m_vertices[right_line ? VX2D_NE : VX2D_NW];
			node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
			node->setIntTag(TagExtended::TAG_INSIDE_AREA, 4);
			// other node
			node = m_vertices[right_line ? VX2D_NW : VX2D_NE];
			double dist = 2*m_dl.x;
			if(!node->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
				node->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist)
			{
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
				int v = (DTriangle2d::det(pt0, pt1, node->coord) > 0.0) ? inside0 : inside1;
				if(v >= 0) node->setIntTag(TagExtended::TAG_INSIDE_AREA, v);
			}
		}
	}
}

void ControlSpace2dQuadTree::QuadLeaf::markInsideNodesX(const DPoint2d& pt0, const DPoint2d& pt1, MeshEdge2d* edge, bool top_line)
{
	if(isSplit()){
		// recursive
		if(top_line){
			m_leaves[VX2D_NW]->markInsideNodesX(pt0, pt1, edge, top_line);
			m_leaves[VX2D_NE]->markInsideNodesX(pt0, pt1, edge, top_line);
		}else{
			m_leaves[VX2D_SW]->markInsideNodesX(pt0, pt1, edge, top_line);
			m_leaves[VX2D_SE]->markInsideNodesX(pt0, pt1, edge, top_line);
		}
	}else{
		const double EPSILON = 1e-2;
		const double MIN_DL = EPSILON * std::max(m_dl.x, m_dl.y);

		MeshElement* el0 = edge->getMeshElement(0);
		MeshElement* el1 = edge->getMeshElement(1);
		int inside0 = -1;
		if(el0) inside0 = (el0->getAreaID() > -1) ? 1 : 0;
		int inside1 = -1;
		if(el1) inside1 = (el1->getAreaID() > -1) ? 1 : 0;

		double x0 = m_middle.x - 2*m_dl.x;
		double x1 = m_middle.x + 2*m_dl.x;
		double xmin = std::min(pt0.x, pt1.x) - MIN_DL;
		double xmax = std::max(pt0.x, pt1.x) + MIN_DL;

		if(xmin < x0 && xmax > x0){ // x0  
			ControlNode2d* node = m_vertices[top_line ? VX2D_NW : VX2D_SW];
			node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
			node->setIntTag(TagExtended::TAG_INSIDE_AREA, 4);
			// other node
			node = m_vertices[top_line ? VX2D_SW : VX2D_NW];
			double dist = 2*m_dl.y;
			if(!node->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
				node->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist)
			{
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
				int v = (DTriangle2d::det(pt0, pt1, node->coord) > 0.0) ? inside0 : inside1;
				if(v >= 0) node->setIntTag(TagExtended::TAG_INSIDE_AREA, v);
			}
		}

		if(xmin < x1 && xmax > x1){ // x1
			ControlNode2d* node = m_vertices[top_line ? VX2D_NE : VX2D_SE];
			node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, 0.0);
			node->setIntTag(TagExtended::TAG_INSIDE_AREA, 4);
			// other node
			node = m_vertices[top_line ? VX2D_SE : VX2D_NE];
			double dist = 2*m_dl.y;
			if(!node->availableTag(TagExtended::TAG_INSIDE_AREA_DIST) ||
				node->getDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST) > dist)
			{
				node->setDoubleTag(TagExtended::TAG_INSIDE_AREA_DIST, dist);
				int v = (DTriangle2d::det(pt0, pt1, node->coord) > 0.0) ? inside0 : inside1;
				if(v >= 0) node->setIntTag(TagExtended::TAG_INSIDE_AREA, v);
			}
		}
	}
}

double ControlSpace2dQuadTree::getLocalResolution(const DPoint2d& pt) const
{
	QuadLeaf* leaf = findLastLeaf(pt);
	return std::min(leaf->m_dl.x, leaf->m_dl.y);
}

/// Propagate control values from within / into leaves
bool ControlSpace2dQuadTree::QuadLeaf::propagateInsideMark()
{
	if(isSplit()){
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			m_leaves[i]->propagateInsideMark();
	}

	int inside_mark[4] = {
		m_vertices[0]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[1]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[2]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1),
		m_vertices[3]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1) };

	int empty_count = 0;
	for(int i = 0; i < 4; i++)
		if(inside_mark[i] < 0) ++empty_count;
	if(empty_count == 0 || empty_count == 4) return (empty_count == 0);

	int forbidden[4] = {
		m_vertices[0]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[1]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[2]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN),
		m_vertices[3]->getIntTag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN) 
	};

//	static int counter = 10;

	// show degug if: at least one "0" (+ at least one forbidden ...)
	//if(inside_mark[0] == 0 || inside_mark[1] == 0 || inside_mark[2] == 0 || inside_mark[3] == 0){
	//	MeshViewSet* set = new MeshViewSet;
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX2D_SW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX2D_SE]->coord, 0.0), 
	//		(m_vertices[VX2D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) ||
	//		 m_vertices[VX2D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX2D_NW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX2D_NE]->coord, 0.0),  
	//		(m_vertices[VX2D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) ||
	//		 m_vertices[VX2D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX2D_SW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX2D_NW]->coord, 0.0),  
	//		(m_vertices[VX2D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) ||
	//		 m_vertices[VX2D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX2D_SE]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX2D_NE]->coord, 0.0),  
	//		(m_vertices[VX2D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) ||
	//		 m_vertices[VX2D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL)) ? 1 : 0);
	//	for(int i = 0; i < 4; i++){
	//		int v = m_vertices[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1);
	//		if(v >= 0) set->addPoint(DPoint3d(m_vertices[i]->coord, 0.0), v, v);
	//	}

	//	++counter;
	//	ostringstream oss;
	//	oss << "Before, counter = " << counter;
	//	SHOW_MESH(oss.str().c_str(), set);
	//}

	while(empty_count > 0){
		// SW - SE
		if( ((inside_mark[VX2D_SW] < 0) != (inside_mark[VX2D_SE] < 0)) &&
			(inside_mark[VX2D_SW] + inside_mark[VX2D_SE] < 1) &&
			(((forbidden[VX2D_SW] & FORBIDDEN_EDGE_XH) == 0) ||
			((forbidden[VX2D_SE] & FORBIDDEN_EDGE_XL) == 0)))
		{
			if(inside_mark[VX2D_SW] == -1){
				m_vertices[VX2D_SW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
					inside_mark[VX2D_SW] = inside_mark[VX2D_SE]);
//				m_vertices[VX2D_SW]->setIntTag(TagExtended::TAG_ID, ++counter);
			}else{
				m_vertices[VX2D_SE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
					inside_mark[VX2D_SE] = inside_mark[VX2D_SW]);
//				m_vertices[VX2D_SE]->setIntTag(TagExtended::TAG_ID, ++counter);
			}
			--empty_count;
			continue;
		}
		// NW - NE
		if( ((inside_mark[VX2D_NW] < 0) != (inside_mark[VX2D_NE] < 0)) &&
			(inside_mark[VX2D_NW] + inside_mark[VX2D_NE] < 1) &&
			(((forbidden[VX2D_NW] & FORBIDDEN_EDGE_XH) == 0) ||
			((forbidden[VX2D_NE] & FORBIDDEN_EDGE_XL) == 0)))
		{
			if(inside_mark[VX2D_NW] == -1){
				m_vertices[VX2D_NW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
					inside_mark[VX2D_NW] = inside_mark[VX2D_NE]);
//				m_vertices[VX2D_NW]->setIntTag(TagExtended::TAG_ID, ++counter);
			}else{
				m_vertices[VX2D_NE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
					inside_mark[VX2D_NE] = inside_mark[VX2D_NW]);
//				m_vertices[VX2D_NE]->setIntTag(TagExtended::TAG_ID, ++counter);
			}
			--empty_count;
			continue;
		}
		// SW - NW
		if( ((inside_mark[VX2D_SW] < 0) != (inside_mark[VX2D_NW] < 0)) &&
			(inside_mark[VX2D_SW] + inside_mark[VX2D_NW] < 1) &&
			(((forbidden[VX2D_SW] & FORBIDDEN_EDGE_YH) == 0) ||
			((forbidden[VX2D_NW] & FORBIDDEN_EDGE_YL) == 0)))
		{
			if(inside_mark[VX2D_SW] == -1){
				m_vertices[VX2D_SW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
					inside_mark[VX2D_SW] = inside_mark[VX2D_NW]);
//				m_vertices[VX2D_SW]->setIntTag(TagExtended::TAG_ID, ++counter);
			}else{
				m_vertices[VX2D_NW]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
					inside_mark[VX2D_NW] = inside_mark[VX2D_SW]);
//				m_vertices[VX2D_NW]->setIntTag(TagExtended::TAG_ID, ++counter);
			}
			--empty_count;
			continue;
		}
		// SE - NE
		if( ((inside_mark[VX2D_SE] < 0) != (inside_mark[VX2D_NE] < 0)) &&
			(inside_mark[VX2D_SE] + inside_mark[VX2D_NE] < 1) &&
			(((forbidden[VX2D_SE] & FORBIDDEN_EDGE_YH) == 0) ||
			((forbidden[VX2D_NE] & FORBIDDEN_EDGE_YL) == 0)))
		{
			if(inside_mark[VX2D_SE] == -1){
				m_vertices[VX2D_SE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
					inside_mark[VX2D_SE] = inside_mark[VX2D_NE]);
//				m_vertices[VX2D_SE]->setIntTag(TagExtended::TAG_ID, ++counter);
			}else{
				m_vertices[VX2D_NE]->setIntTag(TagExtended::TAG_INSIDE_AREA, 
					inside_mark[VX2D_NE] = inside_mark[VX2D_SE]);
//				m_vertices[VX2D_NE]->setIntTag(TagExtended::TAG_ID, ++counter);
			}
			--empty_count;
			continue;
		}
		// no applicable modification ? skip the loop
		break;
	}

	// show degug if: at least one "0" (+ at least one forbidden ...)
	//if(inside_mark[0] == 0 || inside_mark[1] == 0 || inside_mark[2] == 0 || inside_mark[3] == 0){
	//	MeshViewSet* set = new MeshViewSet;
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX2D_SW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX2D_SE]->coord, 0.0), 
	//		(m_vertices[VX2D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) &&
	//		 m_vertices[VX2D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX2D_NW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX2D_NE]->coord, 0.0),  
	//		(m_vertices[VX2D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) &&
	//		 m_vertices[VX2D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX2D_SW]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX2D_NW]->coord, 0.0),  
	//		(m_vertices[VX2D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) &&
	//		 m_vertices[VX2D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL)) ? 1 : 0);
	//	set->addEdge(
	//		DPoint3d(m_vertices[VX2D_SE]->coord, 0.0), 
	//		DPoint3d(m_vertices[VX2D_NE]->coord, 0.0),  
	//		(m_vertices[VX2D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) &&
	//		 m_vertices[VX2D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL)) ? 1 : 0);
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

/// Finish mark tags from within / into leaves
bool ControlSpace2dQuadTree::QuadLeaf::finishInsideMark()
{
	if(isSplit()){
		for(int i = LF_FIRST; i <= LF_LAST; i++)
			m_leaves[i]->finishInsideMark();
	}

	for(int i = 0; i < 4; i++){
		if(!m_vertices[i]->availableTag(TagExtended::TAG_INSIDE_AREA) ||
			m_vertices[i]->getIntTag(TagExtended::TAG_INSIDE_AREA, -1) < 0)
			m_vertices[i]->setIntTag(TagExtended::TAG_INSIDE_AREA, 1);
	}

	return true;
}

/// Propagate inside mark for the main grid
void ControlSpace2dQuadTree::propagateMainGridInsideMark()
{
	// Start from the initialized nodes
	size_t vx = m_nx+1;
	size_t vy = m_ny+1;
	size_t main_count = vx*vy;
	DataVector<size_t> ready_nodes(main_count);
	for(size_t i = 0; i < main_count; i++){
		if(m_grid_vertices.get(i).availableTag(TagExtended::TAG_INSIDE_AREA))
			ready_nodes.add(i);
	}

	size_t start_ready_nodes = ready_nodes.countInt();
	if(start_ready_nodes == 0) return;
	for(size_t i = 0; i < ready_nodes.countInt() && ready_nodes.countInt() < main_count; i++){
		size_t k = ready_nodes[i];
		ControlNode2d& qvk = m_grid_vertices.get(k);
		int inside_tag = qvk.getIntTag(TagExtended::TAG_INSIDE_AREA);
		if(inside_tag > 1) continue; // only 0 and 1 are used to extend info
		size_t ix = k % vx;
		size_t iy = k / vx;
		// check neighbours of 'k'
		size_t nid[]   = { (k-1),    (k+1),       (k-vx),   (k+vx) };
		bool nbrd[] = { (ix > 0), (ix < vx-1), (iy > 0), (iy < vy-1) };
		int nflags[]  = { FORBIDDEN_EDGE_XL, FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_YL, FORBIDDEN_EDGE_YH };
		int nrflags[] = { FORBIDDEN_EDGE_XH, FORBIDDEN_EDGE_XL, FORBIDDEN_EDGE_YH, FORBIDDEN_EDGE_YL };

		for(int j = 0; j < 4; j++){
			if(nbrd[j]){
				ControlNode2d& qv = m_grid_vertices.get(nid[j]);
				if( !qv.availableTag(TagExtended::TAG_INSIDE_AREA) &&
					(!qv.hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, nflags[j]) ||
					 !qvk.hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, nrflags[j])))
				{
					qv.setIntTag(TagExtended::TAG_INSIDE_AREA, inside_tag);
					ready_nodes.add(nid[j]);
				}
			}
		}
	}
}

/// Show CS-grid with inside info for debug
void ControlSpace2dQuadTree::showMarkInside(const MeshContainer2d* boundary_mesh) const
{
	// start with CS-grid edges
	MeshViewSet* set = getViewSet();

	// border edges

	double dx = (m_box.x1 - m_box.x0) / m_nx;
	double dy = (m_box.y1 - m_box.y0) / m_ny;
	const double EPSILON = 0.01;
	double diameter = std::max(dx, dy);
	double z0 = EPSILON * diameter;
	if(boundary_mesh){
		for(IteratorEdge2d it = boundary_mesh->getFirstEdge2d(); it.isValid(); it.nextEdge()){
			const DPoint2d& pt0 = it.getEdge()->getMeshPoint(0)->getCoordinates();
			const DPoint2d& pt1 = it.getEdge()->getMeshPoint(1)->getCoordinates();
			set->addEdge( DPoint3d(pt0, z0), DPoint3d(pt1, z0), 1);
			//set->addPoint( DPoint3d(pt0, z0), 1);
			//set->addPoint( DPoint3d(pt1, z0), 1);
		}
	}

	// add inside-marks for vertices 

	size_t gv_count = m_grid_vertices.countInt();
	for(size_t i = 0; i < gv_count; i++){
		const ControlNode2d& gv = m_grid_vertices.get(i);
		int v = gv.getIntTag(TagExtended::TAG_INSIDE_AREA, -1);
		int id = gv.getIntTag(TagExtended::TAG_ID, -1);
		if(v >= 0) set->addPoint(DPoint3d(gv.coord, 0.0), v, (id >= 0) ? id : v);
	}

	// + forbidden edges

	DataSimpleList<ControlSpace2dQuadTree::QuadLeaf*> leaves;
	int nxy = m_nx * m_ny;
	for(int i = 0; i < nxy; i++)
		leaves.append(m_grid+i);

	while(!leaves.empty()){
		QuadLeaf* leaf = leaves.removeFirst();

		if(leaf->isSplit()){
			for(int i = 0; i < 4; i++)
				leaves.append(leaf->m_leaves[i]);
		}

		if( leaf->m_vertices[QuadLeaf::VX2D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) &&
			leaf->m_vertices[QuadLeaf::VX2D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL))
		{
			set->addEdge(
				DPoint3d(leaf->m_vertices[QuadLeaf::VX2D_SW]->coord, z0),
				DPoint3d(leaf->m_vertices[QuadLeaf::VX2D_SE]->coord, z0), 2);
		}

		if( leaf->m_vertices[QuadLeaf::VX2D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XH) &&
			leaf->m_vertices[QuadLeaf::VX2D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_XL))
		{
			set->addEdge(
				DPoint3d(leaf->m_vertices[QuadLeaf::VX2D_NW]->coord, z0),
				DPoint3d(leaf->m_vertices[QuadLeaf::VX2D_NE]->coord, z0), 2);
		}

		if( leaf->m_vertices[QuadLeaf::VX2D_SW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) &&
			leaf->m_vertices[QuadLeaf::VX2D_NW]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL))
		{
			set->addEdge(
				DPoint3d(leaf->m_vertices[QuadLeaf::VX2D_SW]->coord, z0),
				DPoint3d(leaf->m_vertices[QuadLeaf::VX2D_NW]->coord, z0), 2);
		}

		if( leaf->m_vertices[QuadLeaf::VX2D_SE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YH) &&
			leaf->m_vertices[QuadLeaf::VX2D_NE]->hasIntFlag(TagExtended::TAG_INSIDE_EDGE_FORBIDDEN, FORBIDDEN_EDGE_YL))
		{
			set->addEdge(
				DPoint3d(leaf->m_vertices[QuadLeaf::VX2D_SE]->coord, z0),
				DPoint3d(leaf->m_vertices[QuadLeaf::VX2D_NE]->coord, z0), 2);
		}
	}

	// show
	SHOW_MESH("CS-grid with inside mark", set);
}
