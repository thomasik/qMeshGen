// ControlSpace2dMatrixUniform.cpp: implementation of the ControlSpace2dMatrixUniform class.
//
//////////////////////////////////////////////////////////////////////

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "ControlSpace2dMatrixUniform.h"
#include "common.h"
#include "MeshPoint2d.h"
#include "MeshEdge2d.h"
#include "SurfaceParametric.h"
#include "MeshContainer2d.h"
#include "Curve2dSegment.h"
#include "EPSFile.h"
#include "ControlSpace3dAdaptive.h"

int ControlSpace2dMatrixUniform::param_uniform_nx = 100;

ControlSpace2dMatrixUniform::ControlSpace2dMatrixUniform(
	SurfaceConstPtr surface, const DRect& box, int nx, int ny) 
		: ControlSpace2dAdaptive(surface, box)
{
	LOG4CPLUS_DEBUG(MeshLog::logger_console, "Creating uniform control space.");
	if(nx < 4) nx = 4;
	if(ny < 4) ny = 4;
	if(nx > 200) nx = 200;
	if(ny > 200) ny = 200;
	m_nx = nx;
	m_ny = ny;
	m_dx = box.getDX() / (nx - 1);
	m_dy = box.getDY() / (ny - 1);
	m_space = nullptr;

	int size = m_nx * m_ny;
	m_space = new ControlNode2d[size];
	for(int i = 0; i < size; i++){
		int ix = i % m_nx;
		int iy = i / m_nx; 
		m_space[i].w = 0.0;
		m_space[i].coord.x = box.x0 + ix*m_dx;
		m_space[i].coord.y = box.y0 + iy*m_dy;
	}
}

/// Returns number of control nodes in adaptive control structure
int ControlSpace2dMatrixUniform::getControlNodesCount() const
{
	return m_nx * m_ny;
}

/// Invoke function for all control nodes of this space (read-only)
void ControlSpace2dMatrixUniform::forEachControlNode(const std::function<void(const ControlNode2d & node)>& fg) const {
	for (int i = getControlNodesCount() - 1; i >= 0; i--)
		fg(m_space[i]);
}

/// Invoke function for all control nodes of this space
void ControlSpace2dMatrixUniform::forEachControlNode(const std::function<void(ControlNode2d & node)>& fg) {
	for (int i = getControlNodesCount() - 1; i >= 0; i--)
		fg(m_space[i]);
}

ControlSpace2dMatrixUniform::~ControlSpace2dMatrixUniform()
{
	if(m_space) delete[] m_space;
}

/////////////////////////////////////////////////////////////////
// Wyznaczanie pustych wêz³ów siatki kontrolnej wewn¹trz triangulowanych
//	obszarów
bool ControlSpace2dMatrixUniform::interpolate()
{
	// Start from the initialized nodes
	int main_count = m_nx*m_ny;
	DataVector<int> ready_nodes(main_count);
	int *levels = new int[main_count];
	for(int i = 0; i < main_count; i++){
		ControlNode2d& qv = m_space[i];
		if(qv.w > 0.0){
			qv.control_data /= qv.w;
			qv.w = -1.0;
		}
		if(qv.w < 0.0){
			ready_nodes.add(i);
			levels[i] = 0;
		}else
			levels[i] = main_count;
	}
	int new_nodes[4];
	int nearest_nodes[4];
	for(int i = 0; ready_nodes.countInt() < main_count; i++){
		int j, k = ready_nodes[i];
		int cur_level = levels[k];
		int ix = k % m_nx;
		int iy = k / m_nx;
		int ct = 0;
		// check neighbours of 'k'
		if(ix > 0 && levels[k-1] == main_count)		new_nodes[ct++] = k-1;
		if(ix < m_nx-1 && levels[k+1] == main_count)	new_nodes[ct++] = k+1;
		if(iy > 0 && levels[k-m_nx] == main_count)	new_nodes[ct++] = k-m_nx;
		if(iy < m_ny-1 && levels[k+m_nx] == main_count)	new_nodes[ct++] = k+m_nx;
		// if any new, calculate new value and add to ready_nodes
		for(int l = 0; l < ct; l++){
			k = new_nodes[l];
			ix = k % m_nx;
			iy = k / m_nx;
			int nct = 0;
			for(j = 1; ix-j >= 0 && levels[k-j] > cur_level; j++);	// left
			if(ix-j >= 0) nearest_nodes[nct++] = k-j;
			for(j = 1; ix+j < m_nx && levels[k+j] > cur_level; j++);	// right
			if(ix+j < m_nx) nearest_nodes[nct++] = k+j;
			for(j = 1; iy-j >= 0 && levels[k-m_nx*j] > cur_level; j++);	// down
			if(iy-j >= 0) nearest_nodes[nct++] = k-m_nx*j;
			for(j = 1; iy+j < m_ny && levels[k+m_nx*j] > cur_level; j++);	// right
			if(iy+j < m_ny) nearest_nodes[nct++] = k+m_nx*j;
			assert(nct > 0);
			assert(m_space[k].w == 0.0);
			DMetric2d dmp(base_surface, m_space[k].coord);
			const DPoint2d middle = dmp.transformPStoRS(m_space[k].coord);
			for(int m = 0; m < nct; m++){
				ControlNode2d& nqv = m_space[nearest_nodes[m]];
				double dist = middle.distance2(dmp.transformPStoRS(nqv.coord));
				double nw = (dist < mesh_data.relative_small_number)?mesh_data.relative_infinity:(1.0/dist);
				m_space[k].control_data += nqv.control_data * nw;
				m_space[k].w += nw;
			}
			m_space[k].control_data /= m_space[k].w;
			m_space[k].w = -3.0;
			levels[k] = cur_level+1;
			ready_nodes.add(k);
		}
	}
	delete[] levels;
	m_initialized = 1;
	return true;
}

/// Adds new information (stretching and lengths) for some point within the domain
void ControlSpace2dMatrixUniform::addControlNode(const ControlNode2d& node)
{
	// Ustalenie indeksów odpowiedniego podobszaru
	int k = countLocalCoordinates(node.coord);
	DMetric2d dm(base_surface, node.coord);

	// Uaktualnienie danych
	addControlData(dm, k, node);
	addControlData(dm, k+1, node);
	k += m_nx;
	addControlData(dm, k, node);
	addControlData(dm, k+1, node);
}

double ControlSpace2dMatrixUniform::getMetricGradationRatio(const DPoint2d &pt) const
{
	int i,j;
	countLocalCoordinates(pt, i, j);

	double x = m_box.x0 + i * m_dx;
	double y = m_box.y0 + j * m_dy;
	double tx = (pt.x - x) / m_dx;
	double ty = (pt.y - y) / m_dy;
	double coeff[4] = {(1-tx)*(1-ty), tx *(1-ty), (1-tx)* ty, tx * ty};
	int index[4] = {j*m_nx + i, j*m_nx + i + 1, (j+1)*m_nx + i, (j+1)*m_nx + i + 1};
	double ave = m_space[index[0]].max_gradation_ratio * coeff[0];
	for(i = 1; i < 4; i++)
		ave += m_space[index[i]].max_gradation_ratio * coeff[i];

	return ave;
}

ControlDataMatrix2d ControlSpace2dMatrixUniform::getMetricAtPoint(const DPoint2d& pt) const
{
	int i,j;
	countLocalCoordinates(pt, i, j);

	double x = m_box.x0 + i * m_dx;
	double y = m_box.y0 + j * m_dy;
	double tx = (pt.x - x) / m_dx;
	double ty = (pt.y - y) / m_dy;
	double coeff[4] = {(1-tx)*(1-ty), tx *(1-ty), (1-tx)* ty, tx * ty};
	int index[4] = {j*m_nx + i, j*m_nx + i + 1, (j+1)*m_nx + i, (j+1)*m_nx + i + 1};
	ControlDataMatrix2d ave = m_space[index[0]].control_data * coeff[0];
	for(i = 1; i < 4; i++)
		ave += m_space[index[i]].control_data * coeff[i];

	return ave;
}

/// Return interpolated value of extended tag data from control nodes at some given point;
double ControlSpace2dMatrixUniform::interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const
{
	int i,j;
	countLocalCoordinates(pt, i, j);

	double x = m_box.x0 + i * m_dx;
	double y = m_box.y0 + j * m_dy;
	double tx = (pt.x - x) / m_dx;
	double ty = (pt.y - y) / m_dy;
	double coeff[4] = {(1-tx)*(1-ty), tx *(1-ty), (1-tx)* ty, tx * ty};
	int index[4] = {j*m_nx + i, j*m_nx + i + 1, (j+1)*m_nx + i, (j+1)*m_nx + i + 1};
	double ave = m_space[index[0]].getDoubleTag(type) * coeff[0];
	for(i = 1; i < 4; i++)
		ave += m_space[index[i]].getDoubleTag(type) * coeff[i];

	return ave;
}

void ControlSpace2dMatrixUniform::addControlData(
	const DMetric2d& dmp, int k, const ControlNode2d& node)
{
	double dist = dmp.transformPStoRS(node.coord - m_space[k].coord).length2();
	double nw = (dist < mesh_data.relative_small_number)?mesh_data.relative_infinity:(1.0/dist);
	m_space[k].control_data += node.control_data * nw;
	m_space[k].w += nw;
}

bool ControlSpace2dMatrixUniform::setMinControlData(int k, const ControlDataMatrix2d& data)
{
	if(m_space[k].w < 0)
		return m_space[k].control_data.setMinimum(data);
	else
		m_space[k].control_data = data;

	return true;
}

void ControlSpace2dMatrixUniform::setSurfaceCurvatureControlData()
{
	assert(base_surface);

	double p = m_box.getDiameter();
	double max_len = p * std::min(param_max_diameter_ratio, ControlSpace3dAdaptive::param_max_diameter_ratio);
	double min_len = std::max(p * param_min_diameter_ratio, param_min_length);

	for(int j = 0; j < m_ny; j++){
		for(int i = 0; i < m_nx; i++){
			int k = i+m_nx*j;
			const DPoint2d point(m_box.x0 + i*m_dx, m_box.y0 + j*m_dy);
			double g_ratio;
			const SurfaceCurvature curvature = base_surface->getCurvature(point, g_ratio);
			if(g_ratio > MIN_PARAM_GRATIO){
				const ControlDataStretch2d data = adjustCurvatureData(curvature, param_curvature_ratio,
					p, min_len, max_len, param_stretch_max_ratio);
				// Przeliczenie (lx, ly, alfa) -> (dxx, dxy, dyy)
				m_space[k].w = -1.0;
				m_space[k].control_data = DMetric2d::stretchToMatrix(data);
			}else{
				LOG4CPLUS_WARN(MeshLog::logger_console, "Improper parameterization " << g_ratio);
			}
		}
	}

	interpolate();
}

/// Draws the structure of the control space

void ControlSpace2dMatrixUniform::storeEPS(const char* name, int id)
{
	ostringstream fname;
	fname << name << '-' << id << ".eps";

	DRect rect = m_box;
	rect.inflate(0.05);
	EPSFile eps(fname.str(), rect.x0, rect.x1, rect.y0, rect.y1);

	for(int j = 0; j < m_ny; j++){
		for(int i = 0; i < m_nx; i++){
			int k = i+m_nx*j;
			const DPoint2d pt(m_box.x0 + i*m_dx, m_box.y0 + j*m_dy);
			if(m_space[k].w >= 0) continue;
			DVector2d vx = m_space[k].control_data * DVector2d(0.02, 0.0);
			DVector2d vy = m_space[k].control_data * DVector2d(0.0, 0.02);
			eps.drawLine(pt-vx, pt+vx);
			eps.drawLine(pt-vy, pt+vy);
		}
	}
	for(int j = 0; j < m_ny-1; j++){
		for(int i = 0; i < m_nx-1; i++){
			const DPoint2d pt(m_box.x0 + i*m_dx + 0.5*m_dx, m_box.y0+ j*m_dy + 0.5*m_dy);
			ControlDataMatrix2d dm = getMetricAtPoint(pt);
			const DVector2d vx = dm * DVector2d(0.02, 0.0);
			const DVector2d vy = dm * DVector2d(0.0, 0.02);
			eps.drawLine(pt-vx, pt+vx, true);
			eps.drawLine(pt-vy, pt+vy, true);
		}
	}

	eps.drawLine(DPoint2d(-1.0, -1.0), DPoint2d(-1.0,  1.0));
	eps.drawLine(DPoint2d(-1.0,  1.0), DPoint2d( 1.0,  1.0));
	eps.drawLine(DPoint2d( 1.0,  1.0), DPoint2d( 1.0, -1.0));
	eps.drawLine(DPoint2d( 1.0, -1.0), DPoint2d(-1.0, -1.0));

/*
	int gv_count = m_grid_vertices.countInt();
	for(int i = 0; i < gv_count; i++){
		const ControlNode2d& qv = m_grid_vertices.getDataAt(i);
		double gray_color = 1.0;
		if(qv.w < 0) gray_color = -qv.w*0.2-0.2;
		eps.drawPoint(qv.coord, gray_color);
	}
*/
}

bool ControlSpace2dMatrixUniform::setMinControl(const DPoint2d& pt, const ControlDataMatrix2d& cdm, bool min_value_set)
{
	if(!m_box.contains(pt)) return false;
	assert(m_initialized);
	if(!min_value_set) return false; // can't adapt structure anyway
	double diff = cdm.countDifferenceRR(getMetricAtPoint(pt));
	if(diff < ControlSpace2dAdaptive::param_threshold_diff) return false;

	// Ustalenie indeksów odpowiedniego podobszaru
	int k = countLocalCoordinates(pt);

	bool any_changes = false;
	// Uaktualnienie danych
	any_changes |= setMinControlData(k, cdm);
	any_changes |= setMinControlData(k+1, cdm);
	k += m_nx;
	any_changes |= setMinControlData(k, cdm);
	any_changes |= setMinControlData(k+1, cdm);

	return any_changes;
}

/**
* Calculates the int coordinates of left-bottom box (of matrix) containg the given point
*/
void ControlSpace2dMatrixUniform::countLocalCoordinates(const DPoint2d &pt, int &ix, int &iy) const
{
	ix = (int) ((m_nx-1) * (pt.x - m_box.x0) / (m_box.x1 - m_box.x0));
	iy = (int) ((m_ny-1) * (pt.y - m_box.y0) / (m_box.y1 - m_box.y0));
	if(ix > (m_nx - 2)) ix = m_nx - 2;
	if(iy > (m_ny - 2)) iy = m_ny - 2;
	if(ix < 0) ix = 0;
	if(iy < 0) iy = 0;
}

/**
* Calculates the int coordinates of left-bottom box (of matrix) containg the given point
*/
int ControlSpace2dMatrixUniform::countLocalCoordinates(const DPoint2d &pt) const
{
	int ix = (int) ((m_nx-1) * (pt.x - m_box.x0) / (m_box.x1 - m_box.x0));
	int iy = (int) ((m_ny-1) * (pt.y - m_box.y0) / (m_box.y1 - m_box.y0));
	if(ix > (m_nx - 2)) ix = m_nx - 2;
	if(iy > (m_ny - 2)) iy = m_ny - 2;
	if(ix < 0) ix = 0;
	if(iy < 0) iy = 0;
	return iy*m_nx + ix;
}

void ControlSpace2dMatrixUniform::adaptToParameterization()
{
	assert(m_initialized);
	m_initialized = 2;
}

/// Log basic information about this control space
void ControlSpace2dMatrixUniform::logDescription() const
{
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"MatrixUniform [" << m_nx << "," << m_ny << "], nodes =" << getControlNodesCount());
}
