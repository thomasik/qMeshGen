// ControlSpace3dMatrixUniform.cpp: implementation of the ControlSpace3dMatrixUniform class.
//
//////////////////////////////////////////////////////////////////////

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#include "ControlSpace3dMatrixUniform.h"
#include "ControlSpace2dAdaptive.h"
#include "DataVector.h"
#include "common.h"

int ControlSpace3dMatrixUniform::param_uniform_nx = 30;

ControlSpace3dMatrixUniform::ControlSpace3dMatrixUniform(
	const DBox& box, int nx, int ny, int nz) 
		: ControlSpace3dAdaptive(box)
{
	LOG4CPLUS_INFO(MeshLog::logger_console, "Creating uniform control space 3d.");
	if(nx < 4) nx = 4;
	if(ny < 4) ny = 4;
	if(nz < 4) nz = 4;
	if(nx > 100) nx = 100;
	if(ny > 100) ny = 100;
	if(nz > 100) nz = 100;
	m_box = box;
	m_nx = nx;
	m_ny = ny;
	m_nz = nz;
	m_dx = box.getDX() / (nx - 1);
	m_dy = box.getDY() / (ny - 1);
	m_dz = box.getDZ() / (nz - 1);
	m_space = nullptr;

	int size = m_nx * m_ny * m_nz;
	m_space = new ControlNode3d[size];
	for(int i = 0; i < size; i++){
		int ix, iy, iz;
		countLocalCoordinates(i, ix, iy, iz);
		m_space[i].w = 0.0;
		m_space[i].coord.x = box.x0 + ix*m_dx;
		m_space[i].coord.y = box.y0 + iy*m_dy;
		m_space[i].coord.z = box.z0 + iz*m_dz;
	}
}

/// Returns number of control nodes in adaptive control structure
int ControlSpace3dMatrixUniform::getControlNodesCount() const
{
	return m_nx * m_ny * m_nz;
}

/// Invoke function for all control nodes of this space (read-only)
void ControlSpace3dMatrixUniform::forEachControlNode(const std::function<void(const ControlNode3d&node)>& fg) const {
	for (int i = getControlNodesCount() - 1; i >= 0; i--)
		fg(m_space[i]);
}

/// Invoke function for all control nodes of this space
void ControlSpace3dMatrixUniform::forEachControlNode(const std::function<void(ControlNode3d&node)>& fg) {
	for (int i = getControlNodesCount() - 1; i >= 0; i--)
		fg(m_space[i]);
}

ControlSpace3dMatrixUniform::~ControlSpace3dMatrixUniform()
{
	if(m_space) delete[] m_space;
}

/////////////////////////////////////////////////////////////////
// Wyznaczanie pustych wêz³ów siatki kontrolnej wewn¹trz triangulowanych
//	obszarów
bool ControlSpace3dMatrixUniform::interpolate()
{
	// Start from the initialized nodes
	int main_count = m_nx*m_ny*m_nz;
	DataVector<int> ready_nodes(main_count);
	int *levels = new int[main_count];
	for(int i = 0; i < main_count; i++){
		ControlNode3d& qv = m_space[i];
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
	int new_nodes[6];
	int nearest_nodes[6];
	int nxy = m_nx*m_ny;
	for(int i = 0; ready_nodes.countInt() < main_count; i++){
		int j, k = ready_nodes[i];
		int cur_level = levels[k];
		int ix, iy, iz;
		countLocalCoordinates(k, ix, iy, iz);
		int ct = 0;
		// check neighbours of 'k'
		if(ix > 0 && levels[k-1] == main_count)		new_nodes[ct++] = k-1;
		if(ix < m_nx-1 && levels[k+1] == main_count)	new_nodes[ct++] = k+1;
		if(iy > 0 && levels[k-m_nx] == main_count)	new_nodes[ct++] = k-m_nx;
		if(iy < m_ny-1 && levels[k+m_nx] == main_count)	new_nodes[ct++] = k+m_nx;
		if(iz > 0 && levels[k-nxy] == main_count)	new_nodes[ct++] = k-nxy;
		if(iz < m_nz-1 && levels[k+nxy] == main_count)	new_nodes[ct++] = k+nxy;
		// if any new, calculate new value and add to ready_nodes
		for(int l = 0; l < ct; l++){
			k = new_nodes[l];
			countLocalCoordinates(k, ix, iy, iz);
			int nct = 0;
			for(j = 1; ix-j >= 0 && levels[k-j] > cur_level; j++);	// left
			if(ix-j >= 0) nearest_nodes[nct++] = k-j;
			for(j = 1; ix+j < m_nx && levels[k+j] > cur_level; j++);	// right
			if(ix+j < m_nx) nearest_nodes[nct++] = k+j;
			for(j = 1; iy-j >= 0 && levels[k-m_nx*j] > cur_level; j++);	// down
			if(iy-j >= 0) nearest_nodes[nct++] = k-m_nx*j;
			for(j = 1; iy+j < m_ny && levels[k+m_nx*j] > cur_level; j++);	// right
			if(iy+j < m_ny) nearest_nodes[nct++] = k+m_nx*j;
			for(j = 1; iz-j >= 0 && levels[k-nxy*j] > cur_level; j++);	// near
			if(iz-j >= 0) nearest_nodes[nct++] = k-nxy*j;
			for(j = 1; iz+j < m_nz && levels[k+nxy*j] > cur_level; j++);	// far
			if(iz+j < m_nz) nearest_nodes[nct++] = k+nxy*j;
			assert(nct > 0);
			assert(m_space[k].w == 0.0);
			for(int m = 0; m < nct; m++){
				ControlNode3d& nqv = m_space[nearest_nodes[m]];
				double dist = m_space[k].coord.distance2(nqv.coord);
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
void ControlSpace3dMatrixUniform::addControlNode(const ControlNode3d& node)
{
	// Ustalenie indeksów odpowiedniego podobszaru
	int k = countLocalCoordinates(node.coord);

	// Uaktualnienie danych
	addControlData(k, node);
	addControlData(k+1, node);
	k += m_nx;
	addControlData(k, node);
	addControlData(k+1, node);
	k += m_nx*m_ny;
	addControlData(k, node);
	addControlData(k+1, node);
	k -= m_nx;
	addControlData(k, node);
	addControlData(k+1, node);
}

ControlDataMatrix3d ControlSpace3dMatrixUniform::getMetricAtPoint(const DPoint3d& pt) const
{
	int ix,iy,iz;
	countLocalCoordinates(pt, ix, iy, iz);

	double x = m_box.x0 + ix * m_dx;
	double y = m_box.y0 + iy * m_dy;
	double z = m_box.z0 + iz * m_dz;
	double tx = (pt.x - x) / m_dx;
	double ty = (pt.y - y) / m_dy;
	double tz = (pt.z - z) / m_dz;
	int m_nxy = m_nx*m_ny;

	double coeff[8] = {(1-tx)*(1-ty)*(1-tz), tx*(1-ty)*(1-tz), (1-tx)*ty*(1-tz), tx*ty*(1-tz),
						(1-tx)*(1-ty)*tz, tx*(1-ty)*tz, (1-tx)*ty*tz, tx*ty*tz};
	int index[8] = {iz*m_nxy + iy*m_nx + ix, iz*m_nxy + iy*m_nx + ix + 1, 
					iz*m_nxy + (iy+1)*m_nx + ix, iz*m_nxy + (iy+1)*m_nx + ix + 1,
					(iz+1)*m_nxy + iy*m_nx + ix, (iz+1)*m_nxy + iy*m_nx + ix + 1, 
					(iz+1)*m_nxy + (iy+1)*m_nx + ix, (iz+1)*m_nxy + (iy+1)*m_nx + ix + 1};
	ControlDataMatrix3d ave = m_space[index[0]].control_data * coeff[0];
	for(int i = 1; i < 8; i++)
		ave += m_space[index[i]].control_data * coeff[i];

	return ave;
}

/// Return interpolated value of extended tag data from control nodes at some given point;
double ControlSpace3dMatrixUniform::interpolateDoubleTag(const DPoint3d& pt, TagExtended::TagType type) const
{
	int ix,iy,iz;
	countLocalCoordinates(pt, ix, iy, iz);

	double x = m_box.x0 + ix * m_dx;
	double y = m_box.y0 + iy * m_dy;
	double z = m_box.z0 + iz * m_dz;
	double tx = (pt.x - x) / m_dx;
	double ty = (pt.y - y) / m_dy;
	double tz = (pt.z - z) / m_dz;
	int m_nxy = m_nx*m_ny;

	double coeff[8] = {(1-tx)*(1-ty)*(1-tz), tx*(1-ty)*(1-tz), (1-tx)*ty*(1-tz), tx*ty*(1-tz),
						(1-tx)*(1-ty)*tz, tx*(1-ty)*tz, (1-tx)*ty*tz, tx*ty*tz};
	int index[8] = {iz*m_nxy + iy*m_nx + ix, iz*m_nxy + iy*m_nx + ix + 1, 
					iz*m_nxy + (iy+1)*m_nx + ix, iz*m_nxy + (iy+1)*m_nx + ix + 1,
					(iz+1)*m_nxy + iy*m_nx + ix, (iz+1)*m_nxy + iy*m_nx + ix + 1, 
					(iz+1)*m_nxy + (iy+1)*m_nx + ix, (iz+1)*m_nxy + (iy+1)*m_nx + ix + 1};
	double ave = m_space[index[0]].getDoubleTag(type) * coeff[0];
	for(int i = 1; i < 8; i++)
		ave += m_space[index[i]].getDoubleTag(type) * coeff[i];

	return ave;
}

void ControlSpace3dMatrixUniform::addControlData(int k, const ControlNode3d& node)
{
	assert(k >= 0 && k < m_nx*m_ny*m_nz);
	double dist = node.coord.distance2(m_space[k].coord);
	double nw = (dist < mesh_data.relative_small_number)?mesh_data.relative_infinity:(1.0/dist);
	m_space[k].control_data += node.control_data * nw;
	m_space[k].w += nw;
}

bool ControlSpace3dMatrixUniform::setMinControlData(int k, const ControlDataMatrix3d& data)
{
	assert(k >= 0 && k < m_nx*m_ny*m_nz);
	if(m_space[k].w < 0)
		return m_space[k].control_data.setMinimum(data);
	else
		m_space[k].control_data = data;

	return true;
}

bool ControlSpace3dMatrixUniform::setMinControl(const DPoint3d& pt, 
	const ControlDataMatrix3d& cdm, bool min_value_set)
{
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
	k += m_nx*m_ny;
	any_changes |= setMinControlData(k, cdm);
	any_changes |= setMinControlData(k+1, cdm);
	k -= m_nx;
	any_changes |= setMinControlData(k, cdm);
	any_changes |= setMinControlData(k+1, cdm);

	return any_changes;
}

void ControlSpace3dMatrixUniform::countLocalCoordinates(int k, int &ix, int &iy, int &iz) const
{
	ix = k % m_nx;
	int iyz = k / m_nx; 
	iy = iyz % m_ny;
	iz = iyz / m_ny; 
}

/**
* Calculates the int coordinates of left-bottom box (of matrix) containg the given point
*/
void ControlSpace3dMatrixUniform::countLocalCoordinates(const DPoint3d &pt, int &ix, int &iy, int &iz) const
{
	ix = (int) ((m_nx-1) * (pt.x - m_box.x0) / m_box.getDX());
	iy = (int) ((m_ny-1) * (pt.y - m_box.y0) / m_box.getDY());
	iz = (int) ((m_nz-1) * (pt.z - m_box.z0) / m_box.getDZ());
	if(ix > (m_nx - 2)) ix = m_nx - 2;
	if(iy > (m_ny - 2)) iy = m_ny - 2;
	if(iz > (m_nz - 2)) iz = m_nz - 2;
	if(ix < 0) ix = 0;
	if(iy < 0) iy = 0;
	if(iz < 0) iz = 0;
}

/**
* Calculates the int coordinates of left-bottom box (of matrix) containg the given point
*/
int ControlSpace3dMatrixUniform::countLocalCoordinates(const DPoint3d &pt) const
{
	int ix = (int) ((m_nx-1) * (pt.x - m_box.x0) / m_box.getDX());
	int iy = (int) ((m_ny-1) * (pt.y - m_box.y0) / m_box.getDY());
	int iz = (int) ((m_nz-1) * (pt.z - m_box.z0) / m_box.getDZ());
	if(ix > (m_nx - 2)) ix = m_nx - 2;
	if(iy > (m_ny - 2)) iy = m_ny - 2;
	if(iz > (m_nz - 2)) iz = m_nz - 2;
	if(ix < 0) ix = 0;
	if(iy < 0) iy = 0;
	if(iz < 0) iz = 0;
	return (iz*m_ny + iy)*m_nx + ix;
}

/// Log basic information about this control space
void ControlSpace3dMatrixUniform::logDescription() const
{
	LOG4CPLUS_INFO(MeshLog::logger_mesh, 
		"MatrixUniform3d [" << m_nx << "," << m_ny << "," << m_nz << "], nodes =" << getControlNodesCount());
}

bool ControlSpace3dMatrixUniform::smoothen()
{
	assert(m_initialized == 1);
	
	int ct = m_nx * m_ny * m_nz;

	bool *active = new bool[ct];
	struct NB {
		int x1,x2,y1,y2,z1,z2;
	}	*neighbors = new NB[ct];
	int ix,iy,iz;
	int nxy = m_nx*m_ny;

	for(int k = 0; k < ct; k++){
		active[k] = true;
		countLocalCoordinates(k, ix, iy, iz);
		neighbors[k].x1 = (ix > 0)		? (k-1) : -1;		// x
		neighbors[k].x2 = (ix < m_nx-1) ? (k+1) : -1;
		neighbors[k].y1 = (iy > 0)		? (k-m_nx) : -1;	// y
		neighbors[k].y2 = (iy < m_ny-1) ? (k+m_nx) : -1;
		neighbors[k].z1 = (iz > 0)		? (k-nxy)  : -1;	// z
		neighbors[k].z2 = (iz < m_nz-1) ? (k+nxy)  : -1;
	}

	const DVector3d dvx(m_box.getDX() / (m_nx-1), 0.0, 0.0);
	const DVector3d dvy(0.0, m_box.getDY() / (m_ny-1), 0.0);
	const DVector3d dvz(0.0, 0.0, m_box.getDZ() / (m_nz-1));

	bool any_active = true;
	bool any_change = false;

	int max_steps = m_nx + m_ny + m_nz;
	int steps_counter = 0;

	while(any_active){
		any_active = false;
		if(++steps_counter > max_steps) break;
		for(int k = 0; k < ct; k++)
			if(active[k]){
				any_active = true;
				active[k] = false;
				int k2;
				if((k2=neighbors[k].x1) >= 0){ 
					assert(k2 >= 0 && k2 < ct);
					int r = smoothenMetricForNodes(&m_space[k], &m_space[k2], dvx);
					active[k]  |= ((r&1)>0);
					active[k2] |= ((r&2)>0);
				}
				if((k2=neighbors[k].x2) >= 0){
					assert(k2 >= 0 && k2 < ct);
					int r = smoothenMetricForNodes(&m_space[k], &m_space[k2], dvx);
					active[k]  |= ((r&1)>0);
					active[k2] |= ((r&2)>0);
				}
				//
				if((k2=neighbors[k].y1) >= 0){
					assert(k2 >= 0 && k2 < ct);
					int r = smoothenMetricForNodes(&m_space[k], &m_space[k2], dvy);
					active[k]  |= ((r&1)>0);
					active[k2] |= ((r&2)>0);
				}
				if((k2=neighbors[k].y2) >= 0){
					assert(k2 >= 0 && k2 < ct);
					int r = smoothenMetricForNodes(&m_space[k], &m_space[k2], dvy);
					active[k]  |= ((r&1)>0);
					active[k2] |= ((r&2)>0);
				}
				//
				if((k2=neighbors[k].z1) >= 0){
					assert(k2 >= 0 && k2 < ct);
					int r = smoothenMetricForNodes(&m_space[k], &m_space[k2], dvz);
					active[k]  |= ((r&1)>0);
					active[k2] |= ((r&2)>0);
				}
				if((k2=neighbors[k].z2) >= 0){
					assert(k2 >= 0 && k2 < ct);
					int r = smoothenMetricForNodes(&m_space[k], &m_space[k2], dvz);
					active[k]  |= ((r&1)>0);
					active[k2] |= ((r&2)>0);
				}
			}
		any_change |= any_active;
	}

	delete[] active;
	delete[] neighbors;

	return any_change;
}

double ControlSpace3dMatrixUniform::getMetricGradationRatio(const DPoint3d &pt) const
{
	int ix,iy,iz;
	countLocalCoordinates(pt, ix, iy, iz);

	double x = m_box.x0 + ix * m_dx;
	double y = m_box.y0 + iy * m_dy;
	double z = m_box.z0 + iz * m_dz;
	double tx = (pt.x - x) / m_dx;
	double ty = (pt.y - y) / m_dy;
	double tz = (pt.z - z) / m_dz;
	int m_nxy = m_nx*m_ny;

	double coeff[8] = {(1-tx)*(1-ty)*(1-tz), tx*(1-ty)*(1-tz), (1-tx)*ty*(1-tz), tx*ty*(1-tz),
						(1-tx)*(1-ty)*tz, tx*(1-ty)*tz, (1-tx)*ty*tz, tx*ty*tz};
	int index[8] = {iz*m_nxy + iy*m_nx + ix, iz*m_nxy + iy*m_nx + ix + 1, 
					iz*m_nxy + (iy+1)*m_nx + ix, iz*m_nxy + (iy+1)*m_nx + ix + 1,
					(iz+1)*m_nxy + iy*m_nx + ix, (iz+1)*m_nxy + iy*m_nx + ix + 1, 
					(iz+1)*m_nxy + (iy+1)*m_nx + ix, (iz+1)*m_nxy + (iy+1)*m_nx + ix + 1};
	assert(!m_space[index[0]].gradationUnknown());
	double ave = m_space[index[0]].max_gradation_ratio * coeff[0];
	for (int i = 1; i < 8; i++) {
		assert(!m_space[index[i]].gradationUnknown());
		ave += m_space[index[i]].max_gradation_ratio * coeff[i];
	}

	return ave;
}

