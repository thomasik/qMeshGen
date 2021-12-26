#include <queue>
#include "ControlSpace3dKdTreeLi.h"

ControlDataMatrix3d KdElementLi::getCDM(const DPoint3d & pt, const DBox & box) const 
{
	assert(isLeaf());
	assert(box.contains(pt));

	double tx = (pt.x - box.x0) / box.getDX();
	double ty = (pt.y - box.y0) / box.getDY();
	double tz = (pt.z - box.z0) / box.getDZ();

	double kx = 2 * tx - 1.0; // [-1, 1]
	double ky = 2 * ty - 1.0; // [-1, 1]
	double kz = 2 * tz - 1.0; // [-1, 1]
	double kx2 = kx*kx;
	double ky2 = ky*ky;
	double kz2 = kz*kz;

	double N[28] = { 0.0 }; // shape functions for 27-node hexahedral

	N[1] = 0.125 * (1 - tx)*(1 - ty)*(1 - tz);
	N[2] = 0.125 * (1 + tx)*(1 - ty)*(1 - tz);
	N[3] = 0.125 * (1 + tx)*(1 + ty)*(1 - tz);
	N[4] = 0.125 * (1 - tx)*(1 + ty)*(1 - tz);
	N[5] = 0.125 * (1 - tx)*(1 - ty)*(1 + tz);
	N[6] = 0.125 * (1 + tx)*(1 - ty)*(1 + tz);
	N[7] = 0.125 * (1 + tx)*(1 + ty)*(1 + tz);
	N[8] = 0.125 * (1 - tx)*(1 + ty)*(1 + tz);

	N[9] = 0.25 * (1 - kx2) * (1 - ky)  * (1 - kz);
	N[10] = 0.25 * (1 + kx)  * (1 - ky2) * (1 - kz);
	N[11] = 0.25 * (1 - kx2) * (1 + ky)  * (1 - kz);
	N[12] = 0.25 * (1 - kx)  * (1 - ky2) * (1 - kz);
	N[13] = 0.25 * (1 - kx2) * (1 - ky)  * (1 + kz);
	N[14] = 0.25 * (1 + kx)  * (1 - ky2) * (1 + kz);
	N[15] = 0.25 * (1 - kx2) * (1 + ky)  * (1 + kz);
	N[16] = 0.25 * (1 - kx)  * (1 - ky2) * (1 + kz);
	N[17] = 0.25 * (1 - kx)  * (1 - ky)  * (1 - kz2);
	N[18] = 0.25 * (1 + kx)  * (1 - ky)  * (1 - kz2);
	N[19] = 0.25 * (1 + kx)  * (1 + ky)  * (1 - kz2);
	N[20] = 0.25 * (1 - kx)  * (1 + ky)  * (1 - kz2);

	N[1] -= 0.5 * (N[12] + N[9] + N[17]);
	N[2] -= 0.5 * (N[9] + N[10] + N[18]);
	N[3] -= 0.5 * (N[10] + N[11] + N[19]);
	N[4] -= 0.5 * (N[11] + N[12] + N[20]);
	N[5] -= 0.5 * (N[16] + N[13] + N[17]);
	N[6] -= 0.5 * (N[13] + N[14] + N[18]);
	N[7] -= 0.5 * (N[14] + N[15] + N[19]);
	N[8] -= 0.5 * (N[15] + N[16] + N[20]);

	N[21] = 0.5 * (1 - kx2) * (1 - ky2) * (1 - kz);
	N[22] = 0.5 * (1 - kx2) * (1 - ky2) * (1 + kz);
	N[23] = 0.5 * (1 - kx2) * (1 - ky)  * (1 - kz2);
	N[24] = 0.5 * (1 + kx)  * (1 - ky2) * (1 - kz2);
	N[25] = 0.5 * (1 - kx2) * (1 + ky)  * (1 - kz2);
	N[26] = 0.5 * (1 - kx)  * (1 - ky2) * (1 - kz2);

	N[27] = (1 - kx2)*(1 - ky2)*(1 - kz2);

	N[1] += 0.25 * (N[21] + N[26] + N[23]);
	N[2] += 0.25 * (N[21] + N[23] + N[24]);
	N[3] += 0.25 * (N[21] + N[24] + N[25]);
	N[4] += 0.25 * (N[21] + N[25] + N[26]);
	N[5] += 0.25 * (N[22] + N[26] + N[23]);
	N[6] += 0.25 * (N[22] + N[23] + N[24]);
	N[7] += 0.25 * (N[22] + N[24] + N[25]);
	N[8] += 0.25 * (N[22] + N[25] + N[26]);

	N[9] -= 0.5 * (N[21] + N[23]);
	N[10] -= 0.5 * (N[21] + N[24]);
	N[11] -= 0.5 * (N[21] + N[25]);
	N[12] -= 0.5 * (N[21] + N[26]);
	N[13] -= 0.5 * (N[22] + N[23]);
	N[14] -= 0.5 * (N[22] + N[24]);
	N[15] -= 0.5 * (N[22] + N[25]);
	N[16] -= 0.5 * (N[22] + N[26]);
	N[17] -= 0.5 * (N[26] + N[23]);
	N[18] -= 0.5 * (N[23] + N[24]);
	N[19] -= 0.5 * (N[24] + N[25]);
	N[20] -= 0.5 * (N[25] + N[26]);

	for (int i = 1; i <= 8; i++) N[i] -= 0.125 * N[27];
	for (int i = 9; i <= 20; i++) N[i] += 0.25 * N[27];
	for (int i = 21; i <= 26; i++) N[i] -= 0.5 * N[27];

	const CDM3d& vc = leaf->cnode->control_data;
	const CDM3d& vx0 = leaf->kn[FC3D_WEST] ?
		leaf->kn[FC3D_WEST]->getNearestLeaf(pt)->getCN()->control_data : vc;
	const CDM3d& vx1 = leaf->kn[FC3D_EAST] ?
		leaf->kn[FC3D_EAST]->getNearestLeaf(pt)->getCN()->control_data : vc;
	const CDM3d& vy0 = leaf->kn[FC3D_SOUTH] ?
		leaf->kn[FC3D_SOUTH]->getNearestLeaf(pt)->getCN()->control_data : vc;
	const CDM3d& vy1 = leaf->kn[FC3D_NORTH] ?
		leaf->kn[FC3D_NORTH]->getNearestLeaf(pt)->getCN()->control_data : vc;
	const CDM3d& vz0 = leaf->kn[FC3D_LOW] ?
		leaf->kn[FC3D_LOW]->getNearestLeaf(pt)->getCN()->control_data : vc;
	const CDM3d& vz1 = leaf->kn[FC3D_HIGH] ?
		leaf->kn[FC3D_HIGH]->getNearestLeaf(pt)->getCN()->control_data : vc;

	CDM3d P[28];

	//P[27] = vc;

	P[21] = (vz0 + vc) * 0.5;
	P[22] = (vz1 + vc) * 0.5;
	P[23] = (vy0 + vc) * 0.5;
	P[24] = (vx1 + vc) * 0.5;
	P[25] = (vy1 + vc) * 0.5;
	P[26] = (vx0 + vc) * 0.5;

	P[9] = (P[21] + P[23]) * 0.5;
	P[10] = (P[21] + P[24]) * 0.5;
	P[11] = (P[21] + P[25]) * 0.5;
	P[12] = (P[21] + P[26]) * 0.5;
	P[13] = (P[22] + P[23]) * 0.5;
	P[14] = (P[22] + P[24]) * 0.5;
	P[15] = (P[22] + P[25]) * 0.5;
	P[16] = (P[22] + P[26]) * 0.5;
	P[17] = (P[26] + P[23]) * 0.5;
	P[18] = (P[23] + P[24]) * 0.5;
	P[19] = (P[24] + P[25]) * 0.5;
	P[20] = (P[25] + P[26]) * 0.5;

	P[1] = (P[21] + P[26] + P[23]) * (1.0 / 3.0);
	P[2] = (P[21] + P[23] + P[24]) * (1.0 / 3.0);
	P[3] = (P[21] + P[24] + P[25]) * (1.0 / 3.0);
	P[4] = (P[21] + P[26] + P[23]) * (1.0 / 3.0);
	P[5] = (P[22] + P[26] + P[23]) * (1.0 / 3.0);
	P[6] = (P[22] + P[23] + P[24]) * (1.0 / 3.0);
	P[7] = (P[22] + P[24] + P[25]) * (1.0 / 3.0);
	P[8] = (P[22] + P[25] + P[26]) * (1.0 / 3.0);


	CDM3d result = vc * N[27];
	for (int i = 1; i <= 26; i++)
		result += P[i] * N[i];
	return result;
}

ControlDataMatrix3d KdElementLi::getCDM(const DPoint3d & pt) const {
	assert(false);
	// has to be with box...
	return CDM3d::identity;
}

std::shared_ptr<KdElement> KdElementLi::clone() const {
	assert(isLeaf());
	return std::make_shared<KdElementLi>(leaf->cnode->coord, leaf->cnode->control_data, leaf->kn);
}

int ControlSpace3dKdTreeLi::adaptToField(const std::function<CDM3d(const DPoint3d&pt)>& f, double max_error) 
{
	struct KdDataExtra {
		KdDataExtra(std::shared_ptr<KdElement> _leaf, const DBox& _box, int _level) 
			: leaf((KdElementLi*)(_leaf.get())), box(_box), level(_level) {}
		KdElementLi* leaf;
		DBox box;
		int level;
	};

	DataVector<std::shared_ptr<ControlNode3d>> vertices(10);
	int kd_count = createTopElement(vertices);
	vertices.forEach([&f](std::shared_ptr<ControlNode3d> cn) { cn->control_data = f(cn->coord); });

	std::queue<KdDataExtra> todo_queue;
	KdDataExtra task(m_top, m_box, 0);
	todo_queue.push(task);

	Axis split_axis;
	double split_value;
	DBox box0, box1;

	while (!todo_queue.empty()) {
		task = todo_queue.front();
		todo_queue.pop();

		assert(task.leaf->isLeaf());
		double curr_error = task.leaf->estimateError(task.box, f);
		if (curr_error <= max_error) continue;

		KdTree3dParams::m_cf_splitter->whereToSplit(task.box, f, split_axis, split_value);

		DataVector<std::shared_ptr<ControlNode3d>> split_vertices(10);
		kd_count += task.leaf->splitElementKd(split_vertices, task.box, split_axis, split_value, &box0, &box1);

		split_vertices.forEach([&f](std::shared_ptr<ControlNode3d> node) {
			node->control_data = f(node->coord);
		});

		int next_level = task.level + 1;
		if (next_level < KdTree3dParams::paramMaxLevel) {
			todo_queue.push(KdDataExtra(task.leaf->split->elements[0], box0, next_level));
			todo_queue.push(KdDataExtra(task.leaf->split->elements[1], box1, next_level));
		}

	}
	m_initialized = 1;
	return kd_count;
}

int ControlSpace3dKdTreeLi::createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices)
{
	if (m_top != nullptr) clear();
	KdNeighbors kdn; // empty neighbors
	DPoint3d mid = m_box.getMiddlePoint();
	m_top = std::make_shared<KdElementLi>(mid, CDM3d::identity, kdn);
	new_vertices.add(m_top->getCN());
	return 1;
}
