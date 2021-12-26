#include "ControlSpace3dOctreeLB.h"

int ControlSpace3dKdOctreeLB::createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices)
{
	if (m_top != nullptr) clear();
	KdNeighbors kdn; // empty neighbors
	DPoint3d mid = m_box.getMiddlePoint();
	m_top = std::make_shared<KdOctElementLB>(mid, CDM3d::identity, kdn, m_box, 0);
	new_vertices.add(m_top->getCN());
	return 1;
}

KdOctElementLB::KdOctElementLB(const DPoint3d & mid, const CDM3d & v, const KdNeighbors& kn,
	const DBox& _box, int _level)
	: KdOctElementL(mid, v, kn), leaf_lb(std::make_unique<KdLeafDataLB>(_box, _level)) { }


int KdOctElementLB::splitElementOct(
	DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
	const DBox & box, DBox * boxes,
	DataVector<std::shared_ptr<ControlNode3d>> * kvalues)
{
	if (!isLeaf() || leaf_lb->level >= KdOctree3dParams::paramMaxLevel) return 0;

	int next_level = leaf_lb->level + 1;
	auto split_kn = leaf->kn;

	int result = KdOctElementL::splitElementOct(split_vertices, box, boxes, kvalues);

	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++){
		auto lb = (KdOctElementLB*)(split->elements[vi].get());
		lb->leaf_lb->box = boxes[vi];
		lb->leaf_lb->level = next_level;
	}

	// balance
	DBox sub_boxes[8];
	DPoint3d split_middle = box.getMiddlePoint();
	for (auto kd_e : split_kn.elements) {
		if (kd_e) {
			KdOctElementLB* kd_leaf = (KdOctElementLB*)kd_e->getNearestLeaf(split_middle);
			if (kd_leaf->leaf_lb->level < next_level - KdOctree3dParams::paramMaxBalance)
				result += kd_leaf->splitElementOct(split_vertices, kd_leaf->leaf_lb->box, sub_boxes, kvalues);
		}
	}

	return result;
}

int KdOctElementLB::getTotalBytesExtra() const {
	int total = KdOctElementL::getTotalBytesExtra() + sizeof(KdLeafDataLB*);
	if (leaf != nullptr) total += sizeof(KdLeafDataLB);
	return total;
}

void KdOctElementLB::clear() {
	KdOctElementL::clear();
	leaf_lb.reset();
}

std::shared_ptr<KdElement> KdOctElementLB::clone() const {
	assert(isLeaf());
	return std::make_shared<KdOctElementLB>(leaf->cnode->coord, leaf->cnode->control_data, leaf->kn,
		leaf_lb->box, leaf_lb->level);
}
