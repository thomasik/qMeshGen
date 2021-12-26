#include "ControlSpace3dOctreeL.h"

int ControlSpace3dKdOctreeL::createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices)
{
	if (m_top != nullptr) clear();
	KdNeighbors kdn; // empty neighbors
	DPoint3d mid = m_box.getMiddlePoint();
	m_top = std::make_shared<KdOctElementL>(mid, CDM3d::identity, kdn);
	new_vertices.add(m_top->getCN());
	return 1;
}

int KdOctElementL::adaptToField(
	const std::function<CDM3d(const DPoint3d & pt)> & f,
	const DBox & box, double max_error, int level,
	DataVector<std::shared_ptr<ControlNode3d>> * /* kvalues */)
{
	assert(isLeaf());
	if (level >= KdOctree3dParams::paramMaxLevel) return 0;
	double curr_error = estimateError(box, f);
	if (curr_error <= max_error) return 0;

	DataVector<std::shared_ptr<ControlNode3d>> split_vertices(10);
	DBox boxes[8];
	int result = splitElementOct(split_vertices, box, boxes);

	split_vertices.forEach([&f](std::shared_ptr<ControlNode3d> node) {
		node->control_data = f(node->coord);
	});

	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
		result += split->elements[vi]->adaptToField(f, boxes[vi], max_error, level + 1);
	}
	return result;
}

int KdOctElementL::splitElementOct(
	DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
	const DBox & box, DBox * boxes,
	DataVector<std::shared_ptr<ControlNode3d>> * /* kvalues */)
{
	assert(isLeaf());
	assert(leaf != nullptr && split == nullptr);

	KdNeighbors split_kn = leaf->kn;
	std::shared_ptr<KdElement> leaves[VX3D_COUNT];
	for(int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++)
		leaves[vi] = this->clone();
	clear();

	DPoint3d mid = box.getMiddlePoint();
	split = new KdOctreeSplitElement(mid);
	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
		boxes[vi] = box.splitOct(mid, vi);
		auto oct_leaf = (KdOctElementL*)(leaves[vi].get());
		oct_leaf->leaf->cnode->coord = boxes[vi].getMiddlePoint();
		split_vertices.add(leaves[vi]->getCN());
		split->elements.add(leaves[vi]);
	}

	// update neighbors (outer)
	for (int fi = FC3D_FIRST; fi <= FC3D_LAST; fi++) {
		auto nb = split_kn[fi];
		if (nb != nullptr && !nb->isLeaf()) {
			int ofi = fi ^ 1; // opposite face (xor 1)
			for (int i = 0; i < 4; i++) {
				int vfi = DBox::face_to_vertex[fi][i];
				int ovfi = DBox::face_to_vertex[ofi][i];
				auto sub_okd = nb->split->elements[ovfi];
				auto sub_kd = split->elements[vfi];
				sub_kd->setKN(fi, sub_okd);
				if (sub_okd->isLeaf()) sub_okd->setKN(ofi, sub_kd);
			}
		}
		if (fi % 2 == 0) {
			// inner neighborhood
			int ofi = fi ^ 1; // opposite face (xor 1), also fi+1
			for (int i = 0; i < 4; i++) {
				int vfi = DBox::face_to_vertex[fi][i];
				int ovfi = DBox::face_to_vertex[ofi][i];
				auto sub_okd = split->elements[ovfi];
				auto sub_kd = split->elements[vfi];
				sub_kd->setKN(ofi, sub_okd);
				sub_okd->setKN(fi, sub_kd);
			}
		}
	}

	return VX3D_COUNT;
}

void KdOctElementL::setKN(const KdNeighbors & kn) {
	assert(leaf != nullptr);
	leaf->kn = kn;
}

void KdOctElementL::setKN(int fi, std::shared_ptr<KdElement> const & kde) {
	assert(leaf != nullptr);
	leaf->kn[fi] = kde;
}

int KdOctElementL::smoothenMetricAtNodes() {
	int count = 0;
	if (isLeaf()) {
		// check with neighbors
		auto cn0 = getCN();
		for (auto nb : leaf->kn.elements) {
			if (nb == nullptr) continue;
			auto cn1 = nb->getNearestLeaf(cn0->coord)->getCN();

			int result = ControlSpace3dAdaptive::smoothenMetricForNodes(
				cn0.get(), cn1.get(), cn1->coord - cn0->coord);

			if ((result & 1) != 0) ++count;
			if ((result & 2) != 0) ++count;
		}
	}
	else {
		split->elements.forEach([&](auto& kde) {
			count += kde->smoothenMetricAtNodes();
		});
	}

	return count;
}

void KdOctElementL::gatherControlNodesNeighbors()
{
	if (isLeaf()) {
		auto cn0 = getCN();

		for (auto nb : leaf->kn.elements) {
			if (nb == nullptr) continue;
			auto cn1 = nb->getNearestLeaf(cn0->coord)->getCN();
			cn0->insertNbIfNew(cn1);
			//cn1->insertNbIfNew(cn0); // ???
		}
	}
	else {
		split->elements.forEach([&](auto kde) {
			kde->gatherControlNodesNeighbors();
		});
	}
}

int KdOctElementL::getTotalBytesExtra() const {
	int total = KdElement::getTotalBytesExtra() + sizeof(KdLeafData*);
	if (leaf != nullptr) total += sizeof(KdLeafData);
	return total;
}

ControlDataMatrix3d KdOctElementL::getCDM(
	const DPoint3d & pt, const DBox & box) const
{
	assert(box.contains(pt));
	assert(leaf != nullptr);
	return leaf->cnode->control_data;
}

ControlDataMatrix3d KdOctElementL::getCDM(
	const DPoint3d & pt) const
{
	assert(leaf != nullptr);
	return leaf->cnode->control_data;
}

bool KdOctElementL::setMinimumCDM(const CDM3d & cdm) {
	bool changed = leaf->cnode->control_data.setMinimum(cdm);
	if (changed) leaf->cnode->setGradationUnknown(); // mark as invalid
	return changed;
}

std::shared_ptr<ControlNode3d> KdOctElementL::getCN() { return leaf->cnode; }

std::shared_ptr<const ControlNode3d> KdOctElementL::getCN() const { return leaf->cnode; }

void KdOctElementL::clear() {
	KdElement::clear();
	leaf.reset();
}

std::shared_ptr<KdElement> KdOctElementL::clone() const {
	assert(isLeaf());
	return std::make_shared<KdOctElementL>(leaf->cnode->coord, leaf->cnode->control_data, leaf->kn);
}

KdOctElementL::KdOctElementL(const DPoint3d & mid, const CDM3d & v, const KdNeighbors& kn)
	: leaf( std::make_unique<KdLeafData>(mid, v, kn)) { }
