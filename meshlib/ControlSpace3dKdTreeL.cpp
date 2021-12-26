#include "ControlSpace3dKdTreeL.h"
#include "DataKdTree.h"

int ControlSpace3dKdTreeL::createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices)
{
	if (m_top != nullptr) clear();
	KdNeighbors kdn; // empty neighbors
	DPoint3d mid = m_box.getMiddlePoint();
	m_top = std::make_shared<KdElementL>(mid, CDM3d::identity, kdn);
	new_vertices.add(m_top->getCN());
	return 1;
}

int KdElementL::adaptToField(
	const std::function<CDM3d(const DPoint3d & pt)> & f,
	const DBox & box, double max_error, int level,
	DataVector<std::shared_ptr<ControlNode3d>> * /* kvalues */)
{
	assert(isLeaf());
	if (level >= KdTree3dParams::paramMaxLevel) return 0;
	double curr_error = estimateError(box, f);
	if (curr_error <= max_error) return 0;

	Axis split_axis;
	double split_value;
	KdTree3dParams::m_cf_splitter->whereToSplit(box, f, split_axis, split_value);

	DataVector<std::shared_ptr<ControlNode3d>> split_vertices(10);
	DBox box0, box1;
	splitElementKd(split_vertices, box, split_axis, split_value, &box0, &box1);

	split_vertices.forEach([&f](std::shared_ptr<ControlNode3d> node) {
		node->control_data = f(node->coord);
	});

	return 2 +
		split->elements[0]->adaptToField(f, box0, max_error, level + 1) +
		split->elements[1]->adaptToField(f, box1, max_error, level + 1);
}

int KdElementL::splitElementKd(
	DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
	const DBox & box, Axis split_axis, double split_value,
	DBox * box0, DBox * box1,
	DataVector<std::shared_ptr<ControlNode3d>> * /* kvalues */)
{
	assert(isLeaf());
	assert(leaf != nullptr && split == nullptr);

	auto leaf0 = this->clone();
	auto leaf1 = this->clone();
	clear();

	assert(box0 != nullptr && box1 != nullptr);
	split = new KdTreeSplitElement(split_axis, split_value, leaf0, leaf1);
	box.split(split_axis, split_value, *box0, *box1);

	leaf0->getCN()->coord = box0->getMiddlePoint();
	leaf1->getCN()->coord = box1->getMiddlePoint();

	split_vertices.add(leaf0->getCN());
	split_vertices.add(leaf1->getCN());

	auto lf0 = (KdElementL*)(leaf0.get());
	auto lf1 = (KdElementL*)(leaf1.get());

	switch(split_axis) {
	case Axis::X:
		lf0->setKN(FC3D_EAST, leaf1);
		lf1->setKN(FC3D_WEST, leaf0);
		break;
	case Axis::Y:
		lf0->setKN(FC3D_NORTH, leaf1);
		lf1->setKN(FC3D_SOUTH, leaf0);
		break;
	case Axis::Z:
		lf0->setKN(FC3D_HIGH, leaf1);
		lf1->setKN(FC3D_LOW, leaf0);
		break;
	}

	return 2;
}

void KdElementL::setKN(const KdNeighbors & kn) {
	assert(leaf != nullptr);
	leaf->kn = kn;
}

void KdElementL::setKN(int fi, std::shared_ptr<KdElement> const & kde) {
	assert(leaf != nullptr);
	leaf->kn[fi] = kde;
}

int KdElementL::smoothenMetricAtNodes() {
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
		count += split->elements[0]->smoothenMetricAtNodes();
		count += split->elements[1]->smoothenMetricAtNodes();
	}

	return count;
}

void KdElementL::gatherControlNodesNeighbors()
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
		split->elements[0]->gatherControlNodesNeighbors();
		split->elements[1]->gatherControlNodesNeighbors();
	}
}

int KdElementL::getTotalBytesExtra() const {
	int total = KdElement::getTotalBytesExtra() + sizeof(KdLeafData*);
	if (leaf != nullptr) total += sizeof(KdLeafData);
	return total;
}

ControlDataMatrix3d KdElementL::getCDM(
	const DPoint3d & pt, const DBox & box) const
{
	assert(box.containsEps(pt));
	assert(leaf != nullptr);
	return leaf->cnode->control_data;
}

ControlDataMatrix3d KdElementL::getCDM(
	const DPoint3d & pt) const
{
	assert(leaf != nullptr);
	return leaf->cnode->control_data;
}

bool KdElementL::setMinimumCDM(const CDM3d & cdm) {
	bool changed = leaf->cnode->control_data.setMinimum(cdm);
	if (changed) leaf->cnode->setGradationUnknown(); // mark as invalid
	return changed;
}

std::shared_ptr<ControlNode3d> KdElementL::getCN() { return leaf->cnode; }

std::shared_ptr<const ControlNode3d> KdElementL::getCN() const { return leaf->cnode; }

void KdElementL::clear() {
	KdElement::clear();
	leaf.reset();
}

std::shared_ptr<KdElement> KdElementL::clone() const {
	assert(isLeaf());
	return std::make_shared<KdElementL>(leaf->cnode->coord, leaf->cnode->control_data, leaf->kn);
}

KdElementL::KdElementL(const DPoint3d & mid, const CDM3d & v, const KdNeighbors& kn)
	: leaf( std::make_unique<KdLeafData>(mid, v, kn)) { }

KdElementL::KdElementL(const DPoint3d & mid, KdLeafData * leaf_data)
	: leaf(leaf_data) {
	leaf->cnode->coord = mid;
}
