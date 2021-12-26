#include "ControlSpace3dOctreeV.h"
#include "DataKdTree.h"


/// Standard contructor

ControlSpace3dKdOctreeV::ControlSpace3dKdOctreeV(const DBox & box) : ControlSpace3dKdTree(box) 
{ 
	m_values = std::make_shared<DataVector<std::shared_ptr<ControlNode3d>>>();
}

int ControlSpace3dKdOctreeV::createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices)
{
	if (m_top != nullptr) clear();
	KdNeighbors kn; // empty
	KdVertices kv;
	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
		DPoint3d coord = m_box.getVertex(vi);
		auto scn = std::make_shared<ControlNode3d>(coord, CDM3d::identity);
		kv.vertices[vi] = scn;
		m_values->add(scn);
		new_vertices.add(scn);
	}
	m_top = std::make_shared<KdOctElementV>(kn, kv);
	return 1;
}

KdOctElementV::KdOctElementV(const KdNeighbors& _kn, const KdVertices& _kv) : 
	leaf( std::make_unique<KdLeafData>(_kn, _kv)) {}

int KdOctElementV::adaptToField(
	const std::function<CDM3d(const DPoint3d&pt)>& f,
	const DBox & box, double max_error, int level,
	DataVector<std::shared_ptr<ControlNode3d>> * kvalues)
{
	assert(isLeaf());
	if (level >= KdOctree3dParams::paramMaxLevel) return 0;
	double curr_error = estimateError(box, f);
	if (curr_error <= max_error) return 0;

	int result = 0;
	DataVector<std::shared_ptr<ControlNode3d>> split_vertices(10);
	DBox boxes[VX3D_COUNT];
	result += splitElementOct(split_vertices, box, boxes, kvalues);

	split_vertices.forEach([&f](const std::shared_ptr<ControlNode3d>& node) {
		node->control_data = f(node->coord);
	});

	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
		result += split->elements[vi]->adaptToField(f, boxes[vi], max_error, level+1, kvalues);
	}

	return result;
}

ControlDataMatrix3d KdOctElementV::getCDM(
	const DPoint3d & pt, const DBox & box) const
{
	assert(box.containsEps(pt));

	double tx = (pt.x - box.x0) / box.getDX();
	double ty = (pt.y - box.y0) / box.getDY();
	double tz = (pt.z - box.z0) / box.getDZ();

	double coeff[8] = {
		(1 - tx) * (1 - ty) * (1 - tz),	// VX_LSW
		tx	 * (1 - ty) * (1 - tz),	// VX_LSE
		(1 - tx) *    ty	* (1 - tz),	// VX_LNW
		tx	 *    ty	* (1 - tz),	// VX_LNE
		(1 - tx) * (1 - ty) *    tz ,	// VX_HSW
		tx	 * (1 - ty) *    tz ,	// VX_HSE
		(1 - tx) *    ty	*    tz ,	// VX_HNW
		tx	 *    ty	*    tz 	// VX_HNE
	};
	CDM3d cdm = leaf->kv.vertices[0]->control_data * coeff[0];
	for (int i = 1; i < VX3D_COUNT; i++)
		cdm += leaf->kv.vertices[i]->control_data * coeff[i];
	return cdm;
}

inline ControlDataMatrix3d KdOctElementV::getCDM(const DPoint3d & pt) const {
	assert(false);
	// has to be with box...
	return CDM3d::identity;
}

bool KdOctElementV::setMinimumCDM(const CDM3d & cdm)
{
	bool any_changes = false;
	for (auto cn : leaf->kv.vertices) {
		bool changed = cn->control_data.setMinimum(cdm);
		if (changed) cn->setGradationUnknown(); // mark as invalid
		any_changes |= changed;
	}
	return any_changes;
}

void KdOctElementV::setKV(int vi, std::shared_ptr<ControlNode3d> const & cn) {
	assert(leaf != nullptr);
	leaf->kv.vertices[vi] = cn;
}

void KdOctElementV::setKV(const KdVertices & kv) {
	assert(leaf != nullptr);
	leaf->kv = kv;
}

void KdOctElementV::setKN(const KdNeighbors & kn) {
	assert(leaf != nullptr);
	leaf->kn = kn;
}

void KdOctElementV::setKN(int fi, std::shared_ptr<KdElement> const & kde) {
	assert(leaf != nullptr);
	leaf->kn[fi] = kde;
}



void KdOctElementV::clear() {
	KdElement::clear();
	leaf.reset();
}

int KdOctElementV::splitElementOct(
	DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
	const DBox & box, DBox * boxes,
	DataVector<std::shared_ptr<ControlNode3d>> *kvalues)
{
	assert(isLeaf());
	assert(leaf != nullptr && split == nullptr);

	KdVertices split_kv = leaf->kv;
	KdNeighbors split_kn = leaf->kn;
	std::shared_ptr<KdElement> leaves[VX3D_COUNT];
	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++)
		leaves[vi] = this->clone();
	DPoint3d split_middle = box.getMiddlePoint();

	// -> vertices
	auto mid_box_vertex = std::make_shared<ControlNode3d>(split_middle, getCDM(split_middle, box));
	kvalues->add(mid_box_vertex);
	std::shared_ptr<ControlNode3d> mid_face_vertices[FC3D_COUNT] = { nullptr };
	std::shared_ptr<ControlNode3d> mid_edge_vertices[ED3D_COUNT] = { nullptr };
	for (int fi = FC3D_FIRST; fi <= FC3D_LAST; fi++) {
		DPoint3d face_mid_pt = box.getMiddlePointForFace(fi);
		auto f_kn = (KdOctElementV*)split_kn[fi].get();
		if (f_kn && !f_kn->isLeaf()) {
			mid_face_vertices[fi] = f_kn->getNearestVertexRef(face_mid_pt, box);
		}
		else {
			auto scn = std::make_shared<ControlNode3d>(face_mid_pt, getCDM(face_mid_pt, box));
			mid_face_vertices[fi] = scn;
			kvalues->add(scn);
		}
	}

	for (int ei = ED3D_FIRST; ei <= ED3D_LAST; ei++) {
		DPoint3d edge_mid_pt = box.getMiddlePointForEdge(ei);
		OctFaceWhich fi = DBox::edge_to_face[ei][0];
		auto f_kn = (KdOctElementV*)split_kn[fi].get();
		bool available = (f_kn && !f_kn->isLeaf());
		if (!available) {
			fi = DBox::edge_to_face[ei][1];
			f_kn = (KdOctElementV*)split_kn[fi].get();
			available = (f_kn && !f_kn->isLeaf());
		}
		if (available) {
			mid_edge_vertices[ei] = f_kn->getNearestVertexRef(edge_mid_pt, box);
		}
		else {
			auto scn = std::make_shared<ControlNode3d>(edge_mid_pt, getCDM(edge_mid_pt, box));
			mid_edge_vertices[ei] = scn;
			kvalues->add(scn);
		}
	}


	clear();

	split = new KdOctreeSplitElement(split_middle);
	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
		boxes[vi] = box.splitOct(split_middle, vi);
		split->elements.add(leaves[vi]);
	}

	KdVertices sub_kv[VX3D_COUNT];
	sub_kv[VX3D_LSW][VX3D_LSW] = split_kv[VX3D_LSW];
	sub_kv[VX3D_LSW][VX3D_LSE] = mid_edge_vertices[ED3D_LS];
	sub_kv[VX3D_LSW][VX3D_LNW] = mid_edge_vertices[ED3D_LW];
	sub_kv[VX3D_LSW][VX3D_LNE] = mid_face_vertices[FC3D_LOW];
	sub_kv[VX3D_LSW][VX3D_HSW] = mid_edge_vertices[ED3D_SW];
	sub_kv[VX3D_LSW][VX3D_HSE] = mid_face_vertices[FC3D_SOUTH];
	sub_kv[VX3D_LSW][VX3D_HNW] = mid_face_vertices[FC3D_WEST];
	sub_kv[VX3D_LSW][VX3D_HNE] = mid_box_vertex;

	sub_kv[VX3D_LSE][VX3D_LSW] = mid_edge_vertices[ED3D_LS];
	sub_kv[VX3D_LSE][VX3D_LSE] = split_kv[VX3D_LSE];
	sub_kv[VX3D_LSE][VX3D_LNW] = mid_face_vertices[FC3D_LOW];
	sub_kv[VX3D_LSE][VX3D_LNE] = mid_edge_vertices[ED3D_LE];
	sub_kv[VX3D_LSE][VX3D_HSW] = mid_face_vertices[FC3D_SOUTH];
	sub_kv[VX3D_LSE][VX3D_HSE] = mid_edge_vertices[ED3D_SE];
	sub_kv[VX3D_LSE][VX3D_HNW] = mid_box_vertex;
	sub_kv[VX3D_LSE][VX3D_HNE] = mid_face_vertices[FC3D_EAST];

	sub_kv[VX3D_LNW][VX3D_LSW] = mid_edge_vertices[ED3D_LW];
	sub_kv[VX3D_LNW][VX3D_LSE] = mid_face_vertices[FC3D_LOW];
	sub_kv[VX3D_LNW][VX3D_LNW] = split_kv[VX3D_LNW];
	sub_kv[VX3D_LNW][VX3D_LNE] = mid_edge_vertices[ED3D_LN];
	sub_kv[VX3D_LNW][VX3D_HSW] = mid_face_vertices[FC3D_WEST];
	sub_kv[VX3D_LNW][VX3D_HSE] = mid_box_vertex;
	sub_kv[VX3D_LNW][VX3D_HNW] = mid_edge_vertices[ED3D_NW];
	sub_kv[VX3D_LNW][VX3D_HNE] = mid_face_vertices[FC3D_NORTH];

	sub_kv[VX3D_LNE][VX3D_LSW] = mid_face_vertices[FC3D_LOW];
	sub_kv[VX3D_LNE][VX3D_LSE] = mid_edge_vertices[ED3D_LE];
	sub_kv[VX3D_LNE][VX3D_LNW] = mid_edge_vertices[ED3D_LN];
	sub_kv[VX3D_LNE][VX3D_LNE] = split_kv[VX3D_LNE];
	sub_kv[VX3D_LNE][VX3D_HSW] = mid_box_vertex;
	sub_kv[VX3D_LNE][VX3D_HSE] = mid_face_vertices[FC3D_EAST];
	sub_kv[VX3D_LNE][VX3D_HNW] = mid_face_vertices[FC3D_NORTH];
	sub_kv[VX3D_LNE][VX3D_HNE] = mid_edge_vertices[ED3D_NE];

	sub_kv[VX3D_HSW][VX3D_LSW] = mid_edge_vertices[ED3D_SW];
	sub_kv[VX3D_HSW][VX3D_LSE] = mid_face_vertices[FC3D_SOUTH];
	sub_kv[VX3D_HSW][VX3D_LNW] = mid_face_vertices[FC3D_WEST];
	sub_kv[VX3D_HSW][VX3D_LNE] = mid_box_vertex;
	sub_kv[VX3D_HSW][VX3D_HSW] = split_kv[VX3D_HSW];
	sub_kv[VX3D_HSW][VX3D_HSE] = mid_edge_vertices[ED3D_HS];
	sub_kv[VX3D_HSW][VX3D_HNW] = mid_edge_vertices[ED3D_HW];
	sub_kv[VX3D_HSW][VX3D_HNE] = mid_face_vertices[FC3D_HIGH];

	sub_kv[VX3D_HSE][VX3D_LSW] = mid_face_vertices[FC3D_SOUTH];
	sub_kv[VX3D_HSE][VX3D_LSE] = mid_edge_vertices[ED3D_SE];
	sub_kv[VX3D_HSE][VX3D_LNW] = mid_box_vertex;
	sub_kv[VX3D_HSE][VX3D_LNE] = mid_face_vertices[FC3D_EAST];
	sub_kv[VX3D_HSE][VX3D_HSW] = mid_edge_vertices[ED3D_HS];
	sub_kv[VX3D_HSE][VX3D_HSE] = split_kv[VX3D_HSE];
	sub_kv[VX3D_HSE][VX3D_HNW] = mid_face_vertices[FC3D_HIGH];
	sub_kv[VX3D_HSE][VX3D_HNE] = mid_edge_vertices[ED3D_HE];

	sub_kv[VX3D_HNW][VX3D_LSW] = mid_face_vertices[FC3D_WEST];
	sub_kv[VX3D_HNW][VX3D_LSE] = mid_box_vertex;
	sub_kv[VX3D_HNW][VX3D_LNW] = mid_edge_vertices[ED3D_NW];
	sub_kv[VX3D_HNW][VX3D_LNE] = mid_face_vertices[FC3D_NORTH];
	sub_kv[VX3D_HNW][VX3D_HSW] = mid_edge_vertices[ED3D_HW];
	sub_kv[VX3D_HNW][VX3D_HSE] = mid_face_vertices[FC3D_HIGH];
	sub_kv[VX3D_HNW][VX3D_HNW] = split_kv[VX3D_HNW];
	sub_kv[VX3D_HNW][VX3D_HNE] = mid_edge_vertices[ED3D_HN];

	sub_kv[VX3D_HNE][VX3D_LSW] = mid_box_vertex;
	sub_kv[VX3D_HNE][VX3D_LSE] = mid_face_vertices[FC3D_EAST];
	sub_kv[VX3D_HNE][VX3D_LNW] = mid_face_vertices[FC3D_NORTH];
	sub_kv[VX3D_HNE][VX3D_LNE] = mid_edge_vertices[ED3D_NE];
	sub_kv[VX3D_HNE][VX3D_HSW] = mid_face_vertices[FC3D_HIGH];
	sub_kv[VX3D_HNE][VX3D_HSE] = mid_edge_vertices[ED3D_HE];
	sub_kv[VX3D_HNE][VX3D_HNW] = mid_edge_vertices[ED3D_HN];
	sub_kv[VX3D_HNE][VX3D_HNE] = split_kv[VX3D_HNE];

	KdNeighbors sub_kn[VX3D_COUNT]; // empty for now
	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
		split->elements[vi]->setKV(sub_kv[vi]);
		split->elements[vi]->setKN(sub_kn[vi]);
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

int KdOctElementV::smoothenMetricAtNodes()
{
	int count = 0;
	if (isLeaf()) {
		// check along leaf-edges
		for (int i = ED3D_FIRST; i <= ED3D_LAST; i++) {
			int v0 = DBox::edge_to_vertex[i][0];
			int v1 = DBox::edge_to_vertex[i][1];

			auto cn0 = leaf->kv.vertices[v0];
			auto cn1 = leaf->kv.vertices[v1];

			int result = ControlSpace3dAdaptive::smoothenMetricForNodes(
				cn0.get(), cn1.get(), cn1->coord - cn0->coord);

			if ((result & 1) != 0) ++count;
			if ((result & 2) != 0) ++count;
		}
	}
	else {
		split->elements.forEach([&](auto& kde) {
			count += kde->smoothenMetricAtNodes(); });
	}

	return count;
}

std::shared_ptr<KdElement> KdOctElementV::clone() const {
	assert(isLeaf());
	return std::make_shared<KdOctElementV>(leaf->kn, leaf->kv);
}

std::shared_ptr<ControlNode3d> KdOctElementV::getNearestVertexRef(const DPoint3d & pt, const DBox & box)
{
	DBox vbox = box;
	KdOctElementV* kde = (KdOctElementV*)getNearestLeaf(pt, &vbox);
	int vi = vbox.getMiddlePoint().vertexDirection(pt);
	return kde->leaf->kv[vi];
}

void KdOctElementV::gatherControlNodesNeighbors()
{
	if (isLeaf()) {
		for (int i = ED3D_FIRST; i <= ED3D_LAST; i++) {
			int v0 = DBox::edge_to_vertex[i][0];
			int v1 = DBox::edge_to_vertex[i][1];

			auto cn0 = leaf->kv[v0];
			auto cn1 = leaf->kv[v1];

			assert(cn0);
			assert(cn1);

			cn0->insertNbIfNew(cn1);
			cn1->insertNbIfNew(cn0);
		}
	}
	else {
		split->elements.forEach([&](auto kde) {
			kde->gatherControlNodesNeighbors(); });
	}
}

int KdOctElementV::getTotalBytesExtra() const {
	int total = KdElement::getTotalBytesExtra() + sizeof(KdLeafData*);
	if (leaf != nullptr) total += sizeof(KdLeafData);
	return total;
}

