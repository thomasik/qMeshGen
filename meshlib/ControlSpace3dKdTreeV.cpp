#include "ControlSpace3dKdTreeV.h"
#include "DataKdTree.h"


/// Standard contructor

ControlSpace3dKdTreeV::ControlSpace3dKdTreeV(const DBox & box) 
	: ControlSpace3dKdTree(box)
{
	m_values = std::make_shared<DataVector<std::shared_ptr<ControlNode3d>>>();
}

int ControlSpace3dKdTreeV::createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices)
{
	if (m_top != nullptr) clear();
	KdVertices kv;
	KdNeighbors kdn;

	for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
		DPoint3d coord = m_box.getVertex(vi);
		kv[vi] = std::make_shared<ControlNode3d>(coord, CDM3d::identity);
		m_values->add(kv[vi]);
		new_vertices.add(kv[vi]);
	}
	m_top = std::make_shared<KdElementV>(kv, kdn);
	return 1;
}

int ControlSpace3dKdTreeV::fixCloseControlNodes() 
{
	if (m_top == nullptr) return 0;
	assert(m_values != nullptr && m_values->notEmpty());
	// gather all contol nodes with neighboring nodes
	int cn_count = m_values->countInt();
	DataVector<int> ref_cn(cn_count);
	for (int i = 0; i < cn_count; i++) ref_cn.add(i);
	int close_count = 0;
	for(int i = 1; i < m_values->countInt(); i++) {
		auto cni = m_values->get(i);
		for (int j = 0; j < i; j++) {
			auto cnj = m_values->get(j);
			if (cni->coord.distance2(cnj->coord) < VERY_SMALL_NUMBER) {
				close_count++;
				ref_cn[i] = j;
				assert(ref_cn[j] == j);
				break;
			}
		}
	}
	if (close_count > 0) {

		DataHashTableKeyValue< std::shared_ptr<ControlNode3d>,
			std::shared_ptr<ControlNode3d> > hcnodes(2 * cn_count, nullptr);
		auto new_values = 
			std::make_shared<DataVector<std::shared_ptr<ControlNode3d>>>(cn_count - close_count);
		for (int i = 0; i < cn_count; i++) {
			auto cni = m_values->get(i);
			int j = ref_cn[i];
			if (i == j) {
				ref_cn[i] = new_values->add(cni);
				hcnodes.setValue(cni, new_values->get(ref_cn[i]));
			}
			else {
				hcnodes.setValue(cni, new_values->get(ref_cn[j]));
			}
		}

		m_values = new_values;

		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Fixed " << close_count << " close control nodes.");
		return m_top->replaceControlNodes(hcnodes);
	}
	else
		return 0;
}

KdElementV::KdElementV(const KdVertices & kv, const KdNeighbors & kn) 
	: leaf(std::make_unique<KdLeafData>(kv,kn)) {}

int KdElementV::adaptToField(
	const std::function<CDM3d(const DPoint3d&pt)>& f, 
	const DBox & box, double max_error, int level, 
	DataVector<std::shared_ptr<ControlNode3d>> * kvalues)
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
	int result = splitElementKd(split_vertices, box, split_axis, split_value, &box0, &box1, kvalues);

	split_vertices.forEach([&f](std::shared_ptr<ControlNode3d> node) {
		node->control_data = f(node->coord);
	});

	return result +
		split->elements[0]->adaptToField(f, box0, max_error, level + 1, kvalues) +
		split->elements[1]->adaptToField(f, box1, max_error, level + 1, kvalues);
}

ControlDataMatrix3d KdElementV::getCDM(
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
	CDM3d cdm = leaf->kv[0]->control_data * coeff[0];
	for(int i = 1; i <  VX3D_COUNT; i++)
		cdm += leaf->kv[i]->control_data * coeff[i];
	return cdm;
}

inline ControlDataMatrix3d KdElementV::getCDM(const DPoint3d & pt) const {
	assert(false);
	// has to be with box...
	return CDM3d::identity;
}

bool KdElementV::setMinimumCDM(const CDM3d & cdm) 
{
	bool any_changes = false;
	for (auto cn : leaf->kv.vertices) {
		bool changed = cn->control_data.setMinimum(cdm);
		if (changed) cn->setGradationUnknown(); // mark as invalid
		any_changes |= changed;
	}
	return any_changes;
}

void KdElementV::setKV(int vi, std::shared_ptr<ControlNode3d> const & cn) {
	assert(leaf != nullptr);
	leaf->kv[vi] = cn;
}

void KdElementV::setKN(int fi, std::shared_ptr<KdElement> const & kde) {
	assert(leaf != nullptr);
	leaf->kn[fi] = kde;
}

void KdElementV::clear() {
	KdElement::clear();
	leaf.reset();
}

int KdElementV::splitElementKd(
	DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
	const DBox & box, Axis split_axis, double split_value,
	DBox * box0, DBox * box1,
	DataVector<std::shared_ptr<ControlNode3d>> *kvalues)
{
	assert(isLeaf());
	//assert(valid());

	auto leaf0 = this->clone();
	auto leaf1 = this->clone();

	box.split(split_axis, split_value, *box0, *box1);
	// update values_0 and values_1
	switch (split_axis) {
	case Axis::X:
	{
		const DPoint3d pts[] = {
			DPoint3d(box0->x1, box.y0, box.z0),
			DPoint3d(box0->x1, box.y1, box.z0),
			DPoint3d(box0->x1, box.y0, box.z1),
			DPoint3d(box0->x1, box.y1, box.z1)
		};
		assert(box0->x1 == split_value);
		OctVertexWhich V_LEAF0[] = { VX3D_LSE, VX3D_LNE, VX3D_HSE, VX3D_HNE };
		OctVertexWhich V_LEAF1[] = { VX3D_LSW, VX3D_LNW, VX3D_HSW, VX3D_HNW };
		OctEdgeWhich E_LEAF[] = { ED3D_LS, ED3D_LN, ED3D_HS, ED3D_HN };
		double min_dist2 = LARGE_NUMBER;
		for (int i = 0; i < 4; i++) {
			auto v = getNearestControlNodeForEdge(box, pts[i], E_LEAF[i], &min_dist2);
			if (!v || min_dist2 > SMALL_NUMBER) {
				v = std::make_shared<ControlNode3d>(pts[i], getCDM(pts[i], box));
				kvalues->add(v);
				split_vertices.add(v);
			}
			leaf0->setKV(V_LEAF0[i], v);
			leaf1->setKV(V_LEAF1[i], v);
		}
		break;
	}
	case Axis::Y:
	{
		const DPoint3d pts[] = {
			DPoint3d(box.x0, box0->y1, box.z0),
			DPoint3d(box.x1, box0->y1, box.z0),
			DPoint3d(box.x0, box0->y1, box.z1),
			DPoint3d(box.x1, box0->y1, box.z1)
		};
		assert(box0->y1 == split_value);
		OctVertexWhich V_LEAF0[] = { VX3D_LNW, VX3D_LNE, VX3D_HNW, VX3D_HNE };
		OctVertexWhich V_LEAF1[] = { VX3D_LSW, VX3D_LSE, VX3D_HSW, VX3D_HSE };
		OctEdgeWhich E_LEAF[] = { ED3D_LW, ED3D_LE, ED3D_HW, ED3D_HE };
		double min_dist2 = LARGE_NUMBER;
		for (int i = 0; i < 4; i++){
			auto v = getNearestControlNodeForEdge(box, pts[i], E_LEAF[i], &min_dist2);
			if (!v || min_dist2 > SMALL_NUMBER) {
				v = std::make_shared<ControlNode3d>(pts[i], getCDM(pts[i], box));
				kvalues->add(v);
				split_vertices.add(v);
			}
			leaf0->setKV(V_LEAF0[i], v);
			leaf1->setKV(V_LEAF1[i], v);
		}
		break;
	}
	case Axis::Z:
	{
		const DPoint3d pts[] = {
			DPoint3d(box.x0, box.y0, box0->z1),
			DPoint3d(box.x1, box.y0, box0->z1),
			DPoint3d(box.x0, box.y1, box0->z1),
			DPoint3d(box.x1, box.y1, box0->z1)
		};
		assert(box0->z1 == split_value);
		OctVertexWhich V_LEAF0[] = { VX3D_HSW, VX3D_HSE, VX3D_HNW, VX3D_HNE };
		OctVertexWhich V_LEAF1[] = { VX3D_LSW, VX3D_LSE, VX3D_LNW, VX3D_LNE };
		OctEdgeWhich E_LEAF[] = { ED3D_SW, ED3D_SE, ED3D_NW, ED3D_NE };
		double min_dist2 = LARGE_NUMBER;
		for (int i = 0; i < 4; i++) {
			auto v = getNearestControlNodeForEdge(box, pts[i], E_LEAF[i], &min_dist2);
			if (!v || min_dist2 > SMALL_NUMBER) {
				v = std::make_shared<ControlNode3d>(pts[i], getCDM(pts[i], box));
				kvalues->add(v);
				split_vertices.add(v);
			}
			leaf0->setKV(V_LEAF0[i], v);
			leaf1->setKV(V_LEAF1[i], v);
		}
		break;
	}
	}

	// set neighbors
	switch (split_axis) {
	case Axis::X:
		leaf0->setKN(FC3D_EAST, leaf1);
		leaf1->setKN(FC3D_WEST, leaf0);
		break;
	case Axis::Y:
		leaf0->setKN(FC3D_NORTH, leaf1);
		leaf1->setKN(FC3D_SOUTH, leaf0);
		break;
	case Axis::Z:
		leaf0->setKN(FC3D_HIGH, leaf1);
		leaf1->setKN(FC3D_LOW, leaf0);
		break;
	}

	clear();
	split = new KdTreeSplitElement(split_axis, split_value, leaf0, leaf1);

	return 2;
}

int KdElementV::smoothenMetricAtNodes()
{
	int count = 0;
	if (isLeaf()) {
		// check along leaf-edges
		for (int i = ED3D_FIRST; i <= ED3D_LAST; i++) {
			int v0 = DBox::edge_to_vertex[i][0];
			int v1 = DBox::edge_to_vertex[i][1];

			auto cn0 = leaf->kv[v0];
			auto cn1 = leaf->kv[v1];

			int result = ControlSpace3dAdaptive::smoothenMetricForNodes(
				cn0.get(), cn1.get(), cn1->coord - cn0->coord);

			if ((result & 1) != 0) ++count;
			if ((result & 2) != 0) ++count;
		}
	}else{
		count += split->elements[0]->smoothenMetricAtNodes();
		count += split->elements[1]->smoothenMetricAtNodes();
	}

	return count;
}

std::shared_ptr<KdElement> KdElementV::clone() const {
	assert(isLeaf());
	return std::make_shared<KdElementV>(leaf->kv, leaf->kn);
}

void KdElementV::gatherControlNodesNeighbors()
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
		split->elements[0]->gatherControlNodesNeighbors();
		split->elements[1]->gatherControlNodesNeighbors();
	}
}

int KdElementV::getTotalBytesExtra() const {
	int total = KdElement::getTotalBytesExtra() + sizeof(KdVertices*);
	if (leaf != nullptr) total += sizeof(KdLeafData);
	return total;
}

std::shared_ptr<ControlNode3d> KdElementV::getNearestVertexRef(const DPoint3d & pt, const DBox & box)
{
	DBox vbox = box;
	KdElementV* kde = (KdElementV*)getNearestLeaf(pt, &vbox);
	int vi = vbox.getMiddlePoint().vertexDirection(pt);
	return kde->leaf->kv[vi];
}

std::shared_ptr<ControlNode3d> KdElementV::getNearestControlNodeForEdge(
	const DBox & box, const DPoint3d & edge_pt, int ei, double * ptr_min_dist2)
{
	std::shared_ptr<ControlNode3d> nearest_cn;
	double min_dist2 = LARGE_NUMBER;
	for (int i = 0; i < 2; i++) {
		OctFaceWhich fi = DBox::edge_to_face[ei][i];
		auto f_kn = (KdElementV*)leaf->kn[fi].get();
		if (f_kn && !f_kn->isLeaf()) {
			auto cn = f_kn->getNearestVertexRef(edge_pt, box);
			if (cn == nullptr) continue;
			double dist2 = cn->coord.distance2(edge_pt);
			if (nearest_cn == nullptr || dist2 < min_dist2) {
				nearest_cn = cn;
				min_dist2 = dist2;
			}
		}
	}
	if (ptr_min_dist2 != nullptr) *ptr_min_dist2 = min_dist2;
	return nearest_cn;
}

int KdElementV::replaceControlNodes(
	const DataHashTableKeyValue<std::shared_ptr<ControlNode3d>, std::shared_ptr<ControlNode3d>>& hcnodes)
{
	int count = 0;
	if (isLeaf()) {
		for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
			auto cn = hcnodes.getValue(leaf->kv[vi], nullptr);
			if (cn) {
				leaf->kv[vi] = cn;
				count++;
			}
		}
	}
	else {
		count += split->elements[0]->replaceControlNodes(hcnodes);
		count += split->elements[1]->replaceControlNodes(hcnodes);
	}

	return count;
}
