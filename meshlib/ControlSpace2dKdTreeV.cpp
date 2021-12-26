#include "ControlSpace2dKdTreeV.h"
#include "DataKdTree.h"

#define SAMPLING_SIZE 4
#define MAX_LEVEL 20

int ControlSpace2dKdTreeV::adaptToField(
	const std::function<CDM2d(const DPoint2d&pt)>& f, double max_error) 
{
	if (m_top != nullptr) clear();
	KdElement::KdVertices kv;

	DPoint2d coord(m_box.x0, m_box.y0);
	kv.vx_sw = &m_values.get(m_values.add(ControlNode2d(coord, f(coord))));
	coord = DPoint2d(m_box.x1, m_box.y0);
	kv.vx_se = &m_values.get(m_values.add(ControlNode2d(coord, f(coord))));
	coord = DPoint2d(m_box.x0, m_box.y1);
	kv.vx_nw = &m_values.get(m_values.add(ControlNode2d(coord, f(coord))));
	coord = DPoint2d(m_box.x1, m_box.y1);
	kv.vx_ne = &m_values.get(m_values.add(ControlNode2d(coord, f(coord))));

	m_top = new KdElement(kv);
	return 1 + m_top->adaptToField(f, m_box, max_error, m_values);
}

void ControlSpace2dKdTreeV::clear() {
	if (m_top != nullptr) {
		m_top->clear();
		delete m_top;
		m_top = nullptr;
		m_values.clear();
	}
}

ControlDataMatrix2d ControlSpace2dKdTreeV::getMetricAtPoint(const DPoint2d & pt) const {
	DRect box;
	KdElement* kd_leaf = findLastLeaf(pt, box);
	return kd_leaf->getValue(pt, box);
}

// stat

int ControlSpace2dKdTreeV::getMaxBalance() const {
	int max_balance = 0;
	if (m_top) m_top->getMaxDepthAndBalance(max_balance);
	return max_balance;
}

int ControlSpace2dKdTreeV::getTotalBytes() const {
	return sizeof(*this) // main struct
		+ sizeof(KdElement) * getElementCount() // tree nodes
		+ sizeof(ControlNode2d) * getValueCount(); // values
}

int ControlSpace2dKdTreeV::KdElement::adaptToField(const std::function<CDM2d(const DPoint2d&pt)>& f, 
	const DRect & box, double max_error, DataVector<ControlNode2d>& kvalues, int level)
{

	//static size_t counter = 0;
	//if (++counter % 100 == 0) {
	//	LOG4CPLUS_INFO(MeshLog::logger_console, "** counter", counter);
	//	LOG4CPLUS_INFO(MeshLog::logger_console, "** box-diam", box.getDiameter());
	//}

	assert(is_leaf);
	if (level >= MAX_LEVEL) return 0;
	double curr_error =
		Kd2dRegularErrorEstimator<CDM2d, ControlSpace2dKdTreeV::KdElement, SAMPLING_SIZE>::estimateError(
			box, this, f);
	if (curr_error <= max_error) return 0;

	Axis split_axis;
	double split_value;
	Kd2dMaxGradientSplitter<CDM2d, SAMPLING_SIZE>::whereToSplit(box, f, split_axis, split_value);

	DataVector<ControlNode2d*> split_vertices(10);
	DRect box0, box1;
	split(kvalues, split_vertices, box, split_axis, split_value, box0, box1);

	split_vertices.forEach([&f](ControlNode2d* node) {
		node->control_data = f(node->coord);
	});

	return 2 +
		m_data.split.element_0->adaptToField(f, box0, max_error, kvalues, level + 1) +
		m_data.split.element_1->adaptToField(f, box1, max_error, kvalues, level + 1);
}

ControlDataMatrix2d ControlSpace2dKdTreeV::KdElement::getValue(
	const DPoint2d & pt, const DRect & box) const
{
	assert(box.contains(pt));

	double tx = (pt.x - box.x0) / box.getDX();
	double ty = (pt.y - box.y0) / box.getDY();

	double coeff[4] = {
		(1 - tx) * (1 - ty) ,	// VX_SW
		tx	 * (1 - ty) ,		// VX_SE
		(1 - tx) *    ty ,		// VX_NW
		tx	 *    ty 			// VX_NE
	};
	return
		m_data.values.vx_sw->control_data * coeff[0] +
		m_data.values.vx_se->control_data * coeff[1] +
		m_data.values.vx_nw->control_data * coeff[2] +
		m_data.values.vx_ne->control_data * coeff[3];
}

void ControlSpace2dKdTreeV::KdElement::clear() {
	if (!is_leaf) {
		m_data.split.element_0->clear();
		delete m_data.split.element_0;
		m_data.split.element_1->clear();
		delete m_data.split.element_1;
		is_leaf = true;
	}
}

int ControlSpace2dKdTreeV::KdElement::getElementCount() const {
	return is_leaf ? 1 :
		1 + m_data.split.element_0->getElementCount() + m_data.split.element_1->getElementCount();
}

int ControlSpace2dKdTreeV::KdElement::getMaxDepth() const {
	return is_leaf ? 1 :
		1 + std::max(m_data.split.element_0->getMaxDepth(), m_data.split.element_1->getMaxDepth());
}

int ControlSpace2dKdTreeV::KdElement::getMaxDepthAndBalance(int & max_balance) const {
	if (is_leaf) return 1;
	int left_depth = m_data.split.element_0->getMaxDepthAndBalance(max_balance);
	int right_depth = m_data.split.element_1->getMaxDepthAndBalance(max_balance);
	int balance = std::abs(left_depth - right_depth);
	if (balance > max_balance) max_balance = balance;
	return 1 + std::max(left_depth, right_depth);
}

ControlSpace2dKdTreeV::KdElement * ControlSpace2dKdTreeV::findLastLeaf(
	const DPoint2d & pt, DRect & box) const
{
	box = m_box;
	int level = 0;
	return m_top->findLastLeaf(pt, box, level);
}

ControlSpace2dKdTreeV::KdElement * ControlSpace2dKdTreeV::findLastLeaf(
	const DPoint2d & pt, DRect & box, int & level) const
{
	box = m_box;
	level = 0;
	return m_top->findLastLeaf(pt, box, level);
}

ControlSpace2dKdTreeV::KdElement * ControlSpace2dKdTreeV::KdElement::findLastLeaf(
	const DPoint2d & pt, DRect & box, int & level)
{
	KdElement* kde = this;
	while (!kde->is_leaf) {
		bool split_lower = pt[kde->m_data.split.axis] <= kde->m_data.split.value;
		box.splitAndUpdate(kde->m_data.split.axis, kde->m_data.split.value, split_lower);
		assert(box.contains(pt));
		kde = split_lower ? kde->m_data.split.element_0 : kde->m_data.split.element_1;
		++level;
	}

	return kde;
}

void ControlSpace2dKdTreeV::forEachControlNode(const std::function<void(const ControlNode2d & node)>& fg) const
{
	m_values.forEach(fg);
}

void ControlSpace2dKdTreeV::forEachControlNode(const std::function<void(ControlNode2d & node)>& fg)
{
	m_values.forEach(fg);
}

/// Refines control space at the given point
bool ControlSpace2dKdTreeV::setMinControl(
	const DPoint2d & pt, const ControlDataMatrix2d & cdm, bool min_value_set) 
{
	assert(m_initialized);
	//assert(valid());
	DRect box;
	int level;
	KdElement* kd_leaf = findLastLeaf(pt, box, level);
	return kd_leaf->adaptAndSetControlPoint(m_values, ControlNode2d(pt, cdm), base_surface, box, level, min_value_set);
}

bool ControlSpace2dKdTreeV::KdElement::adaptAndSetControlPoint(
	DataVector<ControlNode2d>& tree_vertices, const ControlNode2d& qv, 
	SurfaceConstPtr surface, DRect& box, int & level, bool min_value_set)
{
	if (!is_leaf) {
		KdElement* kd_leaf = findLastLeaf(qv.coord, box, level);
		return kd_leaf->adaptAndSetControlPoint(tree_vertices, qv, surface, box, level, min_value_set);
	}
	// compare new source node with already existing
	// check with calculated value:
	const ControlDataMatrix2d cdm = getMetricAtPoint(qv.coord, box);
	double diff = qv.control_data.countDifferenceRR(cdm);
	bool any_changes = false;
	if (diff > ControlSpace2dAdaptive::param_threshold_diff) {
		// split may be unnecessary, if qv.control_data is greater than already set for this leaf...
		ControlDataMatrix2d min_cdm = cdm;
		min_cdm.setMinimum(qv.control_data);
		if (cdm.countDifferenceRR(min_cdm) < ControlSpace2dAdaptive::param_threshold_diff) {
			// qv different from current metric, but introduces no changes - so skip
			return false;
		}
		DMetric2d dmp(qv.control_data, surface, qv.coord);
		DVector2d len2(
			dmp.transformPStoMS(DVector2d(box.getDX(), 0.0)).length2(),
			dmp.transformPStoMS(DVector2d(0.0, box.getDY())).length2());
		bool split_needed = (len2.x > METRIC_LENGTH_RATIO2 || len2.y > METRIC_LENGTH_RATIO2);
		if (split_needed && (level < MAX_LEVEL)) {

			Axis split_axis;
			double split_value;
			Kd2dHalfLongestSplitter<CDM2d>::whereToSplit(box, split_axis, split_value);

			DataVector<ControlNode2d*> split_vertices(10);
			DRect box0, box1;
			if (split(tree_vertices, split_vertices, box, split_axis, split_value, box0, box1)) { // split (+ set values for new nodes)

				//split_vertices.forEach([&f](ControlNode2d* node) {
				//	node->control_data = f(node->coord);
				//});

				KdElement* kd_leaf = findLastLeaf(qv.coord, box, level);
				kd_leaf->adaptAndSetControlPoint(tree_vertices, qv, surface, box, level, min_value_set);
				return true;
			}
		}
		if (min_value_set) {
			any_changes |= m_data.values.vx_sw->control_data.setMinimum(qv.control_data);
			any_changes |= m_data.values.vx_se->control_data.setMinimum(qv.control_data);
			any_changes |= m_data.values.vx_nw->control_data.setMinimum(qv.control_data);
			any_changes |= m_data.values.vx_ne->control_data.setMinimum(qv.control_data);
		}
	}

	return any_changes;
}

bool ControlSpace2dKdTreeV::KdElement::split(
	DataVector<ControlNode2d>& tree_vertices,
	DataVector<ControlNode2d*>& split_vertices, 
	const DRect & box, Axis split_axis, double split_value,
	DRect & box0, DRect & box1)
{
	assert(is_leaf);
	//assert(valid());

	KdVertices kv0 = m_data.values;
	KdVertices kv1 = m_data.values;

	is_leaf = false;

	m_data.split.axis = split_axis;
	m_data.split.value = split_value;

	box.split(m_data.split.axis, m_data.split.value, box0, box1);
	// update values_0 and values_1
	switch (m_data.split.axis) {
	case Axis::X:
	{
		const DPoint2d pts[] = {
			DPoint2d(box0.x1, box.y0),
			DPoint2d(box0.x1, box.y1)
		};
		assert(box0.x1 == m_data.split.value);
		ControlNode2d* v[2];
		for (int i = 0; i < 2; i++)
			split_vertices.add(v[i] = &tree_vertices.get(tree_vertices.add(
				ControlNode2d(pts[i], getValue(pts[i], box)))));
		// east for kv0
		kv0.vx_se = v[0];
		kv0.vx_ne = v[1];
		// west for kv1
		kv1.vx_sw = v[0];
		kv1.vx_nw = v[1];
		break;
	}
	case Axis::Y:
	{
		const DPoint2d pts[] = {
			DPoint2d(box.x0, box0.y1),
			DPoint2d(box.x1, box0.y1)
		};
		assert(box0.y1 == m_data.split.value);
		ControlNode2d* v[4];
		for (int i = 0; i < 2; i++)
			split_vertices.add(v[i] = &tree_vertices.get(tree_vertices.add(
				ControlNode2d(pts[i], getValue(pts[i], box)))));
		// north for kv0
		kv0.vx_nw = v[0];
		kv0.vx_ne = v[1];
		// south for kv1
		kv1.vx_sw = v[0];
		kv1.vx_se = v[1];
		break;
	}
	case Axis::Z:
	{
		assert(false);
		break;
	}
	}
	m_data.split.element_0 = new KdElement(kv0);
	m_data.split.element_1 = new KdElement(kv1);

	return true;
}
