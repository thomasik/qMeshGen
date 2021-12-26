#include "ControlSpace2dKdTreeL.h"
#include "DataKdTree.h"

#define SAMPLING_SIZE 4
#define MAX_LEVEL 20

int ControlSpace2dKdTreeL::adaptToField(const std::function<CDM2d(const DPoint2d&pt)>& f, double max_error) {
	DPoint2d mid = m_box.getMiddlePoint();
	if (m_top != nullptr) clear();
	m_top = new KdElement(mid, f(mid));
	return 1 + m_top->adaptToField(f, m_box, max_error);
}

void ControlSpace2dKdTreeL::clear() {
	if (m_top != nullptr) {
		m_top->clear();
		delete m_top;
		m_top = nullptr;
	}
}

int ControlSpace2dKdTreeL::getMaxBalance() const {
	int max_balance = 0;
	if (m_top) m_top->getMaxDepthAndBalance(max_balance);
	return max_balance;
}

int ControlSpace2dKdTreeL::getTotalBytes() const {
	return sizeof(*this) // main struct
		+ sizeof(KdElement) * getElementCount(); // tree nodes
}

/// Get sizing info (matrix mode) at the given point

ControlDataMatrix2d ControlSpace2dKdTreeL::getMetricAtPoint(const DPoint2d & pt) const {
	return findLastLeaf(pt)->getValue(pt);
}

ControlSpace2dKdTreeL::KdElement* ControlSpace2dKdTreeL::findLastLeaf(
	const DPoint2d & pt) const
{
	return m_top->findLastLeaf(pt);
}


ControlSpace2dKdTreeL::KdElement* ControlSpace2dKdTreeL::findLastLeaf(
	const DPoint2d & pt, DRect & box) const
{
	box = m_box;
	return m_top->findLastLeaf(pt, box);
}

ControlSpace2dKdTreeL::KdElement * ControlSpace2dKdTreeL::findLastLeaf(
	const DPoint2d & pt, DRect & box, int & level) const
{
	box = m_box;
	level = 0;
	return m_top->findLastLeaf(pt, box, level);
}

ControlSpace2dKdTreeL::KdElement * ControlSpace2dKdTreeL::KdElement::findLastLeaf(
	const DPoint2d & pt)
{
	KdElement* kde = this;
	while (!kde->is_leaf) {
		kde = (pt[kde->m_data.split.axis] <= kde->m_data.split.value) ?
			kde->m_data.split.element_0 : kde->m_data.split.element_1;
	}
	return kde;
}

ControlSpace2dKdTreeL::KdElement * ControlSpace2dKdTreeL::KdElement::findLastLeaf(
	const DPoint2d & pt, DRect & box)
{
	KdElement* kde = this;
	while (!kde->is_leaf) {
		bool split_lower = pt[kde->m_data.split.axis] <= kde->m_data.split.value;
		box.splitAndUpdate(kde->m_data.split.axis, kde->m_data.split.value, split_lower);
		assert(box.contains(pt));
		kde = split_lower ? kde->m_data.split.element_0 : kde->m_data.split.element_1;
	}
	return kde;
}

ControlSpace2dKdTreeL::KdElement * ControlSpace2dKdTreeL::KdElement::findLastLeaf(
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

void ControlSpace2dKdTreeL::forEachControlNode(const std::function<void(const ControlNode2d & node)>& fg) const
{
	if (m_top) m_top->forEachControlNodeConst(fg);
}

void ControlSpace2dKdTreeL::forEachControlNode(const std::function<void(ControlNode2d & node)>& fg)
{
	if (m_top) m_top->forEachControlNode(fg);
}

void ControlSpace2dKdTreeL::KdElement::forEachControlNodeConst(
	const std::function<void(const ControlNode2d & node)>& fg) const
{
	if (is_leaf) {
		m_data.split.element_0->forEachControlNode(fg);
		m_data.split.element_1->forEachControlNode(fg);
	}
	else
		fg(m_data.value);
}

void ControlSpace2dKdTreeL::KdElement::forEachControlNode(
	const std::function<void(ControlNode2d&node)>& fg)
{
	if (is_leaf) {
		m_data.split.element_0->forEachControlNode(fg);
		m_data.split.element_1->forEachControlNode(fg);
	}
	else
		fg(m_data.value);
}

int ControlSpace2dKdTreeL::KdElement::adaptToField(
	const std::function<CDM2d(const DPoint2d & pt)> & f,
	const DRect & box, double max_error, int level)
{

	//static size_t counter = 0;
	//if (++counter % 100 == 0) {
	//	LOG4CPLUS_INFO(MeshLog::logger_console, "** counter", counter);
	//	LOG4CPLUS_INFO(MeshLog::logger_console, "** box-diam", box.getDiameter());
	//}

	assert(is_leaf);
	if (level >= MAX_LEVEL) return 0;
	double curr_error = 
		Kd2dRegularErrorEstimator<CDM2d,ControlSpace2dKdTreeL::KdElement, SAMPLING_SIZE>::estimateError(box, this, f);
	if (curr_error <= max_error) return 0;

	Axis split_axis;
	double split_value;
	Kd2dMaxGradientSplitter<CDM2d, SAMPLING_SIZE>::whereToSplit(box, f, split_axis, split_value);

	DataVector<ControlNode2d*> split_vertices(10);
	DRect box0, box1;
	split(split_vertices, box, split_axis, split_value, box0, box1);

	split_vertices.forEach([&f](ControlNode2d* node) {
		node->control_data = f(node->coord);
	});

	return 2 +
		m_data.split.element_0->adaptToField(f, box0, max_error, level + 1) +
		m_data.split.element_1->adaptToField(f, box1, max_error, level + 1);
}

bool ControlSpace2dKdTreeL::KdElement::split(
	DataVector<ControlNode2d*>& split_vertices,
	const DRect & box, Axis split_axis, double split_value,
	DRect & box0, DRect & box1)
{
	assert(is_leaf);
	//assert(valid());

	CDM2d cdm_old = m_data.value.control_data;
	is_leaf = false;

	m_data.split.axis = split_axis;
	m_data.split.value = split_value;
	box.split(m_data.split.axis, m_data.split.value, box0, box1);

	DPoint2d mid0 = box0.getMiddlePoint();
	DPoint2d mid1 = box1.getMiddlePoint();
	m_data.split.element_0 = new KdElement(mid0, cdm_old);
	m_data.split.element_1 = new KdElement(mid1, cdm_old);

	split_vertices.add(&m_data.split.element_0->m_data.value);
	split_vertices.add(&m_data.split.element_1->m_data.value);

	return true;
}

const ControlDataMatrix2d& ControlSpace2dKdTreeL::KdElement::getValue(
	const DPoint2d & pt, const DRect & box) const
{
	assert(box.contains(pt));
	return m_data.value.control_data;
}

const ControlDataMatrix2d& ControlSpace2dKdTreeL::KdElement::getValue(
	const DPoint2d & pt) const
{
	return m_data.value.control_data;
}
void ControlSpace2dKdTreeL::KdElement::clear() {
	if (!is_leaf) {
		m_data.split.element_0->clear();
		delete m_data.split.element_0;
		m_data.split.element_1->clear();
		delete m_data.split.element_1;
		is_leaf = true;
	}
}

int ControlSpace2dKdTreeL::KdElement::getElementCount() const {
	return is_leaf ? 1 :
		1 + m_data.split.element_0->getElementCount() + m_data.split.element_1->getElementCount();
}
int ControlSpace2dKdTreeL::KdElement::getValueCount() const {
	return is_leaf ? 1 :
		m_data.split.element_0->getValueCount() + m_data.split.element_1->getValueCount();
}
int ControlSpace2dKdTreeL::KdElement::getMaxDepth() const {
	return is_leaf ? 1 :
		1 + std::max(m_data.split.element_0->getMaxDepth(), m_data.split.element_1->getMaxDepth());
}
int ControlSpace2dKdTreeL::KdElement::getMaxDepthAndBalance(int & max_balance) const {
	if (is_leaf) return 1;
	int left_depth = m_data.split.element_0->getMaxDepthAndBalance(max_balance);
	int right_depth = m_data.split.element_1->getMaxDepthAndBalance(max_balance);
	int balance = std::abs(left_depth - right_depth);
	if (balance > max_balance) max_balance = balance;
	return 1 + std::max(left_depth, right_depth);
}

/// Refines control space at the given point
bool ControlSpace2dKdTreeL::setMinControl(
	const DPoint2d & pt, const ControlDataMatrix2d & cdm, bool min_value_set)
{
	assert(m_initialized);
	//assert(valid());
	DRect box;
	int level;
	KdElement* kd_leaf = findLastLeaf(pt, box, level);
	return kd_leaf->adaptAndSetControlPoint(ControlNode2d(pt, cdm), base_surface, box, level, min_value_set);
}

bool ControlSpace2dKdTreeL::KdElement::adaptAndSetControlPoint(
	const ControlNode2d& qv, SurfaceConstPtr surface,
	DRect& box, int & level, bool min_value_set)
{
	if (!is_leaf) {
		KdElement* kd_leaf = findLastLeaf(qv.coord, box, level);
		return kd_leaf->adaptAndSetControlPoint(qv, surface, box, level, min_value_set);
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
		double len_x = dmp.transformPStoMS(DVector2d(box.getDX(), 0.0)).length2();
		double len_y = dmp.transformPStoMS(DVector2d(0.0, box.getDY())).length2();
		bool split_needed = (len_x > METRIC_LENGTH_RATIO2 || len_y > METRIC_LENGTH_RATIO2);
		if (split_needed && (level < MAX_LEVEL)) {

			Axis split_axis;
			double split_value;
			Kd2dHalfLongestSplitter<CDM2d>::whereToSplit(box, split_axis, split_value);

			DataVector<ControlNode2d*> split_vertices(10);
			DRect box0, box1;
			if (split(split_vertices, box, split_axis, split_value, box0, box1)) { // split (+ set values for new nodes)

				//split_vertices.forEach([&f](ControlNode2d* node) {
				//	node->control_data = f(node->coord);
				//});

				KdElement* kd_leaf = findLastLeaf(qv.coord, box, level);
				kd_leaf->adaptAndSetControlPoint(qv, surface, box, level, min_value_set);
				return true;
			}
		}
		if (min_value_set) {
			any_changes |= m_data.value.control_data.setMinimum(qv.control_data);
		}
	}

	return any_changes;
}
