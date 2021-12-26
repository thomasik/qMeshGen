#include "ControlSpace3dKdTree.h"
#include "DataHashTable.h"
#include "DataStatistics.h"
#include <vector>
#include <memory>
#include <array>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;

#define WITH_TESTS

#ifdef WITH_TESTS
// ========== TESTS BEGIN =======================
#include <regex>

#include "DLine.h"
#include "DMetric3d.h"
#include "DMatrix.h"

//#include "DataKdTree.h"
//#include "DataOctree.h"

#include "MeshGenerator2d.h"
#include "MeshContainer3d.h"
#include "ControlSpace3dAnalytic.h"
#include "MeshGenerator3d.h"
#include "MeshGenerator3dQuality.h"
#include "MeshBRepXML.h"
#include "MeshViewSet.h"

#include "ControlSpace3dKdTreeL.h"
#include "ControlSpace3dKdTreeV.h"
#include "ControlSpace3dKdTreeLi.h"
#include "ControlSpace3dOctreeL.h"
#include "ControlSpace3dOctreeV.h"
#include "ControlSpace3dOctreeLB.h"
#include "ControlSpace3dOctree.h"

#include "DSphere.h"
#include "SurfaceAnalytic.h"
#include "SurfacePlane.h"
#include "DataStatistics.h"
 
#endif  // WITH_TESTS
// ========== TESTS END =======================

/// Returns descriptive label

int KdTree3dParams::paramEstimateErrorSamplingSize = 4;
int KdTree3dParams::paramSplitGradientSamplingSize = 4;
int KdTree3dParams::paramMaxLevel = 20;
int KdTree3dParams::paramGradationMethod = GRADATION_VIC_NB;
double KdTree3dParams::paramMetricLenRatio = 2.0;

int KdOctree3dParams::paramEstimateErrorSamplingSize = 4;
int KdOctree3dParams::paramMaxLevel = 7;
int KdOctree3dParams::paramMaxBalance = 1;

std::shared_ptr<Kd3dSplitterContField>    KdTree3dParams::m_cf_splitter(
	new Kd3dCFMaxGradientSplitter(KdTree3dParams::paramSplitGradientSamplingSize) );

std::shared_ptr<Kd3dSplitterDiscField>    KdTree3dParams::m_df_splitter(
	new Kd3dDFHalfLongestSplitter() );

std::shared_ptr<Kd3dSplitterSingleSource> KdTree3dParams::m_ss_splitter(
	new Kd3dSSHalfLongestSplitter() );


KdElement::~KdElement() { if (split) delete split; }

bool KdElement::adaptAndSetControlPoint(
	std::shared_ptr<ControlNode3d> qv, DBox & box, int & level, bool min_value_set,
	DataVector<std::shared_ptr<ControlNode3d>>* kvalues)
{
	KdElement* kd_leaf = findLastLeaf(qv->coord, &box, &level);

	// compare new source node with already existing
	// check with calculated value:
	const ControlDataMatrix3d cdm = kd_leaf->getCDM(qv->coord, box);
	double diff = qv->control_data.countDifferenceRR(cdm);
	if (diff < ControlSpace2dAdaptive::param_threshold_diff)
		return false;

	// split may be unnecessary, if qv.control_data is greater than already set for this leaf...
	ControlDataMatrix3d min_cdm = cdm;
	min_cdm.setMinimum(qv->control_data);
	if (cdm.countDifferenceRR(min_cdm) < ControlSpace2dAdaptive::param_threshold_diff)
		return false; // qv different from current metric, but introduces no changes - so skip

	DMetric3d dmp(qv->control_data);
	double mlen_ratio2 = KdTree3dParams::paramMetricLenRatio * KdTree3dParams::paramMetricLenRatio;
	Axis split_axis;
	double split_value;
	DataVector<std::shared_ptr<ControlNode3d>> split_vertices;
	DBox box0, box1;

	while (level < KdTree3dParams::paramMaxLevel) {
		DVector3d len2(
			dmp.transformRStoMS(DVector3d(box.getDX(), 0.0, 0.0)).length2(),
			dmp.transformRStoMS(DVector3d(0.0, box.getDY(), 0.0)).length2(),
			dmp.transformRStoMS(DVector3d(0.0, 0.0, box.getDZ())).length2());
		bool split_needed = (len2.x > mlen_ratio2 || len2.y > mlen_ratio2 || len2.z > mlen_ratio2);

		if (split_needed) {
			KdTree3dParams::m_ss_splitter->whereToSplit(box, qv->coord, len2, split_axis, split_value);
			kd_leaf->splitElementKd(split_vertices, box, split_axis, split_value, &box0, &box1, kvalues);
			kd_leaf = kd_leaf->findLastLeaf(qv->coord, &box, &level);
		}
		else break;
	}

	bool any_changes = false;
	if (min_value_set)
		any_changes |= kd_leaf->setMinimumCDM(qv->control_data);

	return any_changes;
}

int KdElement::getTotalBytes() const {
	int total = (int)sizeof(*this) + getTotalBytesExtra();
	if (split != nullptr) {
		split->elements.forEach([&](const auto& kde) {
			total += kde->getTotalBytes(); });
	}
	return total;
}

int KdElement::getTotalBytesExtra() const {
	return (split != nullptr) ? split->getTotalBytes() : 0;
}

void KdElement::clear() {
	if (split != nullptr) {
		delete split; 
		split = nullptr;
	}
}

KdTreeSplitElement::KdTreeSplitElement(Axis _axis, double _value, 
	std::shared_ptr<KdElement> kde0, 
	std::shared_ptr<KdElement> kde1)
	: KdTreeSplitElement(_axis, _value)
{
	elements.add(kde0);
	elements.add(kde1);
}

std::shared_ptr<KdElement> KdTreeSplitElement::selectChild(const DPoint3d& pt, DBox * box) const
{
	bool split_lower = (pt[axis] <= value);
	if (box != nullptr) {
		box->splitAndUpdate(axis, value, split_lower);
		//assert(box->contains(pt)); // not true for nearestLeaf...
	}
	return elements[split_lower ? 0 : 1];
}

std::shared_ptr<KdElement> KdOctreeSplitElement::selectChild(const DPoint3d& pt, DBox * box) const
{
	bool upper_x = pt.x > middle.x;
	bool upper_y = pt.y > middle.y;
	bool upper_z = pt.z > middle.z;

	int vi = 0;
	if (upper_x) vi |= 1;
	if (upper_y) vi |= 2;
	if (upper_z) vi |= 4;

	if (box != nullptr)
		box->splitOctAndUpdate(middle, vi);

	return elements[vi];
}


KdElement* KdElement::getNearestLeaf(const DPoint3d & pt, DBox * box, int * level) {
	KdElement* element = this;
	while (!element->isLeaf()) {
		element = element->split->selectChild(pt, box).get();
		if (level != nullptr) ++(*level);
	}
	return element;
}

KdElement * KdElement::findLastLeaf(const DPoint3d & pt, DBox * box, int * level)
{
	KdElement* kde = this;
	while (!kde->isLeaf()) {
		kde = kde->split->selectChild(pt, box).get();
		if (level != nullptr) ++(*level);
	}

	return kde;
}

void KdElement::forEachControlNodeConst(const std::function<void(std::shared_ptr<const ControlNode3d> node)>& fg) const
{
	if (isLeaf()) {
		auto cn = getCN();
		if (cn) fg(cn);
	}
	else {
		split->elements.forEach([&](const auto& kde){
			kde->forEachControlNodeConst(fg); });
	}
}

void KdElement::forEachControlNode(const std::function<void(std::shared_ptr<ControlNode3d> node)>& fg)
{
	if (isLeaf()) {
		auto cn = getCN();
		if (cn)	fg(cn);
	}
	else {
		split->elements.forEach([&](auto& kde) {
			kde->forEachControlNode(fg); });
	}
}

double KdElement::estimateError(const DBox & box, const std::function<CDM3d(const DPoint3d&pt)>& f)
{
	double error = 0.0;
	double fx = 1.0 / (double)KdTree3dParams::paramEstimateErrorSamplingSize;
	double dx = box.getDX() * fx;
	double dy = box.getDY() * fx;
	double dz = box.getDZ() * fx;

	DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
	for (; pt.x < box.x1; pt.x += dx) {
		for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
			for (pt.z = box.z0 + 0.5*dz; pt.z < box.z1; pt.z += dz) {
				error = std::max(diffKdValue(f(pt), getCDM(pt, box)),
					error);
			}
		}
	}

	return error;
}

int KdElement::getElementCount() const {
	if (isLeaf()) return 1;

	int result = 1;
	split->elements.forEach([&](const auto& kde) {
		result += kde->getElementCount(); });
	return result;
}

int KdElement::getValueCount() const {
	if (isLeaf()) return 1;

	int result = 0;
	split->elements.forEach([&](const auto& kde) {
		result += kde->getValueCount(); });
	return result;
}

int KdElement::getMaxDepth() const {
	if (isLeaf()) return 1;

	int dmax = 0;
	split->elements.forEach([&](const auto& kde) {
		int d = kde->getMaxDepth();
		if (d > dmax) dmax = d;
	});
	return 1 + dmax;
}

int KdElement::getMaxDepthAndBalance(int & max_balance) const {
	if (isLeaf()) return 1;

	int dmax = -1, dmin = -1;
	split->elements.forEach([&](const auto& kde) {
		if (dmin < 0) {
			dmax = dmin = kde->getMaxDepthAndBalance(max_balance);
		}
		else {
			int d = kde->getMaxDepthAndBalance(max_balance);
			if (d > dmax) dmax = d;
			if (d < dmin) dmin = d;
		}
	});

	int balance = dmax - dmin;
	if (balance > max_balance) max_balance = balance;
	return 1 + dmax;
}

string ControlSpace3dKdTree::getLabel() const
{
	string label = getTreeType();
	if (label[0] == 'k'){ // kdtree
		label += "-e";
		label += to_string(KdTree3dParams::paramEstimateErrorSamplingSize);

		label += "-c";
		label += KdTree3dParams::m_cf_splitter->getLabel();
		label += "-d";
		label += KdTree3dParams::m_df_splitter->getLabel();
		label += "-s";
		label += KdTree3dParams::m_ss_splitter->getLabel();
	}else{
		label += "-e";
		label += to_string(KdOctree3dParams::paramEstimateErrorSamplingSize);
	}

	return label;
}

/// Smoothen variance of metric within the control space (returns true if any change)
bool ControlSpace3dKdTree::smoothenViaLeaves()
{
	assert(m_initialized == 1);
	assert(m_top != nullptr);

	// downtree
	int count = m_top->smoothenMetricAtNodes();
	LOG4CPLUS_INFO(MeshLog::logger_mesh, getLabel() << " smoothing (tree-wise) - " << count << " modifications.");
	return count > 0;
}

bool ControlSpace3dKdTree::smoothenViaVicinityPairs()
{
	assert(m_initialized == 1);
	assert(m_top != nullptr);

	// gather all control nodes with neighboring nodes
	createControlNodesNeighbors();
	int cn_count = getControlNodesCount();
	assert(cn_count > 0);

	class HeapCNPairs : public IndexedObject
	{
	public:
		HeapCNPairs(ControlNode3d* _cn0, ControlNode3d* _cn1, double _max_gradation, int _step = 0) 
			: cn{ _cn0, _cn1 }, key(_max_gradation), step(_step)	{ }
		int compareTo(const HeapCNPairs* hcn) const {
			if (key == hcn->key) return 0;
			else return (key < hcn->key) ? 1 : -1;
		}
		void preDeleteAll() {}
		inline bool valid() const { return (step >= cn[0]->wi) && (step >= cn[1]->wi); }
	public:
		ControlNode3d* cn[2];
		double key;
		int step = 0;
	};

	cout << "#cnodes = " << cn_count << endl;

	const double MAX_GRADATION_EPS = ControlSpace2dAdaptive::param_gradation_ratio * 1.01;
	DataContainer<HeapCNPairs> hcnodes(cn_count);
	hcnodes.setHeapOrder(true);

	forEachControlNode([&](std::shared_ptr<ControlNode3d> cn) {
		cn->wi = 0;
		std::shared_ptr<ControlNode3d> cn1;
		for (auto it = cn->nbs->iterator(); it.valid(); it.moveNext()) {
			cn1 = it.item();
			if (cn1.get() > cn.get()) {
				double gr = DMetric3d::getMetricGradation(
					cn->control_data, cn1->control_data,
					cn->coord - cn1->coord);
				if (gr > MAX_GRADATION_EPS)
					hcnodes.addDataItem(new HeapCNPairs(cn.get(), cn1.get(), gr));
			}
		}
	});

	cout << "#hcnodes = " << hcnodes.countInt() << endl;

	int mod_count = 0;
	int step = 0;
	while( !hcnodes.empty() ){
		auto hcnp = hcnodes.pop(); // remove from heap-container
		if (hcnp->valid()) {
			++step;
			if (step % 1000 == 0) {
				cout << "step = " << step << " (#hcnodes = " << hcnodes.countInt() << endl;
			}
			assert(hcnp->key > MAX_GRADATION_EPS);
			//double gr_before = DMetric3d::getMetricGradation(cn0->control_data, cn1->control_data,
			//	cn1->coord - cn0->coord);
			auto result = DMetric3d::adjustMetricGradation(hcnp->cn[0], hcnp->cn[1]);
			//double gr_after = DMetric3d::getMetricGradation(cn0->control_data, cn1->control_data,
			//	cn1->coord - cn0->coord);

			bool res_arr[] = { result.first, result.second };

			for (int i = 0; i < 2; i++) {
				if (!res_arr[i]) continue;
				// add all neighbors of cni but not |cni-cnj|
				auto cni = hcnp->cn[i];
				cni->wi = step;
				std::shared_ptr<ControlNode3d> cn_nb;
				for (auto it = cni->nbs->iterator(); it.valid(); it.moveNext()) {
					cn_nb = it.item();
					if (cn_nb.get() == hcnp->cn[1 - i]) continue;
					double gr = DMetric3d::getMetricGradation(
						cni->control_data, cn_nb->control_data,
						cni->coord - cn_nb->coord);
					if (gr > MAX_GRADATION_EPS) {
						hcnodes.addDataItem(new HeapCNPairs(cni, cn_nb.get(), gr, step));
					}
				}
			}
			if (result.first || result.second) mod_count++;
		}
		delete hcnp;
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh, getLabel() << " smoothing (vicinity-wise-pairs) - " << mod_count << " modifications.");

	double max_gr0 = 0.0;
	forEachControlNode([&](std::shared_ptr<ControlNode3d> cn) {
		std::shared_ptr<ControlNode3d> cn1;
		for (auto it = cn->nbs->iterator(); it.valid(); it.moveNext()) {
			cn1 = it.item();
			if (cn1.get() > cn.get()) {
				double gr = DMetric3d::getMetricGradation(
					cn->control_data, cn1->control_data,
					cn->coord - cn1->coord);
				if (gr > max_gr0) max_gr0 = gr;
			}
		}
	});
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "max_gr, level 0 -> " << max_gr0);

	double max_gr1 = 0.0;
	forEachControlNode([&](std::shared_ptr<ControlNode3d> cn) {
		std::shared_ptr<ControlNode3d> cn1;
		std::shared_ptr<ControlNode3d> cn11;
		for (auto it = cn->nbs->iterator(); it.valid(); it.moveNext()) {
			cn1 = it.item();
			for (auto itj = cn1->nbs->iterator(); itj.valid(); itj.moveNext()) {
				cn11 = itj.item();
				if (cn11.get() > cn.get()) {
					double gr = DMetric3d::getMetricGradation(
						cn->control_data, cn11->control_data,
						cn->coord - cn11->coord);
					if (gr > max_gr1) max_gr1 = gr;
				}
			}
		}
	});

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "max_gr, level 1 -> " << max_gr1);

	clearControlNodesNeighbors();
	return mod_count > 0;
}

bool ControlSpace3dKdTree::smoothenViaVicinityNb()
{
	assert(m_initialized == 1);
	assert(m_top != nullptr);

	// gather all control nodes with neighboring nodes
	createControlNodesNeighbors();
	int cn_count = getControlNodesCount();
	assert(cn_count > 0);

	class HeapCNode : public IndexedObject
	{
	public:
		HeapCNode(ControlNode3d* _cn0, double _max_gradation, int _step = 0)
			: cn0(_cn0), key(_max_gradation), step(_step)	{ }
		int compareTo(const HeapCNode* hcn) const {
			if (key == hcn->key) return 0;
			else return (key < hcn->key) ? 1 : -1;
		}
		void preDeleteAll() {}
		inline bool valid() const { return (step >= cn0->wi); }
	public:
		ControlNode3d* cn0;
		double key;
		int step = 0;
	};

	// invalidate all cn adjacent to modified ones
	DataHashTable<ControlNode3d*> cn_to_invalidate(2*cn_count, nullptr);
	forEachControlNode([&](std::shared_ptr<ControlNode3d> cn) {
		cn->wi = 0;
		if (cn->gradationUnknown()) {
			std::shared_ptr<ControlNode3d> cn1;
			for (auto it = cn->nbs->iterator(); it.valid(); it.moveNext()) {
				cn1 = it.item();
				if (!cn1->gradationUnknown()) cn_to_invalidate.insert(cn1.get());
			}
		}
	});
	for (auto cni : cn_to_invalidate.values())
		cni->setGradationUnknown();

	const double MAX_GRADATION_EPS = ControlSpace2dAdaptive::param_gradation_ratio * 1.01;
	DataContainer<HeapCNode> hcnodes(cn_count);
	hcnodes.setHeapOrder(true);

	DataHashTableKeyValue< t_cn_pair, t_gr_step_pair > gr_cn_pair(6*cn_count, std::make_pair(nullptr, nullptr));
	forEachControlNode([&](std::shared_ptr<ControlNode3d> cn) {
		if (cn->gradationUnknown()) {
			cn->calculateMaxGradationNb(gr_cn_pair, 0);
			hcnodes.addDataItem(new HeapCNode(cn.get(), cn->max_gradation_ratio));
		}
	});

	int mod_count = 0;
	int step = 0;
	while (!hcnodes.empty()) {
		auto hcn = hcnodes.pop(); // remove from heap-container
		if (hcn->valid()) {
			++step;
			assert(hcn->key > MAX_GRADATION_EPS);

			bool changed = DMetric3d::adjustMetricGradationNb(hcn->cn0);
			if (changed) mod_count++;

			std::shared_ptr<ControlNode3d> cn1;
			for (auto it = hcn->cn0->nbs->iterator(); it.valid(); it.moveNext()) {
				cn1 = it.item();
				if (cn1->wi == -1) { // changed
					cn1->wi = step;
					cn1->calculateMaxGradationNb(gr_cn_pair, step);
					hcnodes.addDataItem(new HeapCNode(cn1.get(), cn1->max_gradation_ratio));
				}
			}

			delete hcn;
		}
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh, getLabel() << " smoothing (vicinity-wise-nb) - " << mod_count << " modifications.");

	clearControlNodesNeighbors();
	return mod_count > 0;
}

void ControlSpace3dKdTree::createControlNodesNeighbors()
{
	//auto cn_count = getControlNodesCount();
	forEachControlNode([](std::shared_ptr<ControlNode3d> cn) {
		//cn->nbs.reset();
		cn->nbs = std::make_unique<DataCompoundList<std::shared_ptr<ControlNode3d>>>(6); // max number of neighbours
	});
	m_top->gatherControlNodesNeighbors();
}

void ControlSpace3dKdTree::clearControlNodesNeighbors()
{
	DataStatistics ds;
	forEachControlNode([&ds](std::shared_ptr<ControlNode3d> cn) {
		ds.add(cn->nbs->countInt());
		cn->nbs.reset();
	});
	if (ds.calculate()) {
		LOG4CPLUS_INFO(MeshLog::logger_console, "nbs-min: " << ds.minimum());
		LOG4CPLUS_INFO(MeshLog::logger_console, "nbs-max: " << ds.maximum());
		LOG4CPLUS_INFO(MeshLog::logger_console, "nbs-ave: " << ds.average());
	}
}

bool ControlSpace3dKdTree::smoothenViaVicinityNbMultipass()
{
	assert(m_initialized == 1);
	assert(m_top != nullptr);

//	dumpCNodes("kd-dump-00.txt");

	// gather all control nodes with neighboring nodes
	createControlNodesNeighbors();
	int cn_count = getControlNodesCount();
	assert(cn_count > 0);

//	dumpCNodes("kd-dump-01.txt");

	const double MAX_GRADATION_EPS = ControlSpace2dAdaptive::param_gradation_ratio * 1.05;

	DataVector<std::shared_ptr<ControlNode3d>> current_nodes(cn_count);
	forEachControlNode([MAX_GRADATION_EPS, &current_nodes](std::shared_ptr<ControlNode3d> cn) {
		if (cn->gradationUnknown() || (cn->max_gradation_ratio > MAX_GRADATION_EPS)) {
			current_nodes.add(cn);
		}
	});
	if (current_nodes.empty()) {
		clearControlNodesNeighbors();
		return false; // no change
	}

	const int MAX_STEPS = 10;
	int mod_count = 0;

	std::vector<double> max_grad_ranges = { 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 10.0, GRADATION_UNKNOWN-1.0 };

	{
		ostringstream stat_line;
		stat_line << "STEP\t0.0";
		for (auto d : max_grad_ranges) stat_line << "\t" << d;
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, stat_line.str());
	}

	for (int si = 0; si < MAX_STEPS; si++) {
		//* remove cnodes already set within limit
		//current_nodes.removeIf( [si](ControlNode3d * const & cn) { return cn->wi < si; });
		if (current_nodes.empty()) break;

		for (int i = 0; i < current_nodes.countInt(); i++) {
			auto cn0 = current_nodes[i];
			bool changed = DMetric3d::adjustMetricGradationNb(cn0.get());
			if (changed) mod_count++;
		}

		//string fname = "kd-dump-02-" + to_string(si) + ".txt";
		//dumpCNodes(fname.c_str());

		DataStatistics ds;
		DataHashTableKeyValue< t_cn_pair, double> gr_cn_pair(6 * cn_count, std::make_pair(nullptr, nullptr));
		current_nodes.clear();
		forEachControlNode([&](std::shared_ptr<ControlNode3d> cn) {
			double mg = cn->calculateMaxGradationNb(gr_cn_pair);
			ds.add(mg);
			if (mg > MAX_GRADATION_EPS) {
				cn->wi = si + 1;
				current_nodes.add(cn);
			}
		});

		if (ds.calculate()) {
			auto vec_d = ds.getDataCountInRanges(max_grad_ranges);

			{
				ostringstream stat_line;
				stat_line << si;
				for (auto d : vec_d) stat_line << "\t" << d;
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, stat_line.str());

			}
		}

	}

	//dumpCNodes("kd-dump-03.txt");
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, getLabel() << " smoothing (vicinity-wise-nb-multipass) - " << mod_count << " modifications.");

	clearControlNodesNeighbors();
	return mod_count > 0;
}

ControlDataMatrix3d ControlSpace3dKdTree::getMetricAtPoint(const DPoint3d & pt) const {
	DBox box = m_box;
	return findLastLeaf(pt, &box)->getCDM(pt, box);
}

void ControlSpace3dKdTree::clear() {
	m_top.reset();
	m_values.reset();
}

/// Refines control space at the given point
int ControlSpace3dKdTree::adaptToField(const std::function<CDM3d(const DPoint3d&pt)>& f, double max_error)
{
	DataVector<std::shared_ptr<ControlNode3d>> vertices(10);
	createTopElement(vertices);
	vertices.forEach([&f](std::shared_ptr<ControlNode3d> cn) { cn->control_data = f(cn->coord); });

	m_initialized = 1;
	return 1 + m_top->adaptToField(f, m_box, max_error, 0, m_values.get());
}

bool ControlSpace3dKdTree::setMinControl(const DPoint3d & pt, const ControlDataMatrix3d & cdm, bool min_value_set)
{
	assert(m_initialized);
	assert(!isCompacted());
	//assert(valid());
	DBox box;
	int level;
	KdElement* kd_leaf = findLastLeaf(pt, &box, &level);
	return kd_leaf->adaptAndSetControlPoint(
		std::make_shared<ControlNode3d>(pt, cdm), box, level, min_value_set, m_values.get());
}

/// Initializes the control space with given global metric
void ControlSpace3dKdTree::setGlobalMetric(const ControlDataMatrix3d & cdm)
{
	DataVector<std::shared_ptr<ControlNode3d>> vertices(10);
	createTopElement(vertices);
	vertices.forEach([&cdm](std::shared_ptr<ControlNode3d> cn) { cn->control_data = cdm; });

	m_initialized = 1;
}

int ControlSpace3dKdTree::getValueCount() const {
	if (m_values != nullptr)
		return m_values->countInt();
	else
		return m_top ? m_top->getValueCount() : 0; 
}

int ControlSpace3dKdTree::getMaxBalance() const {
	int max_balance = 0;
	if (m_top) m_top->getMaxDepthAndBalance(max_balance);
	return max_balance;
}

int ControlSpace3dKdTree::getTotalBytes() const 
{
	int total = (int)sizeof(*this);
	if (m_values != nullptr)
		total += (int)sizeof(ControlNode3d) * m_values->countInt() + sizeof(*m_values);
	if (m_top != nullptr)
		total += m_top->getTotalBytes();
	return total;
}

KdElement * ControlSpace3dKdTree::findLastLeaf(const DPoint3d & pt, DBox * box, int * level) const
{
	if (box) *box = m_box;
	if (level) *level = 0;
	assert(m_top != nullptr);
	return m_top->findLastLeaf(pt, box, level);
}

void ControlSpace3dKdTree::forEachControlNode(
	const std::function<void(std::shared_ptr<const ControlNode3d> node)>& fg) const
{
	if(m_values != nullptr)
		m_values->forEach(fg);
	else if (m_top) 
		m_top->forEachControlNodeConst(fg);
}

void ControlSpace3dKdTree::forEachControlNode(
	const std::function<void(std::shared_ptr<ControlNode3d> node)>& fg)
{
	assert(!isCompacted());

	if (m_values != nullptr)
		m_values->forEach(fg);
	else if (m_top) m_top->forEachControlNode(fg);
}

void ControlSpace3dKdTree::dumpCNodes(const char * fname) const
{
	ofstream ofs(fname);

	forEachControlNode([&ofs](std::shared_ptr<const ControlNode3d> cn) {
		ofs << cn->coord << "\t" << cn->control_data << endl;
	});
}

bool KdOctElement::adaptAndSetControlPoint(
	std::shared_ptr<ControlNode3d> qv, DBox & box, int & level, bool min_value_set,
	DataVector<std::shared_ptr<ControlNode3d>>* kvalues)
{
	KdElement* kd_leaf = findLastLeaf(qv->coord, &box, &level);

	// compare new source node with already existing
	// check with calculated value:
	const ControlDataMatrix3d cdm = kd_leaf->getCDM(qv->coord, box);
	double diff = qv->control_data.countDifferenceRR(cdm);
	if (diff < ControlSpace2dAdaptive::param_threshold_diff)
		return false;

	// split may be unnecessary, if qv.control_data is greater than already set for this leaf...
	ControlDataMatrix3d min_cdm = cdm;
	min_cdm.setMinimum(qv->control_data);
	if (cdm.countDifferenceRR(min_cdm) < ControlSpace2dAdaptive::param_threshold_diff)
		return false; // qv different from current metric, but introduces no changes - so skip

	DMetric3d dmp(qv->control_data);
	double mlen_ratio2 = KdTree3dParams::paramMetricLenRatio * KdTree3dParams::paramMetricLenRatio;
	DataVector<std::shared_ptr<ControlNode3d>> split_vertices;
	DBox boxes[8];

	while (level < KdOctree3dParams::paramMaxLevel) {
		DVector3d len2(
			dmp.transformRStoMS(DVector3d(box.getDX(), 0.0, 0.0)).length2(),
			dmp.transformRStoMS(DVector3d(0.0, box.getDY(), 0.0)).length2(),
			dmp.transformRStoMS(DVector3d(0.0, 0.0, box.getDZ())).length2());
		bool split_needed = (len2.x > mlen_ratio2 || len2.y > mlen_ratio2 || len2.z > mlen_ratio2);

		if (split_needed) {
			kd_leaf->splitElementOct(split_vertices, box, boxes, kvalues); // split (+ set values for new nodes)
			kd_leaf = kd_leaf->findLastLeaf(qv->coord, &box, &level);
		}
		else break;
	}

	bool any_changes = false;
	if (min_value_set)
		any_changes |= kd_leaf->setMinimumCDM(qv->control_data);

	return any_changes;
}

#ifdef WITH_TESTS
//============================= TESTS BEGIN ====================================

const int TEST_TREE_RPT = 1;
//const int TEST_TREE_RPT = 10;
const int TEST_ACCESS_RPT = 1;
//	const int TEST_ACCESS_RPT = 1000;

// labbe04 isotropic size-function
CDM3d f_labbe_iso(const DPoint3d& pt) {
	double y = (pt.y + 1.0) * 4.5; // [-1,1] -> [0,9]
	double v;
	if (y < 2) v = 1.0 - (19 * y) / 40;
	else if (y < 4.5) v = std::pow(20, 0.4*y - 1.8);
	else if (y < 7) v = std::pow(5, 1.8 - 0.4*y);
	else v = 0.2 + 0.8*std::pow(0.5*y - 3.5, 4);
	v /= 4.5; //
	if (v < SMALL_NUMBER) v = SMALL_NUMBER;
	return CDM3d(v);
}
// labbe04 anisotropic size-function (extended accordingly for z-axis)
CDM3d f_labbe_aniso(const DPoint3d& pt) {
	double x = (pt.y + 1.0) * 3.5; // [-1,1] -> [0,7]
	double y = (pt.y + 1.0) * 4.5; // [-1,1] -> [0,9]
	double z = (pt.z + 1.0) * 5.5; // [-1,1] -> [0,11]
	double vx;
	if (x < 2) vx = 1.0 - (19 * x) / 40;
	else if (x < 3.5) vx = std::pow(20, 2.0*x / 3.0 - 7.0 / 3.0);
	else if (x < 5) vx = std::pow(5, 7.0 / 3.0 - 2.0*x / 3.0);
	else vx = 0.2 + 0.8*std::pow(0.5*x - 2.5, 4);
	vx /= 3.5;
	if (vx < SMALL_NUMBER) vx = SMALL_NUMBER;
	double vy;
	if (y < 2) vy = 1.0 - (19 * y) / 40;
	else if (y < 4.5) vy = std::pow(20, 0.4*y - 1.8);
	else if (y < 7) vy = std::pow(5, 1.8 - 0.4*y);
	else vy = 0.2 + 0.8*std::pow(0.5*y - 3.5, 4);
	vy /= 4.5;
	if (vy < SMALL_NUMBER) vy = SMALL_NUMBER;
	double vz;
	if (z < 2) vz = 1.0 - (19 * z) / 40;
	else if (z < 5.5) vz = std::pow(20, 2.0*z / 7.0 - 11.0 / 7.0);
	else if (z < 9) vz = std::pow(5, 11.0 / 7.0 - 2.0*z / 7.0);
	else vz = 0.2 + 0.8*std::pow(0.5*z - 4.5, 4);
	vz /= 5.5;
	if (vz < SMALL_NUMBER) vz = SMALL_NUMBER;
	return CDM3d(vx, vy, vz, 0.0, 0.0, 0.0);
}

// a gauss function around the (0,0,0) with 1-sigma equal to 0.1
CDM3d f_gauss(const DPoint3d& pt) {
	double v = std::exp(-(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z) / (2 * 0.1));
	if (v < SMALL_NUMBER) v = SMALL_NUMBER;
	return CDM3d(v);
}

// proportional to (squared) distance from a line
CDM3d f_line_dist(const DPoint3d& pt) {
	static const DPoint3d  line_pt(-1.0, -1.0, -0.5);
	static const DVector3d line_vt = DVector3d(1.5, 2.0, 0.5).normalized();
	double v = DLine3d::distanceToPoint2(line_pt, line_vt, pt, true);
	if (v < SMALL_NUMBER) v = SMALL_NUMBER;
	return CDM3d(v);
}

// proportional to distance from three points
CDM3d f_points_dist(const DPoint3d& pt) {
	DPoint3d a(0.5, 0.5, 0.5);
	DPoint3d b(0.5, 0.5, -0.5);
	DPoint3d c(0.0, 0.0, 0.0);

	double v = (pt.distance(a) + pt.distance(b) + pt.distance(c)) / 3.0;
	if (v < SMALL_NUMBER) v = SMALL_NUMBER;
	return CDM3d(v);
}

CDM3d(*test_functions[5])(const DPoint3d& pt) = {
	f_labbe_iso, f_labbe_aniso, f_gauss, f_line_dist, f_points_dist
};

const char* test_labels[] = {
	"labbe04-iso", "labbe04-aniso", "f-gauss", "f-line-dist2", "f-3pts-dist"
};

struct DiscreteSources
{
public:
	DiscreteSources(const DBox &_box) : src_tree(_box), 
		source_h_ratio_all(1.0), source_h_ratio_normal(1.0), box(_box) {}
public:
	bool matchPredefined(const string & fname) {
		std::regex sphere_regex_sphere("f-sphere-r(\\d+\\.\\d+)");
		std::regex sphere_regex_spheres("f-spheres-r(\\d+\\.\\d+)-r(\\d+\\.\\d+)-n(\\d+)");
		std::regex sphere_regex_spheres_symc("f-symc-spheres-r(\\d+\\.\\d+)-n(\\d+)");
		std::regex sphere_regex_spheres_syml("f-syml-spheres-r(\\d+\\.\\d+)-n(\\d+)");
		std::regex sphere_regex_ellipsoid("f-ellipsoid-r(\\d+\\.\\d+)-r(\\d+\\.\\d+)-r(\\d+\\.\\d+)");
		std::regex sphere_regex_ellipsoids("f-ellipsoids-r(\\d+\\.\\d+)-r(\\d+\\.\\d+)-n(\\d+)");
		std::regex f_regex_sin("f-sin");
		std::regex f_regex_planes("f-planes-a(\\d+\\.\\d+)");
		std::smatch src_match;
		if (std::regex_match(fname, src_match, sphere_regex_sphere)) {
			assert(src_match.size() == 2);
			analytic = true;
			std::ssub_match sub_match = src_match[1];
			istringstream iss(sub_match.str());
			DSphere s(DPoint3d::zero);
			iss >> s.m_radius;
			setSphereSources(s, 0.1);
			DataVector<DSphere> spheres(1);
			spheres.add(s);
			prepareSphereFunc(spheres, 0.1);
			cout << "Created " << sources_list.countInt() << " sources for sphere with r=" << s.m_radius
				<< " and dratio=" << density_ratio << endl;

			//showSources();

			return true;
		}
		if (std::regex_match(fname, src_match, sphere_regex_spheres)) {
			assert(src_match.size() == 4);
			analytic = true;
			std::string data_str = src_match[1].str() + " " + src_match[2].str() + " " + src_match[3].str();
			double r0, r1; // min-max
			int n;
			istringstream iss(data_str);
			iss >> r0 >> r1 >> n;
			RandomGen rg(1);
			assert(box.getDX() > 2 * r1);
			assert(box.getDY() > 2 * r1);
			assert(box.getDZ() > 2 * r1);

			DSphere s(DPoint3d::zero);
			DataVector<DSphere> spheres(n);

			for (int i = 0; i < n; i++) {
				//double rx = rg.doub(r0, r1);
				double t = (n > 1) ? ((double)i / (n - 1)) : 0.5;
				s.m_radius = (1 - t) * r0 + t * r1;
				s.m_center = box.getRandomPoint(s.m_radius);
				spheres.add(s);
				setSphereSources(s);
			}
			prepareSphereFunc(spheres);

			cout << "Created total " << sources_list.countInt() << " sources for " << n
				<< " spheres with r: ";
			spheres.forEach([&](const DSphere& sp) { cout << sp.m_radius << ",";  });
			cout << " and dratio=" << density_ratio << endl;
			return true;
		}
		if (std::regex_match(fname, src_match, sphere_regex_spheres_syml)) {
			assert(src_match.size() == 3);
			analytic = true;
			std::string data_str = src_match[1].str() + " " + src_match[2].str();
			double r;
			int n;
			istringstream iss(data_str);
			iss >> r >> n;
			assert(box.getDX() > 2 * r);
			assert(box.getDY() > 2 * r);
			assert(box.getDZ() > 2 * r);

			DSphere s(DPoint3d::zero);
			DataVector<DSphere> spheres(n);

			double rlen = 1.0 + (n - 1)*(5.0 / 6.0);
			s.m_radius = r / rlen;
			s.m_center.x = -r - s.m_radius;
			for (int i = 0; i < n; i++) {
				s.m_center.x += 2 * s.m_radius;
				spheres.add(s);
				setSphereSources(s);
			}
			prepareSphereFunc(spheres);

			cout << "Created total " << sources_list.countInt() << " sources for " << n
				<< " spheres with r: ";
			spheres.forEach([&](const DSphere& sp) { cout << sp.m_radius << ",";  });
			cout << " and dratio=" << density_ratio << endl;
			return true;
		}		
		if (std::regex_match(fname, src_match, sphere_regex_spheres_symc)) {
			assert(src_match.size() == 3);
			analytic = true;
			std::string data_str = src_match[1].str() + " " + src_match[2].str();
			double r;
			int n;
			istringstream iss(data_str);
			iss >> r >> n;
			assert(box.getDX() > 2 * r);
			assert(box.getDY() > 2 * r);
			assert(box.getDZ() > 2 * r);

			DSphere s(DPoint3d::zero);
			DataVector<DSphere> spheres(n);

			s.m_radius = r  * (11.0/24.0);
			double dr = s.m_radius * (5.0 / 6.0);
			for (int i = 0; i < n; i++) {
				double angle = (2 * PI * i) / n;
				s.m_center.x = dr * cos(angle);
				s.m_center.y = dr * sin(angle);
				spheres.add(s);
				setSphereSources(s);
			}
			prepareSphereFunc(spheres);

			cout << "Created total " << sources_list.countInt() << " sources for " << n
				<< " spheres with r: ";
			spheres.forEach([&](const DSphere& sp) { cout << sp.m_radius << ",";  });
			cout << " and dratio=" << density_ratio << endl;
			return true;
		}
		if (std::regex_match(fname, src_match, sphere_regex_ellipsoid)) {
			assert(src_match.size() == 4);
			analytic = true;
			DEllipsoid ell;
			DataVector<DEllipsoid> ellipsoids(1);
			for (int i = 0; i < 3; i++) {
				std::ssub_match sub_match = src_match[i + 1];
				istringstream iss(sub_match.str());
				iss >> ell.getRadius(i);
			}
			ellipsoids.add(ell);
			setEllipsoidSources(ell);
			prepareEllipsoidFunc(ellipsoids);
			cout << "Created " << sources_list.countInt() << " sources for ellipse with r=["
				<< ell.getRadius(0) << "," << ell.getRadius(1) << "," << ell.getRadius(2)
				<< " and dratio=" << density_ratio << endl;
			return true;
		}
		if (std::regex_match(fname, src_match, sphere_regex_ellipsoids)) {
			assert(src_match.size() == 4);
			analytic = true;
			std::string data_str = src_match[1].str() + " " + src_match[2].str() + " " + src_match[3].str();
			double r0, r1; // min-max
			int n;
			istringstream iss(data_str);
			iss >> r0 >> r1 >> n;
			RandomGen rg(1);
			assert(box.getDX() > 2 * r1);
			assert(box.getDY() > 2 * r1);
			assert(box.getDZ() > 2 * r1);

			DEllipsoid e(DPoint3d::zero);
			DataVector<DEllipsoid> ellipsoids(n);

			for (int i = 0; i < n; i++) {
				//double rx = rg.doub(r0, r1);
				e.getRadius(0) = rg.doub(r0, r1);
				e.getRadius(1) = rg.doub(r0, r1);
				e.getRadius(2) = rg.doub(r0, r1);
				e.m_center = box.getRandomPoint(e.getRadius(0), e.getRadius(1), e.getRadius(2));
				ellipsoids.add(e);
				setEllipsoidSources(e);
			}
			prepareEllipsoidFunc(ellipsoids);

			cout << "Created total " << sources_list.countInt() << " sources for " << n
				<< " ellipsoids with r: ";
			ellipsoids.forEach([&](const DEllipsoid& el) {
				cout << el.getRadius(0) << "x" << el.getRadius(1) << "x" << el.getRadius(2) << ",";  });
			cout << " and dratio=" << density_ratio << endl;
			return true;
		}
		if (std::regex_match(fname, src_match, f_regex_sin)) {
			assert(src_match.size() == 1);
			analytic = true;
			SurfaceAnalytic fsin("u", "sin(PI*u)", "v");
			setFSources(&fsin, 0.1);


			auto plane_sin = std::make_shared<SurfacePlane>(
				DPoint3d::zero, DVector3d::v_ox, DVector3d::v_oz);

			prepareSinDiscreteSourcesFunc(plane_sin);

			LOG4CPLUS_INFO(MeshLog::logger_mesh, "Created " << sources_list.countInt() << " sources for f-sin with dratio=" << density_ratio);

			//showSources();
			return true;
		}
		if (std::regex_match(fname, src_match, f_regex_planes)) {
			assert(src_match.size() == 2);
			analytic = true;
			std::ssub_match sub_match = src_match[1];
			istringstream iss(sub_match.str());
			double a;
			iss >> a;
			SurfaceAnalytic fplane0("-0.5", "u", "v");
			setFSources(&fplane0, 0.1);

			auto cdm0 = sources_list.last()->control_data;

			string y_formula = "-0.2";
			if (a > 0.0) {
				y_formula += "+u*";
				y_formula += to_string(a);
			}
			else if(a < 0.0) {
				y_formula += to_string(a) + "*u";
			}
			SurfaceAnalytic fplane1("u", y_formula.c_str(), "v");
			setFSources(&fplane1, 0.1);
			auto cdm1 = sources_list.last()->control_data;

			//prepareDiscreteSourcesFunc();

			auto plane0 = std::make_shared<SurfacePlane>(
				fplane0.getPoint(DPoint2d::zero), fplane0.getNormalVector(DPoint2d::zero));
			auto plane1 = std::make_shared<SurfacePlane>(
				fplane1.getPoint(DPoint2d::zero), fplane1.getNormalVector(DPoint2d::zero));

			preparePlaneSourcesFunc(plane0, cdm0, plane1, cdm1);

			LOG4CPLUS_INFO(MeshLog::logger_mesh, 
				"Created " << sources_list.countInt() << " sources for f-plane with a=" << a
				<< " and dratio=" << density_ratio);

			showSources();

			return true;
		}
		return false;
	}
	bool readFromFile(const char* fname) {
		analytic = false;

		string fpath = main_dir + fname;

		ifstream ifs(fpath);
		if (!ifs) {
			cout << "Couldn't open file: " << fpath << endl;
			return false;
		}
		DPoint3d pt;
		ControlDataMatrix3d cdm;
		while (ifs) {
			if (pt.readSimple(ifs) && cdm.readSimple(ifs)) {
				sources_list.add(std::make_shared<ControlNode3d>(pt, cdm));
			}
		}
		prepareFullDiscreteSourcesFunc();
		cout << "Read " << sources_list.countInt() << " sources from file: " << fpath << endl;
		return true;
	}
	void setEllipsoidSources(const DEllipsoid& e, double cr = -1.0) {
		double r_min = std::min(e.getRadius(0), std::min(e.getRadius(1), e.getRadius(2)));
		double h_min = ControlSpace2dAdaptive::param_curvature_ratio * ( (cr > 0.0) ? cr : r_min);
		double h_step = h_min * DiscreteSources::density_ratio;

		const double TWO_PI = 2 * PI;
		double arc_len = TWO_PI * r_min; // full circle
		double arc_ratio = arc_len / h_step;
		int nt = 1 + (int)arc_ratio;

		const DPoint3d& emid = e.getCenter();
		string ex = to_string(emid.x) + "+" + to_string(e.getRadius(0)) + "*cos(u)*sin(v)";
		string ey = to_string(emid.y) + "+" + to_string(e.getRadius(1)) + "*sin(u)*sin(v)";
		string ez = to_string(emid.z) + "+" + to_string(e.getRadius(2)) + "*cos(v)";
		SurfaceConstPtr surf_ellipsoid(new SurfaceAnalytic (ex.c_str(), ey.c_str(), ez.c_str()));

		// -> big circle
		double dalpha = PI / nt;
		CDM3d cdm;
		DPoint2d uv(0.0, PI / 2);
		double p = box.getDiameter();

		DVector3d dv1, dv2;
		double hx = h_min * source_h_ratio_all;
		double d[] = { hx * source_h_ratio_normal, hx, hx };

		for (int j = 0; j < nt; j++, uv.x += dalpha) {
			DPoint3d ept = surf_ellipsoid->getPoint(uv);
			if (cr > 0.0) {
				auto dv = surf_ellipsoid->getNormalVector(uv);
				dv.orthonormalVectors(dv1, dv2);
				cdm = ControlDataMatrix3d(dv, dv1, dv2, d);
			}else{
				if (!ControlSpace3dAdaptive::calculateControlFromSurfaceCurvature(
						ept, surf_ellipsoid, p, cdm, &uv))
					continue;
			}
			sources_list.add(std::make_shared<ControlNode3d>(ept, cdm));
			sources_list.add(std::make_shared<ControlNode3d>(emid - (ept - emid), cdm));
		}

		// -> smaller ones (smaller circles, and lower/higher)
		DPoint2d uv_l(0.0, PI / 2);
		DPoint2d uv_h(0.0, PI / 2);
		while (uv_l.y > 0.5*dalpha) {
			double sratio = sin(uv_l.y) * arc_ratio;
			int snt = 1 + (int)sratio;
			double sdalpha = PI / snt;
			uv_l.x = uv_h.x = 0.0;
			for (int j = 0; j < snt; j++, uv_l.x += sdalpha, uv_h.x += sdalpha) {
				DPoint3d ept_l = surf_ellipsoid->getPoint(uv_l);
				if (cr > 0.0) {
					auto dv = surf_ellipsoid->getNormalVector(uv_l);
					dv.orthonormalVectors(dv1, dv2);
					cdm = ControlDataMatrix3d(dv, dv1, dv2, d);
				}
				else {
					if (!ControlSpace3dAdaptive::calculateControlFromSurfaceCurvature(
						ept_l, surf_ellipsoid, p, cdm, &uv_l))
						continue;
				}
				sources_list.add(std::make_shared<ControlNode3d>(ept_l, cdm));
				sources_list.add(std::make_shared<ControlNode3d>(emid - (ept_l - emid), cdm));

				DPoint3d ept_h = surf_ellipsoid->getPoint(uv_h);
				if (cr > 0.0) {
					auto dv = surf_ellipsoid->getNormalVector(uv_h);
					dv.orthonormalVectors(dv1, dv2);
					cdm = ControlDataMatrix3d(dv, dv1, dv2, d);
				}
				else {
					if (!ControlSpace3dAdaptive::calculateControlFromSurfaceCurvature(
						ept_h, surf_ellipsoid, p, cdm, &uv_h))
						continue;
				}
				sources_list.add(std::make_shared<ControlNode3d>(ept_h, cdm));
				sources_list.add(std::make_shared<ControlNode3d>(emid - (ept_h - emid), cdm));
			}
			uv_l.y -= dalpha;
			uv_h.y += dalpha;
		}

		// -> finish
	}
	void setSphereSources(const DSphere& s, double cr = -1.0) {
		double r = s.getRadius();
		if (cr <= 0.0) cr = r;
		const DPoint3d& pt_mid = s.m_center;
		double h = ControlSpace2dAdaptive::param_curvature_ratio * cr;
		//CDM3d cdm(h);

		double h_min = h;

		double h_step = h_min * DiscreteSources::density_ratio;

		const double HALF_PI = 0.5 * PI;
		double arc_len = HALF_PI * r; // one-fourth of the full circle
		double arc_ratio = arc_len / h_step;
		int nt = 1 + (int)arc_ratio;

		DataVector<DVector3d> nvs(nt);

		// -> big circle
		double dalpha = HALF_PI / nt;
		for (int j = 0; j < nt; j++) {
			double sa = sin(j * dalpha);
			double ca = cos(j * dalpha);
			nvs.add(DVector3d::v_ox *  sa + DVector3d::v_oy *  ca);
			nvs.add(DVector3d::v_ox * -sa + DVector3d::v_oy *  ca);
			nvs.add(DVector3d::v_ox *  sa + DVector3d::v_oy * -ca);
			nvs.add(DVector3d::v_ox * -sa + DVector3d::v_oy * -ca);
		}
		// -> smaller ones (smaller circles, and lower/higher)
		for (int i = 1; i < nt; i++) {
			double cos_i_dalpha = cos(i * dalpha);
			double sratio = cos_i_dalpha * arc_ratio;
			int snt = 1 + (int)sratio;
			double sdalpha = HALF_PI / snt;
			for (int j = 0; j < snt; j++) {
				double sa = sin(j * sdalpha) * cos_i_dalpha;
				double ca = cos(j * sdalpha) * cos_i_dalpha;
				DVector3d snvz = DVector3d::v_oz * sin(i * dalpha);
				nvs.add(DVector3d::v_ox *  sa + DVector3d::v_oy *  ca + snvz);
				nvs.add(DVector3d::v_ox * -sa + DVector3d::v_oy *  ca + snvz);
				nvs.add(DVector3d::v_ox *  sa + DVector3d::v_oy * -ca + snvz);
				nvs.add(DVector3d::v_ox * -sa + DVector3d::v_oy * -ca + snvz);
				nvs.add(DVector3d::v_ox *  sa + DVector3d::v_oy *  ca - snvz);
				nvs.add(DVector3d::v_ox * -sa + DVector3d::v_oy *  ca - snvz);
				nvs.add(DVector3d::v_ox *  sa + DVector3d::v_oy * -ca - snvz);
				nvs.add(DVector3d::v_ox * -sa + DVector3d::v_oy * -ca - snvz);
			}
		}
		// -> poles
		nvs.add(DVector3d::v_oz);
		nvs.add(-DVector3d::v_oz);

		DVector3d dv1, dv2;
		double hx = h * source_h_ratio_all;
		double d[] = { hx * source_h_ratio_normal, hx, hx };
		nvs.forEach([&](const DVector3d& dv) {
			dv.orthonormalVectors(dv1, dv2);
			sources_list.add(std::make_shared<ControlNode3d>(pt_mid + dv*r,
				ControlDataMatrix3d(dv, dv1, dv2, d)));
		});

		// -> finish
	}

	void showSources() {
		if (sources_list.empty()) return;

		MeshViewSet * sset = new MeshViewSet;

		sources_list.forEach([&sset](auto node) { 
			sset->addPoint(node->coord); 
			DMatrix3d de;
			double dd[3];
			if (node->control_data.eigensystem(de, dd)) {
				for (int i = 0; i < 3; i++) {
					sset->addEdge(node->coord, node->coord + de.column(i) * dd[i], 2);
				}
			}
		});
		box.drawXYZ(sset, 1);

		SHOW_MESH("sources", sset);
	}

	void setFSources(const SurfaceParametric * sp, 
		double cr,
		double x0 = -1.0, double x1 = 1.0, 
		double y0 = -1.0, double y1 = 1.0) 
	{
		DPoint2d param(x0, y0);

		double h = ControlSpace2dAdaptive::param_curvature_ratio * cr;
		//CDM3d cdm(h);
		double h_min = h;
		double h_step = h_min * DiscreteSources::density_ratio;

		CDM3d cdm;
		SurfaceCurvature sc;
		double p = box.getDiameter();

		auto src_x = std::make_shared<DataVector<double>>();
		auto src_y = std::make_shared<DataVector<double>>();

		disc_src_x.add(src_x);
		disc_src_y.add(src_y);

		src_x->add(param.x);
		src_y->add(param.y);

		double last_y = param.y;
		while (param.y <= y1) {
			double min_dy = (y1-y0);
			while (param.x <= x1) {
				DPoint3d spt = sp->getPoint(param);

				if (param.x > src_x->last())  src_x->add(param.x);
				if (param.y > src_y->last())  src_y->add(param.y);

				//if (!ControlSpace3dAdaptive::calculateControlFromSurfaceCurvature(spt, sp, p, cdm, &param, &sc))
				//	continue;
					
				// sources_list.add(cdm); // oryginal cdm

				//double h_min_xy = cdm.minEigenvalue();
				const DVector3d svn = sp->getNormalVector(param);

				//if (last_y < 0.0 && param.y >= 0.0) {
				//	param.storeSimple(MESHLOG, '\t');
				//	LOG4CPLUS_INFO(MeshLog::logger_mesh, " h_min_xy = " << h_min_xy);
				//}

				DVector3d dv1, dv2;
				double hx = h_min * source_h_ratio_all;
				double d[] = { hx * source_h_ratio_normal, hx, hx };
				svn.orthonormalVectors(dv1, dv2);
				sources_list.add(std::make_shared<ControlNode3d>(spt,	
					ControlDataMatrix3d(svn, dv1, dv2, d)));

				double dy = h_step / sp->getDerivative(DEquation::deriv_dt, param).length();
				if (dy < min_dy) min_dy = dy;

				if (param.x < x1) {
					double dx = h_step / sp->getDerivative(DEquation::deriv_ds, param).length();
					param.x += dx;
					if (param.x > x1) param.x = x1;
				}
				else
					break;
			}
			if (param.y < y1) {
				param.x = x0;
				last_y = param.y;
				param.y += min_dy;
				if (param.y > y1) param.y = y1;
			}
			else
				break;
		}

		// -> finish
	}
	
	void prepareEllipsoidFunc(const DataVector<DEllipsoid>& ellipsoids) {
		prepareFullDiscreteSourcesFunc();
	}

	void prepareSphereFunc(const DataVector<DSphere>& spheres, double cr = -1.0)
	{
		sources_func = [this, spheres, cr](const DPoint3d& pt) -> CDM3d {
			ControlNode3d cn1(DPoint3d::zero, CDM3d(max_metric_len));
			spheres.forEach([&](const DSphere& s) {
				double scr = s.getRadius();
				if (cr > 0.0) scr = cr;
				auto dv = pt - s.getCenter();
				double dist = dv.length();
				double h = ControlSpace2dAdaptive::param_curvature_ratio * scr;
				CDM3d cdm(h);
				if (cr > 0.0 && dist > 0.0) {
					DVector3d dv1, dv2;
					dv /= dist;
					dv.orthonormalVectors(dv1, dv2);
					double hx = h * source_h_ratio_all;
					double d[] = { hx * source_h_ratio_normal, hx, hx };
					cdm = CDM3d(dv, dv1, dv2, d);
				}
				ControlNode3d cn2(DPoint3d(abs(s.getRadius() - dist), 0.0, 0.0), cdm);
				DVector3d cdv(cn2.coord - cn1.coord);
				ControlSpace3dAdaptive::smoothenMetricForNodes(&cn1, &cn2, cdv);
			});
			return cn1.control_data;
		};
	}

	DRect setRectRangeX(double x, double y, const DataVector<double>& tx, const DataVector<double>& ty) {
		DRect rect(tx.first(), tx.last(), ty.first(), ty.last());
		setRange(x, tx, rect.x0, rect.x1);
		if (y - rect.y0 < rect.y1 - y) rect.y0 = y;		else rect.y1 = y;
		return rect;
	}

	DRect setRectRangeY(double x, double y, const DataVector<double>& tx, const DataVector<double>& ty) {
		DRect rect(tx.first(), tx.last(), ty.first(), ty.last());
		setRange(y, ty, rect.y0, rect.y1);
		if (x - rect.x0 < rect.x1 - x) rect.x0 = x;		else rect.x1 = x;
		return rect;
	}

	void setRange(double x, const DataVector<double>& t, double& a, double & b)
	{
		int i0 = 0;
		int ilast = t.countInt() - 1;
		int i1 = ilast;
		while (i1-i0 > 1) {
			int im = (i0 + i1) / 2;
			((t[im] > x) ? i1 : i0) = im;
		}
		if (t[i1] - x < x - t[i0]) i0 = i1; 
		// i0 -> nearest value
		a = ((i0 == 0) ? (t[0] - 0.5*(t[1]-t[0])) : (0.5*(t[i0]+t[i0-1])));
		b = ((i0 == ilast) ? (t[ilast] + 0.5*(t[ilast] - t[ilast-1])) : (0.5*(t[i0] + t[i0 + 1])));
	}

	void preparePlaneSourcesFunc(
		std::shared_ptr<const SurfacePlane> p0, const CDM3d& cdm0,
		std::shared_ptr<const SurfacePlane> p1, const CDM3d& cdm1)
	{
		sources_func = [this, p0, cdm0, p1, cdm1](const DPoint3d& pt) -> CDM3d {
			CDM3d cdm_ret(max_metric_len);

			auto ppt0 = p0->getPoint(p0->getParameters(pt));
			auto ppt1 = p1->getPoint(p1->getParameters(pt));

			ControlSpace3dAdaptive::adjustMetricForNodesLeft(cdm_ret, cdm0, pt - ppt0);
			ControlSpace3dAdaptive::adjustMetricForNodesLeft(cdm_ret, cdm1, pt - ppt1);

			return cdm_ret;
		};

		if (false) {
			std::function<CDM3d(const DPoint3d & pt)> sources_func_full = [this](const DPoint3d& pt) -> CDM3d {
				ControlNode3d cn_ret(pt, CDM3d(max_metric_len));
				double max_h = max_metric_len;
				//int req_ct = 0;
				for (int i = 0; i < sources_list.countInt(); i++) {
					auto cn_i = sources_list[i];
					ControlSpace3dAdaptive::smoothenMetricForNodesLeft(&cn_ret, cn_i.get(), pt - cn_i->coord);
				}
				//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Req: " << req_ct << " out of " << sources_list.countInt());
				return cn_ret.control_data;
			};

			double max_diff = 0.0;
			for (int i = 0; i < 100; i++) {
				auto rnd_pt = box.getRandomPoint();
				auto cdm_planes = sources_func(rnd_pt);
				auto cdm_full = sources_func_full(rnd_pt);

				double diff = cdm_planes.countDifferenceRR(cdm_full);
				if (diff > max_diff) max_diff = diff;
			}
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "max_diff = " << max_diff);
		}
	}

	void prepareSinDiscreteSourcesFunc(std::shared_ptr<const SurfacePlane> psin)
	{
		sources_func = [this,psin](const DPoint3d& pt) -> CDM3d {

			CDM3d cdm_ret(max_metric_len);

			if (use_src_tree) {
				bool is_in_tree = false;
				#pragma omp critical(src_tree_update)
				{
					is_in_tree = src_tree.isInCacheTree(pt, cdm_ret);
				}
				if(is_in_tree) return cdm_ret;
			}

			// select series (v)
			// adjust metric for all (u)

			auto spt = psin->getParameters(pt);

			assert(disc_src_y.countInt() == 1);
			double y0, y1;
			setRange(spt.y, *disc_src_y[0], y0, y1);

			sources_list.forEach([psin,&pt,y0,y1,&cdm_ret](auto cn_i) {
				auto spt_i = psin->getParameters(cn_i->coord);
				if (spt_i.y >= y0 && spt_i.y <= y1)
					ControlSpace3dAdaptive::adjustMetricForNodesLeft(cdm_ret,
						cn_i->control_data, pt - cn_i->coord);
			});

			if (use_src_tree) {
				#pragma omp critical(src_tree_update)
				src_tree.insertIntoCacheTree(pt, cdm_ret);
			}

			return cdm_ret;
		};

		//if (true) {
		//	//std::function<CDM3d(const DPoint3d & pt)> 
		//	auto sources_func_full = [this](const DPoint3d& pt) -> CDM3d {
		//		ControlNode3d cn_ret(pt, CDM3d(max_metric_len));
		//		double max_h = max_metric_len;
		//		for (int i = 0; i < sources_list.countInt(); i++) {
		//			const ControlNode3d& cn_i = sources_list[i];
		//			//ControlSpace3dAdaptive::smoothenMetricForNodesLeft(&cn_ret, &cn_i, pt - cn_i.coord);
		//			ControlSpace3dAdaptive::adjustMetricForNodesLeft(cn_ret.control_data, 
		//				cn_i.control_data, cn_i.coord-pt);
		//		}
		//		return cn_ret.control_data;
		//	};
		//
		//	double max_diff = 0.0;
		//	for (int i = 0; i < 100; i++) {
		//		auto rnd_pt = box.getRandomPoint();
		//		auto ct0 = ControlDataMatrix3d::m_intersection_counter;
		//		auto cdm_filtered = sources_func(rnd_pt);
		//		auto ct1 = ControlDataMatrix3d::m_intersection_counter;
		//		auto cdm_full = sources_func_full(rnd_pt);
		//		auto ct2 = ControlDataMatrix3d::m_intersection_counter;
		//
		//		int ct_diff_filtered = ct1 - ct0;
		//		int ct_diff_full = ct2 - ct1;
		//		double ct_ratio = ct_diff_filtered / (double)ct_diff_full;
		//
		//		double diff = cdm_filtered.countDifferenceRR(cdm_full);
		//		if (diff > max_diff) max_diff = diff;
		//	}
		//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "max_diff = " << max_diff);
		//}
	}



	void prepareOptDiscreteSourcesFunc()
	{
		sources_func = [this](const DPoint3d& pt) -> CDM3d {

			CDM3d cdm_ret(max_metric_len);

			if (use_src_tree) {
				if (src_tree.isInCacheTree(pt, cdm_ret))
					return cdm_ret;
			}


			//int sct = disc_src_x.countInt();
			//assert(sct > 0);
			//assert(sct == disc_src_y.countInt());

			//const SurfacePlane* planes[] = { &p0, &p1 };
			//for (int i = 0; i < sct; i++) {
			//	auto spt = planes[i]->getParameters(pt);

			//	auto rx = setRectRangeX(spt.x, spt.y, *disc_src_x[i], *disc_src_y[i]);
			//	auto ry = setRectRangeY(spt.x, spt.y, *disc_src_x[i], *disc_src_y[i]);

			//	for (const ControlNode3d& cn_i : sources_list) {
			//		auto spt_i = planes[i]->getParameters(cn_i.coord);
			//		if (rx.contains(spt_i) || ry.contains(spt_i)) {
			//			ControlSpace3dAdaptive::adjustMetricForNodesLeft(cdm_ret,
			//				cn_i.control_data, pt - cn_i.coord);
			//		}
			//	}
			//}

			sources_list.forEach([&pt, &cdm_ret](auto cn_i) {
				ControlSpace3dAdaptive::adjustMetricForNodesLeft(cdm_ret,
					cn_i->control_data, pt - cn_i->coord);
			});

			if (use_src_tree)
				src_tree.insertIntoCacheTree(pt, cdm_ret);

			return cdm_ret;
		};

		//if (true) {
		//	//std::function<CDM3d(const DPoint3d & pt)> 
		//	auto sources_func_full = [this](const DPoint3d& pt) -> CDM3d {
		//		ControlNode3d cn_ret(pt, CDM3d(max_metric_len));
		//		double max_h = max_metric_len;
		//		for (int i = 0; i < sources_list.countInt(); i++) {
		//			const ControlNode3d& cn_i = sources_list[i];
		//			//ControlSpace3dAdaptive::smoothenMetricForNodesLeft(&cn_ret, &cn_i, pt - cn_i.coord);
		//			ControlSpace3dAdaptive::adjustMetricForNodesLeft(cn_ret.control_data, 
		//				cn_i.control_data, cn_i.coord-pt);
		//		}
		//		return cn_ret.control_data;
		//	};

		//	double max_diff = 0.0;
		//	for (int i = 0; i < 100; i++) {
		//		auto rnd_pt = box.getRandomPoint();
		//		auto ct0 = ControlDataMatrix3d::m_intersection_counter;
		//		auto cdm_filtered = sources_func(rnd_pt);
		//		auto ct1 = ControlDataMatrix3d::m_intersection_counter;
		//		auto cdm_full = sources_func_full(rnd_pt);
		//		auto ct2 = ControlDataMatrix3d::m_intersection_counter;

		//		int ct_diff_filtered = ct1 - ct0;
		//		int ct_diff_full = ct2 - ct1;
		//		double ct_ratio = ct_diff_filtered / (double)ct_diff_full;

		//		double diff = cdm_filtered.countDifferenceRR(cdm_full);
		//		if (diff > max_diff) max_diff = diff;
		//	}
		//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "max_diff = " << max_diff);
		//}
	}

	void prepareFullDiscreteSourcesFunc() {
		assert(sources_list.notEmpty());

		sources_func = [this](const DPoint3d& pt) -> CDM3d {
			ControlNode3d cn_ret(pt, CDM3d(max_metric_len));
			for (int i = 0; i < sources_list.countInt(); i++) {
				auto cn_i = sources_list[i];
				ControlSpace3dAdaptive::adjustMetricForNodesLeft(cn_ret.control_data, 
					cn_i->control_data, pt - cn_i->coord);
			}
			//LOG4CPLUS_INFO(MeshLog::logger_mesh, "Req: " << req_ct << " out of " << sources_list.countInt());
			return cn_ret.control_data;
		};
	}
public:
	struct SrcTree
	{
		static const int MAX_ITEM_CT = 16;
		static double metric_dist_eps2;
		struct SrcCdmData {
			DPoint3d pt;
			CDM3d cdm;
			CDM3d cdmi;
		};
		struct SrcTreeNode
		{
		public:
			int selectChild(const DPoint3d& pt, DPoint3d& mid, DVector3d& dv) const {
				int vi = 0;
				if (pt.x > mid.x) {	vi |= 1; mid.x += dv.x;	}
				else mid.x -= dv.x;
				if (pt.y > mid.y) {	vi |= 2; mid.y += dv.y;	}
				else mid.y -= dv.y;
				if (pt.z > mid.z) {	vi |= 4; mid.z += dv.z;	}
				else mid.z -= dv.z;

				dv *= 0.5;
				
				return vi;
			}
			int selectChild(const DPoint3d& pt, const DPoint3d& mid) const {
				int vi = 0;
				if (pt.x > mid.x) vi |= 1; 
				if (pt.y > mid.y) vi |= 2; 
				if (pt.z > mid.z) vi |= 4; 

				return vi;
			}
			bool isInCacheTree(const DPoint3d& pt, CDM3d & tree_cdm,
				DPoint3d& mid, DVector3d& dv) const 
			{
				if (leaves) return leaves[selectChild(pt, mid, dv)].isInCacheTree(pt, tree_cdm, mid, dv);
				else {
					for (auto it = items.iterator(); it.valid(); it.moveNext()) {
						const auto & data = it.item();
						double lm2 = (data.cdmi * (pt - data.pt)).length2();
						if (lm2 < metric_dist_eps2) {
							tree_cdm = data.cdm;
							return true;
						}
					}
				}
				
				return false;
			}
			int insertIntoCacheTree(const DPoint3d& pt, const CDM3d & cdm, 
				DPoint3d& mid, DVector3d& dv, int level) 
			{
				if (leaves) return leaves[selectChild(pt, mid, dv)].insertIntoCacheTree(pt, cdm, mid, dv, level+1);
				else if (items_count >= MAX_ITEM_CT) {
					if(level >= 20) return -1;
					// split
					leaves.reset(new SrcTreeNode[8]);
					while (items.notEmpty()) {
						int vi = selectChild(items.getFirst().pt, mid);
						items.moveFirstTo(leaves[vi].items);
						leaves[vi].items_count++;
					}
					items_count = -1;
					// insert
					return leaves[selectChild(pt, mid, dv)].insertIntoCacheTree(pt, cdm, mid, dv, level+1);
				}
				else {
					items.append(SrcCdmData{ pt, cdm, cdm.inverse() });
					items_count++;
					return level;
				}
			}
			void clearAll() {
				leaves.release();
				items.clear();
				items_count = 0;
			}
		public:
			std::unique_ptr<SrcTreeNode[]> leaves;
			DataSimpleList<SrcCdmData> items;
			int items_count = 0;
		};
	public:
		SrcTree(const DBox& bbox, double _metric_dist_eps = 1.0) 
			: middle(bbox.getMiddlePoint()), dxzy(bbox.getDiameterVec()*0.25) 
		{
			metric_dist_eps2 = _metric_dist_eps * _metric_dist_eps;
		}
	public:
		bool isInCacheTree(const DPoint3d& pt, CDM3d & tree_cdm) const {
			DPoint3d middle_copy = middle;
			DVector3d dxzy_copy = dxzy;
			bool res = root.isInCacheTree(pt, tree_cdm, middle_copy, dxzy_copy);
		
			counter_all++;
			if (res) counter_hit++;

			if (counter_all % 1000 == 1) {
				LOG4CPLUS_INFO(MeshLog::logger_mesh, 
					(double)(100.0 * counter_hit / counter_all) << "\t"
					<< "(" << counter_hit << "/" << counter_all << ") tree-nodes: " << counter_items 
					<< "\tmax-depth: " << counter_max_depth);
			}

			return res;
		}
		bool insertIntoCacheTree(const DPoint3d& pt, const CDM3d & cdm) {
			DPoint3d middle_copy = middle;
			DVector3d dxzy_copy = dxzy;
			int depth = root.insertIntoCacheTree(pt, cdm, middle_copy, dxzy_copy, 0);

			if(depth >= 0) counter_items++;
			if(depth > counter_max_depth) counter_max_depth = depth;

			return (depth >= 0);
		}
		void resetWithNewMetricDist(double metric_dist) {
			root.clearAll();
			counter_all = counter_hit = counter_items = counter_max_depth = 0;
			metric_dist_eps2 = metric_dist * metric_dist;
		}
	private:
		SrcTreeNode root;
		DPoint3d middle;
		DVector3d dxzy;
		mutable int counter_all = 0;
		mutable int counter_hit = 0;
		int counter_items = 0;
		int counter_max_depth = 0;
	} src_tree;
public:
	DataVector<std::shared_ptr<ControlNode3d>> sources_list;
	std::function<CDM3d(const DPoint3d & pt)> sources_func;
	bool analytic;
	static double density_ratio;
	double source_h_ratio_all;
	double source_h_ratio_normal;
	double max_metric_len;
	string main_dir;
	DBox box;
	DataVector< std::shared_ptr<DataVector<double>> > disc_src_x;
	DataVector< std::shared_ptr<DataVector<double>> > disc_src_y;
	bool use_src_tree = false;
};
double DiscreteSources::density_ratio = 1.0;
double DiscreteSources::SrcTree::metric_dist_eps2 = 1.0;

bool checkMeshForCS(DataVector<double> & mesh_stats,
	const string& mesh_file, const string& mesh_out_file,
	CS3dPtr cs3d, CS3dPtr ref_cs3d = nullptr,
	DiscreteSources::SrcTree * ds_src_tree = nullptr)
{
	MeshBRepXML desc;
	if (!desc.parseFile(mesh_file) || !desc.validate()) return false;

	MeshContainer3d* domain = desc.createDomainMesh();
	assert(domain);

	assert(domain->getBlocksCount() == 1);
	MeshDomainVolume* mdv = ((MeshDomainVolume*)domain->getBlockAt(0));

	//mdv->clearControlSpace();
	mdv->setControlSpace(cs3d);
	//mdv->createInitialControlSpace(); // should only copy the user-cs ...

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Starting 3d triangulation");

	MeshGenerator2d::show_prediction = false;
	START_CLOCK("3D Triangulation");
	//MeshGenerator2d::autoTriangulate(domain, 2, false, true);
	MeshGenerator2d::triangulateWithCS3d(domain);
	mdv->prepareBoundaryMesh();
	Metric3dContext mc(cs3d);

	MeshDomainVolume::param_min_volume_for_simple_convex = 0.0;
	mdv->discretizeUsingBoundary(mc);

	MeshContainer3d* mesh = mdv->getMesh();

	LOG4CPLUS_DEBUG(MeshLog::logger_console, " - boundary done, adding inner nodes...");

	MeshGenerator3d::addInnerNodes(mc, mesh);

	LOG4CPLUS_DEBUG(MeshLog::logger_console, " - done, now smoothing...");
	MeshGenerator3dQuality::smoothen(mc, mesh, 2, TagExtended::TAG_NONE, 1, true);
	//,	MeshData::SM3_OPT_MIXED | MeshData::SM3_SWAP_COMPLETE | MeshData::SM3_BORDER_PRUNE);
	mesh->setDiscretizationState(2);
	STOP_CLOCK("3D Triangulation");

	LOG4CPLUS_DEBUG(MeshLog::logger_console, "3d triangulation done.");
	//SHOW_MESH("mesh", mesh->getViewSet());

	//mesh->storeTetmesh(mesh_out_file);
	//mesh->storeFacesOFF(mesh_out_file);
	//MeshBRepXML::storeFile(mesh_out_file + ".xml", domain);

	if (true) {

		//SHOW_MESH("mesh", mesh->getViewSet());
		MeshViewSet::ClipPlane cp(FVector3d(0.0f, 0.0f, 1.0f), 0.0f);
		//MeshViewSet::ClipPlane cp(FVector3d(-0.360846f, -0.0662537f, 0.930277f), -0);
		DataVector<MeshViewSet::ClipPlane> clip_planes(1);
		clip_planes.add(cp);

		auto set = mesh->getViewSetWithVisibleFacesOnly(&clip_planes);
		//	auto set = mesh->getViewSetWithVisibleBlocksOnly(&clip_planes);

		float mv_m[] = {
			0.296331f, -2.22229e-08f, 0.399142f, 0.095429f, 
			0.257879f, 0.751049f, -0.191454f, -0.0457738f, 
			0.00603026f, -0.00321179f, -0.00447699f, -0.00107038f, 
			0, 0, 0, 1 };

		float tm_m[] = {
			0.596098f, -4.47035e-08f, 0.802912f, 0.191964f,
			0.315702f, 0.919455f, -0.234383f, -0.0560375f,
			-0.738241f, 0.393196f, 0.548085f, 0.131039f,
			0, 0, 0, 1 };

		////** QtMeshGen *** modelview:
		//float mv_m[] = {
		//	0.49768f,	0.0f,	0.0f,	0.0f,
		//	0.0f,	0.817765f,	0.0f,	0.0f,
		//	0.0f,	0.0f, -0.00817765f, -0.00343403f,
		//	0.0f,	0.0f,	0.0f,	1.0f
		//};
		////** QtMeshGen *** tm :
		//float tm_m[] = {
		//	1,	0,	0,	0,
		//	0,	1,	0,	0,
		//	0,	0,	1,	0.419928f,
		//	0,	0,	0,	1
		//};

		set->setPolygonFillMode(MeshViewSet::FILL_LGRAY);

		string out_fname = mesh_out_file;
		int flen = (int)out_fname.length();
		for (int i = 0; i < flen; i++)
			if (out_fname[i] == '.')
				out_fname[i] = '_';


		//string fname = out_fname + ".eps";
		//set->storeEPSFile(tm_m, mv_m, 15, fname);

		if (true) {
			DMatrix3d mx_id, mx_rot_x, mx_rot_y;
			float id_m[] = {
				1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1
			};
			double ix = -0.9;
			double iy = -0.2;
			//for (double ix = -0.9; ix <= 0.7; ix += 0.1) {
			//	for (double iy = -0.9; iy <= 0.7; iy += 0.1) {
					mx_id.setIdentity(0.7);
					double ay = iy * PI;
					mx_rot_y.setRotationAroundNormalizedVector(DVector3d::v_oy, sin(ay), cos(ay));
					double ax = ix * PI;
					mx_rot_x.setRotationAroundNormalizedVector(DVector3d::v_ox, sin(ax), cos(ax));

					DMatrix3d mx_rot = mx_rot_y * mx_rot_x;
					DMatrix3d mx_rot_s = mx_id * mx_rot;

					int k = 0;
					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 3; j++) {
							tm_m[k] = (float)mx_rot_s(i, j);
							mv_m[k] = (float)mx_rot(i, j);
							k++;
						}
						k++;
					}

					//string fname = "xxx_" + to_string(ix) + "_" + to_string(iy) + ".eps";
					string fname = out_fname + ".eps";
					set->storeEPSFile(tm_m, mv_m, 15, fname);
			//	}
			//}
		}
		delete set;
	}

	//tex_file
	//	<< "\n\\noindent\n"
	//	<< "\\begin{minipage}[c]{0.4\\textwidth}\n"
	//	<< "\\includegraphics[width=\\textwidth]{" << out_fname << ".pdf}\n"
	//	<< "\\end{minipage}\n";
	//tex_file
	//	<< "\\begin{minipage}[c]{0.6\\textwidth}\n"
	//	<< "\\begin{verbatim}\n"
	//	<< mesh_out_file << "\n"
	//	<< "NT=" << mesh->getBlocksCount() << "\n";


	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "=== Mesh-CS ====== " << mesh_out_file << " ===========");
	mesh->logQuality(false);

	if (ref_cs3d != nullptr) {
		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "=== Ref-CS ====== " << mesh_out_file << " ===========");
		CS3dPtr oryg_cs3d = mesh->getControlSpace();
		mesh->setControlSpace(ref_cs3d);
		if (ds_src_tree == nullptr) {
			mesh->logQuality(false);
		}
		else {
			for (double md : {1.5, 1.0, 0.8, 0.5, 0.3, 0.1}) {
				ds_src_tree->resetWithNewMetricDist(md);
				LOG4CPLUS_INFO(MeshLog::logger_mesh, " ==> SRC-TREE: metric distance:" << md);
				LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, " ==> SRC-TREE: metric distance:" << md);
				mesh->logQuality(false);
			}
		}
		mesh->setControlSpace(oryg_cs3d);
	}

	mesh_stats.add(mesh->getBlocksCount());
	DataStatistics estats, qstats;
	Metric3dContext mcs(ref_cs3d ? ref_cs3d : cs3d);
	if (mesh->statMetricEdgeLength(mcs, estats, false) && estats.calculate()) {
		//double e_ave = estats.average();
		//double e_sd = estats.stdDev();
		//mesh_stats.add(e_ave);
		//mesh_stats.add(e_sd);
		int low_ct = estats.getDataCountBottom(0.8);
		int hi_ct = estats.getDataCountTop(1.25);
		int all_ct = estats.countInt();
		int good_ct = all_ct - low_ct - hi_ct;
		double e_good_ratio = ((double)good_ct / all_ct);
		mesh_stats.add(e_good_ratio);
		//tex_file.precision(2);
		//tex_file
		//	<< "ELen_ave=" << e_ave << "\n"
		//	<< "ELen_sdev=" << e_sd << "\n"
		//	<< "Elen_[0.8,1.25]=" << good_ct << " / " << all_ct << " (" << e_good_ratio << ")\n";
	}
	if (mesh->statMetricDifference(mcs, qstats) && qstats.calculate()) {
		//double m_ave = qstats.average();
		//double m_sd = qstats.stdDev();
		//mesh_stats.add(m_ave);
		//mesh_stats.add(m_sd);
		int good_ct = qstats.getDataCountBottom(2.0);
		int all_ct = qstats.countInt();
		double m_good_ratio = ((double)good_ct / all_ct);
		mesh_stats.add(m_good_ratio);
		//tex_file.precision(2);
		//tex_file
		//	<< "MDiff_ave=" << m_ave << "\n"
		//	<< "MDiff_sdev=" << m_sd << "\n"
		//	<< "MDiff_[0,2]=" << good_ct << " / " << all_ct << " (" << m_good_ratio << ")\n";
	}

	delete domain;

	//tex_file
	//	<< "\\end{verbatim}\n"
	//	<< "\\end{minipage}\n";

	return true;
}

void runTestKdTree(
	CS3dPtr tree_cs,
	const DBox& box,
	double accuracy, // double accuracy_range[] = { 4.0, 2.0, 1.0, 0.5, 0.3, 0.1, 0.08, 0.06 };
	int f_type
)
{
	auto tree = (ControlSpace3dKdTree*)(tree_cs.get());
	//	file << "TEST\t ACC\t MAX_DIFF\t ELEM\t BYTES\t DEPTH_MAX\t BAL_MAX\t ADAPT[ms]\t GET[ms]"
	//		<< "\t MESH_NT\t MESH_LE_AVE\t MESH_LE_SD\t MESH_LE_[0.8,1.25]\t MESH_Q_AVE \t MESH_Q_SD\t MESH_Q_[0,2]"
	//		<< endl;

	CDM3d(*f)(const DPoint3d& pt) = test_functions[f_type];

	clock_t cstart = clock();
	int tree_size = tree->adaptToField(f, accuracy);
	clock_t cadapt_ms = (clock() - cstart);
	//cout << " * adapted, tree size:" << tree_size << endl;

	double max_diff = 0.0;

	int GRID_SIZE = 100;
	double fx = 1.0 / (double)GRID_SIZE;
	double dx = box.getDX() * fx;
	double dy = box.getDY() * fx;
	double dz = box.getDZ() * fx;

	DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
	for (; pt.x < box.x1; pt.x += dx) {
		for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
			for (pt.z = box.z0 + 0.5*dz; pt.z < box.z1; pt.z += dz) {
				double diff = diffKdValue(tree->getMetricAtPoint(pt), f(pt));
				if (diff > max_diff) max_diff = diff;
			}
		}
	}

	cstart = clock();
	CDM3d v_total;
	const int RPT = 4;

	dx /= RPT;
	dy /= RPT;
	dz /= RPT;

	pt.x = box.x0 + 0.5*dx;
	pt.y = box.y0;
	pt.z = box.z0;
	for (; pt.x < box.x1; pt.x += dx) {
		for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
			for (pt.z = box.z0 + 0.5*dz; pt.z < box.z1; pt.z += dz) {
				CDM3d v = tree->getMetricAtPoint(pt);
				v_total += v;
				v_total *= 0.5;
			}
		}
	}
	double fget_rpt = (RPT*RPT*RPT);
	double cget_ms = (double)(clock() - cstart) / fget_rpt;

	DataVector<double> mesh_stats;
	// meshing
	string label = tree->getLabel();
	if (false) {
		string mesh_fname = label + "-" + to_string(accuracy);
		CS3dPtr ref_cs(new ControlSpace3dFunctional(f));
		checkMeshForCS(mesh_stats, "cube.xml", mesh_fname, 
			tree_cs, ref_cs);
	}

	// kdtree-Lx4-LONG
	cout << label << "-" << test_labels[f_type] << "\t"
		<< accuracy << "\t"
		<< max_diff << "\t"
		<< tree_size << "\t"
		<< tree->getTotalBytes() << "\t"
		<< tree->getMaxDepth() << "\t"
		<< tree->getMaxBalance() << "\t"
		<< cadapt_ms << "\t"
		<< cget_ms
		<< endl;
}

void runTestKdTreeSS(
	DataVector<CS3dKdPtr> & trees,
	const DBox& box,
	double accuracy, // double accuracy_range[] = { 4.0, 2.0, 1.0, 0.5, 0.3, 0.1, 0.08, 0.06 };
	DiscreteSources& cdm_sources
)
{
	if (trees.empty()) return;

	//	file << "TEST\t ACC\t MAX_DIFF\t ELEM\t BYTES\t DEPTH_MAX\t BAL_MAX\t ADAPT[ms]\t GET[ms]"
	//		<< "\t MESH_NT\t MESH_LE_AVE\t MESH_LE_SD\t MESH_LE_[0.8,1.25]\t MESH_Q_AVE \t MESH_Q_SD\t MESH_Q_[0,2]"
	//		<< endl;

	trees.forEach([](CS3dKdPtr& tree) { tree->setMaxMetric(); });
	KdTree3dParams::paramMetricLenRatio = accuracy;

	if (false) {
		CDM3d cdm_q0 = DMetric3d::stretchToMatrix(CDS3d(4.0, 2.0, 2.0));
		const int ANGLE_STEPS = 10;
		for (int i = 1; i <= ANGLE_STEPS; i++) {
			double angle = PI * i / (double)ANGLE_STEPS;
			CDM3d cdm_q1 = DMetric3d::stretchToMatrix(CDS3d(2.0, 6.0, 2.0, angle, 0.0, 0.0));
			DMatrix3d e; double d[3];
			bool res = cdm_q1.eigensystem(e, d);
			double hx_dir = cdm_q0.getMaxHRatioDirect(cdm_q1);
			double hx_red = cdm_q0.getMaxHRatioReduction(cdm_q1);
			LOG4CPLUS_INFO(MeshLog::logger_mesh, angle << "\t" << hx_dir << "\t" << hx_red);
		}
	}

	clock_t cstart = clock();
	trees.forEach([&](CS3dKdPtr& tree) {
		cdm_sources.sources_list.forEach([&](auto cn) {
			tree->setMinControl(cn->coord, cn->control_data, true);
		});
		LOG4CPLUS_INFO(MeshLog::logger_console, " ... all sources inserted ...");
		//tree->fixCloseControlNodes();
		switch (KdTree3dParams::paramGradationMethod) {
		case GRADATION_NONE: break;
		case GRADATION_LEAVES: tree->smoothenViaLeaves(); break;
		case GRADATION_VIC_PAIRS: tree->smoothenViaVicinityPairs(); break;
		case GRADATION_VIC_NB: tree->smoothenViaVicinityNb(); break;
		case GRADATION_VIC_NB_MULTIPASS: tree->smoothenViaVicinityNbMultipass(); break;
		}
		LOG4CPLUS_INFO(MeshLog::logger_console, " ... gradation smoothed ...");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " tree adapted (+smoothen), tree size:" << tree->getElementCount());
		tree->compact(); // make "read-only"
		DPoint3d cpt = DPoint3d::zero;
		auto cdm = tree->getMetricAtPoint(cpt);
		LOG4CPLUS_INFO(MeshLog::logger_console, cpt <<
			" -> " << cdm.minEigenvalue() << " - " << cdm.maxEigenvalue());
		cpt = DPoint3d(1.0, 1.0, 1.0);
		cdm = tree->getMetricAtPoint(cpt);
		LOG4CPLUS_INFO(MeshLog::logger_console, cpt <<
			" -> " << cdm.minEigenvalue() << " - " << cdm.maxEigenvalue());
		cpt = DPoint3d(0.0, 1.0, 1.0);
		cdm = tree->getMetricAtPoint(cpt);
		LOG4CPLUS_INFO(MeshLog::logger_console, cpt <<
			" -> " << cdm.minEigenvalue() << " - " << cdm.maxEigenvalue());
		cpt = DPoint3d(0.0, 0.0, 1.0);
		cdm = tree->getMetricAtPoint(cpt);
		LOG4CPLUS_INFO(MeshLog::logger_console, cpt <<
			" -> " << cdm.minEigenvalue() << " - " << cdm.maxEigenvalue());
		cpt = DPoint3d(0.0, 0.0, 0.8);
		cdm = tree->getMetricAtPoint(cpt);
		LOG4CPLUS_INFO(MeshLog::logger_console, cpt <<
			" -> " << cdm.minEigenvalue() << " - " << cdm.maxEigenvalue());
		//tree->showMetricLength();
	});
	double cadapt_ms = 1000.0 * (clock() - cstart) / CLOCKS_PER_SEC / trees.countInt();

	auto tree = trees.get(0);

	DataVector<double> grad_data;
	if (false) {
		double dr1 = tree->calculateMaxGradationRatioViaNodes();
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Max Gradation Ratio after smoothing (via nodes) : " << dr1);
		double dr2 = tree->calculateMaxGradationRatioRandom();
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Max Gradation Ratio after smoothing (random)    : " << dr2);
		double dr3 = tree->calculateMaxGradationRatioRegular(11);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Max Gradation Ratio after smoothing (regular-11): " << dr3);
		double dr4 = tree->calculateMaxGradationRatioRegular(21);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "Max Gradation Ratio after smoothing (regular-21): " << dr4);
		grad_data.add(dr1);
		grad_data.add(dr2);
		grad_data.add(dr3);
		grad_data.add(dr4);
	}

	if (false) {
		int cnct = tree->getControlNodesCount();
		DataVector<std::shared_ptr<ControlNode3d>> cnodes(cnct);
		tree->forEachControlNode([&](std::shared_ptr<ControlNode3d> cn) {
			cnodes.add(cn);
		});

		double max_hr_direct = 0.0;
		double max_hr_reduce = 0.0;
		const int TOTAL_RPT = 2;

		DataVector<double> hr_direct(cnct);
		DataVector<double> hr_reduct(cnct);

		cstart = clock();
		for (int i = 0; i < cnct - 1; i++) {
			for (int j = i + 1; j < cnct; j++) {
				double hr = cnodes[i]->control_data.getMaxHRatioDirect(cnodes[j]->control_data);
				if (hr > max_hr_direct) max_hr_direct = hr;
				hr_direct.add(hr);
			}
		}
		double c_max_h_direct_ms = 1e3* (clock() - cstart) / CLOCKS_PER_SEC / TOTAL_RPT;

		cstart = clock();
		for (int i = 0; i < cnct - 1; i++) {
			for (int j = i + 1; j < cnct; j++) {
				double hr = cnodes[i]->control_data.getMaxHRatioReduction(cnodes[j]->control_data);
				if (hr > max_hr_reduce) max_hr_reduce = hr;
				hr_reduct.add(hr);
			}
		}
		double c_max_h_reduce_ms = 1e3* (clock() - cstart) / CLOCKS_PER_SEC / TOTAL_RPT;

		LOG4CPLUS_INFO(MeshLog::logger_mesh, "=====================");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "direct\treduction");
		int close_enough = 0;
		for (int i = 0; i < hr_direct.countInt(); i++) {
			double diff = abs(hr_direct[i] - hr_reduct[i]);
			if (diff > METRIC_SMALL_NUMBER)
				LOG4CPLUS_INFO(MeshLog::logger_mesh, hr_direct[i] << "\t" << hr_reduct[i]);
			else
				++close_enough;
		}
		LOG4CPLUS_INFO(MeshLog::logger_mesh, "=====================");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " CLOSE ENOUGH: " << close_enough);
		LOG4CPLUS_INFO(MeshLog::logger_mesh, " DIFFERENT:    " << (hr_direct.countInt() - close_enough));

		LOG4CPLUS_INFO(MeshLog::logger_console, "max_h_direct [val]" << max_hr_direct);
		LOG4CPLUS_INFO(MeshLog::logger_console, "max_h_direct [ms ]" << c_max_h_direct_ms);
		LOG4CPLUS_INFO(MeshLog::logger_console, "max_h_reduce [val]" << max_hr_reduce);
		LOG4CPLUS_INFO(MeshLog::logger_console, "max_h_reduce [ms ]" << c_max_h_reduce_ms);

	}

	double max_diff_src = 0.0;
	double cget_tree_ms = 0.0;
	double cget_src_ms = 0.0;
	double max_diff_frnd = 0.0;
	double max_diff_freg = 0.0;

	if (true) {

		START_CLOCK("DIFF-S");
		cdm_sources.sources_list.forEach([&](auto cn) {
			double diff = diffKdValue(tree->getMetricAtPoint(cn->coord), cn->control_data);
			if (diff > max_diff_src) max_diff_src = diff;
		});
		STOP_CLOCK("DIFF-S");

		const int REG_PROBE_CT = 30;
		const int RND_PROBE_CT = REG_PROBE_CT * REG_PROBE_CT * REG_PROBE_CT;
		DataVector<DPoint3d> rnd_points(RND_PROBE_CT);
		for (int i = 0; i < RND_PROBE_CT; i++)
			rnd_points.add(box.getRandomPoint());

		//DataVector<CDM3d> tree_val(RND_PROBE_CT);
		//DataVector<CDM3d> src_val(RND_PROBE_CT);

		//START_CLOCK("V-TREE");
		//for (int i = 0; i < RND_PROBE_CT; i++)
		//	tree_val.add(tree->getMetricAtPoint(rnd_points[i]));
		//STOP_CLOCK("V-TREE");

		//START_CLOCK("V-SRC");
		//for (int i = 0; i < RND_PROBE_CT; i++)
		//	src_val.add(cdm_sources.sources_func(rnd_points[i]));
		//STOP_CLOCK("V-SRC");

		cstart = clock();
		CDM3d v = CDM3d::identity;
		const int TOTAL_RPT = TEST_ACCESS_RPT * RND_PROBE_CT;
		for (int i = 0; i < TOTAL_RPT; i++) {
			v += tree->getMetricAtPoint(rnd_points[i%RND_PROBE_CT]);
			v *= 0.5;
		}
		cget_tree_ms = 1e6* (clock() - cstart) / CLOCKS_PER_SEC / TOTAL_RPT;

		//cstart = clock();
		//for (int i = 0; i < TOTAL_RPT; i++) {
		//	v += cdm_sources.sources_func(rnd_points[i%RND_PROBE_CT]);
		//	v *= 0.5;
		//}
		//cget_src_ms = 1e6 * (clock() - cstart) / CLOCKS_PER_SEC / TOTAL_RPT;

		//for (int i = 0; i < RND_PROBE_CT; i++) {
		//	double diff = diffKdValue(tree_val[i], src_val[i]);
		//	if (diff > max_diff_frnd) max_diff_frnd = diff;
		//}

		START_CLOCK("DIFF-R");
		box.forRegularGrid(REG_PROBE_CT, [&](const DPoint3d& pt) {
			double diff = diffKdValue(tree->getMetricAtPoint(pt), cdm_sources.sources_func(pt));
			if (diff > max_diff_freg) max_diff_freg = diff;
		});
		STOP_CLOCK("DIFF-R");

		LOG4CPLUS_INFO(MeshLog::logger_console, "CS tests (with tree depth " << 
			KdTree3dParams::paramMaxLevel << "/" << KdOctree3dParams::paramMaxLevel << ") finished.");
	}

	DataVector<double> mesh_stats;
	// meshing
	string label = tree->getLabel();
	if (true) {
		string mesh_fname = MeshLog::base_log_name;
		if (mesh_fname == "mesh") {
			mesh_fname = label + "-" + to_string(accuracy)
				+ "-" + to_string(cdm_sources.source_h_ratio_all)
				+ "-" + to_string(cdm_sources.source_h_ratio_normal);
		}
		cdm_sources.use_src_tree = true;
		auto ref_cs = std::make_shared<ControlSpace3dFunctional>(cdm_sources.sources_func);
		checkMeshForCS(mesh_stats, cdm_sources.main_dir + "cube.xml", mesh_fname, tree,
			ref_cs, nullptr /* &cdm_sources.src_tree */);
		//checkMeshForCS(mesh_stats, cdm_sources.main_dir + "cube.xml", mesh_fname, tree.ptr(), nullptr);
	}

	// kdtree-Lx4-LONG
	cout << label << "\t" 
		<< cdm_sources.sources_list.countInt() << "\t"
		<< accuracy << "\t"
		<< fixed << max_diff_src << "\t"
		<< max_diff_frnd << "\t"
		<< max_diff_freg << "\t"
		<< tree->getElementCount() << "\t"
		<< tree->getTotalBytes() << "\t"
		<< tree->getMaxDepth() << "\t"
		<< tree->getMaxBalance() << "\t"
		<< cadapt_ms << "\t"
		<< cget_tree_ms << "\t"
		<< cget_src_ms;
	grad_data.forEach([&](double v) { cout << "\t" << v; });
	mesh_stats.forEach([&](double v) { cout << "\t" << fixed << v; });
	cout << endl;
}

void ControlSpace3dKdTree::showMetricLength() const
{
	static const int RES = 30;
	double dx = m_box.getDX() / (RES - 1);
	double dy = m_box.getDY() / (RES - 1);
	double dz = m_box.getDZ() / (RES - 1);

	MeshViewSet* set = new MeshViewSet();

	DPoint3d pt = m_box.getMinPoint();
	double t;

	double min_l = -1.0, max_l;

	for (int iz = 0; iz < RES; iz++) {
		t = (double)iz / (RES-1);
		pt.z = m_box.z0 * (1.0 - t) + m_box.z1 * t;
		for (int iy = 0; iy < RES; iy++) {
			t = (double)iy / (RES - 1);
			pt.y = m_box.y0 * (1.0 - t) + m_box.y1 * t;
			for (int ix = 0; ix < RES; ix++) {
				t = (double)ix / (RES - 1);
				pt.x = m_box.x0 * (1.0 - t) + m_box.x1 * t;

				CDM3d cdm = getMetricAtPoint(pt);
				double mind = cdm.minEigenvalue();

				if (min_l < 0.0) min_l = max_l = mind;
				else {
					if (mind < min_l) min_l = mind;
					if (mind > max_l) max_l = mind;
				}
			}
		}
	}

	for (int iz = 0; iz < RES; iz++) {
		t = (double)iz / (RES - 1);
		pt.z = m_box.z0 * (1.0 - t) + m_box.z1 * t;
		for (int iy = 0; iy < RES; iy++) {
			t = (double)iy / (RES - 1);
			pt.y = m_box.y0 * (1.0 - t) + m_box.y1 * t;
			for (int ix = 0; ix < RES; ix++) {
				t = (double)ix / (RES - 1);
				pt.x = m_box.x0 * (1.0 - t) + m_box.x1 * t;

				CDM3d cdm = getMetricAtPoint(pt);
				double mind = cdm.minEigenvalue();

				double q = (mind - min_l) / (max_l - min_l);
				set->addTetra(pt, (float)dx/2, (float)dy/2, (float)dz/2, q);
			}
		}
	}

	set->setPolygonFillMode(MeshViewSet::FILL_QUALITY);
	SHOW_MESH("metric min length", set);
}


bool parseIntArg(const char* str, int & value, DataHashTableKeyValue<string, int> * harray = nullptr) {
	if (harray != nullptr) {
		string str_key = str;
		if (harray->contains(str_key)) {
			value = harray->getValue(string(str), -1);
			return true;
		}
	}
	istringstream is(str);
	is >> value;
	if (!is) {
		cout << "ERROR parsing arg: " << str << endl;
		return false;
	}
	return true;
}

bool parseDoubleArg(const char* str, double & value) {
	istringstream is(str);
	is >> value;
	if (!is) {
		cout << "ERROR parsing arg: " << str << endl;
		return false;
	}
	return true;
}

bool parseTwoDoubleArgs(const char* str, double & value1, double & value2) {
	istringstream is(str);
	char z;
	is >> value1 >> z >> value2;
	if (!is) {
		cout << "ERROR parsing arg: " << str << endl;
		return false;
	}
	return true;
}
enum TREE_TYPE { KDTREE_L, KDTREE_V, OCTREE_L, OCTREE_V, KDTREE_LI, OCTREE_LB, OCTREE_V_ORYG };
enum CF_SPLIT_TYPE { SPLIT_CF_LONGEST, SPLIT_CF_MAX_GR_SUM, SPLIT_CF_MAX_GR, SPLIT_CF_MAX_DIFF_GR };
enum FUNC_TYPE { F_LABBE_ISO, F_LABBE_ANISO, F_GAUSS, F_LINE_DIST, F_POINTS_DIST };

int ControlSpace3dKdTree::checkKdTreeCF(int argi, int argc, const char* argv[])
{
	//tex_file.open("mesh-plots.tex");

	//tex_file
	//	<< "\\documentclass[a4paper]{article}\n"
	//	<< "\\usepackage{a4wide}\n"
	//	<< "\\usepackage[pdftex]{graphicx}\n"
	//	<< "\\addtolength{\\hoffset}{-1.5cm}\n"
	//	<< "\\addtolength{\\textwidth}{3cm}\n"
	//	<< "\\addtolength{\\voffset}{-1.5cm}\n"
	//	<< "\\addtolength{\\textheight}{3cm}\n"
	//	<< "\\begin{document}\n"
	//	<< "\\pagestyle{empty}\n";

	//ofstream test_file("kdtree-cdm-iso.txt");


	//if (false) {
	//	CS3dPtrHolder fcs0(new ControlSpace3dFunctional(f_labbe_iso));
	//	CS3dPtrHolder fcs1(new ControlSpace3dFunctional(f_labbe_aniso));
	//	DataVector<double> mesh_stats_0, mesh_stats_1;
	//	checkMeshForCS(mesh_stats_0, "cube.xml", "cube-direct-labbe-iso", fcs0.ptr());
	//	double accs[] = { 4.0, 0.06 };
	//	for (double acc : accs) {
	//		test_file << "fdirect-labbe04-iso\t" << acc << "\t"
	//			<< "0\t0\t0\t0\t0\t0\t0\t";
	//		for (double v : mesh_stats_0)
	//			test_file << v << "\t";
	//		test_file << endl;
	//	}
	//	checkMeshForCS(mesh_stats_1, "cube.xml", "cube-direct-labbe-aniso", fcs1.ptr());
	//	for (double acc : accs) {
	//		test_file << "fdirect-labbe04-aniso\t" << acc << "\t"
	//			<< "0\t0\t0\t0\t0\t0\t0\t";
	//		for (double v : mesh_stats_1)
	//			test_file << v << "\t";
	//		test_file << endl;
	//	}
	//	//exit(0);
	//}

	// parse params
	int tree_type = KDTREE_L;		// 0..6
	int cf_split_type = SPLIT_CF_LONGEST; // 0..3
	double accuracy = 1.0;			// acc (0.06 .. 4.0)
	int func_type = F_LABBE_ISO; 	// 0..4

	DataHashTableKeyValue<string, int> h_tree_types(20, "");
	h_tree_types.insert("KDTREE_L", 0);
	h_tree_types.insert("KDTREE_V", 1);
	h_tree_types.insert("OCTREE_L", 2);
	h_tree_types.insert("OCTREE_V", 3);
	h_tree_types.insert("KDTREE_LI", 4);
	h_tree_types.insert("OCTREE_LB", 5);
	h_tree_types.insert("OCTREE_V_ORYG", 6);

	DataHashTableKeyValue<string, int> h_cf_split_types(20, "");
	h_cf_split_types.insert("SPLIT_LONGEST", 0);
	h_cf_split_types.insert("SPLIT_MAX_GR_SUM", 1);
	h_cf_split_types.insert("SPLIT_MAX_GR", 2);
	h_cf_split_types.insert("SPLIT_MAX_DIFF_GR", 3);

	DataHashTableKeyValue<string, int> h_f_types(20, "");
	h_f_types.insert("F_LABBE_ISO", 0);
	h_f_types.insert("F_LABBE_ANISO", 1);
	h_f_types.insert("F_GAUSS", 2);
	h_f_types.insert("F_LINE_DIST", 3);
	h_f_types.insert("F_POINTS_DIST", 4);

	if (++argi < argc && !parseIntArg(argv[argi], tree_type, &h_tree_types)) return -1;
	if (++argi < argc && !parseIntArg(argv[argi], cf_split_type, &h_cf_split_types)) return -1;
	if (++argi < argc && !parseDoubleArg(argv[argi], accuracy)) return -1;
	if (++argi < argc && !parseIntArg(argv[argi], func_type, &h_f_types)) return -1;
	if (++argi < argc && !parseIntArg(argv[argi], KdTree3dParams::paramSplitGradientSamplingSize)) return -1;
	if (++argi < argc && !parseIntArg(argv[argi], KdTree3dParams::paramEstimateErrorSamplingSize)) return -1;


	DBox tree_box(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

	switch (cf_split_type) {
	case SPLIT_CF_LONGEST:
		KdTree3dParams::m_cf_splitter = std::make_shared<Kd3dCFHalfLongestSplitter>();
		break;
	case SPLIT_CF_MAX_GR_SUM:
		KdTree3dParams::m_cf_splitter = std::make_shared<Kd3dCFMaxGradientSumSplitter>(
			KdTree3dParams::paramSplitGradientSamplingSize);
		break;
	case SPLIT_CF_MAX_GR:
		KdTree3dParams::m_cf_splitter = std::make_shared<Kd3dCFMaxGradientSplitter>(
			KdTree3dParams::paramSplitGradientSamplingSize);
		break;
	case SPLIT_CF_MAX_DIFF_GR:
		KdTree3dParams::m_cf_splitter = std::make_shared<Kd3dCFMaxDiffGradientSplitter>(
			KdTree3dParams::paramSplitGradientSamplingSize);
		break;
	default:
		cout << "Error: unknown splitter type: " << cf_split_type << endl;
		return -1;
	}

	CS3dPtr tree;

	switch (tree_type) {
	case KDTREE_L:
		tree = std::make_shared<ControlSpace3dKdTreeL>(tree_box);
		break;
	case KDTREE_LI:
		tree = std::make_shared<ControlSpace3dKdTreeLi>(tree_box);
		break;
	case KDTREE_V:
		tree = std::make_shared<ControlSpace3dKdTreeV>(tree_box);
		break;
	case OCTREE_L:
		tree = std::make_shared<ControlSpace3dKdOctreeL>(tree_box);
		break;
	case OCTREE_LB:
		tree = std::make_shared<ControlSpace3dKdOctreeLB>(tree_box);
		break;
	case OCTREE_V:
		tree = std::make_shared<ControlSpace3dKdOctreeV>(tree_box);
		break;
	case OCTREE_V_ORYG:
		tree = std::make_shared<ControlSpace3dOctree>(tree_box);
		break;
	default:
		cout << "Error: unknown tree type: " << tree_type << endl;
		return -1;
	}

	if (!tree) {
		cout << "Error: not implemented tree type: " << tree_type << endl;
		return -1;
	}

	runTestKdTree(tree, tree_box, accuracy, func_type);

	//if (true) {
	//	//TEST_TREE("kdtree-Lx2-", ControlDataMatrix3d, DataKdTreeL3d, 2);
	//	//TEST_TREE("kdtree-Lx3-", ControlDataMatrix3d, DataKdTreeL3d, 3);
	//	//TEST_TREE("kdtree-Lx6-", ControlDataMatrix3d, DataKdTreeL3d, 6);
	//	//TEST_TREE("kdtree-Lx10-", ControlDataMatrix3d, DataKdTreeL3d, 10);

	//	TEST_TREE("kdtree-Lx4-LONG-", ControlDataMatrix3d, DataKdTreeL3d, 4);
	//	TEST_TREE("kdtree-Vx4-LONG-", ControlDataMatrix3d, DataKdTreeV3d, 4);
	//	TEST_TREE("kdtree-Lix4-LONG-", ControlDataMatrix3d, DataKdTreeLi3d, 4);

	//	TEST_TREE_MAXGRADIENTSUMSPLIT("kdtree-Lx4-MAXSUMGR4-", ControlDataMatrix3d, DataKdTreeL3d, 4, 4);
	//	TEST_TREE_MAXGRADIENTSUMSPLIT("kdtree-Vx4-MAXSUMGR4-", ControlDataMatrix3d, DataKdTreeV3d, 4, 4);
	//	TEST_TREE_MAXGRADIENTSUMSPLIT("kdtree-Lix4-MAXSUMGR4-", ControlDataMatrix3d, DataKdTreeLi3d, 4, 4);

	//	TEST_TREE_MAXGRADIENTSUMSPLIT("kdtree-Lx4-MAXSUMGR8-", ControlDataMatrix3d, DataKdTreeL3d, 4, 8);
	//	TEST_TREE_MAXGRADIENTSUMSPLIT("kdtree-Vx4-MAXSUMGR8-", ControlDataMatrix3d, DataKdTreeV3d, 4, 8);
	//	TEST_TREE_MAXGRADIENTSUMSPLIT("kdtree-Lix4-MAXSUMGR8-", ControlDataMatrix3d, DataKdTreeLi3d, 4, 8);

	//	TEST_TREE_MAXGRADIENTSPLIT("kdtree-Lx4-MAXGRAD4-", ControlDataMatrix3d, DataKdTreeL3d, 4, 4);
	//	TEST_TREE_MAXGRADIENTSPLIT("kdtree-Vx4-MAXGRAD4-", ControlDataMatrix3d, DataKdTreeV3d, 4, 4);
	//	TEST_TREE_MAXGRADIENTSPLIT("kdtree-Lix4-MAXGRAD4-", ControlDataMatrix3d, DataKdTreeLi3d, 4, 4);

	//	//TEST_TREE_MAXDIFFGRADIENTSPLIT("kdtree-Lx4-MAXDIFFGRAD4-", ControlDataMatrix3d, DataKdTreeL3d, 4, 4);
	//	//TEST_TREE_MAXDIFFGRADIENTSPLIT("kdtree-Vx4-MAXDIFFGRAD4-", ControlDataMatrix3d, DataKdTreeV3d, 4, 4);
	//	//TEST_TREE_MAXDIFFGRADIENTSPLIT("kdtree-Lix4-MAXDIFFGRAD4-", ControlDataMatrix3d, DataKdTreeLi3d, 4, 4);
	//}

	//if (true) {
	//	TEST_OCTREE("octree-Lx4-", ControlDataMatrix3d, DataOctreeL, 4);
	//	TEST_OCTREE("octree-LBalx4-", ControlDataMatrix3d, DataOctreeLBalanced, 4);
	//	TEST_OCTREE("octree-Vx4-", ControlDataMatrix3d, DataOctreeV, 4);
	//}

	//tex_file << "\\end{document}\n";
	//tex_file.close();

	return 0;
}

enum SS_SPLIT_TYPE { SPLIT_SS_HALF_LONGEST, SPLIT_SS_GOLDEN_RATIO, SPLIT_SS_MIN_LONGEST };

bool loadCDMSources(const char* fname, DiscreteSources & ds)
{
	if (ds.matchPredefined(fname)) return true;
	else return ds.readFromFile(fname);
}

int ControlSpace3dKdTree::checkKdTreeSS(int argi, int argc, const char* argv[])
{
	// parse params
	int tree_type = KDTREE_L;		// 0..6
	int ss_split_type = SPLIT_SS_HALF_LONGEST; // 0..2
	double accuracy = 1.0;			// acc (0.06 .. 4.0)
	DBox tree_box(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
	DiscreteSources cdm_sources(tree_box);

	cdm_sources.max_metric_len = tree_box.getDiameter() * ControlSpace3dAdaptive::param_max_diameter_ratio;

	cdm_sources.main_dir = argv[0];
	size_t found = cdm_sources.main_dir.find_last_of("/\\");
	if (found != std::string::npos)
		cdm_sources.main_dir = cdm_sources.main_dir.substr(0, found + 1);
	else
		cdm_sources.main_dir = "";

	DataHashTableKeyValue<string, int> h_tree_types(20, "");
	h_tree_types.insert("KDTREE_L", KDTREE_L);
	h_tree_types.insert("KDTREE_V", KDTREE_V);
	h_tree_types.insert("OCTREE_L", OCTREE_L);
	h_tree_types.insert("OCTREE_V", OCTREE_V);
	h_tree_types.insert("KDTREE_LI", KDTREE_LI);
	h_tree_types.insert("OCTREE_LB", OCTREE_LB);
	h_tree_types.insert("OCTREE_V_ORYG", OCTREE_V_ORYG);

	DataHashTableKeyValue<string, int> h_ss_split_types(20, "");
	h_ss_split_types.insert("SPLIT_HALF_LONGEST", SPLIT_SS_HALF_LONGEST);
	h_ss_split_types.insert("SPLIT_GOLDEN_RATIO", SPLIT_SS_GOLDEN_RATIO);
	h_ss_split_types.insert("SPLIT_MIN_LONGEST", SPLIT_SS_MIN_LONGEST);

	DataHashTableKeyValue<string, int> h_gradation_methods(20, "");
	h_gradation_methods.insert("GRADATION_NONE", GRADATION_NONE);
	h_gradation_methods.insert("GRADATION_LEAVES", GRADATION_LEAVES);
	h_gradation_methods.insert("GRADATION_VIC_PAIRS", GRADATION_VIC_PAIRS);
	h_gradation_methods.insert("GRADATION_VIC_NB", GRADATION_VIC_NB);
	h_gradation_methods.insert("GRADATION_VIC_NB_MULTIPASS", GRADATION_VIC_NB_MULTIPASS);

	cout << "Test-kd-ss params: ";
	for (int i = argi + 1; i < argc; i++) cout << " " << argv[i];
	cout << endl;

	const char* cdm_sources_fname = nullptr;
	if (++argi < argc && !parseIntArg(argv[argi], tree_type, &h_tree_types)) return -1;
	if (++argi < argc && !parseIntArg(argv[argi], ss_split_type, &h_ss_split_types)) return -1;
	if (++argi < argc && !parseDoubleArg(argv[argi], accuracy)) return -1;
	if (++argi < argc) cdm_sources_fname = argv[argi];
	if (++argi < argc && !parseDoubleArg(argv[argi], DiscreteSources::density_ratio)) return -1;
	if (++argi < argc && !parseIntArg(argv[argi], KdTree3dParams::paramGradationMethod, &h_gradation_methods)) return -1;
	if (++argi < argc && !parseDoubleArg(argv[argi], ControlSpace2dAdaptive::param_gradation_ratio)) return -1;
	if (++argi < argc && !parseTwoDoubleArgs(argv[argi], 
			cdm_sources.source_h_ratio_all, cdm_sources.source_h_ratio_normal)) return -1;

	if ((cdm_sources_fname != nullptr) && !loadCDMSources(cdm_sources_fname, cdm_sources)) return -1;

	cout << "Params parsed ok." << endl;
	cout << "---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---" << endl;

	switch (ss_split_type) {
	case SPLIT_SS_HALF_LONGEST:
		KdTree3dParams::m_ss_splitter = std::make_shared<Kd3dSSHalfLongestSplitter>();
		break;
	case SPLIT_SS_GOLDEN_RATIO:
		KdTree3dParams::m_ss_splitter = std::make_shared<Kd3dSSGoldenRatioLongestSplitter>();
		break;
	case SPLIT_SS_MIN_LONGEST:
		KdTree3dParams::m_ss_splitter = std::make_shared<Kd3dSSMinLongestSplitter>();
		break;
	default:
		cout << "Error: unknown splitter type: " << ss_split_type << endl;
		return -1;
	}

	DataVector<CS3dKdPtr> trees(TEST_TREE_RPT);
	CS3dKdPtr tree;

	for (int i = 0; i < TEST_TREE_RPT; i++) {
		switch (tree_type) {
		case KDTREE_L:
			tree = std::make_shared<ControlSpace3dKdTreeL>(tree_box);
			break;
		case KDTREE_LI:
			tree = std::make_shared<ControlSpace3dKdTreeLi>(tree_box);
			break;
		case KDTREE_V:
			tree = std::make_shared<ControlSpace3dKdTreeV>(tree_box);
			break;
		case OCTREE_L:
			tree = std::make_shared<ControlSpace3dKdOctreeL>(tree_box);
			break;
		case OCTREE_LB:
			tree = std::make_shared<ControlSpace3dKdOctreeLB>(tree_box);
			break;
		case OCTREE_V:
			tree = std::make_shared<ControlSpace3dKdOctreeV>(tree_box);
			break;
		case OCTREE_V_ORYG:
			tree = std::make_shared<ControlSpace3dOctree>(tree_box);
			break;
		default:
			cout << "Error: unknown tree type: " << tree_type << endl;
			return -1;
		}
		if (tree) trees.add(tree);
	}
	if (tree == nullptr) {
		cout << "Error: not implemented tree type: " << tree_type << endl;
		return -1;
	}

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Tree created OK, starting test...");

	runTestKdTreeSS(trees, tree_box, accuracy, cdm_sources);

	//delete tree;

	LOG4CPLUS_INFO(MeshLog::logger_mesh, "Tree tests finished OK.");

	return 0;
}

#endif // WITH_TESTS
//============================= TESTS END   ====================================
