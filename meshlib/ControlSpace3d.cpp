// ControlSpace3d.cpp: implementation of the ControlSpace2d class.
//
//////////////////////////////////////////////////////////////////////

#include "ControlSpace3d.h"
#include "MeshData.h"

int ControlSpace3d::param_control_type = MeshData::CONTROL_OCTREE_3D;

const ControlSpace3dAdaptive* ControlSpace3d::getAsAdaptive() const
{
	if (!isAdaptive()) return nullptr;
	return (const ControlSpace3dAdaptive*)this;
}

ControlSpace3dAdaptive* ControlSpace3d::getAsAdaptive()
{
	if (!isAdaptive()) return nullptr;
	return (ControlSpace3dAdaptive*)this;
}


double ControlNode3d::calculateMaxGradationNb(
	DataHashTableKeyValue< t_cn_pair, t_gr_step_pair > & gr_cn_pair,
	int step)
{
	auto cn0 = this;
	cn0->max_gradation_ratio = 1.0;
	std::shared_ptr<ControlNode3d> cn1;
	for (auto it = nbs->iterator(); it.valid(); it.moveNext()) {
		auto cn1 = it.item();
		auto cn_pair = (cn1.get() > cn0) ? std::make_pair(cn0, cn1.get()) : std::make_pair(cn1.get(), cn0);
		unsigned int index = 0;
		double gr = 0.0;
		if (gr_cn_pair.contains(cn_pair, index)) {
			gr = gr_cn_pair.slotValue(index).first;
		}else{
			gr = DMetric3d::getMetricGradation(
				cn0->control_data, cn1->control_data,
				cn0->coord - cn1->coord);
			gr_cn_pair.insert(cn_pair, std::make_pair(gr, step));
		}
		if (gr > cn0->max_gradation_ratio)
			cn0->max_gradation_ratio = gr;
	}
	return cn0->max_gradation_ratio;
}

double ControlNode3d::calculateMaxGradationNb(
	DataHashTableKeyValue< t_cn_pair, double> & gr_cn_pair)
{
	auto cn0 = this;
	cn0->max_gradation_ratio = 1.0;
	std::shared_ptr<ControlNode3d> cn1;
	for (auto it = nbs->iterator(); it.valid(); it.moveNext()) {
		auto cn1 = it.item();
		auto cn_pair = (cn1.get() > cn0) ? std::make_pair(cn0, cn1.get()) : std::make_pair(cn1.get(), cn0);
		double gr = gr_cn_pair.getValue(cn_pair, 0.0);
		if(gr == 0.0){
			gr = DMetric3d::getMetricGradation(
				cn0->control_data, cn1->control_data,
				cn0->coord - cn1->coord);
			gr_cn_pair.insert(cn_pair, gr);
		}
		if (gr > cn0->max_gradation_ratio)
			cn0->max_gradation_ratio = gr;
	}
	return cn0->max_gradation_ratio;
}

bool ControlNode3d::insertNbIfNew(const std::shared_ptr<ControlNode3d> & cn1)
{
	if(nbs == nullptr)
		nbs = std::make_unique<DataCompoundList<std::shared_ptr<ControlNode3d>>>(6); // max number of neighbours

	if (nbs->contains(cn1))
		return false;

	nbs->append(cn1);
	return true;
}
