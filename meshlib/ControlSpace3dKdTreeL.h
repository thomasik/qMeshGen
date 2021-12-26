/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dKdTreeL.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DKDTREEL_H__INCLUDED)
#define CONTROLSPACE3DKDTREEL_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace3dKdTree.h"

class KdElementL : public KdElement {
public:
	struct KdLeafData {
		KdLeafData(const DPoint3d& mid, const CDM3d & v, const KdNeighbors& _kn)
			: cnode( std::make_shared<ControlNode3d>(mid, v) ), kn(_kn) {}
		std::shared_ptr<ControlNode3d> cnode;
		KdNeighbors kn;
	};
public:
	KdElementL(const DPoint3d& mid, const CDM3d & v, const KdNeighbors& kn);
	KdElementL(const DPoint3d& mid, KdLeafData* kd_data);
public:
	int adaptToField(const std::function<CDM3d(const DPoint3d & pt)> & f,
		const DBox & box, double max_error, int level = 0,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) override;
	ControlDataMatrix3d getCDM(const DPoint3d& pt, const DBox& box) const override;
	ControlDataMatrix3d getCDM(const DPoint3d& pt) const override;
	bool setMinimumCDM(const CDM3d& cdm) override;
	std::shared_ptr<ControlNode3d> getCN() override;
	std::shared_ptr<const ControlNode3d> getCN() const override;
	void clear() override;
	std::shared_ptr<KdElement> clone() const override;
	/// Split this leaf (+ build relations with neighbours)
	int splitElementKd(DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
		const DBox & box, Axis split_axis, double split_value,
		DBox * box0 = nullptr, DBox * box1 = nullptr,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) override;
	void setKN(const KdNeighbors & kn) override;
	void setKN(int fi, std::shared_ptr<KdElement> const & kde) override;
	int smoothenMetricAtNodes() override;
	void gatherControlNodesNeighbors() override;
	int getTotalBytesExtra() const override;
public:
	std::unique_ptr<KdLeafData> leaf;
};

/**
 * This class implements a kd-tree based control space with nodes in leaves
 */
class ControlSpace3dKdTreeL : public ControlSpace3dKdTree
{
public:
	/// Standard contructor
	ControlSpace3dKdTreeL(const DBox& box) : ControlSpace3dKdTree(box) { }
public:
	int createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices) override;
	string getTreeType() const override { return "kdtree-L"; }
public:
	/// Returns the type of element
	virtual int getType() const override { return MeshData::CONTROL_KDTREE_3D_L; }
};

#endif // !defined(CONTROLSPACE3DKDTREEL_H__INCLUDED)
