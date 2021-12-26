/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dKdTreeV.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DKDTREEV_H__INCLUDED)
#define CONTROLSPACE3DKDTREEV_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace3dKdTree.h"

class KdElementV : public KdElement {
public:
	struct KdLeafData {
		KdLeafData(const KdVertices&_kv, const KdNeighbors& _kn)
			: kv(_kv), kn(_kn) {}
		KdVertices kv;
		KdNeighbors kn;
	};
public:
	KdElementV(const KdVertices & kv, const KdNeighbors& kn);
public:
	int adaptToField(const std::function<CDM3d(const DPoint3d & pt)> & f,
		const DBox & box, double max_error, int level = 0,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) override;
	ControlDataMatrix3d getCDM(const DPoint3d& pt, const DBox& box) const;
	ControlDataMatrix3d getCDM(const DPoint3d& pt) const;
	bool setMinimumCDM(const CDM3d& cdm) override;
	void setKV(int vi, std::shared_ptr<ControlNode3d> const & cn) override;
	void setKN(int fi, std::shared_ptr<KdElement> const & kde) override;
	void clear() override;	/// Adapt treeleaf to new source point and set minimum value
	/// Split this leaf (+ build relations with neighbours)
	int splitElementKd(
		DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
		const DBox & box, Axis split_axis, double split_value,
		DBox * box0 = nullptr, DBox * box1 = nullptr, 
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) override;
	int smoothenMetricAtNodes() override;
	std::shared_ptr<KdElement> clone() const override;
	void gatherControlNodesNeighbors() override;
	int getTotalBytesExtra() const override;
	std::shared_ptr<ControlNode3d> getNearestVertexRef(const DPoint3d & pt, const DBox & box);
	std::shared_ptr<ControlNode3d> getNearestControlNodeForEdge(const DBox& box, const DPoint3d& edge_pt, int ei, double * ptr_min_dist2 = nullptr);
	int replaceControlNodes(const DataHashTableKeyValue< 
		std::shared_ptr<ControlNode3d>, std::shared_ptr<ControlNode3d> > & hcnodes) override;
public:
	std::unique_ptr<KdLeafData> leaf;
};


/**
 * This class implements a kd-tree based control space with metric in vertices
 */
class ControlSpace3dKdTreeV : public ControlSpace3dKdTree
{
public:
	/// Standard contructor
	ControlSpace3dKdTreeV(const DBox& box);
public:
	/// Returns the type of element
	virtual int getType() const override{ return MeshData::CONTROL_KDTREE_3D_V; }
public:
	int createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices) override;
	int fixCloseControlNodes() override;
public: // stat
	string getTreeType() const override { return "kdtree-V"; }
};

#endif // !defined(CONTROLSPACE3DKDTREEV_H__INCLUDED)
