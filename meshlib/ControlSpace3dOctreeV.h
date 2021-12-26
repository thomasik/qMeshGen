/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dKdOctreeV.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DOCTREEV_H__INCLUDED)
#define CONTROLSPACE3DOCTREEV_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace3dKdTree.h"

class KdOctElementV : public KdOctElement {
private:
	struct KdLeafData {
		KdLeafData(const KdNeighbors& _kn, const KdVertices& _kv)
			: kn(_kn), kv(_kv) {}
		KdNeighbors kn;
		KdVertices kv;
	};
public:
	KdOctElementV(const KdNeighbors& _kn, const KdVertices& _kv);
public:
	int adaptToField(const std::function<CDM3d(const DPoint3d & pt)> & f,
		const DBox & box, double max_error, int level = 0,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) override;
	ControlDataMatrix3d getCDM(const DPoint3d& pt, const DBox& box) const;
	ControlDataMatrix3d getCDM(const DPoint3d& pt) const;
	bool setMinimumCDM(const CDM3d & cdm);	
	void setKV(int vi, std::shared_ptr<ControlNode3d> const & cn) override;
	void setKV(const KdVertices& kv) override;
	void setKN(const KdNeighbors & kn) override;
	void setKN(int fi, std::shared_ptr<KdElement> const & kde) override;
	void clear() override;	/// Adapt treeleaf to new source point and set minimum value
	/// Split this leaf (+ build relations with neighbours)
	int splitElementOct(DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
		const DBox & box, DBox * boxes,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) override;
	int smoothenMetricAtNodes() override;
	std::shared_ptr<KdElement> clone() const override;
	std::shared_ptr<ControlNode3d> getNearestVertexRef(const DPoint3d & pt, const DBox & box);
	void gatherControlNodesNeighbors() override;
	int getTotalBytesExtra() const override;
public:
	std::unique_ptr<KdLeafData> leaf;
};


/**
* This class implements a kd-tree based control space with metric in vertices
*/
class ControlSpace3dKdOctreeV : public ControlSpace3dKdTree
{
public:
	/// Standard contructor
	ControlSpace3dKdOctreeV(const DBox& box);
public:
	/// Returns the type of element
	virtual int getType() const override { return MeshData::CONTROL_OCTREE_3D_V; }
public:
	int createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices) override;
public: // stat
	string getTreeType() const override { return "octree-V"; }
};
#endif // !defined(CONTROLSPACE3DKDTREEV_H__INCLUDED)
