/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dOctreeL.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DOCTREEL_H__INCLUDED)
#define CONTROLSPACE3DOCTREEL_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace3dKdTree.h"

/**
 * This class implements a kd-tree based control space with nodes in leaves
 */
class KdOctElementL : public KdOctElement{
private:
	struct KdLeafData {
		KdLeafData(const DPoint3d& mid, const CDM3d & v, const KdNeighbors& _kn)
			: cnode( std::make_shared<ControlNode3d>(mid, v)), kn(_kn) {}
		std::shared_ptr<ControlNode3d> cnode;
		KdNeighbors kn;
	};
public:
	KdOctElementL(const DPoint3d& mid, const CDM3d & v, const KdNeighbors& kn);
public:
	int adaptToField(const std::function<CDM3d(const DPoint3d & pt)> & f,
		const DBox & box, double max_error, int level = 0,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) override;
	ControlDataMatrix3d getCDM(const DPoint3d& pt, const DBox& box) const override;
	ControlDataMatrix3d getCDM(const DPoint3d& pt) const override;
	bool setMinimumCDM(const CDM3d & cdm);
	std::shared_ptr<ControlNode3d> getCN() override;
	std::shared_ptr<const ControlNode3d> getCN() const override;
	void clear() override;
	std::shared_ptr<KdElement> clone() const override;
	/// Split this leaf (+ build relations with neighbours)
	int splitElementOct(DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
		const DBox & box, DBox * boxes = nullptr,
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
class ControlSpace3dKdOctreeL : public ControlSpace3dKdTree
{
public:
	/// Standard contructor
	ControlSpace3dKdOctreeL(const DBox& box) : ControlSpace3dKdTree(box) { }
public:
	int createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices) override;
	// stat	
	string getTreeType() const override { return "octree-L"; }
public:
	/// Returns the type of element
	virtual int getType() const override { return MeshData::CONTROL_OCTREE_3D_L; }
};


#endif // !defined(CONTROLSPACE3DOCTREEL_H__INCLUDED)
