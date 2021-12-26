/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dOctreeLB.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DOCTREELB_H__INCLUDED)
#define CONTROLSPACE3DOCTREELB_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace3dOctreeL.h"

/**
* This class implements a kd-tree based control space with nodes in leaves
*/
class KdOctElementLB : public KdOctElementL {
private:
	struct KdLeafDataLB {
		KdLeafDataLB(const DBox& _box, int _level)
			: box(_box), level(_level) {}
		DBox box;
		int level;
	};
public:
	KdOctElementLB(const DPoint3d& mid, const CDM3d & v, const KdNeighbors& kn,
		const DBox& _box, int _level);
public:
	void clear() override;
	std::shared_ptr<KdElement> clone() const override;
	/// Split this leaf (+ build relations with neighbours)
	int splitElementOct(DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
		const DBox & box, DBox * boxes = nullptr,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) override;
	int getTotalBytesExtra() const override;
public:
	std::unique_ptr<KdLeafDataLB> leaf_lb;
};

/**
* This class implements a kd-tree based control space with nodes in leaves
*/
class ControlSpace3dKdOctreeLB : public ControlSpace3dKdOctreeL
{
public:
	/// Standard contructor
	ControlSpace3dKdOctreeLB(const DBox& box) : ControlSpace3dKdOctreeL(box) { }
	//using ControlSpace3dKdOctreeL::ControlSpace3dKdOctreeL;
public:
	int createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices) override;
	// stat	
	string getTreeType() const override { return "octree-LB"; }
public:
	/// Returns the type of element
	virtual int getType() const override { return MeshData::CONTROL_OCTREE_3D_LB; }
};

#endif // !defined(CONTROLSPACE3DOCTREELB_H__INCLUDED)
