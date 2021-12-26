/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dKdTreeLi.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DKDTREELI_H__INCLUDED)
#define CONTROLSPACE3DKDTREELI_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace3dKdTreeL.h"

class KdElementLi : public KdElementL {
public:
	using KdElementL::KdElementL;
public:
	ControlDataMatrix3d getCDM(const DPoint3d& pt, const DBox& box) const override;
	ControlDataMatrix3d getCDM(const DPoint3d& pt) const override;
	std::shared_ptr<KdElement> clone() const override;
};

/**
 * This class implements a kd-tree based control space with nodes in leaves
 */
class ControlSpace3dKdTreeLi : public ControlSpace3dKdTreeL
{
public:
	/// Standard contructor
	ControlSpace3dKdTreeLi(const DBox& box) : ControlSpace3dKdTreeL(box) { }
	//using ControlSpace3dKdTreeL::ControlSpace3dKdTreeL;
public:
	int adaptToField(const std::function<CDM3d(const DPoint3d & pt)> & f, double max_error) override;
	int createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices) override;
	string getTreeType() const override { return "kdtree-Li"; }
public:
	virtual int getType() const override{ return MeshData::CONTROL_KDTREE_3D_LI; }
};

#endif // !defined(CONTROLSPACE3DKDTREELLI_H__INCLUDED)
