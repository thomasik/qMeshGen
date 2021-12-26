/////////////////////////////////////////////////////////////////////////////
// ControlSpace2dKdTreeL.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACEKDTREEL_H__INCLUDED)
#define CONTROLSPACEKDTREEL_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace2dAdaptive.h"

/**
 * This class implements a kd-tree based control space with nodes in leaves
 */
class ControlSpace2dKdTreeL : public ControlSpace2dAdaptive
{
public:
	/// Standard contructor
	ControlSpace2dKdTreeL(SurfaceConstPtr surface, const DRect& box)
		: ControlSpace2dAdaptive(surface, box) { }
	~ControlSpace2dKdTreeL() { clear(); }
private:
	ControlSpace2dKdTreeL(const ControlSpace2dKdTreeL& ) = delete;
	ControlSpace2dKdTreeL& operator=(const ControlSpace2dKdTreeL& ) = delete;
public:
	int adaptToField(const std::function<CDM2d(const DPoint2d & pt)> & f, double max_error);
	void clear();
// stat	
	int getElementCount() const { return m_top ? m_top->getElementCount() : 0; }
	int getValueCount() const { return m_top ? m_top->getValueCount() : 0; }
	int getMaxDepth() const { return m_top ? m_top->getMaxDepth() : 0; }
	int getMaxBalance() const;
	int getTotalBytes() const;
public:
	/// Returns the type of element
	virtual int getType() const override{ return MeshData::CONTROL_KDTREE_L; }
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt) const override;
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const override { return getValueCount(); } 
	/// Invoke function for all control nodes of this space (read-only)
	virtual void forEachControlNode(const std::function<void(const ControlNode2d& node)>& fg) const override;
	/// Invoke function for all control nodes of this space
	virtual void forEachControlNode(const std::function<void(ControlNode2d& node)>& fg) override;

	/// Adds new information (stretching and lengths) for some point within the domain
	virtual void addControlNode(const ControlNode2d& node) override { assert(false); }
	/// Interpolates (initializes) the control space basing on the previously given sources
	virtual bool interpolate() override { assert(false); return true; }
	/// Refines control space at the given point
	virtual bool setMinControl(const DPoint2d& pt, const ControlDataMatrix2d& cdm, bool min_value_set = true) override;
	/// Returns the "screenshot" of this control space for visualization
	//virtual MeshViewSet* getViewSet(MeshViewSet* set = nullptr, bool /* with_points */ = false) const override {
	//	assert(false); return set; }
	/// Log basic information about this control space
	virtual void logDescription() const override { assert(false); }
	/// Smoothen variance of metric within the control space (returns true if any change)
	virtual bool smoothen() override { assert(false); return false; }
	/// optimize for size (lower memory requirements) - but also from now readonly !
	//virtual void compact() override { assert(false); m_initialized = 10; }
	/// Initializes the control space basing on the curvature of surface
	virtual void setSurfaceCurvatureControlData() override { assert(false); }
	/// Refines control space to parameterization variance
	virtual void adaptToParameterization() override { assert(false); }
	/// Return interpolated value of extended tag data from control nodes at some given point;
	virtual double interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const override { assert(false); return 0.0;  }
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint2d& pt) const override { assert(false); return 1.0; };
private:
	class KdElement {
	public:
		KdElement(const DPoint2d& mid, const CDM2d & v) : m_data(mid, v) { }
		int adaptToField(const std::function<CDM2d(const DPoint2d & pt)> & f,
			const DRect & box, double max_error, int level = 0);
		const ControlDataMatrix2d& getValue(const DPoint2d& pt, const DRect& box) const;
		const ControlDataMatrix2d& getValue(const DPoint2d& pt) const;
		inline ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt, const DRect& box) const { return getValue(pt); }
		void clear();
		/// Returns smallest KdTree leaf for this point
		KdElement* findLastLeaf(const DPoint2d& pt, DRect & box, int & level);
		/// Returns smallest KdTree leaf for this point
		KdElement* findLastLeaf(const DPoint2d& pt, DRect & box);
		/// Returns smallest KdTree leaf for this point
		KdElement* findLastLeaf(const DPoint2d& pt);
		/// Adapt treeleaf to new source point and set minimum value
		bool adaptAndSetControlPoint(const ControlNode2d& qv, SurfaceConstPtr surface,
			DRect & box, int & level, bool min_value_set = true);
		/// Split this leaf (+ build relations with neighbours)
		bool split(DataVector<ControlNode2d*>& split_vertices, 
			const DRect & box, Axis split_axis, double split_value,
			DRect & box0, DRect & box1);
	public:
		int getElementCount() const;
		int getValueCount() const;
		int getMaxDepth() const;
		int getMaxDepthAndBalance(int & max_balance) const;
		/// Invoke function for all control nodes of this space (read-only)
		void forEachControlNodeConst(const std::function<void(const ControlNode2d& node)>& fg) const;
		/// Invoke function for all control nodes of this space
		void forEachControlNode(const std::function<void(ControlNode2d& node)>& fg);
	public:
		bool is_leaf = true;
		union KdData {
			KdData(const DPoint2d& mid, const CDM2d & v) : value(mid, v) {}
			~KdData() {}
			struct KdSplitElement {
				Axis		axis = Axis::X;
				double		value = 0.0;
				KdElement * element_0 = nullptr;
				KdElement * element_1 = nullptr;
			} split;
			ControlNode2d value;
		} m_data;
	};
protected:
	/// Returns smallest KdTree leaf for this point
	KdElement* findLastLeaf(const DPoint2d& pt, DRect & box, int & level) const;
	/// Returns smallest KdTree leaf for this point
	KdElement* findLastLeaf(const DPoint2d& pt, DRect & box) const;
	/// Returns smallest KdTree leaf for this point
	KdElement* findLastLeaf(const DPoint2d& pt) const;
private:
	KdElement* m_top = nullptr;
};

#endif // !defined(CONTROLSPACEKDTREEL_H__INCLUDED)
