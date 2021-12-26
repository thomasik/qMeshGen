/////////////////////////////////////////////////////////////////////////////
// ControlSpace2dKdTreeV.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACEKDTREEV_H__INCLUDED)
#define CONTROLSPACEKDTREEV_H__INCLUDED

#include "MeshData.h"
#include "ControlSpace2dAdaptive.h"

/**
 * This class implements a kd-tree based control space with metric in vertices
 */
class ControlSpace2dKdTreeV : public ControlSpace2dAdaptive
{
public:
	/// Standard contructor
	ControlSpace2dKdTreeV(SurfaceConstPtr surface, const DRect& box)
		: ControlSpace2dAdaptive(surface, box) { }
	~ControlSpace2dKdTreeV() { clear(); }
private:
	ControlSpace2dKdTreeV(const ControlSpace2dKdTreeV&) = delete;
	ControlSpace2dKdTreeV& operator=(const ControlSpace2dKdTreeV&) = delete;
public:
	/// Returns the type of element
	virtual int getType() const override{ return MeshData::CONTROL_KDTREE_V; }
public:
	int adaptToField(const std::function<CDM2d(const DPoint2d & pt)> & f, double max_error);
	void clear();
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
	//	assert(false); return set;
	//}
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
	virtual double interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const override { assert(false); return 0.0; }
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint2d& pt) const override { assert(false); return 1.0; };
public: // stat
	int getElementCount() const { return m_top ? m_top->getElementCount() : 0; }
	int getValueCount() const { return m_values.countInt(); }
	int getMaxDepth() const { return m_top ? m_top->getMaxDepth() : 0; }
	int getMaxBalance() const;
	int getTotalBytes() const;
private:
	class KdElement {
	public:
		struct KdVertices {
			ControlNode2d* vx_sw = nullptr;
			ControlNode2d* vx_se = nullptr;
			ControlNode2d* vx_nw = nullptr;
			ControlNode2d* vx_ne = nullptr;
		};
	public:
		KdElement(const KdVertices& kv) : m_data(kv) { }
		int adaptToField(const std::function<CDM2d(const DPoint2d & pt)> & f,
			const DRect & box, double max_error, DataVector<ControlNode2d> & kvalues, int level = 0);
		ControlDataMatrix2d getValue(const DPoint2d& pt, const DRect& box) const;
		inline ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt, const DRect& box) const { return getValue(pt, box); }
		void clear();
		/// Returns smallest KdTree leaf for this point
		KdElement* findLastLeaf(const DPoint2d& pt, DRect & box, int & level);
		/// Adapt treeleaf to new source point and set minimum value
		bool adaptAndSetControlPoint(DataVector<ControlNode2d>& tree_vertices,
			const ControlNode2d& qv, SurfaceConstPtr surface,
			DRect & box, int & level, bool min_value_set = true);
		/// Split this leaf (+ build relations with neighbours)
		bool split(DataVector<ControlNode2d>& tree_vertices, DataVector<ControlNode2d*>& split_vertices, 
			const DRect & box, Axis split_axis, double split_value,
			DRect & box0, DRect & box1);
	public:
		int getElementCount() const;
		int getMaxDepth() const;
		int getMaxDepthAndBalance(int & max_balance) const;
	public:
		bool is_leaf = true;
		union KdData {
			KdData(const KdVertices & kv) : values(kv) {}
			struct KdSplitElement {
				Axis		axis = Axis::X;
				double		value = 0.0;
				KdElement * element_0 = nullptr;
				KdElement * element_1 = nullptr;
			} split;
			KdVertices values;
		} m_data;
	};
protected:
	/// Returns smallest KdTree leaf for this point
	KdElement* findLastLeaf(const DPoint2d& pt, DRect & box, int & level) const;
	/// Returns smallest KdTree leaf for this point
	KdElement* findLastLeaf(const DPoint2d& pt, DRect & box) const;
private:
	KdElement* m_top = nullptr;
	DataVector<ControlNode2d> m_values;
};


#endif // !defined(CONTROLSPACEKDTREEV_H__INCLUDED)
