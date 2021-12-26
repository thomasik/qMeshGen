/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dKdTree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2017-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DKDTREE_H__INCLUDED)
#define CONTROLSPACE3DKDTREE_H__INCLUDED

#include <vector>
#include <memory>

#include "MeshData.h"
#include "ControlSpace3dAdaptive.h"
#include "DataKdSplit3dContField.h"
#include "DataKdSplit3dDiscField.h"
#include "DataKdSplit3dSingleSource.h"
#include "DataHashTable.h"
#include "DataVector.h"

enum GRADATION_METHOD { 
	GRADATION_NONE, GRADATION_LEAVES, 
	GRADATION_VIC_PAIRS, GRADATION_VIC_NB, 
	GRADATION_VIC_NB_MULTIPASS 
};

typedef DataHashTable<std::shared_ptr<ControlNode3d>> SetCN;
typedef DataHashTableKeyValue< std::shared_ptr<ControlNode3d>, std::shared_ptr<SetCN> > MapSetCN;

class KdTree3dParams
{
public:
	static int paramEstimateErrorSamplingSize;
	static int paramSplitGradientSamplingSize;
	static int paramMaxLevel;
	static int paramGradationMethod;
	static double paramMetricLenRatio;
public:
	static std::shared_ptr<Kd3dSplitterContField>    m_cf_splitter;
	static std::shared_ptr<Kd3dSplitterDiscField>    m_df_splitter;
	static std::shared_ptr<Kd3dSplitterSingleSource> m_ss_splitter;
};

class KdOctree3dParams
{
public:
	static int paramEstimateErrorSamplingSize;
	static int paramMaxLevel;
	static int paramMaxBalance;
};

class KdElement;

class KdSplitElement {
public:
	KdSplitElement(int n) : elements(n) {}
	virtual ~KdSplitElement() {}
public:
	virtual std::shared_ptr<KdElement> selectChild(const DPoint3d& pt, DBox * box = nullptr) const = 0;
	virtual int getTotalBytes() const { return sizeof(elements) + elements.countInt() * 
		sizeof(std::shared_ptr<KdElement>); }
public:
	DataVector<std::shared_ptr<KdElement>> elements;
};

class KdTreeSplitElement : public KdSplitElement{
public:
	KdTreeSplitElement(Axis _axis = Axis::X, double _value = 0.0) 
		: KdSplitElement(2), axis(_axis), value(_value) {}
	KdTreeSplitElement(Axis _axis, double _value, 
		std::shared_ptr<KdElement> kde0, 
		std::shared_ptr<KdElement> kde1);
public:
	std::shared_ptr<KdElement> selectChild(const DPoint3d& pt, DBox * box = nullptr) const override;
	int getTotalBytes() const override { return KdSplitElement::getTotalBytes() + sizeof(Axis) + sizeof(double); }
public:
	Axis		axis = Axis::X;
	double		value = 0.0;
};

class KdOctreeSplitElement : public KdSplitElement{
public:
	KdOctreeSplitElement(const DPoint3d& mid) : KdSplitElement(8), middle(mid) {}
public:
	std::shared_ptr<KdElement> selectChild(const DPoint3d& pt, DBox * box = nullptr) const override;
	int getTotalBytes() const override{ return KdSplitElement::getTotalBytes() + sizeof(DPoint3d); }
public:
	DPoint3d middle;
};

struct KdNeighbors {
public:
	KdNeighbors(std::shared_ptr<KdElement> _low, 
		std::shared_ptr<KdElement> _high, 
		std::shared_ptr<KdElement> _south,
		std::shared_ptr<KdElement> _north, 
		std::shared_ptr<KdElement> _east, 
		std::shared_ptr<KdElement> _west)
	{
		elements[FC3D_LOW] = _low;		elements[FC3D_HIGH] = _high;
		elements[FC3D_SOUTH] = _south;	elements[FC3D_NORTH] = _north;
		elements[FC3D_EAST] = _east;	elements[FC3D_WEST] = _west;
	}
	KdNeighbors() = default;
	std::shared_ptr<KdElement>& operator[](int i) { assert(i < FC3D_COUNT); return elements[i]; }
	std::shared_ptr<KdElement> const & operator[](int i) const { assert(i < FC3D_COUNT); return elements[i]; }
public:
	std::shared_ptr<KdElement> elements[FC3D_COUNT] = { nullptr };
};

struct KdVertices {
public:
	KdVertices(
		std::shared_ptr<ControlNode3d> _lsw, 
		std::shared_ptr<ControlNode3d> _lse,
		std::shared_ptr<ControlNode3d> _lnw, 
		std::shared_ptr<ControlNode3d> _lne,
		std::shared_ptr<ControlNode3d> _hsw, 
		std::shared_ptr<ControlNode3d> _hse,
		std::shared_ptr<ControlNode3d> _hnw, 
		std::shared_ptr<ControlNode3d> _hne)
	{
		vertices[VX3D_LSW] = _lsw;	vertices[VX3D_LSE] = _lse;
		vertices[VX3D_LNW] = _lnw;	vertices[VX3D_LNE] = _lne;
		vertices[VX3D_HSW] = _hsw;	vertices[VX3D_HSE] = _hse;
		vertices[VX3D_HNW] = _hnw;	vertices[VX3D_HNE] = _hne;
	}
	KdVertices() = default;
	std::shared_ptr<ControlNode3d>& operator[](int i) { assert(i < VX3D_COUNT); return vertices[i]; }
	std::shared_ptr<ControlNode3d> const & operator[](int i) const { assert(i < VX3D_COUNT); return vertices[i]; }
public:
	std::shared_ptr<ControlNode3d> vertices[VX3D_COUNT] = { nullptr };
};

class KdElement {
public:
	KdElement() = default;
	virtual ~KdElement();
public:
	virtual CDM3d getCDM(const DPoint3d& pt) const = 0;
	virtual CDM3d getCDM(const DPoint3d& pt, const DBox& box) const = 0;
	virtual bool setMinimumCDM(const CDM3d& cdm) = 0;
	virtual std::shared_ptr<ControlNode3d> getCN() { return nullptr; }
	virtual std::shared_ptr<const ControlNode3d> getCN() const { return nullptr; }
	virtual int adaptToField(
		const std::function<CDM3d(const DPoint3d & pt)> & f,
		const DBox & box, double max_error, int level = 0,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) = 0;
	virtual bool adaptAndSetControlPoint(
		std::shared_ptr<ControlNode3d> qv, DBox & box, int & level, bool min_value_set = true,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr);
	virtual int splitElementKd(DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
		const DBox & box, Axis split_axis, double split_value,
		DBox * box0 = nullptr, DBox * box1 = nullptr,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) { assert(false); return 0; }
	virtual int splitElementOct(DataVector<std::shared_ptr<ControlNode3d>>& split_vertices,
		const DBox & box, DBox * boxes = nullptr,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) { assert(false); return 0; }
	virtual int smoothenMetricAtNodes() = 0;
	virtual void setKN(const KdNeighbors & kn) { assert(false);  }
	virtual void setKN(int fi, std::shared_ptr<KdElement> const & kde) { assert(false);  }
	virtual void setKV(const KdVertices& kv) { assert(false); }
	virtual void setKV(int vi, std::shared_ptr<ControlNode3d> const & cn) { assert(false); }
	virtual std::shared_ptr<KdElement> clone() const = 0;
	virtual void gatherControlNodesNeighbors() = 0;
	virtual int getTotalBytesExtra() const;
	virtual int replaceControlNodes(const DataHashTableKeyValue< 
		std::shared_ptr<ControlNode3d>, std::shared_ptr<ControlNode3d> > & hcnodes) { return 0;  }
	int getTotalBytes() const;
public:
	inline bool isLeaf() const { return split == nullptr; }
	virtual void clear();
	KdElement* getNearestLeaf(const DPoint3d & pt, DBox * box = nullptr, int * level = nullptr);
	/// Returns smallest KdTree leaf for this point
	KdElement* findLastLeaf(const DPoint3d& pt, DBox * box = nullptr, int * level = nullptr);
	void forEachControlNodeConst(const std::function<void(std::shared_ptr<const ControlNode3d> node)>& fg) const;
	void forEachControlNode(const std::function<void(std::shared_ptr<ControlNode3d> node)>& fg);
	double estimateError(const DBox & box, const std::function<CDM3d(const DPoint3d&pt)>& f); 
	int getElementCount() const;
	int getValueCount() const;
	int getMaxDepth() const;
	int getMaxDepthAndBalance(int & max_balance) const;
public:
	KdSplitElement* split = nullptr;
};

class KdOctElement : public KdElement
{
public:
	virtual bool adaptAndSetControlPoint(
		std::shared_ptr<ControlNode3d> qv, DBox & box, int & level, bool min_value_set = true,
		DataVector<std::shared_ptr<ControlNode3d>> *kvalues = nullptr) override;
};

 /**
 * This class implements a general kd-tree based control space
 */
class ControlSpace3dKdTree : public ControlSpace3dAdaptive
{
public:
	/// Standard contructor
	ControlSpace3dKdTree(const DBox& box) : ControlSpace3dAdaptive(box) { }
	~ControlSpace3dKdTree() { clear(); }
private:
	ControlSpace3dKdTree(const ControlSpace3dKdTree&) = delete;
	ControlSpace3dKdTree& operator=(const ControlSpace3dKdTree&) = delete;
public:
	virtual int adaptToField(const std::function<CDM3d(const DPoint3d&pt)>& f, double max_error);
	/// Refines control space at the given point
	virtual bool setMinControl(const DPoint3d & pt, const ControlDataMatrix3d & cdm, bool min_value_set);
	/// Initializes the control space with given global metric
	virtual void setGlobalMetric(const ControlDataMatrix3d & cdm);	// stat	
	virtual int getElementCount() const { return m_top ? m_top->getElementCount() : 0; }
	virtual int getValueCount() const;
	virtual int getMaxDepth() const { return m_top ? m_top->getMaxDepth() : 0; }
	virtual int getMaxBalance() const;
	virtual int getTotalBytes() const;
	virtual int createTopElement(DataVector<std::shared_ptr<ControlNode3d>>& new_vertices) { assert(false); return 0; }
public:
	/// Returns descriptive label
	virtual string getLabel() const;
	virtual string getTreeType() const = 0;
	/// Returns the type of element
	virtual int getType() const override { return MeshData::CONTROL_KDTREE_3D; }
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const override { return getValueCount(); }
	/// Adds new information (stretching and lengths) for some point within the domain
	virtual void addControlNode(const ControlNode3d& /* node */) override { assert(false); }
	/// Interpolates (initializes) the control space basing on the previously given sources
	virtual bool interpolate() override { assert(false); return true; }
	/// Returns the "screenshot" of this control space for visualization
	virtual MeshViewSet* getViewSet(MeshViewSet* set = nullptr, bool /* with_points */ = false) const override {
		assert(false); return set;
	}
	/// Log basic information about this control space
	virtual void logDescription() const override { assert(false); }
	/// optimize for size (lower memory requirements) - but also from now readonly !
	// virtual void compact() override { assert(false); m_initialized = 10; }
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint3d& /* pt */) const override { assert(false); return 1.0; }
	/// Return interpolated value of extended tag data from control nodes at some given point;
	virtual double interpolateDoubleTag(const DPoint3d& /* pt */, TagExtended::TagType /* type */) const override { assert(false); return 0.0; }
	virtual int fixCloseControlNodes() { return 0;  }
public:
	/// Smoothen variance of metric within the control space (returns true if any change)
	bool smoothen() override { return smoothenViaLeaves();  }
	bool smoothenViaLeaves();
	bool smoothenViaVicinity() { return smoothenViaVicinityPairs(); }
	bool smoothenViaVicinityPairs();
	bool smoothenViaVicinityNb();
	bool smoothenViaVicinityNbMultipass();

	void createControlNodesNeighbors();
	void clearControlNodesNeighbors();

	ControlDataMatrix3d getMetricAtPoint(const DPoint3d & pt) const override;
	void clear();
	KdElement * findLastLeaf(const DPoint3d & pt, DBox * box = nullptr, int * level = nullptr) const;
	void forEachControlNode(const std::function<void(std::shared_ptr<const ControlNode3d> node)>& fg) const;
	void forEachControlNode(const std::function<void(std::shared_ptr<ControlNode3d> node)>& fg);

	void dumpCNodes(const char* fname) const;
public:
	void showMetricLength() const;
public:
	static int checkKdTreeSS(int argi, int argc, const char* argv[]);
	static int checkKdTreeCF(int argi, int argc, const char* argv[]);
protected:
	std::shared_ptr<KdElement> m_top;
	std::shared_ptr<DataVector<std::shared_ptr<ControlNode3d>>> m_values;
};
typedef std::shared_ptr<ControlSpace3dKdTree> CS3dKdPtr;

#endif // !defined(CONTROLSPACE3DKDTREE_H__INCLUDED)
