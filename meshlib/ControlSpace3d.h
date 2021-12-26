/////////////////////////////////////////////////////////////////////////////
// ControlSpace3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3D_H__INCLUDED)
#define CONTROLSPACE3D_H__INCLUDED

class ControlSpace3dAdaptive;

#include "DPoint.h"
#include "DMetric3d.h"
#include "MeshData.h"
#include "TagExtended.h"
#include "DataHashTable.h"
#include "DataList.h"

#define GRADATION_UNKNOWN  LARGE_NUMBER

class ControlNode3d;
typedef std::pair<ControlNode3d*, ControlNode3d*> t_cn_pair;
typedef std::pair<double, int> t_gr_step_pair;

namespace std {
	template <>
	struct hash<t_cn_pair>
	{
		std::size_t operator()(const t_cn_pair& cnp) const
		{
			std::size_t h = 0;
	
			h += (size_t)cnp.first;
			h += (h << 10);
			h ^= (h >> 6);

			h += (size_t)cnp.second;
			h += (h << 10);
			h ^= (h >> 6);

			h += (h << 3);
			h ^= (h >> 11);
			h += (h << 15);
			return h;
		}
	};
}

class ControlNode3d : public TagExtended {
public:
	ControlNode3d(const DPoint3d& pt) : coord(pt) {}
	ControlNode3d(const DPoint3d& pt, const ControlDataMatrix3d& cdm)
		: coord(pt), control_data(cdm) {}
	ControlNode3d(double x = 0.0, double y = 0.0, double z = 0.0) : coord(x, y, z) {}
	void reset() { w = 0.0; max_gradation_ratio = GRADATION_UNKNOWN; }
	bool gradationUnknown() const { return max_gradation_ratio == GRADATION_UNKNOWN; }
	void setGradationUnknown() { max_gradation_ratio = GRADATION_UNKNOWN;  }
	double calculateMaxGradationNb(DataHashTableKeyValue< t_cn_pair, t_gr_step_pair > & gr_cn_pair,
		int step);
	double calculateMaxGradationNb(DataHashTableKeyValue< t_cn_pair, double> & gr_cn_pair);
public:
	bool insertNbIfNew(const std::shared_ptr<ControlNode3d> & cn1);
public:
	DPoint3d coord;
	ControlDataMatrix3d control_data;
	union {
		double w = 0.0;
		int wi;
		void* wp;
	};
	double max_gradation_ratio = GRADATION_UNKNOWN;
	std::unique_ptr<DataCompoundList<std::shared_ptr<ControlNode3d>>> nbs; ///list of neighbours
};

/**
 * This class implements an abstract control space
 *  which provides the information about the sizing and/or 
 *  stretching of blocks within the meshed domain.
 */
class ControlSpace3d
{
public:
	/// Standard constructor
	ControlSpace3d() { }
	/// Destructor
	virtual ~ControlSpace3d() {}
public:
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const = 0;
	/// Returns the minimum size of tetrahedra for the given point 3D within the domain
	double getMinSpaceValue() const { return 0.05; }
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix3d getMetricAtPoint(const DPoint3d& pt) const = 0;
	/// If possible to adapt
	virtual bool isAdaptive() const { return false; }
	ControlSpace3dAdaptive* getAsAdaptive();
	const ControlSpace3dAdaptive* getAsAdaptive() const;
	/// Constrain point to the domain of this control space
	virtual DPoint3d fitInPoint(const DPoint3d& pt) const { return pt; }
	/// Whether the given point actually is withing the domain of this space
	virtual bool containsPoint(const DPoint3d& /* pt */) const { return true; }
	/// Invoke function for all control nodes of this space (read-only)
	virtual void forEachControlNode(const std::function<void(const ControlNode3d& node)>& /* fg */) const {}
	/// Invoke function for all control nodes of this space
	virtual void forEachControlNode(const std::function<void(ControlNode3d& node)>& /* fg */) {}
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const { return 0; }
public:
	/// Which type of control space should be used
	static int param_control_type;
};

typedef std::shared_ptr<ControlSpace3d> CS3dPtr;
typedef std::shared_ptr<const ControlSpace3d> CS3dConstPtr;

#endif // !defined(CONTROLSPACE3D_H__INCLUDED)
