/////////////////////////////////////////////////////////////////////////////
// ControlSpace2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE2D_H__INCLUDED)
#define CONTROLSPACE2D_H__INCLUDED

class Curve2dParametric;
class MeshContainer2d;
class MeshViewSet;
class ReparameterizationData;
class ControlSpace2dAdaptive;

#include <memory>

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"
#include "DMetric2d.h"
#include "MeshData.h"
#include "DEquation.h"
#include "SurfaceParametric.h"
#include "ControlSpace2d.h"

class ControlNode2d : public TagExtended {
public:
	ControlNode2d(const DPoint2d& pt) : coord(pt), w(0.0), max_gradation_ratio(1.0) {}
	ControlNode2d(const DPoint2d& pt, const ControlDataMatrix2d& cdm)
		: coord(pt), control_data(cdm), w(-1.0), max_gradation_ratio(1.0) {}
	ControlNode2d(const DPoint2d& pt, const ControlDataMatrix2d& cdm, const ControlDataMatrix2d& pdm)
		: coord(pt), control_data(cdm), param_data(pdm), w(-1.0), max_gradation_ratio(1.0) {}
	ControlNode2d(double x = 0.0, double y = 0.0) : coord(x, y), w(0.0), max_gradation_ratio(1.0) {}
public:
	DPoint2d coord;
	ControlDataMatrix2d control_data;
	ControlDataMatrix2d param_data;
	double w;
	double max_gradation_ratio;
};

/**
 * This class implements an abstract control space
 *  which provides the information about the sizing and/or 
 *  stretching of elements within the meshed domain.
 */
class ControlSpace2d
{
public:
	/// Standard constructor
	ControlSpace2d(SurfaceConstPtr surface) : 
		base_surface(surface) { }
	/// Destructor
	virtual ~ControlSpace2d() { }
public:
	/// Returns stats of comparison of metric in two control spaces
//	virtual MeshData::StatData compareMetricWith(CS2dPtr other_control_space)
//		{ return MeshData::StatData(); }
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const = 0;
	/// If possible to adapt
	virtual bool isAdaptive() const { return false; }
	ControlSpace2dAdaptive* getAsAdaptive();
	const ControlSpace2dAdaptive* getAsAdaptive() const;
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt) const = 0;
	/// Get transformation and parameterization matrices at the given point
	virtual ControlDataMatrix2d getMetricAndParameterizationAtPoint(const DPoint2d& pt, ControlDataMatrix2d& p_cdm) const;
	/// Constrain point to the domain of this control space
	virtual DPoint2d fitInPoint(const DPoint2d& pt) const { return pt; }
	/// Whether reparameterization is used in some subarea (point) of control space
	virtual ReparameterizationData* getReparameterization(const DPoint2d& ) const { return nullptr; }
	/// Invoke function for all control nodes of this space (read-only)
	virtual void forEachControlNode(const std::function<void(const ControlNode2d& node)>& /* fg */) const {}
	/// Invoke function for all control nodes of this space
	virtual void forEachControlNode(const std::function<void(ControlNode2d& node)>& /* fg */) {}
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const { return 0; }
public:
	/// Returns the minimum size of element for the given point 2D within the domain
	static double getMinSpaceValue() { return 0.4; }
public:
	virtual void showMetricLength(const string& label = "") const {}
public:
	/// Sets underlying parameteric surface
	void setBaseSurface(const SurfaceConstPtr& s) { base_surface = s; }
	/// Returns underlying parameteric surface
	const SurfaceConstPtr& getBaseSurface() const { return base_surface; }
public:
	/// Whether the approximated (cached) parameterization matrix should be used
	static int param_cached_parameterization_matrix;
	/// Which type of control space should be used
	static int param_control_type;
	/// Which type of control space should be used
	static double param_inflate_box_factor;
protected:
	/// Pointer to the definition of base surface for the domain of this control space
	SurfaceConstPtr base_surface;
};

typedef std::shared_ptr<ControlSpace2d> CS2dPtr;
typedef std::shared_ptr<const ControlSpace2d> CS2dConstPtr;

#endif // !defined(CONTROLSPACE2D_H__INCLUDED)
