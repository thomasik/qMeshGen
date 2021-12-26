/////////////////////////////////////////////////////////////////////////////
// ControlSpace2dAdaptive.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACEADAPTIVE_H__INCLUDED)
#define CONTROLSPACEADAPTIVE_H__INCLUDED

#include <memory>

#include "ControlSpace2d.h"
#include "ControlSpace3d.h"
#include "DRect.h"
#include "TagExtended.h"
#include "DataVector.h"

class MeshContainer2d;
class MeshEdge2d;
class MeshEdge2dCurve;
class MeshPoint2d;
class Metric2dContext;

/**
 * This class implements an abstract adaptive control space
 *  which provides the information about the sizing and/or 
 *  stretching of elements within the meshed domain.
 */
class ControlSpace2dAdaptive : public ControlSpace2d
{
public:
	/// Standard constructor
	ControlSpace2dAdaptive(SurfaceConstPtr surface, const DRect& box)
		: ControlSpace2d(surface), m_box(box), m_initialized(0) { }
	/// Destructor
	virtual ~ControlSpace2dAdaptive() { }
public:
	/// Refines control space basing on other control space, along given segment
	bool setMinControlLine(Metric2dContext& mc, CS2dPtr space, 
		const DPoint2d& pt0, const DPoint2d& pt1, 
		const DVector2d& dvn, double r);
	/// Refines control space basing on other control space, along given segment
	bool setMinControlPointVicinity(CS2dPtr space, const DPoint2d& pt,  
		const DVector2d& dvn, double r, double dt = 0.5);
	/// Adds new information (stretching and lengths) for some point translated by the vector (properly fitted to bounding box)
	void addControlNodeTranslated(const ControlNode2d& node, const DVector2d& dv);
	/// Refines control space at the given point translated by the vector (properly fitted to bounding box)
	bool setMinControlTranslated(const DPoint2d& pt, const DVector2d& dv, 
		const ControlDataMatrix2d& cdm, bool min_value_set = true);
	/// Adds new information (stretching and lengths) along the curved segment within the domain
	void addControlSegment(std::shared_ptr<const Curve2dParametric> curve, double t0, double t1,
		const ControlDataMatrix2d& data, double r);
	/// Adds new information (stretching and lengths) along the segment within the domain
	void addControlSegment(const DPoint2d& pt1, const DPoint2d& pt2, 
		const ControlDataMatrix2d& data, double r);
	/// Initializes the control space basing on the curvature of boundary
	bool setContourCurvatureControlData(const MeshContainer2d* boundary);
	/// Returns initialization state (0 - no, 1 - control data, 2 - and parameterization data)
	int initializationState() const { return m_initialized; }
	/// Returns the boundary of this control space
	const DRect& getBoundingRect() const { return m_box; }
	/// Refines control space at the given point with the given extended metric (returns number of completely new-set nodes)
	int setMinControlValue(const ControlDataExtMatrix2d& data);
	/// Calculates minimum of two control spaces
	bool applyAsMinimum(CS2dPtr space);
	/// Calculates minimum of two control spaces
	bool applyAsMinimum(CS3dPtr space);
	/// Compares two control spaces (returns max difference)
	double compareMetricWith(CS2dPtr space);
	/// Refines control space at the given point with the given extended metric
	bool setMinControl(const ControlDataExtMatrix2dRadial& data);
	/// Refines control space at the given point with the given extended metric
	bool setMinControl(const ControlDataExtMatrix2dSegment& data);
	/// Refines control space at the given point with the given extended metric
	bool setMinControl(const ControlDataExtMatrix2dCurve& data);
	/// Adjusts control space for segment definition
	bool setSegmentMinControl(const DPoint2d& pt, const DVector2d& dv, double dt,
		const ControlDataMatrix2d& cdm, double r);
	/// Adjusts control space for radial definition
	bool setRadialControlPoint(const DMetric2d& dm, const DPoint2d& pt, double r, 
		const ControlDataMatrix2d& cdm, double a0, double a1, double min_len);
	/// Adjusts control space for radial definition
	void addRadialControlPoint(const DMetric2d& dm, const DPoint2d& pt, double r, 
		const ControlDataMatrix2d& cdm, double a0, double a1, double min_len);
	/// Adds new information - extended! (stretching and lengths) for some point within the domain
	void addControlPoint(const DPoint2d& pt, const ControlDataMatrix2d& data, double r);
	/// Adds new information - extended! (stretching and lengths) for some point within the domain
	void addSegmentControlPoint(const DPoint2d& pt, const DVector2d& dv, double dt,
		const ControlDataMatrix2d& data, double r);
	/// Update control space with boundary triangles
	bool updateForBoundaryTriangle(const MeshEdge2d* edge, const MeshPoint2d* point);
	/// Update control space with length of discretized boundary segment
	bool updateForBoundarySegment(const MeshEdge2d* edge);
	/// Update control space with bad represented shape for edge
	bool updateForBoundaryShape(const MeshEdge2d* edge);
	/// Smoothen variance of metric within the control space (with specified gradation ratio)
	bool smoothenWithGradation(double gradation);
public:
	/// Initializes the control space with maximum metric
	virtual void setMaxMetric(double ratio = -1.0);
	/// Constrain point to the domain of this control space
	virtual DPoint2d fitInPoint(const DPoint2d& pt) const { return m_box.fitInPoint(pt); }
	/// If possible to adapt
	virtual bool isAdaptive() const { return true; }
	/// Adds new information (stretching and lengths) for some point within the domain
	virtual void addControlPoint(const DPoint2d &pt, const ControlDataMatrix2d& data){
		addControlNode(ControlNode2d(pt, data)); 
	}
	/// Invoke function for all control nodes of this space (read-only)
	virtual void forEachControlNode(const std::function<void(const ControlNode2d& node)>& ) const { assert(false); }
	/// Invoke function for all control nodes of this space
	virtual void forEachControlNode(const std::function<void(ControlNode2d& node)>& ) { assert(false); }
public:  // ------- to implement in derived classes --------------
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const = 0;
	/// Adds new information (stretching and lengths) for some point within the domain
	virtual void addControlNode(const ControlNode2d& node) = 0;
	/// Interpolates (initializes) the control space basing on the previously given sources
	virtual bool interpolate() = 0;
	/// Draws the structure of the control space
	virtual void storeEPS(const char* /* name */ = "control", int /* id */ = 0) {}
	/// Draws the metric-size of the nodes within the control space
	virtual void storeSizingEPS(const MeshContainer2d* mesh = nullptr, const char* fname = "control", int id = 0, double vlen = 0.1);
	/// Initializes the control space basing on the curvature of surface
	virtual void setSurfaceCurvatureControlData() = 0;
	/// Refines control space at the given point
	virtual bool setMinControl(const DPoint2d& pt, const ControlDataMatrix2d& cdm, bool min_value_set = true) = 0;
	/// Returns the "screenshot" of this control space for visualization
	virtual MeshViewSet* getViewSet(MeshViewSet* set = nullptr, bool /* with_surface */ = false, 
		TagExtended::TagType /* tag */ = TagExtended::TAG_NONE) const { return set; }
	/// Refines control space to parameterization variance
	virtual void adaptToParameterization() = 0;
	/// Refines control space with respect to reference 3d control space
	virtual bool adaptToControlSpace3d(CS3dConstPtr /* space */) { return false; }
	/// Log basic information about this control space
	virtual void logDescription() const = 0;
	/// Smoothen variance of metric within the control space
	virtual bool smoothen() { return false; }
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint2d& pt) const = 0;
	/// Return interpolated value of extended tag data from control nodes at some given point;
	virtual double interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const = 0;
	/// Set given tag value as minimum at the closest control element
	virtual void setMinDoubleTag(const DPoint2d& /* pt */, TagExtended::TagType /* type */, double /* q */) { assert(false); }
	/// Set given tag value as minimum at the closest control element
	virtual void setMaxDoubleTag(const DPoint2d& /* pt */, TagExtended::TagType /* type */, double /* q */) { assert(false); }
	/// Calculate inside information for nodes/elements
	virtual void markInsideNodes(const MeshContainer2d* /* boundary_mesh */) { }
	/// Returns local resolution
	virtual double getLocalResolution(const DPoint2d& /* pt */) const;
public:
	virtual void showMetricLength(const string& label = "") const override;
public:
	/// Smoothen metric for two nodes
	static int smoothenMetricForNodes(ControlNode2d * cn1, ControlNode2d * cn2, const DVector2d& dv);
	/// Adjusts control data from curvature information
	static ControlDataStretch2d adjustCurvatureData(const SurfaceCurvature& c, double c_ratio, double model_diameter,
		double min_len, double max_len, double max_stretch_ratio);
public:
	/// default starting resolution for quadtree/uniform control
	static int param_control_nxy;
	/// Control metric gradation ratio
	static double param_gradation_ratio;
	/// Curvature ratio for automated sizing
	static double param_curvature_ratio;
	/// Maximum stretching (element length ratio) for automated sizing
	static double param_stretch_max_ratio;
	/// Maximum stretching (element length ratio) for automated sizing (contour curvature)
	static double param_contour_stretch_max_ratio;
	/// Maximum element size as a model diameter ratio for automated sizing
	static double param_max_diameter_ratio;
	/// Maximum element size as a model diameter ratio for automated sizing if the face is an inner boundary
	static double param_inner_boundary_ratio;
	/// Minimum element size as a model diameter ratio for automated sizing
	static double param_min_diameter_ratio;
	/// Minimum element size (length) for automated sizing
	static double param_min_length;
	/// Whether use surface curvature
	static int param_use_surface_curvature;
	/// Whether use contour curvature
	static int param_use_contour_curvature;
	/// Whether metric is too different
	static double param_threshold_diff;
	/// Number of probe points for radial metric definition
	static int param_radial_parts;
protected:
	/// Rectangle specifying the area of control space
	DRect m_box;
	/// Whether this control space is already initialized (has value in all nodes)
	int m_initialized;  // 0 - no, 1 - control data, 2 - and parameterization data
	/// Extended metrics for enforcing minimum size after interpolation
	DataVector<std::shared_ptr<ControlDataExtMatrix2d>> m_ext_source_points;
};

typedef std::shared_ptr<ControlSpace2dAdaptive> CS2dAPtr;
typedef std::shared_ptr<const ControlSpace2dAdaptive> CS2dAConstPtr;

#endif // !defined(CONTROLSPACEADAPTIVE_H__INCLUDED)
