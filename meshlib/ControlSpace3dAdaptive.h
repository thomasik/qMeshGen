/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dAdaptive.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DADAPTIVE_H__INCLUDED)
#define CONTROLSPACE3DADAPTIVE_H__INCLUDED

#include "ControlSpace3d.h"
#include "DRect.h"
#include "DataVector.h"
#include "ControlSpace2dAdaptive.h"

class MeshPoint3d;
class MeshEdge3d;
//class MeshPoint2d;
class MeshContainer3d;
class MeshContainer3dSurface;
class MeshViewSet;
class Metric3dContext;
class MeshFace;
class SurfaceParametric;
class Curve3dParametric;

#define COMPACTED_READ_ONLY 10

/**
 * This class implements an abstract adaptive control space
 *  which provides the information about the sizing and/or 
 *  stretching of elements within the meshed domain.
 */
class ControlSpace3dAdaptive : public ControlSpace3d
{
public:
	/// Standard constructor
	ControlSpace3dAdaptive(const DBox& box)
		: m_box(box), m_initialized(0) { }
	/// Destructor
	virtual ~ControlSpace3dAdaptive() { }
public:
	/// Invoke (test) function in regular grid over this space
	virtual void forRegularGrid(int grid_size, const std::function<void(const DPoint3d& pt)>& fg);
	/// Insert additional info -> about local surface available for the given point
	virtual bool addLocalSurfaceAtPoint(const DPoint3d& /* pt */, SurfaceConstPtr /* local_surface */, 
		DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & /* local_surface_hash */) { assert(false); return false; }
	/// Insert additional info -> about local surface available for the given point
	virtual bool addLocalSurfaceAtBBox(const DBox& /* box */, SurfaceConstPtr /* local_surface */, 
		DataHashTableKeyValue<SurfaceSetPtr, SurfaceSetPtr> & /* local_surface_hash */) { assert(false); return false; }
	/// Returns additional info -> about local surfaces (set) available for the given point
	virtual SurfaceSetConstPtr getLocalSurfaceSetAtPoint(const DPoint3d& /* pt */) { assert(false); return nullptr; }
	/// Refines control space basing on other control space, within given rectangular area
	bool setMinControlRect(Metric3dContext& mc, CS3dPtr space, 
		const DPoint3d& pt00, const DVector3d& dv0, const DVector3d& dv1, 
		const DVector3d& dvn, double r);
	/// Refines control space basing on other control space, along given segment
	bool setMinControlLine(Metric3dContext& mc, CS3dPtr space, 
		const DPoint3d& pt0, const DPoint3d& pt1, 
		const DVector3d& dvn, double r);
	/// Refines control space basing on other control space, along given segment
	bool setMinControlPointVicinity(CS3dPtr space, const DPoint3d& pt,  
		const DVector3d& dvn, double r, double dt = 0.5);
	/// Refines control space at the given point translated by the vector (properly fitted to bounding box)
	bool setMinControlTranslated(const DPoint3d& pt, const DVector3d& dv, 
		const ControlDataMatrix3d& cdm, bool min_value_set = true);
	/// Adds new information (stretching and lengths) along the segment within the domain
//	void addControlSegment(const DPoint3d& pt1, const DPoint3d& pt2, 
//		const ControlDataStretch3d& data1, double r1, 
//		const ControlDataStretch3d& data2, double r2);
	/// Smoothen metric for two nodes
	static int smoothenMetricForNodes(ControlNode3d * cn1, ControlNode3d * cn2, const DVector3d& dv);
	/// Smoothen metric for two nodes
	static bool smoothenMetricForNodesLeft(ControlNode3d * cn1, const ControlNode3d * cn2, const DVector3d& dv);
	/// Smoothen metric for two nodes
	static bool adjustMetricForNodesLeft(ControlDataMatrix3d& cdm0, const ControlDataMatrix3d& cdm1, 
		const DVector3d& dv);
	/// Smoothen metric for two nodes
	static bool adjustMetricForNodesLeft(ControlDataMatrix3d& cdm0, const ControlDataMatrix3d& cdm1, 
		const DVector3d& dv, double& h0_max, double h1_min);
	/// Initializes the control space with maximum metric
	void setMaxMetric(double ratio = -1.0);
	double getMaxMetricLen(double ratio = -1.0) const;
	/// Initializes the control space with given global metric
	virtual void setGlobalMetric(const ControlDataMatrix3d& cdm);
	/// Returns initialization state (0 - no, 1 - control data, 2 - and parameterization data)
	int initializationState() const { return m_initialized; }
	/// Sets initialization state (0 - no, 1 - control data, 2 - and parameterization data)
	void setInitializationState(int state = 1) { m_initialized = state; }
	/// Returns the boundary of this control space
	const DBox& getBoundingBox() const { return m_box; }
	/// Returns the bounding box diameter, containing all points and faces from this container
	double getBoundingBoxDiameter( ) const { 
		return m_box.getDiameter(); }
	/// Refines control space at the given point with the given extended metric (returns number of completely new-set nodes)
	int setMinControlValue(const ControlDataExtMatrix3d& data);
	/// Calculates minimum of two control spaces
	bool applyAsMinimum(CS3dPtr space);
	/// Calculates minimum of two control spaces
	bool applyAsMinimum(CS2dPtr space);
	/// Compares two control spaces (returns max difference)
	double compareMetricWith(CS3dPtr space);
	/// Refines control space at the given point with the given extended metric
	bool setMinControl(const ControlDataExtMatrix3dSphere& data);
	/// Refines control space along segment (two points) with the given extended metric
	bool setMinControl(const ControlDataExtMatrix3dSegment& data);
	/// Refines control space for triangle (three points) with the given extended metric
	bool setMinControl(const ControlDataExtMatrix3dTriangle& data);
	/// Adjusts control space for segment definition
	bool setSegmentMinControl(const DPoint3d& pt, const DVector3d& dv, double dt,
		const ControlDataMatrix3d& cdm, double r);
	/// Adjusts control space for triangle definition
	bool setTriangleMinControl(const DPoint3d& pt, const DVector3d& dv0, const DVector3d& dv1, 
		double dt0, double dt1, const ControlDataMatrix3d& cdm, double r);
	/// Adjusts control space for spherical definition
	bool setSphereMinControl(const DPoint3d& pt, const ControlDataMatrix3d& cdm, double r);
	/// Adjusts control space for radial definition
//	void addSphericalControlPoint(const DMetric2d& dm, const DPoint2d& pt, double r, 
//		const ControlDataMatrix2d& cdm, double a0, double a1, double min_len);
	/// Adds new information - extended! (stretching and lengths) for some point within the domain
//	void addControlPoint(const DPoint3d& pt, const ControlDataStretch3d& data1, double r1, 
//		const ControlDataStretch3d& data2, double r2);
	/// Adds new information - extended! (stretching and lengths) for some point within the domain
//	void addSegmentControlPoint(const DPoint3d& pt, const DPoint3d& dv, double dt,
//		const ControlDataStretch3d& data, double r);
	/// Update control space with boundary triangles
//	bool updateForBoundaryTriangle(const MeshEdge2d* edge, const MeshPoint2d* point);
	/// Update control space with length of discretized boundary segment
	bool updateForBoundarySegment(const MeshEdge3d* edge);
	/// Update control space with bad represented shape for edge
//	bool updateForBoundaryShape(const MeshEdge2d* edge);
	/// Initializes the control space basing on the boundary (using only simple faces)
	virtual void addSimpleBoundaryControlDataIsotroPIc(DataVector<MeshFace*> & bfaces, DataVector<MeshPoint3d*> & bpoints);
	/// Initializes the control space basing on the boundary (using only simple faces)
	virtual void addSimpleBoundaryControlData(DataVector<MeshFace*> & bfaces, DataVector<MeshPoint3d*> & bpoints);
	/// Initializes the control space basing on the boundary (using only simple faces)
	virtual bool addSimpleBoundaryControlData(MeshContainer3dSurface* surface_mesh);
	/// Initializes the control space basing on the mesh (using only blocks)
	virtual bool addSimpleMeshControlData(MeshContainer3d* mesh);
	/// Refines control space at the given point using curvature of the given surface
	virtual bool setMinControlFromSurfaceCurvature(const DPoint3d& pt, SurfaceConstPtr surface, double model_diameter, 
		bool min_value_set = true, ControlDataMatrix3d* surf_cdm = nullptr, SurfaceCurvature *curvature2d = nullptr);
	/// Refines control space at the given face using curvature of the given surface
	virtual bool setMinControlFromSurfaceCurvature(MeshFace* face, SurfaceConstPtr surface, double model_diameter, 
		bool min_value_set = true, ControlDataMatrix3d* surf_cdm = nullptr, SurfaceCurvature *curvature2d = nullptr);
	/// Calculates control space for the given point using curvature of the given surface
	static bool calculateControlFromSurfaceCurvature(const DPoint3d& pt, SurfaceConstPtr surface, 
		double model_diameter, ControlDataMatrix3d& surf_cdm, const DPoint2d * pt2d = nullptr, SurfaceCurvature *curvature2d = nullptr);
	/// Refines control space at the given point using curvature of the given curve
	virtual bool setMinControlFromCurveCurvature(const DPoint3d& pt, 
		std::shared_ptr<const Curve3dParametric> curve, double model_diameter, bool min_value_set = true);
public:
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint3d& pt) const = 0;
	/// If possible to adapt
	virtual bool isAdaptive() const { return true; }
	/// Adds new information (stretching and lengths) for some point within the domain
	virtual void addControlPoint(const DPoint3d &pt, const ControlDataMatrix3d& data){
		addControlNode(ControlNode3d(pt, data)); 
	}
	/// Constrain point to the domain of this control space
	virtual DPoint3d fitInPoint(const DPoint3d& pt) const { return m_box.fitInPoint(pt); }
public:  // ------- to implement in derived classes --------------
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const = 0;
	/// Adds new information (stretching and lengths) for some point within the domain
	virtual void addControlNode(const ControlNode3d& node) = 0;
	/// Interpolates (initializes) the control space basing on the previously given sources
	virtual bool interpolate() = 0;
	/// Refines control space at the given point
	virtual bool setMinControl(const DPoint3d& pt, const ControlDataMatrix3d& cdm, bool min_value_set = true) = 0;
	/// Returns the "screenshot" of this control space for visualization
	virtual MeshViewSet* getViewSet(MeshViewSet* set = nullptr, bool /* with_points */ = false) const { return set; }
	/// Log basic information about this control space
	virtual void logDescription() const = 0;
	/// Smoothen variance of metric within the control space (returns true if any change)
	virtual bool smoothen() { return false; }
	/// optimize for size (lower memory requirements) - but also from now readonly !
	virtual void compact() { m_initialized = COMPACTED_READ_ONLY; }
	virtual bool isCompacted() const { return m_initialized == COMPACTED_READ_ONLY; }
	/// Calculate inside information for nodes/elements
	virtual void markInsideNodes(const MeshContainer3dSurface* /* surface_mesh */) { }
	/// Returns local resolution
	virtual double getLocalResolution(const DPoint3d& /* pt */) const { 
		return std::min(std::min(m_box.getDX(), m_box.getDY()), m_box.getDZ()); }
	/// Return interpolated value of extended tag data from control nodes at some given point;
	virtual double interpolateDoubleTag(const DPoint3d& pt, TagExtended::TagType type) const = 0;
	/// Set given tag value as minimum at the closest control element
	virtual void setMinDoubleTag(const DPoint3d& /* pt */, TagExtended::TagType /* type */, double /* q */) { assert(false); }
	/// Set given tag value as minimum at the closest control element
	virtual void setMaxDoubleTag(const DPoint3d& /* pt */, TagExtended::TagType /* type */, double /* q */) { assert(false); }
	/// Whether the given point actually is withing the domain of this space
	virtual bool containsPoint(const DPoint3d& pt) const { return m_box.contains(pt); }
public:
	/// Smoothen variance of metric within the control space (with specified gradation ratio)
	bool smoothenWithGradation(double gradation);
	/// Check metric lengths along given segment
	bool logMetricAlongSegment(const string& fname, const DPoint3d& pt0, const DPoint3d& pt1, 
		const DVector3d& vn, int n = 100);
	/// Check metric lengths along some segments
	bool testMetric(const string& test_name = "generic");
	/// Metric gradation tests
	double calculateMaxGradationRatioViaNodes() const;
	double calculateMaxGradationRatioRandom(int probe_count = 1000) const;
	double calculateMaxGradationRatioRegular(int sample_count_per_dim = 11) const;
public:
	/// Maximum element size as a model diameter ratio for automated sizing
	static double param_max_diameter_ratio;
	/// Number of probe points for radial metric definition
//	static int param_radial_parts;
protected:
	/// Rectangle specifying the area of control space
	DBox m_box;
	/// Whether this control space is already initialized (has value in all nodes)
	int m_initialized;  // 0 - no, 1 - control data, 2 - and parameterization data, COMPACTED_READ_ONLY - "compacted" read only mode
	/// Extended metrics for enforcing minimum size after interpolation
};

typedef std::shared_ptr<ControlSpace3dAdaptive> CS3dAPtr;
typedef std::shared_ptr<const ControlSpace3dAdaptive> CS3dAConstPtr;

#endif // !defined(CONTROLSPACE3DADAPTIVE_H__INCLUDED)
