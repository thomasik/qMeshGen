/////////////////////////////////////////////////////////////////////////////
// ControlSpace2dMesh.h: 
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACEMESH_H__INCLUDED)
#define CONTROLSPACEMESH_H__INCLUDED

#include "ControlSpace2dAdaptive.h"
#include "DataMatrix.h"
#include "Metric2dContext.h"

class MeshContainer2d;
class MeshTriangle2d;
class MeshViewSet;

/**
 * This class implements a control space given by a traingular mesh ("background mesh").
 */
class ControlSpace2dMesh : public ControlSpace2dAdaptive
{
public:
	/// Standard contructor
	ControlSpace2dMesh(SurfaceConstPtr surface, const DRect& box, int approx);
	/// Destructor
	virtual ~ControlSpace2dMesh();
public:
	/// Returns the type of control space
	virtual int getType() const { return MeshData::CONTROL_MESH; }
	/// Interpolates (initializes) the control space basing on the previously given sizing information in discrete points
	virtual bool interpolate();
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const override;
	/// Invoke function for all control nodes of this space (read-only)
	virtual void forEachControlNode(const std::function<void(const ControlNode2d& node)>& fg) const override;
	/// Invoke function for all control nodes of this space
	virtual void forEachControlNode(const std::function<void(ControlNode2d& node)>& fg) override;
	/// Adds new information (stretching and lengths) for some point within the domain
	virtual void addControlNode(const ControlNode2d& node);
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint2d& pt) const;
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt) const { 
		return getMetricAtPoint(pt, -1); 
	}
	/// Get sizing info (matrix mode) at the given point
	ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt, int method, ControlDataMatrix2d * pcdm = nullptr) const;
	/// Get transformation and parameterization matrices at the given point
	virtual ControlDataMatrix2d getMetricAndParameterizationAtPoint(const DPoint2d& pt, ControlDataMatrix2d& p_cdm) const;
	/// Initializes the control space basing on the curvature of surface and boundary
	virtual void setSurfaceCurvatureControlData();
	/// Log basic information about this control space
	virtual void logDescription() const;
	/// Smoothen variance of metric within the control space
	virtual bool smoothen();
	/// Set maximum metric
	virtual void setMaxMetric(double ratio = -1.0);
	/// Return interpolated value of extended tag data from control nodes at some given point;
	virtual double interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const;
public: // TODO
	/// Refines control space at the given point
	virtual bool setMinControl(const DPoint2d& pt, const ControlDataMatrix2d& cdm, bool min_value_set = true);
	/// Refines control space to parameterization variance
	virtual void adaptToParameterization();
public:
	/// Returns the "screenshot" of this control space for visualization
	virtual MeshViewSet* getViewSet(MeshViewSet* set = nullptr, bool with_surface = false, 
		TagExtended::TagType tag = TagExtended::TAG_NONE) const;
	/// Draws the structure of the control space
	virtual void storeEPS(const char* name = "control", int id = 0);
public:
	/// Check mesh elements type (with respect to the interpolation method
	void checkControlMeshForInterpolationType();
	/// Alters the default approximation method for this control space
	void setApproximationMethod(int method) { m_approximation_method = method; }
	/// Returns number of triangles in the control mesh
	int getTrianglesCount() const;
protected:
	/// Insert new node with mesh refinement (possibly new nodes) in the vicinity
	MeshPoint2d* insertNodeAndRefine(const DPoint2d& point, MeshTriangle2d* containing_triangle);
	/// Count regular grid points
	bool countRegularGridPoints(int count, DataVector<DPoint2d> & points) const;
	/// Add auxiliary nodes into the control mesh
	int addInnerNodes(Metric2dContext& mc, int max_ct = -1);
	/// Check mesh elements type (with respect to the interpolation method
	int getRequiredInterpolationType(const DPoint2d& pt, MeshTriangle2d* triangle, double threshold);
	/// Returns metric from the triangular mesh (simple method - first vertex)
	ControlDataMatrix2d getMetricSimple(const DPoint2d& pt, MeshTriangle2d* triangle) const;
	/// Returns metric from the triangular mesh (interpolation from points - triangle vertices or voronoi)
	ControlDataMatrix2d getMetricInterpolate(const DPoint2d& pt, MeshTriangle2d* triangle, bool use_voronoi_points, bool square_weight) const;
	/// Returns metric from the triangular mesh (interpolation from points - NEM layered)
	//ControlDataMatrix2d getMetricInterpolateNEM(const DPoint2d& pt, MeshTriangle2d* triangle);
	bool isValid() const;
public:
	/// Interpolation method for control mesh
	static int param_interpolation_method;
	/// Whether the calculated weights for the metric value inside the triangle should be squared
	static int param_weights_squared;
	/// Whether inner nodes should be inserted during preparation of this control mesh
	static int param_inner_nodes;
	/// Threshold for metric difference for assesing mesh-control interpolation type
	static double param_mixed_threshold;
private:
	/// Triangular mesh containing the control space data
	MeshContainer2d* m_control_mesh;
	/// Mesh points with control data before triangulation
	MeshContainer2d* m_control_points;
	/// Method of approximating the metric value inside the triangle
	int m_approximation_method;
	/// Nodes
	DataVector<ControlNode2d> m_nodes;
	/// minimum length
	double m_min_length;
};

#endif // !defined(CONTROLSPACEMESH_H__INCLUDED)
