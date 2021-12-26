/////////////////////////////////////////////////////////////////////////////
// ControlSpace2dMatrixUniform.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACEMATRIXUNIFORM_H__INCLUDED)
#define CONTROLSPACEMATRIXUNIFORM_H__INCLUDED

#include "ControlSpace2dAdaptive.h"
#include "DRect.h"

class MeshContainer2d;
class MeshEdge2d;

/**
 * This class implements a regular, uniform control space.
 */
class ControlSpace2dMatrixUniform : public ControlSpace2dAdaptive
{
public:
	/// Standard contructor
	ControlSpace2dMatrixUniform(SurfaceConstPtr surface, 
		const DRect& box, int nx = 10, int ny = 10);
	/// Destructor
	virtual ~ControlSpace2dMatrixUniform();
public:
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const override;
	/// Invoke function for all control nodes of this space (read-only)
	virtual void forEachControlNode(const std::function<void(const ControlNode2d& node)>& fg) const override;
	/// Invoke function for all control nodes of this space
	virtual void forEachControlNode(const std::function<void(ControlNode2d& node)>& fg) override;
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix2d getMetricAtPoint(const DPoint2d& pt) const;
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint2d& pt) const;
	/// Interpolates (initializes) the control space basing on the previously given sizing information in discrete points
	virtual bool interpolate();
	/// Adds new information (stretching and lengths) for some point within the domain
	virtual void addControlNode(const ControlNode2d& node);
	/// Initializes the control space basing on the curvature of surface and boundary
	virtual void setSurfaceCurvatureControlData();
	/// Draws the structure of the control space
	virtual void storeEPS(const char* name = "control", int id = 0);
	// Refines control space at the given point with the given extended metric
	virtual bool setMinControl(const DPoint2d& pt, const ControlDataMatrix2d& cdm, bool min_value_set = true);
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return MeshData::CONTROL_UNIFORM; }
	/// Refines control space to parameterization variance
	virtual void adaptToParameterization();
	/// Log basic information about this control space
	virtual void logDescription() const;
	/// Return interpolated value of extended tag data from control nodes at some given point;
	virtual double interpolateDoubleTag(const DPoint2d& pt, TagExtended::TagType type) const;
protected:
	/// Adds information about size-requirements (as matrix-form)
	void addControlData(const DMetric2d& dm, int k, const ControlNode2d& node);
	/// Sets information about size-requirements (as matrix-form)
	bool setMinControlData(int k, const ControlDataMatrix2d& data);
	/// Translates the point-coordinates into index of control matrix
	void countLocalCoordinates(const DPoint2d& pt, int& ix, int& iy) const;
	/// Translates the point-coordinates into index of control matrix
	int countLocalCoordinates(const DPoint2d& pt) const;
public:
	/// resolution for uniform control grid
	static int param_uniform_nx;
protected:
	/// Number of columns in the matrix
	int m_nx;
	/// Number of rows in the matrix
	int m_ny;
	/// Real-width of the matrix-column
	double m_dx;
	/// Real-height of the matrix-column
	double m_dy;
	/// Main grid
	ControlNode2d* m_space;
};

#endif // !defined(CONTROLSPACEMATRIXUNIFORM_H__INCLUDED)
