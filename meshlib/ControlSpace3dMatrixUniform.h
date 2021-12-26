/////////////////////////////////////////////////////////////////////////////
// ControlSpace3dMatrixUniform.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(CONTROLSPACE3DMATRIXUNIFORM_H__INCLUDED)
#define CONTROLSPACE3DMATRIXUNIFORM_H__INCLUDED

#include "ControlSpace3dAdaptive.h"
#include "DRect.h"
#include "MeshData.h"

class SurfaceParametric;

/**
 * This class implements a regular, uniform control space.
 */
class ControlSpace3dMatrixUniform : public ControlSpace3dAdaptive
{
public:
	/// Standard contructor
	ControlSpace3dMatrixUniform(const DBox& box, int nx = 10, int ny = 10, int nz = 10);
	/// Destructor
	virtual ~ControlSpace3dMatrixUniform();
public:
	/// Return metric gradation ratio at point
	virtual double getMetricGradationRatio(const DPoint3d& pt) const;
	/// Returns number of control nodes in adaptive control structure
	virtual int getControlNodesCount() const override;
	/// Invoke function for all control nodes of this space (read-only)
	virtual void forEachControlNode(const std::function<void(const ControlNode3d& node)>& fg) const override;
	/// Invoke function for all control nodes of this space
	virtual void forEachControlNode(const std::function<void(ControlNode3d& node)>& fg) override;
	/// Get sizing info (matrix mode) at the given point
	virtual ControlDataMatrix3d getMetricAtPoint(const DPoint3d& pt) const;
	/// Return interpolated value of extended tag data from control nodes at some given point;
	virtual double interpolateDoubleTag(const DPoint3d& pt, TagExtended::TagType type) const;
	/// Interpolates (initializes) the control space basing on the previously given sizing information in discrete points
	virtual bool interpolate();
	/// Adds new information (stretching and lengths) for some point within the domain
	virtual void addControlNode(const ControlNode3d& node);
	/// Refines control space at the given point
	virtual bool setMinControl(const DPoint3d& pt, const ControlDataMatrix3d& cdm, bool min_value_set = true);
	/// Returns the type of element (should be reimplemented in derived classes)
	virtual int getType() const { return MeshData::CONTROL_UNIFORM_3D; }
	/// Log basic information about this control space
	virtual void logDescription() const;
	/// Smoothen metric throughout the control domain (returns true if any change)
	virtual bool smoothen();
protected:
	/// Adds information about size-requirements (as matrix-form)
	void addControlData(int k, const ControlNode3d& node);
	/// Sets information about size-requirements (as matrix-form)
	bool setMinControlData(int k, const ControlDataMatrix3d& data);
	/// Translates the point-coordinates into local coordinates
	void countLocalCoordinates(const DPoint3d& pt, int& ix, int& iy, int& iz) const;
	/// Translates the index of control matrix into local coordinates
	void countLocalCoordinates(int k, int& ix, int& iy, int& iz) const;
	/// Translates the point-coordinates into index of control matrix
	int countLocalCoordinates(const DPoint3d& pt) const;
public:
	/// Resolution for uniform control grid
	static int param_uniform_nx;
protected:
	/// Number of columns in the matrix
	int m_nx;
	/// Number of rows in the matrix
	int m_ny;
	/// Number of stores in the matrix
	int m_nz;
	/// Real-width of the matrix-column
	double m_dx;
	/// Real-height of the matrix-column
	double m_dy;
	/// Real-depth of the matrix-column
	double m_dz;
	/// Main grid
	ControlNode3d* m_space;
};

#endif // !defined(CONTROLSPACE3DMATRIXUNIFORM_H__INCLUDED)
