/////////////////////////////////////////////////////////////////////////////
// Metric3dContext.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(METRIC3DCONTEXT_H__INCLUDED)
#define METRIC3DCONTEXT_H__INCLUDED

#include "DMetric3d.h"
#include "DPoint.h"
class ControlSpace3d;
class MeshPoint3d;


class Metric3dContext
{
public:
	/// Standard constructor
	Metric3dContext(std::shared_ptr<const ControlSpace3d> space) :
		m_control_space(space.get()), m_radius(-1.0), m_min_skip_len2(-1.0) { }
	Metric3dContext(const ControlSpace3d* space = nullptr) :
		m_control_space(space), m_radius(-1.0), m_min_skip_len2(-1.0) { }
public:
	/// Return the gradation ratio for given point
	double getMetricGradation(const DPoint3d& pt) const;
	/// Returns coordinates of point transformed from real to metric space
	const DMPoint3d transformRStoMS(const DPoint3d& pt) const {
		return m_metric.transformRStoMS(pt);
	}
	/// Returns coordinates of point transformed from metric to real space
	const DPoint3d transformMStoRS(const DMPoint3d& pt) const {
		return m_metric.transformMStoRS(pt);
	}
	/// Returns coordinates of vector transformed from real to metric space
	const DMVector3d transformRStoMS(const DVector3d& v) const {
		return m_metric.transformRStoMS(v);
	}
	/// Returns coordinates of vector transformed from metric to real space
	const DVector3d transformMStoRS(const DMVector3d& v) const {
		return m_metric.transformMStoRS(v);
	}
	/// Set new metric as special one
	void setSpecialMetric(const ControlDataMatrix3d& cdm);
	/// Count new metric at the given point
	void countMetricAtPoint(const DPoint3d& pt, bool if_needed_only = false);
	/// Count new metric for the given tetrahedron (4 points)
	void countMetricAtPoints(MeshPoint3d** points);
	/// Count new metric for the given edge (2 points)
	void countMetricAtPoints(MeshPoint3d* point0, MeshPoint3d* point1);
	/// Return the control space
	const ControlSpace3d* getControlSpace() const { return m_control_space; }
	/// Returns maximumem length of edge in metric
	double maxLength() const { return m_metric.maxLength(); }
	/// Returns maximumem length of edge in metric
	double minLength() const { return m_metric.minLength(); }
	/// Sets radius
	void setRadius(double r) { m_radius = r; }
	/// Set no radius
	void setNoRadius() { m_radius = -1.0; }
	/// Get radius
	double getRadius() const { return m_radius; }
	/// Check no radius
	bool noRadius() const { return m_radius < 0.0; }
	/// Check point for impact sphere
	bool withinImpactRadius(const DPoint3d& pt) const;
	/// Set minimum skip length
	void setMinSkipLen2(double skip_len2) { m_min_skip_len2 = skip_len2; }
	/** 
	 * Get minimum skip length
	 */
	double getMinSkipLen2() const { return m_min_skip_len2; }
protected:
	/// Control space
	const ControlSpace3d* m_control_space;
	/// Current metric
	DMetric3d m_metric;
	/// Impact sphere radius
	double m_radius;
	/// Center for metric impact sphere (where it was calculated)
	DMPoint3d m_center;
	/// Minimum skip length (squared) for skiPIng validity check for tiny elements during validation by collapse
	double m_min_skip_len2;
};

#endif // !defined(METRIC3DCONTEXT_H__INCLUDED)
