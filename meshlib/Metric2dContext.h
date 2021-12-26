/////////////////////////////////////////////////////////////////////////////
// Metric2dContext.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(METRICCONTEXT_H__INCLUDED)
#define METRICCONTEXT_H__INCLUDED

#include <memory>

#include "DMetric2d.h"
#include "DPoint.h"
#include "DataMatrix.h"

class ControlSpace2d;
class MeshPoint2d;
class ReparameterizationData;

class Metric2dContext
{
public:
	/// Standard constructor
	Metric2dContext(const ControlSpace2d* space = nullptr) : 
			m_control_space(space), m_metric(ControlDataMatrix2d::IDENTITY), 
			m_reparameterization(nullptr) { }
	Metric2dContext(std::shared_ptr<const ControlSpace2d> space) :
		m_control_space(space.get()), m_metric(ControlDataMatrix2d::IDENTITY),
		m_reparameterization(nullptr) { }
public:
	/// Return the gradation ratio for given point
	double getMetricGradation(const DPoint2d& pt) const;
	/// Return metric context id
	unsigned int getCurrentMetricID() const { return m_id; }
	/// Return the control space
	const ControlSpace2d* getControlSpace() const { return m_control_space; }
	/// Returns coordinates of point transformed from parametric to metric space
	const DMPoint2d transformPStoMS(const DPoint2d& pt);
	/// Returns coordinates of point transformed from parametric to real space
	const DPoint2d transformPStoRS(const DPoint2d& pt);
	/// Returns coordinates of point transformed from real to parametric space
	const DPoint2d transformRStoPS(const DPoint2d& pt);
	/// Returns coordinates of point transformed from metric to parametric space
	const DPoint2d transformMStoPS(const DMPoint2d& pt);
	/// Returns coordinates of point transformed from parametric to metric space
	const DMVector2d transformPStoMS(const DVector2d& v);
	/// Returns coordinates of point transformed from parametric to real space
	const DVector2d transformPStoRS(const DVector2d& v);
	/// Returns coordinates of point transformed from real to parametric space
	const DVector2d transformRStoPS(const DVector2d& v);
	/// Returns coordinates of point transformed from metric to parametric space
	const DVector2d transformMStoPS(const DMVector2d& v);
	/// Count new metric at the given point
	void countMetricAtPoint(const DPoint2d& pt);
	/// Count new metric at the given points
	void countMetricAtPoints(const MeshPoint2d* mpt0, const MeshPoint2d* mpt1, const MeshPoint2d* mpt2);
	/// Count new metric at the given points
	void countMetricAtPoints(const MeshPoint2d* mpt0, const MeshPoint2d* mpt1);
	/// Returns metric tensor for current metric
	ControlDataMatrix2d getMetricTensor() const { return m_metric.getMetricTensor(); }
protected:
	/// Unique id, if more than one metric context is used...
	unsigned int m_id;
	/// Control space
	const ControlSpace2d* m_control_space;
	/// Current metric
	DMetric2d m_metric;
	/// Reparameterization data (if available and used)
	ReparameterizationData* m_reparameterization;
};

#endif // !defined(METRICCONTEXT_H__INCLUDED)
