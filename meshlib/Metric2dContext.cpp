// Metric2dContext.cpp: implementation of the Metric2dContext class.
//
//////////////////////////////////////////////////////////////////////

#include "Metric2dContext.h"
#include "ControlSpace2d.h"
#include "ControlSpace2dAdaptive.h"
#include "SurfaceParametric.h"
#include "MeshPoint2d.h"
#include "MeshTriangle2d.h"
#include "ReparameterizationData.h"

void Metric2dContext::countMetricAtPoint(const DPoint2d& pt)
{
	assert(m_control_space);

	const DPoint2d pt_fit = m_control_space->fitInPoint(pt);
	m_reparameterization = m_control_space->getReparameterization(pt_fit);

	if(m_reparameterization){
		const ControlDataMatrix2d p_cdm = m_reparameterization->countParameterizationMatrix(
			m_control_space->getBaseSurface()->getPoint(pt_fit));
		const ControlDataMatrix2d s_cdm = m_control_space->getMetricAtPoint(pt_fit);
		m_metric.setData(s_cdm, p_cdm);
	}else if(ControlSpace2d::param_cached_parameterization_matrix){
		ControlDataMatrix2d p_cdm;
		const ControlDataMatrix2d s_cdm = m_control_space->getMetricAndParameterizationAtPoint(pt_fit, p_cdm);
		m_metric.setData(s_cdm, p_cdm);
	}else{
		double g_ratio;
		const ControlDataMatrix2d p_cdm = 
			m_control_space->getBaseSurface()->countParameterizationMatrix(pt_fit, g_ratio);
		const ControlDataMatrix2d s_cdm = m_control_space->getMetricAtPoint(pt_fit);
		m_metric.setData(s_cdm, p_cdm);
	}
}

void Metric2dContext::countMetricAtPoints(const MeshPoint2d* mpt0, const MeshPoint2d* mpt1, const MeshPoint2d* mpt2)
{
	assert(m_control_space);

	const DPoint2d middle = m_control_space->fitInPoint(
		DPoint2d::average(
			mpt0->getCoordinates(),
			mpt1->getCoordinates(), 
			mpt2->getCoordinates()));
	ControlDataMatrix2d s_cdm, p_cdm;

	m_reparameterization = m_control_space->getReparameterization(middle);

	switch(MeshTriangle2d::param_quality_metric_selection){
	case MeshData::QM_MIDDLE:
		{
			if(m_reparameterization){
				p_cdm = m_reparameterization->countParameterizationMatrix(
					m_control_space->getBaseSurface()->getPoint(middle));
				s_cdm = m_control_space->getMetricAtPoint(middle);
			}else if(ControlSpace2d::param_cached_parameterization_matrix){
				s_cdm = m_control_space->getMetricAndParameterizationAtPoint(middle, p_cdm);
			}else{
				double g_ratio;
				p_cdm = m_control_space->getBaseSurface()->countParameterizationMatrix(middle, g_ratio);
				s_cdm = m_control_space->getMetricAtPoint(middle);
			}
			break;
		}
	case MeshData::QM_MIDDLE_AND_VERTICES_MIN:
	case MeshData::QM_MIDDLE_AND_VERTICES_AVE:
		{
			const DPoint2d pt0_fit = m_control_space->fitInPoint(mpt0->getCoordinates());
			const DPoint2d pt1_fit = m_control_space->fitInPoint(mpt1->getCoordinates());
			const DPoint2d pt2_fit = m_control_space->fitInPoint(mpt2->getCoordinates());

			s_cdm = m_control_space->getMetricAtPoint(pt0_fit);
			if(MeshTriangle2d::param_quality_metric_selection == MeshData::QM_MIDDLE_AND_VERTICES_MIN){
				s_cdm.setMinimum(m_control_space->getMetricAtPoint(pt1_fit));
				s_cdm.setMinimum(m_control_space->getMetricAtPoint(pt2_fit));
			}else{
				s_cdm += m_control_space->getMetricAtPoint(pt1_fit);
				s_cdm += m_control_space->getMetricAtPoint(pt2_fit);
			}

			ControlDataMatrix2d sm_cdm;
			if(m_reparameterization){
				p_cdm = m_reparameterization->countParameterizationMatrix(
					m_control_space->getBaseSurface()->getPoint(middle));
				sm_cdm = m_control_space->getMetricAtPoint(middle);
			}else if(ControlSpace2d::param_cached_parameterization_matrix){
				sm_cdm = m_control_space->getMetricAndParameterizationAtPoint(middle, p_cdm);
			}else{
				double g_ratio;
				SurfaceConstPtr surface = m_control_space->getBaseSurface();
				p_cdm = surface->countParameterizationMatrix(middle, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt0_fit, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt1_fit, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt2_fit, g_ratio);
				sm_cdm = m_control_space->getMetricAtPoint(middle);
			}
			if(MeshTriangle2d::param_quality_metric_selection == MeshData::QM_MIDDLE_AND_VERTICES_MIN){
				s_cdm.setMinimum(sm_cdm);
			}else{
				s_cdm += sm_cdm;
				s_cdm *= 0.25;
			}
			break;
		}
	case MeshData::QM_VERTICES_AVE:
		{
			const DPoint2d pt0_fit = m_control_space->fitInPoint(mpt0->getCoordinates());
			const DPoint2d pt1_fit = m_control_space->fitInPoint(mpt1->getCoordinates());
			const DPoint2d pt2_fit = m_control_space->fitInPoint(mpt2->getCoordinates());

			if(m_reparameterization){
				s_cdm = m_control_space->getMetricAtPoint(pt0_fit);
				s_cdm += m_control_space->getMetricAtPoint(pt1_fit);
				s_cdm += m_control_space->getMetricAtPoint(pt2_fit);
				p_cdm = m_reparameterization->countParameterizationMatrix(
					m_control_space->getBaseSurface()->getPoint(middle));
			}else if(ControlSpace2d::param_cached_parameterization_matrix){
				ControlDataMatrix2d pm_cdm;
				s_cdm = m_control_space->getMetricAndParameterizationAtPoint(pt0_fit, pm_cdm);
				p_cdm = pm_cdm;
				s_cdm += m_control_space->getMetricAndParameterizationAtPoint(pt1_fit, pm_cdm);
				p_cdm += pm_cdm;
				s_cdm += m_control_space->getMetricAndParameterizationAtPoint(pt2_fit, pm_cdm);
				p_cdm += pm_cdm;
				p_cdm *= (1.0/3.0);
			}else{
				s_cdm = m_control_space->getMetricAtPoint(pt0_fit);
				s_cdm += m_control_space->getMetricAtPoint(pt1_fit);
				s_cdm += m_control_space->getMetricAtPoint(pt2_fit);
				double g_ratio;
				SurfaceConstPtr surface = m_control_space->getBaseSurface();
				p_cdm = surface->countParameterizationMatrix(middle, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt0_fit, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt1_fit, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt2_fit, g_ratio);
			}
			s_cdm *= (1.0/3.0);
			break;
		}
	case MeshData::QM_MIDEDGES_MIN:
	case MeshData::QM_MIDEDGES_AVE:
		{
			const DPoint2d pt0_fit = m_control_space->fitInPoint(
				DPoint2d::average(mpt0->getCoordinates(), mpt1->getCoordinates()));
			const DPoint2d pt1_fit = m_control_space->fitInPoint(
				DPoint2d::average(mpt1->getCoordinates(), mpt2->getCoordinates()));
			const DPoint2d pt2_fit = m_control_space->fitInPoint(
				DPoint2d::average(mpt2->getCoordinates(), mpt0->getCoordinates()));

			s_cdm = m_control_space->getMetricAtPoint(pt0_fit);
			if(MeshTriangle2d::param_quality_metric_selection == MeshData::QM_MIDEDGES_MIN){
				s_cdm.setMinimum(m_control_space->getMetricAtPoint(pt1_fit));
				s_cdm.setMinimum(m_control_space->getMetricAtPoint(pt2_fit));
			}else{
				s_cdm += m_control_space->getMetricAtPoint(pt1_fit);
				s_cdm += m_control_space->getMetricAtPoint(pt2_fit);
				s_cdm *= (1.0/3.0);
			}

			if(m_reparameterization){
				p_cdm = m_reparameterization->countParameterizationMatrix(
					m_control_space->getBaseSurface()->getPoint(middle));
			}else if(ControlSpace2d::param_cached_parameterization_matrix){
				m_control_space->getMetricAndParameterizationAtPoint(middle, p_cdm);
			}else{
				double g_ratio;
				SurfaceConstPtr surface = m_control_space->getBaseSurface();
				p_cdm = surface->countParameterizationMatrix(middle, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt0_fit, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt1_fit, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt2_fit, g_ratio);
			}
			break;
		}
	}			

	m_metric.setData(s_cdm, p_cdm);
}

void Metric2dContext::countMetricAtPoints(const MeshPoint2d* mpt0, const MeshPoint2d* mpt1)
{
	assert(m_control_space);

	const DPoint2d middle = m_control_space->fitInPoint(
		DPoint2d::average(mpt0->getCoordinates(), mpt1->getCoordinates()));
	ControlDataMatrix2d s_cdm, p_cdm;

	switch(MeshTriangle2d::param_quality_metric_selection){
	case MeshData::QM_MIDDLE:
		{
			if(m_reparameterization){
				p_cdm = m_reparameterization->countParameterizationMatrix(
					m_control_space->getBaseSurface()->getPoint(middle));
				s_cdm = m_control_space->getMetricAtPoint(middle);
			}else if(ControlSpace2d::param_cached_parameterization_matrix){
				s_cdm = m_control_space->getMetricAndParameterizationAtPoint(middle, p_cdm);
			}else{
				double g_ratio;
				p_cdm = m_control_space->getBaseSurface()->countParameterizationMatrix(middle, g_ratio);
				s_cdm = m_control_space->getMetricAtPoint(middle);
			}
			break;
		}
	case MeshData::QM_MIDDLE_AND_VERTICES_MIN:
	case MeshData::QM_MIDDLE_AND_VERTICES_AVE:
		{
			const DPoint2d pt0_fit = m_control_space->fitInPoint(mpt0->getCoordinates());
			const DPoint2d pt1_fit = m_control_space->fitInPoint(mpt1->getCoordinates());

			s_cdm = m_control_space->getMetricAtPoint(pt0_fit);
			if(MeshTriangle2d::param_quality_metric_selection == MeshData::QM_MIDDLE_AND_VERTICES_MIN){
				s_cdm.setMinimum(m_control_space->getMetricAtPoint(pt1_fit));
			}else{
				s_cdm += m_control_space->getMetricAtPoint(pt1_fit);
			}

			ControlDataMatrix2d sm_cdm;
			if(m_reparameterization){
				p_cdm = m_reparameterization->countParameterizationMatrix(
					m_control_space->getBaseSurface()->getPoint(middle));
				sm_cdm = m_control_space->getMetricAtPoint(middle);
			}else if(ControlSpace2d::param_cached_parameterization_matrix){
				sm_cdm = m_control_space->getMetricAndParameterizationAtPoint(middle, p_cdm);
			}else{
				double g_ratio;
				SurfaceConstPtr surface = m_control_space->getBaseSurface();
				p_cdm = surface->countParameterizationMatrix(middle, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt0_fit, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt1_fit, g_ratio);
				sm_cdm = m_control_space->getMetricAtPoint(middle);
			}
			if(MeshTriangle2d::param_quality_metric_selection == MeshData::QM_MIDDLE_AND_VERTICES_MIN){
				s_cdm.setMinimum(sm_cdm);
			}else{
				s_cdm += sm_cdm;
				s_cdm *= (1.0/3.0);
			}
			break;
		}
	case MeshData::QM_VERTICES_AVE:
		{
			const DPoint2d pt0_fit = m_control_space->fitInPoint(mpt0->getCoordinates());
			const DPoint2d pt1_fit = m_control_space->fitInPoint(mpt1->getCoordinates());

			if(m_reparameterization){
				p_cdm = m_reparameterization->countParameterizationMatrix(
					m_control_space->getBaseSurface()->getPoint(middle));
				s_cdm = m_control_space->getMetricAtPoint(pt0_fit);
				s_cdm += m_control_space->getMetricAtPoint(pt1_fit);
			}else if(ControlSpace2d::param_cached_parameterization_matrix){
				ControlDataMatrix2d pm_cdm;
				s_cdm = m_control_space->getMetricAndParameterizationAtPoint(pt0_fit, pm_cdm);
				p_cdm = pm_cdm;
				s_cdm += m_control_space->getMetricAndParameterizationAtPoint(pt1_fit, pm_cdm);
				p_cdm += pm_cdm;
				p_cdm *= 0.5;
			}else{
				s_cdm = m_control_space->getMetricAtPoint(pt0_fit);
				s_cdm += m_control_space->getMetricAtPoint(pt1_fit);
				double g_ratio;
				SurfaceConstPtr surface = m_control_space->getBaseSurface();
				p_cdm = surface->countParameterizationMatrix(middle, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt0_fit, g_ratio);
				if(g_ratio < MIN_PARAM_GRATIO) p_cdm = surface->countParameterizationMatrix(pt1_fit, g_ratio);
			}
			s_cdm *= 0.5;
			break;
		}
	default:
		assert(false);
		break;
	}

	m_metric.setData(s_cdm, p_cdm);
}

double Metric2dContext::getMetricGradation(const DPoint2d& pt) const
{
	if(!m_control_space || !m_control_space->isAdaptive()) return 100.0;
	else return m_control_space->getAsAdaptive()->getMetricGradationRatio(pt);
}

const DMPoint2d Metric2dContext::transformPStoMS(const DPoint2d& pt) {
	if(m_reparameterization){
		return m_metric.transformPStoMS(
			m_reparameterization->getParameters(
				m_control_space->getBaseSurface()->getPoint(pt)));
	}else 
		return m_metric.transformPStoMS(pt);
}

/// Returns coordinates of point transformed from parametric to real space
const DPoint2d Metric2dContext::transformPStoRS(const DPoint2d& pt) {
	if(m_reparameterization){
		return m_metric.transformPStoRS(
			m_reparameterization->getParameters(
				m_control_space->getBaseSurface()->getPoint(pt)));
	}else 
		return m_metric.transformPStoRS(pt);
}

/// Returns coordinates of point transformed from real to parametric space
const DPoint2d Metric2dContext::transformRStoPS(const DPoint2d& pt) {
	if(m_reparameterization){
		return m_control_space->getBaseSurface()->getParametersNear(
				m_reparameterization->getPoint(m_metric.transformRStoPS(pt)),
				m_reparameterization->getRemoteBasePoint());
	}else 
		return m_metric.transformRStoPS(pt);
}

/// Returns coordinates of point transformed from metric to parametric space
const DPoint2d Metric2dContext::transformMStoPS(const DMPoint2d& pt) {
	if(m_reparameterization){
		return m_control_space->getBaseSurface()->getParametersNear(
				m_reparameterization->getPoint(m_metric.transformMStoPS(pt)),
				m_reparameterization->getRemoteBasePoint());
	}else 
		return m_metric.transformMStoPS(pt);
}

/// Returns coordinates of point transformed from parametric to metric space
const DMVector2d Metric2dContext::transformPStoMS(const DVector2d& v) {
	if(m_reparameterization){
		return m_metric.transformPStoMS(
			m_reparameterization->getParameters(
				m_control_space->getBaseSurface()->getPoint(
					m_reparameterization->getRemoteBasePoint()+v)) 
				- m_reparameterization->getLocalBasePoint());
	}else 
		return m_metric.transformPStoMS(v);
}

/// Returns coordinates of point transformed from parametric to real space
const DVector2d Metric2dContext::transformPStoRS(const DVector2d& v) {
	if(m_reparameterization){
		return m_metric.transformPStoRS(
			m_reparameterization->getParameters(
				m_control_space->getBaseSurface()->getPoint(
					m_reparameterization->getRemoteBasePoint()+v))
				- m_reparameterization->getLocalBasePoint());
	}else 
		return m_metric.transformPStoRS(v);
}

/// Returns coordinates of point transformed from real to parametric space
const DVector2d Metric2dContext::transformRStoPS(const DVector2d& v) {
	if(m_reparameterization){
		return m_control_space->getBaseSurface()->getParametersNear(
				m_reparameterization->getPoint(
					m_reparameterization->getLocalBasePoint() + m_metric.transformRStoPS(v)),
				m_reparameterization->getRemoteBasePoint())
			- m_reparameterization->getRemoteBasePoint();
	}else 
		return m_metric.transformRStoPS(v);
}

/// Returns coordinates of point transformed from metric to parametric space
const DVector2d Metric2dContext::transformMStoPS(const DMVector2d& v) {
	if(m_reparameterization){
		return m_control_space->getBaseSurface()->getParametersNear(
				m_reparameterization->getPoint(
					m_reparameterization->getLocalBasePoint() + m_metric.transformMStoPS(v)),
				m_reparameterization->getRemoteBasePoint())
			- m_reparameterization->getRemoteBasePoint();
	}else 
		return m_metric.transformMStoPS(v);
}
