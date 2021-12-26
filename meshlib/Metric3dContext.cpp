// Metric3dContext.cpp: implementation of the Metric3dContext class.
//
//////////////////////////////////////////////////////////////////////

#include "Metric3dContext.h"
#include "ControlSpace3d.h"
#include "MeshPoint3d.h"
#include "MeshTetrahedron.h"
#include "ControlSpace3dAdaptive.h"

void Metric3dContext::setSpecialMetric(const ControlDataMatrix3d& cdm)
{
	m_metric.setData(cdm);
	m_radius = -1.0;
}

void Metric3dContext::countMetricAtPoint(const DPoint3d& pt, bool if_needed_only)
{
	assert(m_control_space);
	if( if_needed_only && withinImpactRadius( pt ) ) return;

	const DPoint3d pt_fit = m_control_space->fitInPoint(pt);
	m_metric.setData(m_control_space->getMetricAtPoint(pt_fit));

	m_radius = 1.0;
	m_center = transformRStoMS(pt_fit);
}

void Metric3dContext::countMetricAtPoints(MeshPoint3d** points)
{
	assert(m_control_space);

	const DPoint3d middle = m_control_space->fitInPoint(
		DPoint3d::average(
			points[0]->getCoordinates(),
			points[1]->getCoordinates(),
			points[2]->getCoordinates(),
			points[3]->getCoordinates()));
	ControlDataMatrix3d cdm;

	switch(MeshTetrahedron::param_quality_metric_selection){
	case MeshData::QM_MIDDLE:
		{
			cdm = m_control_space->getMetricAtPoint(middle);
			break;
		}
	case MeshData::QM_MIDDLE_AND_VERTICES_MIN:
	case MeshData::QM_MIDDLE_AND_VERTICES_AVE:
		{
			const DPoint3d pts_fit[4] = {
				m_control_space->fitInPoint(points[0]->getCoordinates()),
				m_control_space->fitInPoint(points[1]->getCoordinates()),
				m_control_space->fitInPoint(points[2]->getCoordinates()),
				m_control_space->fitInPoint(points[3]->getCoordinates())};

			cdm = m_control_space->getMetricAtPoint(pts_fit[0]);
			if(MeshTetrahedron::param_quality_metric_selection == MeshData::QM_MIDDLE_AND_VERTICES_MIN){
				cdm.setMinimum(m_control_space->getMetricAtPoint(pts_fit[1]));
				cdm.setMinimum(m_control_space->getMetricAtPoint(pts_fit[2]));
				cdm.setMinimum(m_control_space->getMetricAtPoint(pts_fit[3]));
				cdm.setMinimum(m_control_space->getMetricAtPoint(middle));
			}else{
				cdm += m_control_space->getMetricAtPoint(pts_fit[1]);
				cdm += m_control_space->getMetricAtPoint(pts_fit[2]);
				cdm += m_control_space->getMetricAtPoint(pts_fit[3]);
				cdm += m_control_space->getMetricAtPoint(middle);
				cdm *= 0.2;
			}
			break;
		}
	case MeshData::QM_VERTICES_AVE:
		{
			const DPoint3d pts_fit[4] = {
				m_control_space->fitInPoint(points[0]->getCoordinates()),
				m_control_space->fitInPoint(points[1]->getCoordinates()),
				m_control_space->fitInPoint(points[2]->getCoordinates()),
				m_control_space->fitInPoint(points[3]->getCoordinates())};

			cdm = m_control_space->getMetricAtPoint(pts_fit[0]);
			cdm += m_control_space->getMetricAtPoint(pts_fit[1]);
			cdm += m_control_space->getMetricAtPoint(pts_fit[2]);
			cdm += m_control_space->getMetricAtPoint(pts_fit[3]);
			cdm *= 0.25;
			break;
		}
	case MeshData::QM_MIDEDGES_MIN:
	case MeshData::QM_MIDEDGES_AVE:
		{
			DPoint3d pts_fit[6];
			for(int i = 0, i0 = 0; i0 < 3; i0++){
				for(int i1 = i0+1; i1 < 4; i1++){
					pts_fit[i++] = m_control_space->fitInPoint(
						DPoint3d::average(points[i0]->getCoordinates(), points[i1]->getCoordinates()));
				}
			}

			cdm = m_control_space->getMetricAtPoint(pts_fit[0]);
			if(MeshTetrahedron::param_quality_metric_selection == MeshData::QM_MIDEDGES_MIN){
				for(int i = 1; i < 6; i++)
					cdm.setMinimum(m_control_space->getMetricAtPoint(pts_fit[i]));
			}else{
				for(int i = 1; i < 6; i++)
					cdm += m_control_space->getMetricAtPoint(pts_fit[i]);
				cdm *= (1.0/6.0);
			}
			break;
		}
	}			

	m_metric.setData(cdm);

	m_radius = 1.0;
	m_center = transformRStoMS(middle);
}

void Metric3dContext::countMetricAtPoints(MeshPoint3d* point0, MeshPoint3d* point1)
{
	assert(m_control_space);

	const DPoint3d middle = m_control_space->fitInPoint(
		DPoint3d::average(point0->getCoordinates(), point1->getCoordinates()));
	ControlDataMatrix3d cdm;

	switch(MeshTetrahedron::param_quality_metric_selection){
	case MeshData::QM_MIDDLE:
	case MeshData::QM_MIDEDGES_MIN:
	case MeshData::QM_MIDEDGES_AVE:
		{
			cdm = m_control_space->getMetricAtPoint(middle);
			break;
		}
	case MeshData::QM_MIDDLE_AND_VERTICES_MIN:
	case MeshData::QM_MIDDLE_AND_VERTICES_AVE:
		{
			const DPoint3d pts_fit[2] = {
				m_control_space->fitInPoint(point0->getCoordinates()),
				m_control_space->fitInPoint(point1->getCoordinates())};

			cdm = m_control_space->getMetricAtPoint(pts_fit[0]);
			if(MeshTetrahedron::param_quality_metric_selection == MeshData::QM_MIDDLE_AND_VERTICES_MIN){
				cdm.setMinimum(m_control_space->getMetricAtPoint(pts_fit[1]));
				cdm.setMinimum(m_control_space->getMetricAtPoint(middle));
			}else{
				cdm += m_control_space->getMetricAtPoint(pts_fit[1]);
				cdm += m_control_space->getMetricAtPoint(middle);
				cdm *= (1.0/3.0);
			}
			break;
		}
	case MeshData::QM_VERTICES_AVE:
		{
			const DPoint3d pts_fit[2] = {
				m_control_space->fitInPoint(point0->getCoordinates()),
				m_control_space->fitInPoint(point1->getCoordinates())};

			cdm = m_control_space->getMetricAtPoint(pts_fit[0]);
			cdm += m_control_space->getMetricAtPoint(pts_fit[1]);
			cdm *= 0.5;
			break;
		}
	}			

	m_metric.setData(cdm);

	m_radius = 1.0;
	m_center = transformRStoMS(middle);
}

double Metric3dContext::getMetricGradation(const DPoint3d& pt) const
{
	if(!m_control_space || !m_control_space->isAdaptive()) return 100.0;
	else return m_control_space->getAsAdaptive()->getMetricGradationRatio(pt);
}

/// Check point for impact sphere
bool Metric3dContext::withinImpactRadius(const DPoint3d& pt) const
{
	if(m_radius < 0.0) return false; // no impact sphere set

	return (m_center.distance( transformRStoMS(pt) ) < m_radius);
}
