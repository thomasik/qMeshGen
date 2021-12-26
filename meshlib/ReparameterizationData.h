/////////////////////////////////////////////////////////////////////////////
// ReparameterizationData.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(REPARAMETERIZATIONDATA_H__INCLUDED)
#define REPARAMETERIZATIONDATA_H__INCLUDED

#include "SurfaceParametric.h"

/**
 * Reparameterization info
 */
class ReparameterizationData
{
public:
	/// Standard constructor
	ReparameterizationData(SurfacePtr surface, const DPoint2d& remote_base_point, 
		const DPoint2d& local_base_point, double max_dist_squared) 
		: m_surface(surface), m_remote_base_point(remote_base_point), 
			m_local_base_point(local_base_point), m_max_dist_squared(max_dist_squared) 
	{ 
		assert(surface);
	}
public:
	ControlDataMatrix2d countParameterizationMatrix(const DPoint3d& pt) const {
		double g_ratio;
		return m_surface->countParameterizationMatrix(
			m_surface->getParameters(pt), g_ratio);
	}
	DPoint2d getParameters(const DPoint3d& pt) const {
		const DPoint2d local_pt = m_surface->getParametersNear(pt, m_local_base_point);
//		assert(pt.distance2(m_surface->getPoint(local_pt)) < m_max_dist_squared);
		return local_pt;
	}
	DPoint3d getPoint(const DPoint2d& pt) const {
		return m_surface->getPoint(pt);
	}
	DPoint2d getRemoteBasePoint() const { return m_remote_base_point; }
	DPoint2d getLocalBasePoint() const { return m_local_base_point; }
protected:
	/// Reparameterization surface
	SurfacePtr m_surface;
	/// Base point
	DPoint2d m_remote_base_point;
	/// Base point
	DPoint2d m_local_base_point;
	/// Maximum allowed discrepancy
	double m_max_dist_squared;
};

#endif // !defined(REPARAMETERIZATIONDATA_H__INCLUDED)
