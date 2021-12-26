// ControlSpace2d.cpp: implementation of the ControlSpace2d class.
//
//////////////////////////////////////////////////////////////////////

#include "ControlSpace2d.h"
#include "SurfaceParametric.h"
#include "DEquation.h"
#include "ControlSpace2dAdaptive.h"

int ControlSpace2d::param_control_type = MeshData::CONTROL_QUADTREE;

/// Ratio of inflating bounding box for control space
double ControlSpace2d::param_inflate_box_factor = 0.0;

/// Whether the approximated (cached) parameterization matrix should be used
/// 0 - no, 1 - from leaves, 2 - from nodes,
int ControlSpace2d::param_cached_parameterization_matrix = 2;


const ControlSpace2dAdaptive* ControlSpace2d::getAsAdaptive() const
{
	if (!isAdaptive()) return nullptr;
	return (ControlSpace2dAdaptive*)this;
}

ControlSpace2dAdaptive* ControlSpace2d::getAsAdaptive()
{
	if (!isAdaptive()) return nullptr;
	return (ControlSpace2dAdaptive*)this;
}

/// Get transformation and parameterization matrices at the given point
ControlDataMatrix2d ControlSpace2d::getMetricAndParameterizationAtPoint(
	const DPoint2d& pt, ControlDataMatrix2d& p_cdm) const
{
	if(base_surface){
		double g_ratio;
		p_cdm = base_surface->countParameterizationMatrix(pt, g_ratio);
		// problems if g_ratio < MIN_PARAM_GRATIO, but what solution ???
	}else{
		p_cdm.setIdentity();
	}
	return getMetricAtPoint(pt);
}
