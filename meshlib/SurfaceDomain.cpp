// SurfaceDomain.cpp: implementation of the SurfaceDomain class.
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2014-
//	Generation of unstructured meshes
//////////////////////////////////////////////////////////////////////

#include "SurfaceDomain.h"
#include "MeshViewSet.h"

double SurfaceDomain::w_inside_ratio = 0.0; //  1 -> use for weighted quality
double SurfaceDomain::w_approx_error = 1.0;

bool SurfaceDomain::isInsideOBBox(const DPoint3d& pt) const {
	if(m_box.valid) return m_box.contains(pt);
	else return true;
}

void SurfaceDomain::draw(MeshViewSet* set, const SurfaceParametric * /* surface */, int /* id */) const
{
	set->addInfo("domain type", "SD-true" );
}

