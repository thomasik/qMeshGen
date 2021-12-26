#include "ModelMultiSphere.h"
#include "ControlSpace2dAdaptive.h"
#include "MeshPoint3d.h"

ModelMultiSphere::ModelMultiSphere( int sphere_count )
	: m_spheres( sphere_count )
{
	assert( sphere_count > 0 );

	RandomGen rg( 11 );
	static const double MIN_R = 0.1;
	static const double MAX_R = 1.0;

	m_spheres.add( DSphere( DPoint3d::zero, rg.doub(MIN_R, MAX_R) ) );

	for(int i = 1; i < sphere_count; i++ ){
		int j = (i == 1) ? 0 : (rg.int32() % i);
		const DSphere& ds = m_spheres[j];
		double r = rg.doub(MIN_R, MAX_R);
		DVector3d dv;
		do{
			dv.x = rg.doub(-1.0, 1.0);
			dv.y = rg.doub(-1.0, 1.0);
			dv.z = rg.doub(-1.0, 1.0);
		}while( dv.length2() == 0.0 );
		double ri = ds.getRadius();
		double vmin = std::max(ri, r) - std::min(ri, r);
		double vmax = ri + r;
		double margin = (vmax-vmin) * 0.1;
		double vlen = rg.doub( vmin+margin, vmax-margin);
		m_spheres.add( DSphere( ds.getCenter() + dv.normalized() * vlen, r ) );
	}
}


/// returns bounding box for model geometry
DBox ModelMultiSphere::getBoundingBox() const
{
	DBox box;
	for(size_t i = 0; i < m_spheres.countInt(); i++)
		box.addBox( m_spheres[i].getBoundingBox() );
	return box;
}

/// returns the discretization-grid resolution (min length)
double ModelMultiSphere::getResolution() const
{
	if( m_spheres.empty() ) return 0.0;
	double min_r = m_spheres[0].getRadius();
	for(size_t i = 1; i < m_spheres.countInt(); i++) {
		double r = m_spheres[i].getRadius();
		if( r < min_r ) min_r = r;
	}
	return min_r * ControlSpace2dAdaptive::param_curvature_ratio;
}

/// returns value of the scalar field for the given point
double ModelMultiSphere::isoValue(const DPoint3d& pt, int * sub_surface_id) const
{
	if( m_spheres.empty() ) return 0.0;
	size_t id = 0;
	double min_v = m_spheres[0].implicitValue( pt );
	for(size_t i = 1; i < m_spheres.countInt(); i++) {
		double v = m_spheres[i].implicitValue( pt );
		if( v < min_v ){ min_v = v; id = i; }
	}
	if( sub_surface_id != nullptr ) *sub_surface_id = (int)id;
	return min_v;
}

/// move the givent mesh point onto model surface with the given id
void ModelMultiSphere::moveToSurface( MeshPoint3d* point, int sub_surface_id ) const 
{
	if( sub_surface_id < 0 || (size_t)sub_surface_id >= m_spheres.countInt() ) return;
	const DSphere& ds = m_spheres[sub_surface_id];
	point->setCoordinates( ds.projectToSurface( point->getCoordinates() ) );
}
