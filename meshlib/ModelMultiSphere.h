#pragma once

#ifndef MODELMULTISPHERE_H__INCLUDED
#define MODELMULTISPHERE_H__INCLUDED

#include "MarchingCubes.h"
#include "DSphere.h"

class ModelMultiSphere : public MarchingModel
{
public:
	ModelMultiSphere( int sphere_count );
public:
	/// returns bounding box for model geometry
	virtual DBox getBoundingBox() const override;
	/// returns the discretization-grid resolution (min length)
	virtual double getResolution() const override;
	/// returns value of the scalar field for the given point
	virtual double isoValue(const DPoint3d& pt, int * sub_surface_id = nullptr ) const override;
	/// move the givent mesh point onto model surface with the given id
	virtual void moveToSurface( MeshPoint3d* point, int sub_surface_id ) const override;
private:
	DataVector<DSphere> m_spheres;
};

#endif // MODELMULTISPHERE_H__INCLUDED