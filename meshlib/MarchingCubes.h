#pragma once

#ifndef MARCHINGCUBES_H__INCLUDED
#define MARCHINGCUBES_H__INCLUDED

#include "DPoint.h"
#include "DRect.h"

class MeshContainer3dSurface;
class MeshPoint3d;

class MarchingModel
{
public:
	/// virtual destructor
	virtual ~MarchingModel() {}
	/// returns bounding box for model geometry
	virtual DBox getBoundingBox() const  = 0;
	/// returns the discretization-grid resolution (min length)
	virtual double getResolution() const = 0;
	/// returns value of the scalar field for the given point
	virtual double isoValue(const DPoint3d& pt, int * sub_surface_id = nullptr ) const = 0;
	/// move the givent mesh point onto model surface with the given id
	virtual void moveToSurface( MeshPoint3d* point, int sub_surface_id ) const = 0;
};

class MarchingModelDiscreteLayered
{
public:
	/// virtual destructor
	virtual ~MarchingModelDiscreteLayered() {}
public:
	virtual void getCellCount(int & nx, int & ny) const = 0;
	/// returns bounding box for model geometry
	virtual DRect getLayerBoundingRect() const = 0;
	/// returns bounding box for model geometry
	virtual double getZ0() const = 0;
	/// prepare next layer
	virtual bool prepareNextLayer(double & dz) = 0;
	/// get iso-value
	virtual double isoValue(int iy, int ix, int * sub_surface_id = nullptr) const = 0;
};

class MarchingCubes
{
protected:
	struct VertexData {
		VertexData(const DPoint3d& pt = DPoint3d::zero, const double & iv = 0.0) 
			: coord( pt ), isoValue( iv ), surfId(-1), mpoint( nullptr ) {}
		DPoint3d coord;
		double isoValue;
		int surfId;
		MeshPoint3d* mpoint;
	};
	struct EdgeData {
		EdgeData() : mpoint( nullptr ) {}
		MeshPoint3d* mpoint;
	};
public:
	static MeshContainer3dSurface* createSurfaceMesh( 
		const MarchingModel& model, double isoLevel = 0.0 );
	static MeshContainer3dSurface* createSurfaceMesh(
		MarchingModelDiscreteLayered& model);
protected:
	/// create and insert into mesh triangles for the given cell
	static int polygonise( MarchingCubes::VertexData* cell_vertices[8], 
		MarchingCubes::EdgeData* cell_edges[12],
		MeshContainer3dSurface* mesh, const double & isoLevel = 0.0 );
	/// interpolate point for edge, returns 0 if vd0, 1 if vf1, 2 if between
	static int vertexInterpolate(const VertexData& vd0, const VertexData& vd1, DPoint3d& pt, 
		const double& isoLevel = 0.0, int * surfId = nullptr);
	static MeshPoint3d* getRefMeshPoint(int ci, int i, MarchingCubes::VertexData* cell_vertices[8], 
		MarchingCubes::EdgeData* cell_edges[12] );
	static void markFaces( MeshContainer3dSurface* mesh );
	static void markBorders( MeshContainer3dSurface* mesh, const MarchingModel& model );
protected:
	static int EDGE_TABLE[256];
	static int TRI_TABLE[256][16];
	static int TWO_POWER[12];
	static int EDGE_VERTICES[12][2]; 
};

#endif