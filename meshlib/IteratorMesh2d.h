/////////////////////////////////////////////////////////////////////////////
// IteratorMesh2d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(ITERATORMESH2D_H__INCLUDED)
#define ITERATORMESH2D_H__INCLUDED

class MeshContainer2d;
class MeshContainer3d;
class MeshDomainSurface;
class MeshDomainVolume;

class IteratorMesh2d
{
public:
	/// Constructor
	IteratorMesh2d(MeshContainer3d* domain);
	/// Mesh for current iterator
	MeshContainer2d* getMesh() { return current_mesh; }
	/// Domain-surface for current iterator
	MeshDomainSurface* getDomainSurface() { return current_domain_surface; }
	/// Domain-volume for current iterator
	MeshDomainVolume* getDomainVolume() { return current_domain_volume; }
	/// Get next valid mesh2d
	IteratorMesh2d& nextValidMesh();
	/// for end-checking
	bool isValid() const { return volume_index < volume_count; }
private:
	MeshDomainVolume* current_domain_volume;
	MeshDomainSurface* current_domain_surface;
	MeshContainer2d* current_mesh;
	MeshContainer3d* domain_mesh;
	int volume_index;
	int volume_count;
	int surface_index;
	int surface_count;
};

class IteratorBoundary2d
{
public:
	/// Constructor
	IteratorBoundary2d(MeshContainer3d* domain);
	/// Mesh for current iterator
	MeshContainer2d* getBoundary() { return current_boundary; }
	/// Domain-surface for current iterator
	MeshDomainSurface* getDomainSurface() { return current_domain_surface; }
	/// Domain-volume for current iterator
	MeshDomainVolume* getDomainVolume() { return current_domain_volume; }
	/// Get next valid mesh2d
	IteratorBoundary2d& nextValidBoundary();
	/// for end-checking
	bool isValid() const { return volume_index < volume_count; }
private:
	MeshDomainVolume* current_domain_volume;
	MeshDomainSurface* current_domain_surface;
	MeshContainer2d* current_boundary;
	MeshContainer3d* domain_mesh;
	int volume_index;
	int volume_count;
	int surface_index;
	int surface_count;
};

#endif //!defined(ITERATORMESH2D_H__INCLUDED)
