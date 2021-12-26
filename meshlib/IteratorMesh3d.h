/////////////////////////////////////////////////////////////////////////////
// IteratorMesh3d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(ITERATORMESH3D_H__INCLUDED)
#define ITERATORMESH3D_H__INCLUDED

class MeshContainer3d;
class MeshDomainVolume;

class IteratorMesh3d
{
public:
	/// Constructor
	IteratorMesh3d(MeshContainer3d* domain);
	/// Mesh for current iterator
	MeshContainer3d* getMesh() { return current_mesh; }
	/// Domain-surface for current iterator
	MeshDomainVolume* getDomainVolume() { return current_domain_volume; }
	/// Get next valid mesh3d
	IteratorMesh3d& nextValidMesh();
	/// for end-checking
	bool isValid() const { return volume_index < volume_count; }
private:
	MeshDomainVolume* current_domain_volume;
	MeshContainer3d* current_mesh;
	MeshContainer3d* domain_mesh;
	int volume_index;
	int volume_count;
};

#endif //!defined(ITERATORMESH3D_H__INCLUDED)
