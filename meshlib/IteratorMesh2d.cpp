/////////////////////////////////////////////////////////////////////////////
// IteratorMesh2d.cpp
// Class describing iterator for browsing 2d meshes in the domain
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "IteratorMesh2d.h"
#include "MeshContainer3d.h"
#include "MeshDomainVolume.h"
#include "MeshDomainSurface.h"
#include "common.h"

IteratorMesh2d::IteratorMesh2d(MeshContainer3d* domain) :
	current_domain_volume(nullptr), current_domain_surface(nullptr),
	current_mesh(nullptr), domain_mesh(domain), 
	volume_index(-1), surface_index(0), surface_count(0)
{
	assert(domain_mesh);
	volume_count = domain_mesh->getBlocksCount();
	nextValidMesh();
}

IteratorMesh2d& IteratorMesh2d::nextValidMesh()
{
	do{
		if(++surface_index >= surface_count){
			// end of volume
			if(++volume_index < volume_count){
				current_domain_volume = (MeshDomainVolume*)domain_mesh->getBlockAt(volume_index);
				assert(current_domain_volume && (current_domain_volume->getType() == BLOCK_DOMAIN));
				surface_count = current_domain_volume->getFaceCount();
				if(surface_count > 0){
					current_domain_surface = (MeshDomainSurface*)current_domain_volume->getFace(surface_index=0);
					current_mesh = current_domain_surface->getMesh();
				}
			}else{
				current_domain_volume = nullptr;
				current_domain_surface = nullptr;
				current_mesh = nullptr;	// end of domain
				return *this;
			}
		}else{
			current_domain_volume = (MeshDomainVolume*)domain_mesh->getBlockAt(volume_index);
			current_domain_surface = (MeshDomainSurface*)current_domain_volume->getFace(surface_index);
			current_mesh = current_domain_surface->getMesh();
		}
	}while(!current_mesh);

	return *this;
}

IteratorBoundary2d::IteratorBoundary2d(MeshContainer3d* domain) :
	current_domain_volume(nullptr), current_domain_surface(nullptr), 
	current_boundary(nullptr), domain_mesh(domain), 
	volume_index(-1), surface_index(0), surface_count(0)
{
	assert(domain_mesh);
	volume_count = domain_mesh->getBlocksCount();
	nextValidBoundary();
}

IteratorBoundary2d& IteratorBoundary2d::nextValidBoundary()
{
	do{
		if(++surface_index >= surface_count){
			// end of volume
			if(++volume_index < volume_count){
				current_domain_volume = (MeshDomainVolume*)domain_mesh->getBlockAt(volume_index);
				assert(current_domain_volume && (current_domain_volume->getType() == BLOCK_DOMAIN));
				surface_count = current_domain_volume->getFaceCount();
				if(surface_count > 0){
					current_domain_surface = (MeshDomainSurface*)current_domain_volume->getFace(surface_index=0);
					current_boundary = current_domain_surface->getBoundary();
				}
			}else{
				current_domain_volume = nullptr;
				current_domain_surface = nullptr;
				current_boundary = nullptr;	// end of domain
				return *this;
			}
		}else{
			current_domain_volume = (MeshDomainVolume*)domain_mesh->getBlockAt(volume_index);
			current_domain_surface = (MeshDomainSurface*)current_domain_volume->getFace(surface_index);
			current_boundary = current_domain_surface->getBoundary();
		}
	}while(!current_boundary);

	return *this;
}
