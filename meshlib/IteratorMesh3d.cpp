/////////////////////////////////////////////////////////////////////////////
// IteratorMesh3d.cpp
// Class describing iterator for browsing 3d meshes in the domain
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "IteratorMesh3d.h"
#include "MeshContainer3d.h"
#include "MeshDomainVolume.h"
#include "common.h"

IteratorMesh3d::IteratorMesh3d(MeshContainer3d* domain) :
	current_domain_volume(nullptr), current_mesh(nullptr), domain_mesh(domain), volume_index(-1)
{
	assert(domain_mesh);
	volume_count = domain_mesh->getBlocksCount();
	nextValidMesh();
}

IteratorMesh3d& IteratorMesh3d::nextValidMesh()
{
	do{
		if(++volume_index < volume_count){
			current_domain_volume = (MeshDomainVolume*)domain_mesh->getBlockAt(volume_index);
			assert(current_domain_volume && (current_domain_volume->getType() == BLOCK_DOMAIN));
			current_mesh = current_domain_volume->getMesh();
		}else{
			current_domain_volume = nullptr;	// end of domain
			return *this;
		}
	}while(!current_mesh);

	return *this;
}
