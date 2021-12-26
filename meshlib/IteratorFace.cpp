/////////////////////////////////////////////////////////////////////////////
// IteratorFace.cpp
// Class describing iterator for browsing 3d faces in the mesh
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2007-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "IteratorFace.h"
#include "MeshContainer3d.h"
#include "MeshBlock.h"
#include "MeshFace.h"
#include "MeshDomainVolume.h"

IteratorFace::IteratorFace(const MeshContainer3d* mesh) :
	current_face(nullptr), current_block(nullptr), domain_mesh(mesh), 
   	bindex(-1), fct(0), findex(-1)
{
	assert(domain_mesh);
	bct = mesh->getBlocksCount();
	nextFace();
}

IteratorFace& IteratorFace::nextFace()
{
	bool already_checked = false;
	do{
		if(++findex < fct){	// next face for block
			assert(current_block);
			current_face = current_block->getFace(findex);
		}else{
			if(++bindex < bct){ // next block, first face
				current_block = domain_mesh->getBlockAt(bindex);
				fct = current_block->getFaceCount();
				findex = 0;
				if(fct == 0) continue;
				current_face = current_block->getFace(findex);
			}else{ // end of faces
				current_block = nullptr;
				current_face = nullptr;
				return *this;
			}
		}

		already_checked = current_face->getOtherBlock(current_block) != nullptr &&
			current_face->getBlockIndex(current_block) != 0;
	}while(already_checked);

	return *this;
}
