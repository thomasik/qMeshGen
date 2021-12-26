/////////////////////////////////////////////////////////////////////////////
// MeshBRepRSurf.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHBREPRSURF_H__INCLUDED)
#define MESHBREPRSURF_H__INCLUDED

#include "MeshBRep.h"

class MeshContainer2d;

/**
 * Class responsible for parsing the description of the mesh in ResuSurf format.
 */
class MeshBRepRSurf : public MeshBRep
{
public:
	/// Parse model description data in brep (resu.surf) .txt format
	bool parseFile(const string& fname);
	/// Stores the description of mesh-surface to resu.surf file
	static bool storeFileResuSurf(const string& fname, const MeshContainer2d* mesh);
};

#endif // !defined(MESHBREPRSURF_H__INCLUDED)
