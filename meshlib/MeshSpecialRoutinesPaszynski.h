/////////////////////////////////////////////////////////////////////////////
// MeshSpecialRoutinesPaszynski.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2009-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHSPECIALROUTINESPASZYNSKI_H__INCLUDED)
#define MESHSPECIALROUTINESPASZYNSKI_H__INCLUDED

class MeshContainer3d;

/**
 * This class implements several methods required for the cooperation with
 * adaptive hp-solver (cooperation with Prof. R.Schaefer and M.Paszynski) -> 
 *		txt-file format for mesh, etc.
 */
class MeshSpecialRoutinesPaszynski
{
public:
	static MeshContainer3d* loadTetrahedralMesh(const string& fname);
	static bool storeTetrahedralMesh(const string& fname, const string& header, MeshContainer3d* mesh);
};

#endif // !defined(MESHSPECIALROUTINESPASZYNSKI_H__INCLUDED)
