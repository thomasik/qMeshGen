/////////////////////////////////////////////////////////////////////////////
// MeshSpecialRoutinesUTC.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2010-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHSPECIALROUTINESUTC_H__INCLUDED)
#define MESHSPECIALROUTINESUTC_H__INCLUDED

#include <memory>

#include "DTriangle.h"
#include "DataList.h"
#include "DataVector.h"

class MeshContainer2d;
class MeshContainer3d;
class ControlSpace3d;
class ControlSpace3dAdaptive;
class DBox;
class MeshBlock;
class Metric3dContext;

/**
 * This class implements several methods for decomposition/remeshing with cracks -> 
 *		2010
 */
class MeshSpecialRoutinesUTC  
{
public:
	/// retriangulate given mesh with next-step-CS and smoothing of crack surfaces
	static void remeshCrack(Metric3dContext& mc, MeshContainer3d* mesh, int mode = 0);
	/// identify set of 3d triangles approximating the surface of crack
	static int identifyCrackPlanes(Metric3dContext& mc, 
		const DataVector<MeshBlock*> & crack_blocks, double tolerance,
		DataVector<std::shared_ptr<DTriangle3d>> & crack_planes, 
		DataVector<int> & crack_indices);
	/// prepare (analytically, for now) control space for next simulation step
	static std::shared_ptr<ControlSpace3dAdaptive> prepareNextCS(MeshContainer3d* mesh, 
		const DataCompoundList<DTriangle3d>& crack, const DVector3d& dv);
	/// remove elements within the crack (described by set of (triangle) planes), return number of removed elements
	static int createCrack(Metric3dContext& mc, MeshContainer3d* mesh, const DataCompoundList<DTriangle3d>& crack);
	/// marks elements within the crack (described by box), return number of marked elements
	static int markCrack(MeshContainer3d* mesh, const DBox& box);
	/// collect all blocks from the mesh, crossed by the crack triangles + postprocessing
	static int gatherCrackBlocks(MeshContainer3d* mesh, const DataCompoundList<DTriangle3d>& crack, 
		DataVector<MeshBlock*> & crack_blocks);
	/// store binary file with 3d tetrahedral mesh
	static bool storeBinFile(const string& fname, MeshContainer3d* mesh);
	/// load binary file with 3d tetrahedral mesh
	static MeshContainer3d* readBinFile(const string& fname);
};

#endif // !defined(MESHSPECIALROUTINESUTC_H__INCLUDED)
