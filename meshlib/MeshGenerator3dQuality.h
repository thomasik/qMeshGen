/////////////////////////////////////////////////////////////////////////////
// MeshGenerator3dQuality.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR3DQUALITY_H__INCLUDED)
#define MESHGENERATOR3DQUALITY_H__INCLUDED

#include "MeshData.h"
#include "Metric3dContext.h"
#include "TagExtended.h"

class MeshPoint3d;
class MeshContainer3d;
class MeshTetrahedron;
/**
 * This class gathers several procedures responsible for quality checking/improving
 * for 3D tetrahedral meshes.
 */
class MeshGenerator3dQuality  
{
public:
	/// Insert new point (using Delaunay) for refinement of the given tetrahedra
	static bool insertPointForTetrahedronRefinement(Metric3dContext& mc, MeshContainer3d* mesh, 
		MeshTetrahedron* tetrahedron);
	/// Generate (and log) statistics about mesh
	static void statMesh(Metric3dContext& mc, MeshContainer3d* mesh, const string& caption = "stats");
	/// Attempts to improve tetrahedra with angles below the given threshold, using local modifications
	static int optimizeTetrahedraForMinAngleMod(Metric3dContext& mc, 
		MeshContainer3d* mesh, double threshold = 0.173648); // ~ sin(10o)
	/// Attempts to improve tetrahedra with angles below the given threshold, using simple optimization
	static int optimizeTetrahedraForMinAngle(MeshContainer3d* mesh, double threshold = 0.173648); // ~ sin(10o)
	/// Attempts to improve tetrahedra with quality below the given threshold, using simple optimization
	static int optimizeBadTetrahedra(Metric3dContext& mc, MeshContainer3d* mesh, double threshold = 0.3);
	/// Attempts to improve tetrahedra with point relocation using simple optimization
	static bool movePointUsingSimpleOptimization(Metric3dContext& mc, MeshPoint3d* point);
	/// Attempts to improve tetrahedra with point relocation using simple optimization
	static bool movePointUsingAngleOptimization(MeshPoint3d* point);
	/// Attempts to improve tetrahedra with point relocation using simple optimization
	static bool optimizeNearBoundarySimple(Metric3dContext& mc, MeshContainer3d* mesh);
	/// Attempts to improve the given tetrahedra by cavity retriangulation
	static bool improveTetrahedron(Metric3dContext& mc, MeshContainer3d* mesh, MeshTetrahedron* tetra);
	/// Calls several tetrahedra-smoothing procedures ("steps" times) for the discretization of this volume
	static bool smoothen(Metric3dContext& mc, MeshContainer3d* mesh, int steps = 1, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1, bool with_boundary = false, 
		int method = MeshData::SM3_LAPLACE_MIXED | MeshData::SM3_SWAP_COMPLETE);
	/// Calls the swap-smoothing procedure for all tetrahedra in the given mesh (ready for parallel openmp processing)
	static bool smoothenSwapParallel(Metric3dContext& mc, MeshContainer3d* mesh);
	/// Calls the swap-smoothing procedure for all tetrahedra in the given mesh (complete)
	static bool smoothenSwapComplete(Metric3dContext& mc, MeshContainer3d* mesh, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1,
		bool with_boundary = false);
	/// Calls the swap-smoothing procedure for all tetrahedra in the given mesh
	static bool smoothenSwap(Metric3dContext& mc, MeshContainer3d* mesh,
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Tries to move point to new coordinates (gradually, checking for inverted blocks)
	static bool tryMovingPoint(Metric3dContext& mc, MeshPoint3d *point, const DPoint3d& new_pt, 
		bool no_boundary = true);
	/// Tries to improve the mesh locally by moving the given point to the barycenter of incident points
	static bool movePointByLaplaceForVariableMetric(Metric3dContext& mc, MeshPoint3d *point);
	/// Tries to improve the mesh locally by moving the given (boundary) point to the barycenter of incident points
	static bool moveBoundaryPointByLaplaceForVariableMetric(Metric3dContext& mc, MeshPoint3d* point);
	/// Auxiliary function for calculating minimum angle sine for set of faces and vertex
	static double calculateCavityMinAngleSc(const DPoint3d & dpt, 
		DataVector<DVector3d> tri_fvec[3], DataVector<DPoint3d> tri_fpts[3], 
		const DataVector<DVector3d> & normals, DataVector<int> & min_i);
	/// Auxiliary function for calcuating optimum vector for optimization
	static DVector3d calculateCavityOptVec(const DPoint3d & dpt,
		const DataVector<DPoint3d> & opt_vertices, const DataVector<int> & min_i);
	/// Tries to improve the mesh locally by moving the given point using optimization-based smoothing
	static double movePointOpt(Metric3dContext& mc, MeshPoint3d *point, double threshold = 1.0);
	/// Tries to improve the mesh locally by moving the given point to the barycenter of incident points
	static bool movePointByLaplace(Metric3dContext& mc, MeshPoint3d *point);
	/// Tries to improve the mesh locally by moving the given point to the barycenter of incident points
	static bool moveBoundaryPointByLaplace(Metric3dContext& mc, MeshPoint3d *point);
	/// Calls the optimization-smoothing procedure for all points in the given mesh
	static bool smoothenOpt(Metric3dContext& mc, MeshContainer3d* mesh,
			TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Calls the optimization-smoothing procedure for all points in the given mesh
	static bool smoothenOptMixed(Metric3dContext& mc, MeshContainer3d* mesh,
			TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Calls the Laplacian-smoothing procedure for all points in the given mesh
	static bool smoothenLaplaceMixed(Metric3dContext& mc, MeshContainer3d* mesh, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1,
		bool with_boundary = false);
	/// Calls the Laplacian-smoothing procedure for all points in the given mesh
	static bool smoothenLaplace(Metric3dContext& mc, MeshContainer3d* mesh, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1,
		bool with_boundary = false);
	/// Calls smoothing procedures for all faces in the given description of domain several times
	static bool smoothenBlocks(MeshContainer3d* boundary, int steps = 1, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1,
		int method = MeshData::SM3_LAPLACE_MIXED | MeshData::SM3_SWAP_COMPLETE);
	/// Performs some correctness-checking tests for mesh-consistency
	static void testTetrahedralMesh(const MeshContainer3d* mesh);
	/// Removes slivers on boundary surface
	static int removeBoundarySlivers(Metric3dContext& mc, MeshContainer3d* mesh, 
		double mean_ratio_threshold = 0.2, 
		TagExtended::TagType tag_type = TagExtended::TAG_NONE, int tag_value = 1);
	/// Prints statistical information about minimum inner angles (dihedral)
	static double statMinDihedralAngles(MeshContainer3d* mesh);
protected:
	class PointNode {
	public:
		PointNode(const DPoint3d& _pt = DPoint3d::zero, double q = 1.0) : coord(_pt), quality(q) {}
		bool operator<( const PointNode& rhs ) const { return rhs.quality < quality; } // decreasing order...
		DPoint3d coord;
		double quality;
	};
public:
	/// Criterion used for swapPIng method of quality improvement
	static int param_swap_criterion;
};

#endif // !defined(MESHGENERATOR3DQUALITY_H__INCLUDED)
