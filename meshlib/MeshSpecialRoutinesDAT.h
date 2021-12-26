/////////////////////////////////////////////////////////////////////////////
// MeshTetrahedron.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHSPECIALROUTINESDAT_H__INCLUDED)
#define MESHSPECIALROUTINESDAT_H__INCLUDED

class MeshContainer2d;
class MeshContainer3d;
class MeshContainer3dSurface;

/**
 * This class implements several methods required for the cooperation with
 * fem-solver (developed by Prof. M.PIetrzyk) -> 
 *		adaptation, DAT-file format, trasfer of values, special geometries, etc.
 */
class MeshSpecialRoutinesDAT  
{
public:
	static bool testHessianCurvature(const string& fname, const string& equation, const string& equation2 = "", bool with_quads = false);
	static bool interpolateValues(MeshContainer2d* dest_mesh, MeshContainer2d* source_mesh);
//	static MeshContainer2d* importFromDAT(const string& fname);
//	static MeshContainer2d* importFromDAT(const string& fname_p, const string& fname_i);
//	static bool exportToDAT(MeshContainer2d* mesh, const string& fname, bool mark_boundary_points = false);
//	static MeshContainer3d* createBoundaryForPhaseTransformation(double a, double r, double dr1, double dr2, double dl1, double dl2, bool with_user_control = true, int shape_type = 0);
public:
	static MeshContainer3d* loadSurfaceBoundaryMesh(const string& fname);
	static MeshContainer3d* loadSpecialGrid(const char* fname);
	static double getMaxDiffError(MeshContainer2d* dat_mesh, double tinit);
	static double getMaxOscillationError(MeshContainer2d* dat_mesh);
	static double getMaxSecantError(MeshContainer2d* dat_mesh);
//	static bool checkPhaseCorrectness(MeshContainer2d* mesh);
	static bool setControlSpaceFromDATHesjan(MeshContainer2d* boundary, MeshContainer2d* value_mesh, double h_min_ratio = -1.0, double h_factor = -1.0);
	static double getPhaseArea(double r, double a);
	static double countNextPhaseRadius(double old_r, double a, double ds, int shape_type = 0);
	static double coefficient;
	static bool runSpecialGeometryPhase(int hesjan_control=-1, int from_step=0);
};

#endif // !defined(MESHSPECIALROUTINESDAT_H__INCLUDED)
