/////////////////////////////////////////////////////////////////////////////
// MeshGenerator2dQuad.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR2DQUAD_H__INCLUDED)
#define MESHGENERATOR2DQUAD_H__INCLUDED

#include "MeshData.h"
#include "Metric2dContext.h"

class MeshPoint2d;
class MeshElement;
class MeshTriangle2d;
class FrontLines;
class FrontEdge;
class MeshContainer2d;
class MeshContainer3d;

//#define QUAD_DEBUG_LOG

/**
 * This class gathers several procedures required by the LeeLo algorithm of 
 * triangle-to-quad conversion and several mesh-improvement routines.
 */
class MeshGenerator2dQuad  
{
public:
	struct QuadConversionContext{
		QuadConversionContext() : base_fedge(nullptr), base_element(nullptr), low_ratio(false), 
			qmorph_min_ratio(0.65), qmorph_max_ratio(SQRT3), qmorph_max_ratio2(3.0), 
			merge_method(MeshData::QUADS_UNKNOWN) {}
		bool selectFrontEdge(const FrontLines* front);
		void setMergeMethod(int m) { merge_method = m; }
		bool lowerQMorphRatio() {
			if(low_ratio) return false; // only once
			low_ratio = true;
			qmorph_min_ratio /= 2.0;
			qmorph_max_ratio *= 2.0;
			qmorph_max_ratio2 *= 4.0;
			return true;
		}
		bool init();
		void checkSideEdges();
	public:
		FrontEdge* base_fedge;
		MeshElement *base_element;
		bool low_ratio;
		double qmorph_min_ratio;
		double qmorph_max_ratio;
		double qmorph_max_ratio2;
		MeshPoint2d* mp[4];
		int element_id;
		FrontEdge* side_fedges[3];
		bool side_fedges_available[3];
		MeshPoint2d* side_points[3];
		MeshTriangle2d* side_triangles[2];
		int merge_method;
	};
public:
	/// Improve mesh before front
	static void swapTrianglesAhead(Metric2dContext& mc, QuadConversionContext * qc);
	/// Prepare front edges after forming new quad
	static bool prepareQuadFrontEdges(FrontLines *front, QuadConversionContext * qc);
	/// Update front edges after forming new quad
	static bool updateQuadFrontEdges(Metric2dContext& mc, FrontLines *front, 
		QuadConversionContext * qc);
	/// Selects triangle from the given side for merging
	static MeshTriangle2d* selectSideTriangle(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines *front, 
		QuadConversionContext * qc, int side);
	/// Selects best triangle from the two available sides for merging
	static MeshTriangle2d* selectBestSideTriangle(Metric2dContext& mc,  
		FrontLines *front, QuadConversionContext * qc);
	/// Check for single triangle within front
	static bool checkForSingleFrontTriangle(Metric2dContext& mc, FrontLines *front, 
		QuadConversionContext * qc);
	/// Tries to create single quad for selected front edge
	static bool createLeeloQuad(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines * front, 
		QuadConversionContext * qc);
	/// Split triangles intro quads
	static int splitToQuads(Metric2dContext& mc, MeshContainer2d* mesh);
	/// Improves quadrilaterals with different topological methods
	static void improveQuads(Metric2dContext& mc, MeshContainer2d* mesh, int steps, int method);
	/// Improves quadrilaterals with different topological methods
	static void improveQuadsTopological(Metric2dContext& mc, MeshContainer2d* mesh, int steps = 3);
	/// Improves badly shaped quadrilaterals adjacent to boundary
	static int improveQuadsAtBoundary(Metric2dContext& mc, MeshContainer2d* mesh);
	/// Improves topologically the quadrilateral mesh by removing obsolete nodes of rank 2
	static bool improveQuadsByNodeElimination(MeshContainer2d* mesh);
	/// Improves topologically the quadrilateral mesh by removing poor elements in specific configuration
	static bool improveQuadsByElementElimination(Metric2dContext& mc, MeshContainer2d* mesh);
	/// Improves topologically the quadrilateral mesh by swapPIng of edges (to equalize the ranks of nodes)
	static bool improveQuadsByDiagonalSwapPIng(Metric2dContext& mc, MeshContainer2d* mesh);
	/// Improves topologically the quadrilateral mesh by removing obsolete edges in specific configuration
	static bool improveQuadsByEdgeElimination(MeshContainer2d* mesh);
	/// Call several quad-smoothing procedures several times
	static bool smoothenFaces(MeshContainer3d* boundary, int steps = 1, int method = MeshData::SM_LAPLACE);
	/// Converts all triangular meshes in the boundary into the quardrilateral one, using the selected method
	static int convertFacesToQuads(MeshContainer3d* boundary, int method, int max_ct = -1);
	/// Joins two running-front of triangle-to-quad conversion (if the parity of fronts can be maintained)
	static bool joinFrontConditionally(Metric2dContext& mc, FrontLines* front, MeshPoint2d* point1, MeshPoint2d* point2);
	/// Creates the cycle of front edges starting with the i-th edge of the given element
	static int gatherFrontEdges(FrontLines * front, const MeshElement* element, int index);
	/// Analyses the boundary edges and joins them into "running front"
	static int gatherFrontEdges(MeshContainer2d* mesh, FrontLines* front);
	/// Converts triangular mesh into the quardrilateral one, using the mixed QMorph/LeeLo method
	static int convertToQuadsMixed(Metric2dContext& mc, MeshContainer2d* mesh, 
		int merge_method = MeshData::QUADS_MIXED, int max_ct = -1);
	/// Improves node location using selected method
	static bool smoothenNode(Metric2dContext& mc, MeshPoint2d* point);
public:
	/// Metric gradation threshold for mixed quad conversion
	static double param_metric_gradation_threshold;
	/// Metric edge length (squared) threshold for mixed quad conversion
	static double param_metric_length_threshold;
	/// Minimum front edge level for mixed quad conversion, where LeeLo method is allowed
	static int param_leelo_minimum_level;
	/// Which method of smoothing (tyPIcally Laplace or Laplace-metric) should be used for relocating nodes during qmorph conversion
	static int param_convert_smoothing_method;
};

#endif // !defined(MESHGENERATOR2DQUAD_H__INCLUDED)
