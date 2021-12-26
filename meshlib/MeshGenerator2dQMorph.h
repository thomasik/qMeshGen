/////////////////////////////////////////////////////////////////////////////
// MeshGenerator2dQMorph.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR2DQMORPH_H__INCLUDED)
#define MESHGENERATOR2DQMORPH_H__INCLUDED

#include "MeshData.h"
#include "MeshGenerator2dQuad.h"

class MeshContainer2d;
class FrontLines;
class FrontEdge;
class MeshPoint2d;
class MeshEdge2d;

/**
 * This class gathers several procedures required by the variation of 
 *  the QMorph algorithm of triangle-to-quad conversion.
 */
class MeshGenerator2dQMorph  
{
public:
	/// Tries to create single quad for selected front edge
	static bool createQMorphQuad(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines * front, 
		MeshGenerator2dQuad::QuadConversionContext * qc);
protected:
	struct SeamContext
	{
		SeamContext(Metric2dContext& mc, MeshGenerator2dQuad::QuadConversionContext * qc, int side);
	public:
		FrontEdge *edge1, *edge2;
		MeshEdge2d *mesh_edge1, *mesh_edge2;
		MeshPoint2d *point;
		double angle;
		int level;
		int id;
		double length1, length2;
		double ratio;
	};
	struct InvalidFrontEdge{
		InvalidFrontEdge(FrontEdge* fe = nullptr, MeshPoint2d* _mp1 = nullptr, MeshPoint2d* _mp2 = nullptr) : 
			fedge(fe), mp1(_mp1), mp2(_mp2) {}
	public:
		FrontEdge* fedge;
		MeshPoint2d* mp1;
		MeshPoint2d* mp2;
	};
protected:
	/// Checks the special cases for triangle-to-quad conversion (and fix it)
	static bool checkSeams(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines *front, 
		MeshGenerator2dQuad::QuadConversionContext * qc, int side);
	/// Check seam1 for triangle-to-quad conversion
	static bool trySeam1(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines *front, 
		SeamContext& sc);
	/// Check seam2 for triangle-to-quad conversion
	static bool trySeam2(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines *front, 
		SeamContext& sc);
	/// Check seam3 for triangle-to-quad conversion
	static bool trySeam3(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines *front, 
		SeamContext& sc);
	/// Check if upper front edge can be used for quad forming
	static bool checkQuadUpperEdge(Metric2dContext& mc, 
		MeshGenerator2dQuad::QuadConversionContext * qc, int side);
	/// Form and smoothen new quad
	static bool formNewQuad(Metric2dContext& mc, MeshContainer2d* mesh, 
		MeshGenerator2dQuad::QuadConversionContext * qc);
	/// Set upper edge for quad forming
	static bool setQuadUpperEdge(Metric2dContext& mc, MeshGenerator2dQuad::QuadConversionContext * qc);
	/// Set side edge for quad forming
	static bool setQuadSideEdge(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines *front, 
		MeshGenerator2dQuad::QuadConversionContext * qc, int side);
	/// Check for inner front edges within created quad (true, if found)
	static bool checkForInnerFrontEdges(Metric2dContext& mc, 
		MeshGenerator2dQuad::QuadConversionContext * qc, int side);
protected:
	/// Improves locally the quality (regularity) of triangular/quadrilateral mesh near the "running front"
	static bool smoothenFrontPoint(Metric2dContext& mc, MeshEdge2d* edge, MeshPoint2d* point);
	/// Finds the "front edge" assocciated with the given mesh edge
	static MeshEdge2d* findNeededEdge(Metric2dContext& mc, MeshContainer2d* mesh, FrontLines * front, 
		MeshGenerator2dQuad::QuadConversionContext * qc, MeshEdge2d *edge, MeshPoint2d* point, 
		bool to_left, const DVector2d &dv, double average_length = 1.0);
	/// Removes marked triangles (creating the place for new quadrialteral)
	static int removeInnerTriangles(MeshContainer2d* mesh, MeshTriangle2d* start_triangle,
		const DataVector<MeshEdge2d*> & border_edges);
	/// Improves the local quality (regularity) of triangular mesh ahead of the running front
	static void smoothenAhead(Metric2dContext& mc, MeshContainer2d* mesh, MeshPoint2d* point);
public:
	/// Whether to check size of quads created durign conversion with the control space
	static int param_size_check;
};

#endif // !defined(MESHGENERATOR2DQMORPH_H__INCLUDED)
