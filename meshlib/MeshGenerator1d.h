/////////////////////////////////////////////////////////////////////////////
// MeshGenerator1d.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(MESHGENERATOR1D_H__INCLUDED)
#define MESHGENERATOR1D_H__INCLUDED
#include "MeshData.h"

class MeshContainer3d;
class MeshDomainEdge3d;

/**
 * This class gathers several procedures responsible for discretizing curves into polylines.
 */
class MeshGenerator1d  
{
public:
	/// Discretizes boundary edges of domain using control spaces of all inciden faces (minimum length)
	static int discretizeEdgesMin(MeshContainer3d* boundary);
	/// Discretizes boundary edge of domain using control spaces of all inciden faces (minimum length)
	static int discretizeEdgeMin(MeshDomainEdge3d* edge3d);
public:
	/// Whether to enforce even number of nodes for each disretized boundary edge
	static int param_even_node_count;
	/// Whether to smoothen discretization of edge (for last segment, which may be smaller than required)
	static int param_smoothen_last_node;
	/// Maximum difference for metric-length at edge vertices during discretization
	static double param_boundary_metric_conform_rato;
};

#endif // !defined(MESHGENERATOR1D_H__INCLUDED)
