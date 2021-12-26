/////////////////////////////////////////////////////////////////////////////
// FrontLines.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(FRONTLINES_H__INCLUDED)
#define FRONTLINES_H__INCLUDED

#include "DataContainer.h"
#include "FrontEdge.h"

/**
 * This class represents a set of edges connected into a closed loop - "running front"
 */
class FrontLines  
{
public:
	/// Standard constructor
	FrontLines(int part_size);
	/// Standard destructor
	~FrontLines();
public:
	/// Postpones this edge to be processed later
	void postpone(FrontEdge* fe);
	/// Checks whether this front-cycles are valid
	bool isValid() const;
	/// Finds and returns the front-edge associated with the given mesh-edge (adequate direction given by index)
	FrontEdge* findFrontEdge(const MeshEdge2d* edge, int index) const;
	/// Removes front-edge from the set of edges (no additional connection-updates)
	void removeFrontEdge(FrontEdge* fedge);
	/// Updates the position of this front-edge within the heap-ordered container (based on the rank of this front-edge)
	void updateFrontEdgePosition(FrontEdge* fedge){ m_front_edges->updateDataItemPosition(fedge); }
	/// Classifies (calculates the ranks) of all front-edges
	void classifyAllEdges(Metric2dContext& mc);
	/// Returns the count of front-edges forming this front
	int countInt() const { return m_front_edges->countInt(); }
	/// Returns the pointer to the currently best front-edge (first in the heap order) or nullptr if there is no edges
	FrontEdge* getBestFrontEdge() const;
	/// Returns the pointer to the front-edge at the given index within the container (heap-ordered)
	FrontEdge* getFrontEdge(int i) const { return m_front_edges->getDataAt(i); }
	/// Inserts new front-edge into the container (no classifing or updating of the container)
	void addFrontEdge(FrontEdge* fedge);
protected:
	/// Container of front-edges (parametrized from the template)
	DataContainer<FrontEdge> *m_front_edges;
};

#endif // !defined(FRONTLINES_H__INCLUDED)
