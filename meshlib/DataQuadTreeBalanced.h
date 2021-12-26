/////////////////////////////////////////////////////////////////////////////
// DataQuadTree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DATAQUADTREEBALANCED_H__INCLUDED)
#define DATAQUADTREEBALANCED_H__INCLUDED

#include "DataVector.h"

/**
 * This class implements a quadtree structure
 */
template<class DataItem, int MAXDEPTH=10>
class DataQuadTreeBalanced
{
public:
	class Node {
	public:
		ControlNode(const DPoint2d& pt = DPoint2d::zero) : coord(pt) {}
		ControlNode(const DPoint2d& pt, const DataItem& _item) : coord(pt), data_item(_item) {}
	public:
		DPoint2d coord;
		DataItem item;
	};

public:
	/// Standard contructor
	DataQuadTreeBalanced(const DRect& box, int nxy = 100);
	/// Destructor
	~DataQuadTreeBalanced();
public:
	/// Adds new data for some point
	void addData(const DPoint2d& pt, const DataItem& item);
	/// Return data at a point
	DataItem getData(const DPoint2d& pt) const;
	/// Interpolates (initializes) the structure basing on the previously given data in discrete points
	bool interpolate();
	/// Refines data at the given point
	bool setMinData(const DPoint2d& pt, const DataItem& item, bool min_value_set = true);
public:
	void testQuadTree();
public:
	typedef DataVector<Node*> QuadVertexList;
	typedef DataMatrix<Node> QuadVertexArray;
	// Class QuadLeaf for QuadTree structure
	class QuadLeaf{
	public:
		QuadLeaf();
		QuadLeaf(double mid_x, double mid_y, double dx, double dy, int level = 0);
		~QuadLeaf();
		/// Get sizing info (matrix mode) at the given point
		DataItem getData(const DPoint2d& pt) const;
		/// Set values of vertices indices
		void setVertices(Node* v_sw, Node* v_se, Node* v_nw, Node* v_ne){
			m_vertices[VX_SW] = v_sw;
			m_vertices[VX_SE] = v_se;
			m_vertices[VX_NW] = v_nw;
			m_vertices[VX_NE] = v_ne;
		}
		/// Set pointers to neighbouring leaves
		void setNeighbours(QuadLeaf* nd, QuadLeaf* nr, QuadLeaf* nu, QuadLeaf* nl){
			m_neighbours[NB_DOWN] = nd;
			m_neighbours[NB_RIGHT] = nr;
			m_neighbours[NB_UP] = nu;
			m_neighbours[NB_LEFT] = nl;
		}
		/// Set pointer to neighbouring leaf for this leaf and sub-leaves
		void setNeighbour(int side, QuadLeaf* nb, bool skip_first = false);
		/// Set coordinates (middle point and size)
		void setCoordinates(const DPoint2d& middle, const DVector2d& dl){
			m_middle = middle;
			m_dl = dl;
		}
		/// Adapt quad leaf to new source point and insert this point into proper collection
		void adaptAndInsertNode(QuadVertexArray& grid_vertices, const Node& qv);
		/// Adapt quad leaf to new source point and set minimum value
		bool adaptAndSetNode(QuadVertexArray& grid_vertices, const Node& qv, min_value_set = true);
		/// Split this leaf (+ build relations with neighbours)
		bool split(QuadVertexArray& grid_vertices, QuadVertexList* split_vertices = NULL, bool init_values = false);
		/// Whether leaf is already split
		bool isSplit() const { return m_leaves[0] != NULL; }
		/// Adjust levels of adjacent leaves
		void balance(QuadVertexArray& grid_vertices, QuadVertexList* split_vertices, bool init_values = false);
		/// Find last leaf (not split) for givent point
		QuadLeaf* findLastLeaf(const DPoint2d& pt);
		/// Gather metric value from source nodes within leaf	
		void gatherDataFromSourcePoints();
		/// Propagate control values from within leaves
		int propagateUp();
		/// Propagate control values into leaves
		int propagateDown();
		/// Returns i-th vertex of the leaf
		Node* getVertex(int i) const { return m_vertices[i]; }
	private:
		/// middle point of leaf (parametric)
		DPoint2d m_middle;
		enum QuadVertexWhich {VX_SW = 0, VX_SE = 1, VX_NW = 2, VX_NE = 3,
				VX_FIRST = 0, VX_LAST = 3};
		ControlNode* m_vertices[4]; // only indices, 2-3 / 0-1
		enum QuadLeafWhich {LF_SW = 0, LF_SE = 1, LF_NW = 2, LF_NE = 3,
				LF_FIRST = 0, LF_LAST = 3};
		QuadLeaf* m_leaves[4];
		//-- extra data
		/// depth level
		int m_level;
		/// dimensions of leaf (parametric)
		DVector2d m_dl;
		enum QuadNeighbourWhich {NB_DOWN = 0, NB_RIGHT = 1, NB_UP = 2, NB_LEFT = 3, 
			NB_FIRST = 0, NB_LAST = 3};
		QuadLeaf* m_neighbours[4]; 
		enum QuadMidVertexWhich {MVX_DOWN = 0, MVX_RIGHT = 1, MVX_UP = 2, MVX_LEFT = 3,
			MVX_FIRST = 0, MVX_LAST = 3};
		static QuadVertexWhich mid_to_nodes[4][2];
		static QuadMidVertexWhich nodes_to_mid[4][2];
		Node* m_mid_vertices[4]; 
		//--- point sources before interpolation
		DataVector<Node> *m_source_points;
	};
protected:
	/// Translates the point-coordinates into index of control matrix
	void countLocalCoordinates(const DPoint2d& pt, int& ix, int& iy) const;
	/// Translates the point-coordinates into (consolidated) index of control matrix
	int countLocalCoordinates(const DPoint2d& pt) const;
	/// Returns smalles QuadTree leaf for this point
	QuadLeaf* findLastLeaf(const DPoint2d& pt) const;
	/// Extrapolate control data from nodes
	int extrapolateMainControlNodes();
protected:
	/// Number of columns in the matrix
	int m_nx;
	/// Number of rows in the matrix
	int m_ny;
	/// Values stored in matrix
	QuadLeaf* m_grid;
	/// Grid vertices with defined control data
	QuadVertexArray m_grid_vertices;
};

#endif // !defined(DATAQUADTREEBALANCED_H__INCLUDED)
