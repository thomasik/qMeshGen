/////////////////////////////////////////////////////////////////////////////
// DataOctree.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2016-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef DATAOCTREE_H__INCLUDED
#define DATAOCTREE_H__INCLUDED

#include "common.h"
#include <iostream>
#include <functional>

#include "DPoint.h"
#include "DRect.h"
#include "DataVector.h"

/**
 * This class implements an octree
 *  (values in leaves only) plus some basic operations.
 */
template<typename T,
		size_t SamplingSize = 4,
		size_t MaxLevel = 7>
class DataOctreeL
{
public:
	DataOctreeL(const DBox & bbox) : m_box(bbox) {}
	~DataOctreeL() { clear(); }
private:
	DataOctreeL(const DataOctreeL& tree) = delete;
	DataOctreeL& operator=(const DataOctreeL& tree) = delete;
public:
	int adaptToField(const std::function<T(const DPoint3d & pt)> & f, double max_error) {
		T v = f(m_box.getMiddlePoint());
		if (m_top != nullptr) clear();
		m_top = new KdElement(v);
		return 1 + m_top->adaptToField(f, m_box, max_error);
	}
	void clear() {
		if (m_top != nullptr) {
			m_top->clear();
			delete m_top;
			m_top = nullptr;
		}
	}
	const T& getValue(const DPoint3d& pt) const {
		assert(m_top);
		return m_top->findLeaf(pt)->m_data.value;
	}
	inline const DBox& getBBox() const { return m_box; }
public: // stat
	void forRegularGrid(int grid_size, const std::function<void(const DPoint3d& pt)>& fg) {
		double fx = 1.0 / (double)grid_size;
		double dx = m_box.getDX() * fx;
		double dy = m_box.getDY() * fx;
		double dz = m_box.getDZ() * fx;

		DPoint3d pt(m_box.x0 + 0.5*dx, m_box.y0, m_box.z0);
		for (; pt.x < m_box.x1; pt.x += dx) {
			for (pt.y = m_box.y0 + 0.5*dy; pt.y < m_box.y1; pt.y += dy) {
				for (pt.z = m_box.z0 + 0.5*dz; pt.z < m_box.z1; pt.z += dz) {
					fg(pt);
				}
			}
		}
	}
	int getElementCount() const { return m_top ? m_top->getElementCount() : 0; }
	int getValueCount() const { return m_top ? m_top->getValueCount() : 0; }
	int getMaxDepth() const { return m_top ? m_top->getMaxDepth() : 0; }
	int getMaxBalance() const { 
		int max_balance = 0;
		if(m_top) m_top->getMaxDepthAndBalance(max_balance);
		return max_balance;
	}
	int getTotalBytes() const {
		return sizeof(*this) // main struct
			+ sizeof(KdElement) * getElementCount(); // tree nodes
	}
private:
	class KdElement {
	public:
		KdElement(const T & v) : m_data(v) { }
		int adaptToField(const std::function<T(const DPoint3d & pt)> & f, 
			const DBox & box, double max_error, int level = 0) 
		{

			assert(is_leaf);
			if (level >= MaxLevel) return 0;
			double curr_error = estimateError(box, f);
			if (curr_error <= max_error) return 0;
			is_leaf = false;
			int result = VX3D_COUNT;
			m_data.split.middle = box.getMiddlePoint();
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				DBox vbox = box.splitOct(m_data.split.middle, vi);
				m_data.split.elements[vi] = new KdElement(f(vbox.getMiddlePoint()));
				result += m_data.split.elements[vi]->adaptToField(f,
					vbox, max_error, level + 1);
			}
			return result;
		}
		double estimateError(const DBox &box, const std::function<T(const DPoint3d & pt)> &f)
		{
			double error = 0.0;
			double fx = 1.0 / (double)SamplingSize;
			double dx = box.getDX() * fx;
			double dy = box.getDY() * fx;
			double dz = box.getDZ() * fx;

			DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
			for (; pt.x < box.x1; pt.x += dx) {
				for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
					for (pt.z = box.z0 + 0.5*dz; pt.z < box.z1; pt.z += dz) {
						error = std::max(diffKdValue(f(pt), m_data.value),
							error);
					}
				}
			}

			return error;
		}
		void clear() {
			if (!is_leaf) {
				for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
					m_data.split.elements[vi]->clear();
					delete m_data.split.elements[vi];
				}
				is_leaf = true;
			}
		}
	public:
		KdElement* findLeaf(const DPoint3d& pt) {
			KdElement* kd_leaf = this;
			while (!kd_leaf->is_leaf) {
				if (pt.x < kd_leaf->m_data.split.middle.x)
					if (pt.y < kd_leaf->m_data.split.middle.y)
						if (pt.z < kd_leaf->m_data.split.middle.z)
							kd_leaf = kd_leaf->m_data.split.elements[VX3D_LSW];
						else
							kd_leaf = kd_leaf->m_data.split.elements[VX3D_HSW];
					else
						if (pt.z < kd_leaf->m_data.split.middle.z)
							kd_leaf = kd_leaf->m_data.split.elements[VX3D_LNW];
						else
							kd_leaf = kd_leaf->m_data.split.elements[VX3D_HNW];
				else if (pt.y < kd_leaf->m_data.split.middle.y)
					if (pt.z < kd_leaf->m_data.split.middle.z)
						kd_leaf = kd_leaf->m_data.split.elements[VX3D_LSE];
					else
						kd_leaf = kd_leaf->m_data.split.elements[VX3D_HSE];
				else
					if (pt.z < kd_leaf->m_data.split.middle.z)
						kd_leaf = kd_leaf->m_data.split.elements[VX3D_LNE];
					else
						kd_leaf = kd_leaf->m_data.split.elements[VX3D_HNE];
			}
			return kd_leaf;
		}

		int getElementCount() const {
			if (is_leaf) return 1;
			int result = 1;
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				result += m_data.split.elements[vi]->getElementCount();
			}
			return result;
		}
		int getValueCount() const {
			if (is_leaf) return 1;
			int result = 0;
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				result += m_data.split.elements[vi]->getValueCount();
			}
			return result;
		}
		int getMaxDepth() const {
			if (is_leaf) return 1;
			int dmax = 0;
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				int d = m_data.split.elements[vi]->getMaxDepth();
				if (d > dmax) dmax = d;
			}
			return 1+dmax;
		}
		int getMaxDepthAndBalance(int & max_balance) const {
			if (is_leaf) return 1;
			int dmax,dmin;
			dmax = dmin = m_data.split.elements[VX3D_LAST]->getMaxDepthAndBalance(max_balance);
			for (int vi = VX3D_FIRST; vi < VX3D_LAST; vi++) {
				int d = m_data.split.elements[vi]->getMaxDepthAndBalance(max_balance);
				if (d > dmax) dmax = d;
				if (d < dmin) dmin = d;
			}
			int balance = dmax - dmin;
			if (balance > max_balance) max_balance = balance;
			return 1 + dmax;
		}
	public:
		bool is_leaf = true;
		union KdData {
			KdData(const T & v) : value(v) {}
			struct KdSplitElement{
				DPoint3d middle;
				KdElement * elements[8] = { nullptr };
			} split;
			T value;
		} m_data;
	};
private:
	DBox m_box;
	KdElement* m_top = nullptr;
};


/**
* This class implements a balanced octree 
*  (values in leaves only) plus some basic operations.
*/
template<typename T,
	size_t SamplingSize = 4,
	size_t MaxBalance = 1,
	size_t MaxLevel = 7>
class DataOctreeLBalanced
{
public:
	DataOctreeLBalanced(const DBox & bbox) : m_box(bbox) {}
	~DataOctreeLBalanced() { clear(); }
private:
	DataOctreeLBalanced(const DataOctreeLBalanced& tree) = delete;
	DataOctreeLBalanced& operator=(const DataOctreeLBalanced& tree) = delete;
public:
	int adaptToField(const std::function<T(const DPoint3d & pt)> & f, double max_error) {
		T v = f(m_box.getMiddlePoint());
		if (m_top != nullptr) clear();
		m_top = new KdElement(v, 0, m_box, KdElement::KdNeighbors());
		return 1 + m_top->adaptToField(f, max_error);
	}
	void clear() {
		if (m_top != nullptr) {
			m_top->clear();
			delete m_top;
			m_top = nullptr;
		}
	}
	T getValue(const DPoint3d& pt) const {
		assert(m_top);
		return m_top->getNearestLeaf(pt)->m_data.leaf.value;
	}
	inline const DBox& getBBox() const { return m_box; }
public: // stat
	void forRegularGrid(int grid_size, const std::function<void(const DPoint3d& pt)>& fg) {
		double fx = 1.0 / (double)grid_size;
		double dx = m_box.getDX() * fx;
		double dy = m_box.getDY() * fx;
		double dz = m_box.getDZ() * fx;

		DPoint3d pt(m_box.x0 + 0.5*dx, m_box.y0, m_box.z0);
		for (; pt.x < m_box.x1; pt.x += dx) {
			for (pt.y = m_box.y0 + 0.5*dy; pt.y < m_box.y1; pt.y += dy) {
				for (pt.z = m_box.z0 + 0.5*dz; pt.z < m_box.z1; pt.z += dz) {
					fg(pt);
				}
			}
		}
	}
	int getElementCount() const { return m_top ? m_top->getElementCount() : 0; }
	int getValueCount() const { return m_top ? m_top->getValueCount() : 0; }
	int getMaxDepth() const { return m_top ? m_top->getMaxDepth() : 0; }
	int getMaxBalance() const {
		int max_balance = 0;
		if (m_top) m_top->getMaxDepthAndBalance(max_balance);
		return max_balance;
	}
	int getTotalBytes() const {
		return sizeof(*this) // main struct
			+ sizeof(KdElement) * getElementCount(); // tree nodes
	}
private:
	class KdElement {
	public:
		struct KdNeighbors {
		public:
			KdNeighbors(KdElement* _low, KdElement* _high, KdElement* _south, 
				KdElement* _north, KdElement* _east, KdElement* _west) 
				: low(_low), high(_high), south(_south), north(_north), east(_east), west(_west) {}
			KdNeighbors() = default;
		public:
			KdElement* low = nullptr;
			KdElement* high = nullptr;
			KdElement* south = nullptr;
			KdElement* north = nullptr;
			KdElement* east = nullptr;
			KdElement* west = nullptr;
		};
	public:
		KdElement(const T & v, int lvl, const DBox& bx, const KdNeighbors& kn)
			: m_data(v, lvl, bx, kn) { }
		int split(const std::function<T(const DPoint3d & pt)> & f) {
			if (!is_leaf || m_data.leaf.level >= MaxLevel) return 0;

			int split_level = m_data.leaf.level+1;
			DBox split_box = m_data.leaf.box;
			KdNeighbors split_kn = m_data.leaf.kn;

			is_leaf = false;

			int result = VX3D_COUNT;
			m_data.split.middle = split_box.getMiddlePoint();
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				DBox vbox = split_box.splitOct(m_data.split.middle, vi);
				m_data.split.elements[vi] = new KdElement(f(vbox.getMiddlePoint()),
					split_level, vbox, split_kn);
			}
			// update kn for vi
			m_data.split.elements[VX3D_LSW]->setKN( KdNeighbors( 
				split_kn.low, m_data.split.elements[VX3D_HSW],
				split_kn.south, m_data.split.elements[VX3D_LNW],
				m_data.split.elements[VX3D_LSE], split_kn.west));
			m_data.split.elements[VX3D_LSE]->setKN(KdNeighbors(
				split_kn.low, m_data.split.elements[VX3D_HSE],
				split_kn.south, m_data.split.elements[VX3D_LNE],
				split_kn.east, m_data.split.elements[VX3D_LSW]));
			m_data.split.elements[VX3D_LNW]->setKN(KdNeighbors(
				split_kn.low, m_data.split.elements[VX3D_HNW],
				m_data.split.elements[VX3D_LSW], split_kn.north,
				m_data.split.elements[VX3D_LNE], split_kn.west));
			m_data.split.elements[VX3D_LNE]->setKN(KdNeighbors(
				split_kn.low, m_data.split.elements[VX3D_HNE],
				m_data.split.elements[VX3D_LSE], split_kn.north,
				split_kn.east, m_data.split.elements[VX3D_LNW]));
			m_data.split.elements[VX3D_HSW]->setKN(KdNeighbors(
				m_data.split.elements[VX3D_LSW], split_kn.high,
				split_kn.south, m_data.split.elements[VX3D_HNW],
				m_data.split.elements[VX3D_HSE], split_kn.west));
			m_data.split.elements[VX3D_HSE]->setKN(KdNeighbors(
				m_data.split.elements[VX3D_LSE], split_kn.high,
				split_kn.south, m_data.split.elements[VX3D_HNE],
				split_kn.east, m_data.split.elements[VX3D_HSW]));
			m_data.split.elements[VX3D_HNW]->setKN(KdNeighbors(
				m_data.split.elements[VX3D_LNW], split_kn.high,
				m_data.split.elements[VX3D_HSW], split_kn.north,
				m_data.split.elements[VX3D_HNE], split_kn.west));
			m_data.split.elements[VX3D_HNE]->setKN(KdNeighbors(
				m_data.split.elements[VX3D_LNE], split_kn.high,
				m_data.split.elements[VX3D_HSE], split_kn.north,
				split_kn.east, m_data.split.elements[VX3D_HSW]));
			// balance
			KdElement* kd_t[] = { split_kn.east, split_kn.west, split_kn.south, split_kn.north,
				split_kn.low, split_kn.high };
			for (KdElement *kd_e : kd_t) {
				if (kd_e) {
					KdElement* kd_leaf = kd_e->getNearestLeaf(m_data.split.middle);
					if (kd_leaf->m_data.leaf.level < split_level - MaxBalance) 
						result += kd_leaf->split(f);
				}
			}

			return result;
		}
		int adaptToField(const std::function<T(const DPoint3d & pt)> & f, double max_error)
		{
			int result = 0;

			if (is_leaf) {
				if (m_data.leaf.level >= MaxLevel) return 0;

				double curr_error = estimateError(f);
				if (curr_error <= max_error) return 0;

				result += split(f);
			}

			assert(!is_leaf);
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				result += m_data.split.elements[vi]->adaptToField(f, max_error);
			}

			return result;
		}
		double estimateError(const std::function<T(const DPoint3d & pt)> &f)
		{
			assert(is_leaf);
			double error, max_error = 0.0;
			double fx = 1.0 / (double)SamplingSize;
			double dx = m_data.leaf.box.getDX() * fx;
			double dy = m_data.leaf.box.getDY() * fx;
			double dz = m_data.leaf.box.getDZ() * fx;

			DPoint3d pt(m_data.leaf.box.x0 + 0.5*dx, m_data.leaf.box.y0, m_data.leaf.box.z0);
			for (; pt.x < m_data.leaf.box.x1; pt.x += dx) {
				for (pt.y = m_data.leaf.box.y0 + 0.5*dy; pt.y < m_data.leaf.box.y1; pt.y += dy) {
					for (pt.z = m_data.leaf.box.z0 + 0.5*dz; pt.z < m_data.leaf.box.z1; pt.z += dz) {
						error = diffKdValue(f(pt), m_data.leaf.value);
						if (error > max_error) max_error = error;
					}
				}
			}

			return max_error;
		}
		KdElement* getNearestLeaf(const DPoint3d& pt) {
			KdElement* kd_leaf = this;
			while (!kd_leaf->is_leaf) {
				if (pt.x < kd_leaf->m_data.split.middle.x)
					if (pt.y < kd_leaf->m_data.split.middle.y)
						if (pt.z < kd_leaf->m_data.split.middle.z)
							kd_leaf = kd_leaf->m_data.split.elements[VX3D_LSW];
						else
							kd_leaf = kd_leaf->m_data.split.elements[VX3D_HSW];
					else
						if (pt.z < kd_leaf->m_data.split.middle.z)
							kd_leaf = kd_leaf->m_data.split.elements[VX3D_LNW];
						else
							kd_leaf = kd_leaf->m_data.split.elements[VX3D_HNW];
				else if (pt.y < kd_leaf->m_data.split.middle.y)
					if (pt.z < kd_leaf->m_data.split.middle.z)
						kd_leaf = kd_leaf->m_data.split.elements[VX3D_LSE];
					else
						kd_leaf = kd_leaf->m_data.split.elements[VX3D_HSE];
				else
					if (pt.z < kd_leaf->m_data.split.middle.z)
						kd_leaf = kd_leaf->m_data.split.elements[VX3D_LNE];
					else
						kd_leaf = kd_leaf->m_data.split.elements[VX3D_HNE];
			}
			return kd_leaf;
		}
		void clear() {
			if (!is_leaf) {
				for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
					m_data.split.elements[vi]->clear();
					delete m_data.split.elements[vi];
				}
				is_leaf = true;
			}
		}
		void setKN(const KdNeighbors& kn) {
			assert(is_leaf);
			m_data.leaf.kn = kn;
		}
	public:
		int getElementCount() const {
			if (is_leaf) return 1;
			int result = 1;
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				result += m_data.split.elements[vi]->getElementCount();
			}
			return result;
		}
		int getValueCount() const {
			if (is_leaf) return 1;
			int result = 0;
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				result += m_data.split.elements[vi]->getValueCount();
			}
			return result;
		}
		int getMaxDepth() const {
			if (is_leaf) return 1;
			int dmax = 0;
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				int d = m_data.split.elements[vi]->getMaxDepth();
				if (d > dmax) dmax = d;
			}
			return 1 + dmax;
		}
		int getMaxDepthAndBalance(int & max_balance) const {
			if (is_leaf) return 1;
			int dmax, dmin;
			dmax = dmin = m_data.split.elements[VX3D_LAST]->getMaxDepthAndBalance(max_balance);
			for (int vi = VX3D_FIRST; vi < VX3D_LAST; vi++) {
				int d = m_data.split.elements[vi]->getMaxDepthAndBalance(max_balance);
				if (d > dmax) dmax = d;
				if (d < dmin) dmin = d;
			}
			int balance = dmax - dmin;
			if (balance > max_balance) max_balance = balance;
			return 1 + dmax;
		}
	public:
		bool is_leaf = true;
		union KdData {
			KdData(const T & v, int lvl, const DBox& box, const KdNeighbors& kn)
				: leaf(v, lvl, box, kn) {}
			struct KdSplitElement {
				DPoint3d middle;
				KdElement * elements[VX3D_COUNT] = { nullptr };
			} split;
			struct KdLeaf {
				KdLeaf(const T& v, int lvl, const DBox& bx, const KdNeighbors& kn)
					: value(v), level(lvl), box(bx), kn(kn) {}
				T value;
				int level = 0;
				DBox box;
				KdNeighbors kn;
			} leaf;
		} m_data;
	};
private:
	DBox m_box;
	KdElement* m_top = nullptr;
};

/**
* This class implements an octree
*  (values in vertices) plus some basic operations.
*/
template<typename T,
	size_t SamplingSize = 4,
	size_t MaxLevel = 7>
	class DataOctreeV
{
public:
	DataOctreeV(const DBox & bbox) : m_box(bbox) {}
	~DataOctreeV() { clear(); }
private:
	DataOctreeV(const DataOctreeV& tree) = delete;
	DataOctreeV& operator=(const DataOctreeV& tree) = delete;
public:
	int adaptToField(const std::function<T(const DPoint3d & pt)> & f, double max_error) {
		KdVertices kv;

		for(int v = VX3D_FIRST; v <= VX3D_LAST; v++)
			kv.vertices[v] = &m_values.get(m_values.add(f(m_box.getVertex(v))));

		if (m_top != nullptr) clear();
		m_top = new KdElement(kv, 0, KdElement::KdNeighbors() );
		return 1 + m_top->adaptToField(f, m_box, max_error, m_values);
	}
	void clear() {
		if (m_top != nullptr) {
			m_top->clear();
			delete m_top;
			m_top = nullptr;
			m_values.clear();
		}
	}
	T getValue(const DPoint3d& pt) const {
		assert(m_top);
		DBox vbox;
		KdElement* leaf = m_top->getNearestLeaf(pt, m_box, vbox);
		assert(leaf && leaf->is_leaf);
		return leaf->getValue(pt, vbox);
	}
	inline const DBox& getBBox() const { return m_box; }
public: // stat
	void forRegularGrid(int grid_size, const std::function<void(const DPoint3d& pt)>& fg) {
		double fx = 1.0 / (double)grid_size;
		double dx = m_box.getDX() * fx;
		double dy = m_box.getDY() * fx;
		double dz = m_box.getDZ() * fx;

		DPoint3d pt(m_box.x0 + 0.5*dx, m_box.y0, m_box.z0);
		for (; pt.x < m_box.x1; pt.x += dx) {
			for (pt.y = m_box.y0 + 0.5*dy; pt.y < m_box.y1; pt.y += dy) {
				for (pt.z = m_box.z0 + 0.5*dz; pt.z < m_box.z1; pt.z += dz) {
					fg(pt);
				}
			}
		}
	}
	int getElementCount() const { return m_top ? m_top->getElementCount() : 0; }
	int getValueCount() const { return (int)m_values.countInt(); }
	int getMaxDepth() const { return m_top ? m_top->getMaxDepth() : 0; }
	int getMaxBalance() const {
		int max_balance = 0;
		if (m_top) m_top->getMaxDepthAndBalance(max_balance);
		return max_balance;
	}
	int getTotalBytes() const {
		return sizeof(*this) // main struct
			+ sizeof(KdElement) * getElementCount() // tree nodes
			+ sizeof(T) * getValueCount(); // values
	}
private:
	struct KdVertices {
		KdVertices() = default;
		T*& operator[](int i) { assert(i < VX3D_COUNT); return vertices[i]; }
		T* const & operator[](int i) const { assert(i < VX3D_COUNT); return vertices[i]; }
		T* vertices[VX3D_COUNT] = { nullptr };
	};
	class KdElement {
	public:
		struct KdNeighbors {
		public:
			KdNeighbors(KdElement* _low, KdElement* _high, KdElement* _south,
				KdElement* _north, KdElement* _east, KdElement* _west)
			{
				elements[FC3D_LOW] = _low;		elements[FC3D_HIGH] = _high;
				elements[FC3D_SOUTH] = _south;	elements[FC3D_NORTH] = _north;
				elements[FC3D_EAST] = _east;	elements[FC3D_WEST] = _west;
			}
			KdNeighbors() = default;
			KdElement*& operator[](int i) { assert(i < FC3D_COUNT); return elements[i]; }
		public:
			KdElement* elements[FC3D_COUNT] = { nullptr };
		};
	public:
		KdElement(const KdVertices& kv, int lvl, const KdNeighbors& kn)
			: m_data(kv, lvl, kn) { }
		int split(const std::function<T(const DPoint3d & pt)> & f, const DBox& box, DataVector<T>& kvalues) {
			if (!is_leaf || m_data.leaf.level >= MaxLevel) return 0;

			int split_level = m_data.leaf.level + 1;
			KdVertices split_kv = m_data.leaf.kv;
			KdNeighbors split_kn = m_data.leaf.kn;

			is_leaf = false;

			int result = VX3D_COUNT;
			m_data.split.middle = box.getMiddlePoint();
			// -> vertices
			T* mid_box_vertex = &kvalues.get(kvalues.add(f(m_data.split.middle)));
			T* mid_face_vertices[FC3D_COUNT] = { nullptr };
			T* mid_edge_vertices[ED3D_COUNT] = { nullptr };
			for (int fi = FC3D_FIRST; fi <= FC3D_LAST; fi++) {
				DPoint3d face_mid_pt = box.getMiddlePointForFace(fi);
				KdElement* f_kn = split_kn[fi];
				if (f_kn && !f_kn->is_leaf) {
					mid_face_vertices[fi] = f_kn->getNearestVertexRef(face_mid_pt, box);
				} else {
					mid_face_vertices[fi] = &kvalues.get(kvalues.add(f(face_mid_pt)));
				}
			}

			for (int ei = ED3D_FIRST; ei <= ED3D_LAST; ei++) {
				DPoint3d edge_mid_pt = box.getMiddlePointForEdge(ei);
				OctFaceWhich fi = DBox::edge_to_face[ei][0];
				KdElement* f_kn = split_kn[fi];
				bool available = (f_kn && !f_kn->is_leaf);
				if (!available) {
					fi = DBox::edge_to_face[ei][1];
					f_kn = split_kn[fi];
					available = (f_kn && !f_kn->is_leaf);
				}
				if (available) {
					mid_edge_vertices[ei] = f_kn->getNearestVertexRef(edge_mid_pt, box);
				} else {
					mid_edge_vertices[ei] = &kvalues.get(kvalues.add(f(edge_mid_pt)));
				}
			}

			KdVertices sub_kv[VX3D_COUNT];
			sub_kv[VX3D_LSW][VX3D_LSW] = split_kv[VX3D_LSW];
			sub_kv[VX3D_LSW][VX3D_LSE] = mid_edge_vertices[ED3D_LS];
			sub_kv[VX3D_LSW][VX3D_LNW] = mid_edge_vertices[ED3D_LW];
			sub_kv[VX3D_LSW][VX3D_LNE] = mid_face_vertices[FC3D_LOW];
			sub_kv[VX3D_LSW][VX3D_HSW] = mid_edge_vertices[ED3D_SW];
			sub_kv[VX3D_LSW][VX3D_HSE] = mid_face_vertices[FC3D_SOUTH];
			sub_kv[VX3D_LSW][VX3D_HNW] = mid_face_vertices[FC3D_WEST];
			sub_kv[VX3D_LSW][VX3D_HNE] = mid_box_vertex;

			sub_kv[VX3D_LSE][VX3D_LSW] = mid_edge_vertices[ED3D_LS];
			sub_kv[VX3D_LSE][VX3D_LSE] = split_kv[VX3D_LSE];
			sub_kv[VX3D_LSE][VX3D_LNW] = mid_face_vertices[FC3D_LOW];
			sub_kv[VX3D_LSE][VX3D_LNE] = mid_edge_vertices[ED3D_LE];
			sub_kv[VX3D_LSE][VX3D_HSW] = mid_face_vertices[FC3D_SOUTH];
			sub_kv[VX3D_LSE][VX3D_HSE] = mid_edge_vertices[ED3D_SE];
			sub_kv[VX3D_LSE][VX3D_HNW] = mid_box_vertex;
			sub_kv[VX3D_LSE][VX3D_HNE] = mid_face_vertices[FC3D_EAST];

			sub_kv[VX3D_LNW][VX3D_LSW] = mid_edge_vertices[ED3D_LW];
			sub_kv[VX3D_LNW][VX3D_LSE] = mid_face_vertices[FC3D_LOW];
			sub_kv[VX3D_LNW][VX3D_LNW] = split_kv[VX3D_LNW];
			sub_kv[VX3D_LNW][VX3D_LNE] = mid_edge_vertices[ED3D_LN];
			sub_kv[VX3D_LNW][VX3D_HSW] = mid_face_vertices[FC3D_WEST];
			sub_kv[VX3D_LNW][VX3D_HSE] = mid_box_vertex;
			sub_kv[VX3D_LNW][VX3D_HNW] = mid_edge_vertices[ED3D_NW];
			sub_kv[VX3D_LNW][VX3D_HNE] = mid_face_vertices[FC3D_NORTH];

			sub_kv[VX3D_LNE][VX3D_LSW] = mid_face_vertices[FC3D_LOW];
			sub_kv[VX3D_LNE][VX3D_LSE] = mid_edge_vertices[ED3D_LE];
			sub_kv[VX3D_LNE][VX3D_LNW] = mid_edge_vertices[ED3D_LN];
			sub_kv[VX3D_LNE][VX3D_LNE] = split_kv[VX3D_LNE];
			sub_kv[VX3D_LNE][VX3D_HSW] = mid_box_vertex;
			sub_kv[VX3D_LNE][VX3D_HSE] = mid_face_vertices[FC3D_EAST];
			sub_kv[VX3D_LNE][VX3D_HNW] = mid_face_vertices[FC3D_NORTH];
			sub_kv[VX3D_LNE][VX3D_HNE] = mid_edge_vertices[ED3D_NE];

			sub_kv[VX3D_HSW][VX3D_LSW] = mid_edge_vertices[ED3D_SW];
			sub_kv[VX3D_HSW][VX3D_LSE] = mid_face_vertices[FC3D_SOUTH];
			sub_kv[VX3D_HSW][VX3D_LNW] = mid_face_vertices[FC3D_WEST];
			sub_kv[VX3D_HSW][VX3D_LNE] = mid_box_vertex;
			sub_kv[VX3D_HSW][VX3D_HSW] = split_kv[VX3D_HSW];
			sub_kv[VX3D_HSW][VX3D_HSE] = mid_edge_vertices[ED3D_HS];
			sub_kv[VX3D_HSW][VX3D_HNW] = mid_edge_vertices[ED3D_HW];
			sub_kv[VX3D_HSW][VX3D_HNE] = mid_face_vertices[FC3D_HIGH];

			sub_kv[VX3D_HSE][VX3D_LSW] = mid_face_vertices[FC3D_SOUTH];
			sub_kv[VX3D_HSE][VX3D_LSE] = mid_edge_vertices[ED3D_SE];
			sub_kv[VX3D_HSE][VX3D_LNW] = mid_box_vertex;
			sub_kv[VX3D_HSE][VX3D_LNE] = mid_face_vertices[FC3D_EAST];
			sub_kv[VX3D_HSE][VX3D_HSW] = mid_edge_vertices[ED3D_HS];
			sub_kv[VX3D_HSE][VX3D_HSE] = split_kv[VX3D_HSE];
			sub_kv[VX3D_HSE][VX3D_HNW] = mid_face_vertices[FC3D_HIGH];
			sub_kv[VX3D_HSE][VX3D_HNE] = mid_edge_vertices[ED3D_HE];

			sub_kv[VX3D_HNW][VX3D_LSW] = mid_face_vertices[FC3D_WEST];
			sub_kv[VX3D_HNW][VX3D_LSE] = mid_box_vertex;
			sub_kv[VX3D_HNW][VX3D_LNW] = mid_edge_vertices[ED3D_NW];
			sub_kv[VX3D_HNW][VX3D_LNE] = mid_face_vertices[FC3D_NORTH];
			sub_kv[VX3D_HNW][VX3D_HSW] = mid_edge_vertices[ED3D_HW];
			sub_kv[VX3D_HNW][VX3D_HSE] = mid_face_vertices[FC3D_HIGH];
			sub_kv[VX3D_HNW][VX3D_HNW] = split_kv[VX3D_HNW];
			sub_kv[VX3D_HNW][VX3D_HNE] = mid_edge_vertices[ED3D_HN];

			sub_kv[VX3D_HNE][VX3D_LSW] = mid_box_vertex;
			sub_kv[VX3D_HNE][VX3D_LSE] = mid_face_vertices[FC3D_EAST];
			sub_kv[VX3D_HNE][VX3D_LNW] = mid_face_vertices[FC3D_NORTH];
			sub_kv[VX3D_HNE][VX3D_LNE] = mid_edge_vertices[ED3D_NE];
			sub_kv[VX3D_HNE][VX3D_HSW] = mid_face_vertices[FC3D_HIGH];
			sub_kv[VX3D_HNE][VX3D_HSE] = mid_edge_vertices[ED3D_HE];
			sub_kv[VX3D_HNE][VX3D_HNW] = mid_edge_vertices[ED3D_HN];
			sub_kv[VX3D_HNE][VX3D_HNE] = split_kv[VX3D_HNE];

			KdNeighbors sub_kn[VX3D_COUNT]; // empty for now
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				DBox vbox = box.splitOct(m_data.split.middle, vi);
				m_data.split.elements[vi] = new KdElement(sub_kv[vi],
					split_level, sub_kn[vi]);
			}
			// update neighbors (outer)
			for (int fi = FC3D_FIRST; fi <= FC3D_LAST; fi++) {
				KdElement* nb = split_kn[fi];
				if (nb != nullptr && !nb->is_leaf) {
					int ofi = fi ^ 1; // opposite face (xor 1)
					for (int i = 0; i < 4; i++) {
						int vfi = DBox::face_to_vertex[fi][i];
						int ovfi = DBox::face_to_vertex[ofi][i];
						KdElement* sub_okd = nb->m_data.split.elements[ovfi];
						KdElement* sub_kd  =     m_data.split.elements[vfi];
						sub_kd->setKN(fi, sub_okd);
						if (sub_okd->is_leaf) sub_okd->setKN(ofi, sub_kd);
					}
				}
				if (fi % 2 == 0) {
					// inner neighborhood
					int ofi = fi^1; // opposite face (xor 1), also fi+1
					for (int i = 0; i < 4; i++) {
						int vfi = DBox::face_to_vertex[fi][i];
						int ovfi = DBox::face_to_vertex[ofi][i];
						KdElement* sub_okd = m_data.split.elements[ovfi];
						KdElement* sub_kd  = m_data.split.elements[vfi];
						sub_kd->setKN(ofi, sub_okd);
						sub_okd->setKN(fi, sub_kd);
					}
				}
			}

			return result;
		}
		int adaptToField(const std::function<T(const DPoint3d & pt)> & f, 
			const DBox& box, double max_error, DataVector<T>& kvalues)
		{
			int result = 0;

			if (is_leaf) {
				if (m_data.leaf.level >= MaxLevel) return 0;

				double curr_error = estimateError(box, f);
				if (curr_error <= max_error) return 0;

				result += split(f, box, kvalues);
			}

			assert(!is_leaf);
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				DBox vbox = box.splitOct(m_data.split.middle, vi);
				result += m_data.split.elements[vi]->adaptToField(f, vbox, max_error, kvalues);
			}

			return result;
		}
		double estimateError(const DBox& box, const std::function<T(const DPoint3d & pt)> &f)
		{
			assert(is_leaf);
			double error, max_error = 0.0;
			double fx = 1.0 / (double)SamplingSize;
			double dx = box.getDX() * fx;
			double dy = box.getDY() * fx;
			double dz = box.getDZ() * fx;

			DPoint3d pt(box.x0 + 0.5*dx, box.y0, box.z0);
			for (; pt.x < box.x1; pt.x += dx) {
				for (pt.y = box.y0 + 0.5*dy; pt.y < box.y1; pt.y += dy) {
					for (pt.z = box.z0 + 0.5*dz; pt.z < box.z1; pt.z += dz) {
						error = diffKdValue(f(pt), getValue(pt, box));
						if (error > max_error) max_error = error;
					}
				}
			}

			return max_error;
		}
		T getValue(const DPoint3d& pt, const DBox& box) const
		{
			assert(box.contains(pt));

			double tx = (pt.x - box.x0) / box.getDX();
			double ty = (pt.y - box.y0) / box.getDY();
			double tz = (pt.z - box.z0) / box.getDZ();

			double coeff[8] = {
				(1 - tx) * (1 - ty) * (1 - tz),	// VX_LSW
				tx	 * (1 - ty) * (1 - tz),	// VX_LSE
				(1 - tx) *    ty	* (1 - tz),	// VX_LNW
				tx	 *    ty	* (1 - tz),	// VX_LNE
				(1 - tx) * (1 - ty) *    tz ,	// VX_HSW
				tx	 * (1 - ty) *    tz ,	// VX_HSE
				(1 - tx) *    ty	*    tz ,	// VX_HNW
				tx	 *    ty	*    tz 	// VX_HNE
			};
			T value = *(m_data.leaf.kv[0]) * coeff[0];
			for(int vi = 1; vi <= VX3D_LAST; vi++)
				value +=*(m_data.leaf.kv[vi]) * coeff[vi];
			return value;
		}
		KdElement* getNearestLeaf(const DPoint3d& pt, const DBox& box, DBox& vbox) {
			KdElement* kd_leaf = this;
			vbox = box;
			while (!kd_leaf->is_leaf) {
				int vi = kd_leaf->m_data.split.middle.vertexDirection(pt);
				vbox = vbox.splitOct(kd_leaf->m_data.split.middle, vi);
				kd_leaf = kd_leaf->m_data.split.elements[vi];
			}
			return kd_leaf;
		}
		T* getNearestVertexRef(const DPoint3d& pt, const DBox& box) {
			DBox vbox;
			KdElement* kd_leaf = getNearestLeaf(pt, box, vbox);
			int vi = vbox.getMiddlePoint().vertexDirection(pt);
			return kd_leaf->m_data.leaf.kv[vi];
		}
		void clear() {
			if (!is_leaf) {
				for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
					m_data.split.elements[vi]->clear();
					delete m_data.split.elements[vi];
				}
				is_leaf = true;
			}
		}
		void setKN(const KdNeighbors& kn) {
			assert(is_leaf);
			m_data.leaf.kn = kn;
		}
		void setKN(int fi, KdElement* kde) {
			assert(is_leaf);
			m_data.leaf.kn[fi] = kde;
		}
	public:
		int getElementCount() const {
			if (is_leaf) return 1;
			int result = 1;
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				result += m_data.split.elements[vi]->getElementCount();
			}
			return result;
		}
		int getMaxDepth() const {
			if (is_leaf) return 1;
			int dmax = 0;
			for (int vi = VX3D_FIRST; vi <= VX3D_LAST; vi++) {
				int d = m_data.split.elements[vi]->getMaxDepth();
				if (d > dmax) dmax = d;
			}
			return 1 + dmax;
		}
		int getMaxDepthAndBalance(int & max_balance) const {
			if (is_leaf) return 1;
			int dmax, dmin;
			dmax = dmin = m_data.split.elements[VX3D_LAST]->getMaxDepthAndBalance(max_balance);
			for (int vi = VX3D_FIRST; vi < VX3D_LAST; vi++) {
				int d = m_data.split.elements[vi]->getMaxDepthAndBalance(max_balance);
				if (d > dmax) dmax = d;
				if (d < dmin) dmin = d;
			}
			int balance = dmax - dmin;
			if (balance > max_balance) max_balance = balance;
			return 1 + dmax;
		}
	public:
		bool is_leaf = true;
		union KdData {
			KdData(const KdVertices& kv, int lvl, const KdNeighbors& kn)
				: leaf(kv, lvl, kn) {}
			struct KdSplitElement {
				DPoint3d middle;
				KdElement * elements[VX3D_COUNT] = { nullptr };
			} split;
			struct KdLeaf {
				KdLeaf(const KdVertices& kv, int lvl, const KdNeighbors& kn)
					: kv(kv), level(lvl), kn(kn) {}
				KdVertices kv;
				int level = 0;
				KdNeighbors kn;
			} leaf;
		} m_data;
	};
private:
	DBox m_box;
	KdElement* m_top = nullptr;
	DataVector<T> m_values;
};

#endif // !defined(DATAOCTREE_H__INCLUDED)
