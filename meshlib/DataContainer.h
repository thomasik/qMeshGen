/////////////////////////////////////////////////////////////////////////////
// DataContainer.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(DATACONTAINER_H__INCLUDED)
#define DATACONTAINER_H__INCLUDED

#pragma once

#include "common.h"

/**
 * This template class implements a container of pointers to elements
 *	with coherent indexes.
 * If required can be ordered via "heap - structure"
 * Parameter class "DataItem" must implement methods:
 *	- void setIndex(int i)
 *  - int getIndex()
 *  - int compareTo(DataItem* item)  <- only if heap-mode is used
 */
template<class DataItem>
class DataContainer  
{
public:
	/// Standard constructor (container is an array of arrays[part_size] of pointers)
	DataContainer(int part_size, bool heap_order = false){
		m_part_size = part_size;
		m_data_count = m_array_size = m_array_count = 0;
		m_data = nullptr;
		m_heap_order = heap_order;
	}
	/// Standard destructor (destroys all elements)
	~DataContainer(){ deleteAll(); }
	/// copying - deleted
	DataContainer(const DataContainer& dc) = delete;
	/// copying - deleted
	DataContainer& operator=(const DataContainer& dc) = delete;
public:
	///	Positions the given element according to the heap structure
	void updateDataItemPosition(DataItem* item){
		assert(item != nullptr);
		downheap(item->getIndex()); 
		upheap(item->getIndex()); 
	}
	/// Switches positions of elements at indexes i and j
	void switchDataItems(int i, int j){
		assert(i >= 0 && i < m_data_count);
		assert(j >= 0 && j < m_data_count);
		DataItem* temp_i = m_data[i / m_part_size].data[i % m_part_size];
		DataItem* temp_j = m_data[j / m_part_size].data[j % m_part_size];
		m_data[i / m_part_size].data[i % m_part_size] = temp_j;
		temp_j->setIndex(i);
		m_data[j / m_part_size].data[j % m_part_size] = temp_i;
		temp_i->setIndex(j);
	}
	/// Turns on and off the preservation of the heap structure (with smallest element at the top)
	void setHeapOrder(bool valid = true){
		m_heap_order = valid;
		if(valid){
			// Construction of the heap
			for(int i = m_data_count/2 - 1; i>=0; i--) downheap(i);
		}
	}
	/// Returns whether the container is in "heap mode"
	bool isHeapOrder() const { return m_heap_order; }
	/// Removes the given element from the container
	bool removeDataItemValue(DataItem* item) {
		for(int i = 0; i < m_data_count; i++){
			if(getDataAt(i) == item) {
				removeDataItem(i);
				return true;
			}
		}
		return false;
	}
	/// Returns the element removed (reference only) from container (in heap mode there is additional reorganization)
	DataItem* removeDataItem(int index){
		assert(index >= 0 && index < m_data_count);
		int j = index / m_part_size;
		int i = index % m_part_size;
		DataItem *item = m_data[j].data[i];
		if(index == --m_data_count){	// if last -> no additional operations
			m_data[j].count--;
		}else{
			// else (i.e. the removed element is NOT last)
			//	insert here the last element and decrease the element_count
			int k = m_array_count;
			while(m_data[--k].count == 0);
			m_data[j].data[i] = m_data[k].data[--(m_data[k].count)];
			m_data[j].data[i]->setIndex(index);
			// If in heap order (and not leaf) -> additional reorganization is required
			if(m_heap_order){
				//LOG4CPLUS_INFO(MeshLog::logger_console, "data_count", m_data_count);
				//LOG4CPLUS_INFO(MeshLog::logger_console, "remove, fixing", index);
				DataItem* moved_item = m_data[j].data[i];
				if(index < m_data_count/2) downheap(index);
				index = moved_item->getIndex();
				if(index > 0) upheap(index);
			}
		}
		return item;
	}
	/// Adds new element into the container (in heap mode with additional reorganization)
	int addDataItem(DataItem* item, bool no_heap = false){
		assert(item != nullptr);
		int i;
		for(i = 0; (i < m_array_count) && (m_data[i].count >= m_part_size); i++);
		if(i >= m_array_count){
			// New subarray
			i = m_array_count++;
			if(m_array_count >= m_array_size){
				// The main array is too small -> needs to be enlarged 
				//	all subarrays are transfered to the new one
				m_array_size += m_array_size/2;
				if(m_array_size < m_part_size) m_array_size = m_part_size;
				MData* temp = new MData[m_array_size];
				for(int j = 0; j < m_array_count-1; j++) temp[j] = m_data[j];
				delete[] m_data;
				m_data = temp;
			}
			m_data[i].count = 0;
			m_data[i].data = new DataItem * [m_part_size];
		}
		m_data[i].data[m_data[i].count++] = item;
		m_data_count++;

		i = m_part_size * i + m_data[i].count - 1;
		item->setIndex(i);

		// If in heap mode -> reorganization
		if(m_heap_order && !no_heap) upheap(i);
		return i;
	}
	/// Returns element at the given index (doesn't check for "out of bounds" in release mode)
	DataItem* getDataAt(int i) const { 
		assert(i >= 0 && i < m_data_count); 
		return m_data[i / m_part_size].data[i % m_part_size]; 
	}
	DataItem* top() const {
		assert(!empty());
		return m_data[0].m_data[0];
	}
	inline DataItem* pop() { return removeDataItem(0); }
	DataItem* operator[](int i) const {
		assert(i >= 0 && i < m_data_count);
		return m_data[i / m_part_size].data[i % m_part_size];
	}
	DataItem*& operator[](int i) {
		assert(i >= 0 && i < m_data_count);
		return m_data[i / m_part_size].data[i % m_part_size];
	}
	/// Removes all elements (references AND objects) from container
	void deleteAll(){
		for(int i = 0; i < m_array_count; i++){
			for(int j = 0; j < m_data[i].count; j++){
				m_data[i].data[j]->preDeleteAll();
				delete m_data[i].data[j];
			}
			delete[] m_data[i].data;
		}
		if(m_data) delete[] m_data;
		m_data = nullptr;
		m_array_count = m_array_size = m_data_count = 0;
	}
	/// Returns the number of elements in the container
	inline int countInt() const { return m_data_count; }
	inline bool empty() const { return m_data_count == 0; }
	/// Checks heap order
	bool validHeapOrder() const {
		for(int i = 0; ;i++){
			DataItem* father = getDataAt(i);
			// first son
			int is = 2*i+1;
			if(is >= m_data_count) return true;
			DataItem* son = getDataAt(is);
			if(father->compareTo(son) > 0){
				return false;
			}
			// second son
			if(++is >= m_data_count) return true;
			son = getDataAt(is);
			if(father->compareTo(son) > 0){
				return false;
			}
		}
	}
protected:	
	/// Restores the heap-condition for the k-th element (downward)
	void downheap(int k){
		assert(k >= 0 && k < m_data_count);
		// v - element to order
//		LOG4CPLUS_INFO(MeshLog::logger_console, "downheap k (first)", k);
		DataItem* v = getDataAt(k);
		// StepPIng down the heap
		int ndiv2 = m_data_count / 2;
		while( k < ndiv2 ){
			int j = 2*k+1;	// Next element of k-th in the heap
			DataItem* e1 = getDataAt(j);
			if(j+1 < m_data_count){
				// Choose worse element as the next one
				DataItem* e2 = getDataAt(j+1);
				if(e1->compareTo(e2) > 0){
					++j;
					e1 = e2;
				}
			}
			// if( v >= a[j]) break;
			if( v->compareTo(e1) < 1 ) break;
			m_data[k / m_part_size].data[k % m_part_size] = e1;
//			LOG4CPLUS_INFO(MeshLog::logger_console, "downheap j", j);
			e1->setIndex(k);
			k = j;
		}
//		LOG4CPLUS_INFO(MeshLog::logger_console, "downheap k (last)", k);
		m_data[k / m_part_size].data[k % m_part_size] = v;
		v->setIndex(k);
	}
	/// Restores the heap-condition for the k-th element (upward)
	void upheap(int k){
		// v - element do order
		DataItem* v = getDataAt(k);
		int l = (k-1) / 2;	// Previous element of k-th in the heap
		do{
			DataItem *e = getDataAt(l);
			if(v->compareTo(e) > -1) break;
			m_data[k / m_part_size].data[k % m_part_size] = e;
			e->setIndex(k);
			k = l;
			l = (l-1) / 2;
		}while(l != k);
		m_data[k / m_part_size].data[k % m_part_size] = v;
		v->setIndex(k);
	}
protected:
	/// Structure for subarrays
	struct MData{
		/// Number of elements in the subarray
		int count;
		/// Subarray of pointers to elements
		DataItem**	data;
	};
protected:
	/// Factor of size for the main array and size of subarrays
	int	m_part_size;
	/// Main array (of subarrays)
	MData *m_data;
	/// Total number of elements
	int	m_data_count;
	/// Current number of created subarrays
	int	m_array_count;
	/// Current size of the main array (>= m_array_count)
	int	m_array_size;
	/// Mode of ordering of the elements in the container
	bool m_heap_order;
};

#endif // !defined(DATACONTAINER_H__INCLUDED)
