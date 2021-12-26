/////////////////////////////////////////////////////////////////////////////
// DataHeapVector.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(DATAHEAPVECTOR_H__INCLUDED)
#define DATAHEAPVECTOR_H__INCLUDED

#pragma once

/**
 * This template class implements a vector (growable array) of elements + heap structure
 *	- requires > and == operator
 */
template<class DataItem>
class DataHeapVector
{
public:
	/// Standard constructor
	DataHeapVector(int n = 10){
		if(n < 1) n = 1;
		m_size = n;
		m_step = n;
		m_count = 0;
		m_array = new DataItem[n];
	}
	/// Copying constructor
	DataHeapVector(const DataHeapVector& dv){
		m_size = dv.m_size;
		m_step = dv.m_step;
		m_count = dv.m_count;
		m_array = new DataItem[m_size];
		for(int i = 0; i < m_count; i++)
			m_array[i] = dv.m_array[i];
	}
	/// Standard destructor (destroys all elements)
	~DataHeapVector(){ delete[] m_array; }
public:
	/// Adds new element into the vector
	int add(const DataItem& item){
		if(m_count == m_size){
			DataItem* old_array = m_array;
			int step = std::max(m_step, m_size/4);
			m_array = new DataItem[m_size + step];
			for(int i = 0; i < m_size; i++) m_array[i]=old_array[i];
			delete[] old_array;
			m_size += step;
		}
		m_array[m_count++] = item;
		upheap(m_count-1);
		return m_count-1;
	}
	/// Prepares given amount of free place
	void prepare(int ct){
		if(m_count+ct > m_size){
			DataItem* old_array = m_array;
			m_array = new DataItem[m_count+ct];
			for(int i = 0; i < m_size; i++) m_array[i]=old_array[i];
			delete[] old_array;
			m_size = m_count + ct;
		}
	}
	/// Adds new element into the vector if there is no such element yet
	bool addIfNew(const DataItem& item){
		if(find(item) > -1) 
			return false;
		else {
			add(item);
			return true;
		}
	}
	/// Returns element at the given index (doesn't check for "out of bounds" in release mode)
	DataItem& operator[](int i){ 
		assert(i >= 0 && i < m_count); 
		return m_array[i];
	}
	/// Returns last element
	DataItem& first(){ 
		assert(m_count > 0); 
		return m_array[0];
	}
	/// Assignment operator
	const DataHeapVector& operator=(const DataHeapVector& dv){
		if(&dv != this){
			delete[] m_array;
			m_size = dv.m_size;
			m_step = dv.m_step;
			m_count = dv.m_count;
			m_array = new DataItem[m_size];
			for(int i = 0; i < m_count; i++)
				m_array[i] = dv.m_array[i];
		}
		return dv;
	}
	/// Returns and removes the given element
	DataItem remove(int i = 0){ 
		assert(i >= 0 && i < m_count); 
		DataItem item = m_array[i];
		if(i != m_count-1) m_array[i] = m_array[--m_count];
		else --m_count;
		if(m_count > 0) downheap(0);
		return item;
	}
	/// Returns element at the given index (doesn't check for "out of bounds" in release mode)
	const DataItem& get(int i) const { 
		assert(i >= 0 && i < m_count); 
		return m_array[i];
	}
	/// Removes all elements
	void clear() { m_count = 0; }
	/// Is empty?
	bool empty() const { return m_count == 0; }
	/// Is not empty?
	bool notEmpty() const { return m_count > 0; }
	/// Returns the number of elements in the container
	int count() const { return m_count; }
	/// Returns the number of elements in the container
	int countInt() const { return m_count; }
	/// Finds the given item in the array (returns -1 if not found)
	int find(const DataItem& item) const {
		for(int i = 0; i < m_count; i++)
			if(m_array[i] == item) return i;
		return -1;
	}
	/// updates position of element
	void update(int k) {
		downheap(k);
		upheap(k);
	}
protected:	
	/// Restores the heap-condition for the k-th element (downward)
	void downheap(int k){
		assert(k >= 0 && k < m_count);
		// v - element to order
		DataItem v = m_array[k];
		// StepPIng down the heap
		int ndiv2 = m_count / 2;
		while( k < ndiv2 ){
			int j = 2*k+1;	// Next element of k-th in the heap
			if(j+1 < m_count && (m_array[j+1] > m_array[j])) ++j;
			// if( v >= a[j]) break;
			if(m_array[j] > v) break;
			m_array[k] = m_array[j];
			k = j;
		}
		m_array[k] = v;
	}
	/// Restores the heap-condition for the k-th element (upward)
	void upheap(int k){
		// v - element do order
		DataItem v = m_array[k];
		int l = (k-1) / 2;	// Previous element of k-th in the heap
		do{
			if(m_array[l] > v) break;
			m_array[k] = m_array[l];
			k = l;
			l = (l-1) / 2;
		}while(l != k);
		m_array[k] = v;
	}
protected:
	/// Current size of the array 
	int	m_size;
	/// Factor of size for the array
	int	m_step;
	/// Main array
	DataItem *m_array;
	/// Total number of elements
	int	m_count;
};

#endif // !defined(DATAHEAPVECTOR_H__INCLUDED)
