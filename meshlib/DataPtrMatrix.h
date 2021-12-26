/////////////////////////////////////////////////////////////////////////////
// DataPtrMatrix.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DATAPTRMATRIX_H__INCLUDED)
#define DATAPTRMATRIX_H__INCLUDED

/**
 * This template class implements a container of pointers to elements
 *	with coherent indexes.
 */
template<class DataItem>
class DataPtrMatrix  
{
public:
	/// Standard constructor (matrix is an array of arrays[part_size] of pointers)
	DataPtrMatrix(int part_size){
		m_part_size = part_size;
		m_data_count = m_array_size = m_array_count = 0;
		m_data = nullptr;
	}
	/// Standard destructor (destroys all elements)
	~DataPtrMatrix(){ deleteAll(); }
public:
	/// Returns the element removed (reference only) from container
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
		}
		// free obsolete arrays if necessary
		while(m_array_count > 0 && m_data[m_array_count-1].count == 0){
			delete[] m_data[--m_array_count].data;
		}
		return item;
	}
	/// Adds new element into the container
	DataItem* addDataItem(DataItem* item){
		assert(item != nullptr);
		int j = m_data_count / m_part_size;
		if(j >= m_array_count){
			// New subarray
			if(++m_array_count >= m_array_size){
				// The main array is too small -> needs to be enlarged 
				//	all subarrays are transfered to the new one
				m_array_size += m_array_size/2;
				if(m_array_size < m_part_size) m_array_size = m_part_size;
				MData* temp = new MData[m_array_size];
				for(int i = 0; i < m_array_count-1; i++) temp[i] = m_data[i];
				delete[] m_data;
				m_data = temp;
			}
			m_data[j].count = 0;
			m_data[j].data = new DataItem * [m_part_size];
		}
		m_data[j].data[m_data[j].count++] = item;
		m_data_count++;

		return item;
	}
	/// Returns element at the given index (doesn't check for "out of bounds" in release mode)
	DataItem* getDataAt(int i) const { 
		assert(i >= 0 && i < m_data_count); 
		return m_data[i / m_part_size].data[i % m_part_size]; 
	}
	/// Returns element at the given index (doesn't check for "out of bounds" in release mode)
	DataItem* operator[](int i) const { 
		assert(i >= 0 && i < m_data_count); 
		return m_data[i / m_part_size].data[i % m_part_size]; 
	}
	/// Removes all elements (references AND objects) from container
	void deleteAll(){
		for(int i = 0; i < m_array_count; i++){
			for(int j = 0; j < m_data[i].count; j++){
				delete m_data[i].data[j];
			}
			delete[] m_data[i].data;
		}
		if(m_data) delete[] m_data;
		m_data = nullptr;
		m_array_count = m_array_size = m_data_count = 0;
	}
	/// Returns the number of elements in the container
	int countInt() const { return m_data_count; }
	/// Returns the number of elements in the container
	int count() const { return m_data_count; }
	/// Whether is empty
	bool empty() const { return m_data_count == 0; }
	/// Whether is not empty
	bool notEmpty() const { return m_data_count != 0; }
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
};

#endif // !defined(DATAPTRMATRIX_H__INCLUDED)
