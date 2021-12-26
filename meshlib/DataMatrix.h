/////////////////////////////////////////////////////////////////////////////
// DataMatrix.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2005-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#include <functional>

#if !defined(DATAMATRIX_H__INCLUDED)
#define DATAMATRIX_H__INCLUDED

/**
 * This template class implements a container of elements
 *	with coherent indexes.
 * no reordering !
 * <DataItem> must provide default constructor !
 *		and properly use '=' operator
 */
template<class DataItem>
class DataMatrix
{
public:
	/// Standard constructor (container is an array of arrays[part_size] of pointers)
	DataMatrix(int row_count, int part_size, const DataItem& init_data){
		if(part_size < 1) part_size = 1;
		m_part_size = part_size;
		m_data = new MData[m_array_size = m_array_count = row_count];
		for(int i = 0; i < row_count; i++){
			m_data[i].data = new DataItem[part_size];
			m_data[i].count = part_size;
			for(int j = 0; j < part_size; j++)
				m_data[i].data[j] = init_data;
		}
		m_data_count = row_count * part_size;
	}
	/// Standard destructor (destroys all elements)
	~DataMatrix(){ 
		if(m_data){
			for(int i = 0; i < m_array_count; i++){
				delete[] m_data[i].data;
			}
			delete[] m_data;
		}
	}
private:
	/// operator=
	const DataMatrix& operator=(const DataMatrix& dm){
		if( this == &dm ) return dm;
		clear();

		m_part_size = dm.m_part_size;
		m_array_size = dm.m_array_size;
		m_array_count = dm.m_array_count;
		m_data_count = dm.m_data_count;
		m_data = new MData[m_array_size];

		for(int i = 0; i < m_array_count; i++){
			m_data[i].data = new DataItem[m_part_size];
			int ct = m_data[i].count = dm.m_data[i].count;
			for(int j = 0; j < ct; j++)
				m_data[i].data[j] = dm.m_data[i].data[j];
		}

		return *this;
	}
	/// operator= (from r-value)
	const DataMatrix& operator=(const DataMatrix&& dm){
		assert( this != &dm );
		std::swap( m_part_size, dm.m_part_size );
		std::swap( m_array_size, dm.m_array_size );
		std::swap( m_array_count, dm.m_array_count );
		std::swap( m_data_count, dm.m_data_count );
		std::swap( m_data, dm.m_data );
		return *this;
	}
	/// copying constructor
	DataMatrix(const DataMatrix& dm){
		m_part_size = dm.m_part_size;
		m_array_size = dm.m_array_size;
		m_array_count = dm.m_array_count;
		m_data_count = dm.m_data_count;
		m_data = new MData[m_array_size];

#pragma omp parallel for
		for(int i = 0; i < m_array_count; i++){
			m_data[i].data = new DataItem[m_part_size];
			int ct = m_data[i].count = dm.m_data[i].count;
			for(int j = 0; j < ct; j++)
				m_data[i].data[j] = dm.m_data[i].data[j];
		}
	}
	/// copying constructor (from r-value)
	DataMatrix(const DataMatrix&& dm){
		m_part_size = dm.m_part_size;
		m_array_size = dm.m_array_size;
		m_array_count = dm.m_array_count;
		m_data_count = dm.m_data_count;
		m_data = dm.m_data;
		// prepare dm for cleaning
		dm.m_data = nullptr;
	}
public:
	// copy values from other matrix, with the same dimensions
	bool copyFrom(const DataMatrix& dm) {
		if(m_part_size != dm.m_part_size) return false;
		if(m_array_size != dm.m_array_size) return false;
		if(m_array_count != dm.m_array_count) return false;
		if(m_data_count != dm.m_data_count) return false;

#pragma omp parallel for
		for (int i = 0; i < m_array_count; i++) {
			int ct = m_data[i].count;
			assert(ct == dm.m_data[i].count);
			for (int j = 0; j < ct; j++)
				m_data[i].data[j] = dm.m_data[i].data[j];
		}

		return true;
	}
	/// Returns/access element at the given index in two-dimensional mode (doesn't check for "out of bounds" in release mode)
	DataItem& operator()(int i, int j) { 
		assert(i*m_part_size+j >= 0 && i*m_part_size+j < m_data_count); 
		return m_data[i].data[j]; 
	}
	DataItem& get(int i, int j) { 
		assert(i*m_part_size+j >= 0 && i*m_part_size+j < m_data_count); 
		return m_data[i].data[j]; 
	}
	bool geq(int i, int j, const DataItem& v) {
		assert(i*m_part_size + j >= 0 && i*m_part_size + j < m_data_count);
		return m_data[i].data[j] >= v;
	}
	/// Set value of element at the given index in two-dimensional mode (doesn't check for "out of bounds" in release mode)
	DataItem& set(int i, int j, const DataItem& data) { 
		assert(i*m_part_size+j >= 0 && i*m_part_size+j < m_data_count); 
		return m_data[i].data[j] = data; 
	}
	/// Returns/access element at the given index in two-dimensional mode (doesn't check for "out of bounds" in release mode)
	const DataItem& operator()(int i, int j) const { 
		assert(i*m_part_size+j >= 0 && i*m_part_size+j < m_data_count); 
		return m_data[i].data[j]; 
	}
	/// Returns/access element at the given index in two-dimensional mode (doesn't check for "out of bounds" in release mode)
	const DataItem& get(int i, int j) const { 
		assert(i*m_part_size+j >= 0 && i*m_part_size+j < m_data_count); 
		return m_data[i].data[j]; 
	}
	/// Removes all elements (references AND objects) from container
	void clear(int new_part_size = 0){
		for(int i = 0; i < m_array_count; i++){
			delete[] m_data[i].data;
		}
		if(m_data) delete[] m_data;
		m_data = nullptr;
		m_array_count = m_array_size = m_data_count = 0;
		if(new_part_size > 0) m_part_size = new_part_size;
	}
	/// Returns the number of elements in the container
	int count() const { return m_data_count; }
	/// Returns the number of elements in the container
	int countInt() const { return m_data_count; }
	/// Returns the current number of rows
	int rows() const { return m_array_count; }
	/// Returns the current number of columns
	int colums() const { return m_part_size; }
	/// Whether is empty
	bool empty() const { return m_data_count == 0; }
	/// Whether is not empty
	bool notEmpty() const { return m_data_count != 0; }
	void forEach( std::function <void(DataItem& item)> f ) {
		for(int i = 0; i < m_array_count; i++){
			for(int j = 0; j < m_data[i].count; j++)
				f( m_data[i].data[j] );
		}
	}
	void forEach( std::function <void(const DataItem& item)> f ) const {
		for(int i = 0; i < m_array_count; i++){
			for(int j = 0; j < m_data[i].count; j++)
				f( m_data[i].data[j] );
		}
	}
protected:
	/// Structure for subarrays
	struct MData{
		/// Number of elements in the subarray
		int count;
		/// Subarray of elements
		DataItem*	data;
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

#endif // !defined(DATAMATRIX_H__INCLUDED)
