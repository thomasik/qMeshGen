/////////////////////////////////////////////////////////////////////////////
// DataVector.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2015- (reworked)
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(DATAVECTOR_H__INCLUDED)
#define DATAVECTOR_H__INCLUDED

#pragma once

#include "common.h"
#include <functional>
#include <algorithm>
#include <climits>

/**
 * This template class implements a vector (growable array of arrays) of elements
 */
template<class DataItem, int default_size = 10>
class DataVector
{
public:
	/// Standard constructor (container is an array of arrays[part_size])
	DataVector(size_t part_size = default_size) {
		if(part_size < default_size) part_size = default_size;
		m_part_size = part_size;
	}
	DataVector(size_t part_size, const DataItem& item ) 
		: m_part_size(part_size), m_data(new MData[1]), m_data_count(part_size),
			m_array_count(1), m_array_size(1) {
		assert(part_size > 0);
		m_data[0].count = part_size;
		m_data[0].data = new DataItem[part_size];
		for(size_t i = 0; i < part_size; i++)
			m_data[0].data[i] = item;
	}
	/// Standard destructor (destroys all elements)
	~DataVector(){ 
		if(m_data){
			for(size_t i = 0; i < m_array_count; i++){
				delete[] m_data[i].data;
			}
			delete[] m_data;
		}
	}
public:
	/// operator=
	const DataVector& operator=(const DataVector& dv){
		if( this == &dv ) return dv;
		clear();

		m_part_size = dv.m_part_size;
		m_array_size = dv.m_array_size;
		m_array_count = dv.m_array_count;
		m_data_count = dv.m_data_count;
		m_data = new MData[m_array_size];

		for(size_t i = 0; i < m_array_count; i++){
			m_data[i].data = new DataItem[m_part_size];
			size_t ct = m_data[i].count = dv.m_data[i].count;
			for(size_t j = 0; j < ct; j++)
				m_data[i].data[j] = dv.m_data[i].data[j];
		}

		return *this;
	}
	/// operator= (from r-value)
	const DataVector& operator=(DataVector&& dv){
		assert( this != &dv );
		std::swap( m_part_size, dv.m_part_size );
		std::swap( m_array_size, dv.m_array_size );
		std::swap( m_array_count, dv.m_array_count );
		std::swap( m_data_count, dv.m_data_count );
		std::swap( m_data, dv.m_data );
		return *this;
	}
	/// copying constructor
	DataVector(const DataVector& dv){
		m_part_size = dv.m_part_size;
		m_array_size = dv.m_array_size;
		m_array_count = dv.m_array_count;
		m_data_count = dv.m_data_count;
		m_data = new MData[m_array_size];

		for(size_t i = 0; i < m_array_count; i++){
			m_data[i].data = new DataItem[m_part_size];
			size_t ct = m_data[i].count = dv.m_data[i].count;
			for(size_t j = 0; j < ct; j++)
				m_data[i].data[j] = dv.m_data[i].data[j];
		}
	}
	/// copying constructor (from r-value)
	DataVector(DataVector&& dv){
		m_part_size = dv.m_part_size;
		m_array_size = dv.m_array_size;
		m_array_count = dv.m_array_count;
		m_data_count = dv.m_data_count;
		m_data = dv.m_data;
		// prepare dv for cleaning
		dv.m_data = nullptr;
	}
public:
        typedef size_t size_type;
        class iterator
        {
            public:
                typedef iterator self_type;
                typedef DataItem value_type;
                typedef DataItem& reference;
                typedef DataItem* pointer;
				typedef std::random_access_iterator_tag iterator_category;
                typedef size_t difference_type;

                iterator(DataVector<DataItem> * dv, difference_type idx) : m_dv(dv), m_idx(idx) { }
                self_type& operator++() { m_idx++; return *this; }
                self_type operator++(int junk) { self_type i = *this; m_idx++; return i; }
                self_type& operator+=(difference_type n) { m_idx += n; return *this; }
				friend self_type operator+( const self_type& it, difference_type n ) { 
					return self_type( it.m_dv, it.m_idx+n ); }
				friend self_type operator+( difference_type n, const self_type& it ) { 
					return self_type( it.m_dv, it.m_idx+n ); }
                self_type& operator--() { m_idx--; return *this; }
                self_type operator--(int junk) { self_type i = *this; m_idx--; return i; }
                self_type& operator-=(difference_type n) { m_idx -= n; return *this; }
				friend self_type operator-( const self_type& it, difference_type n ) { 
					return self_type( it.m_dv, it.m_idx-n ); }
				friend difference_type operator-( const self_type& ita, const self_type& itb) { 
					return ita.m_idx - itb.m_idx; }
                reference operator*() { return m_dv->get(m_idx); }
                pointer operator->() { return &m_dv->get(m_idx); }
                reference operator[](difference_type n) { return m_dv->get(m_idx+n); }
                bool operator==(const self_type& rhs) const { return m_idx == rhs.m_idx; }
                bool operator!=(const self_type& rhs) const { return m_idx != rhs.m_idx; }
                bool operator>(const self_type& rhs) const { return m_idx > rhs.m_idx; }
                bool operator>=(const self_type& rhs) const { return m_idx >= rhs.m_idx; }
                bool operator<(const self_type& rhs) const { return m_idx < rhs.m_idx; }
                bool operator<=(const self_type& rhs) const { return m_idx <= rhs.m_idx; }
            private:
				DataVector<DataItem> * m_dv;
                difference_type m_idx;
        };
        class const_iterator
        {
            public:
                typedef const_iterator self_type;
                typedef DataItem value_type;
                typedef DataItem& reference;
                typedef DataItem* pointer;
                typedef int difference_type;
				typedef std::random_access_iterator_tag iterator_category;
                const_iterator(const DataVector<DataItem> * dv, difference_type idx) :  m_dv(dv), m_idx(idx) { }
                self_type& operator++() { m_idx++; return *this; }
                self_type operator++(int junk) { self_type i = *this; m_idx++; return i; }
                self_type& operator+=(difference_type n) { m_idx += n; return *this; }
				friend self_type operator+( const self_type& it, difference_type n ) { 
					return self_type( it.m_dv, it.m_idx+n ); }
				friend self_type operator+( difference_type n, const self_type& it ) { 
					return self_type( it.m_dv, it.m_idx+n ); }
                self_type operator--() { m_idx--; return *this; }
                self_type operator--(int junk) { self_type i = *this; m_idx--; return i; }
                self_type& operator-=(difference_type n) { m_idx -= n; return *this; }
				friend self_type operator-( const self_type& it, difference_type n ) { 
					return self_type( it.m_dv, it.m_idx-n ); }
				friend difference_type operator-( const self_type& ita, const self_type& itb) { 
					return ita.m_idx - itb.m_idx; }
                const value_type& operator*() { return m_dv->get(m_idx); }
                const value_type* operator->() { return &m_dv->get(m_idx); }
                const value_type& operator[](difference_type n) { return m_dv->get(m_idx+n); }
                bool operator==(const self_type& rhs) const { return m_idx == rhs.m_idx; }
                bool operator!=(const self_type& rhs) const { return m_idx != rhs.m_idx; }
                bool operator>(const self_type& rhs) const { return m_idx > rhs.m_idx; }
                bool operator>=(const self_type& rhs) const { return m_idx >= rhs.m_idx; }
                bool operator<(const self_type& rhs) const { return m_idx < rhs.m_idx; }
                bool operator<=(const self_type& rhs) const { return m_idx <= rhs.m_idx; }
            private:
				const DataVector<DataItem> * m_dv;
                difference_type m_idx;
        };
        size_type size() const { return m_data_count; }
        iterator begin() { return iterator(this, 0); }
        iterator end()   { return iterator(this, m_data_count); }
        const_iterator begin() const { return const_iterator(this, 0); }
        const_iterator end() const   { return const_iterator(this, m_data_count); }
public:
	/// Removes all elements (references AND objects) from container
	void clear(size_t new_part_size = 0){
		for(size_t i = 0; i < m_array_count; i++){
			delete[] m_data[i].data;
		}
		if(m_data) delete[] m_data;
		m_data = nullptr;
		m_array_count = m_array_size = m_data_count = 0;
		if(new_part_size > 0) m_part_size = new_part_size;
	}
	size_t prepare(size_t new_count) {
		if(new_count > m_part_size * m_array_size){
			// increase
			size_t required_size = 1 + (new_count-1) / m_part_size;
			//m_array_size = std::max(std::max(required_size, (3*m_array_size)/2), m_part_size);
			m_array_size = std::max(required_size, m_array_size+2 );
			MData* temp = new MData[m_array_size];
			for(size_t j = 0; j < m_array_count; j++) temp[j] = m_data[j];
			delete[] m_data;
			m_data = temp;
		}
		return m_part_size * m_array_size;
	}
	/// Adds new elements into the container
	size_t addItems(size_t ct){
		prepare( m_data_count + ct );

		while(ct > 0){
			size_t i = m_data_count / m_part_size;
			if(i >= m_array_count){
				// New subarray
				i = m_array_count++;
				m_data[i].count = 0;
				m_data[i].data = new DataItem[m_part_size];
			}
			size_t sub_ct = std::min(ct, m_part_size - m_data[i].count);
			m_data[i].count += sub_ct;
			m_data_count += sub_ct;
			ct -= sub_ct;
		}
		return m_data_count;	// current number of items
	}
	/// Adds new elements into the container
	size_t addItems(size_t ct, const DataItem& item){
		prepare( m_data_count + ct );

		while(ct > 0){
			size_t i = m_data_count / m_part_size;
			if(i >= m_array_count){
				// New subarray
				i = m_array_count++;
				m_data[i].count = 0;
				m_data[i].data = new DataItem[m_part_size];
			}
			size_t sub_ct = std::min(ct, m_part_size - m_data[i].count);
			for(size_t j = 0; j < sub_ct; j++)
				m_data[i].data[ m_data[i].count++ ] = item;
			m_data_count += sub_ct;
			ct -= sub_ct;
		}
		return m_data_count;	// current number of items
	}
	/// Adds new element into the container
	int add(const DataItem& item){
		prepare( m_data_count + 1 );
		size_t i = m_data_count / m_part_size;
		if(i >= m_array_count){
			// New subarray
			i = m_array_count++;
			assert( m_array_count <= m_array_size );
			m_data[i].count = 0;
			m_data[i].data = new DataItem[m_part_size];
		}
		m_data[i].data[m_data[i].count++] = item;
		m_data_count++;

		size_t index = m_part_size * i + m_data[i].count - 1;	// index of the added item
		assert(index <= INT_MAX);
		return (int)index;
	}
	/// Adds new element into the vector if there is no such element yet
	bool addIfNew(const DataItem& item){
		if(find(item) != std::string::npos)
			return false;
		else {
			add(item);
			return true;
		}
	}
	size_t addAllIfNew(const DataVector<DataItem> & dv) {
		size_t ct = 0;
		for(size_t i = 0; i < dv.countInt(); i++){
			if(this->addIfNew(dv[i])) ct++;
		}
		return ct;
	}
	/// Returns and removes the last element
	DataItem removeLast(){ 
		assert( notEmpty() );
		--m_data_count;
		size_t row = m_data_count / m_part_size;
		m_data[row].count--; assert( m_data[row].count >= 0 );
		return m_data[row].data[ m_data_count % m_part_size ];
	}
	/// Removes the given element
	void remove(const DataItem& item){ 
		size_t index = find(item);
		if(index >= 0){
			if(index != m_data_count-1) 
				set( index, removeLast() );
			else removeLast();
		}
	}
	/// Returns and removes the given element
	DataItem removeAt(size_t i){
		assert(i >= 0 && i < m_data_count); 
		DataItem item = get(i);
		if(i != m_data_count-1) 
			set( i, removeLast() );
		else removeLast();
		return item;
	}
	/// Returns and removes the given element (keePIng order of elements)
	DataItem removeOrderedAt(size_t i){
		assert(i >= 0 && i < m_data_count); 
		DataItem item = get(i);
		--m_data_count;
		if( i < m_data_count ) {
			size_t row = i / m_part_size;
			size_t col = i % m_part_size;
			for( ; row < m_array_count; row++ ){
				for( ; col < m_data[row].count-1; col++ ){
					m_data[row].data[col] = m_data[row].data[col+1];
				}
				if( row < m_array_count-1 && m_data[row+1].count > 0 ) {
					m_data[row].data[col] = m_data[row+1].data[0];
					col = 0;
				}else{
					m_data[row].count--;
				}
			}
		}
		return item;
	}
	/// Insert new element into the vector
	size_t insertAt(size_t index, const DataItem& item){
		assert( index >= 0 && index < m_data_count );
		size_t i = add( item );
		while( i > index ) {
			set( i, get(i-1) );
			i--;
		}
		assert( i == index );
		set( i, item );
		return m_data_count;
	}
	/// Whether contains the given element
	bool contains(const DataItem& item) const { return find(item) != std::string::npos; }
	/// Whether contains all elements from other container
	bool contains(const DataVector<DataItem> & dv) const {
		for(size_t i = 0; i < dv.countInt(); i++){
			if(this->find(dv[i]) == std::string::npos) return false;
		}
		return true;
	}
	/// Finds the given item in the array (returns npos if not found)
	size_t find(const DataItem& item) const {
		for(size_t row = 0; row < m_array_count; row++)
			for(size_t i = 0; i < m_data[row].count; i++)
				if( m_data[row].data[i] == item ) return row * m_part_size + i;
		return std::string::npos;
	}
	void reverse() { std::reverse( begin(), end() ); }
void switchData(size_t i, size_t j) {
	std::swap(get(i), get(j));
}
//void reverse(){
//	int half = m_data_count / 2;
//	// Reverse points-sequence
//	for(int k = 0; k < half; k++){
//		int k2 = m_data_count - k - 1;
//		std::swap( get(k), get(k2) );
//	}
//}
///// random shuffle of elements
//void shuffle() {
//	int last = m_data_count-1;
//	for(int i = 0; i < last; i++){
//		int ri = rand() % (m_count-i);
//		if(ri > 0) switchData(i, i+ri);
//	}
//}
/// Returns element at the given index (doesn't check for "out of bounds" in release mode)
const DataItem& get(size_t i) const {
	assert(i < m_data_count);
	return m_data[i / m_part_size].data[i % m_part_size];
}
/// Returns/access element at the given index (doesn't check for "out of bounds" in release mode)
DataItem& get(size_t i) {
	assert(i < m_data_count);
	return m_data[i / m_part_size].data[i % m_part_size];
}
/// Returns/access element at the given index (doesn't check for "out of bounds" in release mode)
DataItem& operator[](size_t i) {
	assert(i < m_data_count);
	return m_data[i / m_part_size].data[i % m_part_size];
}
/// Returns/access element at the given index (doesn't check for "out of bounds" in release mode)
const DataItem& operator[](size_t i) const {
	assert(i < m_data_count);
	return m_data[i / m_part_size].data[i % m_part_size];
}
/// Returns/access element at the given index (doesn't check for "out of bounds" in release mode)
DataItem& set(size_t i, const DataItem& data) {
	assert(i < m_data_count);
	return m_data[i / m_part_size].data[i % m_part_size] = data;
}
void setAll(const DataItem& data) {
	for (size_t i = 0; i < m_data_count; i++) set(i, data);
}
/// Returns the number of elements in the container
size_t count() const { return m_data_count; }
/// Returns the number of elements in the container
int countInt() const { assert(m_data_count <= INT_MAX);  return (int)m_data_count; }
/// Removes excess of nodes, leaving only the first ct data
void leaveOnly(size_t new_ct) {
	if (new_ct <= 0) clear();
	else if (new_ct < m_data_count) {
		size_t new_rows = 1 + ((new_ct - 1) / m_part_size);
		while (m_array_count > new_rows)
			delete[] m_data[--m_array_count].data;
		m_data[m_array_count - 1].count = 1 + ((new_ct - 1) % m_part_size);
		m_data_count = new_ct;
	}
}
/// Sets exactly the given number of data
void setExactly(size_t new_ct) {
	if (new_ct < m_data_count)
		leaveOnly(new_ct);
	else if (new_ct > m_data_count)
		addItems(new_ct - m_data_count);
}
/// Returns last element
DataItem& last() {
	assert(m_data_count > 0);
	return get(m_data_count - 1);
}
/// Returns last element
const DataItem& last() const {
	assert(m_data_count > 0);
	return get(m_data_count - 1);
}
/// Returns the first element
DataItem& first() {
	assert(m_data_count > 0);
	return get(0);
}
/// Returns the first element
const DataItem& first() const {
	assert(m_data_count > 0);
	return get(0);
}
/// Whether is empty
bool empty() const { return m_data_count == 0; }
/// Whether is not empty
bool notEmpty() const { return m_data_count != 0; }
public:
	void forEach(std::function <void(DataItem& item)> f) {
		for (size_t i = 0; i < m_data_count; i++) f(get(i));
	}
	void forEach(std::function <void(const DataItem& item)> f) const {
		for (size_t i = 0; i < m_data_count; i++) f(get(i));
	}
	void removeIf(std::function<bool ( const DataItem& item) > f) {
		for (size_t i = 0; i < m_data_count; )
			if (f(get(i)))
				removeAt(i);
			else
				i++;
	}
	void removeOrderedIf(std::function<bool(const DataItem& item) > f) {
		assert(false);
		// to do if needed
	}
public:
	friend std::ostream& operator<<(std::ostream& os, const DataVector<DataItem,default_size> & dv ) {
		os << "DataVector (#" << dv.m_data_count << " in " << dv.m_array_count 
			<< " of " << dv.m_array_size << " x " << dv.m_part_size << ") -> [";
		for(size_t i = 0; i < dv.m_data_count; i++)
			os << dv.get(i) << " ";
		return os << "]";
	}
protected:
	/// Structure for subarrays
	struct MData{
		/// Number of elements in the subarray
		size_t count;
		/// Subarray of elements
		DataItem*	data;
	};
protected:
	/// Factor of size for the main array and size of subarrays
	size_t	m_part_size = 10;
	/// Main array (of subarrays)
	MData *m_data = nullptr;
	/// Total number of elements
	size_t	m_data_count = 0;
	/// Current number of created subarrays
	size_t	m_array_count = 0;
	/// Current size of the main array (>= m_array_count)
	size_t	m_array_size = 0;
};

#endif // !defined(DATAVECTOR_H__INCLUDED)
