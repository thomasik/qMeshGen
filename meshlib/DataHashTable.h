/////////////////////////////////////////////////////////////////////////////
// DataHashTable.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2007-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(DATAHASHTABLE_H__INCLUDED)
#define DATAHASHTABLE_H__INCLUDED

#pragma once

#include <climits>
#include "DataVector.h"
#include "DataMatrix.h"

/**
 * This template class implements a hash table (with linear probing)
 */
template<class DataItem>
class DataHashTable
{
public:
	/// Standard constructor
	DataHashTable(unsigned int n, const DataItem& zero_value) : m_count(0), m_zero(zero_value) {
		// m_size (1.5 * n) rounded up to nearest 2^m
		n += n/2;
		m_size = 1024;
		while(m_size < n) m_size <<= 1;
		m_array = new DataItem[m_size];
		for(unsigned int i = 0; i < m_size; i++) m_array[i] = zero_value;
	}
	/// Standard destructor (destroys all elements)
	~DataHashTable(){ delete[] m_array; }
private:
	/// Copying constructor -- not used
	DataHashTable(const DataHashTable& /* dv */) = delete;
	/// Assignment operator -- not used
	const DataHashTable& operator=(const DataHashTable& /* dv */) = delete;
public:
	/// Adds new element into the hash table, returns false if already is in the hashtable
	bool insert(const DataItem& item){
		unsigned int i = findSlot(item);
		if(m_array[i] == m_zero){
			m_array[i] = item;
			++m_count;
			double ratio = m_count / (double) m_size;
			if(ratio > 0.9){
				LOG4CPLUS_WARN(MeshLog::logger_console, "DataHashTable congestion " << ratio);
				DataItem* old_array = m_array;
				unsigned int old_size = m_size;
				m_size <<= 1;
				m_array = new DataItem[m_size];
				for(i = 0; i < m_size; i++) m_array[i] = m_zero;
				m_count = 0;
				for(i = 0; i < old_size; i++)
					if(old_array[i] != m_zero) insert(old_array[i]);
				delete[] old_array;
			}
			return true;
		}else return false; // already in table
	}
	/// Removes all elements
	void clear(){ 
		for(unsigned int i = 0; i < m_size; i++) m_array[i] = m_zero;
		m_count = 0; 
	}
	/// Returns the number of elements in the container
	unsigned int count() const { return m_count; }
	/// Returns the number of elements in the container
	int countInt() const { assert(m_count <= INT_MAX); return (int)m_count; }
	/// Whether contains the given item
	bool contains(const DataItem& item) const {
		unsigned int i = findSlot(item);
		return m_array[i] == item;
	}
	/// Fills with valid elements
	void getValues(DataVector<DataItem>& data) const {
		if(m_count > 0){
			data.setExactly(m_count);
			unsigned int j = 0;
			for(unsigned int i = 0; i < m_size; i++)
				if(m_array[i] != m_zero)
					data[j++] = m_array[i];
			assert( j == m_count );
		} else {
			data.clear();
		}
	}
	/// Fills with valid elements
	DataVector<DataItem> values() const {
		DataVector<DataItem> data(m_count);
		for (unsigned int i = 0; i < m_size; i++)
			if (m_array[i] != m_zero)
				data.add(m_array[i]);
		return data;
	}
protected:
	unsigned int findSlot(const DataItem& item) const {

		unsigned int i = joaat_hash((const unsigned char*)(&item), sizeof(item)) & (m_size-1);
		while(m_array[i] != m_zero && m_array[i] != item){
			i = (i + 1) & (m_size-1); 
		}
		return i;
	}
	static unsigned int joaat_hash(const unsigned char *key, unsigned int len) {
		 unsigned int hash = 0;
		 for (unsigned int i = 0; i < len; i++) {
			 hash += key[i];
			 hash += (hash << 10);
			 hash ^= (hash >> 6);
		 }
		 hash += (hash << 3);
		 hash ^= (hash >> 11);
		 hash += (hash << 15);
		 return hash;
	 }
protected:
	/// Total number of non-zero elements
	unsigned int	m_count;
	/// Zero value
	DataItem	m_zero;
	/// Current size of the array 
	unsigned int	m_size;
	/// Main array
	DataItem *m_array;
};

/**
 * This template class implements a hash table (with linear probing) - for strings
 */
template<>
class DataHashTable<string>
{
public:
	/// Standard constructor
	DataHashTable(unsigned int n, const string& zero_value) : m_count(0), m_zero(zero_value) {
		// m_size (1.5 * n) rounded up to nearest 2^m
		n += n/2;
		m_size = 1024;
		while(m_size < n) m_size <<= 1;
		m_array = new string[m_size];
		for(unsigned int i = 0; i < m_size; i++) m_array[i] = zero_value;
	}
	/// Standard destructor (destroys all elements)
	~DataHashTable(){ delete[] m_array; }
private:
	/// Copying constructor -- not used
	DataHashTable(const DataHashTable& /* dv */){
		assert(false);
	}
	/// Assignment operator -- not used
	const DataHashTable& operator=(const DataHashTable& /* dv */){
		assert(false);
		return *this;
	}
public:
	/// Adds new element into the hash table
	bool insert(const string& item){
		unsigned int i = findSlot(item);
		if(m_array[i] == m_zero){
			m_array[i] = item;
			++m_count;
			double ratio = m_count / (double)m_size;
			if(ratio > 0.9){
				LOG4CPLUS_WARN(MeshLog::logger_console, "DataHashTable congestion " << ratio);
				string* old_array = m_array;
				unsigned int old_size = m_size;
				m_size <<= 1;
				m_array = new string[m_size];
				for(i = 0; i < m_size; i++) m_array[i] = m_zero;
				m_count = 0;
				for(i = 0; i < old_size; i++)
					if(old_array[i] != m_zero) insert(old_array[i]);
				delete[] old_array;
			}
			return true;
		}else return false; // already in table
	}
	/// Removes all elements
	void clear(){ 
		for(unsigned int i = 0; i < m_size; i++) m_array[i] = m_zero;
		m_count = 0; 
	}
	/// Returns the number of elements in the container
	unsigned int count() const { return m_count; }
	/// Returns the number of elements in the container
	int countInt() const { assert(m_count <= INT_MAX); return (int)m_count; }
	/// Whether contains the given item
	bool contains(const string& item) const {
		unsigned int i = findSlot(item);
		return m_array[i] == item;
	}
	/// Fills with valid elements
	void getValues(DataVector<string>& data){
		if(m_count > 0){
			data.setExactly(m_count);
			unsigned int j = 0;
			for(unsigned int i = 0; i < m_size; i++)
				if(m_array[i] != m_zero)
					data[j++] = m_array[i];
			assert( j == m_count );
		} else {
			data.clear();
		}
	}
protected:
	unsigned int findSlot(const string& item) const {

		unsigned int i = joaat_hash((const unsigned char*)(item.c_str()), (unsigned int)item.length()) & (m_size-1);
		while(m_array[i] != m_zero && m_array[i] != item){
			i = (i + 1) & (m_size-1); 
		}
		return i;
	}
	static unsigned int joaat_hash(const unsigned char *key, unsigned int len) {
		 unsigned int hash = 0;
		 for (unsigned int i = 0; i < len; i++) {
			 hash += key[i];
			 hash += (hash << 10);
			 hash ^= (hash >> 6);
		 }
		 hash += (hash << 3);
		 hash ^= (hash >> 11);
		 hash += (hash << 15);
		 return hash;
	 }
protected:
	/// Total number of elements
	unsigned int	m_count;
	/// Zero value
	string m_zero;
	/// Current size of the array 
	unsigned int	m_size;
	/// Main array
	string *m_array;
};

/**
 * This template class implements a hash table (with linear probing)
 */
template<class DataKey, class DataValue>
class DataHashTableKeyValue
{
public:
	/// Standard constructor
	DataHashTableKeyValue(unsigned int n, const DataKey& zero_value) : m_count(0), m_zero(zero_value) {
		// m_size (1.5 * n) rounded up to nearest 2^m
		n += n/2;
		m_size = 1024;
		while(m_size < n) m_size <<= 1;
		m_key_array   = new DataKey[m_size];
		m_value_array = new DataValue[m_size];
		for(unsigned int i = 0; i < m_size; i++) m_key_array[i] = zero_value;
	}
	/// Standard destructor (destroys all elements)
	~DataHashTableKeyValue(){ delete[] m_key_array; delete[] m_value_array; }
private:
	/// Copying constructor -- not used
	DataHashTableKeyValue(const DataHashTableKeyValue& ) {}
	/// Assignment operator -- not used
	const DataHashTableKeyValue& operator=(const DataHashTableKeyValue& dv) { return dv; }
public:
	/// Adds new element into the hash table
	bool insert(const DataKey& key, const DataValue& value, unsigned int i){
		assert( i < m_size );
		assert(m_key_array[i] == m_zero);
		m_key_array[i] = key;
		m_value_array[i] = value;
		++m_count;
		double ratio = m_count / (double)m_size;
		if(ratio > 0.9){
			LOG4CPLUS_WARN(MeshLog::logger_console, "DataHashTableKeyValue congestion " << ratio);
			DataKey* old_key_array = m_key_array;
			DataValue* old_value_array = m_value_array;
			unsigned int old_size = m_size;
			m_size <<= 1;
			m_key_array = new DataKey[m_size];
			m_value_array = new DataValue[m_size];
			for(i = 0; i < m_size; i++) m_key_array[i] = m_zero;
			m_count = 0;
			for(i = 0; i < old_size; i++)
				if(old_key_array[i] != m_zero) insert(old_key_array[i], old_value_array[i]);
			delete[] old_key_array;
			delete[] old_value_array;
		}
		return true;
	}
	/// Adds new element into the hash table
	bool insert(const DataKey& key, const DataValue& value){
		unsigned int i = findSlot(key);
		if(m_key_array[i] == m_zero){
			return insert( key, value, i );
		}else return false; // already in table
	}
	/// Modify the value for the given key in the hash table
	void setValue(const DataKey& key, const DataValue& value){
		unsigned int i = findSlot(key);
		if(m_key_array[i] == m_zero){
			insert(key, value, i);
		}else{ 
			// already in table
			m_value_array[i] = value;
		}
	}
	/// Modify the value for the given key in the hash table
	void incValue(const DataKey& key, const DataValue& first_value){
		unsigned int i = findSlot(key);
		if(m_key_array[i] == m_zero){
			insert(key, first_value, i);
		}else{ 
			// already in table
			++m_value_array[i];
		}
	}
	/// Removes all elements
	void clear(){ 
		for(unsigned int i = 0; i < m_size; i++) m_key_array[i] = m_zero;
		m_count = 0; 
	}
	/// Returns the number of elements in the container
	unsigned int count() const { return m_count; }
	/// Returns the number of elements in the container
	int countInt() const { assert(m_count <= INT_MAX); return (int)m_count; }
	/// Whether contains the given item
	bool contains(const DataKey& key) const {
		unsigned int i = findSlot(key);
		return m_key_array[i] == key;
	}
	/// Whether contains the given item
	bool contains(const DataKey& key, unsigned int & index) const {
		index = findSlot(key);
		return m_key_array[index] == key;
	}
	/// Return value for given key
	const DataValue& getValue(const DataKey& key, const DataValue& zero_value) const {
		unsigned int i = findSlot(key);
		if(m_key_array[i] == key)
			return m_value_array[i];
		else return zero_value;
	}
	/// Return value for given key
	const DataValue& getValue(const DataKey& key) const {
		unsigned int i = findSlot(key);
		assert(m_key_array[i] == key);
		return m_value_array[i];
	}
	const DataValue& operator[](const DataKey& key) const { return getValue(key); }
	/// Return true if value for the given key exists and is equal to the given value
	bool checkValue(const DataKey& key, const DataValue& value) const {
		unsigned int i = findSlot(key);
		if(m_key_array[i] == key)
			return m_value_array[i] == value;
		else return false;
	}
	/// Fills with valid elements
	void getKeys(DataVector<DataKey>& data) const {
		data.setExactly(m_count);
		unsigned int j = 0;
		for(unsigned int i = 0; i < m_size; i++)
			if(m_key_array[i] != m_zero)
				data[j++] = m_key_array[i];
	}
	DataVector<DataKey> keys() const {
		DataVector<DataKey> data(m_count);
		for (unsigned int i = 0; i < m_size; i++)
			if (m_key_array[i] != m_zero)
				data.add(m_key_array[i]);
		return data;
	}
	/// Fills with valid elements
	void getValues(DataVector<DataValue>& data) const {
		data.setExactly(m_count);
		unsigned int j = 0;
		for(unsigned int i = 0; i < m_size; i++)
			if(m_key_array[i] != m_zero)
				data[j++] = m_value_array[i];
	}
	/// Fills with valid elements
	DataVector<DataValue> values() const {
		DataVector<DataValue> data(m_count);
		for (unsigned int i = 0; i < m_size; i++)
			if (m_key_array[i] != m_zero)
				data.add(m_value_array[i]);
		return data;
	}
	/// Return value for given index (low-level access)
	DataValue slotValue(unsigned int i) const { return m_value_array[i]; }
	/// Return value for given index (low-level access)
	DataKey slotKey(unsigned int i) const { return m_key_array[i]; }
	/// Return data array size (low-level access)
	unsigned int slotCount() const { return m_size; }
protected:
	/// General version
	unsigned int findSlot(const DataKey& key) const {

		unsigned int i = joaat_hash((const unsigned char*)(&key), sizeof(key)) & (m_size-1);
		while(m_key_array[i] != m_zero && m_key_array[i] != key){
			i = (i + 1) & (m_size-1); 
		}
		return i;
	}
	static unsigned int joaat_hash(const unsigned char *key, unsigned int len) {
		 unsigned int hash = 0;
		 for (unsigned int i = 0; i < len; i++) {
			 hash += key[i];
			 hash += (hash << 10);
			 hash ^= (hash >> 6);
		 }
		 hash += (hash << 3);
		 hash ^= (hash >> 11);
		 hash += (hash << 15);
		 return hash;
	 }
protected:
	/// Total number of elements
	unsigned int	m_count;
	/// Zero value
	DataKey m_zero;
	/// Current size of the array 
	unsigned int	m_size;
	/// Key array
	DataKey *m_key_array;
	/// Value array
	DataValue *m_value_array;
};

/**
 * This template class implements a hash table (with linear probing) -> special case for strings
 */
template<class DataValue>
class DataHashTableKeyValue<string, DataValue>
{
public:
	/// Standard constructor
	DataHashTableKeyValue(unsigned int n, const string& zero_value) : m_count(0), m_zero(zero_value) {
		// m_size (1.5 * n) rounded up to nearest 2^m
		n += n/2;
		m_size = 1024;
		while(m_size < n) m_size <<= 1;
		m_key_array   = new string[m_size];
		m_value_array = new DataValue[m_size];
		for(unsigned int i = 0; i < m_size; i++) m_key_array[i] = zero_value;
	}
	/// Standard destructor (destroys all elements)
	~DataHashTableKeyValue(){ delete[] m_key_array; delete[] m_value_array; }
private:
	/// Copying constructor -- not used
	DataHashTableKeyValue(const DataHashTableKeyValue& dv){
	}
	/// Assignment operator -- not used
	const DataHashTableKeyValue& operator=(const DataHashTableKeyValue& dv){
		return *this;
	}
public:
	/// Adds new element into the hash table
	bool insert(const string& key, const DataValue& value){
		unsigned int i = findSlot(key);
		if(m_key_array[i] == m_zero){
			m_key_array[i] = key;
			m_value_array[i] = value;
			++m_count;
			double ratio = m_count / (double)m_size;
			if(ratio > 0.9){
				LOG4CPLUS_WARN(MeshLog::logger_console, "DataHashTableKeyValue congestion " << ratio);
				string* old_key_array = m_key_array;
				DataValue* old_value_array = m_value_array;
				unsigned int old_size = m_size;
				m_size <<= 1;
				m_key_array = new string[m_size];
				m_value_array = new DataValue[m_size];
				for(i = 0; i < m_size; i++) m_key_array[i] = m_zero;
				m_count = 0;
				for(i = 0; i < old_size; i++)
					if(old_key_array[i] != m_zero) insert(old_key_array[i], old_value_array[i]);
				delete[] old_key_array;
				delete[] old_value_array;
			}
			return true;
		}else return false; // already in table
	}
	/// Modify the value for the given key in the hash table
	void setValue(const string& key, const DataValue& value){
		unsigned int i = findSlot(key);
		if(m_key_array[i] == m_zero){
			insert(key, value);
		}else{ 
			// already in table
			m_value_array[i] = value;
		}
	}
	/// Removes all elements
	void clear(){ 
		for(unsigned int i = 0; i < m_size; i++) m_key_array[i] = m_zero;
		m_count = 0; 
	}
	/// Returns the number of elements in the container
	unsigned int count() const { return m_count; }
	/// Returns the number of elements in the container
	int countInt() const { assert(m_count <= INT_MAX); return (int)m_count; }
	/// Whether contains the given item
	bool contains(const string& key) const {
		unsigned int i = findSlot(key);
		return m_key_array[i] == key;
	}
	/// Return value for given key
	DataValue getValue(const string& key, const DataValue& zero_value) const {
		unsigned int i = findSlot(key);
		if(m_key_array[i] == key)
			return m_value_array[i];
		else return zero_value;
	}
	/// Fills with valid elements
	void getKeys(DataVector<string>& data){
		data.setExactly(m_count);
		unsigned int j = 0;
		for(unsigned int i = 0; i < m_size; i++)
			if(m_key_array[i] != m_zero)
				data[j++] = m_key_array[i];
	}
	/// Fills with valid elements
	void getValues(DataVector<DataValue>& data){
		data.setExactly(m_count);
		unsigned int j = 0;
		for(unsigned int i = 0; i < m_size; i++)
			if(m_key_array[i] != m_zero)
				data[j++] = m_value_array[i];
	}
	/// Fills with valid elements
	DataVector<DataValue> values() const {
		DataVector<DataValue> data(m_count);
		for (unsigned int i = 0; i < m_size; i++)
			if (m_key_array[i] != m_zero)
				data.add(m_value_array[i]);
		return data;
	}
	/// Return value for given index (low-level access)
	DataValue slotValue(unsigned int i) const { return m_value_array[i]; }
	/// Return value for given index (low-level access)
	string slotKey(unsigned int i) const { return m_key_array[i]; }
	/// Return data array size (low-level access)
	unsigned int slotCount() const { return m_size; }
protected:
	unsigned int findSlot(const string& key) const {

		unsigned int i = joaat_hash((const unsigned char*)(key.c_str()), (unsigned int)key.length()) & (m_size-1);
		while(m_key_array[i] != m_zero && m_key_array[i] != key){
			i = (i + 1) & (m_size-1); 
		}
		return i;
	}
	static unsigned int joaat_hash(const unsigned char *key, unsigned int len) {
		 unsigned int hash = 0;
		 for (unsigned int i = 0; i < len; i++) {
			 hash += key[i];
			 hash += (hash << 10);
			 hash ^= (hash >> 6);
		 }
		 hash += (hash << 3);
		 hash ^= (hash >> 11);
		 hash += (hash << 15);
		 return hash;
	 }
protected:
	/// Total number of elements
	unsigned int	m_count;
	/// Zero value
	string m_zero;
	/// Current size of the array 
	unsigned int	m_size;
	/// Key array
	string *m_key_array;
	/// Value array
	DataValue *m_value_array;
};

#endif // !defined(DATAHASHTABLE_H__INCLUDED)
