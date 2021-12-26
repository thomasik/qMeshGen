/////////////////////////////////////////////////////////////////////////////
// DataList.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2006-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DATALIST_H__INCLUDED)
#define DATALIST_H__INCLUDED

#include <functional>

/**
 * This template class implements a simple list of elements
 */
template<class DataItem>
class DataSimpleList
{
private:
	struct ListNode
	{
		ListNode(const DataItem& it, ListNode* n = nullptr) : item(it), next(n), marked(false) {}
	public:
		DataItem item;
		ListNode* next;
		bool marked;
	};
public:
	class SimpleListIterator
	{
	public:
		SimpleListIterator(ListNode* head) : current_node(head) { }
		bool moveNext(){
			if(!current_node) return false;
			current_node = current_node->next;
			return (current_node != nullptr);
		}
		void mark() { 
			assert(current_node);
			current_node->marked = true;
		}
		bool valid() const { return current_node != nullptr; }
		DataItem& item() const {
			assert(current_node);
			return current_node->item;
		}
	private:
		ListNode* current_node;
	};
public:
	/// Standard constructor
	DataSimpleList() : head(nullptr), tail(nullptr) { }
	/// Copying constructor
private: // will be implemented if needed...
	DataSimpleList(const DataSimpleList& dl){ }
	/// Assignment operator
	const DataSimpleList& operator=(const DataSimpleList& dv){ return dv; }
public:
	/// Standard destructor (destroys all elements)
	~DataSimpleList(){ 
		while(head){ tail = head->next; delete head; head = tail; }
	}
public:
	/// Remove marked nodes of the list
	void removeMarked() {
		while(head && head->marked) {
			ListNode* t = head->next; delete head; head = t;
		}
		if(head){
			ListNode* p = head;
			while(p->next){
				if(p->next->marked){
					ListNode* t = p->next;
					p->next = p->next->next;
					delete t;
				}else
					p = p->next;
			}
			tail = p;
		}else
			tail = nullptr;
	}
	/// Returns iterator for traversing the simple list
	SimpleListIterator iterator() const { return SimpleListIterator(head); }
	/// Adds new element into the list (at the end)
	DataItem& append(const DataItem& item){
		if(tail) tail = tail->next = new ListNode(item);
		else head = tail = new ListNode(item);
		return tail->item;
	}
	/// Insert new element into the list (at the beginning)
	DataItem& insert(const DataItem& item){
		head = new ListNode(item, head);
		if(!tail) tail = head;
		return head->item;
	}
	/// Returns the first element
	DataItem& getFirst(){ 
		assert(head);
		return head->item;
	}
	/// Returns the first element
	const DataItem& getFirst() const { 
		assert(head);
		return head->item;
	}
	/// Returns the last element
	DataItem& getLast(){ 
		assert(tail);
		return tail->item;
	}
	/// Returns the last element
	const DataItem& getLast() const { 
		assert(tail);
		return tail->item;
	}
	/// Returns and removes the first element
	DataItem removeFirst(){ 
		assert(head);
		ListNode* node = head;
		DataItem item = node->item;
		head = head->next;
		if(!head) tail = nullptr;
		delete node;
		return item;
	}
	void moveFirstTo(DataSimpleList<DataItem>& other_list) {
		assert(head);
		ListNode* node = head;
		head = head->next;
		if (!head) tail = nullptr;

		node->next = other_list.head;
		other_list.head = node;
		if (!other_list.tail) other_list.tail = node;
	}
	/// Removes all elements
	void clear(){ 
		while(head){ tail = head->next; delete head; head = tail; }
		tail = nullptr;
	}
	/// Whether the list is empty
	bool empty() const { return !head; }
	/// Whether the list is not empty
	bool notEmpty() const { return head != nullptr; }
	/// Whether the list contains the given data
	bool contains(const DataItem& item) const { 
		for(ListNode* node = head; node; node = node->next)
			if(node->item == item) return true;
		return false;
	}
	/// Number of elements
	int count() const { 
		int ct = 0;
		for(ListNode* node = head; node; node = node->next) ++ct;
		return ct;
	}
	inline int countInt() const { return count(); }
public:
	/// 
	DataItem* findFirst( std::function <bool(DataItem& item)> f ) const {
		for(ListNode* node = head; node; node = node->next)
			if( f(node->item) ) return & node->item;
		return nullptr;
	}
	void forEach( std::function <void(DataItem& item)> f ) const {
		for(ListNode* node = head; node; node = node->next)
			f(node->item);
	}
protected:
	/// First element
	ListNode* head;
	/// Last element
	ListNode* tail;
};

/**
 * This template class implements a compound list of elements
 */
template<class DataItem>
class DataCompoundList
{
private:
	struct ListNode
	{
		ListNode(const DataItem& it, int s) 
			: first(0), last(0), size(s), next(nullptr), items(new DataItem[s])
		{ items[0] = it; }
		ListNode(const DataItem& it, int s, ListNode* n) 
			: first(s-1), last(s-1), size(s), next(n), items(new DataItem[s]) 
		{ items[s-1] = it; }
		~ListNode() { delete[] items; }
		// at end
		DataItem* append(const DataItem& it) { 
			if(last < size-1) { items[++last] = it; return items+last; }
			else return nullptr;
		}
		// at beginning
		DataItem* insert(const DataItem& it) { 
			if(first > 0) { items[--first] = it; return items+first; }
			else return nullptr;
		}
	public:
		int first, last, size;
		ListNode* next;
		DataItem* items;
	};
public:
	class ListIterator
	{
	public:
		ListIterator(ListNode* head) : current_node(head), current_item(head?(head->first):-1) { }
		bool moveNext(){
			if(!current_node) return false;
			if(++current_item <= current_node->last) return true;
			current_node = current_node->next;
			if(!current_node) return false;
			return (current_item=current_node->first) <= current_node->last;
		}
		bool valid() const { return current_node != nullptr && 
			current_item >= current_node->first &&
			current_item <= current_node->last; }
		DataItem& item() const {
			assert(valid());
			return current_node->items[current_item];
		}
	private:
		ListNode* current_node;
		int current_item;
	};
public:
	/// Standard constructor
	DataCompoundList(int dct = 100) : head(nullptr), tail(nullptr), m_count(0), DCOUNT(dct) { }
	/// Copying constructor
private: // will be implemented if needed...
	DataCompoundList(const DataCompoundList& dcl) = delete;
	/// Assignment operator
	const DataCompoundList& operator=(const DataCompoundList& dcl) = delete;
public:
	/// Standard destructor (destroys all elements)
	~DataCompoundList(){ 
		while(head){ tail = head->next; delete head; head = tail; }
	}
public:
	/// Returns iterator for traversing the compound list
	ListIterator iterator() const { return ListIterator(head); }
	/// Adds new element into the list (at the end)
	DataItem* append(const DataItem& item){
		DataItem* ret = nullptr;
		if(tail){
			ret = tail->append(item);
			if(!ret) tail = tail->next = new ListNode(item, DCOUNT);
		}else head = tail = new ListNode(item, DCOUNT);
		++m_count;
		return ret ? ret : tail->items; // if created new node, item should be at the first position of the tail
	}
	/// Insert new element into the list (at the beginning)
	DataItem* insert(const DataItem& item){
		DataItem* ret = nullptr;
		if(head){
			ret = head->insert(item);
			if(!ret) head = new ListNode(item, DCOUNT, head);
		}else head = tail = new ListNode(item, DCOUNT, nullptr);
		++m_count;
		return ret ? ret : head->items; // if created new node, item should be at the first position of the head
	}
	/// Returns the first element
	DataItem& getFirst(){ 
		assert(head);
		return head->items[head->first];
	}
	/// Returns the last element
	DataItem& getLast(){ 
		assert(tail);
		return tail->items[tail->last];
	}
	/// Returns and removes the first element
	DataItem removeFirst(){ 
		assert(head);
		DataItem item = head->items[head->first++];
		if(head->first > head->last){ // remove
			ListNode* node = head;
			head = head->next;
			if(!head) tail = nullptr;
			delete node;
		}
		--m_count;
		return item;
	}
	/// Removes all elements
	void clear(){ 
		while(head){ tail = head->next; delete head; head = tail; }
		tail = nullptr;
		m_count = 0;
	}
	/// Whether the list is empty
	bool empty() const { return !head; }
	/// Whether the list is not empty
	bool notEmpty() const { return head != nullptr; }
	/// Whether the list contains the given data
	bool contains(const DataItem& item) const { 
		for(ListNode* node = head; node; node = node->next)
			for(int i = node->first; i <= node->last; i++)
				if(node->items[i] == item) return true;
		return false;
	}
	/// Number of elements
	int countInt() const { return m_count; }
	/// Number of elements
	int count() const { return m_count; }
protected:
	/// First element
	ListNode* head;
	/// Last element
	ListNode* tail;
	/// Total count
	unsigned int m_count;
	/// Size of element-table
	const int DCOUNT;
};

/**
 * This template class implements a simple stack of elements on table
 */
template<class DataItem>
class DataTableStack
{
public:
	/// Standard constructor
	DataTableStack(int stack_size = 100) : m_size(stack_size), m_top(0) { 
		assert( stack_size > 0); 
		m_array = new DataItem[stack_size];
	}
	/// Copying constructor
private: // will be implemented if needed...
	DataTableStack(const DataTableStack& dts){ }
	/// Assignment operator
	const DataTableStack& operator=(const DataTableStack& dts){ return dts; }
public:
	/// Standard destructor (destroys table)
	~DataTableStack(){ 
		delete[] m_array;
	}
public:
	/// Push element on the stack
	void push(const DataItem& item){
		assert(m_top < m_size);
		m_array[m_top++] = item;
	}
	/// Pop element from the stack
	DataItem pop(){
		assert(m_top > 0);
		return m_array[--m_top];
	}
	/// Removes all elements
	void clear(){ m_top = 0; }
	/// Whether the stack is empty
	bool empty() const { return m_top == 0; }
	/// Whether the stack is full
	bool full() const { return m_top == m_size; }
	/// Whether the stack is not empty
	bool notEmpty() const { return m_top != 0; }
	/// Number of elements
	int count() const { return m_top; }
protected:
	/// size of stack array
	int m_size;
	/// top of the stack
	int m_top;
	/// array
	DataItem* m_array;
};

#endif // !defined(DATALIST_H__INCLUDED)
