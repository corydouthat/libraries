// *****************************************************************************************************************************
// linkedlist.hpp
// Linked List Template Class
// Double-Linked
// Adapted from class written in 2007 at UMKC
// Author(s): Cory Douthat
// Copyright (c) 2020 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef LINKED_LIST_HPP_
#define LINKED_LIST_HPP_

#include <cstdlib>

template <typename T>
class phLinkedList;

template <typename T>
class phListNode
{
private:
	// DATA
	T data;
	phListNode<T> *next;
	phListNode<T> *prev;
	// CONSTRUCTOR
	phListNode(const T &input) : data(input),next(0),prev(0) {}

public:
	// LINKED LIST FRIEND
	friend class phLinkedList<T>;
};

template <typename T>
class phLinkedList
{
private:
	unsigned int count;
	phListNode<T> *head;						// Pointer to head of list
	phListNode<T> *tail;						// Pointer to tail of list
	phListNode<T> *current;						// Pointer to current item of interest
public:
	// CONSTRUCTORS AND DESTRUCTOR
	phLinkedList() : head(0),tail(0),current(0),count(0) {}		// Constructor: default
	phLinkedList(const phLinkedList<T> &copy);					// Constructor: copy
	~phLinkedList() {clear();}									// Destructor

	// OPERATORS
	const phLinkedList<T>& operator =(const phLinkedList<T> &b);			// Operator =
	T& operator [](unsigned int i) { return getIndex(i); }					// Operator []

	// CHECK FUNCTIONS
	bool isEmpty()const { return head==nullptr; }								// Empty check
	unsigned int size()const { return count; }								// Size

	// DATA CONTROL FUNCTIONS
	T* insertHead(const T &item);											// Add an item at the head
	T* insertTail(const T &item);											// Add an item at the tail
	T* insertCurrent(const T &item);										// Add an item after the current node
	T* insertSorted(const T& item, bool order, bool dup);					// Add an item to a sorted list

	void removeHead();														// Delete item at head and set head to next item
	void removeTail();														// Delete item at tail
	void removeCurrent();													// Delete item at current pointer
	void removeIndex(unsigned int index);									// Delete item at index
	void remove(const phListNode<T> *node);									// Delete specified node

	// TODO: get functions throw exception if list is empty, but not sure how to fix 
	T& getHead() { return head->data; }										// Return reference to data at head
	T& getCurrent() { return current->data; }								// Return reference to data at current
	T& getTail() { return tail->data; }										// Return reference to data at tail
	T& get(const phListNode<T> *node) { return node->data; }				// Return reference to data at specified node
	T& getIndex(unsigned int index);										// Return reference to item at an index from head

	void setHead(const T &d) { head->data = d; }							// Set data at head
	void setCurrent(const T &d) { current->data = d; }						// Set data at current
	void set(const phListNode<T> *node, const T &d) { node->data = d; }		// Set data at specified node
	void setIndex(unsigned int index, const T &d);							// Set data at index from head

	void clear();															// clear list

	// POINTER FUNCTIONS
	bool moveNext();														// Set current pointer to the next item
	bool movePrev();														// Set current pointer to the prev item
	bool moveHead();														// Set current pointer to the head item
	bool moveTail();					                                    // Set current pointer to the tail item
	bool moveIndex(unsigned int index);										// Set current pointer to the index item
	const phListNode<T>* getCurrentPointer()const { return current; }		// Get constant pointer to current node
	const phListNode<T>* getTailPointer()const { return tail; }				// Get constant pointer to tail node
};


// CONSTRUCTORS
template <typename T>
phLinkedList<T>::phLinkedList(const phLinkedList<T> &copy) : head(0),tail(0),current(0),count(0)
{
	*this = copy;
}


// OPERATORS

// Operator = (assignment)
template <typename T>
const phLinkedList<T>& phLinkedList<T>::operator =(const phLinkedList<T> &b)
{
	clear();

	phListNode<T> *temp_node = b.head;
	for (unsigned int i = 0; i < b.count; i++)
	{
		insertTail(temp_node->data);
		temp_node = temp_node->next;
	}
	return *this;
}

// DATA CONTROL
// Add an item at the head
// Returns:	Pointer to new data that was added
template <typename T>
T* phLinkedList<T>::insertHead(const T &item)
{
	if (isEmpty())
	{
		head = tail = current = new phListNode<T>(item);
	}
	else
	{
		phListNode<T> *temp;
		temp = new phListNode<T>(item);
		temp->next = head;
		head->prev = temp;
		head = temp;
	}
	count++;
	return &(head->data);
}

// Add an item at the tail
// Returns:	Pointer to new data that was added
template <typename T>
T* phLinkedList<T>::insertTail(const T &item)
{
	if (isEmpty())
	{
		head = tail = current = new phListNode<T>(item);
	}
	else
	{
		phListNode<T> *temp;
		temp = new phListNode<T>(item);
		tail->next = temp;
		temp->prev = tail;
		tail = temp;
	}
	count++;
	return &(tail->data);
}

// Add an item after the current node
// Returns:	Pointer to new data that was added
template <typename T>
T* phLinkedList<T>::insertCurrent(const T &item)
{
	if (current == tail)
		return insertTail(item);
	else
	{
		phListNode<T>* temp;
		temp = new phListNode<T>(item);
		temp->next = current->next;
		current->next->prev = temp;
		temp->prev = current;
		current->next = temp;
		count++;
		return &(temp->data);
	}
}

// Add an item to a sorted list
// order = false for ascending, true for descending
// dup = false for don't allow duplicates, true for allow
// TODO: could be faster if starting node/index could be provided
// TODO: could be faster with binary-style search
template <typename T>
T* phLinkedList<T>::insertSorted(const T& item, bool order, bool dup)
{
	if (isEmpty())
		return insertTail(item);
	else
	{
		if (order == false) // Ascending
		{
			if (current->data <= item)
				moveTail();
			do
			{
				if (current->data < item)
					return insertCurrent(item);
				else if (current->data == item)
				{
					if (dup)
						return insertCurrent(item);
					else
						return nullptr;
				}
			} while (movePrev());
			// Reaching the end implies head->data is > item
			return insertHead(item);
		}
		// TODO: descending option untested
		else	// Descending
		{
			if (current->data >= item)
				moveTail();
			do
			{
				if (current->data > item)
					return insertCurrent(item);
				else if (current->data == item)
				{
					if (dup)
						return insertCurrent(item);
					else
						return nullptr;
				}
			} while (movePrev());
			// Reaching the end implies head->data is < item
			return insertHead(item);
		}
	}
}

// Delete item at head and set head to next item
template <typename T>
void phLinkedList<T>::removeHead()
{
	if (head == nullptr)
		return;

	if (head == tail)
	{
		delete head;
		head = tail = current = nullptr;
		count = 0;
	}
	else
	{
		if (current == head)
			current = head->next;
		phListNode<T> *temp = head->next;
		delete head;
		head = temp;
		head->prev = nullptr;
		count--;
	}
}

// Delete item at tail
template <typename T>
void phLinkedList<T>::removeTail()
{
	if (head == nullptr)
		return;

	if (head == tail)
	{
		delete tail;
		head = tail = current = nullptr;
		count = 0;
	}
	else
	{
		if (current == tail)
			current = nullptr;
		phListNode<T> *temp = tail->prev;
		delete tail;
		tail = temp;
		tail->next = nullptr;
		count--;
	}
}

// Delete item at current pointer
template <typename T>
void phLinkedList<T>::removeCurrent()
{
	if (head == nullptr)
		return;

	remove(current);
}

// Delete item at index
template <typename T>
void phLinkedList<T>::removeIndex(unsigned int index)
{
	if (index >= count)
		return;

	phListNode<T>* temp = current;

	moveIndex(index);

	if (current == temp)
		removeCurrent()
	else
	{
		removeCurrent();
		current = temp;
	}
}

// Delete specific node
template <typename T>
void phLinkedList<T>::remove(const phListNode<T> *node)
{
	if (head == nullptr)
		return;

	if (node == head)
		removeHead();
	else if (node == tail)
		removeTail();
	else
	{
		if (current == node)
			current = node->next;
		node->prev->next = node->next;
		node->next->prev = node->prev;
		delete node;
		count--;
	}
}

// Return reference to item at an index from head
template <typename T>
T& phLinkedList<T>::getIndex(unsigned int index)
{
	phListNode<T>* temp = current;
	T* temp_data;

	if (moveIndex(index))
	{
		temp_data = &(getCurrent());
		current = temp;
		return *temp_data;
	}
	else
		return tail->data;	// TODO: not ideal
}

// Set data at index from head
template <typename T>
void phLinkedList<T>::setIndex(unsigned int index,const T &d)
{
	phListNode<T>* temp = current;

	if (moveIndex(index))
	{
		setCurrent(d);
		current = temp;
	}
}

// Clear list
template <typename T>
void phLinkedList<T>::clear()
{
	while (!isEmpty())
		removeHead();
}

// POINTER FUNCTIONS
// Set current pointer to the next item
// Returns:	true if next is valid, false if at tail
template <typename T>
bool phLinkedList<T>::moveNext()
{
	if (current == tail)
		return false;
	else
	{
		current = current->next;
		return true;
	}
}

// Set current pointer to the prev item
// Returns:	true if prev is valid, false if at head
template <typename T>
bool phLinkedList<T>::movePrev()
{
	if (current == head)
		return false;
	else
	{
		current = current->prev;
		return true;
	}
}

// Set current pointer to the first item
// Returns:	true if head is valid, false if head is null
template <typename T>
bool phLinkedList<T>::moveHead()
{
	if (head == nullptr)
		return false;
	else
	{
		current = head;
		return true;
	}
}

// Set current pointer to the last item
// Returns:	true if tail is valid, false if tail is null
template <typename T>
bool phLinkedList<T>::moveTail()
{
	if (tail == nullptr)
		return false;
	else
	{
		current = tail;
		return true;
	}
}

// Set current pointer to the index item
// Returns:	true if index is valid, false if not
template <typename T>
bool phLinkedList<T>::moveIndex(unsigned int index)
{
	if (index >= count)
		return false;
	else
	{
		if (index <= count / 2)
		{
			moveHead();
			for (unsigned int i = 0; i < index; i++)
				moveNext();
		}
		else
		{
			moveTail();
			for (unsigned int i = count - 1; i > index; i--)
				movePrev();
		}

		return true;
	}
}

#endif
