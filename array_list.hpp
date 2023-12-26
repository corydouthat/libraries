// *****************************************************************************************************************************
// array_list.hpp
// Array List Template Class
// Self-resizing dynamic array
// Author(s): Cory Douthat
// Copyright (c) 2022 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef ARRAY_LIST_HPP_
#define ARRAY_LIST_HPP_

#include <cstdlib>

template <typename T>
class ArrayList
{
private:
	T* data;
	unsigned int count;
	unsigned int size;
public:
	// CONSTRUCTORS AND DESTRUCTOR
	ArrayList() : data(nullptr), count(0), size(0) {}			// Constructor: default
	ArrayList(unsigned int c) : ArrayList() { allocate(c); }	// Constructor: w/ allocation
	ArrayList(const ArrayList<T>& copy) : ArrayList() { *this = copy; }	// Constructor: copy
	ArrayList(const T in[], unsigned int len);					// Constructor: array
	~ArrayList() { free(); }									// Destructor

	// OPERATORS
	const ArrayList<T>& operator =(const ArrayList<T> &b);		// Operator =
	T& operator [](unsigned int i) { return data[i]; }			// Operator []
	const T& operator [](unsigned int i)const { return data[i]; }	// Operator [] (const)

	// CHECK FUNCTIONS
	bool isEmpty()const { return count == 0; }					// Empty check
	unsigned int getCount()const { return count; }				// Count
	unsigned int getSize()const { return size; }				// Size

	// DATA CONTROL FUNCTIONS
	bool updateCount(unsigned int new_count);					// Update list count
	bool addOneCount() { return updateCount(count + 1); }		// Add one to count (and re-allocate if necessary)
	bool allocate(unsigned int new_size);						// Re-allocate to new size

	bool set(unsigned int index, const T& item);				// Set value at index
	bool setAll(int item);										// Set value of all
	bool assignAll(const T& item);								// Assign value of all (any data type)

	bool insert(unsigned int index, const T& item);				// Insert item at index
	unsigned int insertSorted(const T& item, bool order, bool dup);	// Insert item into sorted list

	void copyData(const ArrayList<T>& b);						// Copy data, only allocate if needed

	unsigned int push(const T& item) { return insert(count, item); }	// Push item on end
	bool pop() { return remove(count - 1); }					// Pop/remove end item
	
	bool remove(unsigned int index);							// Remove item at index

	void clear() { count = 0; }									// Clear list
	void free();												// Clear and de-allocate memory

	// READ FUNCTIONS
	T& get(unsigned int index) { return (*this)[index]; }		// Get reference to item at index
	T& getLast() { return (*this)[count - 1]; }					// Get reference to end item
};

// Constructor: array
template<typename T>
inline ArrayList<T>::ArrayList(const T in[], unsigned int len) : ArrayList()
{
	if (len == 0)
		return;

	updateCount(len);

	memcpy(data, in, len * sizeof(T));
}

// Operator =
template <typename T>
const ArrayList<T>& ArrayList<T>::operator =(const ArrayList<T>& b)
{
	clear();

	allocate(b.size);

	if (!(b.isEmpty()))
	{
		count = b.count;
		memcpy(data, b.data, count * sizeof(T));
	}

	return *this;
}

// Update list count
template <typename T>
bool ArrayList<T>::updateCount(unsigned int new_count)
{
	if (new_count > count)
	{
		allocate(new_count);	// Allocate more memory if needed
		count = new_count;
		return true;
	}
	else
		return false;
}

// Re-allocate to new size
template <typename T>
bool ArrayList<T>::allocate(unsigned int new_size)
{
	if (new_size <= size)
		return false;
	else
	{
		T* temp;
		unsigned int remainder;

		// Re-size
		if (new_size <= 100)
			new_size = 100;
		else if (new_size <= 200)
			new_size = 200;
		else if (new_size <= 500)
			new_size = 500;
		else if (new_size <= 1000)
			new_size = 1000;
		else
		{
			remainder = new_size % 1000;

			if (remainder > 0)
				new_size += 1000 - remainder;
		}

		// Copy data to new memory location
		if (count > 0)
		{
			temp = data;
			data = new T[new_size];
			memcpy(data, temp, count * sizeof(T));
			if (std::is_class<T>::value)	// Is T a class?
				memset(temp, 0, count * sizeof(T));	// Prevent delete from calling class destructors
			delete[] temp;
		}
		else
			data = new T[new_size];

		size = new_size;

		return true;
	}
}

// Set value at index
template <typename T>
bool ArrayList<T>::set(unsigned int index, const T& item)
{
	if (index >= count)
		return false;
	else
	{
		data[index] = item;
		return true;
	}
}

// Set value of all (integer)
template <typename T>
bool ArrayList<T>::setAll(int item)
{
	if (count == 0)
		return false;
	else
	{
		memset(data, item, count * sizeof(T));
		return true;
	}
}

// Assign value of all (any data type)
template <typename T>
bool ArrayList<T>::assignAll(const T& item)
{
	if (count == 0)
		return false;
	else
	{
		for (unsigned int i = 0; i < count; i++)
		{
			data[i] = item;
		}

		return true;
	}
}

// Insert item at index
// Push other items up the list
template <typename T>
bool ArrayList<T>::insert(unsigned int index, const T& item)
{
	// Highest insert allowed is right after last item
	if (index > count)
		return false;

	if (count == 0)
	{
		// Empty
		addOneCount();

		data[0] = item;

		return true;
	}
	else
	{
		addOneCount();

		// Move items up
		for (unsigned int i = count - 2; i >= index; i--)
		{
			data[i] = data[i - 1];
		}

		data[index] = item;

		return true;
	}
}

// Insert item into sorted list
// List must be pre-sorted (TODO: add bool sorted variable)
template <typename T>
unsigned int ArrayList<T>::insertSorted(const T& item, bool order, bool dup)
{
	if (isEmpty())
		insert(0, item);
	else
	{
		// Check if should be inserted at beginning
		if ((dup && data[0] == item) || (!order && data[0] > item) || (order && data[0] < item))
		{
			insert(0, item);
			return 0;
		}

		for (unsigned int i = 0; i < count - 1; i++)
		{
			// Ascending
			if (!order)
			{
				if (data[i + 1] > item && (data[i] < item || dup))
				{
					// Insert between i and i+1
					insert(i + 1, item);
					return i + 1;
				}	
			}
			// Descending
			else
			{
				if (data[i + 1] < item && (data[i] > item || dup))
				{
					// Insert between i and i+1
					insert(i + 1, item);
					return i + 1;
				}
			}
		}

		// End of list
		i = count - 1;
		if (dup || data[i] != item)
		{
			insert(count, item);
			return count;
		}
	}
}

// Copy data, only allocate if needed
template <typename T>
void ArrayList<T>::copyData(const ArrayList<T>& b)
{
	clear();

	if (b.count > size)
		allocate(b.getSize());

	if (!(b.isEmpty()))
	{
		count = b.count;
		memcpy(data, b.data, count * sizeof(T));
	}
}

// Remove item at index
template <typename T>
bool ArrayList<T>::remove(unsigned int index)
{
	if (index >= count)
		return false;	// Invalid index
	else
	{
		for (unsigned int i = index; i < count - 1; i++)
		{
			data[i] = data[i + 1];
		}

		count--;

		return true;
	}
}

// Clear and de-allocate memory
template <typename T>
void ArrayList<T>::free()
{
	clear();

	delete[] data;
}

#endif
