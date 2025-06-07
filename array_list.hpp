// *****************************************************************************************************************************
// array_list.hpp
// Array List Template Class
// Self-resizing dynamic array
// Author(s): Cory Douthat
// Copyright (c) 2022 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************


// TODO: allocate() already zeros out classes to prevent destructor from being called, 
//		 which would cause double deletion for pointers, but what about these cases?
//			1 Structs with destructors
// 		    2 Objects with copy constructors / overloaded assignment operators that make copies of pointer memory
//			3 Nested ArrayList objects - could copy all the data without freeing the old?
// For 2 and 3 - memcpy should take care of that since it bypasses assignment operators and copy constructors?
// But, for 1 - we don't currently have coverage?


#ifndef ARRAY_LIST_HPP_
#define ARRAY_LIST_HPP_

#include <cstdlib>
#include <vector>
#include <utility>

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
	template <typename... Args>
	bool insertEmplace(unsigned int index, Args&&... args);		// Instantiate object in place
	unsigned int insertSorted(const T& item, bool order, bool dup);	// Insert item into sorted list

	int push(const T& item) { return (insert(count, item) ? getCount() - 1 : -1); }	// Push item on end
	template <typename... Args>
	int pushEmplace(Args&&... args);							// Instantiate object in place at back
	bool pop() { return remove(count - 1); }					// Pop/remove end item

	const T* getData()const { return data; }					// Get pointer to data
	void copyData(const ArrayList<T>& b);						// Copy data, allocate only if necessary
	
	bool remove(unsigned int index);							// Remove item at index

	void clear();												// Clear list
	void free();												// Clear and de-allocate memory

	// READ FUNCTIONS
	T& get(unsigned int index) { return (*this)[index]; }		// Get reference to item at index
	T& getLast() { return (*this)[count - 1]; }					// Get reference to end item


	// std::vector functions
	ArrayList(const std::vector<T> &copy) : ArrayList() { *this = copy; }	// Constructor: copy
	const ArrayList<T>& operator =(const std::vector<T>& b);	// Operator =
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
	free();

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
	// Note: be very careful with classes, structs, or pointers
	//	     the copy / delete operations below could cause pointer double-deletion, etc.


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
			if (std::is_destructible<T>::value)		// Does T have a destructor? If so...
				memset(temp, 0, count * sizeof(T));	// prevent delete from calling destructors
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

	addOneCount();

	// TODO: make a shiftData() function since this is done in multiple places?
	// TODO: could make this more efficient by incorporating the data shift into the re-allocation step
	if (count > 1)
	{
		// Move items up
		for (unsigned int i = count - 1; i > index; i--)
		{
			// Use memcpy to avoid calling assignment operator
			memcpy(&data[i], &data[i - 1], sizeof(T));

			if (i == 0)
				break;
		}
	}

	data[index] = item;		// A copy constructor is expected here

	return true;
}

// Instantiate object in place
template <typename T>
template <typename... Args>
bool ArrayList<T>::insertEmplace(unsigned int index, Args&&... args)
{
	// Highest insert allowed is right after last item
	if (index > count)
		return false;

	addOneCount();

	// TODO: make a shiftData() function since this is done in multiple places?
	// TODO: could make this more efficient by incorporating the data shift into the re-allocation step
	if (count > 1)
	{
		// Move items up
		for (unsigned int i = count - 1; i > index; i--)
		{
			// Use memcpy to avoid calling assignment operator
			memcpy(&data[i], &data[i - 1], sizeof(T));

			if (i == 0)
				break;
		}
	}
	
	new (&data[index]) T(std::forward<Args>(args)...);	// Using "placement new" feature

	return true;
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
		unsigned int i = count - 1;
		if (dup || data[i] != item)
		{
			insert(count, item);
			return count;
		}
		else
			return i;
	}
}


// Instantiate object in place at back
template <typename T>
template <typename... Args>
int ArrayList<T>::pushEmplace(Args&&... args) 
{ 
	return (insertEmplace(count, std::forward<Args>(args)...) ? getCount() - 1 : -1);
}

// Empty array and copy data from another array, only allocate if needed
// TODO: how is this different from the assignment operator? Just doesn't re-allocate?
template <typename T>
void ArrayList<T>::copyData(const ArrayList<T>& b)
{
	// Note: This can be dangerous if copying objects with pointers. 
	//	     Most common issues is double deletion.
	//static_assert(!std::is_class_v<T>, "ArrayList<T>::copyData disabled for class/struct types due to pointer duplication risk");

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


// Clear list, but keep memory
template <typename T>
void ArrayList<T>::clear()
{
	for (unsigned int i = 0; i < count; i++)
	{
		if (std::is_destructible<T>::value)	
			data[i].~T();
	}

	memset(data, 0, count * sizeof(T));

	count = 0;
}


// Clear and de-allocate memory
template <typename T>
void ArrayList<T>::free()
{
	if (data)
		delete[] data;	// Calls destructors

	data = nullptr;
	count = 0;
	size = 0;
}


// std::vector functions
template <typename T>
const ArrayList<T>& ArrayList<T>::operator =(const std::vector<T>& b)
{
	free();

	if (b.size() > 0)
	{
		allocate(b.size());
		count = b.size();
		memcpy(data, b.data(), count * sizeof(T));
	}

	return *this;
}


#endif
