// *****************************************************************************************************************************
// array_list.hpp
// Array List Template Class
// Self-resizing dynamic array
// Author(s): Cory Douthat
// Copyright (c) 2025 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************


#ifndef ARRAY_LIST_HPP_
#define ARRAY_LIST_HPP_

#include <cstdlib>
#include <vector>
#include <utility>
#include <type_traits>
#include <new>
#include <stdexcept>

template<typename T>
concept Copyable = std::is_copy_constructible_v<T> && std::is_copy_assignable_v<T>;

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
	ArrayList(const ArrayList<T>& copy) requires Copyable<T> : ArrayList() { *this = copy; }	// Constructor: copy
	ArrayList(ArrayList<T>&& other) noexcept;					// Constructor: move (rvalue)
	ArrayList(unsigned int c) : ArrayList() { allocate(c); }	// Constructor: w/ allocation
	ArrayList(const T in[], unsigned int len) requires Copyable<T>;	// Constructor: array
	~ArrayList() { free(); }									// Destructor

	// OPERATORS
	// Disable assignment operator for non-copyable types
	const ArrayList<T>& operator =(const ArrayList<T>& b) requires Copyable<T>;	// Operator =
	const ArrayList<T>& operator =(ArrayList<T>&& other) noexcept;	// Operator (move rvalue)
	T& operator [](unsigned int i) { return data[i]; }			// Operator []
	const T& operator [](unsigned int i)const { return data[i]; }	// Operator [] (const)

	// CHECK FUNCTIONS
	bool isEmpty()const { return count == 0; }					// Empty check
	unsigned int getCount()const { return count; }				// Count
	unsigned int getSize()const { return size; }				// Size

	// DATA CONTROL FUNCTIONS
	bool updateCount(unsigned int new_count);					// Update list count
	bool addOneCount() { return updateCount(count + 1); }		// Add one to count (and re-allocate if necessary)
	bool allocate(unsigned int new_size, bool downsize = false);	// Re-allocate to new size

	bool set(unsigned int index, const T& item) requires Copyable<T>;	// Set value at index
	bool set(unsigned int index, T&& item) noexcept;			// Set value at index (move operation)
	template <typename... Args>
	bool setEmplace(unsigned int index, Args&&... args);		// Emplace object at index
	bool setAllInt(int item);									// Set value of all
	bool setAll(const T& item) requires Copyable<T>;			// Assign value of all (any data type)
	bool setAll(T&& item) noexcept;								// Assign value of all (move operation)

	bool insert(unsigned int index, const T& item) requires Copyable<T>;	// Insert item at index
	bool insert(unsigned int index, T&& item) noexcept;			// Insert item at index (move operation)
	template <typename... Args>
	bool insertEmplace(unsigned int index, Args&&... args);		// Instantiate object in place
	unsigned int insertSorted(const T& item, bool order, bool dup) requires Copyable<T>;	// Insert item into sorted list
	unsigned int insertSorted(T&& item, bool order, bool dup) noexcept;	// Insert item into sorted list (move operation)

	int push(const T& item) requires Copyable<T> { return (insert(count, item) ? getCount() - 1 : -1); }	// Push item on end
	int push(T&& item) noexcept { return (insert(count, std::move(item)) ? getCount() - 1 : -1); }	// Push item on end (move operation)
	template <typename... Args>
	int pushEmplace(Args&&... args);							// Instantiate object in place at back
	bool pop() { return count ? remove(count - 1) : false; }	// Pop/remove end item (TODO: should this return a copy of the item?)
	bool swap(unsigned int a, unsigned int b) requires Copyable<T>;		// Swap two elements

	const T* getData()const { return data; }					// Get pointer to data
	void copyData(const ArrayList<T>& b) requires Copyable<T>;	// Copy data, allocate only if necessary
	void moveData(ArrayList<T>&& b) noexcept;					// Move data, allocate only if necessary
	
	bool remove(unsigned int index);							// Remove item at index

	void clear();												// Clear list
	void free();												// Clear and de-allocate memory
	
	// READ FUNCTIONS
	T& get(unsigned int index);									// Get reference to item at index
	T& getLast();												// Get reference to end item

	// STD::VECTOR FUNCTIONS (disable for non-copyable types)
	ArrayList(const std::vector<T>& copy) requires Copyable<T> { *this = copy; }
	const ArrayList<T>& operator =(const std::vector<T>& b) requires Copyable<T>;

private:
	void moveDataUp(unsigned int index);	// Shift data up above index
	void moveDataDown(unsigned int index, bool index_initialized = true);	// Shift data down above index
};


// Constructor: move (rvalue-only (&&))
template <typename T>
ArrayList<T>::ArrayList(ArrayList<T>&& other) noexcept
	: data(other.data), count(other.count), size(other.size)
{
	other.data = nullptr;
	other.count = 0;
	other.size = 0;
}

// Constructor: array
template <typename T>
ArrayList<T>::ArrayList(const T in[], unsigned int len) requires Copyable<T> : ArrayList()
{
	if (len == 0)
		return;

	free();

	updateCount(len);

	// TODO: Share this code with operator =
	if (!std::is_trivially_copyable<T>::value)
		for (unsigned int i = 0; i < len; i++)
			new (&data[i]) T(in[i]);
	else
		memcpy(data, in, len * sizeof(T));
}

// Operator =
template <typename T>
const ArrayList<T>& ArrayList<T>::operator =(const ArrayList<T>& b) requires Copyable<T>
{
	free();

	updateCount(b.getCount());

	if (std::is_trivially_copyable<T>::value)
		// Copy raw data
		memcpy(data, b.data, b.getCount() * sizeof(T));
	else
		// Call copy constructors
		for (unsigned int i = 0; i < b.getCount(); i++)
			new (&data[i]) T(b[i]);

	return (const ArrayList<T>&)(*this);
}

// Operator = (move operation with rvalue-only argument (&&))
template <typename T>
const ArrayList<T>& ArrayList<T>::operator =(ArrayList<T>&& other) noexcept
{
	if (this != &other)
	{
		free();
		data = other.data;
		count = other.count;
		size = other.size;

		other.data = nullptr;
		other.count = 0;
		other.size = 0;
	}

	return *this;
}


// Update list count
template <typename T>
bool ArrayList<T>::updateCount(unsigned int new_count)
{
	if (new_count > count)
	{
		allocate(new_count, false);	// Allocate more memory if needed
		
		// Call default constructor for new items
		// TODO: inefficient, but currently necessary to avoid calling destructors on uninitialized memory
		if (!std::is_trivially_copyable<T>::value)
		{
			for (unsigned int i = count; i < new_count; i++)
				new (&data[i]) T();	// Call default constructor
		}

		count = new_count;

		return true;
	}
	else
		return false;
}

// Re-allocate to new size
template <typename T>
bool ArrayList<T>::allocate(unsigned int new_size, bool downsize)
{
	// Note: be very careful with classes, structs, or pointers
	//	     the copy / delete operations below could cause pointer double-deletion, etc.

	if (new_size <= size && !downsize)
		return false;
	else
	{
		if (new_size < count)
			new_size = count;

		// Re-size
		if (new_size <= 10)
			new_size = 10;
		else if (new_size <= 100)
			new_size = 100;
		else if (new_size <= 1000)
			new_size = 1000;
		else if (new_size <= 10000)
			new_size = 10000;
		else
		{
			// TODO: doesn't scale well for very large arrays
			unsigned int remainder = new_size % 10000;

			if (remainder > 0)
				new_size += 10000 - remainder;
		}

		// Copy data to new memory location
		if (count > 0)
		{
			T* temp = data;
			// Manually allocate memory to avoid creating object instances
			data = static_cast<T*>(::operator new(new_size * sizeof(T)));

			// Copy data in used memory (up to count)
			memcpy(data, temp, count * sizeof(T));

			// Manually delete (does not call destructors)
			::operator delete(temp);
		}
		else
			data = static_cast<T*>(::operator new(new_size * sizeof(T)));

		size = new_size;

		return true;
	}
}

// Set value at index
template<typename T>
bool ArrayList<T>::set(unsigned int index, const T& item) requires Copyable<T>
{
	if (index >= count)
		return false;
	else
	{
		if (std::is_trivially_copyable<T>::value)
			data[index] = item;
		else
		{
			data[index].~T();
			new (&data[index]) T(item);
		}
			
		return true;
	}
}

// Set value at index (move operation)
template<typename T>
bool ArrayList<T>::set(unsigned int index, T&& item) noexcept
{
	if (index >= count)
		return false;
	else
	{
		//// TODO: is_trivially_copyable always false for move operation functions?
		//if (!std::is_trivially_copyable<T>::value)
		//	data[index].~T();

		data[index] = std::move(item);

		return true;
	}
}

// Emplace object at index (overwriting previous)
template <typename T>
template <typename... Args>
bool ArrayList<T>::setEmplace(unsigned int index, Args&&... args)
{
	if (index >= count)
		return false;
	else
	{
		if (!std::is_trivially_copyable<T>::value)
			data[index].~T();	// Call destructor

		// Instantiate object in place using placement new
		new (&data[index]) T(std::forward<Args>(args)...);

		return true;
	}
}

// Set value of all (integer)
template <typename T>
bool ArrayList<T>::setAllInt(int item)
{
	static_assert(std::is_trivially_copyable<T>::value,
		"setAllInt is only available for non-class types");

	if (count == 0)
		return false;
	else
	{
		memset(data, item, count * sizeof(T));
		return true;
	}
}

// Assign value of all (any data type)
template<typename T>
bool ArrayList<T>::setAll(const T& item) requires Copyable<T>
{
	if (count == 0)
		return false;
	else
	{
		for (unsigned int i = 0; i < count; i++)
			set(i, item);

		return true;
	}
}

// Assign value of all (move operation, any data type)
template<typename T>
bool ArrayList<T>::setAll(T&& item) noexcept
{
	if (count == 0)
		return false;
	else
	{
		for (unsigned int i = 0; i < count; i++)
			set(i, std::move(item));

		return true;
	}
}

// Insert item at index
// Push other items up the list
template<typename T>
bool ArrayList<T>::insert(unsigned int index, const T& item) requires Copyable<T>
{
	// Highest insert allowed is right after last item (index == count)
	if (index > count)
		return false;
	else if (index == count)
		addOneCount();
	else
		moveDataUp(index);

	return set(index, item);
}

// Insert item at index (move operation)
// Push other items up the list
template<typename T>
bool ArrayList<T>::insert(unsigned int index, T&& item) noexcept
{
	// Highest insert allowed is right after last item (index == count)
	if (index > count)
		return false;
	else if (index == count)
		addOneCount();
	else
		moveDataUp(index);

	return set(index, std::move(item));
}

// Instantiate object in place
template <typename T>
template <typename... Args>
bool ArrayList<T>::insertEmplace(unsigned int index, Args&&... args)
{
	// Highest insert allowed is right after last item (index == count)
	if (index > count)
		return false;
	else if (index == count)
		addOneCount();
	else
		moveDataUp(index);

	return setEmplace(index, std::forward<Args>(args)...);
}


// Insert item into sorted list
// List must be pre-sorted (TODO: add bool "sorted" variable)
template<typename T>
unsigned int ArrayList<T>::insertSorted(const T& item, bool order, bool dup) requires Copyable<T>
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

// Insert item into sorted list (move operation)
// List must be pre-sorted (TODO: add bool "sorted" variable)
template<typename T>
unsigned int ArrayList<T>::insertSorted(T&& item, bool order, bool dup) noexcept
{
	if (isEmpty())
		insert(0, std::move(item));
	else
	{
		// Check if should be inserted at beginning
		if ((dup && data[0] == item) || (!order && data[0] > item) || (order && data[0] < item))
		{
			insert(0, std::move(item));
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
					insert(i + 1, std::move(item));
					return i + 1;
				}
			}
			// Descending
			else
			{
				if (data[i + 1] < item && (data[i] > item || dup))
				{
					// Insert between i and i+1
					insert(i + 1, std::move(item));
					return i + 1;
				}
			}
		}

		// End of list
		unsigned int i = count - 1;
		if (dup || data[i] != item)
		{
			insert(count, std::move(item));
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


// Swap two elements
template <typename T>
bool ArrayList<T>::swap(unsigned int a, unsigned int b) requires Copyable<T>
{
	if (a >= getCount() || b >= getCount())
		return false;

	T temp = get(a);
	set(a, get(b));
	set(b, temp);

	return true;
}


// Empty array and copy data from another array, only allocate if needed
// TODO: how is this different from the assignment operator? Just doesn't re-allocate?
template <typename T>
void ArrayList<T>::copyData(const ArrayList<T>& b) requires Copyable<T>
{
	// Note: This can be dangerous if copying objects with pointers. 
	//	     Most common issues is double deletion.
	//static_assert(!std::is_class_v<T>, "ArrayList<T>::copyData disabled for class/struct types due to pointer duplication risk");

	if (b.isEmpty())
		return;
	
	clear();

	updateCount(b.getCount());

	if (std::is_trivially_copyable<T>::value)
		// Copy raw data
		memcpy(data, b.data, b.getCount() * sizeof(T));
	else
		// Call copy constructors
		for (unsigned int i = 0; i < b.getCount(); i++)
			new (&data[i]) T(b[i]);
}

// Remove item at index
template <typename T>
bool ArrayList<T>::remove(unsigned int index)
{
	if (index > count - 1)
		return false;

	if (!std::is_trivially_copyable<T>::value)
	{
		data[index].~T();
		memset(&data[index], 0, sizeof(T));
	}

	if (index >= count - 1)
		count--;
	else
		moveDataDown(index, false); // Second argument false becasue item at index has already been destructed

	return true;
}


// Clear list, but keep memory
template <typename T>
void ArrayList<T>::clear()
{
	// Call destructors if not basic type
	if (!std::is_trivially_copyable<T>::value)
		for (unsigned int i = 0; i < count; i++)
			data[i].~T();
		
	memset(data, 0, count * sizeof(T));

	count = 0;
}


// Clear and de-allocate memory
template <typename T>
void ArrayList<T>::free()
{
	clear();

	if (data)
		::operator delete(data);  // Raw deallocation only

	data = nullptr;
	count = 0;
	size = 0;
}


// Get reference to item at index
template <typename T>
T& ArrayList<T>::get(unsigned int index)
{
	if (index >= count)
		throw std::out_of_range("ArrayList<T>::get() - index out of range");
	return data[index];
}


// Get reference to end item
template <typename T>
T& ArrayList<T>::getLast()
{
	if (count == 0)
		throw std::out_of_range("ArrayList<T>::getLast() - list is empty");
	get(count - 1);
}


// std::vector functions
template<typename T>
const ArrayList<T>& ArrayList<T>::operator =(const std::vector<T>& b) requires Copyable<T>
{	
	free();

	updateCount(b.size());

	if (std::is_trivially_copyable<T>::value)
		// Copy raw data
		memcpy(data, b.data(), b.size() * sizeof(T));
	else
		// Call copy constructors
		for (unsigned int i = 0; i < b.size(); i++)
			new (&data[i]) T(b[i]);

	return *this;
}


// Shift data up above index
template <typename T>
void ArrayList<T>::moveDataUp(unsigned int index)
{
	if (index >= count)
		return;

	addOneCount();

	// Call destructor on last item (default initialized by addOneCount()) before overwriting memory
	// TODO: Inefficient and only necessary because addOneCount/updateCount currently has to initialize 
	//		 new memory to avoid calling destructors on uninitialized memory. Should re-architect
	if (!std::is_trivially_copyable<T>::value)
		data[count - 1].~T();

	// TODO: could make this more efficient by incorporating the data shift into the re-allocation step?
	// Move memory (even for classes)
	memmove(&data[index + 1], &data[index], (count - 1 - index) * sizeof(T));

	if (std::is_trivially_copyable<T>::value)
		memset(&data[index], 0, sizeof(T));
	else
	{
		// Call default constructor for new item
		// TODO: same as above, this is inefficient and only necessary to avoid calling destructors on uninitialized memory
		new (&data[index]) T();  // Default initialize
	}
}


// Shift data down from index + 1 to index
template <typename T>
void ArrayList<T>::moveDataDown(unsigned int index, bool index_initialized)
{
	if (index >= count - 1)
		return;

	// Call destructor on item at index before overwriting (if index_initialized is true)
	// TODO: May be inefficient in some cases, but calling function has control over this behavior
	if (index_initialized && !std::is_trivially_copyable<T>::value)
		data[count - 1].~T();

	// Move memory (even for classes)
	memmove(&data[index], &data[index + 1], (count - 1 - index) * sizeof(T));

	// Set memory to zero, no need to destruct because object moved
	memset(&data[count - 1], 0, sizeof(T));

	count--;  // Decrease count since we removed an item
}

#endif
