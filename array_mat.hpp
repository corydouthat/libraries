// *****************************************************************************************************************************
// array_mat.hpp
// Linear Math Libraries
// Arbitrary-sized matrix class
// The matrix data is stored in column-major format. [j] returns column j,
// and [j][i] returns the item at column j, row i.
// Author(s): Cory Douthat
// Copyright (c) 2020 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef ARRAY_MAT_HPP_
#define ARRAY_MAT_HPP_

#include <stdlib.h>
#define _USE_MATH_DEFINES	//For PI definition
#include <cmath>

#include "mat.hpp"

template <class T> 
const T& min(const T& a, const T& b)
{
    return !(b < a) ? a : b;
}

template<class T>
const T& max(const T& a, const T& b)
{
    return (a < b) ? b : a;
}

template <typename T>
class ArrayMat
{
    // DATA
    T *v;
    unsigned int rows, cols;
public:
    // FUNCTIONS
    // Constructors and Destructor
    ArrayMat() : v(nullptr), rows(0), cols(0) {}
    ArrayMat(const ArrayMat<T>& b) { *this = b; }
    ArrayMat(unsigned int r, unsigned int c) { allocate(r, c); }
    ArrayMat(unsigned int r, unsigned int c, T s) { allocateValue(r, c, s); }
    ~ArrayMat() { clear(); }

    // Basic Operators
    const ArrayMat<T>& operator =(const ArrayMat<T>& b);    // Operator = (assignment)
    // Matrix Operators
    const ArrayMat<T>& operator +(const ArrayMat<T>& b)const;	// Operator +
    const ArrayMat<T>& operator -(const ArrayMat<T>& b)const;	// Operator -
    const ArrayMat<T>& operator *(const ArrayMat<T>& b)const;	// Operator *
    // Scalar Operators
    const ArrayMat<T>& operator *(T s)const;					// Operator * (scalar)
    const ArrayMat<T>& operator -()const { return *this * -1; }	// Operator - (scalar inverse)
    const ArrayMat<T>& operator *=(T s);                    // Operator *= (scalar)
    // Other Member Functions
    bool isEmpty()const;                                    // Check if empty (null pointer)
    void clear();                                           // Erase and set to null/zero
    void allocate(unsigned int r, unsigned int c);          // Allocate memory
    void allocateZero(unsigned int r, unsigned int c);      // Allocate memory (to zeros)
    void allocateValue(unsigned int r, unsigned int c, T s);// Allocate memory (to s value)
    T get(unsigned int col, unsigned int row)const;		    // Get element
    void set(unsigned int col, unsigned int row, T s);	    // Set element
    const T* getData()const { return v; };                   // Return const pointer to data
    unsigned int getCount()const { return rows * cols; }    // Return number of elements
    //T det()const;                                         // Determinant (square matrix)
    const ArrayMat<T>& transp()const;						// Transpose
    //ArrayMat<T> inv()const;							    // Inverse (square matrix)
    const ArrayMat<T>& invDiag()const;                      // Inverse (diagonal matrix)

    // Friend Functions
    template <typename sT> friend const ArrayMat<sT>& operator *(sT a, const ArrayMat<sT>& b);

    // Static Member Functions
    static T dotProduct(const ArrayMat<T>& a, const ArrayMat<T>& b);

    // Discrete Matrix (Mat2, Mat3, Mat4) Functions
#ifdef MAT_HPP_
    // Operators
    bool operator==(const Mat2<T> b)const;
    bool operator==(const Mat3<T> b)const;
    bool operator==(const Mat4<T> b)const;
    bool operator!=(const Mat2<T> b)const { return !(*this == b); }
    bool operator!=(const Mat3<T> b)const { return !(*this == b); }
    bool operator!=(const Mat4<T> b)const { return !(*this == b); }
#endif
};


// BASIC OPERATORS
// *****************************************************************************************************************************

template <typename T>
const ArrayMat<T>& ArrayMat<T>::operator =(const ArrayMat<T>& b)
{
    // Check for empty b matrix
    if (b.isEmpty() || !(b.v))
    {
        clear();
    }
    else
    {
        // Check for non-null pointer or different dimensions
        if (!v || rows != b.rows || cols != b.cols)
        {
            clear();
            allocate(b.rows, b.cols);
        }

        memcpy(v, b.v, rows * cols * sizeof(T));
    }

    return *this;
}


// MATRIX OPERATORS
// *****************************************************************************************************************************

// + operator
// Will attempt to add matrices even if dimensions don't match
// Output will match dimensions of first addend (a)
template <typename T>
const ArrayMat<T>& ArrayMat<T>::operator +(const ArrayMat<T>& b)const
{
    ArrayMat<T> *y = new ArrayMat<T>;
    unsigned int min_r, min_c;

    // Check for empty b matrix
    if (b.isEmpty())
        return *this;

    // Check for empty this matrix
    if (isEmpty())
        return b;

    y->allocate(rows, cols);
    min_r = min(rows, b.rows);
    min_c = min(cols, b.cols);

    for (unsigned int j = 0; j < cols; j++)
    {
        for (unsigned int i = 0; i < rows; i++)
        {
            if (i < min_r && j < min_c)
                y->set(j, i, get(j, i) + b.get(j, i));
            else
                y->set(j, i, get(j, i));
        }
    }

    return *y;
}

// - operator
// Will attempt to subtract matrices even if dimensions don't match
// Output will match dimensions of first term (a)
template <typename T>
const ArrayMat<T>& ArrayMat<T>::operator -(const ArrayMat<T>& b)const
{
    ArrayMat<T> *y = new ArrayMat<T>;
    unsigned int min_r, min_c;

    // Check for empty b matrix
    if (b.isEmpty())
        return *this;

    // Check for empty this matrix
    if (isEmpty())
        return -b;

    y->allocate(rows, cols);
    min_r = min(rows, b.rows);
    min_c = min(cols, b.cols);

    for (unsigned int j = 0; j < cols; j++)
    {
        for (unsigned int i = 0; i < rows; i++)
        {
            if (i < min_r && j < min_c)
                y->set(j, i, get(j, i) - b.get(j, i));
            else
                y->set(j, i, get(j, i));
        }
    }

    return *y;
}

// operator *
// Only allowed when number of columns in a equals rows in b
// Will return empty matrix if these conditions aren't met
template <typename T>
const ArrayMat<T>& ArrayMat<T>::operator *(const ArrayMat<T>& b)const
{
    ArrayMat<T> *y = new ArrayMat<T>;

    // Check conditions
    if (isEmpty() || b.isEmpty() || cols != b.rows)
        return *y;   // No operation, empty matrix

    y->allocateZero(rows, b.cols);

    for (unsigned int j = 0; j < y->cols; j++)			//j = cols of y/b
    {
        for (unsigned int i = 0; i < y->rows; i++)		//i = rows of y/a
        {
            for (unsigned int k = 0; k < b.rows; k++)	//k = sum index
            {
                y->set(j, i, y->get(j, i) + get(k, i) * b.get(j, k));
            }
        }
    }

    return *y;
}


// SCALAR OPERATORS
// *****************************************************************************************************************************

template <typename T>
const ArrayMat<T>& ArrayMat<T>::operator *(T s)const
{
    ArrayMat<T> *y = new ArrayMat<T>;
    
    if (isEmpty())
        return *this;

    y->allocate(rows, cols);

    for (unsigned int j = 0; j < cols; j++)
    {
        for (unsigned int i = 0; i < rows; i++)
        {
            y->set(j, i, get(j, i) * s);
        }
    }

    return *y;
}

template <typename T>
const ArrayMat<T>& ArrayMat<T>::operator *=(T s)
{
    if (isEmpty())
        return *this;

    for (unsigned int j = 0; j < cols; j++)
    {
        for (unsigned int i = 0; i < rows; i++)
        {
            set(j, i, get(j, i) * s);
        }
    }

    return *this;
}


// OTHER MEMBER FUNCTIONS
// *****************************************************************************************************************************

// Check if empty (null pointer)
template <typename T>
bool ArrayMat<T>::isEmpty()const
{
    if (rows == 0 || cols == 0)
        return true;
    else
        return false;
}

// Erase and set to null/zero
template <typename T>
void ArrayMat<T>::clear()
{
    if (v)
        delete[] v;

    v = nullptr;

    rows = cols = 0;
}


// Allocate memory
template <typename T>
void ArrayMat<T>::allocate(unsigned int r, unsigned int c)
{
    clear();

    rows = r;
    cols = c;
    
    v = new T[cols * rows];  
}

// Allocate memory (to zeros)
template <typename T>
void ArrayMat<T>::allocateZero(unsigned int r, unsigned int c)
{
    allocate(r, c);
    memset(v, 0, rows * cols * sizeof(T));
}

// Allocate memory (to s value)
template <typename T>
void ArrayMat<T>::allocateValue(unsigned int r, unsigned int c, T s)
{
    allocate(r, c);
    memset(v, s, rows * cols * sizeof(T));
}

template <typename T>
T ArrayMat<T>::get(unsigned int col, unsigned int row)const
{
    if (!isEmpty() && col <= cols && row <= rows)
        return v[col * rows + row];
    else
        return T(NAN);
}

template <typename T>
void ArrayMat<T>::set(unsigned int col, unsigned int row, T s)
{
    
    if (!isEmpty() && col <= cols && row <= rows)
        v[col * rows + row] = s;
}

// Transpose
template <typename T>
const ArrayMat<T>& ArrayMat<T>::transp()const
{
    ArrayMat<T> *y = new ArrayMat<T>;
    
    if (isEmpty())
        return *y;

    y->allocate(cols, rows);

    for (unsigned int j = 0; j < cols; j++)
    {
        for (unsigned int i = 0; i < rows; i++)
        {
            y->set(i, j, get(j, i));
        }
    }

    return *y;
}

// Inverse (diagonal matrix only)
template <typename T>
const ArrayMat<T>& ArrayMat<T>::invDiag()const
{
    ArrayMat<T> *y = new ArrayMat<T>;

    // Check for empty or non-square matrix
    if (isEmpty() || rows != cols)
        return *y;

    y->allocateZero(rows, cols);

    for (unsigned int i = 0; i < rows; i++)
    {
        y->set(i, i, T(1) / get(i, i));
    }

    return *y;
}


// FRIEND FUNCTIONS
// *****************************************************************************************************************************
template <typename sT> 
const ArrayMat<sT>& operator *(sT a, const ArrayMat<sT>& b)
{
    return b * a;
}


// STATIC MEMBER FUNCTIONS
// *****************************************************************************************************************************
template <typename T>
T ArrayMat<T>::dotProduct(const ArrayMat<T>& a, const ArrayMat<T>& b)
{
    T temp = 0;
    
    if (a.cols == b.cols == 1 && a.rows == b.rows)
    {
        for (unsigned int i = 0; i < a.rows; i++)
        {
            temp += a.get(0, i) * b.get(0, i);
        }

        return temp;
    }
    else
        return NAN;
}


// DISCRETE MATRIX (MAT2, MAT3, MAT4) FUNCTIONS
// *****************************************************************************************************************************
#ifdef MAT_HPP_
// Operators

template <typename T>
bool ArrayMat<T>::operator==(const Mat2<T> b)const
{
    if (r != 2 || c != 2)
        return false;
    else if (memcmp(v, b.getData(), 2 * 2 * sizeof(T)))
        return false;
    else
        return true;
}

template <typename T>
bool ArrayMat<T>::operator==(const Mat3<T> b)const
{
    if (r != 3 || c != 3)
        return false;
    else if (memcmp(v, b.getData(), 3 * 3 * sizeof(T)))
        return false;
    else
        return true;
}

template <typename T>
bool ArrayMat<T>::operator==(const Mat4<T> b)const
{
    if (rows != 4 || cols != 4)
        return false;
    else if (memcmp(v, b.getData(), 4 * 4 * sizeof(T)))
        return false;
    else
        return true;
}
#endif

#endif
