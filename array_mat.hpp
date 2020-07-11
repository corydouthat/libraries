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

template <typename T = float>
class ArrayMat
{
    // DATA
    T *v;
    unsigned int rows, cols;
public:
    // FUNCTIONS
    // Constructors and Destructor
    ArrayMat() : v(nullptr), rows(0), cols(0) {}			// Constructor: null pointer
    ~ArrayMat() { clear(); }

    // Basic Operators
    const ArrayMat<T>& operator =(const ArrayMat<T>& b);    // Operator = (assignment)
    // Matrix Operators
    ArrayMat<T> operator +(const ArrayMat<T>& b)const;	    // Operator +
    ArrayMat<T> operator -(const ArrayMat<T>& b)const;	    // Operator -
    ArrayMat<T> operator *(const ArrayMat<T>& b)const;		// Operator *
    // Scalar Operators
    ArrayMat<T> operator *(T s)const;						// Operator * (scalar)
    ArrayMat<T> operator -()const { return *this * -1; }	// Operator - (scalar inverse)
    const ArrayMat<T>& operator *=(T s);                    // Operator *= (scalar)
    // Other Member Functions
    bool isEmpty()const;                                    // Check if empty (null pointer)
    void clear();                                           // Erase and set to null/zero
    void allocate(unsigned int r, unsigned int c);          // Allocate memory
    void allocateZero(unsigned int r, unsigned int c);      // Allocate memory (to zeros)
    T get(unsigned int col, unsigned int row)const;		    // Get element
    void set(unsigned int col, unsigned int row, T s);	    // Set element
    //T det()const;                                         // Determinant (square matrix)
    ArrayMat<T> transp()const;								// Transpose
    //ArrayMat<T> inv()const;							    // Inverse (square matrix)
    ArrayMat<T> inv_diag()const;                            // Inverse (diagonal matrix)
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
ArrayMat<T> ArrayMat<T>::operator +(const ArrayMat<T>& b)const
{
    ArrayMat<T> y;
    unsigned int min_r, min_c;

    // Check for empty b matrix
    if (b.isEmpty())
        return *this;

    // Check for empty this matrix
    if (isEmpty())
        return b;

    y.allocate(rows, cols);
    min_r = min(rows, b.rows);
    min_c = min(cols, b.cols);

    for (unsigned int j = 0; j < cols; j++)
    {
        for (unsigned int i = 0; i < rows; i++)
        {
            if (i < min_rows && j < min_cols)
                y.v[j][i] = v[j][i] + b.v[j][i];
            else
                y.v[j][i] = v[j][i];
        }
    }

    return y;
}

// - operator
// Will attempt to subtract matrices even if dimensions don't match
// Output will match dimensions of first term (a)
template <typename T>
ArrayMat<T> ArrayMat<T>::operator -(const ArrayMat<T>& b)const
{
    ArrayMat<T> y;
    unsigned int min_r, min_c;

    // Check for empty b matrix
    if (b.isEmpty())
        return *this;

    // Check for empty this matrix
    if (isEmpty())
        return -b;

    y.allocate(rows, cols);
    min_r = min(rows, b.rows);
    min_c = min(cols, b.cols);

    for (unsigned int j = 0; j < cols; j++)
    {
        for (unsigned int i = 0; i < rows; i++)
        {
            if (i < min_rows && j < min_cols)
                y.v[j][i] = v[j][i] - b.v[j][i];
            else
                y.v[j][i] = v[j][i];
        }
    }

    return y;
}

// operator *
// Only allowed when number of columns in a equals rows in b
// Will return empty matrix if these conditions aren't met
template <typename T>
ArrayMat<T> ArrayMat<T>::operator *(const ArrayMat<T>& b)const
{
    ArrayMat<T> y;

    // Check conditions
    if (isEmpty() || b.isEmpty() || cols != b.rows)
        return y;   // No operation, empty matrix

    y.allocateZero(rows, b.cols);

    for (unsigned int j = 0; j < y.cols; j++)			//j = cols of y/b
    {
        for (unsigned int i = 0; i < y.rows; i++)		//i = rows of y/a
        {
            for (unsigned int k = 0; k < y.cols; k++)	//k = sum index
            {
                y.v[j][i] += v[k][i] * b.v[j][k];
            }
        }
    }

    return y;
}


// SCALAR OPERATORS
// *****************************************************************************************************************************

template <typename T>
ArrayMat<T> ArrayMat<T>::operator *(T s)const
{
    ArrayMat<T> y;
    
    if (isEmpty())
        return *this;

    y.allocate(rows, cols);

    for (unsigned int j = 0; j < cols; j++)
    {
        for (unsigned int i = 0; i < rows; i++)
        {
            y.v[j][i] = v[j][i] * s;
        }
    }

    return y;
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
            v[j][i] *= s;
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
        delete v;

    v = nullptr;

    rows = cols = 0;
}


// Allocate memory
template <typename T>
void ArrayMat<T>::allocate(unsigned int r, unsigned int c)
{
    clear();

    v = new T[cols][rows];

    rows = r;
    cols = c;
}

// Allocate memory (to zeros)
template <typename T>
void ArrayMat<T>::allocateZero(unsigned int rows, unsigned int cols)
{
    clear();

    v = new T[cols][rows];
    memset(v, 0, rows * cols * sizeof(T));

    rows = r;
    cols = c;
}

template <typename T>
T ArrayMat<T>::get(unsigned int col, unsigned int row)const
{
    return v[col][row];
}

template <typename T>
void ArrayMat<T>::set(unsigned int col, unsigned int row, T s)
{
    v[col][row] = s;
}

// Transpose
template <typename T>
ArrayMat<T> ArrayMat<T>::transp()const
{
    ArrayMat<T> y;
    
    if (isEmpty())
        return y;

    y.allocate(cols, rows);

    for (unsigned int j = 0; j < cols; j++)
    {
        for (unsigned int i = 0; i < rows; i++)
        {
            y.[i][j] = v[j][i];
        }
    }
}

// Inverse (diagonal matrix only)
template <typename T>
ArrayMat<T> ArrayMat<T>::inv_diag()const
{
    ArrayMat<T> y;

    // Check for empty or non-square matrix
    if (isEmpty() || rows != cols)
        return y;

    y.allocateZero(rows, cols);

    for (unsigned int i = 0; i < rows; i++)
    {
        y.v[i][i] = T(1) / v[i][i];
    }
}

#endif
