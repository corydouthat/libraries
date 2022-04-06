// *****************************************************************************************************************************
// math_misc.hpp
// Miscellaneous Math Functions
// Author(s): Cory Douthat
// Copyright (c) 2020 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef MATH_MISC_HPP_
#define MATH_MISC_HPP_

#include <stdlib.h>
#define _USE_MATH_DEFINES	//For PI definition
#include <cmath>
#include <limits>

#include "vec.hpp"
#include "mat.hpp"
#include "array_mat.hpp"


// SolveQuadratic()
// Solve quadratic equation using numerically stable method
// q = -1/2(b + sign(b) * sqrt(b^2 - 4ac))
// x1 = q/a		x2 = c/q
// Inputs:	a, b, c - coefficients from ax^2 + bx + c = 0
// Return:	x1, x2 - roots
//			bool - real roots found
template <typename T>
bool SolveQuadratic(T a, T b, T c, T* x1, T* x2)
{
	if (a == 0)		// TODO: check close to zero?
		return false;
	
	T q;
		
	if (b >= 0)
		q = -0.5 * (b + sqrt(b * b - 4 * a * c));
	else
		q = -0.5 * (b - sqrt(b * b - 4 * a * c));

	if (q == 0)		// TODO: check close to zero?
		return false;

	t1* = q / a;
	t2* = c / q;

	return true;
}


// Determinant() - Calculate the determinant of an arbitrary-sized square matrix
// Inputs:  mat = matrix data array
//          size = matrix size (width of square matrix)
// Return:  determinant result
template <typename T>
T Determinant(const T *mat, unsigned int size)
{
    if (size < 2)
        return (T)0;

    // 2x2
    if (size == 2)
        return Mat2<T>(mat).det();

    // 3x3
    if (size == 3)
        return Mat3<T>(mat).det();

    // 4x4
    if (size == 4)
        return Mat4<T>(mat).det();

    // 5x5+
    T det = (T)0;
    T *sub_mat;
    sub_mat = new T[(size - 1)*(size - 1)];
    for (unsigned int i = 0; i < size; i++)
    {
        int col = 0;
        for (unsigned int col_maj = 0; col_maj < size; col_maj++)
        {
            if (col_maj != i)
            {
                for (unsigned int row = 0; row < size - 1; row++)
                {
                    sub_mat[col*(size - 1) + row] = mat[col_maj*(size) + row + 1];
                }
                col++;
            }
        }
        if (i % 2 == 0)
            det += mat[i*size] * Determinant(sub_mat,size - 1);
        else
            det -= mat[i*size] * Determinant(sub_mat,size - 1);
    }

    delete[] sub_mat;

    return det;
}

// CAUTION: Numerically unstable. Use Gaussian Elimination instead unless you know what you're doing.
// SolveCramer() - Solve a system of linear equations using Cramer's Rule.
// Format Ax = b
// Inputs:  A_mat = n x n matrix (array) defining the coefficients of the equations
//          b_vec = n size vector (array) defining right-hand side of each linear equation
//          n = number of equations/variables
// Return:  x_vec = n size vector (array) representing solutions - passed by pointer
//          return boolean specifying whether the operation was successful or not
template <typename T>
bool SolveCramer(const T *A_mat,const T *b_vec,unsigned int n,T *x_vec)
{
    T det_A = Determinant(A_mat,n);
    if (det_A == T(0))
        return false;

    // Matrix for temporary A with column replaced
    T *Ai = new T[n * n];
    memcpy(Ai,A_mat,n * n * sizeof(T));

    for (unsigned int i = 0; i < n; i++)
    {
        // Copy in b vector to column i
        memcpy(&Ai[i * n],b_vec,n * sizeof(T));

        // Solve for xi
        x_vec[i] = Determinant(Ai,n) / det_A;

        // Copy regular values back into column i
        memcpy(&Ai[i * n],&A_mat[i * n],n * sizeof(T));
    }

	delete Ai;
    return true;
}


// TODO: Improve efficiency by checking for special cases or no solution earlier?
// SolveGaussElim() - Solve a system of linear equations using Gaussian Elimination
// Performs Partial Pivoting; does not perform Full Pivoting
// Format: Ax = b
// Inputs:  A_mat = n x n matrix (array) defining the coefficients of the equations
//          b_vec = n size vector (array) defining right-hand side of each linear equation
//          n = number of equations/variables
// Return:  x_vec = n size vector (array) representing solutions - passed by pointer
//          return boolean specifying whether the operation was successful or not
template <typename T>
bool SolveGaussElim(const T *A_mat, const T *b_vec, unsigned int n, T *x_vec)
{
	T *A = new T[n*n];
	T *b = new T[n];
	memcpy(A, A_mat, sizeof(T)*n*n);
	memcpy(b, b_vec, sizeof(T)*n);

	unsigned int *row_index = new unsigned int[n];	// Vector to track row order
	for (unsigned int i = 0; i < n; i++) row_index[i] = i;

	// Reduce matrix to row echelon form using Gaussian Elimination and partial pivoting
	for (unsigned int r = 0; r < n; r++)
	{
		// Pivot row 'r' (check for rows with larger absolute values in the pivot column)
		if (r < n - 1)	// Do not do this when on last row
		{
			for (unsigned int i = r + 1; i < n; i++)
			{
				// Check each successive row to see if absolute value in 'r'th column is larger
				if (abs(A[r*n + row_index[i]]) > abs(A[r*n + row_index[r]]))
				{
					// Swap rows
					unsigned int temp = row_index[r];
					row_index[r] = row_index[i];
					row_index[i] = temp;
				}
			}
		}

		// Check for zero pivot (no solution)
		if (A[r*n + row_index[r]] == (T)0)
		{
			delete[] A; delete[] b; delete[] row_index;
			return false;
		}

		// Reduce successive rows (eliminate columns below pivot of 'r')
		if (r < n - 1)	// Do not do this when on last row
		{
			for (unsigned int i = r + 1; i < n; i++)	// ROWS
			{
				// Multiplication factor
				T temp_mult = A[r*n + row_index[i]] / A[r*n + row_index[r]];
				// First column set 0
				A[r*n + row_index[i]] = (T)0;
				// Addl columns
				for (unsigned int j = r + 1; j < n; j++)
					A[j*n + row_index[i]] -= temp_mult * A[j*n + row_index[r]];
				// b vector
				b[row_index[i]] -= temp_mult * b[row_index[r]];
			}
		}
	}

	// Back-substitue for x values
	for (int x = n - 1; x >= 0; x--)
	{
		// Set 'x' to right side of equation
		x_vec[x] = b[row_index[x]];
		// Subtract known 'x' values from right side of equation
		for (int j = n - 1; j > x; j--)
		{
			x_vec[x] -= x_vec[j] * A[j*n + row_index[x]];
		}
		// Solve for 'x'
		x_vec[x] /= A[x*n + row_index[x]];
	}

	delete[] A; delete[] b; delete[] row_index;
	return true;
}

#endif
