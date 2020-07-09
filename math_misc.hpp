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

// Determinant() - Calculate the determinant of an arbitrary-sized square matrix
// Inputs:  mat = matrix data array
//          size = matrix size (width of square matrix)
// Return:  determinant result
template <typename T>
T Determinant(const T *mat,unsigned int size)
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

    delete sub_mat;

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
			delete A; delete b; delete row_index;
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

	delete A; delete b; delete row_index;
	return true;
}


// TODO: UNTESTED
// SolveMLCP_PGS_Catto() - Solve an MLCP (Mixed Linear Complimentarity Problem) using Projected Gauss-
//	Seidel method, an extension of Gauss-Seidel. In a computer graphics application, the primary use
//	case for this method is to solve systems of inequalities related to constraints such as non-
//	penetrating collisions, joints, etc. This implementation has been written with that particular
//	case in mind and is not be a fully generalized MLCP / PGS algorithm (though functionally the same).
//	This algorithm assumes constraints only exist between pairs of bodies, so J and B are stored in a
//	"sparse" format with fixed number of columns based on body dimensionality (i.e. n_d). This requires
//	J_map to decode.
// References: https://www.toptal.com/game/video-game-physics-part-iii-constrained-rigid-body-simulation
//	Erin Catto GDC 2005 https://code.google.com/archive/p/box2d/downloads
// Format:	JBλ = η or JBx = z
// Note: All matrix arrays stored in column major order
// Inputs:  J_sp = [s by 2*n_d] 
//					Jacobian matrix of constraint functions (sparse)
//			J_map = [s by 2] 
//					Mapping of sparse jacobian matrix bodies (0 indexed)
//					Note: Each element(i,j) is the index of a rigid body. By convention, if a constraint
//					is between a single rigid body and ground, then (i,0) = 0 and the corresponding 
//					J(i,0) is zero.
//			B_sp = [2*n_d by s] 
//					B term: M(inv)*J(transpose) (sparse)
//					M is a mass and inertia diagonal matrix of n*n_d by n*n_d
//					Use J_map for indexing
//			z_vec (η) = [s]
//					z(η) term: (1/dt)*ζ - J((1/dt)*V1 + M(inv)*F_ext)
//					ζ is a constraint bias factor for stabilization
//					V1 is the initial (previous) velocity
//					F_ext is the external force acting on a body
//			x_min_vec (λ-) = [s] 
//					Minimum bound - use std inf (1.0 / 0.0) or max value to represent positive infinity
//			x_max_vec (λ+) = [s]
//					Maximum bound - use std -inf (-1.0 / 0.0) or lowest value to represent negative infinity
//			s = number of constraints
//			n = number of bodies
//			n_d = number of elements per body (i.e. x,y,theta; x,y,z,quat(a,b,c); etc)
//			max_iter = maximum number of iterations
//			x0_vec_in (λ0) [s]
//					Warm start - initial guess at λ, usually result from the previous frame
// Return:  x_vec (λ) = [s]
//					Solution of unknowns
//					Note: x_vec may be initialized with a best-guess (e.g. results from prev. frame) 
//					to warm-start the algorithm.
//          (return) boolean specifying whether the operation was successful or not
template <typename T>
bool SolveMLCP_PGS_Catto(const T *J_sp, const unsigned int *J_map, const T *B_sp, const T *z_vec,
					const T *x_min_vec, const T *x_max_vec, unsigned int s, unsigned int n, 
					unsigned int n_d, unsigned int max_iter, const T *x0_vec_in, T *x_vec)
{
	// TO-DO:
	// Add early termination checks - less than max error, etc
	// Check for valid pointers, ints, etc
	// Check max iterations for <= 0 and do something with it
	// Check for positive semi-definite (if necessary), i.e. the JB matrix

	// Work Variables
	T a* = new T[n * n_d];
	T d* = new T[s];
	T temp1* = new T[n_d];	// Temp array vector
	T temp3* = new T[n_d];	// Temp array vector
	T temp2* = new T[n_d];	// Temp array vector
	T temp4* = new T[n_d];	// Temp array vector
	unsigned int b1;		// Jmap object 1 index
	unsigned int b2;		// Jmap object 2 index
	T x_delta;				// Temporary x(λ) delta
	T a_b1* = new T[n_d];	// Temp array vector
	T a_b2* = new T[n_d];	// Temp array vector

	// Copy x0_vec_in
	T x0_vec* = new T[s];
	if (x0_vec_in)
		memcpy(x0_vec, x0_vec_in, s * sizeof(T));
	else
		memset(x0_vec, 0, s * sizeof(T));

	// Initialize / Warm Start x
	memcpy(x_vec, x0_vec, s * sizeof(T));

	// Initialize work variable a
	a = ArrayMatMult(B, x);

	// Initialize work variable d
	// d(i) = dot(Jsp(i,0),Bsp(0,i)) + dot(Jsp(i,1),Bsp(1,i))
	for (int i = 0; i < s; i++)
	{
		for (int j = 0; j < n_d; j++)
		{
			// Get J_sp(i,0) (i.e. first n_d-length half of row)
			temp1[j] = J_sp[j * s + i];
			// Get J_sp(i,1) (i.e. second n_d-length half of row)
			temp3[j] = J_sp[(j + n_d) * s + i];
			// Get B_sp(0,i) (i.e. first n_d-length half of col)
			temp2[j] = B_sp[i * (2 * n_d) + j];
			// Get B_sp(1,i) (i.e. second n_d-length half of col)
			temp4[j] = B_sp[i * (2 * n_d) + j + n_d];
		}
		d[i] = ArrayVecDot(temp1, temp2) + ArrayVecDot(temp3, temp4);
	}

	// Main Iterative Loop
	for (int iter = 1; iter <= max_iter; iter++)
	{
		for (int i = 0; i < s; i++)
		{
			for (int j = 0; j < n_d; j++)
			{
				// Get J_sp(i,0) (i.e. first n_d-length half of row)
				temp1[j] = J_sp[j * s + i];
				// Get J_sp(i,1) (i.e. second n_d-length half of row)
				temp3[j] = J_sp[(j + n_d) * s + i];
				// Get B_sp(0,i) (i.e. first n_d-length half of col)
				temp2[j] = B_sp[i * (2 * n_d) + j];
				// Get B_sp(1,i) (i.e. second n_d-length half of col)
				temp4[j] = B_sp[i * (2 * n_d) + j + n_d];

				// Get a(b1)
				a_b1[j] = a[b1 * n_d + j];
				// Get a(b2)
				a_b2[j] = a[b2 * n_d + j];
			}

			b1 = J_map[i];		// Object 1 index
			b2 = J_map[s + i];	// Object 2 index

			// Tentative x_delta
			x_delta = (z_vec[i] - ArrayVecDot(temp1, a_b1) - ArrayVecDot(temp3, a_b2)) / d[i];

			// Test against constraints
			x0_vec[i] = x_vec[i];	// Save last iteration result (same on first loop)
			x_vec[i] = max(x_min_vec, min(x0_vec[i] + x_delta, x_max_vec));	// Check min/max
			x_delta = x_vec[i] - x0_vec[i];	// Final delta

			// Update work variable a
			for (int j = 0; j < n_d; j++)
			{
				a[b1 * n_d + j] = a_b1[j] + x_delta[i] * temp2[j];
				a[b2 * n_d + j] = a_b2[j] + x_delta[i] * temp4[j];
			}
		}
	}

	return true;	// TODO - check for issues
}

#endif
