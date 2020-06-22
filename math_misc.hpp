// *****************************************************************************************************************************
// math_misc.hpp
// Miscellaneous Math Functions
// Author(s): Cory Douthat
// Copyright (c) 2017 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef MATH_MISC_HPP_
#define MATH_MISC_HPP_

#include <stdlib.h>
#define _USE_MATH_DEFINES	//For PI definition
#include <cmath>

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

// TODO/NOTE: It looks like this may be a failed experiment. Does not provide correct results when clamping. Problem form may simply be invalid.
// TODO: UNTESTED
// TODO: Add termination criteria (i.e. residual, etc)
// TODO: Handle default min/max for x or individual +/- infinity values
// TODO: TO-DO add checks for convergence issues??
// SolveGaussSeidelClamped() - Solve an MLCP (Mixed Linear Complimentarity Problem) using Gauss-Seidel
//	method with simple clamping constraints on unknowns. This may be equivalent to 
//	Projected Gauss-Seidel method, but it's unclear. No optimization is attempted to
//	target this to rigid body constraints specifically, or take advantage of matrix
//	sparsity by reducing dimensions (as with Erin Catto GDC 2005 presentation)
// References: https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
// Format: Ax = b
// Inputs:  A_mat = [n by n], col-major order, defines the coefficients of the equations
//          b_vec = [n] defines right-hand side of each linear equation
//			x_min = [n] minimum bound of x, null ptr means -infinity (for whole vector)
//			x_max = [n] maximum bound of x, null ptr means +infinity (for whole vector)
//          n = number of equations/variables
//			max_iter = maximum number of iterations
// Return:  x_vec = [n] represents solution variables - passed by pointer
//				Note: x_vec may be initialized with a best-guess (e.g. results from prev.
//				frame) to warm-start the algorithm.
//          return boolean specifying whether the operation was successful or not
template <typename T>
bool SolveGaussSeidelClamped(const T* A_mat, const T* b_vec, const T* x_min, const T* x_max, unsigned int n, unsigned int max_iter, T* x_vec)
{
	// Basic checks
	if (n <= 0)
		return false;

	if (x_min && x_max)
	{
		for (unsigned int c = 0; c < n; c++)
		{
			if (x_min[c] > x_max[c])
				return false;
		}
	}

	T sum1, sum2;

	for (unsigned int iter = 0; iter < max_iter; iter++)
	{
		for (int i = 0; i < n; i++)
		{
			// Calculate next iteration of x (i.e. k+1) using forward-substitution

			// First summation term:
			sum1 = 0;
			for (int j = 0; j < i; j++)
			{
				sum1 += A_mat[i * n + j] * x_vec[j];
			}

			// Second summation term:
			sum2 = 0;
			for (int j = i + 1; j < n; j++)
			{
				sum2 += A_mat[i * n + j] * x_vec[j];
			}

			x_vec[i] = b_vec[i] - sum1 - sum2;

			// Clamp value
			if (x_min && x_vec[i] < x_min[i])
				x_vec[i] = x_min[i];
			else if (x_max && x_vec[i] > x_max[i])
				x_vec[i] = x_max[i];
		}

		// Check termination criteria
		// TO-DO
	}

	return true;
}

// TODO: UNTESTED
// TODO: Check for positive semi-definite (if necessary), i.e. the JB matrix
// TODO: Add alternate methods for termination / error checking?
// TODO: How to handle infinity for min/max?
// SolvePGSCatto() - Solve an MLCP (Mixed Linear Complimentarity Problem) using Projected Gauss-
//	Seidel method, an extension of Gauss-Seidel. In a computer graphics application, the primary use
//	case for this method is to solve systems of inequalities related to constraints such as non-
//	penetrating collisions, joints, etc. This implementation has been written with that particular
//	case in mind and is not be a fully generalized MLCP / PGS algorithm (though functionally the same).
//	This algorithm assumes constraints only exist between pairs of bodies, so J and B are stored in a
//	"sparse" format with fixed number of columns based on body dimensionality (i.e. n_d). This requires
//	J_map to decode.
// References: https://www.toptal.com/game/video-game-physics-part-iii-constrained-rigid-body-simulation
//	Erin Catto GDC 2005 https://code.google.com/archive/p/box2d/downloads
// Format:	JBλ = η or J_mat*B_mat*x_vec = z_vec
// Inputs:  J_sp = [s by 2*n_d] jacobian matrix of constraint functions (sparse)
//			j_map = [s by 2] mapping of sparse jacobian matrix bodies
//					Note: Each element(i,j) is the index of a rigid body. By convention, if a constraint
//					is between a single rigid body and ground, then (i,1) = 0 and the corresponding 
//					J(i,1) is zero.
//			B_mat = [n*n_d by s] M(inv)*J(transpose)
//					M is a mass and inertia diagonal matrix of n*n_d by n*n_d
//			z_vec (η) = [n*n_d] (1/dt)*ζ - J((1/dt)*V1 + M(inv)*F_ext)
//					ζ is a constraint bias factor for stabilization
//					V1 is the initial (previous) velocity
//					F_ext is the external force acting on a body
//			x_min (λ-) = minimum bound
//			x_max (λ+) = maximum bound
//			s = number of constraints
//			n = number of bodies
//			n_d = number of elements per body (i.e. x,y,theta; x,y,z,quat(a,b,c); etc)
// Return:  x_vec (λ) = [s] solution of unknowns
//				Note: x_vec may be initialized with a best-guess (e.g. results from prev.
//				frame) to warm-start the algorithm.
//          return boolean specifying whether the operation was successful or not
template <typename T>
bool SolvePGSCatto()
{
	//TO-DO
	return false;	// TEMP
}

#endif
