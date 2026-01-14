// *****************************************************************************************************************************
// math_misc.hpp
// Miscellaneous Math Functions
// Author(s): Cory Douthat
// Copyright (c) 2025 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef MATH_MISC_HPP_
#define MATH_MISC_HPP_

#include <stdlib.h>
#define _USE_MATH_DEFINES	//For PI definition
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <unordered_set>
#include <algorithm>

#include "vec.h"
#include "mat.h"
#include "array_mat.hpp"

// TODO: replace with C++ 17 std::clamp
template<typename T>
T clamp(T val, T min, T max) 
{
	return std::max(std::min(val, max), min);
}


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
		q = (T)(-0.5 * (b + sqrt(b * b - 4 * a * c)));
	else
		q = (T)(-0.5 * (b - sqrt(b * b - 4 * a * c)));

	if (q == 0)		// TODO: check close to zero?
		return false;

	*x1 = q / a;
	*x2 = c / q;

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

// TODO: Move to phSolvers.hpp?
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

    // Matrix for temporary A with columin replaced
    T *Ai = new T[n * n];
    memcpy(Ai,A_mat,n * n * sizeof(T));

    for (unsigned int i = 0; i < n; i++)
    {
        // Copy in b vector to columin i
        memcpy(&Ai[i * n],b_vec,n * sizeof(T));

        // Solve for xi
        x_vec[i] = Determinant(Ai,n) / det_A;

        // Copy regular values back into columin i
        memcpy(&Ai[i * n],&A_mat[i * n],n * sizeof(T));
    }

	delete Ai;
    return true;
}


// TODO: Improve efficiency by checking for special cases or no solution earlier?
// TODO: Move to phSolvers.hpp?
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
		// Pivot row 'r' (check for rows with larger absolute values in the pivot columin)
		if (r < n - 1)	// Do not do this when on last row
		{
			for (unsigned int i = r + 1; i < n; i++)
			{
				// Check each successive row to see if absolute value in 'r'th columin is larger
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

		// Reduce successive rows (eliminate columins below pivot of 'r')
		if (r < n - 1)	// Do not do this when on last row
		{
			for (unsigned int i = r + 1; i < n; i++)	// ROWS
			{
				// Multiplication factor
				T temp_mult = A[r*n + row_index[i]] / A[r*n + row_index[r]];
				// First columin set 0
				A[r*n + row_index[i]] = (T)0;
				// Addl columins
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

template <typename T>
T CalcTriangleArea(Vec3<T> a, Vec3<T> b, Vec3<T> c)
{
	return ((b - a) % (c - a)).len() / 2;
}

template <typename T>
Vec3<T> CalcTriangleCentroid(Vec3<T> a, Vec3<T> b, Vec3<T> c)
{
	return (a + b + c) / 3;
}

template <typename T>
T CalcTriangleBoundRadius(Vec3<T> centroid, Vec3<T> a, Vec3<T> b, Vec3<T> c)
{
	T ra = (a - centroid).lenSq();
	T rb = (b - centroid).lenSq();
	T rc = (c - centroid).lenSq();
	return sqrt(std::max(std::max(ra, rb), rc));
}

template <typename T>
Vec3<T> ProjectVertexPlane(Vec3<T> v, Vec3<T> a, Vec3<T> b, Vec3<T> c)
{
	// Calculate normal / plane equation coefficients
	Vec4<T> plane = Vec4<T>(Vec3<T>((b - a) % (c - b)).norm(), 0);	// Calculate a, b, c in plane equation
	plane.w = plane.x * a.x + plane.y * a.y + plane.z * a.z;		// Calculate "d" in plane equation

	T translation = plane.w - (plane.xyz() * v);

	return plane.xyz() * translation;
}

template <typename T>
Vec3<T> CalcBarycentricCoefficients(Vec3<T> p, Vec3<T> a, Vec3<T> b, Vec3<T> c)
{
	// Barycentric coordinates can be negative, indicating p is outside the triangle
	// Convention for signed area is counter-clockwise = positive, clockwise = negative
	
	T area_abc = CalcTriangleArea(a, b, c);

	// We assume our abc triangle vertices are listed counter-clockwise
	// TODO: should the normal be passed in?
	// TODO: any other way to optimize the "signed area" checks?
	Vec3<T> n = (b - a) % (c - b);
	
	Vec3<T> temp;

	temp.x = CalcTriangleArea(b, c, p) / area_abc;
	if (((c - b) % (p - c)) * n < 0)
		temp.x = -temp.x;

	temp.y = CalcTriangleArea(c, a, p) / area_abc;
	if (((a - c) % (p - a)) * n < 0)
		temp.y = -temp.y;

	temp.z = CalcTriangleArea(a, b, p) / area_abc;
	if (((b - a) % (p - b)) * n < 0)
		temp.z = -temp.z;

	return temp;
}

// Pack four normalized floats (0-1) in a 32-bit integer
// Typically used for graphics applications - i.e 8-bit color vales
// Input:	v = vector of four floats (0-1)
// Return:  uint32_t - packed 32-bit integer
uint32_t PackFloatInt4x8(Vec4f v)
{
	uint8_t temp = 0;
	uint32_t result = 0;

	for (unsigned int i = 0; i < 4; i++)
	{
		temp = (uint8_t)round(clamp(v[i], 0.f, +1.f) * 255.0);

		result |= static_cast<uint32_t>(temp) << (i * 8);
	}

	return result;
}

// Pack four normalized floats (0-1) in a 32-bit integer
// Typically used for graphics applications - i.e 8-bit color vales
// Input:	v = vector of four floats (0-1)
// Return:  uint32_t - packed 32-bit integer
uint32_t PackFloatInt4x8(Vec4d v)
{
	return PackFloatInt4x8(Vec4f(float(v.x), float(v.y), float(v.z), float(v.w)));
}


// Generate random unique integers from a given range
// Inputs:	count = number of random integers to generate
//          start = start of range (inclusive)
//          end = end of range (inclusive)
//			exclude_count = size of exclude list
// 			exclude_list[] = array of integers to exclude from results
// Return:  results[] = pointer to array of size 'count' to return results
//			Return actual number of unique integers found
unsigned int RandIntsFromRange(unsigned int count, int start, int end, unsigned int exclude_count, int* exclude_list, int* results)
{
	// Validate inputs
	if (results == nullptr || start > end) {
		return 0;
	}

	// Build set of excluded values (only those within valid range)
	std::unordered_set<int> excluded_values;
	if (exclude_list != nullptr && exclude_count > 0) {
		for (unsigned int i = 0; i < exclude_count; ++i) {
			int value = exclude_list[i];
			// Only add if within valid range
			if (value >= start && value <= end) {
				excluded_values.insert(value);
			}
		}
	}

	// Calculate available range size (excluding forbidden values)
	long long range_size = static_cast<long long>(end) - start + 1;
	long long available_size = range_size - excluded_values.size();

	// Can't generate more unique numbers than are available
	if (available_size <= 0) {
		return 0;
	}

	unsigned int max_possible = static_cast<unsigned int>(std::min(static_cast<long long>(count), available_size));

	// Random number generator
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(start, end);

	// Use set to track uniqueness
	std::unordered_set<int> unique_numbers;

	// If we need more than half the available range, use different strategy
	if (max_possible > available_size / 2) {
		// Fill with all valid numbers (excluding forbidden ones), then shuffle
		std::vector<int> all_numbers;
		all_numbers.reserve(available_size);
		for (int i = start; i <= end; ++i) {
			if (excluded_values.find(i) == excluded_values.end()) {
				all_numbers.push_back(i);
			}
		}
		std::shuffle(all_numbers.begin(), all_numbers.end(), gen);

		for (unsigned int i = 0; i < max_possible; ++i) {
			results[i] = all_numbers[i];
		}
	}
	else {
		// Generate random numbers until we have enough unique ones
		while (unique_numbers.size() < max_possible) {
			int num = dist(gen);
			// Check if number is not excluded and not already selected
			if (excluded_values.find(num) == excluded_values.end()) {
				unique_numbers.insert(num);
			}
		}

		// Copy to results array
		unsigned int idx = 0;
		for (int num : unique_numbers) {
			results[idx++] = num;
		}
	}

	return max_possible;
}

#endif
