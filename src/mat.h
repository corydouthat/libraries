// *****************************************************************************************************************************
// Mat.h
// Linear Math Libraries
// Matrix classes 2x2, 3x3, 4x4
// The matrix data is stored in column-major format. [j] returns column j,
// and [j][i] returns the item at column j, row i.
// Author(s): Cory Douthat
// Copyright (c) 2026 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef MAT_HPP_
#define MAT_HPP_

#define _USE_MATH_DEFINES	//For M_PI definition
#include <cmath>
#include <numbers>
#include <cstring>

#include "vec.h"
#include "quat.h"

template <typename T>
class Quat;

template <typename T = float>
class Mat2
{
	// DATA
	T v[4];
public:
	// FUNCTIONS
	// Constructors
	Mat2();												// Constructor: identity matrix
	Mat2(const Mat2<T>& b);								// Constructor: copy
	explicit Mat2(T s);									// Constructor: all elements to s
	Mat2(const T* data_in);								// Constructor: data array
	// Basic Operators
	const Mat2<T>& operator =(const Mat2<T>& b);		// Operator = (assignment)
	template <typename T2>
	const Mat2<T>& operator =(const Mat2<T2>& b);		// Operator = (assignment/conversion)
	bool operator ==(const Mat2<T>& b)const;			// Operator ==
	bool operator !=(const Mat2<T>& b)const;			// Operator !=
	// Matrix Operators
	Mat2<T> operator +(const Mat2<T>& b)const;			// Operator +
	Mat2<T> operator -(const Mat2<T>& b)const;			// Operator -
	Mat2<T> operator -()const;							// Operator - (inverse)
	Mat2<T> operator *(const Mat2<T>& b)const;			// Operator *
	const Mat2<T>& operator +=(const Mat2<T>& b);		// Operator +=
	const Mat2<T>& operator -=(const Mat2<T>& b);		// Operator -=
	const Mat2<T>& operator *=(const Mat2<T>& b);		// Operator *=
	// Scalar Operators
	Mat2<T> operator *(T s)const;						// Operator * (scalar)
	const Mat2<T>& operator *=(T s);					// Operator *= (scalar)
	// Other Member Functions
	T det()const;										// Determinant
	Mat2<T> transp()const;								// Transpose
	Mat2<T> inv()const;									// Inverse
	T get(unsigned int col, unsigned int row)const;		// Get element
	void set(unsigned int col, unsigned int row, T s);	// Set element
	const T* getData()const;							// Get pointer to raw data
	void load(const T* data_in);						// Load in data
	// Static Helpers
	static Mat2<T> zero();								// Generate zero matrix
	static Mat2<T> ident();								// Generate identity matrix
	static Mat2<T> rot(T theta);						// Generate 2D rotation matrix from angle in radians

	// Functions relating to vectors (Vec2)
//#ifdef VEC_HPP_
	// Constructors
	Mat2(const Vec2<T>& v0, const Vec2<T>& v1);			// Constructor: two 2D vectors
	// Operators
	Vec2<T>& operator [](unsigned int col);				// Access column as vector
	const Vec2<T>& operator [](unsigned int col)const;	// Access column as vector
	Vec2<T> operator *(const Vec2<T>& b)const;			// Transform vector
	// 2D Transformation Matrix Generators
	static Mat2<T> scale(const Vec2<T>& t);				// 2D scale matrix from 2D vector
//#endif

	// Friend declarations
	template <class T2>
	friend class Mat2;
	template <class T2>
	friend class Mat3;
};

template <typename T = float>
class Mat3
{
	// DATA
	T v[9];
public:
	// FUNCTIONS
	// Constructors
	Mat3();												// Constructor: identity matrix
	Mat3(const Mat3<T>& b);								// Constructor: copy
	explicit Mat3(T s);									// Constructor: all elements to s
	Mat3(const Mat2<T>& b, T last_element = 1);			// Constructor: from 2x2 matrix, last element usually 1
	Mat3(const T* data_in);								// Constructor: data array
	// Basic Operators
	const Mat3<T>& operator =(const Mat3<T>& b);		// Operator = (assignment)
	template <typename T2>
	const Mat3<T>& operator =(const Mat3<T2>& b);		// Operator = (assignment/conversion)
	bool operator ==(const Mat3<T>& b)const;			// Operator ==
	bool operator !=(const Mat3<T>& b)const;			// Operator !=
	// Matrix Operators
	Mat3<T> operator +(const Mat3<T>& b)const;			// Operator +
	Mat3<T> operator -(const Mat3<T>& b)const;			// Operator -
	Mat3<T> operator -()const;							// Operator - (inverse)
	Mat3<T> operator *(const Mat3<T>& b)const;			// Operator *
	const Mat3<T>& operator +=(const Mat3<T>& b);		// Operator +=
	const Mat3<T>& operator -=(const Mat3<T>& b);		// Operator -=
	const Mat3<T>& operator *=(const Mat3<T>& b);		// Operator *=
	// Scalar Operators
	Mat3<T> operator *(T s)const;						// Operator * (scalar)
	const Mat3<T>& operator *=(T s);					// Operator *= (scalar)
	// Other Member Functions
	T det()const;										// Determinant
	Mat3<T> transp()const;								// Transpose
	Mat3<T> inv()const;									// Inverse
	void decomposeRotScale(Mat3<T>* rot, Vec3<T>* scale)const;	// Decompose
	void decomposeRotScale(Mat3<T>* rot, Mat3<T>* scale)const;	// Decompose
	T get(unsigned int col, unsigned int row)const;		// Get element
	Mat2<T> getSub()const;								// Get 2x2 sub-matrix (upper-left)
	Vec3<T> getTransl()const;							// Get translation vector (last column)
	void set(unsigned int col, unsigned int row, T s);	// Set element
	const T* getData()const;							// Get pointer to raw data
	void load(const T* data_in);						// Load in data
	// Static Helpers
	static Mat3<T> zero();								// Generate zero matrix
	static Mat3<T> ident();								// Generate identity matrix
	static Mat3<T> transf(T theta, T x, T y);			// Generate 2D transformation matrix with rotation and translation
	// Static Functions (other)
	static Mat3<T> chgBasis(Mat3<T> m, Mat3<T> b);		// Change basis of a matrix

	// Functions relating to vectors (Vec3)
//#ifdef VEC_HPP_
	// Constructors
	Mat3(const Vec3<T>& v0, const Vec3<T>& v1, const Vec3<T>& v2);	// Constructor: three 3D vectors
	// Operators
	Vec3<T>& operator [](unsigned int col);							// Access column as vector
	const Vec3<T>& operator [](unsigned int col)const;				// Access column as vector
	Vec3<T> operator *(const Vec3<T>& b)const;						// Transform 3D vector
	Vec2<T> operator *(const Vec2<T>& b)const;						// Transform 2D smaller vector (assume third component is 1)
	// 2D Transformation Matrix Generators
	static Mat3<T> transl(const Vec3<T>& t);						// 2D translation matrix from 3D vector
	static Mat3<T> transl(const Vec2<T>& t);						// 2D translation matrix from 2D vector (3rd element 1)
	// 3D Transformation Matrix Generators
	static Mat3<T> rot(const Vec3<T>& axis);						// 3D rot matrix around axis (magnitude as radians)
	static Mat3<T> rot(const T& theta, const Vec3<T>& axis);		// 3D rot matrix theta radians around normalized axis
	static Mat3<T> rot(const Vec3<T>& a, const Vec3<T>& b);			// 3D rot matrix from vector a to b
	static Mat3<T> scale(const Vec3<T>& t);							// 3D scale matrix from 3D vector
	static Mat3<T> transf(const Mat2<T>& r, const Vec3<T>& t);		// 2D transform matrix from 2D mat + 3D vector
	static Mat3<T> transf(const Mat2<T>& r, const Vec2<T>& t);		// 2D transform matrix from 2D mat + 2D vector
	static Mat3<T> transf(T theta, Vec2<T>& t);						// 2D transform matrix from rotation and translation
	// Other
	static Mat3<T> skewSymCross(Vec3<T> v);							// Generate skew-symmetric (cross product) matrix from vector
//#endif

//#ifdef QUAT_HPP_
	static Mat3<T> rot(const Quat<T> q);				// 3D rot matrix around axis (magnitude as radians)
//#endif

	// Friend declarations
	template <class T2>
	friend class Mat2;
	template <class T2>
	friend class Mat3;
	template <class T2>
	friend class Mat4;
};

template <typename T = float>
class Mat4
{
	// DATA
	T v[16];
public:
	// FUNCTIONS
	// Constructors
	Mat4();												// Constructor: identity matrix
	Mat4(const Mat4<T>& b);								// Constructor: copy
	explicit Mat4(T s);									// Constructor: all elements to s
	Mat4(const Mat3<T>& b, T last_element = 1);			// Constructor: from 3x3 matrix, last element usually 1
	Mat4(const T* data_in);								// Constructor: data array
	// Basic Operators
	const Mat4<T>& operator =(const Mat4<T>& b);		// Operator = (assignment)
	template <typename T2>
	const Mat4<T>& operator =(const Mat4<T2>& b);		// Operator = (assignment/conversion)
	bool operator ==(const Mat4<T>& b)const;			// Operator ==
	bool operator !=(const Mat4<T>& b)const;			// Operator !=
	// Matrix Operators
	Mat4<T> operator +(const Mat4<T>& b)const;			// Operator +
	Mat4<T> operator -(const Mat4<T>& b)const;			// Operator -
	Mat4<T> operator -()const;							// Operator - (inverse)
	Mat4<T> operator *(const Mat4<T>& b)const;			// Operator *
	const Mat4<T>& operator +=(const Mat4<T>& b);		// Operator +=
	const Mat4<T>& operator -=(const Mat4<T>& b);		// Operator -=
	const Mat4<T>& operator *=(const Mat4<T>& b);		// Operator *=
	// Scalar Operators
	Mat4<T> operator *(T s)const;						// Operator * (scalar)
	const Mat4<T>& operator *=(T s);					// Operator *= (scalar)
	// Other Member Functions
	T det()const;										// Determinant
	Mat4<T> transp()const;								// Transpose
	Mat4<T> inv()const;									// Inverse
	void decomposeRotScale(Mat3<T>* rot, Vec3<T>* scale)const;			// Decompose
	void decomposeTransfScale(Mat4<T>* transf, Vec3<T>* scale)const;	// Decompose
	void decomposeTransfScale(Mat4<T>* transf, Mat4<T>* scale)const;	// Decompose
	T get(unsigned int col, unsigned int row)const;		// Get element
	Mat3<T> getSub()const;								// Get 3x3 sub-matrix (upper-left)
	Vec3<T> getTransl()const;							// Get translation vector (last column)
	void set(unsigned int col, unsigned int row, T s);	// Set element
	const T* getData()const;							// Get pointer to raw data
	void load(const T* data_in);						// Load in data
	// Static Helpers
	static Mat4<T> zero();								// Generate zero matrix
	static Mat4<T> ident();								// Generate identity matrix
	// Static Projection Matrix Generators
	static Mat4<T> projOrtho(T left, T right, T bottom, T top, T d_near, T d_far);
	static Mat4<T> projPerspective(T fov_y, T aspect, T d_near, T d_far, bool clip_opengl = false);
	// Static Functions (other)
	static Mat4<T> chgBasis(Mat4<T> m, Mat4<T> b);		// Change basis of a matrix

	// Functions relating to vectors (Vec4/Vec3)
//#ifdef VEC_HPP_
	// Constructors
	Mat4(const Vec4<T>& v0, const Vec4<T>& v1, const Vec4<T>& v2, const Vec4<T>& v3);	// Constructor: Four 4D vectors
	// Operators
	Vec4<T>& operator [](unsigned int col);							// Access column as vector
	const Vec4<T>& operator [](unsigned int col)const;				// Access column as vector
	Vec4<T> operator *(const Vec4<T>& b)const;						// Transform 4D vector
	Vec3<T> operator *(const Vec3<T>& b)const;						// Transform 3D smaller vector (assume fourth component is 1)
	// 2D Transformation Matrix Generators
	static Mat4<T> transl(const Vec4<T>& t);						// 3D translation matrix from 4D vector
	static Mat4<T> transl(const Vec3<T>& t);						// 3D translation matrix from 3D vector (4th element 1)
	// 3D Transformation Matrix Generators
	static Mat4<T> transf(const Mat3<T>& r, const Vec4<T>& t);		// 3D transform matrix from 3D mat + 4D vector
	static Mat4<T> transf(const Mat3<T>& r, const Vec3<T>& t);		// 3D transform matrix from 3D mat + 3D vector
	static Mat4<T> transf(const Mat3<T>& trf);						// 3D transform matrix from 2D transform matrix (affine)
//#endif

//#ifdef QUAT_HPP_
	static Mat4<T> transf(const Quat<T>& q, const Vec4<T>& t);		// 3D transform matrix from 3D mat + 4D vector
	static Mat4<T> transf(const Quat<T>& q, const Vec3<T>& t);		// 3D transform matrix from 3D mat + 3D vector
//#endif

	// Friend declarations
	template <class T2>
	friend class Mat3;
	template <class T2>
	friend class Mat4;
};

typedef Mat2<float> Mat2f;
typedef Mat3<float> Mat3f;
typedef Mat4<float> Mat4f;
typedef Mat2<double> Mat2d;
typedef Mat3<double> Mat3d;
typedef Mat4<double> Mat4d;
typedef Mat2<int> Mat2i;
typedef Mat3<int> Mat3i;
typedef Mat4<int> Mat4i;
typedef Mat2<short> Mat2s;
typedef Mat3<short> Mat3s;
typedef Mat4<short> Mat4s;

#endif