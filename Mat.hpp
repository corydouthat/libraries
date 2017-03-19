// *****************************************************************************************************************************
// Mat.hpp
// Linear Math Libraries
// Matrix classes 2x2, 3x3, //4x4
// The matrix data is stored in column-major format. [j] returns column j,
// and [j][i] returns the item at column j, row i.
// Author(s): Cory Douthat
// Contributor(s): Connor Douthat
// Copyright (c) 2015 Cory Douthat and Connor Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef MAT_HPP_
#define MAT_HPP_

#define _USE_MATH_DEFINES	//For PI definition
#include <cmath>
#include <cstring>

template <typename T = float>
class Mat2
{
	// DATA
	T v[4];
public:
	// FUNCTIONS
	// Constructors
	Mat2();												// Constructor: identity matrix
	Mat2(const Mat2<T> &b);								// Constructor: copy
	explicit Mat2(T s);									// Constructor: all elements to s
	Mat2(const T *data_in);								// Constructor: data array
	// Basic Operators
	const Mat2<T>& operator =(const Mat2<T> &b);		// Operator = (assignment)
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
	T get(unsigned int col,unsigned int row)const;		// Get element
	void set(unsigned int col,unsigned int row,T s);	// Set element
	const T* getData()const;							// Get pointer to raw data
	void load(const T* data_in);						// Load in data
	// Static Helpers
	static Mat2<T> zero();								// Generate zero matrix
	static Mat2<T> ident();								// Generate identity matrix
	static Mat2<T> rot(T theta);						// Generate 2D rotation matrix from angle in radians

	// Functions relating to vectors (Vec2)
#ifdef VEC_HPP_
	// Constructors
	Mat2(const Vec2<T> &v0, const Vec2<T> &v1);			// Constructor: two 2D vectors
	// Operators
	Vec2<T>& operator [](unsigned int col);				// Access column as vector
	const Vec2<T>& operator [](unsigned int col)const;	// Access column as vector
	Vec2<T> operator *(const Vec2<T> &b)const;			// Transform vector
	// 2D Transformation Matrix Generators
	static Mat2<T> scale(const Vec2<T> &t);				// 2D scale matrix from 2D vector
#endif

	// Friend declarations
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
	Mat3(const Mat3<T> &b);								// Constructor: copy
	explicit Mat3(T s);									// Constructor: all elements to s
	Mat3(const Mat2<T> &b, T last_element = 1);			// Constructor: from 2x2 matrix, last element usually 1
	Mat3(const T *data_in);								// Constructor: data array
	// Basic Operators
	const Mat3<T>& operator =(const Mat3<T> &b);		// Operator = (assignment)
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
	T get(unsigned int col,unsigned int row)const;		// Get element
	void set(unsigned int col,unsigned int row,T s);	// Set element
	const T* getData()const;							// Get pointer to raw data
	void load(const T* data_in);						// Load in data
	// Static Helpers
	static Mat3<T> zero();								// Generate zero matrix
	static Mat3<T> ident();								// Generate identity matrix
	static Mat3<T> transf(T theta, T x, T y);			// Generate 3D transformation matrix with rotation and translation

	// Functions relating to vectors (Vec3)
#ifdef VEC_HPP_
	// Constructors
	Mat3(const Vec3<T> &v0,const Vec3<T> &v1,const Vec3<T> &v2);	// Constructor: three 3D vectors
	// Operators
	Vec3<T>& operator [](unsigned int col);							// Access column as vector
	const Vec3<T>& operator [](unsigned int col)const;				// Access column as vector
	Vec3<T> operator *(const Vec3<T> &b)const;						// Transform 3D vector
	Vec2<T> operator *(const Vec2<T> &b)const;						// Transform 2D smaller vector
	// 2D Transformation Matrix Generators
	static Mat3<T> transl(const Vec3<T> &t);						// 2D translation matrix from 3D vector
	static Mat3<T> transl(const Vec2<T> &t);						// 2D translation matrix from 2D vector (3rd element 1)
	// 3D Transformation Matrix Generators
	static Mat3<T> rot(const Vec3<T> &axis);						// 3D rot matrix around axis (magnitude as radians)
	static Mat3<T> rot(const T &theta,const Vec3<T> &axis);			// 3D rot matrix theta radians around normalized axis
	static Mat3<T> scale(const Vec3<T> &t);							// 3D scale matrix from 3D vector
	static Mat3<T> transf(const Mat2<T> &r,const Vec3<T> &t);		// 3D transform matrix from 2D mat + 3D vector
	static Mat3<T> transf(const Mat2<T> &r,const Vec2<T> &t);		// 3D transform matrix from 2D mat + 2D vector
#endif

	// Friend declarations
	template <class T2>
	friend class Mat2;
};

typedef Mat2<float> Mat2f;
typedef Mat3<float> Mat3f;
typedef Mat2<double> Mat2d;
typedef Mat3<double> Mat3d;
typedef Mat2<int> Mat2i;
typedef Mat3<int> Mat3i;
typedef Mat2<short> Mat2s;
typedef Mat3<short> Mat3s;

// CONSTRUCTORS
// **************************************************************************

// Mat2
// Constructor: identity matrix
template <typename T>
Mat2<T>::Mat2()
{
	v[0] = 1;	v[2] = 0; 
	v[1] = 0;	v[3] = 1;
}

// Mat3
// Constructor: identity matrix
template <typename T>
Mat3<T>::Mat3()
{
	v[0] = 1;	v[3] = 0;	v[6] = 0;
	v[1] = 0;	v[4] = 1;	v[7] = 0;
	v[2] = 0;	v[5] = 0;	v[8] = 1;
}

// Mat2
// Constructor: copy
template <typename T>
Mat2<T>::Mat2(const Mat2<T> &b)
{
	memcpy(this,&b,sizeof(*this));
}

// Mat3
// Constructor: copy
template <typename T>
Mat3<T>::Mat3(const Mat3<T> &b)
{
	memcpy(this,&b,sizeof(*this));
}

// Mat2
// Constructor: all elements to s
template <typename T>
Mat2<T>::Mat2(T s)
{
	memset(&v,s,sizeof(v));
}

// Mat3
// Constructor: all elements to s
template <typename T>
Mat3<T>::Mat3(T s)
{
	memset(&v,s,sizeof(v));
}

// Mat2
// Constructor: data array
template <typename T>
Mat2<T>::Mat2(const T *data_in)
{
	memcpy(&v,data_in,sizeof(v));
}

// Mat2
// Constructor: data array
template <typename T>
Mat3<T>::Mat3(const T *data_in)
{
	memcpy(&v,data_in,sizeof(v));
}

// Mat3
// Constructor: from 2x2 matrix, last element usually 1
template <typename T>
Mat3<T>::Mat3(const Mat2<T> &b,T last_element)
{
	v[0] = b.v[0];	v[3] = b.v[2];	v[6] = 0;
	v[1] = b.v[1];	v[4] = b.v[3];	v[7] = 0;
	v[2] = 0;		v[5] = 0;		v[8] = last_element;
}


// BASIC OPERATORS
// **************************************************************************

// Mat2
// Operator = (assignment)
template <typename T>
const Mat2<T>& Mat2<T>::operator =(const Mat2<T> &b)
{
	memcpy(this,&b,sizeof(*this));
	return *this;
}


// Mat3
// Operator = (assignment)
template <typename T>
const Mat3<T>& Mat3<T>::operator =(const Mat3<T> &b)
{
	memcpy(this,&b,sizeof(*this));
	return *this;
}

// Mat2
// Operator ==
template <typename T>
bool Mat2<T>::operator ==(const Mat2<T>& b)const
{
	if (memcmpy(this,&b,sizeof(*this)) == 0)
		return true;
	else
		return false;
}

// Mat3
// Operator ==
template <typename T>
bool Mat3<T>::operator ==(const Mat3<T>& b)const
{
	if (memcmpy(this,&b,sizeof(*this)) == 0)
		return true;
	else
		return false;
}

// Mat2
// Operator !=
template <typename T>
bool Mat2<T>::operator !=(const Mat2<T>& b)const
{
	return !(*this == b);
}

// Mat3
// Operator !=
template <typename T>
bool Mat3<T>::operator !=(const Mat3<T>& b)const
{
	return !(*this == b);
}


// MATRIX OPERATORS
// **************************************************************************

// Mat2
// Operator +
template <typename T>
Mat2<T> Mat2<T>::operator +(const Mat2<T>& b)const
{
	Mat2<T> temp(*this);
	for (int i = 0; i < 4; i++)
		temp.v[i] += b.v[i];
	return temp;
}

// Mat3
// Operator +
template <typename T>
Mat3<T> Mat3<T>::operator +(const Mat3<T>& b)const
{
	Mat3<T> temp(*this);
	for (int i = 0; i < 9; i++)
		temp.v[i] += b.v[i];
	return temp;
}

// Mat2
// Operator -
template <typename T>
Mat2<T> Mat2<T>::operator -(const Mat2<T>& b)const
{
	Mat2<T> temp(*this);
	for (int i = 0; i < 4; i++)
		temp.v[i] -= b.v[i];
	return temp;
}

// Mat3
// Operator -
template <typename T>
Mat3<T> Mat3<T>::operator -(const Mat3<T>& b)const
{
	Mat3<T> temp(*this);
	for (int i = 0; i < 9; i++)
		temp.v[i] -= b.v[i];
	return temp;
}

// Mat2
// Operator - (inverse)
template <typename T>
Mat2<T> Mat2<T>::operator -()const
{
	Mat2<T> temp(*this);
	for (int i = 0; i < 4; i++)
		temp.v[i] = -temp.v[i];
	return temp;
}

// Mat3
// Operator - (inverse)
template <typename T>
Mat3<T> Mat3<T>::operator -()const
{
	Mat3<T> temp(*this);
	for (int i = 0; i < 9; i++)
		temp.v[i] = -temp.v[i];
	return temp;
}

// Mat2
// Operator *
template <typename T>
Mat2<T> Mat2<T>::operator *(const Mat2<T>& b)const
{
	Mat2<T> temp((T)0);

	for (int i = 0; i < 2; i++)				//i = rows of a
	{
		for (int j = 0; j < 2; j++)			//j = columns of b
		{
			for (int k = 0; k < 2; k++)		//k = columns of a, rows of b
			{
				temp.v[2 * j + i] += (*this).v[2 * k + i] * b.v[2 * j + k];
			}
		}
	}

	return temp;
}

// Mat3
// Operator *
template <typename T>
Mat3<T> Mat3<T>::operator *(const Mat3<T>& b)const
{
	Mat3<T> temp((T)0);

	for (int i = 0; i < 3; i++)				//i = rows of a
	{
		for (int j = 0; j < 3; j++)			//j = columns of b
		{
			for (int k = 0; k < 3; k++)		//k = columns of a, rows of b
			{
				temp.v[3 * j + i] += (*this).v[3 * k + i] * b.v[3 * j + k];
			}
		}
	}

	return temp;
}

// Mat2
// Operator +=
template <typename T>
const Mat2<T>& Mat2<T>::operator +=(const Mat2<T>& b)
{
	for (int i = 0; i < 4; i++)
		(*this).v[i] += b.v[i];
	return *this;
}

// Mat3
// Operator +=
template <typename T>
const Mat3<T>& Mat3<T>::operator +=(const Mat3<T>& b)
{
	for (int i = 0; i < 9; i++)
		(*this).v[i] += b.v[i];
	return *this;
}

// Mat2
// Operator -=
template <typename T>
const Mat2<T>& Mat2<T>::operator -=(const Mat2<T>& b)
{
	for (int i = 0; i < 4; i++)
		(*this).v[i] -= b.v[i];
	return *this;
}

// Mat3
// Operator -=
template <typename T>
const Mat3<T>& Mat3<T>::operator -=(const Mat3<T>& b)
{
	for (int i = 0; i < 9; i++)
		(*this).v[i] -= b.v[i];
	return *this;
}

// Mat2
// Operator *=
template <typename T>
const Mat2<T>& Mat2<T>::operator *=(const Mat2<T>& b)
{
	return *this = *this * b;
}

// Mat3
// Operator *=
template <typename T>
const Mat3<T>& Mat3<T>::operator *=(const Mat3<T>& b)
{
	return *this = *this * b;
}


// SCALAR OPERATORS
// **************************************************************************

// Mat2
// Operator * (scalar)
template <typename T>
Mat2<T> Mat2<T>::operator *(T s)const
{
	Mat2<T> temp(*this);
	for (int i = 0; i < 4; i++)
		temp.v[i] *= s;
	return temp;
}

// Mat3
// Operator * (scalar)
template <typename T>
Mat3<T> Mat3<T>::operator *(T s)const
{
	Mat2<T> temp(*this);
	for (int i = 0; i < 9; i++)
		temp.v[i] *= s;
	return temp;
}

// Mat2
// Operator *= (scalar)
template <typename T>
const Mat2<T>& Mat2<T>::operator *=(T s)
{
	for (int i = 0; i < 4; i++)
		(*this).v[i] *= s;
	return *this;
}

// Mat3
// Operator *= (scalar)
template <typename T>
const Mat3<T>& Mat3<T>::operator *=(T s)
{
	for (int i = 0; i < 9; i++)
		(*this).v[i] *= s;
	return *this;
}


// OTHER MEMBER FUNCTIONS
// **************************************************************************

// Mat2
// Determinant
template <typename T>
T Mat2<T>::det()const
{
	return v[0] * v[3] - v[2] * v[1];
}

// Mat3
// Determinant
template <typename T>
T Mat3<T>::det()const
{
	return	v[0] * v[4] * v[8] -
			v[0] * v[7] * v[5] -
			v[3] * v[1] * v[8] +
			v[3] * v[7] * v[2] +
			v[6] * v[1] * v[5] -
			v[6] * v[4] * v[2];
}

// Mat2
// Transpose
template <typename T>
Mat2<T> Mat2<T>::transp()const
{
	Mat2<T> temp(*this);
	temp.v[1] = v[2];
	temp.v[2] = v[1];
	return temp;
}

// Mat3
// Transpose
template <typename T>
Mat3<T> Mat3<T>::transp()const
{
	Mat3<T> temp(*this);
	temp.v[1] = v[3];	temp.v[2] = v[6];
	temp.v[3] = v[1];	temp.v[5] = v[7];
	temp.v[6] = v[2];	temp.v[7] = v[5];
	return temp;
}

// Mat2
// Inverse
template <typename T>
Mat2<T> Mat2<T>::inv()const
{
	Mat2<T> temp((T)0);
	// Check if matrix is inversible by checking that determinant != 0
	if (T d = det())
	{
		// Matrix is inversible
		temp.v[0] = (*this).v[3] / d;
		temp.v[1] = -(*this).v[1] / d;
		temp.v[2] = -(*this).v[2] / d;
		temp.v[3] = (*this).v[0] / d;
	}

	// If matrix is not inversible - return zero matrix to make detection of this issue more easily
	return temp;
}

// Mat3
// Inverse
template <typename T>
Mat3<T> Mat3<T>::inv()const
{
	Mat3<T> temp((T)0);
	// Check if matrix is inversible by checking that determinant != 0
	if (T d = det())
	{
		// Matrix is inversible
		temp.v[0] = ((*this).v[4] * (*this).v[8] - (*this).v[7] * (*this).v[5]) / d;
		temp.v[1] = -((*this).v[1] * (*this).v[8] - (*this).v[7] * (*this).v[2]) / d;
		temp.v[2] = ((*this).v[1] * (*this).v[5] - (*this).v[4] * (*this).v[2]) / d;
		temp.v[3] = -((*this).v[3] * (*this).v[8] - (*this).v[6] * (*this).v[5]) / d;
		temp.v[4] = ((*this).v[0] * (*this).v[8] - (*this).v[6] * (*this).v[2]) / d;
		temp.v[5] = -((*this).v[0] * (*this).v[5] - (*this).v[3] * (*this).v[2]) / d;
		temp.v[6] = ((*this).v[3] * (*this).v[7] - (*this).v[6] * (*this).v[4]) / d;
		temp.v[7] = -((*this).v[0] * (*this).v[7] - (*this).v[6] * (*this).v[1]) / d;
		temp.v[8] = ((*this).v[0] * (*this).v[4] - (*this).v[3] * (*this).v[1]) / d;
	}

	// If matrix is not inversible - return zero matrix to make detection of this issue more easily
	return temp;
}

// Mat2
// Get Element
template <typename T>
T Mat2<T>::get(unsigned int col,unsigned int row)const
{
	if (col > 1 || row > 1) return 0;
	return v[col * 2 + row];
}

// Mat3
// Get Element
template <typename T>
T Mat3<T>::get(unsigned int col,unsigned int row)const
{
	if (col > 2 || row > 2) return 0;
	return v[col * 3 + row];
}

// Mat2
// Set Element
template <typename T>
void Mat2<T>::set(unsigned int col,unsigned int row,T s)
{
	if (col > 1 || row > 1) return;
	v[col * 2 + row] = s;
}

// Mat3
// Set Element
template <typename T>
void Mat3<T>::set(unsigned int col,unsigned int row,T s)
{
	if (col > 2 || row > 2) return;
	v[col * 3 + row] = s;
}

// Mat 2
// Get Data Pointer
template <typename T>
const T* Mat2<T>::getData()const
{
	return (const T*)(&v);
}

// Mat 3
// Get Data Pointer
template <typename T>
const T* Mat3<T>::getData()const
{
	return (const T*)(&v);
}

// Mat 2
// Load Data from Array
template <typename T>
void Mat2<T>::load(const T* data_in)
{
	memcpy(&v,data_in,sizeof(v));
}

// Mat3
// Load Data from Array
template <typename T>
void Mat3<T>::load(const T* data_in)
{
	memcpy(&v,data_in,sizeof(v));
}


// STATIC HELPERS
// **************************************************************************

// Mat2
// Static Zero Matrix Generator
template <typename T>
Mat2<T> Mat2<T>::zero()
{
	return Mat2<T>(0);
}

// Mat3
// Static Zero Matrix Generator
template <typename T>
Mat3<T> Mat3<T>::zero()
{
	return Mat3<T>(0);
}

// Mat2
// Static Identity Matrix Generator
template <typename T>
Mat2<T> Mat2<T>::ident()
{
	return Mat2<T>();
}

// Mat3
// Static Identity Matrix Generator
template <typename T>
Mat3<T> Mat3<T>::ident()
{
	return Mat3<T>();
}

// Mat2
// Static 2D Rotation Matrix Generator
// Inputs:	theta = angle in radians
template <typename T>
Mat2<T> Mat2<T>::rot(T theta)
{
	Mat2<T> temp(0.0);
	temp.v[0] = cos(theta);	temp.v[2] = -sin(theta);
	temp.v[1] = sin(theta);	temp.v[3] = cos(theta);
	return temp;
}

// Mat3
// Static 3D Transformation Matrix Generator
// Inputs:	theta = angle in radians
//			x, y = translation along x and y axis
template <typename T>
Mat3<T> Mat3<T>::transf(T theta,T x,T y)
{
	Mat3<T> temp();	// Identity matrix
	temp.v[0] = cos(theta);	temp.v[3] = -sin(theta);	temp.v[6] = x;
	temp.v[1] = sin(theta);	temp.v[4] = cos(theta);		temp.v[7] = y;
	return temp;
}


// VECTOR RELATED FUNCTIONS
// **************************************************************************

#ifdef VEC_HPP_

// Mat2
// Constructor: two 2D vectors
template <typename T>
Mat2<T>::Mat2(const Vec2<T> &v0,const Vec2<T> &v1)
{
	v[0] = v0[0];	v[2] = v1[0]; 
	v[1] = v0[1];	v[3] = v1[1];
}

// Mat3
// Constructor: three 3D vectors
template <typename T>
Mat3<T>::Mat3(const Vec3<T> &v0,const Vec3<T> &v1,const Vec3<T> &v2)
{
	v[0] = v0[0];	v[3] = v1[0];	v[6] = v2[0];
	v[1] = v0[1];	v[4] = v1[1];	v[7] = v2[1];
	v[2] = v0[2];	v[5] = v1[2];	v[8] = v2[2];
}

// Mat2
// Operator []: Access column as vector (Vec2)
template <typename T>
Vec2<T>& Mat2<T>::operator [](unsigned int col)
{
	return *(reinterpret_cast<Vec2<T>*>((T*)v + col * 2));
}

// Mat3
// Operator []: Access column as vector (Vec3)
template <typename T>
Vec3<T>& Mat3<T>::operator [](unsigned int col)
{
	return *(reinterpret_cast<Vec3<T>*>((T*)v + col * 3));
}

// Mat2
// Operator []: Access column as vector (Vec2) (const)
template <typename T>
const Vec2<T>& Mat2<T>::operator [](unsigned int col)const
{
	return *(reinterpret_cast<const Vec2<T>*>((T*)v + col * 2));
}

// Mat3
// Operator []: Access column as vector (Vec3) (const)
template <typename T>
const Vec3<T>& Mat3<T>::operator [](unsigned int col)const
{
	return *(reinterpret_cast<const Vec3<T>*>((T*)v + col * 3));
}

// Mat2
// Transform a 2D Vector
template <typename T>
Vec2<T> Mat2<T>::operator *(const Vec2<T> &b)const
{
	return Vec2<T> (v[0] * b.x + v[2] * b.y,
					v[1] * b.x + v[3] * b.y);
}

// Mat3
// Transform a 3D Vector
template <typename T>
Vec3<T> Mat3<T>::operator *(const Vec3<T> &b)const
{
	return Vec3<T> (v[0] * b.x + v[3] * b.y + v[6] * b.z,
					v[1] * b.x + v[4] * b.y + v[7] * b.z,
					v[2] * b.x + v[5] * b.y + v[8] * b.z);
}

// Mat3
// Transform a 2D Vector (by 3x3 matrix)
template <typename T>
Vec2<T> Mat3<T>::operator *(const Vec2<T> &b)const
{
	return Vec2<T> (v[0] * b.x + v[3] * b.y + v[6] * 1,
					v[1] * b.x + v[4] * b.y + v[7] * 1);
}

// Mat2
// Static 2D Scale Matrix Generator
template <typename T>
Mat2<T> Mat2<T>::scale(const Vec2<T> &t)
{
	Mat2<T> temp;	// Identity matrix
	temp.v[0] = t.x;
	temp.v[3] = t.y;
	return temp;
}

// Mat3
// Static 3D Scale Matrix Generator
template <typename T>
Mat3<T> Mat3<T>::scale(const Vec3<T> &t)
{
	Mat3<T> temp;	// Identity matrix
	temp.v[0] = t.x;
	temp.v[4] = t.y;
	temp.v[8] = t.z;
	return temp;
}

// Mat3
// Static 2D Translation Matrix Generator: 3D Vector (3rd element usually 1)
template <typename T>
Mat3<T> Mat3<T>::transl(const Vec3<T> &t)
{
	Mat3<T> temp;	// Identity matrix
	temp.v[6] = t.x;
	temp.v[7] = t.y;
	temp.v[8] = t.z;
	return temp;
}

// Mat3
// Static 2D Translation Matrix Generator: 2D Vector (3rd element is 1)
template <typename T>
Mat3<T> Mat3<T>::transl(const Vec2<T> &t)
{
	Mat3<T> temp;	// Identity matrix
	temp.v[6] = t.x;
	temp.v[7] = t.y;
	return temp;
}

// Mat3
// Static 3D Rotation Matrix Generator: Scaled vector
// Axis vector magnitude as rotation angle in radians
template <typename T>
Mat3<T> Mat3<T>::rot(const Vec3<T> &axis)
{
	T length = axis.len();
	if (length == 0) 
		return Mat3<T>();	// Identity matrix
	else
		return Mat3<T>::rot(length, axis / length);
}

// Mat3
// Static 3D Rotation Matrix Generator: Angle, axis
// Angle in radians, axis must be normalized
template <typename T>
Mat3<T> Mat3<T>::rot(const T &theta,const Vec3<T> &axis)
{
	Mat3<T> temp;	// Identity matrix

	if (!theta || (!(axis.x) && !(axis.y) && !(axis.z))) 
		return temp;
	else
	{
		T x = axis.x;
		T y = axis.y;
		T z = axis.z;
		T c = cos(theta);
		T s = sin(theta);
		T t = 1 - cos(theta);

		temp[0] = t * x * x + c;
		temp[1] = t * x * y + s * z;
		temp[2] = t * x * z - s * y;
		temp[3] = t * x * y - s * z;
		temp[4] = t * y * y + c;
		temp[5] = t * y * z + s * x;
		temp[6] = t * x * z + s * y;
		temp[7] = t * y * z - s * x;
		temp[8] = t * z * z + c;

		return temp;
	}
}

// Mat3
// Static 3D Transformation Matrix Generator: 2x2 matrix + 3D vector (last element usually 1)
template <typename T>
Mat3<T> Mat3<T>::transf(const Mat2<T> &r,const Vec3<T> &t)
{
	Mat3<T> temp(r);
	temp.v[6] = t.x;
	temp.v[7] = t.y;
	temp.v[8] = t.z;
	return temp;
}

// Mat3
// Static 3D Transformation Matrix Generator: 2x2 matrix + 2D vector (last element 1)
template <typename T>
Mat3<T> Mat3<T>::transf(const Mat2<T> &r,const Vec2<T> &t)
{
	Mat3<T> temp(r);
	temp.v[6] = t.x;
	temp.v[7] = t.y;
	return temp;
}

#endif

#endif
