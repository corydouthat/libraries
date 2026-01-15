// *****************************************************************************************************************************
// Vec.h
// Linear Math Libraries
// 2D And 3D vector classes
// Template definition for data type
// Author(s): Cory Douthat
// Copyright (c) 2026 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef VEC_H_
#define VEC_H_

#define _USE_MATH_DEFINES	// For PI definition
#include <cmath>
#include <cstring>

template <typename T>
class Vec2;
template <typename T>
class Vec3;
template <typename T>
class Vec4;

template <typename T = float>
class Vec2
{
public:
	// DATA
	union
	{
		struct { T x; T y; };
		T v[2];
	};

	// FUNCTIONS
	// Constructors
	Vec2() : x(0), y(0) {}													// Constructor: initialize all to zeros
	template <typename sT>
	Vec2(const sT& bx, const sT& by) : x((T)bx), y((T)by) {}					// Constructor: x, y
	Vec2(const Vec2<T>& b) { memcpy(this, &b, sizeof(*this)); }				// Constructor: copy
	Vec2(const Vec3<T>& b) { memcpy(this, &b, sizeof(*this)); }				// Constructor: copy (downsize)
	template <typename T2>
	Vec2(const Vec2<T2>& b) { *this = b; }									// Constructor: copy (conversion)
	// Basic Operators
	const Vec2<T>& operator =(const Vec2<T>& b);							// Operator = (assignment)
	template <typename T2>
	const Vec2<T>& operator =(const Vec2<T2>& b);							// Operator = (assignment/conversion)
	T& operator [](unsigned int i) { return v[i]; }							// Operator []
	const T& operator [](unsigned int i)const { return v[i]; }				// Operator [] const
	bool operator ==(const Vec2<T>& b)const;								// Operator ==
	bool operator !=(const Vec2<T>& b)const { return !(*this == b); }		// Operator !=
	// Vector Operators
	Vec2<T> operator +(const Vec2<T>& b)const;								// Operator +
	Vec2<T> operator -(const Vec2<T>& b)const;								// Operator -
	Vec2<T> operator -()const;												// Operator - (inverse)
	T operator %(const Vec2<T>& b)const;								// Cross Product Operator %
	T operator *(const Vec2<T>& b)const;									// Dot Product Operator *
	const Vec2<T>& operator +=(const Vec2<T>& b) { return *this = *this + b; }	// Operator +=
	const Vec2<T>& operator -=(const Vec2<T>& b) { return *this = *this - b; }	// Operator -=
	// Scalar Operators
	Vec2<T> operator *(T s)const;											// Operator * (scalar)
	Vec2<T> operator /(T s)const;											// Operator / (scalar)
	const Vec2<T>& operator *=(T s) { return *this = *this * s; }			// Operator *= (scalar)
	const Vec2<T>& operator /=(T s) { return *this = *this / s; }			// Operator /= (scalar)
	// Other Member Functions
	T len()const { return sqrt(x * x + y * y); }								// Length
	T lenSq()const { return x * x + y * y; }									// Length squared
	Vec2<T> norm()const { return lenSq() != 0 ? *this / len() : *this; }	// Get normalized vector
	const Vec2<T>& normalize() { return lenSq() != 0 ? *this /= len() : *this; }	// Normalize vector
	T angle()const { return atan2(y, x); }									// Angle of vector (from zero)
	T angle(const Vec2<T>& b)const;											// Angle between vectors
	bool similar(const Vec2<T>& b, T margin = 0.01)const;					// Check if vectors are similar
	const T* getData()const { return (const T*)(&v); }						// Get pointer to raw data
	// Static Member Functions
	static T angle(const Vec2<T>& a, const Vec2<T>& b) { return a.angle(b); }	// Angle between vectors (static)
	static T dot(const Vec2<T>& a, const Vec2<T>& b) { return a * b; }			// Dot Product (static)
};

template <typename T = float>
class Vec3
{
public:
	// DATA
	union
	{
		struct { T x; T y; T z; };
		T v[3];
	};

	// FUNCTIONS
	// Constructors
	Vec3() : x(0), y(0), z(0) {}											// Constructor: initialize all to zeros
	template <typename sT>
	Vec3(const sT& bx, const sT& by, const sT& bz) : x((T)bx), y((T)by), z((T)bz) {}	// Constructor: x, y, z
	Vec3(const Vec3<T>& b) { memcpy(this, &b, sizeof(*this)); }				// Constructor: copy
	template <typename T2>
	Vec3(const Vec3<T2>& b) { *this = b; }									// Constructor: copy (conversion)
	Vec3(const Vec2<T>& b, const T& s = 0);									// Constructor: 2D vector
	// Basic Operators
	const Vec3<T>& operator =(const Vec3<T>& b);							// Operator =
	template <typename T2>
	const Vec3<T>& operator =(const Vec3<T2>& b);							// Operator = (assignment/conversion)
	T& operator [](unsigned int i) { return v[i]; }							// Operator []
	const T& operator [](unsigned int i)const { return v[i]; }				// Operator [] const
	bool operator ==(const Vec3<T>& b)const;								// Operator ==
	bool operator !=(const Vec3<T>& b)const { return !(*this == b); }		// Operator !=
	// Vector Operators
	Vec3<T> operator +(const Vec3<T>& b)const;								// Operator +
	Vec3<T> operator -(const Vec3<T>& b)const;								// Operator -
	Vec3<T> operator -()const;												// Operator - (inverse)
	Vec3<T> operator %(const Vec3<T>& b)const;								// Cross Product Operator %
	T operator *(const Vec3<T>& b)const;									// Dot Product Operator *
	const Vec3<T>& operator +=(const Vec3<T>& b) { return *this = *this + b; }	// Operator +=
	const Vec3<T>& operator -=(const Vec3<T>& b) { return *this = *this - b; }	// Operator -=
	const Vec3<T>& operator %=(const Vec3<T>& b) { return *this = *this % b; }	// Cross Product Operator %=
	// Scalar Operators
	Vec3<T> operator *(T s)const;											// Operator * (scalar)
	Vec3<T> operator /(T s)const;											// Operator / (scalar)
	const Vec3<T>& operator *=(T s) { return *this = *this * s; }			// Operator *= (scalar)
	const Vec3<T>& operator /=(T s) { return *this = *this / s; }			// Operator /= (scalar)
	// Other Member Functions
	T len()const { return sqrt(x * x + y * y + z * z); }							// Length
	T lenSq()const { return x * x + y * y + z * z; }								// Length squared
	Vec3<T> norm()const { return lenSq() != 0 ? *this / len() : *this; }	// Get normalized vector
	const Vec3<T>& normalize() { return lenSq() != 0 ? *this /= len() : *this; }// Normalize vector
	T angle(const Vec3<T>& b)const;											// Angle between vectors
	bool similar(const Vec3<T>& b, T margin = 0.01)const;					// Check if vectors are similar
	Vec2<T> xy() { return Vec2<T>(x, y); }									// Get XY components
	const T* getData()const { return (const T*)(&v); }						// Get pointer to raw data
	// Static Member Functions
	static T angle(const Vec3<T>& a, const Vec3<T>& b) { return a.angle(b); }	// Angle between vectors (static)
	static T cross(const Vec3<T>& a, const Vec3<T>& b) { return a % b; }		// Cross Product (static)
	static T dot(const Vec3<T>& a, const Vec3<T>& b) { return a * b; }			// Dot Product (static)
	static Vec3<T> perpendicular(const Vec3<T>& a, const Vec3<T>& b);			// Find perpendicular unit vector
};

template <typename T = float>
class Vec4
{
public:
	// DATA
	union
	{
		struct { T x; T y; T z; T w; };
		T v[4];
	};

	// FUNCTIONS
	// Constructors
	Vec4() : x(0), y(0), z(0), w(0) {}											// Constructor: initialize all to zeros
	template <typename sT>
	Vec4(const sT& bx, const sT& by, const sT& bz, const sT& bw) : x((T)bx), y((T)by), z((T)bz), w((T)bw) {}	    // Constructor: x, y, z, w
	Vec4(const Vec4<T>& b) { memcpy(this, &b, sizeof(*this)); }				// Constructor: copy
	template <typename T2>
	Vec4(const Vec4<T2>& b) { *this = b; }									// Constructor: copy (conversion)
	Vec4(const Vec3<T>& b, const T& s = 0);									// Constructor: 3D vector
	// Basic Operators
	const Vec4<T>& operator =(const Vec4<T>& b);							// Operator =
	template <typename T2>
	const Vec4<T>& operator =(const Vec4<T2>& b);							// Operator = (assignment/conversion)
	T& operator [](unsigned int i) { return v[i]; }							// Operator []
	const T& operator [](unsigned int i)const { return v[i]; }				// Operator [] const
	bool operator ==(const Vec4<T>& b)const;								// Operator ==
	bool operator !=(const Vec4<T>& b)const { return !(*this == b); }		// Operator !=
	// Vector Operators
	Vec4<T> operator +(const Vec4<T>& b)const;								// Operator +
	Vec4<T> operator -(const Vec4<T>& b)const;								// Operator -
	Vec4<T> operator -()const;												// Operator - (inverse)
	T operator *(const Vec4<T>& b)const;									// Dot Product Operator *
	const Vec4<T>& operator +=(const Vec4<T>& b) { return *this = *this + b; }	// Operator +=
	const Vec4<T>& operator -=(const Vec4<T>& b) { return *this = *this - b; }	// Operator -=
	// Scalar Operators
	Vec4<T> operator *(T s)const;											// Operator * (scalar)
	Vec4<T> operator /(T s)const;											// Operator / (scalar)
	const Vec4<T>& operator *=(T s) { return *this = *this * s; }			// Operator *= (scalar)
	const Vec4<T>& operator /=(T s) { return *this = *this / s; }			// Operator /= (scalar)
	// Other Member Functions
	T len()const { return sqrt(x * x + y * y + z * z + w * w); }					// Length
	T lenSq()const { return x * x + y * y + z * z + w * w; }						// Length squared
	Vec4<T> norm()const { return lenSq() != 0 ? *this / len() : *this; }	// Get normalized vector
	const Vec4<T>& normalize() { return lenSq() != 0 ? *this /= len() : *this; }	// Normalize vector
	bool similar(const Vec4<T>& b, T margin = 0.01)const;					// Check if vectors are similar
	Vec3<T> xyz() { return Vec3<T>(x, y, z); }									// Get XYZ components
	const T* getData()const { return (const T*)(&v); }						// Get pointer to raw data
	// Static Member Functions
	static T dot(const Vec4<T>& a, const Vec4<T>& b) { return a * b; }		// Dot Product (static)
};

// Shortcut Types
typedef Vec2<float> Vec2f;
typedef Vec3<float> Vec3f;
typedef Vec4<float> Vec4f;
typedef Vec2<double> Vec2d;
typedef Vec3<double> Vec3d;
typedef Vec4<double> Vec4d;
typedef Vec2<int> Vec2i;
typedef Vec3<int> Vec3i;
typedef Vec4<int> Vec4i;
typedef Vec2<short> Vec2s;
typedef Vec3<short> Vec3s;
typedef Vec4<short> Vec4s;


// Non-member functions
template <typename T, typename sT>
Vec2<T> operator *(sT s, const Vec2<T>& b);		// Operator scalar * Vec2<T>
template <typename T, typename sT>
Vec3<T> operator *(sT s, const Vec3<T>& b);		// Operator scalar * Vec3<T>
template <typename T, typename sT>
Vec4<T> operator *(sT s, const Vec4<T>& b);		// Operator scalar * Vec4<T>

#include "vec.cpp"

#endif