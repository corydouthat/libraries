// *****************************************************************************************************************************
// Vec.hpp
// Linear Math Libraries
// 2D And 3D vector classes
// Template definition for data type
// Author(s): Cory Douthat
// Copyright (c) 2020 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef VEC_HPP_
#define VEC_HPP_

#define _USE_MATH_DEFINES	// For PI definition
#include <cmath>
#include <cstring>

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
    Vec2() : x(0),y(0) {}													// Constructor: initialize all to zeros
    Vec2(const T &bx, const T &by) : x(bx),y(by) {}							// Constructor: x, y
	Vec2(const Vec2<T> &b) { memcpy(this,&b,sizeof(*this)); }				// Constructor: copy
	template <typename T2>
	Vec2(const Vec2<T2>& b) { *this = b; }									// Constructor: copy (conversion)
	// Basic Operators
	const Vec2<T>& operator =(const Vec2<T> &b);							// Operator = (assignment)
	template <typename T2>
	const Vec2<T>& operator =(const Vec2<T2>& b);							// Operator = (assignment/conversion)
	T& operator [](unsigned int i) { return v[i]; }							// Operator []
	const T& operator [](unsigned int i)const { return v[i]; }				// Operator [] const
	bool operator ==(const Vec2<T> &b)const;								// Operator ==
	bool operator !=(const Vec2<T> &b)const { return !(*this == b); }		// Operator !=
	// Vector Operators
	Vec2<T> operator +(const Vec2<T> &b)const;								// Operator +
	Vec2<T> operator -(const Vec2<T> &b)const;								// Operator -
	Vec2<T> operator -()const;												// Operator - (inverse)
	T operator *(const Vec2<T> &b)const;									// Dot Product Operator *
	const Vec2<T>& operator +=(const Vec2<T> &b) { return *this = *this + b; }	// Operator +=
	const Vec2<T>& operator -=(const Vec2<T> &b) { return *this = *this - b; }	// Operator -=
	// Scalar Operators
	Vec2<T> operator *(T s)const;											// Operator * (scalar)
	template <typename sT> friend Vec2<sT> operator *(sT a, const Vec2<sT> &b);	// Operator * (scalar)
	Vec2<T> operator /(T s)const;											// Operator / (scalar)
	const Vec2<T>& operator *=(T s) { return *this = *this * s; }			// Operator *= (scalar)
	const Vec2<T>& operator /=(T s) { return *this = *this / s; }			// Operator /= (scalar)
	// Other Member Functions
	T len()const { return sqrt(x*x + y*y); }								// Length
	T lenSq()const { return x*x + y*y; }									// Length squared
	Vec2<T> norm()const { return lenSq() != 0 ? *this / len() : *this; }	// Get normalized vector
	const Vec2<T>& normalize() { return lenSq() != 0 ? *this /= len() : *this; }	// Normalize vector
	T angle()const { return atan2(y,x); }									// Angle of vector (from zero)
	T angle(const Vec2<T>& b)const;											// Angle between vectors
	bool similar(const Vec2<T> &b,T margin = 0.01)const;					// Check if vectors are similar
	const T* getData()const { return (const T*)(&v); }						// Get pointer to raw data
	// Static Member Functions
	static T angle(const Vec2<T> &a, const Vec2<T> &b) { return a.angle(b); }	// Angle between vectors (static)
	static T dot(const Vec2<T> &a, const Vec2<T> &b) { return a * b; }			// Dot Product (static)
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
    Vec3(const T &bx, const T &by, const T &bz) : x(bx), y(by), z(bz) {}	// Constructor: x, y, z
	Vec3(const Vec3<T> &b) { memcpy(this,&b,sizeof(*this)); }				// Constructor: copy
	template <typename T2>
	Vec3(const Vec3<T2>& b) { *this = b; }									// Constructor: copy (conversion)
    Vec3(const Vec2<T> &b, const T& s = 0);									// Constructor: 2D vector
	// Basic Operators
	const Vec3<T>& operator =(const Vec3<T> &b);							// Operator =
	template <typename T2>
	const Vec3<T>& operator =(const Vec3<T2>& b);							// Operator = (assignment/conversion)
	T& operator [](unsigned int i) { return v[i]; }							// Operator []
	const T& operator [](unsigned int i)const { return v[i]; }				// Operator [] const
	bool operator ==(const Vec3<T> &b)const;								// Operator ==
	bool operator !=(const Vec3<T> &b)const { return !(*this == b); }		// Operator !=
	// Vector Operators
	Vec3<T> operator +(const Vec3<T> &b)const;								// Operator +
	Vec3<T> operator -(const Vec3<T> &b)const;								// Operator -
	Vec3<T> operator -()const;												// Operator - (inverse)
	Vec3<T> operator %(const Vec3<T> &b)const;								// Cross Product Operator %
	T operator *(const Vec3<T> &b)const;									// Dot Product Operator *
	const Vec3<T>& operator +=(const Vec3<T> &b) { return *this = *this + b; }	// Operator +=
	const Vec3<T>& operator -=(const Vec3<T> &b) { return *this = *this - b; }	// Operator -=
	const Vec3<T>& operator %=(const Vec3<T> &b) { return *this = *this % b; }	// Cross Product Operator %=
	// Scalar Operators
	Vec3<T> operator *(T s)const;											// Operator * (scalar)
	template <typename sT> friend Vec3<sT> operator *(sT a, const Vec3<sT> &b);	// Operator * (scalar)
	Vec3<T> operator /(T s)const;											// Operator / (scalar)
	const Vec3<T>& operator *=(T s) { return *this = *this * s; }			// Operator *= (scalar)
	const Vec3<T>& operator /=(T s) { return *this = *this / s; }			// Operator /= (scalar)
	// Other Member Functions
	T len()const { return sqrt(x*x + y*y + z*z); }							// Length
	T lenSq()const { return x*x + y*y + z*z; }								// Length squared
	Vec3<T> norm()const { return lenSq() != 0 ? *this / len() : *this; }	// Get normalized vector
	const Vec3<T>& normalize() { return lenSq() != 0 ? *this /= len() : *this; }// Normalize vector
	T angle(const Vec3<T>& b)const;											// Angle between vectors
	bool similar(const Vec3<T> &b,T margin = 0.01)const;					// Check if vectors are similar
	Vec2<T> xy(){ return Vec2<T>(x, y); }									// Get XY components
	const T* getData()const { return (const T*)(&v); }						// Get pointer to raw data
	// Static Member Functions
	static T angle(const Vec3<T> &a, const Vec3<T> &b) { return a.angle(b); }	// Angle between vectors (static)
	static T cross(const Vec3<T> &a, const Vec3<T> &b) { return a % b; }		// Cross Product (static)
	static T dot(const Vec3<T> &a, const Vec3<T> &b) { return a * b; }			// Dot Product (static)
};

template <typename T = float>
class Vec4
{
public:
    // DATA
    union
    {
        struct { T x; T y; T z; T w;};
        T v[4];
    };

    // FUNCTIONS
    // Constructors
    Vec4() : x(0),y(0),z(0),w(0) {}											// Constructor: initialize all to zeros
    Vec4(const T &bx,const T &by,const T &bz, const T &bw) : x(bx),y(by),z(bz),w(bw) {}	    // Constructor: x, y, z, w
    Vec4(const Vec4<T> &b) { memcpy(this,&b,sizeof(*this)); }				// Constructor: copy
	template <typename T2>
	Vec4(const Vec4<T2>& b) { *this = b; }									// Constructor: copy (conversion)
    Vec4(const Vec3<T> &b,const T& s = 0);									// Constructor: 3D vector
    // Basic Operators
    const Vec4<T>& operator =(const Vec4<T> &b);							// Operator =
	template <typename T2>
	const Vec4<T>& operator =(const Vec4<T2>& b);							// Operator = (assignment/conversion)
    T& operator [](unsigned int i) { return v[i]; }							// Operator []
    const T& operator [](unsigned int i)const { return v[i]; }				// Operator [] const
    bool operator ==(const Vec4<T> &b)const;								// Operator ==
    bool operator !=(const Vec4<T> &b)const { return !(*this == b); }		// Operator !=
    // Vector Operators
    Vec4<T> operator +(const Vec4<T> &b)const;								// Operator +
    Vec4<T> operator -(const Vec4<T> &b)const;								// Operator -
    Vec4<T> operator -()const;												// Operator - (inverse)
    T operator *(const Vec4<T> &b)const;									// Dot Product Operator *
    const Vec4<T>& operator +=(const Vec4<T> &b) { return *this = *this + b; }	// Operator +=
    const Vec4<T>& operator -=(const Vec4<T> &b) { return *this = *this - b; }	// Operator -=
    // Scalar Operators
    Vec4<T> operator *(T s)const;											// Operator * (scalar)
    template <typename sT> friend Vec4<sT> operator *(sT a,const Vec4<sT> &b);	// Operator * (scalar)
    Vec4<T> operator /(T s)const;											// Operator / (scalar)
    const Vec4<T>& operator *=(T s) { return *this = *this * s; }			// Operator *= (scalar)
    const Vec4<T>& operator /=(T s) { return *this = *this / s; }			// Operator /= (scalar)
    // Other Member Functions
    T len()const { return sqrt(x*x + y*y + z*z + w*w); }					// Length
    T lenSq()const { return x*x + y*y + z*z + w*w; }						// Length squared
    Vec4<T> norm()const { return lenSq() != 0 ? *this / len() : *this; }	// Get normalized vector
    const Vec4<T>& normalize() { return lenSq() != 0 ? *this /= len() : *this; }	// Normalize vector
    bool similar(const Vec4<T> &b,T margin = 0.01)const;					// Check if vectors are similar
    Vec3<T> xyz(){ return Vec3<T>(x,y,z); }									// Get XYZ components
    const T* getData()const { return (const T*)(&v); }						// Get pointer to raw data
    // Static Member Functions
    static T dot(const Vec4<T> &a,const Vec4<T> &b) { return a * b; }		// Dot Product (static)
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



// CONSTRUCTORS
// *****************************************************************************************************************************

// Vec3
// Constructor: copy one size smaller vector plus one element
// Inputs	b = input vector
//			s = last element
template <typename T>
Vec3<T>::Vec3(const Vec2<T> &b, const T &s)
{
	memcpy(&v,&(b.v),2 * sizeof(T));
	v[2] = s;
}

// Vec4
// Constructor: copy one size smaller vector plus one element
// Inputs	b = input vector
//			s = last element
template <typename T>
Vec4<T>::Vec4(const Vec3<T> &b,const T& s)
{
    memcpy(&v,&(b.v),3 * sizeof(T));
    v[3] = s;
}


// BASIC OPERATORS
// *****************************************************************************************************************************

// Vec2
// Operator =
template <typename T>
const Vec2<T>& Vec2<T>::operator =(const Vec2<T> &b)
{
	x = b.x;	y = b.y;
	return *this;
}
// Vec2
// Operator = (assignment/conversion)
template <typename T>
template <typename T2>
const Vec2<T>& Vec2<T>::operator =(const Vec2<T2>& b)
{
	x = T(b.x);	y = T(b.y);
	return *this;
}
// Vec3
// Operator =
template <typename T>
const Vec3<T>& Vec3<T>::operator =(const Vec3<T> &b)
{
	x = b.x;	y = b.y;	z = b.z;
	return *this;
}
// Vec3
// Operator = (assignment/conversion)
template <typename T>
template <typename T2>
const Vec3<T>& Vec3<T>::operator =(const Vec3<T2>& b)
{
	x = T(b.x);	y = T(b.y);	z = T(b.z);
	return *this;
}
// Vec4
// Operator =
template <typename T>
const Vec4<T>& Vec4<T>::operator =(const Vec4<T> &b)
{
    x = b.x;	y = b.y;	z = b.z;    w = b.w;
    return *this;
}
// Vec4
// Operator = (assignment/conversion)
template <typename T>
template <typename T2>
const Vec4<T>& Vec4<T>::operator =(const Vec4<T2>& b)
{
	x = T(b.x);	y = T(b.y);	z = T(b.z);	w = T(b.w);
	return *this;
}

// Vec2
// Operator ==
template <typename T>
bool Vec2<T>::operator ==(const Vec2<T> &b)const
{
	return memcmp(this,&b,sizeof(*this)) == 0;
}
// Vec3
// Operator ==
template <typename T>
bool Vec3<T>::operator ==(const Vec3<T> &b)const
{
	return memcmp(this,&b,sizeof(*this)) == 0;
}
// Vec4
// Operator ==
template <typename T>
bool Vec4<T>::operator ==(const Vec4<T> &b)const
{
    return memcmp(this,&b,sizeof(*this)) == 0;
}

// VECTOR OPERATORS
// *****************************************************************************************************************************

// Vec2
// Operator +
template <typename T>
Vec2<T> Vec2<T>::operator +(const Vec2<T> &b)const
{
	return Vec2<T>(x + b.x, y + b.y);
}
// Vec3
// Operator +
template <typename T>
Vec3<T> Vec3<T>::operator +(const Vec3<T> &b)const
{
	return Vec3<T>(x + b.x, y + b.y, z + b.z);
}
// Vec4
// Operator +
template <typename T>
Vec4<T> Vec4<T>::operator +(const Vec4<T> &b)const
{
    return Vec4<T>(x + b.x,y + b.y,z + b.z,w + b.w);
}

// Vec2
// Operator -
template <typename T>
Vec2<T> Vec2<T>::operator -(const Vec2<T> &b)const
{
	return Vec2<T>(x - b.x, y - b.y);
}
// Vec3
// Operator -
template <typename T>
Vec3<T> Vec3<T>::operator -(const Vec3<T> &b)const
{
	return Vec3<T>(x - b.x, y - b.y, z - b.z);
}
// Vec4
// Operator -
template <typename T>
Vec4<T> Vec4<T>::operator -(const Vec4<T> &b)const
{
    return Vec4<T>(x - b.x,y - b.y,z - b.z,w - b.w);
}

// Vec2
// Operator - (inverse)
template <typename T>
Vec2<T> Vec2<T>::operator -()const
{
	return Vec2<T>(-x, -y);
}
// Vec3
// Operator - (inverse)
template <typename T>
Vec3<T> Vec3<T>::operator -()const
{
	return Vec3<T>(-x, -y, -z);
}
// Vec4
// Operator - (inverse)
template <typename T>
Vec4<T> Vec4<T>::operator -()const
{
    return Vec4<T>(-x,-y,-z,-w);
}

// Vec3 (only)
// Cross Product Operator %
template <typename T>
Vec3<T> Vec3<T>::operator %(const Vec3<T> &b)const
{
	return Vec3<T>(
		y * b.z - z * b.y,  //x component result
		z * b.x - x * b.z,  //y component result
		x * b.y - y * b.x); //z component result
}

// Vec2
// Dot Product Operator *
template <typename T>
T Vec2<T>::operator *(const Vec2<T> &b)const
{
	return x * b.x + y * b.y;
}
// Vec3
// Dot Product Operator *
template <typename T>
T Vec3<T>::operator *(const Vec3<T> &b)const
{
	return x * b.x + y * b.y + z * b.z;
}
// Vec4
// Dot Product Operator *
template <typename T>
T Vec4<T>::operator *(const Vec4<T> &b)const
{
    return x * b.x + y * b.y + z * b.z + w * b.w;
}


// SCALAR OPERATORS
// *****************************************************************************************************************************

// Vec2
// Operator * (scalar)
template <typename T>
Vec2<T> Vec2<T>::operator *(T s)const
{
	return Vec2<T>(x * s, y * s);
}
// Vec3
// Operator * (scalar)
template <typename T>
Vec3<T> Vec3<T>::operator *(T s)const
{
    return Vec3<T>(x * s,y * s,z * s);
}
// Vec4
// Operator * (scalar)
template <typename T>
Vec4<T> Vec4<T>::operator *(T s)const
{
    return Vec4<T>(x * s,y * s,z * s,w * s);
}

// Vec2 - Non-Member Friend Operator
// Operator * (scalar * Vec2)
template <typename T>
Vec2<T> operator *(T s,const Vec2<T> &b)
{
	return Vec2<T>(b.x * s,b.y * s);
}
// Vec3 - Non-Member Friend Operator
// Operator * (scalar * Vec3)
template <typename T>
Vec3<T> operator *(T s,const Vec3<T> &b)
{
	return Vec3<T>(b.x * s,b.y * s, b.z * s);
}
// Vec4 - Non-Member Friend Operator
// Operator * (scalar * Vec4)
template <typename T>
Vec4<T> operator *(T s,const Vec4<T> &b)
{
    return Vec4<T>(b.x * s,b.y * s,b.z * s,w * s);
}


// Vec2
// Operator / (scalar)
template <typename T>
Vec2<T> Vec2<T>::operator /(T s)const
{
	return Vec2<T>(x / s, y / s);
}
// Vec3
// Operator / (scalar)
template <typename T>
Vec3<T> Vec3<T>::operator /(T s)const
{
	return Vec3<T>(x / s, y / s, z / s);
}
// Vec4
// Operator / (scalar)
template <typename T>
Vec4<T> Vec4<T>::operator /(T s)const
{
    return Vec4<T>(x / s,y / s,z / s,w / s);
}


// OTHER MEMBER FUNCTIONS
// *****************************************************************************************************************************

// Vec2
// Angle between vectors
template <typename T>
T Vec2<T>::angle(const Vec2<T>& b)const
{
	T den = len() * b.len();
	if (den != T(0))
	{
		T num = *this * b;
		T div = num / den;
		
		if (div <= T(1) && div >= T(-1))
			return acos(div);
		else if (div > T(1))
			// TODO: exception handling
			return acos(1);
		else
			// TODO: exception handling
			return acos(-1);
	}
	else
	{
		// divide-by-zero
		// TODO: exception handling
		return 0;
	}
} 

// Vec3
// Angle between vectors
template <typename T>
T Vec3<T>::angle(const Vec3<T>& b)const
{
	T den = len() * b.len();
	if (den != T(0))
	{
		T num = *this * b;
		T div = num / den;

		if (div <= T(1) && div >= T(-1))
			return acos(div);
		else if (div > T(1))
			// TODO: exception handling
			return acos(1);
		else
			// TODO: exception handling
			return acos(-1);
	}
	else
	{
		// TODO: exception handling
		return 0;
	}
}


// Vec2
// Check if vectors are within a certain margin of being equal
template <typename T>
bool Vec2<T>::similar(const Vec2<T> &b ,T margin)const
{
	return ((x - b.x <= margin) && (y - b.y <= margin));
}
// Vec3
// Check if vectors are within a certain margin of being equal
template <typename T>
bool Vec3<T>::similar(const Vec3<T> &b,T margin)const
{
	return ((x - b.x <= margin) && (y - b.y <= margin) && (z - b.z <= margin));
}
// Vec4
// Check if vectors are within a certain margin of being equal
template <typename T>
bool Vec4<T>::similar(const Vec4<T> &b,T margin)const
{
    return ((x - b.x <= margin) && (y - b.y <= margin) && (z - b.z <= margin) && (w - b.w <= margin));
}

#endif
