// *****************************************************************************************************************************
// quat.h
// Linear Math Libraries
// Quaternion Libraries
// Template definition for data type
// Author(s): Cory Douthat
// Copyright (c) 2026 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef QUAT_HPP_
#define QUAT_HPP_

#define _USE_MATH_DEFINES	// For PI definition
#include <cmath>
#include <cstring>

#include "vec.h"
#include "mat.h"

template <typename T>
class Mat2;
template <typename T>
class Mat3;
template <typename T>
class Mat4;

// TODO: Check / maintain unit quaternion?

template <typename T = float>
class Quat
{
public:
    // DATA
    union
    {
        struct { T w; T x; T y; T z; };
        struct { T a; T b; T c; T d; };
        //struct { T a; T i; T j; T k; };
        T v[4];
    };

    // FUNCTIONS
    // Constructors
    Quat() : a(1), b(0), c(0), d(0) {}											// Constructor: initialize to identity
    Quat(const T& a2, const T& b2, const T& c2, const T& d2) : a(a2), b(b2), c(c2), d(d2) {}	    // Constructor: a, bi, cj, dk
    Quat(const Quat<T>& b2) { memcpy(this, &b2, sizeof(*this)); }			// Constructor: copy
    template <typename T2>
    Quat(const Quat<T2>& b2) { *this = b2; }								// Constructor: copy (conversion)

    // Basic Operators
    const Quat<T>& operator =(const Quat<T>& b);							// Operator =
    template <typename T2>
    const Quat<T>& operator =(const Quat<T2>& b);							// Operator = (assignment/conversion)
    T& operator [](unsigned int i) { return v[i]; }							// Operator []
    const T& operator [](unsigned int i)const { return v[i]; }				// Operator [] const
    bool operator ==(const Quat<T>& b)const;								// Operator ==
    bool operator !=(const Quat<T>& b)const { return !(*this == b); }		// Operator !=
    // Quat Operators
    Quat<T> operator +(const Quat<T>& q2)const;								// Operator +
    Quat<T> operator -(const Quat<T>& q2)const;								// Operator -
    Quat<T> operator *(const Quat<T>& q2)const;								// Multiplication (non-commutative)
    //const Quat<T>& operator +=(const Quat<T> &q2) { return *this = *this + q2; }	// Operator +=
    //const Quat<T>& operator -=(const Quat<T> &a2) { return *this = *this - q2; }	// Operator -=
    // Scalar Operators
    // DANGEROUS - THESE DO NOT SCALE ROTATION
    Quat<T> operator *(T s)const;											// Operator * (scalar)
    Quat<T> operator /(T s)const;											// Operator / (scalar)
    // const Quat<T>& operator *=(T s) { return *this = *this * s; }			// Operator *= (scalar)
    // const Quat<T>& operator /=(T s) { return *this = *this / s; }			// Operator /= (scalar)
    // Other Member Functions
    Quat<T> conj()const;													// Conjugate of Quat
    Quat<T> inv()const;														// Inverse of Quat
    Quat<T> qexp()const;													// Calculate exponential (base e) of Quat
    Quat<T> qln()const;														// Calculate logarithm (base e) of Quat
    Quat<T> qpow(T x)const;													// Raise Quat to a power x
    T norm()const;															// Norm of Quat ||q||
    Quat<T> unit()const;													// Get unit Quat q/||q||
    const Quat<T>& unitize();												// Unitize Quat
    Quat<T> rotate(const Quat<T>& q2)const;									// Rotate vector in Quat format
    bool similar(const Quat<T>& b, T margin = 0.01)const;					// Check if Quats are similar
    const T* getData()const { return (const T*)(&v); }						// Get pointer to raw data

    static Quat<T> slerp(const Quat<T>& q1, const Quat<T>& q2, T t);		// Slerp (interpolation) between two Quats

    // Functions relating to vectors (Vec3)
//#ifdef VEC_HPP_
    Quat(T angle, const Vec3<T>& axis);										// Constructor: Axis-Angle
    Quat(const Vec4<T>& b) { memcpy(this->v, b.getData(), sizeof(*this)); }	// Constructor: 4D vector
    Vec3<T> operator *(const Vec3<T>& vec)const { return rotate(vec); }		// Operator * (Vec3)
    Vec3<T> rotate(const Vec3<T>& p)const;									// Rotate vector by quat
    Vec3<T> toEuler()const;													// To Euler angles (scaled Vec3)
    static Quat<T> fromEuler(const Vec3<T>& v);								// Generate from Euler angles (scaled Vec3)
//#endif

    // Functions relating to Matrices
//#ifdef MAT_HPP_
    //Quat(const Mat3<T> &b);												// Constructor: Mat3
    //Quat(const Mat4<T> &b);												// Constructor: Mat4
    Mat3<T> rotMat()const { return rotMat3(); }							    // Extract rotation matrix
    Mat3<T> rotMat3()const;													// Extract rotation matrix
    Mat4<T> rotMat4()const;													// Extract rotation matrix
//#endif
};

// Shortcut Types
typedef Quat<float> Quatf;
typedef Quat<double> Quatd;
typedef Quat<int> Quati;
typedef Quat<short> Quats;

#endif