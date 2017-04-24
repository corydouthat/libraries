// *****************************************************************************************************************************
// quat.hpp
// Linear Math Libraries
// Quaternion Libraries
// Template definition for data type
// Author(s): Cory Douthat
// Copyright (c) 2017 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef QUAT_HPP_
#define QUAT_HPP_

#define _USE_MATH_DEFINES	// For PI definition
#include <cmath>
#include <cstring>

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
		struct { T a; T i; T j; T k; };
        T v[4];
    };

    // FUNCTIONS
    // Constructors
    Quat() : a(1),b(0),c(0),d(0) {}											// Constructor: initialize to identity
    Quat(const T &a,const T &b,const T &c, const T &d) : a(a),b(b),c(c),d(d) {}	    // Constructor: a, bi, cj, dk
    Quat(const Quat<T> &b) { memcpy(this, &b, sizeof(*this)); }				// Constructor: copy
	
    // Basic Operators
    const Quat<T>& operator =(const Quat<T> &b);							// Operator =
    T& operator [](unsigned int i) { return v[i]; }							// Operator []
    const T& operator [](unsigned int i)const { return v[i]; }				// Operator [] const
    bool operator ==(const Quat<T> &b)const;								// Operator ==
    bool operator !=(const Quat<T> &b)const { return !(*this == b); }		// Operator !=
    // Quat Operators
    Quat<T> operator +(const Quat<T> &q2)const;								// Operator +
    Quat<T> operator -(const Quat<T> &q2)const;								// Operator -
    Quat<T> operator *(const Quat<T> &q2)const;									// Multiplication (non-commutative)
    //const Quat<T>& operator +=(const Quat<T> &q2) { return *this = *this + q2; }	// Operator +=
    //const Quat<T>& operator -=(const Quat<T> &a2) { return *this = *this - q2; }	// Operator -=
    // Scalar Operators
    Quat<T> operator *(T s)const;											// Operator * (scalar)
    template <typename sT> friend Quat<sT> operator *(sT s,const Quat<sT> &q);	// Operator * (scalar)
    Quat<T> operator /(T s)const;											// Operator / (scalar)
    const Quat<T>& operator *=(T s) { return *this = *this * s; }			// Operator *= (scalar)
    const Quat<T>& operator /=(T s) { return *this = *this / s; }			// Operator /= (scalar)
    // Other Member Functions
	Quat<T> conj()const;													// Conjugate of Quat
	Quat<T> inv()const;														// Inverse of Quat
	T norm()const;															// Norm of Quat ||q||
	Quat<T> unit()const { return *this / this->norm(); }					// Get unit Quat q/||q||
	const Quat<T>& unitize() { return *this = *this / this->norm(); }		// Unitize Quat
	Quat<T> rotate(const Quat<T> &q2)const;									// Rotate Quat
    bool similar(const Quat<T> &b,T margin = 0.01)const;					// Check if Quats are similar
    const T* getData()const { return (const T*)(&v); }						// Get pointer to raw data

// Functions relating to vectors (Vec3)
#ifdef VEC_HPP_
	Quat(const T angle, Vec3<T> &axis);										// Constructor: Axis-Angle
	Quat(const Vec4<T> &b) { memcpy(this->v, b.getData(), sizeof(*this)); }	// Constructor: 4D vector
	Vec3<T> rotate(const Vec3<T> &p)const;									// Rotate vector by quat
	Vec3<T> toEuler()const;													// To Euler angles (scaled Vec3)
	static Quat<T> fromEuler(const Vec3<T> &v);								// Generate from Euler angles (scaled Vec3)
#endif

// Functions relating to Matrices
#ifdef MAT_HPP_
	//Quat(const Mat3<T> &b);													// Constructor: Mat3
	//Quat(const Mat4<T> &b);													// Constructor: Mat4
	Mat3<T> rotMat()const { return rotMat3(); }								// Extract rotation matrix
	Mat3<T> rotMat3()const;													// Extract rotation matrix
	Mat4<T> rotMat4()const;													// Extract rotation matrix
#endif
};

// Shortcut Types
typedef Quat<float> Quatf;
typedef Quat<double> Quatd;
typedef Quat<int> Quati;
typedef Quat<short> Quats;

// CONSTRUCTORS
// *****************************************************************************************************************************
#ifdef VEC_HPP_
// Constructor: Axis-Angle
template <typename T>
Quat<T>::Quat(const T angle, Vec3<T> &axis)
{
	if (angle != T(0) && axis != Vec3f(0, 0, 0))
	{
		Vec3<T> r = axis.norm();
		T temp = sin(angle / 2);
		
		a = cos(angle / 2);
		b = temp * r.x;
		c = temp * r.y;
		d = temp * r.z;
	}
	else
		Quat();
}

#endif

#ifdef MAT_HPP_
//// Constructor: Mat3
//template <typename T>
//Quat<T>::Quat(const Mat3<T> &b)
//{
//	// TODO
//}
//
//// Constructor: Mat4
//template <typename T>
//Quat<T>::Quat(const Mat4<T> &b)
//{
//	// TODO
//}
#endif

// BASIC OPERATORS
// *****************************************************************************************************************************
// Operator =
template <typename T>
const Quat<T>& Quat<T>::operator =(const Quat<T> &b)
{
	memcpy(this, &b, sizeof(*this));
	return *this;
}

// Operator ==
template <typename T>
bool Quat<T>::operator ==(const Quat<T> &b)const
{
	return memcmp(this, &b, sizeof(*this)) == 0;
}


// QUATERNION OPERATORS
// *****************************************************************************************************************************
// Operator +
template <typename T>
Quat<T> Quat<T>::operator +(const Quat<T> &q2)const
{
	return Quat<T>(a + q2.a, b + q2.b, c + q2.c, d + q2.d);
}

// Operator -
template <typename T>
Quat<T> Quat<T>::operator -(const Quat<T> &q2)const
{
	return Quat<T>(a - q2.a, b - q2.b, c - q2.c, d - q2.d);
}

// Multiplication (non-commutative)
template <typename T>
Quat<T> Quat<T>::operator *(const Quat<T> &q2)const
{
	Quat<T> temp;

	temp.v[0] = v[0]*q2.v[0] - v[1]*q2.v[1] - v[2]*q2.v[2] - v[3]*q2.v[3];
	temp.v[1] = v[0]*q2.v[1] + v[1]*q2.v[0] + v[2]*q2.v[3] - v[3]*q2.v[2];
	temp.v[2] = v[0]*q2.v[2] - v[1]*q2.v[3] + v[2]*q2.v[0] + v[3]*q2.v[1];
	temp.v[3] = v[0]*q2.v[3] + v[1]*q2.v[2] - v[2]*q2.v[1] + v[3]*q2.v[0];

	return temp;
}

// SCALAR OPERATORS
// *****************************************************************************************************************************
// Operator * (scalar)
template <typename T>
Quat<T> Quat<T>::operator *(T s)const
{
	return Quat<T>(a*s, b*s, c*s, d*s);
}

// Non-Member Operator * (scalar * Quat)
template <typename T> 
Quat<T> operator *(T s, const Quat<T> &q)
{
	return Quat<T>(q.a*s, q.b*s, q.c*s, q.d*s);
}

// Operator / (scalar)
template <typename T>
Quat<T> Quat<T>::operator /(T s)const
{
	return Quat<T>(a/s, b/s, c/s, d/s);
}

// OTHER MEMBER FUNCTIONS
// *****************************************************************************************************************************
// Conjugate of Quat
template <typename T>
Quat<T> Quat<T>::conj()const
{
	return Quat<T>(a, -b, -c, -d);
}

// Inverse of Quat
template <typename T>
Quat<T> Quat<T>::inv()const
{
	// If unit quaternion
	return conj();

	// TODO: (check) - Else non-unit quaternion
	//return Quat(a, -b, -c, -d) / (a*a + b*b + c*c + d*d);
}

// Norm of Quat ||q||
template <typename T>
T Quat<T>::norm()const
{
	// TODO: more efficient way?
	return *this * this->inv();
}

// Rotate Quat
template <typename T>
Quat<T> Quat<T>::rotate(const Quat<T> &q2)const
{
	// TODO: Check if multiplier is a unit quaternion?
	// TODO: Optimize?
	return Quat<T>(*this * q2 * this->inv());
}

// Check if Quats are similar
template <typename T>
bool Quat<T>::similar(const Quat<T> &b, T margin = 0.01)const
{
	return ((a - b.a <= margin) && (b - b.b <= margin) && (c - b.c <= margin) && (d - b.d <= margin));
}

#ifdef VEC_HPP_
// Rotate vector by quat
template <typename T>
Vec3<T> Quat<T>::rotate(const Vec3<T> &p)const 
{ 
	// TODO: Check if multiplier is a unit quaternion
	// TODO: Optimize?
	Quat<T> pq = Quat<T>(0, p.x, p.y, p.z);
	Quat<T> temp = rotate(pq);
	return Vec3<T>(temp.b, temp.c, temp.d);
}

// To Euler angles (scaled Vec3)
template <typename T>
Vec3<T> Quat<T>::toEuler()const
{
	Vec3<T> temp;
	// Based on code from Wikipedia article
	// https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
	T csqr = c * c;

	// roll (x-axis rotation)
	T t0 = +2 * (a * b + c * d);
	T t1 = +1 - 2 * (b * b + csqr);
	temp.x = std::atan2(t0, t1);

	// pitch (y-axis rotation)
	T t2 = +2 * (a * c - d * b);
	t2 = t2 > 1 ? 1 : t2;
	t2 = t2 < -1 ? -1 : t2;
	temp.y = std::asin(t2);

	// yaw (z-axis rotation)
	T t3 = +2 * (a * d + b * c);
	T t4 = +1 - 2 * (csqr + d * d);
	temp.z = std::atan2(t3, t4);

	return temp;
}

// Generate Quat from Euler angles (scaled Vec3)
template <typename T>
Quat<T> Quat<T>::fromEuler(const Vec3<T> &v)
{
	T angle = v.len();
	if (angle != 0)
		return Quat<T>(angle, v / angle);
	else
		return Quat<T>();
}

#endif

#ifdef MAT_HPP_
// Extract rotation matrix (Mat3)
template <typename T>
Mat3<T> Quat<T>::rotMat3()const
{
	return Mat3<T>(
		a*a + b*b - c*c - d*d,	2*b*c - 2*a*d,			2*b*d + 2*a*c,
		2*b*c + 2*a*d,			a*a - b*b + c*c - d*d,	2*c*d - 2*a*b,
		2*b*d - 2*a*c,			2*c*d + 2*a*b,			a*a - b*b - c*c + d*d
		);
}

// Extract rotation matrix (Mat4)
template <typename T>
Mat4<T> Quat<T>::rotMat4()const
{
	return Mat4<T>(
		a*a + b*b - c*c - d*d,	2 * b*c - 2 * a*d,		2 * b*d + 2 * a*c,		0,
		2 * b*c + 2 * a*d,		a*a - b*b + c*c - d*d,	2 * c*d - 2 * a*b,		0,
		2 * b*d - 2 * a*c,		2 * c*d + 2 * a*b,		a*a - b*b - c*c + d*d,	0,
		0,						0,						0,						1
		);
}
#endif

#endif
