// *****************************************************************************************************************************
// quat.cpp
// Linear Math Libraries
// Quaternion Libraries
// Template definition for data type
// Author(s): Cory Douthat
// Copyright (c) 2026 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

// CONSTRUCTORS
// *****************************************************************************************************************************
#ifdef VEC_HPP_
// Constructor: Axis-Angle
template <typename T>
Quat<T>::Quat(T angle, const Vec3<T> &axis)
{
	if (angle != T(0) && axis != Vec3<T>(0, 0, 0))
	{
		Vec3<T> r = axis.norm();
		T temp = sin(angle / 2);

		a = cos(angle / 2);
		b = temp * r.x;
		c = temp * r.y;
		d = temp * r.z;
	}
	else
	{
		a = 1; 
		b = 0; 
		c = 0; 
		d = 0;
	}
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

// Operator = (assignment/conversion)
template <typename T>
template <typename T2>
const Quat<T>& Quat<T>::operator =(const Quat<T2>& qb)
{
	a = T(qb.a);
	b = T(qb.b);
	c = T(qb.c);
	d = T(qb.d);

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
	// TODO: conj() only if unit quaternion:
	// TODO: check for this by tracking in a member variable?
	return conj();

	// TODO: (check) - Else non-unit quaternion
	//return Quat(a, -b, -c, -d) / (a*a + b*b + c*c + d*d);
}

// Calculate exponential (base e) of Quat
template <typename T>
Quat<T> Quat<T>::qexp()const
{
	// https://math.stackexchange.com/questions/939229/unit-quaternion-to-a-scalar-power

	// Check for zero quaternion / denominator == 0
	T len = norm();
	T len_xyz = sqrt(b * b + c * c + d * d);

	if (len == 0 || len_xyz == 0)
		return Quat<T>(1, 0, 0, 0);

	Quat<T> sgn = *this / len;

	return (Quat<T>(cos(len_xyz), 0, 0, 0) * exp(a) + sgn * sin(len_xyz));
}

// Calculate logarithm (base e) of Quat
template <typename T>
Quat<T> Quat<T>::qln()const
{
	// https://math.stackexchange.com/questions/939229/unit-quaternion-to-a-scalar-power
	
	T len = norm();
	T len_xyz = sqrt(x * x + y * y + z * z);

	// ln(0) is undefined, but best we can do is return a zero quat?
	if (len == 0 || len_xyz == 0)
		return Quat<T>(1, 0, 0, 0);

	T qa = log(len);

	Vec3<T> qv = (Vec3<T>(x, y, z) / len_xyz) * acos(a / len);

	return Quat<T>(qa, qv.x, qv.y, qv.z);
}

// Raise Quat to a power x
template <typename T>
Quat<T> Quat<T>::qpow(T x)const
{
	// https://math.stackexchange.com/questions/939229/unit-quaternion-to-a-scalar-power
	// TODO: Is this optimized?

	return (qln() * x).qexp();
}

// Norm of Quat ||q||
template <typename T>
T Quat<T>::norm()const
{
	// TODO: more efficient way?
	return sqrt(a*a + b*b + c*c + d*d);
}

// Get unit Quat q/||q||
template <typename T>
Quat<T> Quat<T>::unit()const 
{ 
	T norm = this->norm();

	if (norm == 0)
		return Quat<T>();	// Identity quaternion
	else
		return *this / norm; 
}

// Unitize this Quat
template <typename T>
const Quat<T>& Quat<T>::unitize() 
{ 
	T norm = this->norm();

	if (norm == 0)
		return *this = Quat<T>();	// Identity quaternion
	else
		return *this = *this / norm; 
}

// Rotate vector in Quat format
template <typename T>
Quat<T> Quat<T>::rotate(const Quat<T> &q2)const
{
	// TODO: Check if multiplier is a unit quaternion?
	// TODO: Optimize?
	return Quat<T>(*this * q2 * this->inv());
}

// Check if Quats are similar
template <typename T>
bool Quat<T>::similar(const Quat<T> &b, T margin)const
{
	return ((a - b.a <= margin) && (b - b.b <= margin) && (c - b.c <= margin) && (d - b.d <= margin));
}

// Slerp (interpolation) between two Quats
template <typename T>
Quat<T> Quat<T>::slerp(const Quat<T>& q1, const Quat<T>& q2, T t)
{
	// TODO: Need to check if they're unit quaternions?
	if (t > 1.0)
		t = 1.0;

	// Math results in zero quat if both are equal - not sure if that's correct
	if (t <= 0 || q1 == q2)
		return q1;

	if (t > 1)
		return q2;

	// https://en.wikipedia.org/wiki/Slerp#Quaternion_Slerp
	// TODO: may need to handle unitization more wholistically (not just this function) for rotation quats
	return (q1 * (q1.inv() * q2).qpow(t)).unitize();
}

#ifdef VEC_HPP_
// Rotate vector by quat
template <typename T>
Vec3<T> Quat<T>::rotate(const Vec3<T> &p)const 
{ 
	// TODO: Check if multiplier is a unit quaternion
	// TODO: Optimize?
	// TODO: Instead use formula from Wikipedia under "Used Methods": vnew = v + 2r x (r x v + wv)
	// https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
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
	T t0 = T(2.0) * (a * b + c * d);
	T t1 = T(1.0) - T(2.0) * (b * b + csqr);
	temp.x = std::atan2(t0, t1);

	// pitch (y-axis rotation)
	T t2 = T(2.0) * (a * c - d * b);
	t2 = t2 > T(1.0) ? T(1.0) : t2;
	t2 = t2 < T(-1.0) ? T(-1.0) : t2;
	temp.y = std::asin(t2);

	// yaw (z-axis rotation)
	T t3 = T(2.0) * (a * d + b * c);
	T t4 = T(1.0) - T(2.0) * (csqr + d * d);
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
