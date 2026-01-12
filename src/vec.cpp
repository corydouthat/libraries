// *****************************************************************************************************************************
// Vec.cpp
// Linear Math Libraries
// 2D And 3D vector classes
// Template definition for data type
// Author(s): Cory Douthat
// Copyright (c) 2026 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************


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

// Vec2 (special)
template <typename T>
T Vec2<T>::operator %(const Vec2<T>& b)const
{
	return x * b.y - y * b.x;	// equivalent of 3D z result
}

// Vec3
// Cross Product Operator %
template <typename T>
Vec3<T> Vec3<T>::operator %(const Vec3<T> &b)const
{
	return Vec3<T>(
		y * b.z - z * b.y,		// x component result
		z * b.x - x * b.z,		// y component result
		x * b.y - y * b.x);		// z component result
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
template <typename T, typename sT>
Vec2<T> operator *(sT s, const Vec2<T> &b)
{
	return Vec2<T>(b.x * s,b.y * s);
}
// Vec3 - Non-Member Friend Operator
// Operator * (scalar * Vec3)
template <typename T, typename sT>
Vec3<T> operator *(sT s, const Vec3<T> &b)
{
	return Vec3<T>(b.x * (T)s, b.y * T(s), b.z * (T)s);
}
// Vec4 - Non-Member Friend Operator
// Operator * (scalar * Vec4)
template <typename T, typename sT>
Vec4<T> operator *(sT s, const Vec4<T> &b)
{
    return Vec4<T>(b.x * s, b.y * s, b.z * s, b.w * s);
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
			return acos((T)1.0);
		else
			// TODO: exception handling
			return acos((T)(-1.0));
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

// Vec3
// Find perpendicular unit vector
template <typename T>
Vec3<T> Vec3<T>::perpendicular(const Vec3<T>& a, const Vec3<T>& b)
{
	Vec3<T> temp = a % b;

	if (temp != Vec3<T>(T(0), T(0), T(0)))
		return temp.norm();
	else
	{
		// a and b are same direction or opposite
		// Any vector on plane of 'a' will work
		if ((temp = a % Vec3<T>(0.0, 1.0, 0.0)) != Vec3<T>(0.0, 0.0, 0.0))
			return temp.norm();
		else
			return (a % Vec3<T>(0.0, 0.0, 1.0)).norm();
	}
}
