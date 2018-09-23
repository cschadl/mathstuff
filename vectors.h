/*
 * vectors.h
 *
 *  Created on: Feb 26, 2012
 *      Author: cds
 */

#ifndef VECTORS_H_
#define VECTORS_H_

#include <valarray>
#include <cmath>
#include <iostream>
#include <limits>

#include <assert.h>

#include "misc.h"

namespace maths
{

/**
 * A vector of N elements
 */
template <typename T, int N>
class n_vector
{
protected:
	std::valarray<T> m_v;	// vector elements

public:
	static const size_t Dim = N;
	typedef T value_type;

public:
	n_vector() : m_v((T)0, N) { }
	n_vector(T x, T y);				// 2-vector
	n_vector(T x, T y, T z);		// 3-vector
	n_vector(T a, T b, T c, T d);	// 4-vector
	n_vector(const std::valarray<T>& va) : m_v(va) { assert(m_v.size() == N); }
	n_vector(const T* data, size_t n) : m_v(data, n) { }

	n_vector(const n_vector<T, N>&) = default;
	n_vector<T, N>& operator=(const n_vector<T, N>&) = default;

	n_vector(n_vector<T, N>&&) = default;
	n_vector<T, N>& operator=(n_vector<T, N>&&) = default;

	~n_vector() { }

	/** The number of elements in the vector **/
	size_t size() const { return Dim; }

	T length() const;								/** length **/
	T length_sq() const;							/** length squared (faster, but less percise than length() (?) **/
	T distance_sq(const n_vector<T, N> & q) const;	/** The squared distance from this point to the point q */

	n_vector<T, N>  make_unit() const;	/** make a unit vector from this vector **/
	n_vector<T, N>& unit();				/** unitize this vector **/

	T norm(unsigned int p) const;
	T max_norm() const;
	T inner_product(const n_vector<T, N>& rhs) const;

	// TODO - bounds checking (std::valarray fails it)
	T& operator[](unsigned int i) { return m_v[i]; }
	const T& operator[](unsigned int i) const { return m_v[i]; }

	const T& x() const 	{ return m_v[0]; }
	      T& x()		{ return m_v[0]; }
	const T& y() const 	{ return m_v[1]; }
		  T& y()		{ return m_v[1]; }
	const T& z() const	{ return m_v[2]; }
		  T& z()		{ return m_v[2]; }

	bool operator==(const n_vector<T, N>& rhs) const;
	bool operator!=(const n_vector<T, N>& rhs) const { return !(*this == rhs); }

	bool is_close(const n_vector<T, N>& v, const T tol) const;
	bool is_null() const;
	bool is_unit() const;

	// Arithmetic operations
	template<typename U, int M> friend n_vector<U, M> operator+(const n_vector<U, M>& v1, const n_vector<U, M>& v2);
	template<typename U, int M> friend n_vector<U, M> operator-(const n_vector<U, M>& v1, const n_vector<U, M>& v2);
	template<typename U, int M> friend U              operator*(const n_vector<U, M>& v1, const n_vector<U, M>& v2);	// dot
	template<typename U, int M> friend n_vector<U, M> operator*(const n_vector<U, M>& v, const U c);
	template<typename U, int M> friend n_vector<U, M> operator*(const U c, const n_vector<U, M>& v);
	template<typename U, int M> friend n_vector<U, M> operator/(const n_vector<U, M>& v, const U d);
	template<typename U, int M> friend n_vector<U, M> operator/(const U d, const n_vector<U, M>& v);
	template<typename U, int M> friend n_vector<U, M> operator-(const n_vector<U, M>& v);	// negation
	template<typename U, int M> friend n_vector<U, M> operator%(const n_vector<U, M>& v1, const n_vector<U, M>& v2);	// cross
	template<typename U, int M> friend n_vector<U, M> outer_product(const n_vector<U, M>& v1, const n_vector<U, M>& v2);
	template<typename U, int M> friend n_vector<U, M> cross(const n_vector<U, M>& v1, const n_vector<U, M>& v2);

	n_vector<T, N>& operator+=(const n_vector<T, N>& v) { m_v += v.m_v; return *this; }
	n_vector<T, N>& operator-=(const n_vector<T, N>& v) { m_v -= v.m_v; return *this; }
	n_vector<T, N>& operator*=(const T c) { m_v *= c; return *this; }
	n_vector<T, N>& operator/=(const T d) { m_v /= d; return *this; }

	const std::valarray<T>& get_valarray() const { return m_v; }

	// unit vectors
	// these don't work very well (need to be called like vector3d::n_vector<double>i()
//	template <typename U> static const n_vector<U, 3> & i();
//	template <typename U> static const n_vector<U, 3> & j();
//	template <typename U> static const n_vector<U, 3> & k();
};

typedef n_vector<double, 2> vector2d;
typedef n_vector<double, 3> vector3d;
typedef n_vector<float,  2> vector2f;
typedef n_vector<float,  3> vector3f;

//template<> template<>
//const n_vector<double, 3>& n_vector<double, 3>::i()
//{
//	static n_vector<double, 3> ux(1.0, 0.0, 0.0);
//	return ux;
//}
//
//template<> template<>
//const n_vector<double, 3>& n_vector<double, 3>::j()
//{
//	static n_vector<double, 3> uy(0.0, 1.0, 0.0);
//	return uy;
//}
//
//template<> template<>
//const n_vector<double, 3>& n_vector<double, 3>::k()
//{
//	static n_vector<double, 3> uz(0.0, 0.0, 1.0);
//	return uz;
//}

template <typename T, int N>
n_vector<T, N>::n_vector(T x, T y)
: m_v((T)0, N)
{
	assert(size() == 2);
	m_v[0] = x;
	m_v[1] = y;
}

template <typename T, int N>
n_vector<T, N>::n_vector(T x, T y, T z)
: m_v((T)0, N)
{
	assert(size() == 3);
	m_v[0] = x;
	m_v[1] = y;
	m_v[2] = z;
}

template <typename T, int N>
n_vector<T, N>::n_vector(T a, T b, T c, T d)
: m_v((T)0, N)
{
	assert(size() == 4);
	m_v[0] = a;
	m_v[1] = b;
	m_v[2] = c;
	m_v[3] = d;
}

template <typename T, int N>
bool n_vector<T, N>::operator==(const n_vector<T, N>& rhs) const
{
	for (size_t i = 0 ; i < Dim ; i++)
		if (m_v[i] != rhs.m_v[i])
			return false;

	return true;
}

template <typename T, int N>
bool n_vector<T, N>::is_close(const n_vector<T, N>& v, const T tol) const
{
	return maths::close(distance_sq(v), 0.0, tol * tol);
}

template <typename T, int N>
bool n_vector<T, N>::is_null() const
{
	const T tol = std::numeric_limits<T>::epsilon();	// should be OK for comparing to 0
	return is_close(n_vector<T, N>(), tol);
}

template <typename T, int N>
bool n_vector<T, N>::is_unit() const
{
	return close(length(), T(1.0), std::numeric_limits<T>::epsilon());
}

template <typename T, int N>
T n_vector<T, N>::length() const
{
	return norm(2);
}

template <typename T, int N>
T n_vector<T, N>::length_sq() const
{
	const n_vector<T, N> & v = *this;
	return v * v;
}

template <typename T, int N>
T n_vector<T, N>::distance_sq(const n_vector<T, N>& q) const
{
	return ((*this) - q).length_sq();
}

template <typename T, int N>
n_vector<T, N> n_vector<T, N>::make_unit() const
{
	return n_vector<T, N>(*this).unit();
}

template <typename T, int N>
n_vector<T, N>& n_vector<T, N>::unit()
{
	m_v /= norm(2);
	return *this;
}

template <typename T, int N>
T n_vector<T, N>::norm(unsigned int p) const
{
	// TODO - p must not be 0

	T p_norm = (T)0;
	for (size_t i = 0 ; i < Dim ; i++)
		p_norm += fabs(pow(m_v[i], (T) p));

	p_norm = pow(p_norm, (T)1 / (T) p);

	return p_norm;
}

template <typename T, int N>
T n_vector<T, N>::max_norm() const
{
	// TODO - inefficient
	std::valarray<T> m_copy = abs(m_v);
	return m_copy.max();
}

template <typename T, int N>
T n_vector<T, N>::inner_product(const n_vector<T, N>& rhs) const
{
	T scalar_product = 0;
	for (int i = 0 ; i < N ; i++)
		scalar_product += (m_v[i] * rhs.m_v[i]);

	return scalar_product;
}

template <typename U, int M>
inline n_vector<U, M> operator+(const n_vector<U, M>& v1, const n_vector<U, M>& v2)
{
	return n_vector<U, M>(v1.m_v + v2.m_v);
}

template <typename U, int M>
inline n_vector<U, M> operator-(const n_vector<U, M>& v1, const n_vector<U, M>& v2)
{
	return n_vector<U, M>(v1.m_v - v2.m_v);
}

template <typename U, int M>
inline U operator*(const n_vector<U, M>& v1, const n_vector<U, M>& v2)
{
	return v1.inner_product(v2);
}

template <typename U, int M>
inline n_vector<U, M>operator*(const n_vector<U, M>& v, const U c)
{
	return n_vector<U, M>(v.m_v * c);
}

template <typename U, int M>
inline n_vector<U, M> operator*(const U c, const n_vector<U, M>& v)
{
	return v * c;
}

template <typename U, int M>
inline n_vector<U, M> operator/(const n_vector<U, M>& v, const U d)
{
	return n_vector<U, M>(v.m_v / d);
}

template <typename U, int M>
inline n_vector<U, M> operator/(const U d, const n_vector<U, M>& v)
{
	return n_vector<U, M>(d / v.m_v);
}

template <typename U, int M>
inline n_vector<U, M> operator-(const n_vector<U, M>& v)
{
	return n_vector<U, M>(-v.m_v);
}

// For now, the outer product only has partial specializations for
// particular dimensions
template <typename U>
inline
n_vector<U, 3> outer_product(const n_vector<U, 3>& v1, const n_vector<U, 3>& v2)
{
	const U n0 = (v1[1] * v2[2]) - (v1[2] * v2[1]);
	const U n1 = (v1[2] * v2[0]) - (v1[0] * v2[2]);
	const U n2 = (v1[0] * v2[1]) - (v1[1] * v2[0]);

	return n_vector<U, 3>(n0, n1, n2);
}

template <typename U>
inline
n_vector<U, 3> cross(const n_vector<U, 3>& v1, const n_vector<U, 3>& v2)
{
	return outer_product(v1, v2);
}

template <typename U>
inline
n_vector<U, 3> operator%(const n_vector<U, 3>& v1, const n_vector<U, 3>& v2)
{
	return outer_product(v1, v2);
}

template <typename U, int M>
inline std::ostream& operator<<(std::ostream& os, const n_vector<U, M>& v)
{
	os << "( ";
	for (size_t i = 0 ; i < v.size() ; i++)
		os << v[i] << " ";
	os << ")";

	return os;
}

template <typename TO, typename FROM, int M>
inline n_vector<TO, M> convert(const n_vector<FROM, M>& src)
{
	n_vector<TO, M> out;

	for (size_t i = 0 ; i < M ; i++)
		out[i] = (TO) src[i];

	return out;
}

// Specialization of maths::abs<T>() for n_vector<double, 3>
// How can we do this for the more general case of n_vector<T, n>?
template <>
inline n_vector<double, 3> abs(n_vector<double, 3> v)
{
	return n_vector<double, 3>(abs(v.x()), abs(v.y()), abs(v.z()));
}

};	// namespace maths

#endif /* VECTORS_H_ */
