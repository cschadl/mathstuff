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

#include <assert.h>

#include "misc.h"

extern void test_vector();

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
	n_vector(const std::valarray<T>& va) : m_v(va) { assert(m_v.size() == N); }

	enum { Dim = N };

public:
	n_vector() : m_v((T)0, N) { }
	n_vector(T x, T y);				// 2-vector
	n_vector(T x, T y, T z);		// 3-vector
	n_vector(T a, T b, T c, T d);	// 4-vector
	virtual ~n_vector() { }

	/** The number of elements in the vector **/
	int size() const { return Dim; }

	virtual T length() const;					/** length **/
	virtual T length_sq() const;				/** length squared **/

	virtual n_vector<T, N>  make_unit() const;	/** make a unit vector from this vector **/
	virtual n_vector<T, N>& unit();				/** unitize this vector **/

	virtual T norm(unsigned int p) const;
	virtual T max_norm() const;
	virtual T inner_product(const n_vector<T, N>& rhs) const;

	T& operator[](unsigned int i) { return m_v[i]; }
	const T& operator[](unsigned int i) const { return m_v[i]; }

	bool operator==(const n_vector<T, N>& rhs) const;
	bool operator!=(const n_vector<T, N>& rhs) const { return !(*this == rhs); }

	bool is_close(const n_vector<T, N>& v, const T tol) const;

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

	virtual n_vector<T, N>& operator+=(const n_vector<T, N>& v) { m_v += v.m_v; return *this; }
	virtual n_vector<T, N>& operator-=(const n_vector<T, N>& v) { m_v -= v.m_v; return *this; }
	virtual n_vector<T, N>& operator*=(const T c) { m_v *= c; return *this; }
	virtual n_vector<T, N>& operator/=(const T d) { m_v /= d; return *this; }

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
	for (int i = 0 ; i < Dim ; i++)
		if (m_v[i] != rhs.m_v[i])
			return false;

	return true;
}

template <typename T, int N>
bool n_vector<T, N>::is_close(const n_vector<T, N>& v, const T tol) const
{
	// I guess this could get hairy if T is not scalar...
	for (int i = 0 ; i < Dim ; i++)
		if (!close(m_v[i], v.m_v[i], tol))
			return false;

	return true;
}

template <typename T, int N>
T n_vector<T, N>::length() const
{
	return norm(2);
}

template <typename T, int N>
T n_vector<T, N>::length_sq() const
{
	T length_sq = 0;
	for (int i = 0 ; i < N ; i++)
		length_sq += (m_v[i] * m_v[i]);

	return length_sq;
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

	T p_norm = 0;
	for (int i = 0 ; i < Dim ; i++)
		p_norm += fabs(pow(m_v[i], (T) p));

	p_norm = pow(p_norm, 1.0 / (T) p);

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
	const U n1 = (v1[3] * v2[0]) - (v1[0] * v2[2]);
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
	for (int i = 0 ; i < v.size() ; i++)
		os << v[i] << " ";
	os << ")";

	return os;
}

};	// namespace maths

#endif /* VECTORS_H_ */