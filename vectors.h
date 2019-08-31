/*
 * vectors.h
 *
 *  Created on: Feb 26, 2012
 *      Author: cds
 */

#ifndef VECTORS_H_
#define VECTORS_H_

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>

#include <assert.h>

#include <stlutil/type_traits.h>
#include <stlutil/make_array.h>

#include "misc.h"

namespace maths
{


/**
 * A vector of N elements
 */
template <typename T, size_t N>
class n_vector
{
protected:
	using elements_t = std::array<T, N>;

	elements_t m_v;	// vector elements

public:
	static const size_t Dim = N;
	typedef T value_type;

public:
	/// This constructor is for lvalues
	template <	typename... Vals,
					std::enable_if_t<stlutil::conjunction<
						std::is_reference<Vals>...,
						std::is_same<T, std::decay_t<Vals>>...>::value, void*> = nullptr>
	n_vector(Vals&&... vals)
		: m_v{std::forward<Vals>(vals)...}
	{
		// We'll still get a compile-time error (too many initializers) if
		// we have more than N arguments, so just check for too few
		static_assert(sizeof...(vals) >= N, "Too few arguments");
	}

	/// And this is for rvalues
	template <	typename... Vals,
					std::enable_if_t<	stlutil::conjunction<
						std::is_arithmetic<Vals>...,
						std::is_convertible<T, Vals>...>::value, void*> = nullptr>
		n_vector(Vals&&... vals)
			: m_v{static_cast<T>(std::forward<Vals>(vals))...}	// static_cast<T> avoids 'narrowing int to double' warning...
		{
			// We'll still get a compile-time error (too many initializers) if
			// we have more than N arguments, so just check for too few
			static_assert(sizeof...(vals) >= N, "Too few arguments");
		}

	n_vector()
		: m_v(stlutil::make_array_val<T, N>(T(0)))
	{

	}

//	n_vector(T const& v)
//		: m_v(make_array_val<T, N>(v))
//	{
//
//	}

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
	template<typename U, size_t M> friend n_vector<U, M> operator+(const n_vector<U, M>& v1, const n_vector<U, M>& v2);
	template<typename U, size_t M> friend n_vector<U, M> operator-(const n_vector<U, M>& v1, const n_vector<U, M>& v2);
	template<typename U, size_t M> friend U              operator*(const n_vector<U, M>& v1, const n_vector<U, M>& v2);	// dot
	template<typename U, size_t M> friend n_vector<U, M> operator*(const n_vector<U, M>& v, U const& c);
	template<typename U, size_t M> friend n_vector<U, M> operator*(U const& c, const n_vector<U, M>& v);
	template<typename U, size_t M> friend n_vector<U, M> operator/(const n_vector<U, M>& v, U const& d);
	template<typename U, size_t M> friend n_vector<U, M> operator/(U const& d, const n_vector<U, M>& v);
	template<typename U, size_t M> friend n_vector<U, M> operator-(const n_vector<U, M>& v);	// negation
	template<typename U, size_t M> friend n_vector<U, M> operator%(const n_vector<U, M>& v1, const n_vector<U, M>& v2);	// cross
	template<typename U, size_t M> friend n_vector<U, M> outer_product(const n_vector<U, M>& v1, const n_vector<U, M>& v2);
	template<typename U, size_t M> friend n_vector<U, M> cross(const n_vector<U, M>& v1, const n_vector<U, M>& v2);

	n_vector<T, N>& operator+=(const n_vector<T, N>& v)
	{
		for (size_t i = 0 ; i < N ; i++)
			m_v[i] += v.m_v[i];

		return *this;
	}
	n_vector<T, N>& operator-=(const n_vector<T, N>& v)
	{
		for (size_t i = 0 ; i < N ; i++)
			m_v[i] -= v.m_v[i];

		return *this;
	}
	n_vector<T, N>& operator*=(T c)
	{
		for (size_t i = 0 ; i < N ; i++)
			m_v[i] *= c;

		return *this;
	}
	n_vector<T, N>& operator/=(T d)
	{
		for (size_t i = 0 ; i < N ; i++)
			m_v[i] /= d;

		return *this;
	}

	const elements_t& get_elements() const { return m_v; }

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

template <typename T, size_t N>
bool n_vector<T, N>::operator==(const n_vector<T, N>& rhs) const
{
	for (size_t i = 0 ; i < Dim ; i++)
		if (m_v[i] != rhs.m_v[i])
			return false;

	return true;
}

template <typename T, size_t N>
bool n_vector<T, N>::is_close(const n_vector<T, N>& v, const T tol) const
{
	return maths::close(distance_sq(v), 0.0, tol * tol);
}

template <typename T, size_t N>
bool n_vector<T, N>::is_null() const
{
	const T tol = std::numeric_limits<T>::epsilon();	// should be OK for comparing to 0
	return is_close(n_vector<T, N>(), tol);
}

template <typename T, size_t N>
bool n_vector<T, N>::is_unit() const
{
	return close(length(), T(1.0), std::numeric_limits<T>::epsilon());
}

template <typename T, size_t N>
T n_vector<T, N>::length() const
{
	return norm(2);
}

template <typename T, size_t N>
T n_vector<T, N>::length_sq() const
{
	const n_vector<T, N> & v = *this;
	return v * v;
}

template <typename T, size_t N>
T n_vector<T, N>::distance_sq(const n_vector<T, N>& q) const
{
	return ((*this) - q).length_sq();
}

template <typename T, size_t N>
n_vector<T, N> n_vector<T, N>::make_unit() const
{
	return n_vector<T, N>(*this).unit();
}

template <typename T, size_t N>
n_vector<T, N>& n_vector<T, N>::unit()
{
	*this /= norm(2);
	return *this;
}

template <typename T, size_t N>
T n_vector<T, N>::norm(unsigned int p) const
{
	// TODO - p must not be 0

	T p_norm = (T)0;
	for (size_t i = 0 ; i < Dim ; i++)
		p_norm += fabs(pow(m_v[i], (T) p));

	p_norm = pow(p_norm, (T)1 / (T) p);

	return p_norm;
}

template <typename T, size_t N>
T n_vector<T, N>::max_norm() const
{
	// TODO - inefficient
	elements_t m_copy = abs(m_v);
	return m_copy.max();
}

template <typename T, size_t N>
T n_vector<T, N>::inner_product(const n_vector<T, N>& rhs) const
{
	T scalar_product = 0;
	for (size_t i = 0 ; i < N ; i++)
		scalar_product += (m_v[i] * rhs.m_v[i]);

	return scalar_product;
}

template <typename U, size_t M>
inline n_vector<U, M> operator+(const n_vector<U, M>& v1, const n_vector<U, M>& v2)
{
	n_vector<U, M> res;
	for (size_t i = 0 ; i < M ; i++)
		res[i] = v1[i] + v2[i];

	return res;
}

template <typename U, size_t M>
inline n_vector<U, M> operator-(const n_vector<U, M>& v1, const n_vector<U, M>& v2)
{
	n_vector<U, M> res;
	for (size_t i = 0 ; i < M ; i++)
		res[i] = v1[i] - v2[i];

	return res;
}

template <typename U, size_t M>
inline U operator*(const n_vector<U, M>& v1, const n_vector<U, M>& v2)
{
	return v1.inner_product(v2);
}

template <typename U, size_t M>
inline n_vector<U, M> operator*(const n_vector<U, M>& v, U const& c)
{
	n_vector<U, M> res;
	for (size_t i = 0 ; i < M ; i++)
		res[i] = v[i] * c;

	return res;
}

template <typename U, size_t M>
inline n_vector<U, M> operator*(U const& c, n_vector<U, M> const& v)
{
	n_vector<U, M> res;
	for (size_t i = 0 ; i < M ; i++)
		res[i] = c * v[i];

	return res;
}

template <typename U, size_t M>
inline n_vector<U, M> operator/(const n_vector<U, M>& v, U const& c)
{
	n_vector<U, M> res;
	for (size_t i = 0 ; i < M ; i++)
		res[i] = v[i] / c;

	return res;
}

template <typename U, size_t M>
inline n_vector<U, M> operator/(U const& c, n_vector<U, M> const& v)
{
	n_vector<U, M> res;
	for (size_t i = 0 ; i < M ; i++)
		res[i] = c / v[i];

	return res;
}

template <typename U, size_t M>
inline n_vector<U, M> operator-(const n_vector<U, M>& v)
{
	n_vector<U, M> res;
	for (size_t i = 0 ; i < M ; i++)
		res[i] = -(res[i]);

	return res;
}

// For now, the outer product only has partial specializations for
// particular dimensions
template <typename U>
inline
n_vector<U, 3> outer_product(const n_vector<U, 3>& v1, const n_vector<U, 3>& v2)
{
	return n_vector<U, 3>((v1[1] * v2[2]) - (v1[2] * v2[1]),
								 (v1[2] * v2[0]) - (v1[0] * v2[2]),
								 (v1[0] * v2[1]) - (v1[1] * v2[0]));
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

template <typename U, size_t M>
inline std::ostream& operator<<(std::ostream& os, const n_vector<U, M>& v)
{
	os << "( ";
	for (size_t i = 0 ; i < v.size() ; i++)
		os << v[i] << " ";
	os << ")";

	return os;
}

template <typename TO, typename FROM, size_t M>
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
