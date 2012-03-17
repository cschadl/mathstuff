/*
 * quaternion.h
 *
 *  Created on: Mar 11, 2012
 *      Author: cds
 */

#ifndef QUATERNION_H_
#define QUATERNION_H_

#include "vectors.h"
#include "matrix.h"

namespace maths
{

template <typename T>
class quaternion
{
protected:
	T				m_scalar;	// scalar part of quaternion
	n_vector<T, 3>	m_vector;	// vector part of quaternion

public:
	explicit quaternion();
	quaternion(const T& a, const T& b, const T& c, const T& d);
	quaternion(const T& r, const n_vector<T, 3>& v);

	bool operator==(const quaternion<T>& rhs) const;
	bool operator!=(const quaternion<T>& rhs) const { return !((*this) == rhs); }
	bool operator==(const T& a) const; // just for convenience
	bool is_close(const quaternion<T>& rhs, const T& tol) const;
	bool is_close(const T& a, const T& tol) const;

	const T& operator[](size_t n) const; // 0 is scalar, 1-3 are i, j, k
	T& operator[](size_t n);

	const n_vector<T, 3>& vector_part() const { return m_vector; }
	const T& scalar_part() const { return m_scalar; }

	static quaternion<T> i() { return quaternion<T>((T)0, (T)1, (T)0, (T)0); }
	static quaternion<T> j() { return quaternion<T>((T)0, (T)0, (T)1, (T)0); }
	static quaternion<T> k() { return quaternion<T>((T)0, (T)0, (T)0, (T)1); }

	const T& a()  const { return m_scalar; }
	const T& b()  const { return m_vector[0]; }
	const T& c()  const { return m_vector[1]; }
	const T& d()  const { return m_vector[2]; }
	const T& bi() const { return m_vector[0] * i(); }
	const T& cj() const { return m_vector[1] * j(); }
	const T& dk() const { return m_vector[2] * k(); }

	quaternion<T>  get_conj() const;
	quaternion<T>& conj();
	quaternion<T>  get_unit() const;
	quaternion<T>& unit();
	quaternion<T>  get_reciprocal() const;
	quaternion<T>& reciprocal();
	T norm() const;
	T norm_sq() const;

	/** get the 4x4 rotation matrix for this quaternion **/
	//matrix<T> rotation_matrix() const;
};

// fwd declarations of binary operators
template <typename _T> quaternion<_T> operator+(const quaternion<_T>& q1, const quaternion<_T>& q2);
template <typename _T> quaternion<_T> operator+(const quaternion<_T>& q, const _T& a);
template <typename _T> quaternion<_T> operator+(const _T& a, const quaternion<_T>& q);
template <typename _T> quaternion<_T> operator-(const quaternion<_T>& q1, const quaternion<_T>& q2);
template <typename _T> quaternion<_T> operator-(const quaternion<_T>& q, const _T& a);
template <typename _T> quaternion<_T> operator-(const _T& a, const quaternion<_T>& q);
template <typename _T> quaternion<_T> operator*(const quaternion<_T>& q1, const quaternion<_T>& q2);
template <typename _T> quaternion<_T> operator*(const quaternion<_T>& q, const _T& c);
template <typename _T> quaternion<_T> operator*(const _T& c, const quaternion<_T>& q);
template <typename _T> quaternion<_T> operator/(const quaternion<_T>& q, const _T& d);
template <typename _T> quaternion<_T> operator-(const quaternion<_T>& q);

typedef quaternion<double> quaternion_d;
typedef quaternion<double> quaternion_f;

};

#include "quaternion.hpp"

#endif /* QUATERNION_H_ */
