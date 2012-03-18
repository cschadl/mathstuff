/*
 * quaternion.hpp
 *
 *  Created on: Mar 11, 2012
 *      Author: cds
 */

namespace maths
{

template <typename T>
quaternion<T>::quaternion()
: m_scalar((T)0)
{

}

template <typename T>
quaternion<T>::quaternion(const T& a, const T& b, const T& c, const T& d)
: m_scalar(a)
, m_vector(b, c, d)
{

}

template <typename T>
quaternion<T>::quaternion(const T& r, const n_vector<T, 3>& v)
: m_scalar(r)
, m_vector(v)
{

}

template <typename T>
bool quaternion<T>::operator==(const quaternion<T>& rhs) const
{
	return m_scalar == rhs.m_scalar && m_vector == rhs.m_vector;
}

template <typename T>
bool quaternion<T>::operator==(const T& c) const
{
	return m_scalar == c && m_vector.is_null();
}

template <typename T>
bool quaternion<T>::is_close(const quaternion<T>& rhs, const T& tol) const
{
	return close(scalar_part(), rhs.scalar_part(), tol) && vector_part().is_close(rhs.vector_part(), tol);
}

template <typename T>
bool quaternion<T>::is_close(const T& a, const T& tol) const
{
	return close(scalar_part(), a, tol) && vector_part().is_null(tol);
}

template <typename T>
const T& quaternion<T>::operator[](size_t n) const
{
	return n == 0 ? m_scalar : m_vector[n - 1];
}

template <typename T>
T& quaternion<T>::operator[](size_t n)
{
	return n == 0 ? m_scalar : m_vector[n - 1];
}

template <typename T>
quaternion<T> quaternion<T>::get_conj() const
{
	return quaternion<T>(m_scalar, -m_vector[0], -m_vector[1], -m_vector[2]);
}

template <typename T>
quaternion<T>& quaternion<T>::conj()
{
	m_vector *= -1.0;
	return *this;
}

template <typename T>
quaternion<T> quaternion<T>::get_unit() const
{
	return (*this) / norm();
}

template <typename T>
quaternion<T>& quaternion<T>::unit()
{
	// xxx - poo
	*this = get_unit();
	return *this;
}

template <typename T>
quaternion<T> quaternion<T>::get_reciprocal() const
{
	return conj() / norm_sq();
}

template <typename T>
quaternion<T>& quaternion<T>::reciprocal()
{
	*this = get_reciprocal();
	return *this;
}

template <typename T>
T quaternion<T>::norm() const
{
	return sqrt((a()*a()) + (b()*b()) + (c()*c()) + (d()*d()));
}

template <typename T>
T quaternion<T>::norm_sq() const
{
	return (a()*a()) + (b()*b()) + (c()*c()) + (d()*d());
}

///////////////
// operators
///////////////

template <typename _T>
quaternion<_T> operator+(const quaternion<_T>& q1, const quaternion<_T>& q2)
{
	return quaternion<_T>(q1.scalar_part() + q2.scalar_part(), q1.vector_part() + q2.vector_part());
}

template <typename _T>
quaternion<_T> operator+(const quaternion<_T>& q1, const _T& a)
{
	return quaternion<_T>(a + q1.scalar_part(), q1.vector_part());
}

template <typename _T>
quaternion<_T> operator+(const _T& a, const quaternion<_T>& q1)
{
	return quaternion<_T>(a + q1.scalar_part(), q1.vector_part());
}

template <typename _T>
quaternion<_T> operator-(const quaternion<_T>& q1, const quaternion<_T>& q2)
{
	return quaternion<_T>(q1.scalar_part() - q2.scalar_part(), q1.vector_part() - q2.vector_part());
}

template <typename _T>
quaternion<_T> operator-(const quaternion<_T>& q1, const _T& a)
{
	return quaternion<_T>(q1.scalar_part() - a, q1.vector_part());
}

template <typename _T>
quaternion<_T> operator-(const _T& a, const quaternion<_T>& q1)
{
	return quaternion<_T>(a - q1.scalar_part(), -q1.vector_part());
}

template <typename _T>
quaternion<_T> operator*(const quaternion<_T>& q1, const quaternion<_T>& q2)
{
	const _T r1r2 = (q1.scalar_part() * q2.scalar_part()) - (q1.vector_part() * q2.vector_part());
	const n_vector<_T, 3> v1v2 =
			  (q1.scalar_part() * q2.vector_part())
			+ (q2.scalar_part() * q1.vector_part())
			+ (q1.vector_part() % q2.vector_part());

	return quaternion<_T>(r1r2, v1v2);
}

template <typename _T>
quaternion<_T> operator*(const quaternion<_T>& q, const _T& c)
{
	return quaternion<_T>(q.scalar_part() * c, q.vector_part() * c);
}

template <typename _T>
quaternion<_T> operator*(const _T&c, const quaternion<_T>& q)
{
	return quaternion<_T>(q.scalar_part() * c, q.vector_part() * c);
}

template <typename _T>
quaternion<_T> operator*(const quaternion<_T>& q, const n_vector<_T, 3>& v)
{
	return quaternion<_T>(q.scalar_part(), q.vector_part() * v);
}

template <typename _T>
quaternion<_T> operator*(const n_vector<_T, 3>& v, const quaternion<_T>& q)
{
	return q * v;
}

template <typename _T>
quaternion<_T> operator/(const quaternion<_T>& q, const _T& c)
{
	return quaternion<_T>(q.scalar_part() / c, q.vector_part() / c);
}

template <typename _T>
quaternion<_T> operator-(const quaternion<_T>& q)
{
	return quaternion<_T>(-q.scalar_part(), -q.vector_part());
}

};
