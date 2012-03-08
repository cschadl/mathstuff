/*
 * matrix.h
 *
 *  Created on: Feb 26, 2012
 *      Author: cds
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <valarray>
#include <iostream>
#include <string>
#include <sstream>
#include <memory>

#include <assert.h>

#include "misc.h"

extern void test_matrix();

namespace maths
{

class matrix_exception;
template <typename T> class matrix_index_exception;
template <typename T> class matrix_dimension_mismatch_exception;
class invalid_matrix_exception;

template <typename T>
class matrix
{
protected:
	// Index helper interface
	// Allows 1-based indexing and instant transposition
	class indexer
	{
	protected:
		const matrix<T>* mp;
	public:
		indexer(const matrix<T>* m) : mp(m) { }
		const matrix<T>* get_mp() const { return mp; }

		virtual size_t index(size_t i, size_t j) const = 0;
		virtual size_t row_index_start() const = 0;
		virtual size_t row_index_end() const = 0;
		virtual size_t col_index_start() const = 0;
		virtual size_t col_index_end() const = 0;
		virtual indexer* clone(const matrix<T>* m) const = 0;	// clone this indexer for the passed in matrix
		virtual indexer* make_transpose() const = 0;	// returns the transpose indexer of this indexer type
		virtual indexer* make_zero() const = 0;			// returns the zero-based indexer of this indexer type
		virtual indexer* make_one() const = 0;			// returns the one-based indexer of this index type
		virtual ~indexer() { }
	};

protected:
	const size_t					m_n_rows;
	const size_t					m_n_cols;
	std::valarray<T>				m_A;	// matrix entries
	mutable std::auto_ptr<indexer>	m_idx;	// should be a unique_ptr.  never NULL.

	virtual size_t _idx(size_t i, size_t j) const {  return m_idx->index(i, j); }

public:
	explicit matrix(size_t m, size_t n);

	matrix(const matrix<T>& rhs);
	matrix& operator=(const matrix<T>& rhs);

	size_t rows() const { return m_n_rows; }
	size_t cols() const { return m_n_cols; }
	size_t n_entries() const { return m_n_rows * m_n_cols; }
	size_t r_begin() const { return m_idx->row_index_start(); }
	size_t r_end() const { return m_idx->row_index_end(); }
	size_t c_begin() const { return m_idx->col_index_start(); }
	size_t c_end() const { return m_idx->col_index_end(); }

	bool operator==(const matrix<T>& rhs) const;
	bool operator!=(const matrix<T>& rhs) const { return !(*this == rhs); }
	bool is_close(const matrix<T>& rhs, const T& tol) const;

	matrix<T>& transpose();
	matrix<T>  get_transpose() const;

	/** returns a nxn identity matrix **/
	static matrix<T> I(size_t n);

	/** create a nxn matrix whose diagonal entries are w **/
	static matrix<T> diag(const std::valarray<T>& w);

	/**
	 * Compute the singular value decomposition A = U*W*Vt of this matrix.
	 * The entries of this matrix are replaced with U, w is a vector whose
	 * entries are the diagonal elements of the m x m matrix W.
	 * V is returned rather than Vt.
	 * (Adapted from Numerical Recipes in C, sec. 2.6)
	 */
	matrix<T>& svd(std::valarray<T>& w, matrix<T>& V);

	const T& operator()(size_t i, size_t j) const;
	T& operator()(size_t i, size_t j);

	template <typename _T>
	friend matrix<_T> operator*(const matrix<_T>& a, const matrix<_T>& b);

	template <typename _T>
	friend matrix<_T> operator*(const matrix<_T>& m, const _T& c);

	template <typename _T>
	friend matrix<_T> operator*(const _T& c, const matrix<_T>& m);

	template <typename _T>
	friend matrix<_T> operator+(const matrix<_T>& a, const matrix<_T>& b);

	template <typename _T>
	friend matrix<_T> operator-(const matrix<_T>& a, const matrix<_T>& b);

	template <typename _T>
	friend matrix<_T> operator-(const matrix<_T>& m);	// negation

protected:
	// indexers
	class row_major_indexer : public indexer
	{
	public:
		// We need to use this-> when referencing indexer members because of how G++
		// does class member name lookups (see http://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html)
		row_major_indexer(const matrix<T>* m) : indexer(m) { }
		virtual size_t index(size_t i, size_t j) const { return this->mp->m_n_cols * i + j; }
		virtual size_t row_index_start() const { return 0; }
		virtual size_t row_index_end() const { return this->mp->m_n_rows; }
		virtual size_t col_index_start() const { return 0; }
		virtual size_t col_index_end() const { return this->mp->m_n_cols; }
		virtual indexer* clone(const matrix<T>* m) const { return new row_major_indexer(m); }
		virtual indexer* make_transpose() const { return new col_major_indexer(this->mp); }
		virtual indexer* make_zero() const { return clone(this->mp); }
		virtual indexer* make_one() const { return new row_major_1_indexer(this->mp); }
	};

	class col_major_indexer : public indexer
	{
	public:
		col_major_indexer(const matrix<T>* m) : indexer(m) { }
		virtual size_t index(size_t i, size_t j) const { return i + (this->mp->m_n_rows * j); }
		virtual size_t row_index_start() const { return 0; }
		virtual size_t row_index_end() const { return this->mp->m_n_rows; }
		virtual size_t col_index_start() const { return 0; }
		virtual size_t col_index_end() const { return this->mp->m_n_cols; }
		virtual indexer* clone(const matrix<T>* m) const { return new col_major_indexer(m); }
		virtual indexer* make_transpose() const { return new row_major_indexer(this->mp); }
		virtual indexer* make_zero() const { return clone(this->mp); }
		virtual indexer* make_one() const { return new col_major_1_indexer(this->mp); }
	};

	class row_major_1_indexer : public indexer
	{
	public:
		row_major_1_indexer(const matrix<T>* m) : indexer(m) { }
		virtual size_t index(size_t i, size_t j) const { return this->mp->m_n_cols * (i - 1) + (j - 1); }
		virtual size_t row_index_start() const { return 1; }
		virtual size_t row_index_end() const { return this->mp->m_n_rows + 1; }
		virtual size_t col_index_start() const { return 1; }
		virtual size_t col_index_end() const { return this->mp->m_n_cols + 1; }
		virtual indexer* clone(const matrix<T>* m) const { return new row_major_1_indexer(m); }
		virtual indexer* make_transpose() const { return new col_major_1_indexer(this->mp); }
		virtual indexer* make_zero() const { return new row_major_indexer(this->mp); }
		virtual indexer* make_one() const { return clone(this->mp); }
	};

	class col_major_1_indexer : public indexer
	{
	public:
		col_major_1_indexer(const matrix<T>* m) : indexer(m) { }
		virtual size_t index(size_t i, size_t j) const { return (i - 1) + (this->mp->m_n_rows * (j - 1)); }
		virtual size_t row_index_start() const { return 1; }
		virtual size_t row_index_end() const { return this->mp->m_n_rows + 1; }
		virtual size_t col_index_start() const { return 1; }
		virtual size_t col_index_end() const { return this->mp->m_n_cols + 1; }
		virtual indexer* clone(const matrix<T>* m) const { return new col_major_1_indexer(m); }
		virtual indexer* make_transpose() const { return new row_major_1_indexer(this->mp); }
		virtual indexer* make_zero() const { return new col_major_indexer(this->mp); }
		virtual indexer* make_one() const { return clone(this->mp); }
	};
public:
	// for testing
	enum IndexType
	{
		RowMajor,
		ColMajor,
		RowMajorOne,
		ColMajorOne,
	};

	void set_index_type(const enum IndexType t)
	{
		switch (t)
		{
		case RowMajor:
			m_idx = std::auto_ptr<indexer>(new row_major_indexer(this));
			break;
		case ColMajor:
			m_idx = std::auto_ptr<indexer>(new col_major_indexer(this));
			break;
		case RowMajorOne:
			m_idx = std::auto_ptr<indexer>(new row_major_1_indexer(this));
			break;
		case ColMajorOne:
			m_idx = std::auto_ptr<indexer>(new col_major_1_indexer(this));
			break;
		default:
			assert(false);
		}
	}

protected:

	class scoped_index_change
	{
	protected:
		const matrix<T>& 		mp;
		std::auto_ptr<indexer>	orig_index;

	public:
		scoped_index_change(const matrix<T>& m, indexer* new_index)
		: mp(m)
		, orig_index(m.m_idx.release())
		{
			assert(new_index->get_mp() == &mp);
			mp.m_idx = std::auto_ptr<indexer>(new_index);
		}

		~scoped_index_change()
		{
			mp.m_idx = std::auto_ptr<indexer>(orig_index.release());
		}
	};
};

template <typename T>
matrix<T>::matrix(size_t m, size_t n)
: m_n_rows(m)
, m_n_cols(n)
, m_A((T)0, m * n)
, m_idx(new row_major_indexer(this))
{
	// we'll probably run into problems before we get here
	// when we try to allocate m_A, but oh well
	if (m_n_rows < 1 || m_n_cols < 1)
		throw invalid_matrix_exception();
}

template <typename T>
matrix<T>::matrix(const matrix<T>& rhs)
: m_n_rows(rhs.m_n_rows)
, m_n_cols(rhs.m_n_cols)
, m_A(rhs.m_A)
, m_idx(rhs.m_idx->clone(this))	// inherit rhs index type
{
	assert(m_idx.get() != NULL);	// wtf doesnt !! work?
	assert(m_n_rows >= 1 && m_n_cols >= 1);
}

template <typename T>
matrix<T>& matrix<T>::operator=(const matrix<T>& rhs)
{
	if (m_n_rows != rhs.m_n_rows || m_n_cols != rhs.m_n_cols)
		throw matrix_dimension_mismatch_exception<T>(*this, rhs);

	m_A = rhs.m_A;

	return *this;
}

template <typename T>
bool matrix<T>::operator==(const matrix<T>& rhs) const
{
	if (m_n_rows != rhs.m_n_rows || m_n_cols != rhs.m_n_cols)
		return false;

	for (size_t i = 0 ; i < m_n_rows ; i++)
		for (size_t j = 0 ; j < m_n_cols ; j++)
			if ((*this)(i, j) != rhs(i, j))
				return false;

	return true;
}

template <typename T>
bool matrix<T>::is_close(const matrix<T>& rhs, const T& tol) const
{
	if (m_n_rows != rhs.m_n_rows || m_n_cols != rhs.m_n_cols)
		throw matrix_dimension_mismatch_exception<T>(*this, rhs);

	// we can't just compare m_A - what if we have two nxn
	// matrices with the same entries of m_A, but one is
	// row major, and the other is column major?

	// make sure we're using 0-based indices
	scoped_index_change ai(*this, m_idx->make_zero());
	scoped_index_change bi(rhs, rhs.m_idx->make_zero());

	for (size_t i = 0 ; i < m_n_rows ; i++)
		for (size_t j = 0 ; j < m_n_cols ; j++)
			if (!close((*this)(i, j), rhs(i, j), tol))
				return false;

	return true;
}

template <typename T>
matrix<T>& matrix<T>::transpose()
{
	// just set the indexer type to the transpose of the current index type
	const size_t m_orig_rows = m_n_rows;
	const size_t m_orig_cols = m_n_cols;
	const_cast<size_t &>(m_n_rows) = m_orig_cols;
	const_cast<size_t &>(m_n_cols) = m_orig_rows;
	m_idx = std::auto_ptr<typename matrix<T>::indexer>(m_idx->make_transpose());

	return *this;
}

template <typename T>
matrix<T> matrix<T>::get_transpose() const
{
	matrix<T> mt = *this;

	return mt.transpose();
}

template <typename T>
matrix<T> matrix<T>::I(size_t n)
{
	matrix<T> ident(n, n);
	for (size_t i = 0 ; i < n ; i++)
		ident(i, i) = 1.0;

	return ident;
}

template <typename T>
matrix<T> matrix<T>::diag(const std::valarray<T>& w)
{
	const size_t n = w.size();	assert(n > 0);
	matrix<T> d(n, n);

	for (size_t i = 0 ; i < n ; i++)
		d(i, i) = w[i];

	return d;
}

template <typename T>
matrix<T>& matrix<T>::svd(std::valarray<T>& w, matrix<T>& V)
{
	scoped_index_change this_index(*this, m_idx->make_one());	// used one-based indexing for this algorithm

	const size_t m = rows();
	const size_t n = cols();

	w = std::valarray<T>((T)0, n);
	V = matrix<T>(n, n);	// not transposed until the very end of the algorithm
	scoped_index_change Vt_index(V, V.m_idx->make_one());	// make sure V is indexed 1..n

	matrix<T> rv1(1, n);
	rv1.set_index_type(RowMajorOne);

	matrix<T>& a = *this;	// little bit less ugly. TODO - add static svd(), pass *this, or copy

	size_t i, its, j, jj, k, l, nm;
	T anorm, c, f, g, h, s, scale, x, y, z;

	for (i = 1 ; i <= n ; i++)
	{
		l = i + 1;
		rv1(1, i) = scale * g;

		// Householder reduction to bidagonal form
		s = 0.0;
		g = 0.0;
		scale = 0.0;

		if (i <= m)
		{
			for (k = i ; k <= m ; k++)
				scale += abs(a(k, i));

			if (scale != 0.0)	// epsilon?
			{
				for (k = i ; k <= m ; k++)
				{
					a(k, i) /= scale;
					s += a(k, i) * a(k, i);
				}
				f = a(i, i);
				g = -sign(sqrt(s), f);
				h = f * g - s;
				a(i, i) = f - g;
				for (j = l ; j <=n ; j++)
				{
					s = 0.0;
					for (k = i ; k <= m ; k++)
						s += a(k, i) * a(k, j);

					f = s / h;

					for (k = i ; k <=m ; k++)
						a(k, j) += f * a(k, i);
				}
				for (k = i ; k <= m ; k++)
					a(k, i) *= scale;
			}
		}
		w[i - 1] = scale * g;

		g = 0.0;
		s = 0.0;
		scale = 0.0;

		if (i <= m && i != n)
		{
			for (k = l ; k <= n ; k++)
				scale += abs(a(i, k));

			if (scale != 0.0)	// epsilon?
			{
				for (k = l ; k <= n ; k++)
				{
					a(i, k) /= scale;
					s += a(i, k) * a(i, k);
				}
				f = a(i, l);
				g = -sign(sqrt(s), f);
				h = f * g - s;
				a(i, l) = f - g;

				for (k = l ; k <=n ; k++)
					rv1(1, k) = a(i, k) / h;

				for (j = l ; j <= m ; j++)
				{
					s = 0.0;
					for (k = l ; k <= n ; k++)
						s += a(j, k) * a(i, k);

					for (k = l ; k <=n ; k++)
						a(j, k) += s * rv1(1, k);
				}
				for (k = l ; k <= n ; k++)
					a(i, k) *= scale;
			}
		}
		anorm = max(anorm, (abs(w[i - 1]) + abs(rv1(1, i))));
	}

	assert(l == n + 1);

	// Accumulation of right-hand transformations
	for (i = n ; i >= 1 ; i--)
	{
		if (i < n)
		{
			if (g != 0)
			{
				for (j = l ; j <= n ; j++)
					V(j, i) = (a(i, j) / a(i, l)) / g;

				for (j= l ; j <= n ; j++)
				{
					s = 0.0;
					for (k = l ; k <= n ; k++)
						s+= a(i, k) * V(k, j);

					for (k = l ; k <= n ; k++)
						V(k, j) += s * V(k, i);
				}
			}
			for (j = l ; j <= n ; j++)
				V(i, j) = V(j, i) = 0.0;
		}
		V(i, i) = 1.0;
		g = rv1(1, i);
		l = i;
	}

	// Accumulation of left-hand transformations
	for (i = min(m, n) ; i >= 1 ; i--)
	{
		l = i + 1;
		g = w[i - 1];

		for (j = l ; j <= n ; j++)
			a(i, j) = 0.0;

		if (g != 0.0)
		{
			g = 1.0 / g;
			for (j = l ; j <= n ; j++)
			{
				s = 0.0;
				for (k = l ; k <= m ; k++)
					s += a(k, i) * a(k, j);

				f = (s / a(i, i)) * g;

				for (k = i ; k <= m ; k++)
					a(k, j) += f * a(k, i);
			}
			for (j = i ; j <= m ; j++)
				a(j, i) *= g;
		}
		else
		{
			for (j = i ; j <= m ; j++)
				a(j, i) = 0.0;
		}
		a(i, i) += 1.0; // ++a[i][i] in NRC...
	}

	// Diagonalization of the bidiagonal form:
	// Loop over singular values, and over allowed iterations
	for (k = n ; k >= 1 ; k--)
	{
		const size_t MAX_SVD_ITS = 30;	// xxx
		for (its = 1 ; its <= MAX_SVD_ITS ; its++)
		{
			bool flag = true;
			for (l = k ; l >= 1 ; l--)
			{
				// Test for splitting.
				// Note that rv1[1] is always 0.
				nm = l - 1;

				// NRC typecasts this sum before comparing for some reason
				if ((abs(rv1(1, l)) + anorm) == anorm)
				{
					flag = false;
					break;
				}

				// Supposedly, we shouldn't ever hit this assert because
				// of the condition that rv[1] is always 0 (we'll break
				// out of the loop above)
				assert(nm >= 1);
				if ((abs(w[nm - 1]) + anorm) == anorm)
					break;
			}

			if (flag)
			{
				c = 0.0;	// Cancellation of rv1[l] if l > 1
				s = 1.0;

				for (i = l ; i <= k ; i++)
				{
					f = s * rv1(1, i);
					rv1(1, i) = c * rv1(1, i);
					if ((abs(f) + anorm) == anorm)
						break;

					g = w[i - 1];
					h = pythag(f, g);
					w[i - 1] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 1 ; j <= m ; j++)
					{
						y = a(j, nm);
						z = a(j, i);
						a(j, nm) = y * c + z * s;
						a(j, i) = z * c - y * s;
					}
				}
			}
			z = w[k - 1];
			if (l == k)	// Convergence
			{
				if (z < 0.0)
				{
					w[k - 1] = -z;	// Singular value is made nonnegative
					for (j = 1 ; j <= n ; j++)
						V(j, k) = -V(j, k);
				}

				break;
			}

			if (its == MAX_SVD_ITS)
			{
				// no convergence
				assert(false);	// xxx - fix this
			}

			nm = k - 1;
			assert(nm >= 1);

			x = w[l - 1];
			y = w[nm - 1];
			g = rv1(1, nm);
			h = rv1(1, k);
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, (T)1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;

			// Next QR transformation
			c = 1.0;
			s = 1.0;
			for (j = l ; j <= nm ; j++)
			{
				i = j + 1;
				g = rv1(1, i);
				y = w[i - 1];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1(1, j) = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;

				for (jj = 1 ; jj <= n ; jj++)
				{
					x = V(jj, j);
					z = V(jj, i);
					V(jj, j) = x * c + z * s;
					V(jj, i) = z * c - x * s;
				}

				z = pythag(f, h);
				w[j - 1] = z;
				if (z != 0.0)
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;

				for (jj = 1 ; jj <= m ; jj++)
				{
					y = a(jj, j);
					z = a(jj, i);
					a(jj, j) = y * c + z * s;
					a(jj, i) = z * c - y * s;
				}
			}
			rv1(1, l) = 0.0;
			rv1(1, k) = f;
			w[k - 1] = x;
		}
	}

	return a;
}

template <typename T>
T& matrix<T>::operator()(size_t i, size_t j)
{
	if (i < r_begin() || i > r_end() || j < c_begin() || j > c_end())
		throw matrix_index_exception<T>(this, i, j);

	return m_A[_idx(i, j)];
}

template <typename T>
const T& matrix<T>::operator()(size_t i, size_t j) const
{
	if (i < r_begin() || i > r_end() || j < c_begin() || j > c_end())
		throw matrix_index_exception<T>(this, i, j);

	return m_A[_idx(i, j)];
}

template <typename _T>
inline
matrix<_T> operator*(const matrix<_T>& a, const matrix<_T>& b)
{
	if (a.cols() != b.rows())
		throw matrix_dimension_mismatch_exception<_T>(a, b);

	matrix<_T> c(a.rows(), b.cols());

	// let's make sure that matrix::scoped_index_change works...
	typename matrix<_T>::indexer* aip = a.m_idx.get();
	typename matrix<_T>::indexer* bip = b.m_idx.get();
	{
		// use the zero-based index for each of our matrices
		typename matrix<_T>::scoped_index_change a_index(a, a.m_idx->make_zero());
		typename matrix<_T>::scoped_index_change b_index(b, b.m_idx->make_zero());

		const size_t m = a.rows();
		const size_t n = b.cols();
		const size_t p = a.cols();	// == b.rows();

		for (size_t i = 0 ; i < m ; i++)
		{
			for (size_t j = 0 ; j < n ; j++)
			{
				for (size_t k = 0 ; k < p ; k++)
				{
					const _T& c_ij = c(i, j);
					const _T& a_ik = a(i, k);
					const _T& b_kj = b(k, j);

					c(i, j) = c_ij + a_ik * b_kj;
				}
			}
		}
	}
	assert(a.m_idx.get() == aip);
	assert(b.m_idx.get() == bip);

	return c;
}

template <typename _T>
inline
matrix<_T> operator*(const matrix<_T>& m, const _T& c)
{
	matrix<_T> mc(m);
	mc.m_A = mc.m_A * c;

	return mc;
}

template <typename _T>
inline
matrix<_T> operator*(const _T& c, const matrix<_T>& m)
{
	matrix<_T> mc(m);
	mc.m_A = c * mc.m_A;

	return mc;
}

template <typename _T>
inline
matrix<_T> operator+(const matrix<_T>& a, const matrix<_T>& b)
{
	if (a.rows() != b.rows() || a.cols() != b.cols())
		throw matrix_dimension_mismatch_exception<_T>(a, b);

	// xxx for this and operator-, what if the matrices
	// are the same size, but one is transposed?
	// Better just add item by item
	matrix<_T> c(a.rows(), a.cols());
	c.m_A = a.m_A + b.m_A;

	return c;
}

template <typename _T>
inline
matrix<_T> operator-(const matrix<_T>& a, const matrix<_T>& b)
{
	if (a.rows() != b.rows() || a.cols() != b.cols())
		throw matrix_dimension_mismatch_exception<_T>(a, b);

	matrix<_T> c(a.rows(), a.cols());
	c.m_A = a.m_A - b.m_A;

	return c;
}

template <typename _T>
inline
matrix<_T> operator-(const matrix<_T>& m)
{
	matrix<_T> mn(m.rows(), m.cols());
	mn.m_A = -m.m_A;

	return mn;
}

template <typename _T>
inline
std::ostream& operator<<(std::ostream& os, const matrix<_T>& m)
{
	for (size_t i = m.r_begin() ; i < m.r_end() ; i++)
	{
		for (size_t j = m.c_begin() ; j < m.c_end() ; j++)
		{
			os << m(i, j) << " ";
		}

		os << std::endl;
	}

	return os;
}

class matrix_exception
{
public:
	virtual std::string what() const = 0;
	virtual ~matrix_exception() { }
};

template <typename T>
class matrix_index_exception : public matrix_exception
{
protected:
	const size_t num_rows;
	const size_t num_cols;
	const size_t req_i;
	const size_t req_j;
public:
	// TODO - busted
	matrix_index_exception(const matrix<T>* m, size_t i, size_t j)
	: num_rows(m->rows())
	, num_cols(m->cols())
	, req_i(i)
	, req_j(j)
	{
		assert(req_i >= num_rows || req_j >= num_cols);
	}

	virtual std::string what() const
	{
		std::stringstream ss;
		if (req_i >= num_rows && req_j < num_cols)
		{
			ss << "Element (" << req_i << ", " << req_j << " )"
				<< "has invalid row (max row is " << num_rows - 1;
		}
		else if (req_i < num_rows && req_j >= num_cols)
		{
			ss << "Element (" << req_i << ", " << req_j << " )"
				<< "has invalid column (max column is " << num_cols - 1;
		}
		else if (req_i >= num_rows && req_j >= num_cols)
		{
			ss << "Element (" << req_i << ", " << req_j << " )"
				<< "invalid (max row " << num_rows - 1
				<< ", max col " << num_cols - 1;
		}
		else
		{
			ss << "INVALID EXCEPTION!"
				<< " (i " << req_i
				<< ", j " << req_j
				<< ", num_rows " << num_rows
				<< ", num_cols " << num_cols
				<< " )";
		}

		return ss.str();
	}
};

template <typename T>
class matrix_dimension_mismatch_exception : public matrix_exception
{
public:
	const size_t m1_rows;
	const size_t m1_cols;
	const size_t m2_rows;
	const size_t m2_cols;

	matrix_dimension_mismatch_exception(const matrix<T>& m1, const matrix<T>& m2)
	: m1_rows(m1.rows())
	, m1_cols(m1.cols())
	, m2_rows(m2.rows())
	, m2_cols(m2.cols())
	{

	}

	virtual std::string what() const { return std::string("dimension mismatch");	 }// bleh
};

class invalid_matrix_exception : public matrix_exception
{
	virtual std::string what() const { return std::string("invalid matrix");	}// bleh
};

};

#endif /* MATRIX_H_ */
