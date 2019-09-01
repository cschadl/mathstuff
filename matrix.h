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
#include <algorithm>
#include <limits>
#include <exception>

#include <assert.h>

#include "misc.h"
#include "vectors.h"

namespace maths
{

class matrix_exception;
template <typename T> class matrix_index_exception;
template <typename T> class matrix_dimension_mismatch_exception;
class invalid_matrix_exception;

template <typename T> class diag_matrix;

/// \remark It was kind of fun to write this, but it's probably best to
///			just use Eigen instead :(
template <typename T>
class matrix
{
protected:
	// Index helper interface
	// Allows 1-based indexing and instant transposition
	class indexer
	{
	protected:
		const matrix<T> * const mp;
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
	const size_t						m_n_rows;
	const size_t						m_n_cols;
	std::valarray<T>					m_A;	// matrix entries
	mutable std::unique_ptr<indexer>	m_idx;	// should be a unique_ptr

	/** Protected constructor for subclasses (e.g. diag_matrix) */
	matrix(size_t m, size_t n, const std::valarray<T>& A, indexer* idx);

	virtual size_t _idx(size_t i, size_t j) const {  return m_idx->index(i, j); }

public:
	explicit matrix(size_t m, size_t n);
	explicit matrix(const T* cols, size_t m, size_t n);
	matrix(const matrix<T>& rhs);
	matrix<T>& operator=(const matrix<T>& rhs);

	virtual ~matrix() { }

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
	bool is_null() const;	// No tolerance, returns true if all entries are exactly 0

	matrix<T>& transpose();
	matrix<T>  get_transpose() const;

	//virtual matrix<T>& inverse();
	//matrix<T>& pseudo_inverse();	// wouldn't be any more efficient than get_pseudo_inverse()...
	matrix<T>  get_pseudo_inverse() const;	// gets the pseudo-inverse of the matrix from the SVD

	virtual matrix<T>& resize(size_t m, size_t n);
	virtual matrix<T>& fill(const T& val);

	virtual matrix<T>& row_swap(size_t i, size_t j);
	virtual matrix<T>& col_swap(size_t i, size_t j);

	/** returns a nxn identity matrix **/
	static matrix<T> I(size_t n);

	/**
	 *  create an nxn diagonal matrix
	 *
	 *  \details We return a pointer because we need to enforce that the object returned is const
	 */
	static const diag_matrix<T>* diag(const std::valarray<T>& d);
	static const diag_matrix<T>* diag(const size_t n, const T& x);

	/**
	 *  Compute the singular value decomposition A = U*W*Vt of this matrix.
	 *
	 *  \details	The entries of a are replaced with U, w is a vector whose
	 *  			entries are the diagonal elements of the n x n matrix W.
	 *  			V is returned rather than Vt.
	 *
	 *  \param a	An m x n matrix, whose entries are replaced with U
	 *  \param w	The diagonal entries of the n x n singular value matrix
	 *  \param V	an n x n orthogonal matrix
	 */
	static bool svd(matrix<T>& a, std::valarray<T>& w, matrix<T>& V);

	/**
	 * Calls static svd() function with *this as A matrix.
	 */
	bool svd(std::valarray<T>& w, matrix<T>& V);

	/**
	 * Determine the rank of A given matrix from its singular values
	 */
	static size_t sv_rank(const std::valarray<T>& w);

	/**
	 * Determine the rank of this matrix using the SVD
	 */
	size_t rank() const;

	/**
	 * 	Solves Ax = b using SVD back-substitution algorithm
	 *
	 * 	\returns	n x 1 solution vector x
	 *  \bug		BROKEN!
	 */
	std::valarray<T> svd_solve(const std::valarray<T>& b) const;

	/**
	 * 	Solves Ax = b using SVD back-substitution algorithm
	 *
	 * 	\details 	Re-uses already calculated SVD.  Assumes "small" entries in w have already been set to 0
	 * 	\returns	the n x 1 solution vector x
	 *  \ bug		BROKEN!
	 */
	static std::valarray<T> svd_solve(const matrix<T>& U, const std::valarray<T>& w, const matrix<T>& V, const std::valarray<T>& b);

	/**
	 * Loads this matrix into the specified array, in row major format
	 */
	void as_array(T* array) const;

	virtual const T& operator()(size_t i, size_t j) const;
	virtual T& operator()(size_t i, size_t j);

	// TODO - operators that take matrix rhs are inefficient
	virtual matrix<T>& operator+=(const matrix<T>& rhs) { *this = *this + rhs; return *this; }
	virtual matrix<T>& operator+=(const T& c) { m_A += c; return *this; }
	virtual matrix<T>& operator-=(const matrix<T>& rhs) { *this = *this - rhs; return *this; }
	virtual matrix<T>& operator-=(const T& c) { m_A -= c; return *this; }
	virtual matrix<T>& operator*=(const matrix<T>& rhs) { *this = *this * rhs; return *this; }
	virtual matrix<T>& operator*=(const T& c) { m_A *= c; return *this; }
	virtual matrix<T>& operator/=(const T& d) { m_A /= d; return *this; }

	template <typename _T>
	friend matrix<_T> operator*(const matrix<_T>& a, const matrix<_T>& b);

	template <typename _T>
	friend matrix<_T> operator*(const matrix<_T>& m, const _T& c);

	template <typename _T>
	friend matrix<_T> operator*(const _T& c, const matrix<_T>& m);

	template <typename _T>
	friend std::valarray<_T> operator*(const matrix<_T>& m, const std::valarray<_T>& v);

	template <typename _T>
	friend std::valarray<_T> operator*(const std::valarray<_T>& v, const matrix<_T>& m);

	template <typename _T>
	friend matrix<_T> operator+(const matrix<_T>& a, const matrix<_T>& b);

	template <typename _T>
	friend matrix<_T> operator-(const matrix<_T>& a, const matrix<_T>& b);

	template <typename _T>
	friend matrix<_T> operator-(const matrix<_T>& m);	// negation

	// misc

	/** Too lazy to write a general matrix inversion (e.g. with LU decomposition)
	 *  but a 4x4 matrix inversion is often useful.
	 */
	static bool invert4x4(const matrix<T>& m, matrix<T>& m_inv);

	/**
	 *  Construct the 4x4 rotation matrix for rotating a 3-vector about the given axis.
	 */
	static matrix<T> rotation(const maths::vector3f& axis, T angle_deg);
	static matrix<T> translation(const maths::vector3f& v);

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

protected:

	class scoped_index_change
	{
	protected:
		const matrix<T> * const 	mp;
		std::unique_ptr<indexer>	orig_index;

	public:
		scoped_index_change(const matrix<T>* m, indexer* new_index)
		: mp(m)
		, orig_index(std::move(m->m_idx))
		{
			assert(new_index->get_mp() == mp);
			mp->m_idx = std::unique_ptr<indexer>(new_index);
		}

		~scoped_index_change()
		{
			mp->m_idx = std::move(orig_index);
		}
	};
};

/**
 * An n x n diagonal matrix
 */
template <typename T>
class diag_matrix : public matrix<T>
{
protected:
	const T m_off_diag_element;	// probably 0
	T m_quit_bitching;

	class diag_matrix_indexer : public matrix<T>::indexer
	{
	public:
		diag_matrix_indexer(const diag_matrix<T>* m) : matrix<T>::indexer(m) { }
		virtual size_t index(size_t i, size_t j) const { return (i == j) ? i : -1; }
		virtual size_t row_index_start() const { return 0; }
		virtual size_t row_index_end() const { return this->mp->rows(); }
		virtual size_t col_index_start() const { return 0; }
		virtual size_t col_index_end() const { return this->mp->cols(); }
		virtual typename matrix<T>::indexer* clone(const matrix<T>* m) const;
		virtual typename matrix<T>::indexer* make_transpose() const { return clone(this->mp); }
		virtual typename matrix<T>::indexer* make_zero() const { return clone(this->mp); }
		virtual typename matrix<T>::indexer* make_one() const;
	};

	class diag_matrix_1_indexer : public matrix<T>::indexer
	{
	public:
		diag_matrix_1_indexer(const diag_matrix<T>* m) : matrix<T>::indexer(m) { }
		virtual size_t index(size_t i, size_t j) const { assert(i > 0 && j > 0); return (i == j) ? i - 1 : -1; }
		virtual size_t row_index_start() const { return 1; }
		virtual size_t row_index_end() const { return this->mp->rows() + 1; }
		virtual size_t col_index_start() const { return 1; }
		virtual size_t col_index_end() const { return this->mp->cols() + 1; }
		virtual typename matrix<T>::indexer* clone(const matrix<T>* m) const;
		virtual typename matrix<T>::indexer* make_transpose() const { return clone(this->mp); }
		virtual typename matrix<T>::indexer* make_zero() const;
		virtual typename matrix<T>::indexer* make_one() const { return clone(this->mp); }
	};

	diag_matrix(const std::valarray<T>& d);	// new diagonal matrix with entries of D as elements

private:
	// This class should not be used in a non-const context
	virtual T& operator()(size_t m, size_t n) { assert(false); return m_quit_bitching; }
	diag_matrix<T>& operator=(const diag_matrix<T>& rhs) { assert(false); return *this; }
	diag_matrix(const diag_matrix<T>& rhs) { }
	virtual matrix<T>& operator+=(const matrix<T>& rhs) { assert(false); return *this; }
	virtual matrix<T>& operator+=(const T& c) { assert(false); return *this; }
	virtual matrix<T>& operator-=(const matrix<T>& rhs) { assert(false); return *this; }
	virtual matrix<T>& operator-=(const T& c) { assert(false); return *this; }
	virtual matrix<T>& operator*=(const matrix<T>& rhs) { assert(false); return *this; }
	virtual matrix<T>& operator*=(const T& c) { assert(false); return *this; }
	virtual matrix<T>& operator/=(const T& d) { assert(false); return *this; }
	virtual matrix<T>& resize(size_t m, size_t n) { assert(false); return *this; }
	virtual matrix<T>& fill(const T& val) { assert(false); return *this; }

public:
	virtual const T& operator()(size_t i, size_t j) const;
	virtual diag_matrix<T>& inverse();

	friend class matrix<T>;
};

////////////////////////////////////
/// exception types

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
	const size_t row_start;
	const size_t row_end;
	const size_t col_start;
	const size_t col_end;
public:
	matrix_index_exception(const matrix<T>* m, size_t i, size_t j)
	: num_rows(m->rows())
	, num_cols(m->cols())
	, req_i(i)
	, req_j(j)
	, row_start(m->r_begin())
	, row_end(m->r_end())
	, col_start(m->c_begin())
	, col_end(m->c_end())
	{

	}

	virtual std::string what() const
	{
		std::stringstream ss;
		ss << "Bad index: (" << req_i << ", " << req_j << ")";
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
#include "matrix.hpp"

#endif /* MATRIX_H_ */
