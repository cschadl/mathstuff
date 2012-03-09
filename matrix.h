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
	matrix<T>& operator=(const matrix<T>& rhs);

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
	matrix<T>& resize(size_t m, size_t n);
	matrix<T>& fill(const T& val);

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

};
#include "matrix.hpp"

#endif /* MATRIX_H_ */
