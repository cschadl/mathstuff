// see matrix.h

namespace maths
{

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
matrix<T>::matrix(size_t m, size_t n, const std::valarray<T>& A, matrix<T>::indexer* idx)
: m_n_rows(m)
, m_n_cols(n)
, m_A(A)
, m_idx(idx)
{

}

/**
 * Constructor that takes a pointer to an array as elements
 * Entries are considered column-major (e.g. for OpenGL)
 * The number of elements in cols[] must be m * n
 */
template <typename T>
matrix<T>::matrix(const T* cols, size_t m, size_t n)
: m_n_rows(m)
, m_n_cols(n)
, m_A(cols, m * n)
, m_idx(new col_major_indexer(this))
{
	if (m_n_rows < 1 || m_n_cols < 1)
		throw invalid_matrix_exception();

	assert(m_A.size() > 0);
}

template <typename T>
matrix<T>& matrix<T>::operator=(const matrix<T>& rhs)
{
	// We used to check that rhs had the same dimension as *this
	// and throw an exception if that wasn't the case.  But, this
	// broke A *= B, and also seemed awkward in some situations,
	// so we no longer do that.
	// I'm still a bit iffy on this.
	const_cast<size_t&>(m_n_rows) = rhs.m_n_rows;
	const_cast<size_t&>(m_n_cols) = rhs.m_n_cols;
	m_A = rhs.m_A;
	m_idx.reset(rhs.m_idx->clone(this));

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
	scoped_index_change ai(this, m_idx->make_zero());
	scoped_index_change bi(&rhs, rhs.m_idx->make_zero());

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
matrix<T>& matrix<T>::resize(size_t m, size_t n)
{
	// xxx - throw if m == 0 || n == 0
	const_cast<size_t &>(m_n_rows) = m;
	const_cast<size_t &>(m_n_cols) = n;

	// resizing the matrix re-initializes to row-major
	// check for 1-based index (kinda icky)
	std::auto_ptr<indexer> new_indexer;
	if (c_begin() > 0)
	{
		assert(c_begin() == 1);
		assert(r_begin() == 1);
		new_indexer.reset(new row_major_1_indexer(this));
	}
	else
	{
		assert(c_begin() == 0);
		assert(r_begin() == 0);
		new_indexer.reset(new row_major_indexer(this));
	}
	m_idx = new_indexer;
	m_A = std::valarray<T>((T)0, m * n);

	return *this;
}

template <typename T>
matrix<T>& matrix<T>::fill(const T& val)
{
	m_A = val;
	return *this;
}

template <typename T>
matrix<T>& matrix<T>::row_swap(size_t i, size_t j)
{
	// xxx - need to add some sort of "row / col" iterator
	for (size_t ci = c_begin() ; ci < r_end() ; ci++)
	{
		std::swap((*this)(i, ci), (*this)(j, ci));
	}

	return *this;
}

template <typename T>
matrix<T>& matrix<T>::col_swap(size_t i, size_t j)
{
	for (size_t ri = r_begin() ; ri < r_end() ; ri++)
	{
		std::swap((*this)(ri, i), (*this)(ri, j));
	}

	return *this;
}

//static
template <typename T>
matrix<T> matrix<T>::I(size_t n)
{
	matrix<T> ident(n, n);
	for (size_t i = 0 ; i < n ; i++)
		ident(i, i) = 1.0;

	return ident;
}

//static
template <typename T>
const diag_matrix<T>* matrix<T>::diag(const std::valarray<T>& w)
{
	return new diag_matrix<T>(w);
}

// static
template <typename T>
bool matrix<T>::svd(matrix<T>& a, std::valarray<T>& w, matrix<T>& V)
{
	scoped_index_change this_index(&a, a.m_idx->make_one());	// used one-based indexing for this algorithm

	const size_t m = a.rows();
	const size_t n = a.cols();

	// re-initialize w and V, make sure they're the correct size.
	// (maybe just fill() here if their size is compatible)
	w = std::valarray<T>((T)0, n);
	V = matrix<T>(n, n);

	scoped_index_change Vt_index(&V, V.m_idx->make_one());	// make sure V is indexed 1..n

	matrix<T> rv1(1, n);
	rv1.m_idx = std::auto_ptr<indexer>(new row_major_1_indexer(&rv1));

	size_t i, its, j, jj, k, l, nm;
	T anorm, c, f, g, h, s, scale, x, y, z;

	scale = 0.0;
	s = 0.0;
	g = 0.0;

	for (i = 1 ; i <= n ; i++)
	{
		l = i + 1;
		rv1(1, i) = scale * g;

		g = 0.0;
		s = 0.0;
		scale = 0.0;

		// Householder reduction to bidagonal form
		if (i <= m)
		{
			for (k = i ; k <= m ; k++)
				scale += abs(a(k, i));

			if (scale != 0.0)
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

			if (scale != 0.0)
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
					// I think these tests for 0.0 + anorm can
					// be improved by taking a small epsilon value
					// into account (improved convergence)
					flag = false;
					break;
				}

				// Supposedly, we shouldn't ever hit this assert because
				// of the condition that rv[1] is always 0 (we'll break
				// out of the loop above)
				if (nm < 1)
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

			if (its == MAX_SVD_ITS)	// didn't converge :(
				return false;

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

	// sort the singular values in descending order
	static bool sort_sv = true;
	if (sort_sv)
	{
		for (i = 0 ; i < n - 1 ; i++)
		{
			T sv_largest = w[i];
			size_t sv_largest_idx = i;
			for (j = i + 1 ; j < n ; j++)
			{
				if (w[j] > sv_largest)
				{
					sv_largest = w[j];
					sv_largest_idx = j;
				}
			}

			if (sv_largest_idx != i)
			{
				std::swap(w[i], w[sv_largest_idx]);
				a.col_swap(i + 1, sv_largest_idx + 1);
				V.col_swap(i + 1, sv_largest_idx + 1);
			}
		}
	}

	return true;
}

template <typename T>
bool matrix<T>::svd(std::valarray<T>& w, matrix<T>& V)
{
	return svd(*this, w, V);
}

template <typename T>
T& matrix<T>::operator()(size_t i, size_t j)
{
	if (i < r_begin() || i > r_end() || j < c_begin() || j > c_end())
		throw matrix_index_exception<T>(this, i, j);

	const size_t n = _idx(i, j);
	if (n >= m_A.size())
		throw matrix_index_exception<T>(this, i, j);

	return m_A[n];
}

template <typename T>
const T& matrix<T>::operator()(size_t i, size_t j) const
{
	if (i < r_begin() || i > r_end() || j < c_begin() || j > c_end())
		throw matrix_index_exception<T>(this, i, j);

	const size_t n = _idx(i, j);
	if (n >= m_A.size())
		throw matrix_index_exception<T>(this, i, j);

	return m_A[n];
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
		typename matrix<_T>::scoped_index_change a_index(&a, a.m_idx->make_zero());
		typename matrix<_T>::scoped_index_change b_index(&b, b.m_idx->make_zero());

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

	matrix<_T> c(a.rows(), a.cols());

	typename matrix<_T>::scoped_index_change ai(&a, a.m_idx->make_zero());
	typename matrix<_T>::scoped_index_change bi(&b, b.m_idx->make_zero());
	for (size_t i = 0 ; i < a.rows() ; i++)
		for (size_t j = 0 ; j < b.cols() ; j++)
			c(i, j) = a(i, j) + b(i, j);

	return c;
}

template <typename _T>
inline
matrix<_T> operator-(const matrix<_T>& a, const matrix<_T>& b)
{
	if (a.rows() != b.rows() || a.cols() != b.cols())
		throw matrix_dimension_mismatch_exception<_T>(a, b);

	matrix<_T> c(a.rows(), a.cols());

	typename matrix<_T>::scoped_index_change ai(&a, a.m_idx->make_zero());
	typename matrix<_T>::scoped_index_change bi(&b, b.m_idx->make_zero());
	for (size_t i = 0 ; i < a.rows() ; i++)
		for (size_t j = 0 ; j < b.cols() ; j++)
			c(i, j) = a(i, j) - b(i, j);

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

//////////////////////////////////////////////
// diag_matrix
template <typename T>
diag_matrix<T>::diag_matrix(const std::valarray<T>& d)
: matrix<T>(d.size(), d.size(), d, new typename diag_matrix<T>::diag_matrix_indexer(this))
, m_off_diag_element((T)0)
, m_quit_bitching((T)0)
{

}

template <typename T>
const T& diag_matrix<T>::operator()(size_t i, size_t j) const
{
	const size_t _idx = this->m_idx->index(i, j);
	if (_idx == (size_t)-1)
		return m_off_diag_element;

	return this->m_A[_idx];
}

template <typename T>
typename matrix<T>::indexer* diag_matrix<T>::diag_matrix_indexer::clone(const matrix<T>* m) const
{
	const diag_matrix<T>* dm = dynamic_cast< const diag_matrix<T>* >(m);
	assert(dm);

	return new diag_matrix_indexer(dm);
}

template <typename T>
typename matrix<T>::indexer* diag_matrix<T>::diag_matrix_indexer::make_one() const
{
	const diag_matrix<T>* dm = dynamic_cast< const diag_matrix<T>* >(this->mp);
	assert(dm);

	return new diag_matrix_1_indexer(dm);
}

template <typename T>
typename matrix<T>::indexer* diag_matrix<T>::diag_matrix_1_indexer::clone(const matrix<T>* m) const
{
	const diag_matrix<T>* dm = dynamic_cast< const diag_matrix<T>* >(m);
	assert(dm);

	return new diag_matrix_1_indexer(dm);
}

template <typename T>
typename matrix<T>::indexer* diag_matrix<T>::diag_matrix_1_indexer::make_zero() const
{
	const diag_matrix<T>* dm = dynamic_cast< const diag_matrix<T>* >(this->mp);
	assert(dm);

	return new diag_matrix_indexer(dm);
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
