#ifndef GEOM_H
#define GEOM_H

#include "vectors.h"
#include "misc.h"

namespace maths
{

template <typename T>
class triangle_3
{
private:
	const n_vector<T, 3>	m_verts[3];
	const n_vector<T, 3>	m_normal;

public:
	/** Computes the normal from ABC (right-hand rule) */
	triangle_3(const n_vector<T, 3>& a, const n_vector<T, 3>& b, const n_vector<T, 3>& c);
	/** Explicitly specify the normal */
	triangle_3(const n_vector<T, 3>& a, const n_vector<T, 3>& b, const n_vector<T, 3>& c, const n_vector<T, 3>& n);

	/** Two triangles T1 and T2 are equal iff there is a cyclic permutation from the
	 *  vertices of T1 onto the vertices of T2.
	 */
	bool operator==(const triangle_3<T>& t2) const;
	bool operator!=(const triangle_3<T>& t2) const;

	/** A triangle is degenerate if it has two or more equal vertices */
	bool is_degenerate() const;

	/** Gets the i-th vertex of the triangle modulo 3 */
	const n_vector<T, 3>& vert(int i) const { return m_verts[i % 3]; }
	const n_vector<T, 3>& operator[](int i) const { return m_verts[i % 3]; }
	const n_vector<T, 3>& normal() const { return m_normal; }

	/** The same triangle, with the normal reversed */
	triangle_3<T> flipped() const;

	// TODO -
	//plane<T> supporting_plane() const;
	//bbox_3<T> bbox() const;

	//bool is_point_on(const n_vector<T, 3>& point) const;
	//bool line_intersection(const line<T>& line, n_vector<T, 3>& intersect_pt) const;
	//bool barycentric_point(const n_vector<T, 3>& point, n_vector<T, 2>& barycentric_pt) const;
};

typedef triangle_3<double> triangle3d;
typedef triangle_3<float>  triangle3f;

template <typename T>
triangle_3<T>::triangle_3(const n_vector<T, 3>& a, const n_vector<T, 3>& b, const n_vector<T, 3>& c)
{
	const_cast< n_vector<T, 3>& >(m_verts[0]) = a;
	const_cast< n_vector<T, 3>& >(m_verts[1]) = b;
	const_cast< n_vector<T, 3>& >(m_verts[2]) = c;

	// compute the normal
	const_cast< n_vector<T, 3>& >(m_normal) = ((b - a) % (c - a)).unit();
}

template<typename T>
triangle_3<T>::triangle_3(const n_vector<T, 3>& a, const n_vector<T, 3>& b, const n_vector<T, 3>& c, const n_vector<T, 3>& n)
{
	const_cast< n_vector<T, 3>& >(m_verts[0]) = a;
	const_cast< n_vector<T, 3>& >(m_verts[1]) = b;
	const_cast< n_vector<T, 3>& >(m_verts[2]) = c;
	const_cast< n_vector<T, 3>& >(m_normal) = n;
}

template <typename T>
bool triangle_3<T>::operator==(const triangle_3<T>& t2) const
{
	// For now, we'll just compare the vertices in order
	// This isn't sufficient, but it should do for now.
	const triangle_3<T>& t1 = *this;
	return (t1[0] == t2[0] && t1[1] == t2[1] && t1[2] == t2[2]);
}

template <typename T>
bool triangle_3<T>::operator!=(const triangle_3<T>& t2) const
{
	return !((*this) == t2);
}

template <typename T>
bool triangle_3<T>::is_degenerate() const
{
	const T eps = std::numeric_limits<T>::epsilon();
	return (m_verts[0].is_close(m_verts[1], eps * eps) ||
			m_verts[0].is_close(m_verts[2], eps * eps) ||
			m_verts[1].is_close(m_verts[2], eps * eps));
}

template <typename T>
triangle_3<T> triangle_3<T>::flipped() const
{
	return triangle_3<T>((*this)[0], (*this)[1], (*this)[2], -normal());
}

};

#endif
