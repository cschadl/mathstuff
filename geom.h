#ifndef GEOM_H
#define GEOM_H

#include "vectors.h"
#include "misc.h"

namespace maths
{

/// An axis aligned bounding box in 3 dimensions
template <typename T>
class bbox_3
{
private:
	n_vector<T, 3>	m_min;
	n_vector<T, 3>	m_max;

public:
	/** Construct an empty bounding box.
	 *  Note that an empty bounding box is different than
	 *  a bounding box consisting of a single point.
	 *  Use is_empty() and is_point() to distinguish
	 *  these cases.
	 */
	bbox_3();

	/** Construct a bounding box with the given points as min, max */
	bbox_3(const n_vector<T, 3>& min, const n_vector<T, 3>& max) : m_min(min), m_max(max) { }

	// Queries
	const n_vector<T, 3>& min() const { return m_min; }
	const n_vector<T, 3>& max() const { return m_max; }
	n_vector<T, 3>  center() const;

	bool operator==(const bbox_3<T>& bbox) const;
	bool operator!=(const bbox_3<T>& bbox) const { return !((*this) == bbox); }

	bool is_empty() const;
	bool is_point() const;

	T extent_x() const;
	T extent_y() const;
	T extent_z() const;
	T max_extent() const;

	bool contains(const n_vector<T, 3>& p) const;

	/** Corners are output CCW top / bottom */
	template <typename OutputIter> void get_corners(OutputIter it) const;

	/** Add a single point
	 *  Returns true if adding the point increased
	 *  the size of the bounding box.
	 */
	bool add_point(const n_vector<T, 3>& pt);

	/** Add a bunch of points.
	 *  After adding the points, the bbox will be the minimum
	 *  axis-aligned bbox that encloses the given points.
	 *  Returns true if any of the points added increased the
	 *  size of the bounding box.
	 */
	template <typename Iter> bool add_points(Iter begin, Iter end);
};

template <typename T>
bbox_3<T>::bbox_3()
: m_min(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max())
, m_max(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max())
{
	// box is empty if min / max is std::numeric_limits<T>::max()
}

template <typename T>
bool bbox_3<T>::is_empty() const
{
	return 	m_min == n_vector<T, 3>(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max()) &&
			m_max == n_vector<T, 3>(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max());
}

template <typename T>
bool bbox_3<T>::is_point() const
{
	return !is_empty() && (m_min == m_max);
}

template <typename T>
n_vector<T, 3> bbox_3<T>::center() const
{
	return (m_min + m_max) / 2.0;
}

template <typename T>
T bbox_3<T>::extent_x() const
{
	return m_max.x() - m_min.x();
}

template <typename T>
T bbox_3<T>::extent_y() const
{
	return m_max.y() - m_min.y();
}

template <typename T>
T bbox_3<T>::extent_z() const
{
	return m_max.z() - m_min.z();
}

template <typename T>
T bbox_3<T>::max_extent() const
{
	return std::max(std::max(extent_x(), extent_y()), extent_z());
}

template <typename T>
bool bbox_3<T>::contains(const n_vector<T, 3>& p) const
{
	return 	p.x() >= m_min.x() && p.x() <= m_max.x() &&
			p.y() >= m_min.y() && p.y() <= m_min.y() &&
			p.z() >= m_min.z() && p.z() <= m_min.z();
}

template <typename T> template <typename OutputIter>
void bbox_3<T>::get_corners(OutputIter it) const
{
	n_vector<T, 3> v = m_min;
	*it++ = v;			// lower bottom left
	v.x() = m_max.x();
	*it++ = v;			// lower bottom right
	v.y() = m_max.y();
	*it++ = v;			// lower top right
	v.x() = m_min.x();
	*it++ = v;			// lower top left
	v.y() = m_min.y();

	v.z() = m_max.z();
	*it++ = v;			// upper bottom left
	v.x() = m_max.x();
	*it++ = v;			// upper bottom right
	v.y() = m_max.y();
	*it++ = v;			// upper top left
	v.x() = m_min.x();
	*it++ = v;			// upper top right
}

template <typename T>
bool bbox_3<T>::add_point(const n_vector<T, 3>& p)
{
	if (is_empty())
	{
		m_min = p;
		m_max = p;

		return true;
	}

	if (contains(p))
		return false;

	bool changed = false;
	if (p.x() < m_min.x())
	{	m_min.x() = p.x();	changed = true; }
	else if (p.x() > m_max.x())
	{	m_max.x() = p.x();	changed = true; }

	if (p.y() < m_min.y())
	{	m_min.y() = p.y();	changed = true; }
	else if (p.y() > m_max.y())
	{	m_max.y() = p.y();	changed = true; }

	if (p.z() < m_min.z())
	{	m_min.z() = p.z();	changed = true; }
	else if (p.z() > m_max.z())
	{	m_max.z() = p.z();	changed = true; }

	return changed;
}

template <typename T> template <typename Iter>
bool bbox_3<T>::add_points(Iter begin, Iter end)
{
	bool added = false;
	for (Iter i = begin ; i != end ; ++i)
	{
		bool b = add_point(*i);
		added = added || b;
	}

	return added;
}

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
	bbox_3<T> bbox() const;

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

template <typename T>
bbox_3<T> triangle_3<T>::bbox() const
{
	bbox_3<T> bbox;
	bbox.add_points(m_verts, m_verts + 3);

	return bbox;
}

};

#endif
