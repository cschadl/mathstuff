#ifndef _INCLUDE_PRIMITIVE_FITTING_H
#define _INCLUDE_PRIMITIVE_FITTING_H

#include <type_traits>
#include <numeric>
#include <iterator>
#include <utility>

#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/Dense>

#include "geom.h"

namespace maths
{

namespace primitive_fitting
{
	template <typename PT>
	using pt_traits = traits::point_traits<PT>;

	template <typename InputIterator, typename PointType, size_t Dim = pt_traits<PointType>::dimension()>
	Eigen::JacobiSVD< Eigen::Matrix<typename pt_traits<PointType>::value_type, Dim, Eigen::Dynamic> >
	pointssvd(InputIterator begin, InputIterator end, const PointType & origin)
	{
		static_assert(std::is_same<typename std::decay<decltype(*begin)>::type, PointType>::value, "Type mismatch");

		typedef typename PointType::value_type scalar;
		typedef Eigen::Matrix<scalar, Dim, Eigen::Dynamic> MatrixType;

		const size_t n = std::distance(begin, end);

		MatrixType A(Dim, n);

		auto it = begin;
		for (size_t i = 0 ; i < n ; i++)
		{
			PointType p_i = *it - origin;

			for (size_t j = 0 ; j < Dim ; j++)
			{
				A(j, i) = p_i[j];
			}

			++it;
		}

		Eigen::JacobiSVD<MatrixType> svd(A, Eigen::ComputeThinU);

		return svd;
	}

	/** Best fit plane through a cloud of 3d points. */
	template <typename InputIterator, typename PointType>
	bool plane(	InputIterator begin,
				InputIterator end,
				PointType & out_point,
				PointType & out_normal)
	{
		static_assert(pt_traits<PointType>::dimension() > 2, "Invalid point dimension");

		const size_t n = std::distance(begin, end);
		if (n < 3)
			return false;

		PointType origin = centroid(begin, end);
		auto svd = pointssvd<InputIterator, PointType, 3>(begin, end, origin);

		if (svd.rank() < 2)
			return false;

		auto u2 = svd.matrixU().col(2);

		out_point = origin;

		out_normal.x() = u2(0);
		out_normal.y() = u2(1);
		out_normal.z() = u2(2);

		return true;
	}

	/** Best fit line through a point cloud.
	 */
	template <typename InputIterator, typename LineType>
	bool line(InputIterator begin, InputIterator end, LineType & out_line)
	{
		typedef typename traits::halfspace_traits<LineType>::point_type point_type;
		typedef typename traits::halfspace_traits<LineType>::vector_type vector_type;

		constexpr size_t Dim = traits::halfspace_traits<LineType>::dimension();

		using create_line = adapters::create_line<LineType, Dim>;

		const size_t n = std::distance(begin, end);
		if (n < 2)
			return false;

		if (n == 2)
		{
			out_line = create_line::from_points(*begin, *std::next(begin));

			return true;
		}

		point_type origin = centroid(begin, end);
		auto svd = pointssvd(begin, end, origin);

		auto u0 = svd.matrixU().col(0);

		vector_type dir;

		for (size_t j = 0 ; j < Dim ; j++)
			dir[j] = u0(j);

		out_line = create_line::from_point_dir(std::move(origin), std::move(dir));

		return true;
	}

	/** Best fit sphere around a 3d point cloud
	 *  Taken from http://math.stackexchange.com/questions/820815/about-sphere-equation-z-abxcydx2-ey2
	 *  which was taken from http://fr.scribd.com/doc/14819165/Regressions-coniques-quadriques-circulaire-spherique
	 */
	template <typename InputIterator, typename PointType>
	bool sphere(InputIterator begin, InputIterator end,
				PointType & out_sphere_center,
				typename PointType::value_type & out_sphere_radius)
	{
		typedef typename PointType::value_type scalar;
		typedef Eigen::Matrix<scalar, 4, 4>	Matrix_t;
		typedef Eigen::Matrix<scalar, 4, 1> Vector_t;

		auto const n = (scalar) std::distance(begin, end);

		scalar sum_x = scalar(0);
		scalar sum_y = scalar(0);
		scalar sum_z = scalar(0);
		scalar sum_xy = scalar(0);
		scalar sum_xz = scalar(0);
		scalar sum_yz = scalar(0);
		scalar sum_rho_sq = scalar(0);
		scalar sum_rho_sq_x = scalar(0);
		scalar sum_rho_sq_y	= scalar(0);
		scalar sum_rho_sq_z = scalar(0);
		scalar sum_x_sq = scalar(0);
		scalar sum_y_sq = scalar(0);
		scalar sum_z_sq = scalar(0);

		for (auto it = begin ; it != end ; ++it)
		{
			const PointType & p = *it;

			sum_x += p.x();
			sum_y += p.y();
			sum_z += p.z();

			sum_xy += p.x() * p.y();
			sum_xz += p.x() * p.z();
			sum_yz += p.y() * p.z();

			scalar rho_sq = (p.x() * p.x()) + (p.y() * p.y()) + (p.z() * p.z());
			sum_rho_sq += rho_sq;
			sum_rho_sq_x += rho_sq * p.x();
			sum_rho_sq_y += rho_sq * p.y();
			sum_rho_sq_z += rho_sq * p.z();

			sum_x_sq += p.x() * p.x();
			sum_y_sq += p.y() * p.y();
			sum_z_sq += p.z() * p.z();
		}

		Matrix_t M;
		M(0, 0) = n;		M(0, 1) = sum_x;	M(0, 2) = sum_y;	M(0, 3) = sum_z;
		M(1, 0) = sum_x;	M(1, 1) = sum_x_sq;	M(1, 2) = sum_xy;	M(1, 3) = sum_xz;
		M(2, 0) = sum_y;	M(2, 1) = sum_xy;	M(2, 2) = sum_y_sq;	M(2, 3) = sum_yz;
		M(3, 0) = sum_z;	M(3, 1) = sum_xz;	M(3, 2) = sum_yz;	M(3, 3) = sum_z_sq;

		Vector_t q;
		q[0] = sum_rho_sq;
		q[1] = sum_rho_sq_x;
		q[2] = sum_rho_sq_y;
		q[3] = sum_rho_sq_z;

		Eigen::HouseholderQR<Matrix_t> M_qr(M);
		Vector_t A = M_qr.solve(q);

		scalar const a = 0.5 * A[1];
		scalar const b = 0.5 * A[2];
		scalar const c = 0.5 * A[3];

		scalar const abc_sq = (a * a) + (b * b) + (c *  c);

		scalar const r = ::sqrt(abc_sq + A[0]);	// sources say this is sqrt(abc_sq - A[0]), but that's wrong...

		out_sphere_center.x() = a;
		out_sphere_center.y() = b;
		out_sphere_center.z() = c;

		out_sphere_radius = r;

		return true;
	}

};

};

#endif
































