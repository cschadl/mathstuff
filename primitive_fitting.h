#ifndef _INCLUDE_PRIMITIVE_FITTING_H
#define _INCLUDE_PRIMITIVE_FITTING_H

#include <type_traits>
#include <numeric>

#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/Dense>

#include "geom.h"

namespace maths
{

namespace primitive_fitting
{
	template <typename InputIterator, typename PointType>
	Eigen::JacobiSVD< Eigen::Matrix<typename PointType::value_type, 3, Eigen::Dynamic> >
	point3svd(InputIterator begin, InputIterator end, const PointType & origin)
	{
		static_assert(std::is_same<typename std::decay<decltype(*begin)>::type, PointType>::value, "Type mismatch");

		typedef typename PointType::value_type value_type;
		typedef Eigen::Matrix<value_type, 3, Eigen::Dynamic> MatrixType;

		const size_t n = std::distance(begin, end);

		MatrixType A(3, n);
		for (size_t i = 0 ; i < n ; i++)
		{
			PointType p_i = (*(begin + i)) - origin;

			A(0, i) = p_i.x();
			A(1, i) = p_i.y();
			A(2, i) = p_i.z();
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
		const size_t n = std::distance(begin, end);
		if (n < 3)
			return false;

		PointType origin = centroid(begin, end);
		auto svd = point3svd(begin, end, origin);

		auto u2 = svd.matrixU().col(2);

		out_point = origin;

		out_normal.x() = u2(0);
		out_normal.y() = u2(1);
		out_normal.z() = u2(2);

		return true;
	}

	/** Best fit line through a 3d point cloud.
	 *  TODO - template this to work with 2d as well.
	 */
	template <typename InputIterator, typename PointType>
	bool line(	InputIterator begin, InputIterator end,
				PointType & out_point,
				PointType & out_dir	)
	{
		const size_t n = std::distance(begin, end);
		if (n < 2)
			return false;

		PointType origin = centroid(begin, end);
		auto svd = point3svd(begin, end, origin);

		auto u0 = svd.matrixU().col(0);

		out_point = origin;

		out_dir.x() = u0(0);
		out_dir.y() = u0(1);
		out_dir.z() = u0(2);

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
		typedef typename PointType::value_type value_type;
		typedef Eigen::Matrix<value_type, 4, 4>	Matrix_t;
		typedef Eigen::Matrix<value_type, 4, 1> Vector_t;

		auto const n = (value_type) std::distance(begin, end);

		value_type sum_x = value_type(0);
		value_type sum_y = value_type(0);
		value_type sum_z = value_type(0);
		value_type sum_xy = value_type(0);
		value_type sum_xz = value_type(0);
		value_type sum_yz = value_type(0);
		value_type sum_rho_sq = value_type(0);
		value_type sum_rho_sq_x = value_type(0);
		value_type sum_rho_sq_y	= value_type(0);
		value_type sum_rho_sq_z = value_type(0);
		value_type sum_x_sq = value_type(0);
		value_type sum_y_sq = value_type(0);
		value_type sum_z_sq = value_type(0);

		for (auto it = begin ; it != end ; ++it)
		{
			const PointType & p = *it;

			sum_x += p.x();
			sum_y += p.y();
			sum_z += p.z();

			sum_xy += p.x() * p.y();
			sum_xz += p.x() * p.z();
			sum_yz += p.y() * p.z();

			value_type rho_sq = (p.x() * p.x()) + (p.y() * p.y()) + (p.z() * p.z());
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

		value_type const a = 0.5 * A[1];
		value_type const b = 0.5 * A[2];
		value_type const c = 0.5 * A[3];

		value_type const abc_sq = (a * a) + (b * b) + (c *  c);

		value_type const r = ::sqrt(abc_sq + A[0]);

		out_sphere_center.x() = a;
		out_sphere_center.y() = b;
		out_sphere_center.z() = c;

		out_sphere_radius = r;

		return true;
	}

};

};

#endif
































