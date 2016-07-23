#ifndef _INCLUDE_PRIMITIVE_FITTING_H
#define _INCLUDE_PRIMITIVE_FITTING_H

#include <type_traits>
#include <Eigen/SVD>
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

	/** Best fit line through a 3d point cloud. */
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
};

};

#endif
