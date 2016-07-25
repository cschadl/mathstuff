/*
 * testgeom.cpp
 *
 *  Created on: Mar 3, 2013
 *      Author: cds
 */

#include "geom.h"
#include "primitive_fitting.h"
#include "tut.h"

#include <stdio.h>
#include <unistd.h>
#include <libgen.h>
#include <sys/param.h>

#include <stdexcept>
#include <fstream>

using namespace maths;

namespace tut
{
	struct test_geom_data
	{
		const vector3d x;
		const vector3d y;
		const vector3d z;

		test_geom_data()
		: x(1.0, 0.0, 0.0)
		, y(0.0, 1.0, 0.0)
		, z(0.0, 0.0, 1.0)
		{

		}

		struct point_wrapper
		{
			vector2d point;
			point_wrapper(double x, double y) : point(x, y) { }

			const vector2d& p() const { return point; }
		};

		static std::string test_data_path()
		{
			// Get path of TUT test data.
			// We'll assume that it's one directory up in test_data/
			// Note - this isn't portable at all
			char tmp[32];
			char path_buf[MAXPATHLEN];

			::sprintf(tmp, "/proc/%d/exe", ::getpid());
			const int bytes = std::min(::readlink(tmp, path_buf, MAXPATHLEN), (ssize_t)(MAXPATHLEN - 1));
			if (bytes < 0)
				throw std::runtime_error("Couldn't readlink() path");	// shouldn't happen...

			const std::string cur_path = ::dirname(path_buf);

			char cur_path_buf[MAXPATHLEN];	// because dirname() doesn't take a const char*...
			std::size_t cur_path_len = cur_path.copy(cur_path_buf, cur_path.length());
			cur_path_buf[cur_path_len] = '\0';

			const std::string base_path = ::dirname(cur_path_buf);

			std::string path = base_path + "/test_data";

			return path;
		}
	};

	typedef test_group<test_geom_data> test_geom_data_t;
	test_geom_data_t geom_tests("Geometry tests");

	template <> template <>
	void test_geom_data_t::object::test<1>()
	{
		set_test_name("determine normal");
		triangle3d t(y, -x, x);	// CCW orientation
		ensure(t.normal().is_close(z, std::numeric_limits<double>::epsilon()));
	}

	template <> template <>
	void test_geom_data_t::object::test<2>()
	{
		set_test_name("flip normal");
		triangle3d t(y, x, -x);	// CW orientation
		ensure(t.normal().is_close(-z, std::numeric_limits<double>::epsilon()));
		ensure(t.flipped().normal().is_close(z, std::numeric_limits<double>::epsilon()));
	}

	template <> template <>
	void test_geom_data_t::object::test<3>()
	{
		set_test_name("is_degenerate()");

		triangle3d dg(y, vector3d(-std::numeric_limits<double>::min() * 10.0, 0.0, 0.0), vector3d(0.0, 0.0, 0.0));
		ensure(dg.is_degenerate());

		triangle3d t(y, -x, x);
		ensure(!t.is_degenerate());
	}

	template <> template <>
	void test_geom_data_t::object::test<4>()
	{
		set_test_name("bbox()");

		triangle3d t(z, (-y - x), (-y + x));
		const bbox_3<double> bbox = t.bbox();
		ensure(maths::close(bbox.extent_x(), 2.0, std::numeric_limits<double>::epsilon()));
		ensure(maths::close(bbox.extent_y(), 1.0, std::numeric_limits<double>::epsilon()));
		ensure(maths::close(bbox.extent_z(), 1.0, std::numeric_limits<double>::epsilon()));
	}

	template <> template <>
	void test_geom_data_t::object::test<5>()
	{
		set_test_name("Centroid");

		std::vector<vector2d> square;
		square.push_back(vector2d(-1, -1));
		square.push_back(vector2d(1, -1));
		square.push_back(vector2d(1, 1));
		square.push_back(vector2d(-1, 1));

		vector2d c = centroid(square.begin(), square.end());
		ensure(close(c.x(), 0.0, std::numeric_limits<double>::epsilon()));
		ensure(close(c.y(), 0.0, std::numeric_limits<double>::epsilon()));
	}

	template <> template <>
	void test_geom_data_t::object::test<6>()
	{
		set_test_name("Centroid functor");

		std::vector<test_geom_data::point_wrapper> points;
		points.push_back(test_geom_data::point_wrapper(-1, -1));
		points.push_back(test_geom_data::point_wrapper(1, -1));
		points.push_back(test_geom_data::point_wrapper(1, 1));
		points.push_back(test_geom_data::point_wrapper(-1, 1));

		vector2d c = centroid(points.begin(), points.end(),
			[](const test_geom_data::point_wrapper & pw) { return pw.p(); });

		ensure(close(c.x(), 0.0, std::numeric_limits<double>::epsilon()));
		ensure(close(c.y(), 0.0, std::numeric_limits<double>::epsilon()));
	}

	// This should probably be moved to its own file
	template <> template <>
	void test_geom_data_t::object::test<7>()
	{
		set_test_name("Primitive fitting: plane");

		// A bunch of planar points on a plane centered at the origin and orthogonal to the Z axis
		std::vector<vector3d> plane_points;
		plane_points.push_back(vector3d(1, 0, 0));
		plane_points.push_back(vector3d(0, 1, 0));
		plane_points.push_back(vector3d(-1, 0, 0));
		plane_points.push_back(vector3d(0, -1, 0));
		plane_points.push_back(vector3d(-0.25, 0.5, 0));
		plane_points.push_back(vector3d(0.1, -0.7, 0));

		vector3d plane_point(5, 5, 5);
		vector3d plane_normal(0, 0, 0);
		bool success = primitive_fitting::plane(plane_points.begin(), plane_points.end(), plane_point, plane_normal);
		ensure(success);

		ensure(plane_point.is_close(vector3d(-0.1 * 0.25, -0.1 * (1.0 / 3.0), 0), std::numeric_limits<double>::epsilon() * 20));
		ensure(plane_normal.is_close(vector3d(0, 0, 1), std::numeric_limits<double>::epsilon() * 20));
	}

	template <> template <>
	void test_geom_data_t::object::test<8>()
	{
		set_test_name("Fit point cloud to plane");

		std::ifstream ifs(test_data_path() + "/buddha-statue_pts.txt");
		std::vector<vector3d> points;

		for (std::string line ; std::getline(ifs, line) ; )
		{
			double v[3];
			std::istringstream line_ss(line);
			for (int i = 0 ; i < 3 ; i++)
			{
				std::string tok;
				std::getline(line_ss, tok, ' ');
				v[i] = std::strtod(tok.c_str(), nullptr);
			}

			points.push_back(vector3d(v[0], v[1], v[2]));
		}

		vector3d plane_point(5, 5, 5);
		vector3d plane_normal(0, 0, 0);
		ensure(primitive_fitting::plane(points.begin(), points.end(), plane_point, plane_normal));

		const double tol = 1.0e-6;	// from input
		ensure_distance(plane_point.x(), 0.0, tol);
		ensure_distance(plane_point.y(), 0.0, tol);
		ensure_distance(plane_point.z(), 0.0, tol);

		// Verified with GNU Octave
		ensure_distance(plane_normal.x(), -0.155759, tol);
		ensure_distance(plane_normal.y(),  0.907574, tol);
		ensure_distance(plane_normal.z(), -0.389935, tol);
	}
};





































