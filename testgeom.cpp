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
#include <random>
#include <chrono>
#include <iomanip>	// for setprecision
#include <list>

using namespace maths;
using namespace std;

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

		class stupid_2d_vector : public maths::vector2d
		{
		public:
			stupid_2d_vector() : maths::vector2d(numeric_limits<double>::max()) { }
			stupid_2d_vector(double x, double y) : maths::vector2d(x, y) { }
		};

		class stupid_3d_vector : public maths::vector3d
		{
		public:
			stupid_3d_vector() : maths::vector3d(numeric_limits<double>::max()) { }
			stupid_3d_vector(double x, double y, double z) : maths::vector3d(x, y, z) { }
		};

		static string test_data_path()
		{
			// Get path of TUT test data.
			// We'll assume that it's one directory up in test_data/
			// Note - this isn't portable at all
			char tmp[32];
			char path_buf[MAXPATHLEN];

			::sprintf(tmp, "/proc/%d/exe", ::getpid());
			const int bytes = std::min(::readlink(tmp, path_buf, MAXPATHLEN), (ssize_t)(MAXPATHLEN - 1));
			if (bytes < 0)
				throw runtime_error("Couldn't readlink() path");	// shouldn't happen...

			path_buf[bytes] = '\0';

			const string cur_path = ::dirname(path_buf);

			char cur_path_buf[MAXPATHLEN];	// because dirname() doesn't take a const char*...
			size_t cur_path_len = cur_path.copy(cur_path_buf, cur_path.length());
			cur_path_buf[cur_path_len] = '\0';

			const string base_path = ::dirname(cur_path_buf);

			string path = base_path + "/test_data";

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
		ensure(t.normal().is_close(z, numeric_limits<double>::epsilon()));
	}

	template <> template <>
	void test_geom_data_t::object::test<2>()
	{
		set_test_name("flip normal");
		triangle3d t(y, x, -x);	// CW orientation
		ensure(t.normal().is_close(-z, numeric_limits<double>::epsilon()));
		ensure(t.flipped().normal().is_close(z, numeric_limits<double>::epsilon()));
	}

	template <> template <>
	void test_geom_data_t::object::test<3>()
	{
		set_test_name("is_degenerate()");

		triangle3d dg(y, vector3d(-numeric_limits<double>::min() * 10.0, 0.0, 0.0), vector3d(0.0, 0.0, 0.0));
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
		ensure(maths::close(bbox.extent_x(), 2.0, numeric_limits<double>::epsilon()));
		ensure(maths::close(bbox.extent_y(), 1.0, numeric_limits<double>::epsilon()));
		ensure(maths::close(bbox.extent_z(), 1.0, numeric_limits<double>::epsilon()));
	}

	template <> template <>
	void test_geom_data_t::object::test<5>()
	{
		set_test_name("Centroid");

		vector<vector2d> square;
		square.push_back(vector2d(-1, -1));
		square.push_back(vector2d(1, -1));
		square.push_back(vector2d(1, 1));
		square.push_back(vector2d(-1, 1));

		vector2d c = centroid(square.begin(), square.end());
		ensure(close(c.x(), 0.0, numeric_limits<double>::epsilon()));
		ensure(close(c.y(), 0.0, numeric_limits<double>::epsilon()));
	}

	template <> template <>
	void test_geom_data_t::object::test<6>()
	{
		set_test_name("Centroid functor");

		vector<test_geom_data::point_wrapper> points;
		points.push_back(test_geom_data::point_wrapper(-1, -1));
		points.push_back(test_geom_data::point_wrapper(1, -1));
		points.push_back(test_geom_data::point_wrapper(1, 1));
		points.push_back(test_geom_data::point_wrapper(-1, 1));

		vector2d c = centroid(points.begin(), points.end(),
			[](const test_geom_data::point_wrapper & pw) { return pw.p(); });

		ensure(close(c.x(), 0.0, numeric_limits<double>::epsilon()));
		ensure(close(c.y(), 0.0, numeric_limits<double>::epsilon()));
	}

	// This should probably be moved to its own file
	template <> template <>
	void test_geom_data_t::object::test<7>()
	{
		set_test_name("Primitive fitting: plane");

		// A bunch of planar points on a plane centered at the origin and orthogonal to the Z axis
		vector<vector3d> plane_points;
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

		ensure(plane_point.is_close(vector3d(-0.1 * 0.25, -0.1 * (1.0 / 3.0), 0), numeric_limits<double>::epsilon() * 20));
		ensure(plane_normal.is_close(vector3d(0, 0, 1), numeric_limits<double>::epsilon() * 20));
	}

	template <> template <>
	void test_geom_data_t::object::test<8>()
	{
		set_test_name("Fit point cloud to plane");

		ifstream ifs(test_data_path() + "/buddha-statue_pts.txt");

		// Use a std::list to make sure that primitive_fitting::plane
		// works with non random-access iterators.
		list<vector3d> points;

		for (string line ; getline(ifs, line) ; )
		{
			double v[3];
			istringstream line_ss(line);
			for (int i = 0 ; i < 3 ; i++)
			{
				string tok;
				getline(line_ss, tok, ' ');
				v[i] = strtod(tok.c_str(), nullptr);
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

	template <> template <>
	void test_geom_data_t::object::test<9>()
	{
		set_test_name("Best-fit sphere 1");

		auto test_file_path = test_data_path() + "/sphere-points.txt";

		ifstream ifs(test_file_path);
		ensure(!ifs.fail());

		vector<vector3d> points;
		for (string line ; getline(ifs, line) ; )
		{
			double v[3];
			istringstream line_ss(line);
			for (int i = 0 ; i < 3 ; i++)
			{
				string tok;
				getline(line_ss, tok, ' ');
				v[i] = strtod(tok.c_str(), nullptr);
			}

			vector3d p(v[0], v[1], v[2]);
			vector3d p_t = p + vector3d(5, 4, 3);
			points.push_back(p_t);
		}

		ensure(points.size() == 146);

		double const blah = numeric_limits<double>::max();
		vector3d sphere_center(blah);
		double sphere_radius = blah;

		ensure(primitive_fitting::sphere(points.begin(), points.end(), sphere_center, sphere_radius));

		double tol = 1.0e-6;
		ensure_distance(sphere_center.x(), 5.0, tol);
		ensure_distance(sphere_center.y(), 4.0, tol);
		ensure_distance(sphere_center.z(), 3.0, tol);
		ensure_distance(sphere_radius, 1.5, tol);
	}

	template <> template <>
	void test_geom_data_t::object::test<10>()
	{
		set_test_name("Best-fit sphere 2");

		auto test_file_path = test_data_path() + "/unit_sphere-ascii-points.txt";

		ifstream ifs(test_file_path);
		ensure(!ifs.fail());

		vector<vector3d> points;
		for (string line ; getline(ifs, line) ; )
		{
			double v[3];
			istringstream line_ss(line);
			for (int i = 0 ; i < 3 ; i++)
			{
				string tok;
				getline(line_ss, tok, ' ');
				v[i] = strtod(tok.c_str(), nullptr);
			}

			vector3d p(v[0], v[1], v[2]);
			points.push_back(p);
		}

		ensure(points.size() == 2692);

		double const blah = numeric_limits<double>::max();
		vector3d sphere_center(blah);
		double sphere_radius = blah;

		ensure(primitive_fitting::sphere(points.begin(), points.end(), sphere_center, sphere_radius));

		double const radius_tol = 1.0e-6;
		double const center_tol_sq = 2.0e-2;	// quite a bit of noise here...
		ensure(sphere_center.distance_sq(vector3d(0, 0, 0)) < center_tol_sq);
		ensure_distance(sphere_radius, 1.0, radius_tol);
	}

	template <> template <>
	void test_geom_data_t::object::test<11>()
	{
		set_test_name("2d line fit");

		// Define a 2d line from a point (origin) and a direction
		vector2d const line_pt(0.0, 0.0);
		vector2d const line_dir(1.0 / ::sqrt(2), 1.0 / ::sqrt(2));

		// Collect a random sample of points on this line, evaluated via
		// parameters in the interval (-5, 5)
		//auto line_pt_seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		size_t line_pt_seed = 0x147574b0584aa4cf;
		mt19937_64 line_pt_generator(line_pt_seed);	// maybe have separate generators for x and y
		uniform_real_distribution<double> line_pt_dist(-5.0, 5.0);

		//auto pt_noise_seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		size_t pt_noise_seed = 0x147574b0584ab42f;
		mt19937_64 pt_noise_generator(pt_noise_seed);
		uniform_real_distribution<double> line_noise_dist(-0.025, 0.025);

		const size_t num_pts = 5000;
		vector<maths::vector2d> line_points(num_pts);

		for (size_t i = 0 ; i < num_pts ; i++)
		{
			vector2d pt = line_pt + line_dir * line_pt_dist(line_pt_generator);
			pt.y() += line_noise_dist(pt_noise_generator);

			line_points.push_back(pt);
		}

		maths::line<double, 2> line2d;
		ensure(primitive_fitting::line(line_points.begin(), line_points.end(), line2d));

		ensure_equals(signbit(line2d.dir().x()), signbit(line2d.dir().y()));

		ensure_distance(line2d.point().x(), 0.0, 1.0e-2);
		ensure_distance(line2d.point().y(), 0.0, 1.0e-2);
		ensure_distance(abs(line2d.dir().x()), 1.0 / ::sqrt(2), 1.0e-3);
		ensure_distance(abs(line2d.dir().y()), 1.0 / ::sqrt(2), 1.0e-3);

		double rmsd = accumulate(line_points.begin(), line_points.end(), 0.0,
			[&line2d](double err_sq, const maths::vector2d & p)
			{
				const double dist = line2d.distance(p);
				err_sq += dist * dist;

				return err_sq;
			});

		rmsd /= (double) line_points.size();

		// What's the average expected value of a uniform distribution [-0.025, 0.025] with n = 5000?
		// Oh well, let's just give it a magic number for now.
		ensure(rmsd < 1.0e-4);
	}

	template <> template <>
	void test_geom_data_t::object::test<12>()
	{
		set_test_name("maths::adapters::create_point");

		maths::adapters::create_point<maths::vector2d> point2d_factory;
		maths::adapters::create_point<maths::vector3d> point3d_factory;

		auto point2d = point2d_factory(0.0, 1.0);
		auto point3d = point3d_factory(1.0, 2.0, 3.0);

		ensure(point2d.x() == 0.0);
		ensure(point2d.y() == 1.0);

		ensure(point3d.x() == 1.0);
		ensure(point3d.y() == 2.0);
		ensure(point3d.z() == 3.0);
	}

	template<> template<>
	void test_geom_data_t::object::test<13>()
	{
		set_test_name("maths::traits::point_traits<T>::origin()");

		auto v3d = maths::traits::point_traits<stupid_3d_vector>::origin();
		ensure(v3d.x() == 0.0);
		ensure(v3d.y() == 0.0);
		ensure(v3d.z() == 0.0);

		auto v2d = maths::traits::point_traits<stupid_2d_vector>::origin();
		ensure(v2d.x() == 0.0);
		ensure(v2d.y() == 0.0);
	}

	template<> template<>
	void test_geom_data_t::object::test<14>()
	{
		set_test_name("fit degenerate set of points to plane");

		vector3d const line_pt(1.0, 2.0, 3.0);
		vector3d const line_dir = vector3d(0.5, 0.5, 0.5).unit();

		maths::line<double, 3> line(line_pt, line_dir);

		//auto line_pt_seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		auto line_pt_seed = 0xdeadbeefdeadbeef;
		mt19937_64 line_pt_generator(line_pt_seed);
		uniform_real_distribution<double> line_pt_dist(-1.0, 1.0);

		// There isn't anything very scientific about this,
		// but adding more points tends to cause primitive_fitting::plane()
		// to actually succeed (with very small singular values in the SVD)
		size_t const num_pts = 10;

		vector<maths::vector3d> line_pts(num_pts);

		for (size_t i = 0 ; i < num_pts ; i++)
			line_pts[i] = line_pt + line_dir * line_pt_dist(line_pt_generator);

		vector3d plane_origin, plane_normal;
		ensure(!primitive_fitting::plane(line_pts.begin(), line_pts.end(), plane_origin, plane_normal));

		maths::line<double, 3> line_lsq;
		ensure(primitive_fitting::line(line_pts.begin(), line_pts.end(), line_lsq));

		double const line_lsq_pt_dist = line.distance(line_lsq.point());
		double const line_dir_dist_sq = line.dir().distance_sq(maths::abs(line_lsq.dir()));

		double const rmsd = accumulate(line_pts.begin(), line_pts.end(), 0.0,
			[&line_lsq](double err_sq, const maths::vector3d & p)
			{
				const double dist = line_lsq.distance(p);
				err_sq += dist * dist;

				return err_sq;
			});

		// B.S. tolerances
		ensure(rmsd < 1.0e-16);
		ensure(line_lsq_pt_dist < 1.0e-12);
		ensure(line_dir_dist_sq < 1.0e-16);
	}
};





































