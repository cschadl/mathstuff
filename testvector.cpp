/*
 * testvector.cpp
 *
 *  Created on: Feb 26, 2012
 *      Author: cds
 */

#include "vectors.h"
#include "assert.h"

#include <iostream>
#include <iomanip>

#include <tut.h>

using maths::n_vector;
using maths::vector3d;
using maths::vector3f;
using maths::close;

namespace tut
{
	struct vector_test_data
	{
		const vector3d i_;
		const vector3d j_;
		const vector3d k_;

		vector_test_data()
		: i_(1.0, 0.0, 0.0)
		, j_(0.0, 1.0, 0.0)
		, k_(0.0, 0.0, 1.0)
		{

		}

		vector3d const& i() { return i_; }
		vector3d const& j() { return j_; }
		vector3d const& k() { return k_; }
	};
	typedef test_group<vector_test_data> vector_tests;
	vector_tests vector_test_group("mathstuff::vector tests");

	template <> template <>
	void vector_tests::object::test<1>()
	{
		set_test_name("Constructors");

		vector3d v0;
		ensure(v0[0] == 0);
		ensure(v0[1] == 0);
		ensure(v0[2] == 0);

		vector3d v1(3.5);
		ensure(v1[0] == 3.5);
		ensure(v1[1] == 3.5);
		ensure(v1[2] == 3.5);

		vector3d v2(1);
		ensure(v2[0] == 1);
		ensure(v2[1] == 1);
		ensure(v2[2] == 1);

		vector3d v1_2(v1);
		ensure(v1_2 == v1);
		ensure(v1_2 != v0);
	}

	template <> template <>
	void vector_tests::object::test<2>()
	{
		set_test_name("inner product");

		const double tol = std::numeric_limits<double>::epsilon();

		ensure(close(i().inner_product(j()), 0.0, tol));
		ensure(close(j().inner_product(i()), 0.0, tol));
		ensure(close(i().inner_product(k()), 0.0, tol));
		ensure(close(k().inner_product(i()), 0.0, tol));
		ensure(close(j().inner_product(k()), 0.0, tol));
		ensure(close(k().inner_product(j()), 0.0, tol));
		ensure(outer_product(i(), j()).is_close(k(),  tol));
		ensure(outer_product(j(), i()).is_close(-k(), tol));	// shouldn't need to do this...
		ensure((i() % j()).is_close(k(), tol));
		ensure((j() % i()).is_close(-k(), tol));
	}

	template <> template<>
	void vector_tests::object::test<3>()
	{
		set_test_name("vector length");

		const double tol = std::numeric_limits<double>::epsilon();

		vector3d v = i() + j() + k();
		ensure(v[0] == 1.0 && v[1] == 1.0 && v[1] == 1.0);
		ensure(v.length() > 1.0);
		const vector3d vu = v.make_unit();
		ensure(vu.length() == v.unit().length());	// unitize v
		ensure(close(vu.length(), 1.0, tol));
		ensure(close(vu.length_sq(), 1.0, tol * 2.0));
		ensure(close(v.length_sq(), 1.0, tol * 2.0));
		ensure(close(v.length(), 1.0, tol));
	}

	template <> template <>
	void vector_tests::object::test<4>()
	{
		set_test_name("arithmetic");

		vector3d v1(0.0, 0.5, 0.0);
		vector3d v2(0.5, 0.0, 0.0);
		vector3d v3(0.0, 0.0, 0.5);

		const double tol = std::numeric_limits<double>::epsilon();
		ensure((v1 + v2).is_close(vector3d(0.5, 0.5, 0.0), tol));
		ensure((v2 + v2).is_close(vector3d(1.0, 0.0, 0.0), tol));
		ensure(((v1 + v2) - v2).is_close(v1, tol));
	}

	template <> template <>
	void vector_tests::object::test<5>()
	{
		set_test_name("initializers");

		// These should fail to compile if there are problems

		n_vector<double, 3> v3d{0.5, 0.2, 1.0};
		ensure(v3d[0] == 0.5);
		ensure(v3d[1] == 0.2);
		ensure(v3d[2] == 1.0);

		n_vector<double, 3> v3d_2(i().x(), i().y(), i().z());
		ensure(v3d_2 == i());

		n_vector<int, 4> v4d = { 0, 1, 2, 3 };
		ensure(v4d == n_vector<int, 4>{0, 1, 2, 3 });

		n_vector<double, 3> v3d_3({0.0, j().y(), 0.0});	// have to construct from std::array
		ensure(v3d_3 == j());
	}
};
