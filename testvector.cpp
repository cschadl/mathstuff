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
		const vector3d i;
		const vector3d j;
		const vector3d k;

		vector_test_data()
		: i(1.0, 0.0, 0.0)
		, j(0.0, 1.0, 0.0)
		, k(0.0, 0.0, 1.0)
		{

		}
	};
	typedef test_group<vector_test_data> vector_tests;
	vector_tests vector_test_group("mathstuff::vector tests");

	template <> template <>
	void vector_tests::object::test<1>()
	{
		set_test_name("inner product");

		const double tol = std::numeric_limits<double>::epsilon();

		ensure(close(i.inner_product(j), 0.0, tol));
		ensure(close(j.inner_product(i), 0.0, tol));
		ensure(close(i.inner_product(k), 0.0, tol));
		ensure(close(k.inner_product(i), 0.0, tol));
		ensure(close(j.inner_product(k), 0.0, tol));
		ensure(close(k.inner_product(j), 0.0, tol));
		ensure(outer_product(i, j).is_close(k,  tol));
		ensure(outer_product(j, i).is_close(-k, tol));	// shouldn't need to do this...
		ensure((i % j).is_close(k, tol));
		ensure((j % i).is_close(-k, tol));
	}

	template <> template<>
	void vector_tests::object::test<2>()
	{
		set_test_name("vector length");

		const double tol = std::numeric_limits<double>::epsilon();

		vector3d v = i + j + k;
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
	void vector_tests::object::test<3>()
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
};
