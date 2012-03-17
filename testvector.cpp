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

		ensure(i.inner_product(j) == 0);
		ensure(j.inner_product(i) == 0);
		ensure(i.inner_product(k) == 0);
		ensure(k.inner_product(i) == 0);
		ensure(j.inner_product(k) == 0);
		ensure(k.inner_product(j) == 0);
		ensure(outer_product(i, j).is_close(k,  1.0e-64));
		ensure(outer_product(j, i).is_close(-k, 1.0e-64));	// shouldn't need to do this...
		ensure((i % j) == k);
	}

	template <> template<>
	void vector_tests::object::test<2>()
	{
		set_test_name("vector length");

		vector3d v = i + j + k;
		ensure(v[0] == 1.0 && v[1] == 1.0 && v[1] == 1.0);
		ensure(v.length() > 1.0);
		const vector3d vu = v.make_unit();
		ensure(vu.length() == v.unit().length());
		ensure(vu.length() == 1.0);
		ensure(maths::close(vu.length_sq(), 1.0, 1.0e-15));
		ensure(maths::close(v.length_sq(), 1.0, 1.0e-15));
		ensure(v.length() == 1.0);
	}

	template <> template <>
	void vector_tests::object::test<3>()
	{
		set_test_name("arithmetic");

		vector3d v1(0.0, 0.5, 0.0);
		vector3d v2(0.5, 0.0, 0.0);
		vector3d v3(0.0, 0.0, 0.5);

		const double tol = 1.0e-64;
		ensure((v1 + v2).is_close(vector3d(0.5, 0.5, 0.0), tol));
		ensure((v2 + v2).is_close(vector3d(1.0, 0.0, 0.0), tol));
		ensure(((v1 + v2) - v2).is_close(v1, tol));
	}
};
