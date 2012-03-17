/*
 * testquaternion.cpp
 *
 *  Created on: Mar 12, 2012
 *      Author: cds
 */

#include "quaternion.h"

#include <tut.h>

using maths::quaternion;

namespace tut
{
	struct quaternion_test_data {};
	typedef test_group<quaternion_test_data> quaternion_tests;
	quaternion_tests quaternion_test_group("maths::quaternion tests");

	template <> template <>
	void quaternion_tests::object::test<1>()
	{
		set_test_name("basic arithmetic properties");

		const quaternion<double> i = quaternion<double>::i();
		const quaternion<double> j = quaternion<double>::j();
		const quaternion<double> k = quaternion<double>::k();

		const quaternion<double> isq = i * i;
		const quaternion<double> jsq = j * j;
		const quaternion<double> ksq = k * k;
		const quaternion<double> ijk = i * j * k;
		const quaternion<double> ij = i * j;
		const quaternion<double> ji = j * i;
		const quaternion<double> jk = j * k;
		const quaternion<double> kj = k * j;
		const quaternion<double> ki = k * i;
		const quaternion<double> ik = i * k;

		const double tol = 1.0e-64;
		ensure(isq.is_close(-1.0, tol));
		ensure(jsq.is_close(-1.0, tol));
		ensure(ksq.is_close(-1.0, tol));
		ensure(ijk.is_close(-1.0, tol));

		ensure(ij.is_close(k, tol));
		ensure(ji.is_close(-k, tol));
		ensure(jk.is_close(i, tol));
		ensure(kj.is_close(-i, tol));
		ensure(ki.is_close(j, tol));
		ensure(ik.is_close(-j, tol));
	}
};
