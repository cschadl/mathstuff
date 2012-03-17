/*
 * testquaternion.cpp
 *
 *  Created on: Mar 12, 2012
 *      Author: cds
 */

#include "quaternion.h"

#include <tut.h>

using maths::quaternion;
using maths::quatd;

namespace tut
{
	struct quaternion_test_data {};
	typedef test_group<quaternion_test_data> quaternion_tests;
	quaternion_tests quaternion_test_group("maths::quaternion tests");

	template <> template <>
	void quaternion_tests::object::test<1>()
	{
		set_test_name("basic arithmetic properties");

		const quatd i = quatd::i();
		const quatd j = quatd::j();
		const quatd k = quatd::k();

		const quatd isq = i * i;
		const quatd jsq = j * j;
		const quatd ksq = k * k;
		const quatd ijk = i * j * k;
		const quatd ij = i * j;
		const quatd ji = j * i;
		const quatd jk = j * k;
		const quatd kj = k * j;
		const quatd ki = k * i;
		const quatd ik = i * k;

		const double tol = 1.0e-128;
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

	template <> template <>
	void quaternion_tests::object::test<2>()
	{
		set_test_name("conjugate");

		const quatd i = quatd::i();
		const quatd j = quatd::j();
		const quatd k = quatd::k();

		const quatd q(1.0 / 3.0, 0.5, 1.0, 2.0 / 3.0);	// random crap
		const quatd q_star = -0.5 * (q + (i*q*i) + (j*q*j) + (k*q*k));

		ensure(q.get_conj().is_close(q_star, 1.0e-15));
		ensure(q.get_conj().is_close(q.a() - q.bi() - q.cj() - q.dk(), 1.0e-15));
	}
};
