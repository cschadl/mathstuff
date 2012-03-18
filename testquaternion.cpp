/*
 * testquaternion.cpp
 *
 *  Created on: Mar 12, 2012
 *      Author: cds
 */

#include "quaternion.h"
#include "misc.h"

#include <tut.h>

using maths::quaternion;
using maths::quatd;
using maths::close;
using maths::sqrt;

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

		const double tol = std::numeric_limits<double>::min();
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

		const double tol = std::numeric_limits<double>::epsilon();
		ensure(q.get_conj().is_close(-0.5 * (q + (i*q*i) + (j*q*j) + (k*q*k)), tol));
		ensure(q.get_conj().is_close(q.a() - q.bi() - q.cj() - q.dk(), tol));
	}

	template <> template <>
	void quaternion_tests::object::test<3>()
	{
		set_test_name("norm");

		const quatd q(0.37865, 1.28747, -3.2847005, 1.0);	// more random crap...

		const double qn = q.norm();
		const quatd q_star_q = q * q.get_conj();
		const quatd q_q_star = q.get_conj() * q;

		// q_star_q, q_q_star should be scalar
		ensure(q_star_q.is_scalar());
		ensure(q_q_star.is_scalar());

		const double tol = std::numeric_limits<double>::epsilon();
		ensure(close(qn, sqrt(q_star_q.scalar_part()), tol));
		ensure(close(qn, sqrt(q_q_star.scalar_part()), tol));
		ensure(close(qn, sqrt((q.a() * q.a()) + (q.b() * q.b()) + (q.c() * q.c()) + (q.d() * q.d())), tol));

		const double alpha = 0.5;
		ensure(close((alpha * q).norm(), (alpha * q.norm()), tol));

		const quatd q1(8.27843, -2.359783, -1.84394, 5.28497);
		ensure(close((q * q1).norm(), q.norm() * q1.norm(), tol * 100.0));	// tol a little looser due to sqrt()
	}
};
