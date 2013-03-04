/*
 * testgeom.cpp
 *
 *  Created on: Mar 3, 2013
 *      Author: cds
 */

#include "geom.h"
#include "tut.h"

using maths::triangle3d;
using maths::vector3d;
using maths::bbox_3;

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
};
