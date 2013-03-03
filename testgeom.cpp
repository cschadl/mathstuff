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

namespace tut
{
	struct test_geom_data
	{

	};

	typedef test_group<test_geom_data> test_geom_data_t;
	test_geom_data_t geom_tests("Geometry tests");

	template <> template <>
	void test_geom_data_t::object::test<1>()
	{
		set_test_name("determine normal");
		triangle3d t1(vector3d(0.0, 1.0, 0.0), vector3d(-1.0, 0.0, 0.0), vector3d(1.0, 0.0, 0.0));
		ensure(t1.normal().is_close(vector3d(0.0, 0.0, 1.0), std::numeric_limits<double>::epsilon()));
	}
};
