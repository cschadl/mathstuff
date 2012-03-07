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

using maths::n_vector;
using maths::vector3d;
using maths::vector3f;

void test_vector()
{
	const vector3d& i = vector3d(1.0, 0.0, 0.0);
	const vector3d& j = vector3d(0.0, 1.0, 0.0);
	const vector3d& k = vector3d(0.0, 0.0, 1.0);

	assert(i.inner_product(j) == 0);
	assert(j.inner_product(i) == 0);
	assert(i.inner_product(k) == 0);
	assert(k.inner_product(i) == 0);
	assert(j.inner_product(k) == 0);
	assert(k.inner_product(j) == 0);
	assert(outer_product(i, j).is_close(k,  1.0e-64));
	assert(outer_product(j, i).is_close(-k, 1.0e-64));	// shouldn't need to do this...
	assert((i % j) == k);

	vector3d v = i + j + k;
	assert(v[0] == 1.0 && v[1] == 1.0 && v[1] == 1.0);
	assert(v.length() > 1.0);
	const vector3d vu = v.make_unit();
	std::cout << std::setprecision(3) << vu << std::endl;
	assert(vu.length() == v.unit().length());
	assert(vu.length() == 1.0);
	assert(maths::close(vu.length_sq(), 1.0, 1.0e-15));
	assert(maths::close(v.length_sq(), 1.0, 1.0e-15));
	assert(v.length() == 1.0);
}
