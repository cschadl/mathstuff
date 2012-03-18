/*
 * testmatrix.cpp
 *
 *  Created on: Feb 26, 2012
 *      Author: cds
 */

#include "matrix.h"

#include <iomanip>
#include <memory>
#include <limits>

#include <tut.h>

using maths::matrix;
using maths::diag_matrix;
using maths::matrix_exception;
using std::cout;
using std::setprecision;
using std::endl;

namespace tut
{
	struct test_matrix_data { };
	typedef test_group<test_matrix_data> matrix_tests;
	matrix_tests matrix_test_group("mathstuff::matrix tests");

	template <> template <>
	void matrix_tests::object::test<1>()
	{
		set_test_name("basic matrix arithmetic");

		matrix<double> A(4, 3);
		matrix<double> B(3, 2);

		A(0, 0) = 14.0;
		A(0, 1) = 9.0;
		A(0, 2) = 3.0;

		A(1, 0) = 2.0;
		A(1, 1) = 11.0;
		A(1, 2) = 15.0;

		A(2, 0) = 0.0;
		A(2, 1) = 12.0;
		A(2, 2) = 17.0;

		A(3, 0) = 5.0;
		A(3, 1) = 2.0;
		A(3, 2) = 3.0;

		B(0, 0) = 12.0;
		B(0, 1) = 25.0;

		B(1, 0) = 9.0;
		B(1, 1) = 10.0;

		B(2, 0) = 8.0;
		B(2, 1) = 5.0;

		matrix<double> C(4, 2);
		C(0, 0) = 273.0;
		C(0, 1) = 455.0;
		C(1, 0) = 243.0;
		C(1, 1) = 235.0;
		C(2, 0) = 244.0;
		C(2, 1) = 205.0;
		C(3, 0) = 102.0;
		C(3, 1) = 160.0;

		ensure(A + A == A * 2.0);
		ensure(A - A == matrix<double>(4, 3).fill(0.0));

		bool caught_matrix_exception = false;
		try
		{
			matrix<double> nop = A + B;
		}
		catch (matrix_exception&)
		{
			caught_matrix_exception = true;
		}
		ensure(caught_matrix_exception);

		ensure(C == A * B);

		bool caught_matrix_exception2 = false;
		try
		{
			matrix<double> nop = A * A;
		}
		catch (matrix_exception&)
		{
			caught_matrix_exception2 = true;
		}
		ensure(caught_matrix_exception2);
	}

	template <> template <>
	void matrix_tests::object::test<2>()
	{
		set_test_name("+=, -= and *= operators");

		matrix<double> A(4, 3);
		A(0, 0) = 14.0;
		A(0, 1) = 9.0;
		A(0, 2) = 3.0;

		A(1, 0) = 2.0;
		A(1, 1) = 11.0;
		A(1, 2) = 15.0;

		A(2, 0) = 0.0;
		A(2, 1) = 12.0;
		A(2, 2) = 17.0;

		A(3, 0) = 5.0;
		A(3, 1) = 2.0;
		A(3, 2) = 3.0;

		matrix<double> B(3, 2);
		B(0, 0) = 12.0;
		B(0, 1) = 25.0;

		B(1, 0) = 9.0;
		B(1, 1) = 10.0;

		B(2, 0) = 8.0;
		B(2, 1) = 5.0;

		matrix<double> A1t = A;
		A1t += A;
		matrix<double> A2t = A;
		A2t *= 2.0;
		ensure(A1t == A2t);
		matrix<double> A3t = A;
		ensure((A3t *= B) == (A * B));
	}

	template <> template <>
	void matrix_tests::object::test<3>()
	{
		set_test_name("matrix transpose tests");

		matrix<double> E1(3, 3);
		E1(0, 0) = 1.0;
		E1(0, 1) = 2.0;
		E1(0, 2) = 3.0;
		E1(1, 0) = 4.0;
		E1(1, 1) = 5.0;
		E1(1, 2) = 6.0;
		E1(2, 0) = 7.0;
		E1(2, 1) = 8.0;
		E1(2, 2) = 9.0;
		matrix<double> E2 = E1;
		ensure(E1 == E2);
		ensure(E2 == E1);
		ensure(E1 != E2.transpose());

		double u = 0.0;
		matrix<double> M7(3, 3);
		for (size_t i = 0 ; i < 3 ; i++)
			for (size_t j = 0 ; j < 3 ; j++)
				M7(i, j) = ++u;
		matrix<double> M7sum = M7 + M7;
		ensure(M7sum == M7 * 2.0);
		matrix<double> M7Tsum = M7 + M7.get_transpose();
		matrix<double> eM7Tsum(3, 3);
		eM7Tsum(0, 0) = 2.0;
		eM7Tsum(0, 1) = 6.0;
		eM7Tsum(0, 2) = 10.0;
		eM7Tsum(1, 0) = 6.0;
		eM7Tsum(1, 1) = 10.0;
		eM7Tsum(1, 2) = 14.0;
		eM7Tsum(2, 0) = 10.0;
		eM7Tsum(2, 1) = 14.0;
		eM7Tsum(2, 2) = 18.0;
		ensure(M7Tsum == eM7Tsum);

		matrix<double> C(4, 2);
		C(0, 0) = 273.0;
		C(0, 1) = 455.0;
		C(1, 0) = 243.0;
		C(1, 1) = 235.0;
		C(2, 0) = 244.0;
		C(2, 1) = 205.0;
		C(3, 0) = 102.0;
		C(3, 1) = 160.0;

		// we can't multiply C by itself, but we should be able
		// to multiply it by it's transpose
		matrix<double> D = C * C.get_transpose();	// tut test fails if this throws
	}

	template <> template <>
	void matrix_tests::object::test<4>()
	{
		set_test_name("SVD test 1");

		matrix<double> A(4, 5);
		matrix<double> V(5, 5);
		std::valarray<double> w;
		A(0,0) = 1.0;
		A(0,4) = 2.0;
		A(1,2) = 3.0;
		A(3,1) = 4.0;
		matrix<double> M_ = A;

		ensure(A.svd(w, V));
		ensure((A * A.get_transpose()).is_close(matrix<double>::I(A.rows()), 1.0e-12));
		ensure((V * V.get_transpose()).is_close(matrix<double>::I(A.cols()), 1.0e-12));

		std::auto_ptr< const diag_matrix<double> > pDw(matrix<double>::diag(w));
		const diag_matrix<double>& Dw = *pDw;

		matrix<double> R1 = A * Dw * V.get_transpose();

		static bool show_me = false;
		if (show_me)
		{
		cout << "U: " << A << endl << "w: " << endl << Dw << endl << "V: " << endl << V << endl;
		cout << "M_:" << endl << M_ << endl;
		cout << "R1:" << endl << R1 << endl;
		}

		ensure(M_.is_close(R1, 1.0e-15));
	}

	template <> template <>
	void matrix_tests::object::test<5>()
	{
		set_test_name("SVD test 2");

		// Note - sometimes this test will fail, but if it's run again, it will work fine.
		// I don't know why.

		matrix<float> A(4, 2);
		A(0, 0) = 2.0;
		A(0, 1) = 4.0;
		A(1, 0) = 1.0;
		A(1, 1) = 3.0;
		A(2, 0) = 2.5;	// these should work, too...
		A(2, 1) = 3.6;
		A(3, 0) = 1.7;
		A(3, 1) = 2.8;
		matrix<float> A1 = A;

		matrix<float> V(2, 2);
		std::valarray<float> w;
		ensure(A.svd(w, V));
		// The "compact" SVD that .svd() returns apparently does not necessarily
		// have the attribute that S * St = I or V * Vt = I...
		// GNU Octave seems to agree with the results that matrix::svd() returns,
		// and we can still always reconstruct the matrix...

		//ensure((A * A.get_transpose()).is_close(matrix<float>::I(A.rows()), 1.0e-8));
		//ensure((V * V.get_transpose()).is_close(matrix<float>::I(A.cols()), 1.0e-8));
		std::auto_ptr< const diag_matrix<float > > pDw(matrix<float>::diag(w));
		const diag_matrix<float>& Dw = *pDw;

		matrix<float> R1 = A * Dw * V.get_transpose();

		static bool show_me = false;
		if (show_me)
		{
		cout << "U: " << endl << A << endl <<
				"w: " << endl << Dw << endl <<
				"V: " << endl << V.get_transpose() << endl;
		cout << "A1:" << endl << A1 << endl << "R1:" << endl << R1 << endl;
		}

		if (!A1.is_close(R1, 1.0e-5))
		{
			cout << "Uh oh, A1 test failed!" << endl;
			cout << "A1: " << endl << A1 << endl;
			cout << "R1: " << endl << R1 << endl;
		}

		ensure(A1.is_close(R1, 1.0e-5));
	}

	template <> template <>
	void matrix_tests::object::test<6>()
	{
		set_test_name("row_swap(), col_swap()");

		double u = 0.0;
		matrix<double> E1(5, 3);
		for (size_t i = E1.r_begin() ; i < E1.r_end() ; i++)
			for (size_t j = E1.c_begin() ; j < E1.c_end() ; j++)
				E1(i, j) = ++u;

		matrix<double> E2 = E1;
		E2.row_swap(0, 1);
		for (size_t i = E1.c_begin() ; i < E1.c_end() ; i++)
		{
			ensure(E2(0, i) == E1(1, i));
		}

		matrix<double> E3 = E1;
		E3.col_swap(0, 2);
		for (size_t i = E1.r_begin() ; i < E1.r_end() ; i++)
		{
			ensure(E3(i, 0) == E1(i, 2));
		}

		matrix<double> E4 = E1.get_transpose();
		E1.transpose();
		E4.col_swap(3, 4);
		for (size_t i = E1.r_begin() ; i < E1.r_end() ; i++)
		{
			ensure(E4(i, 3) == E1(i, 4));
		}
	}
};
