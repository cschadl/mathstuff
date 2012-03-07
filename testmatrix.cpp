/*
 * testmatrix.cpp
 *
 *  Created on: Feb 26, 2012
 *      Author: cds
 */

#include "matrix.h"

#include <iomanip>

using maths::matrix;
using maths::matrix_exception;
using std::cout;
using std::setprecision;
using std::endl;

void test_matrix()
{
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

	matrix<double> M1 = A * B;
	assert(C == M1);
	matrix<double> M2(4, 2);
	M2 = M1;
	assert(M2 == C);
	assert(A * B == M1);
	assert(A * B == M2);

	bool caught_exception = false;
	try
	{
		matrix<double> M3 = A + B;
	}
	catch (matrix_exception &e)
	{
		caught_exception = true;
		cout << "Caught matrix exception: " << e.what() << endl;
	}
	assert(caught_exception);

	matrix<double> M4 = M1 + M2;
	assert(M4 == M1 * 2.0);

	matrix<double> M5 = -A;
	assert(M5 == A * -1.0);
	assert(M5 == -1.0 * A);

	matrix<double> M6 = M1 - M2;
	assert(M6 == matrix<double>(M1.rows(), M1.cols()));

	// We can't multiply C by itself
	bool caught_exception2 = false;
	try
	{
		matrix<double> nop = C * C;
	}
	catch (matrix_exception &e)
	{
		caught_exception2 = true;
	}
	assert(caught_exception2);

	// But, we should be able to multiply it by it's transpose
	matrix<double> C1 = C;
	matrix<double> Ct = C1.get_transpose();
	C1.set_index_type(matrix<double>::RowMajorOne);	// make things a little more difficult...
	std::cout << "C: " << endl << C << endl;
	matrix<double> D = C1 * Ct;
	cout << "D: " << endl << setprecision(20) << D << endl;

	matrix<double> Msvd(4, 5);
	matrix<double> VTsvd(5, 5);
	std::valarray<double> w;
	Msvd(0,0) = 1.0;
	Msvd(0,4) = 2.0;
	Msvd(1,2) = 3.0;
	Msvd(3,1) = 4.0;
	matrix<double> M_ = Msvd;
	Msvd.svd(w, VTsvd);
	cout << "Msvd: " << endl << Msvd << endl;
	cout << "w: " << endl;
	for (size_t i = 0 ; i < w.size() ; i++) std::cout << w[i] << " ";
	cout << endl << endl
	     << setprecision(12) << "VTsvd: " << endl << VTsvd << endl
	     << "M * Mt" << endl << Msvd * Msvd.get_transpose() << endl
	     << "V * Vt" << endl << VTsvd * VTsvd.get_transpose() << endl;
	assert((Msvd * Msvd.get_transpose()).is_close(matrix<double>::I(Msvd.rows()), 1.0e-12));
	assert((VTsvd * VTsvd.get_transpose()).is_close(matrix<double>::I(Msvd.cols()), 1.0e-12));
	matrix<double> R1 = Msvd * matrix<double>::diag(w) * VTsvd.get_transpose();
	cout << "R1: " << endl << setprecision(12) << R1 << endl;
	assert(M_.is_close(R1, 1.0e-15));

	matrix<double> Ssvd(4, 2);
	Ssvd(0, 0) = 2.0;
	Ssvd(0, 1) = 4.0;
	Ssvd(1, 0) = 1.0;
	Ssvd(1, 1) = 3.0;
//	Ssvd(2, 0) = 2.5;
//	Ssvd(2, 1) = 3.6;
//	Ssvd(3, 0) = 1.7;
//	Ssvd(3, 1) = 2.8;
	matrix<double> S_ = Ssvd;
	cout << "S: " << endl << Ssvd << endl;

	matrix<double> V(2, 2);
	std::valarray<double> w2;
	Ssvd.svd(w2, V);
	cout << "Ssvd: " << endl << Ssvd << endl;
	cout << "w: " << endl;
	for (size_t i = 0 ; i < w2.size() ; i++) std::cout << w2[i] << " ";
	cout << endl << endl
	     << setprecision(12) << "VTsvd: " << endl << V << endl
	     << "S * St" << endl << Ssvd * Ssvd.get_transpose() << endl
	     << "V * Vt" << endl << V * V.get_transpose() << endl;
	// The "compact" SVD that .svd() returns apparently does not necessarily
	// have the attribute that S * St = I or V * Vt = I...
	// GNU Octave seems to agree with the results that matrix::svd() returns,
	// and we can still always reconstruct the matrix...

	//assert((Ssvd * Ssvd.get_transpose()).is_close(matrix<double>::I(Ssvd.rows()), 1.0e-8));
	//assert((V * V.get_transpose()).is_close(matrix<double>::I(Ssvd.cols()), 1.0e-8));
	matrix<double> R2 = Ssvd * matrix<double>::diag(w2) * V.get_transpose();
	cout << "R2: " << endl << R2 << endl;
	assert(S_.is_close(R2, 1.0e-14));

	cout << "All done!" << endl;
}
