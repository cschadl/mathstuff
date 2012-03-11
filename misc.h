/*
 * misc.h
 *
 *  Created on: Mar 4, 2012
 *      Author: cds
 */

#ifndef MISC_H_
#define MISC_H_

#include <cmath>

namespace maths
{
/**
 * templated abs function.
 */
template <typename T>
inline
T abs(T x)
{
	return x < 0 ? -x : x;
}

// do I really have to use template specialization for this?
template<>
inline
double abs(double x)
{
	// I think fabs() is faster than straight comparison
	return fabs(x);
}

template<>
inline
float abs(float x)
{
	return fabs(x);
}

template <typename T>
T sqrt(T x);

template <>
inline
double sqrt(double x)
{
	return ::sqrt(x);
}

template <>
inline
float sqrt(float x)
{
	return ::sqrtf(x);
}

/**
 *  returns true if a is close to b by a distance of tol
 */
template <typename T>
inline
bool close(T a, T b, T tol)
{
	return abs(a - b) < tol;
}

/**
 * Computes (a^2 + b^2)^(1/2) without destructive overflow
 * or underflow
 * (Numerical recipes in C sec. 2.6)
 */
template <typename T>
inline
T pythag(const T& a, const T& b)
{
	T absa = abs(a);
	T absb = abs(b);
	if (absa > absb)
		return absa * sqrt(1.0 + (absb / absa) * (absb / absa));
	else
		return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb) * (absa / absb)));
}

template <typename T>
inline
T sign(const T& a, const T& b)
{
	return b >= 0.0 ? abs(a) : -abs(b);
}

template <typename T>
inline
T max(const T& a, const T& b)
{
	return a > b ? a : b;
}

template <typename T>
inline
T min(const T& a, const T& b)
{
	return a < b ? a : b;
}

};

#endif /* MISC_H_ */
