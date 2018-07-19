//
//  Model.h
//  Quadrature
//
//  Created by Stephen Weston on 05/07/2018.
//  Copyright © 2018 Stephen Weston. All rights reserved.
//

#ifndef Model_h
#define Model_h

#include <stdio.h>

// Defined variables
#define EPS 1.0e-6

// 1 / sqrt(2 * pi)
#define OOSQR2PI 0.398942280401433

#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; });
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _b : _a; });

typedef int (*ObjectFunc) (float x, void * para, float *f);

// Used with LogNormInt
static float Kk, c1k, c2k, c3k, c4k, c5k, c6k, c7k;

/// Cumulative normal distribution, used for single underlier, lognormally distributed European exercise options
/// correct to 8 decimal places
float cn(float x);

/// Payoff function for standard spread option - assumes no indivisibles in either of the underlyings
static float fSpread(float S1);

/// Spread option pricing via lognormal integration
/// Uses a combination of Simpson and trapezoidal rules
float SpreadOptLNI(long *nSteps, float FwdPrice1,float FwdPrice2, float Strike, float Vol1, float Vol2, float Rho, float T);

/// Generalised integration of the lognormal using a modified version of Haselgrove's algorithm
/// ("A Method for Numerical Integration", Haselgrove, 1960
float LogNormInt(float *Delta, float *Gamma, float *Vega, long *nSteps, float (*fn)(float), float L, float U, float q, float S0);


#endif /* Model_h */
