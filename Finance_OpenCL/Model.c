//
//  Model.c
//  Quadrature
//
//  Created by Stephen Weston on 05/07/2018.
//  Copyright © 2018 Stephen Weston. All rights reserved.
//

#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include "Model.h"


float test(float x){
    return x * 2;
}

/// Cumulative normal distribution, used for single underlier, lognormally distributed European exercise options
/// correct to 8 decimal places
float cn(float x)
{
    float t;
    
    if (x == 0) return .5;
    t = 1 / (1 + .2316419 * fabs(x));
    t *= OOSQR2PI * exp(-.5 * x * x) * (.31938153 + t * (-.356563782 +
                                                         t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    return (x >= 0? 1 - t: t);
}


/// Spread option pricing via numerical integration - assumes the two assets behave lognormally.
/// Integration is done analytically in one direction
/// Pass a parameterisable object function to the lognormal integrator - this will make it easier to add
/// further option models in future - e.g. chooser options and compound options (options on options)
/// Numerical integration of the lognormal, i.e.
/// I = \int_L^U pd(S) dS func(S)
/// where pd(S) is the measure exp(-(ln(S/S0)+Vol^2*T/2)^2/(2*Vol^2*T))/(S*Vol*sqrt(2*pi*T))
/// which applies to a random walk of volatility after time T such that <S> = S0.
/// Uses Simpson's rule evaluation, after a simple change of variables.
/// func(S)    The function of S to be evaluated
/// L          The lower limit (minimum 0)
/// U          The upper limit: -ve value means infinity
/// q          volatility * sqrt(T)
/// S0         The expectation value of the underlying at time T
/// nSteps     Max no. of ordinates; if zero or negative no limit;
/// The actual number used is returned here.

/// Returns:   the value & (optionall) sensitivities:
/// Delta = {\partial I \over \partial S_0};
/// Gamma = {\partial^2 I \over \partial S_0^2}  and
/// Vega = {\partial I \over \partial q}             */

/// Payoff function for standard spread option - assumes no indivisibles in either of the underlyings
float fSpread(float S1)
{
    float h, l1, l2, SK;
    
    SK = S1 - Kk;
    l1 = log(S1);
    l2 = log(SK);
    h = c1k * l1 + c2k * l2 + c3k;
    return SK * cn(h) - c4k * exp(c6k * l1) * cn(h - c5k);
}

/// Spread option pricing via lognormal integration
/// Uses a combination of Simpson and trapezoidal rules
float SpreadOptLNI(long *nSteps, float FwdPrice1,float FwdPrice2, float Strike, float Vol1, float Vol2, float Rho, float T)
{
    float Delta, Gamma, Vega, sqrt_t, X1, X2, q1, q2, s;
    
    if (T < 0) return 0;
    if (T == 0) return max(FwdPrice1 - FwdPrice2 - Strike, 0);
    if (fabs(Rho) >= 1) return 0;
    
    sqrt_t = sqrt(T);
    q1 = Vol1 * sqrt_t;
    q2 = Vol2 * sqrt_t;
    X1 = FwdPrice1 * exp(-q1 * q1 / 2);
    X2 = FwdPrice2 * exp(-q2 * q2 / 2);
    s = 1 / sqrt(1 - Rho * Rho);
    Kk = Strike;
    
    c1k = -Rho * s / q1;
    c2k = s / q2;
    c3k = (Rho * log(X1) / q1 - log(X2) / q2) * s;
    c4k = FwdPrice2 * exp(-Rho * Rho * q2 * q2 / 2) * pow(X1, -Rho * Vol2 / Vol1);
    c5k = q2 / s;
    c6k = Rho * q2 / q1;
    return LogNormInt(&Delta, &Gamma, &Vega, nSteps, fSpread, Strike, -1, q1, FwdPrice1);
}

/// Generalised integration of the lognormal using a modified version of Haselgrove's algorithm
/// ("A Method for Numerical Integration", Haselgrove, 1960
float LogNormInt(float *Delta, float *Gamma, float *Vega, long *nSteps, float (*fn)(float), float L, float U, float q, float S0)
{
    float q2, l, u, lnS1, h, s, s1, s2, sm, sm1, sm2, os, os1, os2,
    osm, v, w, y, z, p, sum, sum1, sum2, xx, yy;
    long n_max = *nSteps, n_used = 2;
    
    if (S0 <= 0 || q <= 0) return DBL_MAX;
    q2 = q * q;
    lnS1 = log(S0) - 0.5 * q2;
    
    // Change limits to new variables for convenience ...
    if (L > 0)
    {
        v = (log(L) - lnS1) / q;
        l = (sqrt(0.25 + v * v) - 0.5) / v;
    }
    else l = -1;
    
    if (U < 0) u = 1;
    else if (U == 0) u = -1;
    else
    {
        v = (log(U) - lnS1) / q;
        u = (sqrt(0.25 + v * v) - 0.5) / v;
    }
    
    if (u < l) {v = u; u = l; l = v;}
    else if (u == l) {*nSteps = 0; *Delta = *Gamma = *Vega = 0; return 0;}
    
    if (l > -1)
    {
        xx = exp(v = q * l / (1 - l * l) + lnS1);
        yy = (*fn)(xx);
        s = exp(-0.5 * pow(l / (1 - l * l), 2.)) * (1 + l * l) /
        pow(1 - l * l, 2.0) * yy;
        s1 = s * v;
        s2 = s1 * v;
    }
    else s = s1 = s2 = 0;
    
    if (u < 1)
    {
        xx = exp(v = q * u / (1 - u * u) + lnS1);
        yy = (*fn)(xx);
        s += (w = exp(-0.5 * pow(u / (1 - u * u), 2.)) * (1 + u * u) /
              pow(1 - u * u, 2.0) * yy);
        s1 += (w *= v);
        s2 += w * v;
    }
    
    s *= (h = u - l) / 2;
    s1 *= h / 2;
    s2 *= h / 2;
    sm = 0;
    sm1 = 0;
    sm2 = 0;
    
    for (n_used = 2; 2 * n_used - 1 <= n_max || n_max <= 0; n_used += n_used - 1)
    {
        os = s;
        os1 = s1;
        os2 = s2;
        osm = sm;
        for (z = l + h / 2, sum = sum1 = sum2 = 0; z < u; z += h)
        {
            v = z * (y = 1 / (1 - z * z));
            w = v * v;
            p = q * v + lnS1;
            xx = exp(p);
            yy = (*fn)(xx);
            sum += (v = exp(-.5 * w) * (y * y + w) * yy);
            sum1 += (v *= p);
            sum2 += v * p;
        }
        s = s / 2 + (h /= 2) * sum;  // include midpoints under trapezoid rule
        s1 = s1 / 2 + h * sum1;
        s2 = s2 / 2 + h * sum2;
        sm = (4 * s - os) / 3;       // convert to Simpson's rule
        sm1 = (4 * s1 - os1) / 3;
        sm2 = (4 * s2 - os2) / 3;
        if (fabs(sm - osm) <  1e-9 + EPS * fabs(osm) && n_used >= 33) break;
    }
    
    *nSteps = n_used;
    sm *= OOSQR2PI;
    sm1 *= OOSQR2PI;
    sm2 *= OOSQR2PI;
    w = lnS1 / q2;
    v = ((w + 1) * lnS1 - 1) * sm - (1 + 2 * w) * sm1 + sm2 / q2;
    y = S0 * q2;
    *Delta = (sm1 - lnS1 * sm) / y;
    *Gamma = v / (S0 * y);
    *Vega = v / q;
    return sm;
}
