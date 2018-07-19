//
//  main.c
//  Quadrature
//
//  Created by Stephen Weston on 05/07/2018.
//  Copyright © 2018 Stephen Weston. All rights reserved.
//

#include <stdio.h>
#include "Model.h"
#include <time.h>

int main()
{
    int i, numOpts;
    long numSteps;
    float spresult;
    float FwdPrice1, FwdPrice2, Strike, Vol1, Vol2,Rho, T;
    
    printf("Hello, World!\n");
    clock_t begin = clock();
    /// Spread option pricing via lognormal integration
    /// Uses a combination of Simpson and trapezoidal rules
    FwdPrice1 = 100;
    FwdPrice2 = 105;
    Strike = 8.0;
    Vol1 = 0.20;
    Vol2 = 0.3;
    Rho = 0.5;
    T = 2.0;
    spresult = SpreadOptLNI(&numSteps, FwdPrice1, FwdPrice2, Strike, Vol1, Vol2, Rho, T);
    printf("Spread option price = %15.15f.\n", spresult);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("time spent =%f\n",time_spent);
    // Spread option value for the above settings should be: 9.204206973451248
    
    // Now loop over more option prices to verify and test validity
    // numOpts =
    // for (i = 0; i < numOpts; i++)
    //{
    //   spresult = 0.0;
    //  spresult = SpreadOptLNI(&numSteps, FwdPrice1, FwdPrice2, Strike, Vol1, Vol2, Rho, T);
    // printf("Spread option price = %15.15f.\n", spresult);
    //}
    
    
    return 0;
}
