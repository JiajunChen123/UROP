#define OOSQR2PI 0.398942280401433

float helper (float input, float q1, float lnS1)
{
    float expArg = q1 * (input / (1 - input * input)) + lnS1;
    return expArg;
}

float cn(float x)
{
    float t;
    if (x == 0) return .5;
    t = 1 / (1 + .2316419 * fabs(x));
    t *= OOSQR2PI * exp((float)(-.5 * x * x)) * (.31938153 + t * (-.356563782 +
                                                           t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    return (x >= 0? 1 - t: t);
}

float fspread(float input,
             const float Kk, const float c1k, const float c2k,
               const float c3k, const float c4k, const float c5k,
               const float c6k)
{
    printf("input is: %f, Kk is: %f, c1k: %f, c2k: %f, c3k: %f, c4k: %f, c5k: %f, c6k: %f, \n", input, Kk, c1k, c2k, c3k, c4k,c5k,c6k);
    float l1 = log((float)input);
    float l2 = log((float)(input - Kk));
    float h = c1k * l1 + c2k * l2 + c3k;
    return (input - Kk) * cn(h) - c4k * exp((float)(c6k * l1)) * cn(h - c5k);
}

__kernel void fx(global float* input,
                 global float* output,
                 global float* output2,
                 global float* output3,
                 const float q1,
                 const float lnS1,
                 const float Kk, const float c1k,
                 const float c2k, const float c3k,
                 const float c4k, const float c5k,
                 const float c6k)
{
    
    size_t i = get_global_id(0);
//    printf("input is: %f\n", input[i]);
    float expArg = helper(input[i], q1, lnS1);
    float fspreadOut = fspread(exp(expArg), Kk, c1k, c2k, c3k, c4k ,c5k, c6k );
    printf("fspreadout = %f\n", fspreadOut);
    float expr = exp(-0.5 * (pow(input[i] / pow((1 - input[i]), 2), 2)));
    output[i] = expr * ((1 + pow(input[i],2))/pow(1-input[i]*input[i],2)) * fspreadOut;
    output2[i] = output[i] * expArg;
    output3[i] = output[i] * expArg * expArg;
    printf("output at %i is: %f\n", i, output);
    
//    float oneMinInput = 1 - input[i];
//    float oneDivOMI = input[i] / (oneMinInput * oneMinInput);
//    float exp = 0.5 * oneDivOMI * oneDivOMI;
//    float inputSq = input[i] * input[i];
//    float onePlusInputSq = 1 + inputSq;
//    float oneMinInputSq = 1 - inputSq;
//    float div = onePlusInputSq / (oneMinInputSq * oneMinInputSq);
//    output[i] = exp * div * fspreadOut;
    //exp(-0.5*pow(input[i]/(1-input[i]^2),2)) * ((1 + pow(input[i],2))/pow(1-input[i]*input[i],2)) * fspread[i]
}
