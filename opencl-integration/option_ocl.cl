#pragma OPENCL EXTENSION cl_khr_fp64: enable


#define OOSQR2PI 0.398942280401433
#define q1 0.282843
#define lnS1 4.565170
#define Kk 8.000000
#define c1k -2.041241
#define c2k 2.721655
#define c3k -3.102912
#define c4k 3.345388
#define c6k 0.750000
#define c5k 0.367423

float g (float input)
{
    return  q1 * (input / (1 - input * input)) + lnS1;
    
}

float cn(float x)
{
    double t;
    if (x == 0) return .5;
    t = 1 / (1 + .2316419 * fabs(x));
    t *= OOSQR2PI * exp((float)(-.5 * x * x)) * (.31938153 + t * (-.356563782 +
                                                                  t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    return (x >= 0? 1 - t: t);
}

float fspread(float input)
{
    float l1 = log((float)input);
    float l2 = log((float)(input - Kk));
    float h = c1k * l1 + c2k * l2 + c3k;
    return (input - Kk) * cn(h) - c4k * exp((float)(c6k * l1)) * cn(h - c5k);
}

 float fx(float input)
{
    float expArg = g(input);
    float fspreadOut = fspread(exp(expArg));
    return exp(-0.5f * pow(input / (1-pow(input, 2)), 2)) * ((1.0f + pow(input, 2)) / pow(1 - input * input, 2)) * fspreadOut;
}


void reduce(
            __local  float*,
            __global float*);


__kernel void option(
                 const int          niters,
                 const float        step_size,
                 __local  float*    local_sums,
                 __global float*    partial_sums)
{
    int num_wrk_items  = get_local_size(0);
    int local_id       = get_local_id(0);
    int group_id       = get_group_id(0);
    
    float x, accum = 0.0f;
    int i,istart,iend;
    
    istart = (group_id * num_wrk_items + local_id) * niters;
    iend   = istart+niters;
    
    for(i= istart; i<iend; i++){
        x = (i+0.5f)*step_size;
        printf("x=%f\n",x);
        accum += fx(x);

    }
    local_sums[local_id] = accum;
    printf("local_sum=%f\n",accum);
    barrier(CLK_LOCAL_MEM_FENCE);
    
    reduce(local_sums, partial_sums);
}

//------------------------------------------------------------------------------
//
// OpenCL function:  reduction
//
// Purpose: reduce across all the work-items in a work-group
//
// input: local float* an array to hold sums from each work item
//
// output: global float* partial_sums   float vector of partial sums
//

void reduce(
            __local  float*    local_sums,
            __global float*    partial_sums)
{
    int num_wrk_items  = get_local_size(0);
    int local_id       = get_local_id(0);
    int group_id       = get_group_id(0);
    
    float sum;
    int i;
    
    if (local_id == 0) {
        sum = 0.0f;
        
        for (i=0; i<num_wrk_items; i++) {
            sum += local_sums[i];
        }
        printf("sums = %f\n",sum);

        partial_sums[group_id] = sum;

    }
}

