
#include "calculate_numbers.h"
#include <time.h>




/* Create program from a file and compile it */
cl_program build_program(cl_context ctx, cl_device_id dev, const char* filename) {
    
    cl_program program;
    FILE *program_handle;
    char *program_buffer, *program_log;
    size_t program_size, log_size;
    int err;
    
    /* Read program file and place content into buffer */
    program_handle = fopen(filename, "r");
    if(program_handle == NULL) {
        perror("Couldn't find the program file");
        exit(1);
    }
    fseek(program_handle, 0, SEEK_END);
    program_size = ftell(program_handle);
    rewind(program_handle);
    program_buffer = (char*)malloc(program_size + 1);
    program_buffer[program_size] = '\0';
    fread(program_buffer, sizeof(char), program_size, program_handle);
    fclose(program_handle);
    
    /* Create program from file
     Read the file's content into the program_buffer, and then call clCreateProgramWithSource.
     */
    
    program = clCreateProgramWithSource(ctx, 1,
                                        (const char**)&program_buffer, &program_size, &err);
    if(err < 0) {
        perror("Couldn't create the program");
        exit(1);
    }
    free(program_buffer);
    
    char arg[] = "-I ";
    char dir[FILENAME_MAX];
    getcwd(dir, FILENAME_MAX);
    
    char* linkerArg;
    linkerArg = malloc(strlen(dir) + strlen(arg));
    strcpy(linkerArg, arg);
    strcat(linkerArg, dir);
    printf("dir is: %s\n", linkerArg);
    
    /* Build program */
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if(err < 0) {
        
        /* Find size of log and print to std output */
        clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG,
                              0, NULL, &log_size);
        program_log = (char*) malloc(log_size + 1);
        program_log[log_size] = '\0';
        clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG,
                              log_size + 1, program_log, NULL);
        printf("%s\n", program_log);
        free(program_log);
        exit(1);
    }
    
    return program;
}

static void createDeviceStructures(cl_command_queue *commands, cl_context *context, cl_device_id *device_id, int *err) {
    *err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, device_id, NULL);
    if (*err != CL_SUCCESS)
    {
        printf("Error: Failed to create a device group!\n");
        exit(1);
    }
    
    // Create a compute context
    *context = clCreateContext(0, 1, device_id, NULL, NULL, err);
    if (!*context)
    {
        printf("Error: Failed to create a compute context!\n");
        exit(1);
    }
    
    // Create a command commands
    *commands = clCreateCommandQueue(*context, *device_id, 0, err);
    if (!*commands)
    {
        printf("Error: Failed to create a command commands!\n");
        exit(1);
    }
}

static void setupKernelWorkGroups(cl_command_queue commands, cl_device_id device_id, int *err, cl_kernel kernel, long n_used) {
    size_t global;
    size_t local;
    
    // Setting local: Get the maximum work group size for executing the kernel on the device
    *err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
    if (*err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve kernel work group info! %d\n", *err);
        exit(1);
    }

    local = min(local, n_used);

    // Execute the kernel over the entire range of our 1d input data set
    // using the maximum number of work group items for this device
    global = n_used;
    *err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
    if (*err)
    {
        printf("Error: Failed to execute kernel!\n");
        exit(1);
    }
}

int main(int argc, char** argv)
{
    //----------------------------- Program setup -------------------------------//
    int err;                            // error code returned from api calls
    cl_device_id device_id;             // compute device id
    cl_context context;                 // compute context
    cl_command_queue commands;          // compute command queue
    cl_program program;                 // compute program
    cl_kernel kernel;                   // compute kernel
    cl_mem m_inputPoints;               // device memory used for the input array
    cl_mem m_output;                    // device memory used for the output array
    cl_mem m_output2;                   // device memory used for the output array
    cl_mem m_output3;                   // device memory used for the output array

    
    
    // Connect to a compute device
    createDeviceStructures(&commands, &context, &device_id, &err);
    
    // Build program
    program = build_program(context, device_id, PROGRAM_FILE);
    
    //----------------------------- Program setup -------------------------------//
    clock_t begin = clock();

    //----------------------------- Quadrature setup -------------------------------//
    
    float * inputPoints;
    float  Gamma ;
    float  Delta;
    float  Vega;
    float FwdPrice1, FwdPrice2, Strike, Vol1, Vol2, Rho, T, L, u, l, v, osm, sm, U, h;
    float osm2,sm2,osm3,sm3,sum,sum2,sum3;
    long num_points;
    FwdPrice1 = 100.0;
    FwdPrice2 = 105.0;
    Strike = 8.0;
    Vol1 = 0.2;
    Vol2 = 0.3;
    Rho = 0.5;
    T = 2.0;
    
    float sqrt_t, X1, X2, q1, q2, s;
    sqrt_t = sqrt(T);
    q1 = Vol1 * sqrt_t;
    q2 = Vol2 * sqrt_t;
    X1 = FwdPrice1 * exp(-q1 * q1 / 2);
    X2 = FwdPrice2 * exp(-q2 * q2 / 2);
    s = 1 / sqrt(1 - Rho * Rho);
    
    float lnS1 = log(FwdPrice1) - 0.5 * pow(q1, 2);
    printf("q1=%f\n,lnS1=%f\n",q1,lnS1);
    Kk = Strike;
    c1k = -Rho * s / q1;
    c2k = s / q2;
    c3k = (Rho * log(X1) / q1 - log(X2) / q2) * s;
    c4k = FwdPrice2 * exp(-Rho * Rho * q2 * q2 / 2) * pow(X1, -Rho * Vol2 / Vol1);
    c5k = q2 / s;
    c6k = Rho * q2 / q1;
    
    printf("Kk = %f\n, c1k =%f\n,c2k=%f\n,c3k=%f\n,c4k=%f\n,c5k=%f\n,c6k=%f\n",Kk,c1k,c2k,c3k,c4k,c5k,c6k);
    
    long n_max = 1000, n_used;
    
    L = Strike;
    U = -1;
    u = 1;
    if (L > 0)
    {
        v = (log(L) - lnS1) / q1;
//        printf("L is: %f, lns1 is :%f, q1 is %f, v is %f\n", L, lnS1, q1, v);
        l = (sqrt(0.25 + v * v) - 0.5) / v;
    }
    else l = -1;
    if (u < l) {v = u; u = l; l = v;}
    osm = 0;
    osm2 = 0;
    osm3 = 0;
    sm = 0;
    sm2 = 0;
    sm3 = 0;
    h = (u-l);
    
    for (n_used = 3; 2 * n_used - 1 <= n_max || n_max <= 0; n_used += n_used - 1)
    {
        num_points = (n_used - 1) / 2;
        inputPoints = malloc(sizeof(float) * num_points);     // original data set given to device
        osm = sm;
        osm2 = sm2;
        osm3 = sm3;
        sum=sum2=sum3=0;
//        printf("h is %f\n",h);
//        printf("there are %ld points needed to be computed \n",n_used-1);
        //create the input array with size = n_used
        int i = 0;
        for (float z = l+h/2; z< u; z+=h){
            inputPoints[i] = z;
//            printf("input point is %f\n", inputPoints[i]);
            i++;

        }
        //-------------------------Kernel and input/output buffer setup---------------------------------//
        
        float results[num_points];          // results returned from device
        float results2[num_points];          // results returned from device
        float results3[num_points];          // results returned from device
        // Create the input and output arrays in device memory for our calculation
        m_inputPoints = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(float) * (num_points), NULL, NULL);
        m_output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * (num_points), NULL, NULL);
        m_output2 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * (num_points), NULL, NULL);
        m_output3 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * (num_points), NULL, NULL);
        
        if (!m_inputPoints || !m_output)
        {
            printf("Error: Failed to allocate device memory!\n");
            exit(1);
        }
        
        // Write the data set into the input array in device memory
        err = 0;
        err = clEnqueueWriteBuffer(commands, m_inputPoints, CL_TRUE, 0, sizeof(float) * (num_points), inputPoints, 0, NULL, NULL);
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to write to source array!\n");
            exit(1);
        }
        
        // Create the compute kernel in the program we wish to run
        kernel = clCreateKernel(program, "fx", &err);
        if (!kernel || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel!\n");
            exit(1);
        }
        
        // Set the arguments to our compute kernel
        err = 0;
        err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &m_inputPoints);
        err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &m_output);
        err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &m_output2);
        err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &m_output3);
        err |= clSetKernelArg(kernel, 4, sizeof(float), &q1);
        err |= clSetKernelArg(kernel, 5, sizeof(float), &lnS1);
        err |= clSetKernelArg(kernel, 6, sizeof(float), &Kk);
        err |= clSetKernelArg(kernel, 7, sizeof(float), &c1k);
        err |= clSetKernelArg(kernel, 8, sizeof(float), &c2k);
        err |= clSetKernelArg(kernel, 9, sizeof(float), &c3k);
        err |= clSetKernelArg(kernel, 10, sizeof(float), &c4k);
        err |= clSetKernelArg(kernel, 11, sizeof(float), &c5k);
        err |= clSetKernelArg(kernel, 12, sizeof(float), &c6k);
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to set kernel arguments! %d\n", err);
            exit(1);
        }
    
        setupKernelWorkGroups(commands, device_id, &err, kernel, num_points);
        
        // Execute the kernel with the given task
        err = clEnqueueTask(commands, kernel, 0, NULL, NULL);
        if (err)
        {
            printf("Error: Failed to execute kernel!\n");
            exit(1);
        }
        
        // Wait for the commands to finish
        clFlush(commands);
        clFinish(commands);
        
        // Read back the results from the device to verify the output
        err = clEnqueueReadBuffer(commands, m_output, CL_TRUE, 0, sizeof(float) * (num_points), results, 0, NULL, NULL );
        err = clEnqueueReadBuffer(commands, m_output2, CL_TRUE, 0, sizeof(float) * (num_points), results2, 0, NULL, NULL );
        err = clEnqueueReadBuffer(commands, m_output3, CL_TRUE, 0, sizeof(float) * (num_points), results3, 0, NULL, NULL );

        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to read output array! %d\n", err);
            exit(1);
        }
        
        clReleaseMemObject(m_inputPoints);
        clReleaseMemObject(m_output);
        clReleaseMemObject(m_output2);
        clReleaseMemObject(m_output3);

        clReleaseKernel(kernel);
        
        // ----------------------- Kernel end ----------------------------//
        
        for(int i = 0; i < num_points; i++){
                sum += results[i];
//                printf("output point is %f\n", results[i]);
            
                sum2 += results2[i];
                sum3 += results3[i];
            }
        sm = sm / 2 + (h /= 2) * sum;
        sm2 = sm2 / 2 + h * sum2;
        sm3 = sm3 / 2 + h * sum3;
        
        sm = (4 * sm - osm) / 3;

        sm2 = (4 * sm2 - osm2) / 3;
        sm3 = (4 * sm3 - osm3) / 3;
//        printf("sm is %f, osm is %f, nused is %lu, fabs is %f\n",sm, osm, n_used, fabs(sm-osm));
        if(fabs(sm - osm) < (1e-9 + EPS * fabs(osm)) && n_used >= 33) break;
        
    }
    printf("n_used is %lu\n", n_used);

    sm *= OOSQR2PI;
    sm2 *= OOSQR2PI;

    sm3 *= OOSQR2PI;

    float w = lnS1 / (q1*q1);
    v = ((w + 1) * lnS1 - 1) * sm - (1 + 2 * w) * sm2 + sm3 / (q1*q1);
    float y = FwdPrice1 * (q1*q1);
    Delta = (sm2 - lnS1 * sm) / y;
    Gamma = v / (FwdPrice1 * y);
    Vega = v / q1;
    printf("The option price is %f\n", sm);
    printf("The Greeks Delta is %f\n", Delta);
    printf("The Greeks Gamma is %f\n", Gamma);
    printf("The Greeks Vega is %f\n", Vega);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("time spent =%f\n",time_spent);

    // Shutdown and cleanup
    clReleaseProgram(program);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    
    return 0;
}
