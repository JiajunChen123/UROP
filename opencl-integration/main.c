//
//  main.c
//  opencl-device
//
//  Created by Macintosh HD on 2018/7/17.
//  Copyright © 2018年 Macintosh HD. All rights reserved.
//

#include "main.h"



void clPrintDevInfo(cl_device_id device) {
    char device_string[1024];
    
    // CL_DEVICE_NAME
    clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_string), &device_string, NULL);
    printf("  CL_DEVICE_NAME: \t\t\t%s\n", device_string);
    
    // CL_DEVICE_VENDOR
    clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(device_string), &device_string, NULL);
    printf("  CL_DEVICE_VENDOR: \t\t\t%s\n", device_string);
    
    // CL_DRIVER_VERSION
    clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(device_string), &device_string, NULL);
    printf("  CL_DRIVER_VERSION: \t\t\t%s\n", device_string);
    
    // CL_DEVICE_INFO
    cl_device_type type;
    clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(type), &type, NULL);
    if( type & CL_DEVICE_TYPE_CPU )
        printf("  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_CPU");
    if( type & CL_DEVICE_TYPE_GPU )
        printf("  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_GPU");
    if( type & CL_DEVICE_TYPE_ACCELERATOR )
        printf("  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_ACCELERATOR");
    if( type & CL_DEVICE_TYPE_DEFAULT )
        printf("  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_DEFAULT");
    
    // CL_DEVICE_MAX_COMPUTE_UNITS
    cl_uint compute_units;
    clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(compute_units), &compute_units, NULL);
    printf("  CL_DEVICE_MAX_COMPUTE_UNITS:\t\t%u\n", compute_units);
    
    // CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
    size_t workitem_dims;
    clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(workitem_dims), &workitem_dims, NULL);
    printf("  CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS:\t%zu\n", workitem_dims);
    
    // CL_DEVICE_MAX_WORK_ITEM_SIZES
    size_t workitem_size[3];
    clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(workitem_size), &workitem_size, NULL);
    printf("  CL_DEVICE_MAX_WORK_ITEM_SIZES:\t%zu / %zu / %zu \n", workitem_size[0], workitem_size[1], workitem_size[2]);
    
    // CL_DEVICE_MAX_WORK_GROUP_SIZE
    size_t workgroup_size;
    clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(workgroup_size), &workgroup_size, NULL);
    printf("  CL_DEVICE_MAX_WORK_GROUP_SIZE:\t%zu\n", workgroup_size);
    
    // CL_DEVICE_MAX_CLOCK_FREQUENCY
    cl_uint clock_frequency;
    clGetDeviceInfo(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clock_frequency), &clock_frequency, NULL);
    printf("  CL_DEVICE_MAX_CLOCK_FREQUENCY:\t%u MHz\n", clock_frequency);
    
    // CL_DEVICE_ADDRESS_BITS
    cl_uint addr_bits;
    clGetDeviceInfo(device, CL_DEVICE_ADDRESS_BITS, sizeof(addr_bits), &addr_bits, NULL);
    printf("  CL_DEVICE_ADDRESS_BITS:\t\t%u\n", addr_bits);
    
    // CL_DEVICE_MAX_MEM_ALLOC_SIZE
    cl_ulong max_mem_alloc_size;
    clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(max_mem_alloc_size), &max_mem_alloc_size, NULL);
    printf("  CL_DEVICE_MAX_MEM_ALLOC_SIZE:\t\t%u MByte\n", (unsigned int)(max_mem_alloc_size / (1024 * 1024)));
    
    // CL_DEVICE_GLOBAL_MEM_SIZE
    cl_ulong mem_size;
    clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
    printf("  CL_DEVICE_GLOBAL_MEM_SIZE:\t\t%u MByte\n", (unsigned int)(mem_size / (1024 * 1024)));
    
    // CL_DEVICE_ERROR_CORRECTION_SUPPORT
    cl_bool error_correction_support;
    clGetDeviceInfo(device, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(error_correction_support), &error_correction_support, NULL);
    printf("  CL_DEVICE_ERROR_CORRECTION_SUPPORT:\t%s\n", error_correction_support == CL_TRUE ? "yes" : "no");
    
    // CL_DEVICE_LOCAL_MEM_TYPE
    cl_device_local_mem_type local_mem_type;
    clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(local_mem_type), &local_mem_type, NULL);
    printf("  CL_DEVICE_LOCAL_MEM_TYPE:\t\t%s\n", local_mem_type == 1 ? "local" : "global");
    
    // CL_DEVICE_LOCAL_MEM_SIZE
    clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
    printf("  CL_DEVICE_LOCAL_MEM_SIZE:\t\t%u KByte\n", (unsigned int)(mem_size / 1024));
    
    // CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
    clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(mem_size), &mem_size, NULL);
    printf("  CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE:\t%u KByte\n", (unsigned int)(mem_size / 1024));
    
    // CL_DEVICE_QUEUE_PROPERTIES
    cl_command_queue_properties queue_properties;
    clGetDeviceInfo(device, CL_DEVICE_QUEUE_PROPERTIES, sizeof(queue_properties), &queue_properties, NULL);
    if( queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE )
        printf("  CL_DEVICE_QUEUE_PROPERTIES:\t\t%s\n", "CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE");
    if( queue_properties & CL_QUEUE_PROFILING_ENABLE )
        printf("  CL_DEVICE_QUEUE_PROPERTIES:\t\t%s\n", "CL_QUEUE_PROFILING_ENABLE");
    
    // CL_DEVICE_IMAGE_SUPPORT
    cl_bool image_support;
    clGetDeviceInfo(device, CL_DEVICE_IMAGE_SUPPORT, sizeof(image_support), &image_support, NULL);
    printf("  CL_DEVICE_IMAGE_SUPPORT:\t\t%u\n", image_support);
    
    // CL_DEVICE_MAX_READ_IMAGE_ARGS
    cl_uint max_read_image_args;
    clGetDeviceInfo(device, CL_DEVICE_MAX_READ_IMAGE_ARGS, sizeof(max_read_image_args), &max_read_image_args, NULL);
    printf("  CL_DEVICE_MAX_READ_IMAGE_ARGS:\t%u\n", max_read_image_args);
    
    // CL_DEVICE_MAX_WRITE_IMAGE_ARGS
    cl_uint max_write_image_args;
    clGetDeviceInfo(device, CL_DEVICE_MAX_WRITE_IMAGE_ARGS, sizeof(max_write_image_args), &max_write_image_args, NULL);
    printf("  CL_DEVICE_MAX_WRITE_IMAGE_ARGS:\t%u\n", max_write_image_args);
    
    // CL_DEVICE_IMAGE2D_MAX_WIDTH, CL_DEVICE_IMAGE2D_MAX_HEIGHT, CL_DEVICE_IMAGE3D_MAX_WIDTH, CL_DEVICE_IMAGE3D_MAX_HEIGHT, CL_DEVICE_IMAGE3D_MAX_DEPTH
    size_t szMaxDims[5];
    printf("\n  CL_DEVICE_IMAGE <dim>");
    clGetDeviceInfo(device, CL_DEVICE_IMAGE2D_MAX_WIDTH, sizeof(size_t), &szMaxDims[0], NULL);
    printf("\t\t\t2D_MAX_WIDTH\t %zu\n", szMaxDims[0]);
    clGetDeviceInfo(device, CL_DEVICE_IMAGE2D_MAX_HEIGHT, sizeof(size_t), &szMaxDims[1], NULL);
    printf("\t\t\t\t\t2D_MAX_HEIGHT\t %zu\n", szMaxDims[1]);
    clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_WIDTH, sizeof(size_t), &szMaxDims[2], NULL);
    printf("\t\t\t\t\t3D_MAX_WIDTH\t %zu\n", szMaxDims[2]);
    clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_HEIGHT, sizeof(size_t), &szMaxDims[3], NULL);
    printf("\t\t\t\t\t3D_MAX_HEIGHT\t %zu\n", szMaxDims[3]);
    clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_DEPTH, sizeof(size_t), &szMaxDims[4], NULL);
    printf("\t\t\t\t\t3D_MAX_DEPTH\t %zu\n", szMaxDims[4]);
    
    // CL_DEVICE_PREFERRED_VECTOR_WIDTH_<type>
    printf("  CL_DEVICE_PREFERRED_VECTOR_WIDTH_<t>\t");
    cl_uint vec_width [6];
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, sizeof(cl_uint), &vec_width[0], NULL);
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, sizeof(cl_uint), &vec_width[1], NULL);
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, sizeof(cl_uint), &vec_width[2], NULL);
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, sizeof(cl_uint), &vec_width[3], NULL);
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, sizeof(cl_uint), &vec_width[4], NULL);
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(cl_uint), &vec_width[5], NULL);
    printf("CHAR %u, SHORT %u, INT %u, FLOAT %u, DOUBLE %u\n\n\n",
           vec_width[0], vec_width[1], vec_width[2], vec_width[3], vec_width[4]);
}

//------------------------------------------------------------------------------
char * getKernelSource(char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        fprintf(stderr, "Error: Could not open kernel source file\n");
        exit(EXIT_FAILURE);
    }
    fseek(file, 0, SEEK_END);
    long len = ftell(file) + 1;
    rewind(file);
    
    char *source = (char *)calloc(sizeof(char), len);
    if (!source)
    {
        fprintf(stderr, "Error: Could not allocate memory for source string\n");
        exit(EXIT_FAILURE);
    }
    fread(source, sizeof(char), len, file);
    fclose(file);
    return source;
}
//------------------------------------------------------------------------------
#define INSTEPS (8192)
#define ITERS (2)


int main(int argc, const char** argv) {
//-----------------------------------------------------------------------------------------------------
    // start logs
    printf("clDeviceQuery Starting...\n\n");
    bool bPassed = true;
    
    // Get OpenCL platform ID for NVIDIA if avaiable, otherwise default
    char cBuffer[1024];
    cl_platform_id clSelectedPlatformID = NULL;
    cl_platform_id* clPlatformIDs;
    
    cl_uint num_platforms;
    cl_int ciErrNum = clGetPlatformIDs(0, NULL, &num_platforms);
    if (ciErrNum != CL_SUCCESS) {
        printf(" Error %i in clGetPlatformIDs Call!\n\n", ciErrNum);
        bPassed = false;
    }
    
    if (num_platforms == 0) {
        printf("No OpenCL platform found! \n");
        bPassed = false;
        } else {
            // if there's one platform or more, make space for ID's
            if ((clPlatformIDs = (cl_platform_id*)malloc(num_platforms * sizeof(cl_platform_id))) == NULL) {
                printf("Failed to allocate memory for cl_platform ID's!\n\n");
                bPassed = false;
            }
            
        printf("%d OpenCL Platforms found\n\n", num_platforms);
        // get platform info for each platform
        ciErrNum = clGetPlatformIDs (num_platforms, clPlatformIDs, NULL);
        for(cl_uint i = 0; i < num_platforms; ++i) {
            ciErrNum = clGetPlatformInfo (clPlatformIDs[i], CL_PLATFORM_NAME, 1024, &cBuffer, NULL);
            if(ciErrNum == CL_SUCCESS) {
                clSelectedPlatformID = clPlatformIDs[i];
                // Get OpenCL platform name and version
                ciErrNum = clGetPlatformInfo (clSelectedPlatformID, CL_PLATFORM_NAME, sizeof(cBuffer), cBuffer, NULL);
                    if (ciErrNum == CL_SUCCESS) {
                        printf(" CL_PLATFORM_NAME: \t%s\n", cBuffer);
                    } else {
                        printf(" Error %i in clGetPlatformInfo Call !!!\n\n", ciErrNum);
                        bPassed = false;
                    }
                    
                ciErrNum = clGetPlatformInfo (clSelectedPlatformID, CL_PLATFORM_VERSION, sizeof(cBuffer), cBuffer, NULL);
                if (ciErrNum == CL_SUCCESS) {
                    printf(" CL_PLATFORM_VERSION: \t%s\n\n", cBuffer);
                } else {
                    printf(" Error %i in clGetPlatformInfo Call !!!\n\n", ciErrNum);
                    bPassed = false;
                }
                    
                // Log OpenCL SDK Version # (for convenience:  not specific to OpenCL)
                    
                // Get and log OpenCL device info
                cl_uint ciDeviceCount;
                cl_device_id *devices;
                printf("OpenCL Device Info:\n\n");
                ciErrNum = clGetDeviceIDs (clSelectedPlatformID, CL_DEVICE_TYPE_ALL, 0, NULL, &ciDeviceCount);
                    
                // check for 0 devices found or errors...
                if (ciDeviceCount == 0) {
                    printf(" No devices found supporting OpenCL (return code %i)\n\n", ciErrNum);
                    bPassed = false;
                } else if (ciErrNum != CL_SUCCESS) {
                    printf(" Error %i in clGetDeviceIDs call !!!\n\n", ciErrNum);
                    bPassed = false;
                } else {
                    // Get and log the OpenCL device ID's
                    ciErrNum = clGetPlatformInfo (clSelectedPlatformID, CL_PLATFORM_NAME, sizeof(cBuffer), cBuffer, NULL);
                    printf(" %u devices found supporting OpenCL on: %s\n\n", ciDeviceCount, cBuffer);
                    char cTemp[2];
                    sprintf(cTemp, "%u", ciDeviceCount);
                    if ((devices = (cl_device_id*)malloc(sizeof(cl_device_id) * ciDeviceCount)) == NULL) {
                        printf(" Failed to allocate memory for devices !!!\n\n");
                        bPassed = false;
                    }
                    ciErrNum = clGetDeviceIDs (clSelectedPlatformID, CL_DEVICE_TYPE_ALL, ciDeviceCount, devices, &ciDeviceCount);
                    if (ciErrNum == CL_SUCCESS) {
                        for(unsigned int i = 0; i < ciDeviceCount; ++i )  {
                            printf(" ----------------------------------\n");
                            clGetDeviceInfo(devices[i], CL_DEVICE_NAME, sizeof(cBuffer), &cBuffer, NULL);
                            printf(" Device %s\n", cBuffer);
                            printf(" ---------------------------------\n");
                            clPrintDevInfo(devices[i]);
                            }
                    } else {
                        printf(" Error %i in clGetDeviceIDs call !!!\n\n", ciErrNum);
                        bPassed = false;
                    }
                }
                    
            }
        }
    }
//------------------------------------------------------------------------------------------------------
    clock_t begin = clock();

    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    cl_kernel        kernel_option;     // compute kernel
    cl_device_id        device_id;     // compute device id
    cl_int err;
    
    char *kernelsource = getKernelSource("./option_ocl.cl");             // Kernel source
    
    err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &device_id, NULL);

    // Create a compute context
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    //    checkError(err, "Creating context");
    // Create a command queue
    commands = clCreateCommandQueue(context, device_id, 0, &err);
    //    checkError(err, "Creating command queue");
    // Create the compute program from the source buffer
    program = clCreateProgramWithSource(context, 1, (const char **) & kernelsource, NULL, &err);
    //    checkError(err, "Creating program");
    // Build the program
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];
        
        printf("Error: Failed to build program executable!");
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        return EXIT_FAILURE;
    }
    
    // Create the compute kernel from the program
    kernel_option = clCreateKernel(program, "option", &err);
    
   
    float *h_psum;              // vector to hold partial sum
    int in_nsteps = INSTEPS;    // default number of steps (updated later to device preferable)
    int niters = ITERS;         // number of iterations
    u_long nsteps;
    float step_size;
    size_t nwork_groups;
    size_t max_size, work_group_size = 8;
    float option_res;
    
    cl_mem d_partial_sums;

    // Find kernel work-group size
    err = clGetKernelWorkGroupInfo (kernel_option, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    // Now that we know the size of the work-groups, we can set the number of
    // work-groups, the actual number of steps, and the step size
    nwork_groups = in_nsteps/(work_group_size*niters);
    printf("number of work group is %zu\n",nwork_groups);

    if (nwork_groups < 1)
    {
        err = clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(size_t), &nwork_groups, NULL);
        //        checkError(err, "Getting device compute unit info");
        work_group_size = in_nsteps / (nwork_groups * niters);
    }
    
    nsteps = work_group_size * niters * nwork_groups;
    printf("step size is %lu\n",nsteps);

    step_size = 1.0f/(float)nsteps;
    printf("step size is %f\n",step_size);
    h_psum = calloc(sizeof(float), nwork_groups);
    if (!h_psum)
    {
        printf("Error: could not allocate host memory for h_psum\n");
        return EXIT_FAILURE;
    }
    printf(" %ld work-groups of size %ld, %lu Integration steps\n",
           nwork_groups,
           work_group_size,
           nsteps);
    
    d_partial_sums = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * nwork_groups, NULL, &err);
    //    checkError(err, "Creating buffer d_partial_sums");
    
    // Set kernel arguments
    err  = clSetKernelArg(kernel_option, 0, sizeof(int), &niters);
    err |= clSetKernelArg(kernel_option, 1, sizeof(float), &step_size);
    err |= clSetKernelArg(kernel_option, 2, sizeof(float) * work_group_size, NULL);
    err |= clSetKernelArg(kernel_option, 3, sizeof(cl_mem), &d_partial_sums);
    //    checkError(err, "Settin kernel args");
    
    // Execute the kernel over the entire range of our 1D input data set
    // using the maximum number of work items for this device
    size_t global = nsteps / niters;
    size_t local = work_group_size;
    err = clEnqueueNDRangeKernel(commands,kernel_option,1, NULL,&global,&local,0, NULL, NULL);
    //    checkError(err, "Enqueueing kernel");
    
    err = clEnqueueReadBuffer(commands,d_partial_sums,CL_TRUE,0,sizeof(float) * nwork_groups,h_psum,0, NULL, NULL);
    //    checkError(err, "Reading back d_partial_sums");
    
    // complete the sum and compute the final integral value on the host
    option_res = 0.0f;
    for (unsigned int i = 0; i < nwork_groups; i++)
    {
        option_res += h_psum[i];
    }
    option_res *= step_size;
    
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nThe calculation ran in %lf seconds\n", time_spent);
    printf(" option = %f for %lu steps\n", option_res, nsteps);
    
    // clean up
    clReleaseMemObject(d_partial_sums);
    clReleaseProgram(program);
    clReleaseKernel(kernel_option);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    free(kernelsource);
    free(h_psum);
}

