#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <cuComplex.h>
#include <cublas_v2.h>
#include <cusolverDn.h>


int nDevices = -1;

cublasHandle_t handle_cublas;

cusolverDnHandle_t handle_cusolver;


typedef size_t devptr_t;
  // device pointer arrays

devptr_t *dev_ptrs;

// Fortran interface macro funx_
#define F90(name) name##_

#if CUDART_VERSION >= 12000
static const char *
cublasGetErrorString(cublasStatus_t err) {
  return cublasGetStatusString(err);
}
#else
static const char *
cublasGetErrorString(cublasStatus_t err) {
  switch (err) {
    case CUBLAS_STATUS_SUCCESS:
      return "CUBLAS_STATUS_SUCCESS";
    case CUBLAS_STATUS_NOT_INITIALIZED:
      return "CUBLAS_STATUS_NOT_INITIALIZED";
    case CUBLAS_STATUS_ALLOC_FAILED:
      return "CUBLAS_STATUS_ALLOC_FAILED";
    case CUBLAS_STATUS_INVALID_VALUE:
      return "CUBLAS_STATUS_INVALID_VALUE";
    case CUBLAS_STATUS_ARCH_MISMATCH:
      return "CUBLAS_STATUS_ARCH_MISMATCH";
    case CUBLAS_STATUS_MAPPING_ERROR:
      return "CUBLAS_STATUS_MAPPING_ERROR";
    case CUBLAS_STATUS_EXECUTION_FAILED:
      return "CUBLAS_STATUS_EXECUTION_FAILED";
    case CUBLAS_STATUS_INTERNAL_ERROR:
      return "CUBLAS_STATUS_INTERNAL_ERROR";
    case CUBLAS_STATUS_NOT_SUPPORTED:
      return "CUBLAS_STATUS_NOT_SUPPORTED";
    default:
      return "UNKNOWN_ERROR";
  }
}
#endif  // CUDART_VERSION >= 12000

static const char *cusolverGetErrorString(cusolverStatus_t status) {
  switch (status) {
    case CUSOLVER_STATUS_SUCCESS:
      return "CUSOLVER_STATUS_SUCCESS";
    case CUSOLVER_STATUS_NOT_INITIALIZED:
      return "CUSOLVER_STATUS_NOT_INITIALIZED";
    case CUSOLVER_STATUS_ALLOC_FAILED:
      return "CUSOLVER_STATUS_ALLOC_FAILED";
    case CUSOLVER_STATUS_INVALID_VALUE:
      return "CUSOLVER_STATUS_INVALID_VALUE";
    case CUSOLVER_STATUS_ARCH_MISMATCH:
      return "CUSOLVER_STATUS_ARCH_MISMATCH";
    case CUSOLVER_STATUS_EXECUTION_FAILED:
      return "CUSOLVER_STATUS_EXECUTION_FAILED";
    case CUSOLVER_STATUS_INTERNAL_ERROR:
      return "CUSOLVER_STATUS_INTERNAL_ERROR";
    case CUSOLVER_STATUS_MAPPING_ERROR:
      return "CUSOLVER_STATUS_MAPPING_ERROR";
    default:
      return "UNKNOWN_ERROR";
  }
}

inline void cudaError(cudaError_t err, const char *file, const int line) {
  if (err != cudaSuccess) {
    printf("*** CUDA Error in %s at line %d : %s -code: %i \n", file, line,
           cudaGetErrorString(err), err);
    cudaDeviceReset();
    exit(-1);
  }
}

inline void checkCudaError(const char *file, const int line)
{
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    printf("*** CUDA Error in %s, line %i\n", file, line);
    printf("%s\n", cudaGetErrorString(error));
    cudaDeviceReset();
    exit(-1);
  }
}

inline void cublasError(cublasStatus_t err, const char *file, const int line) {
  if (err != CUBLAS_STATUS_SUCCESS) {
    printf("*** CUBLAS Error in %s at line %d : %s -code: %i\n", file, line,
           cublasGetErrorString(err), err);
    cudaDeviceReset();
    exit(-1);
  }
}

inline void cusolverError(cusolverStatus_t err, const char *file,
                          const int line) {
  if (err != CUSOLVER_STATUS_SUCCESS) {
    printf("*** CUSOLVER Error in %s at line %d : %s -code: %i \n", file, line,
           cusolverGetErrorString(err), err);
    exit(-1);
  }
}


#define CHECK_FOR_ERROR() (checkCudaError(__FILE__, __LINE__))
#define HANDLE_CUDA( err ) (cudaError( err, __FILE__, __LINE__ ))
#define HANDLE_CUBLAS( err ) (cublasError( err, __FILE__, __LINE__ ))
#define HANDLE_CUSOLVER( err ) (cusolverError( err, __FILE__, __LINE__ ))

// conj complex
__global__ void cmplxConj_gpu( cuDoubleComplex *vec1, cuDoubleComplex *vec2, const int n )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  while( tid < n){
    vec2[tid] = cuConj(vec1[tid]);
    tid += stride;
  }
}
// init cuda and cublas cusolver
//extern "C" void F90(cudalib_init)(int *dev_id) {
extern "C" void F90(cudalib_init)(void) {
  //HANDLE_CUDA(cudaSetDevice(*dev_id));
  //HANDLE_CUDA(cudaDeviceReset());
  HANDLE_CUBLAS(cublasCreate(&handle_cublas));
  HANDLE_CUSOLVER(cusolverDnCreate(&handle_cusolver));
}

extern "C" void F90(cudalib_free)(void) {
  if (handle_cublas) cublasDestroy(handle_cublas);
  if (handle_cusolver) cusolverDnDestroy(handle_cusolver);
}

extern "C" void F90(set_gpu)(int *dev_id) {
  HANDLE_CUDA(cudaSetDevice(*dev_id));
}

extern "C" void F90(alloc_cublas)(const int *n, const int *elemSize,
                                  devptr_t *devicePtr) {
  void *ptr;
  HANDLE_CUDA(cudaMalloc((void **)&ptr, (*n) * (*elemSize)));
  *devicePtr = (devptr_t)(ptr);
}

extern "C" void F90(free_cublas)(devptr_t *devicePtr) {
  HANDLE_CUDA(cudaFree(devicePtr));
}

extern "C" void F90(getgpusavailable)(int *result){
  HANDLE_CUDA(cudaGetDeviceCount(result));
  nDevices = *result;
}

extern "C" void F90(malloc_gpu)(long *ptr, int *length, int *elemSize) {
  double *temp;
  long l = *length;
  long e = *elemSize;
  long s = l * e;
  HANDLE_CUDA(cudaMalloc((void **)&temp, s));
  *ptr = (long)(temp);
}

extern "C" void F90(memcpy_gpu)(double *host_ptr, long *device_ptr,
                                 int *length, int *elemSize, int *flag) {
  long l = *length;
  long e = *elemSize;
  long size = l * e;
  double *ptr_d = (double *)(*device_ptr);
  if (*flag == 1) {
    // copy from CPU to GPU
    HANDLE_CUDA(cudaMemcpy(ptr_d, host_ptr, size, cudaMemcpyHostToDevice));
  } else if (*flag == 2) {
    // copy from GPU to CPU
    HANDLE_CUDA(cudaMemcpy(host_ptr, ptr_d, size, cudaMemcpyDeviceToHost));
  } else {
    printf(" Wrong parameter.... \n  Copy flag : %d \n ", *flag);
    exit(-1);
  }
}

extern "C" void F90(free_gpu)(long *ptr) {
  double *temp;
  temp = (double *)(*ptr);
  CHECK(cudaFree(temp) == cudaSuccess);
}


extern "C" void F90(zhegvd_gpu)(int *type, char *job, char *uplo_h, int *len,
                                long *dev_A, long *dev_B, long *dev_eigs,
                                int *info) {
  cuDoubleComplex *A = (cuDoubleComplex *)(*dev_A);
  cuDoubleComplex *B = (cuDoubleComplex *)(*dev_B);
  double *eigs = (double *)(*dev_eigs);

  cuDoubleComplex *d_work = NULL;
  int *devInfo = NULL;

  cusolverEigType_t itype;
  cusolverEigMode_t jobz;
  cublasFillMode_t uplo;

  char zjob = *job;
  char uplo_c = *uplo_h;
  int n = *len;
  int lda = n;
  int ldb = n;

  if (zjob == 'v' || zjob == 'V')
    jobz = CUSOLVER_EIG_MODE_VECTOR;
  else
    jobz = CUSOLVER_EIG_MODE_NOVECTOR;

  if (uplo_c == 'U' || uplo_c == 'u')
    uplo = CUBLAS_FILL_MODE_UPPER;
  else
    uplo = CUBLAS_FILL_MODE_LOWER;

  if (*type == 1) itype = CUSOLVER_EIG_TYPE_1;
  if (*type == 2) itype = CUSOLVER_EIG_TYPE_2;
  if (*type == 3) itype = CUSOLVER_EIG_TYPE_3;

  int lwork = -1;
  HANDLE_CUSOLVER(cusolverDnZhegvd_bufferSize(handle_cusolver,
                                              itype,
                                              jobz,
                                              uplo,
                                              n,
                                              A,
                                              lda,
                                              B,
                                              ldb,
                                              eigs,
                                              &lwork));

  HANDLE_CUDA(cudaMalloc((void **)&d_work, sizeof(cuDoubleComplex) * lwork));

  HANDLE_CUDA(cudaMalloc((void **)&devInfo, sizeof(int)));

  HANDLE_CUSOLVER(cusolverDnZhegvd(handle_cusolver,
                                   itype,
                                   jobz,
                                   uplo,
                                   n,
                                   A,
                                   lda,
                                   B,
                                   ldb,
                                   eigs,
                                   d_work,
                                   lwork,
                                   devInfo));

  HANDLE_CUDA(cudaMemcpy(info, devInfo, sizeof(int), cudaMemcpyDeviceToHost));

  HANDLE_CUDA(cudaFree(d_work));
  HANDLE_CUDA(cudaFree(devInfo));
}
