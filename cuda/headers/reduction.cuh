#ifndef REDUCTION_CUH_INCLUDED
#define REDUCTION_CUH_INCLUDED 

#include "cuda.cuh"


template <typename T>
struct Identity
{
    __device__ __forceinline__
    T operator()(const T &a) const {
        return a;
    }
};

template <typename T>
struct Sum
{
    __device__ __forceinline__
    T operator()(const T &a, const T &b) const {
        return a + b;
    }
};


// Array reduction, adaption from: 
//   Mark Harris. "Optimizing CUDA." SC07: High Performance Computing With CUDA (2007).


/*!
 * \brief Kernel that performs the reduction.
 *
 * \see reduce_array
 */
template <unsigned int num_threads, typename T_in, typename T_out, typename T_op, typename T_init>
__global__ void reduce_kernel(
        const T_in * __restrict__ data,
        T_out * __restrict__ out,
        const int n,
        const T_op &op,
        const T_init &init,
        T_out null
        )
{
    #define REDUCE_SYNC(n_, op) \
        if (num_threads >= 2 * n_) { \
            if (tid < n_) { \
                s_data[tid] = op(s_data[tid], s_data[tid + n_]); \
            } \
            __syncthreads(); \
        }

    #define REDUCE_SIMD(n_, op) \
        if (num_threads >= 2 * n_) { \
            s_data[tid] = op(s_data[tid], s_data[tid + n_]); \
            __syncthreads(); /* otherwise, the shared data must be volatile */ \
        }

    __shared__ T_out s_data[num_threads];
    const int tid = threadIdx.x;
    int i = blockIdx.x * num_threads + tid;

    if (i < n) {
        s_data[tid] = init(data[i]);
    }
    else {
        s_data[tid] = null;
    }

    __syncthreads();
    REDUCE_SYNC(512, op)
    REDUCE_SYNC(256, op)
    REDUCE_SYNC(128, op)
    REDUCE_SYNC(64, op)
    if (tid < 32) {
        REDUCE_SIMD(32, op)
        REDUCE_SIMD(16, op)
        REDUCE_SIMD(8, op)
        REDUCE_SIMD(4, op)
        REDUCE_SIMD(2, op)
    }
    if (tid == 0) {
        out[blockIdx.x] = op(s_data[0], s_data[1]);
    }

    #undef REDUCE_SYNC
    #undef REDUCE_SIMD
}


/*!
 * \brief Run the reduction kernel.
 *
 * \see reduce_array
 */
template <typename T_in, typename T_out, typename T_op, typename T_init>
__inline__
static void run_reduce_kernel(
        const T_in * __restrict__ data,
        T_out * __restrict__ out,
        const int n,
        const int block_size,
        const int grid_size,
        const T_op &op,
        const T_init &init,
        const T_out null)
{
    #define CASE(num_threads) \
        case num_threads: \
            reduce_kernel<num_threads, T_in, T_out, T_op, T_init> \
                         <<<grid_size, num_threads>>> \
                         (data, out, n, op, init, null); \
            break;

    switch (block_size) {
        CASE(1024)
        CASE(512)
        CASE(256)
        CASE(128)
        CASE(64)
        CASE(32)
        CASE(16)
        CASE(8)
        CASE(4)
        CASE(2)
        CASE(1)
    }

    #undef CASE
}


/*!
 * \brief Reduction of an array.
 *
 * Execute the reduction of a unidimensional array. Support distinct
 * input and output types, allowing an operation to transform the input
 * type objects to the output type objects.
 *
 * \tparam num_thread Size of a block.
 * \tparam T_in       Input type.
 * \tparam T_out      Output type.
 * \tparam T_op       Functor with a member T_out oprerator(const T_out, const T_out)
 * \tparam T_init     Functor with a member T_out oprerator(const T_in)
 *
 * \param data       Input array.
 * \param out        Location where the result will be stored.
 * \param n          Length of the array.
 * \param block_size Number of threads in a block.
 * \param op         Reduction operator.
 * \param init       Operator to transform from input to output type.
 * \param null       Neuter element w.r.t. the operator op.
 */
template <typename T_in, typename T_out, typename T_op, typename T_init>
cudaError_t reduce_array(
        const T_in * __restrict__ data,
        T_out * __restrict__ out,
        const int n,
        const int block_size,
        const T_op &op,
        const T_init &init,
        const T_out null)
{
    Identity<T_out> id;
    int block_count = (n + block_size - 1) / block_size;

    T_out *d_data, *d_out;
    cuda_return_error(cudaMalloc(&d_data, block_count * sizeof (T_out)));
    cuda_return_error(cudaMalloc(&d_out, block_count * sizeof (T_out)));

    // Reduce volume to one value per block
    run_reduce_kernel<T_in, T_out, T_op, T_init>
                     (data, d_out, n, block_size, block_count, op, init, null);

    // Reduce the number of blocks 
    while (block_count > 1) {
        const int blocks = (block_count + block_size - 1) / block_size;
        cuda_return_error(cudaMemcpy(d_data, d_out, block_count * sizeof (T_out), cudaMemcpyDeviceToDevice));
        run_reduce_kernel<T_out, T_out, T_op, Identity<T_out>>
                         (d_data, d_out, block_count, block_size, blocks, op, id, null);
        block_count = blocks;
    }

    cuda_return_error(cudaMemcpy(out, d_out, sizeof (T_out), cudaMemcpyDeviceToHost));

    cuda_return_error(cudaFree(d_out));
    cuda_return_error(cudaFree(d_data));

    return cudaSuccess;
}

#endif /* ifndef REDUCTION_CUH_INCLUDED */
