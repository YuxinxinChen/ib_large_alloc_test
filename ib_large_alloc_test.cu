#include<stdio.h>
#include <iostream>
#include <assert.h>
#include <unistd.h>
#include "mpi.h"
#include "nvshmem.h"
#include "nvshmemx.h"

#include "../util/error_util.cuh"
#include "../util/nvshmem_util.cuh"

template<typename RECV_T, typename COUNTER_T, int PADDING_SIZE>
class Queues {
public:
	int n_pes;
	int my_pe;
	int nodes_size;
	int node_id;
	int group_size;
	int group_id;

	COUNTER_T recv_capacity;
		
	RECV_T *recv_queues;

	COUNTER_T *counters;
	int num_counters=4;
	COUNTER_T *start, *start_alloc, *end, *end_alloc;

	Queues() {}
    ~Queues() {}

    __host__ void baseInit(int _n_pes, int _my_pe, int _group_id, int _group_size, int local_id, int local_size,
    COUNTER_T r_capacity, bool PRINT_INFO = false)
    {
        n_pes = _n_pes;
        my_pe = _my_pe;
        nodes_size= _group_size;
        node_id = _group_id;
        group_size = local_size;
        group_id = local_id;

        recv_capacity = r_capacity;

        alloc(PRINT_INFO);
    }

private:
	void alloc(bool PRINT_INFO = false)
    {
        if(recv_capacity <= 0) return;
        if(PRINT_INFO)
                std::cout << "pe "<< my_pe << " called distributed queue base allocator\n";

        recv_queues = (RECV_T *)nvshmem_malloc(sizeof(RECV_T)*recv_capacity*n_pes);
        CUDA_CHECK(cudaMemset(recv_queues, 0xffffffff, sizeof(RECV_T)*recv_capacity*n_pes));

        counters = (COUNTER_T *)nvshmem_malloc(sizeof(COUNTER_T)*num_counters*PADDING_SIZE*n_pes);
        start = counters;
        start_alloc = (counters+1*PADDING_SIZE*n_pes);
        end_alloc = (counters+2*PADDING_SIZE*n_pes);
        end = (counters+3*PADDING_SIZE*n_pes);
    }

	void release(bool PRINT_INFO = false) {
        if(recv_capacity <= 0) return;
        if(PRINT_INFO)
            std::cout << "pe "<< my_pe << " call distributed queue base destructor\n";

		// if uncommend following lines, got errors: src/util/cs.cpp:26: non-zero status: 22: No such file or directory, exiting... mutex lock failed
        //nvshmem_free(agg_queues);
        //nvshmem_free(recv_queues);
        //nvshmem_free(counters);
    }
};

template<typename RECV_T, typename COUNTER_T, int PADDING_SIZE>
__global__ void set(Queues<RECV_T, COUNTER_T, PADDING_SIZE> queue) {
	if(threadIdx.x == 0)
	nvshmem_uint32_p((uint32_t *)(queue.end), 12314, (queue.my_pe^1));
}

int main(int argc, char** argv)
{
	uint32_t recv_size = 600000000;
	if(argc > 1)
        for(int i=1; i<argc; i++) {
            if(std::string(argv[i]) == "-size")
                recv_size = std::stoi(argv[i+1]);
        }

	printf("recv_size %lld\n", size_t(recv_size));

    int n_pes, my_pe, group_id, group_size, local_id, local_size; 
    nvshm_mpi_init(my_pe, n_pes, group_id, group_size, local_id, local_size, &argc, &argv);

	Queues<uint32_t, uint32_t, 32> queue;
	queue.baseInit(n_pes, my_pe, group_id, group_size, local_id, local_size, recv_size);
	
	cudaStream_t stream;
    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
	nvshmem_barrier_all();
	
	set<<<1,32, 0, stream>>>(queue);

	CUDA_CHECK(cudaStreamSynchronize(stream));
	nvshmem_barrier_all();

	uint32_t end;
	CUDA_CHECK(cudaMemcpy(&end, queue.end, sizeof(uint32_t), cudaMemcpyDeviceToHost));
	printf("PE %d, end %d\n", my_pe, end);

	nvshmem_barrier_all();
    nvshm_mpi_finalize();
	return 0;
}
