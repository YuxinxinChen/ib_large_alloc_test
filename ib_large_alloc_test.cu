#include<stdio.h>
#include <iostream>
#include <assert.h>
#include <unistd.h>
#include "mpi.h"
#include "nvshmem.h"
#include "nvshmemx.h"

#include "../util/error_util.cuh"
#include "../util/nvshmem_util.cuh"

namespace dev {
	template<typename LOCAL_T, typename RECV_T, typename COUNTER_T, int PADDING_SIZE>
	class Queues {
	public:
		int n_pes;
		int my_pe;
		int group_size;
	
		COUNTER_T recv_capacity;
		COUNTER_T local_capacity;
			
		LOCAL_T *local_queues;
		RECV_T *recv_queues;
		RECV_T *agg_queues;	

		int num_agg_queues;
		int num_local_queues;
		int total_num_queues;

		volatile COUNTER_T *start, *start_alloc, *end, *end_alloc, *end_max, *end_count;
	    volatile int *stop;
		Queues(int _n_pes, int _my_pe, int _group_size, LOCAL_T *l_q, RECV_T *r_q, RECV_T *a_q, COUNTER_T l_capacity, COUNTER_T r_capacity,
    	       COUNTER_T *_start, COUNTER_T *_start_alloc, COUNTER_T *_end_alloc, COUNTER_T *_end, COUNTER_T *_end_max, COUNTER_T *_end_count, int *_stop, int _num_l_q,
    	       int _num_a_q): n_pes(_n_pes), my_pe(_my_pe), group_size(_group_size), local_queues(l_q), recv_queues(r_q), agg_queues(a_q), 
    	       num_local_queues(_num_l_q), num_agg_queues(_num_a_q), local_capacity(l_capacity), recv_capacity(r_capacity) 
    	{
    	    start = (volatile COUNTER_T *)_start;
    	    start_alloc = (volatile COUNTER_T *)_start_alloc;
    	    end_alloc = (volatile COUNTER_T *)_end_alloc;
    	    end = (volatile COUNTER_T *)_end;
    	    end_max = (volatile COUNTER_T *)_end_max;
    	    end_count = (volatile COUNTER_T *)_end_count;
    	    stop = (volatile int *)_stop;

    	    total_num_queues=num_local_queues+n_pes-1+num_agg_queues;
    	}
	};
}

template<typename LOCAL_T, typename RECV_T, typename COUNTER_T, int PADDING_SIZE>
class Queues {
public:
	int n_pes;
	int my_pe;
	int nodes_size;
	int node_id;
	int group_size;
	int group_id;

	COUNTER_T recv_capacity;
	COUNTER_T local_capacity;
		
	LOCAL_T *local_queues;
	RECV_T *recv_queues;
	RECV_T *agg_queues;	

	int num_agg_queues;
	int num_local_queues;
	int total_num_queues;

	COUNTER_T *counters;
	int num_counters=7;
	COUNTER_T *start, *start_alloc, *end, *end_alloc, *end_max, *end_count;
    int *stop;

	Queues() {}
    ~Queues() { release(); }

    __host__ void baseInit(int _n_pes, int _my_pe, int _group_id, int _group_size, int local_id, int local_size,
    COUNTER_T l_capacity, COUNTER_T r_capacity, int l_queues = 1, bool PRINT_INFO = false)
    {
        n_pes = _n_pes;
        my_pe = _my_pe;
        nodes_size= _group_size;
        node_id = _group_id;
        group_size = local_size;
        group_id = local_id;

        local_capacity = l_capacity;
        recv_capacity = r_capacity;

        num_local_queues = l_queues;
        num_agg_queues = (nodes_size-1)*group_size;
        total_num_queues = num_local_queues+n_pes-1+num_agg_queues;

        alloc(PRINT_INFO);
    }

	dev::Queues<LOCAL_T, RECV_T, COUNTER_T, PADDING_SIZE>
    deviceObject() const {
    	return dev::Queues<LOCAL_T, RECV_T, COUNTER_T, PADDING_SIZE>
        (n_pes, my_pe, group_size, local_queues,  recv_queues, agg_queues, local_capacity, recv_capacity,
        start, start_alloc, end_alloc, end, end_max, end_count, stop, num_local_queues, num_agg_queues);
    }

private:
	void alloc(bool PRINT_INFO = false)
    {
        if(local_capacity+recv_capacity <= 0) return;
        if(PRINT_INFO)
                std::cout << "pe "<< my_pe << " called distributed queue base allocator\n";
        CUDA_CHECK(cudaMalloc(&local_queues, sizeof(LOCAL_T)*local_capacity*num_local_queues));
        CUDA_CHECK(cudaMemset(local_queues, 0xffffffff, sizeof(LOCAL_T)*local_capacity*num_local_queues));

        recv_queues = (RECV_T *)nvshmem_malloc(sizeof(RECV_T)*recv_capacity*(n_pes-1));
        CUDA_CHECK(cudaMemset(recv_queues, 0xffffffff, sizeof(RECV_T)*recv_capacity*(n_pes-1)));

        agg_queues = (RECV_T *)nvshmem_malloc(sizeof(RECV_T)*recv_capacity*num_agg_queues);
        CUDA_CHECK(cudaMemset(agg_queues, 0xffffffff, sizeof(RECV_T)*recv_capacity*num_agg_queues));

        counters = (COUNTER_T *)nvshmem_malloc(sizeof(COUNTER_T)*num_counters*PADDING_SIZE*total_num_queues);
        start = counters;
        start_alloc = (counters+1*PADDING_SIZE*total_num_queues);
        end_alloc = (counters+2*PADDING_SIZE*total_num_queues);
        end = (counters+3*PADDING_SIZE*total_num_queues);
        end_max = (counters+4*PADDING_SIZE*total_num_queues);
        end_count = (counters+5*PADDING_SIZE*total_num_queues);

        stop = (int *)(counters+6*PADDING_SIZE*total_num_queues);
    }

	void release(bool PRINT_INFO = false) {
        if(local_capacity+recv_capacity <= 0) return;
        if(PRINT_INFO)
            std::cout << "pe "<< my_pe << " call distributed queue base destructor\n";
        if(local_queues!=NULL)
        CUDA_CHECK(cudaFree(local_queues));

		// if uncommend following lines, got errors: src/util/cs.cpp:26: non-zero status: 22: No such file or directory, exiting... mutex lock failed
        //nvshmem_free(agg_queues);
        //nvshmem_free(recv_queues);
        //nvshmem_free(counters);
    }
};

template<typename T, typename Y>
struct Pack {
	T id;
	Y value;
	__device__ __host__ Pack() {}
	__device__ __host__ Pack(T _id, Y _value): id(_id), value(_value) {}
};

template<typename LOCAL_T, typename RECV_T, typename COUNTER_T, int PADDING_SIZE>
__global__ void set(dev::Queues<LOCAL_T, RECV_T, COUNTER_T, PADDING_SIZE> queue) {
	if(threadIdx.x == 0)
	nvshmem_uint32_p((uint32_t *)(queue.end+queue.num_local_queues*PADDING_SIZE), 12314, (queue.my_pe^1));
}

int main(int argc, char** argv)
{

    int n_pes, my_pe, group_id, group_size, local_id, local_size; 
    nvshm_mpi_init(my_pe, n_pes, group_id, group_size, local_id, local_size, &argc, &argv);

	Queues<int, Pack<int, float>, uint32_t, 32> queue;
	queue.baseInit(n_pes, my_pe, group_id, group_size, local_id, local_size, 1024, 300000000);
	
	cudaStream_t stream;
    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
	nvshmem_barrier_all();
	
	set<<<1,32, 0, stream>>>(queue.deviceObject());

	CUDA_CHECK(cudaStreamSynchronize(stream));
	nvshmem_barrier_all();

	uint32_t end;
	CUDA_CHECK(cudaMemcpy(&end, queue.end+32, sizeof(uint32_t), cudaMemcpyDeviceToHost));
	printf("PE %d, end %d\n", my_pe, end);

	nvshmem_barrier_all();
    nvshm_mpi_finalize();
	return 0;
}
