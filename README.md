This test record a nvshmem error.

In this example, a class is created. In side the class, recv\_queue and agg\_queue are allocated by nvshmem\_malloc. After the queues, a bunch of counters are allocated by nvshmem\_malloc. However, when the size of queue is beyond certain size, nvshmem\_uint32\_p fails to update counters on remote GPU.
