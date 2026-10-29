/* CHUI opt-in MPI collective call-sequence tracer (diagnostic only).
 * Build: mpicc -std=c11 -O2 -fPIC -shared -o /tmp/chui_mpi_trace.so tools/diagnostics/mpi_collective_trace.c
 * Use:   COSMOSIM_TRACE_MPI=1 LD_PRELOAD=/tmp/chui_mpi_trace.so ctest --test-dir build/mpi-hdf5-fftw-debug -R '^integration_reference_workflow_distributed_treepm_mpi_two_rank$' --output-on-failure -j1
 * Prints each Allreduce/Bcast before entry, including world rank and communicator-local rank,
 * count and datatype size; can reveal rank-divergent MPI collective ordering and sizes.
 * Logs are intentionally verbose; never enable during production simulations.
 */
#include <mpi.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>

static atomic_ullong chui_collective_sequence = 0;

static int chui_trace_enabled(void) {
    const char *value = getenv("COSMOSIM_TRACE_MPI");
    return value && value[0] == '1' && value[1] == '\0';
}

static void chui_trace(const char *operation, MPI_Comm comm, int count, MPI_Datatype datatype) {
    if (!chui_trace_enabled()) return;
    int initialized = 0;
    PMPI_Initialized(&initialized);
    if (!initialized) return;
    int rank = -1, world_rank = -1, size = -1, typesize = -1;
    PMPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    PMPI_Comm_rank(comm, &rank);
    PMPI_Comm_size(comm, &size);
    if (datatype != MPI_DATATYPE_NULL) PMPI_Type_size(datatype, &typesize);
    const unsigned long long seq = atomic_fetch_add(&chui_collective_sequence, 1) + 1;
    fprintf(stderr, "CHUI_MPI_TRACE seq=%llu world_rank=%d comm_rank=%d comm_size=%d op=%s count=%d datatype_bytes=%d\n",
            seq, world_rank, rank, size, operation, count, typesize);
    fflush(stderr);
}

int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
    chui_trace("Allreduce", comm, count, datatype);
    return PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    chui_trace("Bcast", comm, count, datatype);
    return PMPI_Bcast(buffer, count, datatype, root, comm);
}

int MPI_Barrier(MPI_Comm comm) {
    chui_trace("Barrier", comm, 0, MPI_DATATYPE_NULL);
    return PMPI_Barrier(comm);
}

int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    chui_trace("Allgather", comm, sendcount, sendtype);
    return PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}

int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    chui_trace("Alltoall", comm, sendcount, sendtype);
    return PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}
