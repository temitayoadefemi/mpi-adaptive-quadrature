#include <mpi.h>
#include "queue.c"

MPI_Datatype create_interval_type() {
    Interval interval;
    MPI_Datatype interval_type;
    int blocklengths[6] = {1, 1, 1, 1, 1, 1};
    MPI_Aint displacements[6];
    MPI_Datatype types[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

    MPI_Get_address(&interval.left, &displacements[0]);
    MPI_Get_address(&interval.right, &displacements[1]);
    MPI_Get_address(&interval.tol, &displacements[2]);
    MPI_Get_address(&interval.f_left, &displacements[3]);
    MPI_Get_address(&interval.f_right, &displacements[4]);
    MPI_Get_address(&interval.f_mid, &displacements[5]);

    for (int i = 5; i > 0; i--) {
        displacements[i] -= displacements[0];
    }
    displacements[0] = 0;

    MPI_Type_create_struct(6, blocklengths, displacements, types, &interval_type);
    MPI_Type_commit(&interval_type);

    return interval_type;
}

int main(int argc, char **argv) {
    struct Queue queue;
    struct Interval whole;

    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Datatype interval_type = create_interval_type();

    if (world_rank == 0) {  // Master process
        init(&queue);
        whole.left = 0.0;
        whole.right = 10.0;
        whole.tol = 1e-6;
        whole.f_left = func1(whole.left);
        whole.f_right = func1(whole.right);
        whole.f_mid = func1((whole.left + whole.right) / 2.0);

        // Divide the interval among available workers
        for (int i = 1; i < world_size; i++) {
            MPI_Send(&whole, 1, interval_type, i, 0, MPI_COMM_WORLD);  // Send the same or different data
        }

        // Collect results from workers
        double result = 0.0, worker_result;
        for (int i = 1; i < world_size; i++) {
            MPI_Recv(&worker_result, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            result += worker_result;  // Accumulate results
        }

        printf("Result = %e\n", result);
    } else {  // Worker processes
        Interval local_interval;
        MPI_Recv(&local_interval, 1, interval_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Perform local computation, e.g., a local integration
        double local_result = simpson(func1, local_interval.left, local_interval.right, local_interval.tol);

        // Send the result back to the master
        MPI_Send(&local_result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Type_free(&interval_type);
    MPI_Finalize();
    return 0;
}
