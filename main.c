#include <mpi.h>
#include "queue.c"


MPI_Datatype create_interval_type() {
    Interval interval;
    MPI_Datatype interval_type;
    int blocklengths[6] = {1, 1, 1, 1, 1, 1};
    MPI_Aint displacements[6];
    MPI_Datatype types[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}

    MPI_Get_address(&interval.left, &displacements[0])
    MPI_Get_address(&interval.right, &displacements[1])
    MPI_Get_address(&interval.left, &displacements[2])
    MPI_Get_address(&interval.tol, &displacements[3])
    MPI_Get_address(&interval.f_right, &displacements[4])
    MPI_Get_address(&interval.f_mid, &displacements[5])

    for (int i = 5; i > 0; i--) {
        displacements[i] -= displacements[0]
    }
    displacements[0] = 0

    MPI_Type_create_struct(6, blocklengths, displacements, types, &interval_type);
    MPI_Type_commit(&interval_type);

    return interval_type;

}

int main(void) {
    struct Queue queue;
    struct Interval whole;

    MPI_Init(&argc, &argv);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 
    if (rank == 0) {

        init(&queue);
        whole.left = 0.0;
        whole.right = 10.0;
        whole.tol = 1e-6;
        whole.f_left = func1(whole.left);
        whole.f_right = func1(whole.right);
        whole.f_mid = func1((whole.left + whole.right) / 2.0);

        enqueue(whole, &queue);

        while (!isempty(queue)) {

            int active_processes = 0;
            int tasks_sent = 0;
            struct Interval interval = dequeue(queue);
            for (i = 0; i < world_size; i++) {
                
            }

        }

    }
    // Perform integration
    double result = simpson(func1, &queue);
    printf("Result = %e\n", result);

    return 0;
}
