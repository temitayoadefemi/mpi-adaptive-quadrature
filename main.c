#include <mpi.h>
#include "queue.c"

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

    if (world_rank == 0) {

        init(&queue);
        whole.left = 0.0;
        whole.right = 10.0;
        whole.tol = 1e-6;
        whole.f_left = func1(whole.left);
        whole.f_right = func1(whole.right);
        whole.f_mid = func1((whole.left + whole.right) / 2.0);

        enqueue(whole, &queue);

        }

    // Perform integration
    double result = simpson(func1, &queue);
    printf("Result = %e\n", result);

    return 0;
}
