#include <mpi.h>
#include "queue.c"
#include <stdbool.h>

MPI_Datatype create_interval_type() {
    struct Interval interval;
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

    for (int i = 1; i < 6; i++) {
        displacements[i] -= displacements[0];
    }
    displacements[0] = 0;

    MPI_Type_create_struct(6, blocklengths, displacements, types, &interval_type);
    MPI_Type_commit(&interval_type);

    return interval_type;
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Datatype interval_type = create_interval_type();

    struct Queue queue;  // Define the queue
    init(&queue);  // Initialize the queue

    if (world_rank == 0) {
        struct Interval whole = {.left = 0.0, .right = 10.0, .tol = 1e-6, .f_left = func1(0.0), .f_right = func1(10.0), .f_mid = func1((0.0 + 10.0) / 2.0)};
        enqueue(whole, &queue);

        while (!isempty(&queue)) {
            int active_processes = 0;
            for (int i = 1; i < world_size && !isempty(&queue); i++) {
                struct Interval interval = dequeue(&queue);
                MPI_Send(&interval, 1, interval_type, i, 0, MPI_COMM_WORLD);
                active_processes++;
            }

            for (int i = 1; i <= active_processes; i++) {
                struct Interval results[2];
                MPI_Recv(results, 2, interval_type, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                enqueue(results[0], &queue);
                enqueue(results[1], &queue);
            }
        }
    } else {
        struct Interval local_interval;
        while (true) {  // Continue processing until there are no more intervals to process
            MPI_Status status;
            MPI_Recv(&local_interval, 1, interval_type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // Check if the master has sent a termination signal
            if (status.MPI_TAG == 1) break;

            double result;
            struct Interval sub_intervals[2];
            simpson(func1, local_interval, &sub_intervals[0], &sub_intervals[1], &result);
            MPI_Send(sub_intervals, 2, interval_type, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Type_free(&interval_type);  // Clean up the MPI data type
    MPI_Finalize();
    return 0;
}
