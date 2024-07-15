#include <mpi.h>
#include "queue.h"
#include "function.h"
#include <stdbool.h>
#include <mpi.h>
#include <stdio.h>


MPI_Datatype create_interval_type() {
    MPI_Datatype interval_type;
    int blocklengths[] = {1, 1, 1, 1, 1, 1};
    MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[6];

    offsets[0] = offsetof(struct Interval, left);
    offsets[1] = offsetof(struct Interval, right);
    offsets[2] = offsetof(struct Interval, tol);
    offsets[3] = offsetof(struct Interval, f_left);
    offsets[4] = offsetof(struct Interval, f_right);
    offsets[5] = offsetof(struct Interval, f_mid);

    MPI_Type_create_struct(6, blocklengths, offsets, types, &interval_type);
    MPI_Type_commit(&interval_type);

    return interval_type;
}

int main(int argc, char **argv) {
    int mpi_result = MPI_Init(&argc, &argv);
    if (mpi_result != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Init failed\n");
        return 1;
    }

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Datatype interval_type = create_interval_type();

    if (world_rank == 0) {
        struct Queue queue;
        init(&queue);  // Initialize the queue

        // Define the whole interval to process
        struct Interval whole = {.left = 0.0, .right = 10.0, .tol = 1e-6, .f_left = func1(0.0), .f_right = func1(10.0), .f_mid = func1(5.0)};
        enqueue(whole, &queue);

        double total_integral = 0.0;

        while (!isempty(&queue)) {
            int active_processes = 0;
            for (int i = 1; i < world_size && !isempty(&queue); i++) {
                struct Interval interval = dequeue(&queue);
                mpi_result = MPI_Send(&interval, 1, interval_type, i, 0, MPI_COMM_WORLD);
                if (mpi_result != MPI_SUCCESS) {
                    fprintf(stderr, "MPI_Send failed\n");
                    MPI_Abort(MPI_COMM_WORLD, mpi_result);
                }
                active_processes++;
            }

            for (int i = 1; i <= active_processes; i++) {
                struct Interval results[2];
                double partial_result;
                mpi_result = MPI_Recv(&partial_result, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (mpi_result != MPI_SUCCESS) {
                    fprintf(stderr, "MPI_Recv failed\n");
                    MPI_Abort(MPI_COMM_WORLD, mpi_result);
                }
                
                if (partial_result != 0.0) {
                    total_integral += partial_result;
                } else {
                    mpi_result = MPI_Recv(results, 2, interval_type, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (mpi_result != MPI_SUCCESS) {
                        fprintf(stderr, "MPI_Recv failed\n");
                        MPI_Abort(MPI_COMM_WORLD, mpi_result);
                    }
                    enqueue(results[0], &queue);
                    enqueue(results[1], &queue);
                }
            }
        }

        // Send termination message to all workers
        for (int i = 1; i < world_size; i++) {
            struct Interval terminate = {.left = 0, .right = 0, .tol = 0};
            mpi_result = MPI_Send(&terminate, 1, interval_type, i, 1, MPI_COMM_WORLD);
            if (mpi_result != MPI_SUCCESS) {
                fprintf(stderr, "MPI_Send failed\n");
                MPI_Abort(MPI_COMM_WORLD, mpi_result);
            }
        }

        printf("Integration complete. Total integral: %f\n", total_integral);
    } else {
        struct Interval local_interval;
        double partial_integral = 0.0;
        while (true) {
            MPI_Status status;
            mpi_result = MPI_Recv(&local_interval, 1, interval_type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (mpi_result != MPI_SUCCESS) {
                fprintf(stderr, "MPI_Recv failed\n");
                MPI_Abort(MPI_COMM_WORLD, mpi_result);
            }
            
            if (status.MPI_TAG == 1) {  // Termination message
                break;
            }

            double result;
            struct Interval sub_intervals[2];
            simpson(func1, local_interval, &sub_intervals[0], &sub_intervals[1], &result);

            if (result != 0) {
                partial_integral += result;
                mpi_result = MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                if (mpi_result != MPI_SUCCESS) {
                    fprintf(stderr, "MPI_Send failed\n");
                    MPI_Abort(MPI_COMM_WORLD, mpi_result);
                }
            } else {
                double zero = 0.0;
                mpi_result = MPI_Send(&zero, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                if (mpi_result != MPI_SUCCESS) {
                    fprintf(stderr, "MPI_Send failed\n");
                    MPI_Abort(MPI_COMM_WORLD, mpi_result);
                }
                mpi_result = MPI_Send(sub_intervals, 2, interval_type, 0, 0, MPI_COMM_WORLD);
                if (mpi_result != MPI_SUCCESS) {
                    fprintf(stderr, "MPI_Send failed\n");
                    MPI_Abort(MPI_COMM_WORLD, mpi_result);
                }
            }
        }

        // Send partial result back to master
        mpi_result = MPI_Send(&partial_integral, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        if (mpi_result != MPI_SUCCESS) {
            fprintf(stderr, "MPI_Send failed\n");
            MPI_Abort(MPI_COMM_WORLD, mpi_result);
        }
    }

    MPI_Type_free(&interval_type);
    MPI_Finalize();
    return 0;
}