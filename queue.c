#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "function.h" // Ensure this header defines func1 correctly

#define MAXQUEUE 10000

struct Interval {
    double left;   // left boundary
    double right;  // right boundary
    double tol;    // tolerance
    double f_left; // function value at left boundary
    double f_mid;  // function value at midpoint
    double f_right; // function value at right boundary
};

struct Queue {
    struct Interval entry[MAXQUEUE]; // array of queue entries
    int top; // index of last entry
};

// Initialize queue
void init(struct Queue *queue_p) {
    queue_p->top = -1;
}

// Check if queue is empty
int isempty(struct Queue *queue_p) {
    return (queue_p->top == -1);
}

// Add an interval to the queue
void enqueue(struct Interval interval, struct Queue *queue_p) {
    if (queue_p->top == MAXQUEUE - 1) {
        printf("Maximum queue size exceeded - exiting\n");
        exit(1);
    }
    queue_p->top++;
    queue_p->entry[queue_p->top] = interval;
}

// Extract last interval from queue
struct Interval dequeue(struct Queue *queue_p) {
    if (queue_p->top == -1) {
        printf("Attempt to extract from empty queue - exiting\n");
        exit(1);
    }
    return queue_p->entry[queue_p->top--];
}

// Simpson's rule integration using a queue
double simpson(double (*func)(double), struct Queue *queue_p) {
    double quad = 0.0;

    while (!isempty(queue_p)) {
        struct Interval interval = dequeue(queue_p);
        double h = (interval.right - interval.left);
        double c = (interval.left + interval.right) / 2.0;
        double d = (interval.left + c) / 2.0;
        double e = (c + interval.right) / 2.0;
        double fd = func(d);
        double fe = func(e);

        double q1 = h / 6.0 * (interval.f_left + 4.0 * interval.f_mid + interval.f_right);
        double q2 = h / 12.0 * (interval.f_left + 4.0 * fd + 2.0 * interval.f_mid + 4.0 * fe + interval.f_right);

        if (fabs(q2 - q1) < interval.tol || h < 1.0e-12) {
            quad += q2 + (q2 - q1) / 15.0;
        } else {
            struct Interval i1 = {interval.left, c, interval.tol / 2, interval.f_left, fd, interval.f_mid};
            struct Interval i2 = {c, interval.right, interval.tol / 2, interval.f_mid, fe, interval.f_right};
            enqueue(i1, queue_p);
            enqueue(i2, queue_p);
        }
    }

    return quad;
}

