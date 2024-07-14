#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "function.h" // Ensure this header defines func1 correctly
#include "queue.h"

void init(struct Queue *queue_p) {
    queue_p->top = -1;
}

int isempty(struct Queue *queue_p) {
    return queue_p->top == -1;
}

int size(struct Queue *queue_p) {
    return queue_p->top + 1;
}

void enqueue(struct Interval interval, struct Queue *queue_p) {
    if (queue_p->top == MAXQUEUE - 1) {
        printf("Maximum queue size exceeded - exiting\n");
        exit(1);
    }
    queue_p->top++;
    queue_p->entry[queue_p->top] = interval;
}

struct Interval dequeue(struct Queue *queue_p) {
    if (queue_p->top == -1) {
        printf("Attempt to extract from empty queue - exiting\n");
        exit(1);
    }
    return queue_p->entry[queue_p->top--];
}
// Simpson's rule integration using a queue
void simpson(double (*func)(double), struct Interval interval, struct Interval *i1, struct Interval *i2, double *result) {
    double quad = 0.0;
    double h = (interval.right - interval.left);
    double c = (interval.left + interval.right) / 2.0;
    double d = (interval.left + c) / 2.0;
    double e = (c + interval.right) / 2.0;
    double fd = func(d);
    double fe = func(e);

    double q1 = h / 6.0 * (interval.f_left + 4.0 * interval.f_mid + interval.f_right);
    double q2 = h / 12.0 * (interval.f_left + 4.0 * fd + 2.0 * interval.f_mid + 4.0 * fe + interval.f_right);

    if (fabs(q2 - q1) < interval.tol || h < 1.0e-12) {
        quad = q2 + (q2 - q1) / 15.0;
        *result = quad;  // Store the result via pointer
    } else {
        // Populate the interval pointers for further processing
        *i1 = (struct Interval){interval.left, c, interval.tol / 2, interval.f_left, fd, interval.f_mid};
        *i2 = (struct Interval){c, interval.right, interval.tol / 2, interval.f_mid, fe, interval.f_right};
        *result = 0;  // Indicative that further processing is needed
    }
}
