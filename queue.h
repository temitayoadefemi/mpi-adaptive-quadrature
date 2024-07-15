#ifndef QUEUE_H
#define QUEUE_H

#define MAXQUEUE 10000

struct Interval {
    double left;
    double right;
    double tol;
    double f_left;
    double f_mid;
    double f_right;
};

struct Queue {
    struct Interval entry[MAXQUEUE];
    int top;
};

void init(struct Queue *queue_p);
void enqueue(struct Interval interval, struct Queue *queue_p);
struct Interval dequeue(struct Queue *queue_p);
int isempty(struct Queue *queue_p);
int size(struct Queue *queue_p);
void simpson(double (*func)(double), struct Interval interval, struct Interval *i1, struct Interval *i2, double *result);

#endif // QUEUE_H
