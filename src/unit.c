#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "unit.h"

void stack_init(TStack* stack) {
    (*stack) = NULL;
}

void stack_push(TStack* stack, int val) {
    TStack st = *stack;
    if ( NULL == st ) {     // Stack is empty
        st = (TStack)malloc(sizeof(struct elems));
        if ( NULL == st ) {
            fprintf(stderr, "error [stack init]: memory allocation error\n");
            exit(2);
        }
        st->val  = val;
        st->next = NULL;
        *stack   = st;
    } else {                // st is head
        st = (TStack)malloc(sizeof(struct elems));
        if ( NULL == st ) {
            fprintf(stderr, "error [stack push]: memory allocation error\n");
            exit(2);
        }
        st->next = *stack;
        st->val  = val;
        *stack   = st;
    }
}

int stack_pop (TStack* stack, int *val_ptr) {
    TStack st = *stack;
    if ( NULL == st ) return 0;
    *val_ptr = st->val;
    *stack   = st->next;

    free(st);
    return 1;
}

void queue_init(TQueue* queue) {
    (*queue) = NULL;
}

void queue_push(TQueue* queue, int val) {
    TQueue qu = (TQueue)malloc(sizeof(struct elems));
    if ( NULL == qu ) {
        fprintf(stderr, "error [queue push]: memory allocation error\n");
        exit(2);
    }
    qu->val  = val;
    qu->next = NULL;

    TQueue run = *queue;
    if ( NULL == run ) {
        *queue = qu;
        return;
    }

    while (run->next != NULL) run=run->next;
    run->next = qu;
}

int queue_pop (TQueue* queue, int *val_ptr) {
    TQueue qu = *queue;
    if ( NULL == qu ) return 0;
    *val_ptr = qu->val;
    *queue   = qu->next;

    free(qu);
    return 1;
}
