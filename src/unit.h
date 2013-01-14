#ifndef _UNIT_H_
#define _UNIT_H_

struct elems {
	int val;
	struct elems *next;
};

typedef struct elems* TStack;
typedef struct elems* TQueue;

void stack_init(TStack* /*stack*/                  );
void stack_push(TStack* /*stack*/, int  /*val*/    );
int  stack_pop (TStack* /*stack*/, int* /*val_ptr*/);

void queue_init(TQueue* /*queue*/                  );
void queue_push(TQueue* /*queue*/, int  /*val*/    );
int  queue_pop (TQueue* /*queue*/, int* /*val_ptr*/);

#endif
