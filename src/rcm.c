#include "rcm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "unit.h"

#define ERROR_MESSAGE(c,e) do { printf("%s (err_code: %d)\n", c, e); return e; } while(0)

int find_periphery(TWGraph *gr, int *per) {
	int i;
	int level;
	int max_level, max_vertex, min_neighbour;
	int glob_max_level, glob_max_vertex;
	int viewed;
	int xadj_start, xadj_end; 
	int *vertex_level = (int*)malloc(sizeof(int)*(gr->size));

	if ( NULL == vertex_level ) {
		fprintf(stderr, "error [find_periphery]: memory allocation error\n");
		return ERROR_MEMORY_ALLOCATION;
	}

	int g, s;
	int *xadj, *adjncy;
	xadj   = gr->xadj;
	adjncy = gr->adjncy;

	TQueue queue;
	queue_init(&queue);

	max_level  = 1;
	max_vertex = 0;

	do {
		memset(vertex_level, 0, sizeof(int)*(gr->size));
//		for (i = 0; i < gr->size; i++) vertex_level[i] = 0;

		glob_max_level  = max_level;
		glob_max_vertex = max_vertex;

		// Set level for first vertex equals 1
		queue_push(&queue, glob_max_vertex);
		vertex_level[glob_max_vertex] = 1;

		// Set maxmin variables
		max_level     = 1;
		max_vertex    = glob_max_vertex;
		min_neighbour = gr->size;
		viewed = 1;

		while (queue_pop(&queue, &g)) {
			level = vertex_level[g];				// parent level

			xadj_start = xadj[g];
			xadj_end   = xadj[g+1];

			for (i = xadj_start; i < xadj_end; i++) {
				s = adjncy[i];
				if (!vertex_level[s]) {
					viewed++;						// one more vertex is viewed now
					vertex_level[s] = level + 1;	// increasing level
					queue_push(&queue, s);			// will 
				}
			}

			if (viewed == gr->size) {				// smart thing to decrease number of operations
				if (level > max_level) {			// search for vertex with max level
					max_level     = level;
					max_vertex    = g;
					min_neighbour = xadj_end- xadj_start;
				}
													// but with minimal number of neighbours
				if ( (level == max_level) && (xadj_end - xadj_start < min_neighbour) ) {
					max_vertex = g;
					min_neighbour = xadj_end - xadj_start;
				}
			}
		}
#ifdef _DEBUG_LEVEL_1
		printf("[debug(1)]{find_periphery} current max vertex: %d [level:%d]\n", max_vertex, max_level);
#endif
	} while (max_level > glob_max_level);			// repeat until glob_max_level isn't changed

#ifdef _DEBUG_LEVEL_1
	printf("[debug(1)]{find_periphery} glob_max_level: %d, periphery: %d\n", glob_max_level, glob_max_vertex);
#endif

	free(vertex_level);

	*per = glob_max_vertex;
	return ERROR_NO_ERROR;
}

int find_permutation(TWGraph *gr, int root, int** _perm, int** _invp) {
	// Making permututation vector
	//
	// v-th vertex comes perm[v] 
	// e.g. perm[0] <-- root
	// 		invp[root] <-- 0
	int v, g, i, s;
	TQueue queue;
	queue_init(&queue);

	int *perm = (int*)malloc(sizeof(int)*(gr->size));	// permutation vector
	int *invp = (int*)malloc(sizeof(int)*(gr->size));	// but also we have to use inverted permutation

	if ( NULL == perm || NULL == invp ) {
		if (perm) free(perm);
		if (invp) free(invp);
		fprintf(stderr, "error [find_permutation]: memory allocation error\n");
		return ERROR_MEMORY_ALLOCATION;
	}

	//	-2   -- not yet seen
	//	-1   -- already in queue
	//	0..n -- already numbered 
	for (i = 0; i < gr->size; i++) perm[i] = invp[i] = -2;

#ifdef _DEBUG_LEVEL_1
printf("[debug(1)]{find_permutation}: root = %d\n", root);
#endif

	// actually it is simple bs
	invp[root] = -1;
	queue_push(&queue, root);

	s = 0;
	while (queue_pop(&queue, &v)) {

		for (i = gr->xadj[v]; i < gr->xadj[v+1]; i++) {
			g = gr->adjncy[i];
			if ( invp[g] < 0 && FABS(gr->wvert[g]) < EPS_THRESHOLD ) {
				stack_push(&queue, v);
				v = g;
				break;
			}
		}

#ifdef DIRECT_CM
		perm[s] = v;
		invp[v] = s++;				// put connected vertices closer
#else
		perm[gr->size - 1 - s] = v;
		invp[v] = gr->size - 1 - s++;		// put connected vertices closer
#endif
		for (i = gr->xadj[v]; i < gr->xadj[v+1]; i++) {
			g = gr->adjncy[i];
			if ( invp[g] == -2 ) {
				invp[g] = -1;
				queue_push(&queue, g);	// put unindexed vertices to queue
			}
		}
	}

	*_perm = perm;
	*_invp = invp;

#ifdef _DEBUG_LEVEL_1
int __i;
printf("[debug(1)]{find_permutation}:");
printf("\n   i:");
for (__i = 0; __i < gr->size; __i++) printf("%4d ", __i);
printf("\nperm:");
for (__i = 0; __i < gr->size; __i++) printf("%4d ", perm[__i]);
printf("\ninvp:");
for (__i = 0; __i < gr->size; __i++) printf("%4d ", invp[__i]);
printf("\n");
#endif

	return ERROR_NO_ERROR;
}

int graph_reorder(TWGraph *gr, int *perm, int *invp) {
	real *wvert, *wedge;
	int  *adjncy, *xadj;
	int  i, j, ci, g;

#ifdef _DEBUG_LEVEL_1
printf("[debug(1)]{graph_reorder} INPUT graph\n");
graph_show(gr);
#endif

	wedge  = (real*)malloc(sizeof(real)*(gr->nonz));
	adjncy = ( int*)malloc(sizeof( int)*(gr->nonz));
	xadj   = ( int*)malloc(sizeof( int)*(gr->size+1));
	wvert  = (real*)malloc(sizeof(real)*(gr->size));

	if ( NULL == wedge || NULL == adjncy || NULL == xadj || NULL == wvert) {
		if (  wedge ) free( wedge);
		if ( adjncy ) free(adjncy);
		if (   xadj ) free(  xadj);
		if (  wvert ) free( wvert);

		fprintf(stderr, "error [graph_reorder]: memory allocation error\n");
		return ERROR_MEMORY_ALLOCATION;
	}	

	// let's reorder
	ci = 0;
	for (i = 0; i < gr->size; i++) {					// run over all rows
		wvert[i] = gr->wvert[perm[i]];					// copy vertices weight
		xadj[i] = ci;
		for (j=gr->xadj[perm[i]]; j<gr->xadj[perm[i]+1]; j++) {		// run over all neighbours
			g = gr->adjncy[j];
			adjncy[ci] = invp[g];					// make edge
			wedge[ci]  = gr->wedge[j];				// set  edge weight
			ci++;
		}
	}
	xadj[i] = ci;								// xadj[size] = nonz

#ifdef _DEBUG_LEVEL_1
printf("[debug(1)]{graph_reorder} REORDERED graph\n");
graph_show(gr);
#endif

	free(gr->wedge);
	free(gr->adjncy);
	free(gr->xadj);
	free(gr->wvert);

	gr->wedge  = wedge;
	gr->adjncy = adjncy;
	gr->xadj   = xadj;
	gr->wvert  = wvert;

	return ERROR_NO_ERROR;
}

//	REORDERING: RCM. Reverse Cuthill-McKey
int rcm(TMatrix_DCSR *matr) {
	TWGraph gr;
	int *perm, *invp;
	int err = 0;

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rcm}: graph_builder\n");
#endif
	err = build_graph(&gr, matr);
	if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rcm: graph_builder failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rcm} find_periphery\n");
#endif
	int root;
	err = find_periphery(&gr, &root);	 				// RCM starts here (root vertex)
	if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rcm: find_periphery failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rcm} find_permutation\n");
#endif
	err = find_permutation(&gr, root, &perm, &invp);	// <-- reordering
	if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rcm: find_permutation failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rcm} graph_reoder\n");
#endif
	err = graph_reorder(&gr, perm, invp);
	if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rcm: graph_reorder failed", err);

	free(perm);
	free(invp);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rcm} matrix_builder\n");
#endif
	err = build_matrix(&gr, matr, 0);
	if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rcm: matrix_builder failed", err);

	free(gr.wedge);
	free(gr.adjncy);
	free(gr.xadj);
	free(gr.wvert);

	return ERROR_NO_ERROR;
}

