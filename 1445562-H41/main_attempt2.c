#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab4_IO.h"
#include "timer.h"
#include "mpi.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

int main (int argc, char* argv[]){
    struct node *nodehead;
    int nodecount, sub_node;
    int *num_in_links, *num_out_links;
    double *r, *r_pre, *sub_r;
    int i, j;
    double damp_const;
    double cst_addapted_threshold;
    double error;
    double start, end, time = 0;
    int num_procs, my_id;

    MPI_Init(NULL, NULL);

    /* Get the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    /* Get my rank among all the processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    get_node_stat(&nodecount, &num_in_links, &num_out_links);

    int sub_start = nodecount*my_id/num_procs;
    int sub_end = nodecount*(my_id+1)/num_procs;
    sub_node = nodecount/num_procs;

    // initialize nodes to have same ranks
    r = malloc(nodecount * sizeof(double));
    sub_r = malloc(sub_node * sizeof(double)); //local r for each process

    for (i = 0; i < nodecount; ++i) {
        r[i] = 1.0 / nodecount;
    }
    GET_TIME(start);

    // initialize the node
    if (node_init(&nodehead, num_in_links, num_out_links, sub_start, sub_end)) return 254;

    //printf("%d\n", nodecount);
    r_pre = malloc(nodecount * sizeof(double));

    int iterationcount = 0;
    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;
    // CORE CALCULATION
    do{
        ++iterationcount;
        // make a copy of the vector.
        vec_cp(r, r_pre, nodecount);

        for ( i = 0; i < sub_node; ++i){
            sub_r[i] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j) {
                sub_r[i] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            }
            sub_r[i] *= DAMPING_FACTOR;
            sub_r[i] += damp_const;
        }

        MPI_Allgather(sub_r, sub_node, MPI_DOUBLE, r, sub_node, MPI_DOUBLE, MPI_COMM_WORLD);
    } while(rel_error(r, r_pre, nodecount) >= EPSILON);
    	printf("Program converges at %d th iteration.\n", iterationcount);

    GET_TIME(end);
    time = end-start;

    MPI_Finalize();

    if (my_id == 0) {
        Lab4_saveoutput(r, nodecount, time);
    }

    // post processing
    node_destroy(nodehead, sub_node);

    // free memory
    free(r);
    free(sub_r);
    free(r_pre);

}
