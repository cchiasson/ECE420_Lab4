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
    double *r, *r_pre, *sub_r, *sub_r_pre;
    int i, j;
    double damp_const;
    int iterationcount = 0;
    double error;
    double start, end, time = 0;

    int num_procs, my_id;

    MPI_Init(NULL, NULL);

    /* Get the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    /* Get my rank among all the processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    get_node_stat(&nodecount, &num_in_links, &num_out_links);
    sub_node = nodecount/num_procs; //278 for 4 processes, 556 for 2
    int sub_start = nodecount*my_id/num_procs;
    int sub_end = nodecount*(my_id+1)/num_procs;

    // initialize nodes to have same ranks
    r = malloc(nodecount * sizeof(double));
    sub_r = malloc(sub_node * sizeof(double)); //local r for each process
//    sub_r_pre = malloc(sub_node * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));
    for ( i = 0; i < nodecount; ++i) {
        r[i] = 1.0 / nodecount;
//        if (i<=sub_node) {
//            sub_r[i] = 1.0 / sub_node;
//        }
    }
    GET_TIME(start);
    // Calculate the result

    if (node_init(&nodehead, num_in_links, num_out_links,sub_start, sub_end)) return 254;

    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;
    // CORE CALCULATION
    do{
        ++iterationcount;
        // make a copy of the vector.
        vec_cp(r, r_pre, nodecount);
//        vec_cp(sub_r,sub_r_pre,sub_node);

        for ( i = 0; i < sub_node; ++i){
            sub_r[i] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                sub_r[i] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            sub_r[i] *= DAMPING_FACTOR;
            sub_r[i] += damp_const;
        }


        MPI_Scatter(sub_r, sub_node, MPI_DOUBLE, r, sub_node, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//        MPI_Gather(sub_r, sub_node, MPI_DOUBLE, r, sub_node, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    } while(rel_error(r, r_pre, nodecount) >= EPSILON);
    printf("Program converges at %d th iteration for rank %d\n", iterationcount, my_id);
    if (my_id==0) {
        MPI_Gather(r, sub_node, MPI_DOUBLE, sub_r, sub_node, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Gather(r, sub_node, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

//    if (my_id==0) {
//        MPI_Gather(sub_r, sub_node, MPI_DOUBLE, r, sub_node, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    }
//
//    else {
//        MPI_Gather(sub_r, sub_node, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    }
    GET_TIME(end);
    time = end-start;

    MPI_Finalize();

    if (my_id == 0) {
    Lab4_saveoutput(r, nodecount, time);
    }
    // post processing
    node_destroy(nodehead, nodecount);
    free(r);
    free(sub_r);
//    free(sub_r_pre);
    free(r_pre);/*
    free(num_in_links);
    free(num_out_links);*/
//    printf("done\n");

}

