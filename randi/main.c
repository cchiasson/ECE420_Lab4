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


    // initialize nodes to have same ranks
    r = malloc(nodecount * sizeof(double));
    sub_r = malloc(sub_node * sizeof(double)); //local r for each process
    sub_r_pre = malloc(sub_node * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));
    for ( i = 0; i < nodecount; ++i) {
        r[i] = 1.0 / nodecount;
	sub_r[i]=1.0/nodecount;
    }
    GET_TIME(start);
    // Calculate the result
    
    printf("before if %d\n", my_id);


    if (node_init(&nodehead, num_in_links, num_out_links, 0, nodecount)) return 254;
    printf("before scatter %d\n", my_id);
    MPI_Scatter(r, sub_node, MPI_DOUBLE, sub_r, sub_node, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    printf("after scatter %d\n", my_id);


    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;
    // CORE CALCULATION
    do{
        ++iterationcount;
        // make a copy of the vector.
        //vec_cp(r, r_pre, nodecount);
	vec_cp(sub_r,sub_r_pre,sub_node);

        for ( i = 0; i < sub_node; ++i){
            sub_r[i] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                sub_r[i] += sub_r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            sub_r[i] *= DAMPING_FACTOR;
            sub_r[i] += damp_const;

	    /*MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	    if (my_id==1) {
		printf("%f\n",sub_r[i]);
	    }*/        
	}

	//printf("Before gather\n");
	//MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	//printf("ID: %d\n",my_id);
        /*if (my_id==0) {
                printf("attempting gather for 0\n");
        	MPI_Gather(sub_r, sub_node, MPI_DOUBLE, r, sub_node, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	else {
        	MPI_Gather(sub_r, sub_node, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}*/
	//printf("after gather\n");
    } while(rel_error(sub_r, sub_r_pre, sub_node) >= EPSILON);
    	printf("Program converges at %d th iteration for rank %d\n", iterationcount, my_id);

    
    if (my_id==0) {
	printf("attempting gather for 0\n");
       	MPI_Gather(sub_r, sub_node, MPI_DOUBLE, r, sub_node, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    else {
	printf("attempting gather for others\n");
       	MPI_Gather(sub_r, sub_node, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    GET_TIME(end);
    time = end-start;
    
	printf("attempting finalize\n");
    MPI_Finalize();
    
	printf("after finalize\n");


    //if (my_id == 0) {
        Lab4_saveoutput(r, nodecount, time);
    //}
    // post processing
    printf("destroying\n");
    node_destroy(nodehead, nodecount);
    printf("destroying1\n");/*
    free(r);
    printf("destroying2\n");
    free(sub_r);
    printf("destroying3\n");
    free(sub_r_pre);
    printf("destroying4\n");
    free(r_pre);/*
    printf("destroying5\n");
    free(num_in_links);
    printf("destroying6\n");
    free(num_out_links);*/
    printf("done\n");
    
}
