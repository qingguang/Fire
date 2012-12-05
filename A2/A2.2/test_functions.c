/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "metis.h"
#include "util_read_files.h"
//#include "initialization.h"
#include "mpi.h"
#include "util_write_files.h"

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index, int num_elems,
                      double *cgup, int* epart, int* npart, int* objval) {
    int i;
    int my_rank, num_procs;
    /** Simulation parameters parsed from the input datasets */
    
    int nintci, nintcf;    /// internal cells start and end index
    int nextci, nextcf;
    int **lcc;    /// link cell-to-cell array - stores neighboring information
    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs, *be, *bn, *bw, *bl, *bh;
    double *bp;    /// Pole coefficient
    double *su;    /// Source values
    /** Additional vectors required for the computation */
    //double *distr, *oc, *cnorm;

    /** Geometry data */
    int points_count;    /// total number of points that define the geometry
    int** points;    /// coordinates of the points that define the cells - size [points_cnt][3]
    int* elems;    /// definition of the cells using their nodes (points) - each cell has 8 points
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processes
    

    int f_status = read_binary_geo(file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc, &bs,
                                   &be, &bn, &bw, &bl, &bh, &bp, &su, &points_count,
                                   &points, &elems);
    if ( f_status != 0 ) {
        fprintf(stderr, "Failed to initialize data!\n");
        MPI_Abort(MPI_COMM_WORLD, my_rank);
    }

    double *distr = (double*) calloc(sizeof(double), (nextcf + 1));

    for( i = (nintci); i <= (nintcf); i++ ) { 
       (distr)[i] = 0.0;
    }
    int npro = num_elems / num_procs;
    int exter = nextcf -nextci +1;
    int remain = 0;
    if (my_rank == (num_procs-1) ) {
        remain = num_elems % num_procs;
    }
    for ( i= 0; i < npro+remain+exter ;i++) {
          int k = local_global_index[i];
          (distr)[k] = cgup[i]; 
    }
    vtk_write_unstr_grid_header(file_in, file_vtk_out, nintci, nintcf, points_count, points, elems);
    vtk_append_double(file_vtk_out, "CGUP", nintci, nintcf, distr);    
    
    // Return an error if not implemented
    return -1;
    }

    int test_communication(char *file_in, char *file_vtk_out, int *local_global_index, int num_elems,
                     int neighbors_count, int* send_count, int** send_list, int* recv_count,
                   int** recv_list) {
   int i=0, j=0, kp=0, k=0;
    int my_rank, num_procs;
    /** Simulation parameters parsed from the input datasets */

    int nintci, nintcf;    /// internal cells start and end index
    int nextci, nextcf;
    int **lcc;    /// link cell-to-cell array - stores neighboring information
    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs, *be, *bn, *bw, *bl, *bh;
    double *bp;    /// Pole coefficient
    double *su;    /// Source values
    /** Additional vectors required for the computation */
    //double *distr, *oc, *cnorm;

    /** Geometry data */
    int points_count;    /// total number of points that define the geometry
    int** points;    /// coordinates of the points that define the cells - size [points_cnt][3]
    int* elems;    /// definition of the cells using their nodes (points) - each cell has 8 points
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processes


    int f_status = read_binary_geo(file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc, &bs,
                                   &be, &bn, &bw, &bl, &bh, &bp, &su, &points_count,
                                   &points, &elems);
    if ( f_status != 0 ) {
        fprintf(stderr, "Failed to initialize data!\n");
        MPI_Abort(MPI_COMM_WORLD, my_rank);
    }
    int npro = num_elems / num_procs;
     int remain = 0;
     if (my_rank == (num_procs-1) ) {
        remain = num_elems % num_procs;
    }    
    int exter = nextcf -nextci +1; 
    int *commlist = (int*) calloc(sizeof(int), (nextcf + 1));
     for ( i= 0; i < npro+remain+exter ;i++) {
          k = local_global_index[i];
          (commlist)[k] = 15;
    }
 
// for (i=0; i<nintcf; i++){
  //      (commlist)[i]=15;
    //    } 
    int *neighbors = (int*) calloc(sizeof(int), (6));
    printf("recv_list%d\n",recv_list[0][0]);
    for ( i = 0; i < num_procs; i++){
    //printf("i is:%d",i);
      // printf("send_count%d\n",send_count[i]);
    for ( j = 0; j < recv_count[i] ;j++) {
    //for ( i = 0; i < num_procs; i++){    
         k = recv_list[i][j];
      //   printf("i is:%d, j is %d\n",kp,j);
    //printf("recv_list%d\n",recv_list[0][j]);
          (commlist)[k] = 5;
    //if (k != 0) printf("k is:%d\n",k);
    } 
    for (kp =0; kp< send_count[i]; kp++){
         k=send_list[i][kp];
        (commlist)[k]=10;
    }
    }
   
    vtk_append_integer(file_vtk_out, "commlist", nintci, nintcf, commlist);
    
    // Return an error if not implemented
      return -1;
    }

