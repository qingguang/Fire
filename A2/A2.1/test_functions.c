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
                      double *cgup_local, int* epart, int* npart, int* objval) {
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
    for (i=0; i<=num_elems/num_procs ;i++) {
    int k = local_global_index[i];
    (distr)[k] = cgup_local[i]; 
    }
    printf("processor 1 cgup[0] is%f ", cgup_local[0]);
    vtk_write_unstr_grid_header(file_in, file_vtk_out, nintci, nintcf, points_count, points, elems);
    vtk_append_double(file_vtk_out, "CGUP", nintci, nintcf, cgup_local);    
    
    // Return an error if not implemented
    return -1;
}

    //int test_communication(char *file_in, char *file_vtk_out, int *local_global_index, int *num_elems,
      //                 int neighbors_count, int* send_count, int** send_list, int* recv_count,
        //               int** recv_list) {
    // Return an error if not implemented
  //  return -1;
//}

