/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include "util_write_files.h"
#include "mpi.h"
void finalization(char* file_in, char* out_prefix, int total_iters, double residual_ratio,
                  int nintci, int nintcf, int points_count, int** points, int* elems, double* var,
                  double* cgup, double* su,  int* local_global_index, int num_elems_local) {
    char file_out[100];
    sprintf(file_out, "%s_summary.out", out_prefix);
    int my_rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processe
    int num_elems = nintcf-nintci+1;
    int nc_global= 0;
    int i;
    double *su_global = (double*) calloc(sizeof(double), (num_elems));
    double *var_global = (double*) calloc(sizeof(double), (num_elems));
    double *cgup_global = (double*) calloc(sizeof(double), (num_elems));

    for (i = 0; i < num_elems_local; i++) {
      nc_global = local_global_index[i];
      var_global[nc_global] = var[i];
      su_global[nc_global] = su[i];
      cgup_global[nc_global] = cgup[i];
    }
    double *su_global_sum = (double*) calloc(sizeof(double), (num_elems));
    double *var_global_sum = (double*) calloc(sizeof(double), (num_elems));
    double *cgup_global_sum = (double*) calloc(sizeof(double), (num_elems));
    MPI_Reduce(&var_global[0], &var_global_sum[0], num_elems,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&su_global[0], &su_global_sum[0], num_elems,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&cgup_global[0], &cgup_global_sum[0], num_elems,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_rank == 0) {
        int status = store_simulation_stats(file_in, file_out, nintci, nintcf,
                                            var_global, total_iters, residual_ratio);
        sprintf(file_out, "%s_data.vtk", out_prefix);
        vtk_write_unstr_grid_header(file_in, file_out, nintci, nintcf, points_count, points, elems);
        vtk_append_double(file_out, "CGUP", nintci, nintcf, cgup_global_sum);
    }

    free(su_global);
    free(var_global);
    free(cgup_global);
    free(su_global_sum);
    free(var_global_sum);
    free(cgup_global_sum);
    // if (status != 0) {
       // fprintf(stderr, "Error when trying to write to file %s\n", file_out);
       // }
    }

