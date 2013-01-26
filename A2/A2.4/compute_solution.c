/**
 * Computational loop
 *
 * @file compute_solution.c
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int compute_solution(const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                     double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                     double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                     int* local_global_index, int* global_local_index, int neighbors_count,
                     int* send_count, int** send_list, int* recv_count, int** recv_list,
                     int num_elems_local, int* epart) {
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int my_rank, num_procs;
    int i = 0, j = 0, k = 0;
    double residual_ratio_sum;
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    // Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    // get number of processe
    MPI_Datatype send_type[num_procs], recv_type[num_procs];
    // allocate arrays used in gccg
    int nomax = 3;

    /** the reference residual*/
    double resref = 0.0;
    double resref_sum = 0.0;
    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (num_elems_local));
    //printf("num_elems_local is:%d\n",num_elems_local); 
    // initialize the reference residual
    for ( nc = 0; nc < num_elems_local; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }
    MPI_Allreduce(&resref, &resref_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    resref_sum = sqrt(resref_sum);
    if ( resref_sum < 1.0e-15 ) {
        fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
        return 0;
    }
 
   
  /** the computation vectors */
    double *direc1 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (num_elems_local));
    double *adxor2 = (double *) calloc(sizeof(double), (num_elems_local));
    double *dxor1 = (double *) calloc(sizeof(double), (num_elems_local));
    double *dxor2 = (double *) calloc(sizeof(double), (num_elems_local));
    double *cnorm_sum = (double *) calloc(sizeof(double), 4);
    int block_lens[nextcf+1];
    for (i = 0; i < 4; i++){
      cnorm_sum[i] = 1;}
    int nc_global = 0;
    int  *block_len_send;
    int *block_len_recv;
    // Define MPI Datatype
    for (i = 0; i < num_procs; i++) {
        block_len_send = (int* ) calloc(sizeof(int), (send_count[i]));
        block_len_recv = (int* ) calloc(sizeof(int), (recv_count[i]));
        if (send_count[i] != 0) {
            for (k = 0 ; k < send_count[i]; k++) {
                 block_len_send[k] = 1;
            }
            MPI_Type_indexed(send_count[i], block_len_send,
                             send_list[i], MPI_DOUBLE, &send_type[i]);
            MPI_Type_commit(&send_type[i]);
            free(block_len_send);
        }
        if (recv_count[i] != 0) {
             for (k = 0; k < recv_count[i]; k++) { 
                  block_len_recv[k] = 1;
             }
        MPI_Type_indexed(recv_count[i], block_len_recv,
                         recv_list[i], MPI_DOUBLE, &recv_type[i]);
        MPI_Type_commit(&recv_type[i]);
        free(block_len_recv);
        }
    }

    while ( iter <  max_iters ) {
    /**********  START COMP PHASE 1 **********/
    // update the old values of direc
    for ( nc = 0; nc < num_elems_local; nc++ ) {
	  nc_global= local_global_index[nc];   
          direc1[nc_global] = direc1[nc_global] + resvec[nc] * cgup[nc];
        }

    // communicate the direc1

    for (i = 0; i < num_procs; i++) {
         if (send_count[i] != 0) {
             MPI_Sendrecv(direc1,1, send_type[i], i, i, 
			  direc1,1, recv_type[i], i, my_rank, MPI_COMM_WORLD, &status);
	 }
    }

    // compute new guess (approximation) for direc
    for ( nc = 0; nc < num_elems_local; nc++ ) {  
	  nc_global = local_global_index[nc];   
          direc2[nc] = bp[nc] * direc1[nc_global] - bs[nc] * direc1[lcc[nc][0]]
                       - bw[nc] * direc1[lcc[nc][3]] - bl[nc] * direc1[lcc[nc][4]]
                       - bn[nc] * direc1[lcc[nc][2]] - be[nc] * direc1[lcc[nc][1]]
                       - bh[nc] * direc1[lcc[nc][5]];
    }
       
    /********** END COMP PHASE 1 **********/

    /********** START COMP PHASE 2 **********/
    // execute normalization steps
    double oc1, oc2, occ;
    double occ_sum;
            
    if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

        for ( nc = 0; nc < num_elems_local; nc++ ) {
              occ = occ + adxor1[nc] * direc2[nc];
        }
           
        MPI_Allreduce(&occ, &occ_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        oc1 = occ_sum / cnorm_sum[1];
           
        for ( nc = 0; nc < num_elems_local; nc++ ) {
              direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
              nc_global=local_global_index[nc]; 
              direc1[nc_global] = direc1[nc_global] - oc1 * dxor1[nc];
        }
        if1++;
        } else {
          if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;

          for ( nc = 0; nc < num_elems_local; nc++ ) {
                occ = occ + adxor1[nc] * direc2[nc];
          }
                MPI_Allreduce(&occ, &occ_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                oc1 = occ_sum / cnorm_sum[1];
                oc2 = 0;
                occ = 0;
          for ( nc = 0; nc <num_elems_local; nc++ ) {
                occ = occ + adxor2[nc] * direc2[nc];
          }
                MPI_Allreduce(&occ, &occ_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                oc2 = occ_sum / cnorm_sum[2];
          for ( nc = 0; nc < num_elems_local; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                nc_global=local_global_index[nc];      
                direc1[nc_global] = direc1[nc_global] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
          }

          if2++;
          }
        }
        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
        double omega_sum = 0;
        for ( nc = 0; nc < num_elems_local; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }
        MPI_Allreduce(&cnorm[nor], &cnorm_sum[nor], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&omega, &omega_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        omega = omega_sum / cnorm_sum[nor];
  
        double res_updated = 0.0;
	double res_updated_sum;
        for ( nc = 0; nc < num_elems_local; nc++ ) {
            nc_global=local_global_index[nc];

            var[nc] = var[nc] + omega * direc1[nc_global];
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
        }

        MPI_Allreduce(&res_updated, &res_updated_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        res_updated_sum = sqrt(res_updated_sum);
   
        *residual_ratio = res_updated_sum / resref_sum;
 
 
        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 ) break;

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( nc = 0; nc <num_elems_local; nc++ ) {
                    nc_global=local_global_index[nc];
                    dxor1[nc] = direc1[nc_global];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = 0; nc < num_elems_local; nc++ ) {
                  nc_global=local_global_index[nc];                  
                        dxor2[nc] = direc1[nc_global];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
    }
           
    // printf("myrank%d,residual_sum is %e\n",my_rank,residual_ratio_sum);
    // if (my_rank == 0){
    // printf("rank is: %d, iter is: %d, residual ration is: %e\n",my_rank,iter,*residual_ratio);
    // printf("pent iter is: %d,residual ration is: %e\n",546,9.695255e-11);
    // }
    free(resvec);
    free(direc1);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);

    return iter;
}


