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
                     int* send_count, int** send_list, int* recv_count, int** recv_list, int num_elems_local, int* epart) {
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int my_rank, num_procs;
    int i = 0, j=0, k=0;
    double residual_ratio_sum; 
    MPI_Status status;
    int rank;
    MPI_Datatype send_type;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    // Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    // get number of processe
    // allocate arrays used in gccg
    int nomax = 3;
    int blocklengths[1] = {send_count[0]};
    int displacements[1] = {send_count[0]};

    /** the reference residual*/
    double resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (num_elems_local));
    printf("num_elems_local is:%d\n",num_elems_local); 
    // initialize the reference residual
    for ( nc = 0; nc < num_elems_local; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }
//   MPI_Allreduce(&resref, &resref, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    resref = sqrt(resref);
    if ( resref < 1.0e-15 ) {
        fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
        return 0;
    }
   //printf("myrank is: %d,send_count is: %d:%d:%d, recv_count is: %d:%d:%d\n",
   // my_rank,send_count[0],send_count[1],send_count[2],recv_count[0],recv_count[1],recv_count[2]);
   //printf("num_elems%d,global_local_index[num_elems_pro+1]:%d\n",nintcf-nintci,global_local_index[nintcf-nintci]);
  /** the computation vectors */
    double *direc1 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (num_elems_local));
    double *adxor2 = (double *) calloc(sizeof(double), (num_elems_local));
    double *dxor1 = (double *) calloc(sizeof(double), (num_elems_local));
    double *dxor2 = (double *) calloc(sizeof(double), (num_elems_local));
    int    nc_global = 0;
    double *direc1_sum = (double *) calloc(sizeof(double), (nextcf + 1));
    int send_sum = 0;
    for (i=0; i < num_procs; i++){
         send_sum =send_sum + send_count[i]; 
        } 
   double *direc1_send = (double *) calloc(sizeof(double), (send_sum));
   printf("blocklengths is:%d\n",blocklengths[0]);
//MPI_Barrier(MPI_COMM_WORLD);

    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/
        // update the old values of direc
        for ( nc = 0; nc < num_elems_local; nc++ ) {
           nc_global= local_global_index[nc];
          // direc1[nc_global] = direc1[nc_global] + resvec[nc] * cgup[nc]; 
           direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
           direc1_sum[nc_global]=direc1[nc];
//           MPI_Bcast(&direc1_sum[nc_global],1,MPI_DOUBLE,my_rank,MPI_COMM_WORLD);
        }

//MPI_Barrier(MPI_COMM_WORLD);
if (iter == 10)      printf("direc1[0] is : %e\n", direc1[0]);
int index=0;
for (i=0;i<num_procs;i++){
for (j=0;j<send_count[i];j++){
k=global_local_index[send_list[i][j]];
direc1_send[index]=direc1[k];
index=index+1;
}
}

if (my_rank==0){

}
//printf("direc1_send[send_sum-1] is:%e\n",direc1_send[send_sum-1]);

/*
for (k = 0; k < num_procs; k++){
if (k != my_rank){
{
for (i=0; i < send_count[k]; i++)
{
   MPI_Recv(&values[100/np-1], 1, MPI_INT, rnbr, 10, MPI_COMM_WORLD, &status);
    MPI_Send(&buf, 1, MPI_INT, lnbr, 10, MPI_COMM_WORLD);
 j = global_local_index(send_list[k][i]);
MPI_Send(&direc1[j],1,MPI_DOUBLE,k,10,MPI_COMM_WORLD);
MPI_Recv(&buffer[i],1,MPI_DOUBLE,20,MPI_COMM_WORLD,&status);
}
for (j=0; j<recv_count[k]; j++)
{
 MPI_recv(buffer[k][j],)
}
}
}
}*/
//printf("send_count is: %d\n",send_count[0]);
//MPI_Barrier(MPI_COMM_WORLD);
 
       // compute new guess (approximation) for direc
       
        for ( nc = 0; nc < num_elems_local; nc++ ) {
///lcc_1=lcc[nc][3]

        /*for (i = 0; i < 6; i++) {
           if (lcc[nc][i] >= nintcf)
               lcc[nc][i] =   
         
        }*/
          direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[nc][0]]
                         - bw[nc] * direc1[lcc[nc][3]] - bl[nc] * direc1[lcc[nc][4]]
                         - bn[nc] * direc1[lcc[nc][2]] - be[nc] * direc1[lcc[nc][1]]
                         - bh[nc] * direc1[lcc[nc][5]];
        }
       
        /********** END COMP PHASE 1 **********/

        /********** START COMP PHASE 2 **********/
        // execute normalization steps
        double oc1, oc2, occ;
        double oc1_sum=0, oc2_sum, occ_sum;
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

            for ( nc = 0; nc < num_elems_local; nc++ ) {
                occ = occ + adxor1[nc] * direc2[nc];
            }
            
            oc1 = occ / cnorm[1];
            //printf("oc1 is:%e\n",oc1);
        //    MPI_Barrier(MPI_COMM_WORLD);

          // MPI_Allreduce(&oc1, &oc1_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
             //printf("oc1_sum is:%e\n",oc1_sum);
            for ( nc = 0; nc < num_elems_local; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }

            if1++;
        } else {
            if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;

                for ( nc = 0; nc < num_elems_local; nc++ ) {
                    occ = occ + adxor1[nc] * direc2[nc];
                }

                oc1 = occ / cnorm[1];
        //        MPI_Allreduce(&oc1, &oc1_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                oc2 = 0;
                occ = 0;
                for ( nc = 0; nc <num_elems_local; nc++ ) {
                    occ = occ + adxor2[nc] * direc2[nc];
                }

                oc2 = occ / cnorm[2];
        ///        MPI_Allreduce(&oc2, &oc2_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for ( nc = 0; nc < num_elems_local; nc++ ) {
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                    direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                }

                if2++;
            }
        }
//        MPI_Barrier(MPI_COMM_WORLD);
        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
        double omega_sum = 0;
        for ( nc = 0; nc < num_elems_local; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }

        omega = omega / cnorm[nor];
   //     MPI_Allreduce(&omega, &omega_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double res_updated = 0.0;
        for ( nc = 0; nc < num_elems_local; nc++ ) {
            var[nc] = var[nc] + omega * direc1[nc];
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
        }
//        MPI_Allreduce(&res_updated, &res_updated, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        res_updated = sqrt(res_updated);
        *residual_ratio = res_updated / resref;
        //MPI_Allreduce(residual_ratio, &residual_ratio_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 
        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 ) break;

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( nc = 0; nc <num_elems_local; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = 0; nc < num_elems_local; nc++ ) {
                        dxor2[nc] = direc1[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
    }
           
    //printf("myrank%d,residual_sum is %e\n",my_rank,residual_ratio_sum);
    printf("rank is: %d, iter is: %d, residual ration is: %e\n",my_rank,iter,*residual_ratio);
    printf("pent iter is: %d,residual ration is: %e\n",546,9.695255e-11);
    free(resvec);
    free(direc1);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);

    return iter;
}


