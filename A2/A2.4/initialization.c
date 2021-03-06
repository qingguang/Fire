/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdlib.h>
#include "metis.h"
#include "util_read_files.h"
#include "initialization.h"
#include "mpi.h"
int initialization(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index, int** global_local_index,
                   int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
                   int*** recv_list, int** epart, int** npart, int** objval, int* num_elems_local) {
    /********** START INITIALIZATION **********/
    int i = 0;
    int j = 0;
    int num_elems_pro;  // number of elements in each processor
    int my_rank, num_procs;

    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs_a, *be_a, *bn_a, *bw_a, *bl_a, *bh_a;
    double *bp_a;    /// Pole coefficient
    double *su_a;    /// Source values
    int** lcc_a;    /// link cell-to-cell array - stores neighboring information
    int** lcc_b;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processe

    // read-in the input file by one processor
    if ( my_rank == 0 ) {
         int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf,
                                        &lcc_a, &bs_a, &be_a, &bn_a, &bw_a, &bl_a, &bh_a,
                                        &bp_a, &su_a, &*points_count, &*points, &*elems);
         if ( f_status != 0 ) return f_status;
    }

    // Send the common information to other processors
    MPI_Bcast(nintci, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nintcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nextci, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nextcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(points_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // local arrays and share parameters
    int num_elems = *nintcf - *nintci + 1;
    if (my_rank != 0) {
        *elems = (int*) calloc(sizeof(int), num_elems * 8);
    }
    MPI_Bcast(*elems, num_elems * 8, MPI_INT, 0, MPI_COMM_WORLD);
    int points_num = *points_count;
    int npro = num_elems / num_procs;
    int exter = *nextcf - *nextci + 1;
    int remain = 0;
    int *k = (int*) calloc(sizeof(int), num_procs);
    int *k_sum = (int*) calloc(sizeof(int), num_procs);
    int last_proc = num_procs - 1;
    if (my_rank == last_proc) {
        remain = num_elems % num_procs;
    }
    int local_array_size = npro + remain + exter;
    *local_global_index = (int*) calloc(sizeof(int), npro + remain + exter);
    *global_local_index = (int*) calloc(sizeof(int), num_elems);
    *bs = (double*) calloc(sizeof(double), (local_array_size));
    *bn = (double*) calloc(sizeof(double), (local_array_size));
    *bw = (double*) calloc(sizeof(double), (local_array_size));
    *be = (double*) calloc(sizeof(double), (local_array_size));
    *bl = (double*) calloc(sizeof(double), (local_array_size));
    *bh = (double*) calloc(sizeof(double), (local_array_size));
    *bp = (double*) calloc(sizeof(double), (local_array_size));
    *su = (double*) calloc(sizeof(double), (local_array_size));
    *var = (double*) calloc(sizeof(double), (local_array_size));
    *cgup = (double*) calloc(sizeof(double), (local_array_size));
    *oc = (double*) calloc(sizeof(double), (npro + remain));
    *cnorm = (double*) calloc(sizeof(double), (npro + remain));
    *lcc = (int**) calloc(sizeof(int*), (local_array_size));
    for ( i = 0; i < local_array_size; i++ ) {
         (*lcc)[i] = (int *) calloc(sizeof(int), (6));
    }
    int *data = (int *) calloc(sizeof(int), (num_elems*6));
    lcc_b = (int **) calloc(sizeof(int*), (num_elems));
    for (i = 0; i < num_elems; i++) {
        lcc_b[i] = &(data[6 * i]);
    }
    if ( my_rank == 0 ) {
         for ( i = 0; i< num_elems; i++ ) {
         for ( j = 0; j < 6; j++ ) {
               lcc_b[i][j] = lcc_a[i][j];
    }
    }
    }

    MPI_Bcast(&(lcc_b[0][0]), num_elems*6, MPI_INT, 0, MPI_COMM_WORLD);
    // choose part type
    if (strcmp(part_type, "classical") == 0) {
        k[my_rank] = npro + remain;
        (num_elems_pro) = npro + remain;
         int p = 0;
    for (p = 0; p < num_procs; p++) {
        MPI_Bcast(&k[p], 1, MPI_INT, p, MPI_COMM_WORLD);
    }

    if (my_rank == 0) {
         for (i = 1; i < num_procs; i++) {
               k_sum[i] = k_sum[i-1] + k[i-1];
    }
    }

    // ditribute all B* array
    MPI_Scatterv(bs_a, k, k_sum, MPI_DOUBLE, *bs, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bn_a, k, k_sum, MPI_DOUBLE, *bn, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bw_a, k, k_sum, MPI_DOUBLE, *bw, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(be_a, k, k_sum, MPI_DOUBLE, *be, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bl_a, k, k_sum, MPI_DOUBLE, *bl, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bh_a, k, k_sum, MPI_DOUBLE, *bh, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bp_a, k, k_sum, MPI_DOUBLE, *bp, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(su_a, k, k_sum, MPI_DOUBLE, *su, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i = 0; i < num_elems_pro; i++) {
          (*local_global_index)[i] = my_rank * npro + i;
       for (j = 0; j < 6; j++) {
            (*lcc)[i][j] = lcc_b[my_rank*npro+i][j];
       }
    }
    for (i = 0; i < num_elems_pro; i++) {
          if (i > npro) {
              (*global_local_index)[my_rank*npro+i] = (my_rank*npro+i) % npro + npro;
          } else {
            (*global_local_index)[my_rank*npro+i] = (my_rank*npro+i) % npro;
          }
    }
    // part type is not classics but metis
    } else {
         *epart = (int*) calloc(sizeof(int), num_elems);
         *npart = (int*) calloc(sizeof(int), num_elems*8);
    // if ( my_rank == 0 ) {
         // parametes and array for metis partition libary
         idx_t ne = (idx_t) num_elems;
         idx_t nn = (idx_t) points_num;
         idx_t ncommon = 4;
         idx_t nparts = num_procs;
         int node_num = ne * 8;
         idx_t *eptr = (idx_t*) calloc(sizeof(idx_t), num_elems + 1);
         idx_t *eind = (idx_t*) calloc(sizeof(idx_t), node_num);
         idx_t objval_METIS;
         idx_t *epart_METIS = (idx_t*) calloc(sizeof(idx_t), num_elems);
         idx_t *npart_METIS = (idx_t*) calloc(sizeof(idx_t), node_num);
         int metis_final;
    for (i = (*nintci); i <= (*nintcf + 1) ; i++) {
          eptr[i] = (idx_t) i * 8;
    }

    for (i = 0; i < node_num; i++) {
         eind[i] = (idx_t) (*elems)[i];
    }

    if (strcmp(part_type, "dual") == 0) {
         metis_final = METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL,
                                          &ncommon, &nparts, NULL, NULL, &objval_METIS,
                                          epart_METIS, npart_METIS);
    } else if (strcmp(part_type, "noda") == 0) {
                metis_final = METIS_PartMeshNodal(&ne, &nn, eptr, eind, NULL, NULL,
                                                  &nparts, NULL, NULL, &objval_METIS,
                                                  epart_METIS, npart_METIS);
    }

    if (metis_final != METIS_OK) {
         printf("Metis part fails\n");
    }
    (*objval) = (int*) calloc(sizeof(int), 1);
    (*objval)[0] = (int) objval_METIS;
    for (i = 0; i < num_elems; i++) {
          (*epart)[i] = (int) epart_METIS[i];
    }

    for (i = 0; i < node_num; i++) {
          (*npart)[i] = (int) npart_METIS[i];
    }

    // ditribute data according to METIS Partition
    int p = 0;
    // store local to global mapping
    for (p = 0; p < num_procs; p++) {
          if (my_rank == p) {
              for (j = 0; j < num_elems; j++) {
                   if ((*epart)[j] == my_rank) {
                        (*local_global_index)[k[my_rank]] = j;
                       for (i = 0; i < 6; i++) {
                           (*lcc)[k[my_rank]][i] = lcc_b[j][i];
                       }
                       (*global_local_index)[j] = k[my_rank];
                       k[my_rank] = k[my_rank] + 1;
                   }
              }
          }
          MPI_Bcast(&k[p], 1, MPI_INT, p, MPI_COMM_WORLD);    /// send k[p] to other processors
    }    /// finish storing local to global mapping
    (num_elems_pro) = k[my_rank];
    int *local_global_index_sum = (int*) calloc(sizeof(int), num_elems);
    if (my_rank == 0) {
         for (i = 1; i < num_procs; i++) {
               k_sum[i] = k_sum[i-1] + k[i-1];
         }
    }
    MPI_Gatherv(*local_global_index, k[my_rank], MPI_INT,
                 local_global_index_sum, k, k_sum,
                 MPI_INT, 0, MPI_COMM_WORLD);
    // copy B* array into new array accoring to order from metis partition
    double *bs_b = (double*) calloc(sizeof(double), (num_elems));
    double *bn_b = (double*) calloc(sizeof(double), (num_elems));
    double *bw_b = (double*) calloc(sizeof(double), (num_elems));
    double *be_b = (double*) calloc(sizeof(double), (num_elems));
    double *bl_b = (double*) calloc(sizeof(double), (num_elems));
    double *bh_b = (double*) calloc(sizeof(double), (num_elems));
    double *bp_b = (double*) calloc(sizeof(double), (num_elems));
    double *su_b = (double*) calloc(sizeof(double), (num_elems));
    if (my_rank == 0) {
        for (i= 0; i < num_elems; i++) {
             j = local_global_index_sum[i];
             bs_b[i] = bs_a[j];
             bn_b[i] = bn_a[j];
             bw_b[i] = bw_a[j];
             be_b[i] = be_a[j];
             bl_b[i] = bl_a[j];
             bh_b[i] = bh_a[j];
             bp_b[i] = bp_a[j];
             su_b[i] = su_a[j];
         }
    }
    MPI_Scatterv(bs_b, k, k_sum , MPI_DOUBLE, *bs, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bn_b, k, k_sum , MPI_DOUBLE, *bn, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bw_b, k, k_sum , MPI_DOUBLE, *bw, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(be_b, k, k_sum , MPI_DOUBLE, *be, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bl_b, k, k_sum , MPI_DOUBLE, *bl, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bh_b, k, k_sum , MPI_DOUBLE, *bh, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(bp_b, k, k_sum , MPI_DOUBLE, *bp, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(su_b, k, k_sum , MPI_DOUBLE, *su, k[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(bp_b);
    free(bh_b);
    free(bl_b);
    free(bw_b);
    free(bn_b);
    free(be_b);
    free(bs_b);
    free(su_b);
    free(local_global_index_sum);
    }  // finish choose part type section and all local array are stored

    // initialization computational array
    for (i = 0; i <= 10; i++) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    for (i = 0; i < num_elems_pro; i++) {
        (*cgup)[i] = 0.0;
        (*var)[i] = 0.0;
    }

    for (i = num_elems_pro; i < local_array_size; i++) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bl)[i] = 0.0;
        (*bh)[i] = 0.0;
    }

    for (i = 0; i < num_elems_pro; i++) {
       (*cgup)[i] = 1.0 / ((*bp)[i]);
    }

    // ************Comunication List********* //
    *num_elems_local = num_elems_pro;
    *neighbors_count = num_procs-1;
    *send_count = (int*) calloc(sizeof(int), (num_procs));
    *recv_count = (int*) calloc(sizeof(int), (num_procs));
    *send_list = (int **) calloc(sizeof(int*), (*neighbors_count + 1));

    for (i = 0; i < *neighbors_count + 1; i++) {
        (*send_list)[i] = (int *) calloc(sizeof(int), (6 * num_elems_pro));
    }

    int num_elems_global = 0;
    int* rank = (int*) calloc(sizeof(int), (num_procs));
    int m = 0;
    int** count_time_local = (int**) calloc(sizeof(int*), (num_procs));
    for (i = 0; i < num_procs; i++) {
        count_time_local[i] = (int *) calloc(sizeof(int), (num_elems));
    }
    int* count_time = (int*) calloc(sizeof(int), (num_elems));
    for (i = 0; i < num_elems_pro; i++) {
    for (j = 0; j < 6; j++) {
         num_elems_global = (*lcc)[i][j];
    // choose only ghost cell
    if (num_elems_global < num_elems) {
    // choose part type
    if (strcmp(part_type, "classical") == 0) {
        if (num_elems_global >= npro * num_procs) {
            rank[my_rank]= num_elems_global / npro - 1;
        } else {
            rank[my_rank]= num_elems_global / npro;
        }
    } else {
          rank[my_rank]=(*epart)[num_elems_global];
    }   /// end choosing part type
    // record times of this elems occur
    if (rank[my_rank] != my_rank) {
    count_time_local[rank[my_rank]][i] = count_time_local[rank[my_rank]][i] + 1;
    if (count_time_local[rank[my_rank]][i] == 1) {
    (*send_list)[rank[my_rank]][(*send_count)[rank[my_rank]]] = (*local_global_index)[i];
    (*send_count)[rank[my_rank]]=(*send_count)[rank[my_rank]] + 1;
    }
    }
    }    /// choose ghost cell
    }    /// end j for loop
    }    /// end i for loop
    // Set the order in send_list and recv_list same
    for (i = 0; i < num_procs; i++) {
        MPI_Sendrecv(&((*send_count)[i]), 1, MPI_INT, i, i * 1000,
                     &((*recv_count)[i]), 1, MPI_INT, i,
                     my_rank * 1000, MPI_COMM_WORLD, &status);
    }
    *recv_list = (int **) calloc(sizeof(int*), (*neighbors_count+1));
    for (i = 0; i < *neighbors_count+1; i++) {
         (*recv_list)[i] = (int *) calloc(sizeof(int), ((*recv_count)[i]));
    }

    for (i = 0; i < num_procs; i++) {
         if (my_rank == i) {
                for (j = 0; j < num_procs; j++) {
                    if (j != my_rank) {
                        MPI_Send((*send_list)[j], (*send_count)[j], MPI_INT, j, 100,
                                 MPI_COMM_WORLD);
                    }
                }
         } else {
                MPI_Recv((*recv_list)[i], (*recv_count)[i], MPI_INT, i, 100,
                         MPI_COMM_WORLD, &status);
           }
    }

    free(lcc_b);
    free(count_time_local);
    return 0;
    }
