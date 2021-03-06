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
                   int*** recv_list, int** epart, int** npart, int** objval, double** cgup_local) {
    /********** START INITIALIZATION **********/
    int i = 0;
    int my_rank, num_procs;
    //int **lcc;    /// link cell-to-cell array - stores neighboring information
    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs_a, *be_a, *bn_a, *bw_a, *bl_a, *bh_a;
    double *bp_a;    /// Pole coefficient
    double *su_a;    /// Source values
    //int** points_a;    /// coordinates of the points that define the cells - size [points_cnt][3]
    //int* elems_a;    /// definition of the cells using their nodes (points) - each cell has 8 points
    int** lcc_a;    /// link cell-to-cell array - stores neighboring information
    
    //MPI_Status status[2];
    //MPI_Request request[2];
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processe
    // read-in the input file by one processor
    
    if ( my_rank == 0 ) {
         int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &lcc_a, &bs_a,
                                        &be_a, &bn_a, &bw_a, &bl_a, &bh_a, &bp_a, &su_a, &*points_count,
                                        &*points, &*elems);
         if ( f_status != 0 ) return f_status;
    }   

    //Send the common information to other processors
    MPI_Bcast (nintci,1, MPI_INT, 0, MPI_COMM_WORLD);    
    MPI_Bcast (nintcf,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (nextci,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (nextcf,1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast (points_count,1, MPI_INT, 0, MPI_COMM_WORLD);

    //local arrays and share parameters
    int num_elems = *nintcf-*nintci+1;
    int points_num = *points_count;
    int npro = num_elems/num_procs;
    int exter = *nextcf - *nextci + 1;
    int remain = 0;
    if (my_rank == (num_procs-1) ) {
        remain = num_elems % num_procs;
    }
    int local_array_size = npro + remain + exter;
    *local_global_index = (int*) calloc(sizeof(int), npro+exter);
    *epart = (int*) calloc(sizeof(int), num_elems);
    *npart = (int*) calloc(sizeof(int), num_elems*8);
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
    *oc = (double*) calloc(sizeof(double), (npro+remain));
    *cnorm = (double*) calloc(sizeof(double), (npro+remain));

    //choose part type 
    if (strcmp(part_type,"classical") == 0) {

        //ditribute all B* array
        MPI_Scatter(bs_a, npro, MPI_DOUBLE, *bs, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
        MPI_Scatter(bn_a, npro, MPI_DOUBLE, *bn, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
        MPI_Scatter(bw_a, npro, MPI_DOUBLE, *bw, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
        MPI_Scatter(be_a, npro, MPI_DOUBLE, *be, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
        MPI_Scatter(bl_a, npro, MPI_DOUBLE, *bl, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
        MPI_Scatter(bh_a, npro, MPI_DOUBLE, *bh, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
        MPI_Scatter(bp_a, npro, MPI_DOUBLE, *bp, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
        MPI_Scatter(su_a, npro, MPI_DOUBLE, *su, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
    //initialization of computational array 
    for ( i = 0; i <= 10; i++ ) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    for ( i = 0; i < npro; i++ ) {
        (*cgup)[i] = 0.0;
        (*var)[i] = 0.0;
    }

    for ( i = npro; i < npro+exter; i++ ) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bl)[i] = 0.0;
        (*bh)[i] = 0.0;
    }
    
    for ( i = 0; i < npro; i++ ) {
       (*cgup)[i] = 1.0 / ((*bp)[i]);
    }
    
    for ( i = 0; i < npro; i++ ) { 
        (*local_global_index)[i] = my_rank * npro + i; 
    }
    
    //part type is not classics but metis   
    }else{
    
    if ( my_rank == 0 ) {

         //parametes and array for metis partition libary
         idx_t ne = (idx_t) num_elems;
         idx_t nn = (idx_t) points_num;
         idx_t ncommon = 4;
         idx_t nparts = num_procs;
         int node_num = ne*8;
         idx_t *eptr = (idx_t*) calloc(sizeof(idx_t), num_elems + 1);
         idx_t *eind = (idx_t*) calloc(sizeof(idx_t), node_num);
         idx_t objval_METIS;
         idx_t *epart_METIS = (idx_t*) calloc(sizeof(idx_t), num_elems);
         idx_t *npart_METIS = (idx_t*) calloc(sizeof(idx_t), node_num);
         int metis_final;

    for ( i = (*nintci); i <= (*nintcf + 1) ; i++ ) {
          eptr[i] = (idx_t) i*8;
    }

    for( i = 0; i < node_num; i++ ) {
         eind[i] = (idx_t) (*elems)[i];
    }
   
    if ( strcmp(part_type,"dual") == 0 ) {
         metis_final = METIS_PartMeshDual(&ne,&nn,eptr, eind, NULL, NULL, 
                                          &ncommon, &nparts, NULL,NULL, &objval_METIS, epart_METIS, npart_METIS); 
    } else if ( strcmp(part_type,"noda") == 0 ) {
                metis_final = METIS_PartMeshNodal(&ne,&nn,eptr, eind, NULL, NULL,
                                                  &nparts, NULL,NULL, &objval_METIS, epart_METIS, npart_METIS);
    }
    
    if ( metis_final != METIS_OK ) {
         printf("Metis part fails\n");
    }
    (*objval)=(int*) calloc(sizeof(int), 1);
    (*objval)[0]=(int) objval_METIS;   
    for ( i = 0; i < num_elems; i++ ) {
          (*epart)[i] = (int) epart_METIS[i];
    }
    
    for ( i = 0; i < node_num; i++ ) {
          (*npart)[i] = (int) npart_METIS[i]; 
    }
    }//single processor 
    
    //Full METIS arrary should be avaible for every processor
    MPI_Bcast(*epart,num_elems,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(*npart,num_elems*8,MPI_INT,0,MPI_COMM_WORLD);
    //ditribute data according to METIS Partition
    MPI_Barrier(MPI_COMM_WORLD);
    int p = 0;
    int j= 0;
    int *k = (int*) calloc(sizeof(int), num_procs);

    //store local to global mapping
    for ( p = 0; p < num_procs; p++ ) {
    
    if (my_rank == p){
    
    for (j = 0; j < num_elems; j++){
         if ( (*epart)[j] == my_rank ) {
              (*local_global_index)[k[my_rank]] = j ;
              k[my_rank] = k[my_rank] + 1;
    }
    }
    }
    // send k[p] to other processors 
    MPI_Bcast(&k[p],1,MPI_INT,p,MPI_COMM_WORLD);
    }//finish storing local to global mapping

    //copy B* array into new array accoring to order from metis partition 
    int *k_sum = (int*) calloc(sizeof(int), num_procs);
    int *local_global_index_sum = (int*) calloc(sizeof(int), num_elems);
    if ( my_rank == 0 ) {
         for (i = 1; i < num_procs; i++ ) {
               k_sum[i]=k_sum[i-1]+k[i-1];
    }
    }  
    MPI_Gatherv( *local_global_index, k[my_rank], MPI_INT,
                 local_global_index_sum, k, k_sum,
                 MPI_INT, 0, MPI_COMM_WORLD);
    double *bs_b = (double*) calloc(sizeof(double), (num_elems));
    double *bn_b = (double*) calloc(sizeof(double), (num_elems));
    double *bw_b = (double*) calloc(sizeof(double), (num_elems));
    double *be_b = (double*) calloc(sizeof(double), (num_elems));
    double *bl_b = (double*) calloc(sizeof(double), (num_elems));
    double *bh_b = (double*) calloc(sizeof(double), (num_elems));
    double *bp_b = (double*) calloc(sizeof(double), (num_elems));
    double *su_b = (double*) calloc(sizeof(double), (num_elems));

    if (my_rank==0){ 
    for (i= 0; i<num_elems; i++){
        j=local_global_index_sum[i]; 
        bs_b[i]=bs_a[j];
        bn_b[i]=bn_a[j];
        bw_b[i]=bw_a[j];
        be_b[i]=be_a[j];
        bl_b[i]=bl_a[j];
        bh_b[i]=bh_a[j];
        bp_b[i]=bp_a[j];
        su_b[i]=su_a[j]; 
    }
    }
    MPI_Scatterv( bs_b, k, k_sum , MPI_DOUBLE, *bs, k[my_rank],MPI_DOUBLE,0, MPI_COMM_WORLD);  
    MPI_Scatterv( bn_b, k, k_sum , MPI_DOUBLE, *bn, k[my_rank],MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatterv( bw_b, k, k_sum , MPI_DOUBLE, *bw, k[my_rank],MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatterv( be_b, k, k_sum , MPI_DOUBLE, *be, k[my_rank],MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatterv( bl_b, k, k_sum , MPI_DOUBLE, *bl, k[my_rank],MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatterv( bh_b, k, k_sum , MPI_DOUBLE, *bh, k[my_rank],MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatterv( bp_b, k, k_sum , MPI_DOUBLE, *bp, k[my_rank],MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatterv( su_b, k, k_sum , MPI_DOUBLE, *su, k[my_rank],MPI_DOUBLE,0, MPI_COMM_WORLD);

    //initialization computational array
    for ( i = 0; i <= 10; i++ ) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    for ( i = 0; i < k[my_rank]; i++ ) {
        (*cgup)[i] = 0.0;
        (*var)[i] = 0.0;
    }

    for ( i = k[my_rank]; i < npro+exter; i++ ) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bl)[i] = 0.0;
        (*bh)[i] = 0.0;
    }
 
    for ( i = 0; i < k[my_rank]; i++ )
          (*cgup)[i] = 1.0 / ((*bp)[i]);

    MPI_Barrier(MPI_COMM_WORLD);

    free(bp_b);
    free(bh_b);
    free(bl_b);
    free(bw_b);
    free(bn_b);
    free(be_b);
    free(bs_b);  
    free(su_b);
    free(local_global_index_sum);
    
    }//finish choose part type section and all local array are stored

    MPI_Barrier(MPI_COMM_WORLD);

    //printf("processor %d executes sucessfully\n", my_rank);

    return 0;
    }
