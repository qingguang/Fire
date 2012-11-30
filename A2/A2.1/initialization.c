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
    MPI_Status status[2];
    MPI_Request request[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processe
    // read-in the input file by one processor
    printf("processor number is%d\n", my_rank);
    
    if (my_rank==0) {      
     //

  int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &lcc_a, &bs_a,
                             &be_a, &bn_a, &bw_a, &bl_a, &bh_a, &bp_a, &su_a, &*points_count,
                               &*points, &*elems);
//int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                    //                  &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
                      //                &*points, &*elems);

       if ( f_status != 0 ) return f_status;
   printf("lcc %d,bp %f,bp_npro\n", (lcc_a)[0][0],bp_a[0]);
    }   

    //Send the information to other processors
    MPI_Bcast (nintci,1, MPI_INT, 0, MPI_COMM_WORLD);    
    MPI_Bcast (nintcf,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (nextci,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (nextcf,1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast (points_count,1, MPI_INT, 0, MPI_COMM_WORLD);

    int num_elems = *nintcf-*nintci+1;
    printf("number of elemts%d\n",num_elems);
    int points_num = *points_count;
    int npro = num_elems/num_procs;
    int remain = *nextcf - *nextci + 1;
    *local_global_index = (int*) calloc(sizeof(int), num_elems);
    *cgup_local = (double*) calloc(sizeof(double), npro);
    *epart = (int*) calloc(sizeof(int), num_elems);
    *npart = (int*) calloc(sizeof(int), num_elems*8);
    
    *bs = (double*) calloc(sizeof(double), (npro+remain));
    *bn = (double*) calloc(sizeof(double), (npro+remain));
    *bw = (double*) calloc(sizeof(double), (npro+remain));
    *be = (double*) calloc(sizeof(double), (npro+remain)); 
    *bl = (double*) calloc(sizeof(double), (npro+remain));
    *bh = (double*) calloc(sizeof(double), (npro+remain));
    *bp = (double*) calloc(sizeof(double), (npro+remain));
    *su = (double*) calloc(sizeof(double), (npro+remain));
    *var = (double*) calloc(sizeof(double), (npro+remain));
    *cgup = (double*) calloc(sizeof(double), (npro+remain));
    *oc = (double*) calloc(sizeof(double), (npro));
    *cnorm = (double*) calloc(sizeof(double), (npro));
    if (my_rank==0){ 
    printf("bp %f,bp_npro%f\n", bp_a[0], bp_a[npro]);
    }
    if (strcmp(part_type,"classical") == 0) {
    MPI_Scatter(bs_a, npro, MPI_DOUBLE, *bs, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatter(bn_a, npro, MPI_DOUBLE, *bn, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatter(bw_a, npro, MPI_DOUBLE, *bw, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatter(be_a, npro, MPI_DOUBLE, *be, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatter(bl_a, npro, MPI_DOUBLE, *bl, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatter(bh_a, npro, MPI_DOUBLE, *bh, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatter(bp_a, npro, MPI_DOUBLE, *bp, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatter(su_a, npro, MPI_DOUBLE, *su, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
    //MPI_Scatter(lcc_a, npro, MPI_DOUBLE, *lcc, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);

    printf("bs after Scan is: %f\n", (*bp)[0]);
    //Data distribution
    //classical data distribution from one processor to all another 
    for ( i = 0; i <= 10; i++ ) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    for ( i = 0; i < npro; i++ ) {
        (*cgup)[i] = 0.0;
        (*var)[i] = 0.0;
    }

    for ( i = npro; i < npro+remain; i++ ) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bl)[i] = 0.0;
        (*bh)[i] = 0.0;
    }
    
    for ( i = 0; i < npro; i++ )
       (*cgup)[i] = 1.0 / ((*bp)[i]);
       //MPI_Scatter(*cgup, npro, MPI_DOUBLE, *cgup_local, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);

    for ( i = 0; i < npro; i++ ) { 
        (*local_global_index)[i] = my_rank * npro + i; 
    }
    printf("local_global_index %d\n",(*local_global_index)[0]);    
   
    }else{
    
    if (my_rank==0){
    //Metis Dual
    idx_t ne = (idx_t) num_elems;
    idx_t nn = (idx_t) points_num;
    idx_t ncommon = 4;
    idx_t nparts = num_procs;
    int node_num = ne*8;
    idx_t *eptr = (idx_t*) calloc(sizeof(idx_t), num_elems + 1);

    for ( i = (*nintci); i <= (*nintcf + 1) ; i++ ) {
        eptr[i]=(idx_t) i*8;
    }
    idx_t *eind = (idx_t*) calloc(sizeof(idx_t), node_num);
    for(i = 0; i < node_num; i++ )
        eind[i] = (idx_t) (*elems)[i];
    // int* options[METIS_NOPTIONS];
    //options[METIS_OPTION_NUMBERING]=0;
    idx_t objval_METIS;
    idx_t *epart_METIS = (idx_t*) calloc(sizeof(idx_t), num_elems);
    idx_t *npart_METIS = (idx_t*) calloc(sizeof(idx_t), node_num);
    int metis_final;
    if (strcmp(part_type,"dual") == 0) {
       metis_final = METIS_PartMeshDual(&ne,&nn,eptr, eind, NULL, NULL, 
                                        &ncommon, &nparts, NULL,NULL, &objval_METIS, epart_METIS, npart_METIS); 
    } else if (strcmp(part_type,"nodal") == 0) {
                metis_final = METIS_PartMeshNodal(&ne,&nn,eptr, eind, NULL, NULL,
                                                  &nparts, NULL,NULL, &objval_METIS, epart_METIS, npart_METIS);
    }
    if (metis_final != METIS_OK) {
       printf("Metis part fails\n");
    }
    (*objval)=(int*) calloc(sizeof(int), 1);
    (*objval)[0]=(int) objval_METIS;
    for(i = 0; i < num_elems; i++ )
        (*epart)[i] = (int) epart_METIS[i];
    for(i = 0; i < node_num; i++ )
        (*npart)[i] = (int) npart_METIS[i]; 
    printf("epart is %d,%d, %d\n",(*epart)[0],(*epart)[4],(*epart)[num_elems]);
    //Full METIS arrary should be avaible for every processor
    //int j = 0;
    }//single processor 
    MPI_Bcast(*epart,num_elems,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(*npart,num_elems*8,MPI_INT,0,MPI_COMM_WORLD);
    //ditribute data according to METIS Partition
    MPI_Barrier(MPI_COMM_WORLD);
   // if (my_rank==0){
    //int j = 1;
   // int k = 1;
    //MPI_Sendrecv ( &k, 1, MPI_INT, 0, 10, &j, 1, MPI_INT, 1, MPI_INT, MPI_COMM_WORLD, &status[0]); 
    //  MPI_Isend(&j,1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &request[0]);
     // MPI_Irecv(&k,1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &request[1]);
      //MPI_Waitall(2, request, status);
   // (*bp)[0]=k;
    //if (my_rank==0)
    //printf("k is %f\n", (*bp)[0]);
    int p = 0;
    int j= 0;
    //int k = 0;
    int *k = (int*) calloc(sizeof(int), num_procs);
    for (p=0; p< num_procs; p++){
    //k=0;
   if (my_rank == p){
    //for ( i = 0; i < npro; i++ ) {
     //   (*local_global_index)[i] = my_rank * npro + i;
   // }
    //while  ( k < npro) {
    for (j =0; j < num_elems; j++){
        //double buf = (*cgup)[i];
     //if ( my_rank == (*epart)[j]){
        // (*cgup)[k] = (*cgup)[j];
          //k=k+1;
         if ( (*epart)[j] == my_rank){
         (*local_global_index)[k[my_rank]]= j ;
             k[my_rank]=k[my_rank]+1;
      //   (*bp)[k[my_rank]]= 
    }
    //printf("j is,%d,%d,%d,%f\n",num_elems,(*epart)[j],my_rank,(*local_global_index)[k],(*cgup_local)[k]);
    }
  // printf("pro is : %d and k0 is: %d k1 is %d\n", my_rank, (*local_global_index)[k[0]-1],(*local_global_index)[k[1]-1]);
   // } 
  }

   MPI_Bcast(&k[p],1,MPI_INT,p,MPI_COMM_WORLD);

   }
  //put new array store vector 
  int *k_sum = (int*) calloc(sizeof(int), num_procs);
 if ( my_rank == 0 ) {
 //double *bp_b = (double*) calloc(sizeof(double), (num_elems+num_procs));
 //for (i= 0; i<num_elems; i++){
 //*bp_b[i]= 
 //}
// int *k_sum = (int*) calloc(sizeof(int), num_procs);
 for (i =1; i<num_procs; i++)
 {
   k_sum[i]=k_sum[i-1]+k[i];
  }
}  
  //MPI_Sendrecv_replace( &k[1], 1, MPI_INT, 0, 10, &j, 1, MPI_INT, 1, MPI_INT, MPI_COMM_WORLD, &status[0]);  
  // MPI_Sendrecv_replace( &k[1], 1, MPI_INT, 0, 10, 1, 10, MPI_COMM_WORLD, &status[0]);  
  
   printf("pro is : %d and k0 is: %d k1 is %d\n", my_rank, k[0],k[1]);

  //MPI_Scatter(bp_a, npro, MPI_DOUBLE, *bp, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);
  MPI_Scatterv( bp_a, k, k_sum , MPI_DOUBLE, *bp, k[my_rank],MPI_DOUBLE,0, MPI_COMM_WORLD);  
  for ( i = 0; i < k[my_rank]; i++ )
       (*cgup)[i] = 1.0 / ((*bp)[i]);

    //printf("local processor%d,%d\n",my_rank, (*local_global_index)[npro]); 
   // MPI_Wait (&request, &status);
   // MPI_Irecv (recv_buffer, BUFSIZ, MPI_DOUBLE, 0, 77, MPI_COMM_WORLD,
    //           &request);
   // MPI_Wait (&request, &status);
   
    }//finish metis

    MPI_Barrier(MPI_COMM_WORLD);

    if ( my_rank == 0 ) {
       printf("processor 0 npro is%d,cgup%f \n", npro,(*cgup)[0]);
    }
    if (my_rank==1)
        printf("processor 1 cgup_local is: %f,  epart is%d \n",(*cgup)[0],(*local_global_index)[npro]);

    return 0;
    }
