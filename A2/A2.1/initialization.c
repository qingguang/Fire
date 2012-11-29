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
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processe
    // read-in the input file by one processor
    printf("processor number is%d\n", my_rank);
    if (my_rank==0) {      
       int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                                      &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
                                      &*points, &*elems);

       if ( f_status != 0 ) return f_status;
    }    
    //Send the information to other processors
    MPI_Bcast (nintci,1, MPI_INT, 0, MPI_COMM_WORLD);    
    MPI_Bcast (nintcf,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (nextci,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (nextcf,1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast (points_count,1, MPI_INT, 0, MPI_COMM_WORLD);
    *var = (double*) calloc(sizeof(double), (*nextcf + 1));
    *cgup = (double*) calloc(sizeof(double), (*nextcf + 1));
    *oc = (double*) calloc(sizeof(double), (*nintcf + 1));
    *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));
    int num_elems = *nintcf-*nintci+1;
    int points_num = *points_count;
    int npro = num_elems/num_procs;
    *local_global_index = (int*) calloc(sizeof(int), npro);
    *cgup_local = (double*) calloc(sizeof(double), npro);
    *epart = (int*) calloc(sizeof(int), num_elems);
    *npart = (int*) calloc(sizeof(int), num_elems*8);
    if (my_rank==0)
    {

    // initialize the arrays
    for ( i = 0; i <= 10; i++ ) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ ) {
        (*cgup)[i] = 0.0;
        (*var)[i] = 0.0;
    }

    for ( i = (*nextci); i <= (*nextcf); i++ ) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bl)[i] = 0.0;
        (*bh)[i] = 0.0;
    }
    
    for ( i = (*nintci); i <= (*nintcf); i++ )
       (*cgup)[i] = 1.0 / ((*bp)[i]);
    //End of reading and initialize by one processor
    } 
 
    //Data distribution
    //classical data distribution from one processor to all another 
    if (strcmp(part_type,"classical") == 0) {
       MPI_Scatter(*cgup, npro, MPI_DOUBLE, *cgup_local, npro,MPI_DOUBLE,0, MPI_COMM_WORLD);

    for ( i = 0; i <= npro; i++ ) { 
        (*local_global_index)[i] = my_rank * npro + i; 
    }
    printf("local_global_index %d\n",(*local_global_index)[0]);
    //MPI_Bcast(*cgup, *nintcf+1, MPI_DOUBLE,0, MPI_COMM_WORLD);  
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
    printf("epart is %d,%d, %d\n",(*epart)[0],(*epart)[1],(*epart)[2]);
    //Full METIS arrary should be avaible for every processor
    MPI_Bcast(*epart,num_elems,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(*npart,num_elems*8,MPI_INT,0,MPI_COMM_WORLD);
    //ditribute data according to METIS Partition
        
    }//finish sinle 
    }//finish metis

    MPI_Barrier(MPI_COMM_WORLD);

    if ( my_rank == 0 ) {
       printf("processor 0 npro is%d,cgup%f \n", npro,(*cgup_local)[0]);
    }
    if (my_rank==1)
        printf("processor 1 after MBI_scatter cgup_local is,%fepart is%d \n",(*cgup_local)[0],(*epart)[0]);

    return 0;
    }
