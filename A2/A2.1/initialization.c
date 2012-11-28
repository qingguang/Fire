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
                   int*** recv_list, int** epart, int** npart, int** objval) {
    /********** START INITIALIZATION **********/
    int i = 0;
    int my_rank, num_procs;
    //double *cgupart = (double*) calloc(sizeof(double), (*nextcf + 1)/2);
    //MPI_Init(&argc, &argv);    /// Start MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processe
    // read-in the input file by one processor
        printf("processor number is%d \n", my_rank);
    if (my_rank==0) {      
    int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                                   &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
                                   &*points, &*elems);
    
    if ( f_status != 0 ) return f_status;
    }
    MPI_Bcast (nintci,1, MPI_INT, 0, MPI_COMM_WORLD);    
    MPI_Bcast (nintcf,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (nextci,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (nextcf,1, MPI_INT, 0, MPI_COMM_WORLD);
    if (my_rank==1){
    printf("nintci is : %d\n", *nextcf);
    }
    *var = (double*) calloc(sizeof(double), (*nextcf + 1));
    *cgup = (double*) calloc(sizeof(double), (*nextcf + 1));
    *oc = (double*) calloc(sizeof(double), (*nintcf + 1));
    *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));
    
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
    //if (my_rank==0)
   // {
    for ( i = (*nintci); i <= (*nintcf); i++ )
        (*cgup)[i] = 1.0 / ((*bp)[i]);
    //End of reading and initialize by one processor
    
    //test
    int t=1;
    printf("values cgup:%f\n", (*cgup)[1]);
    printf("nintcf is:%d,%d \n", *nintci,*nintcf);
    //t=1;
    
    //Data distribution
    //classical data distribution from one processor to all another 
    //double* rebu=(double*) calloc(sizeof(double), (*nintcf + 1)/num_procs);
    //int npartpro=(*nextcf + 1)/num_procs;
    //printf("npartpro is :%d\n",t);
    } 
    //MPI_Scatter (cgup, 47312, MPI_DOUBLE, rebu, npartpro, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Bcast (&*cgup,*nextcf + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Metis Dual
    int ne = *nintcf-*nintci+1;
    int nn = *points_count;
    int ncommon= 4;
    int nparts = 2;
     //int* eptr = *element
    //int *ncommon = (int*) malloc(sizeof(int));
    //*ncommon = 4;
    //int *nparts =  (int*) malloc(sizeof(int));
   // *nparts = 6;
    int *eptr = (int*) calloc((ne + 1),sizeof(int));;
    //*eptr[0]=0;
    for ( i = (*nintci); i <= (*nintcf + 1) ; i++ ) {
        eptr[i]=i*8;
    }
    printf("numberelement and node%d,%d,%d\n",(*elems)[(*nintcf+1)*8-1], *points_count,(eptr)[*nintcf+1]);
    //int *tpwgts;
    //int* options; 
    int* options[METIS_NOPTIONS];
    options[METIS_OPTION_NUMBERING]=0;
    int metis_final = METIS_PartMeshDual(&ne,&nn,eptr, *elems, NULL, NULL, 
                                       &ncommon, &nparts, NULL,options, *objval, *epart, *npart);
    if (metis_final != METIS_OK){
         printf("Metis part Dual fails\n");
       }
     
    /*//Metis Node
    METIS_PartMeshDual(ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgts, options, objval, epart, npart);
    */
    return 0;
}

