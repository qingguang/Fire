#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <string.h>
#include "papi.h"
#include "xread.h"
#include "xwrite.h"
#define NUM_EVENTS 4
#define PeakPer  7970.50
int main(int argc, char *argv[]) {

    if (argc < 4) {
	   printf("Usage: %s data_type(text or bin) input_file output_file\n", argv[0]);
	   return EXIT_FAILURE;
    }

    char *file_type =argv[1]; 
    char *file_in = argv[2];
    char *file_out = argv[3];
    char *str1="SU.vtk";
    char *str2="VAR.vtk";
    char *str3="CGUP.vtk";
    char *file_perf="pstats.dat"    
    printf("file_type:%s\n",file_type); 
    //strcat(str2,file_out); 
    //strcat(str3,file_out); 
    int status = 0;
    /** internal cells start and end index*/
    int nintci, nintcf;
    /** external cells start and end index. The external cells are only ghost cells. 
    They are accessed only through internal cells*/
    int nextci, nextcf;
    /** link cell-to-cell array. Stores topology information*/
    int **lcc;
    /** red-black colouring of the cells*/
    int *nboard;
    /** boundary coefficients for each volume cell */
    double *bs, *be, *bn, *bw, *bl, *bh, *bp, *su;
    /**parameter used for volmesh and reading binary input file */
    int* nodeCnt;// = (int *) calloc(sizeof(int), (nintcf + 1)); 
    int*** points;// = (int ***) calloc(3 * sizeof(int*),nintcf + 1);
    int*** elems;// = (int ***) calloc(8 * sizeof(int*),nintcf + 1);
    /**Papi parameters*/
    long long *values_i[NUM_EVENTS] = (long long *) calloc(sizeof(long long), 3);
    long long *values_c[NUM_EVENTS] = (long long *) calloc(sizeof(long long), 3);
    long long *values_o[NUM_EVENTS] = (long long *) calloc(sizeof(long long), 3);  
    double *mflops = (double *) calloc(sizeof(double), 3);///, mflops_i,mflops_c,mflops_o;
    double *L1mira = (double *) calloc(sizeof(double), 3);
    double *Lmirate = (double *) calloc(sizeof(double), 3);
    double *util = (double *) calloc(sizeof(double), 3);
    long long *et = (long long *) calloc(sizeof(long long), 3);
    long long start_cycles, start_usec,end_cycles_1, end_usec_1, end_cycles_2, end_cycles_3, end_usec_2, end_usec_3;
    /**In cluster mpp_inter L1 and L2 events can not computed at the same time, 
    so set into two groups*/
    int Events[NUM_EVENTS]={PAPI_L2_TCM,PAPI_L2_TCA,PAPI_FP_INS,PAPI_TOT_CYC};
   // int Events[NUM_EVENTS]={PAPI_L1_TCM,PAPI_L1_TCA,PAPI_FP_INS,PAPI_TOT_CYC};
    /* initialization  */
    if ( PAPI_start_counters( Events, NUM_EVENTS ) != PAPI_OK )
    printf("Fail to start PAPI counter\n");    
    start_cycles = PAPI_get_real_cyc(); // Gets the starting time in clock cycles
    start_usec = PAPI_get_real_usec(); // Gets the starting time in microseconds
    
    // read-in the input file
    int f_status;
    if (strcmp(file_type,"text") == 0) {

             f_status = read_formatted(file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
		        &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &nboard);
    } else if (strcmp(file_type,"bin") == 0) { 

           f_status = read_formatted_bin(file_in, &nintci, &nintcf, &nextci,
                      &nextcf, &lcc, &bs, &be, &bn, &bw,
                      &bl, &bh, &bp, &su,&nboard);
    } else { 
            printf ("Input file format is nor correct\n");
             return EXIT_FAILURE;
    }
    if (f_status != 0){
		
        printf("failed to initialize data!\n");
	return EXIT_FAILURE;
    }
    // allocate arrays used in gccg
    int nomax = 3;
    /** the reference residual*/
    double resref = 0.0;
    /** the ratio between the reference and the current residual*/
    double ratio;
    /** array storing residuals */
    double* resvec = (double *) calloc(sizeof(double), (nintcf + 1));
    /** the variation vector -> keeps the result in the end */
    double* var = (double *) calloc(sizeof(double), (nextcf + 1));
    /** the computation vectors */
    double* direc1 = (double *) calloc(sizeof(double), (nextcf + 1));
    double* direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    /** additional vectors */
    double* cgup = (double *) calloc(sizeof(double), (nextcf + 1));
    double* oc = (double *) calloc(sizeof(double), (nintcf + 1));
    double* cnorm = (double *) calloc(sizeof(double), (nintcf + 1));
    double* adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double* adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    double* dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double* dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    /**array used for writing measured performance */
    
    /**store volume information*/
    int nc=0;
    // initialize the reference residual
    for ( nc = nintci; nc <= nintcf; nc++) {
	  resvec[nc] = su[nc];
	  resref = resref + resvec[nc] * resvec[nc];
    }
    resref = sqrt(resref);
    if (resref < 1.0e-15){

	printf("i/o - error: residue sum less than 1.e-15 - %lf\n", resref);
	return EXIT_FAILURE;
    }

    // initialize the arrays
    for (nc = 0; nc <= 10; nc++){
	oc[nc] = 0.0;
	cnorm[nc] = 1.0;
    }

    for (nc = nintci; nc <= nintcf; nc++){
	cgup[nc] = 0.0;
	var[nc] = 0.0;
    }

    for (nc = nextci; nc <= nextcf; nc++){
	var[nc] = 0.0;
	cgup[nc] = 0.0;
	direc1[nc] = 0.0;
	bs[nc] = 0.0;
	be[nc] = 0.0;
	bn[nc] = 0.0;
	bw[nc] = 0.0;
	bl[nc] = 0.0;
	bh[nc] = 0.0;
    }

    for (nc = nintci; nc <= nintcf; nc++){
	cgup[nc] = 1.0 / bp[nc];
    }
    int if1 = 0;
    int if2 = 0;
    int iter = 1;
    int nor = 1;
    int nor1 = nor - 1;
	
    /* finished initalization */
    //read PAPI event counters 
        if ( PAPI_read_counters( values_i, NUM_EVENTS ) != PAPI_OK ){ 
    printf("fail to stop papi counter");
    }else{ 
    printf("****Input phase*********\n");
    printf("EVENTS :%lld,%lld,%lld,%lld\n" , values_i[0],values_i[1],values_i[2],values_i[3]);}
    Lmirate[0]=(double)values_i[0]/values_i[1];
    printf("cache miss rate: %4.6f\n",Lmirate[0]);
    end_usec_1 = PAPI_get_real_usec(); 
    
    mflops[0] = (double) values_i[2] / (end_usec_1-start_usec);
    printf("Mflops: %4.4f\n",mflops[0]);
    util[0]=mflops/PeakPer;
    /* start computation loop */
    while (iter < 10000){

    /* start phase 1 */
    // update the old values of direc
    for (nc = nintci; nc <= nintcf; nc++){
	direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
    }

    // compute new guess (approximation) for direc
    for (nc = nintci; nc <= nintcf; nc++){
	direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[0][nc]]
			- bw[nc] * direc1[lcc[3][nc]] - bl[nc] * direc1[lcc[4][nc]]
			- bn[nc] * direc1[lcc[2][nc]] - be[nc] * direc1[lcc[1][nc]]
			- bh[nc] * direc1[lcc[5][nc]];
    } /* end phase 1 */
	
    /*  start phase 2 */
    // execute normalization steps
    double oc1, oc2, occ;
    if (nor1 == 1){
    oc1 = 0;
    occ = 0;
    for (nc = nintci; nc <= nintcf; nc++){	
	occ = occ + adxor1[nc] * direc2[nc];
    }
    oc1 = occ / cnorm[1];
    for (nc = nintci; nc <= nintcf; nc++){
	direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
	direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
    }
    if1++;
    }else if (nor1 == 2){
	oc1 = 0;
	occ = 0;
    for (nc = nintci; nc <= nintcf; nc++){
	occ = occ + adxor1[nc] * direc2[nc];
    }
	oc1 = occ / cnorm[1];
	oc2 = 0;
	occ = 0;
    for (nc = nintci; nc <= nintcf; nc++){
	occ = occ + adxor2[nc] * direc2[nc];
    }
    oc2 = occ / cnorm[2];
    for (nc = nintci; nc <= nintcf; nc++){
	direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
	direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
    }

    if2++;
    }

    cnorm[nor] = 0;
    double omega = 0;

    // compute the new residual
    for (nc = nintci; nc <= nintcf; nc++){
	cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
	omega = omega + resvec[nc] * direc2[nc];
    }
    omega = omega / cnorm[nor];
    double resnew = 0.0;
    for (nc = nintci; nc <= nintcf; nc++){
	var[nc] = var[nc] + omega * direc1[nc];
	resvec[nc] = resvec[nc] - omega * direc2[nc];
	resnew = resnew + resvec[nc] * resvec[nc];
    }
    resnew = sqrt(resnew);
    ratio = resnew / resref;

    // exit on no improvements of residual
    if (ratio <= 1.0e-10){
	break;
    }
    iter++;

    // prepare additional arrays for the next iteration step
    if (nor == nomax){
	nor = 1;
    }else{
    if (nor == 1){
    for (nc = nintci; nc <= nintcf; nc++){
	dxor1[nc] = direc1[nc];	
	adxor1[nc] = direc2[nc];
    }
    } else if (nor == 2){
    for (nc = nintci; nc <= nintcf; nc++){
	 dxor2[nc] = direc1[nc];
  	 adxor2[nc] = direc2[nc];
    }
    }
    nor++;
    }
    nor1 = nor - 1;

    }/* end phase 2 */

    /* finished computation loop */
    end_cycles_2 = PAPI_get_real_cyc(); // Gets the ending time in clock cycles
    end_usec_2 = PAPI_get_real_usec(); // Gets the ending time in microseconds
    if ( PAPI_read_counters( values_c, NUM_EVENTS ) != PAPI_OK ){ 
    printf("fail to stop papi counter");
    }else{
    printf("****Computation phase*********\n");
    printf("EVENTS: %lld,%lld,%lld,%lld\n" , values_c[0],values_c[1],values_c[2],values_c[3]);}

    Lmirate[1]=(double)values_c[0]/values_c[1];
    printf("cache miss rate: %4.6f\n",Lmirate[1]);

    //Caculate Mflops using total floating caculation and total caculation time
    mflops[1] = (double) ( values_c[2] ) / ( end_usec_2-end_usec_1 );
    printf("Mflops: %4.6f\n",mflops_c);	
    util[1]=mflops[1]/PeakPer;
    /* write output file  */
    /**
    if ( write_result(file_in, file_out, nintci, nintcf, var, iter, ratio) != 0 )
    printf("error when trying to write to file %s\n", file_out);*/
    
    //transfer volume to mesh
    // if (strcmp(file_type,"text")==0){    
    if (vol2mesh(nintci, nintcf, lcc, &nodeCnt, &points, &elems) != 0 ){ 
        printf("error when trying to converge topology to volume");
    }   
    //write output to vtk file    
    if (write_result_vtk(str1, nintci, nintcf, nodeCnt, points, elems, su) != 0){
       printf("error when write SU to vtk file");
    }
    if (write_result_vtk(str2, nintci, nintcf, nodeCnt, points, elems, var) != 0){
       printf("error when write VAR to vtk file");
    }
    if (write_result_vtk(str3, nintci, nintcf, nodeCnt, points, elems, cgup) != 0){
       printf("error when write CGUP to vtk file");
    }
    if (write_result_dat(file_perf, values_i,values_c,*values_o,Lmirate, et, mlops, util) != 0 ){
    printf("error when write measured performance to data file");
    }
       
    if ( PAPI_stop_counters( values_o, NUM_EVENTS ) != PAPI_OK ){ 
    printf("fail to stop papi counter");
    }else{ 
    printf("****Output phase*********\n");
    printf("EVENTS : %lld,%lld,%lld,%lld\n" , values_o[0],values_o[1],values_o[2],values_o[3]);}
    
    Lmirate[2]=(double) values_o[0]/values_o[1];
    printf("cache miss rate :%4.6f\n",Lmirate[2]);
    end_cycles_3 = PAPI_get_real_cyc(); // Gets the ending time in clock cycles
    end_usec_3 = PAPI_get_real_usec(); // Gets the ending time in microseconds 
    
    mflops[2] = (float) (values_o[2])/(end_usec_3-end_usec_2);
    printf("Mflops: %4.6f\n",mflops[2]);
    util[2] = mflops[2]/PeakPer;
    et[0] = end_usec_1-start_usec;
    et[1] = end_usec_2-end_usec_1;
    et[2] = end_usec_3-end_usec_2;
    printf("Input time:%lld,Computation time:%lld, output time: %lld\n",end_usec_1-start_usec,end_usec_2-end_usec_1,end_usec_3-end_usec_2); 
    /* Free all the dynamically allocated memory */
    free(direc2); free(direc1); free(dxor2); free(dxor1); free(adxor2); free(adxor1);
    free(cnorm); free(oc); free(var); free(cgup); free(resvec); free(su); free(bp);
    free(bh); free(bl); free(bw); free(bn); free(be); free(bs);
    printf("Simulation completed successfully!\n");
    return EXIT_SUCCESS;
    }
