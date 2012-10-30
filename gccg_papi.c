#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "xread.h"
#include "xwrite.h"
#include "papi.h"
#define NUM_EVENTS 4
char format_string[] = { "Real_time: %f Proc_time: %f Total flpins: %lld MFLOPS: %f\n" };
int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		printf("Usage: %s input_file output_file\n", argv[0]);
		return EXIT_FAILURE;
	}

	char *file_in = argv[1];
	char *file_out = argv[2];

	int status = 0;

	/** internal cells start and end index*/
	int nintci, nintcf;
	/** external cells start and end index. The external cells are only ghost cells. They are accessed only through internal cells*/
	int nextci, nextcf;
	/** link cell-to-cell array. Stores topology information*/
	int **lcc;
	/** red-black colouring of the cells*/
	int *nboard;

	/** boundary coefficients for each volume cell */
	double *bs, *be, *bn, *bw, *bl, *bh, *bp, *su;
    /**Papi parameters*/
    long long values[NUM_EVENTS];

    int retval;
    float real_time, proc_time, mflops;
    long long flops;
 
    int Events[NUM_EVENTS]={PAPI_L2_TCM,PAPI_L2_TCA,PAPI_FP_INS,PAPI_TOT_CYC};
    //int Events[NUM_EVENTS]={PAPI_L1_TCM,PAPI_L1_TCA,PAPI_FP_INS,PAPI_TOT_CYC};  
 	/* initialization  */
	// read-in the input file

	int f_status = read_formatted(file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
			&bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &nboard);

	if (f_status != 0)
	{
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
    /**useful value*/
    int nc=0;
	// initialize the reference residual
	for ( nc = nintci; nc <= nintcf; nc++)
	{
		resvec[nc] = su[nc];
		resref = resref + resvec[nc] * resvec[nc];
	}
	resref = sqrt(resref);
	if (resref < 1.0e-15)
	{
		printf("i/o - error: residue sum less than 1.e-15 - %lf\n", resref);
		return EXIT_FAILURE;
	}

	// initialize the arrays
	for ( nc = 0; nc <= 10; nc++)
	{
		oc[nc] = 0.0;
		cnorm[nc] = 1.0;
	}

	for (nc = nintci; nc <= nintcf; nc++)
	{
		cgup[nc] = 0.0;
		var[nc] = 0.0;
	}

	for (nc = nextci; nc <= nextcf; nc++)
	{
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

	for (nc = nintci; nc <= nintcf; nc++)
		cgup[nc] = 1.0 / bp[nc];

	int if1 = 0;
	int if2 = 0;
	int iter = 1;
	int nor = 1;
	int nor1 = nor - 1;
	/* finished initalization*/

//   if (PAPI_flops(&real_time, &proc_time, &flops, &mflops) != PAPI_OK) 
  //   printf("Fail to start flops\n ");

    if ( PAPI_start_counters( Events, NUM_EVENTS ) != PAPI_OK ) 
    printf("Fail to start PAPI counter\n");

//if (PAPI_flops(&real_time, &proc_time, &flops, &mflops) != PAPI_OK) 
  //   printf("Fail to start flops\n ");

	/* start computation loop */
	while (iter < 10000)
	{

		/* start phase 1 */

		// update the old values of direc
		for (nc = nintci; nc <= nintcf; nc++)
		{
			direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
		}

		// compute new guess (approximation) for direc
		for (nc = nintci; nc <= nintcf; nc++)
		{
			direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[0][nc]]
					- bw[nc] * direc1[lcc[3][nc]] - bl[nc] * direc1[lcc[4][nc]]
					- bn[nc] * direc1[lcc[2][nc]] - be[nc] * direc1[lcc[1][nc]]
					- bh[nc] * direc1[lcc[5][nc]];
		} /* end phase 1 */

		/*  start phase 2 */
		// execute normalization steps
		double oc1, oc2, occ;
		if (nor1 == 1)
		{
			oc1 = 0;
			occ = 0;
			for (nc = nintci; nc <= nintcf; nc++)
			{
				occ = occ + adxor1[nc] * direc2[nc];
			}
			oc1 = occ / cnorm[1];
			for (nc = nintci; nc <= nintcf; nc++)
			{
				direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
				direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
			}
			if1++;

		}
		else if (nor1 == 2)
		{
			oc1 = 0;
			occ = 0;
			for (nc = nintci; nc <= nintcf; nc++)
				occ = occ + adxor1[nc] * direc2[nc];

			oc1 = occ / cnorm[1];
			oc2 = 0;
			occ = 0;
			for (nc = nintci; nc <= nintcf; nc++)
				occ = occ + adxor2[nc] * direc2[nc];

			oc2 = occ / cnorm[2];
			for (nc = nintci; nc <= nintcf; nc++)
			{
				direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
				direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
			}

			if2++;
		}

		cnorm[nor] = 0;
		double omega = 0;

		// compute the new residual
		for (nc = nintci; nc <= nintcf; nc++)
		{
			cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
			omega = omega + resvec[nc] * direc2[nc];
		}
		omega = omega / cnorm[nor];

		double resnew = 0.0;
		for (nc = nintci; nc <= nintcf; nc++)
		{
			var[nc] = var[nc] + omega * direc1[nc];
			resvec[nc] = resvec[nc] - omega * direc2[nc];
			resnew = resnew + resvec[nc] * resvec[nc];
		}
		resnew = sqrt(resnew);
		ratio = resnew / resref;

		// exit on no improvements of residual
		if (ratio <= 1.0e-10)
			break;

		iter++;

		// prepare additional arrays for the next iteration step
		if (nor == nomax)
			nor = 1;
		else
		{
			if (nor == 1)
			{
				for (nc = nintci; nc <= nintcf; nc++)
				{
					dxor1[nc] = direc1[nc];
					adxor1[nc] = direc2[nc];
				}

			} else if (nor == 2)
			{
				for (nc = nintci; nc <= nintcf; nc++)
				{
					dxor2[nc] = direc1[nc];
					adxor2[nc] = direc2[nc];
				}
			}
			nor++;
		}
		nor1 = nor - 1;

	}/* end phase 2 */

	/* finished computation loop */
//if (PAPI_flops(&real_time, &proc_time, &flops, &mflops) != PAPI_OK) {
  //    printf("fail to read flogs\n");
  // }


   //printf(format_string, real_time, proc_time, flops, mflops);
if ( PAPI_stop_counters( values, NUM_EVENTS ) != PAPI_OK ) 
printf("fail to stop papi counter");
else 
printf("%lld,%lld,%lld,%lld,%lld\n" , values[0],values[1],values[2],values[3]);
	
/* write output file  */


	if ( write_result(file_in, file_out, nintci, nintcf, var, iter, ratio) != 0 )
		printf("error when trying to write to file %s\n", file_out);

	/* Free all the dynamically allocated memory */
	free(direc2); free(direc1); free(dxor2); free(dxor1); free(adxor2); free(adxor1);
	free(cnorm); free(oc); free(var); free(cgup); free(resvec); free(su); free(bp);
	free(bh); free(bl); free(bw); free(bn); free(be); free(bs);

	printf("Simulation completed successfully!\n");
	return EXIT_SUCCESS;
}
