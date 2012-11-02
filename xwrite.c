#include <stdio.h>
#include <stdlib.h>
#include "xwrite.h"

/**
 * Export the results to a file
 *
 * @param inFileName
 * @param outFileName
 * @param NINTCI
 * @param NINTCF
 * @param VAR
 * @param ITER
 * @param RATIO
 * @return
 */
int write_result(char *inFileName, char *outFileName, int NINTCI, int NINTCF,
		double *VAR, int ITER, double RATIO)
{
	double *IPOINT = (double *) malloc((NINTCF + 1) * sizeof(double));
	int I1, I2, I3, I4, I5;

	int NC;
	for (NC = NINTCI; NC <= NINTCF; NC++)
		IPOINT[NC] = NC;

	if (NINTCF <= 1)
	{
		printf("Error: NINTCF <= 1\n");
		return -1;
	}

	I1 = NINTCF + 1;

	while (I1 != 0)
	{
		I1 = I1 / 2;
		I2 = NINTCF + 1 - I1;
		I4 = 1;

		do
		{
			I3 = I4;
			do
			{
				I5 = I3 + I1;
				if (VAR[I3] <= VAR[I5])
					break;

				double 	ZDUM = VAR[I3],
								IDUM  = IPOINT[I3];

				VAR[I3] = VAR[I5];
				IPOINT[I3] = IPOINT[I5];
				VAR[I5] = ZDUM;
				IPOINT[I5] = IDUM;
				I3 = I3 - I1;
			} while (I3 >= 1);
			I4++;
		} while (I4 < I2);
	}

	FILE *fp = fopen(outFileName, "w");
	if (fp == NULL)
	{
		printf("Error opening file %s for writing\n", outFileName);
		return -1;
	}

	printf("Writing data to %s: ", outFileName);
	fprintf(fp, "               ----------------------------------------\n");
	fprintf(fp, "               - AVL -  Linear Equation Solver - GCCG -\n");
	fprintf(fp, "               ----------------------------------------\n");
	fprintf(fp, "                     diagonal scaling\n");
	fprintf(fp, "\n\n");
	fprintf(fp, "     Input File:  %s\n", inFileName);
	fprintf(fp, "     ===========\n\n");
	fprintf(fp, "     Output File:  %s\n", outFileName);
	fprintf(fp, "     ============\n\n");
	fprintf(fp, "     No. of Active Cells:  %d\n", NINTCF);
	fprintf(fp, "     ====================\n\n");
	fprintf(fp, "     Iteration Count - Residual Ratio\n");
	fprintf(fp, "     ================================\n\n");
	fprintf(fp, "%d %e\n", ITER, RATIO);
	fprintf(fp, "     Addresses Solution (Minima)           Addresses Solution (Maxima)\n");
	fprintf(fp, "     ===========================           ===========================\n\n");

	int N;
	for (N = 1; N <= 10; N++)
		fprintf(fp, "%lf %lf %lf %lf\n", IPOINT[N], VAR[N],	IPOINT[NINTCF - N + 1], VAR[NINTCF - N + 1]);

	fclose(fp);
	printf("done!\n");

	free(IPOINT);
	return 0;
}

int write_result_vtk(char *outFileName, int startIdx, int endIdx, int nodeCnt, int **points, int **elems, double *vector)
{
	int i,j;
	int cellCnt = endIdx - startIdx + 1;
	FILE *fp = fopen(outFileName, "w");
	if(fp == NULL)
	{
		printf("Error opening file %s for writing\n", outFileName);
		return -1;
	}

	fprintf(fp, "# vtk DataFile Version 3.0\n"); // file version and identifier
	fprintf(fp,"vtk output\n");					 // header
	fprintf(fp,"ASCII\n");						 // file format
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");	 // dataset structure

	fprintf(fp,"POINTS %d float\n",nodeCnt);

	//printf("nodeCnt = %d, startInd= %d, endInd = %d, cellCnt=%d\n",nodeCnt, startIdx, endIdx, cellCnt);
	for (i = 1; i<=nodeCnt/2; i++){
		fprintf(fp,"%d %d %d ",points[0][2*i],points[1][2*i],points[2][2*i]);
		fprintf(fp,"%d %d %d\n",points[0][2*i+1],points[1][2*i+1],points[2][2*i+1]);
	}

	fprintf(fp,"\nCELLS %d %d\n",cellCnt, 9*cellCnt);
	for (i = startIdx; i<= endIdx; i++){
		fprintf(fp,"8 %d %d %d %d %d %d %d %d\n",elems[0][i],elems[1][i],elems[2][i],
			elems[3][i],elems[4][i],elems[5][i],elems[6][i],elems[7][i]);
	}

	fprintf(fp,"\nCELL_TYPES %d\n",cellCnt);
	for (i = startIdx; i<= endIdx; i++){
		fprintf(fp,"11\n");//cell type vtk voxel
	}

	fprintf(fp,"\nCELL_DATA %d\n",cellCnt);
	fprintf(fp,"SCALARS scalars float\nLOOKUP_TABLE default\n");

	for (i = startIdx; i<= endIdx; i++){
		fprintf(fp,"%f\n ",vector[i]);
	}

	fclose(fp);
	return 0;
}

int write_result_dat(char *outFileName, long long *values_i,long long *values_c,long long *values_o, 
                     double *Lmirate, long long *et, double *mlops, double *util)
{
        
        FILE *fp = fopen(outFileName, "w");
        if(fp == NULL)
        {
                printf("Error opening file %s for writing\n", outFileName);
                return -1;
        }

       
        fprintf(fp,"Measured Performance\n");                            // header
        fprintf(fp,"ASCII\n");                    // file format
        fprintf(fp,"\n");       // 
 fprintf(fp,"*****Input phase*****\n",); 
 fprintf(fp,"INPUT_PAPI_L2_TCM=%lld\n",values_i[0]); 
 fprintf(fp,"INPUT_PAPI_L2_TCA=%lld\n",values_i[1]); 
 fprintf(fp,"INPUT_L2_MissRate=%4.4f\n",Lmira[0]); 
 fprintf(fp,"INPUT_PAPI_FP_INS=%lld\n",values_i[2]); 
 fprintf(fp,"INPUT_PAPI_TOT_CYC=%lld\n",values_i[3]);
 fprintf(fp,"INPUT_Execution_Time=%lld\n",et[0]);i
 fprintf(fp,"Mflop=%f\n",mlops[0]); 
 fprintf(fp,"INPUT_Util=%f\n",util[0]);
    fprintf(fp,"\n");
 fprintf(fp,"*****Computation phase*****\n",);
 fprintf(fp,"COM_PAPI_L2_TCM=%lld\n",values_c[0]);
 fprintf(fp,"COM_PAPI_L2_TCA=%lld\n",values_c[1]);
 fprintf(fp,"COM_L2_MissRate=%f\n",Lmira[1]);
 fprintf(fp,"COM_PAPI_FP_INS=%lld\n",values_c[2]);
 fprintf(fp,"COM_PAPI_TOT_CYC=%lld\n",values_c[3]);
 fprintf(fp,"COM_Execution_Time=%lld\n",et[1]);
  fprintf(fp,"Mflop=%f\n",mlops[1]);
 fprintf(fp,"COM_Util=%f\n",util[1]);

   fprintf(fp,"\n");
 fprintf(fp,"*****Output phase*****\n",);
 fprintf(fp,"OUT_PAPI_L2_TCM=%lld\n",values_o[0]);
 fprintf(fp,"OUT_PAPI_L2_TCA=%lld\n",values_o[1]);
 fprintf(fp,"OUT_L2_MissRate=%f\n",Lmira[2]);
 fprintf(fp,"OUT_PAPI_FP_INS=%lld\n",values_o[2]);
 fprintf(fp,"OUT_PAPI_TOT_CYC=%lld\n",values_o[3]);
 fprintf(fp,"OUT_Execution_Time=%lld\n",et[2]);
 fprintf(fp,"Mflop=%f\n",mlops[2]); 
 fprintf(fp,"OUT_Util=%f\n",util[2]);
  

        fclose(fp);
        return 0;
}

