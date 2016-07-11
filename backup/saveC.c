# include "math.h"
# include "omp.h"
# include "mex.h"

# include "stdio.h"

// beamFormDAQ(int * rawData, int tTot)
// mex beamFormAngle.c COMPFLAGS="/openmp /Ot /fp:fast $COMPFLAGS"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int m, n;
    
    float * anaDataR, * anaDataI;
    FILE * fp;
    
    anaDataR   = (float *) mxGetPr(prhs[0]);
    anaDataI   = (float *) mxGetPi(prhs[0]);
    m          = (int) mxGetM(prhs[0]);
    n          = (int) mxGetN(prhs[0]);
    
    //printf("%d\n", m);
    //printf("%d\n", n);
    
    // printf("%f\n", anaDataR[644 * 125 * 200 + 10]);
    
    //fp = fopen("C:\\DaqData\\data.dat", "w+");
    fp = fopen("D:\\data.dat", "w+");
    
    if (fp==NULL)
    {
        printf("ha\n");
        return;
    }
    else
    {
        fwrite(anaDataR, sizeof(float), m*n, fp);
        fwrite(anaDataI, sizeof(float), m*n, fp);
        fclose(fp);
        return;
    }
    
    //return;
}