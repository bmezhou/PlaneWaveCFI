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
    
    float * anaDataR;
    short * anaDataI;
    
    int i;
    
    int a[10];
    
    for (i = 0; i < 10; i++)
        a[i] = i;

    
    anaDataR = (float *) a;
    anaDataI = a;
    for (i = 0; i < 10; i++)
        printf("%d\n", a + i);
    
    printf("%d\n", anaDataR + 2);
    
    printf("%d\n", anaDataI + 6);
    
    
    return;
}