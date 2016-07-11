# include "math.h"
# include "omp.h"
# include "mex.h"
# include "malloc.h"

#include <emmintrin.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    float * fp, delta;
    //float alphaR;
    __m128 f;
    
    int i;
    char * c;
            
    delta = 2.0f;
    fp = (float *)malloc(sizeof(float) * 4);
    
    f = _mm_set1_ps(3.0f);
    // f = _mm_setzero_ps();
        
    _mm_storeu_ps(fp, f);
    
    memset(fp, -3.0, sizeof(float) * 4);
    
    c  = (char *) fp;
    
    printf("%d\n", sizeof(char));
    
    for (i = 0; i < 16; i++)
        printf("%d\n", *(c + i));
    
    for (i = 0; i < 4; i++)
        printf("%f\n", *(fp + i));

    free(fp);
    
    fp = NULL;
    
    return;
}