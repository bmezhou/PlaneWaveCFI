// autoData = sum(wallFilt(:, :, filterLen + 10: emsenbel    ) .* ...
//           conj(wallFilt(:, :, filterLen +  9: emsenbel - 1) ), 3);

# include "math.h"
# include "mex.h"
# include "malloc.h"

# include "omp.h"
//  mex auto1Lag.c COMPFLAGS="/Ot /fp:fast $COMPFLAGS"
// In Matlab, use the syntax
// [autoDataR, autoDataI] = auto1Lag(xq, yq, begin, avgWin);

// #define inteLeng 40

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    const mwSize * dimArray;
    int m, n, l;
    float * xq, * yq;
        
    float * xTabR, * xTabI; //, * xI, * yI;
    float * autoDataR, * autoDataI;
    float * autoDataR1, * autoDataI1;
    int i, j, k, ja;
    int * jLabel;
    double temp;
    
    int inteLeng;
    
    int begin, avgWin, stepSize;
    
    m = (int) (*mxGetPr(prhs[0]));
    n = (int) (*mxGetPr(prhs[1]));
    
    xq       = (float *)mxGetPr(prhs[2]);
    yq       = (float *)mxGetPr(prhs[3]);
    begin    = (int) (*(mxGetPr(prhs[4])));
    
    inteLeng = (int)(*(mxGetPr(prhs[5])));
//     mexPrintf("%d\n", begin);
    
    dimArray = mxGetDimensions(prhs[2]);
    
//     mexPrintf("%d\n", dimArray[0]);
//     mexPrintf("%d\n", dimArray[1]);
//     mexPrintf("%d\n", dimArray[2]);
    
    plhs[0] = mxCreateNumericMatrix(n, m, 
                                    mxSINGLE_CLASS, 0);
    plhs[1] = mxCreateNumericMatrix(n, m,
                                    mxSINGLE_CLASS, 0);
//     plhs[2] = mxCreateNumericMatrix(1, 
//                                     1, 
//                                     mxINT32_CLASS, 0);
//     plhs[2] = mxCreateNumericMatrix(dimArray[0], 
//                                     dimArray[1], 
//                                     mxSINGLE_CLASS, 0);
//     plhs[3] = mxCreateNumericMatrix(dimArray[0], 
//                                     dimArray[1], 
//                                     mxSINGLE_CLASS, 0);     
    
    autoDataR = (float *)mxGetPr(plhs[0]);
    autoDataI = (float *)mxGetPr(plhs[1]);
    
//     jLabel = (int *)mxGetPr(plhs[2]);
        
    xTabR = (float *)malloc(sizeof(float) * m * n * inteLeng);
    xTabI = (float *)malloc(sizeof(float) * m * n * inteLeng);
// 
    memset(xTabR, 0, sizeof(float) * m * n * inteLeng);
    memset(xTabI, 0, sizeof(float) * m * n * inteLeng);
    
    // New table method
    // Note that avgWin should be always larger than 1.
    
    #pragma omp parallel for shared(xq, yq, xTabR, xTabI, begin, m, n, inteLeng) \
                             private(k, i)
    
    for (k = begin + 0; k < dimArray[1]; k ++)
        for (i = 0; i < m * n * inteLeng; i ++) //dimArray[2]
        {            
            xTabR[i] +=
                        xq[ k    * m * n * inteLeng + i]
                      * xq[(k-1) * m * n * inteLeng + i]
                     +  yq[ k    * m * n * inteLeng + i]
                      * yq[(k-1) * m * n * inteLeng + i];
                //
            xTabI[i] +=
                        yq[ k    * m * n * inteLeng + i]
                      * xq[(k-1) * m * n * inteLeng + i] 
                     -  xq[ k    * m * n * inteLeng + i]
                      * yq[(k-1) * m * n * inteLeng + i];
        }
    
//     printf("\n\n\n");
    
    ////////
    for (i = 0; i < m; i ++)
    {
        for (j = 0; j < n; j ++)
        {
            for (ja = 0; ja < inteLeng; ja ++)
            {
                autoDataR[i * n + j] += xTabR[i * n * inteLeng + j * inteLeng + ja];                       
                autoDataI[i * n + j] += xTabI[i * n * inteLeng + j * inteLeng + ja];
            }
        }
    }    
    
    free(xTabR);
    free(xTabI);
    
    xTabR = NULL;
    xTabI = NULL;
    
    return;
}     

// auto1D(xq, yq, 1, 1);
/*

101.104309 + 15.148970 -32.648422
-34.804688 - 0.561438 + 14.631359
-127.197601 + 0.953510 -30.636353
-79.956238 + 22.893644 -81.643326
57.468513 + 2.998265 -11.712244
77.000076 + 1.707104 -3.053586
80.222290 + 14.023602 -36.558235
 */