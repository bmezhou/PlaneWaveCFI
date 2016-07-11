// autoData = sum(wallFilt(:, :, filterLen + 10: emsenbel    ) .* ...
//           conj(wallFilt(:, :, filterLen +  9: emsenbel - 1) ), 3);

# include "math.h"
# include "mex.h"
# include "malloc.h"
//  mex auto1Lag.c COMPFLAGS="/Ot /fp:fast $COMPFLAGS"
// In Matlab, use the syntax
// [autoDataR, autoDataI] = auto1Lag(xq, yq, begin, avgWin);

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
    
    int begin, avgWin, stepSize;
    
    xq       = (float *)mxGetPr(prhs[0]);
    yq       = (float *)mxGetPr(prhs[1]);
    begin    = (int) (*(mxGetPr(prhs[2])));
    avgWin   = (int) (*(mxGetPr(prhs[3])));
    stepSize = (int) (*(mxGetPr(prhs[4])));
//     mexPrintf("%d\n", begin);
    
    dimArray = mxGetDimensions(prhs[0]);
    
//     mexPrintf("%d\n", dimArray[0]);
//     mexPrintf("%d\n", dimArray[1]);
//     mexPrintf("%d\n", dimArray[2]);
    
    plhs[0] = mxCreateNumericMatrix(dimArray[0], 
                                    dimArray[1], 
                                    mxSINGLE_CLASS, 0);
    plhs[1] = mxCreateNumericMatrix(dimArray[0], 
                                    dimArray[1], 
                                    mxSINGLE_CLASS, 0);
    plhs[2] = mxCreateNumericMatrix(1, 
                                    1, 
                                    mxINT32_CLASS, 0);
//     plhs[2] = mxCreateNumericMatrix(dimArray[0], 
//                                     dimArray[1], 
//                                     mxSINGLE_CLASS, 0);
//     plhs[3] = mxCreateNumericMatrix(dimArray[0], 
//                                     dimArray[1], 
//                                     mxSINGLE_CLASS, 0);     
    
    autoDataR = (float *)mxGetPr(plhs[0]);
    autoDataI = (float *)mxGetPr(plhs[1]);
    
    jLabel = (int *)mxGetPr(plhs[2]);
        
    xTabR = (float *)malloc(sizeof(float) * dimArray[1] * dimArray[0]);
    xTabI = (float *)malloc(sizeof(float) * dimArray[1] * dimArray[0]);
// 
    memset(xTabR, 0, sizeof(float) * dimArray[1] * dimArray[0]);
    memset(xTabI, 0, sizeof(float) * dimArray[1] * dimArray[0]);    
    
    // Original
//     for (i = 0; i < dimArray[1]; i ++)
//     {
//         for (j = avgWin; j < dimArray[0] - avgWin; j ++)
//         {
//             for (ja = -avgWin; ja < avgWin + 1; ja ++)
//             {
//                 for (k = begin - 1; k < dimArray[2]; k ++) // dimArray[2]
//                 {
//                     autoDataR[i*dimArray[0] + j]  += 
//                             xq[ k    * dimArray[0] * dimArray[1] + i * dimArray[0] + j + ja]
//                           * xq[(k-1) * dimArray[0] * dimArray[1] + i * dimArray[0] + j + ja]
//                          +  yq[ k    * dimArray[0] * dimArray[1] + i * dimArray[0] + j + ja]
//                           * yq[(k-1) * dimArray[0] * dimArray[1] + i * dimArray[0] + j + ja];
//                     
//                     autoDataI[i*dimArray[0] + j]  +=
//                             yq[ k    * dimArray[0] * dimArray[1] + i * dimArray[0] + j + ja]
//                           * xq[(k-1) * dimArray[0] * dimArray[1] + i * dimArray[0] + j + ja] 
//                          -  xq[ k    * dimArray[0] * dimArray[1] + i * dimArray[0] + j + ja]
//                           * yq[(k-1) * dimArray[0] * dimArray[1] + i * dimArray[0] + j + ja];
//                 }
//             }   
//         }
//     }
    // New table method
    // Note that avgWin should be larger than 1.

//   
    for (i = 0; i < dimArray[1]; i ++)
    {
        for (j = 0; j < dimArray[0]; j++)
        {                    
            for (k = begin - 1; k < dimArray[2]; k ++) //dimArray[2]
            {
                xTabR[i*dimArray[0] + j] +=
                            xq[ k    * dimArray[0] * dimArray[1] + i * dimArray[0] + j]
                          * xq[(k-1) * dimArray[0] * dimArray[1] + i * dimArray[0] + j]
                         +  yq[ k    * dimArray[0] * dimArray[1] + i * dimArray[0] + j]
                          * yq[(k-1) * dimArray[0] * dimArray[1] + i * dimArray[0] + j];
                //
                xTabI[i*dimArray[0] + j] +=
                            yq[ k    * dimArray[0] * dimArray[1] + i * dimArray[0] + j]
                          * xq[(k-1) * dimArray[0] * dimArray[1] + i * dimArray[0] + j] 
                         -  xq[ k    * dimArray[0] * dimArray[1] + i * dimArray[0] + j]
                          * yq[(k-1) * dimArray[0] * dimArray[1] + i * dimArray[0] + j];
            }
        }
    }
    
//     printf("\n\n\n");
    
    ////////
    for (i = 0; i < dimArray[1]; i ++)
    {
        *jLabel = 0;
        for (j = avgWin; j < dimArray[0] - avgWin; j += stepSize)
        {
            for (ja = -avgWin; ja < avgWin + 1; ja ++)
            {
                autoDataR[i*dimArray[0] + *jLabel] += xTabR[i * dimArray[0] + j + ja];                       
                autoDataI[i*dimArray[0] + *jLabel] += xTabI[i * dimArray[0] + j + ja];
            }
            *jLabel += 1;
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