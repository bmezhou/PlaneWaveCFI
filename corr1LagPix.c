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
    
    float tempD0R[40], tempD0I[40],
          tempD1R[40], tempD1I[40],
          tempD2R[40], tempD2I[40];
    
    float corrX[21], maxCorr;
    
    int i, j, k, ii;
    int m, n;
    
    int inteLeng;
    
    float * xq, * yq, * vel;
    
    float * cmpR, * cmpI, * refR, * refI;
    
    float  displacement;
    int movLabel, movWin;
    float num, den0, den1;
    
    int begin;
    
    for (i = 0; i < m; i ++)
    {
        for (j = 0; j < n; j ++)
        {
            displacement = 0;
            for (k = begin; k < dimArray[1]; k ++)
            {
                cmpR = xq + k    * m * n * inteLeng \
                                 + i * n * inteLeng + j * inteLeng + 10;
                cmpI = yq + k    * m * n * inteLeng \
                                 + i * n * inteLeng + j * inteLeng + 10;                
                
                refR = xq + (k-1) * m * n * inteLeng \
                                 + i * n * inteLeng + j * inteLeng +  0;
                refI = yq + (k-1) * m * n * inteLeng \
                                 + i * n * inteLeng + j * inteLeng +  0;
                
                maxCorr = 0;
                movLabel = 0;                
                
                for (movWin = 0; movWin < 21; movWin ++)
                {
                    num = 0;
                    den0 = 0;
                    den1 = 0;                    

                    for (ii = 0; ii < 20; ii ++)
                    {
                        num +=    *(cmpR + ii) * *(refR + ii + movWin)\
                                + *(cmpI + ii) * *(refI + ii + movWin)\
                                + *(cmpI + ii) * *(refR + ii + movWin)\
                                - *(cmpR + ii) * *(refI + ii + movWin);

                        den0 +=   *(cmpR + ii) * *(cmpR + ii) \
                                + *(cmpI + ii) * *(cmpI + ii);
                        
                        den1 +=   *(refR + ii + movWin) * *(refR + ii + movWin) \
                                + *(refI + ii + movWin) * *(refI + ii + movWin);  
                    }
                    
                    corrX[movWin] = num/sqrtf(den0 * den1);
                    
                    if (maxCorr < corrX[movWin])
                    {
                        maxCorr  = corrX[movWin];
                        movLabel = movWin;
                    }
                }
                
                displacement += (corrX[movLabel + 1] - corrX[movLabel - 1]) / 2 /
                        (corrX[movLabel + 1] + corrX[movLabel - 1] - 2 * corrX[movLabel]);
            }
            
            vel[i * n + j] = (movLabel - 11 - displacement) / (dimArray[1] - begin);
        }
    }
    
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