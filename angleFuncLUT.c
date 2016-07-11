# include "math.h"
# include "omp.h"
# include "mex.h"

// beamFormDAQ(int * rawData, int tTot)
// mex beamFormAngle.c COMPFLAGS="/openmp /Ot /fp:fast $COMPFLAGS"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    float fs, c, eleSpace, spaSpac;
    float alphaR;
    float eleSpac, sapSpac;
    

    int m, n;
    float tTot;
    float * rawData, * dasData;
    float deltaZ, deltaX;

    int i, j, k;
    float pi = 3.1415926535897932384626f;
    float z, x, xi, ttx, trx;
    float delayTemp;

    int apeSize;
    int iIntep;  
    
    ///////////
    float z0, x0;
    int m0, m1, n1;
    ///////////
    int * ape;
    float * delay;
    float * apeWin;
    /////////////////////////
    fs = 40e6;
    c  = 1540.0;
    //////////////////////////
    rawData    = (float *) mxGetPr(prhs[0]);
    m          = (int) mxGetM(prhs[0]);
    n          = (int) mxGetN(prhs[0]);
    
    // The delay between the transmission of the high-voltage pulse
    // and the start of the acquisition of the SonixDAQ
    tTot       = (float) * (mxGetPr(prhs[1]));
    
    // Note that the unit of transmit delay determined in Texo
    // between the adjacent element is 12.5 ns (or 80 MHz sampling rate)
    // refer the note in {Apr. 22, 2014}
    
    /////////////////////////
    //mexPrintf("%d\n", m);
    //mexPrintf("%d\n", n);
    //mexPrintf("%f\n", *tTot);
    
    /////////////////////////
    // Window functions
    //////////////////////////////
    

//     plhs[1] = mxCreateNumericMatrix(n, 1, mxINT32_CLASS, 0);
//     delayData = mxGetPr(plhs[1]);
    
    // apeSize = 24;
    

    
//     deltaZ = eleSpac / 4 / tan(alphaR) / 20;
//     deltaX = eleSpac / 4 / 20 * *steerLabel;
    // x = 0;

    // alphaR = 12/180 * pi;

    // Initial parameters:
    // z0, the z-axis corrdinate of the begining point
    // x0, ...
    // m0, the index of centeral position for the begining point
    // alphaR, the steered angle for receiving.
    // m1, width of the ROI
    // n1, height of the ROI
    alphaR = (float)*(mxGetPr(prhs[2]));
    x0     = (float)*(mxGetPr(prhs[3]));
    z0     = (float)*(mxGetPr(prhs[4]));

    m0     = (int)*(mxGetPr(prhs[5]));     // 
    
    m1     = (int)*(mxGetPr(prhs[6]));
    n1     = (int)*(mxGetPr(prhs[7]));
    //
    ape    = (int *) mxGetPr(prhs[8]);
    delay  = (float *) mxGetPr(prhs[9]);
    apeWin = (float *) mxGetPr(prhs[10]);
    ///////////////////////////////
    // Parameter of the L14-38 Xducer.
    eleSpac = 0.3048e-3f;
    sapSpac = eleSpac / 4.0f / sinf(fabsf(alphaR)) / 20.0f;
    
    deltaX = sapSpac * sinf(alphaR);
    deltaZ = sapSpac * cosf(alphaR);
    
    plhs[0] = mxCreateNumericMatrix(n1, m1,  mxSINGLE_CLASS, 0);
    dasData = (float *) mxGetPr(plhs[0]);
    
    #pragma omp parallel for shared(dasData, rawData, n, m, m1, n1, m0, deltaZ, deltaX, pi, x0, z0, eleSpac, c, fs, tTot, apeSize, delay, apeWin, ape) \
                             private(i, j, k, iIntep, z, x, xi, ttx, trx)
    
    for (i = 0; i < m1; i ++)      // x-direction
    {
        //x = x0 + i * eleSpac/2 - 30 * deltaX;
        //z = z0                 - 30 * deltaZ;
               
        iIntep = (int) ((m0 + i)/2);
        
        for (j = 0; j < n1; j ++)  // z-direction
        {
             // ttx = z;
             // Unit in [m]
             // apeSize = (int) (z/2.8/0.3048*1e3);
             // if (apeSize < 4)
//                  apeSize = 4;
//              if (apeSize > 32)
//                  apeSize = 32;
             
             apeSize = ape[i * n1 + j];
             
             for (k = -apeSize; k < apeSize + 1; k ++)
             {
                                 
                 if ((k + iIntep) < 0 || (k + iIntep) >= n)
                     continue;

                 //xi = (k + iIntep) * eleSpac;
                 // Unit in [m]
                 //trx = sqrtf( (x - xi) * (x - xi) + z * z);
                 // Unit in [m]
                 //delayTemp = (z + trx)/c * fs + tTot;
                 // Uint in [s]
                 delayTemp = delay[(i * n1 + j) * 65 + k + apeSize];
                 //mexPrintf("%d\n", delay);
                 if (delayTemp >= (m - 1))
                     continue;
                 
               //   dasData[i * n1 + j] = j;
//                  dasData[i * n1 + j] += (0.54f - 0.46f * cosf(pi/apeSize*(k + apeSize)))  \
//                                  * ( (delayTemp - (int) (delayTemp))                                \
//                                   * (   rawData[(k + iIntep) * m + (int) (delayTemp + 1) ]      \
//                                       - rawData[(k + iIntep) * m + (int) (delayTemp)       ])   \
//                                      + rawData[(k + iIntep) * m + (int) (delayTemp) ] );
                 dasData[i * n1 + j] += apeWin[(i * n1 + j) * 65 + k + apeSize]  \
                                 * ( (delayTemp - (int) (delayTemp))                                \
                                  * (   rawData[(k + iIntep) * m + (int) (delayTemp + 1) ]      \
                                      - rawData[(k + iIntep) * m + (int) (delayTemp)       ])   \
                                     + rawData[(k + iIntep) * m + (int) (delayTemp) ] );                 
             }
             //z += deltaZ;
             //x += deltaX;
         }
     }   
    
    return;
}
