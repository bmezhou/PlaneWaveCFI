# include "math.h"
# include "omp.h"
# include "mex.h"

#define sampleSize 50

// beamFormDAQ(int * rawData, int tTot)
// mex beamFormAngle.c COMPFLAGS="/openmp /Ot /fp:fast $COMPFLAGS"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    float fs, c, eleSpace, spaSpac;
    float alphaR;
    float eleSpac, sapSpac;
    

    // int m, n;
    float tTot;
    float * rawData, * dasData;
    float deltaZ, deltaX;

    int i, j, j1, k;
    float pi = 3.1415926535897932384626f;
    float z, x, xi, ttx, trx;
    // float delay;

    int apeSize;
    int iIntep;  
    
    ///////////
    float z0, x0, zc, xc, r;
    int m1, n1;
    //////////
    int * ape;
    float * delay0, * apeWin; 
    short * delay1;
    /////////////////////////
    fs = 40e6;
    c  = 1540.0;
    //////////////////////////
    // rawData    = (float *) mxGetPr(prhs[0]);
    //m          = (int) mxGetM(prhs[0]);
    //n          = (int) mxGetN(prhs[0]);
    
    // The delay between the transmission of the high-voltage pulse
    // and the start of the acquisition of the SonixDAQ
    
    
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
    tTot       = (float) * (mxGetPr(prhs[0]));
    alphaR = (float)*(mxGetPr(prhs[1]));
    x0     = (float)*(mxGetPr(prhs[2]));
    z0     = (float)*(mxGetPr(prhs[3]));

    //m0     = (int)*(mxGetPr(prhs[4]));     // 
    
    m1     = (int)*(mxGetPr(prhs[4]));
    n1     = (int)*(mxGetPr(prhs[5]));
    ///////////////////////////////
    // Parameter of the L14-38 Xducer.
    eleSpac = 0.3048e-3f;
    sapSpac = c / 2.0f / fs;
    
    deltaX = sapSpac * sinf(alphaR);  // negative or positive ...
    deltaZ = sapSpac * cosf(alphaR);
    
    //
    plhs[0] = mxCreateNumericMatrix(n1 * m1,  1, mxINT32_CLASS, 0);
    ape = (int *) mxGetPr(plhs[0]);    
    plhs[1] = mxCreateNumericMatrix(n1 * m1 * 65, 1,  mxSINGLE_CLASS, 0);
    apeWin = (float *)mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateNumericMatrix(n1 * m1 * 65 * sampleSize, 1,  mxSINGLE_CLASS, 0);
    delay0 = (float *)mxGetPr(plhs[2]);
    plhs[3] = mxCreateNumericMatrix(n1 * m1 * 65 * sampleSize, 1,  mxINT16_CLASS, 0);
    delay1 = (short *)mxGetPr(plhs[3]);
    
//     #pragma omp parallel for shared(dasData, rawData, n, m, m1, n1, m0, deltaZ, deltaX, pi, x0, z0, eleSpac, c, fs, tTot) \
//                              private(i, j, k, iIntep, z, x, xi, ttx, trx, delay, apeSize)
    
    for (i = 0; i < m1; i ++)      // x-direction
    {
        xc = x0 + i * 0.2e-3f;
        // zc = z0;
                       
        for (j = 0; j < n1; j ++)  // z-direction
        {
            zc = z0 + j * 0.2e-3f;
             
            iIntep = (int) ((xc - zc*tanf(alphaR)) / eleSpac + 0.5);
             // ttx = z;
             // Unit in [m]
            r = zc/cosf(alphaR);
            apeSize = (int) (r/2.8f/eleSpac);
            if (apeSize < 4)
                apeSize = 4;
            if (apeSize > 32)
                apeSize = 32;
             
            ape[i * n1 + j] = apeSize;
            
            for (k = -apeSize; k < apeSize + 1; k ++)
                apeWin[i * n1 * 65 + j * 65 + k + apeSize] = 
                               0.54f - 0.46f * cosf(pi/32*(k + apeSize));
             
            x = xc - 30 * deltaX;
            z = zc - 30 * deltaZ;
             
            for (j1 = 0; j1 < sampleSize; j1 ++)
            {   
                for (k = -apeSize; k < apeSize + 1; k ++)
                {
//                      if ((k + iIntep) < 0 || (k + iIntep) >= n)
//                          continue;

                    xi = (k + iIntep) * eleSpac;
                     // Unit in [m]
                    trx = sqrtf( (x - xi) * (x - xi) + z * z);
                     // Unit in [m]
                     // delay = ;
                     // Uint in [s]
                    delay0[i * n1 * 65 * sampleSize + j * 65 * sampleSize + j1 * 65 + k + apeSize] = 
                                 (z + trx)/c * fs + tTot;
                     
                    delay1[i * n1 * 65 * sampleSize + j * 65 * sampleSize + j1 * 65 + k + apeSize] = 
                                 k + iIntep;
                     
//                      //mexPrintf("%d\n", delay);
//                      if (*delay >= (m - 1))
//                          continue;

                    //   dasData[i * n1 + j] = j;
//                  dasData[i * n1 + j] += (0.54f - 0.46f * cosf(pi/apeSize*(k + apeSize)))  \
//                                  * ( (delay - (int) (delay))                                \
//                                   * (   rawData[(k + iIntep) * m + (int) (delay + 1) ]      \
//                                       - rawData[(k + iIntep) * m + (int) (delay)       ])   \
//                                      + rawData[(k + iIntep) * m + (int) (delay) ] );
                }
                z += deltaZ;
                x += deltaX;
            }
        }
    }   
    
    return;
}
