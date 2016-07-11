# include "math.h"
# include "omp.h"
# include "mex.h"

// beamFormDAQ(int * rawData, int tTot)

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double fs, c, eleSpace, spaSpac;
    double alphaT, alphaR;
    double eleSpac, sapSpac;
    
    float * rawData;
    int m, n;
    double *tTot, * delayValueT, * delayValueR;
    float * dasData;
    double deltaZ, deltaX;
    
    double * steerLabel;
    
    int * delayData;

    int i, j, k;
    float pi = 3.1415926535897932384626;
    double z, x, xi, ttx, trx;
    float delay;

    int apeSize;
    int iIntep;
    
    int i1, j1, k1, apeSize1, iIntep1, delay1;
    double z1, x1, xi1, ttx1, trx1;    
    
    ///////////
    double z0, x0;
    int m0, m1, n1;
    /////////////////////////
    fs = 40e6;
    c  = 1540.0;
    //////////////////////////
    rawData    = mxGetPr(prhs[0]);
    m          = mxGetM(prhs[0]);
    n          = mxGetN(prhs[0]);
    
    // The delay between the transmission of the high-voltage pulse
    // and the start of the acquisition of the SonixDAQ
    tTot       = mxGetPr(prhs[1]);
    
    
    
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
    alphaR = *(mxGetPr(prhs[2]));    
    x0     = *(mxGetPr(prhs[3]));
    z0     = *(mxGetPr(prhs[4]));

    m0     = (int)*(mxGetPr(prhs[5]));           // 
    
    m1     = (int)*(mxGetPr(prhs[6]));
    n1     = (int)*(mxGetPr(prhs[7]));
    ///////////////////////////////
    // Parameter of the L14-38 Xducer.
    eleSpac = 0.3048e-3;
    sapSpac = eleSpac/4/sin(fabs(alphaR))/20;
    
    deltaX = sapSpac * sin(alphaR);
    deltaZ = sapSpac * cos(alphaR);
    
    plhs[0] = mxCreateNumericMatrix(n1, m1,  mxSINGLE_CLASS, 0);
    dasData = mxGetPr(plhs[0]);
    
    #pragma omp parallel for shared(dasData, rawData, n, m, m1, n1, m0, deltaZ, deltaX, pi, x0, z0, eleSpac, c, fs, tTot) \
                             private(i, j, k, iIntep, z, x, xi, ttx, trx, delay, apeSize)
    
    for (i = 0; i < m1; i ++)      // x-direction
    {
        x = x0 + i * eleSpac/2 - 30 * deltaX;
        z = z0                 - 30 * deltaZ;
               
        iIntep = (int) ((m0 + i)/2);
        
        for (j = 0; j < n1; j ++)  // z-direction
        {
             ttx = z;
             // Unit in [m]
             apeSize = (int) (z/2.8/0.30480*1e3);
             if (apeSize < 4)
                 apeSize = 4;
             if (apeSize > 32)
                 apeSize = 32;
            
             for (k = -apeSize; k < apeSize + 1; k ++)
             {
                                 
                 if ((k + iIntep) < 0 || (k + iIntep) >= n)
                     continue;

                 xi = (k + iIntep) * eleSpac;
                 // Unit in [m]
                 trx = sqrt( (x - xi) * (x - xi) + z * z);
                 // Unit in [m]
                 delay = (ttx + trx)/c * fs + *tTot;
                 // Uint in [s]
                
                 //mexPrintf("%d\n", delay);
                 if (delay >= (m - 1))
                     continue;
                 
               //   dasData[i * n1 + j] = j;
                 dasData[i * n1 + j] += (0.54f - 0.46f * cosf(pi*(k + apeSize + 0.0f)/apeSize))     \
                                 * ( (delay - (int) (delay))                                  \
                                  * (   rawData[(k + iIntep) * m + (int) (delay + 1) ]  \
                                      - rawData[(k + iIntep) * m + (int) (delay)       ]) \
                                     + rawData[(k + iIntep) * m + (int) (delay) ] );
             }
             z += deltaZ;
             x += deltaX;
         }
     }   
    
    return;
}
