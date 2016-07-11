# include "math.h"
# include "omp.h"
# include "mex.h"
// mex beamFormRect.c COMPFLAGS="/openmp /O2 /fp:fast $COMPFLAGS"
// 
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double fs, c, width, kerf, apeWid, eleSpace, spaSpac;
    double alpha, alphaR;
    double eleSpac, sapSpac;
    
    float * rawData;
    int m, n;
    double *tTot, * delayValue, * delayValueR;
    float * dasData;
    double deltaZ, deltaX;
    
    double * steerLabel;
    
    int * delayData;

    int i, j, k;
    double pi = 3.1415926535897932384626;
    double z, x, xi, ttx, trx;
    float delayTemp;

    int apeSize;
    int iIntep;
    
    //**********************//
    int *ape;
    float *delay;
    float *apeWin;
    //**********************//
    
    int i1, j1, k1, apeSize1, iIntep1, delay1;
    double z1, x1, xi1, ttx1, trx1;    
    /////////////////////////
    fs = 40e6;
    c  = 1540.0;
    // Parameter of the L14-38 Xducer.
    width = 0.2798e-3;
    kerf  = 0.0250e-3;
    eleSpac = width + kerf;
    
    apeWid = 127 * eleSpac;
    
    sapSpac = c/fs/2;
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
    // delayValueT = mxGetPr(prhs[2]);
    // alphaT = atan(12.5e-9 * *delayValueT * c/eleSpac);
    
    delayValue = mxGetPr(prhs[2]);
    // alpha = atan(12.5e-9 * *delayValue * c/eleSpac);
    
    // steerLabel = mxGetPr(prhs[3]);
    /////////////////////////
    // Window functions
    
//     for (i = 0; i < 66; i ++)
//         printf("%f\n", win[i]);
    //////////////////////////////
    
    plhs[0] = mxCreateNumericMatrix(m, (255), mxINT32_CLASS, 0);
    ape = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericMatrix(m * 255 * 65, 1, mxSINGLE_CLASS, 0);
    delay = mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericMatrix(m * 255 * 65, 1, mxSINGLE_CLASS, 0);
    apeWin = mxGetPr(plhs[2]);
    // apeSize = 24;
    deltaZ = c/fs/2; //eleSpac / 4 / 20;
    // x = 0;
    
    #pragma omp parallel for shared(dasData, rawData, alpha, eleSpac, steerLabel, n, m, deltaZ, deltaX, pi) \
                             private(i, j, k, iIntep, z, x, xi, ttx, trx, delay, apeSize)
    //////
    for (i = 0; i < 255; i ++)      // x-direction
    {
        z = 1e-3;
        x = i * eleSpac/2;
        
        // iIntep = (int) (i/2);
        iIntep = (int)( x / eleSpac + 0.5);  // Position of the central element        
        
        for (j = 0; j < m; j ++)  // z-direction
        { 
            ttx = z;
            // Unit in [m]
            apeSize = (int) (z/2.0/0.3048*1e3);
            if (apeSize < 8)
                apeSize = 8;
            
            if (apeSize > 32)
                apeSize = 32;
            
            ape[i * m + j] = apeSize;
            
            for (k = 0; k < 2 * apeSize + 1; k ++)
            {
                if ((k - apeSize + iIntep) < 0 || (k - apeSize + iIntep) >= n)
                     continue;
                
                xi = (iIntep + k - apeSize) * eleSpac;
                // Unit in [m]
                trx = sqrt( (x - xi) * (x - xi) + z * z);
                // Unit in [m]
                delayTemp = (ttx + trx)/c * fs + *tTot;
                
                delay[(i * m + j) * 65 + k] = delayTemp;
                apeWin[(i * m + j) * 65 + k] = (0.54 - 0.46 * cos(pi*(k + 0.0)/apeSize));
                // Uint in [s]     
            }
            z += deltaZ;
        }
    }   
    
    return;
}
