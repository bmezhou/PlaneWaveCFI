# include "math.h"
# include "omp.h"
# include "mex.h"
# include "malloc.h"

#include <emmintrin.h>

#define sampleSize 30
#define filtSize 20

void filter(float * _src, float * coeff,  float * _dst, int width, int length);
// beamFormDAQ(int * rawData, int tTot)
// mex beamFormAngle.c COMPFLAGS="/openmp /Ot /fp:fast $COMPFLAGS"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    float fs, c;
    //float alphaR;
    

    int m, n;
    //float tTot;
    float * rawData, * dasDataR, * dasDataI;
    //float deltaZ, deltaX;

    int i, j, j1, k;
    float pi = 3.1415926535897932384626f;
    float delayTemp0;
    int   delayTemp1;

    int apeSize;
    
    mxArray *iniData;
    float * filtData0, * filtData1;
    
    ///////////
    //float z0, x0;
    int m1, n1;
    ///////////
    int * ape;
    float * delay0;
    short * delay1;
    float * apeWin;
    
    float * filt0, * filt1;
    
    int lenFilt;
    /////////////////////////////////
    // __m128 d4;
    
    /////////////////////////
    fs = 40e6;
    c  = 1540.0;
    //////////////////////////
 
    
    // The delay between the transmission of the high-voltage pulse
    // and the start of the acquisition of the SonixDAQ
    //tTot       = (float) * (mxGetPr(prhs[1]));
    
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
    rawData    = (float *) mxGetPr(prhs[0]);
    m          = (int) mxGetM(prhs[0]);
    n          = (int) mxGetN(prhs[0]);
    //alphaR = (float)*(mxGetPr(prhs[2]));
    //x0     = (float)*(mxGetPr(prhs[3]));
    //z0     = (float)*(mxGetPr(prhs[4]));

    //m0     = (int)*(mxGetPr(prhs[5]));     // 
    
    m1     = (int)*(mxGetPr(prhs[1]));
    n1     = (int)*(mxGetPr(prhs[2]));
    //
    ape    = (int *) mxGetPr(prhs[3]);
    apeWin = (float *) mxGetPr(prhs[4]);
    delay0 = (float *) mxGetPr(prhs[5]);
    delay1 = (short *) mxGetPr(prhs[6]);
    
    filt0 = (float *) mxGetPr(prhs[7]);
    filt1 = (float *) mxGetPr(prhs[8]);
    
    lenFilt = (int) mxGetN(prhs[8]);
    ///////////////////////////////
    // Parameter of the L14-38 Xducer.
    //eleSpac = 0.3048e-3f;
    //sapSpac = c/2.0f/fs;
    
    //deltaX = sapSpac * sinf(alphaR);
    //deltaZ = sapSpac * cosf(alphaR);
    
    plhs[0] = mxCreateNumericMatrix(n1 * m1 * sampleSize, 1,  mxSINGLE_CLASS, 0);
    dasDataR = (float *) mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericMatrix(n1 * m1 * (sampleSize - lenFilt + 1), 1,  mxSINGLE_CLASS, 0);
    dasDataI = (float *) mxGetPr(plhs[1]);
//     
    iniData   = mxCreateNumericMatrix(sampleSize, 1,  mxSINGLE_CLASS, 0);
    
    filtData0 = (float *)malloc(sizeof(float) * m1 * n1 * sampleSize);
//     filtData1
    
    // (float *)malloc(sizeof(float) * sampleSize);
//     filtData0 = (float *)malloc(sizeof(float) * sampleSize);
//     filtData1 = (float *)malloc(sizeof(float) * sampleSize);
    
    #pragma omp parallel for shared(dasDataR, dasDataI, rawData, n, m, m1, n1,  delay0, delay1, apeWin, ape, filtData0) \
                             private(i, j, k, delayTemp0, delayTemp1, apeSize)
    
    for (i = 0; i < m1; i ++)      // x-direction
    {
        //x = x0 + i * eleSpac/2 - 30 * deltaX;
        //z = z0                 - 30 * deltaZ;
               
        // iIntep = (int) ((m0 + i)/2);
        
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
            
            // memset(iniData, 0, sizeof(float) * sampleSize);
            for (j1 = 0; j1 < sampleSize; j1 ++)
            {
                
                for (k = -apeSize; k < apeSize + 1; k ++)
                {
                    delayTemp0 = 
                             delay0[i * n1 * 65 * sampleSize + j * 65 * sampleSize + j1 * 65 + k + apeSize];
                    //mexPrintf("%d\n", delay);
                    
                    if (delayTemp0 >= (m - 1))
                         continue;
                    
                    delayTemp1 = 
                            delay1[i * n1 * 65 * sampleSize + j * 65 * sampleSize + j1 * 65 + k + apeSize];
                    
                    if (delayTemp1 < 0 || delayTemp1 >= n)
                        continue;

                    dasDataR[i * n1 * sampleSize + j * sampleSize + j1] += 
                              apeWin[i * n1 * 65 + j * 65 + k + apeSize]  \
                                 * ( (delayTemp0 - (int) (delayTemp0))                         \
                                  * (  rawData[delayTemp1 * m + (int) (delayTemp0 + 1) ]      \
                                     - rawData[delayTemp1 * m + (int) (delayTemp0)     ])     \
                                   + rawData[delayTemp1 * m + (int) (delayTemp0) ] );
                }
            }
            // Analitical filtering.
            // float data;
            // for 
            
//             filter(filtData0 + i * n1 * sampleSize + j * sampleSize,
//                            filt0, dasDataI + i*filtSize*n1 + j*filtSize, sampleSize - lenFilt + 1, lenFilt);
//             filter(filtData0 + i * n1 * sampleSize + j * sampleSize,
//                            filt1, dasDataR + i*filtSize*n1 + j*filtSize, sampleSize - lenFilt + 1, lenFilt);
        }
    }
    
//     free(iniData);
    free(filtData0);
//     free(filtData1);
//     
//     iniData   = NULL;
    filtData0 = NULL;
//     filtData1 = NULL;
    

    
    return;
}


void filter(float * _src, float * coeff,  float * _dst, int width, int length)
{
        const float * kf   = coeff;
        float * src = _src;
        float * dst = _dst;
        int i = 0, k, nz = length;
        
        // float delta = 0.000001f;
        __m128 d4 = _mm_setzero_ps();
        
        float * S;
        
        __m128 s0, s1, s2, s3, 
               t0, t1, t2, t3;
        __m128 f;
        
        for(i = 0; i <= width - 16; i += 16 )
        {
            s0 = d4, s1 = d4, s2 = d4, s3 = d4;

            for( k = 0; k < nz; k++ )
            {
                f = _mm_load_ss(kf + k);
                f = _mm_shuffle_ps(f, f, 0);  // (__m128 f, __m128 f, unsigned int imm8)
                S = src + i + k;

                t0 = _mm_loadu_ps(S);
                t1 = _mm_loadu_ps(S + 4);
                s0 = _mm_add_ps(s0, _mm_mul_ps(t0, f));
                s1 = _mm_add_ps(s1, _mm_mul_ps(t1, f));

                t0 = _mm_loadu_ps(S + 8);
                t1 = _mm_loadu_ps(S + 12);
                s2 = _mm_add_ps(s2, _mm_mul_ps(t0, f));
                s3 = _mm_add_ps(s3, _mm_mul_ps(t1, f));
            }

            _mm_storeu_ps(dst + i, s0);
            _mm_storeu_ps(dst + i + 4, s1);
            _mm_storeu_ps(dst + i + 8, s2);
            _mm_storeu_ps(dst + i + 12, s3);
        }
// 
        for( ; i <= width - 4; i += 4 )
        {
            s0 = d4;

            for( k = 0; k < nz; k++ )
            {
                f = _mm_load_ss(kf + k);
                f = _mm_shuffle_ps(f, f, 0);
                t0 = _mm_loadu_ps(src + k + i);
                s0 = _mm_add_ps(s0, _mm_mul_ps(t0, f));
            }
            _mm_storeu_ps(dst + i, s0);
        }
        
        for (; i < width; i++)
        {
            for( k = 0; k < nz; k++ )
            {
                *(dst + i) += *(src + i + k) * *(kf + k); 
            }
        }

        return;
}
