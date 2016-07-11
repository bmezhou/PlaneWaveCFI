void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
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