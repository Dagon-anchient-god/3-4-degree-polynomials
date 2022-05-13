#include "polynomials.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <corecrt_math_defines.h>
#include <complex>
using namespace std;

const float one_third = 1.0f / 3.0f;

inline float algorith(const float b, const float c, const float d)
{
    float e = b * one_third;
    float f = c - b * e;
    float g = c * e - d - 2 * (e * e * e);
    if (f == 0)
    {
        return cbrt(g) - e;
    }
    else
    {
        float h = sqrt(4 * abs(f) * one_third);
        float i = 4 * g / (h * h * h);
        if (f > 0)
        {
            return h * sinh(one_third * asinh(i)) - e;
        }
        else
        {
            if (abs(i) > 1)
            {
                return (h * i / abs(i)) * cosh(one_third * acosh(abs(i))) - e;
            }
            else
            {

                float j = acos(i) * one_third;
                if ((h * cos(j) - e) < 0)
                    return  h * cos(j) - e;
                else if (h * cos(((float((2 * M_PI) * one_third)) + j)) - e)
                    return h * cos(((float((2 * M_PI) * one_third)) + j)) - e;
                else
                    return h * cos(((float((2 * M_PI) * one_third)) - j)) - e;
                
            }
        }

    }
}

void A_Salzer(const vector<float> coef, vector<complex<float>>& solution)
{
    float b = -coef[2],c= coef[3]* coef[1]-4* coef[0],
        d= coef[0]*(4* coef[2]- coef[3]* coef[3])- coef[1]* coef[1];
    float x1 = algorith(b,  c,  d);
    float f = 0.25f * coef[3] * coef[3] - coef[2] + x1;
    complex<float> m = 0.25f * coef[3] * coef[3] - coef[2] + x1;
    m =  sqrt(m);
    complex<float> n;
    if (real(m) == 0 && imag(m) == 0)
    {
        n = sqrt(0.25f *x1*x1-coef[0]);
    }
    else
    {
        n = (coef[3] * x1 - 2 * coef[1]) * 0.25f / m;
    }
    if (imag(m) == 0 )
    {
        complex <float> al = 0.5f * coef[3] * coef[3] - x1 - coef[2], be = 4 * real(n) - coef[3] * real(m);
        complex <float> gm = sqrt(al + be), dl = sqrt(al - be);
        solution[0] = (-0.5f * coef[3] + m+ gm) * 0.5f;
        solution[2] = (-0.5f * coef[3] + m - gm) * 0.5f;
        solution[1] = (-0.5f * coef[3] - m + dl) * 0.5f;
        solution[3] = (-0.5f * coef[3] - m - dl) * 0.5f;
    }
    else
    {
        m = complex <float>(imag(m), -real(m));
        n = complex <float>(imag(n), -real(n));
        complex <float> al = 0.5f * coef[3] * coef[3] - x1 - coef[2], be = 4.0f * n - coef[3] * m;
        complex <float> p = sqrt(al * al + be * be), gm = sqrt(0.5f * (al + p));
        complex <float> dl;
        if (real(gm) == 0 && imag(gm) == 0)
        {
            al = (-real(al), -imag(al));
            dl = sqrt(al);
        }
        else
        {
            dl = 0.5f * be / gm;
        }
        complex <float> t = 0.5f * (m + dl);
        solution[0] = complex <float>(0.5f * (-0.5f * coef[3] + real(gm)) - imag(t), (real(t) + imag(gm) * 0.5f));
        solution[1] = conj(solution[0]);
        t = 0.5f * (m - dl);
        solution[2] = complex <float>(0.5f * (-0.5f * coef[3] - real(gm)) - imag(t), real(t) - imag(gm) * 0.5f);
        solution[3] = conj(solution[2]);
       
        cout<<'\n';
    }
}

int main()
{
    vector<float> coef = { -1.5f,0,2.5f,0 };
    vector<complex<float>> solution(4);
    A_Salzer(coef, solution);
    for (int i = 0; i < solution.size(); i++)
    {
        cout << solution[i] << "\n";
        //cout << answer[i] << "\n";
    }
}

