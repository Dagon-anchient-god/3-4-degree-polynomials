#include <vector>
#include <iostream>
#include <cmath>
#include <corecrt_math_defines.h>
#include <complex>
#include "polynomials.h"
using namespace std;

const float one_third = 1.0f / 3.0f;


inline float sq(const float& x) { 
    return x * x;
}

inline int eqn_quadratic(const vector<float> a, vector<complex<float>>& x, int i, int j) {
    float p = -0.5f * a[1],
        d = sq(p) - a[0];
    if (d >= 0.0) {
        d = sqrt(d);
        x[i] = p - d;
        x[j] = p + d;
        return 2;
    }
    else
    {
        complex<float> v = d;
        x[i] = p - sqrt(v);
        x[j] = p + sqrt(v);
        return 3;
    }
}

inline int gornor(const vector<float> coef, vector<complex<float>>& solution)
{
    vector<float> coef1(2);
    coef1[1] = real((solution[0])) + coef[2];
    coef1[0] = real((solution[0])) * coef1[1] + coef[1];
    return eqn_quadratic(coef1, solution, 1, 2);
}

int algorith(const vector<float> coef, vector<complex<float>>& solution)
{
    float e = coef[2] * one_third;
    float f = coef[1] - coef[2] * e;
    float g = coef[1] * e - coef[0] - 2 * (e * e * e);
    if (f == 0)
    {
        solution[0] =  cbrt(g) - e;
        return gornor(coef, solution);          
    }
    else
    {
        float h = sqrt(4*abs(f)*one_third);
        float i = 4 * g / (h * h * h);
        if (f > 0)
        {
            solution[0] = h * sinh(one_third * asinh(i)) - e;
            return gornor(coef, solution);
        }
        else
        {
            if (abs(i) > 1)
            {
                solution[0] = (h * i / abs(i)) * cosh(one_third *  acosh(abs(i))) - e;
                return gornor(coef, solution);
            }
            else
            {
               
                float j = acos(i) * one_third;
                solution[0] = h * cos(j) - e;
                solution[1] = h * cos(((float((2 * M_PI)*one_third))+j))-e;
                if (j==0 || j == 2 * M_PI ||j == M_PI * one_third)
                {
                    vector<float> coef1;
                    coef1[0] = real((solution[0])) + coef[0];
                    coef1[1] = real((solution[0])) * coef1[0] + coef[1];
                    solution[2] = -(solution[1] + coef1[0]);
                    return 1;
                }
                else
                {
                    solution[2] = h * cos(((float((2 * M_PI) * one_third)) - j) ) - e;
                    return 0;
                }
            }
        }

    }
}


int main()
{
    vector<float> answer = { 0.7f , 13 , -1.98f };
    third_degree_polynomial<float> A(answer);
    vector<float> coef = A.get_coefs();
    coef = {0.75f,-3.33f,12.7f};
    vector<complex<float>> solution(3);
    algorith(coef, solution);
    for (int i = 0; i < solution.size(); i++)
    {
        cout << solution[i] << "\n";
        //cout << answer[i] << "\n";
    }
   
}

