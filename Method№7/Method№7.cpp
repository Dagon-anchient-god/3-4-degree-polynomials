
#include "polynomials.h"
#include <vector>
#include <iostream>
using namespace std;

inline float sq(const float& x) { // Remark #1
    return x * x;
}

inline int signum(const float& x) {
    return (0.0 < x) - (x < 0.0);
}

inline float newton1(const float* a, float x) {
    float y, y1;
    y = a[2] + x;
    y1 = 2 * y + x;
    y1 = x * y1 + a[1];
    y = (x * y + a[1]) * x + a[0];
    if (y1 != 0.0) x -= y / y1;
    return x;
}

inline int eqn_quadratic(const float* a, float* x, int i, int j) {
    float p = -0.5f * a[1],
        d = sq(p) - a[0];
    if (d >= 0.0) {
        d = sqrt(d);
        x[i] = p - d;
        x[j] = p + d;
        return 2;
    }
    return 0;
}

inline int eqn_quadratic(const float* a, const float* o, float* x,
    int i, int j) {
    float p = -0.5f * a[1],
        d = sq(p) - a[0];
    if (d >= 0.0) {
        d = sqrt(d);
        if (p < 0.0) {
            x[i] = newton1(o, p - d);
            x[j] = p + d;
        }
        else {
            x[i] = p - d;
            x[j] = newton1(o, p + d);
        }
            return 2;
    }
    return 0;
}

int eqn_cubic(const float* a, float* x) { // Remark #2
    int i_slope, i_loc;
    float w, xh, y, y1, y2, dx, b[3], c[2], d, * pb;
    const float* pa,
        prec = 1.0e-6f; // termination criterion, Remark #3
    if (a[3] == 0.0) { // a less-than-cubic problem?
        if (a[2] == 0.0) {
            if (a[1] == 0.0) // ill-defined problem
                return 0;
            else { // linear problem
                x[0] = -a[0] / a[1];
                return 1;
            }
        }
        else { // quadratic problem
            w = 1.0f / a[2];
            c[0] = a[0] * w;
            c[1] = a[1] * w;
            return eqn_quadratic(c, x, 0, 1);
        }
    }
    w = (float)1.0f / a[3]; // normalize
    pa = a;
    pb = b;
    while (pb <= b + 2) *pb++ = *pa++ * w;
    if (b[0] == 0.0) { // root at zero? Remark #4
        x[0] = 0.0;
        if (eqn_quadratic(b + 1, x, 1, 2) == 0)
            return 1;
        else { // sort results
            if (x[2] < 0.0) {
                x[0] = x[1];
                x[1] = x[2];
                x[2] = 0.0;
            }
            else if (x[1] < 0.0) {
                x[0] = x[1];
                x[1] = 0.0;
            }
            return 3;
        }
    }

    xh = -1.0f / 3.0f * b[2]; // inflection point, Remark #5
    y = b[0] + xh * (b[1] + xh * (b[2] + xh));
    if (y == 0.0) { // is inflection point a root?
        x[0] = x[1] = xh;
        c[1] = xh + b[2]; // deflation
        c[0] = c[1] * xh + b[1];
        return 1 + eqn_quadratic(c, x, 0, 2);
    }
    i_loc = (y >= 0.0);
    d = sq(b[2]) - 3 * b[1];
    if ((i_slope = signum(d)) == 1) // Laguerre-Nair-Samuelson bounds
        xh += ((i_loc) ? -2.0f / 3.0f : 2.0f / 3.0f) * sqrt(d);
    else if (i_slope == 0) { // saddle point?
        x[0] = xh - cbrt(y);
        return 1;
    }
    do { // iteration (Halley’s method)
        y = b[2] + xh;
        y1 = 2 * y + xh;
        y2 = y1 + 3 * xh;
        y1 = xh * y1 + b[1];
        y = (xh * y + b[1]) * xh + b[0];
        dx = y * y1 / (sq(y1) - 0.5f * y * y2);
        xh -= dx;
    } while (fabs(dx) > prec * fabs(xh)); // Remark #6
    x[0] = x[2] = xh;
    if (i_slope == 1) {
        c[1] = xh + b[2]; // deflation
        c[0] = c[1] * xh + b[1];
        return 1 + eqn_quadratic(c, b, x, i_loc, i_loc + 1);
    }
    return 1;
}



int main()
{
    vector<float> answer = { -1 ,3 , 4 };
    third_degree_polynomial<float> A(answer);
    vector<float> coef = A.get_coefs();
    float* a = new float[4]{ 0,0,0,1 };
    float* x = new float[3];
   for (int i =0; i < 3; i++)
    {
        a[i] = coef[i];
        cout << a[i] << "\n";
        cout << coef[i] << "\n";
    }
    eqn_cubic(a, x);
    for (int i = 0; i < 3; i++)
    {
        cout << x[i] << "\n";
    }
}

