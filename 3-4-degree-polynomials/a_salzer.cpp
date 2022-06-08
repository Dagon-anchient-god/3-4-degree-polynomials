#include "polynomials.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <math.h>
#include <corecrt_math_defines.h>
#include <complex>
using namespace std;

const float one_third = 1.0f / 3.0f; /*���������, ������������ ����� 1\3 */

inline float algorith(const float b, const float c, const float d)/*������� ��������, ���������� ����� float, ��������� ��� ��������� float, b ��� ���������� ��� x^2,c ��� x � d - ��� ��������� ���� */
{/*������� algorith ����� ��� ���������� �������������� ����� � ��������� 3 �������
 ���� �� ���������, ������� ���, ������� ��� ������ 0, ��� 0, ��� ���������*/

    float e = b * one_third;
    float f = c - b * e;/*�� �������� ������ ���������� ������� ������� ����� ������������ ������ � ������� �� ����������*/
    float g = e * (c - 2 * e * e) - d;
    /*���������� e, g, h, j, ����� ��� ���������� ���������� ������ � ������ ������� */
    if (f == 0)
    {
        return cbrt(g) - e;/*����������� ������������� ������������� �����*/
    }
    else
    {
        float h = sqrt(4 * abs(f) * one_third);
        float i = 4 * g / (h * h * h);
        if (f > 0)
        {
            return h * sinh(one_third * asinh(i)) - e;/*����������� ������������� ������������� �����*/
        }
        else
        {
            if (abs(i) > 1)
            {
                return (h * signbit(i) * cosh(one_third * acosh(abs(i))) - e);/*����������� ������������� ������������� �����*/
            }
            else/*���� ���������� f ������ ���� � ������ ���������� i ������ 1, �� � ������� ���������� 3 ������������ �����*/
            {

                float j = acos(i) * one_third;
                if ((h * cos(j) - e) <= 0)
                    return  h * cos(j) - e;
                else if ((h * cos(((float((2 * M_PI) * one_third)) + j)) - e) <= 0)
                    return h * cos(((float((2 * M_PI) * one_third)) + j)) - e;
                else
                    return h * cos(((float((2 * M_PI) * one_third)) - j)) - e;

            }
        }

    }
}

void A_Salzer(const float z, const float x, const float y, const float v, vector<complex<float>>& solution)/*������� �_������ ������ �� ����������, ��������� 4 ���������, ������������ ����������� ����������, � ������ solution, ���� ����� �������� ������*/
{
    //z-���������� ��� x^3
    //x-���������� ��� x^2
    //y-���������� ��� x
    //v-��������� ����
    float z2 = z * z;//����������, ��������� ����� �� ��������� ������ ��������� ��������� ���
    float b = -x, c = z * y - 4 * v,
        d = v * (4 * x - z2) - y * y;
    /*���������� b, c, d ��� ����������� ������ ���������� 3 �������, ������������ ������ �������� ����� �����*/
    float x1 = algorith(b, c, d);/*������ ���������� 3 �������*/
    complex<float> m = 0.25f * z2 - x + x1;
    m = sqrt(m);
    complex<float> n;
    if (real(m) == 0 && imag(m) == 0)/* ���� ���������� m ����� ����, �� n ����������� �� ������ �������*/
    {
        n = sqrt(0.25f * x1 * x1 - v);
    }
    else
    {
        n = (z * x1 - 2 * y) * 0.25f / m;
    }
    if (imag(m) == 0)/* ���� ���������� m ������������, �� ����� ���������� ������������ ����� ���������*/
    {
        complex <float> al = 0.5f * z2 - x1 - x, be = 4 * real(n) - z * real(m);
        complex <float> gm = sqrt(al + be), dl = sqrt(al - be);
        /*���������� ����� ������ ����������*/
        complex<float> temp = (-0.5f * z + m) * 0.5f;//��������� ���������� ��� ����� ���������� ����������
        solution[0] = temp + gm * 0.5f;
        solution[2] = temp - gm * 0.5f;
        temp = temp - m;//��������� �������� ��������� ����������
        solution[1] = temp + dl * 0.5f;
        solution[3] = temp - dl * 0.5f;
    }
    else/* ���� m �����������, �� ������*/
    {
        //���� m ���������� �����, �� ��� ����� �������� �����������
        m = complex <float>(imag(m), -real(m));
        n = complex <float>(imag(n), -real(n));
        complex <float> al = 0.5f * z2 - x1 - x, be = 4.0f * n - z * m;
        complex <float> p = sqrt(al * al + be * be), gm = sqrt(0.5f * (al + p));
        complex <float> dl;
        if (real(gm) == 0 && imag(gm) == 0)/* ���� ���������� gm ����� 0, �� ���������� dl ����������� �� �������, � ���������� al ������ ���� �������� */
        {
            al = (-real(al), -imag(al));
            dl = sqrt(al);
        }
        else
        {
            dl = 0.5f * be / gm;
        }
        complex <float> t = 0.5f * (m + dl);
        /*���������� ����� ������ ����������*/
        float temp = -0.25f * z - imag(t);//��������� ���������� ��� ����� ���������� ����������
        solution[0] = complex <float>(temp + real(gm) * 0.5f, (real(t) + imag(gm) * 0.5f));
        solution[1] = conj(solution[0]);//���������� ���������� ������������ �����
        t = 0.5f * (m - dl);
        solution[2] = complex <float>(temp - real(gm) * 0.5f, real(t) - imag(gm) * 0.5f);
        solution[3] = conj(solution[2]);//���������� ���������� ������������ �����

    }
}