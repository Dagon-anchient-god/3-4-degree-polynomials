#pragma once
#define _USE_MATH_DEFINES

#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>

template<typename T>
class polynomial //abstract class for polynomial
{
public:
    polynomial() : coefs(0), roots(0) {}
    polynomial(int degree) : coefs(degree, 0), roots(degree, 0) {}
    virtual ~polynomial() {}

    void info()
    {
        std::cout << "~~~~~~~~~~~~~~~~" << std::endl << "Degree of the polynomial: " << this->roots.size();
        std::cout << std::endl << "The polynomial: x^" << this->coefs.size();
        for (int i = this->coefs.size() - 1; i > 0; --i)
        {
            if (coefs[i] < 0) std::cout << " - " << std::abs(this->coefs[i]);
            else std::cout << " + " << this->coefs[i];
            std::cout << "*x^" << i;
        }
        if (this->coefs[0] < 0) std::cout << " - " << std::abs(this->coefs[0]);
        else std::cout << " + " << this->coefs[0];
        std::cout << std::endl << "Roots of the polynomial: ";
        for (auto v : this->roots)
        {
            std::cout << v << ", ";
        }
        std::cout << std::endl << "~~~~~~~~~~~~~~~~" << std::endl;
    }
    std::vector<T> get_coefs()
    {
        return this->coefs;
    }
    std::vector<T> get_roots()
    {
        return this->roots;
    }
    float error_est_sum(std::vector<T> est_roots)  //error estimation
    {
        //if (est_roots.size()!=roots.size())
        //  exception
        //else
        float est = 0;
        for (int i = 0; i < this->roots.size(); ++i)
            est += abs(est_roots[i] - this->roots[i]);
        return est;
    }
    float error_est_max(std::vector<T> est_roots)
    {
        //if (est_roots.size()!=roots.size())
        //  exception
        //else
        float est = 0;
        for (int i = 0; i < this->roots.size(); ++i)
            est += abs(est_roots[i] - this->roots[i]) / std::max(est_roots[i], this->roots[i]);
        return est;
    }
protected:
    virtual void count_coefs() {}
    std::vector<T> roots;
    std::vector<float> coefs;   //coefficients of a polynomial
};

template<typename T>
class third_degree_polynomial : public polynomial<T>
    ///
    ///         x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    third_degree_polynomial() : polynomial<T>(3) {}
    third_degree_polynomial(float a, float b, float c) : polynomial<T>(3) //constructor that receives coefficients
    {
        this->coefs[2] = a; this->coefs[1] = b; this->coefs[0] = c;
    }
    third_degree_polynomial(std::vector<T> rec_roots) : polynomial<T>(3) //constructor that receives roots
    {
        for (int i = 0; i < 3; ++i) this->roots[i] = rec_roots[i];
        std::sort(this->roots.begin(), this->roots.end());
        count_coefs();
    }
    ~third_degree_polynomial() {}
protected:
    void count_coefs()
    {
        this->coefs[0] = -(this->roots[0] * this->roots[1] * this->roots[2]);
        this->coefs[1] = this->roots[0] * this->roots[1] + this->roots[0] * this->roots[2] + this->roots[1] * this->roots[2];
        this->coefs[2] = -(this->roots[0] + this->roots[1] + this->roots[2]);
    }
};

template<typename T>
class fourth_degree_polynomial : public polynomial<T>
    ///
    ///         x^4 + coefs[3]*x^3 + coefs[2]*x^2 + coefs[1] * x + coefs[0]
    ///
{
public:
    fourth_degree_polynomial() : polynomial<T>(4)
    {
    }
    fourth_degree_polynomial(float a, float b, float c, float d) : polynomial<T>(4) //constructor that receives coefficients
    {
        this->coefs[3] = a; this->coefs[2] = b; this->coefs[1] = c; this->coefs[0] = d;
    }
    fourth_degree_polynomial(std::vector<T> rec_roots) : polynomial<T>(4) //constructor that receives roots
    {
        for (int i = 0; i < 4; ++i) this->roots[i] = rec_roots[i];
        std::sort(this->roots.begin(), this->roots.end());
        count_coefs();
    }
    ~fourth_degree_polynomial() {}
protected:
    void count_coefs()
    {
        this->coefs[0] = this->roots[0] * this->roots[1] * this->roots[2] * this->roots[3];
        this->coefs[1] = -(this->roots[0] * this->roots[1] * this->roots[2] + this->roots[0] * this->roots[1] * this->roots[3] + this->roots[0] * this->roots[2] *
            this->roots[3] + this->roots[1] * this->roots[2] * this->roots[3]);
        this->coefs[2] = this->roots[0] * this->roots[1] + this->roots[0] * this->roots[2] + this->roots[0] * this->roots[3] + this->roots[1] * this->roots[2] +
            this->roots[1] * this->roots[3] + this->roots[2] * this->roots[3];
        this->coefs[3] = -(this->roots[0] + this->roots[1] + this->roots[2] + this->roots[3]);
    }
};