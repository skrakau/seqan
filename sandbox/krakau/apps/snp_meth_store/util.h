/*==========================================================================

Functors and functions for likelihood function evaluation

==========================================================================*/

#ifndef SANDBOX_KRAKAU_APPS_SNP_METH_STORE_UTIL_H_
#define SANDBOX_KRAKAU_APPS_SNP_METH_STORE_UTIL_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

//#include <boost/math/distributions.hpp>
//#include <boost/math/tools/tuple.hpp>
//#include <boost/math/tools/roots.hpp>

//#include "bs_one_calling.h"

using namespace std;
using namespace seqan;

// Prob. function
struct Naive_;
typedef Tag<Naive_> Naive;

struct Polynomial_;
typedef Tag<Polynomial_> Polynomial;

struct LogFunction_;
typedef Tag<LogFunction_> LogFunction;


// Optimize
struct Sampling_;
typedef Tag<Sampling_> Sampling;

struct Newton_;
typedef Tag<Newton_> Newton;

// Polynomial evaluation
struct Horner_;
typedef Tag<Horner_> Horner;

struct HornerScaling_;
typedef Tag<HornerScaling_> HornerScaling;

struct Classic_;
typedef Tag<Classic_> Classic;

struct LogProbValue
{
    LogProb<long double> c;
    bool sign;
    bool zero;

    LogProbValue()
    {
        c = 1.0;
        sign = true;
        zero = false;
    }
};

// Functor for extended Horner function
// Returns values for f'(x) and f''(x)
template <typename TValue>
struct Fct_12N
{   
    
    Fct_12N(String<TValue> const& coeffs) : coeffs(coeffs)
    { // Constructor 
    }
    boost::math::tuple<TValue, TValue> operator()(TValue const&beta)
    { // z is estimate so far.
        //std::cout << "start calculate horner " << std::endl;

        /*
        TValue f_0 = back(coeffs);
        TValue f_1 = 0.0;
        TValue f_2 = 0.0;

        for (int i = length(coeffs)-2; i >= 0; --i) 
        {
            f_2 = f_1 + beta*f_2;
            f_1 = f_0 + beta*f_1;
            f_0 = coeffs[i] + beta*f_0;
        }

        f_2 = 2*f_2;

        instead:
        Compute f'(x) independent of f(x)
        Makes a big difference for error bounds
        at positions where f'(x) is small but f(x) large
        error cause newton to jump away near root, if f'(x) slightly wrong
        */
        
        /* only f'(independent)
        TValue f_1 = back(coeffs)*(length(coeffs)-1);
        TValue f_2 = 0.0;

        for (int i = length(coeffs)-2; i >= 1; --i)
        {
            f_2 = f_1 + beta*f_2;
            f_1 = coeffs[i]*i + beta*f_1;
        }
        */

        // f' and f'' independent computed 
        TValue f_1 = back(coeffs)*(length(coeffs)-1);
        TValue f_2 = back(coeffs)*(length(coeffs)-1)*(length(coeffs)-2);

        for (int i = length(coeffs)-2; i >= 2; --i)
        {
            f_2 = coeffs[i]*i*(i-1) + beta*f_2;
            f_1 = coeffs[i]*i + beta*f_1;
        }
        f_1 = coeffs[1]*1 + beta*f_1;
        
        //std::cout << std::setprecision (25) << "Beta NSpace: " << beta << "  f': " << f_1 << "  f'': " << f_2 << std::endl;
        return boost::math::make_tuple(f_1, f_2);
    }
private:
    String<TValue> coeffs;
};

// NSpace with error bound
template <typename TValue>
struct Fct_12eBN
{
    Fct_12eBN(String<TValue> const& coeffs) : coeffs(coeffs)
    { // Constructor 
    }
    boost::math::tuple<TValue, TValue> operator()(bool &eBViolated, TValue const&beta)
    { // z is estimate so far.
        //std::cout << "start calculate horner " << std::endl;
        TValue f_1 = back(coeffs)*(length(coeffs)-1);
        TValue f_2 = back(coeffs)*(length(coeffs)-1)*(length(coeffs)-2);

        TValue eB_1 = 0.0;
        TValue eB_2 = 0.0;

        for (int i = length(coeffs)-2; i >= 2; --i)
        {
            // 1
            eB_2 =  abs(coeffs[i]*i*(i-1)) + abs(beta*f_2) + beta*eB_2;
            f_2 = coeffs[i]*i*(i-1) + beta*f_2;
            eB_2 += abs(f_2);
            // 0
            eB_1 = abs(coeffs[i]*i) + abs(beta*f_1) + beta*eB_1; 
            f_1 = coeffs[i]*i + beta*f_1;
            eB_1 += abs(f_1);
            //TODO: since f' independent of f: maybe simple eB can be used again?
        } 
        eB_1 = abs(coeffs[1]*1) + abs(beta*f_1) + beta*eB_1;
        f_1 = coeffs[1]*1 + beta*f_1;
        eB_1 += abs(f_1);


        eB_1 *= 6*pow(10, -18);
        eB_2 *= 6*pow(10, -18);

        //std::cout << "beta: " << beta <<  "  f_0: " << f_0 << " eB_0: " << eB_0 << "   f_1: " << f_1 << " eB_1: " << eB_1 << "   f_2: " << f_2 << " eB_2: " << eB_2 << std::endl;

        if (f_1 <= eB_1 || f_2 <= eB_2) eBViolated = true;
        else eBViolated = false;
 
        //  std::cout << "Beta NSpace: " << z << "  f': " << f_1 << "  f'': " << f_2 << std::endl;
        return boost::math::make_tuple(f_1, f_2);
    }
private:
    String<TValue> coeffs;
};

// NSpace with error bound, only for f (sampling horner)
template <typename TValue>
struct Fct_0eBN
{
    Fct_0eBN(String<TValue> const& coeffs) : coeffs(coeffs)
    { // Constructor 
    }
    TValue operator()(bool &eBViolated, TValue const&beta)
    { // z is estimate so far.
        //std::cout << "start calculate horner " << std::endl;
        TValue f_0 = back(coeffs);
        TValue eB_0 = 0.0;

        for (int i = length(coeffs)-2; i >= 0; --i)
        {
            eB_0 = abs(beta*f_0) + beta*eB_0;
            f_0 = coeffs[i] + beta*f_0;
            eB_0 += abs(f_0);
        }
        eB_0 *= 6*pow(10, -18);


        if (f_0 <= eB_0) eBViolated = true;
        else eBViolated = false;

        return f_0;
    }
private:
    String<TValue> coeffs;
};

// LogSpace
// Functor for Horner function for Sampling
// Returns values for f(x)
template <typename TValue>
struct Fct_0eBL
{
    Fct_0eBL(String<LogProbValue> const& coeffs) : coeffs(coeffs)     
    { // Constructor 
    }
    LogProbValue operator()(bool &eBViolated, TValue const & beta)
    { // beta is estimate so far.
        LogProbValue f;
        LogProb<TValue> eB;

        if (back(coeffs).zero)
        {
            f.zero = true;
        }
        else
        {
            f.zero = false;
            f.c = back(coeffs).c;
            f.sign = back(coeffs).sign;
            eB = f.c/2.0;
        }
        for (int i = length(coeffs)-2; i >= 0; --i) 
        {
            // Compute f
            if (coeffs[i].zero && f.zero)
                f.zero = true;
            else if (!coeffs[i].zero && f.zero)
            {
                f.zero = false;
                f.c = coeffs[i].c;        
                f.sign = coeffs[i].sign;
            }
            else if (coeffs[i].zero && !f.zero)
            {
                f.zero = false;
                f.c = f.c * beta;       // non-LogProb value must be on right hand side!
                f.sign = f.sign;        // nothing to do, stays the same
            }
            else if (coeffs[i].sign && f.sign)
            {
                f.zero = false;
                f.c = coeffs[i].c + f.c*beta;
                f.sign = true;        // nothing to do, stays the same
            }
            else if (!coeffs[i].sign && !f.sign)
            {
                f.zero = false;
                f.c = coeffs[i].c + f.c*beta;
                f.sign = false; 
            }
            else if (coeffs[i].sign && !f.sign)
            {
                f.zero = false;
                if (coeffs[i].c < f.c*beta)
                {
                    f.c = f.c*beta - coeffs[i].c;
                    f.sign = false;
                }
                else  
                {
                    f.c = coeffs[i].c - f.c*beta;
                    f.sign = true;
                }
            }
            else if (!coeffs[i].sign && f.sign)
            {
                f.zero = false;
                if (coeffs[i].c < f.c*beta)
                {
                    f.c = f.c*beta - coeffs[i].c;
                    f.sign = true;
                }
                else
                {
                    f.c = coeffs[i].c - f.c*beta;
                    f.sign = false;
                }
            }
            eB = eB*beta + f.c;

        } 
        eB = (eB*2.0 - f.c)*(6*pow(10, -18));    // TODO u ?

        //SEQAN_ASSERT_EQ(f.sign, true);
        SEQAN_ASSERT_EQ(f.zero, false);

        if (!f.sign || f.c <= eB) eBViolated = true;
        else eBViolated = false;
        return f;
    }
private:
    String<LogProbValue> coeffs;
};


// LogSpace
// Functor for Horner function
// Returns values for f'(x) and f''(x)
template <typename TValue>
struct Fct_12L
{
    Fct_12L(String<LogProbValue> const& coeffs) : coeffs(coeffs)     
    { // Constructor 
    }
    boost::math::tuple<TValue, TValue> operator()(TValue const & beta)
    { // beta is estimate so far.
        // TODO: name f, f_1 -> f_1, f_2
        LogProbValue f;
        LogProbValue f_1;

        if (back(coeffs).zero)
        {
            f.zero = true;
            f_1.zero = true;
        }
        else
        {
            f.zero = false;
            f.c = back(coeffs).c*(length(coeffs)-1);
            f.sign = back(coeffs).sign;
            f_1.zero = false;
            f_1.c = back(coeffs).c*(length(coeffs)-1)*(length(coeffs)-2);
            f_1.sign = back(coeffs).sign;
        }
        for (int i = length(coeffs)-2; i >= 2; --i) // Ignore 1st coeff -> derivative!
        {
            // Compute f_1 (idependently of f
            if (coeffs[i].zero && f_1.zero)
                f_1.zero = true;
            else if (coeffs[i].zero && !f_1.zero)
            {
                f_1.zero = false;
                f_1.c = f_1.c*beta; 
                f_1.sign = f_1.sign;
            }
            else if (!coeffs[i].zero && f_1.zero)
            {
                f_1.zero = false;
                f_1.c = coeffs[i].c*i*(i-1);       
                f_1.sign = f.sign;        // nothing to do, stays the same
            }
            else if (coeffs[i].sign && f_1.sign)             
            {
                f_1.zero = false;
                f_1.c = coeffs[i].c*i*(i-1) + f_1.c*beta;
                f_1.sign = true;        // nothing to do, stays the same
            }
            else if (!coeffs[i].sign && !f_1.sign)
            {
                f_1.zero = false;
                f_1.c = coeffs[i].c*i*(i-1) + f_1.c*beta;
                f_1.sign = false; 
            }
            else if (coeffs[i].sign && !f_1.sign)
            {
                f_1.zero = false;
                if (coeffs[i].c*i*(i-1) < f_1.c*beta) // TODO ab hier weiter
                {
                    f_1.c = f_1.c*beta - coeffs[i].c*i*(i-1);
                    f_1.sign = false;
                }
                else
                {
                    f_1.c = coeffs[i].c*i*(i-1) - f_1.c*beta;
                    f_1.sign = true;
                }
            }
            else if (!coeffs[i].sign && f_1.sign)
            {
                f_1.zero = false;
                if (coeffs[i].c*i*(i-1) < f_1.c*beta)
                {
                    f_1.c = f_1.c*beta - coeffs[i].c*i*(i-1);
                    f_1.sign = true;
                }
                else
                {
                    f_1.c = coeffs[i].c*i*(i-1) - f_1.c*beta;
                    f_1.sign = false;
                }
            }
            // Compute f
            if (coeffs[i].zero && f.zero)
                f.zero = true;
            else if (!coeffs[i].zero && f.zero)
            {
                f.zero = false;
                f.c = coeffs[i].c*i;        // derivative -> *i
                f.sign = coeffs[i].sign;
            }
            else if (coeffs[i].zero && !f.zero)
            {
                f.zero = false;
                f.c = f.c * beta;       // non-LogProb value must be on right hand side!
                f.sign = f.sign;        // nothing to do, stays the same
            }
            else if (coeffs[i].sign && f.sign)
            {
                f.zero = false;
                f.c = coeffs[i].c*i + f.c*beta;
                f.sign = true;        // nothing to do, stays the same
            }
            else if (!coeffs[i].sign && !f.sign)
            {
                f.zero = false;
                f.c = coeffs[i].c*i + f.c*beta;
                f.sign = false; 
            }
            else if (coeffs[i].sign && !f.sign)
            {
                f.zero = false;
                if (coeffs[i].c*i < f.c*beta)
                {
                    f.c = f.c*beta - coeffs[i].c*i;
                    f.sign = false;
                }
                else
                {
                    f.c = coeffs[i].c*i - f.c*beta;
                    f.sign = true;
                }
            }
            else if (!coeffs[i].sign && f.sign)
            {
                f.zero = false;
                if (coeffs[i].c*i < f.c*beta)
                {
                    f.c = f.c*beta - coeffs[i].c*i;
                    f.sign = true;
                }
                else
                {
                    f.c = coeffs[i].c*i - f.c*beta;
                    f.sign = false;
                }
           }
        } 
        if (coeffs[1].zero && f.zero)
            f.zero = true;
        else if (!coeffs[1].zero && f.zero)
        {
            f.zero = false;
            f.c = coeffs[1].c;        // derivative -> *i
            f.sign = coeffs[1].sign;
        }
        else if (coeffs[1].zero && !f.zero)
        {
            f.zero = false;
            f.c = f.c * beta;       // non-LogProb value must be on right hand side!
            f.sign = f.sign;        // nothing to do, stays the same
        }
        else if (coeffs[1].sign && f.sign)
        {
            f.zero = false;
            f.c = coeffs[1].c + f.c*beta;
            f.sign = true;        // nothing to do, stays the same
        }
        else if (!coeffs[1].sign && !f.sign)
        {
            f.zero = false;
            f.c = coeffs[1].c + f.c*beta;
            f.sign = false; 
        }
        else if (coeffs[1].sign && !f.sign)
        {
            f.zero = false;
            if (coeffs[1].c < f.c*beta)
            {
                f.c = f.c*beta - coeffs[1].c;
                f.sign = false;
            }
            else
            {
                f.c = coeffs[1].c - f.c*beta;
                f.sign = true;
            }
        }
        else if (!coeffs[1].sign && f.sign)
        {
            f.zero = false;
            if (coeffs[1].c < f.c*beta)
            {
                f.c = f.c*beta - coeffs[1].c;
                f.sign = true;
            }
            else
            {
                f.c = coeffs[1].c - f.c*beta;
                f.sign = false;
            }
       }

        long double r1;
        long double r2;
        if (f.zero) r1 = 0.0;
        else if (!f.sign) r1 = - static_cast<long double>(f.c);
        else r1 = f.c;
        if (f_1.zero) r2 = 0.0;
        else if (!f_1.sign) r2 = - static_cast<long double>(f_1.c);
        else r2 = f_1.c;

        //std::cout << "Beta logspace: " << beta << "  f': " << r1 << "  f'': " << r2 << std::endl;
        return boost::math::make_tuple(r1, r2);
    }
private:
    String<LogProbValue> coeffs;
};

// LogSpace
// Functor for Horner function with error bound
// Returns values for f'(x) and f''(x) and eB
template <typename TValue>
struct Fct_12eBL
{
    Fct_12eBL(String<LogProbValue> const& coeffs) : coeffs(coeffs)     // TODO: calculate derivatie coeffs of polynome in advance, only once 
    { // Constructor 
    }
    boost::math::tuple<TValue, TValue> operator()(bool &eBViolated, TValue const & beta)
    { // beta is estimate so far.
        LogProbValue f;
        LogProbValue f_1;
        LogProb<TValue> eB;             // constuctor with std::log(0.0);
        LogProb<TValue> eB_1;           // class LogProb should be able to handle zeros in calculations automatically... ??

        if (back(coeffs).zero)
        {
            f.zero = true;
            f_1.zero = true;
        }
        else
        {
            f.zero = false;
            f.c = back(coeffs).c*(length(coeffs)-1);
            f.sign = back(coeffs).sign;
            f_1.zero = false;
            f_1.c = back(coeffs).c*(length(coeffs)-1)*(length(coeffs)-2);
            f_1.sign = back(coeffs).sign;
        }
        for (int i = length(coeffs)-2; i >= 2; --i) // Ignore 1st coeff -> derivative!
        {
            // Compute f_1
            if (coeffs[i].zero && f_1.zero)
            {
                f_1.zero = true;
            }
            else if (coeffs[i].zero && !f_1.zero)
           {
                eB_1 = f_1.c*beta + eB_1*beta;
                f_1.zero = false;
                f_1.c = f_1.c*beta;
                f_1.sign = f_1.sign;
                eB_1 += f_1.c;
            }
            else if (!coeffs[i].zero && f_1.zero)
           {
                f_1.zero = false;
                f_1.c = coeffs[i].c*i*(i-1);      
                f_1.sign = f.sign;        // nothing to do, stays the same
                eB_1 = f_1.c;
            }
            else 
            {
                eB_1 = coeffs[i].c*i*(i-1) + f_1.c*beta + eB_1*beta;
                if (coeffs[i].sign && f_1.sign)             
                {
                    f_1.zero = false;
                    f_1.c = coeffs[i].c*i*(i-1) + f_1.c*beta;
                    f_1.sign = true;        // nothing to do, stays the same
                }
                else if (!coeffs[i].sign && !f_1.sign)
                {
                    f_1.zero = false;
                    f_1.c = coeffs[i].c*i*(i-1) + f_1.c*beta;
                    f_1.sign = false; 
                }
                else if (coeffs[i].sign && !f_1.sign)
                {
                    f_1.zero = false;
                    if (coeffs[i].c*i*(i-1) < f_1.c*beta) // TODO ab hier weiter
                    {
                        f_1.c = f_1.c*beta - coeffs[i].c*i*(i-1);
                        f_1.sign = false;
                    }
                    else
                    {
                        f_1.c = coeffs[i].c*i*(i-1) - f_1.c*beta;
                        f_1.sign = true;
                    }
                }
                else if (!coeffs[i].sign && f_1.sign)
                {
                    f_1.zero = false;
                    if (coeffs[i].c*i*(i-1) < f_1.c*beta)
                    {
                        f_1.c = f_1.c*beta - coeffs[i].c*i*(i-1);
                        f_1.sign = true;
                    }
                    else
                    {
                        f_1.c = coeffs[i].c*i*(i-1) - f_1.c*beta;
                        f_1.sign = false;
                    }
                }
               eB_1 += f_1.c;
            }
            
            // Compute f
            if (coeffs[i].zero && f.zero)
            {
                f.zero = true;
            }
            else if (!coeffs[i].zero && f.zero)
            {
                f.zero = false;
                f.c = coeffs[i].c*i;        // derivative -> *i
                f.sign = coeffs[i].sign;
                eB = f.c; 
            }
            else if (coeffs[i].zero && !f.zero)
            {
                f.zero = false;
                f.c = f.c * beta;       // non-LogProb value must be on right hand side!
                f.sign = f.sign;        // nothing to do, stays the same
                eB = f.c + f.c*beta + eB*beta; 
            }
            else                            // f != zero
            {
                eB = coeffs[i].c*i + f.c*beta + eB*beta; 
                if (coeffs[i].sign && f.sign)
                {
                    f.zero = false;
                    f.c = coeffs[i].c*i + f.c*beta;
                    f.sign = true;        // nothing to do, stays the same
                }
                else if (!coeffs[i].sign && !f.sign)
                {
                    f.zero = false;
                    f.c = coeffs[i].c*i + f.c*beta;
                    f.sign = false; 
                }
                else if (coeffs[i].sign && !f.sign)
                {
                    f.zero = false;
                    if (coeffs[i].c*i < f.c*beta)
                    {
                        f.c = f.c*beta - coeffs[i].c*i;
                        f.sign = false;
                    }
                    else
                    {
                        f.c = coeffs[i].c*i - f.c*beta;
                        f.sign = true;
                    }
                }
                else if (!coeffs[i].sign && f.sign)
                {
                    f.zero = false;
                    if (coeffs[i].c*i < f.c*beta)
                    {
                        f.c = f.c*beta - coeffs[i].c*i;
                        f.sign = true;
                    }
                    else
                    {
                        f.c = coeffs[i].c*i - f.c*beta;
                        f.sign = false;
                    }
                }
                eB += f.c; 
            }
        }
        //
        if (coeffs[1].zero && f.zero)
        {
            f.zero = true;
        }
        else if (!coeffs[1].zero && f.zero)
        {
            f.zero = false;
            f.c = coeffs[1].c;        // derivative -> *i
            f.sign = coeffs[1].sign;
        }
        else if (coeffs[1].zero && !f.zero)
        {
            f.zero = false;
            f.c = f.c * beta;       // non-LogProb value must be on right hand side!
            f.sign = f.sign;        // nothing to do, stays the same
            eB = f.c + f.c*beta + eB*beta; 
        }
        else                            // f != zero
        {
            eB = f.c*beta + eB*beta; 
            if (coeffs[1].sign && f.sign)
            {
                f.zero = false;
                f.c = coeffs[1].c + f.c*beta;
                f.sign = true;        // nothing to do, stays the same
            }
            else if (!coeffs[1].sign && !f.sign)
            {
                f.zero = false;
                f.c = coeffs[1].c + f.c*beta;
                f.sign = false; 
            }
            else if (coeffs[1].sign && !f.sign)
            {
                f.zero = false;
                if (coeffs[1].c < f.c*beta)
                {
                    f.c = f.c*beta - coeffs[1].c;
                    f.sign = false;
                }
                else
                {
                    f.c = coeffs[1].c - f.c*beta;
                    f.sign = true;
                }
            }
            else if (!coeffs[1].sign && f.sign)
            {
                f.zero = false;
                if (coeffs[1].c < f.c*beta)
                {
                    f.c = f.c*beta - coeffs[1].c;
                    f.sign = true;
                }
                else
                {
                    f.c = coeffs[1].c - f.c*beta;
                    f.sign = false;
                }
            }
            eB += f.c; 
        }

        eB = eB*(6*pow(10, -18));    // TODO u? 
        eB_1 = eB_1*(6*pow(10, -18)); 

        long double r1;
        long double r2;
        if (f.zero) r1 = 0.0;
        else if (!f.sign) r1 = - static_cast<long double>(f.c);
        else r1 = f.c;
        if (f_1.zero) r2 = 0.0;
        else if (!f_1.sign) r2 = - static_cast<long double>(f_1.c);
        else r2 = f_1.c;

        //std::cout << "Beta logspace eB: " << beta << "  f': " << r1 << "  f'': " << r2 << "  eB: " << eB << " eB_1: " << eB_1 << std::endl;
        if (f.c <= eB || f_1.c <= eB_1) eBViolated = true;
        else eBViolated = false;
        return boost::math::make_tuple(r1, r2);
    }
private:
    String<LogProbValue> coeffs;
};


// Classical polynomial evaluation
template<typename TLHood, typename TBeta, typename TCoeffs>
inline TLHood 
getProbForBeta(TBeta &beta, TCoeffs &coeffs, Classic const&)
{
    TLHood f = 0.0;
    TBeta b = 1.0;
    for (int i = 0; i < length(coeffs); ++i)
    {
        f += coeffs[i] * b;
        b *= beta;
    }
    return f;
}

// Horner
template<typename TValue, typename TCoeffs>
inline TValue 
getProbForBeta(TValue &beta, TCoeffs &coeffs, Horner const&)
{
    TValue f = back(coeffs);
    TValue eB = abs(f)/2;
    for (int i = length(coeffs)-2; i >= 0; --i)
    {
        f = coeffs[i] + beta*f;
        eB = beta*eB + abs(f);
    }
    eB = (6*pow(10, -18))*(2*eB - abs(f));

    std::cout << " eB " << eB << std::endl;
    if (f <= eB) return 0.0;
    return f;
}


// Horner in LogSpace
template<typename TValue>
inline LogProbValue 
getProbForBeta(TValue &beta, String<LogProbValue> &coeffs, Horner const&)
{
    LogProbValue f;
    LogProb<TValue> eB;
    if (back(coeffs).zero)
    {
        f.zero = true;
    }
    else
    {
        f.zero = false;
        f.c = back(coeffs).c;
        f.sign = back(coeffs).sign;
        eB = f.c/2.0;
    }
    for (int i = length(coeffs)-2; i >= 0; --i)
    {
        // Check all different cases !!!!
        if (coeffs[i].zero && f.zero)
            f.zero = true;
        else if (!coeffs[i].zero && f.zero)
        {
            f.zero = false;
            f.c = coeffs[i].c;
            f.sign = coeffs[i].sign;
        }
        else if (coeffs[i].zero && !f.zero)
        {
            f.zero = false;
            f.c = f.c * beta;       // non-LogProb value must be on right hand side!
            f.sign = f.sign;        // nothing to do, stays the same
        }
        else if (coeffs[i].sign && f.sign)
        {
            f.zero = false;
            f.c = coeffs[i].c + f.c*beta;
            f.sign = true;        // nothing to do, stays the same
        }
        else if (!coeffs[i].sign && !f.sign)
        {
            f.zero = false;
            f.c = coeffs[i].c + f.c*beta;
            f.sign = false; 
        }
        else if (coeffs[i].sign && !f.sign)
        {
            f.zero = false;
            if (coeffs[i].c < f.c*beta)
            {
                f.c = f.c*beta - coeffs[i].c;
                f.sign = false;
            }
            else
            {
                f.c = coeffs[i].c - f.c*beta;
                f.sign = true;
            }
        }
        else if (!coeffs[i].sign && f.sign)
        {
            f.zero = false;
            if (coeffs[i].c < f.c*beta)
            {
                f.c = f.c*beta - coeffs[i].c;
                f.sign = true;
            }
            else
            {
                f.c = coeffs[i].c - f.c*beta;
                f.sign = false;
            }
       }
        eB = eB*beta + f.c;
    }
    eB = (eB*2.0 - f.c)*(6*pow(10, -18));    // TODO u ?
    //std::cout << "eB: " << eB << std::endl;
    if (!f.sign || f.c <= eB) f.zero = true;
    else SEQAN_ASSERT_EQ(f.sign, true);
    return f;
}


// Horner with scaling
template<typename TValue, typename TCoeffs>
inline TValue 
getProbForBeta(TValue &beta, TCoeffs &coeffs, HornerScaling const&)
{
    int a = -10;
    int b = 10;
    int g = min(-floor((a+1)/2), floor(b/2));

    TValue f = back(coeffs);
    long double s;
    int r = 0;
    for (int i = length(coeffs)-2; i >= 0; --i)
    {
        int exp1 = floor(log10(coeffs[i])); // TODO what if = 0 check!
        int exp2 = floor(log10(beta*f));
        if (abs(exp1) > g && abs(exp2) > g)   // Both exponents are out of range -> scale 
        {
            int r_i;
            if (abs(coeffs[i]) >= abs(beta*f))    // Larger one will have exponent g
                r_i = g - exp1;
            else 
                r_i = g - exp2;
            s = pow(10, r_i);    // Scaling factor
            r += r_i;
            //std::cout << " exp1: " << exp1 << " exp2: " << exp2 << " r_i: " << r_i << " g: " << g << "  s: " << s <<  std::endl;
            //std::cout << " s*coeffs[i]: " << coeffs[i] << "  " <<  s*coeffs[i] << " s*beta*f: " << beta*f << "  " <<  s*(beta*f) << "     " << g <<  std::endl;
        }
        else
            s = 1;
        f = s*(coeffs[i] + beta*f);
    }
    //std::cout << " f: " << f << " r: " << r << std::endl;
    s = pow(10, -r);
    return (f/s);
}


// Naive
template<typename TValue, typename TConstants>
inline TValue 
getProbForBeta(TValue &beta, TConstants &constants, Naive const&)
{
    // TODO LogProb
    TValue f = 1.0;
    for (unsigned i = 0; i < length(constants[0]); ++i)
    {
        f *= constants[0][i] + beta*constants[1][i];
    }
    return f;
}

// NSpace
// Functor for  naive evaliuation
// Returns values for f(x) 
template <typename TValue>
struct FctNaive_0N
{    
    FctNaive_0N(String<String<TValue> > const& constants) : constants(constants)
    { // Constructor 
    }
    TValue operator()(TValue const&beta)
    { // beta is estimate so far.
        TValue f = 1.0;
        for (unsigned i = 0; i < length(constants[0]); ++i)
        {
            f *= (constants[0][i] + beta*constants[1][i]);
        }

        return f;
    }
private:
    String<String<TValue> > constants;
};



// Naive LogSpace
template<typename TValue>
inline LogProbValue 
getProbForBeta(TValue &beta, String<String<LogProbValue> > &constants, Naive const&)
{
    LogProbValue f;
    f.c = 1.0;
    for (unsigned i = 0; i < length(constants[0]); ++i)
    {
        if (constants[1][i].zero)       // b = 0
            f.c *= constants[0][i].c;
        else if (constants[1][i].sign) 
            f.c *= constants[0][i].c + constants[1][i].c*beta;
        else
        {
            // TODO check:  constants[0][i].c == constants[1][i]*beta
            if (constants[0][i].c < constants[1][i].c*beta)
            {
                f.c *= constants[1][i].c*beta - constants[0][i].c;
                f.sign = false;
            }
            else
                f.c *= constants[0][i].c - constants[1][i].c*beta;
        }
    }
    return f;
}

// LogSpace
// Functor for  naive evaluation
// Returns values for f(x) 
template <typename TValue>
struct FctNaive_0L
{    
    FctNaive_0L(String<String<LogProbValue> > const& constants) : constants(constants)
    { // Constructor 
    }
    LogProbValue operator()(TValue const&beta)
    { // beta is estimate so far.
        LogProbValue f;
        f.c = 1.0;
        for (int i = 0; i < length(constants[0]); ++i)
        {
            if (constants[1][i].zero)       // b = 0
                f.c *= constants[0][i].c;
            else if (constants[1][i].sign) 
                f.c *= constants[0][i].c + constants[1][i].c*beta;
            else
            {
                // TODO check:  constants[0][i].c == constants[1][i]*beta
                if (constants[0][i].c < constants[1][i].c*beta)
                {
                    f.c *= constants[1][i].c*beta - constants[0][i].c;
                    f.sign = false;
                }
                else
                    f.c *= constants[0][i].c - constants[1][i].c*beta;
            }
        }

        return f;
    }
private:
    String<String<LogProbValue> > constants;
};

// Functor for log function
// Returns values for f(x)
template <typename TValue>
struct FctLog_02N
{   
    FctLog_02N(String<String<TValue> > const& constants) : constants(constants)
    { // Constructor 
    }
    ::boost::math::tuple<TValue, TValue> operator()(TValue const&beta)
    { 
        TValue f_0 = 0;
        TValue f_2 = 0;
        for (unsigned i = 0; i < length(constants[0]); ++i)
        {
            f_0 += std::log10(constants[0][i] + constants[1][i]*beta);  
            f_2 += -pow(constants[1][i],2)/( pow(constants[0][i] , 2) + 2*constants[0][i]*constants[1][i]*beta + pow(constants[1][i], 2)*pow(beta, 2) );
        }
        f_2 *= 1.0/(std::log(10));  

        return  boost::math::make_tuple(f_0, f_2);
    }
private:
    String<String<TValue> > constants;
};

// Functor for log function
// Returns values for f'(x) and f''(x)
template <typename TValue>
struct FctLog_12N
{   
    FctLog_12N(String<String<TValue> > const& constants) : constants(constants)
    { // Constructor 
    }
    ::boost::math::tuple<TValue, TValue> operator()(TValue const&beta)
    { 
        TValue f_1 = 0;
        TValue f_2 = 0;
        for (unsigned i = 0; i < length(constants[0]); ++i)
        {
            f_1 += constants[1][i]/(constants[0][i] + constants[1][i]*beta);
            f_2 += -pow(constants[1][i],2)/( pow(constants[0][i], 2) + 2*constants[0][i]*constants[1][i]*beta + pow(constants[1][i], 2)*pow(beta, 2) );
        }

        f_1 *= 1.0/(std::log(10)) ;  // 1/(ln*b)
        f_2 *= 1.0/(std::log(10)) ;  
        
        return boost::math::make_tuple(f_1, f_2);
    }
private:
    String<String<TValue> > constants;
};


/////////////////////////////////////////////////////////////////
// Bisektions
/////////////////////////////////////////////////////////////////



// x1 valid, x2 unvalid 
// Looking for new x2 between x1 and x2_old, where error bound is not violated and sign difference holds
// Bisektion until no error bound is violated anymore
template<typename TValue, typename TFunctor, typename TFunctor2>
inline TValue
bisektionR(TValue x1, TValue x2, TValue f_1x1, TValue f_2x1, TFunctor &functor, TFunctor2 &functorN)
{
    boost::uintmax_t maxIter = 10;
    //TValue beta;
    TValue x = x1 + (x2 - x1)/2.0;
    
    bool eBViolated = false;
    boost::math::tuple<TValue, TValue> tuple = functor(eBViolated, x);

    if (eBViolated)     // Find new x
    {
        if (x - x1 < 0.01)
            return x1;
        else
            return bisektionR(x1, x, f_1x1, f_2x1, functor, functorN);
    }
    else                // Valid
    {
        long double f_1x = boost::math::get<0>(tuple);
        long double f_2x = boost::math::get<1>(tuple);

        // if sign diff: Root between x1 and x
        if ((f_1x1 < 0 && f_1x > 0) || (f_1x1 > 0 && f_1x < 0))
        {
             return boost::math::tools::newton_raphson_iterate(functorN, x, x1, x, std::numeric_limits<long double>::digits, maxIter);   // TODO: start with x? -> if no good value, at least automatically bisektion in newton is started directly
        }
        else   // else: root between x and x2
        {
            if (x - x1 < 0.01)
                return x1;
            else
                return bisektionR(x, x2, f_1x, f_2x, functor, functorN);
        }
    }
}
// x1 unvalid, x2 valid 
template<typename TValue, typename TFunctor, typename TFunctor2>
inline TValue
bisektionL(TValue x1, TValue x2, TValue f_1x2, TValue f_2x2, TFunctor &functor, TFunctor2 &functorN)
{
    boost::uintmax_t maxIter = 10;
    TValue x = x1 + (x2 - x1)/2.0;
    
    bool eBViolated = false;
    boost::math::tuple<TValue, TValue> tuple = functor(eBViolated, x);

    if (eBViolated)     // Find new x
    {
        if (x2 - x < 0.01)      // if no big difference anymore and still error bound violated, stop!
            return x2;
        else
            return bisektionL(x, x2, f_1x2, f_2x2, functor, functorN);
    }
    else                // Valid
    {
        long double f_1x = boost::math::get<0>(tuple);
        long double f_2x = boost::math::get<1>(tuple);
        // if sign diff: Root between x and x2
        if ( ((f_1x2 < 0 && f_1x > 0) || (f_1x2 > 0 && f_1x < 0)) )
        {
             return boost::math::tools::newton_raphson_iterate(functorN, x, x, x2, std::numeric_limits<long double>::digits, maxIter);   // TODO: start with x? -> if no good value, at least automatically bisektion in newton is started directly
        }
        else // else: root between x1 and x
        {
            if (x2 - x < 0.01)      // if no big difference anymore and still no sign diff, stop!
                return x2;
            else
                return bisektionL(x1, x, f_1x, f_2x, functor, functorN);
        }
    }
    return 0;
}

#endif
