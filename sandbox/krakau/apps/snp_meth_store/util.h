/*==========================================================================

Functors and functions for likelihood function evaluation

==========================================================================*/

#ifndef SANDBOX_KRAKAU_APPS_SNP_METH_STORE_UTIL_H_
#define SANDBOX_KRAKAU_APPS_SNP_METH_STORE_UTIL_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

using namespace std;
using namespace seqan;

// Prob. function
struct Naive_;
typedef Tag<Naive_> Naive;

struct LogFunction_;
typedef Tag<LogFunction_> LogFunction;


// Optimize
struct Sampling_;
typedef Tag<Sampling_> Sampling;

struct Newton_;
typedef Tag<Newton_> Newton;


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


#endif
