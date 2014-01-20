/*==========================================================================

Methylation Calling

==========================================================================*/

#ifndef __SANDBOX_KRAKAU_APPS_SNP_METH_STORE_METHS_H__
#define __SANDBOX_KRAKAU_APPS_SNP_METH_STORE_METHS_H__

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

#include "bs_alphabets.h"
#include "meths.h"

using namespace std;
using namespace seqan;


template<typename TDna>
bool
isTransition(TDna base1, TDna base2)  // why not error with & ?
{
    if (base1 == 'A' && base2 == 'G') return true;
    else if (base1 == 'G' && base2 == 'A') return true;
    else if (base1 == 'C' && base2 == 'T') return true;
    else if (base1 == 'T' && base2 == 'C') return true;
    else return false;
}

template<typename TMethOptions, typename TOptions>
void
computeGenotypePriors(TMethOptions &methOptions, TOptions &options)
{
    resize(methOptions.genPriors, 4*4*4, -1.0, Exact());
    long double pHetSnp = 0.001;
    long double pHomoSnp = 0.0005;
    long double r_j;
    long double r_k;

    if(options._debugLevel > 1)
    {
        std::cout << "Genotype prior probabilities: "<< std::endl;
        std::cout << "Ref." << '\t' << "Allel1" << '\t' << "Allele2" << std::endl;
    }

    if (!methOptions.uniformGenPriors)
    {
        for (unsigned i = 0; i < 4; ++i)
        {
            for (unsigned j = 0; j < 4; ++j)
            {
                for (unsigned k = 0; k < 4; ++k)
                {
                    r_j = (isTransition((Dna)i, (Dna)j)? (double)(4.0/6.0):(double)(1.0/6.0));
                    r_k = (isTransition((Dna)i, (Dna)k)? (double)(4.0/6.0):(double)(1.0/6.0));
                    
                    if (j == i && k == i)
                        methOptions.genPriors[i<<4| j<<2| k] = (1.0-pHetSnp -pHomoSnp);
                    else if (j == i)
                        methOptions.genPriors[i<<4| j<<2| k] = pHetSnp*r_k;
                    else if (k == i)
                        methOptions.genPriors[i<<4| j<<2| k] = pHetSnp*r_j;
                    else if (j == k)
                        methOptions.genPriors[i<<4| j<<2| k] = pHomoSnp*r_j;    // + pError ?
                    else
                        methOptions.genPriors[i<<4| j<<2| k] = pHetSnp * pHomoSnp;  // pHetSnp*(4/6)*0.0005*(4/6) + pRefError*pHetSnp*(4/6)

                    if(options._debugLevel > 1)
                    {              
                        std::cout << (Dna)i << '\t' << (Dna)j << '\t' << (Dna)k << "\t\t" << (i<<4| j<<2| k) << '\t' << std::setprecision (25) << methOptions.genPriors[i<<4| j<<2| k] << std::endl;
                    }
                }
            }
        }
    }
    else
    {
        for (unsigned i = 0; i < 4; ++i)
        {
            for (unsigned j = 0; j < 4; ++j)
            {
                for (unsigned k = 0; k < 4; ++k)
                {
                    if (j == i && k == i)
                        methOptions.genPriors[i<<4| j<<2| k] = (1.0-pHetSnp -pHomoSnp);
                    else if (j == i)
                        methOptions.genPriors[i<<4| j<<2| k] = pHetSnp;
                    else if (k == i)
                        methOptions.genPriors[i<<4| j<<2| k] = pHetSnp;
                    else if (j == k)
                        methOptions.genPriors[i<<4| j<<2| k] = pHomoSnp;    
                    else
                        methOptions.genPriors[i<<4| j<<2| k] = pHetSnp * pHomoSnp;  // pHetSnp*(4/6)*0.0005*(4/6) + pRefError*pHetSnp*(4/6)
               }
            }
        }
    }
}


#endif

