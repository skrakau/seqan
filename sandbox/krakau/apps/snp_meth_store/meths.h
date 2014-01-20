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

// Compute meth. state prior probabilities for all genotypes (and possible meth. states) and for all possible contexts or context combinations
template<typename TMethOptions, typename TOptions>
void
computeBsPriors(TMethOptions &methOptions, TOptions &/*options*/)
{
    // access with bsPriors[h2][h1][context], while h1 <= h2

    // Prepare table
    clear(methOptions.bsPriors);
    resize(methOptions.bsPriors, 6);
    // for each haplotype h2
    for (unsigned h2 = 0; h2 < 6; ++h2)
    {
       resize(methOptions.bsPriors[h2], h2+1);
       // for each haplotype h1 with h1 <= h2
       for (unsigned h1 = 0; h1 <= h2; ++h1)
       {
            // TODO: check if (Dna)testDnaM works properly !!!
            if ((Dna)((DnaM)h1) == 'C' && (Dna)((DnaM)h2) == 'G') // for genotype CG we need to take different contexts for both strands into account
                resize(methOptions.bsPriors[h2][h1], 9); 
            else
                resize(methOptions.bsPriors[h2][h1], 3);    
       }
    }

    // 
    for (unsigned h2 = 0; h2 < 6; ++h2)
    {
        for (unsigned h1 = 0; h1 <=h2 ; ++h1)
        {
            // For homo cases 
            //long double bsPrior = 1.0;
            if ((DnaM)h1 == 'C' && (DnaM)h2 == 'C')     // CC
            {
                methOptions.bsPriors[h2][h1][0] = methOptions.homoNonMethPriorCG; 
                methOptions.bsPriors[h2][h1][1] = methOptions.homoNonMethPriorCHG; 
                methOptions.bsPriors[h2][h1][2] = methOptions.homoNonMethPriorCHH; 
            }
            else if ((DnaM)h1 == 'D' && (DnaM)h2 == 'D')     // CmCm
            {
                methOptions.bsPriors[h2][h1][0] = methOptions.homoMethPriorCG; 
                methOptions.bsPriors[h2][h1][1] = methOptions.homoMethPriorCHG; 
                methOptions.bsPriors[h2][h1][2] = methOptions.homoMethPriorCHH; 
            }
            else if ((DnaM)h1 == 'C' && (DnaM)h2 == 'D')     // CCm
            {
                methOptions.bsPriors[h2][h1][0] = methOptions.heteroMethPriorCG; 
                methOptions.bsPriors[h2][h1][1] = methOptions.heteroMethPriorCHG; 
                methOptions.bsPriors[h2][h1][2] = methOptions.heteroMethPriorCHH; 
            }
            if ((DnaM)h1 == 'G' && (DnaM)h2 == 'G')     // GG
            {
                methOptions.bsPriors[h2][h1][0] = methOptions.homoNonMethPriorCG; 
                methOptions.bsPriors[h2][h1][1] = methOptions.homoNonMethPriorCHG; 
                methOptions.bsPriors[h2][h1][2] = methOptions.homoNonMethPriorCHH; 
            }
            else if ((DnaM)h1 == 'H' && (DnaM)h2 == 'H')     // GmGm
            {
                methOptions.bsPriors[h2][h1][0] = methOptions.homoMethPriorCG; 
                methOptions.bsPriors[h2][h1][0] = methOptions.homoMethPriorCHG; 
                methOptions.bsPriors[h2][h1][0] = methOptions.homoMethPriorCHH; 
            }
            else if ((DnaM)h1 == 'G' && (DnaM)h2 == 'H')     // GGm
            {
                methOptions.bsPriors[h2][h1][0] = methOptions.heteroMethPriorCG; 
                methOptions.bsPriors[h2][h1][0] = methOptions.heteroMethPriorCHG; 
                methOptions.bsPriors[h2][h1][0] = methOptions.heteroMethPriorCHH; 
            }
            else 
            {                               
                // Snp: mixed cases 
                // (Maybe later: take into account, if reAllele is C or not)
                if ((DnaM)h1 == 'C' && (DnaM)h2 == 'G')         // CG
                {
                    methOptions.bsPriors[h2][h1][0] = (1-methOptions.cxRefCMethPriorCG)*(1-methOptions.cxRefCMethPriorCG);          // 00 forward context CG, reverse CG
                    methOptions.bsPriors[h2][h1][1] = (1-methOptions.cxRefCMethPriorCG)*(1-methOptions.cxRefCMethPriorCHG);         // 01
                    methOptions.bsPriors[h2][h1][2] = (1-methOptions.cxRefCMethPriorCG)*(1-methOptions.cxRefCMethPriorCHH);         // 02
                    methOptions.bsPriors[h2][h1][3] = (1-methOptions.cxRefCMethPriorCHG)*(1-methOptions.cxRefCMethPriorCG);          // 10 double, not neccessary !?
                    methOptions.bsPriors[h2][h1][4] = (1-methOptions.cxRefCMethPriorCHG)*(1-methOptions.cxRefCMethPriorCHG);         // 11
                    methOptions.bsPriors[h2][h1][5] = (1-methOptions.cxRefCMethPriorCHG)*(1-methOptions.cxRefCMethPriorCHH);         // 12
                    methOptions.bsPriors[h2][h1][6] = (1-methOptions.cxRefCMethPriorCHH)*(1-methOptions.cxRefCMethPriorCG);          // 20 ...
                    methOptions.bsPriors[h2][h1][7] = (1-methOptions.cxRefCMethPriorCHH)*(1-methOptions.cxRefCMethPriorCHG);         // 21 ...
                    methOptions.bsPriors[h2][h1][8] = (1-methOptions.cxRefCMethPriorCHH)*(1-methOptions.cxRefCMethPriorCHH);         // 22
                }
                else if ((DnaM)h1 == 'D' && (DnaM)h2 == 'H')    // CmGm
                {
                    methOptions.bsPriors[h2][h1][0] = methOptions.cxRefCMethPriorCG*methOptions.cxRefCMethPriorCG;          // 00 forward context CG, reverse CG
                    methOptions.bsPriors[h2][h1][1] = methOptions.cxRefCMethPriorCG*methOptions.cxRefCMethPriorCHG;         // 01
                    methOptions.bsPriors[h2][h1][2] = methOptions.cxRefCMethPriorCG*methOptions.cxRefCMethPriorCHH;         // 02
                    methOptions.bsPriors[h2][h1][3] = methOptions.cxRefCMethPriorCHG*methOptions.cxRefCMethPriorCG;          // 10 
                    methOptions.bsPriors[h2][h1][4] = methOptions.cxRefCMethPriorCHG*methOptions.cxRefCMethPriorCHG;         // 11
                    methOptions.bsPriors[h2][h1][5] = methOptions.cxRefCMethPriorCHG*methOptions.cxRefCMethPriorCHH;         // 12
                    methOptions.bsPriors[h2][h1][6] = methOptions.cxRefCMethPriorCHH*methOptions.cxRefCMethPriorCG;          // 20 
                    methOptions.bsPriors[h2][h1][7] = methOptions.cxRefCMethPriorCHH*methOptions.cxRefCMethPriorCHG;         // 21 
                    methOptions.bsPriors[h2][h1][8] = methOptions.cxRefCMethPriorCHH*methOptions.cxRefCMethPriorCHH;        // 22 
                }
                else if ((DnaM)h1 == 'D' && (DnaM)h2 == 'G')    // CmG
                {
                    methOptions.bsPriors[h2][h1][0] = (methOptions.cxRefCMethPriorCG)*(1-methOptions.cxRefCMethPriorCG);          // 00 forward context CG, reverse CG
                    methOptions.bsPriors[h2][h1][1] = (methOptions.cxRefCMethPriorCG)*(1-methOptions.cxRefCMethPriorCHG);         // 01
                    methOptions.bsPriors[h2][h1][2] = (methOptions.cxRefCMethPriorCG)*(1-methOptions.cxRefCMethPriorCHH);         // 02
                    methOptions.bsPriors[h2][h1][3] = (methOptions.cxRefCMethPriorCHG)*(1-methOptions.cxRefCMethPriorCG);          // 10 
                    methOptions.bsPriors[h2][h1][4] = (methOptions.cxRefCMethPriorCHG)*(1-methOptions.cxRefCMethPriorCHG);         // 11
                    methOptions.bsPriors[h2][h1][5] = (methOptions.cxRefCMethPriorCHG)*(1-methOptions.cxRefCMethPriorCHH);         // 12
                    methOptions.bsPriors[h2][h1][6] = (methOptions.cxRefCMethPriorCHH)*(1-methOptions.cxRefCMethPriorCG);          // 20 
                    methOptions.bsPriors[h2][h1][7] = (methOptions.cxRefCMethPriorCHH)*(1-methOptions.cxRefCMethPriorCHG);         // 21 
                    methOptions.bsPriors[h2][h1][8] = (methOptions.cxRefCMethPriorCHH)*(1-methOptions.cxRefCMethPriorCHH);        // 22 
                }
                else if ((DnaM)h1 == 'C' && (DnaM)h2 == 'H')    // CGm
                {
                    methOptions.bsPriors[h2][h1][0] = (1-methOptions.cxRefCMethPriorCG)*methOptions.cxRefCMethPriorCG;          // 00 forward context CG, reverse CG
                    methOptions.bsPriors[h2][h1][1] = (1-methOptions.cxRefCMethPriorCG)*methOptions.cxRefCMethPriorCHG;         // 01
                    methOptions.bsPriors[h2][h1][2] = (1-methOptions.cxRefCMethPriorCG)*methOptions.cxRefCMethPriorCHH;         // 02
                    methOptions.bsPriors[h2][h1][3] = (1-methOptions.cxRefCMethPriorCHG)*methOptions.cxRefCMethPriorCG;          // 10 double !?
                    methOptions.bsPriors[h2][h1][4] = (1-methOptions.cxRefCMethPriorCHG)*methOptions.cxRefCMethPriorCHG;         // 11
                    methOptions.bsPriors[h2][h1][5] = (1-methOptions.cxRefCMethPriorCHG)*methOptions.cxRefCMethPriorCHH;         // 12
                    methOptions.bsPriors[h2][h1][6] = (1-methOptions.cxRefCMethPriorCHH)*methOptions.cxRefCMethPriorCG;          // 20 ...
                    methOptions.bsPriors[h2][h1][7] = (1-methOptions.cxRefCMethPriorCHH)*methOptions.cxRefCMethPriorCHG;         // 21 ...
                    methOptions.bsPriors[h2][h1][8] = (1-methOptions.cxRefCMethPriorCHH)*methOptions.cxRefCMethPriorCHH;        // 22 
                }
                else if ((DnaM)h1 == 'C' || (DnaM)h2 == 'C')             // CX
                {
                    methOptions.bsPriors[h2][h1][0] = 1-methOptions.cxRefCMethPriorCG;
                    methOptions.bsPriors[h2][h1][1] = 1-methOptions.cxRefCMethPriorCHG;
                    methOptions.bsPriors[h2][h1][2] = 1-methOptions.cxRefCMethPriorCHH;
                }
                else if ((DnaM)h1 == 'D' || (DnaM)h2 == 'D')        // CmX
                {
                    methOptions.bsPriors[h2][h1][0] = methOptions.cxRefCMethPriorCG;
                    methOptions.bsPriors[h2][h1][1] = methOptions.cxRefCMethPriorCHG;
                    methOptions.bsPriors[h2][h1][2] = methOptions.cxRefCMethPriorCHH;
                }
                else if ((DnaM)h1 == 'G' || (DnaM)h2 == 'G')                 // GX
                {
                    methOptions.bsPriors[h2][h1][0] = 1-methOptions.cxRefCMethPriorCG;
                    methOptions.bsPriors[h2][h1][1] = 1-methOptions.cxRefCMethPriorCHG;
                    methOptions.bsPriors[h2][h1][2] = 1-methOptions.cxRefCMethPriorCHH;
                }
                else if ((DnaM)h1 == 'H' || (DnaM)h2 == 'H')             // GmX
                {
                    methOptions.bsPriors[h2][h1][0] = methOptions.cxRefCMethPriorCG;
                    methOptions.bsPriors[h2][h1][1] = methOptions.cxRefCMethPriorCHG;
                    methOptions.bsPriors[h2][h1][2] = methOptions.cxRefCMethPriorCHH;
                }
                else
                {
                    methOptions.bsPriors[h2][h1][0] = 1.0;              // For non C/G genotypes
                    methOptions.bsPriors[h2][h1][1] = 1.0; 
                    methOptions.bsPriors[h2][h1][2] = 1.0; 
                }
            }
        }
    }
}









#endif

