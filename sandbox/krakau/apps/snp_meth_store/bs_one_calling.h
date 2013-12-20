/*==========================================================================

Methylation Calling: Version 2

Calculate snp and methylation probability in one step
        Pr(Dj | G = A:T)
    Pr(Dj | G = C^(1-beta)C_M^(beta):T)
    ...

Add additional cases: e.g. prob. converted and sequencing error
Add underconversion rates...

==========================================================================*/

#ifndef __SANDBOX_KRAKAU_APPS_SNP_METH_STORE_BS_ONE_CALLING_H__
#define __SANDBOX_KRAKAU_APPS_SNP_METH_STORE_BS_ONE_CALLING_H__

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

//#include <boost/math/distributions.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/roots.hpp>


#include "bs_alphabets.h"
#include "util.h"
//#include "snp_meth_store.h"
#include "meths.h"


using namespace std;
using namespace seqan;



// TODO:
// precalculate probs and get rid of thousand different cases to check
// Compute prob. to observe given base i under the assumption that the underlying haplotype h
// For reads mapped to the forward strand
template<typename TProb, typename TErrorProb, typename TMethOptions>
inline void
getSingleBaseProbHaploF(TProb &singleProb, Dna i, DnaM h, TErrorProb &e, bool & origin, TMethOptions &methOptions)
{
    if (methOptions.uniformSeqErrorsCalling) 
    {
        if ( h == 'C')                              // Haplotype C
        {                     
            if (i == 'C')                                               // correct and not converted + converted and error  
                singleProb = (1.0-e)*(1-methOptions.convRate) ; //+ methOptions.convRate*(e/3.0);  
            else if (i == 'T')                                          // error + correct and converted 
                singleProb = e/3.0 + (1.0-e)*(methOptions.convRate);          
            else                                                        // error            
                singleProb = e/3.0;
        } 
        else if ( h == 'D')                         // Haplotype Cm
        {                     
            if (i == 'C')                                               // correct and not converted + converted and error    
                singleProb = (1.0-e)*(1.0-methOptions.methConvRate); // + methOptions.methConvRate*(e/3.0);
            else if (i == 'T')                                          // error + correct and converted
                singleProb = e/3.0 + (1.0-e)*methOptions.methConvRate;
            else                                                        // error 
                singleProb = e/3.0;
        } 
        else if ( h == 'G')                         // Haplotype G
        {
            if (i == 'G')          
                singleProb = 1.0 - e;
            else                                  
                singleProb = e/3.0;
        } 
        else if ( h == 'H')                         // Haplotype Gm
        {
            if (i == 'G')         
                singleProb = 1.0 - e;
            else                                    
                singleProb = e/3.0;
        } 
        else                                        // Haplotype T, A: no bs 
        {
            if (i == h)
                singleProb = 1.0 - e; 
            else
                singleProb = e/3.0; 
        }
    }
    else
    {
        // temporary, just to try
        // 
        TProb const *seqErrorFreqs = SeqErrorFreqsN<TProb, BsNonSimple>::getData();  
        //TProb const *seqErrorFreqsTo = SeqErrorFreqsTo<TProb, BsNonSimple>::getData();  

        if (origin)
        {
            if ( h == 'C')                              // Haplotype C
            {                     
                if (i == 'C')                                                
                    singleProb = (1.0-e)*(1-methOptions.convRate) ; 
                else if (i == 'T')                                       
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)] + (1.0-e)*(methOptions.convRate);          
                else                                                
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)];
            } 
            else if ( h == 'D')                         // Haplotype Cm
            {                     
                if (i == 'C')                                                  
                    singleProb = (1.0-e)*(1.0-methOptions.methConvRate); 
                else if (i == 'T')                          
                    singleProb = e*seqErrorFreqs[1*5 + ordValue(i)]  + (1.0-e)*methOptions.methConvRate;
                else                                                        
                    singleProb = e*seqErrorFreqs[1*5 + ordValue(i)] ;
            } 
            else if ( h == 'G')                         // Haplotype G
            {
                if (i == 'G')          
                    singleProb = 1.0 - e;
                else                                  
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)] ;
            } 
            else if ( h == 'H')                         // Haplotype Gm
            {
                if (i == 'G')         
                    singleProb = 1.0 - e;
                else                                    
                    singleProb = e*seqErrorFreqs[2*5 + ordValue(i)] ;
            } 
            else                                        // Haplotype T, A: no bs 
            {
                if (i == h)
                    singleProb = 1.0 - e; 
                else
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)] ; 
            }
        }
        else    // use sequencing error probs from RC
        {
            FunctorDna5OrdValueComplement<int> fCompl; 
             if ( h == 'C')                              // Haplotype C
            {                     
                if (i == 'C')                                               
                    singleProb = (1.0-e)*(1-methOptions.convRate) ; 
                else if (i == 'T')                                           
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))] + (1.0-e)*(methOptions.convRate);          
                else                                                               
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))];
            } 
            else if ( h == 'D')                         // Haplotype Cm
            {                     
                if (i == 'C')                                        
                    singleProb = (1.0-e)*(1.0-methOptions.methConvRate); 
                else if (i == 'T')                                      
                    singleProb = e*seqErrorFreqs[fCompl(1)*5 + fCompl(ordValue(i))]  + (1.0-e)*methOptions.methConvRate;
                else                                            
                    singleProb = e*seqErrorFreqs[fCompl(1)*5 + fCompl(ordValue(i))] ;
            } 
            else if ( h == 'G')                         // Haplotype G
            {
                if (i == 'G')          
                    singleProb = 1.0 - e;
                else                                  
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))] ;
            } 
            else if ( h == 'H')                         // Haplotype Gm
            {
                if (i == 'G')         
                    singleProb = 1.0 - e;
                else                                    
                    singleProb = e*seqErrorFreqs[fCompl(2)*5 + fCompl(ordValue(i))] ;
            } 
            else                                        // Haplotype T, A: no bs 
            {
                if (i == h)
                    singleProb = 1.0 - e; 
                else
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))] ; 
            }       
        }
    }
}
// For reads mapped to the reverse strand
template<typename TProb, typename TErrorProb, typename TMethOptions>
inline void
getSingleBaseProbHaploR(TProb &singleProb, Dna i, DnaM h, TErrorProb &e, bool & origin, TMethOptions &methOptions)
{
    if (methOptions.uniformSeqErrorsCalling)
    {
        // i is read base corresponding to forward strand
        if ( h == 'C')                              // Haplotype C
        {                     
            if (i == 'C')           
                singleProb = 1.0-e;  
            else                                   
                singleProb = e/3.0;
        } 
        else if ( h == 'D')                         // Haplotype Cm
        {                     
            if (i == 'C')          
                singleProb = 1.0-e;
            else                                   
                singleProb = e/3.0;
        } 
        else if ( h == 'G')                         // Haplotype G
        {
            if (i == 'G')          
                singleProb = (1.0-e)*(1-methOptions.convRate); //+ methOptions.convRate*(e/3.0);    // TODO: why not included???
            else if (i == 'A')
                singleProb = e/3.0 + (1.0-e)*methOptions.convRate;
            else                                  
                singleProb = e/3.0;
        } 
        else if ( h == 'H')                         // Haplotype Gm
        {
            if (i == 'G')         
                singleProb = (1.0-e)*(1-methOptions.methConvRate); // + methOptions.methConvRate*(e/3.0);
            else if (i == 'A')
                singleProb = e/3.0 + (1.0-e)*methOptions.methConvRate;
            else                                    
                singleProb = e/3.0;
         } 
        else                                        // Haplotype T, A: no bs 
        {
            if (i == h)
                singleProb = 1.0 - e; 
            else
                singleProb = e/3.0; 
        }
    }
    else
    {
        // temporary, just to try
        TProb const *seqErrorFreqs = SeqErrorFreqsN<TProb, BsNonSimple>::getData(); 
        if (origin)
        {
            FunctorDna5OrdValueComplement<int> fCompl; 
           // i is read base corresponding to forward strand
            if ( h == 'C')                              // Haplotype C
            {                     
                if (i == 'C')           
                    singleProb = 1.0-e;  
                else                                   
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))];
            } 
            else if ( h == 'D')                         // Haplotype Cm
            {                     
                if (i == 'C')          
                    singleProb = 1.0-e;
                else                                   
                    singleProb = e*seqErrorFreqs[fCompl(1)*5 + fCompl(ordValue(i))];
            } 
            else if ( h == 'G')                         // Haplotype G
            {
                if (i == 'G')          
                    singleProb = (1.0-e)*(1-methOptions.convRate); //+ methOptions.convRate*(e/3.0);    // TODO: why not included???
                else if (i == 'A')
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))] + (1.0-e)*methOptions.convRate;
                else                                  
                    singleProb = e/3.0;
            } 
            else if ( h == 'H')                         // Haplotype Gm
            {
                if (i == 'G')         
                    singleProb = (1.0-e)*(1-methOptions.methConvRate); // + methOptions.methConvRate*(e/3.0);
                else if (i == 'A')
                    singleProb = e*seqErrorFreqs[fCompl(2)*5 + fCompl(ordValue(i))] + (1.0-e)*methOptions.methConvRate;
                else                                    
                    singleProb = e*seqErrorFreqs[fCompl(2)*5 + fCompl(ordValue(i))];
             } 
            else                                        // Haplotype T, A: no bs 
            {
                if (i == h)
                    singleProb = 1.0 - e; 
                else
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))]; 
            }
        }
        else
        {
            // This read was projected to original BS strand, thus this is the original base
            if ( h == 'C')                              // Haplotype C
            {                     
                if (i == 'C')           
                    singleProb = 1.0-e;  
                else                                   
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)];
            } 
            else if ( h == 'D')                         // Haplotype Cm
            {                     
                if (i == 'C')          
                    singleProb = 1.0-e;
                else                                   
                    singleProb = e*seqErrorFreqs[1*5 + ordValue(i)];
            } 
            else if ( h == 'G')                         // Haplotype G
            {
                if (i == 'G')          
                    singleProb = (1.0-e)*(1-methOptions.convRate); //+ methOptions.convRate*(e/3.0);    // TODO: why not included???
                else if (i == 'A')
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)] + (1.0-e)*methOptions.convRate;
                else                                  
                    singleProb = e/3.0;
            } 
            else if ( h == 'H')                         // Haplotype Gm
            {
                if (i == 'G')         
                    singleProb = (1.0-e)*(1-methOptions.methConvRate); // + methOptions.methConvRate*(e/3.0);
                else if (i == 'A')
                    singleProb = e*seqErrorFreqs[2*5 + ordValue(i)] + (1.0-e)*methOptions.methConvRate;
                else                                    
                    singleProb = e*seqErrorFreqs[2*5 + ordValue(i)];
             } 
            else                                        // Haplotype T, A: no bs 
            {
                if (i == h)
                    singleProb = 1.0 - e; 
                else
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)]; 
            }
        }
    }
}


// For Naive multiplication of lHoods
template<typename TConstants, typename TValue>
inline void
addFactorToConstants(TConstants &constants, TValue &pC, TValue &pCm, TValue &pOther, unsigned &r, Naive const &)
{
    SEQAN_ASSERT_NEQ(pC, 0.0);
    SEQAN_ASSERT_NEQ(pCm, 0.0);

    // Compute constant values:
    if (pOther > 0)
    {
        constants[0][r] = 0.5*(pC + pOther);
        constants[1][r] = 0.5*(-pC + pCm);
    }
    else
    {
        constants[0][r] = pC;
        constants[1][r] = -pC + pCm;
    }
}

// For Naive multiplication of lHoods
template< typename TValue>
inline void
addFactorToConstants(String<String<LogProbValue> > &constants, TValue &pC, TValue &pCm, TValue &pOther, unsigned &r, Naive const &)
{
    SEQAN_ASSERT_NEQ(pC, 0.0);
    SEQAN_ASSERT_NEQ(pCm, 0.0);
    // Compute constant values:
    if (pOther > 0)
    {
        constants[0][r].c = 0.5*(pC + pOther);
        if (pC == pCm)
            constants[1][r].zero = true;
        else if (pC < pCm)
            constants[1][r].c = 0.5*(pCm - pC);
        else
        {
            constants[1][r].c = 0.5*(pC - pCm);
            constants[1][r].sign = false;
        }
    }
    else
    {
        constants[0][r].c = pC;
        if (pC == pCm)
            constants[1][r].zero = true;
        else if (pC < pCm)
            constants[1][r].c = pCm - pC;
        else
        {
            constants[1][r].c = pC - pCm;
            constants[1][r].sign = false;
        }
    }
}


template<typename TCoeffs, typename TValue>
inline void
addFactorToConstants(TCoeffs &coeffs, TValue &pC, TValue &pCm, TValue &pOther, unsigned &r, Polynomial const &)
{
    SEQAN_ASSERT_NEQ(pC, 0.0);
    SEQAN_ASSERT_NEQ(pCm, 0.0);
    SEQAN_ASSERT_LEQ(r+1, length(coeffs)-1);
    // Compute constant values:
    TValue a;
    TValue b;
    if (pOther > 0)
    {
        a = 0.5*(pC + pOther);
        b = 0.5*(-pC + pCm);
    }
    else
    {
        a = pC;
        b = -pC + pCm;
    }

    SEQAN_ASSERT_LT(r, length(coeffs));
    if (r == 0) 
    {
        coeffs[0] = a;
        coeffs[1] = b;
    }
    else 
    {
        coeffs[r+1] = b*coeffs[r];
        for (unsigned ii = r; ii > 0; --ii)
        {
            coeffs[ii] = a*coeffs[ii] + b*coeffs[ii-1] ;
        }
        coeffs[0] = a*coeffs[0];
    }
}


template<typename TValue>
inline void
addFactorToConstants(String<LogProbValue> &coeffs, TValue &pC, TValue &pCm, TValue &pOther, unsigned &r, Polynomial const &)
{
    SEQAN_ASSERT_NEQ(pC, 0.0);
    SEQAN_ASSERT_NEQ(pCm, 0.0);


    // Compute constant values:
    LogProbValue a;
    LogProbValue b;
    if (pOther > 0.0)         
    {
        a.c = 0.5*(pC + pOther);        // TODO check if > 0, otherwise a.zero true!
        if (pC == pCm)
            b.zero = true;
        else if (pC < pCm) 
            b.c = 0.5*(pCm - pC);
        else
        {
            b.c = 0.5*(pC - pCm);
            b.sign = false;
        }
    }
    else
    {
        if (pC == 0)               // TODO in the moment assume single prob are coming in as non-log values
            a.zero = true;
        else
            a.c = pC;
        if (pC == pCm)
            b.zero = true;
        else if (pC < pCm)
            b.c = pCm - pC;
        else
        {
            b.c = pC - pCm;
            b.sign = false;
        }
    }

    SEQAN_ASSERT_LT(r, length(coeffs));
    // Recompute each coeff:
    if (r == 0) 
    {
        if (a.zero)
            coeffs[0].zero = true;
        else
        {
            coeffs[0].zero = false;
            coeffs[0].c = a.c;
            coeffs[0].sign = a.sign;    // nothing to do, is true
        }
        if (b.zero)
            coeffs[1].zero = true;
        else
        {
            coeffs[1].zero = false;
            coeffs[1].c = b.c;
            coeffs[1].sign = b.sign;
        }
    }
    else 
    {
        if (b.zero || coeffs[r].zero)
            coeffs[r+1].zero = true;
        else
        {
            coeffs[r+1].zero = false;
            coeffs[r+1].c = b.c*coeffs[r].c;
            if ((b.sign && !coeffs[r].sign) || (!b.sign && coeffs[r].sign)) 
               coeffs[r+1].sign = false;
            else
                coeffs[r+1].sign = true;
        }
        for (unsigned ii = r; ii > 0; --ii)
        {
            //std::cout << " ii: " << ii << "  a: " << a.sign << " " << a.c << "  coeffs[ii]: " << coeffs[ii].sign << " " << coeffs[ii].c << "  +  b: " << b.sign << " " << b.c << "  coeffs[ii-1]: "<< coeffs[ii-1].sign << " " <<  coeffs[ii-1].c << std::endl;
            if ((a.zero || coeffs[ii].zero) && (b.zero || coeffs[ii-1].zero))
                coeffs[ii].zero = true;
            else if (a.zero || coeffs[ii].zero)
            {
                coeffs[ii].c = b.c*coeffs[ii-1].c;
                if ((b.sign && !coeffs[ii-1].sign) || (!b.sign && coeffs[ii-1].sign)) 
                    coeffs[ii].sign = false;
                else
                   coeffs[ii].sign = true;

                coeffs[ii].zero = false;
            }
            else if (b.zero || coeffs[ii-1].zero)
            {
                coeffs[ii].c = a.c*coeffs[ii].c;
                if ((a.sign && !coeffs[ii].sign) || (!a.sign && coeffs[ii].sign)) 
                    coeffs[ii].sign = false;
                else
                   coeffs[ii].sign = true;

                coeffs[ii].zero = false;
            }
            else
            {
                coeffs[ii].zero = false;
                LogProb<TValue> part1 = a.c*coeffs[ii].c;
                LogProb<TValue> part2 = b.c*coeffs[ii-1].c;
                bool sign1 = true;
                bool sign2 = true;
                if (!coeffs[ii].sign)
                    sign1 = false;
                if ((!b.sign && coeffs[ii-1].sign) || (b.sign && !coeffs[ii-1].sign))
                    sign2 = false;

                if (sign1 && sign2)
                {
                    coeffs[ii].c = part1 + part2;
                    coeffs[ii].sign = true; 

                }
                else if (!sign1 && !sign2)
                {
                    coeffs[ii].c = part1 + part2;
                    coeffs[ii].sign = false;
                }
                else if (sign1 && !sign2)
                {
                    if (part1 > part2)
                    {
                        coeffs[ii].c = part1 - part2;
                        coeffs[ii].sign = true;
                    }
                    else
                    {
                        coeffs[ii].c = part2 - part1;
                        coeffs[ii].sign = false;
                    }
                }
                else if (!sign1 && sign2)
                {
                    if (part1 > part2)
                    {
                        coeffs[ii].c = part1 - part2;
                        coeffs[ii].sign = false;
                    }
                    else
                    {
                        coeffs[ii].c = part2 - part1;
                        coeffs[ii].sign = true;
                    }
                }
           }
        //std::cout << "            coeffs[ii]: " << coeffs[ii].sign << " " << coeffs[ii].c <<  std::endl;
        }
        if (a.zero || coeffs[0].zero)
            coeffs[0].zero = true;
        else
        {
            coeffs[0].zero = false;
            coeffs[0].c = a.c*coeffs[0].c;
            if ((a.sign && !coeffs[0].sign) || (!a.sign && coeffs[0].sign)) 
                coeffs[0].sign = false;
            else
                coeffs[0].sign = true;
        }
    }
}


// For Log function
template<typename TConstants, typename TValue>
inline void
addFactorToConstants(TConstants &constants, TValue &pC, TValue &pCm, TValue &pOther, unsigned &r, LogFunction const &)
{
    addFactorToConstants(constants, pC, pCm, pOther, r, Naive());
}

template<typename TConstantSet, typename TStrand>
inline void
adjustConstantsSize(TConstantSet &constantSet, TStrand strand, Naive const &)
{
    if (strand == 'F')
    {
        for (unsigned i = 0; i <= 1; ++i)
        {
            eraseBack(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')][i]);
            eraseBack(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')][i]);
            eraseBack(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')][i]);
            eraseBack(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')][i]);
        }
    }
    else 
    {
        for (unsigned i = 0; i <= 1; ++i)
        {
            eraseBack(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')][i]);
            eraseBack(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')][i]);
            eraseBack(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')][i]);
            eraseBack(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')][i]);
        }
    }
}

template<typename TConstantSet, typename TStrand>
inline void
adjustConstantsSize(TConstantSet &constantSet, TStrand strand, Polynomial const &)
{
    if (strand == 'F')
    {
        eraseBack(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')]);
        eraseBack(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')]);
        eraseBack(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')]);
        eraseBack(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')]);
    }
    else 
    {
        eraseBack(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')]);
        eraseBack(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')]);
        eraseBack(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')]);
        eraseBack(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')]);
    }
}

template<typename TConstantSet, typename TStrand>
inline void
adjustConstantsSize(TConstantSet &constantSet, TStrand strand, LogFunction const &)
{
    adjustConstantsSize(constantSet, strand, Naive());
}

template<typename TCoeffs>
inline TCoeffs polynomialMultipl(TCoeffs &coeffs1, TCoeffs &coeffs2)
{
    int n = length(coeffs1)-1;
    int m = length(coeffs2)-1;
    TCoeffs result;
    resize(result, n+m+1, 0.0, Exact());
    for (int k = 0; k <= m+n; ++k)
    {
        for (int i = 0; (i <= n && i <=k); ++i)
        {
            if (k-i <= m)                                   // clean up
                result[k] += coeffs1[i]*coeffs2[k-i];
        }
    }
    return result;
}


inline void
assignLHood(LogProbValue &lHood1, LogProbValue const &lHood2)
{
    lHood1 = lHood2;
}
template<typename TValue>
inline void
assignLHood(LogProbValue &lHood1, TValue const &lHood2)
{
    if (lHood2 == 0.0)
        lHood1.zero = true;
    else
        lHood1.c = lHood2;
}
template<typename TValue>
inline void
assignLHood(TValue &lHood1, TValue const &lHood2)
{
    lHood1 = lHood2;
}

template<typename TValue>
inline void
multiplyLHoods(TValue &lHood1, TValue const &lHood2)
{
    lHood1 *= lHood2;
}
inline void
multiplyLHoods(LogProbValue &lHood1, LogProbValue const &lHood2)
{
    if (!lHood1.zero && !lHood2.zero)
        lHood1.c *= lHood2.c;
    else
        lHood1.zero = true;
}
template<typename TValue>
inline void
multiplyLHoods(LogProbValue &lHood1, TValue const &lHood2)
{
    if (lHood2 == 0.0)
        lHood1.zero = true;
    else
        lHood1.c *= lHood2;
}




// Get coeffs for meth cases and (partial) likelihoods for other cases
template<typename TConstantSet, typename TLHoods, typename TQStrings, typename TMapqs, typename TOriginString, typename TCounts, typename TMethOptions, typename TMethod>
inline void 
constructConstantsAndLHoods(TConstantSet &constantSet, 
                            TLHoods &lHoods, 
                            TQStrings &qualF, TQStrings &qualR, 
                            TMapqs &mapqsF, TMapqs &mapqsR, 
                            TOriginString & originStringF, TOriginString & originStringR, 
                            TCounts &countF, TCounts &countR, 
                            TMethOptions &methOptions, 
                            TMethod const &)
{
    typedef typename Value<TLHoods>::Type   TLHood;
    unsigned minCountCT = 1;
    unsigned countF_CT = countF[ordValue((Dna)'C')] + countF[ordValue((Dna)'T')];
    unsigned countR_CT = countR[ordValue((Dna)'G')] + countR[ordValue((Dna)'A')];
    unsigned rF = 0;
    unsigned rR = 0;

    // ATTENTIONE: values for C are stored in constantSet[CG] and betas[CG]; values for G are stored in constantSet[GC] and betas[GC]

    // Get constants for reads on corresponding strand
    // For reads on other strand: caclulate already likelihood, since this is independent on beta
    // Multiply later
    // for each observed base type
    for (unsigned i = 0; i < 4; ++i)  
    {
        // for all reads mapped on forward strand
        for (unsigned j = 0; j < length(qualF[i]); ++j)
        {
            long double qual =  static_cast<long double>(ordValue(qualF[i][j])-33);
            /*
            if (candidatePos + startCoord == 908640 || candidatePos + startCoord == 985089 )  
            {
                std::cout << qual << std::endl;
            }
            */

            // If quality is below threshold, ignore read for all further calculations
            if (qual < 1 || (methOptions.useMapq && mapqsF[i][j] < 1))   
            {
                adjustConstantsSize(constantSet, 'F', TMethod());
                continue;
            }
            long double e = pow(10.0, (long double)(-qual/10.0));
            // for each possible candidate genotype
            // likelihood to observe single base under assumption of given genotype
            String<long double> singleProbs;
            resize(singleProbs, 6);
            for (unsigned h = 0; h < 6; ++h)
            {
                getSingleBaseProbHaploF(singleProbs[h], (Dna)i, (DnaM)h, e, originStringF[i][j], methOptions);
                if (methOptions.useMapq) 
                {
                    //singleProbs[h] *= mapqsF[i][j];
                    //std::cerr << "Use mapq : " << mapqsF[i][j]  << std::endl;
                    singleProbs[h] = mapqsF[i][j] * singleProbs[h];
                }
            }
         
            // Calculate likelihood for diploid genotypes
            for (unsigned h1 = 0; h1 < 4; ++h1)
            {
                for (unsigned h2 = h1; h2 < 4; ++h2)
                {
                    //std::cout << " test F" << (Dna)h1 << (Dna)h2 << " rF: " << rF << std::endl;
                    // Build up polynom coeffs
                    if ((Dna)h1 == 'C' && (Dna)h2 == 'C')   // CC
                    {
                        if (countF_CT >= minCountCT)
                        {
                            long double pC = singleProbs[ordValue((Dna)'C')];
                            long double pCm = singleProbs[ordValue((DnaM)'D')]; 
                            long double pOther = 0.0;
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rF, TMethod());
                        }
                    }
                    else if ( (Dna)h1 == 'C' && (Dna)h2 == 'G') // CG TODO same as CX
                    {
                        if ((countF_CT >= minCountCT) && (countR_CT >= minCountCT))
                        {
                            long double pC = singleProbs[ordValue((Dna)'C')];  // Beta independent of meth Level of Cs on other strand, hence to work with G is enough (?)
                            long double pCm = singleProbs[ordValue((DnaM)'D')];
                            long double pOther = singleProbs[ordValue((Dna)'G')];
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rF, TMethod());
                        }
                    }
                    else if ((Dna)h1 == 'C' || (Dna)h2 == 'C')  // CX
                    { 
                        if (countF_CT >= minCountCT)
                        {
                            long double pC = singleProbs[ordValue((Dna)'C')];
                            long double pCm = singleProbs[ordValue((DnaM)'D')];
                            long double pOther = singleProbs[((Dna)h1=='C')? h2:h1];
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rF, TMethod());
                        }
                    }
                    else            // XX, no beta to maximize
                    {
                        long double p = 0.5*singleProbs[h1] + 0.5*singleProbs[h2];
                        multiplyLHoods(lHoods[(h1<<2)|h2], p);
                    }
                }
            }
            ++rF;
        }
        // for all reads mapped on reverse strand
        for (unsigned j = 0; j < length(qualR[i]); ++j)
        {
            long double qual =  static_cast<long double>(ordValue(qualR[i][j])-33);
            // If quality is below threshold, ignore read for all further calculations
            if (qual < 1 || (methOptions.useMapq && mapqsR[i][j] < 1))   
            {
                adjustConstantsSize(constantSet, 'R', TMethod());
                continue;
            }

            long double e = pow(10.0, (long double)(-qual/10.0));  
            // likelihood to observe single base under assumption of given genotype
            String<long double> singleProbs;
            resize(singleProbs, 6);
            // for each possible candidate genotype
            for (unsigned h = 0; h < 6; ++h)
            {
                getSingleBaseProbHaploR(singleProbs[h], (Dna)i, (DnaM)h, e, originStringR[i][j], methOptions);
                if (methOptions.useMapq)
                    singleProbs[h] = mapqsR[i][j] * singleProbs[h];
            }
            
            // Calculate likelihood for diploid genotypes
            for (int h1 = 0; h1 < 4; ++h1)
            {
                for (int h2 = h1; h2 < 4; ++h2)
                {
                    // Build up polynom coeffs
                    if ( (Dna)h1 == 'C' && (Dna)h2 == 'G') // CG
                    {
                        if ( (countF_CT >= minCountCT) && (countR_CT >= minCountCT))
                        {
                            // G
                            long double pC = singleProbs[ordValue((Dna)'G')];  
                            long double pCm = singleProbs[ordValue((DnaM)'H')];
                            long double pOther = singleProbs[ordValue((DnaM)'C')];
                            addFactorToConstants(constantSet[(h2<<2)|h1], pC, pCm, pOther, rR, TMethod());
                        }
                    } 
                    else if ((Dna)h1 == 'G' && (Dna)h2 == 'G')  // GG
                    {
                        if (countR_CT >= minCountCT)
                        {
                            long double pC = singleProbs[ordValue((Dna)'G')];
                            long double pCm = singleProbs[ordValue((DnaM)'H')];
                            long double pOther = 0.0;
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rR, TMethod());
                        }
                    } 
                    else if ((Dna)h1 == 'G' || (Dna)h2 == 'G')   // GX
                    {
                        if (countR_CT >= minCountCT)
                        {
                            long double pC = singleProbs[ordValue((Dna)'G')];
                            long double pCm = singleProbs[ordValue((DnaM)'H')];
                            long double pOther = singleProbs[((Dna)h1=='G')? h2:h1];
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rR, TMethod());
                        }
                    } 
                    else            // XX, no beta to maximize
                    {
                        long double p = 0.5*singleProbs[h1] + 0.5*singleProbs[h2];
                        multiplyLHoods(lHoods[(h1<<2)|h2], p);
                    }
                }
            }
            ++rR;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////


// Check if borders are better (for newton method)
template<typename TLHood, typename TBeta, typename TFunctor, typename TEvalMethod>
inline TBeta 
verifyBeta(TLHood &lHood, TBeta &beta, TFunctor &functor, TEvalMethod const &)
{
    bool eBViolated = false;
    TLHood fBeta = functor(eBViolated, beta);
    if (eBViolated) fBeta = 0.0;

    TBeta b_0 = 0.0;
    TLHood f_0 = functor(eBViolated, b_0);
    if (eBViolated) f_0 = 0.0; 

    TBeta b_1 = 1.0;
    TLHood f_1 = functor(eBViolated, b_1);
    if (eBViolated) f_1 = 0.0;

    if (fBeta >= f_0 && fBeta >= f_1)
    {
        lHood = fBeta;
        return beta;
    }
    else if (f_1 >= f_0)
    {
        lHood = f_1;
        return 1.0;
    }
    else
    {
        lHood = f_0;
        return 0.0;
    }
}
// Check if borders are better (for newton method)
template<typename TBeta, typename TFunctor, typename TEvalMethod>
inline TBeta
verifyBeta(LogProbValue &lHood, TBeta &beta, TFunctor &functor, TEvalMethod const &)
{
    bool eBViolated = false;
    LogProbValue fBeta = functor(eBViolated, beta);
    if (eBViolated) fBeta.zero = true;

    TBeta b_0 = 0.0;
    LogProbValue f_0 = functor(eBViolated, b_0);
    if (eBViolated) f_0.zero = true; 

    TBeta b_1 = 1.0;
    LogProbValue f_1 = functor(eBViolated, b_1);
    if (eBViolated) f_1.zero = true;

    if (fBeta.zero && f_0.zero && f_1.zero)
    {
        lHood.zero = true;
        return 666.0;
    }
    else if (!fBeta.zero && (f_0.zero || fBeta.c >= f_0.c) && (f_1.zero || fBeta.c >= f_1.c))
    {
        lHood = fBeta;
        return beta;
    }
    else if (!f_1.zero && (f_0.zero || f_1.c >= f_0.c))
    {
        lHood = f_1;
        return 1.0;
    }
    else
    {
        lHood = f_0;
        return 0.0;
    }
}

// NSpace
template<typename TBeta, typename TLHood, typename TConstants, typename TMethOptions>
inline void
getMaximizingBeta(TBeta &beta, TLHood &lHood, TConstants &constants, TBeta /*&guess*/, TMethOptions &/*methOptions*/, Sampling const &, Naive const &)
{
#ifdef CALL_PROFILE 
    double timeStamp = sysTime();
#endif 
    FctNaive_0N<long double> fctNaive(constants);

    TLHood maxlHood = 0.0;
    TBeta betaMax = 0.0;
    for (TBeta currBeta = 0.0; currBeta <= 1.0; currBeta+=0.01)       // TODO: choose betas dependend guess and not equally distributed
    {
        TLHood currlHood =  fctNaive(currBeta);
        //std::cout << "curr Beta: " << currBeta << "  CurrlHood: " << currlHood << std::endl;
        if (currlHood > maxlHood)
        {
            betaMax = currBeta;
            maxlHood = currlHood;
        }
    } 
    beta = betaMax;
    lHood = maxlHood;
#ifdef CALL_PROFILE
    Times::instance().time_optimization += (sysTime() - timeStamp);
#endif
}
template<typename TBeta, typename TLHood, typename TConstants, typename TMethOptions>
inline void
getMaximizingBeta(TBeta &beta, TLHood &lHood, TConstants &constants, TBeta /*&guess*/, TMethOptions &/*methOptions*/, Sampling const &, Horner const &)
{
#ifdef CALL_PROFILE 
    double timeStamp = sysTime();
#endif 
    Fct_0eBN<long double> fctHorner(constants);
    bool eBViolated = false;
    TLHood maxlHood = 0.0;
    TBeta betaMax = 0.0;
    for (TBeta currBeta = 0.0; currBeta <= 1.0; currBeta+=0.01)       // TODO: choose betas dependend guess and not equally distributed
    {
        TLHood currlHood =  fctHorner(eBViolated, currBeta);
        if (!eBViolated && currlHood > maxlHood)
        {
            betaMax = currBeta;
            maxlHood = currlHood;
        }
    } 
    beta = betaMax;
    lHood = maxlHood;
#ifdef CALL_PROFILE
    Times::instance().time_optimization += (sysTime() - timeStamp);
#endif
}

// LogSpace
template<typename TBeta, typename TConstants, typename TMethOptions>
inline void
getMaximizingBeta(TBeta &beta, LogProbValue &lHood, TConstants &constants, TBeta /*&guess*/, TMethOptions &/*methOptions*/, Sampling const &, Naive const &)
{
#ifdef CALL_PROFILE 
    double timeStamp = sysTime();
#endif 

    FctNaive_0L<long double> fctNaive(constants);

    LogProbValue maxlHood;
    maxlHood.zero = true;
    TBeta betaMax = 0.0;
    for (TBeta currBeta = 0.0; currBeta <= 1.0; currBeta+=0.01)       // TODO: choose betas dependend guess and not equally distributed
    {
        LogProbValue currlHood = fctNaive(currBeta); 
        if ((!currlHood.zero && maxlHood.zero) || (!currlHood.zero && currlHood.c > maxlHood.c))
        {
            betaMax = currBeta;
            maxlHood = currlHood;
        }
    } 
    beta = betaMax;
    lHood = maxlHood;
#ifdef CALL_PROFILE
    Times::instance().time_optimization += (sysTime() - timeStamp);
#endif
}
template<typename TBeta, typename TConstants, typename TMethOptions>
inline void
getMaximizingBeta(TBeta &beta, LogProbValue &lHood, TConstants &constants, TBeta /*&guess*/, TMethOptions &/*methOptions*/, Sampling const &, Horner const &)
{
#ifdef CALL_PROFILE 
    double timeStamp = sysTime();
#endif 

    Fct_0eBL<long double> fctHorner(constants);
    bool eBViolated = false;
    LogProbValue maxlHood;
    maxlHood.zero = true;
    TBeta betaMax = 0.0;
    for (TBeta currBeta = 0.0; currBeta <= 1.0; currBeta+=0.01)       // TODO: choose betas dependend guess and not equally distributed
    {
        LogProbValue currlHood = fctHorner(eBViolated, currBeta); 
        if (!eBViolated && ((!currlHood.zero && maxlHood.zero) || (!currlHood.zero && currlHood.c > maxlHood.c)))
        {
            betaMax = currBeta;
            maxlHood = currlHood;
        }
    } 
    beta = betaMax;
    lHood = maxlHood;
#ifdef CALL_PROFILE
    Times::instance().time_optimization += (sysTime() - timeStamp);
#endif
}


// NSpace
template<typename TBeta, typename TLHood, typename TCoeffs, typename TMethOptions>
inline void
getMaximizingBeta(TBeta &beta, TLHood &lHood, TCoeffs &coeffs, TBeta &guess, TMethOptions &methOptions, Newton const &, Horner const &)
{
#ifdef CALL_PROFILE 
    double timeStamp = sysTime();
#endif 

    boost::uintmax_t maxIter = 10;
    //std::cout << " Guess: " << guess << std::endl;
    int digits = std::numeric_limits<double>::digits/2;
    Fct_12N<long double> myFN(coeffs);
    beta = boost::math::tools::newton_raphson_iterate(myFN, guess, (TBeta)0.0, (TBeta)1.0, digits, maxIter);
    // Check error bounds
    typedef typename boost::math::tuple<TBeta, TBeta> TTuple;   // get rid of boost tuple...
    Fct_12eBN<long double> myF(coeffs);
    bool eBViolatedR = false;
    TTuple tupleR = myF(eBViolatedR, beta);
    long double f_2R = boost::math::get<1>(tupleR); 
    if (eBViolatedR || f_2R >= 0)     // if eB of result is violated or if f'' not < 0
    {
        /*
        if (eBViolated)
            ++methOptions.counteBViolated;

        TTuple tupleG = myF(eBViolated, guess);
        if (eBViolated)         // if guess has already violated error bound, then stop! (probably we are curretnly not looking at the real genotype)
        {
            beta = 0.0;     //
            ++methOptions.countNoPlanB;
        }
        else
        {
            TBeta test_b = guess;
            TLHood l = 0.0;
            Fct_0eBN<long double> myF1(coeffs);
            verifyBeta(l, test_b, myF1, Horner());
            if (test_b > 0.1 || test_b < 0.9)     // if guess not at border (1 border: f'' always > 0)
            {
                // Do Sampling:
                //getMaximizingBeta(beta, lHood, coeffs, guess, methOptions, Sampling(), Horner());
                ++methOptions.countPlanB;
            }
        }
*/
        
        //std::cout << " Try bisektion: " << std::endl;
        bool eBViolatedG = false;
        TTuple tupleG = myF(eBViolatedG, guess);
        if (eBViolatedG)         // if guess is already wrong, then stop!
        {
            beta = 0.0;     // 
        }
        else
        {
            // TODO: if f'' > 0 check if borders have higher likelihood than beta and guess 
            if (!eBViolatedR && f_2R >= 0)
            {
                Fct_0eBN<long double> myF1(coeffs);
                TBeta beta1 = verifyBeta(lHood, beta, myF1, Horner());       // TODO just to test, clean up!
                TLHood lHood2;
                TBeta beta2 = verifyBeta(lHood2, guess, myF1, Horner());
                if (lHood < lHood2) // if guess is better than result, take guess
                {
                    beta = beta2;
                }
                else beta = beta1;

                if (beta < 0.01 || beta > 0.09)
                {
                    // its the border
                }
                //tupleR = myF(eBViolatedR, beta);
                //f_2R = boost::math::get<1>(tupleR);     // New f'', not: could be different than before
            }
            if ( !(!eBViolatedR && f_2R >= 0 && (beta < 0.01 || beta > 0.09)))  // Only apply bisektion if verified beta not at border
            {
                ++methOptions.countPlanB;

                long double f_1G = boost::math::get<0>(tupleG);
                long double f_2G = boost::math::get<1>(tupleG);
                if (f_1G < 0)       // Root between 0.0 and guess  
                {
                    if (beta < guess)   // eB violated at left border...
                    {
                        // Look for position between beta and guess, where eB is not violated, but sign diff still holds!
                        // Bisektion
                        beta = bisektionL(beta, guess, f_1G, f_2G, myF, myFN);
                    }
                    else 
                    {
                        // Check 0.0 border
                        // Bisektion
                       bool eBViolated = false;
                        /*TTuple tuple0 =*/ myF(eBViolated, (TBeta)0.0);
                        if (eBViolated)
                        {
                            beta = bisektionL((TBeta)0.0, guess, f_1G, f_2G, myF, myFN);
                        }
                        else
                            beta = boost::math::tools::newton_raphson_iterate(myFN, (TBeta)0.0, (TBeta)0.0, guess, digits, maxIter);   // Start with 0.0, with guess we tried already...
                        // Error bound should not be violated here, since no violation on left or right side
                        // Check anyway?
                    }     
                }
                else                // Root between guess and 1.0
                {
                    if (guess < beta)
                    {
                         // Look for position between guess and beta, where eB is not violated, but sign diff still holds!
                         // Bisektion
                         beta = bisektionR(guess, beta, f_1G, f_2G, myF, myFN);
                    }
                    else
                    {
                        // check 1.0 border
                        // Bisektion
                        bool eBViolated = false;
                        /*TTuple tuple1 =*/ myF(eBViolated, (TBeta)1.0);
                        if (eBViolated)
                        {
                            beta = bisektionR(guess, (TBeta)1.0, f_1G, f_2G, myF, myFN);
                        }
                        else
                            beta = boost::math::tools::newton_raphson_iterate(myFN, (TBeta)1.0, guess, (TBeta)1.0, digits, maxIter);   // Start with 0.0, with guess we tried already...
                        // Error bound should not be violated here, since no violation on left or right side
                        // Check anyway?
                    }
                }
            }
       } 
    }
    else ++methOptions.countNoPlanB;

    // TODO only if beta not already at border!
    // necessary here at all?
    Fct_0eBN<long double> myF1(coeffs);
    beta = verifyBeta(lHood, beta, myF1, Horner());

    if (methOptions.helpPrint)   // certain pos
    {
        std::cout << "Polynomial: " << std::endl;
        for (TBeta b =0.0; b <= 1.0; b+=0.05)
        {
            TTuple t12 = myFN(b);    // 12
            std::cout << "beta: " << b << "\t f: " << myF1(eBViolatedR, b) << "\t f': " << boost::math::get<0>(t12) << "\t f'': " <<  boost::math::get<1>(t12)  << std::endl;
        }
        std::cout << "beta: " << beta << " lHood: " << lHood << std::endl;
    }
#ifdef CALL_PROFILE
    Times::instance().time_optimization += (sysTime() - timeStamp);
#endif
}



// LogSpace
template<typename TBeta, typename TCoeffs, typename TMethOptions>
inline void
getMaximizingBeta(TBeta &beta, LogProbValue &lHood, TCoeffs &coeffs, TBeta &guess, TMethOptions &methOptions, Newton const &, Horner const &)
{
#ifdef CALL_PROFILE 
    double timeStamp = sysTime();
#endif 

    boost::uintmax_t maxIter = 10;
    //std::cout << " Guess: " << guess << std::endl;
    int digits = std::numeric_limits<long double>::digits;
    Fct_12L<long double> myFN(coeffs);
    beta = boost::math::tools::newton_raphson_iterate(myFN, guess, (TBeta)0.0, (TBeta)1.0, digits, maxIter);

    // Check error bounds
    typedef typename boost::math::tuple<TBeta, TBeta> TTuple;   // get rid of boost tuple...
    Fct_12eBL<long double> myF(coeffs);
    bool eBViolatedR = false;
    TTuple tupleR = myF(eBViolatedR, beta);
    long double f_2R = boost::math::get<1>(tupleR); 
    if (eBViolatedR || f_2R >= 0)     // if eB of result is violated or if f'' not < 0
    {
        bool eBViolatedG = false;
        TTuple tupleG = myF(eBViolatedG, guess);
        if (eBViolatedG)         // if guess is already wrong, then stop!
        {
            beta = 0.0;     // 
        }
        else
        {
            // TODO: if f'' > 0 check if borders have higher likelihood than beta and guess 
            if (!eBViolatedR && f_2R >= 0)
            {
                Fct_0eBL<long double> myF1(coeffs);
                TBeta beta1 = verifyBeta(lHood, beta, myF1, Horner());       // TODO just to test, clean up!
                LogProbValue lHood2;
                TBeta beta2 = verifyBeta(lHood2, guess, myF1, Horner());
                if (!lHood2.zero && !lHood.zero && lHood2.c > lHood.c) // if guess is better than result, take guess
                {
                    beta = beta2;
                }
                else beta = beta1;

                if (beta < 0.01 || beta > 0.09)
                {
                    // its the border
                }
                //tupleR = myF(eBViolatedR, beta);
                //f_2R = boost::math::get<1>(tupleR);     // New f'', not: could be different than before
            }
            if ( !(!eBViolatedR && f_2R >= 0 && (beta < 0.01 || beta > 0.09)))  // Only apply bisektion if verified beta not at border
            {
                ++methOptions.countPlanB;

                long double f_1G = boost::math::get<0>(tupleG);
                long double f_2G = boost::math::get<1>(tupleG);
                if (f_1G < 0)       // Root between 0.0 and guess  
                {
                    if (beta < guess)   // eB violated at left border...
                    {
                        // Look for position between beta and guess, where eB is not violated, but sign diff still holds!
                        // Bisektion
                        beta = bisektionL(beta, guess, f_1G, f_2G, myF, myFN);
                    }
                    else 
                    {
                        // Check 0.0 border
                        // Bisektion
                       bool eBViolated = false;
                        /*TTuple tuple0 = */ myF(eBViolated, (TBeta)0.0);
                        if (eBViolated)
                        {
                            beta = bisektionL((TBeta)0.0, guess, f_1G, f_2G, myF, myFN);
                        }
                        else
                            beta = boost::math::tools::newton_raphson_iterate(myFN, (TBeta)0.0, (TBeta)0.0, guess, digits, maxIter);   // Start with 0.0, with guess we tried already...
                        // Error bound should not be violated here, since no violation on left or right side
                        // Check anyway?
                    }     
                }
                else                // Root between guess and 1.0
                {
                    if (guess < beta)
                    {
                         // Look for position between guess and beta, where eB is not violated, but sign diff still holds!
                         // Bisektion
                         beta = bisektionR(guess, beta, f_1G, f_2G, myF, myFN);
                    }
                    else
                    {
                        // check 1.0 border
                        // Bisektion
                        bool eBViolated = false;
                        /*TTuple tuple1 =*/ myF(eBViolated, (TBeta)1.0);
                        if (eBViolated)
                        {
                            beta = bisektionR(guess, (TBeta)1.0, f_1G, f_2G, myF, myFN);
                        }
                        else
                            beta = boost::math::tools::newton_raphson_iterate(myFN, (TBeta)1.0, guess, (TBeta)1.0, digits, maxIter);   // Start with 0.0, with guess we tried already...
                        // Error bound should not be violated here, since no violation on left or right side
                        // Check anyway?
                    }
                }
            }
       } 
    }
    else ++methOptions.countNoPlanB;

    Fct_0eBL<long double> myF1(coeffs);
    beta = verifyBeta(lHood, beta, myF1, Horner());
#ifdef CALL_PROFILE
    Times::instance().time_optimization += (sysTime() - timeStamp);
#endif
}

// NSpace, Log function
template<typename TBeta, typename TLHood, typename TConstants, typename TMethOptions>
inline void
getMaximizingBeta(TBeta &beta, TLHood &lHood, TConstants &constants, TBeta &guess, TMethOptions &methOptions, Newton const &, LogFunction const &)
{
#ifdef CALL_PROFILE 
    double timeStamp = sysTime();
#endif 
    typedef typename boost::math::tuple<TBeta, TBeta> TTuple; 

    boost::uintmax_t maxIter = 10;
    //std::cout << " Guess: " << guess << std::endl;
    int digits = std::numeric_limits<long double>::digits;
    FctLog_12N<long double> myF12(constants);
    beta = boost::math::tools::newton_raphson_iterate(myF12, guess, (TBeta)0.0, (TBeta)1.0, digits, maxIter);
    // Check if f2 < 0 and get lHood
    FctLog_02N<long double> myF02(constants);
    TTuple tupleR = myF02(beta);
    long double f_2R = boost::math::get<1>(tupleR);
    if (f_2R >= 0)     // TODO Bisektion ?
    {
        getMaximizingBeta(beta, lHood, constants, guess, methOptions, Sampling(), Naive()); // TODO other plan b 
        ++methOptions.countPlanB;
    }
    else
    {
        TLHood logLHood = boost::math::get<0>(tupleR);
        lHood = pow(10, logLHood);
        ++methOptions.countNoPlanB;
    }

    if (methOptions.helpPrint)   // certain pos
    { 
        std::cout << "Log function: " << std::endl;
        for (TBeta b =0.0; b <= 1.0; b+=0.05)
        {
            TTuple t02 = myF02(b);
            TTuple t12 = myF12(b);
            std::cout << "beta: " << b << "\t f: " << boost::math::get<0>(t02) << "\t f': " << boost::math::get<0>(t12) << "\t f'': " <<  boost::math::get<1>(t12)  << std::endl;
        }
        std::cout << "beta: " << beta << " lHood: " << lHood << std::endl;
    }
#ifdef CALL_PROFILE
    Times::instance().time_optimization += (sysTime() - timeStamp);
#endif

}



template<typename TBetas, typename TLHoods, typename TConstantSet, typename TCounts, typename TMethOptions, typename TMethod, typename TEvalMethod>
inline void 
getBetasAndLHoods(TBetas &betas, TLHoods &lHoods, TConstantSet &constantSet, TCounts &countF, TCounts &countR, TMethOptions &methOptions, TMethod const &, TEvalMethod const &, unsigned &pos)
{
    typedef typename Value<TLHoods>::Type   TLHood;

    unsigned minCountCT = 1;    // We do not have to limit this here, since this is taken into account in prob. calculations later anyway
    resize(betas, 4*4, 0.666, Exact());   
    long double guess;   
    unsigned countF_CT = countF[ordValue((Dna)'C')] + countF[ordValue((Dna)'T')];
    unsigned countR_CT = countR[ordValue((Dna)'G')] + countR[ordValue((Dna)'A')];

    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            // Compute beta values and assign corresponding probs to lHoods
            if ( ((Dna)h1 == 'C' && (Dna)h2 == 'G')) // CG
            {
                if ((countF_CT >= minCountCT) && (countR_CT >= minCountCT))                // TODO C/T threshold to take meth level into account ?? problem: could bias score, result  
                {
                    // C: lHood from froward strand reads 
                    guess = (long double)countF[ordValue((Dna)'C')]/(long double)(countF_CT);
                    TLHood lHood;
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    assignLHood(lHoods[(h1<<2)|h2], lHood);
                    // G: lHood from reverse strand reads 
                    guess = (long double)countR[ordValue((Dna)'G')]/(long double)(countR_CT);
                    getMaximizingBeta(betas[h2<<2|h1], lHood, constantSet[(h2<<2)|h1], guess, methOptions, TMethod(), TEvalMethod());
                    multiplyLHoods(lHoods[(h1<<2)|h2], lHood);
                }
                else        // if not enough Cs and Ts: do not iterate and set prob to 0
                {
                    betas[h1<<2|h2] = 666.0;
                    betas[h2<<2|h1] = 666.0;
                    assignLHood(lHoods[(h1<<2)|h2], (long double)0.0);
               }
            }
            else if ( ((Dna)h1 == 'C' && (Dna)h2 == 'T'))
            {
                if (countF_CT >= minCountCT)
                {
                    int diff = (long double)countF[ordValue((Dna)'T')] - (long double)(countF_CT)/2;
                    if (diff > 0) guess = (long double)countF[ordValue((Dna)'C')]/(long double)(countF[ordValue((Dna)'C')] + diff);
                    else guess = 1.0;   // Set to 0.9 to avoid algorithm to get stuck?
                    TLHood lHood;
                    if (pos == 104)
                    {
                        methOptions.helpPrint = true; 
                        std::cout << (Dna)h1 << (Dna)h2 << std::endl;
                    }
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    methOptions.helpPrint = false; 
                    multiplyLHoods(lHoods[(h1<<2)|h2], lHood);
               }
                else        // if not enough Cs and Ts: do not iterate and set prob to 0
                {
                    betas[h1<<2|h2] = 666.0;
                    assignLHood(lHoods[(h1<<2)|h2], (long double)0.0);
                }
            }
            else if ( ((Dna)h1 == 'A' && (Dna)h2 == 'G'))
            {
                if (countR_CT >= minCountCT)
                {
                    int diff = (long double)countR[ordValue((Dna)'A')] - (long double)(countR_CT)/2;
                    if (diff > 0) guess = (long double)countR[ordValue((Dna)'G')]/(long double)(countR[ordValue((Dna)'G')] + diff);
                    else guess = 1.0;   // Set to 0.9 to avoid algorithm to get stuck?
                    TLHood lHood;
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    multiplyLHoods(lHoods[(h1<<2)|h2], lHood);
               }
                else
                {
                    betas[h1<<2|h2] = 666.0;
                    assignLHood(lHoods[(h1<<2)|h2], (long double)0.0);        
                }
            }
            else if ((Dna)h1 == 'C' || (Dna)h2 == 'C')
            {
                if (countF_CT >= minCountCT)
                {
                    guess = (long double)countF[ordValue((Dna)'C')]/(long double)(countF_CT);
                    TLHood lHood;
                    if (pos == 104 && (Dna)h1 == 'C' && (Dna)h2 == 'C')
                    {
                        methOptions.helpPrint = true; 
                        std::cout << (Dna)h1 << (Dna)h2 << std::endl;
                    }
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    methOptions.helpPrint = false; 
                    multiplyLHoods(lHoods[(h1<<2)|h2], lHood);
               }
                else
                {
                    betas[h1<<2|h2] = 666.0;
                    assignLHood(lHoods[(h1<<2)|h2], (long double)0.0);
                }
            }
            else if ((Dna)h1 == 'G' || (Dna)h2 == 'G')
            {
                if (countR_CT >= minCountCT)
                {
                    guess = (long double)countR[ordValue((Dna)'G')]/(long double)(countR_CT);
                    TLHood lHood;
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    multiplyLHoods(lHoods[(h1<<2)|h2], lHood);
              }
                else
                {
                    betas[h1<<2|h2] = 666.0;
                    assignLHood(lHoods[(h1<<2)|h2], (long double)0.0);
               }
            } 
        }
    }
}

template<typename TConstantSet, typename TLHoods, typename TCounts>
inline void
setUpConstants(TConstantSet &constantSet, TLHoods &lHoods, TCounts &countF, TCounts &countR, Polynomial const &)
{
    clear(lHoods);
    resize(lHoods, 4*4, 1.0, Exact()); 

    clear(constantSet);
    resize(constantSet, 4*4, Exact());
    unsigned covF = countF[0] + countF[1] + countF[2] + countF[3];
    unsigned covR = countR[0] + countR[1] + countR[2] + countR[3];

    // ATTENTIONE: values for C are stored in constantSet[CG] and betas[CG]; values for G are stored in constantSet[GC] and betas[GC]

    // Resize coeffs corresponding to read number taking into account
    // Forward strand methylations
    resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')], covF+1, 0.0, Exact());  
    resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')], covF+1, 0.0, Exact());  
    resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')], covF+1, 0.0, Exact());  
    resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')], covF+1, 0.0, Exact());  
    // Reverse strand methylations
    resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')], covR+1, 0.0, Exact()); 
    resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')], covR+1, 0.0, Exact());   
    resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')], covR+1, 0.0, Exact()); 
    resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')], covR+1, 0.0, Exact()); 
}
template<typename TCounts>
inline void
setUpConstants(String<String<LogProbValue> > &constantSet, String<LogProbValue> &lHoods, TCounts &countF, TCounts &countR, Polynomial const &)
{
    clear(lHoods);
    resize(lHoods, 4*4, Exact()); 

    clear(constantSet);
    resize(constantSet, 4*4, Exact());
    unsigned covF = countF[0] + countF[1] + countF[2] + countF[3];
    unsigned covR = countR[0] + countR[1] + countR[2] + countR[3];

    // ATTENTIONE: values for C are stored in constantSet[CG] and betas[CG]; values for G are stored in constantSet[GC] and betas[GC]

    // Resize coeffs corresponding to read number taking into account
    // Forward strand methylations
    resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')], covF+1, Exact());  
    resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')], covF+1, Exact());  
    resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')], covF+1, Exact());  
    resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')], covF+1, Exact());  
    // Reverse strand methylations
    resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')], covR+1, Exact()); 
    resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')], covR+1, Exact());   
    resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')], covR+1, Exact()); 
    resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')], covR+1, Exact()); 
}

template<typename TConstantSet, typename TLHoods, typename TCounts>
inline void
setUpConstants(TConstantSet &constantSet, TLHoods &lHoods, TCounts &countF, TCounts &countR, Naive const &)
{
    clear(lHoods);
    resize(lHoods, 4*4, 1.0, Exact()); 

    clear(constantSet);
    resize(constantSet, 4*4, Exact());
    unsigned covF = countF[0] + countF[1] + countF[2] + countF[3];
    unsigned covR = countR[0] + countR[1] + countR[2] + countR[3];

    // ATTENTIONE: values for C are stored in constantSet[CG] and betas[CG]; values for G are stored in constantSet[GC] and betas[GC]

    // resize for as and bs
    for (unsigned i = 0; i <= 1; ++i)
    {
        // Forward strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')], 2, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')], 2, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')], 2, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')], 2, Exact());  
        // Reverse strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')], 2, Exact()); 
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')], 2, Exact());   
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')], 2, Exact()); 
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')], 2, Exact()); 
    }
    // resize for number of reads
    for (unsigned j = 0; j <= 1; ++j)
    {
        // Forward strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')][j], covF, 0.0, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')][j], covF, 0.0, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')][j], covF, 0.0, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')][j], covF, 0.0, Exact());  
        // Reverse strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')][j], covR, 0.0, Exact()); 
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')][j], covR, 0.0, Exact());   
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')][j], covR, 0.0, Exact()); 
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')][j], covR, 0.0, Exact()); 
    }
}
template< typename TCounts>
inline void
setUpConstants(String<String<String<LogProbValue> > > &constantSet, String<LogProbValue> &lHoods, TCounts &countF, TCounts &countR, Naive const &)
{
    clear(lHoods);
    resize(lHoods, 4*4, Exact()); 

    clear(constantSet);
    resize(constantSet, 4*4, Exact());
    unsigned covF = countF[0] + countF[1] + countF[2] + countF[3];
    unsigned covR = countR[0] + countR[1] + countR[2] + countR[3];

    // ATTENTIONE: values for C are stored in constantSet[CG] and betas[CG]; values for G are stored in constantSet[GC] and betas[GC]

    // resize for as and bs
    for (unsigned i = 0; i <= 1; ++i)
    {
        // Forward strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')], 2, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')], 2, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')], 2, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')], 2, Exact());  
        // Reverse strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')], 2, Exact()); 
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')], 2, Exact());   
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')], 2, Exact()); 
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')], 2, Exact()); 
    }
    // resize for number of reads
    for (unsigned j = 0; j <= 1; ++j)
    {
        // Forward strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')][j], covF, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')][j], covF, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')][j], covF, Exact());  
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')][j], covF, Exact());  
        // Reverse strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')][j], covR, Exact()); 
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')][j], covR, Exact());   
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')][j], covR, Exact()); 
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')][j], covR, Exact()); 
    }
}

template<typename TConstantSet, typename TLHoods, typename TCounts>
inline void
setUpConstants(TConstantSet &constantSet, TLHoods &lHoods, TCounts &countF, TCounts &countR, LogFunction const &)
{
    setUpConstants(constantSet, lHoods, countF, countR, Naive());  
}

template<typename TPostProbs, typename TLHoods, typename TRefContext, typename TOptions, typename TMethOptions>
inline void
computePostProbs(TPostProbs &postProbs, TLHoods &lHoods, TRefContext &refContext, TOptions &options, TMethOptions &methOptions)
{
    // Calculate Pr(D)
    long double obsBasesProb = 0.0;
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            obsBasesProb += methOptions.genPriors[refContext.refAllele<<4| ordValue((Dna)h1)<<2| ordValue((Dna)h2)] * lHoods[(h1<<2)| h2];
         }
    }
    // Calculate posterior probs. for possible genotypes
    // Pr(G|D) = pi(G) * Pr(D|G) / Pr(D9)
    // For bs genotypes take prior methylation/non-methylation probability into account, dependent on context
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            if(options._debugLevel > 1)
                std::cout << (Dna)refContext.refAllele << (Dna)h1 << (Dna)h2 << "  " << std::setprecision (25) << "genPrior: " << std::setprecision (25) << (long double)methOptions.genPriors[ refContext.refAllele<<4| ordValue((Dna)((DnaM)h1))<<2 | ordValue((Dna)((DnaM)h2))] << "  lHood  " << (long double)lHoods[(h1<<2)| h2] << std::endl;
            
            // TODO take bsPriors into account, take context into account !!!
            postProbs[(h1<<2)| h2] = methOptions.genPriors[ refContext.refAllele<<4| h1<<2 | h2] * lHoods[(h1<<2)| h2] / obsBasesProb;   

            if(options._debugLevel > 1)
                 std::cout << std::setprecision (25) << "candidateProb: " << (Dna)refContext.refAllele << (Dna)h1 << (Dna)h2 << "  " << postProbs[(h1<<2)| h2] << "  context" << refContext.contextF << refContext.contextR << std::endl;
        }
    }
}

template<typename TPostProbs, typename TRefContext, typename TOptions, typename TMethOptions>
inline void
computePostProbs(TPostProbs &postProbs, String<LogProbValue> &lHoods, TRefContext &refContext, TOptions &options, TMethOptions &methOptions)
{
    // Calculate Pr(D)
    LogProbValue obsBasesProb;
    obsBasesProb.zero = true;
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            if (obsBasesProb.zero && !lHoods[(h1<<2)| h2].zero)
            {
                 obsBasesProb.c = lHoods[(h1<<2)| h2].c * methOptions.genPriors[refContext.refAllele<<4| ordValue((Dna)h1)<<2| ordValue((Dna)h2)];
                 obsBasesProb.zero = false;
            }
            else if (!lHoods[(h1<<2)| h2].zero)
                obsBasesProb.c += lHoods[(h1<<2)| h2].c * methOptions.genPriors[refContext.refAllele<<4| ordValue((Dna)h1)<<2| ordValue((Dna)h2)];
         }
    }
    SEQAN_ASSERT_NEQ(obsBasesProb.c, 0.0);
    // Calculate posterior probs. for possible genotypes
    // Pr(G|D) = pi(G) * Pr(D|G) / Pr(D9)
    // For bs genotypes take prior methylation/non-methylation probability into account, dependent on context -> how to deal with beta?
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            if(options._debugLevel > 1)
                std::cout << (Dna)refContext.refAllele << (Dna)h1 << (Dna)h2 << "  " << std::setprecision (25) << "genPrior: " << std::setprecision (25) << (long double)methOptions.genPriors[ refContext.refAllele<<4| ordValue((Dna)((DnaM)h1))<<2 | ordValue((Dna)((DnaM)h2))] << "  lHood  " << /*(long double)lHoods[(h1<<2)| h2] <<*/ std::endl;
            
            // TODO take bsPriors into account, take context into account !!!
            if (!lHoods[(h1<<2)| h2].zero)
                postProbs[(h1<<2)| h2] = (lHoods[(h1<<2)| h2].c * methOptions.genPriors[ refContext.refAllele<<4| h1<<2 | h2]) / obsBasesProb.c; 
            else
                postProbs[(h1<<2)| h2] = 0.0;

            if(options._debugLevel > 1)
                 std::cout << std::setprecision (25) << "candidateProb: " << (Dna)refContext.refAllele << (Dna)h1 << (Dna)h2 << "  " << postProbs[(h1<<2)| h2] << "  context" << refContext.contextF << refContext.contextR << std::endl;
        }
    }
}
 


template<typename TProbs, typename TBetas, typename TMethOptions, typename TOptions, typename TQStrings, typename TMapqs, typename TOriginString, typename TCounts, typename TRefContext>
inline void 
getCandidateProbs(TProbs &postProbs, TBetas &betas, 
                  TMethOptions &methOptions, TOptions &options, 
                  TQStrings &qualF, TQStrings &qualR, 
                  TMapqs &mapqsF, TMapqs &mapqsR, 
                  TOriginString & originStringF, 
                  TOriginString & originStringR,  
                  TCounts &countF, TCounts &countR, 
                  TRefContext &refContext)
{
    if (methOptions.ignoreBs) // Snp calling without bs conversions
        methOptions.convRate = 0.0;
 
    String<long double> lHoods;  // likelihoods to observe observed data under assumption of given genotypes 

    // Naive and Sampling Method
    if (methOptions.betaSampling && !methOptions.polynomialProbFunction && !methOptions.logSpace)           // sampling, naive, NSpace
    {
        String<long double> lHoods; 
        String<String<String<long double> > > constantSet;   
        setUpConstants(constantSet, lHoods, countF, countR, Naive());
        constructConstantsAndLHoods(constantSet, lHoods, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, methOptions, Naive());
        getBetasAndLHoods(betas, lHoods, constantSet, countF, countR, methOptions, Sampling(), Naive(), refContext.pos);
        computePostProbs(postProbs, lHoods, refContext, options, methOptions);
    }
    else if (methOptions.betaSampling && methOptions.polynomialProbFunction && !methOptions.logSpace)       // Sampling, polynomial, NSpace
    {
        String<long double> lHoods; 
        String<String<long double> > constantSet;    
        setUpConstants(constantSet, lHoods, countF, countR, Polynomial());
        constructConstantsAndLHoods(constantSet, lHoods, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, methOptions, Polynomial());
        getBetasAndLHoods(betas, lHoods, constantSet, countF, countR, methOptions, Sampling(), Horner(), refContext.pos);
        computePostProbs(postProbs, lHoods, refContext, options, methOptions);
    }
    else if (!methOptions.betaSampling && methOptions.polynomialProbFunction && methOptions.logSpace)   // Newton, polynomial, LogSpace
    {
        String<LogProbValue> lHoods;
        String<String<LogProbValue> > constantSet; 
        setUpConstants(constantSet, lHoods, countF, countR, Polynomial());
        constructConstantsAndLHoods(constantSet, lHoods, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, methOptions, Polynomial());
        getBetasAndLHoods(betas, lHoods, constantSet, countF, countR, methOptions, Newton(), Horner(), refContext.pos);
        computePostProbs(postProbs, lHoods, refContext, options, methOptions);
    }
    else if (!methOptions.betaSampling && methOptions.polynomialProbFunction && !methOptions.logSpace)   // Newton, polynomial, NSpace
    {
        String<long double> lHoods;
        String<String<long double> > constantSet; 
        setUpConstants(constantSet, lHoods, countF, countR, Polynomial());
        constructConstantsAndLHoods(constantSet, lHoods, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, methOptions, Polynomial());
        getBetasAndLHoods(betas, lHoods, constantSet, countF, countR, methOptions, Newton(), Horner(), refContext.pos);
        computePostProbs(postProbs, lHoods, refContext, options, methOptions);
    }
    else if (methOptions.betaSampling && !methOptions.polynomialProbFunction && methOptions.logSpace)   // sampling, naive, LogSpace
    {
        String<LogProbValue> lHoods;
        String<String<LogProbValue> > constantSet; 
        setUpConstants(constantSet, lHoods, countF, countR, Polynomial());
        constructConstantsAndLHoods(constantSet, lHoods, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, methOptions, Polynomial());
        getBetasAndLHoods(betas, lHoods, constantSet, countF, countR, methOptions, Sampling(), Horner(), refContext.pos);
        computePostProbs(postProbs, lHoods, refContext, options, methOptions);
    }
    else if (!methOptions.betaSampling && !methOptions.polynomialProbFunction && !methOptions.logSpace)   // Newton, logFunction, NSpace
    {
        //std::cout << " Run with LogFunction ...................................................!" << std::endl;
        String<long double> lHoods;
        String<String<String<long double> > > constantSet;  
        setUpConstants(constantSet, lHoods, countF, countR, LogFunction());
        constructConstantsAndLHoods(constantSet, lHoods, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, methOptions, LogFunction());
        getBetasAndLHoods(betas, lHoods, constantSet, countF, countR, methOptions, Newton(), LogFunction(), refContext.pos);
        computePostProbs(postProbs, lHoods, refContext, options, methOptions);
    }
}

///////////////////////////////////////////////////////////////////////
// SNP and meth calling in one step
template<typename TCounts, typename TQualities, typename TMapqs, typename TOriginString, typename TRefContext, typename TMethOptions, typename TOptions, typename TMethylVariant>
inline bool
doBsCalling(TCounts & countF,
          TCounts & countR,
          TQualities & qualF,
          TQualities & qualR,
          TMapqs & mapqsF,
          TMapqs & mapqsR,
          TOriginString & originStringF,
          TOriginString & originStringR,
          TRefContext & refContext,
          TMethOptions &methOptions,
          TOptions & options,
          TMethylVariant &meth
          )
{
    int genotypeRef = (refContext.refAllele<<2) | refContext.refAllele;

    String<long double> candidateProbs;
    String<long double> betas; 
    resize(candidateProbs, 4*4); // for simplicity; not all are used

    getCandidateProbs(candidateProbs, betas, methOptions, options, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, refContext);
    //Choose genotype which maximizes the posterior prob.

    int genotype1 = genotypeRef;
    int allele1 = 666;
    int allele2 = 666;
    int genotype2 = genotypeRef;
    long double maxProb1 = 0.0;
    long double maxProb2 = 0.0;
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            /*if (refContext.pos == 30032)
            {
                std::cout << std::setprecision (50) << "current genotype: " << (Dna)h1 << (Dna)h2 <<  "probs: " << candidateProbs[(h1<<2)| h2] << std::endl;
            }*/
            if (candidateProbs[(h1<<2)| h2] >= maxProb1) 
            {
                maxProb2 = maxProb1;
                genotype2 = (allele1<<2)|allele2;
                maxProb1 = candidateProbs[(h1<<2)| h2];
                allele1 = h1;
                allele2 = h2;
            }
            else if (candidateProbs[(h1<<2)| h2] >= maxProb2)
            {
                maxProb2 = candidateProbs[(h1<<2)| h2];
                genotype2 = (h1<<2)|h2;
            }
        }
    }
    genotype1 = (allele1<<2)|allele2;

    //unsigned totalCoverage = countF[0] + countF[1] +countF[2] +countF[3] +countF[4]
                           //+ countR[0] + countR[1] +countR[2] +countR[3] +countR[4];
    meth.genotype = genotype1;
    meth.score = log( (long double)candidateProbs[genotype1]/ (long double)candidateProbs[genotype2]);     // genotype calling score (only underlying genotypes with best beta, no bs types)?

    if (meth.score <= methOptions.minScoreToCallSnp)
    {
        meth.genotypeCalled = false;
        ++methOptions.countScoreTooLow;
    }
    else 
        meth.genotypeCalled = true;
    
    /*if (refContext.pos == 693)
    {
        std::cout << " Pos 693" << std::endl;
        std::cout << "Score: " << meth.score << " msc: " << methOptions.minScoreToCallSnp << std::endl;
        std::cout << " Called: " << (meth.genotypeCalled ? "true":"false") << std::endl;
    }*/

    meth.genotypeProb = candidateProbs[genotype1];
    meth.methLevel1 = betas[genotype1];
    if ((Dna)allele1 == 'C' && (Dna)allele2 == 'G')
        meth.methLevel2 = betas[(ordValue((Dna)'G')<<2)|ordValue((Dna)'C')];
    
    // If bs case:
    // Genotype was called = score was good enough to call
    if ( (meth.genotypeCalled && ((Dna)allele1 == 'C' || (Dna)allele1 == 'G' || (Dna)allele2 == 'C' || (Dna)allele2 == 'G')) )     // TODO: think how to hand call threshold
    {
        meth.bsCalled = true;
    }
    else 
        meth.bsCalled = false;


    /*if (refContext.pos == 30032 )
    {
        std::cout << "Pos: " << refContext.pos << std::endl;
        std::cout << std::setprecision (25) << " prob genotype1.." << (long double)candidateProbs[genotype1] << " prob genotype2.." << (long double)candidateProbs[genotype2]  <<  std::endl;
        std::cout << " genotype1  allele1: "<< (Dna)(genotype1>>2) << "  allele2: " << (Dna)(genotype1 % 4) << "  beta: " << betas[genotype1] << std::endl;
        std::cout << " genotype2  allele1: "<< (Dna)(genotype2>>2) << "  allele2: " << (Dna)(genotype2 % 4) << "  beta: " << betas[genotype2] << std::endl;
        std::cout << std::setprecision (50) << " prob genotype CG.." << (long double)candidateProbs[ordValue((Dna)'C')<<2|ordValue((Dna)'G')]  << "  beta: " << betas[ordValue((Dna)'C')<<2|ordValue((Dna)'G')] <<  std::endl;
        //std::cout << std::setprecision (25) << " prob genotype3.." << (long double)candidateProbs[0>>2|2]  <<  std::endl;
    }*/
 
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Output
////////////////////////////////////////////////////////////////////////////////////////////

// write to file
// Merge with other version?
template<typename TFile, typename TMethylVariant, typename TQualities, typename TRefContext, typename TMethOptions, typename TOptions>
inline bool
writeMeth(TFile &file,  
       TMethylVariant &meth,
       TQualities &qualityStringF, 
       TQualities &qualityStringR,
       TRefContext &refContext, 
       unsigned realCoverage,
       TMethOptions &methOptions,
       TOptions &options)
{
#ifdef CALL_PROFILE 
    double timeStamp = sysTime();
#endif 
//IOREV _nodoc_ what kind of format is this?
    if (!file.is_open()) 
    {
        ::std::cerr << "SNP/Meth output file is not open" << std::endl;
        return false;
    }
    //chromosome
    file << refContext.genomeID << '\t';
    //file << candPos + options.positionFormat<< '\t';
    file << refContext.pos << '\t';        // 0-based

    file << (Dna)refContext.refAllele <<'\t';

    if(options.showQualityStrings)
    {
        file << "["<<qualityStringF[0] <<"]\t";
        file << "["<<qualityStringF[1] <<"]\t";
        file << "["<<qualityStringF[2] <<"]\t";
        file << "["<<qualityStringF[3] <<"]\t";
        file << "["<<qualityStringR[0] <<"]\t";
        file << "["<<qualityStringR[1] <<"]\t";
        file << "["<<qualityStringR[2] <<"]\t";
        file << "["<<qualityStringR[3] <<"]\t";
    }

    file << realCoverage;

    // GenotypeCalled to string
    if(meth.genotypeCalled && ((meth.genotype>>2) != refContext.refAllele || (meth.genotype%4) != refContext.refAllele) )
        file << '\t' << (Dna)(meth.genotype>>2) << (Dna)(meth.genotype%4); 
    else
        file << "\t."; 
 
    // Methylation level
    if (meth.bsCalled)
    {
        file << '\t' << meth.methLevel1;
        if ((Dna)(meth.genotype>>2)  == 'C' && (Dna)(meth.genotype%4) == 'G')   // we don't know here anymore know if bad prob is maybe only cause by one strand
        {
            file << ':' << meth.methLevel2;
            // Some stats...
            if (refContext.contextF == 0) 
            {
                ++methOptions.countCG;
                methOptions.statsCGMethylated += meth.methLevel1; 
            }
            else if (refContext.contextF == 1) 
            {
                ++methOptions.countCHG;
                methOptions.statsCHGMethylated += meth.methLevel1; 
            }
            else
            {
                ++methOptions.countCHH;
                methOptions.statsCHHMethylated += meth.methLevel1; 
            }
            if (refContext.contextR == 0) 
            {
                ++methOptions.countCG;
                methOptions.statsCGMethylated += meth.methLevel2; 
            }
            else if (refContext.contextR == 1) 
            {
                ++methOptions.countCHG;
                methOptions.statsCHGMethylated += meth.methLevel2; 
            }
            else
            {
                ++methOptions.countCHH;
                methOptions.statsCHHMethylated += meth.methLevel2; 
            }
        }
        else if ((Dna)(meth.genotype>>2)  == 'C')
        {
            if (refContext.contextF == 0) 
            {
                ++methOptions.countCG;
                methOptions.statsCGMethylated += meth.methLevel1; 
            }
            else if (refContext.contextF == 1) 
            {
                ++methOptions.countCHG;
                methOptions.statsCHGMethylated += meth.methLevel1; 
            }
            else
            {
                ++methOptions.countCHH;
                methOptions.statsCHHMethylated += meth.methLevel1; 
            }
        } 
        else if ((Dna)(meth.genotype>>2)  == 'G')
        {
            if (refContext.contextR == 0) 
            {
                ++methOptions.countCG;
                methOptions.statsCGMethylated += meth.methLevel1; 
            }
            else if (refContext.contextR == 1) 
            {
                ++methOptions.countCHG;
                methOptions.statsCHGMethylated += meth.methLevel1; 
            }
            else
            {
                ++methOptions.countCHH;
                methOptions.statsCHHMethylated += meth.methLevel1; 
            }
        }
    }
    else
        file << "\t."; 

    file <<  '\t' << meth.genotypeProb << ':' << meth.score;

    if (meth.bsCalled)
    {
        file << ':';
        if (meth.methLevel1 >= 0.75)
            file << 'M';
        else if (meth.methLevel1 > 0.25)
            file << 'F';
        else
            file << 'U';
        if ((Dna)(meth.genotype>>2)  == 'C' && (Dna)(meth.genotype%4) == 'G')
        {
            if (meth.methLevel2 >= 0.75)
                file << 'M';
            else if (meth.methLevel2 > 0.25)
                file << 'F';
            else
                file << 'U';
        }
    }
    file << std::endl;
#ifdef CALL_PROFILE
    Times::instance().time_IO += (sysTime() - timeStamp);
#endif
    return true;
}


/*
template<typename TValue, typename F>
TValue newtonIterate(TValue prob, F functor, TValue guess, TValue min, TValue max, TValue thresholdErr, unsigned thresholdIt)
{
    //F p_functor(coeffs);
    std::cout << "start newtonIterate " << std::endl;

    TValue err = 1.0;
    unsigned countIt= 0;
    TValue x = guess;
    TValue x1;
    Tuple<TValue, 3> tuple;
    while( err > thresholdErr && countIt < thresholdIt )
    {
            tuple = functor(x);
            x1 = x - tuple[1]/tuple[2];             // We want to maximize f(x), so we deal with f'(x) and f''(x)
            err = fabs( x1 - x );
            x = x1;
            ++countIt;
    }

    // Check quick and dirty(!) if f(x) is maximum
    Tuple<TValue, 3> tuple0 = functor(0.0);
    Tuple<TValue, 3> tuple1 = functor(1.0);

    if (tuple[0] >= tuple0[0] && tuple[0] >= tuple1[0])
    {
        prob = tuple[0];
        return x;
    }
    else if (tuple0[0] >= tuple[0] && tuple0[0] >= tuple1[0])
    {
        prob = tuple0[0];
        return 0.0;
    }
    prob = tuple1[0];

    std::cout << "End newtonIterate " << std::endl;

    return 1.0;
}
// Functor for extended Horner function
// Returns values for f(x), f'(x) and f''(x)
template <typename TValue>
struct P_functor
{
    P_functor(String<TValue> const& coeffs) : coeffs(coeffs)
    { // Constructor 
    }
    Tuple<TValue, 3> operator()(TValue const&z)
    { // z is estimate so far.
        std::cout << "start calculate horner " << std::endl;

        TValue f_0 = back(coeffs);
        TValue f_1 = 0.0;
        TValue f_2 = 0.0;
        std::cout << "   1 " << std::endl;
        for (int i = length(coeffs)-2; i >= 0; --i)
        {
            std::cout << "   i " << i<< "  " <<  length(coeffs) << std::endl;

            f_2 = f_1 + z*f_2;
            f_1 = f_0 + z*f_1;
            f_0 = coeffs[i] + z*f_0;
        }
        std::cout << "   2" << std::endl;
        Tuple<TValue, 3> tuple;
        assignValue(tuple, 0, f_0);     // Returns value for Pr(Dj|G=...); no prior probs taken into account for maximization!
        assignValue(tuple, 1, f_1);
        assignValue(tuple, 2, f_2);
         std::cout << "end calculate horner " << std::endl;

        return tuple;
    }
private:
    String<TValue> coeffs;
};
*/

/*
 *template<typename TValue, typename TCoeffs>
TValue getProbForBeta_Der(TValue &beta, TCoeffs &coeffs, Horner const&)
{
        TValue f_0 = back(coeffs);
        TValue f_1 = 0.0;

        for (int i = length(coeffs)-2; i >= 0; --i)
        {
            f_1 = f_0 + beta*f_1;
            f_0 = coeffs[i] + beta*f_0;
        }
        return f_1;
}
*/


#endif  // #ifndef SANDBOX_KRAKAU_APPS_SNP_METH_STORE_BS_ONE_CALLING_H_
