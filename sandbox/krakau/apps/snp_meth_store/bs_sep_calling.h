/*==========================================================================

Methylation Calling: Version 2

==========================================================================*/

#ifndef SANDBOX_KRAKAU_APPS_SNP_METH_STORE_SNP_METH_STORE_H
#define SANDBOX_KRAKAU_APPS_SNP_METH_STORE_SNP_METH_STORE_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

#include "bs_alphabets.h"
#include "util.h"
//#include "snp_meth_store.h"
//#include "meths.h"


using namespace std;
using namespace seqan;

template<typename TProbs, typename TMethOptions, typename TOptions, typename TQStrings, typename TCounts>
void
getCandidatePostProbs(TProbs &candidatePostProbs, TMethOptions &methOptions, TOptions &options, TQStrings &qualF, TQStrings &qualR, TCounts &countF, TCounts &countR, int &refAllele)
{
    // methylation probability (locus-specific for the beginning)
    long double betaF = 0.0;
    if (countF[ordValue((Dna)'C')] + countF[ordValue((Dna)'T')] > 0) betaF = (long double)countF[ordValue((Dna)'C')] /(long double)(countF[ordValue((Dna)'C')] + countF[ordValue((Dna)'T')]);
    long double betaR = 0.0;
    if (countR[ordValue((Dna)'G')] + countF[ordValue((Dna)'A')] > 0) betaR = (long double)countR[ordValue((Dna)'G')] /(long double)(countR[ordValue((Dna)'G')] + countR[ordValue((Dna)'A')]);

    if (methOptions.ignoreBs) // Snp calling without bs conversations
    {
        betaF = 1.0;
        betaR = 1.0;
        methOptions.convRate = 0.0;
    }


    if(options._debugLevel > 1)
    {
        std::cout << std::setprecision (25) << "betaR: "<< betaR << " countR G " << countR[ordValue((Dna)'G')] << "   countR A " << countR[ordValue((Dna)'A')] << std::endl;
        std::cout << std::setprecision (25) << "betaR: "<< betaR << " countR G " << countR[2] << "   countR A " << countR[0] << std::endl;
    }
	long double e; // error probability represented by quals

    // likelihood to observe single base under assumption of given genotype
    String<long double> singleProbs;
    resize(singleProbs, 4);

    String<long double> lHoods;  // likelihoods to observe observed data under assumption of given genotypes 
    resize(lHoods, 4*4, 1.0, Exact());
    //FunctorComplement<Dna> f;
    // for each possible observed base type
    for (unsigned i = 0; i < 4; ++i)
    {
        if(options._debugLevel > 1)
            std::cout << " base: "<< (Dna)i << std::endl;  

        // for all reads mapped on forward strand
        for (unsigned j = 0; j < length(qualF[i]); ++j)
        {
            double qual =  static_cast<double>(ordValue(qualF[i][j])-33);
            e = pow(10.0, (double)(-qual/10.0));
            if(options._debugLevel > 1)
                std::cout <<  "  e  " << e << " qual " << (double)getValue(qualF[i], j) <<  "  betaF "  << betaF  <<std::endl;   
            // for each possible candidate genotype
            for (unsigned g = 0; g < 4; ++g)
            {
                if ( g == ordValue((Dna)'C'))          // Genotype C
                {                     
                    if (i == ordValue((Dna)'C'))           // forward strand C (correct and meth. and not converted or unmeth. and not converted)
                        singleProbs[g] = (1.0-e)*(betaF + (1.0-betaF)*(1.0-methOptions.convRate));
                    else if (i == ordValue((Dna)'T'))      // forward strand T (base call error or correct and unmeth. and converted)
                        singleProbs[g] = e/3.0 + (1.0 - e)*((1.0-betaF)*methOptions.convRate);
                    else                                    // base call error
                        singleProbs[g] = e/3.0;
                } 
                else if ( g == ordValue((Dna)'G'))     // Genotype G
                {
                    if (i == ordValue((Dna)'G'))           // forward strand correct ('G')
                        singleProbs[g] = 1.0 - e;
                    else                                    // forward strand base call error    
                        singleProbs[g] = e/3.0;
                } 
                else                                    // Genotype T, A: no bs 
                {
                    if (i != g)
                        singleProbs[g] = e/3.0;
                    else
                        singleProbs[g] = 1.0 - e;
                }
            }
            // Calculate likelihood for diploid genotypes
            for (int g1 = 0; g1 < 4; ++g1)
            {
                for (int g2 = g1; g2 < 4; ++g2)
                {
                    lHoods[(g1<<2)|g2] *= (long double)(0.5*singleProbs[g1] + 0.5*singleProbs[g2]);
                    if(options._debugLevel > 1)
                    {    
                        if ((Dna)g1 == 'G' && (Dna)g2 == 'G')
                        {
                            std::cout << "Forward:"<< std::endl;
                            //std::cout << std::setprecision (25) << "p(g1): "<< (Dna)g1 << " " << singleProbs[g1] << "   p(g2): " << (Dna)g2 << " " << singleProbs[g1] << std::endl;
                            std::cout << std::setprecision (25) <<  "  my prob: " << (Dna)g1 << (Dna)g2 << " " << (0.5*singleProbs[g1] + 0.5*singleProbs[g2])  << std::endl;
                        }
                    }
                }
            }
        }
        // for all reads mapped on reverse strand
        for (unsigned j = 0; j < length(qualR[i]); ++j)
        {
            double qual =  static_cast<double>(ordValue(qualR[i][j])-33);
            e = pow(10.0, (double)(-qual/10.0));
            if(options._debugLevel > 1)
                std::cout <<  "  e  " << e  << " qual " << (double)getValue(qualR[i], j) <<  "  betaR "  << betaR  <<std::endl;   
            // for each possible candidate genotype
            for (unsigned g = 0; g < 4; ++g)
            {
                if ( g == ordValue((Dna)'C'))          // Genotype C
                { 
                    if (i == ordValue((Dna)'C'))           // reverse strand correct ('G')
                        singleProbs[g] = 1.0 - e; 
                    else                                    // reverse strand base call error                    
                        singleProbs[g] = e/3.0;                    
                }
                else if ( g == ordValue((Dna)'G'))     // Genotype G
                {
                    if (i == ordValue((Dna)'G'))           // reverse strand C (correct and meth. and not converted or unmeth. and not converted)
                        singleProbs[g] = (1.0 - e)*(betaR + (1.0-betaR)*(1.0-methOptions.convRate)); 
                    else if (i == ordValue((Dna)'A'))      // reverse strand T (base call error or correct and unmeth. and converted)
                         singleProbs[g] = e/3.0 + (1.0 - e)*((1.0-betaR)*methOptions.convRate);  
                    else                                    // base call error
                        singleProbs[g] = e/3.0; 
                } 
                else                                    // Genotype T, A: no bs
                {
                    if (i != g)
                        singleProbs[g] = e/3.0;    
                    else
                        singleProbs[g] = 1.0 - e;    
                }
            }
              
            // Calculate likelihood for diploid genotypes
            for (int g1 = 0; g1 < 4; ++g1)
            {
                for (int g2 = g1; g2 < 4; ++g2)
                {
                    lHoods[(g1<<2)|g2] *=(long double)(0.5*singleProbs[g1] + 0.5*singleProbs[g2]);
                    if(options._debugLevel > 1)
                    {
                        if ((Dna)g1 == 'G' && (Dna)g2 == 'G')
                        {
                            std::cout << "Reverse:" << std::endl;
                            //std::cout << std::setprecision (25) << "p(g1): "<< (Dna)g1 << " " << singleProbs[g1] << "   p(g2): " << (Dna)g2 << " " << singleProbs[g1] << std::endl;
                            std::cout << std::setprecision (25) << "  my prob: " << (Dna)g1 << (Dna)g2 << " " << (long double)(0.5*singleProbs[g1] + 0.5*singleProbs[g2]) << std::endl;
                        }
                    }                    
                }
            }
        }
    }
  
    // Calculate Pr(D)
    long double obsBasesProb = 0.0;
    for (int g1 = 0; g1 < 4; ++g1)
    {
        for (int g2 = g1; g2 < 4; ++g2)
        {
            obsBasesProb += methOptions.genPriors[ refAllele<<4| g1<<2| g2] * lHoods[(g1<<2)| g2];
        }
    }

    // Calculate posterior probs. for possible genotypes
    // Pr(G|D) = pi(G) * Pr(D|G) / Pr(D)
    for (int g1 = 0; g1 < 4; ++g1)
    {
        for (int g2 = g1; g2 < 4; ++g2)
        {
            if(options._debugLevel > 1)
                std::cout << (Dna)refAllele << (Dna)g1 << (Dna)g2 << "  " << std::setprecision (25) << "genPrior: " << std::setprecision (25) << (long double)methOptions.genPriors[ refAllele<<4| g1<<2| g2] << "  lHood  " << (long double)lHoods[(g1<<2)| g2] << std::endl;
            candidatePostProbs[(g1<<2)| g2] = methOptions.genPriors[ refAllele<<4| g1<<2| g2] * lHoods[(g1<<2)| g2] / obsBasesProb; 
            if(options._debugLevel > 1)
                 std::cout << std::setprecision (25) << "candidateProb: " << (Dna)refAllele << (Dna)g1 << (Dna)g2 << "  " << candidatePostProbs[(g1<<2)| g2] << std::endl;
        }
    }

}

// Calculate most reasonable state
// by looking at C-T level, Threshold
template<typename TMethylVariant, typename TCounts, typename TMethOptions>
bool
getMethState(TMethylVariant &meth, TCounts &countF, TCounts &countR, int &refAllele, TMethOptions &methOptions)
{ 
    int minCoverageSingle = 3.0;
    double offset = 0.25; 

    // Compute C-T ratios (methylation level)
    if ( ((Dna)refAllele == 'C' && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'C' || (Dna)(meth.genotype%4) == 'C')) )             // Methylation on 'F' ?
    {
        if (meth.coverageF >= minCoverageSingle && meth.coverageR >= minCoverageSingle
                   && countF[ordValue((Dna)'C')] > 0 )
        {
            meth.methLevelF = (double)countF[ordValue((Dna)'C')] / (double)(countF[ordValue((Dna)'C')]+countF[ordValue((Dna)'T')]); // oder so ähnlich...  
            //std::cout << "methLevelF: " << meth.methLevelF << " count C " << countF[ordValue((Dna)'C')] << " count T " << countF[ordValue((Dna)'T')] << std::endl;
        }
    }
    if ( ((Dna)refAllele == 'G' && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' || (Dna)(meth.genotype%4) == 'G')) )          // Methylation on 'R' ?
    {
         if (meth.coverageF >= minCoverageSingle && meth.coverageR >= minCoverageSingle
                    && countR[ordValue((Dna)'G')] > 0 )
        {
            meth.methLevelR = (double)countR[ordValue((Dna)'G')] / (double)(countR[ordValue((Dna)'G')]+countR[ordValue((Dna)'A')]); // oder so ähnlich...                
        }
    }

    double heteroR = 0.5;       // better way? Heinricht et. all? or approximate 

    // Call methylation state dependent on C-T ratios for different scenarios
    // if hete snp is called, with CG, state1 represents methylation state on F, state2 on R
    // same with called1, called2
    //
    // make sure that ranges do not overlap
    // TODO: check for actual number of reads instead of doubles...?
    // CC, GG, CH, GH
    SEQAN_ASSERT_LT((1-methOptions.convRate)+offset, heteroR*(1-offset));
    SEQAN_ASSERT_LT(heteroR*(1-offset), 1 - heteroR*(methOptions.convRate-offset));
    SEQAN_ASSERT_LT(1 - heteroR*(methOptions.convRate-offset), 1 - offset);
    // CT, AG
    SEQAN_ASSERT_LT(heteroR*(1-methOptions.convRate) + offset, heteroR*(1-offset));

    //std::cout << " ranges to call methylation state: " << std::endl;
    //std::cout << (1-methOptions.convRate)+offset << '\t' << heteroR*(1-offset) << '\t' << 1 - heteroR*(methOptions.convRate-offset) << '\t' <<  1 - offset << std::endl;

    if ( (meth.snpCalled && ((meth.genotype>>2) == ordValue((Dna)'C')) && ((meth.genotype % 4) == ordValue((Dna)'C')) ) ||    // CC
         (!meth.snpCalled && (refAllele == ordValue((Dna)'C')) ) )          
    {
        if (meth.methLevelF <= (1-methOptions.convRate) + offset )
        {
            meth.state1 = false;
            meth.state2 = false;
            meth.called1 = true;
            meth.called2 = true;
        }
        else if (meth.methLevelF >= 1 - offset)  // both Cs are methylated -> not converted
        {
            meth.state1 = true;
            meth.state2 = true;
            meth.called1 = true;
            meth.called2 = true;
        }
        else if ( (meth.methLevelF >= heteroR*(1-offset) ) &&   // one C is methylated
                  (meth.methLevelF <= 1 - heteroR*(methOptions.convRate-offset) )  )
        {
            meth.state1 = true;
            meth.state2 = false;
            meth.called1 = true;
            meth.called2 = true;
        }
    }
    else if ( (meth.snpCalled && ((meth.genotype>>2) == ordValue((Dna)'G')) && ((meth.genotype % 4) == ordValue((Dna)'G')) ) ||    // GG
         (!meth.snpCalled && (refAllele == ordValue((Dna)'G')) ) )  
    {
        if (meth.methLevelR <= (1-methOptions.convRate) + offset )     // both Cs are not methylated -> converted
        {                               // if no. of Cs is equal or smaller than the non-conversion rate + offset 
            meth.state1 = false;
            meth.state2 = false;
            meth.called1 = true;
            meth.called2 = true;
        }
        else if (meth.methLevelR >= 1 - offset)  // both Cs are methylated -> not converted
        {
            meth.state1 = true;
            meth.state2 = true;
            meth.called1 = true;
            meth.called2 = true;
        }
        else if ( (meth.methLevelR >= heteroR*(1-offset) ) &&   // one C is methylated
                  (meth.methLevelR <= 1 - heteroR*(methOptions.convRate-offset) )  )
        {
            meth.state1 = true;
            meth.state2 = false;
            meth.called1 = true;
            meth.called2 = true;
        }
    }
    else 
    {
        if ( (meth.genotype>>2 == ordValue((Dna)'C') && (meth.genotype % 4 == ordValue((Dna)'T')) ) ||      // CT || TC  
             (meth.genotype>>2 == ordValue((Dna)'T') && (meth.genotype % 4 == ordValue((Dna)'C')) ) )       // should only be CT and not TC anyway ...
        {
            if (meth.methLevelF <= heteroR*(1-methOptions.convRate) + offset)   // C is not methylated -> converted    // + error Prob. that allele2 is C ?? 
            {
                meth.state1 = false;
                meth.called1 = true;
            }
            else if (meth.methLevelF >= heteroR*(1-offset) )        // C is methylated -> not converted
            {
                meth.state1 = true;
                meth.called1 = true;
            }
        }
        else if ( (meth.genotype>>2 == ordValue((Dna)'C') || (meth.genotype % 4 == ordValue((Dna)'C')) ) )       // CH || HC
        {                                                                                                       // merging with case CC possible?
            if (meth.methLevelF <= (1-methOptions.convRate) + offset ) //
            {
                meth.state1 = false;
                meth.called1 = true;
            }
            else if (meth.methLevelF >= 1 - offset)
            {
                meth.state1 = true;
                meth.called1 = true;
            }
        }
        if ( (meth.genotype>>2 == ordValue((Dna)'G') && (meth.genotype % 4 == ordValue((Dna)'A')) ) ||      // AG || GA   
             (meth.genotype>>2 == ordValue((Dna)'A') && (meth.genotype % 4 == ordValue((Dna)'G')) ) )   
        {
            if (meth.methLevelR <= heteroR*(1-methOptions.convRate) + offset)   // C is not methylated -> converted    // + error Prob. that allele2 is C ?? 
            {
                meth.state2 = false;
                meth.called2 = true;
            }
            else if (meth.methLevelR >= heteroR*(1-offset) )        // C is methylated -> not converted
            {
                meth.state2 = true;
                meth.called2 = true;
            }
        }
        else if ( (meth.genotype>>2 == ordValue((Dna)'G') || (meth.genotype % 4 == ordValue((Dna)'G')) ) ) // GH || HG
        {
            if (meth.methLevelR <= (1-methOptions.convRate) + offset ) //
            {
                meth.state2 = false;
                meth.called2 = true;
            }
            else if (meth.methLevelR >= 1 - offset)
            {
                meth.state2 = true;
                meth.called2 = true;
            }
        }  
    }
    return 0;
}


// Calculate meth state probs/likelihoods (Bayesian likelihood model)
template<typename TMethylVariant, typename TQualities, typename TCounts, typename TMethOptions, typename TOptions>
bool
getMethProb(TMethylVariant &meth, TQualities &qualF, TQualities &qualR, TCounts &/*countF*/, TCounts &/*countR*/, int &refAllele, unsigned &contextF, unsigned &contextR, TMethOptions &methOptions, TOptions &options)
{
	long double e; // error probability represented by quals

    // likelihood to observe single base under assumption of given genotype
    String<long double> singleProbs;
    resize(singleProbs, 2);

    String<long double> lHoods;  // likelihoods to observe observed data under assumption of given genotypes 
    resize(lHoods, 3, 1.0, Exact());

   
    if ( (refAllele == ordValue((Dna)'C') && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'C' && (Dna)(meth.genotype%4) == 'C')) )  // *** CC ***
    {
        // for each possible observed read base
        for (unsigned i = 0; i < 4; ++i) 
        {
            // for all reads mapped on forward strand
            for (unsigned j = 0; j < length(qualF[i]); ++j)
            {
                double qual =  static_cast<double>(ordValue(qualF[i][j])-33);
                e = pow(10.0, (double)(-qual/10.0));
                // for methylated genotype            
                if (i == ordValue((Dna)'C'))           // forward strand C                
                     singleProbs[0] = (1.0-e);
                else                                    // base call error
                    singleProbs[0] = e/3.0;    
                // for non-methylated genotype
                if (i == ordValue((Dna)'C'))           // forward strand C                
                    singleProbs[1] = (1.0-e)*(1-methOptions.convRate);
                else if (i == ordValue((Dna)'T'))      // forward strand T 
                      singleProbs[1] = e/3.0 + (1.0-e)*methOptions.convRate;
                else                                    // base call error
                    singleProbs[1] = e/3.0;

                // Calculate likelihood for diploid genotypes
                lHoods[0] *= (long double)(singleProbs[0]);    // homoMeth
                lHoods[1] *= (long double)(singleProbs[1]);    // homoNonMeth
                lHoods[2] *= (long double)(0.5*singleProbs[0] + 0.5*singleProbs[1]);    // heteroMeth
           }
            // for all reads mapped on reverse strand
            // TODO: maybe better ignore? because here we assume genotype is called right, and just look at methylation probability...???
            // errors on other side are represented in genotype score ? or discard position maybe?
            for (unsigned j = 0; j < length(qualR[i]); ++j)
           {
                double qual =  static_cast<double>(ordValue(qualR[i][j])-33);
                e = pow(10.0, (double)(-qual/10.0));                 
                // for methylated genotype            
                if (i == ordValue((Dna)'G'))           // reverse strand C                
                    singleProbs[0] = (1.0-e);
                else                                    // base call error
                    singleProbs[0] = e/3.0;
                // for non-methylated genotype
                if (i == ordValue((Dna)'G'))           // reverse strand C                
                    singleProbs[1] = (1.0-e);
                else                                    // base call error
                    singleProbs[1] = e/3.0;

                // Calculate likelihood for diploid genotypes
                lHoods[0] *= (long double)singleProbs[0];                               // homoMeth
                lHoods[1] *= (long double)singleProbs[1];                               // homoNonMeth
                lHoods[2] *= (long double)(0.5*singleProbs[0] + 0.5*singleProbs[1]);    // heteroMeth
            }
        }

        long double obsBasesProb;
        if (contextF == 0) 
        {
            // Calculate Pr(D)
            obsBasesProb = methOptions.homoMethPriorCG*lHoods[0] + methOptions.homoNonMethPriorCG*lHoods[1] + methOptions.heteroMethPriorCG*lHoods[2];    
           // Calculate posterior probs. for possible genotypes
           // Pr(G|D) = pi(G) * Pr(D|G) / Pr(D)
            meth.homoMethProbF = methOptions.homoMethPriorCG*lHoods[0] / obsBasesProb;
            meth.homoNonMethProbF = methOptions.homoNonMethPriorCG*lHoods[1] / obsBasesProb;
            meth.heteroMethProbF = methOptions.heteroMethPriorCG*lHoods[2] / obsBasesProb;
        }
        else if (contextF == 1)
        {
            obsBasesProb = methOptions.homoMethPriorCHG*lHoods[0] + methOptions.homoNonMethPriorCHG*lHoods[1] + methOptions.heteroMethPriorCHG*lHoods[2];    
            meth.homoMethProbF = methOptions.homoMethPriorCHG*lHoods[0] / obsBasesProb;
            meth.homoNonMethProbF = methOptions.homoNonMethPriorCHG*lHoods[1] / obsBasesProb;
            meth.heteroMethProbF = methOptions.heteroMethPriorCHG*lHoods[2] / obsBasesProb;
        }
        else 
        {
            obsBasesProb = methOptions.homoMethPriorCHH*lHoods[0] + methOptions.homoNonMethPriorCHH*lHoods[1] + methOptions.heteroMethPriorCHH*lHoods[2];    
            meth.homoMethProbF = methOptions.homoMethPriorCHH*lHoods[0] / obsBasesProb;
            meth.homoNonMethProbF = methOptions.homoNonMethPriorCHH*lHoods[1] / obsBasesProb;
            meth.heteroMethProbF = methOptions.heteroMethPriorCHH*lHoods[2] / obsBasesProb;
        }
        if(options._debugLevel > 1)
        {
            std::cout << std::endl;
            std::cout << "obsBasesProb: " << obsBasesProb << std::endl;
            std::cout << "Methylated: lHoods  " << lHoods[0] << "  posterior prob:  " << meth.homoMethProbF << std::endl;
            std::cout << "Non-Methylated: lHoods  " << lHoods[1] << "  posterior prob:  " << meth.homoNonMethProbF << std::endl;
            std::cout << "Hetero-Methylated: lHoods  " << lHoods[2] << "  posterior prob:  " << meth.heteroMethProbF << std::endl;
            std::cout << std::endl;
        }
    } 
    else if ( (refAllele == ordValue((Dna)'G') && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' && (Dna)(meth.genotype%4) == 'G')) )  // *** GG ***
    {
        // for each possible observed read base
        for (unsigned i = 0; i < 4; ++i)
        {
            // for all reads mapped on forward strand
            for (unsigned j = 0; j < length(qualF[i]); ++j)
            {
                double qual =  static_cast<double>(ordValue(qualF[i][j])-33);
                e = pow(10.0, (double)(-qual/10.0));
                // for methylated genotype            
                if (i == ordValue((Dna)'G'))           // forward strand G                
                     singleProbs[0] = (1.0-e);
                else                                    // base call error
                    singleProbs[0] = e/3.0;    
                // for non-methylated genotype
                if (i == ordValue((Dna)'G'))           // forward strand G                
                    singleProbs[1] = (1.0-e);
                else                                    // base call error
                    singleProbs[1] = e/3.0;

                // Calculate likelihood for diploid genotypes
                lHoods[0] *= (long double)singleProbs[0];    // homoMeth
                lHoods[1] *= (long double)singleProbs[1];    // homoNonMeth
                lHoods[2] *= (long double)(0.5*singleProbs[0] + 0.5*singleProbs[1]);    // heteroMeth
           }
            // for all reads mapped on reverse strand
            for (unsigned j = 0; j < length(qualR[i]); ++j)
           {
                e = pow(10.0, (double)(-getValue(qualR[i], j)/10.0));                 
                // for methylated genotype            
                if (i == ordValue((Dna)'G'))           // reverse strand C                
                    singleProbs[0] = (1.0-e);
                else                                    // base call error
                    singleProbs[0] = e/3.0;
                // for non-methylated genotype
                if (i == ordValue((Dna)'G'))           // reverse strand C                
                    singleProbs[1] = (1.0-e)*(1.0-methOptions.convRate);
                else if (i == ordValue((Dna)'A'))      // reverse strand T
                    singleProbs[1] = e/3.0 + (1.0-e)*methOptions.convRate;
                else                                    // base call error
                    singleProbs[1] = e/3.0;

                // Calculate likelihood for diploid genotypes
                lHoods[0] *= (long double)singleProbs[0];    // homoMeth
                lHoods[1] *= (long double)singleProbs[1];    // homoNonMeth
                lHoods[2] *= (long double)(0.5*singleProbs[0] + 0.5*singleProbs[1]);    // heteroMeth
            }
        }

        long double obsBasesProb;
        if (contextR == 0) 
        {
            // Calculate Pr(D)
            obsBasesProb = methOptions.homoMethPriorCG*lHoods[0] + methOptions.homoNonMethPriorCG*lHoods[1] + methOptions.heteroMethPriorCG*lHoods[2];    
           // Calculate posterior probs. for possible genotypes
           // Pr(G|D) = pi(G) * Pr(D|G) / Pr(D)
            meth.homoMethProbR = methOptions.homoMethPriorCG*lHoods[0] / obsBasesProb;
            meth.homoNonMethProbR = methOptions.homoNonMethPriorCG*lHoods[1] / obsBasesProb;
            meth.heteroMethProbR = methOptions.heteroMethPriorCG*lHoods[2] / obsBasesProb;
        }
        else if (contextR == 1)
        {
            obsBasesProb = methOptions.homoMethPriorCHG*lHoods[0] + methOptions.homoNonMethPriorCHG*lHoods[1] + methOptions.heteroMethPriorCHG*lHoods[2];    
            meth.homoMethProbR = methOptions.homoMethPriorCHG*lHoods[0] / obsBasesProb;
            meth.homoNonMethProbR = methOptions.homoNonMethPriorCHG*lHoods[1] / obsBasesProb;
            meth.heteroMethProbR = methOptions.heteroMethPriorCHG*lHoods[2] / obsBasesProb;
        }
        else 
        {
            obsBasesProb = methOptions.homoMethPriorCHH*lHoods[0] + methOptions.homoNonMethPriorCHH*lHoods[1] + methOptions.heteroMethPriorCHH*lHoods[2];    
            meth.homoMethProbR = methOptions.homoMethPriorCHH*lHoods[0] / obsBasesProb;
            meth.homoNonMethProbR = methOptions.homoNonMethPriorCHH*lHoods[1] / obsBasesProb;
            meth.heteroMethProbR = methOptions.heteroMethPriorCHH*lHoods[2] / obsBasesProb;
        }
    }
    else
    {
        if ( meth.snpCalled && ((Dna)(meth.genotype>>2) == 'C' || (Dna)(meth.genotype%4) == 'C') )     // *** CH || HC ***
        {
            unsigned allele2 = ( ((Dna)(meth.genotype>>2) == 'C')? (meth.genotype%4):(meth.genotype>>2) );
            long double singleProbAllele2;
            // for each possible observed read base
            for (unsigned i = 0; i < 4; ++i)
            {
                if(options._debugLevel > 1)
                    std::cout << " base: "<< (Dna)i << std::endl; 
                // for all reads mapped on forward strand
                for (unsigned j = 0; j < length(qualF[i]); ++j)
                {
                    double qual =  static_cast<double>(ordValue(qualF[i][j])-33);
                    e = pow(10.0, (double)(-qual/10.0));
                    // for methylated genotype            
                    if (i == ordValue((Dna)'C'))           // forward strand C                
                         singleProbs[0] = (1.0-e);
                    else                                    // base call error
                        singleProbs[0] = e/3.0;    
                    // for non-methylated genotype
                    if (i == ordValue((Dna)'C'))           // forward strand C                
                        singleProbs[1] = (1.0-e)*(1-methOptions.convRate);
                    else if (i == ordValue((Dna)'T'))      // forward strand T 
                        singleProbs[1] = e/3.0 + (1.0-e)*methOptions.convRate;
                    else                                    // base call error
                        singleProbs[1] = e/3.0;
                    // for non-C genotype
                    if (i == allele2)                       // correct                
                        singleProbAllele2 = (1.0-e);
                    else                                    // base call error
                        singleProbAllele2 = e/3.0;

                    // Calculate likelihood for diploid genotypes
                    lHoods[0] *= (long double)(0.5*singleProbs[0] + 0.5*singleProbAllele2);    // in this case: only one C, methylated
                    lHoods[1] *= (long double)(0.5*singleProbs[1] + 0.5*singleProbAllele2);    // non-methylated          
 
                    if(options._debugLevel > 1)
                    {
                        std::cout << "Forward:"<< std::endl;
                        std::cout << std::setprecision (25) <<  "  my prob: methylated" << " " << (0.5*singleProbs[0] + 0.5*singleProbAllele2)  << std::endl;
                        std::cout << std::setprecision (25) <<  "  my prob: non-methylated" << " " << (0.5*singleProbs[1] + 0.5*singleProbAllele2)  << std::endl;
                    }
                }
                // for all reads mapped on reverse strand
                for (unsigned j = 0; j < length(qualR[i]); ++j)
                {   
                    double qual =  static_cast<double>(ordValue(qualR[i][j])-33);
                    e = pow(10.0, (double)(-qual/10.0));                 
                    // for methylated genotype            
                    if (i == ordValue((Dna)'G'))           // reverse strand C                
                        singleProbs[0] = (1.0-e);
                    else                                    // base call error
                        singleProbs[0] = e/3.0;
                    // for non-methylated genotype
                    if (i == ordValue((Dna)'G'))           // reverse strand C                
                        singleProbs[1] = (1.0-e);
                    else                                    // base call error
                        singleProbs[1] = e/3.0;
                    // for non-C genotype
                    if (i == allele2)                       // correct                
                        singleProbAllele2 = (1.0-e);
                    else                                    // base call error
                        singleProbAllele2 = e/3.0;
    
                    // Calculate likelihood for diploid genotypes
                    lHoods[0] *= (long double)(0.5*singleProbs[0] + 0.5*singleProbAllele2);    // in this case: only one C, methylated
                    lHoods[1] *= (long double)(0.5*singleProbs[1] + 0.5*singleProbAllele2);    // non-methylated
                    
                    if(options._debugLevel > 1)
                    {
                        std::cout << "Reverse:"<< std::endl;
                        std::cout << std::setprecision (25) <<  "  my prob: methylated" << " " << (0.5*singleProbs[0] + 0.5*singleProbAllele2)  << std::endl;
                        std::cout << std::setprecision (25) <<  "  my prob: non-methylated" << " " << (0.5*singleProbs[1] + 0.5*singleProbAllele2)  << std::endl;
                    }
                } 
            }
            
            long double obsBasesProb;
            // TODO: extra case, if ref. base is non-C
            if (contextF == 0) 
            {
                // Calculate Pr(D)
                obsBasesProb = methOptions.cxRefCMethPriorCG*lHoods[0] + (1-methOptions.cxRefCMethPriorCG)*lHoods[1];    
                // Calculate posterior probs. for possible genotypes
                // Pr(G|D) = pi(G) * Pr(D|G) / Pr(D)
                meth.homoMethProbF = methOptions.cxRefCMethPriorCG*lHoods[0] / obsBasesProb;
                meth.homoNonMethProbF = (1-methOptions.cxRefCMethPriorCG)*lHoods[1] / obsBasesProb;
            }
            else if (contextF == 1)
            {
                obsBasesProb = methOptions.cxRefCMethPriorCHG*lHoods[0] + (1-methOptions.cxRefCMethPriorCHG)*lHoods[1];   
                meth.homoMethProbF = methOptions.cxRefCMethPriorCHG*lHoods[0] / obsBasesProb;
                meth.homoNonMethProbF = (1-methOptions.cxRefCMethPriorCHG)*lHoods[1] / obsBasesProb;
            }
            else 
            {
                obsBasesProb = methOptions.cxRefCMethPriorCHH*lHoods[0] + (1-methOptions.cxRefCMethPriorCHH)*lHoods[1];   
                meth.homoMethProbF = methOptions.cxRefCMethPriorCHH*lHoods[0] / obsBasesProb;
                meth.homoNonMethProbF = (1-methOptions.cxRefCMethPriorCHH)*lHoods[1] / obsBasesProb;
            }
        
            if(options._debugLevel > 1)
            {
                std::cout << std::endl;
                std::cout << "obsBasesProb: " << obsBasesProb << std::endl;
                std::cout << "Methylated: lHoods  " << lHoods[0] << "  posterior prob:  " << meth.homoMethProbF << std::endl;
                std::cout << "Non-Methylated: lHoods  " << lHoods[1] << "  posterior prob:  " << meth.homoNonMethProbF << std::endl;
                std::cout << std::endl;
            }
        }
        if ( meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' || (Dna)(meth.genotype%4) == 'G') )     // *** GH || GC ***
        {
            clear(lHoods);
            resize(lHoods, 3, 1.0, Exact());            
            unsigned allele2 = ( ((Dna)(meth.genotype>>2) == 'G')? (meth.genotype%4):(meth.genotype>>2) );
            long double singleProbAllele2;
            // for each possible observed read base
            for (unsigned i = 0; i < 4; ++i)
            {
                // for all reads mapped on forward strand
                for (unsigned j = 0; j < length(qualF[i]); ++j)
                {
                    double qual =  static_cast<double>(ordValue(qualF[i][j])-33);
                    e = pow(10.0, (double)(-qual/10.0));
                    // for methylated genotype            
                    if (i == ordValue((Dna)'G'))           // forward strand G                
                         singleProbs[0] = (1.0-e);
                    else                                    // base call error
                        singleProbs[0] = e/3.0;    
                    // for non-methylated genotype
                    if (i == ordValue((Dna)'G'))           // forward strand G                
                        singleProbs[1] = (1.0-e);
                    else                                    // base call error
                        singleProbs[1] = e/3.0;
                    // for non-C genotype
                    if (i == allele2)                       // correct                
                        singleProbAllele2 = (1.0-e);
                    else                                    // base call error
                        singleProbAllele2 = e/3.0;

                    // Calculate likelihood for diploid genotypes
                    lHoods[0] *= (long double)(0.5*singleProbs[0] + 0.5*singleProbAllele2);    // in this case: only one C, methylated
                    lHoods[1] *= (long double)(0.5*singleProbs[1] + 0.5*singleProbAllele2);    // non-methylated      
                }
                // for all reads mapped on reverse strand
                for (unsigned j = 0; j < length(qualR[i]); ++j)
                {   
                    double qual =  static_cast<double>(ordValue(qualR[i][j])-33);
                    e = pow(10.0, (double)(-qual/10.0));                 
                    // for methylated genotype            
                    if (i == ordValue((Dna)'G'))           // reverse strand C                
                        singleProbs[0] = (1.0-e);
                    else                                    // base call error
                        singleProbs[0] = e/3.0;
                    // for non-methylated genotype
                    if (i == ordValue((Dna)'G'))           // reverse strand C                
                        singleProbs[1] = (1.0-e)*(1-methOptions.convRate);
                    else if (i == ordValue((Dna)'A'))      // reverse strand T                
                        singleProbs[1] = e/3.0 + (1.0-e)*methOptions.convRate;
                    else                                    // base call error
                        singleProbs[1] = e/3.0;
                    // for non-C genotype
                    if (i == allele2)                       // correct                
                        singleProbAllele2 = (1.0-e);
                    else                                    // base call error
                        singleProbAllele2 = e/3.0;
    
                    // Calculate likelihood for diploid genotypes
                    lHoods[0] *= (long double)(0.5*singleProbs[0] + 0.5*singleProbAllele2);    // in this case: only one C, methylated
                    lHoods[1] *= (long double)(0.5*singleProbs[1] + 0.5*singleProbAllele2);    // non-methylated
                } 
            }

            long double obsBasesProb;
            // TODO: extra case, if ref. base is non-C
            if (contextR == 0) 
            {
                // Calculate Pr(D)
                obsBasesProb = methOptions.cxRefCMethPriorCG*lHoods[0] + (1-methOptions.cxRefCMethPriorCG)*lHoods[1];    
               // Calculate posterior probs. for possible genotypes
               // Pr(G|D) = pi(G) * Pr(D|G) / Pr(D)
                meth.homoMethProbR = methOptions.cxRefCMethPriorCG*lHoods[0] / obsBasesProb;
                meth.homoNonMethProbR = (1-methOptions.cxRefCMethPriorCG)*lHoods[1] / obsBasesProb;
            }
            else if (contextR == 1)
            {
                obsBasesProb = methOptions.cxRefCMethPriorCHG*lHoods[0] + (1-methOptions.cxRefCMethPriorCHG)*lHoods[1];   
                meth.homoMethProbR = methOptions.cxRefCMethPriorCHG*lHoods[0] / obsBasesProb;
                meth.homoNonMethProbR = (1-methOptions.cxRefCMethPriorCHG)*lHoods[1] / obsBasesProb;
           }
           else 
            {
                obsBasesProb = methOptions.cxRefCMethPriorCHH*lHoods[0] + (1-methOptions.cxRefCMethPriorCHH)*lHoods[1];   
                meth.homoMethProbR = methOptions.cxRefCMethPriorCHH*lHoods[0] / obsBasesProb;
                meth.homoNonMethProbR = (1-methOptions.cxRefCMethPriorCHH)*lHoods[1] / obsBasesProb;
           }
        }
    }
    return 0;
}


///////////////////////////////////////////////////////////////////////
// snp calling and meth calling separate
template<typename TCounts, typename TQualities, typename TMapqs, typename TMethOptions, typename TOptions, typename TMethylVariant>
inline bool
doBsCalling(TCounts & countF,
          TCounts & countR,
          TQualities & qualF,
          TQualities & qualR,
          TMapqs &mapqsF,
          TMapqs &mapqsR,
          int &refAllele,
          TMethOptions &methOptions,
          TOptions & options,
          TMethylVariant &meth,
          unsigned &contextF,
          unsigned &contextR,
          int pos
          )
{

    // the diploid reference genotype
    int genotypeRef = (refAllele<<2) | refAllele; 

    String<long double> candidateProbs;
    resize(candidateProbs, 4*4); // for simplicity; only 10 are used
    
    getCandidatePostProbs(candidateProbs, methOptions, options, qualF, qualR, countF, countR, refAllele);
   
    // Choose genotype which maximizes the posterior prob.
    int genotype1 = genotypeRef;
    int a1 = 666;
    int a2 = 666;
    int genotype2 = genotypeRef;
    long double maxProb1 = 0.0;
    long double maxProb2 = 0.0;
    for (int g1 = 0; g1 < 4; ++g1)
    {
        for (int g2 = 0; g2 < 4; ++g2)
        {
            if (g1 <= g2)
            {
                if (candidateProbs[(g1<<2)| g2] >= maxProb1) 
                {
                    maxProb2 = maxProb1;
                    genotype2 = genotype1;
                    maxProb1 = candidateProbs[(g1<<2)| g2];
                    genotype1 = (g1<<2)|g2;
                    a1 = g1;
                    a2 = g2;
                }
                else if (candidateProbs[(g1<<2)| g2] >= maxProb2)
                {
                    maxProb2 = candidateProbs[(g1<<2)| g2];
                    genotype2 = (g1<<2)|g2;
                }
            }
        }
    }
    
    unsigned totalCoverage = countF[0] + countF[1] +countF[2] +countF[3] +countF[4]
                           + countR[0] + countR[1] +countR[2] +countR[3] +countR[4];
  
    meth.genotype   = genotype1;
    //std::cout << "a1 " <<  (Dna)a1 << "  a2  " <<  (Dna)a2 << std::endl;
    int consideredCount = countF[a1] + countR[a1] + countF[a2] + countR[a2];                    // TODO: necessary? other checkpoints...
    if ((Dna)a1 == 'C' || (Dna)a2 == 'C') consideredCount += countF[ordValue((Dna)'T')];
    if ((Dna)a1 == 'G' || (Dna)a2 == 'G') consideredCount += countR[ordValue((Dna)'A')];

    if (genotype1 == genotypeRef) meth.snpCalled = false;
    else meth.snpCalled = true;
    if((double)consideredCount/totalCoverage < options.minExplainedColumn)
    {
        meth.snpCalled = false;
    }
    meth.coverageF = countF[0] + countF[1] +countF[2] +countF[3] +countF[4];
    meth.coverageR = countR[0] + countR[1] +countR[2] +countR[3] +countR[4];

    meth.score = log( (long double)candidateProbs[genotype1]/ (long double)candidateProbs[genotype2]);     // genotype calling score
        if (meth.score <= methOptions.minScoreToCallSnp)
        meth.snpCalled = false;

    //if(options._debugLevel > 1)
    if (pos == 6501 || pos == 9945 || pos == 16770 || pos == 20832)
    {
        std::cout << std::setprecision (25) << " prob genotype1.." << (long double)candidateProbs[genotype1] << " prob genotype2.." << (long double)candidateProbs[genotype2]  <<  std::endl;
        std::cout << " genotype1  allele1: "<< (Dna)(genotype1>>2) << "  allele2: " << (Dna)(genotype1 % 4) << std::endl;
        std::cout << " genotype2  allele1: "<< (Dna)(genotype2>>2) << "  allele2: " << (Dna)(genotype2 % 4) << std::endl;

        std::cout << " snp called: "<< (meth.snpCalled? "true":"False") << "  genotypeRef: " << (Dna)(genotypeRef>>2) << " " << (Dna)(genotypeRef % 4) << std::endl;
        std::cout << " genotype1: "<< genotype1 << "  genotypeRef " << genotypeRef << std::endl;
    }
    // Get type of qualities
    typedef typename Value<TQualities>::Type    TQuals;
    typedef typename Value<TQuals>::Type        TQual;

    bool qualOK = false;
    // check avg. quality of C's or C's and T's ...
    if ( (refAllele == ordValue((Dna)'C') && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'C' || (Dna)(meth.genotype%4) == 'C')))
    {
        double qualC = 0.0;
        for(unsigned i = 0; i < length(qualF[ordValue((Dna)'C')]); ++i)
        {
            qualC += static_cast<double>(ordValue(qualF[ordValue((Dna)'C')][i])-33);
            //qualC += getValue(qualF[ordValue((Dna)'C')], i);
        }
        for(unsigned i = 0; i < length(qualR[ordValue((Dna)'C')]); ++i)
        {
            qualC += static_cast<double>(ordValue(qualR[ordValue((Dna)'C')][i])-33);
            //qualC += getValue(qualR[ordValue((Dna)'C')], i);
        }
        double qualT = 0.0;
        for(unsigned i = 0; i < length(qualF[ordValue((Dna)'T')]); ++i)
        {
            qualT += static_cast<double>(ordValue(qualF[ordValue((Dna)'T')][i])-33);
            //qualT += getValue(qualF[ordValue((Dna)'T')], i);
        }
        if (qualC/(countF[ordValue((Dna)'C')]+countR[ordValue((Dna)'C')]) >= options.avgQualT || 
                (qualC + qualT)/(countF[ordValue((Dna)'C')]+countR[ordValue((Dna)'C')]+countF[ordValue((Dna)'T')]) >= options.avgQualT)
            qualOK = true;
    }
    if ( (refAllele == ordValue((Dna)'G') && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' || (Dna)(meth.genotype%4) == 'G'))) // TODO: zuasammenfassen irgendwo
    {
        double qualG = 0.0;
        for(unsigned i = 0; i < length(qualF[ordValue((Dna)'G')]); ++i)
        {
            qualG += static_cast<double>(ordValue(qualF[ordValue((Dna)'G')][i])-33);
            //qualG += getValue(qualF[ordValue((Dna)'G')], i);
        }
        for(unsigned i = 0; i < length(qualR[ordValue((Dna)'G')]); ++i)
        {
            qualG += static_cast<double>(ordValue(qualR[ordValue((Dna)'G')][i])-33);
            //qualG += getValue(qualR[ordValue((Dna)'G')], i);
        }
        double qualA = 0.0;
        for(unsigned i = 0; i < length(qualR[ordValue((Dna)'A')]); ++i)
        {
            qualA += static_cast<double>(ordValue(qualR[ordValue((Dna)'A')][i])-33);
            //qualA += getValue(qualR[ordValue((Dna)'A')], i);
        }
        if (qualG/(countF[ordValue((Dna)'G')]+countR[ordValue((Dna)'G')])  >= options.avgQualT || 
                (qualG + qualA)/(countF[ordValue((Dna)'G')]+countR[ordValue((Dna)'G')]+countR[ordValue((Dna)'A')])  >= options.avgQualT)
            qualOK = true;
    }

    // output if outputCandidates ?
    if (methOptions.outputMethStates &&  qualOK)
        getMethState(meth, countF, countR, refAllele, methOptions);

    if (methOptions.outputMethProbs && qualOK)
        getMethProb(meth, qualF, qualR, countF, countR, refAllele, contextF, contextR, methOptions, options);

    return true;
}


// write to file
// Merge with other version?
template<typename TFile, typename TMethylVariant, typename TQualities, typename TString, typename TPos, typename TMethOptions, typename TOptions>
inline bool
writeMeth(TFile &file, 
       TMethylVariant &meth,
       TQualities &qualityStringF, 
       TQualities &qualityStringR,
       int refAllele, 
       TString &genomeID, 
       TPos candPos, 
       unsigned realCoverage,
       TMethOptions &methOptions,
       TOptions &options)
{
//IOREV _nodoc_ what kind of format is this?
    if (!file.is_open()) 
    {
        ::std::cerr << "SNP/Meth output file is not open" << std::endl;
        return false;
    }

    //chromosome
    file << genomeID << '\t';
    //file << candPos + options.positionFormat<< '\t';
    file << candPos << '\t';        // 0-based:w

    file << (Dna)refAllele <<'\t';

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
    

    //genotypeCalled to string
    if(meth.snpCalled)//genotypeCalled != genotypeRef) 

        file << '\t' << (Dna)(meth.genotype>>2) << (Dna)(meth.genotype % 4) << '\t'<< meth.score; // << '\t' << snp.quality << '\t' << snp.snpQuality;
    else
        file << "\t.\t."; 
  

    // Methylation Calling
    
    if (methOptions.outputMethStates)
    {
        if ( (refAllele == ordValue((Dna)'C') && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'C' && (Dna)(meth.genotype%4) == 'C')) )  // CC
        {
            file << '\t' << meth.methLevelF;
            if (meth.called1 == true)
                file << '\t' << (meth.state1 ? 'M':'X') << (meth.state2 ? 'M':'X');
            else
                file << '\t' << "..";
        }
        else if ( (refAllele == ordValue((Dna)'G') && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' && (Dna)(meth.genotype%4) == 'G')) )  // GG
        {
            file << '\t' << meth.methLevelR;
            if (meth.called2 == true)
                file << '\t' << (meth.state1 ? 'M':'X') << (meth.state2 ? 'M':'X');
            else
                file << '\t' << "..";
        }
        else if (refAllele == ordValue((Dna)'C') || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'C' || (Dna)(meth.genotype%4) == 'C')) )     // CH || HC
        {
            file << '\t' << meth.methLevelF;

            // and G...
            if (refAllele == ordValue((Dna)'G') || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' || (Dna)(meth.genotype%4) == 'G')) )     // GC
                file << ':' << meth.methLevelR;
         
            if (meth.called1 == true)
                file << '\t' << (meth.state1 ? 'M':'X') ;
            else
                file << '\t' << ".";
            
            // and G...
            if (refAllele == ordValue((Dna)'G') || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' || (Dna)(meth.genotype%4) == 'G')) )     // GC
            {
                if (meth.called2 == true)
                    file << ':' << (meth.state2 ? 'M':'X');
                else
                    file << ':' << ".";
            }


        }
        else if (refAllele == ordValue((Dna)'G') || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' || (Dna)(meth.genotype%4) == 'G')) )     // GH
        {
            file << '\t' << meth.methLevelR;
            if (meth.called2 == true)
                file << '\t' << (meth.state2 ? 'M':'X');
            else
                file << '\t' << ".";

        }
    }

    if (methOptions.outputMethProbs)
    {
        file << '\t';
        if ( (refAllele == ordValue((Dna)'C') && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'C' && (Dna)(meth.genotype%4) == 'C')) )  // CC
        {
            file << meth.homoMethProbF << ',' << meth.homoNonMethProbF << ',' << meth.heteroMethProbF ;
        }
        else if (refAllele == ordValue((Dna)'C') || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'C' || (Dna)(meth.genotype%4) == 'C')) )     // CH || HC
        {
            file << meth.homoMethProbF << ',' << meth.homoNonMethProbF ;
            file << ":";
        }

        if ( (refAllele == ordValue((Dna)'G') && !meth.snpCalled) || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' && (Dna)(meth.genotype%4) == 'G')) )  // GG
        {
            file << meth.homoMethProbR << ',' << meth.homoNonMethProbR << ',' << meth.heteroMethProbR ;
        }
        else if (refAllele == ordValue((Dna)'G') || (meth.snpCalled && ((Dna)(meth.genotype>>2) == 'G' || (Dna)(meth.genotype%4) == 'G')) )     // GH || GC
        {
            file << meth.homoMethProbR << ',' << meth.homoNonMethProbR ;
        }
    }
    file << std::endl;
    return true;
}


#endif  // #ifndef SANDBOX_KRAKAU_APPS_SNP_METH_STORE_BS_SEP_CALLING_H_
