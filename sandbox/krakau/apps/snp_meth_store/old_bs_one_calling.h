/*==========================================================================

Methylation Calling: Version 2

Calculate snp and methylation probability in one step
    Pr(Dj | H1 = T)
    Pr(Dj | H1 = C_M)
    Pr(Dj | H1 = C)
    ...

Add additional cases: e.g. prob. converted and sequencing error
Add underconversion rates...

==========================================================================*/

#ifndef SANDBOX_KRAKAU_APPS_SNP_METH_STORE_BS_ONE_CALLING_H_
#define SANDBOX_KRAKAU_APPS_SNP_METH_STORE_BS_ONE_CALLING_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

#include "bs_alphabets.h"
//#include "snp_meth_store.h"
//#include "meths.h"

using namespace std;
using namespace seqan;

// TODO:
// precalculate probs and get rid of thausend different cases to check
// Compute prob. to observe given base i under the assumption that the underlying haplotype h
// For reads mapped to the forward strand
template<typename TProb, typename TErrorProb, typename TMethOptions>
inline void
getSingleBaseProbHaploF(TProb &singleProb, Dna i, DnaM h, TErrorProb &e, TMethOptions &methOptions)
{
    if ( h == 'C')                              // Haplotype C
    {                     
        if (i == 'C')           
            singleProb = (1.0-e)*(1-methOptions.convRate);
        else if (i == 'T')     
            singleProb = e/3.0 + (1.0-e)*(methOptions.convRate);          
        else                                   
            singleProb = e/3.0;
    } 
    else if ( h == 'D')                         // Haplotype Cm
    {                     
        if (i == 'C')          
            singleProb = (1.0-e)*(1.0-methOptions.methConvRate);
        else if (i == 'T')   
            singleProb = e/3.0 + (1.0-e)*methOptions.methConvRate;
        else                                   
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
// For reads mapped to the reverse strand
template<typename TProb, typename TErrorProb,  typename TMethOptions>
inline void
getSingleBaseProbHaploR(TProb &singleProb, Dna i, DnaM h, TErrorProb &e, TMethOptions &methOptions)
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
            singleProb = (1.0-e)*(1-methOptions.convRate);
        else if (i == 'A')
            singleProb = e/3.0 + (1.0-e)*methOptions.convRate;
        else                                  
            singleProb = e/3.0;
    } 
    else if ( h == 'H')                         // Haplotype Gm
    {
        if (i == 'G')         
            singleProb = (1.0-e)*(1-methOptions.methConvRate);
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


template<typename TProbs, typename TMethOptions, typename TOptions, typename TQStrings, typename TCounts>
void 
getCandidatePostProbs(TProbs &candidatePostProbs, TMethOptions &methOptions, TOptions &options, TQStrings &qualF, TQStrings &qualR, TCounts &countF, TCounts &countR, int &refAllele,
                     unsigned &contextF, unsigned &contextR )
{
    if (methOptions.ignoreBs) // Snp calling without bs conversions
        methOptions.convRate = 0.0;
 
    // likelihood to observe single base under assumption of given genotype
    String<long double> singleProbs;
    resize(singleProbs, 6);

    String<long double> lHoods;  // likelihoods to observe observed data under assumption of given genotypes 
    resize(lHoods, 8*8, 1.0, Exact()); //??

    //std::cout << "Test 0 " << std::endl;
    // for each observed base type
    for (unsigned i = 0; i < 4; ++i)        
    {
        //std::cout << "Test 0.0 " << std::endl;
        // for all reads mapped on forward strand
        for (unsigned j = 0; j < length(qualF[i]); ++j)
        {
            double qual =  static_cast<double>(ordValue(qualF[i][j])-33);
            long double e = pow(10.0, (double)(-qual/10.0));
            // for each possible candidate genotype
            for (unsigned h = 0; h < 6; ++h)
                getSingleBaseProbHaploF(singleProbs[h], (Dna)i, (DnaM)h, e, methOptions);
            
            // Calculate likelihood for diploid genotypes
            for (unsigned h1 = 0; h1 < 6; ++h1)
            {
                for (unsigned h2 = h1; h2 < 6; ++h2)
                {
                    //std::cout << "Test genotype in int" << ((h1<<3)|h2) << std::endl;
                    lHoods[(h1<<3)|h2] *= (long double)(0.5*singleProbs[h1] + 0.5*singleProbs[h2]);
                    if(options._debugLevel > 1)
                    {    
                        if ((Dna)h1 == 'G' && (Dna)h2 == 'G')
                        {
                            std::cout << "Forward:"<< std::endl;
                            std::cout << std::setprecision (25) <<  "  my prob: " << (Dna)h1 << (Dna)h2 << " " << (0.5*singleProbs[h1] + 0.5*singleProbs[h2])  << std::endl;
                        }
                    }
                }
            }
        }
        //std::cout << "Test 0.1 " << std::endl;

        // for all reads mapped on reverse strand
        for (unsigned j = 0; j < length(qualR[i]); ++j)
        {
            double qual =  static_cast<double>(ordValue(qualR[i][j])-33);
            long double e = pow(10.0, (double)(-qual/10.0));  
            // for each possible candidate genotype
            for (unsigned h = 0; h < 6; ++h)
                getSingleBaseProbHaploR(singleProbs[h], (Dna)i, (DnaM)h, e, methOptions);
            
            // Calculate likelihood for diploid genotypes
            for (int h1 = 0; h1 < 6; ++h1)
            {
                for (int h2 = h1; h2 < 6; ++h2)
                {
                    lHoods[(h1<<3)|h2] *=(long double)(0.5*singleProbs[h1] + 0.5*singleProbs[h2]);
                    if(options._debugLevel > 1)
                    {
                        if ((Dna)h1 == 'G' && (Dna)h2 == 'G')
                        {
                            std::cout << "Reverse:" << std::endl;
                            std::cout << std::setprecision (25) << "  my prob: " << (Dna)h1 << (Dna)h2 << " " << (long double)(0.5*singleProbs[h1] + 0.5*singleProbs[h2]) << std::endl;
                        }
                    }                    
                }
            }
        }
    }

    // Calculate Pr(D)
    long double obsBasesProb = 0.0;
    for (int h1 = 0; h1 < 6; ++h1)
    {
        for (int h2 = h1; h2 < 6; ++h2)
        {
            obsBasesProb += methOptions.genPriors[refAllele<<4| ordValue((Dna)((DnaM)h1))<<2| ordValue((Dna)((DnaM)h2))] * lHoods[(h1<<3)| h2];
        }
    }

    // Calculate posterior probs. for possible genotypes
    // Pr(G|D) = pi(G) * Pr(D|G) / Pr(D9)
    // For methylated/nonmethylated genotypes take prior methylation/non-methylation probability into account, dependent on context
    for (int h1 = 0; h1 < 6; ++h1)
    {
        for (int h2 = h1; h2 < 6; ++h2)
        {
            if(options._debugLevel > 1)
                std::cout << (Dna)refAllele << (Dna)h1 << (Dna)h2 << "  " << std::setprecision (25) << "genPrior: " << std::setprecision (25) << (long double)methOptions.genPriors[ refAllele<<4| ordValue((Dna)((DnaM)h1))<<2 | ordValue((Dna)((DnaM)h2))] << "  lHood  " << (long double)lHoods[(h1<<3)| h2] << std::endl;
                  
            // for bsPriors: maybe later: take into account, if reAllele is C or not
            if ((Dna)((DnaM)h1) == 'C' && (Dna)((DnaM)h2) == 'G')
                candidatePostProbs[(h1<<3)| h2] = methOptions.bsPriors[h2][h1][contextF*3+contextR] * methOptions.genPriors[ refAllele<<4| ordValue((Dna)((DnaM)h1))<<2| ordValue((Dna)((DnaM)h2))] * lHoods[(h1<<3)| h2] / obsBasesProb;   
            else if ((Dna)((DnaM)h1) == 'C' || (Dna)((DnaM)h2) == 'C' )
                candidatePostProbs[(h1<<3)| h2] = methOptions.bsPriors[h2][h1][contextF] * methOptions.genPriors[ refAllele<<4| ordValue((Dna)((DnaM)h1))<<2| ordValue((Dna)((DnaM)h2))] * lHoods[(h1<<3)| h2] / obsBasesProb; 
            else if ((Dna)((DnaM)h1) == 'G' || (Dna)((DnaM)h2) == 'G' )
                candidatePostProbs[(h1<<3)| h2] = methOptions.bsPriors[h2][h1][contextR] * methOptions.genPriors[ refAllele<<4| ordValue((Dna)((DnaM)h1))<<2| ordValue((Dna)((DnaM)h2))] * lHoods[(h1<<3)| h2] / obsBasesProb; 
            else
                candidatePostProbs[(h1<<3)| h2] = methOptions.genPriors[ refAllele<<4| ordValue((Dna)((DnaM)h1))<<2| ordValue((Dna)((DnaM)h2))] * lHoods[(h1<<3)| h2] / obsBasesProb;     // No methylation possible, only A or T

            if(options._debugLevel > 1)
                 std::cout << std::setprecision (25) << "candidateProb: " << (Dna)refAllele << (Dna)h1 << (Dna)h2 << "  " << candidatePostProbs[(h1<<3)| h2] << std::endl;
        }
    }

}

///////////////////////////////////////////////////////////////////////
// SNP and meth calling in one step
template<typename TCounts, typename TQualities, typename TMethOptions, typename TOptions>
inline bool
doBsCalling(TCounts & countF,
          TCounts & countR,
          TQualities & qualF,
          TQualities & qualR,
          int &refAllele,
          TMethOptions &methOptions,
          TOptions & options,
          MethylVariant<OneCallMethod> &meth,
          unsigned &contextF,
          unsigned &contextR,
          int pos
          )
{
    int genotypeRef = (refAllele<<3) | refAllele;
    int qCall1 = 0;  // genotype call quality
    int qSnp = 0; // SNP call quality

    String<long double> candidateProbs;
    resize(candidateProbs, 8*8); // for simplicity; not all are used
    getCandidatePostProbs(candidateProbs, methOptions, options, qualF, qualR, countF, countR, refAllele, contextF, contextR);

    // Choose genotype which maximizes the posterior prob.
    int genotype1 = genotypeRef;
    int bsAllele1 = 666;
    int bsAllele2 = 666;
    int genotype2 = genotypeRef;
    long double maxProb1 = 0.0;
    long double maxProb2 = 0.0;
    for (int h1 = 0; h1 < 6; ++h1)
    {
        for (int h2 = h1; h2 < 6; ++h2)
        {
            if (candidateProbs[(h1<<3)| h2] >= maxProb1) 
            {
                maxProb2 = maxProb1;
                genotype2 = (bsAllele1<<3)|bsAllele2;
                maxProb1 = candidateProbs[(h1<<3)| h2];
                bsAllele1 = h1;
                bsAllele2 = h2;
            }
            else if (candidateProbs[(h1<<3)| h2] >= maxProb2)
            {
                maxProb2 = candidateProbs[(h1<<3)| h2];
                genotype2 = (h1<<3)|h2;
            }
        }
    }
    genotype1 = (bsAllele1<<3)|bsAllele2;
    int allele1 = ordValue((Dna)((DnaM)bsAllele1));
    int allele2 = ordValue((Dna)((DnaM)bsAllele2));

    unsigned totalCoverage = countF[0] + countF[1] +countF[2] +countF[3] +countF[4]
                           + countR[0] + countR[1] +countR[2] +countR[3] +countR[4];
  
    meth.genotype = genotype1;
    meth.score = log( (long double)candidateProbs[genotype1]/ (long double)candidateProbs[genotype2]);     // genotype calling score (only undrlying genotypes, no bs types)?

    if ( ((allele1<<3)|allele2) == genotypeRef) 
        meth.snpCalled = false;      // Convert DnaM ordValues into Dna ordValues and check if genotype called is the as the reference
    else if (meth.score <= methOptions.minScoreToCallSnp)
        meth.snpCalled = false;
    else 
        meth.snpCalled = true;

    meth.snpProb = candidateProbs[genotype1];

    // TODO: get clear cases, if snp not called, score to low, meth states called, 
    // what if snp score to low, do we want to call bs states?
    // maybe output different probs (could be the reason for to low score, in snps with C or T) ???

    // If bs case:
    if ((Dna)refAllele == 'C' || (Dna)refAllele == 'G' || 
       (meth.snpCalled && ((Dna)allele1 == 'C' || (Dna)allele1 == 'G' || (Dna)allele2 == 'C' || (Dna)allele2 == 'G')) )
    {
        meth.bsCalled = true;
    }


    //if(options._debugLevel > 1)
    if (pos == 33  || pos == 389 || pos == 5113 )
    {
        std::cout << "Pos: " << pos << std::endl;
        std::cout << std::setprecision (25) << " prob genotype1.." << (long double)candidateProbs[genotype1] << " prob genotype2.." << (long double)candidateProbs[genotype2]  <<  std::endl;
        std::cout << " genotype1  allele1: "<< (DnaM)(genotype1>>3) << "  allele2: " << (DnaM)(genotype1 % 8) << std::endl;
        std::cout << " genotype2  allele1: "<< (DnaM)(genotype2>>3) << "  allele2: " << (DnaM)(genotype2 % 8) << std::endl;
        std::cout << std::setprecision (25) << " prob genotype AG.." << (long double)candidateProbs[ordValue((DnaM)'A')>>3|ordValue((DnaM)'G')]  <<  std::endl;
        std::cout << std::setprecision (25) << " prob genotype3.." << (long double)candidateProbs[0>>3|2]  <<  std::endl;

    }

    // calculate average quality of genotype quality ?
    //std::cout << " genotype1  allele1: "<< (DnaM)(genotype1>>3) << "  allele2: " << (DnaM)(genotype1 % 8) << std::endl;

    if (methOptions.outputAllBsStateProbs && meth.bsCalled)
    {
        clear(meth.bsStateProbs);
        // Store all possible meth. states for underlying genotype
        // Following order holds for diff. genotypes:
            // CC, CD, DD
            // GG, GH, HH
            // CG, DG, CH, DH
            // CX, DX
            // GX, HX
        if ( (Dna)allele1 == 'C' && (Dna)allele2 == 'C')
        {
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'D')>>3|ordValue((DnaM)'D')]);
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'C')>>3|ordValue((DnaM)'D')]);
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'D')>>3|ordValue((DnaM)'D')]);
        }
        else if ((Dna)allele1 == 'G' && (Dna)allele2 == 'G')
        {
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'G')>>3|ordValue((DnaM)'G')]);
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'G')>>3|ordValue((DnaM)'H')]);
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'H')>>3|ordValue((DnaM)'H')]);
        }
        else if ((Dna)allele1 == 'C' && (Dna)allele2 == 'G')
        {
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'C')>>3|ordValue((DnaM)'G')]);
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'D')>>3|ordValue((DnaM)'G')]);
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'C')>>3|ordValue((DnaM)'H')]);
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'D')>>3|ordValue((DnaM)'H')]);
        }
        else if ((Dna)allele1 == 'C')
        {
            //std::cout << " 1 " << ordValue((DnaM)'C')>>3<< std::endl;
            //std::cout << "bsAllele2 " << bsAllele2 << std::endl;
            //std::cout << "prob " << candidateProbs[ordValue((DnaM)'C')>>3|bsAllele2] << std::endl;
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'C')>>3|bsAllele2]);
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'D')>>3|bsAllele2]);
        }
        else if ((Dna)allele2 == 'C')
        {
            appendValue(meth.bsStateProbs, candidateProbs[bsAllele1|ordValue((DnaM)'C')]);
            appendValue(meth.bsStateProbs, candidateProbs[bsAllele1|ordValue((DnaM)'D')]);
        }
        else if ((Dna)allele1 == 'G')
        {
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'G')>>3|bsAllele2]);
            appendValue(meth.bsStateProbs, candidateProbs[ordValue((DnaM)'H')>>3|bsAllele2]);
        }
        else if ((Dna)allele2 == 'G')
        {
            appendValue(meth.bsStateProbs, candidateProbs[bsAllele1|ordValue((DnaM)'G')]);
            appendValue(meth.bsStateProbs, candidateProbs[bsAllele1|ordValue((DnaM)'H')]);
        }
    }
 

  
    return true;
}



// write to file
// Merge with other version?
template<typename TFile, typename TQualities, typename TString, typename TPos, typename TMethOptions, typename TOptions>
inline bool
writeMeth(TFile &file,  
       MethylVariant<OneCallMethod> &meth,
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
        file << '\t' << (Dna)((DnaM)(meth.genotype>>3)) << (Dna)((DnaM)(meth.genotype%8)) << '\t' << meth.snpProb << '\t' << meth.score; 
    else
        file << "\t.\t."; 
 

    // Methylation Calling
    // testing
    file << '\t' << (DnaM)(meth.genotype>>3) <<  (DnaM)(meth.genotype%8) << '\t';

    // output called meth state...
    if ( (DnaM)(meth.genotype>>3) == 'C' || (DnaM)(meth.genotype>>3) == 'G' )
        file << 'X'; 
    if ( (DnaM)(meth.genotype>>3) == 'D' || (DnaM)(meth.genotype>>3) == 'H' )
        file << 'M';
    if ( (DnaM)(meth.genotype%8) == 'C' || (DnaM)(meth.genotype%8) == 'G' )
        file << 'X'; 
    if ( (DnaM)(meth.genotype%8) == 'D' || (DnaM)(meth.genotype%8) == 'H' )
        file << 'M';


    if (methOptions.outputAllBsStateProbs && meth.bsCalled)
    {
        file << '\t';
        // write probs of underlying genotype for other possible meth states 
        for (unsigned i = 0; i < length(meth.bsStateProbs); ++i)
        {
            //std::cout << "test 4 " << meth.bsStateProbs[i] << std::endl;
            file << meth.bsStateProbs[i] << ',';
        }
    }
    //std::cout << "test 5 " << std::endl;
    file << std::endl;
    return true;
}




#endif  // #ifndef SANDBOX_KRAKAU_APPS_SNP_METH_STORE_BS_ONE_CALLING_H_
