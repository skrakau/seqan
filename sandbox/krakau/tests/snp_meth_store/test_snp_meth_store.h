// ==========================================================================
//                               snp_meth_store
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef SANDBOX_KRAKAU_TESTS_SNP_METH_STORE_TEST_SNP_METH_STORE_H_
#define SANDBOX_KRAKAU_TESTS_SNP_METH_STORE_TEST_SNP_METH_STORE_H_

#include <seqan/basic.h>

#include <seqan/platform.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <math.h>
#include <cmath>

#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/assert.hpp>
#include <boost/math/tools/rational.hpp>
#include <boost/math/tools/real_cast.hpp>
#include <boost/math/special_functions/binomial.hpp>


#include "../../apps/snp_meth_store/snp_meth_store.h"
#include "../../apps/snp_meth_store/meths.h"
#include "../../apps/snp_meth_store/bs_alphabets.h"

#include "../../apps/snp_meth_store/bs_one_calling.h"

using namespace seqan;

template<typename TQStrings, typename TCounts>
void 
setUp_reads1(TQStrings &qualF, TQStrings &qualR, TCounts &countF, TCounts &countR)
{
    resize(qualF, 4, Exact());
    resize(qualR, 4, Exact());
   /* 
    qualF[0] = "";
    qualF[1] = "";
    qualF[2] = "JEU>CMKJHFLJKKKESFPJMEKEMIDIJH?LMHIDIFIOLLKGFGKKIKIIJGEFFFKHHKGHIHIHGJIIHIHHI";
    qualF[3] = "";
    qualR[0] = "IHHFFEGIIGIMEIGLHIKEKIEMBGJBMELGIGEJ>PBCDIOLLKGFGKKIKIIJGEFFFKHHKGHIHIHGJIIHIHHI";
    qualR[1] = "";
    qualR[2] = "IIHHIHH";
    qualR[3] = "";
    */
    qualF[0] = "";
    qualF[1] = "";
    qualF[2] = "ICH@EOHNGIGFHHHHIII"; //CPKCRFMHJKJHHJGI";
    qualF[3] = "";
    qualR[0] = "IGJDNK"; //MDMKHBD";
    qualR[1] = "";
    qualR[2] = "IJGBEHA"; //IHIFGIEGMI";
    qualR[3] = "";
  
    /*
    qualF[0] = "";
    qualF[1] = "";
    qualF[2] = "PFNHR";
    qualF[3] = "";
    qualR[0] = "HHHIGGGGHHHHKFFJEHDKJKBEHKIMJGGGGHHHHK";
    qualR[1] = "";
    qualR[2] = "HHHIGGG";
    qualR[3] = "";
*/

    resize(countF, 4, Exact());
    resize(countR, 4, Exact());
    for (int i = 0; i < 4; ++i)
    {
        countF[i] = length(qualF[i]);
        countR[i] = length(qualR[i]);
    }
}

// Test polynomial evaluation
SEQAN_DEFINE_TEST(test_snp_meth_store_polynom_evaluation)
{
    SNPCallingOptions<>     options;
    MethCallingOptions      methOptions;
    String<CharString>      qualF;
    String<CharString>      qualR;
    String<int>             countF;
    String<int>             countR;
    setUp_reads1(qualF, qualR, countF, countR);

    int pos = 123;

    std::cout << " NSpace, Newton: " << std::endl;
    String<long double> betas1; 
    String<long double> lHoods1; 
    String<String<long double> > constantSet1;    
    setUpConstants(constantSet1, lHoods1, countF, countR, Polynomial());
    constructConstantsAndLHoods(constantSet1, lHoods1, qualF, qualR, countF, countR, methOptions, Polynomial());
    getBetasAndLHoods(betas1, lHoods1, constantSet1, countF, countR, methOptions, Newton(), Horner(), pos);
    std::cout << " End  NSpace, Newton: " << std::endl;

    std::cout << " NSpace, Naive sampling: " << std::endl;
    String<long double> betas2;
    String<long double> lHoods2; 
    String<String<String<long double> > > constantSet2;   
    setUpConstants(constantSet2, lHoods2, countF, countR, Naive());
    constructConstantsAndLHoods(constantSet2, lHoods2, qualF, qualR, countF, countR, methOptions, Naive());
    getBetasAndLHoods(betas2, lHoods2, constantSet2, countF, countR, methOptions, Sampling(), Naive(), pos);
   
    std::cout << " LogSpace, Newton: " << std::endl;
    String<long double> betas3;
    String<LogProbValue> lHoods3; 
    String<String<LogProbValue> > constantSet3;   
    setUpConstants(constantSet3, lHoods3, countF, countR, Polynomial());
    constructConstantsAndLHoods(constantSet3, lHoods3, qualF, qualR, countF, countR, methOptions, Polynomial());
    getBetasAndLHoods(betas3, lHoods3, constantSet3, countF, countR, methOptions, Newton(), Horner(), pos);

    std::cout << " LogSpace, Naive Sampling: " << std::endl;
    String<long double> betas4;
    String<LogProbValue> lHoods4; 
    String<String<String<LogProbValue> > > constantSet4;   
    setUpConstants(constantSet4, lHoods4, countF, countR, Naive());
    constructConstantsAndLHoods(constantSet4, lHoods4, qualF, qualR, countF, countR, methOptions, Naive());
    getBetasAndLHoods(betas4, lHoods4, constantSet4, countF, countR, methOptions, Sampling(), Naive(), pos);


    /////////////////////////////////////////////
    // Compare
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            std::cout << (Dna)h1 << (Dna)h2 << std::endl;
            std::cout << "  beta1: " << betas1[(h1<<2)|h2] << "   lHood1: " << lHoods1[(h1<<2)|h2] << std::endl;
            std::cout << "  beta2: " << betas2[(h1<<2)|h2] << "   lHood2: " << lHoods2[(h1<<2)|h2]  <<  std::endl;
            std::cout << "  beta3: " << betas3[(h1<<2)|h2] << "   lHood3: " <<  lHoods3[(h1<<2)|h2].c  << "  " << lHoods3[(h1<<2)|h2].zero << std::endl;
            std::cout << "  beta4: " << betas4[(h1<<2)|h2] << "   lHood4: " << lHoods4[(h1<<2)|h2].c  << "  " << lHoods4[(h1<<2)|h2].zero <<  std::endl;


            if ((Dna)h1 == 'C' && (Dna)h2 == 'G')
            {
                std::cout << "  beta1: " << betas1[(h2<<2)|h1] << "   lHood1: " << lHoods1[(h1<<2)|h2] << std::endl;
                std::cout << "  beta2: " << betas2[(h2<<2)|h1] << "   lHood2: " << lHoods2[(h1<<2)|h2]  <<  std::endl;
                std::cout << "  beta3: " << betas3[(h2<<2)|h1] << "   lHood3: " << lHoods3[(h1<<2)|h2].c  <<  std::endl;
            }
            //SEQAN_ASSERT_LT(abs(lHoods1[(h1<<2)|h2] - lHoods2[(h1<<2)|h2]), 0.001);
        }
    }

    /*
    // test : scale coeffs directly 
    for (int i = length(constantSet3[2<<2|2])-1; i >= 0; --i)
    {
        constantSet3[2<<2|2][i] *= pow(10, -14);
    }
*/

    /*
    // GG
    for (long double beta = 0.0; beta <= 1.00001; beta += 0.01)
    {
            std::cout << "Beta: " << beta << '\t' << "NS Newton: " << getProbForBeta(beta, constantSet1[2<<2|2], Horner()) << std::endl;
    }
    
    std::cout << "coeffs:   normal    logspace " << std::endl;
    for (int i = 0; i < length(constantSet1[2<<2|3])-1; ++i)
    {
        std::cout << '\t' << constantSet1[2<<2|3][i] << "\t\t"   << "  " << constantSet3[2<<2|3][i].c << " " << constantSet3[2<<2|3][i].sign << " " << constantSet3[2<<2|3][i].zero  << std::endl; 
    }

    std::cout << "coeffs:   " << std::endl;
    for (int i = 0; i < length(constantSet1[2<<2|3]); ++i)
    {
        std::cout << '\t' << constantSet1[2<<2|3][i] << "\t\t"   << "  " << constantSet3[2<<2|3][i].c << " " << constantSet3[2<<2|3][i].sign << " " << constantSet3[2<<2|3][i].zero  << std::endl; 
    }
*/

    /*
    String<LogProbValue> co;
    resize(co, 4);
    LogProbValue c0;
    c0.zero = false;
    c0.sign = true;
    c0.c = 4;
    co[0] = c0;
    LogProbValue c1;
    c1.zero = false;
    c1.sign = true;
    c1.c = 2;
    co[1] = c1;
    LogProbValue c2;
    c2.zero = false;
    c2.sign = true;
    c2.c = 7;
    co[2] = c2;
    LogProbValue c3;
    c3.zero = false;
    c3.sign = true;
    c3.c = 5;
    co[3] = c3;
    // GT
    std::cout << "LogSpace: " << std::endl;
    P_functorLeB<long double> myFL(co);
    for (long double beta = 0.0; beta <= 1.00001; beta += 0.1)
    {
             boost::math::tuple<long double, long double, long double> tupleL = myFL(beta);      
    }
    // GT
    String<long double> coN;
    resize(coN, 4);
    coN[0] = 4;
    coN[1] = 2;
    coN[2] = 7;
    coN[3] = 5;
    std::cout << "NSpace: " << std::endl;
    P_functor<long double> myFN(coN);
    for (long double beta = 0.0; beta <= 1.0001; beta += 0.1)
    {
             boost::math::tuple<long double, long double> tupleN = myFN(beta);      
    }
    */
    // Check coeffs construction:
/*
    String<long double> testSet;
    resize(testSet, 5);
    
    long double a;
    long double b;
    unsigned r;
    r = 0;
    a = 2;
    b = 3;
    addFactorToConstants(testSet, a, b, r, Polynomial(), NSpace()); // TODO: interface changed
    r = 1;
    a = 4;
    b = 5;
    addFactorToConstants(testSet, a, b, r, Polynomial(), NSpace());
    r = 2;
    a = 6;
    b = 7;
    addFactorToConstants(testSet, a, b, r, Polynomial(), NSpace());
    r = 3;
    a = 8;
    b = 9;
    addFactorToConstants(testSet, a, b, r, Polynomial(), NSpace());

    std::cout << "test coeffs: " << std::endl;
    for (int i = 0; i < length(testSet); ++i)
    {
        std::cout << "  " << testSet[i] << std::endl; 
    }
    // should be 384, 1936, 3644, 3036, 945
    */

    long double bla1 = 260723.786688;
    long double bla2 = 392.3678835;

    long double res = bla1 - bla2;

    LogProb<long double> log1 = LogProb<long double>(bla1);
    LogProb<long double> log2 = LogProb<long double>(bla2);

    LogProb<long double> logRes = log1 - log2;

    std::cout << " bla1: " << bla1 << "  bla2: " << bla2 << "  res: " << res << std::endl;

    std::cout << " log1: " << log1 << "  log: " << log2 << "  logRes: " << logRes << "  " <<  static_cast<long double>(logRes)  <<  std::endl;

    std::cout << " Precision: " << std::numeric_limits<long double>::digits10 << std::endl;
    std::cout << " std::numeric_limits<T>::digits " << std::numeric_limits<long double>::digits << std::endl;
}


// Test polynomial evaluation
SEQAN_DEFINE_TEST(test_snp_meth_store_sampling_vs_newton)
{
    SNPCallingOptions<>     options;
    MethCallingOptions      methOptions;
    String<CharString>      qualF;
    String<CharString>      qualR;
    String<int>             countF;
    String<int>             countR;

    setUp_reads1(qualF, qualR, countF, countR);
/*
    // Sampling
    String<String<long double> > coeffSet;    
    String<long double> lHoodsS;
    constructCoeffsAndLHoods(coeffSet, lHoodsS, qualF, qualR, countF, countR, methOptions);
    String<long double> betasS;
    int pos = 123;
    getBetasAndLHoods(betasS, lHoodsS, coeffSet, countF, countR, Sampling(), pos);

    
    // Newton
    clear(coeffSet);
    String<long double> lHoodsN;
    constructCoeffsAndLHoods(coeffSet, lHoodsN, qualF, qualR, countF, countR, methOptions);
    String<long double> betasN;
    getBetasAndLHoods(betasN, lHoodsN, coeffSet, countF, countR, Newton(), pos);

    
    /////////////////////////////////////////////
    // Compare
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            std::cout << (Dna)h1 << (Dna)h2 << std::endl;
            std::cout << "  betaS: " << betasS[(h1<<2)|h2] << "   lHoodS: " << lHoodsS[(h1<<2)|h2] << std::endl;
            std::cout << "  betaN: " << betasN[(h1<<2)|h2] << "   lHoodN: " << lHoodsN[(h1<<2)|h2]  <<  std::endl;
            if ((Dna)h1 == 'C' && (Dna)h2 == 'G')
            {
                std::cout << "  betaS: " << betasS[(h2<<2)|h1] << "   lHoodS: " << lHoodsS[(h1<<2)|h2] << std::endl;
                std::cout << "  betaN: " << betasN[(h2<<2)|h1] << "   lHoodN: " << lHoodsN[(h1<<2)|h2]  <<  std::endl;
            }
             //SEQAN_ASSERT_LT(abs(lHoods1[(h1<<2)|h2] - lHoods2[(h1<<2)|h2]), 0.001);
        }
    }
    */
}

/*
template <class T>
T chebyshev_coefficient(unsigned n, unsigned m)
{
   BOOST_MATH_STD_USING
   if(m > n)
      return 0;
   if((n & 1) != (m & 1))
      return 0;
   if(n == 0)
      return 1;
   T result = T(n) / 2;
   unsigned r = n - m;
   r /= 2;

   BOOST_ASSERT(n - 2 * r == m);

   if(r & 1)
      result = -result;
   result /= n - r;
   result *= boost::math::binomial_coefficient<T>(n - r, r);
   result *= ldexp(1.0f, m);
   return result;
}

template <class Seq>
Seq polynomial_to_chebyshev(const Seq& s)
{
   // Converts a Polynomial into Chebyshev form:
   //typedef typename Seq::value_type value_type;
   //typedef typename Seq::difference_type difference_type;
   Seq result(s);
   int order = length(s) - 1;
   int even_order = order & 1 ? order - 1 : order;
   int odd_order = order & 1 ? order : order - 1;

   for(int i = even_order; i >= 0; i -= 2)
   {
      long double val = s[i];
      for(int k = even_order; k > i; k -= 2)
      {
         val -= result[k] * chebyshev_coefficient<long double>(static_cast<unsigned>(k), static_cast<unsigned>(i));
      }
      val /= chebyshev_coefficient<long double>(static_cast<unsigned>(i), static_cast<unsigned>(i));
      result[i] = val;
   }
   result[0] *= 2;

   for(int i = odd_order; i >= 0; i -= 2)
   {
      long double val = s[i];
      for(int k = odd_order; k > i; k -= 2)
      {
         val -= result[k] * chebyshev_coefficient<long double>(static_cast<unsigned>(k), static_cast<unsigned>(i));
      }
      val /= chebyshev_coefficient<long double>(static_cast<unsigned>(i), static_cast<unsigned>(i));
      result[i] = val;
   }
   return result;
}

template <class Seq, class T>
T evaluate_chebyshev(const Seq& a, const T& x)
{
   // Clenshaw's formula:
   T yk2 = 0;
   T yk1 = 0;
   T yk = 0;
   for(int i = length(a) - 1; i >= 1; --i)
   {
      yk2 = yk1;
      yk1 = yk;
      yk = 2 * x * yk1 - yk2 + a[i];
   }
   return a[0] / 2 + yk * x - yk1;
}

SEQAN_DEFINE_TEST(test_snp_meth_store_polynom_calculation)
{

    // set up test case
    // Options & methOptions
    SNPCallingOptions<>     options;
    MethCallingOptions      methOptions;
    
    // Quality strings
    String<CharString> qualF;
    String<CharString> qualR;

    resize(qualF, 4, Exact());
    resize(qualR, 4, Exact());

    //qualF[2] = "KCJ]QLEGJDEL=HDRKGVKKIAFJ=O?EKKILEOHENHBFHGLGDHGGHIIHHGGFIIHJHHIIHHHH";
    qualR[0] = "HIIIIJGHIJHJGIMLJJHPHKJHKKJK??KBRC>OIU";
    //qualR[0] = "HIIIIJGHIJHJGIMLJJHPHKJHKKJKKBRCOIU";
    qualR[2] = "HHHHIHIHKIIFGKGJKCKFJDLJCPHOKALQFLHTIE";

    //qualF[2] = "HHEGH@FHGCNIHJIHHHHHI";
    //qualR[0] = "FJJLIOIGMJI";
    //qualR[2] = "HGHICJ";


    // Calculate candidate post probs with polynomial: TODO use function instead of copy and paste
    String<long double> singleProbs;
    resize(singleProbs, 6);

    String<long double> lHoods;  
    resize(lHoods, 4*4, 1.0, Exact()); 

    unsigned R = length(qualR[0]) + length(qualR[2]);  // number of mapped reads (Ns not counted)
    String<String<long double> > coeffSet;      // get polynom coeffs for each possible genotype polynom
    String<long double> coeffs;
    resize(coeffs, R+1, 0.0, Exact());
    resize(coeffSet, 4*4, coeffs, Exact());
    long double a;
    long double b;
    unsigned r = 0;

    long double lHood1 = 1.0;
    long double beta = 0.5;
    // for each observed base type
    for (unsigned i = 0; i < 4; ++i)        
    {
        //std::cout << "Test 0.0 " << std::endl;
        // for all reads mapped on forward strand
        for (unsigned j = 0; j < length(qualF[i]); ++j)
        {
            long double qual =  static_cast<long double>(ordValue(qualF[i][j])-33);
            if (qual < 1)   
            {
                continue;
            }
            long double e = pow(10.0, (long double)(-qual/10.0));
            // for each possible candidate genotype
            for (unsigned h = 0; h < 6; ++h)
                getSingleBaseProbHaploF(singleProbs[h], (Dna)i, (DnaM)h, e, methOptions);
           
            lHood1 *= (1-beta)*singleProbs[ordValue((Dna)'G')] + beta*singleProbs[ordValue((DnaM)'H')];
        }
        // for all reads mapped on reverse strand
        for (unsigned j = 0; j < length(qualR[i]); ++j)
        {
            long double qual =  static_cast<long double>(ordValue(qualR[i][j])-33);
            if (qual < 1)   
            {
                eraseBack(coeffSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')]);
                continue;
            }

            long double e = pow(10.0, (long double)(-qual/10.0));  
            // for each possible candidate genotype
            for (unsigned h = 0; h < 6; ++h)
                getSingleBaseProbHaploR(singleProbs[h], (Dna)i, (DnaM)h, e, methOptions);
             
            a = singleProbs[ordValue((Dna)'G')];
            b = -singleProbs[ordValue((Dna)'G')] + singleProbs[ordValue((DnaM)'H')];
            addFactorToPolynom(coeffSet[(ordValue((Dna)'G')<<2)|ordValue((Dna)'G')], a, b, r);

            ++r;
        }
    }
    std::cout << "coeffs: " << std::endl;
    for (unsigned i = 0; i < length(coeffSet[(ordValue((Dna)'G')<<2)|ordValue((Dna)'G')]); ++i)
    {
        std::cout << coeffSet[(ordValue((Dna)'G')<<2)|ordValue((Dna)'G')][i] << "  ,  ";
    }
    std::cout << std::endl;
    std::cout << "lHood1 for forward strand: " <<  lHood1 << std::endl;
    std::cout << "lHood1 for reverse strand: (Horner)" <<  getPolynomialValue(beta, coeffSet[(ordValue((Dna)'G')<<2)|ordValue((Dna)'G')]) << std::endl;
    std::cout << "lHood1 for reverse strand: (Clenshaw) " <<  getPolynomialValueC(beta, coeffSet[(ordValue((Dna)'G')<<2)|ordValue((Dna)'G')]) << std::endl;

    String<double> chebys;
    chebys = polynomial_to_chebyshev(coeffSet[(ordValue((Dna)'G')<<2)|ordValue((Dna)'G')]);
    std::cout << "lHood1 for reverse strand: (chebyshe): " <<  evaluate_chebyshev(chebys, beta) << std::endl;


    lHood1 *= getPolynomialValue(beta, coeffSet[(ordValue((Dna)'G')<<2)|ordValue((Dna)'G')]); 

    ///////////////////////////////////////////////////////////////////
    // Compute likelihoods for given beta values without polynom
    ///////////////////////////////////////////////////////////////////
    long double lHood2 = 1.0;
    long double lHood2_F = 1.0;
    long double lHood2_R = 1.0;
    // for each observed base type
    for (unsigned i = 0; i < 4; ++i)        
    {
        //std::cout << "Test 0.0 " << std::endl;
        // for all reads mapped on forward strand
        for (unsigned j = 0; j < length(qualF[i]); ++j)
        {
            long double qual =  static_cast<long double>(ordValue(qualF[i][j])-33);
            if (qual < 1)   
            {
                continue;
            }
            long double e = pow(10.0, (long double)(-qual/10.0));
            // for each possible candidate genotype
            for (unsigned h = 0; h < 6; ++h)
                getSingleBaseProbHaploF(singleProbs[h], (Dna)i, (DnaM)h, e, methOptions);
            
            lHood2 *= (1-beta)*singleProbs[ordValue((Dna)'G')] + beta*singleProbs[ordValue((DnaM)'H')];
            lHood2_F *= (1-beta)*singleProbs[ordValue((Dna)'G')] + beta*singleProbs[ordValue((DnaM)'H')];

        }
        // for all reads mapped on reverse strand
        for (unsigned j = 0; j < length(qualR[i]); ++j)
        {
            long double qual =  static_cast<long double>(ordValue(qualR[i][j])-33);
            if (qual < 1)   
            {
                continue;
            }
            long double e = pow(10.0, (long double)(-qual/10.0));  
            // for each possible candidate genotype
            for (unsigned h = 0; h < 6; ++h)
                getSingleBaseProbHaploR(singleProbs[h], (Dna)i, (DnaM)h, e, methOptions);
             
            std::cout << std::setprecision (50) << " probs G: " <<  singleProbs[ordValue((Dna)'G')] << std::endl;
            std::cout << std::setprecision (50) << " probs H: " <<  singleProbs[ordValue((DnaM)'H')] << std::endl;
            std::cout << std::setprecision (10) << " probs G: " <<  (1.0-beta)*singleProbs[ordValue((Dna)'G')] << std::endl;
            std::cout << std::setprecision (10) << " probs H: " <<  beta*singleProbs[ordValue((DnaM)'H')] << std::endl;
            std::cout << " result: " <<  ((long double)0.42)*singleProbs[ordValue((Dna)'G')] + ((long double)0.58)*singleProbs[ordValue((DnaM)'H')] << std::endl;
            std::cout << ((long double)0.42)*singleProbs[ordValue((Dna)'G')] + ((long double)0.58)*singleProbs[ordValue((DnaM)'H')] << std::cout << std::endl;
            

            lHood2 *= (1.0-beta)*singleProbs[ordValue((Dna)'G')] + beta*singleProbs[ordValue((DnaM)'H')];
            lHood2_R *= (1.0-beta)*singleProbs[ordValue((Dna)'G')] + beta*singleProbs[ordValue((DnaM)'H')];
        }
    }
    std::cout << "lHood2 for forward strand: " <<  lHood2_F << std::endl;
    std::cout << "lHood2 for reverse strand: " <<  lHood2_R << std::endl;


    std::cout << std::setprecision (25) << "lHood1: " << lHood1 << "  lHood2: " << lHood2 << std::endl;
    // Check, if likelihoods are more or less the same
    SEQAN_ASSERT_LT(abs(lHood1 - lHood2), 0.0001);

   // Test clenshaw:
   clear(coeffs);
   appendValue(coeffs, 3.0);
   appendValue(coeffs, 4.0);
   appendValue(coeffs, 2.0);
   appendValue(coeffs, 1.0);

   long double z = 1.0;
    std::cout << "result: (Horner)" <<  getPolynomialValue(z, coeffs) << std::endl;
    std::cout << "result: (Clenshaw) " <<  getPolynomialValueC(z, coeffs) << std::endl;

    
    chebys = polynomial_to_chebyshev(coeffs);
    std::cout << "result: (chebyshe): " <<  evaluate_chebyshev(chebys, z) << std::endl;

}
*/

#endif  // SANDBOX_KRAKAU_TESTS_SNP_METH_STORE_TEST_SNP_METH_STORE_H_
