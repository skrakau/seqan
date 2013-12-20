#ifndef SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_VS_SNPSTORE_H
#define SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_VS_SNPSTORE_H

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>


using namespace seqan;

// Return genotype base for pseudo(?) iupac code used in snpStore
void castGenotype(Dna5 &a1, Dna5 &a2, Iupac &iu)
{
    if (iu == 'A')
    {
        a1 = 'A';
        a2 = 'A';
    }
    else if (iu == 'M') 
    {
        a1 = 'A';
        a2 = 'C';
    }
    else if (iu == 'R') 
    {
        a1 = 'A';
        a2 = 'G';
    }
    else if (iu == 'W') 
    {
        a1 = 'A';
        a2 = 'T';
    }
    else if (iu == 'C') 
    {
        a1 = 'C';
        a2 = 'C';
    }
    else if (iu == 'S') 
    {
        a1 = 'C';
        a2 = 'G';
    }
    else if (iu == 'Y') 
    {
        a1 = 'C';
        a2 = 'T';
    }
    else if (iu == 'G') 
    {
        a1 = 'G';
        a2 = 'G';
    }
    else if (iu == 'K') 
    {
        a1 = 'G';
        a2 = 'T';
    }
    else if (iu == 'T') 
    {
        a1 = 'T';
        a2 = 'T';
    }
}

template<typename TOptions>
bool benchmark_vs_snpstore(TOptions &options)
{
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;

    std::fstream simFile(toCString(options.simSnpsFile), std::ios::binary | std::ios::in);
    TRecordReader simReader(simFile);

    std::fstream callFile(toCString(options.calledSnpsFile), std::ios::binary | std::ios::in);
    TRecordReader callReader(callFile);

    // snpStore - ss
    std::fstream snpStoreFile(toCString(options.snpStoreFile), std::ios::binary | std::ios::in);
    TRecordReader ssReader(snpStoreFile);

    // First skip lines with metainforamation
    CharString helper = "#"; 
    while(!atEnd(simReader) && helper[0] == '#')
    {
       clear(helper);
       readNChars(helper, simReader, 1);
       skipLine(simReader);
    }
    helper = "#"; 
    while(!atEnd(callReader) && helper[0] == '#')
    {
       clear(helper);
       readNChars(helper, callReader, 1);
       skipLine(callReader);
    }
    helper = "#"; 
    while(!atEnd(ssReader) && helper[0] == '#')
    {
       clear(helper);
       readNChars(helper, ssReader, 1);
       skipLine(ssReader);
    }

    
    CharString simContig;
    CharString callContig;
    CharString ssContig;


    CharString simPos;
    CharString callPos;
    CharString ssPos;
    CharString callCoverage;
    CharString ssCoverage;

    
    Dna5 simRef;
    Dna5 callRef;
    Dna5 ssRef;


    Dna5 simA1;
    Dna5 simA2;
    Dna5 callA1;
    Dna5 callA2;
    Dna5 ssA1;
    Dna5 ssA2;

    unsigned countSimSnps = 0;
    unsigned countCalledSnps = 0;
    unsigned countNotListed = 0;    // Not found snps
    unsigned countNotCalled = 0;    // 
    unsigned countWrongCalled = 0;
    unsigned countRightCalled = 0;
    unsigned countFalsePositive = 0;    // False found snps

    unsigned countssNotListed = 0;
    unsigned countssNotCalled = 0;
    unsigned countssWrongCalled = 0;
    unsigned countssRightCalled = 0;
    unsigned countssFalsePositive = 0;

    unsigned countBothNotListed = 0;
    unsigned countBothNotCalled = 0;
    unsigned countBothWrongCalled = 0;
    unsigned countBothRightCalled = 0;
    unsigned countBothFalsePositive = 0;

    unsigned countSameWrongCall = 0;


    clear(simContig);
    clear(callContig);
    clear(ssContig);
    clear(simPos);
    clear(callPos);
    clear(ssPos);

    readUntilWhitespace(simContig, simReader);
    skipWhitespaces(simReader);
    readUntilWhitespace(simPos, simReader);

    readUntilWhitespace(callContig, callReader);
    skipWhitespaces(callReader);
    readUntilWhitespace(callPos, callReader);

    readUntilWhitespace(ssContig, ssReader);
    skipWhitespaces(ssReader);
    readUntilWhitespace(ssPos, ssReader);


    // get rid of record reader, doesnt make sense, only for sequence files etc...
    // TODO go till end of other files...
    while(!atEnd(simReader) && !atEnd(callReader) && !atEnd(ssReader))
    {
        //std::cout << "SimPos: " << simPos << "  callPos: " << callPos << std::endl;
        if (simContig != callContig || simContig != ssContig)
        {
            // ja was?
            // skip and go to position as in other cases
            // continue
        }

        if (lexicalCast<size_t>(simPos) < lexicalCast<size_t>(callPos) && lexicalCast<size_t>(simPos) < lexicalCast<size_t>(ssPos)) // simSnp not found in both
        {
            ++countNotListed;
            ++countSimSnps;
            ++countssNotListed;
            ++countBothNotListed;

            //if (simPos > 13000000 && simPos < 13000000)
            //{
                //std::cout << "Not listed at simPos: " << simPos <<  std::endl;
            //}

            skipLine(simReader);                        // skip
            clear(simContig);
            readUntilWhitespace(simContig, simReader);  // chrom
            skipWhitespaces(simReader);
            clear(simPos);
            readUntilWhitespace(simPos, simReader);     // pos
            continue;
        }
        else if (lexicalCast<size_t>(callPos) < lexicalCast<size_t>(simPos) && lexicalCast<size_t>(callPos) < lexicalCast<size_t>(ssPos))   // false positive snp in snp_meth_store, but not in snpStore
        {
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // ref
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // cov
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // genotype called
            if (helper != '.')
            {
                callA1 = helper[0];
                callA2 = helper[1];
                ++countFalsePositive;
                ++countCalledSnps;

                //if (countFalsePositive > 17000)
                    //std::cout << "False positive at callPos: " << callPos << std::endl;
            }

            skipLine(callReader);                           // skip
            clear(callContig);
            readUntilWhitespace(callContig, callReader);    // chrom
            skipWhitespaces(callReader);
            clear(callPos);
            readUntilWhitespace(callPos, callReader);       // pos
            continue;
        }
        else if (lexicalCast<size_t>(ssPos) < lexicalCast<size_t>(simPos) && lexicalCast<size_t>(ssPos) < lexicalCast<size_t>(callPos))     // false postive snp in snpStore, but not in snp_meth_store
        {
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // ref
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualA
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualC
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualG
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualT
            skipWhitespaces(ssReader);
            clear(ssCoverage);
            readUntilWhitespace(ssCoverage, ssReader);        // cov
            readNChars(helper, ssReader, 1);
            clear(helper);
            readNChars(helper, ssReader, 1);        // genotype called, if called
            if (isalpha(helper[0]))
            {
                Iupac i = helper[0];
                castGenotype(ssA1, ssA2, i);
                ++countssFalsePositive;
                if (countssFalsePositive < 10)
                    std::cout << " FP in snpStore: " << ssA1 << " " << ssA2 << " coverage " << ssCoverage <<  std::endl;
            }

            skipLine(ssReader);                           // skip
            clear(ssContig);
            readUntilWhitespace(ssContig, ssReader);    // chrom
            skipWhitespaces(ssReader);
            clear(ssPos);
            readUntilWhitespace(ssPos, ssReader);       // pos
            continue;
        }
        else if (lexicalCast<size_t>(ssPos) < lexicalCast<size_t>(simPos) && lexicalCast<size_t>(ssPos) == lexicalCast<size_t>(callPos))     // false postive snp in both
        {
            // callReader
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // ref
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // cov
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // genotype called
            bool bothFP = false;
            if (helper != '.')
            {
                callA1 = helper[0];
                callA2 = helper[1];
                ++countFalsePositive;
                ++countCalledSnps;
                bothFP = true;
            }

            skipLine(callReader);                           // skip
            clear(callContig);
            readUntilWhitespace(callContig, callReader);    // chrom
            skipWhitespaces(callReader);
            clear(callPos);
            readUntilWhitespace(callPos, callReader);       // pos

            // ssReader 
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // ref
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualA
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualC
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualG
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualT
            skipWhitespaces(ssReader);
            clear(ssCoverage);
            readUntilWhitespace(ssCoverage, ssReader);        // cov
            readNChars(helper, ssReader, 1);
            clear(helper);
            readNChars(helper, ssReader, 1);        // genotype called, if called
            if (isalpha(helper[0]))
            {
                Iupac i = helper[0];
                castGenotype(ssA1, ssA2, i);
                ++countssFalsePositive;
                if (bothFP)
                    ++countBothFalsePositive;
            }

            skipLine(ssReader);                           // skip
            clear(ssContig);
            readUntilWhitespace(ssContig, ssReader);    // chrom
            skipWhitespaces(ssReader);
            clear(ssPos);
            readUntilWhitespace(ssPos, ssReader);       // pos

            continue;
        }
        
        // snp found somewhere
        if (lexicalCast<size_t>(simPos) == lexicalCast<size_t>(callPos) && lexicalCast<size_t>(simPos) < lexicalCast<size_t>(ssPos) )       // simSnp only found in snp_meth_store, but not in snpStore
        {
            ++countssNotListed;
            // simReader
            skipWhitespaces(simReader);
            skipChar(simReader, '.');
            skipWhitespaces(simReader);
            clear(helper);
            readUntilWhitespace(helper, simReader);     // ref
            simRef = helper[0];
            skipWhitespaces(simReader);
            clear(helper);
            readUntilWhitespace(helper, simReader);     // genotype called
            if (length(helper) == 1)
            {
                simA1 = helper[0];
                // Check if hetero or homo snp
                skipWhitespaces(simReader);
                skipChar(simReader, '.');
                skipWhitespaces(simReader);
                skipUntilWhitespace(simReader);
                skipWhitespaces(simReader);
                clear(helper);
                readUntilWhitespace(helper, simReader);
                if (prefix(helper, 3) == "AF=")
                {
                    if (suffix(helper, 3) == "50")          // Create hetero genotype 
                    {
                        if (simA1 < simRef)
                            simA2 = simRef;
                        else
                        {
                            simA2 = simA1;
                            simA1 = simRef;
                        }
                    }
                    else if (suffix(helper, 3) == "100")    // Create homo genotype
                        simA2 = simA1;                   
                    else
                    {
                        std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simPos << std::endl;
                        return 1;
                    }
                }
                else
                {
                    std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simPos << std::endl;
                    return 1;
                }
            }
            else if (length(helper) == 3)
            {
                simA1 = helper[0];
                simA2 = helper[2];
            }
            else if (length(helper) != 0)
            {
                std::cerr << "ERROR: Something wrong in simulated SNPs file with ALT entry at position: " << simPos << std::endl;
                return 1;
            }

            ++countSimSnps;
            skipLine(simReader);                        // skip
            clear(simContig);
            readUntilWhitespace(simContig, simReader);  // chrom
            skipWhitespaces(simReader);
            clear(simPos);
            readUntilWhitespace(simPos, simReader);     // pos

            // callReader
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // ref
            skipWhitespaces(callReader);
            clear(callCoverage);
            readUntilWhitespace(callCoverage, callReader);        // cov
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // genotype called 
            if (helper != '.')
            {
                callA1 = helper[0];
                callA2 = helper[1];
                ++countCalledSnps;
                // compare called genotypes
                if (simA1 != callA1 || simA2 != callA2)
                {
                    ++countWrongCalled;
                    //if (options.verbosity >= 4)
                    if (countWrongCalled < 70)
                    {
                        std::cout << "Wrong called at snp_meth position (only found in snp_meth): " << callPos << std::endl;
                        std::cout << "simA1: " << simA1 << " simA2: " << simA2 << " callA1: " << callA1 << " callA2: " << callA2 << std::endl;
                    }
                }
                else
                {
                    ++countRightCalled;
                   // std::cout << "Right called at callPos: " << callPos << std::endl;
                }
            }
            else
            {
                ++countNotCalled;
                if (options.verbosity >= 4)
                {
                    std::cout << "Not called at position: " << callPos << " coverage: " << callCoverage << std::endl;
                }
            }
            skipLine(callReader);
            clear(callContig);
            readUntilWhitespace(callContig, callReader);    // chrom
            skipWhitespaces(callReader);
            clear(callPos);
            readUntilWhitespace(callPos, callReader);       // pos

            continue;
        }
        else if (lexicalCast<size_t>(simPos) < lexicalCast<size_t>(callPos) && lexicalCast<size_t>(simPos) == lexicalCast<size_t>(ssPos) )  // simSnp only found in SnpStore, but not in snp_meth_store 
        {
            ++countNotListed; // in snp_meth_store
            // simReader
            skipWhitespaces(simReader);
            skipChar(simReader, '.');
            skipWhitespaces(simReader);
            clear(helper);
            readUntilWhitespace(helper, simReader);     // ref
            simRef = helper[0];
            skipWhitespaces(simReader);
            clear(helper);
            readUntilWhitespace(helper, simReader);     // genotype called
            if (length(helper) == 1)
            {
                simA1 = helper[0];
                // Check if hetero or homo snp
                skipWhitespaces(simReader);
                skipChar(simReader, '.');
                skipWhitespaces(simReader);
                skipUntilWhitespace(simReader);
                skipWhitespaces(simReader);
                clear(helper);
                readUntilWhitespace(helper, simReader);
                if (prefix(helper, 3) == "AF=")
                {
                    if (suffix(helper, 3) == "50")          // Create hetero genotype 
                    {
                        if (simA1 < simRef)
                            simA2 = simRef;
                        else
                        {
                            simA2 = simA1;
                            simA1 = simRef;
                        }
                    }
                    else if (suffix(helper, 3) == "100")    // Create homo genotype
                        simA2 = simA1;                   
                    else
                    {
                        std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simPos << std::endl;
                        return 1;
                    }
                }
                else
                {
                    std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simPos << std::endl;
                    return 1;
                }
            }
            else if (length(helper) == 3)
            {
                simA1 = helper[0];
                simA2 = helper[2];
            }
            else if (length(helper) != 0)
            {
                std::cerr << "ERROR: Something wrong in simulated SNPs file with ALT entry at position: " << simPos << std::endl;
                return 1;
            }

            ++countSimSnps;
            skipLine(simReader);                        // skip
            clear(simContig);
            readUntilWhitespace(simContig, simReader);  // chrom
            skipWhitespaces(simReader);
            clear(simPos);
            readUntilWhitespace(simPos, simReader);     // pos

            // ssReader
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // ref
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualA
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualC
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualG
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualT
            skipWhitespaces(ssReader);
            clear(ssCoverage);
            readUntilWhitespace(ssCoverage, ssReader);        // cov
            readNChars(helper, ssReader, 1);
            clear(helper);
            readNChars(helper, ssReader, 1);        // genotype called, if called
            if (isalpha(helper[0]))
            {
                Iupac i = helper[0];
                castGenotype(ssA1, ssA2, i);
                if (simA1 != ssA1 || simA2 != ssA2)
                {
                    ++countssWrongCalled;
                    //if (options.verbosity >= 5)
                    if (countssWrongCalled < 1)
                    {
                        std::cout << "Wrong called (only) in SnpStore at position: " << ssPos << std::endl;
                        std::cout << "simA1: " << simA1 << " simA2: " << simA2 << " ssA1: " << ssA1 << " ssA2: " << ssA2 << " i " << i <<  " helper " << helper[0] << std::endl;
                    }
               }
                else
                {
                    ++countssRightCalled;
                }
            }
            else
                ++countssNotCalled;


            skipLine(ssReader);                           // skip
            clear(ssContig);
            readUntilWhitespace(ssContig, ssReader);    // chrom
            skipWhitespaces(ssReader);
            clear(ssPos);
            readUntilWhitespace(ssPos, ssReader);       // pos
            continue;
        }
        else if (lexicalCast<size_t>(simPos) == lexicalCast<size_t>(callPos) && lexicalCast<size_t>(simPos) == lexicalCast<size_t>(ssPos) )     // simSnp found in both
        {
            bool bothWrong = false;
            bool bothNotCalled = false;

            // simReader
            skipWhitespaces(simReader);
            skipChar(simReader, '.');
            skipWhitespaces(simReader);
            clear(helper);
            readUntilWhitespace(helper, simReader);     // ref
            simRef = helper[0];
            skipWhitespaces(simReader);
            clear(helper);
            readUntilWhitespace(helper, simReader);     // genotype called
            if (length(helper) == 1)
            {
                simA1 = helper[0];
                // Check if hetero or homo snp
                skipWhitespaces(simReader);
                skipChar(simReader, '.');
                skipWhitespaces(simReader);
                skipUntilWhitespace(simReader);
                skipWhitespaces(simReader);
                clear(helper);
                readUntilWhitespace(helper, simReader);
                if (prefix(helper, 3) == "AF=")
                {
                    if (suffix(helper, 3) == "50")          // Create hetero genotype 
                    {
                        if (simA1 < simRef)
                            simA2 = simRef;
                        else
                        {
                            simA2 = simA1;
                            simA1 = simRef;
                        }
                    }
                    else if (suffix(helper, 3) == "100")    // Create homo genotype
                        simA2 = simA1;                   
                    else
                    {
                        std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simPos << std::endl;
                        return 1;
                    }
                }
                else
                {
                    std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simPos << std::endl;
                    return 1;
                }
            }
            else if (length(helper) == 3)
            {
                simA1 = helper[0];
                simA2 = helper[2];
            }
            else if (length(helper) != 0)
            {
                std::cerr << "ERROR: Something wrong in simulated SNPs file with ALT entry at position: " << simPos << std::endl;
                return 1;
            }

            ++countSimSnps;
            skipLine(simReader);                        // skip
            clear(simContig);
            readUntilWhitespace(simContig, simReader);  // chrom
            skipWhitespaces(simReader);
            clear(simPos);
            readUntilWhitespace(simPos, simReader);     // pos

            // callReader
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // ref
            skipWhitespaces(callReader);
            clear(callCoverage);
            readUntilWhitespace(callCoverage, callReader);        // cov
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // genotype called 
            if (helper != '.')
            {
                callA1 = helper[0];
                callA2 = helper[1];
                ++countCalledSnps;
                // compare called genotypes
                if (simA1 != callA1 || simA2 != callA2)
                {
                    ++countWrongCalled;
                    bothWrong = true;
                    //if (options.verbosity >= 4)
                    if (countWrongCalled < 100)
                    {
                        std::cout << "Wrong called at snp_meth position (also something found in snpStore): " << callPos << std::endl;
                        std::cout << "simA1: " << simA1 << " simA2: " << simA2 << " callA1: " << callA1 << " callA2: " << callA2 << std::endl;
                    }
                }
                else
                {
                    ++countRightCalled;
                   // std::cout << "Right called at callPos: " << callPos << std::endl;
                }
            }
            else
            {
                ++countNotCalled;
                bothNotCalled = true;
                if (options.verbosity >= 4)
                {
                    std::cout << "Not called at position: " << callPos << " coverage: " << callCoverage << std::endl;
                }
            }
            skipLine(callReader);
            clear(callContig);
            readUntilWhitespace(callContig, callReader);    // chrom
            skipWhitespaces(callReader);
            clear(callPos);
            readUntilWhitespace(callPos, callReader);       // pos

            // ssReader
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // ref
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualA
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualC
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualG
            skipWhitespaces(ssReader);
            readUntilWhitespace(helper, ssReader);        // qualT
            skipWhitespaces(ssReader);
            clear(ssCoverage);
            readUntilWhitespace(ssCoverage, ssReader);        // cov
            readNChars(helper, ssReader, 1);
            clear(helper);
            readNChars(helper, ssReader, 1);        // genotype called, if called
            if (isalpha(helper[0]))
            {
                Iupac i = helper[0];
                // attention: not normal iupac coding..
                // toIupac = "AMRWMCSYRSGKWYKT";
                castGenotype(ssA1, ssA2, i);
                //std::cout << " ssA1: " << ssA1 << " ssA2: " << ssA2 <<  std::endl;

                if (simA1 != ssA1 || simA2 != ssA2)
                {
                    ++countssWrongCalled;
                    if (!bothNotCalled && ssA1 == callA1 && ssA2 == callA2)
                        ++countSameWrongCall;
 
                    if (countWrongCalled < 100)
                    {
                        if (!bothNotCalled && ssA1 == callA1 && ssA2 == callA2)
                        {
                            std::cout << "Same wrong call in both " << ssPos << std::endl;
                            std::cout << " ssA1: " << ssA1 << " ssA2: " << ssA2 <<  std::endl;
                        }
                    }

                    if (bothWrong)
                        ++countBothWrongCalled;
                }
                else
                {
                    ++countssRightCalled;
                    if (ssA1 == callA1 && ssA2 == callA2)
                        ++countBothRightCalled;
                }
            }
            else
            {
                ++countssNotCalled;
                if (bothNotCalled) 
                    ++countBothNotCalled;
            }

            skipLine(ssReader);                           // skip
            clear(ssContig);
            readUntilWhitespace(ssContig, ssReader);    // chrom
            skipWhitespaces(ssReader);
            clear(ssPos);
            readUntilWhitespace(ssPos, ssReader);       // pos
            continue;
        }


    }

    std::cout << "Simulated snps: \t" << countSimSnps << std::endl;
    //std::cout << "Called snps: \t" << countCalledSnps << std::endl;
    std::cout << std::endl;
    std::cout << "Snp_meth_store: " << std::endl;
    std::cout << "Not listed: \t" << countNotListed << std::endl;
    std::cout << "Not called: \t" << countNotCalled << std::endl;
    std::cout << "False positive: \t" << countFalsePositive << std::endl;
    std::cout << "Wrong called genotype: \t" << countWrongCalled << std::endl; 
    std::cout << "Right called genotype: \t" << countRightCalled << std::endl;
    std::cout << std::endl;
    std::cout << "SnpStore:" << std::endl;
    std::cout << "Not listed: \t" << countssNotListed << std::endl;
    std::cout << "Not called: \t" << countssNotCalled << std::endl;
    std::cout << "False positive: \t" << countssFalsePositive << std::endl;
    std::cout << "Wrong called genotype: \t" << countssWrongCalled << std::endl; 
    std::cout << "Right called genotype: \t" << countssRightCalled << std::endl;
    std::cout << std::endl;
    std::cout << "Both:" << std::endl;
    std::cout << "Not listed: \t" << countBothNotListed << std::endl;
    std::cout << "Not called: \t" << countBothNotCalled << std::endl;
    std::cout << "False positive: \t" << countBothFalsePositive << std::endl;
    std::cout << "Wrong called genotype: \t" << countBothWrongCalled << std::endl; 
    std::cout << "Right called genotype: \t" << countBothRightCalled << std::endl;

    std::cout << std::endl;
    std::cout << "Same Wrong calls in both: " << countSameWrongCall << std::endl;
    return 0;
}













    

 
#endif
