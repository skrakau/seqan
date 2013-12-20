// ==========================================================================
// 
//                           simulate_illumina_bs.h
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
// Author: Sabrina Krakau <sabrina.krakau@fu-berlin.de>
// ==========================================================================
//
//
// SIMULATE METHYLATION PATTERNS DEPENDING ON CONTEXT 
// OR
// BASED ON GIVEN METHYLATION RATES
//

#ifndef SANDBOX_KRAKAU_APPS_BS_MASON_SIMULATE_ILLUMINA_BS_H_
#define SANDBOX_KRAKAU_APPS_BS_MASON_SIMULATE_ILLUMINA_BS_H_

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/beta.hpp>

#include <seqan/store.h>

#include "bs_mason.h"

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================
 
struct IlluminaReadsBS_;
typedef Tag<IlluminaReadsBS_> IlluminaReadsBS;

template<>
struct Options<IlluminaReadsBS> : public Options<IlluminaReads>
{
    double probabilityMethylationCG;
    double probabilityMethylationCHG;
    double probabilityMethylationCHH;
    // Protocol: Lister (0), Cokus (1) 
    unsigned sequencingProtocol;
    // Rate of unmethylated Cs converted to Ts
    double conversionRate;
    // For haplotype creation: 
    // rate of unmethylated Cs becoming methylated
    double haplotypeMethylInsRate;
    // rate of methylated Cs becoming unmethylated
    double haplotypeMethylDelRate;

    // For import of given methylation rates
    bool useMethylRates;
    CharString methylRatesFile;   

    // For simulation of methylation level
    // p, q values for beta distributions
    double pCG;
    double qCG;
    double pCHG;
    double qCHG;
    double pCHH;
    double qCHH;
    double thresholdMethLevel;


    bool writeVCFFile;  // Output simulated snps, only MISMATCHES, for max. 2 haplotypes, internal option at the moment
    bool writeMethFile;
    // TODO: check, if vcf file given
    CharString vcfOutputFile;
    CharString methOutputFile;

    Options():
                // Context dependent methylation probabilities:
                probabilityMethylationCG(0.8),
                probabilityMethylationCHG(0.1),
                probabilityMethylationCHH(0.01),
                // 
                sequencingProtocol(0),
                //
                conversionRate(0.98),
                //
                haplotypeMethylInsRate(0.001),  // TODO: normaly 0.0001 ? contex dependent?
                haplotypeMethylDelRate(0.001),  
                //
                useMethylRates(false),
                //
                pCG(0.22),
                qCG(0.715),
                pCHG(0.5),
                qCHG(1.0),
                pCHH(0.4),
                qCHH(1.0),
                thresholdMethLevel(0.5),
                //
                writeVCFFile(true),
                writeMethFile(true)
    {}                
};

template<>
struct ModelParameters<IlluminaReadsBS> : public ModelParameters<IlluminaReads>
{
    double thresholdYCG;
    double thresholdYCHG;
    double thresholdYCHH;

    // For import of given methylation rates:
    StringSet<String<double> > topMethylRates; 
    StringSet<String<double> > bottomMethylRates;
    
    // StringSet: corresponding to contigs in fragmentStore
    StringSet<String<double, Journaled<Alloc<> > > > topMethylRatesHaplotype;            
    StringSet<String<double, Journaled<Alloc<> > > > bottomMethylRatesHaplotype;

};

template <>
struct ReadSimulationInstruction<IlluminaReadsBS> : public ReadSimulationInstruction<IlluminaReads> 
{
    // Bs conversion String: regarding original read before base call errors (editString) are simulated
    // At the beginning: contains FALSE where methylations in reference prevent bs conversions
    // After applySimulationInstructions: contains TRUE where bs conversions where applied (only C or G positions)
    String<bool> bsConversionString;
    // Wheather or not the read is coming from the original top strand 
    bool isFromTopStrand;

};

// ============================================================================
// Metafunctions.
// ============================================================================

// ============================================================================
// Functions.
// ============================================================================

// Extract Options<IlluminaReads> out of Options<IlluminaReadsBS>
void extractOptionsBase(Options<IlluminaReads> & options, Options<IlluminaReadsBS> const & optionsBS)
{
    assign(options.readLength, optionsBS.readLength);
    assign(options.probabilityInsert, optionsBS.probabilityInsert);
    assign(options.probabilityDelete, optionsBS.probabilityDelete);
    assign(options.probabilityMismatchFromFile, optionsBS.probabilityMismatchFromFile);
    assign(options.probabilityMismatchScale, optionsBS.probabilityMismatchScale);
    assign(options.probabilityMismatch, optionsBS.probabilityMismatch);
    assign(options.probabilityMismatchBegin, optionsBS.probabilityMismatchBegin);
    assign(options.probabilityMismatchEnd, optionsBS.probabilityMismatchEnd);
    assign(options.positionRaise, optionsBS.positionRaise);
    assign(options.illuminaNoN, optionsBS.illuminaNoN);
    assign(options.meanQualityBegin, optionsBS.meanQualityBegin);
    assign(options.meanQualityEnd, optionsBS.meanQualityEnd);
    assign(options.stdDevQualityBegin, optionsBS.stdDevQualityBegin);
    assign(options.stdDevQualityEnd, optionsBS.stdDevQualityEnd);
    assign(options.meanMismatchQualityBegin, optionsBS.meanMismatchQualityBegin);
    assign(options.meanMismatchQualityEnd, optionsBS.meanMismatchQualityEnd);
    assign(options.stdDevMismatchQualityBegin, optionsBS.stdDevMismatchQualityBegin);
    assign(options.stdDevMismatchQualityEnd, optionsBS.stdDevMismatchQualityEnd);
}
// Extract ModelParameters <IlluminaReads> out of ModelParameters <IlluminaReadsBS>
void extractModelParametersBase(ModelParameters<IlluminaReads> & parameters, ModelParameters<IlluminaReadsBS> const & parametersBS)
{
    assign(parameters.mismatchProbabilities, parametersBS.mismatchProbabilities);
    assign(parameters.mismatchQualityMeans, parametersBS.mismatchQualityMeans);
    assign(parameters.mismatchQualityStdDevs, parametersBS.mismatchQualityStdDevs);
    assign(parameters.qualityMeans, parametersBS.qualityMeans);
    assign(parameters.qualityStdDevs, parametersBS.qualityStdDevs);
}
// Merge Options<IlluminaReads> into Options<IlluminaReadsBS>
void mergeBaseIntoBSOptions(Options<IlluminaReadsBS> & optionsBS, Options<IlluminaReads> const & options)
{
    assign(optionsBS.readLength, options.readLength);
    assign(optionsBS.probabilityInsert, options.probabilityInsert);
    assign(optionsBS.probabilityDelete, options.probabilityDelete);
    assign(optionsBS.probabilityMismatchFromFile, options.probabilityMismatchFromFile);
    assign(optionsBS.probabilityMismatchScale, options.probabilityMismatchScale);
    assign(optionsBS.probabilityMismatch, options.probabilityMismatch);
    assign(optionsBS.probabilityMismatchBegin, options.probabilityMismatchBegin);
    assign(optionsBS.probabilityMismatchEnd, options.probabilityMismatchEnd);
    assign(optionsBS.positionRaise, options.positionRaise);
    assign(optionsBS.illuminaNoN, options.illuminaNoN);
    assign(optionsBS.meanQualityBegin, options.meanQualityBegin);
    assign(optionsBS.meanQualityEnd, options.meanQualityEnd);
    assign(optionsBS.stdDevQualityBegin, options.stdDevQualityBegin);
    assign(optionsBS.stdDevQualityEnd, options.stdDevQualityEnd);
    assign(optionsBS.meanMismatchQualityBegin, options.meanMismatchQualityBegin);
    assign(optionsBS.meanMismatchQualityEnd, options.meanMismatchQualityEnd);
    assign(optionsBS.stdDevMismatchQualityBegin, options.stdDevMismatchQualityBegin);
    assign(optionsBS.stdDevMismatchQualityEnd, options.stdDevMismatchQualityEnd);
}
// Merge ModelParameters<IlluminaReads> into ModelParameters<IlluminaReadsBS>
void mergeBaseIntoBSModelParameters(ModelParameters<IlluminaReadsBS> & parametersBS, ModelParameters<IlluminaReads> & parameters)
{
    assign(parametersBS.mismatchProbabilities, parameters.mismatchProbabilities);
    assign(parametersBS.mismatchQualityMeans, parameters.mismatchQualityMeans);
    assign(parametersBS.mismatchQualityStdDevs, parameters.mismatchQualityStdDevs);
    assign(parametersBS.qualityMeans, parameters.qualityMeans);
    assign(parametersBS.qualityStdDevs, parameters.qualityStdDevs);
}

template <typename TStream>
TStream & operator<<(TStream & stream, Options<IlluminaReadsBS> const & options) 
{
    stream << static_cast<Options<IlluminaReads> >(options);
    stream << "bisulfite-options {" << std::endl
           << "  probabilityMethylationCG;               " << options.probabilityMethylationCG << std::endl
           << "  probabilityMethylationCHG;              " << options.probabilityMethylationCHG << std::endl           
           << "  probabilityMethylationCHH;              " << options.probabilityMethylationCHH<< std::endl
           << "  sequencingProtocol;                     " << options.sequencingProtocol << std::endl
           << "  conversionRate;                         " << options.conversionRate << std::endl
           << "  haplotypeMethylInsRate;                 " << options.haplotypeMethylInsRate << std::endl
           << "  haplotypeMethylDelRate;                 " << options.haplotypeMethylDelRate << std::endl
           << "  useMethylRates:                    " << (options.useMethylRates ? "true" : "false") << std::endl
           << "  methylRatesFile:                 \"" << options.methylRatesFile << "\"" << std::endl
           << "  writeVCFFile:                          " << options.writeVCFFile << std::endl
           << "  writeMethFile:                          " << options.writeMethFile << std::endl
           << "  vcfOutputFile:                          " << options.vcfOutputFile << std::endl
           << "  methOutputFile:                          " << options.methOutputFile << std::endl
           << "}" << std::endl;
    return stream;
}

void setUpArgumentParser(ArgumentParser & parser,
                            IlluminaReadsBS const &)
{
    setUpArgumentParser(parser, IlluminaReads());

    addSection(parser, "Illumina bisulfite parameters");

    addOption(parser, ArgParseOption("pcg",  "prob-cg", "The probability of an 'C' beeing methylated in CG-context.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "pcg", "0.25");
    addOption(parser, ArgParseOption("pchg",  "prob-chg", "The probability of an 'C' beeing methylated in CHG-context.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "pchg", "0.06");
    addOption(parser, ArgParseOption("pchh",  "prob-chh", "The probability of an 'C' beeing methylation in CHH-context.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "pchh", "0.01");

    addOption(parser, ArgParseOption("p",  "sequence-protocol", "The sequence protocol used (Lister or Cokus).  Default: Cokus.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "p", "0");
    addOption(parser, ArgParseOption("c", "conversion-rate", "The C-T conversion rate.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "c", "0.98");
    addOption(parser, ArgParseOption("pmi", "prob-meth-ins", "The probability that a methylation is inserted into haplotype C, which was not methylated before (not if -umr is set).", ArgParseOption::DOUBLE));
    addOption(parser, ArgParseOption("pmd", "prob-meth-del", "The probability that a methylation is deleted at haplotye C (not if -umr is set).", ArgParseOption::DOUBLE)); 

    addOption(parser, ArgParseOption("umr", "use-methyl-rates", "Use given methylation rates.  Default: false."));
    setDefaultValue(parser, "umr", "false");
    addOption(parser, ArgParseOption("mr",  "methyl-rates-file", "Input file with given methylation rates.  Default: \"\".", ArgParseOption::STRING));

}

ArgumentParser::ParseResult
parseArgumentsAndCheckModelSpecific(Options<IlluminaReadsBS> & options,
                                          ArgumentParser & parser)
{
    Options<IlluminaReads> optionsBase;
    extractOptionsBase(optionsBase, options);  

    ArgumentParser::ParseResult res = parseArgumentsAndCheckModelSpecific(optionsBase, parser);
    mergeBaseIntoBSOptions(options, optionsBase);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Lazy....
    //options.haplotypeNoN = true;
    options.simulateQualities = true;
   // options.numHaplotypes = 2;
    options.readNamePrefix = "read";
    //options.readLength = 60;

    //options.haplotypeIndelRate = 0.0;

    getOptionValue(options.probabilityMethylationCG, parser, "prob-cg");
    getOptionValue(options.probabilityMethylationCHG, parser, "prob-chg");
    getOptionValue(options.probabilityMethylationCHH, parser, "prob-chh");
    
    getOptionValue(options.sequencingProtocol, parser, "sequence-protocol");
    getOptionValue(options.conversionRate, parser, "conversion-rate");
    getOptionValue(options.haplotypeMethylInsRate, parser, "prob-meth-ins");
    getOptionValue(options.conversionRate, parser, "conversion-rate");
    getOptionValue(options.haplotypeMethylDelRate, parser, "prob-meth-del");


    if (isSet(parser, "use-methyl-rates"))
    {
        options.useMethylRates = true;
        getOptionValue(options.methylRatesFile, parser, "methyl-rates-file");
    }
    else
    {
        options.methOutputFile = options.outputFile;
        append(options.methOutputFile, ".meths");
    }

    if (options.writeVCFFile)
    {
        options.vcfOutputFile = options.outputFile;
        append(options.vcfOutputFile, ".vcf");
    }
    /*
    if (options.writeMethFile)
    {
        options.methOutputFile = options.outputFile;
        append(options.methOutputFile, ".meths");
    }
    */

    return ArgumentParser::PARSE_OK;
}

// Load given methylation rates from a file.
int loadMethylRates(ModelParameters<IlluminaReadsBS> & parameters, 
                         Options<IlluminaReadsBS> const & options, 
                         FragmentStore<MyFragmentStoreConfig> & fragmentStore)
{
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;

    std::fstream inputFile(toCString(options.methylRatesFile), std::ios::binary | std::ios::in);
    TRecordReader reader(inputFile);

    // Resize methylation rate strings corresponding contig lengths   
    clear(parameters.topMethylRates);
    clear(parameters.bottomMethylRates);
    resize(parameters.topMethylRates, length(fragmentStore.contigStore), Exact());
    resize(parameters.bottomMethylRates, length(fragmentStore.contigStore), Exact());
    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i){
        resize(value(parameters.topMethylRates, i), length(fragmentStore.contigStore[i].seq), 0.0, Exact());      
        resize(value(parameters.bottomMethylRates, i), length(fragmentStore.contigStore[i].seq), 0.0, Exact());     
    }

    CharString contigName;
    CharString contigName_Prev = "";
    CharString position;
    CharString strand;
    CharString methylationRate;

    size_t contigId = 0;
    refresh(fragmentStore.contigNameStoreCache);
    while(!atEnd(reader))
    {
        clear(contigName);
        clear(position);
        clear(strand);
        clear(methylationRate);
        // Read contigName
        readUntilWhitespace(contigName, reader);
        skipWhitespaces(reader);
        // Read position
        readUntilWhitespace(position, reader);
        skipWhitespaces(reader);
        // Read strand
        readUntilWhitespace(strand, reader);
        skipWhitespaces(reader);
        // Read methylationRate
        readUntilWhitespace(methylationRate, reader);
        skipLine(reader);

        // If 1st occurence of contigName:
        if (contigName != contigName_Prev){
            assign(contigName_Prev, contigName);
            contigId = 0;
            // Get contigId for current contigName out of contigNameStore 
            if (!getIdByName(fragmentStore.contigNameStore, contigName, contigId, fragmentStore.contigNameStoreCache)) {
                        std::cerr << "ERROR: Could not find contig with name \"" << contigName << "\" (from methylation rates file) in contigs." << std::endl;
                        return 1;
            }         
        }
        if (lexicalCast<size_t>(position) >= length(fragmentStore.contigStore[contigId].seq))
        {
            std::cerr << "ERROR: Positions in file with methylation rate do not fit to length of reference sequence! (Reference sequence has length " << length(fragmentStore.contigStore[contigId].seq) << ", position with methylation rate was " << position << ")" << std::endl;
            return 1;
        }
        // Assign methylation rate to top/bottomMethylRate strings
        if (strand == '+'){
            SEQAN_ASSERT_EQ(fragmentStore.contigStore[contigId].seq[lexicalCast<size_t>(position)], 'C'); 
            assignValue(value(parameters.topMethylRates, contigId), lexicalCast<size_t>(position), lexicalCast<double>(methylationRate));   
        } else {
            SEQAN_ASSERT_EQ(fragmentStore.contigStore[contigId].seq[lexicalCast<size_t>(position)], 'G'); 
            assignValue(value(parameters.bottomMethylRates, contigId), lexicalCast<size_t>(position), lexicalCast<double>(methylationRate)); 
        }
    }
    return 0;
}


template <typename TMethylRates, typename TContigId, typename TPos, typename TRNG, typename TContig, typename TParameters>
inline void
assignMethLevelTop(TMethylRates &topMethylRates, 
                 TContigId &contigId, 
                 TPos &pos, 
                 TContig &contig, 
                 TRNG & rng,
                 TParameters &parameters,
                 Options<IlluminaReadsBS> const & options)
{
    SEQAN_ASSERT_LEQ(pos, length(topMethylRates[contigId]));
    SEQAN_ASSERT_LEQ(pos, length(contig));

    double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
    double y;
    if (pos < length(contig)-1 && contig[pos+1] == 'G')         // CG context
    {
        // Decide if C is 'methylated' or 'unmethylated'
        if (x < options.probabilityMethylationCG)                                           // methylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(parameters.thresholdYCG, 1));
        else                                                                                // unmethylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(0, parameters.thresholdYCG));
        // Pick methylation level corresponding to given distribution for choosen 'methylation state'
        topMethylRates[contigId][pos] = boost::math::ibeta_inv(options.pCG, options.qCG, y);
    }
    else if (pos < length(contig)-2 && contig[pos+2] == 'G')    // CHG context
    {
        // Decide if C is 'methylated' or 'unmethylated'
        if (x < options.probabilityMethylationCHG)                                           // methylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(parameters.thresholdYCHG, 1));
        else                                                                                // unmethylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(0, parameters.thresholdYCHG));
        // Pick methylation level corresponding to given distribution for choosen 'methylation state'
        topMethylRates[contigId][pos] = boost::math::ibeta_inv(options.pCHG, options.qCHG, y);
    }
    else                                                    // CHH context
    {
        // Decide if C is 'methylated' or 'unmethylated'
        if (x < options.probabilityMethylationCHH)                                           // methylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(parameters.thresholdYCHH, 1));
        else                                                                                // unmethylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(0, parameters.thresholdYCHH));
        // Pick methylation level corresponding to given distribution for choosen 'methylation state'
        topMethylRates[contigId][pos] = boost::math::ibeta_inv(options.pCHH, options.qCHH, y);
    }
}

template <typename TMethylRates, typename TContigId, typename TPos, typename TRNG, typename TContig, typename TParameters>
inline void
assignMethLevelBottom(TMethylRates & bottomMethylRates, 
                 TContigId &contigId, 
                 TPos &pos, 
                 TContig &contig, 
                 TRNG & rng, 
                 TParameters &parameters,
                 Options<IlluminaReadsBS> const & options)
{
    SEQAN_ASSERT_LEQ(pos, length(bottomMethylRates[contigId]));
    SEQAN_ASSERT_LEQ(pos, length(contig));

    double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
    double y;
    if (pos >= 1 && contig[pos+1] == 'C')         // CG context
    {
        // Decide if C is 'methylated' or 'unmethylated'
        if (x < options.probabilityMethylationCG)                                           // methylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(parameters.thresholdYCG, 1));
        else                                                                                // unmethylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(0, parameters.thresholdYCG));
        // Pick methylation level corresponding to given distribution for choosen 'methylation state'
        bottomMethylRates[contigId][pos] = boost::math::ibeta_inv(options.pCG, options.qCG, y);
    }
    else if (pos >= 2 && contig[pos+2] == 'C')    // CHG context
    {
        // Decide if C is 'methylated' or 'unmethylated'
        if (x < options.probabilityMethylationCHG)                                           // methylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(parameters.thresholdYCHG, 1));
        else                                                                                // unmethylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(0, parameters.thresholdYCHG));
        // Pick methylation level corresponding to given distribution for choosen 'methylation state'
        bottomMethylRates[contigId][pos] = boost::math::ibeta_inv(options.pCHG, options.qCHG, y);
    }
    else                                        // CHH context
    {
        // Decide if C is 'methylated' or 'unmethylated'
        if (x < options.probabilityMethylationCHH)                                           // methylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(parameters.thresholdYCHH, 1));
        else                                                                                // unmethylated
            y = pickRandomNumber(rng, Pdf<Uniform<double> >(0, parameters.thresholdYCHH));
        // Pick methylation level corresponding to given distribution for choosen 'methylation state'
        bottomMethylRates[contigId][pos] = boost::math::ibeta_inv(options.pCHH, options.qCHH, y);
    }
}


// Simulate methylation rates
int simulateMethylRates(ModelParameters<IlluminaReadsBS> & parameters, 
                         Options<IlluminaReadsBS> const & options, 
                         FragmentStore<MyFragmentStoreConfig> & fragmentStore)
{
    Rng<MersenneTwister> rng(options.seed);

    // Resize methylation rate strings corresponding contig lengths   
    clear(parameters.topMethylRates);
    clear(parameters.bottomMethylRates);
    resize(parameters.topMethylRates, length(fragmentStore.contigStore), Exact());
    resize(parameters.bottomMethylRates, length(fragmentStore.contigStore), Exact());

    // For each contig
    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i){
        resize(value(parameters.topMethylRates, i), length(fragmentStore.contigStore[i].seq), 0.0, Exact());      
        resize(value(parameters.bottomMethylRates, i), length(fragmentStore.contigStore[i].seq), 0.0, Exact());    

        for (unsigned j = 0; j < length(fragmentStore.contigStore[i].seq); ++j)
        {
            if (fragmentStore.contigStore[i].seq[j] == 'C') assignMethLevelTop(parameters.topMethylRates, i, j, fragmentStore.contigStore[i].seq, rng, parameters, options);
            else if (fragmentStore.contigStore[i].seq[j] == 'G') assignMethLevelBottom(parameters.bottomMethylRates, i, j, fragmentStore.contigStore[i].seq, rng, parameters, options);
        }
    } 
    return 0;
}


// Write simulated methylation level
// For reference, not adjusted for different haplotypes!!
// something goes wrong for output if only one haplotype 
template<typename TOptions>
bool writeSimMethLevel(FragmentStore<MyFragmentStoreConfig> &fragmentStore, ModelParameters<IlluminaReadsBS> & parameters, TOptions const & options)
{
    std::fstream file(toCString(options.methOutputFile), std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Could not open meth output file \"" << options.methOutputFile << "\"" << std::endl;
        return 1;
    }
    file << "#CHROM\tPOS\tSTRAND\tMETHLEVEL\n"; 

    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore), length(parameters.topMethylRates));
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore), length(parameters.bottomMethylRates));

    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
    {
        SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[i].seq), length(parameters.topMethylRates[i]));
        SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[i].seq), length(parameters.bottomMethylRates[i]));

        for (unsigned j = 0; j < length(fragmentStore.contigStore[i].seq); ++j)
        {
            if (fragmentStore.contigStore[i].seq[j] == 'C')
            {
                streamPut(file, fragmentStore.contigNameStore[i]);
                streamPut(file, '\t');
                streamPut(file, j);
                streamPut(file, '\t');
                streamPut(file, '+');
                streamPut(file, '\t');
                streamPut(file, parameters.topMethylRates[i][j]);
                streamPut(file, '\n');
            }
            else if (fragmentStore.contigStore[i].seq[j] == 'G')
            {
                streamPut(file, fragmentStore.contigNameStore[i]);
                streamPut(file, '\t');
                streamPut(file, j);
                streamPut(file, '\t');
                streamPut(file, '-');
                streamPut(file, '\t');
                streamPut(file, parameters.bottomMethylRates[i][j]);
                streamPut(file, '\n');
            }
        }
    }
    return 0;
}


// Calls original functions for non-BS tags
template<typename TModelParameters, typename TOptions>
int simulateReadsSetupModelSpecificData(TModelParameters & parameters,
                                        TOptions const & options,
                                        FragmentStore<MyFragmentStoreConfig> & /*NOP*/)
{
    int ret = simulateReadsSetupModelSpecificData(parameters, options);
    if (ret != 0)
        return ret;

    return 0;
}

// Bs_change: need fragmentStore for contig size -> size of methylationRates strings
int simulateReadsSetupModelSpecificData(ModelParameters<IlluminaReadsBS> & parameters,
                                        Options<IlluminaReadsBS> const & options,
                                        FragmentStore<MyFragmentStoreConfig> & fragmentStore)
{
    // Build standard illumina parameters:
    ModelParameters<IlluminaReads> parametersBase;
    // Extract IlluminaReads parameters out of IlluminaReadsBS parameters 
    extractModelParametersBase(parametersBase, parameters); 
    int ret = simulateReadsSetupModelSpecificData(parametersBase, static_cast<Options<IlluminaReads> >(options) );
    // Merge IlluminaReads parameters again into IlluminaReadsBS parameters
    mergeBaseIntoBSModelParameters(parameters, parametersBase);
    if (ret != 0)
        return ret;
   
    // Read methylation rates if given
    if (options.useMethylRates){
        std::cerr << "Loading methylation rates from \"" << options.methylRatesFile << "\"" << std::endl;
        ret = loadMethylRates(parameters, options, fragmentStore);
        if (ret != 0)
            return ret;
    }
    else
    {
        // Get y values for meth level thresholds...
        parameters.thresholdYCG = boost::math::ibeta(options.pCG, options.qCG, options.thresholdMethLevel);
        parameters.thresholdYCHG = boost::math::ibeta(options.pCHG, options.qCHG, options.thresholdMethLevel);
        parameters.thresholdYCHH = boost::math::ibeta(options.pCHH, options.qCHH, options.thresholdMethLevel);

        std::cout << "Simulate methylation rates." << std::endl;
        simulateMethylRates(parameters, options, fragmentStore);
        writeSimMethLevel(fragmentStore, parameters, options);
    }

    return 0;
}


template <typename TRNG>
unsigned pickReadLength(TRNG const &, Options<IlluminaReadsBS> const & options)
{
    return options.readLength;
}


// Build a haplotype, based on the contigs from the given fragment store.
// And adjust methylation patterns 
template <typename TRNG>
void buildHaplotype(StringSet<String<Dna5, Journaled<Alloc<> > > > & haplotype,
                    ModelParameters<IlluminaReadsBS> & parameters,
                    FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                    String<String<Snp> > & snpSet,
                    TRNG & rng,
                    Options<IlluminaReadsBS> const & options) {

    resize(haplotype, length(fragmentStore.contigStore), Exact());
    String<Dna5> buffer;
    reserve(buffer, options.haplotypeIndelRangeMax);
    Snp snp;
  
    // Bs_change:
    // Journaled strings: store only differences in haplotypes
    clear(parameters.topMethylRatesHaplotype);
    clear(parameters.bottomMethylRatesHaplotype);
    resize(parameters.topMethylRatesHaplotype, length(parameters.topMethylRates));
    resize(parameters.bottomMethylRatesHaplotype, length(parameters.bottomMethylRates));

    for (unsigned i = 0; i < length(parameters.topMethylRates); ++i){
        setHost(parameters.topMethylRatesHaplotype[i], parameters.topMethylRates[i]);
        setHost(parameters.bottomMethylRatesHaplotype[i], parameters.bottomMethylRates[i]);
    }

    String<double> bs_bufferTop;
    String<double> bs_bufferBottom;

    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
    {
        // statistics
        unsigned numSNPs = 0;
        unsigned numIndels = 0;
        unsigned indelLenSum = 0;

        std::cout << "    contig # " << i+1 << "/" << length(fragmentStore.contigStore) << '\t' << fragmentStore.contigNameStore[i] << std::flush;

        clear(haplotype[i]);
        setHost(haplotype[i], fragmentStore.contigStore[i].seq);
        String<Dna5> const & contig = fragmentStore.contigStore[i].seq;
        String<Dna5, Journaled<Alloc<> > > & haplotypeContig = haplotype[i];

        // Only generate Ns in the haplotype if allowed.
        int maxOrdValue = options.haplotypeNoN ? 3 : 4;

        // j is position in original sequence, k is position in haplotype
        for (size_t j = 0, k = 0; j < length(contig);)
        {
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            snp.virtualPos = k;
            if (x < options.haplotypeSnpRate) {
                // SNP
                ++numSNPs;
                snp.length = 1;

                Dna5 c = Dna5(pickRandomNumber(rng, Pdf<Uniform<int> >(0, maxOrdValue - 1)));
                if (c == contig[j])
                    c = Dna5(ordValue(c) + 1);
                if (options.haplotypeNoN)
                    SEQAN_ASSERT(c != Dna5('N'));
                assignValue(haplotypeContig, k, c);
                snp.altBase = c;
                // Bs_change: adjust methylation state if C is changed to different base
                if (contig[j] == 'C')
                    parameters.topMethylRatesHaplotype[i][k] = 0.0;
                else if (contig[j] == 'G')
                    parameters.bottomMethylRatesHaplotype[i][k] = 0.0;
                // If conversion into 'C': maybe methylate it, depending on probabilities
                // Assign meth level dependend on context and prob. distribution
                // ATTENTIONE: Dirty, using context from current haplotype sequence for meth. level generation (context behind this position could still change, but chance is anyway low)
                if (c == 'C') assignMethLevelTop(parameters.topMethylRatesHaplotype, i, k, haplotypeContig, rng, parameters, options);
                else if (c == 'G') assignMethLevelBottom(parameters.bottomMethylRatesHaplotype, i, k, haplotypeContig, rng, parameters, options);

                ++j;
                ++k;
                // HAPTYPE  ----GATTACA---- 
                // REF      ----CACACAC----
                snp.type = ERROR_TYPE_MISMATCH;
                appendValue(snpSet[i], snp);
            }
            else if (x < options.haplotypeSnpRate + options.haplotypeIndelRate) 
            {
                // Indel of random length.
                ++numIndels;
                unsigned rangeLen = options.haplotypeIndelRangeMax - options.haplotypeIndelRangeMin;
                unsigned indelLen = options.haplotypeIndelRangeMin + static_cast<unsigned>(pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1)) * rangeLen);
                snp.length = indelLen;
                indelLenSum += indelLen;
                if (pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1)) < 0.5) 
                {
                    // Insertion.
                    clear(buffer);
                    for (unsigned ii = 0; ii < indelLen; ++ii)
                        appendValue(buffer, Dna5(pickRandomNumber(rng, Pdf<Uniform<int> >(0, maxOrdValue))));
                    insert(haplotypeContig, k, buffer);
                    // Bs_change: adjust length of methylation string
                    clear(bs_bufferTop);
                    clear(bs_bufferBottom);
                    resize(bs_bufferTop, indelLen, 0.0, Exact());
                    resize(bs_bufferBottom, indelLen, 0.0, Exact());
                    // Insert methylations at new indel Cs dependent on context and prob. distribution
                    // ATTENTIONE: Dirty, using context from current haplotype sequence for meth. level generation (context behind this position could still change, but chance is anyway low)
                    insert(value(parameters.topMethylRatesHaplotype, i), k, bs_bufferTop);
                    insert(value(parameters.bottomMethylRatesHaplotype, i), k, bs_bufferBottom);
                    for (unsigned ii = 0; ii < indelLen; ++ii)
                    {
                        size_t vPos = k+ii;
                        if (buffer[ii] == 'C') assignMethLevelTop(parameters.topMethylRatesHaplotype, i, vPos, haplotypeContig, rng, parameters, options);
                        else if (buffer[ii] == 'G') assignMethLevelBottom(parameters.bottomMethylRatesHaplotype, i, vPos, haplotypeContig, rng, parameters, options);
                    }

                    k += indelLen;

                    // HAPTYPE  ----GATTACA---- 
                    // REF      ----       ----
                    snp.type = ERROR_TYPE_INSERT;
                }
                else
                {
                    // Deletion.
                    indelLen = _min(indelLen, length(haplotypeContig) - k);
                    erase(haplotypeContig, k, k + indelLen);
                    // Bs_change: adjust length of methylation string
                    erase(value(parameters.topMethylRatesHaplotype, i), k, k + indelLen);
                    erase(value(parameters.bottomMethylRatesHaplotype, i), k, k + indelLen); 

                    j += indelLen;

                    // HAPTYPE  ----       ----
                    // REF      ----GATTACA---- 
                    snp.type = ERROR_TYPE_DELETE;
                }
                appendValue(snpSet[i], snp);
            }
            else 
            {
                // Match.
                j += 1;
                k += 1;
            }
        }
        std::cout << "\tSNPs:" << numSNPs << "\tindels:" << numIndels << "\tindel len sum:\t" << indelLenSum << "\tvrate:" << (numSNPs+indelLenSum)/(double)length(contig) << std::endl;
      
        SEQAN_ASSERT_EQ(length(haplotype[i]), length(parameters.topMethylRatesHaplotype[i]));
    }

}

template <typename TRNG, typename TContig>
void buildSimulationInstructions(ReadSimulationInstruction<IlluminaReadsBS> & inst, TRNG & rng, unsigned readLength, TContig const & contig, ModelParameters<IlluminaReadsBS> const & parameters, Options<IlluminaReadsBS> const & options) {
    clear(inst.editString);
    reserve(inst.editString, static_cast<size_t>(1.2 * readLength), Generous());
    inst.delCount = 0;
    inst.insCount = 0;

    //
    // Build Edit String.    //
    for (unsigned i = 0; i < readLength; /*NOP*/) {
        double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
        double pMismatch = parameters.mismatchProbabilities[i];
        double pInsert   = options.probabilityInsert;
        double pDelete   = options.probabilityDelete;
        double pMatch    = 1.0 - pMismatch - pInsert - pDelete;
        if (x < pMatch) {
            // match
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MATCH);
        } else if (x < pMatch + pMismatch) {
            // mismatch
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MISMATCH);
        } else if (x < pMatch + pMismatch + pInsert) {
            // insert
            if (length(inst.editString) > 0 && back(inst.editString == ERROR_TYPE_DELETE)) {
                inst.delCount -= 1;
                eraseBack(inst.editString);
            } else {
                i += 1;
                inst.insCount += 1;
                appendValue(inst.editString, ERROR_TYPE_INSERT);
            }
        } else {
            // Decrement string size, do not add a delete if string is
            // too short, possibly remove insert from edit string.
            if (length(inst.editString) > 0) {
                if (back(inst.editString == ERROR_TYPE_INSERT)) {
                    i -= 1;
                    inst.insCount -= 1;
                    eraseBack(inst.editString);
                } else {
                    inst.delCount += 1;
                    appendValue(inst.editString, ERROR_TYPE_DELETE);
                }
            }
        }
    }
    SEQAN_ASSERT_EQ(readLength, length(inst.editString) - inst.delCount);

    //
    // Adjust Positions.
    //

    // If the number of deletions does not equal the number of inserts
    // then we have to adjust the read positions.
    if (inst.delCount != inst.insCount) {
        int delta = static_cast<int>(inst.delCount) - static_cast<int>(inst.insCount);
        inst.endPos += delta;
        if (inst.endPos > length(contig)) {
            delta = inst.endPos - length(contig);
            inst.endPos -= delta;
            inst.beginPos -= delta;
        }
        SEQAN_ASSERT_EQ(inst.endPos - inst.beginPos + inst.insCount - inst.delCount,
                        readLength);
    }

    //
    // Build bsConversion String (regarding adjusted positions in forward haplotype contig sequence)
    // Apart of begin and end positions, this is independently of the edit string, 
    // because later this bs conversions will be applied before the edit string 
    //
    clear(inst.bsConversionString);
    resize(inst.bsConversionString, inst.endPos - inst.beginPos, false, Exact());
    if (inst.isFromTopStrand) {
        for (unsigned i = 0; i < (inst.endPos - inst.beginPos); ++i) {
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            if (x > getValue(parameters.topMethylRatesHaplotype[inst.contigId], inst.beginPos + i)) {  
                // Only methylated positions are protected from conversion
                // remember later to check if position is C or G
                x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                if (x < options.conversionRate) {
                    // Bs conversion
                    assignValue(inst.bsConversionString, i, true);
                }
            } else
                SEQAN_ASSERT_EQ(inst.bsConversionString[i], false);                
        }
    } else {
        for (unsigned i = 0; i < (inst.endPos - inst.beginPos); ++i) {
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            if (x > getValue(parameters.bottomMethylRatesHaplotype[inst.contigId], inst.beginPos + i)) {
                x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                if (x < options.conversionRate) {
                    // Bs conversion
                    assignValue(inst.bsConversionString, i, true);
                }
            } else
                SEQAN_ASSERT_EQ(inst.bsConversionString[i], false);
        }
    }
    SEQAN_ASSERT_EQ((inst.endPos - inst.beginPos), length(inst.bsConversionString));

    //
    // Quality Simulation.
    //
    SEQAN_ASSERT_GT(length(inst.editString), 0u);
    if (options.simulateQualities) {
        clear(inst.qualities);
        resize(inst.qualities, length(inst.editString), 0, Exact());

        for (unsigned i = 0, j = 0; i < length(inst.editString); i++) {
            SEQAN_ASSERT_LEQ(j, inst.endPos - inst.beginPos + inst.delCount);
            if (inst.editString[i] == ERROR_TYPE_MISMATCH) {
                // std::cout << "i == " << i << ", j == " << j << ", parameters.mismatchQualityMeans[j] == " << parameters.mismatchQualityMeans[j] << ", parameters.mismatchQualityStdDevs[j] == " << parameters.mismatchQualityStdDevs[j] << std::endl;
                Pdf<Normal> pdf(parameters.mismatchQualityMeans[j], parameters.mismatchQualityStdDevs[j]);
                inst.qualities[i] = static_cast<int>(pickRandomNumber(rng, pdf));
            } else {
                Pdf<Normal> pdf(parameters.qualityMeans[j], parameters.qualityStdDevs[j]);
                inst.qualities[i] = static_cast<int>(pickRandomNumber(rng, pdf));
            }

            if (inst.editString[i] == ERROR_TYPE_MISMATCH || inst.editString[i] == ERROR_TYPE_MATCH)
                j += 1;
        }
    }

}

template<typename TRNG, typename TReadSimulationInstruction>
void pickOriginalStrand(TRNG & /*NOP*/, TReadSimulationInstruction & /*NOP*/){}

template<typename TRNG>
void pickOriginalStrand(TRNG & rng, ReadSimulationInstruction<IlluminaReadsBS> & inst){
    inst.isFromTopStrand = pickRandomNumber(rng, Pdf<Uniform<int> >(0, 1));
}

// Build read simulation instructions for a haplotype.
//
// pick a contig, probability is proportional to the length
// pick a start position, end position = start position + read length
// pick whether to match on the forward or reverse strand
// simulate edit string
// build quality values
// possibly adjust mate if left read has insert at the beginning or right read has insert at the right
//
// Note that we set the isForward flag of the simulation instructions here.  Notably, the technology specific simulation
// parts do not interpret this flag but only provide simulation instructions on the forward strand.  Later, when actually
// cutting out the sequences and modifying the sampled read, we have to take this in mind.  This means, there we have to
// reverse the edit string etc. there.
template <typename TRNG, typename THaplotypeSpec>
int buildReadSimulationInstruction(
        String<ReadSimulationInstruction<IlluminaReadsBS> > & instructions,
        TRNG & rng,
        unsigned const & haplotypeId,
        StringSet<String<Dna5, THaplotypeSpec> > const & haplotype,
        String<double> const & relativeContigLengths,
        size_t const & contigId,
        bool fixedContigId,
        ModelParameters<IlluminaReadsBS> const & parameters,
        Options<IlluminaReadsBS> const & options)
{

    ReadSimulationInstruction<IlluminaReadsBS> inst;
    inst.haplotype = haplotypeId;

    // We have to retry simulation if the mate pair did not fit in.
    bool invalid = false;
    do {
        clear(instructions);
        invalid = false;  // By default, we do not want to repeat.
        if (fixedContigId) {
            // Use precomputed contig id.
            inst.contigId = contigId;
        } else {
            // Pick contig id, probability is proportional to the length.
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            for (unsigned i = 0; i < length(relativeContigLengths); ++i) {
                if (x < relativeContigLengths[i]) {
                    inst.contigId = i - 1;
                    break;
                }
            }
        }

        // bs_change:
        // does nothing for non-bs instructions
        // is for matepairs the same
        pickOriginalStrand(rng, inst);

        // Pick whether on forward or reverse strand.
        if (options.sequencingProtocol == 0)    // Lister, directional
        {
            if (inst.isFromTopStrand)
                inst.isForward  = true;
            else
                inst.isForward = false;
        }
        else
            inst.isForward = pickRandomNumber(rng, Pdf<Uniform<int> >(0, 1));



        // Pick the length in the haplotype infix the read comes from, possibly randomly.
        unsigned readLength = pickReadLength(rng, options);
        // This cannot work if the haplotype is shorter than the length of the read to simulate.
        if (length(haplotype[inst.contigId]) < readLength) {
            std::cerr << "ERROR: haplotype (== " << length(haplotype[inst.contigId]) << ") < read length!" << std::endl;
            return 1;
        }
        // Pick a start and end position.
        inst.beginPos = pickRandomNumber(rng, Pdf<Uniform<size_t> >(0, length(haplotype[inst.contigId]) - readLength - 1));
        inst.endPos = inst.beginPos + readLength;
        // Simulate the read with these parameters.
        buildSimulationInstructions(inst, rng, readLength, haplotype[inst.contigId], parameters, options);
        // Append read to result list.
        appendValue(instructions, inst);
        //std::cout << " 2 beginPoss: " << inst.beginPos << std::endl;


        ReadSimulationInstruction<IlluminaReadsBS> inst2(inst);

        // Maybe create a mate for this read.
        if (options.generateMatePairs) {
            // Pick a read length, possibly randomly.
            unsigned readLength = pickReadLength(rng, options);
            // Pick a library length, according to the options.
            size_t libraryLength = pickLibraryLength(rng, options);
            // Compute start and end position.
            if (inst.isForward)
            {
                inst.endPos = inst.beginPos + libraryLength;
                inst.beginPos = inst.endPos - readLength;
            }
            else
            {
                if (inst.endPos < libraryLength)
                    invalid = true;
                inst.beginPos = inst.endPos - libraryLength;
                inst.endPos = inst.beginPos + readLength;
            }
            // Set orientation of second mate.
            inst.isForward = !back(instructions).isForward;
            // Verify that the mate fits right of the originally simulated read.
            size_t contigLength = length(haplotype[inst.contigId]);
            if ((inst.beginPos > contigLength) || (inst.endPos > contigLength)) {
                // Mate did not fit!  Remove previously added read and set
                // invalid to true so we repeat this simulation.
                SEQAN_ASSERT_GT(length(instructions), 0u);
                eraseBack(instructions);
                invalid = true;

                if (options.verbose) {
                    std::cerr << "INFO: Mate did not fit! Repeating..." << std::endl;
                    std::cerr << "      inst2.beginPos == " << inst2.beginPos << ", inst2.endPos == " << inst2.endPos << std::endl;
                    std::cerr << "       inst.beginPos == " << inst.beginPos << ",  inst.endPos == " << inst.endPos << std::endl;
                }
                continue;
            }
            // Simulate the read with these parameters.
            buildSimulationInstructions(inst, rng, readLength, haplotype[inst.contigId], parameters, options);
            // Append read to result list.
            appendValue(instructions, inst);
        }

        // Check whether there are Ns in the selected areas.
        if (!options.allowNFromGenome) {
            for (unsigned i = 0; i < length(instructions); ++i) {
                String<Dna5, THaplotypeSpec> const & contig = haplotype[instructions[i].contigId];
                typedef typename Position<Dna5String>::Type TPosition;
                TPosition beginPos = instructions[i].beginPos;
                TPosition endPos = instructions[i].endPos;
                if (beginPos > endPos)
                    std::swap(beginPos, endPos);
                for (unsigned i = beginPos; i != endPos; ++i) {
                    if (contig[i] == Dna5('N')) {
                        invalid = true;
                        break;
                    }
                }
            }
        }
    } while (invalid);
	
	if (options.generateMatePairs)
		SEQAN_ASSERT_EQ(length(instructions), 2u);
	else
		SEQAN_ASSERT_EQ(length(instructions), 1u);
	
    return 0;
}



template <typename TRNG, typename TString>
void applySimulationInstructions(TString & read, TRNG & rng, ReadSimulationInstruction<IlluminaReadsBS> & inst, Options<IlluminaReadsBS> const & options)
{
    typedef typename Value<TString>::Type TAlphabet;

    if (options.simulateQualities)
        SEQAN_ASSERT_EQ(length(inst.qualities), length(inst.editString));
    
    //
    // apply bs conversions
    //    

    // Top strand, forward strand sequence of read is edited
    // For bottom strand and/or reverse strand: modifying will happen later
    if (inst.isFromTopStrand){
        for (unsigned i = 0; i < length(inst.bsConversionString); ++i){
            if (inst.bsConversionString[i] == true && read[i] == 'C'){
                assignValue(read, i, 'T');
            } else {
                // Modify bsConversionString on the fly to a indicator string, if bs conversion was applied or not
                // TRUE only if current base is C and this C is converted
                assignValue(inst.bsConversionString, i, false);
            }
        }
    }else{
        for (unsigned i = 0; i < length(inst.bsConversionString); ++i){
            if (inst.bsConversionString[i] == true && read[i] == 'G'){
                assignValue(read, i, 'A');
            } else {
                assignValue(inst.bsConversionString, i, false);
            }
        }
    }
   
    //
    // apply errors
    //
    TString tmp;
    reserve(tmp, length(read) + inst.insCount - inst.delCount);
    unsigned j = 0;
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        SEQAN_ASSERT_LEQ(j, i);

        TAlphabet c;
        //int x, xold;
        switch (inst.editString[i]) {
            case ERROR_TYPE_MATCH:
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                appendValue(tmp, read[j]);
                if (options.simulateQualities)
                    assignQualityValue(back(tmp), inst.qualities[i]);
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " match" << std::endl;
                //std::cout << back(tmp) << " " << read[j] << " " << inst.qualities[i] << std::endl;
                j += 1;
                break;
            case ERROR_TYPE_MISMATCH:
                if (options.illuminaNoN) {
                    c = TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 3)));  // -3, N not allowed
                } else {
                    c = TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 2)));  // -2, N allowed
                }
                //xold = ordValue(c);
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                if (ordValue(c) >= ordValue(read[j]))
                    c = TAlphabet(ordValue(c) + 1);
                if (options.illuminaNoN)
                    SEQAN_ASSERT(c != TAlphabet('N'));
                //x = ordValue(c);
                appendValue(tmp, c);
                if (options.simulateQualities) {
                    if (options.illuminaNoN)  // Ns can be introduced through quality, too.
                        assignQualityValue(back(tmp), _max(1, inst.qualities[i]));
                    else
                        assignQualityValue(back(tmp), inst.qualities[i]);
                }
                // std::cout << i << " q(q_i)=" << getQualityValue(back(tmp)) << " q(i)=" << inst.qualities[i] << " char=" << convert<char>(back(tmp)) << " c_old=" << xold << " c=" << x << " r_j=" << ordValue(read[j]) << std::endl;
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " mismatch" << std::endl;
                //std::cout << "MM " << c << " " << back(tmp) << " " << inst.qualities[i] << std::endl;
                j += 1;
                break;
            case ERROR_TYPE_INSERT:
                if (options.illuminaNoN)
                    appendValue(tmp, TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 2))));  // -2 == no N
                else
                    appendValue(tmp, TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 1))));  // -1 == N allowed
                if (options.simulateQualities) {
                    if (options.illuminaNoN)  // Ns can be introduced through quality, too.
                        assignQualityValue(back(tmp), _max(1, inst.qualities[i]));
                    else
                        assignQualityValue(back(tmp), inst.qualities[i]);
                }
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " insertion" << std::endl;
                break;
            case ERROR_TYPE_DELETE:
                j += 1;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid error type.");
        }
    }
    SEQAN_ASSERT_EQ(j, length(read));
    SEQAN_ASSERT_GEQ(length(tmp), options.readLength);

    //std::cout << "tmp == " << tmp << std::endl;
    resize(tmp, options.readLength, Exact());
    move(read, tmp);
}

template<typename TStream, typename TContigName, typename TPos, typename TValue>
inline void
streamPutMethLevel(TStream &file, TContigName &contigName, TPos &pos, char const strand, TValue const methLevel)
{
    streamPut(file, contigName);
    streamPut(file, '\t');
    streamPut(file, pos);
    streamPut(file, '\t');
    streamPut(file, strand);
    streamPut(file, '\t');
    streamPut(file, methLevel);
    streamPut(file, '\n');
}

template<typename TStream, typename TContigName, typename TPos, typename TValue>
inline void
streamPutMethLevel(TStream &file, TContigName &contigName, TPos &pos, TValue const methLevelTop, TValue const methLevelBottom)
{
    streamPut(file, contigName);
    streamPut(file, '\t');
    streamPut(file, pos);
    streamPut(file, '\t');
    streamPut(file, "+:-");
    streamPut(file, '\t');
    streamPut(file, methLevelTop);
    streamPut(file, ':');
    streamPut(file, methLevelBottom);
    streamPut(file, '\n');
}


// Write vcf file for snps in simulated haplotypes
// here: only write out MISMATCHES, no insertions, deletions ..!!!
template <typename TAllHaploSnpSets, typename TMethylRates, typename TParameters, typename TOptions>
bool writeVCFOutput(FragmentStore<MyFragmentStoreConfig> &fragmentStore, 
                    TAllHaploSnpSets &allSnpSets, 
                    TMethylRates &topMethylRatesHaplotype0, 
                    TMethylRates &bottomMethylRatesHaplotype0,
                    TMethylRates &topMethylRatesHaplotype1,
                    TMethylRates &bottomMethylRatesHaplotype1,
                    TParameters &parameters,
                    TOptions const & options)
{
    std::fstream snpFile(toCString(options.vcfOutputFile), std::ios_base::out);
    if (!snpFile.is_open()) {
        std::cerr << "Could not open snp output VCF file \"" << options.vcfOutputFile << "\"" << std::endl;
        return 1;
    }
    snpFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n";
    std::fstream methFile(toCString(options.methOutputFile), std::ios_base::out);
    if (!methFile.is_open()) {
        std::cerr << "Could not open simulated meth output file \"" << options.methOutputFile << "\"" << std::endl;
        return 1;
    }
    methFile << "#CHROM\tPOS\tSTRAND\tMETHLEVEL\n"; 

    // Max no. of haplotypes: 2 !!!
    std::cout << "Write snp and meth level output" << std::endl;
       
    typedef String<Snp> TContigSnpSet;

    SEQAN_ASSERT_EQ(length(allSnpSets[0]), length(fragmentStore.contigStore));
    // For each contig
    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
    { 
        CharString &contigName = fragmentStore.contigNameStore[i];
        // For 2 haplotypes
        if (options.numHaplotypes == 2)
        {
            int diff0 = 0;    // diff. between virtual position and real contig position
            int diff1 = 0;    // TODO: not checked for simulation with indels...!
            unsigned p0 = 0;
            unsigned p1 = 0;
            unsigned virtualPos0;
            unsigned virtualPos1;
            unsigned methPos = 0;
            // TODO use iterators, code ugly
            // Iterate over snps in both haplotypes parallel
            while (p0 < length(allSnpSets[0][i]) || p1 < length(allSnpSets[1][i]))    
            {
                unsigned contigPos0;      
                if (p0 < length(allSnpSets[0][i]))
                {
                    //Snp &currSnp0 = allSnpSets[0][i][p0];
                    // Real contig position is calculated by virtual position + length of deletion - length of insertions up to this position
                    virtualPos0 = allSnpSets[0][i][p0].virtualPos;
                    contigPos0 = virtualPos0 + diff0;
                    if (allSnpSets[0][i][p0].type == ERROR_TYPE_INSERT)  // If insertion, subtract length, iterate to next snp entry and continue

                    {
                        diff0 -= allSnpSets[0][i][p0].length;
                        ++p0;
                        continue;
                    }
                    else if (allSnpSets[0][i][p0].type == ERROR_TYPE_DELETE)  // If deletion, add length, iterate to next snp entry and continue

                    {
                        diff0 += allSnpSets[0][i][p0].length;
                        ++p0;
                        continue;
                    }
                }
                unsigned contigPos1;
                if (p1 < length(allSnpSets[1][i]))
                {
                    //Snp &currSnp1 =  allSnpSets[1][i][p1];
                    virtualPos1 = allSnpSets[1][i][p1].virtualPos;
                    contigPos1 = virtualPos1 + diff1;
                    if (allSnpSets[1][i][p1].type == ERROR_TYPE_INSERT)
                    {
                        diff1 -= allSnpSets[1][i][p1].length;
                        ++p1;
                        continue;
                    }
                    else if (allSnpSets[1][i][p1].type == ERROR_TYPE_DELETE)
                    {
                        diff1 += allSnpSets[1][i][p1].length;
                        ++p1;
                        continue;
                    }
                }
                // Between snps: output meth level for C positions in original contig from original meth levels
                while (methPos < contigPos0 && methPos < contigPos1)
                {
                    if (fragmentStore.contigStore[i].seq[methPos] == 'C')
                    {
                        streamPutMethLevel(methFile, contigName, methPos, '+', parameters.topMethylRates[i][methPos]);
                    }
                    else if (fragmentStore.contigStore[i].seq[methPos] == 'G')
                    {
                        streamPutMethLevel(methFile, contigName, methPos, '-', parameters.bottomMethylRates[i][methPos]);
                    }
                    ++methPos;
                }

                if (p0 < length(allSnpSets[0][i]) && p1 < length(allSnpSets[1][i])  &&                  // not at end of snpSet
                    allSnpSets[0][i][p0].type == ERROR_TYPE_MISMATCH && allSnpSets[1][i][p1].type == ERROR_TYPE_MISMATCH &&     // if both MISMATCH
                   (contigPos0 == contigPos1) )                                                         // at the same contig position
                {
                    // Snp output
                    streamPut(snpFile, contigName);
                    streamPut(snpFile, '\t');
                    streamPut(snpFile, contigPos0);
                    streamPut(snpFile, '\t');
                    streamPut(snpFile, '.');
                    streamPut(snpFile, '\t');
                    // REF
                    Dna5 const &ref = fragmentStore.contigStore[i].seq[contigPos0];
                    streamPut(snpFile, ref);
                    streamPut(snpFile, '\t');
                    // ALT
                    if (allSnpSets[0][i][p0].altBase == allSnpSets[1][i][p1].altBase)       // homo snp
                    {
                        streamPut(snpFile, allSnpSets[0][i][p0].altBase);
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, '.');
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, "PASS");
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, "AF=100");
                        streamPut(snpFile, '\n');
                    }
                    else                                            // hetero snp
                    {
                        streamPut(snpFile, allSnpSets[0][i][p0].altBase);
                        streamPut(snpFile, ',');
                        streamPut(snpFile, allSnpSets[1][i][p1].altBase);
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, '.');
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, "PASS");
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, "AF=50,50");
                        streamPut(snpFile, '\n');
                    }
                    // Meth output
                    if (allSnpSets[0][i][p0].altBase == allSnpSets[1][i][p1].altBase && allSnpSets[0][i][p0].altBase == 'C') // Haplotype CC
                    {
                        double methRate = double(topMethylRatesHaplotype0[i][virtualPos0] + topMethylRatesHaplotype1[i][virtualPos1])/2.0;  // different simulated meth level are not possible anymore to distinguish -> build average
                        streamPutMethLevel(methFile, contigName, methPos, '+', methRate);
                    }
                    else if (allSnpSets[0][i][p0].altBase == allSnpSets[1][i][p1].altBase && allSnpSets[0][i][p0].altBase == 'G') // Haplotype GG
                    {
                        double methRate = double(bottomMethylRatesHaplotype0[i][virtualPos0] + bottomMethylRatesHaplotype1[i][virtualPos1])/2.0; 
                        streamPutMethLevel(methFile, contigName, methPos, '-', methRate);
                    }
                    else if ( allSnpSets[0][i][p0].altBase == 'C' && allSnpSets[1][i][p1].altBase == 'G') // Haplotype CG
                        streamPutMethLevel(methFile, contigName, methPos, (double)topMethylRatesHaplotype0[i][virtualPos0], (double)bottomMethylRatesHaplotype1[i][virtualPos1]);
                    else if ( allSnpSets[0][i][p0].altBase == 'G' && allSnpSets[1][i][p1].altBase == 'C') // Haplotype GC
                        streamPutMethLevel(methFile, contigName, methPos, (double)topMethylRatesHaplotype1[i][virtualPos1], (double)bottomMethylRatesHaplotype0[i][virtualPos0]);
                    else if (allSnpSets[0][i][p0].altBase == 'C') // CX
                        streamPutMethLevel(methFile, contigName, methPos, '+', topMethylRatesHaplotype0[i][virtualPos0]);
                    else if (allSnpSets[1][i][p1].altBase == 'C') // XC
                        streamPutMethLevel(methFile, contigName, methPos, '+', topMethylRatesHaplotype1[i][virtualPos1]);
                    else if (allSnpSets[0][i][p0].altBase == 'G') // GX
                        streamPutMethLevel(methFile, contigName, methPos, '-', bottomMethylRatesHaplotype0[i][virtualPos0]);
                    else if (allSnpSets[1][i][p1].altBase  == 'G') // XG
                        streamPutMethLevel(methFile, contigName, methPos, '-', bottomMethylRatesHaplotype1[i][virtualPos1]);
                   
                    ++p0;
                    ++p1;
                    ++methPos;
                }
                else
                {
                    if (((p1 == length(allSnpSets[1][i]))  ||
                         (p0 < length(allSnpSets[0][i]) && contigPos0 <= contigPos1) ) 
                       && allSnpSets[0][i][p0].type == ERROR_TYPE_MISMATCH)
                    {
                        // Snp output
                        streamPut(snpFile, contigName);
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, contigPos0);
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, '.');
                        streamPut(snpFile, '\t');
                        // REF
                        Dna5 const &ref = fragmentStore.contigStore[i].seq[contigPos0];
                        streamPut(snpFile, ref);
                        streamPut(snpFile, '\t');
                        // ALT
                        streamPut(snpFile, allSnpSets[0][i][p0].altBase);
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, '.');
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, "PASS");
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, "AF=50\tH=0");
                        streamPut(snpFile, '\n');

                        // Meth output
                        if (ref == 'C') SEQAN_ASSERT_NEQ(allSnpSets[0][i][p0].altBase, 'C');
                        else if (ref == 'G') SEQAN_ASSERT_NEQ(allSnpSets[0][i][p0].altBase, 'G');
                        if (allSnpSets[0][i][p0].altBase == 'G' && ref == 'C') // Haplotype GC
                            streamPutMethLevel(methFile, contigName, methPos, (double)parameters.topMethylRates[i][methPos], (double)bottomMethylRatesHaplotype0[i][virtualPos0]);
                        else if (ref == 'C') // Haplotype XC
                            streamPutMethLevel(methFile, contigName, methPos, '+', parameters.topMethylRates[i][methPos]);
                        else if (allSnpSets[0][i][p0].altBase == 'C' && ref == 'G') // Haplotype CG
                            streamPutMethLevel(methFile, contigName, methPos, (double)topMethylRatesHaplotype0[i][virtualPos0], (double)parameters.bottomMethylRates[i][methPos]);
                        else if (ref == 'G') // Haplotype XG
                            streamPutMethLevel(methFile, contigName, methPos, '-', parameters.bottomMethylRates[i][methPos]);
                        else if (allSnpSets[0][i][p0].altBase == 'C' && ref != 'C' && ref != 'G') // CX
                            streamPutMethLevel(methFile, contigName, methPos, '+', topMethylRatesHaplotype0[i][virtualPos0]);
                        else if (allSnpSets[0][i][p0].altBase == 'G' && ref != 'C' && ref != 'G') // GX
                            streamPutMethLevel(methFile, contigName, methPos, '-', bottomMethylRatesHaplotype0[i][virtualPos0]);
                       
                        ++p0; 
                       ++methPos;
                    }
                    else if (((p0 == length(allSnpSets[0][i])) ||  
                              (p1 < length(allSnpSets[1][i]) && contigPos1 <= contigPos0) ) 
                            && allSnpSets[1][i][p1].type == ERROR_TYPE_MISMATCH)
                    {
                        // Snp output
                        streamPut(snpFile, contigName);
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, contigPos1);
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, '.');
                        streamPut(snpFile, '\t');
                        // REF
                        Dna5 const &ref = fragmentStore.contigStore[i].seq[contigPos1];
                        streamPut(snpFile, ref);
                        streamPut(snpFile, '\t');
                        // ALT
                        streamPut(snpFile, allSnpSets[1][i][p1].altBase);
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, '.');
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, "PASS");
                        streamPut(snpFile, '\t');
                        streamPut(snpFile, "AF=50\tH=1");
                        streamPut(snpFile, '\n');

                        // Meth output
                        if (ref == 'C') SEQAN_ASSERT_NEQ(allSnpSets[1][i][p1].altBase, 'C');
                        else if (ref == 'G') SEQAN_ASSERT_NEQ(allSnpSets[1][i][p1].altBase, 'G');
                        if (ref == 'C' && allSnpSets[1][i][p1].altBase == 'G') // Haplotype CG
                            streamPutMethLevel(methFile, contigName, methPos, (double)parameters.topMethylRates[i][methPos], (double)bottomMethylRatesHaplotype0[i][virtualPos1]);
                        else if (ref == 'C') // Haplotype CX
                            streamPutMethLevel(methFile, contigName, methPos, '+', parameters.topMethylRates[i][methPos]);
                        else if (ref == 'G' && allSnpSets[1][i][p1].altBase == 'C') // Haplotype GC
                            streamPutMethLevel(methFile, contigName, methPos, (double)topMethylRatesHaplotype1[i][virtualPos1], (double)parameters.bottomMethylRates[i][methPos]);
                        else if (ref == 'G') // Haplotype GX
                            streamPutMethLevel(methFile, contigName, methPos, '-', parameters.bottomMethylRates[i][methPos]);
                        else if (ref != 'C' && ref != 'G' && allSnpSets[1][i][p1].altBase == 'C') // XC
                            streamPutMethLevel(methFile, contigName, methPos, '+', topMethylRatesHaplotype1[i][virtualPos1]);
                        else if (ref != 'C' && ref != 'G' && allSnpSets[1][i][p1].altBase == 'G') // XG
                            streamPutMethLevel(methFile, contigName, methPos, '-', bottomMethylRatesHaplotype1[i][virtualPos1]);
                       
                        ++p1;
                        ++methPos;
                    }
                } 
            }
        }
        // For 1 haplotypes
        if (options.numHaplotypes == 1)
        {
            int diff0 = 0;    
            unsigned p0 = 0;

            // Iterate over snps in the one haplotype 
            while (p0 < length(allSnpSets[0][i]) )   
            {
                unsigned contigPos0;      
                Snp &currSnp0 = allSnpSets[0][i][p0];
                // Real contig position is calculated by virtual position + length of deletion - length of insertions up to this position
                contigPos0 = currSnp0.virtualPos + diff0;
                if (currSnp0.type == ERROR_TYPE_INSERT)  // If insertion, subtract length, iterate to next snp entry and continue
                {
                    diff0 -= currSnp0.length;
                    ++p0;
                    continue;
                }
                else if (currSnp0.type == ERROR_TYPE_DELETE)  // If deletion, add length, iterate to next snp entry and continue
                {
                    diff0 += currSnp0.length;
                    ++p0;
                    continue;
                }
                
                if (currSnp0.type == ERROR_TYPE_MISMATCH)                                                        
                {
                    streamPut(snpFile, contigName);
                    streamPut(snpFile, '\t');
                    streamPut(snpFile, contigPos0);
                    streamPut(snpFile, '\t');
                    streamPut(snpFile, '.');
                    streamPut(snpFile, '\t');
                    // REF
                    Dna5 const &ref = fragmentStore.contigStore[i].seq[contigPos0];
                    streamPut(snpFile, ref);
                    streamPut(snpFile, '\t');
                    // ALT
                    streamPut(snpFile, currSnp0.altBase);
                    streamPut(snpFile, '\t');
                    streamPut(snpFile, '.');
                    streamPut(snpFile, '\t');
                    streamPut(snpFile, "PASS");
                    streamPut(snpFile, '\t');
                    streamPut(snpFile, "AF=100");
                    streamPut(snpFile, '\n');

                    // Meth output
                    if (ref == 'C') SEQAN_ASSERT_NEQ(currSnp0.altBase, 'C');
                    else if (ref == 'G') SEQAN_ASSERT_NEQ(currSnp0.altBase, 'G');
                    if (currSnp0.altBase == 'C') // Haplotype C
                        streamPutMethLevel(methFile, contigName, contigPos0, '+', parameters.topMethylRatesHaplotype[i][currSnp0.virtualPos]);
                    else if (currSnp0.altBase == 'G') // Haplotype G
                        streamPutMethLevel(methFile, contigName, contigPos0, '-', parameters.bottomMethylRatesHaplotype[i][currSnp0.virtualPos]);
                    
                    ++p0;
                }
            }
        }
    } 
    return 0;
}


// Performs the actual read simulation.
template <typename TRNG, typename TOptions>
int simulateReadsMain(FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                      TRNG & rng,
                      TOptions const & options,
                      ModelParameters<IlluminaReadsBS> & parameters) {        // bs_change: not const, need to adjust a few parameters
    typedef FragmentStore<MyFragmentStoreConfig> TFragmentStore;
    typedef Value<TFragmentStore::TMatePairStore>::Type TMatePairStoreElement;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type	TAlignedElement;
    typedef typename TAlignedElement::TGapAnchors								TReadGapAnchors;

    if (options.verbose)
        std::cerr << "Simulating reads..." << std::endl;

    typedef Position<CharString>::Type TPos;

    // Number of reads comes from command line by default.  If sample
    // counts are given, we compute it from this instead.
    size_t numReads = options.numReads;
    if (length(parameters.sampleCounts) > 0u) {
        numReads = 0;
        for (unsigned i = 0; i < length(parameters.sampleCounts); ++i)
            numReads += parameters.sampleCounts[i];
    }

    // First, we randomly pick the haplotype for each read to be
    // simulated or read it from the sample counts in parameters.
    String<unsigned> haplotypeIds;
    // Pick random haplotype origin.
    reserve(haplotypeIds, numReads);
    for (size_t i = 0; i < numReads; ++i)
        appendValue(haplotypeIds, pickRandomNumber(rng, Pdf<Uniform<unsigned> >(0, options.numHaplotypes - 1)));

    // Maybe pick contig ids to sample from.
    String<unsigned> contigIds;
    if (length(parameters.sampleCounts) > 0) {
        for (unsigned i = 0; i < length(fragmentStore.contigNameStore); ++i) {
            resize(contigIds, length(contigIds) + parameters.sampleCounts[i], i);
            if (options.veryVerbose)
                std::cerr << parameters.sampleCounts[i] << " reads from haplotype " << i << "..." << std::endl;
        }
        shuffle(contigIds, rng);
    }

    // We do not build all haplotypes at once since this could cost a
    // lot of memory.
    //
    // TODO(holtgrew): Would only have to switch pointers to journals which is possible.
    //
    // for each haplotype id
    //   simulate haplotype
    //   for each simulation instruction for this haplotype:
    //     build simulated read
    reserve(fragmentStore.readSeqStore, numReads, Exact());
    reserve(fragmentStore.readNameStore, numReads, Exact());
    CharString readNameBuf;
    char outFileName[151];
    snprintf(outFileName, 150, "%s", toCString(options.readNamePrefix));

    // bs_change
    // we need snpSets for all haplotypes together to write vcf file 
    typedef String<Snp> TContigSnpSet;          // Snps of one contig
    typedef String<TContigSnpSet> TSnpSet;     // Snps of all contigs of one haplotype
    typedef String<TSnpSet> TAllHaploSnpSets;    // Snps of all haplotypes
    TAllHaploSnpSets allSnpSets;
    resize(allSnpSets, options.numHaplotypes);
    // We need both methylation haplotypes to write meth states to output
    // Store first haplotype methyl positions
    StringSet<String<double, Journaled<Alloc<> > > > topMethylRatesHaplotype0;
    StringSet<String<double, Journaled<Alloc<> > > > bottomMethylRatesHaplotype0;

    for (unsigned haplotypeId = 0; haplotypeId < options.numHaplotypes; ++haplotypeId) {

        std::cerr << "Simulating for haplotype #" << haplotypeId << "..." << std::endl;
        std::cout << "  Building haplotype..." << std::endl;
        StringSet<String<Dna5, Journaled<Alloc<> > > > haplotypeContigs;

        // bs_change
        String<TContigSnpSet> & snpSet = allSnpSets[haplotypeId];
        
        resize(snpSet, length(fragmentStore.contigNameStore));

        double buildStart = sysTime();

        //if (empty(options.vcfFile))
            // bs_change
            // need of parameters, adjust corresponding to haplotype changes
            buildHaplotype(haplotypeContigs, parameters, fragmentStore, snpSet, rng, options);
        
        // bs_change: Create copy, needed later to write output with methylation states
        if (options.writeVCFFile && haplotypeId == 0) // 
        {
            topMethylRatesHaplotype0 = parameters.topMethylRatesHaplotype; 
            bottomMethylRatesHaplotype0 = parameters.bottomMethylRatesHaplotype;
        }

        //else
        //    loadHaplotype(haplotypeContigs, fragmentStore, snpSet, options);  // for bs: do not load haplotype for the beginning
        std::cout << "  Finished haplotype creation in " << (sysTime() - buildStart) << 's' << std::endl;

        // TODO(holtgrew): Assigning of string set with compatible string should be possible.
        StringSet<String<Dna5> > haplotypeContigsCopy;
        for (unsigned i = 0; i < length(haplotypeContigs); ++i)
            appendValue(haplotypeContigsCopy, haplotypeContigs[i]);

        // Build partial sums over relative contig lengths so we can pick the contigs later on.
        size_t totalLength = 0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
            totalLength += length(fragmentStore.contigStore[i].seq);
        String<double> relativeContigLengths;
        resize(relativeContigLengths, length(fragmentStore.contigStore) + 1, Exact());
        front(relativeContigLengths) = 0.0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
            double l = static_cast<double>(length(fragmentStore.contigStore[i].seq));
            relativeContigLengths[i + 1] = l / totalLength;
        }
        std::partial_sum(begin(relativeContigLengths), end(relativeContigLengths), begin(relativeContigLengths));
        back(relativeContigLengths) = 1.0;

        // Simulate the reads...
        std::cerr << "  Simulating reads for haplotype #" << haplotypeId << "..." << std::endl;

//         std::cerr << "Journal: " << haplotypeContigs[0]._journalEntries << std::endl;

        Snp searchSnp;
        char tagBuffer[40];

        for (unsigned j = 0; j < length(haplotypeIds); ++j) {
            if (haplotypeIds[j] != haplotypeId)
                continue;  // Guard against instructions on wrong haplotype.

            // Build simulation instructions.
            String<ReadSimulationInstruction<IlluminaReadsBS> > instructions;
            // TODO(holtgrew): Pick contig id outside of instructions.
            size_t contigId = 0;
            bool fixedContigId = false;
            if (length(parameters.sampleCounts) > 0) {
                fixedContigId = true;
                if (options.generateMatePairs)
                    contigId = contigIds[length(fragmentStore.readSeqStore) / 2];
                else
                    contigId = contigIds[length(fragmentStore.readSeqStore)];
            }
            int res = buildReadSimulationInstruction(instructions, rng, haplotypeId, haplotypeContigsCopy, relativeContigLengths, contigId, fixedContigId, parameters, options);
            if (res != 0)
                return res;

            for (unsigned k = 0; k < length(instructions); ++k) {
                ReadSimulationInstruction<IlluminaReadsBS> & inst = instructions[k];

                // Reverse edit string and qualities string if this instruction simulates from reverse strand.
                if (!inst.isForward)
                {
                    reverse(inst.editString);
                    reverse(inst.qualities);
                }
                // Apply simulation instructions.
                SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore), length(fragmentStore.readNameStore));
                // Cut out segment from haplotype.
                String<Dna5Q> read = infix(haplotypeContigsCopy[inst.contigId], inst.beginPos, inst.endPos);
                String<Dna5Q> haplotypeInfix = read;  // Copy for printing later on.
                applySimulationInstructions(read, rng, inst, options);
                if (!inst.isForward)
                {
                    reverseComplement(read);
                    // Reconstruct edit string and qualities in sequencing direction.
                    reverse(inst.editString);
                    reverse(inst.qualities);
                }
                // Append read sequence to read seq store and mate pair to read name store.  This also yields the read
                // id.  We will generate and append the read name below, depending on the read id.
                unsigned readId;
                if (options.generateMatePairs)
                    readId = appendRead(fragmentStore, read, length(fragmentStore.matePairStore));
                else
                    readId = appendRead(fragmentStore, read);

                // Get expected begin/end position in the original sequence.
                TPos origBeginPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.beginPos);
                TPos origEndPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.endPos);
               
                unsigned numSNPs = 0;
                unsigned numIndels = 0;
                
                if (options.includeReadInformation && inst.contigId < length(snpSet))
                {
                    typedef Iterator<TContigSnpSet, Standard>::Type TIter;                    
                    TIter allBeg = begin(snpSet[inst.contigId], Standard());
                    TIter allEnd = end(snpSet[inst.contigId], Standard());

//                    searchSnp.contigId = inst.contigId;
                    searchSnp.virtualPos = inst.beginPos;
                    TIter snpBeg = std::lower_bound(allBeg, allEnd, searchSnp, SnpLess());
                    searchSnp.virtualPos = inst.endPos;
                    TIter snpEnd = std::upper_bound(allBeg, allEnd, searchSnp, SnpLess());
                    
                    // find first SNP that overlaps the read
                    while (snpBeg != allBeg && (snpBeg - 1)->virtualPos + (snpBeg - 1)->length > inst.beginPos)
                        --snpBeg;
                    
//                    if (readId == 88704)
//                    {
//                        std::cout << "begPos:\t" << inst.beginPos << std::endl;
//                        std::cout << "endPos:\t" << inst.endPos << std::endl;
//                        for (TIter it = snpBeg; it != snpEnd; ++it)
//                            std::cout << it->virtualPos << '\t' << (int)it->type << '\t' << it->length << std::endl;
//                    }

                
                    for (TIter it = snpBeg; it != snpEnd; ++it)
                    {
                        int len = it->length;

                        // cut at the right end of the read
                        if (it->virtualPos + len > inst.endPos)
                            len = inst.endPos - it->virtualPos;

                        // cut at the left end of the read
                        if (it->virtualPos <= inst.beginPos)
                        {
                            // deletions must occur in the read to affect it
                            if (it->type == ERROR_TYPE_DELETE)
                                continue;
                            len -= inst.beginPos - it->virtualPos;
                        }
                        
                        if (len <= 0) continue;
                        
                        if (it->type == ERROR_TYPE_MISMATCH)
                            numSNPs += len;
                        else if (it->type == ERROR_TYPE_INSERT)
                            numIndels += len;
                        else // if (it->type == ERROR_TYPE_DELETE)
                            // at least one base left and right of the deletion is required
                            if (it->virtualPos < inst.endPos)
                                numIndels += it->length;
                    }
                }
                

                // Generate read name.
                // TODO(holtgrew): Remove mateNum, not required?
                if (options.generateMatePairs)
                {
                    // Generate the mate num \in {1, 2}, randomly but consistent so two entries belonging together have
                    // different nums.
                    if (options.includeReadInformation)
                    {
                        resize(readNameBuf, 1024 + length(haplotypeInfix) + length(inst.editString) + length(inst.bsConversionString));
                        sprintf(&readNameBuf[0], "%s.%09u contig=%s haplotype=%u length=%lu orig_begin=%lu orig_end=%lu snps=%u indels=%u haplotype_infix=%s edit_string=", outFileName, readId / 2, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, static_cast<long unsigned>(length(read)), static_cast<long unsigned>(origBeginPos), static_cast<long unsigned>(origEndPos),numSNPs, numIndels, toCString(CharString(haplotypeInfix)));
                    }
                    else
                    {
                        resize(readNameBuf, 1024);
                        sprintf(&readNameBuf[0], "%s.%09u", outFileName, readId / 2);
                    }
                } else {
                    if (options.includeReadInformation)
                    {
                        resize(readNameBuf, 1024 + length(haplotypeInfix) + length(inst.editString) + length(inst.bsConversionString));
                        sprintf(&readNameBuf[0], "%s.%09u contig=%s haplotype=%u length=%lu orig_begin=%lu orig_end=%lu snps=%u indels=%u haplotype_infix=%s edit_string=", outFileName, readId, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, static_cast<long unsigned>(length(read)), static_cast<long unsigned>(origBeginPos), static_cast<long unsigned>(origEndPos), numSNPs, numIndels, toCString(CharString(haplotypeInfix)));
                    }
                    else
                    {
                        resize(readNameBuf, 1024);
                        sprintf(&readNameBuf[0], "%s.%09u", outFileName, readId);
                    }
                }
                if (options.includeReadInformation) 
                {
                    for (unsigned i = 0; i < length(inst.editString); ++i) {
                        char buffer[2] = "*";
                        buffer[0] = "MEID"[static_cast<int>(inst.editString[i])];
                        strcat(&readNameBuf[0], buffer);
                    }
                    // bs_change:
                    strcat(&readNameBuf[0], " bs_conversion_string=");
                    for (unsigned i = 0; i < length(inst.bsConversionString); ++i) {
                        char buffer[2] = "*";
                        buffer[0] = "01"[static_cast<int>(inst.bsConversionString[i])];
                        strcat(&readNameBuf[0], buffer);
                    }

                }
                // bs_change:
                CharString readName(&readNameBuf[0]);
                //if (options.includeReadInformation) {
                    
                if (inst.isFromTopStrand)
                    append(readName, " original_strand=top");
                else
                    append(readName, " original_strand=bottom");

                if (inst.isForward)
                    append(readName, " strand=forward");
                else
                    append(readName, " strand=reverse");
                //}
                appendValue(fragmentStore.readNameStore, readName);


                // Print info about read and haplotype.
                if (options.veryVerbose) {
                    std::cout << ",-- Read #" << readId << std::endl
                              << "| inst.beginPos    " << inst.beginPos << std::endl
                              << "| inst.endPos      " << inst.endPos << std::endl
                              << "| origBeginPos     " << origBeginPos << std::endl
                              << "| origEndPos       " << origEndPos << std::endl
                              << "| numSNPs          " << numSNPs << std::endl
                              << "| numIndels        " << numIndels << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos-1) << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos) << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos+1) << std::endl
                              << "| name:            " << readName << std::endl
                              << "| original infix:  " << infix(fragmentStore.contigStore[inst.contigId].seq, origBeginPos, origEndPos) << std::endl
                              << "| haplotype infix: " << infix(haplotypeContigsCopy[inst.contigId], inst.beginPos, inst.endPos) << std::endl
                              << "| read:            " << read << std::endl
                              << "`-- " << std::endl;
                }

                // Swap original begin and end position if from reverse strand.
                if (!inst.isForward)
                    std::swap(origBeginPos, origEndPos);

                // Add matches to aligned read store.
                if (options.generateMatePairs)
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos, length(fragmentStore.matePairStore));
                else
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos);
                
                sprintf(tagBuffer, "XE:i:%u\tXS:i:%u\tXI:i:%u", (inst.mismatchCount + inst.insCount + inst.delCount), numSNPs, numIndels);
                appendValue(fragmentStore.alignedReadTagStore, tagBuffer);

                // Adding mate pair information.
                if (options.generateMatePairs && readId % 2 == 1) {  // Only append mate pair info after simulating second mate.
                    // Append mate pair element to fragment store's mate pair store.
                    TMatePairStoreElement matePair;
                    matePair.readId[0] = readId - 1;
                    matePair.readId[1] = readId;
                    appendValue(fragmentStore.matePairStore, matePair);
                }
            }
			if (options.generateMatePairs)
            {
                // When generating mate pairs, an even number of reads is generated in each step.
 				SEQAN_ASSERT_EQ(length(fragmentStore.alignedReadStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readNameStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore) % 2, 0u);
			}
		}
        // bs_change
        // Only store snpSets etc. of all haplotypes, if vcf file output is required
        if (!options.writeVCFFile)
        {
            clear(allSnpSets[haplotypeId]);
            clear(topMethylRatesHaplotype0);
            clear(bottomMethylRatesHaplotype0);
        }
    }
    std::cout << "Write vcf output..." << std::endl;
    // bs_change
    if (options.writeVCFFile && options.numHaplotypes <= 2)
        writeVCFOutput(fragmentStore, allSnpSets, topMethylRatesHaplotype0, bottomMethylRatesHaplotype0, parameters.topMethylRatesHaplotype, parameters.bottomMethylRatesHaplotype, parameters, options);

    std::cout << "... Done with vcf output..." << std::endl;


    // Last but not least, convert the matches collected before to a global alignment.
    convertMatchesToGlobalAlignment(fragmentStore, Score<int, EditDistance>(), True());
	
	// AlignedReadLayout layout;
	// layoutAlignment(layout, fragmentStore);
	// printAlignment(std::cout, Raw(), layout, fragmentStore, 0, 0, 300, 0, 100);
    
    if (options.verbose)
        std::cerr << "Simulated " << length(fragmentStore.readSeqStore) << " reads" << std::endl;

    return 0;
}



#endif  // #ifndef SANDBOX_KRAKAU_APPS_BS_MASON_SIMULATE_ILLUMINA_BS_H_
