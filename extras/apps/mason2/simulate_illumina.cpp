// ==========================================================================
//                         Mason - A Read Simulator
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "sequencing.h"
#include "simulate_illumina_error_frequencies.h"

// ===========================================================================
// Class IlluminaSequencingOptions
// ===========================================================================

class IlluminaModel
{
public:
    // Probabilities for a mismatch at a given position.
    seqan::String<double> mismatchProbabilities;

    // Standard deviations for the normal distributions of base qualities for the mismatch case.
    seqan::String<double> mismatchQualityMeans;
    // Standard deviations for the normal distributions of base qualities for the mismatch case.
    seqan::String<double> mismatchQualityStdDevs;

    // Standard deviations for the normal distributions of base qualities for the non-mismatch case.
    seqan::String<double> qualityMeans;
    // Standard deviations for the normal distributions of base qualities for the non-mismatch case.
    seqan::String<double> qualityStdDevs;

    double const * baseErrorFreqsFrom;
    double scalingFactorFrom;
    double const * seqErrorFreqs;

    double const * insErrorFreqsM;
    double const * delErrorFreqsM;
    double scalingFactorDel;

    IlluminaModel()
    {
        (*this).baseErrorFreqsFrom = seqan::BaseErrorFreqsFrom<double, seqan::BsSimple>::getData();
        (*this).seqErrorFreqs = seqan::SeqErrorFreqs<double, seqan::BsSimple>::getData();
        (*this).insErrorFreqsM = seqan::InsErrorFreqsM<double, seqan::BsSimple>::getData();
        (*this).delErrorFreqsM = seqan::DelErrorFreqsM<double, seqan::BsSimple>::getData();
        scalingFactorFrom = 5.0;
        scalingFactorDel = 5.0;
    }
};

// ===========================================================================
// Class IlluminaSequencingSimulator
// ===========================================================================

// ---------------------------------------------------------------------------
// Constructor IlluminaSequencingSimulator::IlluminaSequencingSimulator
// ---------------------------------------------------------------------------

IlluminaSequencingSimulator::IlluminaSequencingSimulator(TRng & rng,
                                                         TRng & methRng,
                                                         SequencingOptions const & seqOptions,
                                                         IlluminaSequencingOptions const & illuminaOptions) :
        SequencingSimulator(rng, methRng, seqOptions), illuminaOptions(illuminaOptions),
        model(new IlluminaModel())
{
    this->_initModel();
}

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::_initModel()
// ---------------------------------------------------------------------------

void IlluminaSequencingSimulator::_initModel()
{
    if (illuminaOptions.nonSimpleSubstErrorsFrom) 
    {
        model->baseErrorFreqsFrom = seqan::BaseErrorFreqsFrom<double, seqan::BsNonSimple>::getData();
        model->scalingFactorFrom = 3.5;
    }
    if (illuminaOptions.nonSimpleSubstErrors) model->seqErrorFreqs = seqan::SeqErrorFreqs<double, seqan::BsNonSimple>::getData();
    if (illuminaOptions.nonSimpleInsErrors) model->insErrorFreqsM = seqan::InsErrorFreqsM<double, seqan::BsNonSimple>::getData();
    if (illuminaOptions.nonSimpleDelErrors)
    {
        model->delErrorFreqsM = seqan::DelErrorFreqsM<double, seqan::BsNonSimple>::getData();
        model->scalingFactorDel = 3.5;
    }

    std::cout << "illuminaOptions.nonSimpleSubstErrorsFrom: " << (illuminaOptions.nonSimpleSubstErrorsFrom ? "true":"false") << std::endl;
    std::cout << "illuminaOptions.nonSimpleSubstErrors: " << (illuminaOptions.nonSimpleSubstErrors ? "true":"false") << std::endl;
    std::cout << "illuminaOptions.nonSimpleInsErrors: " << (illuminaOptions.nonSimpleInsErrors ? "true":"false") << std::endl;
    std::cout << "illuminaOptions.nonSimpleDelErrors: " << (illuminaOptions.nonSimpleDelErrors ? "true":"false") << std::endl;

    // Compute mismatch probabilities, piecewise linear function.
    resize(model->mismatchProbabilities, illuminaOptions.readLength);
    // Compute probability at raise point.
    double y_r = 2 * illuminaOptions.probabilityMismatch - illuminaOptions.positionRaise * illuminaOptions.probabilityMismatchBegin - illuminaOptions.probabilityMismatchEnd + illuminaOptions.probabilityMismatchEnd * illuminaOptions.positionRaise;
    if (illuminaOptions.verbosity >= 2)
    {
        std::cerr << "Illumina error curve:\n"
                  << "  (0, " << illuminaOptions.probabilityMismatchBegin << ") -- (" << illuminaOptions.positionRaise << ", " << y_r << ") -- (1, " << illuminaOptions.probabilityMismatchEnd << ")\n";
    }
    // std::cout << "y_r = " << y_r << std::endl;
    // Compute mismatch probability at each base.
    if (!empty(illuminaOptions.probabilityMismatchFile))
    {
        // Open file.
        std::fstream file;
        file.open(toCString(illuminaOptions.probabilityMismatchFile), std::ios_base::in);
        if (!file.is_open())
        {
            std::cerr << "Failed to load mismatch probabilities from " << illuminaOptions.probabilityMismatchFile << std::endl;
            // return 1;
        }
        // Load probabilities.
        double x;
        file >> x;
        unsigned i;
        for (i = 0; i < illuminaOptions.readLength && !file.eof(); ++i) {
            model->mismatchProbabilities[i] = x;
            file >> x;
        }
        if (i != illuminaOptions.readLength)
        {
            std::cerr << "Not enough mismatch probabilites in " << illuminaOptions.probabilityMismatchFile << " (" << i << " < " << illuminaOptions.readLength << ")!" << std::endl;
            // return 1;
        }
    } else {
        // Use piecewise linear function for mismatch probability simulation.
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
            double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
            if (x < illuminaOptions.positionRaise) {
                double b = illuminaOptions.probabilityMismatchBegin;
                double m = (y_r - illuminaOptions.probabilityMismatchBegin) / illuminaOptions.positionRaise;
                model->mismatchProbabilities[i] = m * x + b;
                // std::cout << "model->mismatchProbabilities[" << i << "] = " << model->mismatchProbabilities[i] << std::endl;
            } else {
                double b = y_r;
                double m = (illuminaOptions.probabilityMismatchEnd - y_r) / (1 - illuminaOptions.positionRaise);
                x -= illuminaOptions.positionRaise;
                model->mismatchProbabilities[i] = m * x + b;
                // std::cout << "model->mismatchProbabilities[" << i << "] = " << model->mismatchProbabilities[i] << std::endl;
            }
        }
    }
    if (illuminaOptions.probabilityMismatchScale != 1.0) {
        for (unsigned i = 0; i < illuminaOptions.readLength; ++i)
            model->mismatchProbabilities[i] *= illuminaOptions.probabilityMismatchScale;
    }

    // Compute match/mismatch means and standard deviations.
    resize(model->mismatchQualityMeans, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.meanMismatchQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.meanMismatchQualityEnd - illuminaOptions.meanMismatchQualityBegin);
        model->mismatchQualityMeans[i] = m * x + b;
        // std::cout << "model->mismatchQualityMeans[" << i << "] = " << model->mismatchQualityMeans[i] << std::endl;
    }
    resize(model->mismatchQualityStdDevs, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.stdDevMismatchQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.stdDevMismatchQualityEnd - illuminaOptions.stdDevMismatchQualityBegin);
        model->mismatchQualityStdDevs[i] = m * x + b;
        // std::cout << "model->mismatchQualityStdDevs[" << i << "] = " << model->mismatchQualityStdDevs[i] << std::endl;
    }
    resize(model->qualityMeans, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.meanQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.meanQualityEnd - illuminaOptions.meanQualityBegin);
        model->qualityMeans[i] = m * x + b;
        // std::cout << "model->qualityMeans[" << i << "] = " << model->qualityMeans[i] << std::endl;
    }
    resize(model->qualityStdDevs, illuminaOptions.readLength);
    for (unsigned i = 0; i < illuminaOptions.readLength; ++i) {
        double b = illuminaOptions.stdDevQualityBegin;
        double x = static_cast<double>(i) / (illuminaOptions.readLength - 1);
        double m = (illuminaOptions.stdDevQualityEnd - illuminaOptions.stdDevQualityBegin);
        model->qualityStdDevs[i] = m * x + b;
        // std::cout << "model->qualityStdDevs[" << i << "] = " << model->qualityStdDevs[i] << std::endl;
    }
}

// ---------------------------------------------------------------------------
// Function _simulateRead()
// ---------------------------------------------------------------------------

namespace {

// Simulate the characters that polymorphisms turn into and inserted characters.
//
// Through the usage of ModifiedString, we will always go from the left to the right end.
template <typename TFrag>
void _simulateSequence(TRead & read, TRng & rng, TFrag const & frag,
                       TCigarString const & cigar)
{
    clear(read);

    typedef typename seqan::Iterator<TFrag>::Type TFragIter;
    TFragIter it = begin(frag, seqan::Standard());

    for (unsigned i = 0; i < length(cigar); ++i)
    {
        //unsigned numSimulate = 0;
        if (cigar[i].operation == 'M')
        {
            for (unsigned j = 0; j < cigar[i].count; ++j, ++it)
                appendValue(read, *it);
            continue;
        }
        else if (cigar[i].operation == 'D')
        {
            it += cigar[i].count;
            continue;
        }

        // Otherwise, we have insertions or mismatches.
        for (unsigned j = 0; j < cigar[i].count; ++j)
        {
            // Pick a value between 0 and 1.
            double x = 1.0;
            while (x == 1.0)
                x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
            int num = static_cast<int>(x / 0.25);

            // NOTE: We can only insert CGAT, but we can have a polymorphism to N.

            if (cigar[i].operation == 'I')
                appendValue(read, seqan::Dna5(num));
            else
                appendValue(read, seqan::Dna5(num + (num == ordValue(*it))));
        }

        if (cigar[i].operation == 'X')
            it += cigar[i].count;
    }
}

template<typename TFrag, typename TModel>
void _simulateSequence(TRead & read, TRng & rng, TFrag const & frag,
                       TCigarString const & cigar, TModel const & model, SequencingSimulationInfo & info)
{
    clear(read);

    if (info.debugRead)
    {
        std::cout << "cigar: ";
        for (unsigned i = 0; i < length(cigar); ++i)
        {
            if (cigar[i].operation == 'D') std::cout << 'D' << cigar[i].count;
            else if (cigar[i].operation == 'I')  std::cout << 'I' << cigar[i].count;
            else if (cigar[i].operation == 'M')  std::cout << 'M' << cigar[i].count;
            else if (cigar[i].operation == 'X')  std::cout << 'X' << cigar[i].count;
        }
        std::cout << std::endl;
    }

    typedef typename seqan::Iterator<TFrag>::Type TFragIter;
    TFragIter it = begin(frag, seqan::Standard());

    for (unsigned i = 0; i < length(cigar); ++i)
    {
        //unsigned numSimulate = 0;
        if (cigar[i].operation == 'M')
        {
            for (unsigned j = 0; j < cigar[i].count; ++j, ++it)
                appendValue(read, *it);
            continue;
        }
        else if (cigar[i].operation == 'D')
        {
            it += cigar[i].count;
            info.countD += cigar[i].count;
            continue;
        }

        // Otherwise, we have insertions or mismatches.
        if (cigar[i].operation == 'I')
        {
            for (unsigned j = 0; j < cigar[i].count; ++j)
            {
                double pA = model->insErrorFreqsM[0];
                double pC = model->insErrorFreqsM[1];
                double pG = model->insErrorFreqsM[2];
                double pT = model->insErrorFreqsM[3];

                int num;
                double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));

                if (x < pA) num = 0;
                else if (x < pA + pC) num = 1;
                else if (x < pA + pC + pG) num = 2;
                else if (x < pA + pC + pG + pT) num = 3;
                else  num = 4;

                appendValue(read, seqan::Dna5(num));
                ++info.countI;
            }
            continue;
        }

        if (cigar[i].operation == 'X')
        {
            for (unsigned j = 0; j < cigar[i].count; ++j)
            {
                // in case of gold: doesn't matter, which base, only cigar counts
                double pA = model->seqErrorFreqs[ordValue(*it)*5 + 0];
                double pC = model->seqErrorFreqs[ordValue(*it)*5 + 1];
                double pG = model->seqErrorFreqs[ordValue(*it)*5 + 2];
                double pT = model->seqErrorFreqs[ordValue(*it)*5 + 3];

                int num;
                double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));

                if (x < pA) num = 0;
                else if (x < pA + pC) num = 1;
                else if (x < pA + pC + pG) num = 2;
                else if (x < pA + pC + pG + pT) num = 3;
                else num = 4;
                
                //std::cout << (*it) << " pA: " << pA << " pC:" << pC << " pG: " << pG << " pT: " << (1-pA - pC - pG) << std::endl;
                appendValue(read, seqan::Dna5(num));
                ++info.count;
            }
            it += cigar[i].count;
        }
    }
    //std::cout << "read: " << read << std::endl;

}

template<typename TFrag, typename TModel>
void _simulateCigar2(TFrag const & frag, TRng & rng, TCigarString & cigar, 
                    unsigned len, seqan::String<bool> bsEditString, TModel const & model, IlluminaSequencingOptions const & illuminaOptions,  SequencingSimulationInfo & info)  // TODO rm info
{
    clear(cigar);
    if (info.debugRead)
    {
        std::cout << "frag: " << frag << std::endl;
        std::cout << "bsEdit: " << bsEditString << std::endl;
    }
    typedef typename seqan::Iterator<TFrag>::Type TFragIter;
    typedef typename seqan::Iterator<seqan::String<bool> >::Type TGoldIter;
    TFragIter it = begin(frag, seqan::Standard());
    TGoldIter itG = begin(bsEditString, seqan::Standard());    // For build gold sam, true if bs conversion in final bs read
    
    for (int i = 0; i < (int)len;)
    {
        double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
        double pInsert = illuminaOptions.probabilityInsert;
        // Use theoretical bs base for error freqs 
        double pMismatch;
        if (*itG)  // take care that error is different than original? shouldn't matter  
        {
            if ((*it) == 'C') pMismatch = illuminaOptions.probabilityMismatch * model->baseErrorFreqsFrom[3] * model->scalingFactorFrom;     // Use del rate for T
            else if ((*it) == 'G') pMismatch = illuminaOptions.probabilityMismatch * model->baseErrorFreqsFrom[0] * model->scalingFactorFrom;     
            else pMismatch = illuminaOptions.probabilityMismatch * model->baseErrorFreqsFrom[ordValue(*it)] * model->scalingFactorFrom;
        } 
        else pMismatch = illuminaOptions.probabilityMismatch * model->baseErrorFreqsFrom[ordValue(*it)] * model->scalingFactorFrom; 
        double pDelete;
        if (*itG)   
        {
            if ((*it) == 'C') pDelete = illuminaOptions.probabilityDelete * model->delErrorFreqsM[3] * model->scalingFactorDel;     // Use del rate for T
            else if ((*it) == 'G') pDelete = illuminaOptions.probabilityDelete * model->delErrorFreqsM[0] *  model->scalingFactorDel;      
            else pDelete = illuminaOptions.probabilityDelete * model->delErrorFreqsM[ordValue(*it)] *  model->scalingFactorDel;
        }   
        else  pDelete = illuminaOptions.probabilityDelete * model->delErrorFreqsM[ordValue(*it)] *  model->scalingFactorDel;     

        double pMatch = 1.0 - pMismatch - pInsert - pDelete;

        if (info.debugRead)
        {
            std::cout << "x: " << x << std::endl;
            std::cout << "pMismatch: " << pMismatch << std::endl;
            std::cout << "pInsert: " << pInsert << std::endl;
            std::cout << "pDelete: " << pDelete << std::endl;
            std::cout << "pMatch: " << pMatch << std::endl;
        }


        if (x < pMatch)  // match
        {
            i += appendOperation(cigar, 'M').first;
            ++it;
            ++itG;
        }
        else if (x < pMatch + pMismatch)  // point polymorphism
        {
            i += appendOperation(cigar, 'X').first;
            ++it;
            ++itG;
        }
        else if (x < pMatch + pMismatch + pInsert) // insertion
        {
            i += appendOperation(cigar, 'I').first;
        }
        else  // deletion
        {
            i += appendOperation(cigar, 'D').first;
            ++it;
            ++itG;
        }
    }
    //std::cout << std::endl;
}

}  // namespace (anonymous)

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::simulateRead()
// ---------------------------------------------------------------------------

// Actually simulate read and qualities from fragment and direction forward/reverse strand.
void IlluminaSequencingSimulator::simulateRead(TRead & seq, TQualities & quals, SequencingSimulationInfo & info,
                                               TFragment const & frag, Direction dir, Strand strand)
{
    // std::cerr << "simulateRead(" << (char const *)(dir == LEFT ? "L" : "R") << ", " << (char const *)(strand == FORWARD ? "-->" : "<--") << ")\n";
    // Simulate sequencing operations.
    TCigarString cigar;
    typedef seqan::ModifiedString<seqan::ModifiedString<TFragment, seqan::ModView<seqan::FunctorComplement<seqan::Dna5> > >, seqan::ModReverse> TRevCompFrag;
    typedef seqan::ModifiedString<seqan::String<bool>, seqan::ModReverse>                                                                       TRevBoolString;

    if (illuminaOptions.uniformSequencingErrors) _simulateCigar(cigar);
    else 
    {   // Currently only works for R1 --> <-- R2, since we don't know read length in advance
        if ((dir == LEFT) && (strand == FORWARD))
            _simulateCigar2(frag, rng, cigar, this->readLength(), info.bsEditString, model, illuminaOptions, info);
        else //if ((dir == RIGHT) && (strand == REVERSE))
            _simulateCigar2(TRevCompFrag(frag), rng, cigar, this->readLength(), TRevBoolString(info.bsEditString), model, illuminaOptions, info);
        /*else
        {
            std::cerr << "Dir: " << ((dir == LEFT) ? "LEFT":"RIGHT") << std::endl;
            std::cerr << "strand: " << ((strand == FORWARD) ? "FORWARD":"REVERSE") << std::endl;

            throw std::runtime_error("This should not happend!");
        }*/

    }
    unsigned lenInRef = 0;
    _getLengthInRef(cigar, lenInRef);

    if (lenInRef > length(frag))
    {
        throw std::runtime_error("Illumina read is too long, increase fragment length");
    }

    // Simulate sequence (materialize mismatches and insertions).
    if (illuminaOptions.uniformSequencingErrors) 
    {
        if ((dir == LEFT) && (strand == FORWARD))
            _simulateSequence(seq, rng, prefix(frag, lenInRef), cigar);
        else if ((dir == LEFT) && (strand == REVERSE))
            _simulateSequence(seq, rng, TRevCompFrag(prefix(frag, lenInRef)), cigar);
        else if ((dir == RIGHT) && (strand == FORWARD))
            _simulateSequence(seq, rng, suffix(frag, length(frag) - lenInRef), cigar);
        else  // ((dir == RIGHT) && (strand == REVERSE))
            _simulateSequence(seq, rng, TRevCompFrag(suffix(frag, length(frag) - lenInRef)), cigar);
    }
    else
    {
        if ((dir == LEFT) && (strand == FORWARD))
            _simulateSequence(seq, rng, prefix(frag, lenInRef), cigar, model, info);
        else  
            _simulateSequence(seq, rng, TRevCompFrag(suffix(frag, length(frag) - lenInRef)), cigar, model, info);
    }

    // Simulate qualities.
    _simulateQualities(quals, cigar);
    SEQAN_ASSERT_EQ(length(seq), length(quals));

    // // Reverse qualities if necessary.
    // if (strand == REVERSE)
    //     reverse(quals);

    // Write out extended sequencing information info if configured to do so.  We always write out the sample position
    // and alignment information.
    info.cigar = cigar;
    unsigned len = 0;
    _getLengthInRef(cigar, len);
    info.beginPos = (dir == LEFT) ? beginPosition(frag) : (beginPosition(frag) + length(frag) - len);
    info.isForward = (strand == FORWARD);
    // std::cerr << "  beginPos=" << info.beginPos - beginPosition(frag) << "\n";

    if (seqOptions->embedReadInfo)
    {
        if (dir == LEFT)
            info.sampleSequence = prefix(frag, len);
        else
            info.sampleSequence = suffix(frag, length(frag) - len);
        if (strand == REVERSE)
            reverseComplement(info.sampleSequence);
    }
    // std::cerr << "  sampleSequence  =" << info.sampleSequence << " beginPos: " << info.beginPos << std::endl;
    // std::cerr << "  fragRC=" << TRevCompFrag(frag) << "\n";
    // std::cerr << "  seq=" << seq << "\tquals=" << quals << "\n";
}

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::_simulateQualities()
// ---------------------------------------------------------------------------

// Simulate PHRED qualities from the CIGAR string.
void IlluminaSequencingSimulator::_simulateQualities(TQualities & quals, TCigarString const & cigar)
{
    clear(quals);

    unsigned pos = 0;
    for (unsigned i = 0; i < length(cigar); ++i)
    {
        for (unsigned j = 0; j < cigar[i].count; ++j)
        {
            int q = 0;
            if (cigar[i].operation == 'M')
            {
                seqan::Pdf<seqan::Normal> pdf(model->qualityMeans[pos], model->qualityStdDevs[pos]);
                q = static_cast<int>(pickRandomNumber(rng, pdf));
                ++pos;
            }
            else if (cigar[i].operation == 'I' || cigar[i].operation == 'X')
            {
                seqan::Pdf<seqan::Normal> pdf(model->mismatchQualityMeans[pos], model->mismatchQualityStdDevs[pos]);
                q = static_cast<int>(pickRandomNumber(rng, pdf));
                ++pos;
            }
            else
            {
                // Deletion/padding, no quality required.
                continue;
            }
            q = std::max(0, std::min(q, 40));  // limit quality to 0..40
            appendValue(quals, (char)('!' + q));
        }
    }
}

// ---------------------------------------------------------------------------
// Function IlluminaSequencingSimulator::_simulateCigar()
// ---------------------------------------------------------------------------

// Simulate CIGAR string.  We can do this with position specific parameters only and thus independent of any
// context.
void IlluminaSequencingSimulator::_simulateCigar(TCigarString & cigar)
{
    clear(cigar);
    unsigned len = this->readLength();

    for (int i = 0; i < (int)len;)
    {
        double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double> >(0, 1));
        double pMismatch = model->mismatchProbabilities[i];
        double pInsert   = illuminaOptions.probabilityInsert;
        double pDelete   = illuminaOptions.probabilityDelete;
        double pMatch    = 1.0 - pMismatch - pInsert - pDelete;

        // Simulate mutation/insertion/deletion events.  If possible we reuse the last CIGAR entry.  Adjacent
        // insertion/deletion pairs cancel each other out.

        // TODO(holtgrew): No indel at beginning or ending! Same for other simulators!

        if (x < pMatch)  // match
            i += appendOperation(cigar, 'M').first;
        else if (x < pMatch + pMismatch)  // point polymorphism
            i += appendOperation(cigar, 'X').first;
        else if (x < pMatch + pMismatch + pInsert) // insertion
            i += appendOperation(cigar, 'I').first;
        else  // deletion
            i += appendOperation(cigar, 'D').first;
    }
}

// ============================================================================
// Class SequencingSimulatorFactory
// ============================================================================

// ----------------------------------------------------------------------------
// Function SequencingSimulatorFactory::make()
// ----------------------------------------------------------------------------

std::SEQAN_AUTO_PTR_NAME<SequencingSimulator> SequencingSimulatorFactory::make()
{
    std::SEQAN_AUTO_PTR_NAME<SequencingSimulator> res;

    switch (seqOptions.sequencingTechnology)
    {
        case SequencingOptions::ILLUMINA:
            res.reset(new IlluminaSequencingSimulator(rng, methRng, seqOptions, illuminaOptions));
            break;
        case SequencingOptions::SANGER:
            res.reset(new SangerSequencingSimulator(rng, methRng, seqOptions, sangerOptions));
            break;
        case SequencingOptions::ROCHE_454:
            res.reset(new Roche454SequencingSimulator(rng, methRng, seqOptions, roche454Options));
            break;
    }

    return res;
}
