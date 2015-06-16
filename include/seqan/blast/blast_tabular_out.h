// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// This file contains routines to generate BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_BLAST_BLAST_TABULAR_WRITE_H_
#define SEQAN_BLAST_BLAST_TABULAR_WRITE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Type BlastTabularFileOut
// ----------------------------------------------------------------------------

/*!
 * @class BlastTabularFileOut
 * @signature template <typename TBlastIOContext>
 * using BlastTabularFileOut = FormattedFile<BlastTabular, Output, TBlastIOContext>;
 * @extends FormattedFileOut
 * @headerfile <seqan/blast.h>
 * @brief FormattedFileOut abstraction for @link BlastTabular @endlink
 *
 * This is a @link FormattedFile @endlink specialization for writing @link BlastTabular @endlink formats. For details
 * on how to influence the writing of files and how to differentiate between the tabular format without headers and the
 * one with headers, see @link BlastIOContext @endlink.
 * Please note that you have specify the type of the context as a template parameter to BlastTabularFileOut, see the example
 * below.
 *
 * @section Example
 *
 * The following short program creates the pairwise alignments between three query sequences and two database sequences,
 * it computes the e-values, sorts the matches and prints all results that score above a threshold. <i>The same example
 * is used for @link BlastTabularFileOut @endlink and @link BlastReportFileOut @endlink, you only need to change one line.</i>
 *
 * @include demos/blast/blast_out_example.cpp
 *
 * The file generated in /tmp/output.blast looks like this:
 *
 * @include demos/blast/blast_out_example.tabular
 *
 * @see BlastRecord
 */

template <typename TBlastIOContext = BlastIOContext<>>
using BlastTabularFileOut = FormattedFile<BlastTabular, Output, TBlastIOContext>;

// ============================================================================
// Metafunctions and global const-expressions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function guessFormat()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool guessFormat(FormattedFile<BlastTabular, Output, TSpec> &)
{
    return true;
}

// ----------------------------------------------------------------------------
// Function _writeFieldLabels()
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeFieldLabels(TFwdIterator & stream,
                  BlastIOContext<TScore, p, h> & context,
                  BlastTabular const &)
{
    if (!empty(context.fieldsAsStrings)) // give preference to string labels
    {
         write(stream, context.fieldsAsStrings, ", ");
    }
    else
    {
        for (auto it = seqan::begin(context.fields),
                  itB = it,
                  itEnd = seqan::end(context.fields);
            it != itEnd;
            ++it)
        {
            if (it != itB)
                write(stream, ", ");

            write(stream, BlastMatchField<>::columnLabels[static_cast<uint8_t>(*it)]);
        }
    }

    writeValue(stream, '\n');
}

// ----------------------------------------------------------------------------
// Function _writeCommentLines()
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename TScore,
          typename TMatch,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeCommentLinesWithoutColumnLabels(TFwdIterator & stream,
                                      BlastIOContext<TScore, p, h> & context,
                                      BlastRecord<TMatch> const & r,
                                      BlastTabular const & /*tag*/)
{
    write(stream, "# ");
    if (empty(context.versionString))
        context._setDefaultVersionString();
    write(stream, context.versionString);

    write(stream, "\n# Query: ");
    write(stream, r.qId);
    write(stream, "\n# Database: ");
    write(stream, context.dbName);
    write(stream, '\n');
}

//NOTE(h-2): dox disabled to clean-up interface
/*
 * @fn BlastTabular#_writeCommentLines
 * @headerfile seqan/blast.h
 * @brief write the header of a @link BlastRecord @endlink to file
 * @signature void _writeCommentLines(stream, context, blastRecord, blastTabular)
 *
 * @param[in,out] stream       The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in,out] context      A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastRecord  The @link BlastRecord @endlink whose header you want to print.
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 *
 * @section Remarks
 *
 * This function writes the header of a record unless the context's tabularSpec
 * (@link BlastIOContext#getBlastTabularSpec @endlink) is set to BlastTabularSpec::NO_COMMENTS (in which case this is
 * a NOOP).
 *
 * If context.@link BlastIOContext::versionString @endlink is set, this will be written,
 * otherwise one is generated. If either context.@link BlastIOContext::fields @endlink or
 * context.@link BlastIOContext::fieldsAsStrings @endlink is specified these will be printed
 * as column labels. If both are specified than @link BlastIOContext::fieldsAsStrings @endlink
 * are given preference. Please note that it is recommended to use
 * @link BlastIOContext::fields @endlink and not @link BlastIOContext::fieldsAsStrings @endlink to stay
 * "standards"-compliant. Also only @link BlastIOContext::fields @endlink has an influence on
 * the values printed by @link BlastTabular#_writeMatch @endlink.
 *
 * @throw IOError On low-level I/O errors.
 *
 * @see BlastRecord
 * @see BlastIOContext
 */

template <typename TFwdIterator,
          typename TScore,
          typename TMatch,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeCommentLines(TFwdIterator & stream,
                  BlastIOContext<TScore, p, h> & context,
                  BlastRecord<TMatch> const & r,
                  BlastTabular const & /*tag*/)
{
    ++context.numberOfRecords;

    if (context.tabularSpec == BlastTabularSpec::NO_COMMENTS)
        return;

    _writeCommentLinesWithoutColumnLabels(stream, context, r, BlastTabular());

    if (SEQAN_LIKELY(!context.legacyFormat))
    {
        // only write fields line if matches will follow
        if (length(r.matches) > 0)
        {
            write(stream, "# Fields: ");
            _writeFieldLabels(stream, context, BlastTabular());
        }
        // write # hits line
        write(stream, "# ");
        write(stream, length(r.matches));
        write(stream, " hits found\n");
    }
    else
    {
        write(stream, "# Fields: ");
        write(stream, BlastMatchField<>::legacyColumnLabels);
        write(stream, '\n');

        #if SEQAN_ENABLE_DEBUG
        if ((length(context.fields) != 1) || (context.fields[0] != BlastMatchField<>::Enum::STD))
            std::cerr << "Warning: custom fields set, but will be ignored, because legacyFormat is also set.\n";
        #endif
    }
}

// ----------------------------------------------------------------------------
// Function _writeField() [match object given]
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename TScore,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeField(TFwdIterator & s,
            BlastIOContext<TScore, p, h> & context,
            BlastMatch<TAlign, TPos, TQId, TSId> const & match,
            typename BlastMatchField<>::Enum const fieldId,
            BlastTabular const &)
{
    switch (fieldId)
    {
        case BlastMatchField<>::Enum::STD:
             // STD is handled from the calling function
            break;
        case BlastMatchField<>::Enum::Q_SEQ_ID:
            write(s, prefix(match.qId,
                            std::find(begin(match.qId, Standard()), end(match.qId, Standard()), ' ')
                            - begin(match.qId, Standard()))); // truncate at first ' '
            break;
//         case ENUM::Q_GI: write(s,  * ); break;
//         case ENUM::Q_ACC: write(s,  * ); break;
//         case ENUM::Q_ACCVER: write(s,  * ); break;
        case BlastMatchField<>::Enum::Q_LEN:
            write(s, match.qLength);
            break;
        case BlastMatchField<>::Enum::S_SEQ_ID:
            write(s, prefix(match.sId,
                            std::find(begin(match.sId, Standard()), end(match.sId, Standard()), ' ')
                            - begin(match.sId, Standard()))); // truncate at first ' '
            break;
//         case ENUM::S_ALL_SEQ_ID: write(s,  * ); break;
//         case ENUM::S_GI: write(s,  * ); break;
//         case ENUM::S_ALL_GI: write(s,  * ); break;
//         case ENUM::S_ACC: write(s,  * ); break;
//         case ENUM::S_ACCVER: write(s,  * ); break;
//         case ENUM::S_ALLACC: write(s,  * ); break;
        case BlastMatchField<>::Enum::S_LEN:
            write(s, match.sLength);
            break;
        case BlastMatchField<>::Enum::Q_START:
        {
            TPos effectiveQStart    = match.qStart;
            TPos effectiveQEnd      = match.qEnd;
            _untranslateQPositions(effectiveQStart, effectiveQEnd, match.qFrameShift, match.qLength,
                                   context.blastProgram);
            write(s, effectiveQStart);
        } break;
        case BlastMatchField<>::Enum::Q_END:
        {
            TPos effectiveQStart    = match.qStart;
            TPos effectiveQEnd      = match.qEnd;
            _untranslateQPositions(effectiveQStart, effectiveQEnd, match.qFrameShift, match.qLength,
                                   context.blastProgram);
            write(s, effectiveQEnd);
        } break;
        case BlastMatchField<>::Enum::S_START:
        {
            TPos effectiveSStart    = match.sStart;
            TPos effectiveSEnd      = match.sEnd;
            _untranslateSPositions(effectiveSStart, effectiveSEnd, match.sFrameShift, match.sLength,
                                   context.blastProgram);
            write(s, effectiveSStart);
        } break;
        case BlastMatchField<>::Enum::S_END:
        {
            TPos effectiveSStart    = match.sStart;
            TPos effectiveSEnd      = match.sEnd;
            _untranslateSPositions(effectiveSStart, effectiveSEnd, match.sFrameShift, match.sLength,
                                   context.blastProgram);
            write(s, effectiveSEnd);
        } break;
//         case ENUM::Q_SEQ: write(s,  * ); break;
//         case ENUM::S_SEQ: write(s,  * ); break;
        case BlastMatchField<>::Enum::E_VALUE:
        {
            std::string formatString;
            // imported from NCBI code
            if (match.eValue < 1.0e-180)
                formatString = "%3.1lf";
            else if (match.eValue < 1.0e-99)
                formatString = "%2.0le";
            else if (match.eValue < 0.0009)
                formatString = "%3.0le";
            else if (match.eValue < 0.1)
                formatString = "%4.3lf";
            else if (match.eValue < 1.0)
                formatString = "%3.2lf";
            else if (match.eValue < 10.0)
                formatString = "%2.1lf";
            else
                formatString = "%5.0lf";

            write(s, FormattedNumber<double>(formatString.c_str(), match.eValue));
        } break;
        case BlastMatchField<>::Enum::BIT_SCORE:
        {
            std::string formatString;
            // imported from NCBI code
            if (match.bitScore > 9999)
                formatString = "%4.3le";
            else if (match.bitScore > 99.9)
                formatString = "%4.0ld";
            else
                formatString = "%4.1lf";

            write(s, FormattedNumber<double>(formatString.c_str(), match.bitScore));
        } break;
        case BlastMatchField<>::Enum::SCORE:
            write(s, match.alignStats.alignmentScore);
            break;
        case BlastMatchField<>::Enum::LENGTH:
            write(s, match.alignStats.alignmentLength);
            break;
        case BlastMatchField<>::Enum::P_IDENT:
            write(s, FormattedNumber<float>("%.2f", match.alignStats.alignmentIdentity));
            break;
        case BlastMatchField<>::Enum::N_IDENT:
            write(s, match.alignStats.numMatches);
            break;
        case BlastMatchField<>::Enum::MISMATCH:
            if (context.legacyFormat) // legacy format includes gaps in mismatches
                write(s, match.alignStats.numMismatches +
                         match.alignStats.numGaps);
            else
                write(s, match.alignStats.numMismatches);
            break;
        case BlastMatchField<>::Enum::POSITIVE:
            write(s, match.alignStats.numPositiveScores);
            break;
        case BlastMatchField<>::Enum::GAP_OPEN:
            write(s, match.alignStats.numGapOpens);
            break;
        case BlastMatchField<>::Enum::GAPS:
            write(s, match.alignStats.numGaps);
            break;
        case BlastMatchField<>::Enum::P_POS:
            write(s, FormattedNumber<double>("%.2f", match.alignStats.alignmentSimilarity));
            break;
        case BlastMatchField<>::Enum::FRAMES:
            // for formats that don't have frames, blast says 0 instead of +1
            if (qNumFrames(context.blastProgram) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.qFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            write(s, '/');
            if (sNumFrames(context.blastProgram) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.sFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            break;
        case BlastMatchField<>::Enum::Q_FRAME:
            if (qNumFrames(context.blastProgram) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.qFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            break;
        case BlastMatchField<>::Enum::S_FRAME:
            if (sNumFrames(context.blastProgram) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.sFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            break;
//         case ENUM::BTOP: write( * ); break;
//         case ENUM::S_TAX_IDS: write( * ); break;
//         case ENUM::S_SCI_NAMES: write( * ); break;
//         case ENUM::S_COM_NAMES: write( * ); break;
//         case ENUM::S_BLAST_NAMES: write( * ); break;
//         case ENUM::S_S_KINGDOMS: write( * ); break;
//         case ENUM::S_TITLE: write( * ); break;
//         case ENUM::S_ALL_TITLES: write( * ); break;
//         case ENUM::S_STRAND: write( * ); break;
//         case ENUM::Q_COV_S: write( * ); break;
//         case ENUM::Q_COV_HSP:
        default:
            write(s, "n/i"); // not implemented
    };
}

// ----------------------------------------------------------------------------
// Function _writeFields() [match object given]
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename TScore,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeFields(TFwdIterator & stream,
             BlastIOContext<TScore, p, h> & context,
             BlastMatch<TAlign, TPos, TQId, TSId> const & match,
             BlastTabular const &)
{
    if (SEQAN_LIKELY(!context.legacyFormat))
    {
        for (auto it = std::begin(context.fields), itB = it, itEnd = std::end(context.fields); it != itEnd; ++it)
        {
            if (it != itB)
                write(stream, '\t');

            if (*it != BlastMatchField<>::Enum::STD)
            {
                _writeField(stream, context, match, *it, BlastTabular());
            }
            else // STD is placeholder for multiple fields
            {
                for (auto it2 = std::begin(BlastMatchField<>::defaults), it2B = it2,
                     it2End = std::end(BlastMatchField<>::defaults); it2 != it2End; ++it2)
                {
                    if (it2 != it2B)
                        write(stream, '\t');

                    _writeField(stream, context, match, *it2, BlastTabular());
                }
            }
        }
    }
    else
    {
        for (auto it = std::begin(BlastMatchField<>::defaults), itB = it,
             itEnd = std::end(BlastMatchField<>::defaults); it != itEnd; ++it)
        {
            if (it != itB)
                write(stream, '\t');

            _writeField(stream, context, match, *it, BlastTabular());
        }

        #if SEQAN_ENABLE_DEBUG
        if ((length(context.fields) != 1) || (context.fields[0] != BlastMatchField<>::Enum::STD))
            std::cerr << "Warning: custom fields set, but will be ignored, because legacyFormat is also set.\n";
        #endif
    }

    write(stream, '\n');
}

// ----------------------------------------------------------------------------
// Function _writeMatch()
// ----------------------------------------------------------------------------

//NOTE(h-2): dox disabled to clean-up interface
/*
 * @fn BlastTabular#_writeMatch
 * @headerfile seqan/blast.h
 * @brief write a @link BlastMatch @endlink to file
 * @signature void _writeMatch(stream, context, blastMatch, blastTabular)
 *
 * @param[in,out] stream       The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in,out] context      A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastMatch   The @link BlastMatch @endlink you wish to print.
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 *
 * @section Remarks
 *
 * Please note that BLAST is 1-indexed and considers the last position
 * to be the back, not the end, i.e. last one included in a match/sequence/...,
 * not the one behind it (as SeqAn does); this functions corrects for both of
 * these bahaviours, so you don't have to. Additionally, based on the context's
 * @link BlastProgram @endlink, positions are transformed back to DNA space, if
 * translation has taken place.
 * Please note also that query and subject IDs are truncated at the first space
 * character in NCBI BLAST, this is also done by default here.
 *
 * By setting context.@link BlastIOContext::fields @endlink you can specify which columns you
 * wish to print (if you don't want defaults); the same conversions mentioned above will me made. See
 * @link BlastMatchField::Enum @endlink for a list of fields available. If the context's legacy-flag is set
 * (@link BlastIOContext::legacyFormat @endlink) the @link BlastIOContext::fields @endlink
 * variable is ignored.
 *
 * Many guides recommend always printing the default 12 columns and using only
 * additional columns with additional (custom) data.
 *
 * Please see @link BlastTabular#_writeMatch0 @endlink for an implementation that
 * does not require a @link BlastMatch @endlink object.
 *
 * @throw IOError On low-level I/O errors.
 *
 * @see BlastRecord
 * @see BlastIOContext
 */

template <typename TQId,
          typename TSId,
          typename TFwdIterator,
          typename TScore,
          typename TPos,
          typename TAlign,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_writeMatch(TFwdIterator & stream,
           BlastIOContext<TScore, p, h> & context,
           BlastMatch<TAlign, TPos, TQId, TSId> const & match,
           BlastTabular const & /*tag*/)
{
    _writeFields(stream, context, match, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularFileOut#writeRecord
 * @headerfile seqan/blast.h
 * @brief Write a @link BlastRecord @endlink including it's @link BlastMatch @endlinkes and possible comment lines to a file.
 * @signature void writeRecord(blastTabularOut, blastRecord);
 *
 * @param[in,out] blastTabularOut A @link BlastTabularFileOut @endlink formattedFile.
 * @param[in]     blastRecord     The @link BlastRecord @endlink you wish to print.
 *
 * @section Remarks
 *
 * This function first writes the comment lines of the record and then writes match lines for every match in it. If the
 * context's @link BlastIOContext::tabularSpec @endlink is set to BlastTabularSpec::NO_COMMENTS, no
 * comment lines will be written and this function immediately prints the match lines.
 *
 * @subsection Comment lines
 *
 * The contents of the comment lines is largely defined by the members of the context.
 *
 * If the context's @link BlastIOContext::versionString @endlink is set, this will be written,
 * otherwise a stanadard version string is generated.
 *
 * The context's @link BlastIOContext::dbName @endlink is printed as the database name / path.
 *
 * If either the context's @link BlastIOContext::fields @endlink or
 * @link BlastIOContext::fieldsAsStrings @endlink member are specified these will be printed
 * as column labels. If both are specified than @link BlastIOContext::fieldsAsStrings @endlink
 * is given preference. Please note that it is recommended to use
 * @link BlastIOContext::fields @endlink and not @link BlastIOContext::fieldsAsStrings @endlink to stay
 * "standards"-compliant. Also only @link BlastIOContext::fields @endlink has an influence on
 * the values printed in the match lines.
 *
 * @subsection Matches
 *
 * Please note that BLAST is 1-indexed and considers the last position
 * to be the back, not the end, i.e. last one included in a match/sequence/...,
 * not the one behind it (as SeqAn does); this functions corrects for both of
 * these bahaviours, so you don't have to. Additionally, based on the context's
 * @link BlastIOContext::blastProgram @endlink, positions are transformed back to DNA space, if
 * translation has taken place.
 * Please note also that query and subject IDs are truncated at the first space
 * character in NCBI BLAST, this is also done by default here.
 *
 * By setting the context's @link BlastIOContext::fields @endlink meber you can specify which columns you
 * wish to print (if you don't want defaults); the same conversions mentioned above will me made. See
 * @link BlastMatchField::Enum @endlink for a list of fields available. If the context's legacy-flag is set
 * (@link BlastIOContext::legacyFormat @endlink) the @link BlastIOContext::fields @endlink
 * variable is ignored.
 *
 * Many guides recommend always printing the default 12 columns and using only
 * additional columns with additional (custom) data.
 *
 * @throw IOError On low-level I/O errors.
 *
 * @see BlastRecord
 * @see BlastIOContext
 */

template <typename TFwdIterator,
          typename TScore,
          typename TMatch,
          BlastProgram p,
          BlastTabularSpec h>
inline void
writeRecord(TFwdIterator & stream,
            BlastIOContext<TScore, p, h> & context,
            BlastRecord<TMatch> const & r,
            BlastTabular const & /*tag*/)
{
    //TODO if debug, do lots of sanity checks on record

    //NOOP for TABULAR
    _writeCommentLines(stream, context, r, BlastTabular());
    for (auto it = r.matches.begin(); it != r.matches.end(); ++it)
    {
        //SOME SANITY CHECKS
        SEQAN_ASSERT(startsWith(r.qId, it->qId));

        _writeMatch(stream, context, *it, BlastTabular());
    }
}

template <typename TContext,
          typename TMatch>
inline void
writeRecord(BlastTabularFileOut<TContext> & formattedFile,
            BlastRecord<TMatch> const & r)
{
    writeRecord(formattedFile.iter, context(formattedFile), r, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function writeHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularFileOut#writeHeader
 * @headerfile seqan/blast.h
 * @brief Write the header (top-most section) of a BlastTabular file (this is a NOOP).
 * @signature void writeHeader(blastTabularOut);
 *
 * @param[in,out] blastTabularOut A @link BlastTabularFileOut @endlink formattedFile.
 *
 * @throw IOError On low-level I/O errors.
 *
 * @see BlastRecord
 * @see BlastIOContext
 */

template <typename TFwdIterator,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
writeHeader(TFwdIterator & ,
            BlastIOContext<TScore, p, h> &,
            BlastTabular const & /*tag*/)
{
}

template <typename TContext>
inline void
writeHeader(BlastTabularFileOut<TContext> & formattedFile)
{
    writeHeader(formattedFile.iter, context(formattedFile), BlastTabular());
}

// ----------------------------------------------------------------------------
// Function writeFooter()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularFileOut#writeFooter
 * @headerfile seqan/blast.h
 * @brief Write the footer of a BlastTabular file.
 * @signature void writeFooter(blastTabularOut);
 *
 * @param[in,out] blastTabularOut A @link BlastTabularFileOut @endlink formattedFile.
 *
 * @throw IOError On low-level I/O errors.
 *
 * @see BlastRecord
 * @see BlastIOContext
 */

template <typename TFwdIterator,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
writeFooter(TFwdIterator & stream,
            BlastIOContext<TScore, p, h> & context,
            BlastTabular const & /*tag*/)
{
    if ((!context.legacyFormat) && (context.tabularSpec != BlastTabularSpec::NO_COMMENTS))
    {
        write(stream, "# BLAST processed ");
        write(stream, context.numberOfRecords); // number of records equals number of queries
        write(stream, " queries\n");
    }
}

template <typename TContext>
inline void
writeFooter(BlastTabularFileOut<TContext> & formattedFile)
{
    writeFooter(formattedFile.iter, context(formattedFile), BlastTabular());
}

} // namespace seqan
#endif // header guard
