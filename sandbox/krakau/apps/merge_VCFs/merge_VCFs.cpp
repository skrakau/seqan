// ==========================================================================
//                                 merge_VCFs
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

using namespace std;
using namespace seqan;

struct AppOptions
{
    int verbosity;

    CharString inputName1;
    CharString inputName2;
    CharString outputName;

    AppOptions() :
        verbosity(1),
        inputName1(""),
        inputName2(""),
        outputName("outputWithoutName.vcf")
    {}
};

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("merge_VCFs");
    // Set short description, version, and date.
    setShortDescription(parser, "Put a Short Description Here");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, 0, "vcf");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, 1, "vcf");


    addOption(parser, ArgParseOption("o", "output-file", "Name of output file.", ArgParseArgument::OUTPUTFILE));
    setValidValues(parser, "o", "vcf");

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBmerge_VCFs\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getArgumentValue(options.inputName1, parser, 0);
    getArgumentValue(options.inputName2, parser, 1);
    getOptionValue(options.outputName, parser, "output-file");

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    std::cout << "EXAMPLE PROGRAM\n"
              << "===============\n\n";

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n';
    }

    // Open input stream
    VcfStream vcfIn1(toCString(options.inputName1));
    VcfStream vcfIn2(toCString(options.inputName2));

    if (!isGood(vcfIn1) || !isGood(vcfIn2))
    {
        std::cerr << "ERROR: Could not open vcf files \n";
        return 1;
    }
    // Open output stream, filename "-" means stdout.
    VcfStream vcfOut(toCString(options.outputName), VcfStream::WRITE);
    // Copy over header.
    vcfOut.header = vcfIn1.header;
    // Read the file record by record.
    VcfRecord record;
    VcfRecord record1;
    VcfRecord record2;
    __int32 prev_rID = -1;
    if (readRecord(record1, vcfIn1) != 0 || readRecord(record2, vcfIn2) != 0 )
    {
        std::cerr << "ERROR: Problem reading first entry from vcf file\n";
        return 1;
    }
    while (!atEnd(vcfIn1) || !atEnd(vcfIn2) )
    {
        clear(record);
/*
CHROM -> rID
    Name of the chromosome/reference sequence that the variant lies on. 
POS -> beginPos
    The 1-based position of the variant. 
ID -> id
    A name of the variant. . is used if no name is available. 
REF
    The value of the reference allele. 
ALT
    The alternate allele values (multiple values are comma-separated). 
QUAL
    Quality value of the call (float). 
FILTER
    A value for the filter result (given in a FILTER meta information line). 
INFO
    Information about a variant. 
FORMAT
    Colon-separated list of entries that are found for each variant. 
*/
        
        if (record1.rID == record2.rID)
        {
            record.rID = record1.rID;
            if (record1.beginPos == record2.beginPos) // SNP in both haplotypes
            {
                record.beginPos = record1.beginPos;
                if (record1.alt == record2.alt)
                {
                    record.ref = record1.ref;
                    record.alt = record1.alt;
                    clear(record1.genotypeInfos);
                    appendValue(record1.genotypeInfos, "1|1");
                }
                else if (record1.alt < record2.alt)
                {
                    record.ref = record1.ref;
                    record.alt = record1.alt;
                    appendValue(record.alt, ',');
                    append(record.alt, record2.alt);
                    appendValue(record.genotypeInfos, "1|2");
                }
                else
                {
                    record.ref = record2.ref;
                    record.alt= record2.alt;
                    appendValue(record.alt, ',');
                    append(record.alt, record1.alt);
                    appendValue(record.genotypeInfos, "1|2");
                }
                if (!atEnd(vcfIn1) && !atEnd(vcfIn2))
                {
                    if (readRecord(record1, vcfIn1) != 0 || readRecord(record2, vcfIn2) != 0 ) 
                    {
                        std::cerr << "ERROR: Problem reading from 1. or 2. vcf file\n";
                        return 1;
                    }
                }
                else if (!atEnd(vcfIn1) && readRecord(record1, vcfIn1) != 0 )
                {
                    std::cerr << "ERROR: Problem reading from 1. vcf file\n";
                    return 1;
                }
                else if (!atEnd(vcfIn2) && readRecord(record2, vcfIn2) != 0 )
                {
                    std::cerr << "ERROR: Problem reading from 2. vcf file\n";
                    return 1;
                }
                if (atEnd(vcfIn1)) record1.rID = -1;
                if (atEnd(vcfIn2)) record2.rID = -1;
            }
            else if (record1.beginPos < record2.beginPos)   // SNP only in haplotype1
            {
                record.beginPos = record1.beginPos;
                record.ref = record1.ref;
                record.alt = record1.alt;
                appendValue(record.genotypeInfos, "1|0");
                if (!atEnd(vcfIn1) && readRecord(record1, vcfIn1) != 0)
                {
                    std::cerr << "ERROR: Problem reading from 1. vcf file\n";
                    return 1;
                }
                else if (atEnd(vcfIn1))  record1.rID = -1;
            }
            else                                            // SNP only in haplotype2
            {
                record.beginPos = record2.beginPos;
                record.ref = record2.ref;
                record.alt = record2.alt;
                appendValue(record.genotypeInfos, "0|1");
                if (!atEnd(vcfIn2) && readRecord(record2, vcfIn2) != 0)
                {
                    std::cerr << "ERROR: Problem reading from 2. vcf file\n";
                    return 1;
                }
                else if (atEnd(vcfIn2)) record2.rID = -1;
            }
        }
        else if (record1.rID == prev_rID)                   // SNP only in haplotype1, reached already next contig in hapltoype2
        {
            record.rID = record1.rID;
            record.beginPos = record1.beginPos;
            record.ref = record1.ref;
            record.alt = record1.alt;
            appendValue(record.genotypeInfos, "1|0");
            if (!atEnd(vcfIn1) && readRecord(record1, vcfIn1) != 0)
            {
                std::cerr << "ERROR: Problem reading from 1. vcf file\n";
                return 1;
            }
            else if (atEnd(vcfIn1)) record1.rID = -1;
        }
        else if (record2.rID == prev_rID)                 // SNP only in haplotype2, reached already next contig in hapltoype1
        {
            record.rID = record2.rID;
            record.beginPos = record2.beginPos;
            record.ref = record2.ref;
            record.alt = record2.alt;
            appendValue(record.genotypeInfos, "0|1");
            if (!atEnd(vcfIn2) && readRecord(record2, vcfIn2) != 0)
            {
                std::cerr << "ERROR: Problem reading from 2. vcf file\n";
                return 1;
            }
            else if (atEnd(vcfIn2)) record2.rID = -1;
        }

        // Write
        --record.beginPos;  // TODO did it now, since other modules are 0-based still, change to official stupid style !!!!
        if (writeRecord(vcfOut, record) != 0)
        {
            std::cerr << "ERROR: Problem writing to output.\n";
            return 1;
        }
        prev_rID = record.rID;
    }

    return 0;
}



