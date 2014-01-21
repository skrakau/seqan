// ==========================================================================
//                              benchmark_calls
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

#include <cmath>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

#include "benchmark_calls.h"
#include "benchmark_calls_vs_bissnp.h"
#include "parse_meth_level.h"

struct Options
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The first (and only) argument of the program is stored here.

    CharString simSnpsFile;
    CharString simMethsFile;
    CharString calledSnpsFile;
    CharString bisSnpFile;
    CharString bedFile;
    unsigned minCovMeth;    // min coverage to benchmark methylation level
    bool methFASTA;
    int simulatedAtRefN;

    CharString outputFile;

    Options() :
        verbosity(1),
        simMethsFile(""),
        bisSnpFile(""),
        bedFile(""),
        minCovMeth(2),
        methFASTA(false),
        simulatedAtRefN(0),
        outputFile("")
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("benchmark_calls");
    // Set short description, version, and date.
    setShortDescription(parser, "Put a Short Description Here");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addOption(parser, ArgParseOption("s", "sim-snps-file", "Path to simulated snps file.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("m", "sim-meths-file", "Path to simulated meth states file.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("c", "called-snps-file", "Path to called snps file.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("bis", "bis-snp-file", "Path to BisSNP output file.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("bed", "bed-file", "Path to BED output file.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("mcm", "min-cov-meth", "Minimal coverage to benchmark meth level.", ArgParseOption::INTEGER));
    addOption(parser, ArgParseOption("o", "output-file", "Path to output file.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("mf", "meth-fasta", "Simulated methylation levels are given in FASTA format."));
    addOption(parser, ArgParseOption("sn", "sim-at-n", "Number of simulated SNPs at reference N positions (not counted as false negatives).", ArgParseOption::INTEGER));


    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBbenchmark_calls\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

        // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

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
    //getArgumentValue(options.text, parser, 0);

    getOptionValue(options.simSnpsFile, parser, "sim-snps-file");
    getOptionValue(options.simMethsFile, parser, "sim-meths-file");
    getOptionValue(options.calledSnpsFile, parser, "called-snps-file");
    getOptionValue(options.bisSnpFile, parser, "bis-snp-file");
    getOptionValue(options.bedFile, parser, "bed-file");
    getOptionValue(options.minCovMeth, parser, "min-cov-meth");
    getOptionValue(options.outputFile, parser, "output-file");
    if (isSet(parser, "meth-fasta"))
        options.methFASTA = true;
    getOptionValue(options.simulatedAtRefN, parser, "sim-at-n");

    return ArgumentParser::PARSE_OK;
}


int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "BENCHMARK CALLS\n"
              << "===============\n\n";
    
    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "sim-snps-file: \t" << options.simSnpsFile << '\n'
                  << "sim-meths-file: \t" << options.simMethsFile << '\n'
                  << "called-snps-file: \t" << options.calledSnpsFile << '\n'
                  << "bis-snp-file: \t" << options.bisSnpFile << '\n'
                  << "bed-file: \t" << options.bedFile << '\n'
                  << "min-cov-meth: \t" << options.minCovMeth << '\n'
                  << "output-file: \t" << options.outputFile << '\n'
                  << "sim-at-N: \t" << options.simulatedAtRefN << '\n';
    }
    
    if (options.bisSnpFile != "") benchmark_bisSNP(options);
    else benchmark(options, SnpMason2(), MethMason2());

    // Write out simulated vs. called meth. levels for plot
    parse_methLevel_My(options);
    if (options.bedFile != "") parse_methLevel_BisSNP(options);

    return 0;
}
