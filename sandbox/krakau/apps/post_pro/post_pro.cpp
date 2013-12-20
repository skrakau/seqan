// ==========================================================================
//                                  post_pro
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

#define POST_PRO_PROFILE

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/bam_io.h>
#include <seqan/score.h> 

#include "bs_score_data.h"
#include "bs_score.h"
#include "post_pro_base.h"
#include "post_pro.h"


//using namespace std;
using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    CharString readFileName;
    CharString readFileName2;
    CharString samFileName;
    CharString refFileName;
    CharString outputFileName;
    unsigned intervalOffset;
    double minMapq;
    unsigned max4Error;        // max. allowed real error rate
    unsigned max3Error;        // max. allowed error rate in 3 letter alphabet (corresponding to mapper settings)
    double maxScore;
    unsigned maxBasePenalty;    // limit the penalty for a single base 

    int minScore;
    bool outputSingleMates;

    double scoreMatch;
    double scoreMismatch;

    bool simpleScore;
    bool nonSimpleSubstErrors;
    bool nonSimpleInsErrors;
    bool nonSimpleDelErrors;
    double delErrorRate;

    double lambda;
    double gapOpenScore;
    double gapExtendScore;
    double scalingFactorDelErrors;
    double scalingFactorInsErrors;

    double bsConversionRate;
    double globalMethRate;

    AppOptions() :
        verbosity(1),
        intervalOffset(3),
        minMapq(10),
        max4Error(4),
        max3Error(3),
        maxScore(1000000),  // TODO: what would be reasonable?
        maxBasePenalty(-3),  // scaled to single penalties
        minScore(0),
        outputSingleMates(true),    // Output also read whose mate didn't map & if no match mate pair found, output mates single 
        scoreMatch(10.0),
        scoreMismatch(0.1),
        simpleScore(true),
        nonSimpleSubstErrors(false),
        nonSimpleInsErrors(false),
        nonSimpleDelErrors(false),
        delErrorRate(0.001),
        lambda(1.0),
        gapOpenScore(-4.0),
        gapExtendScore(-1.0),
        scalingFactorDelErrors(5.0),
        scalingFactorInsErrors(5.0),
        bsConversionRate(0.98),
        globalMethRate(0.2)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("post_pro");
    // Set short description, version, and date.
    setShortDescription(parser, "Put a Short Description Here");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    // We require ... arguments.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, 0, "sam");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, 1, "fa fasta FA FASTA");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "IN", true));
    setValidValues(parser, 2, "fastq fq FASTQ FQ");

    addOption(parser, ArgParseOption("o", "output-file", "SAM output file.", ArgParseArgument::OUTPUTFILE));
    setValidValues(parser, "o", "sam");

    addOption(parser, ArgParseOption("e3", "max3-error", "Max. error rate in 3-letter alphabet.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("e4", "max4-error", "Max. error rate in 4-letter alphabet.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("mq", "min-mapq", "Min required mapping quality.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("ns", "non-simple", "Use non-simple scoring method for all distributions (Non unifrom distributions)."));
    addOption(parser, ArgParseOption("nse", "ns-subst-errors", "Use nonuniform substitution error frequency."));
    addOption(parser, ArgParseOption("nsi", "ns-ins-errors", "Use nonuniform insertion error frequency."));
    addOption(parser, ArgParseOption("nsd", "ns-del-errors", "Use nonuniform deletion error frequency."));
    addOption(parser, ArgParseOption("der", "del-error-rate", "Deletion error rate.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("gas", "gap-open-score", "Gap open score (Original, must be proportional to mismatch scores) Dafault: -4.0.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("ges", "gap-extend-score", "Gap extend score.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("bsc", "bs-conversion-rate", "Bs conversion rate.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("gmr", "global-meth-rate", "Global methylation rate for background frequencies.", ArgParseArgument::DOUBLE));
 
    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBpost_pro\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.samFileName, parser, 0);
    getArgumentValue(options.refFileName, parser, 1);
 
    if (1 == getArgumentValueCount(parser, 2))
        getArgumentValue(options.readFileName, parser, 2, 0);
    else if (2 == getArgumentValueCount(parser, 2))
    {
        getArgumentValue(options.readFileName, parser, 2, 0);
        getArgumentValue(options.readFileName2, parser, 2, 1);  
    }

    getOptionValue(options.outputFileName, parser, "output-file");
    getOptionValue(options.max3Error, parser, "max3-error");
    getOptionValue(options.max4Error, parser, "max4-error");
    getOptionValue(options.minMapq, parser, "min-mapq");
    if (isSet(parser, "non-simple"))
        options.simpleScore = false;
    if (isSet(parser, "ns-subst-errors"))
        options.nonSimpleSubstErrors = true;
    if (isSet(parser, "ns-ins-errors"))
        options.nonSimpleInsErrors = true;
    if (isSet(parser, "ns-del-errors"))
        options.nonSimpleDelErrors = true;
    getOptionValue(options.delErrorRate, parser, "del-error-rate");
    getOptionValue(options.gapOpenScore, parser, "gap-open-score");
    getOptionValue(options.gapExtendScore, parser, "gap-extend-score");
    getOptionValue(options.bsConversionRate, parser, "bs-conversion-rate");
    getOptionValue(options.globalMethRate, parser, "global-meth-rate");

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    return ArgumentParser::PARSE_OK;
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

    std::cout << "POST PROCESSING\n"
              << "===============\n\n";
    
    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "readfileName\t" << options.readFileName << '\n'
                  << "readfileName\t" << options.readFileName2 << '\n'
                  << "samFileName\t" << options.samFileName << '\n'
                  << "refFileName     \t" << options.refFileName << '\n'
                  << "outputFileName     \t" << options.outputFileName << '\n'
                  << "\n";
    }

#ifdef POST_PRO_PROFILE 
    double timeStamp = sysTime();
#endif
    if (!options.simpleScore)
        postProcessMain(options, BsNonSimple());
    else
        postProcessMain(options, BsSimple());

#ifdef POST_PRO_PROFILE
    Times::instance().time_all = sysTime() - timeStamp;
    std::cout << "  Time needed for all: " << Times::instance().time_all/60.0 << "min" << std::endl;
    std::cout << "  Time needed for globalAlignment: " << Times::instance().time_globalAlignment/60.0 << "min" << std::endl;
    std::cout << "  Time needed for writeBsAlignment: " << Times::instance().time_writeBsAlignment/60.0 << "min" << std::endl;
#endif

    return 0;
}
