// ==========================================================================
//                                  pre_pro
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

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/arg_parse.h>
#include "pre_pro.h"

using namespace std;
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

    // The first (and only) argument of the program is stored here.
    CharString inputFileName;
    bool ctConversion;

    CharString outputFileName;

    AppOptions() :
        verbosity(1),
        ctConversion(true)
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
    ArgumentParser parser("pre_pro");
    // Set short description, version, and date.
    setShortDescription(parser, "Put a Short Description Here");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, 0, "fa fasta FASTA fastq fq FASTQ");

    addOption(parser, ArgParseOption("o", "output-file", "Name of output file.", ArgParseArgument::OUTPUTFILE));
    setValidValues(parser, "o", "fa fasta FASTA fastq fq FASTQ");
    addOption(parser, ArgParseOption("ga", "ga-conversion", "Convert Gs to As, instead of Cs to Ts."));

    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBpre_pro\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;
            
    getArgumentValue(options.inputFileName, parser, 0);
    getOptionValue(options.outputFileName, parser, "output-file");

    if (isSet(parser, "ga-conversion"))
        options.ctConversion = false;


    CharString tmp1 = options.inputFileName;
    toLower(tmp1);
    CharString tmp2 = options.outputFileName;
    toLower(tmp2);

    if ( ( (endsWith(tmp1, ".fa") || endsWith(tmp1, ".fasta")) &&  (endsWith(tmp2, ".fastq") || endsWith(tmp2, ".fq")) ) ||
         ( (endsWith(tmp2, ".fa") || endsWith(tmp2, ".fasta")) &&  (endsWith(tmp1, ".fastq") || endsWith(tmp2, ".fq")) )  )
    {
        std::cerr << "ERROR: Output file must have the same format as input file!" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

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

    std::cout << "PRE PROCESSING\n"
              << "===============\n\n";
    
    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "inputFileName     \t" << options.inputFileName << "\n\n"
                  << "outputFileName     \t" << options.outputFileName << "\n\n";
    }

    preProcess(options);

    return 0;
}

/*
#include <iostream>
#include <fstream>
#include <sstream>
#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/parallel.h>
#include <seqan/stream.h>

using namespace seqan;
using namespace std;



int main()
{
    // Readers
    std::fstream readFileStream; 
    readFileStream.open(toCString("testReadInput.txt"), std::ios::binary | std::ios_base::in);
    if(!readFileStream.is_open())
    {
        std::cerr << "Failed to open read input " << "outFile.vcf" << std::endl;
        return 1;
    }

    StringSet<CharString> nameStore;
    appendValue(nameStore, "20");
    appendValue(nameStore, "21");
    appendValue(nameStore, "22");
    appendValue(nameStore, "23");
    appendValue(nameStore, "24");
    appendValue(nameStore, "25");
    appendValue(nameStore, "26");
    appendValue(nameStore, "27");
    appendValue(nameStore, "28");
    appendValue(nameStore, "29");
    appendValue(nameStore, "30");

    NameStoreCache<StringSet<CharString> > nameStoreCache(nameStore);
    BamIOContext<StringSet<CharString>, NameStoreCache<StringSet<CharString> > > context(nameStore, nameStoreCache);

    BamHeader header;
    String<RecordReader<std::fstream, SinglePass< > >* > recordReaders;
    resize(recordReaders, 1);
    recordReaders[0] = new RecordReader<std::fstream,SinglePass< > >(readFileStream);
    readRecord(header, context, *recordReaders[0], Sam()); 
    //BamAlignmentRecord record; 

    std::cout << "Parallel: " << std::endl;

    String<CharString> contigFilenames;
    resize(contigFilenames, length(nameStore));

    // Read parallel
    // Assign task by hand to contigs
    unsigned threads = omp_get_max_threads();
    std::cout << "max threads: " << threads << std::endl; 
    unsigned chunkSize = floor((double)(length(nameStore))/(double)threads);
    unsigned rest = length(nameStore) % threads;
    SEQAN_OMP_PRAGMA(parallel for)
    for (unsigned t = 0; t < threads; ++t)   // Each thread should have it's own recordReader and parse a chunk of contigs
    {
        std::cout << "Curr. thread: " << t << std::endl; 

        std::fstream testFileStream; 
        testFileStream.open(toCString("testReadInput.txt"), std::ios::binary | std::ios_base::in);
        if(!testFileStream.is_open())
        {
            std::cerr << "Failed to open read input " << "testReadInput.txt" << std::endl;
        }

        String<RecordReader<std::fstream, SinglePass< > >* > recordReaders2;
        resize(recordReaders2, 1);
        recordReaders2[0] = new RecordReader<std::fstream,SinglePass< > >(testFileStream);
        BamHeader header2;
        BamIOContext<StringSet<CharString>, NameStoreCache<StringSet<CharString> > > context2(nameStore, nameStoreCache);
        readRecord(header2, context2, *recordReaders2[0], Sam()); 

        BamAlignmentRecord record;
        clear(record);
        for (unsigned i = 0; (i < chunkSize || ((t < rest) && i <= (chunkSize))); ++i)
        {
            // ContigId in current chunk
            unsigned id;
            if(t < rest) id = t*chunkSize + i + t;  // chunk shifted by additional rest contigs distributed to previous threads already (t)
            else id = t*chunkSize + i + rest;       // ...

            std::cout << "Curr. id: " << id << std::endl; 


            CharString currFN = "/tmp/SEQAN.XXXXXXXXXXXXXX";
            stringstream ss;
            ss << id;
            append(currFN, ss.str());
            contigFilenames[id] = currFN;
            std::cout << "tmp name: " << currFN  <<  "  id: " << id << std::endl; // "  length(recordReaders): " << length(recordReaders2) << std::endl;
            std::ofstream tempFileStream; 
            tempFileStream.open(toCString(currFN), std::ios::binary | std::ios_base::out);
            if(!tempFileStream.is_open())
            {
                std::cerr << "Failed to open read file 0" << currFN  << std::endl;
            }

            int res = 0;
            while (!atEnd(*recordReaders2[0]) && res == 0)
            {
                std::cout << "Test: " << id  << std::endl; 
                if (empty(record.qName))
                {
                    res = readRecord(record, context2, *recordReaders2[0], Sam());
                    std::cout << "Test: id" << id << "  record.rID: " << record.rID << "  record.qName: " << record.qName << std::endl; 
                }
                if (record.rID < (int)id) 
                {
                    clear(record);
                    continue;
                }
                if (record.rID > (int)id) break;
                tempFileStream << "Test: " << id << "  " << record.rID << '\t' << record.qName << std::endl;
                clear(record);
            }

            tempFileStream.close();
        }
    }
    

    std::fstream outFileStream; 
    outFileStream.open(toCString("outFile.vcf"), std::ios::binary | std::ios_base::out);
    if(!outFileStream.is_open())
    {
        std::cerr << "Failed to open read file 1 " << "outFile.vcf" << std::endl;
        return 1;
    }

    std::cout << "Append temp files: " << std::endl;

    CharString buffer;
    resize(buffer, 1000);
    for (unsigned i = 0; i < length(nameStore); ++i)
    {
        std::cout << "tmp name: " << contigFilenames[i] << std::endl;
        std::fstream tempFileStream; 
        tempFileStream.open(toCString(contigFilenames[i]), std::ios::binary | std::ios_base::in);
        if(!tempFileStream.is_open())
        {
            std::cerr << "Failed to open read file 2 " << contigFilenames[i] << std::endl;
        }

        while (!streamEof(tempFileStream) && seqan::streamError(tempFileStream) == 0)
        {
            int num = streamReadBlock(&buffer[0], tempFileStream, length(buffer));
            streamWriteBlock(outFileStream, &buffer[0], num);
        }
        remove(toCString(contigFilenames[i]));
    }
    return 0;
}
*/


