// ==========================================================================
//                                getMethLevel
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
/*
#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <string>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/refinement.h>

#include "getMethLevel.h"

using namespace seqan;

// Program entry point
int main(int argc, char const ** argv)
{
    // Setup command line parser.
    CommandLineParser parser;
    Options options;
    setupCommandLineParser(parser, options);
    
    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    int ret = parseCommandLineAndCheck(options, parser, argc, argv);
    if (ret != 0)
    {
        std::cerr << "Invalid usage!" << std::endl;
        return ret;
    }
    if (options.showHelp || options.showVersion)
        return 0;


    CharString inputFileName = options.inputFileName;
    CharString outputFileName = options.outputFileName;
 
    getMethLevel(outputFileName, inputFileName);
    
    // Finally, launch the program.
    ret = mainWithOptions(options);
    return ret;
}
*/

#include <seqan/platform.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>

#include <boost/regex.hpp>
#include <iostream>
#include <string>

//#include <boost/tuple/tuple.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/roots.hpp>


using namespace std;
using namespace seqan;

//using boost::math::tools::newton_raphson_iterate;


// Functor for extended Horner function
// Returns values for f'(x) and f''(x)
template <typename TValue>
struct Fct_bla
{   
    
    Fct_bla(String<TValue> const& coeffs) : coeffs(coeffs)
    { // Constructor 
    }
    boost::math::tuple<TValue, TValue> operator()(TValue const&beta)
    //double operator()(TValue const&beta)
    { // z is estimate so far.
        //std::cout << "start calculate horner " << std::endl;
        //std::cout << std::setprecision (25) << "Beta NSpace: " << beta << "  f': " << f_1 << "  f'': " << f_2 << std::endl;
        std::cout << "Beta: " << beta << std::endl;
        return boost::math::make_tuple(0.33, 0.66);
    }
private:
    String<TValue> coeffs;
};


int main()
{

    std::string line;
    boost::regex pat( "^Subject: (Re: |Aw: )*(.*)" );

    while (std::cin)
    {
        std::getline(std::cin, line);
        boost::smatch matches;
        if (boost::regex_match(line, matches, pat))
            std::cout << matches[2] << std::endl;
    
    }

    //typedef typename boost::tuples::tuple<int, int> TTuple;

    //boost::tuples::tuple<int, int> myTuple = boost::tuples::make_tuple(2, 3); 
    

    String<double> coeffs;
    resize(coeffs, 5, 0.0);

    Fct_bla<double> myFN(coeffs);

    boost::uintmax_t maxIter = 3;
    int digits = std::numeric_limits<double>::digits;

     boost::math::tools::newton_raphson_iterate(myFN, (double)0.5, (double)0.0, (double)1.0, digits, maxIter);

    return 0;
}

/*
#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>

int main()
{
    using namespace boost::lambda;
    typedef std::istream_iterator<int> in;

    std::for_each(
        in(std::cin), in(), std::cout << (_1 * 3) << " " );
}
*/
