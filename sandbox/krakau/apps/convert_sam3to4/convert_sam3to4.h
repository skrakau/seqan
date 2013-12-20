#ifndef SANDBOX_KRAKAU_APPS_CONVERT_SAM3TO4_CONVERT_SAM3TO4_H_
#define SANDBOX_KRAKAU_APPS_CONVERT_SAM3TO4_CONVERT_SAM3TO4_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <fstream>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/bam_io.h>
#include <seqan/score.h> 

using namespace seqan;

template <typename TFSSpec, typename TFSConfig, typename TFileName>
bool loadReadsCroppedId(FragmentStore<TFSSpec, TFSConfig> &store, TFileName &fileName)
{
	MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(fileName), OPEN_RDONLY))
		return false;

	// guess file format and split into sequence fractions
	AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);

	// reserve space in fragment store
	unsigned seqOfs = length(store.readStore);
	unsigned seqCount = length(multiSeqFile);
	reserve(store.readStore, seqOfs + seqCount);
	reserve(store.readSeqStore, seqOfs + seqCount);
	reserve(store.readNameStore, seqOfs + seqCount);

	// read sequences
	String<Dna5Q> seq;
	CharString qual;
	CharString _id;
	for (unsigned i = 0; i < seqCount; ++i)
	{
		assignSeq(seq, multiSeqFile[i], format);    // read sequence
		assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
        assignCroppedSeqId(_id, multiSeqFile[i], format);  // read sequence id up to the first whitespace
        
		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
		// TODO: support different ASCII represenations of quality values
		assignQualities(seq, qual);
		appendRead(store, seq, _id);
	}
    return true;
}

template <typename TFSSpec, typename TFSConfig, typename TReadInfos, typename TFileName>
bool loadReadsCroppedId(FragmentStore<TFSSpec, TFSConfig> &store, TReadInfos &readInfos, TFileName &fileName)
{
	MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(fileName), OPEN_RDONLY))
		return false;

	// guess file format and split into sequence fractions
	AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);

	// reserve space in fragment store
	unsigned seqOfs = length(store.readStore);
	unsigned seqCount = length(multiSeqFile);
	reserve(store.readStore, seqOfs + seqCount);
	reserve(store.readSeqStore, seqOfs + seqCount);
	reserve(store.readNameStore, seqOfs + seqCount);

	// read sequences
	String<Dna5Q> seq;
	CharString qual;
	CharString _id;
	CharString id_long;

	for (unsigned i = 0; i < seqCount; ++i)
	{
		assignSeq(seq, multiSeqFile[i], format);    // read sequence
		assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
        assignCroppedSeqId(_id, multiSeqFile[i], format);  // read sequence id up to the first whitespace
        assignSeqId(id_long, multiSeqFile[i], format);  // long read sequence id 
        
		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
		// TODO: support different ASCII represenations of quality values
		assignQualities(seq, qual);
		appendRead(store, seq, _id);
		appendValue(readInfos, suffix(id_long, length(id_long)-1));
	}
    return true;
}

template <typename TFSSpec, typename TFSConfig, typename TFileName>
bool loadReadsCroppedId(FragmentStore<TFSSpec, TFSConfig> & store, TFileName & fileNameL, TFileName & fileNameR)
{
	MultiSeqFile multiSeqFileL, multiSeqFileR;
	if (!open(multiSeqFileL.concat, toCString(fileNameL), OPEN_RDONLY))
		return false;
	if (!open(multiSeqFileR.concat, toCString(fileNameR), OPEN_RDONLY))
		return false;

	// Guess file format and split into sequence fractions
	AutoSeqFormat formatL, formatR;
	guessFormat(multiSeqFileL.concat, formatL);
	split(multiSeqFileL, formatL);
	guessFormat(multiSeqFileR.concat, formatR);
	split(multiSeqFileR, formatR);

    // Check that both files have the same number of reads
	SEQAN_ASSERT_EQ(length(multiSeqFileL), length(multiSeqFileR));

	// Reserve space in fragment store
	unsigned seqOfs = length(store.readStore);
	unsigned seqCountL = length(multiSeqFileL);
	unsigned seqCountR = length(multiSeqFileR);
	reserve(store.readStore, seqOfs + seqCountL + seqCountR);
	reserve(store.readSeqStore, seqOfs + seqCountL + seqCountR);
	reserve(store.readNameStore, seqOfs + seqCountL + seqCountR);

	// Read in sequences
	String<Dna5Q> seq[2];
	CharString qual[2];
	CharString _id[2];
	for (unsigned i = 0; i < seqCountL; ++i) {
		assignSeq(seq[0], multiSeqFileL[i], formatL);    // read sequence
		assignQual(qual[0], multiSeqFileL[i], formatL);  // read ascii quality values
        assignCroppedSeqId(_id[0], multiSeqFileL[i], formatL);  // read sequence id up to the first whitespace 

		assignSeq(seq[1], multiSeqFileR[i], formatR);    // read sequence
		assignQual(qual[1], multiSeqFileR[i], formatR);  // read ascii quality values
		assignCroppedSeqId(_id[1], multiSeqFileR[i], formatR);  // read sequence id up to the first whitespace

		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
		// TODO: support different ASCII represenations of quality values
		for (int j = 0; j < 2; ++j)
			assignQualities(seq[j], qual[j]);
		
		appendMatePair(store, seq[0], seq[1], _id[0], _id[1]);
	}
	return true;
}


template <typename TFSSpec, typename TFSConfig, typename TReadInfos, typename TFileName>
bool loadReadsCroppedId(FragmentStore<TFSSpec, TFSConfig> & store, TReadInfos &readInfos, TFileName & fileNameL, TFileName & fileNameR)
{
	MultiSeqFile multiSeqFileL, multiSeqFileR;
	if (!open(multiSeqFileL.concat, toCString(fileNameL), OPEN_RDONLY))
		return false;
	if (!open(multiSeqFileR.concat, toCString(fileNameR), OPEN_RDONLY))
		return false;

	// Guess file format and split into sequence fractions
	AutoSeqFormat formatL, formatR;
	guessFormat(multiSeqFileL.concat, formatL);
	split(multiSeqFileL, formatL);
	guessFormat(multiSeqFileR.concat, formatR);
	split(multiSeqFileR, formatR);

    // Check that both files have the same number of reads
	SEQAN_ASSERT_EQ(length(multiSeqFileL), length(multiSeqFileR));

	// Reserve space in fragment store
	unsigned seqOfs = length(store.readStore);
	unsigned seqCountL = length(multiSeqFileL);
	unsigned seqCountR = length(multiSeqFileR);
	reserve(store.readStore, seqOfs + seqCountL + seqCountR);
	reserve(store.readSeqStore, seqOfs + seqCountL + seqCountR);
	reserve(store.readNameStore, seqOfs + seqCountL + seqCountR);

	// Read in sequences
	String<Dna5Q> seq[2];
	CharString qual[2];
	CharString _id[2];
	CharString id_long[2];

	for (unsigned i = 0; i < seqCountL; ++i) {
		assignSeq(seq[0], multiSeqFileL[i], formatL);    // read sequence
		assignQual(qual[0], multiSeqFileL[i], formatL);  // read ascii quality values
        assignCroppedSeqId(_id[0], multiSeqFileL[i], formatL);  // read sequence id up to the first whitespace 
        assignSeqId(id_long[0], multiSeqFileL[i], formatR);  // long read sequence id 
		appendValue(readInfos, suffix(id_long[0], length(id_long[0])-1));

		assignSeq(seq[1], multiSeqFileR[i], formatR);    // read sequence
		assignQual(qual[1], multiSeqFileR[i], formatR);  // read ascii quality values
		assignCroppedSeqId(_id[1], multiSeqFileR[i], formatR);  // read sequence id up to the first whitespace
        assignSeqId(id_long[1], multiSeqFileR[i], formatR);  // long read sequence id 
		appendValue(readInfos, suffix(id_long[1], length(id_long[1])-1));

		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
		// TODO: support different ASCII represenations of quality values
		for (int j = 0; j < 2; ++j)
			assignQualities(seq[j], qual[j]);
		
		appendMatePair(store, seq[0], seq[1], _id[0], _id[1]);
	}
	return true;
}

template<typename TOptions>
bool
convertSam(TOptions &options)
{
#if SEQAN_ENABLE_PARALLELISM
    std::cout << "NOTE: parellelism is enabled. Set number of threads to 1 for this app. " << std::endl;
#endif 

    typedef FragmentStore<>                                                         TFragmentStore;
    typedef typename Value<typename TFragmentStore::TMatePairStore>::Type           TMatePairStoreElement;
    typedef typename TFragmentStore::TReadNameStore                                 TReadNameStore;
    typedef typename TFragmentStore::TContigNameStore                               TContigNameStore;
    typedef NameStoreCache<TReadNameStore, CharString>                              TReadNameStoreCache;
    typedef NameStoreCache<TContigNameStore, CharString>                            TContigNameStoreCache;
    typedef BamIOContext<TContigNameStore, TContigNameStoreCache>                   TBamIOContext;
    
    // Load all reads 
    TFragmentStore store;
    String<CharString> readInfos; 
    // Load original reads into fragmentStore
    if (options.readFileName2 == "")
    {
        if (options.includeStrandInfo) loadReadsCroppedId(store, readInfos, options.readFileName);
        else loadReadsCroppedId(store, options.readFileName);
    }
    else
    {
        if (options.includeStrandInfo) loadReadsCroppedId(store, readInfos, options.readFileName, options.readFileName2);
        else loadReadsCroppedId(store, options.readFileName, options.readFileName2);
    }
    
    TReadNameStoreCache readNameCache(store.readNameStore); // 
    refresh(readNameCache);
    std::cout << "No. of reads in readStore: " << length(store.readStore) << std::endl;
  
    // Parse SAM file ....
    std::fstream inStream(toCString(options.samFileName), std::ios::binary | std::ios::in);
    if (!inStream.good())
    {
        std::cerr << "ERROR: Could not open " << toCString(options.samFileName) << " for reading.\n";
        return 1;
    }
    RecordReader<std::fstream, seqan::SinglePass<> > reader(inStream);
    std::fstream outStream(toCString(options.outputFileName), std::ios::binary | std::ios::out);
    if (!outStream.good())
    {
        std::cerr << "ERROR: Could not open " << toCString(options.outputFileName) << " for writing.\n";
        return 1;
    }

    BamHeader header;
    BamAlignmentRecord record;
    TBamIOContext bamIOContext(store.contigNameStore, store.contigNameStoreCache);
    // Read header
    if (readRecord(header, bamIOContext, reader, Sam()) != 0)  
    {
        std::cerr << "ERROR: Could not read SAM header record!\n";
        return 1;
    }  
    // Write out header again. Maybe add information about BS mapping ?
    if (write2(outStream, header, bamIOContext, Sam()) != 0)
    {
        std::cerr << "ERROR: Could not write header to SAM file "  << "\n";
        return 1;
    }
    std::cout << "Header writen. " << std::endl;
    unsigned count = 0;
    while (!atEnd(reader))
    {
        ++count;
        if (readRecord(record, bamIOContext, reader, Sam()) != 0)  
        {
            std::cerr << "ERROR: Could not read sam record!\n";
            return 1;
        }     

        unsigned readId;
        if (!getIdByName(store.readNameStore, record.qName, readId, readNameCache)) {
            std::cout << "NOTE: Read not found in nameStore: " << record.qName << std::endl;
            continue;
        }
        /*if ((record.flag & 0x1) == 1)    // If paired: Get readId for current mate
        {
            int inPair = 1 - ((record.flag & 0x40) >> 6);
            unsigned matePairId = store.readStore[readId].matePairId;
            if (matePairId != TMatePairStoreElement::INVALID_ID)
            {    
                readId = store.matePairStore[matePairId].readId[inPair];
                if (readId == TMatePairStoreElement::INVALID_ID) continue; 
            }
        }*/

        record.seq = store.readSeqStore[readId];
        if (/*!options.includeStrandInfo &&*/ hasFlagRC(record))
        {
            reverseComplement(record.seq);
        }
        /*else if (options.includeStrandInfo)
        {
            if(readInfos[readId] == "R")
            {
                record.flag = 0x0010;
                if(hasFlagRC(record)) reverseComplement(record.seq);
            }
            else
            {
                record.flag = 0x0000;
                if(hasFlagRC(record)) reverseComplement(record.seq);
           }
        }*/

        if (write2(outStream, record, bamIOContext, Sam()) != 0)
        {
     std::cout << "No. of reads in readStore: " << length(store.readStore) << std::endl;
           return 1;
        }
        if (count % 500000 == 0 ) std::cout << "Parsed " << count << " reads." << std::endl;
    }


    return 0;
}


#endif
