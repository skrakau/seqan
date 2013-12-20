/*==========================================================================

Get methylation stats for each annotation

- for each context CG, CHG, CHH:

    - methylated
    - fuzzy methylated
    - unmethylated

Counts only for reference Cs (Gs).

==========================================================================*/

#ifndef SANDBOX_KRAKAU_APPS_ANNOMETHS_H_
#define SANDBOX_KRAKAU_APPS_ANNOMETHS_H_

using namespace seqan;

struct AnnoMethStats
{
    unsigned methsCG;
    unsigned methsCHG;
    unsigned methsCHH;

    unsigned fuzzyCG;
    unsigned fuzzyCHG;
    unsigned fuzzyCHH;

    unsigned unmethsCG;
    unsigned unmethsCHG;
    unsigned unmethsCHH;

    AnnoMethStats() 
    {
        methsCG = 0;
        methsCHG = 0;
        methsCHH = 0;

        fuzzyCG = 0;
        fuzzyCHG = 0;
        fuzzyCHH = 0;

        unmethsCG = 0;
        unmethsCHG = 0;
        unmethsCHH = 0;
    }
};

template<typename TInterval, typename TStore>
inline 
void 
extractGeneIntervals(String<String<TInterval> > &intervalsF, String<String<TInterval> > &intervalsR, TStore const & store)
{
    typedef typename TStore::TAnnotationStore           TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotationStoreElement;
    typedef typename TAnnotationStoreElement::TId       TId;
    typedef typename TAnnotationStoreElement::TPos      TPos;
 	static const TPos INVALID_POS = TAnnotationStoreElement::INVALID_POS; 

    // extract intervals from gene annotations (grouped by contigId)
    resize(intervalsF, length(store.contigStore));
    resize(intervalsR, length(store.contigStore));
    Iterator<FragmentStore<> const, AnnotationTree<> >::Type it = begin(store, AnnotationTree<>());

    while (!atEnd(it))
    {
        if (getAnnotation(it).beginPos != INVALID_POS)
        {
            TPos beginPos = getAnnotation(it).beginPos;
            TPos endPos = getAnnotation(it).endPos;
            TId contigId = getAnnotation(it).contigId;
            if (beginPos < endPos)
                appendValue(intervalsF[contigId], TInterval(beginPos, endPos, value(it)));
            else 
            {
                std::swap(beginPos, endPos);
                appendValue(intervalsR[contigId], TInterval(beginPos, endPos, value(it)));
            }
        }
        goNext(it);
    }
}

template<typename TIntervalTree, typename TInterval>
inline void 
constructIntervalTrees(String<TIntervalTree> & intervalTreesF,
                        String<TIntervalTree> & intervalTreesR,
                        String<String<TInterval> > & intervalsF,
                        String<String<TInterval> > & intervalsR)
{
    int numContigs = length(intervalsF);
    resize(intervalTreesF, numContigs);
    resize(intervalTreesR, numContigs);

    SEQAN_OMP_PRAGMA(parallel for)
    for (int i = 0; i < numContigs; ++i)
    {
        createIntervalTree(intervalTreesF[i], intervalsF[i]);
        createIntervalTree(intervalTreesR[i], intervalsR[i]);
    }
}

// Counts also Snp Cs, independent if only one haplotype is C
template<typename TIntervalTree, typename TStore, typename TOptions>
inline void 
parseAndCountMeths(String<AnnoMethStats> &stringAnnoStats, String<AnnoMethStats> &stringTypeStats, 
                   String<TIntervalTree> const & intervalTreesF, String<TIntervalTree> const & intervalTreesR, 
                   TStore const & store, 
                   TOptions &options)
{
    typedef typename TStore::TAnnotationStore           TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotationStoreElement;
    typedef typename TAnnotationStoreElement::TId       TId;
    typedef typename TAnnotationStoreElement::TPos      TPos;
    typedef Stream<std::fstream>                        TStream;
    typedef RecordReader<std::fstream, SinglePass<> >   TRecordReader;

    std::fstream methsFile(toCString(options.methsFileName), std::ios::binary | std::ios::in);
    TRecordReader methsReader(methsFile);

    resize(stringAnnoStats, length(store.annotationStore));
    if (options.showTypeStats) resize(stringTypeStats, length(store.annotationTypeStore));

    while(!atEnd(methsReader))
    {
        CharString contigName;
        TPos pos;
        CharString refAllele;
        double methLevel;
        CharString helper;
        bool topStrand;

        readUntilWhitespace(contigName, methsReader);           // chrom
        if (contigName[0] == '#') 
        {
            skipLine(methsReader);
            continue;
        }
        skipWhitespaces(methsReader);
        clear(helper);
        readUntilWhitespace(helper, methsReader);           // pos
        pos = lexicalCast<size_t>(helper);
        skipWhitespaces(methsReader);
        readUntilWhitespace(refAllele, methsReader);        // refAllele
        if (refAllele == 'C') topStrand = true;
        else if (refAllele == 'G') topStrand = false;
        else 
        {
            skipLine(methsReader);
            continue;
        }
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // qA
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // qC
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // qG
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // qT
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // qA
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // qC
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // qG
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // qT
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // cov
        skipWhitespaces(methsReader);
        readUntilWhitespace(helper, methsReader);           // genotype called
        skipWhitespaces(methsReader);
        clear(helper);
        readUntilWhitespace(helper, methsReader);           // MethLevel
        if (helper == '.')
        {
            skipLine(methsReader);
            continue;
        }
        methLevel = lexicalCast<double>(helper);
        // Get corresponding contigId
        TId contigId;
        getIdByName(store.contigNameStore, contigName, contigId, store.contigNameStoreCache);
        // Get context for both strands CG:0, CHG:1, CHH:2
        unsigned context = 2;
        if (options.refFileName!="")
        {
            if (topStrand)
            {
                if (pos < (int)length(store.contigStore[contigId].seq)-1 && store.contigStore[contigId].seq[pos+1] == 'G')
                    context = 0;
                else if (pos < (int)length(store.contigStore[contigId].seq)-2 && store.contigStore[contigId].seq[pos+2] == 'G')
                    context = 1;
            }
            else
            {
                // Check if not first position
                if (pos > 0 && store.contigStore[contigId].seq[pos-1] == 'C')
                    context = 0;
                else if (pos > 1 && store.contigStore[contigId].seq[pos-2] == 'C')
                    context = 1;
            }
        }

        // Get annotation ids which contain the current position 
        // (could be much more efficient, but forget about it for the moment)
        String<TId> result;
        if (topStrand) findIntervals(intervalTreesF[contigId], pos, result);
        else findIntervals(intervalTreesR[contigId], pos, result);
        // For each anno id: increase corresponding counts
        for (unsigned j = 0; j < length(result); ++j)
        {
            TId annoId = result[j];
            TId typeId = store.annotationStore[annoId].typeId;

            if (context == 0)
            {
                if (methLevel < options.thresholdUnmeth)
                {
                    ++stringAnnoStats[annoId].unmethsCG;
                    if (options.showTypeStats) ++stringTypeStats[typeId].unmethsCG;
                }
                else if (methLevel > options.thresholdMeth)
                {
                    ++stringAnnoStats[annoId].methsCG;
                    if (options.showTypeStats) ++stringTypeStats[typeId].methsCG;
                }
                else
                {
                    ++stringAnnoStats[annoId].fuzzyCG;
                    if (options.showTypeStats) ++stringTypeStats[typeId].fuzzyCG;
                }            }
            else if (context == 1)
            {
                if (methLevel < options.thresholdUnmeth) 
                {
                    ++stringAnnoStats[annoId].unmethsCHG;
                    if (options.showTypeStats) ++stringTypeStats[typeId].unmethsCHG;
                }
                else if (methLevel > options.thresholdMeth)
                {
                    ++stringAnnoStats[annoId].methsCHG;
                    if (options.showTypeStats) ++stringTypeStats[typeId].methsCHG;
                }
                else 
                {
                    ++stringAnnoStats[annoId].fuzzyCHG;
                    if (options.showTypeStats) ++stringTypeStats[typeId].fuzzyCHG;
                }
            }
            else
            {
                if (methLevel < options.thresholdUnmeth)
                {
                    ++stringAnnoStats[annoId].unmethsCHH;
                    if (options.showTypeStats) ++stringTypeStats[typeId].unmethsCHH;
                }
                else if (methLevel > options.thresholdMeth)
                {
                    ++stringAnnoStats[annoId].methsCHH;
                    if (options.showTypeStats) ++stringTypeStats[typeId].methsCHH;
                }
                else
                {
                    ++stringAnnoStats[annoId].fuzzyCHH;
                    if (options.showTypeStats) ++stringTypeStats[typeId].fuzzyCHH;
                }
            }
        }
        skipLine(methsReader); 
    }
    /*
    std::cout << " unmethCG " << stringAnnoStats[3].unmethsCG << std::endl;
    std::cout << " unmethCHG " << stringAnnoStats[3].unmethsCHG << std::endl;
    std::cout << " unmethCHH " << stringAnnoStats[3].unmethsCHH << std::endl;
    */
}

// Assign stats to annotation key-value pairs
template<typename TStore, typename TOptions>
inline void
assignToKeyValuePairs(TStore &store, String<AnnoMethStats> &stringAnnoStats, TOptions &options)
{
    // Create iterator
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it  = begin(store, AnnotationTree<>());
    // Iterate over annotation tree and count different elements and compute exon lengths
    while (!atEnd(it)){
        
        // For parents: count too if possible
        if (options.refFileName!="")
        {
            // CG
            // methylated
            unsigned countCGs = stringAnnoStats[value(it)].methsCG + stringAnnoStats[value(it)].fuzzyCG + stringAnnoStats[value(it)].unmethsCG;
            /*if (value(it) == 6)
            {
                std::cout << "countCGs:" << countCGs << "  " << stringAnnoStats[value(it)].methsCG << " " << stringAnnoStats[value(it)].fuzzyCG << " " << stringAnnoStats[value(it)].unmethsCG << std::endl;
            }*/
            if (countCGs > 0)
            {
                double methRateCG = (double)stringAnnoStats[value(it)].methsCG/(double)countCGs;
                std::stringstream tmp_methRateCG;
                tmp_methRateCG << methRateCG;
                assignValueByKey(it, "meth_rate_cg", tmp_methRateCG.str());
                // fuzzy methylated
                double fuzzyRateCG = (double)stringAnnoStats[value(it)].fuzzyCG/(double)countCGs;
                std::stringstream tmp_fuzzyRateCG;
                tmp_fuzzyRateCG << fuzzyRateCG;
                assignValueByKey(it, "fuzzy_rate_cg", tmp_fuzzyRateCG.str());
            }
            else
            {
                assignValueByKey(it, "meth_rate_cg", "-1");
                assignValueByKey(it, "fuzzy_rate_cg", "-1");
            }

            // CHH
            // methylated
            unsigned countCHGs = stringAnnoStats[value(it)].methsCHG + stringAnnoStats[value(it)].fuzzyCHG + stringAnnoStats[value(it)].unmethsCHG;

            if (countCHGs > 0)
            {
                double methRateCHG = (double)stringAnnoStats[value(it)].methsCHG/(double)countCHGs;
                std::stringstream tmp_methRateCHG;
                tmp_methRateCHG << methRateCHG;
                assignValueByKey(it, "meth_rate_chg", tmp_methRateCHG.str());
                // fuzzy methylated
                double fuzzyRateCHG = (double)stringAnnoStats[value(it)].fuzzyCHG/(double)countCHGs;
                std::stringstream tmp_fuzzyRateCHG;
                tmp_fuzzyRateCHG << fuzzyRateCHG;
                assignValueByKey(it, "fuzzy_rate_chg", tmp_fuzzyRateCHG.str());
            }
            else 
            {
                assignValueByKey(it, "meth_rate_chg", "-1");
                assignValueByKey(it, "fuzzy_rate_chg", "-1");
            }
        }
        // CHH
        // methylated
        unsigned countCHHs = stringAnnoStats[value(it)].methsCHH + stringAnnoStats[value(it)].fuzzyCHH + stringAnnoStats[value(it)].unmethsCHH;
        if (countCHHs > 0)
        {
            double methRateCHH = (double)stringAnnoStats[value(it)].methsCHH/(double)countCHHs;
            std::stringstream tmp_methRateCHH;
            tmp_methRateCHH << methRateCHH;
            assignValueByKey(it, "meth_rate_chh", tmp_methRateCHH.str());
            // fuzzy methylated
            double fuzzyRateCHH = (double)stringAnnoStats[value(it)].fuzzyCHG/(double)countCHHs;
            std::stringstream tmp_fuzzyRateCHH;
            tmp_fuzzyRateCHH << fuzzyRateCHH;
            assignValueByKey(it, "fuzzy_rate_chh", tmp_fuzzyRateCHH.str());
        }
        else
        {
            assignValueByKey(it, "meth_rate_chh", "-1");
            assignValueByKey(it, "fuzzy_rate_chh", "-1");
        }
        goNext(it);
    }
}

template<typename TStore, typename TOptions>
inline void
outputTypeStats(String<AnnoMethStats> &stringTypeStats, TStore &store, TOptions &options)
{
    std::cout << "Annotation type specific methylation rates: " << std::endl;
    std::cout << "Context: \t" << "Methylated: \t" << "Fuzzy methylated: \t" << "Unmethylated: " << std::endl;
    for (unsigned i = 0; i < length(store.annotationTypeStore); ++i)
    {
        if (store.annotationTypeStore[i] == "<root>" || store.annotationTypeStore[i] == "<deleted>") continue;
        std::cout << "-----------------------------------------------------------" << std::endl;
        std::cout << "| " << store.annotationTypeStore[i] << " |" << std::endl;
        std::cout << "-----------------------------------------------------------" << std::endl;
        if (options.refFileName!="")
        {
            unsigned countCG = stringTypeStats[i].methsCG + stringTypeStats[i].fuzzyCG + stringTypeStats[i].unmethsCG;
            if (countCG > 0) std::cout << "CG: \t" << (double)stringTypeStats[i].methsCG/(double)countCG << '\t' 
                                                    << (double)stringTypeStats[i].fuzzyCG/(double)countCG << '\t'
                                                    << (double)stringTypeStats[i].unmethsCG/(double)countCG << '\n';
            else std::cout << "CG: \t" << "-1" << '\t' << "-1" << '\t' << "-1" << '\n';

            unsigned countCHG = stringTypeStats[i].methsCHG + stringTypeStats[i].fuzzyCHG + stringTypeStats[i].unmethsCHG;
            if (countCHG > 0) std::cout << "CHG: \t" << (double)stringTypeStats[i].methsCHG/(double)countCHG << '\t' 
                                                     << (double)stringTypeStats[i].fuzzyCHG/(double)countCHG << '\t'
                                                     << (double)stringTypeStats[i].unmethsCHG/(double)countCHG << '\n';
            else std::cout << "CHG: \t" << "-1" << '\t' << "-1" << '\t' << "-1" << '\n';
        }
        unsigned countCHH = stringTypeStats[i].methsCHH + stringTypeStats[i].fuzzyCHH + stringTypeStats[i].unmethsCHH;
        if (countCHH > 0) std::cout << "CHH: \t" << (double)stringTypeStats[i].methsCHH/(double)countCHH << '\t' 
                                                 << (double)stringTypeStats[i].fuzzyCHH/(double)countCHH << '\t'
                                                 << (double)stringTypeStats[i].unmethsCHH/(double)countCHH << '\n';
        else std::cout << "CHH: \t" << "-1" << '\t' << "-1" << '\t' << "-1" << '\n';
        std::cout << '\n';
    }
}

template<typename TOptions>
bool annoMeths(TOptions & options)
{
    typedef FragmentStore<>                         TStore;
    typedef TStore::TAnnotationStore                TAnnotationStore;
    typedef Value<TAnnotationStore>::Type           TAnnotationStoreElement;
    typedef TAnnotationStoreElement::TId            TId;
    typedef TAnnotationStoreElement::TPos           TPos;
    typedef IntervalAndCargo<TPos, TId>             TInterval;
    typedef IntervalTree<TPos, TId>                 TIntervalTree;


    TStore store;
    // Load reference if given
    if (options.refFileName != "")
         loadContigs(store, options.refFileName);
    // Read annotations from gff
    std::ifstream gffFile(toCString(options.gffFileName), std::ios_base::in | std::ios_base::binary);
    if (!gffFile.good())
    {
        std::cerr << "Couldn't open annotation file" << options.gffFileName << std::endl;
        return false;
    }
    std::cerr << "Loading genome annotation ... " << std::endl;
    read(gffFile, store, Gtf());
    refresh(store.contigNameStoreCache);

    // Create intervalTrees
    std::cerr << "Create interval trees ... " << std::endl;
    String<String<TInterval> > intervalsF;
    String<String<TInterval> > intervalsR;
    extractGeneIntervals(intervalsF, intervalsR, store);
    String<TIntervalTree> intervalTreesF;
    String<TIntervalTree> intervalTreesR;
    constructIntervalTrees(intervalTreesF, intervalTreesR, intervalsF, intervalsR);

    // Read methylation level and count for each annotation which contains the position
    // For each annotation: store counts for each state (unmethylated, fuzzy, methylated) and context (CG, CHG, CHH)
    std::cerr << "Read and count meths ... " << std::endl;
    String<AnnoMethStats> stringAnnoStats;
    String<AnnoMethStats> stringTypeStats;
    parseAndCountMeths(stringAnnoStats, stringTypeStats, intervalTreesF, intervalTreesR, store, options);

    std::cerr << "Compute rates and write output ... " << std::endl;
    // Assign stats to annotation key-value pairs
    assignToKeyValuePairs(store, stringAnnoStats, options);

    // Write annotations (including key-value pairs containing information about methylation rates)
    std::ofstream fileOut(toCString(options.outputFileName), std::ios_base::out | std::ios_base::binary);
    write(fileOut, store, Gtf());

    // Output stats for annotation types
    if (options.showTypeStats) 
        outputTypeStats(stringTypeStats, store, options);

    return 0;
}






#endif  
