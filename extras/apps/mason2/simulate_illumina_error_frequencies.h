

namespace seqan {

struct BsNonSimple_;        // Tag for assuming most simple model: uniform error distributions
typedef seqan::Tag<BsNonSimple_>   BsNonSimple; 

struct BsSimple_;
typedef seqan::Tag<BsSimple_>      BsSimple; 

template <typename TValue = double, typename TSpec = BsSimple>
struct BaseErrorFreqsFrom;

template <typename TValue>
struct BaseErrorFreqsFrom<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE
    };

    static inline TValue const * getData() {
        SEQAN_CHECKPOINT; 
                                                // A,  C,    G,    T     N
        static TValue const _data[VALUE_SIZE] =  {0.25, 0.19, 0.24, 0.33, 0.0};      // siehe Dohm et.al. 2008
        return _data;                                                                // N  stays N
    }
};

template <typename TValue>
struct BaseErrorFreqsFrom<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE
    };

    static inline TValue const * getData() {
        SEQAN_CHECKPOINT; 
        TValue f = 1.0/5.0;
                                                // A,  C,    G,    T     N
        static TValue const _data[VALUE_SIZE] = {f, f, f, f, f};      // siehe Dohm et.al. 2008
        return _data;
    }
};



template <typename TValue = double, typename TSpec = BsNonSimple>
struct SeqErrorFreqs;

// Default sequencing error frequencies (substitutions only)
template <typename TValue>
struct SeqErrorFreqs<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline TValue const * getData() {
        SEQAN_CHECKPOINT;
        
        // Row sum must be 1, since error is given  
        static TValue const _data[TAB_SIZE] = {
                // To  A,             C,            G,             T
                        0,            (0.14/0.25),   (0.05/0.25),    (0.05/0.25),    0,    // From A
                       (0.13/0.19),    0,            (0.02/0.19),    (0.04/0.19),    0,    // C
                       (0.04/0.24),   (0.08/0.24),    0,             (0.12/0.24),    0,
                       (0.08/0.33),   (0.15/0.33),   (0.09/0.33),     0,       0,    
                       0.25,           0.25,          0.25,           0.25,    0 
        };
        return _data;
    }
};

// Simple sequencing error frequencies (substitutions only)
// [ordValue(realBase) * 5 + ordValue(observed base)]
template <typename TValue>
struct SeqErrorFreqs<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline TValue const * getData() {
        SEQAN_CHECKPOINT;
        
        TValue fE = 1.0/3.0;                        // Column sum must be 1    
        static TValue const _data[TAB_SIZE] = {
                // To A,    C,    G,    T,   N
                      0,    fE,   fE,   fE,  0,     // From A
                      fE,   0,    fE,   fE,  0,     // C
                      fE,   fE,   0,    fE,  0,
                      fE,   fE,   fE,   0,   0,
                      0.25,   0.25,   0.25,   0.25,   0
        };
        return _data;
    }
};


template <typename TValue = double, typename TSpec = BsNonSimple>
struct InsErrorFreqsM;

template <typename TValue>
struct InsErrorFreqsM<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
        SEQAN_CHECKPOINT; 
                                                // A,    C,     G,     T,    N   
        static TValue const _data[VALUE_SIZE] = {0.43, 0.065, 0.065, 0.43, 0.01};      // siehe Dohm et.al. 2011 // TODO N?

        return _data;
    }
};

template <typename TValue>
struct InsErrorFreqsM<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
        SEQAN_CHECKPOINT; 
        TValue f = 1.0/4.0;                   
        static TValue const _data[VALUE_SIZE] = {f, f, f, f, 0};   
        return _data;
    }
};


template <typename TValue = double, typename TSpec = BsNonSimple>
struct DelErrorFreqsM;

template <typename TValue>
struct DelErrorFreqsM<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
        SEQAN_CHECKPOINT; 
                                                // A,    C,     G,     T,    N   
        static TValue const _data[VALUE_SIZE] = {0.42, 0.075, 0.075, 0.42, 0.00};      // siehe Dohm et.al. 2011  N?

        return _data;
    }
};

template <typename TValue>
struct DelErrorFreqsM<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
        SEQAN_CHECKPOINT; 
        TValue f = 1.0/5.0;
        static TValue const _data[VALUE_SIZE] = {f, f, f, f, f};
        return _data;
    }
};



}


