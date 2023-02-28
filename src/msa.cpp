#include <msa.hpp>
#include <hostdevicefunctions.cuh>
#include <alignmentorientation.hpp>

#include <map>
#include <bitset>
namespace care{

void MultipleSequenceAlignment::build(const InputData& args){

    assert(args.anchorLength > 0);
    assert(args.anchor != nullptr);

    inputData = args;
    addedSequences = 0;

    nCandidates = args.nCandidates;

    //determine number of columns in pileup image
    int startindex = 0;
    int endindex = args.anchorLength;

    for(int i = 0; i < nCandidates; ++i){
        const int shift = args.candidateShifts[i];
        const int candidateEndsAt = args.candidateLengths[i] + shift;
        startindex = std::min(shift, startindex);
        endindex = std::max(candidateEndsAt, endindex);
    }

    nColumns = endindex - startindex;

    anchorColumnsBegin_incl = std::max(-startindex,0);
    anchorColumnsEnd_excl = anchorColumnsBegin_incl + args.anchorLength;

    resize(nColumns);

    fillzero();

    addSequence(args.useQualityScores, args.anchor, args.anchorQualities, args.anchorLength, 0, 1.0f);

    for(int candidateIndex = 0; candidateIndex < nCandidates; candidateIndex++){
        const char* ptr = args.candidates + candidateIndex * args.candidatesPitch;
        const char* qptr = args.candidateQualities + candidateIndex * args.candidateQualitiesPitch;
        const int candidateLength = args.candidateLengths[candidateIndex];
        const int shift = args.candidateShifts[candidateIndex];
        const float defaultWeightFactor = args.candidateDefaultWeightFactors[candidateIndex];

        addSequence(args.useQualityScores, ptr, qptr, candidateLength, shift, defaultWeightFactor);
    }

    findConsensus();

    findOrigWeightAndCoverage(args.anchor);
}

void MultipleSequenceAlignment::resize(int cols){

    consensus.resize(cols);
    support.resize(cols);
    coverage.resize(cols);
    origWeights.resize(cols);
    origCoverages.resize(cols);
    countsA.resize(cols);
    countsC.resize(cols);
    countsG.resize(cols);
    countsT.resize(cols);
    weightsA.resize(cols);
    weightsC.resize(cols);
    weightsG.resize(cols);
    weightsT.resize(cols);
}

void MultipleSequenceAlignment::fillzero(){
    auto zero = [](auto& vec){
        std::fill(vec.begin(), vec.end(), 0);
    };

    zero(consensus);
    zero(support);
    zero(coverage);
    zero(origWeights);
    zero(origCoverages);
    zero(countsA);
    zero(countsC);
    zero(countsG);
    zero(countsT);
    zero(weightsA);
    zero(weightsC);
    zero(weightsG);
    zero(weightsT);
}

void MultipleSequenceAlignment::addSequence(bool useQualityScores, const char* sequence, const char* quality, int length, int shift, float defaultWeightFactor){
    assert(sequence != nullptr);
    assert(!useQualityScores || quality != nullptr);

    for(int i = 0; i < length; i++){
        const int globalIndex = anchorColumnsBegin_incl + shift + i;
        const char base = sequence[i];
        const float weight = defaultWeightFactor * (useQualityScores ? qualityConversion->getWeight(quality[i]) : 1.0f);
        switch(base){
            case 'A': countsA[globalIndex]++; weightsA[globalIndex] += weight;break;
            case 'C': countsC[globalIndex]++; weightsC[globalIndex] += weight;break;
            case 'G': countsG[globalIndex]++; weightsG[globalIndex] += weight;break;
            case 'T': countsT[globalIndex]++; weightsT[globalIndex] += weight;break;
            default: assert(false); break;
        }
        coverage[globalIndex]++;
    }

    addedSequences++;
}

/*
void MultipleSequenceAlignment::removeSequence(bool useQualityScores, const char* sequence, const char* quality, int length, int shift, float defaultWeightFactor){
    assert(sequence != nullptr);
    assert(!useQualityScores || quality != nullptr);

    for(int i = 0; i < length; i++){
        const int globalIndex = anchorColumnsBegin_incl + shift + i;
        const char base = sequence[i];
        const float weight = defaultWeightFactor * (useQualityScores ? qualityConversion.getWeight(quality[i]) : 1.0f);
        switch(base){
            case 'A': countsA[globalIndex]--; weightsA[globalIndex] -= weight;break;
            case 'C': countsC[globalIndex]--; weightsC[globalIndex] -= weight;break;
            case 'G': countsG[globalIndex]--; weightsG[globalIndex] -= weight;break;
            case 'T': countsT[globalIndex]--; weightsT[globalIndex] -= weight;break;
            default: assert(false); break;
        }
        coverage[globalIndex]--;
    }
}
*/

void MultipleSequenceAlignment::findConsensus(){
    for(int column = 0; column < nColumns; ++column){
        char cons = 'A';
        float consWeight = weightsA[column];
        if(weightsC[column] > consWeight){
            cons = 'C';
            consWeight = weightsC[column];
        }
        if(weightsG[column] > consWeight){
            cons = 'G';
            consWeight = weightsG[column];
        }
        if(weightsT[column] > consWeight){
            cons = 'T';
            consWeight = weightsT[column];
        }
        consensus[column] = cons;

        const float columnWeight = weightsA[column] + weightsC[column] + weightsG[column] + weightsT[column];
        support[column] = consWeight / columnWeight;
        assert(!std::isnan(support[column]));
    }
}

void MultipleSequenceAlignment::findOrigWeightAndCoverage(const char* anchor){
    for(int column = anchorColumnsBegin_incl; column < anchorColumnsEnd_excl; ++column){

        const int localIndex = column - anchorColumnsBegin_incl;
        const char anchorBase = anchor[localIndex];
        switch(anchorBase){
            case 'A':origWeights[column] = weightsA[column]; origCoverages[column] = countsA[column]; break;
            case 'C':origWeights[column] = weightsC[column]; origCoverages[column] = countsC[column]; break;
            case 'G':origWeights[column] = weightsG[column]; origCoverages[column] = countsG[column]; break;
            case 'T':origWeights[column] = weightsT[column]; origCoverages[column] = countsT[column]; break;
            default: assert(false); break;
        }

    }
}


void MultipleSequenceAlignment::print(std::ostream& os) const{
    std::vector<int> indices(nCandidates+1);
    std::iota(indices.begin(), indices.end(), 0);

    auto get_shift_of_row = [&](int k){
        if(k == 0) return 0;
        else return inputData.candidateShifts[k-1];
    };

    std::sort(indices.begin(), indices.end(),
              [&](int l, int r){return get_shift_of_row(l) < get_shift_of_row(r);});

    for(int row = 0; row < nCandidates+1; row++) {
        int sortedrow = indices[row];

        if(sortedrow == 0){
            os << ">> ";

            for(int i = 0; i < anchorColumnsBegin_incl; i++){
                os << "0";
            }

            for(int i = 0; i < inputData.anchorLength; i++){
                os << inputData.anchor[i];
            }

            for(int i = anchorColumnsEnd_excl; i < nColumns; i++){
                os << "0";
            }

            os << " <<";
        }else{
            os << "   ";
            int written = 0;
            for(int i = 0; i < anchorColumnsBegin_incl + get_shift_of_row(sortedrow); i++){
                os << "0";
                written++;
            }

            for(int i = 0; i < inputData.candidateLengths[sortedrow-1]; i++){
                os << inputData.candidates[(sortedrow-1) * inputData.candidatesPitch + i];
                written++;
            }

            for(int i = anchorColumnsBegin_incl + get_shift_of_row(sortedrow) 
                        + inputData.candidateLengths[sortedrow-1]; 
                    i < nColumns; i++){
                os << "0";
                written++;
            }

            assert(written == nColumns);

            os << "   " << inputData.candidateLengths[sortedrow-1] << " " << get_shift_of_row(sortedrow);
        }

        os << '\n';
    }
}


void MultipleSequenceAlignment::printWithDiffToConsensus(std::ostream& os) const{
    std::vector<int> indices(nCandidates+1);
    std::iota(indices.begin(), indices.end(), 0);

    auto get_shift_of_row = [&](int k){
        if(k == 0) return 0;
        else return inputData.candidateShifts[k-1];
    };

    std::sort(indices.begin(), indices.end(),
              [&](int l, int r){return get_shift_of_row(l) < get_shift_of_row(r);});

    for(int row = 0; row < nCandidates+1; row++) {
        int sortedrow = indices[row];

        if(sortedrow == 0){
            os << ">> ";

            for(int i = 0; i < anchorColumnsBegin_incl; i++){
                os << "0";
            }

            for(int i = 0; i < inputData.anchorLength; i++){
                const int globalIndex = anchorColumnsBegin_incl + i;
                const char c = consensus[globalIndex] == inputData.anchor[i] ? '=' : inputData.anchor[i];
                os << c;
            }

            for(int i = anchorColumnsEnd_excl; i < nColumns; i++){
                os << "0";
            }

            os << " <<";
        }else{
            os << "   ";
            int written = 0;
            for(int i = 0; i < anchorColumnsBegin_incl + get_shift_of_row(sortedrow); i++){
                os << "0";
                written++;
            }

            for(int i = 0; i < inputData.candidateLengths[sortedrow-1]; i++){
                const int globalIndex = anchorColumnsBegin_incl + get_shift_of_row(sortedrow) + i;
                const char base = inputData.candidates[(sortedrow-1) * inputData.candidatesPitch + i];
                const char c = consensus[globalIndex] == base ? '=' : base;

                os << c;
                written++;
            }

            for(int i = anchorColumnsBegin_incl + get_shift_of_row(sortedrow) 
                        + inputData.candidateLengths[sortedrow-1]; 
                    i < nColumns; i++){
                os << "0";
                written++;
            }

            assert(written == nColumns);

            os << "   " << inputData.candidateLengths[sortedrow-1] << " " << get_shift_of_row(sortedrow);
        }

        os << '\n';
    }
}


MSAProperties MultipleSequenceAlignment::getMSAProperties(
    int firstCol,
    int lastCol, //exclusive
    float estimatedErrorrate,
    float estimatedCoverage,
    float m_coverage
) const {

    assert(firstCol <= lastCol);

    const float avg_support_threshold = 1.0f-1.0f*estimatedErrorrate;
    const float min_support_threshold = 1.0f-3.0f*estimatedErrorrate;
    const float min_coverage_threshold = m_coverage / 6.0f * estimatedCoverage;

    //const int firstCol = anchorColumnsBegin_incl; //0;
    //const int lastCol = anchorColumnsEnd_excl; //nColumns; //exclusive
    const int distance = lastCol - firstCol;

    MSAProperties msaProperties;

    msaProperties.min_support = *std::min_element(support.data() + firstCol, support.data() + lastCol);

    const float supportsum = std::accumulate(support.data() + firstCol, support.data() + lastCol, 0.0f);
    msaProperties.avg_support = supportsum / distance;

    auto minmax = std::minmax_element(coverage.data() + firstCol, coverage.data() + lastCol);

    msaProperties.min_coverage = *minmax.first;
    msaProperties.max_coverage = *minmax.second;

    auto isGoodAvgSupport = [=](float avgsupport){
        return fgeq(avgsupport, avg_support_threshold);
    };
    auto isGoodMinSupport = [=](float minsupport){
        return fgeq(minsupport, min_support_threshold);
    };
    auto isGoodMinCoverage = [=](float mincoverage){
        return fgeq(mincoverage, min_coverage_threshold);
    };

    msaProperties.failedAvgSupport = !isGoodAvgSupport(msaProperties.avg_support);
    msaProperties.failedMinSupport = !isGoodMinSupport(msaProperties.min_support);
    msaProperties.failedMinCoverage = !isGoodMinCoverage(msaProperties.min_coverage);


    const float avg_support = msaProperties.avg_support;
    const float min_support = msaProperties.min_support;
    const int min_coverage = msaProperties.min_coverage;

    msaProperties.isHQ = false;

    const bool allGood = isGoodAvgSupport(avg_support) 
                                        && isGoodMinSupport(min_support) 
                                        && isGoodMinCoverage(min_coverage);
    if(allGood){
        int smallestErrorrateThatWouldMakeHQ = 100;

        const int estimatedErrorratePercent = ceil(estimatedErrorrate * 100.0f);
        for(int percent = estimatedErrorratePercent; percent >= 0; percent--){
            const float factor = percent / 100.0f;
            const float avg_threshold = 1.0f - 1.0f * factor;
            const float min_threshold = 1.0f - 3.0f * factor;
            if(fgeq(avg_support, avg_threshold) && fgeq(min_support, min_threshold)){
                smallestErrorrateThatWouldMakeHQ = percent;
            }
        }

        msaProperties.isHQ = isGoodMinCoverage(min_coverage)
                            && fleq(smallestErrorrateThatWouldMakeHQ, estimatedErrorratePercent * 0.5f);
    }

    return msaProperties;
}


CorrectionResult MultipleSequenceAlignment::getCorrectedAnchor(
    MSAProperties msaProperties,
    float estimatedErrorrate,
    float estimatedCoverage,
    float m_coverage,
    int /*neighborRegionSize*/,
    read_number /* readId*/
) const {

    if(nCandidates == 0){
        //cannot be corrected without candidates

        CorrectionResult result{};
        result.isCorrected = false;
        return result;
    }

    auto canBeCorrectedByConsensus = [&](){

        const float avg_support_threshold = 1.0f-1.0f*estimatedErrorrate;
        const float min_support_threshold = 1.0f-3.0f*estimatedErrorrate;
        const float min_coverage_threshold = m_coverage / 6.0f * estimatedCoverage;

        auto isGoodAvgSupport = [=](float avgsupport){
            return fgeq(avgsupport, avg_support_threshold);
        };
        auto isGoodMinSupport = [=](float minsupport){
            return fgeq(minsupport, min_support_threshold);
        };
        auto isGoodMinCoverage = [=](float mincoverage){
            return fgeq(mincoverage, min_coverage_threshold);
        };

        const float avg_support = msaProperties.avg_support;
        const float min_support = msaProperties.min_support;
        const int min_coverage = msaProperties.min_coverage;

        return isGoodAvgSupport(avg_support) 
                                        && isGoodMinSupport(min_support) 
                                        && isGoodMinCoverage(min_coverage);
    };

    CorrectionResult result{};
    result.isCorrected = false;
    result.correctedSequence.resize(inputData.anchorLength);

    result.isCorrected = true;

    
    int flag = 0;    

    if(canBeCorrectedByConsensus()){
        flag = msaProperties.isHQ ? 2 : 1;
    }

    if(flag > 0){
        std::copy_n(
            consensus.data() + anchorColumnsBegin_incl,
            inputData.anchorLength,
            result.correctedSequence.begin()
        );
    }else{
        //correct only positions with high support to consensus, else leave position unchanged.
        for(int i = 0; i < inputData.anchorLength; i += 1){
            //assert(consensus[i] == 'A' || consensus[i] == 'C' || consensus[i] == 'G' || consensus[i] == 'T');
            if(support[anchorColumnsBegin_incl + i] > 0.90f && origCoverages[anchorColumnsBegin_incl + i] <= 2){
                result.correctedSequence[i] = consensus[anchorColumnsBegin_incl + i];
            }else{
                result.correctedSequence[i] = inputData.anchor[i];
            }
        }
    }

    return result;
}



std::vector<CorrectedCandidate> MultipleSequenceAlignment::getCorrectedCandidates(
    float estimatedErrorrate,
    float estimatedCoverage,
    float m_coverage,
    int new_columns_to_correct
) const {

    //const float avg_support_threshold = 1.0f-1.0f*estimatedErrorrate;
    const float min_support_threshold = 1.0f-3.0f*estimatedErrorrate;
    const float min_coverage_threshold = m_coverage / 6.0f * estimatedCoverage;

    std::vector<CorrectedCandidate> result;
    result.reserve(inputData.nCandidates);

    for(int candidate_index = 0; candidate_index < inputData.nCandidates; ++candidate_index){

        const int queryColumnsBegin_incl = anchorColumnsBegin_incl + inputData.candidateShifts[candidate_index];
        const int candidateLength = inputData.candidateLengths[candidate_index];
        const int queryColumnsEnd_excl = queryColumnsBegin_incl + candidateLength;

        bool candidateShouldBeCorrected = false;

        //check range condition and length condition
        if(anchorColumnsBegin_incl - new_columns_to_correct <= queryColumnsBegin_incl
            && queryColumnsBegin_incl <= anchorColumnsBegin_incl + new_columns_to_correct
            && queryColumnsEnd_excl <= anchorColumnsEnd_excl + new_columns_to_correct){

            float newColMinSupport = 1.0f;
            int newColMinCov = std::numeric_limits<int>::max();

            //check new columns left of anchor
            for(int columnindex = anchorColumnsBegin_incl - new_columns_to_correct;
                columnindex < anchorColumnsBegin_incl;
                columnindex++){

                assert(columnindex < nColumns);

                if(queryColumnsBegin_incl <= columnindex){
                    newColMinSupport = support[columnindex] < newColMinSupport ? support[columnindex] : newColMinSupport;
                    newColMinCov = coverage[columnindex] < newColMinCov ? coverage[columnindex] : newColMinCov;
                }
            }
            //check new columns right of anchor
            for(int columnindex = anchorColumnsEnd_excl;
                columnindex < anchorColumnsEnd_excl + new_columns_to_correct
                && columnindex < nColumns;
                columnindex++){

                newColMinSupport = support[columnindex] < newColMinSupport ? support[columnindex] : newColMinSupport;
                newColMinCov = coverage[columnindex] < newColMinCov ? coverage[columnindex] : newColMinCov;
            }

            candidateShouldBeCorrected = fgeq(newColMinSupport, min_support_threshold)
                            && fgeq(newColMinCov, min_coverage_threshold);

            candidateShouldBeCorrected = true;
        }

        if(candidateShouldBeCorrected){
            std::string correctedString(&consensus[queryColumnsBegin_incl], &consensus[queryColumnsEnd_excl]);
            result.emplace_back(candidate_index, inputData.candidateShifts[candidate_index], std::move(correctedString));
        }
    }

    return result;
}





//remove all candidate reads from alignment which are assumed to originate from a different genomic region
//the indices of remaining candidates are returned in MinimizationResult::remaining_candidates
//candidates in vector must be in the same order as they were inserted into the msa!!!

RegionSelectionResult MultipleSequenceAlignment::findCandidatesOfDifferentRegion(
    int dataset_coverage
) const {

    auto is_significant_count = [&](int count, int dataset_coverage){
        if(int(dataset_coverage * 0.3f) <= count)
            return true;
        return false;
    };

    constexpr std::array<char, 4> index_to_base{'A','C','G','T'};

    //find column with a non-consensus base with significant coverage
    int col = 0;
    bool foundColumn = false;
    char foundBase = 'F';
    int foundBaseIndex = 0;
    int consindex = 0;

    //if anchor has no mismatch to consensus, don't minimize
    auto pair = std::mismatch(inputData.anchor,
                                inputData.anchor + inputData.anchorLength,
                                consensus.data() + anchorColumnsBegin_incl);

    if(pair.first == inputData.anchor + inputData.anchorLength){
        RegionSelectionResult result;
        result.performedMinimization = false;
        return result;
    }

    for(int columnindex = anchorColumnsBegin_incl; columnindex < anchorColumnsEnd_excl && !foundColumn; columnindex++){
        std::array<int,4> counts;
        //std::array<float,4> weights;

        counts[0] = countsA[columnindex];
        counts[1] = countsC[columnindex];
        counts[2] = countsG[columnindex];
        counts[3] = countsT[columnindex];

        /*weights[0] = weightsA[columnindex];
        weights[1] = weightsC[columnindex];
        weights[2] = weightsG[columnindex];
        weights[3] = weightsT[columnindex];*/

        char cons = consensus[columnindex];
        consindex = -1;

        switch(cons){
            case 'A': consindex = 0;break;
            case 'C': consindex = 1;break;
            case 'G': consindex = 2;break;
            case 'T': consindex = 3;break;
        }

        //const char originalbase = anchor[columnindex - columnProperties.anchorColumnsBegin_incl];

        //find out if there is a non-consensus base with significant coverage
        int significantBaseIndex = -1;
        //int maxcount = 0;
        for(int i = 0; i < 4; i++){
            if(i != consindex){
                bool significant = is_significant_count(counts[i], dataset_coverage);

                bool process = significant; //maxcount < counts[i] && significant && (cons == originalbase || index_to_base[i] == originalbase);

                significantBaseIndex = process ? i : significantBaseIndex;

                //maxcount = process ? std::max(maxcount, counts[i]) : maxcount;
            }
        }

        if(significantBaseIndex != -1){
            foundColumn = true;
            col = columnindex;
            foundBase = index_to_base[significantBaseIndex];
            foundBaseIndex = significantBaseIndex;

            //printf("found col %d, baseIndex %d\n", col, foundBaseIndex);
        }
    }



    RegionSelectionResult result;
    result.performedMinimization = foundColumn;
    result.column = col;

    if(foundColumn){

        result.differentRegionCandidate.resize(nCandidates);

        //compare found base to original base
        const char originalbase = inputData.anchor[col - anchorColumnsBegin_incl];

        result.significantBase = foundBase;
        result.originalBase = originalbase;
        result.consensusBase = consensus[col];

        std::array<int,4> counts;

        counts[0] = countsA[col];
        counts[1] = countsC[col];
        counts[2] = countsG[col];
        counts[3] = countsT[col];

        result.significantCount = counts[foundBaseIndex];
        result.consensuscount = counts[consindex];

        auto discard_rows = [&](bool keepMatching){

            std::array<int, 4> seenCounts{0,0,0,0};

            for(int candidateIndex = 0; candidateIndex < nCandidates; candidateIndex++){
                //check if row is affected by column col
                const int row_begin_incl = anchorColumnsBegin_incl + inputData.candidateShifts[candidateIndex];
                const int row_end_excl = row_begin_incl + inputData.candidateLengths[candidateIndex];
                const bool notAffected = (col < row_begin_incl || row_end_excl <= col);
                const char base = notAffected ? 'F' : inputData.candidates[candidateIndex * inputData.candidatesPitch + (col - row_begin_incl)];

                /*printf("k %d, candidateIndex %d, row_begin_incl %d, row_end_excl %d, notAffected %d, base %c\n", candidateIndex, candidateIndex,
                row_begin_incl, row_end_excl, notAffected, base);
                for(int i = 0; i < row_end_excl - row_begin_incl; i++){
                    if(i == (col - row_begin_incl))
                        printf("_");
                    printf("%c", candidates[candidateIndex * candidatesPitch + i]);
                    if(i == (col - row_begin_incl))
                        printf("_");
                }
                printf("\n");*/

                if(base == 'A') seenCounts[0]++;
                if(base == 'C') seenCounts[1]++;
                if(base == 'G') seenCounts[2]++;
                if(base == 'T') seenCounts[3]++;

                if(notAffected){
                    result.differentRegionCandidate[candidateIndex] = false;
                }else if(keepMatching && (base == foundBase)){
                    //keep candidates which match the found base
                    result.differentRegionCandidate[candidateIndex] = false;
                }else if(!keepMatching && (base != foundBase)){
                    //keep candidates which do not match the found base
                    result.differentRegionCandidate[candidateIndex] = false;
                }else{
                    result.differentRegionCandidate[candidateIndex] = true;
                }
            }

            if(originalbase == 'A') seenCounts[0]++;
            if(originalbase == 'C') seenCounts[1]++;
            if(originalbase == 'G') seenCounts[2]++;
            if(originalbase == 'T') seenCounts[3]++;

            assert(seenCounts[0] == countsA[col]);
            assert(seenCounts[1] == countsC[col]);
            assert(seenCounts[2] == countsG[col]);
            assert(seenCounts[3] == countsT[col]);


#if 1
            //check that no candidate which should be removed has very good alignment.
            //if there is such a candidate, none of the candidates will be removed.
            bool veryGoodAlignment = false;
            for(int candidateIndex = 0; candidateIndex < nCandidates; candidateIndex++){
                if(result.differentRegionCandidate[candidateIndex]){
                    const float overlapweight = inputData.candidateDefaultWeightFactors[candidateIndex];
                    assert(overlapweight <= 1.0f);
                    assert(overlapweight >= 0.0f);

                    if(overlapweight >= 0.90f){
                        veryGoodAlignment = true;
                    }
                }
            }

            if(veryGoodAlignment){
                for(int candidateIndex = 0; candidateIndex < nCandidates; candidateIndex++){
                    result.differentRegionCandidate[candidateIndex] = false;
                }
            }
#endif
        };



        if(originalbase == foundBase){
            //discard all candidates whose base in column col differs from foundBase
            discard_rows(true);
        }else{
            //discard all candidates whose base in column col matches foundBase
            discard_rows(false);
        }

        //if(result.num_discarded_candidates > 0){
        //    find_consensus();
        //}

        return result;
    }else{

        return result;
    }
}


MultipleSequenceAlignment::PossibleMsaSplits MultipleSequenceAlignment::inspectColumnsRegionSplit(int firstColumn) const{
    return inspectColumnsRegionSplit(firstColumn, nColumns);
}

MultipleSequenceAlignment::PossibleMsaSplits MultipleSequenceAlignment::inspectColumnsRegionSplit(int firstColumn, int lastColumnExcl) const{
    assert(lastColumnExcl >= 0);
    assert(lastColumnExcl <= nColumns);

    std::vector<PossibleSplitColumn> possibleColumns;

    //find columns in which two different nucleotides each make up approx 50% (40% - 60%) of the column
    for(int col = firstColumn; col < lastColumnExcl; col++){
        std::array<PossibleSplitColumn, 4> array{};
        int numPossibleNucs = 0;

        auto checkNuc = [&](const auto& counts, const char nuc){
            const float ratio = float(counts[col]) / float(coverage[col]);
            //if((counts[col] == 2 && fgeq(ratio, 0.4f) && fleq(ratio, 0.6f)) || counts[col] > 2){
            if(counts[col] >= 2 && fgeq(ratio, 0.4f) && fleq(ratio, 0.6f)){
                array[numPossibleNucs] = {nuc, col, ratio};
                numPossibleNucs++;
            }
        };

        checkNuc(countsA, 'A');
        checkNuc(countsC, 'C');
        checkNuc(countsG, 'G');
        checkNuc(countsT, 'T');

        
        if(numPossibleNucs == 2){
            possibleColumns.insert(possibleColumns.end(), array.begin(), array.begin() + numPossibleNucs);
        }
    }

    assert(possibleColumns.size() % 2 == 0);

    if(possibleColumns.size() >= 2 && possibleColumns.size() <= 32){

        PossibleMsaSplits result;

        //calculate proper results
        {
            

            //for each candidate check if it appears in the selected columns
            //there can be at most 16 selected columns.the result is bit-encoded in 32 bits, 2 bits per column
            //0b10 - candidate base matches first selected nuc in column
            //0b11 - candidate base matches second selected nuc in column
            //0b00 - candidate base does not match, or candidate does not occupy column

            //map groups candidates with the same bit pattern
            std::map<unsigned int, std::vector<int>> map;
            
            for(int l = 0; l < nCandidates; l++){
                const int candidateRow = l;

                constexpr std::size_t numPossibleColumnsPerFlag = sizeof(unsigned int) * CHAR_BIT / 2;

                unsigned int flags = 0;

                const int numPossibleColumns = std::min(numPossibleColumnsPerFlag, possibleColumns.size() / 2);
                const char* const candidateString = &inputData.candidates[candidateRow * inputData.candidatesPitch];

                for(int k = 0; k < numPossibleColumns; k++){
                    flags <<= 2;

                    const PossibleSplitColumn psc0 = possibleColumns[2*k+0];
                    const PossibleSplitColumn psc1 = possibleColumns[2*k+1];
                    assert(psc0.column == psc1.column);

                    const int candidateColumnsBegin_incl = inputData.candidateShifts[candidateRow] + anchorColumnsBegin_incl;
                    const int candidateColumnsEnd_excl = inputData.candidateLengths[candidateRow] + candidateColumnsBegin_incl;
                    
                    //column range check for row
                    if(candidateColumnsBegin_incl <= psc0.column && psc0.column < candidateColumnsEnd_excl){
                        const int positionInCandidate = psc0.column - candidateColumnsBegin_incl;

                        if(candidateString[positionInCandidate] == psc0.letter){
                            flags = flags | 0b10;
                        }else if(candidateString[positionInCandidate] == psc1.letter){
                            flags = flags | 0b11;
                        }else{
                            flags = flags | 0b00;
                        } 

                    }else{
                        flags = flags | 0b00;
                    } 

                }

                map[flags].emplace_back(l);
            }

            std::vector<std::pair<unsigned int, std::vector<int>>> flatmap(map.begin(), map.end());

            std::map<unsigned int, std::vector<int>> finalMap;

            const int flatmapsize = flatmap.size();
            for(int i = 0; i < flatmapsize; i++){
                //try to merge flatmap[i] with flatmap[k], i < k, if possible
                const unsigned int flagsToSearch = flatmap[i].first;
                unsigned int mask = 0;
                for(int s = 0; s < 16; s++){
                    if(flagsToSearch >> (2*s+1) & 1){
                        mask = mask | (0x03 << (2*s));
                    }
                }

                //std::cerr << "i = " << i << ", flags = " << std::bitset<32>(flagsToSearch) << ", mask = " << std::bitset<32>(mask) << "\n";

                bool merged = false;
                for(int k = i+1; k < flatmapsize; k++){
                    //if both columns are identical not including wildcard columns
                    if((mask & flatmap[k].first) == flagsToSearch){
                        //std::cerr << "k = " << k << ", flags = " << std::bitset<32>(flatmap[k].first) << " equal" << "\n";
                        flatmap[k].second.insert(
                            flatmap[k].second.end(),
                            flatmap[i].second.begin(),
                            flatmap[i].second.end()
                        );

                        std::sort(flatmap[k].second.begin(), flatmap[k].second.end());

                        flatmap[k].second.erase(
                            std::unique(flatmap[k].second.begin(), flatmap[k].second.end()),
                            flatmap[k].second.end()
                        );

                        merged = true;
                    }else{
                        //std::cerr << "k = " << k << ", flags = " << std::bitset<32>(flatmap[k].first) << " not equal" << "\n";
                    }
                }

                if(!merged){
                    finalMap[flatmap[i].first] = std::move(flatmap[i].second);
                }
            }

            for(auto& pair : finalMap){

                std::vector<int> listOfCandidates = std::move(pair.second);
                std::vector<PossibleSplitColumn> columnInfo;

                const unsigned int flag = pair.first;
                const int num = possibleColumns.size() / 2;
                for(int i = 0; i < num; i++){
                    const unsigned int cur = (flag >> (num - i - 1) * 2) & 0b11;
                    const bool match = (cur & 0b10) == 0b10;
                    if(match){
                        const int which = cur & 1;
                        columnInfo.emplace_back(possibleColumns[2*i + which]);
                    }
                }

                result.splits.emplace_back(std::move(columnInfo), std::move(listOfCandidates));                
            }
        }

#if 0
        //calculate sorted and debug prints
        {

            std::map<unsigned int, std::vector<int>> debugprintmap;
            std::vector<int> sortedindices(nCandidates);
            std::iota(sortedindices.begin(), sortedindices.end(), 0);

            auto get_shift_of_row = [&](int k){
                return inputData.candidateShifts[k];
            };

            std::sort(sortedindices.begin(), sortedindices.end(),
                    [&](int l, int r){return get_shift_of_row(l) < get_shift_of_row(r);});

            for(int l = 0; l < nCandidates; l++){
                const int candidateRow = sortedindices[l];

                constexpr std::size_t numPossibleColumnsPerFlag = sizeof(unsigned int) * CHAR_BIT / 2;

                unsigned int flags = 0;

                const int numPossibleColumns = std::min(numPossibleColumnsPerFlag, possibleColumns.size() / 2);
                const char* const candidateString = &inputData.candidates[candidateRow * inputData.candidatesPitch];

                for(int k = 0; k < numPossibleColumns; k++){
                    flags <<= 2;

                    const PossibleSplitColumn psc0 = possibleColumns[2*k+0];
                    const PossibleSplitColumn psc1 = possibleColumns[2*k+1];
                    assert(psc0.column == psc1.column);

                    const int candidateColumnsBegin_incl = inputData.candidateShifts[candidateRow] + anchorColumnsBegin_incl;
                    const int candidateColumnsEnd_excl = inputData.candidateLengths[candidateRow] + candidateColumnsBegin_incl;
                    
                    //column range check for row
                    if(candidateColumnsBegin_incl <= psc0.column && psc0.column < candidateColumnsEnd_excl){
                        const int positionInCandidate = psc0.column - candidateColumnsBegin_incl;

                        if(candidateString[positionInCandidate] == psc0.letter){
                            flags = flags | 0b10;
                        }else if(candidateString[positionInCandidate] == psc1.letter){
                            flags = flags | 0b11;
                        }else{
                            flags = flags | 0b00;
                        } 

                    }else{
                        flags = flags | 0b00;
                    } 

                }

                debugprintmap[flags].emplace_back(l);
            }

            std::vector<std::pair<unsigned int, std::vector<int>>> flatmap(debugprintmap.begin(), debugprintmap.end());

            std::map<unsigned int, std::vector<int>> finalMap;

            const int flatmapsize = flatmap.size();
            for(int i = 0; i < flatmapsize; i++){
                //try to merge flatmap[i] with flatmap[k], i < k, if possible
                const unsigned int flagsToSearch = flatmap[i].first;
                unsigned int mask = 0;
                for(int s = 0; s < 16; s++){
                    if(flagsToSearch >> (2*s+1) & 1){
                        mask = mask | (0x03 << (2*s));
                    }
                }

                //std::cerr << "i = " << i << ", flags = " << std::bitset<32>(flagsToSearch) << ", mask = " << std::bitset<32>(mask) << "\n";

                bool merged = false;
                for(int k = i+1; k < flatmapsize; k++){
                    //if both columns are identical not including wildcard columns
                    if((mask & flatmap[k].first) == flagsToSearch){
                        //std::cerr << "k = " << k << ", flags = " << std::bitset<32>(flatmap[k].first) << " equal" << "\n";
                        flatmap[k].second.insert(
                            flatmap[k].second.end(),
                            flatmap[i].second.begin(),
                            flatmap[i].second.end()
                        );

                        std::sort(flatmap[k].second.begin(), flatmap[k].second.end());

                        flatmap[k].second.erase(
                            std::unique(flatmap[k].second.begin(), flatmap[k].second.end()),
                            flatmap[k].second.end()
                        );

                        merged = true;
                    }else{
                        //std::cerr << "k = " << k << ", flags = " << std::bitset<32>(flatmap[k].first) << " not equal" << "\n";
                    }
                }

                if(!merged){
                    finalMap[flatmap[i].first] = std::move(flatmap[i].second);
                }
            }

            auto printMap = [&](const auto& map){
                for(const auto& pair : map){
                    const unsigned int flag = pair.first;
                    //convert flag to position and nuc
                    const int num = possibleColumns.size() / 2;

                    std::cerr << "flag " << flag << " : ";
                    for(int i = 0; i < num; i++){
                        const unsigned int cur = (flag >> (num - i - 1) * 2) & 0b11;
                        const bool match = (cur & 0b10) == 0b10;
                        const int column = possibleColumns[2*i].column;
                        char nuc = '-';
                        if(match){
                            int which = cur & 1;
                            nuc = possibleColumns[2*i + which].letter;
                        }
                        std::cerr << "(" << nuc << ", " << column << ") ";
                    }
                    std::cerr << ": ";

                    for(int c : pair.second){
                        std::cerr << c << " ";
                    }
                    std::cerr << "\n";
                }
            };

            if(possibleColumns.size() > 0){
                std::cerr << possibleColumns.size() << "\n";
                for(const auto& p : possibleColumns){
                    std::cerr << "{" << p.letter << ", " << p.column << ", " << p.ratio << "} ";
                }
                std::cerr << "\n";
                
                printMap(debugprintmap);

                std::cerr << "final map: \n";

                printMap(finalMap);
                std::cerr << "\n";

                // print(std::cerr);
                // std::cerr << "\n";
            }
        }
#endif    
        return result;
    }else{
        // single split with all candidates
        std::vector<int> listOfCandidates(nCandidates);
        std::iota(listOfCandidates.begin(), listOfCandidates.end(), 0);

        PossibleMsaSplits result;

        std::vector<PossibleSplitColumn> columnInfo;


        result.splits.emplace_back(std::move(columnInfo), std::move(listOfCandidates));
        
        return result;
    }
    // for(int which = 0; which < 4; which++){

    //     const std::vector<int>* vec = nullptr;
    //     if(which == 0) vec = &countsMatrixA;
    //     if(which == 1) vec = &countsMatrixC;
    //     if(which == 2) vec = &countsMatrixG;
    //     if(which == 3) vec = &countsMatrixT;


    //     for(int row = 0; row < addedSequences; row++){
    //         const int* const data = vec->data() + row * nColumns;
    //         for(int col = 0; col < nColumns; col++){
    //             os << data[col] << " ";
    //         }
    //         os << "\n";
    //     }
    // }

}

std::vector<MultipleSequenceAlignment::PossibleSplitColumn> computePossibleSplitColumns(
    int firstColumn, 
    int lastColumnExcl,
    const int* countsA,
    const int* countsC,
    const int* countsG,
    const int* countsT,
    const int* coverages
){
    std::vector<MultipleSequenceAlignment::PossibleSplitColumn> possibleColumns;

    //find columns in which two different nucleotides each make up approx 50% (40% - 60%) of the column
    for(int col = firstColumn; col < lastColumnExcl; col++){
        std::array<MultipleSequenceAlignment::PossibleSplitColumn, 4> array{};
        int numPossibleNucs = 0;

        auto checkNuc = [&](const auto& counts, const char nuc){
            const float ratio = float(counts[col]) / float(coverages[col]);
            //if((counts[col] == 2 && fgeq(ratio, 0.4f) && fleq(ratio, 0.6f)) || counts[col] > 2){
            if(counts[col] >= 2 && fgeq(ratio, 0.4f) && fleq(ratio, 0.6f)){
                array[numPossibleNucs] = {nuc, col, ratio};
                numPossibleNucs++;
            }
        };

        checkNuc(countsA, 'A');
        checkNuc(countsC, 'C');
        checkNuc(countsG, 'G');
        checkNuc(countsT, 'T');

        
        if(numPossibleNucs == 2){
            possibleColumns.insert(possibleColumns.end(), array.begin(), array.begin() + numPossibleNucs);
        }
    }

    assert(possibleColumns.size() % 2 == 0);

    return possibleColumns;
}

MultipleSequenceAlignment::PossibleMsaSplits inspectColumnsRegionSplit(
    const MultipleSequenceAlignment::PossibleSplitColumn* possibleColumns,
    int numPossibleColumns,
    int /*firstColumn*/, 
    int /*lastColumnExcl*/,
    int anchorColumnsBegin_incl,
    int numCandidates,
    const char* candidates,
    int decodedSequencePitchBytes,
    const int* candidateShifts,
    const int* candidateLengths
){   

    if(2 <= numPossibleColumns && numPossibleColumns <= 32 && numPossibleColumns % 2 == 0){

        MultipleSequenceAlignment::PossibleMsaSplits result;

        //calculate proper results
        {
            

            //for each candidate check if it appears in the selected columns
            //there can be at most 16 selected columns.the result is bit-encoded in 32 bits, 2 bits per column
            //0b10 - candidate base matches first selected nuc in column
            //0b11 - candidate base matches second selected nuc in column
            //0b00 - candidate base does not match, or candidate does not occupy column

            //map groups candidates with the same bit pattern
            std::map<unsigned int, std::vector<int>> map;
            
            for(int l = 0; l < numCandidates; l++){
                const int candidateRow = l;

                constexpr int numPossibleColumnsPerFlag = sizeof(unsigned int) * CHAR_BIT / 2;

                unsigned int flags = 0;

                const int numColumnsToCheck = std::min(numPossibleColumnsPerFlag, numPossibleColumns / 2);
                const char* const candidateString = &candidates[candidateRow * decodedSequencePitchBytes];

                for(int k = 0; k < numColumnsToCheck; k++){
                    flags <<= 2;

                    const MultipleSequenceAlignment::PossibleSplitColumn psc0 = possibleColumns[2*k+0];
                    const MultipleSequenceAlignment::PossibleSplitColumn psc1 = possibleColumns[2*k+1];
                    assert(psc0.column == psc1.column);

                    const int candidateColumnsBegin_incl = candidateShifts[candidateRow] + anchorColumnsBegin_incl;
                    const int candidateColumnsEnd_excl = candidateLengths[candidateRow] + candidateColumnsBegin_incl;
                    
                    //column range check for row
                    if(candidateColumnsBegin_incl <= psc0.column && psc0.column < candidateColumnsEnd_excl){
                        const int positionInCandidate = psc0.column - candidateColumnsBegin_incl;

                        if(candidateString[positionInCandidate] == psc0.letter){
                            flags = flags | 0b10;
                        }else if(candidateString[positionInCandidate] == psc1.letter){
                            flags = flags | 0b11;
                        }else{
                            flags = flags | 0b00;
                        } 

                    }else{
                        flags = flags | 0b00;
                    } 

                }

                map[flags].emplace_back(l);
            }

            std::vector<std::pair<unsigned int, std::vector<int>>> flatmap(map.begin(), map.end());

            std::map<unsigned int, std::vector<int>> finalMap;

            const int flatmapsize = flatmap.size();
            for(int i = 0; i < flatmapsize; i++){
                //try to merge flatmap[i] with flatmap[k], i < k, if possible
                const unsigned int flagsToSearch = flatmap[i].first;
                unsigned int mask = 0;
                for(int s = 0; s < 16; s++){
                    if(flagsToSearch >> (2*s+1) & 1){
                        mask = mask | (0x03 << (2*s));
                    }
                }

                //std::cerr << "i = " << i << ", flags = " << std::bitset<32>(flagsToSearch) << ", mask = " << std::bitset<32>(mask) << "\n";

                bool merged = false;
                for(int k = i+1; k < flatmapsize; k++){
                    //if both columns are identical not including wildcard columns
                    if((mask & flatmap[k].first) == flagsToSearch){
                        //std::cerr << "k = " << k << ", flags = " << std::bitset<32>(flatmap[k].first) << " equal" << "\n";
                        flatmap[k].second.insert(
                            flatmap[k].second.end(),
                            flatmap[i].second.begin(),
                            flatmap[i].second.end()
                        );

                        std::sort(flatmap[k].second.begin(), flatmap[k].second.end());

                        flatmap[k].second.erase(
                            std::unique(flatmap[k].second.begin(), flatmap[k].second.end()),
                            flatmap[k].second.end()
                        );

                        merged = true;
                    }else{
                        //std::cerr << "k = " << k << ", flags = " << std::bitset<32>(flatmap[k].first) << " not equal" << "\n";
                    }
                }

                if(!merged){
                    finalMap[flatmap[i].first] = std::move(flatmap[i].second);
                }
            }

            for(auto& pair : finalMap){

                std::vector<int> listOfCandidates = std::move(pair.second);
                std::vector<MultipleSequenceAlignment::PossibleSplitColumn> columnInfo;

                const unsigned int flag = pair.first;
                const int num = numPossibleColumns / 2;
                for(int i = 0; i < num; i++){
                    const unsigned int cur = (flag >> (num - i - 1) * 2) & 0b11;
                    const bool match = (cur & 0b10) == 0b10;
                    if(match){
                        const int which = cur & 1;
                        columnInfo.emplace_back(possibleColumns[2*i + which]);
                    }
                }

                result.splits.emplace_back(std::move(columnInfo), std::move(listOfCandidates));                
            }
        }

#if 0
        //calculate sorted and debug prints
        {

            std::map<unsigned int, std::vector<int>> debugprintmap;
            std::vector<int> sortedindices(numCandidates);
            std::iota(sortedindices.begin(), sortedindices.end(), 0);

            auto get_shift_of_row = [&](int k){
                return candidateShifts[k];
            };

            std::sort(sortedindices.begin(), sortedindices.end(),
                    [&](int l, int r){return get_shift_of_row(l) < get_shift_of_row(r);});

            for(int l = 0; l < numCandidates; l++){
                const int candidateRow = sortedindices[l];

                constexpr std::size_t numPossibleColumnsPerFlag = sizeof(unsigned int) * CHAR_BIT / 2;

                unsigned int flags = 0;

                const int numColumnsToCheck = std::min(numPossibleColumnsPerFlag, possibleColumns.size() / 2);
                const char* const candidateString = &candidates[candidateRow * decodedSequencePitchBytes];

                for(int k = 0; k < numColumnsToCheck; k++){
                    flags <<= 2;

                    const PossibleSplitColumn psc0 = possibleColumns[2*k+0];
                    const PossibleSplitColumn psc1 = possibleColumns[2*k+1];
                    assert(psc0.column == psc1.column);

                    const int candidateColumnsBegin_incl = candidateShifts[candidateRow] + anchorColumnsBegin_incl;
                    const int candidateColumnsEnd_excl = candidateLengths[candidateRow] + candidateColumnsBegin_incl;
                    
                    //column range check for row
                    if(candidateColumnsBegin_incl <= psc0.column && psc0.column < candidateColumnsEnd_excl){
                        const int positionInCandidate = psc0.column - candidateColumnsBegin_incl;

                        if(candidateString[positionInCandidate] == psc0.letter){
                            flags = flags | 0b10;
                        }else if(candidateString[positionInCandidate] == psc1.letter){
                            flags = flags | 0b11;
                        }else{
                            flags = flags | 0b00;
                        } 

                    }else{
                        flags = flags | 0b00;
                    } 

                }

                debugprintmap[flags].emplace_back(l);
            }

            std::vector<std::pair<unsigned int, std::vector<int>>> flatmap(debugprintmap.begin(), debugprintmap.end());

            std::map<unsigned int, std::vector<int>> finalMap;

            const int flatmapsize = flatmap.size();
            for(int i = 0; i < flatmapsize; i++){
                //try to merge flatmap[i] with flatmap[k], i < k, if possible
                const unsigned int flagsToSearch = flatmap[i].first;
                unsigned int mask = 0;
                for(int s = 0; s < 16; s++){
                    if(flagsToSearch >> (2*s+1) & 1){
                        mask = mask | (0x03 << (2*s));
                    }
                }

                //std::cerr << "i = " << i << ", flags = " << std::bitset<32>(flagsToSearch) << ", mask = " << std::bitset<32>(mask) << "\n";

                bool merged = false;
                for(int k = i+1; k < flatmapsize; k++){
                    //if both columns are identical not including wildcard columns
                    if((mask & flatmap[k].first) == flagsToSearch){
                        //std::cerr << "k = " << k << ", flags = " << std::bitset<32>(flatmap[k].first) << " equal" << "\n";
                        flatmap[k].second.insert(
                            flatmap[k].second.end(),
                            flatmap[i].second.begin(),
                            flatmap[i].second.end()
                        );

                        std::sort(flatmap[k].second.begin(), flatmap[k].second.end());

                        flatmap[k].second.erase(
                            std::unique(flatmap[k].second.begin(), flatmap[k].second.end()),
                            flatmap[k].second.end()
                        );

                        merged = true;
                    }else{
                        //std::cerr << "k = " << k << ", flags = " << std::bitset<32>(flatmap[k].first) << " not equal" << "\n";
                    }
                }

                if(!merged){
                    finalMap[flatmap[i].first] = std::move(flatmap[i].second);
                }
            }

            auto printMap = [&](const auto& map){
                for(const auto& pair : map){
                    const unsigned int flag = pair.first;
                    //convert flag to position and nuc
                    const int num = numPossibleColumns / 2;

                    std::cerr << "flag " << flag << " : ";
                    for(int i = 0; i < num; i++){
                        const unsigned int cur = (flag >> (num - i - 1) * 2) & 0b11;
                        const bool match = (cur & 0b10) == 0b10;
                        const int column = possibleColumns[2*i].column;
                        char nuc = '-';
                        if(match){
                            int which = cur & 1;
                            nuc = possibleColumns[2*i + which].letter;
                        }
                        std::cerr << "(" << nuc << ", " << column << ") ";
                    }
                    std::cerr << ": ";

                    for(int c : pair.second){
                        std::cerr << c << " ";
                    }
                    std::cerr << "\n";
                }
            };

            if(numPossibleColumns > 0){
                std::cerr << numPossibleColumns << "\n";
                for(const auto& p : possibleColumns){
                    std::cerr << "{" << p.letter << ", " << p.column << ", " << p.ratio << "} ";
                }
                std::cerr << "\n";
                
                printMap(debugprintmap);

                std::cerr << "final map: \n";

                printMap(finalMap);
                std::cerr << "\n";

                // print(std::cerr);
                // std::cerr << "\n";
            }
        }
#endif    
        return result;
    }else{
        // single split with all candidates
        std::vector<int> listOfCandidates(numCandidates);
        std::iota(listOfCandidates.begin(), listOfCandidates.end(), 0);

        MultipleSequenceAlignment::PossibleMsaSplits result;

        std::vector<MultipleSequenceAlignment::PossibleSplitColumn> columnInfo;


        result.splits.emplace_back(std::move(columnInfo), std::move(listOfCandidates));
        
        return result;
    }
}


}