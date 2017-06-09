#include <iostream>
#include <seqan/file.h>
#include <seqan/bam_io.h>
#include <set>
#include <map>
#include <string>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>

using namespace std;
using namespace seqan;

struct MateEditInfo {
    int beginPosShift = 0;
    int fragLenChange = 0;
    bool mateRemoved = false;
    bool matePrinted = false;
} ;

struct DeletionStats {
    int nSoftClippedBp = 0;
    int nQualityClippedBp = 0;
    int nAdapterClippedBp = 0;
    int nMatchRemovedReads = 0;
    int nAdapterReads = 0;
    unsigned nTotalReads = 0;
    unsigned nCoverageFiltered = 0;
} ;

DeletionStats delStats;

void removeBeginReads(map<unsigned, String<BamAlignmentRecord> >& beginPosToReads, map<CharString, Pair<MateEditInfo> >& mateEditMap, unsigned readyPos, unsigned start, unsigned end)
{
    bool sameOrientation;
    map<unsigned, String<BamAlignmentRecord> >::const_iterator it = beginPosToReads.begin();
    while (it->first<=readyPos && it->first < start && it!=beginPosToReads.end())
    {
        String<BamAlignmentRecord> & readsToCheck = beginPosToReads[it->first];
        String<BamAlignmentRecord> readsToKeep;
        for (unsigned i=0; i<length(readsToCheck); ++i)
        {
        	BamAlignmentRecord & record = readsToCheck[i];
            sameOrientation = hasFlagRC(record) == hasFlagNextRC(record);
            if (hasFlagNextUnmapped(record) || hasFlagUnmapped(record))
                continue;
            MateEditInfo editInfo;
            if (sameOrientation)
                continue;
            else
            {
                if (!hasFlagRC(record))
                    editInfo = mateEditMap[record.qName].i2;
                else
                    editInfo = mateEditMap[record.qName].i1;
            }
            int fragLen, pNextAlnLen;
            if (record.tLen>0)
            {
                fragLen = editInfo.fragLenChange - record.beginPos;
                pNextAlnLen = fragLen - (record.pNext-record.beginPos);
            }
            else
                pNextAlnLen = getAlignmentLengthInRef(record);
            if ((record.pNext + pNextAlnLen < start && record.beginPos + getAlignmentLengthInRef(record) < start) || editInfo.mateRemoved || (record.pNext > end && record.beginPos > end))
            {
                mateEditMap[record.qName].i2.mateRemoved = true;
                mateEditMap[record.qName].i1.mateRemoved = true;
                continue;
            }
            else
                appendValue(readsToKeep, record);
        }
        beginPosToReads[it->first] = readsToKeep;
        ++it;
    }
}

void removeUnPairReads(map<unsigned, String<BamAlignmentRecord> >& beginPosToReads, map<CharString, Pair<MateEditInfo> >& mateEditMap)
{
    bool sameOrientation;
    map<unsigned, String<BamAlignmentRecord> >::const_iterator itEnd = beginPosToReads.end();
    for (map<unsigned, String<BamAlignmentRecord> >::const_iterator it = beginPosToReads.begin(); it != itEnd; ++it)
    {
        String<BamAlignmentRecord> & readsToCheck = beginPosToReads[it->first];
        String<BamAlignmentRecord> readsToKeep;
        for (unsigned i=0; i<length(readsToCheck); ++i)
        {
            BamAlignmentRecord & record = readsToCheck[i];
            if (hasFlagNextUnmapped(record) || hasFlagUnmapped(record))
                continue;
            sameOrientation = hasFlagRC(record) == hasFlagNextRC(record);
            MateEditInfo editInfo;
            if (sameOrientation)
                continue;
            else
            {
                if (!hasFlagRC(record))
                    editInfo = mateEditMap[record.qName].i2;
                else
                    editInfo = mateEditMap[record.qName].i1;
            }
            if (editInfo.matePrinted)
                appendValue(readsToKeep, record);
        }
        beginPosToReads[it->first] = readsToKeep;
    }
    return;
}

void removeTags(BamAlignmentRecord& record, bool keepMapQual)
{
    BamTagsDict tagsDict(record.tags);
    unsigned numOfTags = length(tagsDict);
    String<CharString> keysToErase;
    for (unsigned i=0; i<numOfTags; ++i)
    {
        CharString key = getTagKey(tagsDict, i);
        if (!(key == "RG" || (key == "MQ" && keepMapQual)))
            appendValue(keysToErase, key);
    }
    for (unsigned i=0; i<length(keysToErase); ++i)
        eraseTag(tagsDict, keysToErase[i]);
    return;
}

void makeUnpaired(BamAlignmentRecord& record, bool keepMapQual)
{
    record.tLen = 0;
    record.pNext = record.INVALID_POS;
    record.rNextId = record.INVALID_REFID;
    //unset: FlagNextUnmapped, FlagAllProper,FlagMultiple, FlagNextRC:
    record.flag &= ~BAM_FLAG_NEXT_UNMAPPED;
    record.flag &= ~BAM_FLAG_ALL_PROPER;
    record.flag &= ~BAM_FLAG_MULTIPLE;
    record.flag &= ~BAM_FLAG_NEXT_RC;
}

void printReadyReads(map<CharString, Pair<MateEditInfo> >& mateEditMap, map<unsigned, String<BamAlignmentRecord> >& beginPosToReads, unsigned readyPos, BamFileOut& bamFileOut, map<CharString, unsigned>& readNameToNum, unsigned& currIdx, bool keepMapQual, Pair<unsigned> start_end)
{
    bool sameOrientation;
    map<unsigned, String<BamAlignmentRecord> >::const_iterator it = beginPosToReads.begin();
    while (it->first<=readyPos && it!=beginPosToReads.end())
    {
        String<BamAlignmentRecord> & readsToWrite = beginPosToReads[it->first];
        for (unsigned i=0; i<length(readsToWrite); ++i)
        {
        	BamAlignmentRecord & record = readsToWrite[i];
            MateEditInfo editInfo;
            sameOrientation = hasFlagRC(record) == hasFlagNextRC(record);
            if (hasFlagMultiple(record) && !sameOrientation)
            {
                if (!hasFlagRC(record))
                {
                    editInfo = mateEditMap[record.qName].i2;
                    if (editInfo.matePrinted)
                        mateEditMap.erase(record.qName);
                    else
                        mateEditMap[record.qName].i1.matePrinted = true;
                }
                else
                {
                    editInfo = mateEditMap[record.qName].i1;
                    if (editInfo.matePrinted)
                        mateEditMap.erase(record.qName);
                    else
                        mateEditMap[record.qName].i2.matePrinted = true;
                }
                if (!editInfo.matePrinted && record.beginPos > record.pNext)
                    continue;
                if (editInfo.mateRemoved)
                {
                    makeUnpaired(record, keepMapQual);
                    mateEditMap.erase(record.qName);
                    if (hasFlagUnmapped(record))
                        continue;
                }
                else
                {
                    record.pNext += editInfo.beginPosShift;
                    if (!hasFlagNextUnmapped(record) && !hasFlagUnmapped(record))
                    {
                        if (record.tLen>0)
                            record.tLen = editInfo.fragLenChange - record.beginPos;
                        else
                            record.tLen = -1 * (record.beginPos + getAlignmentLengthInRef(record) - editInfo.fragLenChange);
                    }
                    else
                        record.tLen = 0;
                }
            }
            else
            {
                mateEditMap.erase(record.qName);
                makeUnpaired(record, keepMapQual);
            }
            if (readNameToNum.count(record.qName) == 0)
            {
                if (hasFlagMultiple(record))
                    readNameToNum[record.qName] = currIdx;
                stringstream ss;
                ss << currIdx;
                CharString str = ss.str();
                record.qName = str;
                ++currIdx;
            }
            else
            {
                stringstream ss;
                ss << readNameToNum[record.qName];
                CharString str = ss.str();
                record.qName = str;
                readNameToNum.erase(record.qName);
            }
            removeTags(record, keepMapQual);
            writeRecord(bamFileOut, record);
        }
        beginPosToReads.erase(it->first);
        it = beginPosToReads.begin();
    }
    return;
}

unsigned countMatchingBases(String<CigarElement<> >& cigarString)
{
    unsigned numOfMatches = 0;
    for (unsigned i=0; i<length(cigarString); ++i)
    {
        CharString cigarOperation = cigarString[i].operation;
        string cigarOperationStr = toCString(cigarOperation);
        if (cigarOperationStr.compare("M")==0)
            numOfMatches += cigarString[i].count;
    }
    return numOfMatches;
}

bool cigarAndSeqMatch(BamAlignmentRecord& record)
{
    String<CigarElement<> >& cigarString = record.cigar;
    unsigned counter = 0;
    for (unsigned i=0; i<length(cigarString); ++i)
    {
        CharString cigarOperation = cigarString[i].operation;
        string cigarOperationStr = toCString(cigarOperation);
        if (cigarOperationStr.compare("D")!=0)
            counter += cigarString[i].count;
    }
    if (length(record.seq)!=counter)
        return false;
    return true;
}

void resetCigarStringEnd(BamAlignmentRecord& record, unsigned nRemoved, String<CigarElement<> >& removed)
{
    String<CigarElement<> >& cigarString = record.cigar;
    CharString cigarOperation = cigarString[length(cigarString)-1].operation;
    string cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("D")==0)
    {
        appendValue(removed, cigarString[length(cigarString)-1]);
        erase(cigarString, length(cigarString)-1);
    }
    if (cigarString[length(cigarString)-1].count > nRemoved)
    {
        CigarElement<> t = cigarString[length(cigarString)-1];
        t.count = nRemoved;
        appendValue(removed, t);
        cigarString[length(cigarString)-1].count -= nRemoved;
        return;
    }
    else
    {
        if (cigarString[length(cigarString)-1].count == nRemoved)
        {
            appendValue(removed, cigarString[length(cigarString)-1]);
            erase(cigarString, length(cigarString)-1);
            cigarOperation = cigarString[length(cigarString)-1].operation;
            cigarOperationStr = toCString(cigarOperation);
            if (cigarOperationStr.compare("D")==0)
            {
                appendValue(removed, cigarString[length(cigarString)-1]);
                erase(cigarString, length(cigarString)-1);
            }
            return;
        }
        else
        {
            unsigned nLeft = nRemoved - cigarString[length(cigarString)-1].count;
            appendValue(removed, cigarString[length(cigarString)-1]);
            erase(cigarString, length(cigarString)-1);
            resetCigarStringEnd(record, nLeft, removed);
        }
    }
}

bool qualityClipEnd(BamAlignmentRecord& record, int windowSize, map<CharString, Pair<MateEditInfo> >& mateEditMap)
{
    bool sameOrientation = hasFlagRC(record) == hasFlagNextRC(record);
    int qualSum = 0;
    for (unsigned i = length(record.qual)-1; i>=length(record.qual)-windowSize; --i)
        qualSum += (record.qual[i]-33);
    int averageQual = round((double)qualSum/(double)windowSize);
    if (averageQual>=25)
        return true;
    int index = 0;
    while (averageQual<25 && index < length(record.qual))
    {
        ++index;
        qualSum -= (record.qual[length(record.qual)-index]-33);
        qualSum += (record.qual[length(record.qual)-windowSize-index]-33);
        averageQual = round((double)qualSum/(double)windowSize);
    }
    delStats.nQualityClippedBp += index;
    erase(record.seq, length(record.seq)-index, length(record.seq));
    erase(record.qual,length(record.qual)-index, length(record.qual));
    if (hasFlagRC(record) && !sameOrientation)
    {
        if (record.tLen > 0)
            record.tLen -= index;
        else
            record.tLen += index;
    }
    String<CigarElement<> > removedCigar;
    resetCigarStringEnd(record, index, removedCigar);
    if (length(record.seq)>=50)
    {
        if (!sameOrientation)
        {
            if (hasFlagRC(record))
                mateEditMap[record.qName].i2.fragLenChange = record.beginPos + getAlignmentLengthInRef(record);
        }
        return true;
    }
    else
    {
        if (!sameOrientation)
        {
            if (hasFlagRC(record))
                mateEditMap[record.qName].i2.mateRemoved = true;
            else
                mateEditMap[record.qName].i1.mateRemoved = true;
        }
        return false;
    }
}

void resetCigarStringBegin(BamAlignmentRecord& record, unsigned nRemoved, String<CigarElement<> >& removed)
{
    String<CigarElement<> >& cigarString = record.cigar;
    CharString cigarOperation = cigarString[0].operation;
    string cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("D")==0)
    {
        appendValue(removed, cigarString[0]);
        erase(cigarString, 0);
    }
    if (cigarString[0].count > nRemoved)
    {
        CigarElement<> t = cigarString[0];
        t.count = nRemoved;
        appendValue(removed, t);
        cigarString[0].count -= nRemoved;
        return;
    }
    else
    {
        if (cigarString[0].count == nRemoved)
        {
            appendValue(removed, cigarString[0]);
            erase(cigarString, 0);
            cigarOperation = cigarString[0].operation;
            cigarOperationStr = toCString(cigarOperation);
            if (cigarOperationStr.compare("D")==0)
            {
                appendValue(removed, cigarString[0]);
                erase(cigarString, 0);
            }
            return;
        }
        else
        {
            unsigned nLeft = nRemoved - cigarString[0].count;
            appendValue(removed, cigarString[0]);
            erase(cigarString, 0);
            resetCigarStringBegin(record, nLeft, removed);
        }
    }
}

bool qualityClipBegin(BamAlignmentRecord& record, int windowSize, map<CharString, Pair<MateEditInfo> >& mateEditMap)
{
    bool sameOrientation = hasFlagRC(record) == hasFlagNextRC(record);
    int qualSum = 0;
    for (unsigned i = 0; i<windowSize; ++i)
        qualSum += (record.qual[i]-33);
    int averageQual = round((double)qualSum/(double)windowSize);
    if (averageQual>=25)
        return true;
    unsigned index = 0;
    while (averageQual<25 && index+windowSize < length(record.qual))
    {
        qualSum -= (record.qual[index]-33);
        qualSum += (record.qual[index+windowSize]-33);
        averageQual = round((double)qualSum/(double)windowSize);
        ++index;
    }
    delStats.nQualityClippedBp += index;
    erase(record.seq, 0, index);
    erase(record.qual, 0, index);
    CharString cigarOperation = record.cigar[0].operation;
    string cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("I")==0)
    {
        int shift = index - record.cigar[0].count;
        if (shift > 0)
            record.beginPos += shift;
    }
    else
        record.beginPos += index;
    if (!hasFlagRC(record) && !sameOrientation)
    {
        if (record.tLen > 0)
            record.tLen -= index;
        else
            record.tLen += index;
    }
    String<CigarElement<> > removedCigar;
    resetCigarStringBegin(record, index, removedCigar);
    if (length(record.seq)>=50)
    {
        if (!sameOrientation)
        {
            if (hasFlagRC(record))
            {
                mateEditMap[record.qName].i2.beginPosShift += index;
            }
            else
            {
                mateEditMap[record.qName].i1.beginPosShift += index;
                mateEditMap[record.qName].i1.fragLenChange = record.beginPos;
            }
        }
        return true;
    }
    else
    {
        if (!sameOrientation)
        {
            if (hasFlagRC(record))
                mateEditMap[record.qName].i2.mateRemoved = true;
            else
                mateEditMap[record.qName].i1.mateRemoved = true;
        }
        return false;
    }
}

bool removeSoftClipped(BamAlignmentRecord& record, map<CharString, Pair<MateEditInfo> >& mateEditMap, unsigned minMatchingBases)
{
    bool sameOrientation = hasFlagRC(record) == hasFlagNextRC(record);
    String<CigarElement<> > cigarString = record.cigar;
    CharString readName = record.qName;
    CharString cigarOperation = cigarString[0].operation;
    string cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("S")==0)
    {
        delStats.nSoftClippedBp += cigarString[0].count;
        erase(record.seq, 0, cigarString[0].count);
        erase(record.qual, 0, cigarString[0].count);
        erase(record.cigar, 0);
    }
    cigarOperation = cigarString[length(cigarString)-1].operation;
    cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("S")==0)
    {
        delStats.nSoftClippedBp += cigarString[length(cigarString)-1].count;
        for (unsigned i = 0; i<cigarString[length(cigarString)-1].count; ++i)
        {
            eraseBack(record.seq);
            eraseBack(record.qual);
        }
        eraseBack(record.cigar);
    }
    if (length(record.seq)>=minMatchingBases)
        return true;
    else
    {
        if (!sameOrientation)
        {
            if (hasFlagRC(record))
                mateEditMap[readName].i2.mateRemoved = true;
            else
                mateEditMap[readName].i1.mateRemoved = true;
        }
        return false;
    }
}

bool removeNsAtEnds(BamAlignmentRecord& record, map<CharString, Pair<MateEditInfo> >& mateEditMap, unsigned minMatchingBases)
{
    CharString firstBaseC = record.seq[0];
    string firstBase = toCString(firstBaseC);
    CharString lastBaseC = record.seq[length(record.seq)-1];
    string lastBase = toCString(lastBaseC);
    int nOfNs = 0;
    if (firstBase.compare("N")==0)
    {
        ++nOfNs;
        int idx =1;
        CharString base2checkC = record.seq[idx];
        string base2check = toCString(base2checkC);
        while (base2check.compare("N")==0 && idx < length(record.seq)-1)
        {
            ++nOfNs;
            ++idx;
            base2checkC = record.seq[idx];
            base2check = toCString(base2checkC);
        }
        //Remove ns from beginning of sequence and qual fields:
        erase(record.seq, 0, nOfNs);
        erase(record.qual, 0, nOfNs);
        if (!hasFlagUnmapped(record)) //Only have to fix CIGAR, beginPos and fragLen if the read is mapped.
        {
            String<CigarElement<> > removedCigar;
            resetCigarStringBegin(record, nOfNs, removedCigar);
            int shift = 0;
            for (unsigned i=0; i<length(removedCigar);++i)
            {
                CharString cigarOperation = removedCigar[i].operation;
                string cigarOperationStr = toCString(cigarOperation);
                if (cigarOperationStr.compare("M")==0 || cigarOperationStr.compare("D")==0)
                    shift += removedCigar[i].count;
            }
            record.beginPos += shift;
            if (!hasFlagRC(record))
            {
                mateEditMap[record.qName].i1.beginPosShift += shift;
                mateEditMap[record.qName].i1.fragLenChange = record.beginPos;
            }
            else
                mateEditMap[record.qName].i2.beginPosShift += shift;
            //cout << "Removed: " << nOfNs << " Ns from beginning of " << record.qName << endl;
        }
    }
    if (length(record.seq) < minMatchingBases)
        return false;
    nOfNs = 0;
    if (lastBase.compare("N")==0)
    {
        ++nOfNs;
        int idx =length(record.seq)-2;
        CharString base2checkC = record.seq[idx];
        string base2check = toCString(base2checkC);
        while (base2check.compare("N")==0 && idx > 0)
        {
            ++nOfNs;
            --idx;
            base2checkC = record.seq[idx];
            base2check = toCString(base2checkC);
        }
        //Remove ns from end of sequence and qual fields:
        erase(record.seq, length(record.seq)-nOfNs, length(record.seq));
        erase(record.qual,length(record.qual)-nOfNs, length(record.qual));
        if (!hasFlagUnmapped(record)) //Only have to fix CIGAR, beginPos and fragLen if the read is mapped.
        {
            String<CigarElement<> > removedCigar;
            resetCigarStringEnd(record, nOfNs, removedCigar);
            if (hasFlagRC(record))
                mateEditMap[record.qName].i2.fragLenChange = record.beginPos + getAlignmentLengthInRef(record);
            //cout << "Removed: " << nOfNs << " Ns from end of " << record.qName << endl;
        }
    }
    if (length(record.seq) < minMatchingBases)
        return false;
    return true;
}

bool qualityFilterLevel2(BamAlignmentRecord& record, map<CharString, Pair<MateEditInfo> >& mateEditMap, unsigned minMatchingBases)
{
    /*if (!removeSoftClipped(record, mateEditMap))
        return false;
    if (!qualityClipBegin(record, 5, mateEditMap))
        return false;
    if (!qualityClipEnd(record, 5, mateEditMap))
        return false;*/
    if (!removeNsAtEnds(record, mateEditMap, minMatchingBases))
    {
        ++delStats.nMatchRemovedReads;
        mateEditMap[record.qName].i2.mateRemoved = true;
        mateEditMap[record.qName].i1.mateRemoved = true;
        return false;
    }
    if (!hasFlagUnmapped(record))
    {
        if (!cigarAndSeqMatch(record))
            cout << "THE CIGAR STRING AND READ-LENGTH DON'T MATCH: " << record.qName << endl;
        unsigned matchingBases = countMatchingBases(record.cigar);
        if (matchingBases < minMatchingBases)
        {
            ++delStats.nMatchRemovedReads;
            mateEditMap[record.qName].i2.mateRemoved = true;
            mateEditMap[record.qName].i1.mateRemoved = true;
            return false;
        }
    }
    return true;
}

Pair<int> findNum2Clip(BamAlignmentRecord& recordReverse, int forwardStartPos)
{
    int num2clip = 0, num2shift = 0;
    unsigned cigarIndex = 0, reverseStartPos = recordReverse.beginPos, n;
    CharString cigarOperation = recordReverse.cigar[cigarIndex].operation;
    string cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("S")==0)
    {
        num2clip = recordReverse.cigar[cigarIndex].count;
        ++cigarIndex;
    }
    for (unsigned i=cigarIndex; i<length(recordReverse.cigar); ++i)
    {
        cigarOperation = recordReverse.cigar[cigarIndex].operation;
        cigarOperationStr = toCString(cigarOperation);
        n=1;
        while (reverseStartPos < forwardStartPos && n <= recordReverse.cigar[cigarIndex].count)
        {
            if (cigarOperationStr.compare("D")!=0)
                ++num2clip;
            if (cigarOperationStr.compare("I")!=0)
                ++reverseStartPos;
            ++n;
        }
        if (reverseStartPos == forwardStartPos)
            break;
        if (reverseStartPos>forwardStartPos)
            cout << recordReverse.qName << " Reverse startpos has become bigger than forward, reverse beginPos: " << reverseStartPos << " forward beginPos is: " << forwardStartPos << endl;
        ++cigarIndex;
    }
    if (cigarOperationStr.compare("D")==0)
    {
        num2shift = recordReverse.cigar[cigarIndex].count - n + 1;
        if (num2shift > 0)
            cout << "Will shift reverse read because it starts with a deletion: " << recordReverse.qName << endl;
    }
    return Pair<int>(num2clip,num2shift);
}

bool removeAdapters(BamAlignmentRecord& recordForward, BamAlignmentRecord& recordReverse, map<CharString, Pair<MateEditInfo> >& mateEditMap, unsigned minMatchingBases)
{
    //Check for soft clipped bases at beginning of forward record.
    CharString cigarOperation = recordForward.cigar[0].operation;
    string cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("S")==0)
    {
        if (recordForward.beginPos - recordForward.cigar[0].count <= recordReverse.beginPos)
        {
            return true;
        }
    }
    if (!removeSoftClipped(recordForward, mateEditMap, minMatchingBases) || !removeSoftClipped(recordReverse, mateEditMap, minMatchingBases))
        return false;
    delStats.nAdapterReads += 2;
    int startPosDiff = recordForward.beginPos - recordReverse.beginPos;
    if (startPosDiff<0)
        return true;
    Pair<int> clipAndShift = findNum2Clip(recordReverse, recordForward.beginPos);
    int index = clipAndShift.i1;
    int shift = clipAndShift.i2;
    //erase from reverse read bases 0 to index
    erase(recordReverse.seq, 0, index);
    erase(recordReverse.qual, 0, index);
    String<CigarElement<> > removedCigar;
    resetCigarStringBegin(recordReverse, index, removedCigar);
    //erase from forward read bases from length(reverse.seq) to end
    delStats.nAdapterClippedBp += index;
    int forwardClip = index;
    if (length(recordForward.seq)>length(recordReverse.seq) && index>0)
    {
        forwardClip = length(recordForward.seq) - length(recordReverse.seq);
        erase(recordForward.seq, length(recordReverse.seq), length(recordForward.seq));
        erase(recordForward.qual, length(recordReverse.qual), length(recordForward.qual));
        delStats.nAdapterClippedBp += forwardClip;
        String<CigarElement<> > removedCigar;
        resetCigarStringEnd(recordForward, forwardClip, removedCigar);
    }
    // if (forwardClip != index)
    // {
    //     cout << "Did not clip the same number of bases from both reads! " << recordReverse.qName << endl;
    //     cout << "Clipped " << index << "bp from beginning of reverse read and " << forwardClip << " from the end of the forward read." << endl;
    // }
    recordReverse.beginPos = recordForward.beginPos;
    if (shift > 0)
        recordReverse.beginPos += shift;
    recordForward.pNext = recordReverse.beginPos;
    mateEditMap[recordReverse.qName].i2.fragLenChange = recordReverse.beginPos + getAlignmentLengthInRef(recordReverse);
    mateEditMap[recordForward.qName].i1.fragLenChange = recordForward.beginPos;
    if (!cigarAndSeqMatch(recordForward))
        cout << "The cigar string and sequence length don't match for the forward read!!" << endl;
    if (!cigarAndSeqMatch(recordReverse))
        cout << "The cigar string and sequence length don't match for the reverse read!!" << endl;
    if (length(recordForward.seq)>=minMatchingBases)
        return true;
    else
        return false;
}

bool qualityFilter(BamAlignmentRecord& record, map<CharString, Pair<MateEditInfo> >& mateEditMap, int maxFragmentLength, map<CharString, BamAlignmentRecord>& adapterMap, unsigned minMatchingBases, bool keepMapQual, map<unsigned, String<BamAlignmentRecord> >& beginPosToReads)
{
    MateEditInfo editInfo;
    BamAlignmentRecord reverseRecord;
    bool sameOrientation = hasFlagRC(record) == hasFlagNextRC(record);
    if (hasFlagNextUnmapped(record) && sameOrientation)
    {
        //If mate is unmapped and both have same orientation, need to flip mate orientation
        if (hasFlagRC(record))
            record.flag |= ~BAM_FLAG_NEXT_RC;
        else
            record.flag |= BAM_FLAG_NEXT_RC;
        sameOrientation = false;
    }
    if (hasFlagUnmapped(record) && sameOrientation)
    {
        //If read is unmapped and both have same orientation, need to flip read orientation
        if (hasFlagRC(record))
            record.flag |= ~BAM_FLAG_RC;
        else
            record.flag |= BAM_FLAG_RC;
        sameOrientation = false;
    }
    if (hasFlagRC(record))
    {
        if (mateEditMap.count(record.qName)==0)
            mateEditMap[record.qName].i2 = editInfo;
        mateEditMap[record.qName].i2.fragLenChange = record.beginPos + getAlignmentLengthInRef(record);
    }
    else
    {
        if (mateEditMap.count(record.qName)==0)
            mateEditMap[record.qName].i1 = editInfo;
        mateEditMap[record.qName].i1.fragLenChange = record.beginPos;
    }
    if (sameOrientation)
    {
        makeUnpaired(record, keepMapQual);
        mateEditMap[record.qName].i2.mateRemoved = true;
        mateEditMap[record.qName].i1.mateRemoved = true;
    }
    //if (hasFlagDuplicate(record) || hasFlagUnmapped(record))
    if (hasFlagDuplicate(record))
    {
        if (hasFlagRC(record))
            mateEditMap[record.qName].i2.mateRemoved = true;
        else
            mateEditMap[record.qName].i1.mateRemoved = true;
        return false;
    }
    //if (hasFlagNextUnmapped(record) || abs(record.tLen) > maxFragmentLength || record.rID != record.rNextId)
    if (abs(record.tLen) > maxFragmentLength || record.rID != record.rNextId)
    {
        makeUnpaired(record, keepMapQual);
        mateEditMap[record.qName].i2.mateRemoved = true;
        mateEditMap[record.qName].i1.mateRemoved = true;
    }
    if (abs(record.tLen)<length(record.seq) && record.rID == record.rNextId && adapterMap.count(record.qName)==0 && !hasFlagNextUnmapped(record) && !hasFlagUnmapped(record))
    {
            adapterMap[record.qName] = record;
            return false;
    }
    if (abs(record.tLen)<length(record.seq) && record.rID == record.rNextId && adapterMap.count(record.qName)!=0)
    {
        if (!sameOrientation)
        {
            if (!hasFlagRC(record))
                reverseRecord = adapterMap[record.qName];
            else
            {
                reverseRecord = record;
                record = adapterMap[record.qName];
            }
            adapterMap.erase(record.qName);
            if (!removeAdapters(record,reverseRecord, mateEditMap, minMatchingBases))
                return false;
        }
        else
            reverseRecord = adapterMap[record.qName];
        if (qualityFilterLevel2(reverseRecord, mateEditMap, minMatchingBases))
        {
            CharString binary_qual = "";
            for (unsigned x=0; x<length(reverseRecord.qual); ++x)
            {
                if (reverseRecord.qual[x]-33 >= 25)
                    append(binary_qual,"I");
                else
                    append(binary_qual,"!");
            }
            reverseRecord.qual = binary_qual;
            appendValue(beginPosToReads[reverseRecord.beginPos], reverseRecord);
        }
        else
            mateEditMap[record.qName].i2.mateRemoved = true;
    }
    if (!qualityFilterLevel2(record, mateEditMap, minMatchingBases))
        return false;
    return true;
}

void removeHardClipped(BamAlignmentRecord& record)
{
    //cout << "Working on record: " << record.qName << endl;
    if (hasFlagUnmapped(record))
        return;
    CharString cigarOperation = record.cigar[0].operation;
    string cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("H")==0)
        erase(record.cigar, 0);
    cigarOperation = record.cigar[length(record.cigar)-1].operation;
    cigarOperationStr = toCString(cigarOperation);
    if (cigarOperationStr.compare("H")==0)
        erase(record.cigar, length(record.cigar)-1);
}

int qualityFilterSlice(Triple<CharString, int, int >& chr_start_end, CharString baiPathIn, BamFileIn& bamFileIn, BamFileOut& bamFileOut, map<CharString, Pair<MateEditInfo> >& mateEditMap, bool keepMapQual, int maxFragLen, map<unsigned, String<BamAlignmentRecord> >& beginPosToReads, map<CharString, BamAlignmentRecord>& adapterMap, unsigned minMatchingBases, double avgCovByReadLen, map<CharString, unsigned>& readNameToNum, unsigned& currIdx)
{
    double maxQueSum = avgCovByReadLen*(double)50.0*(double)3.0;
    unsigned currBeginPos = 0;
    deque<unsigned> myQue;
    BamIndex<Bai> baiIndex;
    if (!open(baiIndex, toCString(baiPathIn)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << baiPathIn << "\n";
        return 1;
    }
    int rID = 0;
    if (!getIdByName(rID, contigNamesCache(context(bamFileIn)), chr_start_end.i1))
    {
        std::cerr << "ERROR: Reference sequence named " << chr_start_end.i1 << " not known.\n";
        return 1;
    }
    bool hasAlignments = false;
    if (!jumpToRegion(bamFileIn, hasAlignments, rID, std::max((int)0,(int)chr_start_end.i2-maxFragLen), chr_start_end.i3+maxFragLen, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << chr_start_end.i2 << ":" << chr_start_end.i3 << "\n";
        return 1;
    }
    if (!hasAlignments)
    {
        cout << "No alignments found in the interval: " << std::max((int)0,(int)chr_start_end.i2-maxFragLen) << " to " << chr_start_end.i3+maxFragLen << "\n";
        return 0;
    }
    BamAlignmentRecord record;
    CharString binary_qual;
    while (!atEnd(bamFileIn) )
    {
        readRecord(record, bamFileIn);
        //cout << "Processing read: " << record.qName << " at:" << record.beginPos << endl;
        if (record.rID == -1 || record.rID > rID || record.beginPos > chr_start_end.i3+maxFragLen)
            break;
        if (record.beginPos < std::max((int)0,chr_start_end.i2-maxFragLen))
            continue;
        removeHardClipped(record);
        ++delStats.nTotalReads;

        if (record.beginPos != currBeginPos)
        {
            if (record.beginPos-currBeginPos > 1)
            {
                for (unsigned i=1; i<record.beginPos-currBeginPos; ++i)
                {
                    myQue.push_front(0);
                    while (myQue.size()>50)
                        myQue.pop_back();
                }
            }
            myQue.push_front(1);
            while (myQue.size()>50)
                myQue.pop_back();
            currBeginPos = record.beginPos;
        }
        else
            ++myQue.front();
        if (std::accumulate(myQue.begin(),myQue.end(),0) > maxQueSum)
        {
            ++delStats.nCoverageFiltered;
            --myQue.front();
            if (hasFlagRC(record))
                mateEditMap[record.qName].i2.mateRemoved = true;
            else
                mateEditMap[record.qName].i1.mateRemoved = true;
            continue;
        }
        if (qualityFilter(record, mateEditMap, maxFragLen, adapterMap, minMatchingBases, keepMapQual, beginPosToReads))
        {
            binary_qual = "";
            for (unsigned x=0; x<length(record.qual); ++x)
            {
                if (record.qual[x]-33 >= 25)
                    append(binary_qual,"I");
                else
                    append(binary_qual,"!");
            }
            record.qual = binary_qual;
            appendValue(beginPosToReads[record.beginPos], record);
            if (record.beginPos-maxFragLen >=0)
            {
                //If I am printing reads infront of the interval I need to remove ones that don't have a mate in the interval first.
                if (beginPosToReads.begin()->first < chr_start_end.i2)
                    removeBeginReads(beginPosToReads, mateEditMap, record.beginPos-maxFragLen, chr_start_end.i2, chr_start_end.i3);
                printReadyReads(mateEditMap, beginPosToReads, record.beginPos-maxFragLen, bamFileOut, readNameToNum, currIdx, keepMapQual, Pair<unsigned>(chr_start_end.i2, chr_start_end.i3));
            }
        }
        else
            --myQue.front();
    }
    printReadyReads(mateEditMap, beginPosToReads, chr_start_end.i3, bamFileOut, readNameToNum, currIdx, keepMapQual, Pair<unsigned>(chr_start_end.i2, chr_start_end.i3));
    removeUnPairReads(beginPosToReads, mateEditMap);
    printReadyReads(mateEditMap, beginPosToReads, record.beginPos, bamFileOut, readNameToNum, currIdx, keepMapQual, Pair<unsigned>(chr_start_end.i2, chr_start_end.i3));
    if (beginPosToReads.size()>0)
    {
        for (map<unsigned, String<BamAlignmentRecord> >::const_iterator it = beginPosToReads.begin(); it != beginPosToReads.end(); ++it)
        {
            if (length(beginPosToReads[it->first])>0)
                cout << "There are " << length(beginPosToReads[it->first]) << " reads left at: " << it->first << endl;
        }
    }
    return 0;
}

String<Triple<CharString, int, int > > readIntervals(CharString& intervalFile, int maxFragLen)
{
    //String of intervals to return
    String<Triple<CharString, int, int > > intervalString;
    string chrStr;
    Triple<CharString, int, int > chr_start_end;
    ifstream intFile(toCString(intervalFile));
    if(intFile.fail())
    {
        cout << "Unable to locate interval file at: " << intervalFile << endl;
        return intervalString;
    }
    //Read first interval and add to string
    intFile >> chrStr;
    chr_start_end.i1 = chrStr;
    intFile >> chr_start_end.i2;
    --chr_start_end.i2;
    intFile >> chr_start_end.i3;
    --chr_start_end.i3;
    //cout << "Added: " << chr_start_end.i1 << ":" << chr_start_end.i2 << "-" << chr_start_end.i3 << " to interval string.\n";
    append(intervalString, chr_start_end);
    //If file only contains one interval, return.
    if (intFile.eof())
        return intervalString;
    while (!intFile.eof())
    {
        intFile >> chrStr;
        chr_start_end.i1 = chrStr;
        intFile >> chr_start_end.i2;
        --chr_start_end.i2;
        intFile >> chr_start_end.i3;
        --chr_start_end.i3;
        if (intFile.eof())
            break;
        //If beginning of interval is closer than 2*maxFragLen bases to the previous interval we merge them. Otherwise we cannot ensure sorting of reads.
        if (chr_start_end.i2-intervalString[length(intervalString)-1].i3 <= 2*maxFragLen && chr_start_end.i1 == intervalString[length(intervalString)-1].i1)
        {
            //cout << "Merging: " << intervalString[length(intervalString)-1].i2 << "-" << intervalString[length(intervalString)-1].i3 << " and " << chr_start_end.i2 << "-" << chr_start_end.i3 << " on " << chr_start_end.i1 << "\n";
            intervalString[length(intervalString)-1].i3 = chr_start_end.i3;
        }
        else
        {
            //cout << "Added: " << chr_start_end.i1 << ":" << chr_start_end.i2 << "-" << chr_start_end.i3 << " to interval string.\n";
            append(intervalString, chr_start_end);
        }
    }
    intFile.close();
    return intervalString;
}

int main(int argc, char const ** argv)
{
    if (argc != 7 && argc != 9)
    {
        cerr << "USAGE: " << argv[0] << " IN.bam OUT.bam maxFragmentLength keepMapQuality(Y/N) minNumMatches avgCovByReadLen.sh [baiFile intervalFile]\n";
        return 1;
    }
    cout<< "File to filter: " << argv[1] << endl;
    double avgCovByReadLen = lexicalCast<double>(argv[6]);
    CharString bamPathIn = argv[1], baiPathIn, intervalFile;
    int maxFragLen = lexicalCast<unsigned>(argv[3]), minMatchingBases = lexicalCast<unsigned>(argv[5]);
    bool keepMapQual = false, readBamSlice = false;
    string keepMapQualStr = argv[4];
    String<Triple<CharString, int, int > > intervalString;
    if (keepMapQualStr.compare("Y")==0)
        keepMapQual = true;
    if (argc == 9)
    {
        baiPathIn = argv[7];
        intervalFile = argv[8];
        intervalString = readIntervals(intervalFile, maxFragLen);
        // cout << "Listing intervals: " << endl;
        // for (unsigned i=0; i<length(intervalString); ++i)
        //     cout << intervalString[i].i2 << "-" << intervalString[i].i3 << endl;
        if (length(intervalString)==0)
        {
            std::cerr << "The interval file contained no intervals!" << endl;
            return 1;
        }
        readBamSlice = true;
    }
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamPathIn)))
    {
        std::cerr << "ERROR: Could not open " << bamPathIn << std::endl;
        return 1;
    }
    map<CharString, Pair<MateEditInfo> > mateEditMap;
    map<unsigned, String<BamAlignmentRecord> > beginPosToReads;
    map<CharString, BamAlignmentRecord> adapterMap;
    map<CharString, unsigned> readNameToNum;
    unsigned currIdx = 0;
    BamFileOut bamFileOut(context(bamFileIn), argv[2]);
    BamAlignmentRecord record;
    try
    {
        BamHeader header;
        try
        {
            readHeader(header, bamFileIn);
        } catch (...){
            std::cerr<<"Failed to read the header from the BAM file"<<endl;
            return 1;
        }
        writeHeader(bamFileOut, header);
        if (readBamSlice)
        {
            for (unsigned i=0; i<length(intervalString); ++i)
            {
                //cout << "Quality filtering interval: " << i << ", which is: " << intervalString[i].i1 << ":" << intervalString[i].i2 << "-" << intervalString[i].i3 << endl;
                int returnValue = qualityFilterSlice(intervalString[i], baiPathIn, bamFileIn, bamFileOut, mateEditMap, keepMapQual, maxFragLen, beginPosToReads, adapterMap, minMatchingBases, avgCovByReadLen, readNameToNum, currIdx);
                beginPosToReads.clear();
                if (returnValue != 0)
                {
                    std::cerr << "Something went wrong in filtering:" << intervalString[i].i1 << ":" << intervalString[i].i2 << "-" << intervalString[i].i3 << endl;
                    return 1;
                }
            }
        }
        else
        {
            double maxQueSum = avgCovByReadLen*(double)50.0*(double)3.0;
            unsigned currBeginPos = 0;
            deque<unsigned> myQue;
            BamAlignmentRecord record;
            CharString binary_qual;
            while (!atEnd(bamFileIn))
            {
                ++delStats.nTotalReads;
                readRecord(record, bamFileIn);
                removeHardClipped(record);
                if (record.beginPos != currBeginPos)
                {
                    if (record.beginPos-currBeginPos > 1)
                    {
                        for (unsigned i=1; i<record.beginPos-currBeginPos; ++i)
                        {
                            myQue.push_front(0);
                            if (myQue.size()>50)
                                myQue.pop_back();
                        }
                    }
                    myQue.push_front(1);
                    if (myQue.size()>50)
                        myQue.pop_back();
                    currBeginPos = record.beginPos;
                }
                else
                    ++myQue.front();
                if (std::accumulate(myQue.begin(),myQue.end(),0) > maxQueSum)
                {
                    ++delStats.nCoverageFiltered;
                    --myQue.front();
                    if (hasFlagRC(record))
                        mateEditMap[record.qName].i2.mateRemoved = true;
                    else
                        mateEditMap[record.qName].i1.mateRemoved = true;
                    continue;
                }
                if (qualityFilter(record, mateEditMap, maxFragLen, adapterMap, minMatchingBases, keepMapQual, beginPosToReads))
                {
                    binary_qual = "";
                    for (unsigned x=0; x<length(record.qual); ++x)
                    {
                        if (record.qual[x]-33 >= 25)
                            append(binary_qual,"I");
                        else
                            append(binary_qual,"!");
                    }
                    record.qual = binary_qual;
                    appendValue(beginPosToReads[record.beginPos], record);
                    if (record.beginPos-maxFragLen >=0)
                        printReadyReads(mateEditMap, beginPosToReads, record.beginPos-maxFragLen, bamFileOut, readNameToNum, currIdx, keepMapQual, Pair<unsigned>(0,260000000));
                }
                else
                    --myQue.front();
            }
            if (beginPosToReads.size()>0)
            {
                cout << "Smallest position in map before last call to printReadyReads: " << beginPosToReads.begin()->first << endl;
                int safePos = (beginPosToReads.rbegin()->first)+1;
                printReadyReads(mateEditMap, beginPosToReads, safePos, bamFileOut, readNameToNum, currIdx, keepMapQual, Pair<unsigned>(0,260000000));
            }
            if (beginPosToReads.size()>0)
                cout << beginPosToReads.size() << " positions were not printed, something is wrong. " << endl;
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    close(bamFileOut);
    close(bamFileIn);
    cout << "Soft clipped bp: " << delStats.nSoftClippedBp << " Number of coverage filtered reads: "<< delStats.nCoverageFiltered << " Quality clipped bp: " << delStats.nQualityClippedBp << " Not enough matches reads: " << delStats.nMatchRemovedReads << " Adapter removed bp: " << delStats.nAdapterClippedBp << " Number of adapter trimmed reads: " << delStats.nAdapterReads << " Total number of reads: " << delStats.nTotalReads << " Fragment of adapter reads: " << (double)delStats.nAdapterReads/(double)delStats.nTotalReads << endl;
    return 0;
}
