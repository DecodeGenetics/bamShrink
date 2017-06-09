# bamShrink

BamShrink extracts sequence reads of a region and reduces the output file size by binarizing base qualities values, removing unused BAM tags, removing unaligned read pairs, duplicate reads and reads that have less than a user specified number of matching bases in their alignment.
In addition to this, bamShrink performs coverage filtering in regions where the coverage is more than 3 times the average coverage, removes Ns if present on either end of a read and removes hard clipped entries from CIGAR strings.
Last, bamShrink performs adapter removal by clipping overhanging ends of read pairs where the reverse read has been aligned in front of the forward read and their alignments overlap.

## Installation
dependencies: SeqAn (zip-file provided with source code)
Compiler: Version 4.8.2 of g++

```sh
make bamShrink
```

## Usage
```sh
bamShrink IN.bam OUT.bam maxFramgentLength keepMapQuality(Y/N) minNumMatches avgCovByReadLen.sh [baiFile intervalFile]
```

## Things that bamShrink does
1. Fetches reads in a region or list of regions provided by user and their mates if they are within a user specified distance from each end of the region.
2. Unpairs reads that at further apart than a user specified distance or have the same orientation.
3. Filters read where both reads in the pair are unaligned, duplicates or with less than a user specified number of matching bp in BWA alignment.
4. Coverage filters in areas where coverage > 3*averageCoverage
5. Binarizes the quality string (qual >= 25 --> qual = max , qual < 25 --> qual = min)
6. Removes hard clipped entries from CIGAR strings
7. Removes N(s) at ends of reads where they are present
