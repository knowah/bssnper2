# bssnper2

## Background

bssnper2 is an updated version of [BS-Snper](https://github.com/hellbelly/BS-Snper).

## Changes from BS-Snper

* Runtime efficiency updates:
    + `Factorial` (`BS-Snper.pl`) has been optimized (`genotyping.c:logFactorial_quotient()`).
    + Base-transition probability constants from `Bayes` (`BS-Snper.pl`) are now hard coded (`genotyping.h`).
* Memory efficiency updates:
    + Reference sequence FASTA file now read in 10 kbp chunks (`RefCache`) instead of the entire chromosome at once.
    + Chromosome-length arrays for storing basecall counts and quality scores have been abandoned in favor of a novel data structure (`gt_buffer`) with a default size of 1000 bp.
    + The memory usage for human chr1 (assuming 32-bit integers) is reduced from the order of 10+ GB to less than 1 MB.
* CIGAR parsing bugs corrected (see related [BS-Snper issue](https://github.com/hellbelly/BS-Snper/issues/10)).
* Overlapping bases of paired-end reads are now only counted once using the strategy from htslib (`bssnper2.c:tweak_overlap_quality()`)
* Support for outputting homozygous reference genotypes has been added (`--homf`, `--homrefInVCF`, `--assumeHomref`).
* Greater reliance on htslib library functions.
* Methylation calling functionality has been removed.

## Installation

bssnper2 requires an installation of htslib.
If the latest version does not work, the version used for bssnper2 development is commit `246c146f3f46d184b1dc3877ca35b16d13ee220a` of samtools/htslib, which can be downloaded [directly from github](https://github.com/samtools/htslib/archive/246c146f3f46d184b1dc3877ca35b16d13ee220a.zip).

After verifying htslib is installed, replace `HTSLIB` in `build.sh` with its directory and build bssnper2 using that script.