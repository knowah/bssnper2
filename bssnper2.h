#ifndef BSSNPER2_H
#define BSSNPER2_H

#include <htslib/sam.h>
#include "genotyping.h"

#define BAMF_REV_R2 144
#define mate_same_chrom(b) (bool)((b)->core.mtid>=0 && (b)->core.tid==(b)->core.mtid)
#define await_mate(b) (bool)((b)->core.mpos>(b)->core.pos && (b)->core.mpos<bam_endpos(b));

enum read_calls { W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G };

void genotype_bam(struct GenotypingOptions *opts);

#endif

