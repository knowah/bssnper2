#ifndef BSSNPER2_GENOTYPING_H
#define BSSNPER2_GENOTYPING_H

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <htslib/hts.h> // for hts_pos_t typedef

/*
log constants of genotype probability priors
adapted from BS-Snper.pl
source: 
Li et al. (2009) SNP detection for massively parallel whole-genome resequencing.
  Genome Res. doi:10.1101/gr.088013.108
assumptions:
- hom SNP rate = 0.0005
- het SNP rate = 0.001
- transition:transversion ratio = 4:1
*/
static const double REF_HOM = -0.001501126126;    // log(0.9985)
static const double TSIT_HOM = -8.006367567650;   // log(0.0005*2/3)
static const double TVER_HOM = -9.392661928770;   // log(0.0005/6)
static const double REF_TSIT = -7.313220387090;   // log(0.001*2/3)
static const double REF_TVER = -8.699514748210;   // log(0.001/6)
static const double TSIT_TVER = -16.012735135300; // log(0.001*2/3)+log(0.001/6)
static const double TVER_TVER = -17.399029496420; // 2*log(0.001/6)

enum Genotype { gt_AA, gt_AT, gt_AC, gt_AG, gt_TT, gt_TC, gt_TG, gt_CC, gt_CG, gt_GG, gt_NN=-1 };

struct GenotypingOptions
{
    char *bam_fname;
    char *ref_fname;
    char *vcf_fname;
    char *homref_fname;
    uint8_t min_base_qual;
    uint32_t min_depth;
    uint32_t max_depth;
    double min_hom_freq;
    double min_het_freq;
    double error_rate;
    uint8_t min_mapq;
    uint32_t min_alt_count;
    uint16_t buffer_size;
    bool assume_homref;
    bool homref_in_vcf;
    char *sample_name;
    int cmd_argc;
    char **cmd_argv;
};

double logFactorial_quotient(uint32_t a, uint32_t t, uint32_t c, uint32_t g);
void bayesian_genotype_inference(
    enum Genotype* Gt, int* Qual, char ref, 
    uint32_t w_A, uint32_t w_T, uint32_t w_C, uint32_t w_G, uint32_t c_A, uint32_t c_T, uint32_t c_C, uint32_t c_G,
    uint32_t w_Aq, uint32_t w_Tq, uint32_t w_Cq, uint32_t w_Gq, uint32_t c_Aq, uint32_t c_Tq, uint32_t c_Cq, uint32_t c_Gq
);
void genotype(
    FILE* vcf_fptr, FILE* homref_fptr, const char* chrom, hts_pos_t pos, char ref,
    struct GenotypingOptions *opts,
    //uint32_t min_cover, uint32_t min_qual, uint32_t min_cover_alt, double min_hom_freq, double min_het_freq,
    uint32_t w_A, uint32_t w_T, uint32_t w_C, uint32_t w_G, uint32_t c_A, uint32_t c_T, uint32_t c_C, uint32_t c_G,
    uint32_t w_Aq, uint32_t w_Tq, uint32_t w_Cq, uint32_t w_Gq, uint32_t c_Aq, uint32_t c_Tq, uint32_t c_Cq, uint32_t c_Gq
);
void write_uint_array(FILE* fptr, uint32_t *arr, int len, char sep);
void write_dbl_array(FILE* fptr, double *arr, int len, char sep);

#endif
