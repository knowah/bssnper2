#include "genotyping.h"

// calculates log of [(a+t+c+g)! / (a!*t!*c!*g!)]
double logFactorial_quotient(uint32_t a, uint32_t t, uint32_t c, uint32_t g)
{
    uint32_t total = a + t + c + g;
    double logfac_bases = 0;
    double logfac_total = 0;

    for (uint32_t i = 2; i <= total; i++)
    {
        logfac_total += log((double)i);
        if (i == a) logfac_bases += logfac_total;
        if (i == t) logfac_bases += logfac_total;
        if (i == c) logfac_bases += logfac_total;
        if (i == g) logfac_bases += logfac_total;
    }
    return logfac_total - logfac_bases;
}

// guess genotype
void bayesian_genotype_inference(
    enum Genotype* Gt, int* Qual, char ref, 
    uint32_t w_A, uint32_t w_T, uint32_t w_C, uint32_t w_G, uint32_t c_A, uint32_t c_T, uint32_t c_C, uint32_t c_G,
    uint32_t w_Aq, uint32_t w_Tq, uint32_t w_Cq, uint32_t w_Gq, uint32_t c_Aq, uint32_t c_Tq, uint32_t c_Cq, uint32_t c_Gq
)
{
    int i;

    // convert quality scores
    double baseq_A = pow(0.1, (w_Aq>c_Aq ? w_Aq : c_Aq)/10);
    double baseq_T = pow(0.1, (w_Tq>c_Tq ? w_Tq : c_Tq)/10);
    double baseq_C = pow(0.1, (w_Cq>c_Cq ? w_Cq : c_Cq)/10);
    double baseq_G = pow(0.1, (w_Gq>c_Gq ? w_Gq : c_Gq)/10);

    // get total base counts (adjusted for A/T)
    uint32_t totals_A = ref=='A' ? (uint32_t)ceil(w_A+0.8*c_A) : w_A+c_A;
    uint32_t totals_T = ref=='T' ? (uint32_t)ceil(0.8*w_T+c_T) : w_T+c_T;
    uint32_t totals_C = w_C+c_C;
    uint32_t totals_G = w_G+c_G;
    
    double NN = logFactorial_quotient(totals_A, totals_T, totals_C, totals_G);
    double gt_probs[10] = {0};
 
    // some concerns with the following section:
    // why is the final term general baseq_X/3 ? should it not be /2 ?
    // why is it different for CG? and why is AT not 1-(bqA+bqT)/4 instead of 1/2-(bqA+bqT/4) ?
    if (baseq_A < 1) // some A reads with nonzero qual
    {
        gt_probs[gt_AA] = NN + totals_A*log(1-baseq_A) + (totals_T+totals_C+totals_G)*log(baseq_A/3);
        if (totals_T>0) gt_probs[gt_AT] = NN + (totals_A+totals_T)*log((1-(baseq_A+baseq_T)/2)/2) + (totals_C+totals_G)*log((baseq_A+baseq_T)/4);
        if (totals_C>0) gt_probs[gt_AC] = NN + (totals_A+totals_C)*log((1-baseq_A)/2) + (totals_T+totals_G)*log(baseq_A/3);
        if (totals_G>0) gt_probs[gt_AG] = NN + (totals_A+totals_G)*log((1-baseq_A)/2) + (totals_T+totals_C)*log(baseq_A/3);
    }
    if (baseq_T < 1) // some T reads with nonzero qual
    {
        gt_probs[gt_TT] = NN + totals_T*log(1-baseq_T) + (totals_A+totals_C+totals_G)*log(baseq_T/3);
        if (totals_C>0) gt_probs[gt_TC] = NN + (totals_T+totals_C)*log((1-baseq_T)/2) + (totals_A+totals_G)*log(baseq_T/3);
        if (totals_G>0) gt_probs[gt_TG] = NN + (totals_T+totals_G)*log((1-baseq_T)/2) + (totals_A+totals_C)*log(baseq_T/3);
    }
    if (baseq_C < 1) // some C reads with nonzero qual
    {
        gt_probs[gt_CC] = NN + totals_C*log(1-baseq_C) + (totals_A+totals_T+totals_G)*log(baseq_C/3);
        if (totals_G>0) gt_probs[gt_CG] = NN + (totals_C+totals_G)*log((1-baseq_C)/2) + (totals_A+totals_T)*log(baseq_C/2);
    }
    if (baseq_G < 1) // some G reads with nonzero qual
    {
        gt_probs[gt_GG] = NN + totals_G*log(1-baseq_G) + (totals_A+totals_T+totals_C)*log(baseq_G/3);
    }

    // identify which genotypes were supported
    int supported_gt[10] = {0};
    int support_count = 0;
    for (i = 0; i < 10; i++) if (gt_probs[i] != 0) { supported_gt[i] = 1; support_count += 1; }
    
    // if no useful basecalls, end here
    if (support_count == 0)
    {
        *Gt = gt_NN;
        *Qual = 0;
        return;
    }
    
    // add priors
    switch (ref)
    {
        case 'A':
            gt_probs[gt_AA] += REF_HOM;
            gt_probs[gt_AT] += REF_TVER;
            gt_probs[gt_AC] += REF_TVER;
            gt_probs[gt_AG] += REF_TSIT;
            gt_probs[gt_TT] += TVER_HOM;
            gt_probs[gt_TC] += TVER_TVER;
            gt_probs[gt_TG] += TSIT_TVER;
            gt_probs[gt_CC] += TVER_HOM;
            gt_probs[gt_CG] += TSIT_TVER;
            gt_probs[gt_GG] += TSIT_HOM;
            break;
        case 'T':
            gt_probs[gt_AA] += TVER_HOM;
            gt_probs[gt_AT] += REF_TVER;
            gt_probs[gt_AC] += TSIT_TVER;
            gt_probs[gt_AG] += TVER_TVER;
            gt_probs[gt_TT] += REF_HOM;
            gt_probs[gt_TC] += REF_TSIT;
            gt_probs[gt_TG] += REF_TVER;
            gt_probs[gt_CC] += TSIT_HOM;
            gt_probs[gt_CG] += TSIT_TVER;
            gt_probs[gt_GG] += TVER_HOM;
            break;
        case 'C':
            gt_probs[gt_AA] += TVER_HOM;
            gt_probs[gt_AT] += TSIT_TVER;
            gt_probs[gt_AC] += REF_TVER;
            gt_probs[gt_AG] += TVER_TVER;
            gt_probs[gt_TT] += TSIT_HOM;
            gt_probs[gt_TC] += REF_TSIT;
            gt_probs[gt_TG] += TSIT_TVER;
            gt_probs[gt_CC] += REF_HOM;
            gt_probs[gt_CG] += REF_TVER;
            gt_probs[gt_GG] += TVER_HOM;
            break;
        case 'G':
            gt_probs[gt_AA] += TSIT_HOM;
            gt_probs[gt_AT] += TSIT_TVER;
            gt_probs[gt_AC] += TSIT_TVER;
            gt_probs[gt_AG] += REF_TSIT;
            gt_probs[gt_TT] += TVER_HOM;
            gt_probs[gt_TC] += TVER_TVER;
            gt_probs[gt_TG] += REF_TVER;
            gt_probs[gt_CC] += TVER_HOM;
            gt_probs[gt_CG] += REF_TVER;
            gt_probs[gt_GG] += REF_HOM;
            break;
    }
    
    // iterate through supported genotypes, identifying
    // the genotype with the highest probability while
    // calculating the sum of all probabilities
    double fenmu = 0;
    double this_val, max_val;
    enum Genotype best_gt;
    bool seen = false;
    for (i = 0; i < 10; i++) {
        if (supported_gt[i]) {
            this_val = exp(gt_probs[i]);
            fenmu += this_val;

            if (!seen || this_val > max_val) {
                max_val = this_val;
                best_gt = i;
                seen = true;
            }
        }
    }
    
    // determine genotype quality 
    double gt_prob = 0;
    double qual;
    if (support_count > 1) {   
        if (fenmu == 0) qual = 0;
        else {
            gt_prob = 1-max_val/fenmu;
            qual = gt_prob==0 ? 1000 : -10*log10(gt_prob);
        }
    }
    else if (support_count == 1) {
        switch (best_gt) {
            case gt_AA: gt_prob = 1-1/(1+pow(0.5,w_A+c_A)); break;
            case gt_TT: gt_prob = 1-1/(1+pow(0.5,w_T+c_T)); break;
            case gt_CC: gt_prob = 1-1/(1+pow(0.5,w_C+c_C+w_T)); break;
            case gt_GG: gt_prob = 1-1/(1+pow(0.5,w_G+c_G+c_A)); break;
            default: gt_prob = 1; //any het genotype
            // all hets get prob=1 -> qual = 0 ?!?!
        }
        qual = gt_prob==0 ? 1000 : -10*log10(gt_prob);
    }

    // adjust genotype guesses

    // account for unevenly-stranded libraries
    // e.g. target capture
    // basically, if only evidence for G>A is from w_A, then assume no change
    // same goes for C>T and no c_T
    // maybe this should be what the 'assume_homref' does (with a new name,
    // e.g. assume_novar)
    if (ref=='G' && w_A==0) {
        if (best_gt == gt_AT) best_gt = gt_TG;
        if (best_gt == gt_AC) best_gt = gt_CG;
        if (best_gt == gt_AA || best_gt == gt_AG)
            best_gt = gt_GG;
    }
    else if (ref=='C' && c_T==0) {
        if (best_gt == gt_AT) best_gt = gt_AC;
        if (best_gt == gt_TG) best_gt = gt_CG;
        if (best_gt == gt_TT || best_gt == gt_TC)
            best_gt = gt_CC;
    }

    // adjustments from BS-Snper.pl
    if (best_gt == gt_TC && c_T+c_C > 0) {
        double t_ratio = (double)c_T/(c_T+c_C);
        if (t_ratio <= 0.05) best_gt = gt_CC;
        else if (t_ratio < 0.15 || t_ratio > 0.85) best_gt = gt_NN;
    }
    else if (best_gt == gt_AG && w_A+w_G > 0) {
        double a_ratio = (double)w_A/(w_A+w_G);
        if (a_ratio <= 0.05) best_gt = gt_GG;
        else if (a_ratio < 0.15 || a_ratio > 0.85) best_gt = gt_NN;
    }
    else if (best_gt == gt_TT && c_C > 1) {
        if ((double)c_C/(c_T+c_C) > 0.05) best_gt = gt_NN;
    }
    else if (best_gt == gt_AA && w_G > 1) {
        if ((double)w_G/(w_A+w_G) > 0.05) best_gt = gt_NN;
    }
    
    // assign genotype and quality score to passed variables
    *Gt = best_gt;   
    *Qual = (int)qual;

    return;
}

#define ddiv(a,b) ((double)(a)/(b))
void genotype(
    FILE* vcf_fptr, FILE* homref_fptr,
    const char* chrom, hts_pos_t pos, char ref,
    struct GenotypingOptions *opts,
    //uint32_t min_cover, uint32_t min_qual, uint32_t min_cover_alt, double min_hom_freq, double min_het_freq,
    uint32_t w_A, uint32_t w_T, uint32_t w_C, uint32_t w_G, uint32_t c_A, uint32_t c_T, uint32_t c_C, uint32_t c_G,
    uint32_t w_Aq, uint32_t w_Tq, uint32_t w_Cq, uint32_t w_Gq, uint32_t c_Aq, uint32_t c_Tq, uint32_t c_Cq, uint32_t c_Gq
)
{
    enum Genotype gt;
    int GT_qual;
    
    uint32_t total = w_A+w_T+w_C+w_G+c_A+c_T+c_C+c_G;

    // filter by read depth
    if (total < opts->min_depth || total > opts->max_depth) return;
    // TODO: implement option for saving 'uncalled' positions with >0 reads in a separate file?
    // TODO: allow different read depths for calling homref vs variants?
    // e.g. require 10x to report variant but only 5x to call homref?

    bayesian_genotype_inference(&gt, &GT_qual, ref, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq);
    //fprintf(stderr, "{%c %d}", ref, gt);
    if (gt == gt_NN) return;
    
    // sometimes the amount to compare to min_alt_cov
    // isn't equal to what's in ad[] ? according to perl script...
    // check this out 
    char *alts, *gt_str;
    double min_freq, variant_freq, alfr[3];
    uint32_t N_GT, qval, alt_cov, depth, adF[3], adR[3], ad[3];
    bool passed, homref = false;
    switch (ref) {
        case 'A': switch (gt) {
            case gt_AA:  // A>AA
                homref = true; if (!opts->homref_in_vcf) break;
                N_GT = 1; alts = "."; gt_str = "0/0";
                qval = w_Aq>c_Aq ? w_Aq : c_Aq;
                adF[0]=w_A; adR[0]=c_A; ad[0]=adF[0]+adR[0];
                depth = ad[0];
                alfr[0]=ddiv(ad[0],total);
                break;
            case gt_AT:  // A>AT
                N_GT = 2; alts = "T"; gt_str = "0/1";
                qval = w_Tq>c_Tq ? w_Tq : c_Tq;
                adF[0]=w_A; adR[0]=c_A; ad[0]=adF[0]+adR[0];
                adF[1]=w_T; adR[1]=c_T; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],total);
                break;
            case gt_AC:  // A>AC
                N_GT = 2; alts = "C"; gt_str = "0/1";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                adF[0]=w_A;     adR[0]=c_A; ad[0]=adF[0]+adR[0];
                adF[1]=w_C+w_T; adR[1]=c_C; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1]-w_T;
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],total);
                break;
            case gt_AG:  // A>AG
                N_GT = 2; alts = "G"; gt_str = "0/1";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                adF[0]=w_A; adR[0]=0; ad[0]=adF[0];
                adF[1]=w_G; adR[1]=0; ad[1]=adF[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_TT:  // A>TT
                N_GT = 2; alts = "T"; gt_str = "1/1";
                qval = c_Tq;
                adF[0]=w_A; adR[0]=c_A; ad[0]=adF[0]+adR[0];
                adF[1]=w_T; adR[1]=c_T; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_TC:  // A>TC 
                N_GT = 3; alts = "C,T"; gt_str = "1/2";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                qval = qval>c_Tq ? c_Tq : qval; // qval = min(c_Tq, max(w_Cq,c_Cq))
                adF[0]=0; adR[0]=c_A; ad[0]=adR[0];
                adF[1]=0; adR[1]=c_C; ad[1]=adR[1];
                adF[2]=0; adR[2]=c_T; ad[2]=adR[2];
                depth = c_C+c_T+w_A+c_A;
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_TG:  // A>TG
                N_GT = 3; alts = "G,T"; gt_str = "1/2";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                qval = qval>c_Tq ? c_Tq : qval; // qval = min(c_Tq, max(w_Cq,c_Cq))
                adF[0]=w_A; adR[0]=0;       ad[0]=adF[0];
                adF[1]=w_G; adR[1]=c_G+c_A; ad[1]=adF[1]+adR[1];
                adF[2]=w_T; adR[2]=c_T;     ad[2]=adF[2]+adR[2];
                uint32_t varG = w_G+c_G+w_A; // why?? shouldn't the w_A be c_A?
                depth = varG+ad[2];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(varG,depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_CC:  // A>CC
                N_GT = 2; alts = "C"; gt_str = "1/1";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                adF[0]=w_A;     adR[0]=c_A; ad[0]=adF[0]+adR[0];
                adF[1]=w_C+w_T; adR[1]=c_C; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_CG:  // A>CG
                N_GT = 3; alts = "C,G"; gt_str = "1/2";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                uint32_t qvalG = w_Gq>c_Gq ? w_Gq : c_Gq;
                qval = qval<qvalG ? qval : qvalG; // qval = min(max(w_Cq,c_Cq),max(w_Gq,c_Gq))
                adF[0]=w_A;     adR[0]=0;       ad[0]=adF[0];
                adF[1]=w_C+w_T; adR[1]=c_C;     ad[1]=adF[1]+adR[1];
                adF[2]=w_G;     adR[2]=c_G+c_A; ad[2]=adF[2]+adR[2];
                depth = ad[1]+ad[2];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_GG:  // A>GG
                N_GT = 2; alts = "G"; gt_str = "1/1";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                adF[0]=w_A; adR[0]=0;       ad[0]=adF[0];
                adF[1]=w_G; adR[1]=c_G+c_A; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            }; break;
        case 'T': switch (gt) {
            case gt_AA:  // T>AA
                N_GT = 2; alts = "A"; gt_str = "1/1";
                qval = w_Aq>c_Aq ? w_Aq : c_Aq;
                adF[0]=w_T; adR[0]=c_T; ad[0]=adF[0]+adR[0];
                adF[1]=w_A; adR[1]=c_A; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],total);
                break;
            case gt_AT:  // T>AG  
                N_GT = 2; alts = "A"; gt_str = "0/1";
                qval = w_Aq>c_Aq ? w_Aq : c_Aq;
                adF[0]=w_T; adR[0]=c_T; ad[0]=adF[0]+adR[0];
                adF[1]=w_A; adR[1]=c_A; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],total);
                break;
            case gt_AC:  // T>AC
                N_GT = 3; alts = "A,C"; gt_str = "1/2";
                qval = w_Aq>c_Aq ? w_Aq : c_Aq;
                uint32_t qvalC = w_Cq>c_Cq ? w_Cq : c_Cq;
                qval = qval<qvalC ? qval : qvalC; // qval = min(max(w_Aq,c_Aq),max(w_Cq,c_Cq))
                adF[0]=w_T;     adR[0]=c_T; ad[0]=adR[0]; // don't count w_T in AD
                adF[1]=w_A;     adR[1]=c_A; ad[1]=adF[1]+adR[1];
                adF[2]=w_C+w_T; adR[2]=c_C; ad[2]=adF[2]+adR[2];
                depth = ad[1]+ad[2]-w_T;
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_AG:  // T>AG
                N_GT = 3; alts = "A,G"; gt_str = "1/2";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                qval = qval<w_Aq ? qval : w_Aq; // qval = min(w_Aq,max(w_G,c_Gq))
                adF[0]=w_T; adR[0]=0; ad[0]=adF[0]+c_T; // not sure why c_T added here...
                adF[1]=w_A; adR[1]=0; ad[1]=adF[1];
                adF[2]=w_G; adR[2]=0; ad[2]=adF[2];
                depth = ad[1]+ad[2]; // T calls excluded from depth, but depth used as denominator for T ALFR, also strange...
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_TT:  // T>TT
                homref = true; if (!opts->homref_in_vcf) break;
                N_GT = 1; alts = "."; gt_str = "0/0";
                qval = w_Tq>c_Tq ? w_Tq : c_Tq;
                adF[0]=w_T; adR[0]=c_T; ad[0]=adF[0]+adR[0];
                depth = ad[0];
                alfr[0]=ddiv(ad[0],total);
                break;
            case gt_TC:  // T>TC
                N_GT = 2; alts = "C"; gt_str = "0/1";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                adF[0]=0; adR[0]=c_T; ad[0]=adR[0];
                adF[1]=0; adR[1]=c_C; ad[1]=adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_TG:  // T>TG
                N_GT = 2; alts = "G"; gt_str = "0/1";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                adF[0]=w_T; adR[0]=c_T;     ad[0]=adF[0]+adR[0];
                adF[1]=w_G; adR[1]=c_G+c_A; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],depth); // why one total and one depth?
                break;
            case gt_CC:  // T>CC
                N_GT = 2; alts = "C"; gt_str = "1/1";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                adF[0]=0;       adR[0]=c_T; ad[0]=adR[0];
                adF[1]=w_C+w_T; adR[1]=c_C; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_CG:  // T>CG
                N_GT = 3; alts = "C,G"; gt_str = "1/2";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                uint32_t qvalG = w_Gq>c_Gq ? w_Gq : c_Gq;
                qval = qval<qvalG ? qval : qvalG; // qval = min(max(w_Cq,c_Cq),max(w_Gq,c_Gq))
                adF[0]=0;       adR[0]=c_T;     ad[0]=adR[0];
                adF[1]=w_C+w_T; adR[1]=c_C;     ad[1]=adF[1]+adR[1];
                adF[2]=w_G;     adR[2]=c_G+c_A; ad[2]=adF[2]+adR[2];
                depth = ad[1]+ad[2];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_GG:  // T>GG
                N_GT = 2; alts = "G"; gt_str = "1/1";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                adF[0]=w_T; adR[0]=c_T;     ad[0]=adF[0]+adR[0];
                adF[1]=w_G; adR[1]=c_G+c_A; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            }; break;
        case 'C': switch (gt) {
            case gt_AA:  // C>AA
                N_GT = 2; alts = "A"; gt_str = "1/1";
                qval = w_Aq>c_Aq ? w_Aq : c_Aq;
                adF[0]=w_C; adR[0]=c_C; ad[0]=adF[0]+adR[0];
                adF[1]=w_A; adR[1]=c_A; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],total);
                break;
            case gt_AT:  // C>AG
                N_GT = 3; alts = "A,T"; gt_str = "1/2";
                qval = w_Aq>c_Aq ? w_Aq : c_Aq;
                qval = qval>c_Tq ? c_Tq : qval; // qval = min(c_Tq, max(w_Aq, c_Aq))
                adF[0]=w_C; adR[0]=c_C; ad[0]=adF[0]+adR[0];
                adF[1]=w_A; adR[1]=c_A; ad[1]=adF[1]+adR[1];
                adF[2]=w_T; adR[2]=c_T; ad[2]=adF[2]+adR[2];
                depth = ad[1]+ad[2];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],total); alfr[2]=ddiv(ad[2],total);
                break;
            case gt_AC:  // C>AC
                N_GT = 2; alts = "A"; gt_str = "0/1";
                qval = w_Aq>c_Aq ? w_Aq : c_Aq;
                adF[0]=w_C+w_T; adR[0]=c_C; ad[0]=adF[0]+adR[0]; // messed up in perl version...
                adF[1]=w_A;     adR[1]=c_A; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1]-w_T;
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],total);
                break;
            case gt_AG:  // C>AG
                N_GT = 3; alts = "A,G"; gt_str = "1/2";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                qval = qval<w_Aq ? qval : w_Aq; // qval = min(w_Aq,max(w_G,c_Gq))
                adF[0]=w_C; adR[0]=0; ad[0]=adF[0]; // perl version gives ad[0]=w_T...
                adF[1]=w_A; adR[1]=0; ad[1]=adF[1];
                adF[2]=w_G; adR[2]=0; ad[2]=adF[2];
                depth = ad[1]+ad[2];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_TT:  // C>TT
                /*
                if (c_T==0) {
                    // no c_T reads -> probably still TT genotype
                    if (opts->assume_homref || (w_C+c_C)>0) {
                        homref=1; break;
                    } else return;
                }
                */
                N_GT = 2; alts = "T"; gt_str = "1/1";
                qval = c_Tq;
                adF[0]=w_C; adR[0]=c_C; ad[0]=adF[0]+adR[0];
                adF[1]=w_T; adR[1]=c_T; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_TC:  // C>TC
                /*
                if (c_T==0) {
                    // no c_T reads -> probably still TT genotype
                    if (opts->assume_homref || (w_C+c_C)>0) {
                        homref=1; break;
                    } else return;
                }
                */
                N_GT = 2; alts = "T"; gt_str = "0/1";
                qval = c_Tq;
                adF[0]=0; adR[0]=c_C; ad[0]=adR[0];
                adF[1]=0; adR[1]=c_T; ad[1]=adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_TG:  // C>TG
                N_GT = 3; alts = "G,T"; gt_str = "1/2";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                qval = qval<c_Tq ? qval : c_Tq; // qval = min(c_Tq,max(w_Gq,c_Gq))
                adF[0]=w_C; adR[0]=c_C;     ad[0]=adF[0]+adR[0];
                adF[1]=w_G; adR[1]=c_G+c_A; ad[1]=adF[1]+adR[1];
                adF[2]=w_T; adR[2]=c_T;     ad[2]=adF[2]+adR[2];
                depth = ad[1]+ad[2];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_CC:  // C>CC
                homref = true; if (!opts->homref_in_vcf) break;
                N_GT = 1; alts = "."; gt_str = "0/0";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                adF[0]=w_C+w_T; adR[0]=c_C; ad[0]=adF[0]+adR[0];
                depth = ad[0];
                alfr[0]=ddiv(ad[0],total);
                break;
            case gt_CG:  // C>CG
                N_GT = 2; alts = "G"; gt_str = "0/1";
                qval = w_Gq;
                adF[0]=w_C+w_T; adR[0]=c_C;     ad[0]=adF[0]+adR[0];
                adF[1]=w_G;     adR[1]=c_G+c_A; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_GG:  // C>GG
                N_GT = 2; alts = "G"; gt_str = "1/1";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                adF[0]=w_C; adR[0]=c_C;     ad[0]=adF[0]+adR[0];
                adF[1]=w_G; adR[1]=c_G+c_A; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1]+w_T; // perl gives c_T instead of c_A...
                alfr[0]=ddiv(ad[0],ad[0]+ad[1]); alfr[1]=ddiv(ad[1],depth);
                break;
            }; break;
        case 'G': switch (gt) {
            case gt_AA:  // G>AA
                /*
                if (w_A==0) {
                    // no w_A reads -> probably still GG genotype
                    if (opts->assume_homref || (w_G+c_G)>0) {
                        homref=1; break;
                    } else return;
                }
                */
                N_GT = 2; alts = "A"; gt_str = "1/1";
                qval = w_Aq;
                adF[0]=w_G; adR[0]=c_G; ad[0]=adF[0]+adR[0];
                adF[1]=w_A; adR[1]=c_A; ad[1]=adF[1]+adR[1];
                
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth); // C>AA uses 'total' as denom. why different here?
                break;
            case gt_AT:  // G>AT
                N_GT = 3; alts = "A,T"; gt_str = "1/2";
                qval = w_Tq>c_Tq ? w_Tq : c_Tq;
                qval = qval>w_Aq ? w_Aq : qval; // qval = min(w_Aq, max(w_Tq, c_Tq))
                adF[0]=w_G; adR[0]=c_G; ad[0]=adF[0]+adR[0];
                adF[1]=w_A; adR[1]=c_A; ad[1]=adF[1]+adR[1];
                adF[2]=w_T; adR[2]=c_T; ad[2]=adF[2]+adR[2];
                depth = ad[1]+ad[2];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],total); alfr[2]=ddiv(ad[2],total);
                break;
            case gt_AC:  // G>AC
                N_GT = 3; alts = "A,C"; gt_str = "1/2";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                qval = qval<w_Aq ? qval : w_Aq; // qval = min(w_Aq,max(w_Cq,c_Cq))
                adF[0]=w_G;     adR[0]=c_G; ad[0]=adF[0]+adR[0];
                adF[1]=w_A;     adR[1]=c_A; ad[1]=adF[1]+adR[1];
                adF[2]=w_C+w_T; adR[2]=c_C; ad[2]=adF[2]+adR[2];
                depth = ad[1]+ad[2]-w_T;
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_AG:  // G>AG
                /*
                if (w_A==0) {
                    // no w_A reads -> probably still GG genotype
                    if (opts->assume_homref || (w_G+c_G)>0) {
                        homref=1; break;
                    } else return;
                }
                */
                N_GT = 2; alts = "A"; gt_str = "0/1";
                qval = w_Aq;
                adF[0]=w_G; adR[0]=0; ad[0]=adF[0];
                adF[1]=w_A; adR[1]=0; ad[1]=adF[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_TT:  // G>TT
                N_GT = 2; alts = "T"; gt_str = "1/1";
                qval = c_Tq;
                adF[0]=w_G; adR[0]=c_G; ad[0]=adF[0]+adR[0];
                adF[1]=w_T; adR[1]=c_T; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],total); // why depth for one and total for the other?
                break;
            case gt_TC:  // G>TC
                N_GT = 3; alts = "C,T"; gt_str = "1/2";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                qval = qval>c_Tq ? c_Tq : qval; // qval = min(c_Tq, max(w_Cq,c_Cq))
                adF[0]=0; adR[0]=c_G; ad[0]=adR[0];
                adF[1]=0; adR[1]=c_C; ad[1]=adR[1];
                adF[2]=0; adR[2]=c_T; ad[2]=adR[2];
                depth = c_C+c_T+w_G+c_G;
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth); alfr[2]=ddiv(ad[2],depth);
                break;
            case gt_TG:  // G>TG
                N_GT = 2; alts = "T"; gt_str = "0/1";
                qval = c_Tq;
                adF[0]=w_G; adR[0]=c_G+c_A; ad[0]=adF[0]+adR[0];
                adF[1]=w_T; adR[1]=c_T;     ad[1]=adF[1]+adF[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],total); alfr[1]=ddiv(ad[1],depth); // why total for one and depth for the other?
                break;
            case gt_CC:  // G>CC
                N_GT = 2; alts = "C"; gt_str = "1/1";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                adF[0]=w_G;     adR[0]=c_G; ad[0]=adF[0]+adR[0];
                adF[1]=w_C+w_T; adR[1]=c_C; ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_CG:  // G>CG
                N_GT = 2; alts = "C"; gt_str = "0/1";
                qval = w_Cq>c_Cq ? w_Cq : c_Cq;
                adF[0]=w_G;     adR[0]=c_G+c_A; ad[0]=adF[0]+adR[0];
                adF[1]=w_C+w_T; adR[1]=c_C;     ad[1]=adF[1]+adR[1];
                depth = ad[0]+ad[1];
                alfr[0]=ddiv(ad[0],depth); alfr[1]=ddiv(ad[1],depth);
                break;
            case gt_GG:  // G>GG
                homref = true; if (!opts->homref_in_vcf) break;
                N_GT = 1; alts = "."; gt_str = "0/0";
                qval = w_Gq>c_Gq ? w_Gq : c_Gq;
                adF[0]=w_G; adR[0]=c_G+c_A; ad[0]=adF[0]+adR[0];
                depth = ad[0];
                alfr[0]=ddiv(ad[0],total);
                break;
            }; break;
        default: return;
    }
    
        if (homref) {
        if (opts->homref_in_vcf) {
            passed = depth >= opts->min_depth && qval >= opts->min_base_qual && alfr[0] >= opts->min_hom_freq;
        } else {
            // write entry to homref file
            fprintf(homref_fptr,
                "%s\t%ld\t%c\t%u,%u,%u,%u,%u,%u,%u,%u\n",
                chrom, pos+1, ref, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G
            );
            fflush(homref_fptr);
            return;
        }
    } else {
        min_freq = (gt==gt_AA||gt==gt_TT||gt==gt_CC||gt==gt_GG) ? opts->min_hom_freq : opts->min_het_freq;
        alt_cov = N_GT<3 ? ad[1] : (ad[1]<ad[2] ? ad[1] : ad[2]);
        variant_freq = N_GT<3 ? alfr[1] : (alfr[1]<alfr[2] ? alfr[1] : alfr[2]);
        passed = depth >= opts->min_depth && qval >= opts->min_base_qual && alt_cov >= opts->min_alt_count && variant_freq >= min_freq;
    }
    
    if (depth == 0) return;
    // write VCF entry
    // initial mandatory fields
    fprintf(vcf_fptr, "%s\t%ld\t.\t%c\t%s\t%d\t%s\t", chrom, pos+1, ref, alts, GT_qual, passed ? "PASS" : "Low");
    //summary fields
    fprintf(vcf_fptr, "DP=%u;ADF=", total); write_uint_array(vcf_fptr,adF,N_GT,',');
    fprintf(vcf_fptr, ";ADR="); write_uint_array(vcf_fptr,adR,N_GT,',');
    fprintf(vcf_fptr, ";AD="); write_uint_array(vcf_fptr,ad,N_GT,',');
    //sample fields
    fprintf(vcf_fptr, "\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t%s:%u:", gt_str, depth); //GT:DP:
    write_uint_array(vcf_fptr,adF,N_GT,','); fprintf(vcf_fptr, ":"); //ADF:
    write_uint_array(vcf_fptr,adR,N_GT,','); fprintf(vcf_fptr, ":"); //ADR:
    write_uint_array(vcf_fptr,ad,N_GT,','); //AD
    fprintf(vcf_fptr,
        ":%u,%u,%u,%u,%u,%u,%u,%u:%u,%u,%u,%u,%u,%u,%u,%u:",
        w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G,
        w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq
    ); //:BSD:BSQ:
    write_dbl_array(vcf_fptr,alfr,N_GT,','); //ALFR
    fprintf(vcf_fptr, "\n");
    fflush(vcf_fptr);
}

void write_uint_array(FILE* fptr, uint32_t* arr, int len, char sep)
{
    for (int i = 0; i < len-1; i++)
		fprintf(fptr, "%u%c", arr[i], sep);
    fprintf(fptr, "%u", arr[len-1]);
    return;
}

void write_dbl_array(FILE* fptr, double* arr, int len, char sep)
{
    for (int i = 0; i < len-1; i++)
		fprintf(fptr, "%0.3f%c", arr[i], sep);
    fprintf(fptr, "%0.3f", arr[len-1]);
    return;
}
