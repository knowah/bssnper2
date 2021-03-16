#include "bssnper2.h"
#include "uncalled_pos.h"
#include "genotyping.h"
#include <ctype.h>
#include <stdio.h>
#include <htslib/faidx.h>
#include <time.h>

/**
 *  cigar_iref2iseq_set()  - find the first CMATCH setting the ref and the read index
 *  cigar_iref2iseq_next() - get the next CMATCH base
 *  @cigar:       pointer to current cigar block (rw)
 *  @cigar_max:   pointer just beyond the last cigar block
 *  @icig:        position within the current cigar block (rw)
 *  @iseq:        position in the sequence (rw)
 *  @iref:        position with respect to the beginning of the read (iref_pos - b->core.pos) (rw)
 *
 *  Returns BAM_CMATCH, -1 when there is no more cigar to process or the requested position is not covered,
 *  or -2 on error.
 */
static inline int cigar_iref2iseq_set(uint32_t **cigar, uint32_t *cigar_max, hts_pos_t *icig, hts_pos_t *iseq, hts_pos_t *iref)
{
    hts_pos_t pos = *iref;
    if ( pos < 0 ) return -1;
    *icig = 0;
    *iseq = 0;
    *iref = 0;
    while ( *cigar<cigar_max )
    {
        int cig  = (**cigar) & BAM_CIGAR_MASK;
        int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

        if ( cig==BAM_CSOFT_CLIP ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CHARD_CLIP || cig==BAM_CPAD ) { (*cigar)++; *icig = 0; continue; }
        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            pos -= ncig;
            if ( pos < 0 ) { *icig = ncig + pos; *iseq += *icig; *iref += *icig; return BAM_CMATCH; }
            (*cigar)++; *iseq += ncig; *icig = 0; *iref += ncig;
            continue;
        }
        if ( cig==BAM_CINS ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CDEL || cig==BAM_CREF_SKIP )
        {
            pos -= ncig;
            if ( pos<0 ) pos = 0;
            (*cigar)++; *icig = 0; *iref += ncig;
            continue;
        }
        hts_log_error("Unexpected cigar %d", cig);
        return -2;
    }
    *iseq = -1;
    return -1;
}
static inline int cigar_iref2iseq_next(uint32_t **cigar, uint32_t *cigar_max, hts_pos_t *icig, hts_pos_t *iseq, hts_pos_t *iref)
{
    while ( *cigar < cigar_max )
    {
        int cig  = (**cigar) & BAM_CIGAR_MASK;
        int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            //if ( *icig >= ncig - 1 ) { *icig = 0;  (*cigar)++; continue; }
            if ( *icig >= ncig - 1) { 
				*icig = 0; (*cigar)++;
				if (*cigar < cigar_max) continue;
			} else
				(*icig)++;
            (*iseq)++; (*iref)++;
            return BAM_CMATCH;
        }
        if ( cig==BAM_CDEL || cig==BAM_CREF_SKIP ) { (*cigar)++; (*iref) += ncig; *icig = 0; continue; }
        if ( cig==BAM_CINS ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CSOFT_CLIP ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CHARD_CLIP || cig==BAM_CPAD ) { (*cigar)++; *icig = 0; continue; }
        hts_log_error("Unexpected cigar %d", cig);
        return -2;
    }
    *iseq = -1;
    *iref = -1;
    return -1;
}

static int tweak_overlap_quality(bam1_t *a, bam1_t *b)
{
    uint32_t *a_cigar = bam_get_cigar(a), *a_cigar_max = a_cigar + a->core.n_cigar;
    uint32_t *b_cigar = bam_get_cigar(b), *b_cigar_max = b_cigar + b->core.n_cigar;
    hts_pos_t a_icig = 0, a_iseq = 0;
    hts_pos_t b_icig = 0, b_iseq = 0;
    uint8_t *a_qual = bam_get_qual(a), *b_qual = bam_get_qual(b);
    uint8_t *a_seq  = bam_get_seq(a), *b_seq = bam_get_seq(b);

    hts_pos_t iref   = b->core.pos;
    hts_pos_t a_iref = iref - a->core.pos;
    hts_pos_t b_iref = iref - b->core.pos;
    int a_ret = cigar_iref2iseq_set(&a_cigar, a_cigar_max, &a_icig, &a_iseq, &a_iref);
    if ( a_ret<0 ) return a_ret<-1 ? -1:0;  // no overlap or error
    int b_ret = cigar_iref2iseq_set(&b_cigar, b_cigar_max, &b_icig, &b_iseq, &b_iref);
    if ( b_ret<0 ) return b_ret<-1 ? -1:0;  // no overlap or error

    #if DBG
        fprintf(stderr,"tweak %s  n_cigar=%d %d  .. %d-%d vs %"PRIhts_pos"-%"PRIhts_pos"\n", bam_get_qname(a), a->core.n_cigar, b->core.n_cigar,
            a->core.pos+1,a->core.pos+bam_cigar2rlen(a->core.n_cigar,bam_get_cigar(a)), b->core.pos+1, b->core.pos+bam_cigar2rlen(b->core.n_cigar,bam_get_cigar(b)));
    #endif

    int err = 0;
    while ( 1 )
    {
        // Increment reference position
        while ( a_ret >= 0 && a_iref>=0 && a_iref < iref - a->core.pos )
            a_ret = cigar_iref2iseq_next(&a_cigar, a_cigar_max, &a_icig, &a_iseq, &a_iref);
        if ( a_ret<0 || a_iseq == a->core.l_qseq) { err = a_ret<-1?-1:0; break; }   // done
        if ( iref < a_iref + a->core.pos ) iref = a_iref + a->core.pos;

        while ( b_ret >= 0 && b_iref>=0 && b_iref < iref - b->core.pos )
            b_ret = cigar_iref2iseq_next(&b_cigar, b_cigar_max, &b_icig, &b_iseq, &b_iref);
        if ( b_ret<0 || b_iseq == b->core.l_qseq) { err = b_ret<-1?-1:0; break; }  // done
        if ( iref < b_iref + b->core.pos ) iref = b_iref + b->core.pos;

        iref++;
        if ( a_iref+a->core.pos != b_iref+b->core.pos ) continue;   // only CMATCH positions, don't know what to do with indels

        if (a_iseq > a->core.l_qseq || b_iseq > b->core.l_qseq)
            return -1;  // Fell off end of sequence, bad CIGAR?

        if ( bam_seqi(a_seq,a_iseq) == bam_seqi(b_seq,b_iseq) )
        {
            #if DBG
                fprintf(stderr,"%c",seq_nt16_str[bam_seqi(a_seq,a_iseq)]);
            #endif
            // we are very confident about this base
            int qual = a_qual[a_iseq] + b_qual[b_iseq];
            a_qual[a_iseq] = qual>200 ? 200 : qual;
            b_qual[b_iseq] = 0;
        }
        else
        {
            if ( a_qual[a_iseq] >= b_qual[b_iseq] )
            {
                #if DBG
                    fprintf(stderr,"[%c/%c]",seq_nt16_str[bam_seqi(a_seq,a_iseq)],tolower_c(seq_nt16_str[bam_seqi(b_seq,b_iseq)]));
                #endif
                a_qual[a_iseq] = 0.8 * a_qual[a_iseq];  // not so confident about a_qual anymore given the mismatch
                b_qual[b_iseq] = 0;
            }
            else
            {
                #if DBG
                    fprintf(stderr,"[%c/%c]",tolower_c(seq_nt16_str[bam_seqi(a_seq,a_iseq)]),seq_nt16_str[bam_seqi(b_seq,b_iseq)]);
                #endif
                b_qual[b_iseq] = 0.8 * b_qual[b_iseq];
                a_qual[a_iseq] = 0;
            }
        }
    }
    #if DBG
        fprintf(stderr,"\n");
    #endif
    return err;
}

#define udiv(n,d) ((uint32_t)((d)>0?((double)(n)/(d)+0.5):0))
void call_genotype(FILE *vcf_fptr, FILE *homref_fptr, const char* chrom, gt_pos* gp, char ref, struct GenotypingOptions *opts)
{
	/*
	fprintf(vcf_out, "%s\t%ld\t", chrom, gp->pos+1);
	uint32_t total = 0;
	for (int i = 0; i < 8; i++) {
		total += gp->call_counts[i];
		fprintf(vcf_out, "%u,%u;", gp->call_counts[i], gp->qual_sums[i]);
	}
	fprintf(vcf_out, "\t%u\n", total);
	return;
	*/
	uint32_t quals[8];
	uint32_t *calls = gp->call_counts;
	for (int i = 0; i < 8; i++)
		quals[i] = udiv(gp->qual_sums[i],calls[i]);

	genotype(
	  vcf_fptr, homref_fptr, chrom, gp->pos, ref, opts,
	  calls[W_A], calls[W_T], calls[W_C], calls[W_G], calls[C_A], calls[C_T], calls[C_C], calls[C_G],
	  quals[W_A], quals[W_T], quals[W_C], quals[W_G], quals[C_A], quals[C_T], quals[C_C], quals[C_G]
	);
}

void update_pos_buffer(gt_buffer *gb, hts_pos_t pos, uint8_t basecall, uint8_t qual, bool is_W)
{
	gt_pos *gp = retrieve_gt_pos(gb, pos, true);
	int index = -1;
	if (is_W) {
		switch (basecall) {
			case 'A': index = W_A; break;
			case 'T': index = W_T; break;
			case 'C': index = W_C; break;
			case 'G': index = W_G; break;
		}
	} else {
		switch (basecall) {
			case 'A': index = C_A; break;
			case 'T': index = C_T; break;
			case 'C': index = C_C; break;
			case 'G': index = C_G; break;
		}
	}
	if (index < W_A || index > C_G) return;
	gp->call_counts[index] += 1;
	gp->qual_sums[index] += qual;
}

void evaluate_read(gt_buffer *gb, bam1_t *b, hts_pos_t eval_start, hts_pos_t eval_end, uint8_t min_qual)
{
	const bam1_core_t *read = &b->core;
	hts_pos_t refpos = read->pos;
	if (refpos < 0) return;
	
	bool past_start = eval_start <= 0 ? true : (refpos >= eval_start);
	bool watson = (read->flag&BAMF_REV_R2)==BAMF_REV_R2 || (read->flag&BAMF_REV_R2)==0;

	//uint8_t *seq = bam_get_seq(b);
	uint8_t *quals = bam_get_qual(b);
	
	uint32_t op_len;
	uint16_t op_num;
	uint8_t cigar_op;
	uint32_t readpos = 0;
	for (op_num = 0; op_num < read->n_cigar; op_num++) {
		op_len = bam_get_cigar(b)[op_num]>>BAM_CIGAR_SHIFT;
		cigar_op = BAM_CIGAR_STR[bam_get_cigar(b)[op_num]&BAM_CIGAR_MASK];
		switch(cigar_op) {
			case 'M':
			case '=':
			case 'X':
				while (op_len-- > 0) { // check this is correct number of operations
					if (!past_start && refpos >= eval_start)
						past_start = true;
					if (eval_end > 0 && refpos >= eval_end)
						return;
					
					uint8_t this_qual = quals[readpos];
					if (past_start && this_qual >= min_qual)
						update_pos_buffer(gb, refpos, seq_nt16_str[bam_seqi(bam_get_seq(b), readpos)], this_qual, watson);
					readpos += 1;
					refpos += 1;
				}
				break;
			case 'D':
				// register_deletion(refpos+1, op_len);
				refpos += op_len;
				break;
			case 'I':
				// register_insertion(refpos+1, op_len, b, seq+readpos);
			case 'S':
				readpos += op_len;
				break;
			case 'H': break;
			default:
				fprintf(stderr, "ERROR: Encountered malformed/unsupported CIGAR operation in read %s\n", bam_get_qname(b));
				exit(1);
		}
	}
}

typedef struct 
{
	faidx_t *fai;
	char *seq;
	hts_pos_t size;
	hts_pos_t start, end;
	hts_pos_t ret;
} RefCache;

void reset_refcache(RefCache *r)
{
	if (r->seq != NULL) {
		free(r->seq);
		r->seq = NULL;
	}
	r->start = -1;
	r->end = -1;
	r->ret = 0;
}

char get_ref_base(const char *chrom, hts_pos_t pos, RefCache *ref)
{
	if (ref->seq != NULL) {
		if (pos >= ref->start && pos <= ref->end)
			return ref->seq[pos-(ref->start)];
		free(ref->seq);
		ref->seq = NULL;
	}
	
	ref->start = pos;
	ref->seq = faidx_fetch_seq64(ref->fai, chrom, ref->start, pos + ref->size - 1, &(ref->ret));
	if (ref->ret == -2) {
		fprintf(stderr, "ERROR: Chromosome %s missing from FastA file.\n", chrom);
		exit(-3);
	}
	if (ref->ret < 0 || ref->seq == NULL) {
		fprintf(stderr, "ERROR: Unable to fetch sequence.\n");
		exit(-3);
	}
	ref->end = ref->start + ref->ret - 1;

	// convert to capital letters
	for (hts_pos_t i = 0; i < ref->ret; i++)
		ref->seq[i] = toupper(ref->seq[i]);

	return ref->seq[0];
}

void write_vcf_header(FILE *vcf_fptr, sam_hdr_t *bam_hdr, struct GenotypingOptions *opts)
{
	int i;

	// mandatory fileformat line
	fputs("##fileformat=VCFv4.3\n", vcf_fptr);

	// file creation date
	time_t tt = time(NULL);
  	struct tm lt = *localtime(&tt);
	fprintf(vcf_fptr, "##fileDate=%4d%02d%02d\n", lt.tm_year+1900, lt.tm_mon+1, lt.tm_mday);
	
	// source info
	fputs("##bssnper2Version=0.1\n", vcf_fptr);
	fputs("##bssnper2Command=", vcf_fptr);
	for (i = 1; i < opts->cmd_argc-1; i++)
		fprintf(vcf_fptr, "%s ", opts->cmd_argv[i]);
	fprintf(vcf_fptr, "%s\n", opts->cmd_argv[i]);
	// chromosome reference location
	fprintf(vcf_fptr, "##reference=file:///%s\n", opts->ref_fname);
	
	// chromosome lengths
	for (i = 0; i < bam_hdr->n_targets; i++) {
		fputs("##contig=<ID=", vcf_fptr);
		fputs(bam_hdr->target_name[i], vcf_fptr);
		fputs(",length=", vcf_fptr);
		fprintf(vcf_fptr, "%d>\n", bam_hdr->target_len[i]);
	}

	// INFO fields
	fputs("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n", vcf_fptr);
	fputs("##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand\">\n", vcf_fptr);
	fputs("##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand\">\n", vcf_fptr);
	fputs("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">\n", vcf_fptr);
	
	// FILTER field
	fputs("##FILTER=<ID=Low,Description=\"Low Quality\">\n", vcf_fptr); // TODO make more descriptive

	// FORMAT fields
	fputs("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", vcf_fptr);
	fputs("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">\n", vcf_fptr);
	fputs("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">\n", vcf_fptr);
	fputs("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">\n", vcf_fptr);
	fputs("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n", vcf_fptr);
	fputs("##FORMAT=<ID=BSD,Number=8,Type=Integer,Description=\"Depth, ATCG in Watson strand and Crick strand\">\n", vcf_fptr);
	fputs("##FORMAT=<ID=BSQ,Number=8,Type=Integer,Description=\"Average Base Quality, ATCG in Watson strand and Crick strand\">\n", vcf_fptr);
	fputs("##FORMAT=<ID=ALFR,Number=R,Type=Float,Description=\"Allele frequency\">\n", vcf_fptr);

	// mandatory header line
	fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", vcf_fptr);
	fputs(opts->sample_name, vcf_fptr);
	fputs("\n", vcf_fptr);

	fflush(vcf_fptr);
}

void genotype_bam(struct GenotypingOptions *opts)
{
	int i;
	char *bam_filename = opts->bam_fname;

	// TODO implement vSnpPerBase / error_rate?
	// not sure how important this is unless the mapper is unreliable?
	
	// instantiate geno_buf 'hash' array;
	gt_buffer *geno_buf = gt_buffer_init(opts->buffer_size);
	gt_pos *geno_pos;
	
	// initialize bam file handle and header
	samFile *bamf = (samFile*)sam_open(bam_filename, "r");
    sam_hdr_t *header = sam_hdr_read(bamf);

	// initialize fasta index and reference seq cache
	RefCache ref_cache;
	ref_cache.fai = fai_load(opts->ref_fname);
	if (ref_cache.fai == NULL) {
		fprintf(stderr, "Unable to load FastA index for reference.\n");
		exit(-2);
	}
	ref_cache.size = 10000;
	ref_cache.seq = NULL;
	reset_refcache(&ref_cache);
	
	// keep track of whether a chrom has been seen yet
	bool seen_chrom[header->n_targets];
	for (i = 0; i < header->n_targets; i++) seen_chrom[i] = false;

	// open VCF and homref files for output
	FILE *vcf_fptr = fopen(opts->vcf_fname, "w");
	write_vcf_header(vcf_fptr, header, opts);
	FILE *homref_fptr = fopen(opts->homref_fname, "w");

	int32_t curr_chrom = -1;
	hts_pos_t curr_pos = -1;
	char *read_name;
	bam1_t *old_read, *b;
	b = bam_init1();
	bool new_chrom;
	char* chrom_name;

	// walk through alignment file
	while (sam_read1(bamf, header, b) >= 0) {
		// ignore unmapped reads or low-MAPQ reads
		if (b->core.flag&BAM_FUNMAP || b->core.tid < 0 || b->core.pos < 0 || b->core.qual < opts->min_mapq)
			continue; 

		read_name = bam_get_qname(b);
		new_chrom = curr_chrom != b->core.tid;

		// verify sorting by position
		if (!new_chrom && b->core.pos < curr_pos) {
			fprintf(stderr, "ERROR: BAM file not sorted [failing read: %s]\n", read_name);
			exit(1);
		}
		
		// write genotypes of positions already passed
		while (geno_buf->occupancy > 0 && (curr_pos < b->core.pos || new_chrom)) {
			geno_pos = retrieve_gt_pos(geno_buf, curr_pos, false);
			if (geno_pos != NULL) {
				while (geno_pos->num_reads > 0) {
					old_read = pop_read_pointer(geno_pos, NULL);
					fprintf(stderr, "Warning: Read %s mate not found\n", bam_get_qname(old_read));
					evaluate_read(geno_buf, old_read, curr_pos, -1, opts->min_base_qual); // start at curr_pos as the positions before will already have been counted
					bam_destroy1(old_read);
				}
				call_genotype(vcf_fptr, homref_fptr, chrom_name, geno_pos, get_ref_base(chrom_name, curr_pos, &ref_cache), opts);
				gt_pos_destroy(pop_gt_pos(geno_buf, curr_pos)); // destroy geno_pos and update geno_buf->occupancy if necessary
			}
			curr_pos++;
		}
		curr_pos = b->core.pos; // update curr_pos to where we are now
		
		// verify sorting by chrom
		if (new_chrom) {
			if (seen_chrom[curr_chrom]) {
				fprintf(stderr, "ERROR: BAM file not sorted [failing read: %s]\n", read_name);
				exit(1);
			}
			curr_chrom = b->core.tid;
			chrom_name = header->target_name[curr_chrom];
			seen_chrom[curr_chrom] = true;
			reset_refcache(&ref_cache);
		}

		if (b->core.flag&BAM_FPROPER_PAIR && mate_same_chrom(b)) {
			bool mate_later = await_mate(b);
			if (mate_later) {
				// evaluate non-overlapping bases and hash read
				evaluate_read(geno_buf, b, -1, b->core.mpos, opts->min_base_qual);
				insert_read(geno_buf, b->core.mpos, b);
				/*if (readname was already hashed) {
					fprintf(stderr, "ERROR: Read %s mates unsorted/malformed\n", read_name);
					exit(2);
				}*/
			} else {
				old_read = pop_read_pointer_gb(geno_buf, curr_pos, bam_get_qname(b));
				if (old_read != NULL) {
					// handle overlap and evaluate both mates
					tweak_overlap_quality(old_read, b);
					evaluate_read(geno_buf, b, -1, -1, opts->min_base_qual);
					evaluate_read(geno_buf, old_read, curr_pos, -1, opts->min_base_qual); 
					bam_destroy1(old_read);
				} else {
					// mate not hashed
					if (b->core.mpos == b->core.pos) {
						// await mate at same pos
						insert_read(geno_buf, b->core.mpos, b);
					} else {
						// mate already seen or doesn't overlap
						evaluate_read(geno_buf, b, -1, -1, opts->min_base_qual);
					}
				}
			}
		} else {
			// mate unmapped or different chrom
			evaluate_read(geno_buf, b, -1, -1, opts->min_base_qual);
		}

	}

	// finish genotyping any remaining uncalled positions
	while (geno_buf->occupancy > 0) {
		geno_pos = retrieve_gt_pos(geno_buf, curr_pos, false);
		//fprintf(stderr, "pos: %ld\toccupancy: %u\tgp: %p\n", curr_pos, geno_buf->occupancy, geno_pos);
		if (geno_pos != NULL) {
			while (geno_pos->num_reads > 0) {
				old_read = pop_read_pointer(geno_pos, NULL);
				fprintf(stderr, "Warning: Read %s mate not found\n", bam_get_qname(old_read));
				evaluate_read(geno_buf, old_read, curr_pos, -1, opts->min_base_qual); // start at curr_pos as the positions before will already have been counted
				bam_destroy1(old_read);
			}
			call_genotype(vcf_fptr, homref_fptr, chrom_name, geno_pos, get_ref_base(chrom_name, curr_pos, &ref_cache), opts);
			gt_pos_destroy(pop_gt_pos(geno_buf, curr_pos)); // destroy geno_pos and update geno_buf->occupancy if necessary
		}
		curr_pos++;
	}
	
	// close files and free memory
	fclose(vcf_fptr);
	fclose(homref_fptr);

	fai_destroy(ref_cache.fai);
	reset_refcache(&ref_cache);
	gt_buffer_destroy(geno_buf);
}	