#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "genotyping.h"
#include "bssnper2.h"
#include <libgen.h>

int main(int argc, char **argv)
{
    // set defaults
    struct GenotypingOptions options;
    options.homref_fname = "/dev/null";
    options.indel_fname = "/dev/null";
    options.min_base_qual = 13;
    options.min_depth = 1;
    options.max_depth = UINT32_MAX;
    options.min_hom_freq = 0.85;
    options.min_het_freq = 0.1;
    options.error_rate = 0.02;
    options.min_mapq = 13;
    options.min_alt_count = 2;
    options.buffer_size = 1000;
    options.assume_homref = false;
    options.homref_in_vcf = false;
    options.sample_name = NULL;
    options.cmd_argc = argc;
    options.cmd_argv = argv;

    // parse arguments
    int opt, optIdx;
    struct option longOpts[] = 
    {
        {"ref", required_argument, NULL, 'r'},
        {"vcf", required_argument, NULL, 'v'},
        {"homf", required_argument, NULL, 'H'},
        {"indelf", required_argument, NULL, 'I'},
        {"minCoverage", required_argument, NULL, 'c'},
        {"maxCoverage", required_argument, NULL, 'C'},
        {"minBaseQ", required_argument, NULL, 'b'},
        {"minMapQ", required_argument, NULL, 'q'},
        {"minHomFreq", required_argument, NULL, 'f'},
        {"minHetFreq", required_argument, NULL, 'F'},
        {"errorRate", required_argument, NULL, 'e'},
        {"minAltCoverage", required_argument, NULL, 'a'},
        {"bufferSize", required_argument, NULL, 'i'},
        {"assumeHomref", no_argument, NULL, 'A'},
        {"homrefInVCF", no_argument, NULL, 'V'},
        {"sampleName", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}
    };
    
    while ((opt = getopt_long(argc, argv, "r:v:H:I:c:C:b:q:f:F:e:a:i:AVs:", longOpts, &optIdx)) != -1)
    {
        switch (opt)
        {
            case 'r': options.ref_fname = optarg; break;
            case 'v': options.vcf_fname = optarg; break;
            case 'H': options.homref_fname = optarg; break;
            case 'I': options.indel_fname = optarg; break;
            case 'c': options.min_depth = atoi(optarg); break;
            case 'C': options.max_depth = atoi(optarg); break;
            case 'b': options.min_base_qual = atoi(optarg); break;
            case 'q': options.min_mapq = atoi(optarg); break;
            case 'f': options.min_hom_freq = atof(optarg); break;
            case 'F': options.min_het_freq = atof(optarg); break;
            case 'e': options.error_rate = atof(optarg); break;
            case 'a': options.min_alt_count = atoi(optarg); break;
            case 'i': options.buffer_size = atoi(optarg); break;
            case 'A': options.assume_homref = true; break;
            case 'V': options.homref_in_vcf = true; break;
            case 's': options.sample_name = optarg; break;
            case '?':
                fprintf(stderr, "ERROR: Invalid argument.\n");
                exit(-1);
        }
    }
 
    // only expecting one positional argument
    if (optind == argc)
    {
        fprintf(stderr, "ERROR: No BAM file specified.\n");
        exit(-1);
    }
    if (optind != argc-1)
    {
        fprintf(stderr, "ERROR: Too many arguments provided:");
        for (int i = optind+1; i < argc; i++) fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n");
        exit(-1);
    }
    options.bam_fname = argv[optind];
    
    if (options.ref_fname == NULL) { fprintf(stderr, "ERROR: Missing required argument -r/--ref (reference FastA file)\n"); exit(-1); }
    if (options.vcf_fname == NULL) { fprintf(stderr, "ERROR: Missing required argument -v/--vcf (output VCF file)\n"); exit(-1); }
    //if (options.homref_fname == NULL) { fprintf(stderr, "ERROR: Missing required argument -H/--homf (output homref file)\n"); exit(-1); }
    //if (options.homref_fname == NULL) { options.homref_fname = "/dev/null"; } // TODO: handle this better
    char *dummy = NULL;
    if (options.sample_name == NULL) {
        dummy = strdup(options.bam_fname);
        options.sample_name = basename(dummy);
    }
    fprintf(stderr, "bam: %s\n", options.bam_fname);
    fprintf(stderr, "ref: %s\n", options.ref_fname);
    fprintf(stderr, "vcf: %s\n", options.vcf_fname);
    fprintf(stderr, "homref: %s\n", options.homref_fname);
    fprintf(stderr, "indelf: %s\n", options.indel_fname);
    fprintf(stderr, "minBaseQual: %u\n", options.min_base_qual);
    fprintf(stderr, "minCount: %u\n", options.min_depth);
    fprintf(stderr, "maxCount: %u\n", options.max_depth);
    fprintf(stderr, "minHomFreq: %f\n", options.min_hom_freq);
    fprintf(stderr, "minHetFreq: %f\n", options.min_het_freq);
    fprintf(stderr, "errorRate: %f\n", options.error_rate);
    fprintf(stderr, "minMapQual: %u\n", options.min_mapq);
    fprintf(stderr, "minAltCount: %u\n", options.min_alt_count);
    fprintf(stderr, "bufferSize: %u\n", options.buffer_size);
    fprintf(stderr, "assumeHomref: %s\n", options.assume_homref ? "true" : "false");
    fprintf(stderr, "homrefInVCF: %s\n", options.homref_in_vcf ? "true" : "false");
    
    // do genotyping
    genotype_bam(&options);
    
    if (dummy != NULL) free(dummy);
    return 0;
}
