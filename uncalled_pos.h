
#ifndef BSSNPER2_UNCALLED_POS_H
#define BSSNPER2_UNCALLED_POS_H

#include <htslib/sam.h>
#include <stdbool.h>

#define INDEL_HASHSIZE 10

typedef struct read_ptr_node read_ptr_node;
struct read_ptr_node
{
    char *qname;
    bam1_t *read_ptr;
    read_ptr_node *next;
};

typedef struct del_entry del_entry;
struct del_entry
{
    uint32_t length;
    int count;
    del_entry *next;
};

typedef struct 
{
    uint32_t call_counts[8];
    uint32_t qual_sums[8];
} basecalls;

typedef struct ins_entry ins_entry;
struct ins_entry
{
    uint32_t length;
    int count;
    basecalls *bc_list;
    ins_entry *next;
};

typedef struct
{
    del_entry *entries[INDEL_HASHSIZE];
    int total_reads;
    int occupancy;
} del_dict;

typedef struct
{
    ins_entry *entries[INDEL_HASHSIZE];
    int total_reads;
    int occupancy;
} ins_dict;


typedef struct gt_pos gt_pos;
struct gt_pos
{
    hts_pos_t pos;
    //uint32_t call_counts[8];
    //uint32_t qual_sums[8];
    basecalls bc;
    int num_reads;
    ins_dict *ins;
    del_dict *del;
    read_ptr_node *reads;
    gt_pos *next;
};

typedef uint32_t gt_hash_t;

typedef struct
{
    gt_hash_t size;
    gt_hash_t occupancy;
    gt_pos **buf;
} gt_buffer;

#define gt_buffer_offset(gb, p) ((gt_hash_t)(p)%((gb)->size))

bam1_t *pop_read_pointer(gt_pos *gp, const char *qname);
gt_pos *retrieve_gt_pos(gt_buffer *gb, hts_pos_t, bool create);
#define pop_read_pointer_gb(g,p,q) ((bam1_t*)pop_read_pointer(retrieve_gt_pos((g),(p),false),(q)))
void insert_read(gt_buffer *gb, hts_pos_t pos, bam1_t *read_ptr);
gt_buffer *gt_buffer_init(gt_hash_t size);
gt_pos *pop_gt_pos(gt_buffer *gb, hts_pos_t pos);
void gt_pos_destroy(gt_pos *gp);
void gt_buffer_destroy(gt_buffer *gb);
void increment_del(gt_pos *gp, uint32_t del_len);
ins_entry* retrieve_ins_entry(gt_pos *gp, uint32_t ins_len);
del_entry* retrieve_del_entry(del_dict *dd, uint32_t del_len);
int gp_count(gt_pos *gp);

#endif 