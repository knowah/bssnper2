#include "uncalled_pos.h"

void free_read_ptr_node_chain(read_ptr_node *rpn)
{
    free(rpn->qname);
    bam_destroy1(rpn->read_ptr);
    if (rpn->next != NULL) free_read_ptr_node_chain(rpn->next);
    free(rpn);
    return;
}

void read_ptr_node_destroy(read_ptr_node *rpn, bool destroy_read_ptr)
{
    free(rpn->qname);
    if (destroy_read_ptr) bam_destroy1(rpn->read_ptr);
    rpn->next = NULL;
    free(rpn);
}

read_ptr_node *read_ptr_node_init(bam1_t *read_ptr)
{
    read_ptr_node *rpn = malloc(sizeof(read_ptr_node));
    if (rpn == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for read_ptr_node\n");
        exit(3);
    }
    rpn->qname = malloc((read_ptr->core.l_qname+1)*sizeof(char));
    strcpy(rpn->qname, bam_get_qname(read_ptr));
    rpn->read_ptr = bam_dup1(read_ptr);
    rpn->next = NULL;
    return rpn;
}

gt_pos *gt_pos_init(hts_pos_t pos)
{
    gt_pos *gp = malloc(sizeof *gp);
    gp->pos = pos;
    for (int i = 0; i < 8; i++) {
        gp->call_counts[i] = 0;
        gp->qual_sums[i] = 0;
    }
    gp->num_reads = 0;
    gp->reads = NULL;
    gp->next = NULL;
    return gp;
}

// NB: qname==NULL will pop first read pointer regardless of content
bam1_t *pop_read_pointer(gt_pos *gp, const char *qname)
{
    if (gp == NULL || gp->reads == NULL) return NULL;

    read_ptr_node *rpn = gp->reads;
    bam1_t *read = NULL;
    if (qname == NULL || strcmp(rpn->qname, qname) == 0) {
        gp->reads = rpn->next;
        read = rpn->read_ptr;
        read_ptr_node_destroy(rpn, false);
        gp->num_reads -= 1;
        return read;
    }
    while (rpn->next != NULL) {
        if (strcmp(rpn->next->qname, qname) == 0) {
            read = rpn->next->read_ptr;
            read_ptr_node *tmp = rpn->next->next;
            read_ptr_node_destroy(rpn->next, false);
            rpn->next = tmp;
            gp->num_reads -= 1;
            return read;
        }
        rpn = rpn->next;
    }
    return NULL;
}


gt_pos *retrieve_gt_pos(gt_buffer *gb, hts_pos_t pos, bool create)
{
    gt_hash_t offset = gt_buffer_offset(gb, pos);
    gt_pos *gp = gb->buf[offset];

    // no entry in hash table
    if (gp == NULL) { 
        if (create) {
            gb->buf[offset] = gt_pos_init(pos);      
            gb->occupancy += 1;
            return gb->buf[offset];
        } else return NULL;
    }

    // first entry is the target pos
    if (gp->pos == pos)
        return gp;

    // first entry is after the target pos
    if (gp->pos > pos) {
        if (create) {
            gp = gt_pos_init(pos);
            gp->next = gb->buf[offset];
            gb->buf[offset] = gp;
            return gp;
        } else return NULL;
    }

    // navigate through list of entries, stopping at the entry
    // immediatlely prior to where the target pos belongs
    while (gp != NULL && gp->next != NULL && gp->next->pos <= pos)
        gp = gp->next;
    
    // current entry is the target pos
    if (gp->pos == pos)
        return gp;

    // target pos should be after current entry
    if (create) {            
        gt_pos *newgt = gt_pos_init(pos);
        newgt->next = gp->next;
        gp->next = newgt;
        return newgt;
    }

    return NULL;
}

void insert_read(gt_buffer *gb, hts_pos_t pos, bam1_t *read_ptr)
{
    gt_pos *gp = retrieve_gt_pos(gb, pos, true);
    read_ptr_node *rpn = read_ptr_node_init(read_ptr);
    read_ptr_node *prev = gp->reads;
    if (prev == NULL) {
        gp->reads = rpn;
    } else {
        while (prev->next != NULL) prev = prev->next;
        prev->next = rpn;
    }
    gp->num_reads += 1;
    return;
}

void gt_pos_destroy(gt_pos *gp)
{
    if (gp->reads != NULL) {
        fprintf(stderr, "Warning: destroying gt_pos with remaining reads.\n");
        free_read_ptr_node_chain(gp->reads);
    }
    gp->next = NULL;
    free(gp);
    return;
}



gt_buffer *gt_buffer_init(gt_hash_t size)
{
    gt_buffer *gb;
    gb = malloc(sizeof *gb);
    //gb = malloc((sizeof *gb) + size * sizeof(gt_pos*));
    gb->size = size;
    gb->occupancy = 0;
    //for (gt_hash_t i = 0; i < size; i++)
    //    gb->buf[i] = NULL;
    gb->buf = calloc(size, sizeof(gt_pos*)); // initialise buffer with 0 / NULL pointer
    return gb;
}

gt_pos *pop_gt_pos(gt_buffer *gb, hts_pos_t pos)
{
    // this function only returns the desired gt_pos if
    // the gt_buffer enforced ordering when inserting 
    // gt_pos pointers
    gt_hash_t offset = gt_buffer_offset(gb, pos);
    if (gb->buf[offset] == NULL) return NULL;
    gt_pos *ret = gb->buf[offset];
    if (ret->pos != pos) return NULL;
    gb->buf[offset] = ret->next;
    if (gb->buf[offset] == NULL) gb->occupancy -= 1;
    return ret;
}

void free_gt_pos_chain(gt_pos *gp)
{   
    if (gp == NULL) return;
    if (gp->next != NULL) free_gt_pos_chain(gp->next);
    gt_pos_destroy(gp);
    return;
}

void gt_buffer_destroy(gt_buffer *gb)
{
    for (gt_hash_t i = 0; i < gb->size; i++)
        free_gt_pos_chain(gb->buf[i]);
    free(gb->buf);
    free(gb);
    return;
}