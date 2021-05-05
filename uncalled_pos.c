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
void bc_zero(basecalls *p)
{
    for (int i = 0; i < 8; i++) {
        p->call_counts[i] = 0;
        p->qual_sums[i] = 0;
    }
}

gt_pos *gt_pos_init(hts_pos_t pos)
{
    gt_pos *gp = malloc(sizeof *gp);
    gp->pos = pos;
    /*
    for (int i = 0; i < 8; i++) {
        gp->bc.call_counts[i] = 0;
        gp->bc.qual_sums[i] = 0;
    }
    */
    bc_zero(&(gp->bc));
    gp->num_reads = 0;
    gp->ins = NULL;
    gp->del = NULL;
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

del_dict* del_dict_init()
{
    del_dict *dd = malloc(sizeof(del_dict));
    dd->total_reads = 0;
    for (int i = 0; i < INDEL_HASHSIZE; i++) dd->entries[i] = NULL;
    return dd;
}

del_entry* del_entry_init(uint32_t del_len)
{
    del_entry *e = malloc(sizeof(del_entry));
    e->count = 0;
    e->length = del_len;
    e->next = NULL;
    return e;
}

del_entry* retrieve_del_entry(del_dict *dd, uint32_t del_len)
{
    del_entry *temp, *e;
    uint32_t index = del_len % INDEL_HASHSIZE;
    e = dd->entries[index];

    // if entry non-existent, create it
    if (e == NULL) {
        dd->entries[index] = del_entry_init(del_len);
        return dd->entries[index];
    }

    // entry exists and is the correct one
    if (e->length == del_len) return e;
    
    if (e->length > del_len) {
        // existing entry in hash is for a larger del_len
        temp = e;
        e = del_entry_init(del_len);
        e->next = temp;
    } else {
        // initial entry is for a smaller del_len
        while (e->next != NULL && e->next->length < del_len) e = e->next;
        if (e->length != del_len) {
            del_entry *temp = e->next;
            e->next = del_entry_init(del_len);
            e = e->next;
            e->next = temp;
        }
    }

    return e;
}

void increment_del(gt_pos *gp, uint32_t del_len)
{
    if (gp->del == NULL)
        gp->del = del_dict_init();
    del_entry *e = retrieve_del_entry(gp->del, del_len);
    e->count += 1;
    gp->del->total_reads += 1;
}

ins_entry* ins_entry_init(uint32_t ins_len)
{
    ins_entry *e = malloc(sizeof(ins_entry));
    e->count = 0;
    e->length = ins_len;
    e->bc_list = malloc(ins_len * sizeof(basecalls));
    for (uint32_t i = 0; i < ins_len; i++)
        bc_zero(&(e->bc_list[i]));
    e->next = NULL;
    return e;
}

ins_dict* ins_dict_init()
{
    ins_dict *id = malloc(sizeof(ins_dict));
    id->total_reads = 0;
    for (int i = 0; i < INDEL_HASHSIZE; i++) id->entries[i] = NULL;
    return id;
}


ins_entry* retrieve_ins_entry(gt_pos *gp, uint32_t ins_len)
{
    ins_entry *temp, *e;
    uint32_t index = ins_len % INDEL_HASHSIZE;
    if (gp->ins == NULL) {
        gp->ins = ins_dict_init();
    }
    
    e = gp->ins->entries[index];

    // if entry non-existent, create it
    if (e == NULL) {
        gp->ins->entries[index] = ins_entry_init(ins_len);
        return gp->ins->entries[index];
    }

    // entry exists and is the correct one
    if (e->length == ins_len) return e;
    
    if (e->length > ins_len) {
        // existing entry in hash is for a larger ins_len
        temp = e;
        e = ins_entry_init(ins_len);
        e->next = temp;
    } else {
        // initial entry is for a smaller ins_len
        while (e->next != NULL && e->next->length < ins_len) e = e->next;
        if (e->length != ins_len) {
            ins_entry *temp = e->next;
            e->next = ins_entry_init(ins_len);
            e = e->next;
            e->next = temp;
        }
    }

    return e;
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

void ins_entry_destroy(ins_entry *e)
{
    free(e->bc_list);
    free(e);
}

void ins_dict_destroy(ins_dict *id)
{
    ins_entry *e, *temp;
    for (int i = 0; i < INDEL_HASHSIZE; i++) {
        e = id->entries[i];
        while (e != NULL) {
            temp = e->next;
            ins_entry_destroy(e);
            e = temp;
        }
    }
    free(id);
}

void del_dict_destroy(del_dict *dd)
{
    del_entry *e, *temp;
    for (int i = 0; i < INDEL_HASHSIZE; i++) {
        e = dd->entries[i];
        while (e != NULL) {
            temp = e->next;
            free(e);
            e = temp;
        }
    }
    free(dd);

}

void gt_pos_destroy(gt_pos *gp)
{
    if (gp->reads != NULL) {
        fprintf(stderr, "Warning: destroying gt_pos with remaining reads.\n");
        free_read_ptr_node_chain(gp->reads);
    }
    gp->next = NULL;
    if (gp->ins != NULL) ins_dict_destroy(gp->ins);
    if (gp->del != NULL) del_dict_destroy(gp->del);
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