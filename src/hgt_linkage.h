//
// Created by LinMingzhi on 7/12/15.
//

#ifndef HGT_HGT_LINKAGE_H
#define HGT_HGT_LINKAGE_H
typedef struct _hgt_linkage hgt_linkage;
struct _hgt_linkage {
    unsigned  int numChildren;
    unsigned long birthTime;
    hgt_linkage * parent;
};
hgt_linkage * hgt_linkage_alloc();
hgt_linkage * hgt_linkage_new(hgt_linkage * parent, unsigned long birthTime);
int hgt_linkage_free(hgt_linkage *l);
int hgt_linkage_free_more(hgt_linkage **l, unsigned int size);
int hgt_linkage_prune(hgt_linkage *l);
int hgt_linkage_prune_more(hgt_linkage **l, unsigned int size);
typedef unsigned long hgt_linkage_find_time_func(hgt_linkage ** linkages, int size);
unsigned long hgt_linkage_find_most_rescent_ancestor_time(hgt_linkage ** linkages, int size);
unsigned long hgt_linkage_find_most_rescent_coalescence_time(hgt_linkage ** linkages, int size) ;
#endif //HGT_HGT_LINKAGE_H
