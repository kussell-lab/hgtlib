//
// Created by LinMingzhi on 7/12/15.
//
#include "hgt_linkage.h"
#include <stdlib.h>

hgt_linkage * hgt_linkage_alloc() {
    hgt_linkage *l = (hgt_linkage *) malloc(sizeof(hgt_linkage));
    l->birthTime = 0;
    l->numChildren = 0;
    return l;
}

hgt_linkage * hgt_linkage_new(hgt_linkage * parent, unsigned long birthTime) {
    hgt_linkage *l = hgt_linkage_alloc();
    l->birthTime = birthTime;
    l->parent = parent;
    if (parent) {
        parent->numChildren++;
    }
    return l;
}

int hgt_linkage_free(hgt_linkage *l) {
    hgt_linkage * parent;
    parent = l->parent;
    if (l->numChildren < 1) {
        if (parent) {
            // reduce the number of children of the parent by 1,
            // and try to free the parent.
            parent->numChildren--;
            hgt_linkage_free(parent);
        }
        free(l);
    }

    return EXIT_SUCCESS;
}

int hgt_linkage_free_more(hgt_linkage **l, unsigned int size) {
    int i;
    for (i = 0; i < size; i++) {
        hgt_linkage_free(l[i]);
    }
    return  EXIT_SUCCESS;
}
unsigned long hgt_linkage_find_most_rescent_ancestor_time(hgt_linkage ** linkages, int size) {
    int bad = 0, found = 1;
    hgt_linkage * parent = linkages[0]->parent;
    int i;
    for (i = 0; i < size; ++i) {
        if (!linkages[i]->parent) {
            bad = 1;
            break;
        } else {
            if (parent != linkages[i]->parent ) {
                found = 0;
            }
        }
    }

    unsigned long maxBirthTime, birthTime;
    unsigned int maxIndex;
    if (bad == 1) {
        return 0;
    } else {
        if (found == 1) {
            birthTime = linkages[0]->birthTime;
            return birthTime;
        } else {
            // find the one that has max birth time.
            maxBirthTime = linkages[0]->birthTime;
            maxIndex = 0;
            for (i = 0; i < size; i++) {
                birthTime = linkages[i]->birthTime;
                if (maxBirthTime < birthTime) {
                    maxBirthTime = birthTime;
                    maxIndex = i;
                }
            }
            linkages[maxIndex] = linkages[maxIndex]->parent;
            return hgt_linkage_find_most_rescent_ancestor_time(linkages, size);
        }
    }

}

unsigned long hgt_linkage_find_most_rescent_coalescence_time(hgt_linkage ** linkages, int size) {
    int i, j;
    int bad, maxIndex;
    unsigned long birthTime, maxBirthTime;

    bad = 0;
    for (i = 0; i < size; ++i) {
        if (!linkages[i]->parent) {
            bad = 1;
            break;
        }
    }

    if (bad == 1) {
        return 0;
    } else {
        maxBirthTime = 0;
        for (i = 0; i < size; i ++) {
            for (j = i+1; j < size; j++) {
                if (linkages[i]->parent == linkages[j]->parent) {
                    birthTime = linkages[i]->birthTime;
                    if (maxBirthTime < birthTime) {
                        maxBirthTime = birthTime;
                    }
                }
            }
        }
        if (maxBirthTime > 0) {
            return maxBirthTime;
        } else {
            // find the one that has max birth time.
            maxBirthTime = linkages[0]->birthTime;
            maxIndex = 0;
            for (i = 0; i < size; i++) {
                if (maxBirthTime < linkages[i]->birthTime) {
                    maxBirthTime = linkages[i]->birthTime;
                    maxIndex = i;
                }
            }
            linkages[maxIndex] = linkages[maxIndex]->parent;
            return hgt_linkage_find_most_rescent_coalescence_time(linkages, size);
        }
    }
}

int hgt_linkage_prune(hgt_linkage *l) {
    hgt_linkage *parent;
    parent = l->parent;

    if (l && parent) {
        if (parent->numChildren <= 1) {
            l->parent = parent->parent;
            free(parent);
            hgt_linkage_prune(l);
        }
    }
    return EXIT_SUCCESS;
}

int hgt_linkage_prune_more(hgt_linkage **l, unsigned int size) {
    int i;
    for (i = 0; i < size; i++) {
        hgt_linkage_prune(l[i]);
    }
    return EXIT_SUCCESS;
}
