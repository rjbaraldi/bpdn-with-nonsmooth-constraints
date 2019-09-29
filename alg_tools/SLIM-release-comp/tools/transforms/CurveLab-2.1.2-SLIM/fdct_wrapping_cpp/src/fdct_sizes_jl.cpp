#include <cstdio>

#include "fdct_sizes.hpp"

static bool dbg = false;

extern "C" {
int jl_fdct_sizes_map_size(int, int, int);
void jl_fdct_sizes(int, int, int, int, int, int*, size_t*);
}

int jl_fdct_sizes_map_size(int nbscales, int nbangles_coarse, int all_crvlts){
    int map_size = 2+(nbscales-2+all_crvlts)*5+(1-all_crvlts)*2;
    return map_size;
}

void jl_fdct_sizes(int nbscales, int nbangles_coarse, int all_crvlts, int n1, int n2, int *pmap, size_t *totalcoeffs){
    int map_size = 2+(nbscales-2+all_crvlts)*5+(1-all_crvlts)*2;
    vector<int> map;
    fdct_sizes( map, nbscales, nbangles_coarse, all_crvlts, n1, n2, *totalcoeffs);
    for (int i=0;i<map_size;i++) {
        pmap[i]=map[i];
        if(dbg) printf("\t%d %d, %d\n",i,map[i],pmap[i]);
    }
}
