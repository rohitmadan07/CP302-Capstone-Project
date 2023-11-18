#include "cache.h"
#include "uncore.h"
#define _BSD_SOURCE
#include <getopt.h>
#include "ooo_cpu.h"
#include <fstream>
// #include "cache.cc"

long long int PREV_SUM = 0;
long long int counter = 0;
int flag = 1;
float till_prev_epoch_total_hits = 0;
CACHE* cache = &uncore.LLC;

void CACHE::llc_prefetcher_initialize() 
{
    cout << "LLC Next Line Prefetcher" << endl;
}

uint32_t CACHE::llc_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in)
{
    // long long int total_hits = cache->pf_useful;
    // counter = (ooo_cpu[0].num_retired - ooo_cpu[0].begin_sim_instr);

    // if(counter % 10000 == 0){ 
    //     long long int epoch_hits = (total_hits - till_prev_epoch_total_hits);
    //     if(epoch_hits >= 5000 ){
    //         flag = 2;
    //     }
    //     else{
    //         flag = 1;
    //     }
    //     till_prev_epoch_total_hits = total_hits;
    // }
    // int blocks_to_prefetch = num_of_blocks_to_prefetch();
    for(int k=1; k<=BLOCKS_TO_PREFETCH_A; k++){
        uint64_t pf_addr = ((addr>>LOG2_BLOCK_SIZE)+k) << LOG2_BLOCK_SIZE;
        prefetch_line(ip, addr, pf_addr, FILL_LLC, 0);
    }
    return metadata_in;
}

uint32_t CACHE::llc_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in)
{
  return metadata_in;
}

void CACHE::llc_prefetcher_final_stats()
{
  cout << "LLC Next Line Prefetcher Final Stats: none" << endl;
  // cout << "Writebacks to below level" << writesbacks_for_below_level << endl;
}