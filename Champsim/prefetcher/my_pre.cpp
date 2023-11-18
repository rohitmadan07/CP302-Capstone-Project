#include "cache.h"
#include "uncore.h"
#define _BSD_SOURCE
#include <getopt.h>
#include "ooo_cpu.h"
#include <fstream>

long long int PREV_SUM = 0;
long long int counter = 0;
int flag = 1;
float epoch_sum_bwd = 0;

void CACHE::llc_prefetcher_initialize() 
{
    cout << "LLC Next Line Prefetcher" << endl;
}

uint32_t CACHE::llc_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in)
{
    long long int SUM = 0;
    long long int CURR_SUM = 0;
    for (uint32_t i=0; i<DRAM_CHANNELS; i++) {
        SUM = SUM + uncore.DRAM.RQ[i].ROW_BUFFER_HIT + uncore.DRAM.WQ[i].ROW_BUFFER_HIT + uncore.DRAM.RQ[i].ROW_BUFFER_MISS + uncore.DRAM.WQ[i].ROW_BUFFER_MISS; 
    }
    CURR_SUM = SUM - PREV_SUM;
    // float crr_bdw = ((CURR_SUM*64)/(32*1000000*64))*100; //in percentage
    double crr_bdw = (double)(SUM*64)/(double)(current_core_cycle[0] - ooo_cpu[0].begin_sim_cycle)*100;
    epoch_sum_bwd = epoch_sum_bwd + crr_bdw;
    counter = (ooo_cpu[0].num_retired - ooo_cpu[0].begin_sim_instr);
    PREV_SUM = SUM;

    cout<<"Curr bwd = "<<crr_bdw<< " SUM "<< SUM << " cyecles " << (current_core_cycle[0] - ooo_cpu[0].begin_sim_cycle)<< endl;

    if(counter % 10000 == 0){ 
        float epoch_bwd = epoch_sum_bwd / 10000;
        cout<<"epoch bwd = "<<epoch_bwd<<endl;
        if(epoch_bwd <= 20){
            flag = 2;
        }
        else{
            flag = 1;
        }
        epoch_bwd = 0;
        epoch_sum_bwd = 0;
    }
    for(int k=1; k<=flag; k++){
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
}