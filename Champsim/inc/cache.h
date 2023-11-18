#ifndef CACHE_H
#define CACHE_H
#include "memory_class.h"
#include <vector>
#include<map>

// PAGE
extern uint32_t PAGE_TABLE_LATENCY, SWAP_LATENCY;

// CACHE TYPE
#define IS_ITLB 0
#define IS_DTLB 1
#define IS_STLB 2
#define IS_L1I  3
#define IS_L1D  4
#define IS_L2C  5
#define IS_LLC  6

// INSTRUCTION TLB
#define ITLB_SET 16
#define ITLB_WAY 4
#define ITLB_RQ_SIZE 16
#define ITLB_WQ_SIZE 16
#define ITLB_PQ_SIZE 0
#define ITLB_MSHR_SIZE 8
#define ITLB_LATENCY 1

// DATA TLB
#define DTLB_SET 16
#define DTLB_WAY 4
#define DTLB_RQ_SIZE 16
#define DTLB_WQ_SIZE 16
#define DTLB_PQ_SIZE 0
#define DTLB_MSHR_SIZE 8
#define DTLB_LATENCY 1

// SECOND LEVEL TLB
#define STLB_SET 128
#define STLB_WAY 12
#define STLB_RQ_SIZE 32
#define STLB_WQ_SIZE 32
#define STLB_PQ_SIZE 0
#define STLB_MSHR_SIZE 16
#define STLB_LATENCY 8

// L1 INSTRUCTION CACHE
#define L1I_SET 64
#define L1I_WAY 8
#define L1I_RQ_SIZE 64
#define L1I_WQ_SIZE 64 
#define L1I_PQ_SIZE 32
#define L1I_MSHR_SIZE 8
#define L1I_LATENCY 4

// L1 DATA CACHE
#define L1D_SET 64
#define L1D_WAY 12
#define L1D_RQ_SIZE 64
#define L1D_WQ_SIZE 64 
#define L1D_PQ_SIZE 8
#define L1D_MSHR_SIZE 16
#define L1D_LATENCY 5 

// L2 CACHE
#define L2C_SET 1024
#define L2C_WAY 8
#define L2C_RQ_SIZE 32
#define L2C_WQ_SIZE 32
#define L2C_PQ_SIZE 16
#define L2C_MSHR_SIZE 32
#define L2C_LATENCY 10  // 4/5 (L1I or L1D) + 10 = 14/15 cycles

// LAST LEVEL CACHE
#define LLC_SET NUM_CPUS*2048
#define LLC_WAY 16
#define LLC_RQ_SIZE NUM_CPUS*L2C_MSHR_SIZE //48
#define LLC_WQ_SIZE NUM_CPUS*L2C_MSHR_SIZE //48
#define LLC_PQ_SIZE NUM_CPUS*32
#define LLC_MSHR_SIZE NUM_CPUS*64
#define LLC_LATENCY 20  // 4/5 (L1I or L1D) + 10 + 20 = 34/35 cycles

extern uint32_t count_load_hit, count_load_miss, count_store_hit, count_store_miss, count_else_hit, count_else_miss;//saurabh

class CACHE : public MEMORY {
  public:
    uint32_t cpu;
    const string NAME;
    const uint32_t NUM_SET, NUM_WAY, NUM_LINE, WQ_SIZE, RQ_SIZE, PQ_SIZE, MSHR_SIZE;
    uint32_t LATENCY;
    BLOCK **block;
    int fill_level;
    uint32_t MAX_READ, MAX_FILL;
    uint32_t reads_available_this_cycle;
    uint8_t cache_type;
    // Nupur, Antara
    int BLOCKS_TO_PREFETCH_A = 1;
    int READ_REQ_DELAY = 0;
    
    vector<vector <uint64_t>> return_pkt_mshr; //saurabh
    map<unsigned,unsigned> evected_block;

    // prefetch stats
    uint64_t pf_requested,
             pf_issued,
             pf_useful,
             pf_useless,
             //Nupur, Antara
             POLLUTION,
             WRITEBACK_POLLUTION,
             OUR_MSHR_COUNT,
             OUR_RQ_COUNT,
             DIRTY_EVICTS,
             
             pf_fill,
             pf_usefull_on_time, pf_usefull_late,
             pf_late,
             PQ_over_flow,
             outside_baundary,
             counter_0, counter_1, counter_2, counter_3, counter_4, counter_5, counter_6, counter_7, counter_8, counter_9, counter_10, counter_11;
    uint32_t count_load_hit, count_load_miss, count_store_hit, count_store_miss, count_else_hit, count_else_miss;//saurabh
    

    // queues
    PACKET_QUEUE WQ{NAME + "_WQ", WQ_SIZE}, // write queue
                 RQ{NAME + "_RQ", RQ_SIZE}, // read queue
                 PQ{NAME + "_PQ", PQ_SIZE}, // prefetch queue
                 MSHR{NAME + "_MSHR", MSHR_SIZE}, // MSHR
                 PROCESSED{NAME + "_PROCESSED", ROB_SIZE}; // processed queue

    float Avg[11];
    uint64_t sim_access[NUM_CPUS][NUM_TYPES],
             sim_hit[NUM_CPUS][NUM_TYPES],
             sim_miss[NUM_CPUS][NUM_TYPES],
             my_stat[NUM_CPUS][NUM_TYPES],
             //sim_hit_rate[NUM_CPUS][NUM_TYPES],   //saurabh
             //AMAT[NUM_CPUS][NUM_TYPES],
             LateBy[11],counter[11],
             roi_access[NUM_CPUS][NUM_TYPES],
             roi_hit[NUM_CPUS][NUM_TYPES],
             roi_miss[NUM_CPUS][NUM_TYPES],
             roi_my_stat[NUM_CPUS][NUM_TYPES];
             //roi_hit_rate[NUM_CPUS][NUM_TYPES],   //saurabh
             //roi_AMAT[NUM_CPUS][NUM_TYPES];

    uint64_t total_miss_latency;
    //uint64_t temp2; //saurabh
    //uint64_t New_Total = 0;
    // constructor
    CACHE(string v1, uint32_t v2, int v3, uint32_t v4, uint32_t v5, uint32_t v6, uint32_t v7, uint32_t v8) 
        : NAME(v1), NUM_SET(v2), NUM_WAY(v3), NUM_LINE(v4), WQ_SIZE(v5), RQ_SIZE(v6), PQ_SIZE(v7), MSHR_SIZE(v8) {

        LATENCY = 0;

        // cache block
        block = new BLOCK* [NUM_SET];
        for (uint32_t i=0; i<NUM_SET; i++) {
             // last element keeps the cache pollution counter.
            block[i] = new BLOCK[NUM_WAY+4*NUM_WAY + 1]; 

            for (uint32_t j=0; j<NUM_WAY; j++) {
                block[i][j].lru = j;
            }
        }

        for (uint32_t i=0; i<=10; i++)
        {
            LateBy[i] = 0;
            counter[i] = 0;
            Avg[i] = 0;
        }

        for (uint32_t i=0; i<NUM_CPUS; i++) {
            upper_level_icache[i] = NULL;
            upper_level_dcache[i] = NULL;

            for (uint32_t j=0; j<NUM_TYPES; j++) {
                sim_access[i][j] = 0;
                sim_hit[i][j] = 0;
                sim_miss[i][j] = 0;
                //sim_hit_rate[i][j] = 0;  //saurabh
                //AMAT[i][j] = 0; //saurabh
                roi_access[i][j] = 0;
                roi_hit[i][j] = 0;
                roi_miss[i][j] = 0;
                my_stat[i][j] = 0;
                roi_my_stat[i][j] = 0;
                //roi_hit_rate[i][j] = 0;   //saurabh
                //roi_AMAT[i][j] = 0;
            }
        }

	total_miss_latency = 0;
    
        lower_level = NULL;
        extra_interface = NULL;
        fill_level = -1;
        MAX_READ = 1;
        MAX_FILL = 1;

        pf_requested = 0;
        //Nupur, Antara
        POLLUTION= 0;
        WRITEBACK_POLLUTION = 0;
        OUR_MSHR_COUNT = 0;
        OUR_RQ_COUNT = 0;
        DIRTY_EVICTS = 0;

        pf_issued = 0;
        pf_useful = 0;
        pf_useless = 0;
        pf_fill = 0;
        PQ_over_flow=0;
        outside_baundary=0;
        counter_0 = 0;
        counter_1 = 0;
        counter_2 = 0;
        counter_3 = 0;
        counter_4 = 0;
        counter_5 = 0;
        counter_6 = 0;
        counter_7 = 0;
        counter_8 = 0;
        counter_9 = 0;
        counter_10 = 0;
        counter_11 = 0;
        pf_usefull_on_time = 0;
        pf_usefull_late = 0;
        pf_late = 0;        
    };

    // destructor
    ~CACHE() {
        for (uint32_t i=0; i<NUM_SET; i++)
            delete[] block[i];
        delete[] block;
    };

    // functions
    int  add_rq(PACKET *packet),
         add_wq(PACKET *packet),
         add_pq(PACKET *packet);

    void return_data(PACKET *packet),
         operate(),
         increment_WQ_FULL(uint64_t address);

    uint32_t get_occupancy(uint8_t queue_type, uint64_t address),
             get_size(uint8_t queue_type, uint64_t address);

    int  check_hit(PACKET *packet),
         invalidate_entry(uint64_t inval_addr),
         check_mshr(PACKET *packet),
         prefetch_line(uint64_t ip, uint64_t base_addr, uint64_t pf_addr, int prefetch_fill_level, uint32_t prefetch_metadata),
         kpc_prefetch_line(uint64_t base_addr, uint64_t pf_addr, int prefetch_fill_level, int delta, int depth, int signature, int confidence, uint32_t prefetch_metadata);

    void handle_fill(),
         handle_writeback(),
         handle_read(),
         handle_prefetch();

    void add_mshr(PACKET *packet),
         update_fill_cycle(),         
         llc_initialize_replacement(),
         update_replacement_state(uint32_t cpu, uint32_t set, uint32_t way, uint64_t full_addr, uint64_t ip, uint64_t victim_addr, uint32_t type, uint8_t hit),
         llc_update_replacement_state(uint32_t cpu, uint32_t set, uint32_t way, uint64_t full_addr, uint64_t ip, uint64_t victim_addr, uint32_t type, uint8_t hit),
         lru_update(uint32_t set, uint32_t way),
         fill_cache(uint32_t set, uint32_t way, PACKET *packet),
         replacement_final_stats(),
         llc_replacement_final_stats(),
         //prefetcher_initialize(),
         l1d_prefetcher_initialize(),
         l2c_prefetcher_initialize(),
         llc_prefetcher_initialize(uint32_t cpu),
         prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type),
         l1d_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type),
         prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr),
         l1d_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in),
         //prefetcher_final_stats(),
         l1d_prefetcher_final_stats(),
         l2c_prefetcher_final_stats(),
         llc_prefetcher_final_stats();
    void (*l1i_prefetcher_cache_operate)(uint32_t, uint64_t, uint8_t, uint8_t);
    void (*l1i_prefetcher_cache_fill)(uint32_t, uint64_t, uint32_t, uint32_t, uint8_t, uint64_t);

    uint32_t l2c_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in),
        //  llc_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t metadata_in),
        l2c_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in),
        llc_prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t metadata_in);
    
    uint32_t get_set(uint64_t address),
             get_way(uint64_t address, uint32_t set),
             find_victim(uint32_t cpu, uint64_t instr_id, uint32_t set, const BLOCK *current_set, uint64_t ip, uint64_t full_addr, uint32_t type),
             llc_find_victim(uint32_t cpu, uint64_t instr_id, uint32_t set, const BLOCK *current_set, uint64_t ip, uint64_t full_addr, uint32_t type),
             lru_victim(uint32_t cpu, uint64_t instr_id, uint32_t set, const BLOCK *current_set, uint64_t ip, uint64_t full_addr, uint32_t type);

    // for vldp - Nupur, Antara
    void llc_prefetcher_cache_fill_(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t cpu);
    void llc_prefetcher_inform_warmup_complete_();
    void llc_prefetcher_inform_roi_complete_(uint32_t cpu);
    void llc_prefetcher_roi_stats_(uint32_t cpu);
    void llc_prefetcher_final_stats_(uint32_t cpu);
    void llc_prefetcher_initialize_(uint32_t cpu);
    uint32_t llc_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t cpu);

    // for vldp l2c - Nupur, Antara
    void l2c_prefetcher_cache_fill_(uint64_t addr, uint32_t set, uint32_t way, uint8_t prefetch, uint64_t evicted_addr, uint32_t cpu);
    void l2c_prefetcher_inform_warmup_complete_();
    void l2c_prefetcher_inform_roi_complete_(uint32_t cpu);
    void l2c_prefetcher_roi_stats_(uint32_t cpu);
    void l2c_prefetcher_final_stats_(uint32_t cpu);
    void l2c_prefetcher_initialize_(uint32_t cpu);
    // uint32_t l2c_prefetcher_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit, uint8_t type, uint32_t cpu);
    //for bandwidth - Nupur, Antara
    void bandwidth_calc(uint32_t cpu, CACHE *cache, uint64_t DRAM_RESPONSES);
    // void pol_based_bandwidth_calc(uint32_t cpu, CACHE *cache, uint64_t DRAM_RESPONSES);
    void pol_based_bandwidth_calc(uint32_t cpu, CACHE *cache, CACHE *cache_l2c, uint64_t DRAM_RESPONSES);
    void poll_based_pref_degree(uint32_t cpu, CACHE *cache);
    void pf_hit_rate_based(uint32_t cpu, CACHE *cache);
    void timeliness_based_pref_degree();
};

#endif
