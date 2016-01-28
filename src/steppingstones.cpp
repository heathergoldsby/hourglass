//
//  steppingstones.cpp
//  hourglass
//
//  Created by Heather Goldsby on 1/11/16.
//  Copyright © 2016 Michigan State University. All rights reserved.
//


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>

#include <ea/mkv/markov_network_evolution.h>
#include <ea/generational_models/moran_process.h>
#include <ea/fitness_function.h>
#include <ea/cmdline_interface.h>
#include <ea/datafiles/fitness.h>
#include "markov_movie.h"
#include "hourglass.h"

using namespace std;
using namespace boost::accumulators;
using namespace ealib;
using namespace mkv;


LIBEA_MD_DECL(W1, "ea.hourglass.steppingstones.w1", double); //
LIBEA_MD_DECL(W2, "ea.hourglass.steppingstones.w2", double); //
LIBEA_MD_DECL(W3, "ea.hourglass.steppingstones.w3", double); //
LIBEA_MD_DECL(W4, "ea.hourglass.steppingstones.w4", double); //




template <typename EA>
double outline(int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {

    
    double f1_01 = 1.0;
    double f2_11 = 1.0;
    
    
    // Compute fitness.
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        if (p == -1) {
            continue;
        }
        
        if ((agent_x ==0) || (agent_y == 0) || (agent_x == (max_x -1)) || (agent_y == (max_y -1))) {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                ++f1_01;
            }
        } else {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                ++f2_11;
            }
        }
        
    }
    
    //f = f1_01 * f2_11;
    double e = f1_01 + f2_11 - 2;
    double f = pow(1.5, e);
    
    return f;
}


template <typename EA>
double outline_split(int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {

    
    double f1_01 = 1.0;
    double f2_11 = 1.0;
    double f3_10 = 1.0;
    
    
    // Compute fitness.
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        if (p == -1) {
            continue;
        }
        
        
        
        if (agent_x < (floor(max_x) / 2 )) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                ++f3_10;
            }
            
        } else if (agent_x >= (floor(max_x) / 2 )) {
            if ((agent_x ==0) || (agent_y == 0) || (agent_x == (max_x -1)) || (agent_y == (max_y -1))) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f1_01;
                }
            } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                ++f2_11;
            }
        }
        
        
        
        
    }
    
    //f = f1_01 * f2_11 * f3_10;
    double e = f1_01 + f2_11 + f3_10 - 3;
    double f = pow(1.5, e);
    
    return f;
}

template <typename EA>
double outline_split2(int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {

        double f1_01 = 1.0;
        double f2_11 = 1.0;
        double f3_10 = 1.0;
        double f4_00 = 1.0;
        
        // Compute fitness.
        for (int xy = 0; xy<grid_size; xy++) {
            
            // set the input states...
            int agent_x = floor(xy / max_x);
            int agent_y = xy % max_x;
            
            int p = agent_pos[xy];
            
            if (p == -1) {
                continue;
            }
            
            
            
            if (agent_x < (floor(max_x) / 2 )) {
                if ((agent_x ==0) || (agent_y == 0) || (agent_x == (max_x -1)) || (agent_y == (max_y -1))) {
                    if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                        ++f4_00;
                    }
                } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f3_10;
                }
                
            } else if (agent_x >= (floor(max_x) / 2 )) {
                if ((agent_x ==0) || (agent_y == 0) || (agent_x == (max_x -1)) || (agent_y == (max_y -1))) {
                    if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                        ++f1_01;
                    }
                } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f2_11;
                }
            }
            
            
            
            
            
        }
        
        
        //f = f1_01 * f2_11 * f3_10 * f4_00;
        double e = f1_01 + f2_11 + f3_10 + f4_00 - 4;
        double f = pow(1.5, e);
        
        return f;
}

template <typename EA>
double outline_quad(int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {

    
    double f1_01 = 0.0;
    double f2_11 = 0.0;
    double f3_10 = 0.0;
    double f4_00 = 0.0;
    double f5_11 = 0.0;
    double f6_10 = 0.0;
    
    // Compute fitness.
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        if (p == -1) {
            continue;
        }
        
        
        
        if (agent_x < (floor(max_x) / 2 )) {
            if ((agent_x ==0) || (agent_y == 0) || (agent_x == (max_x -1)) || (agent_y == (max_y -1))) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                    ++f4_00;
                }
            } else {
                if (agent_y < (floor(max_y) / 2 )) {
                    if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                        ++f3_10;
                    }
                } else {
                    if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                        ++f5_11;
                    }
                }
            }
            
        } else if (agent_x >= (floor(max_x) / 2 )) {
            if ((agent_x ==0) || (agent_y == 0) || (agent_x == (max_x -1)) || (agent_y == (max_y -1))) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f1_01;
                }
            } else {
                if (agent_y < (floor(max_y) / 2 )) {
                    
                    if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                        ++f2_11;
                    }
                } else {
                    if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                        ++f6_10;
                    }
                    
                }
                
            }
        }
        
        
        
        
        
    }
    
    
    double e = f1_01 + f2_11 + f3_10 + f4_00 + f5_11 + f6_10;
    double f = pow(1.5, e);
    
    return f;

}



/*! Sample fitness function for Markov networks.
 */
struct steppingstones : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    template <typename Individual, typename RNG, typename EA>
    double operator()(Individual& ind, RNG& rng, EA& ea) {
        
        int max_x = get<X_SIZE>(ea,10);
        int max_y = get<Y_SIZE>(ea,10);
        int grid_size = max_x * max_y;
        vector<int> agent_pos (grid_size, -1);
        vector<int> exec_order (grid_size);
        double f=0.0;
        
        
        // get the "prototype" phenotype (markov network):
        typename EA::phenotype_type &N = ealib::phenotype(ind, ea);
        vector<typename EA::phenotype_type> as; //
        
        
        for (int q=0; q<get<NUM_START_AGENTS>(ea,1); q++) {
            as.push_back(N); //grid_size, N); // my agents or networks
            agent_pos[q] = q;
            as[q].reset(rng.seed());
        }
        
        for (int i=0; i<grid_size; ++i) {
            exec_order[i] = i;
        }
        
        
        std::vector< std::vector<int> > cell_color(grid_size, std::vector<int>(2, 0));
        update_world_stigmergic_communication_N(get<WORLD_UPDATES>(ea,10), agent_pos, exec_order, as, cell_color, ind, rng, ea);
        
        accumulator_set<double, stats<tag::max> > fs;
        

        // given the final state, evaluate on a variety of different patterns
        fs(get<W1>(ea) * outline(grid_size, max_x, max_y, agent_pos, as, ea));
        fs(get<W2>(ea) * outline_split(grid_size, max_x, max_y, agent_pos, as, ea));
        fs(get<W3>(ea) * outline_split2(grid_size, max_x, max_y, agent_pos, as, ea));
        fs(get<W4>(ea) * outline_quad(grid_size, max_x, max_y, agent_pos, as, ea));
        
        f = boost::accumulators::max(fs);
        return (f);
        
        
        
 
    }
};

// Evolutionary algorithm definition.
typedef markov_network_evolution
< steppingstones
, recombination::asexual
, generational_models::moran_process< >
> ea_type;

/*! Define the EA's command-line interface.
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        add_mkv_options(this);
        
        add_option<POPULATION_SIZE>(this);
        add_option<MORAN_REPLACEMENT_RATE_P>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        
        add_option<X_SIZE>(this);
        add_option<Y_SIZE>(this);
        add_option<BRAIN_UPDATES>(this);
        add_option<WORLD_UPDATES>(this);
        add_option<INPUT_BIT_ERROR_PROB>(this);
        add_option<CAPABILITIES_OFF>(this);
        add_option<AGENT_DEATH_PROB>(this);
        add_option<NUM_START_AGENTS>(this);
        add_option<W1>(this);
        add_option<W2>(this);
        add_option<W3>(this);
        add_option<W4>(this);
        
    }
    
    
    virtual void gather_tools() {
        add_tool<analysis::dominant_genetic_graph>(this);
        add_tool<analysis::dominant_causal_graph>(this);
        add_tool<analysis::dominant_reduced_graph>(this);
        
        add_tool<ealib::analysis::markov_movie>(this);
        
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);



