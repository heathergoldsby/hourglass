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
#include "body_plans.h"

using namespace std;
using namespace boost::accumulators;
using namespace ealib;
using namespace mkv;





/*! Sample fitness function for Markov networks.
 */
struct body_plans : fitness_function<unary_fitness<double>, constantS, stochasticS> {
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
        std::vector<int> apop_count(grid_size,0);
        update_world_stigmergic_communication_apop_N(get<WORLD_UPDATES>(ea,10), agent_pos, exec_order, as, cell_color, apop_count, ind, rng, ea);
        
        
        int fit_func = get<BODYPLAN>(ea,0) ;
        switch(fit_func) {
            case 0:
                f = body_plan0(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 1:
                f = body_plan1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 2:
                f = body_plan2(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 3:
                f = body_plan3(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 4:
                f = body_plan4(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 5:
                f = body_plan5(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 6:
                f = body_plan6(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 7:
                f = body_plan7(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 8:
                f = body_plan8(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 9:
                f = body_plan9(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 10:
                f = body_plan10(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 11:
                f = body_plan11(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 12:
                f = body_plan12(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 13:
                f = body_plan13(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 14:
                f = body_plan14(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 15:
                f = body_plan15(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 16:
                f = body_plan16(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 17:
                f = body_plan17(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 18:
                f = body_plan18(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 19:
                f = body_plan19(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
        }

        if (f == 0) {
            f = 1; 
        }
        return (f);
        
        
        
        
    }
};

// Evolutionary algorithm definition.
typedef markov_network_evolution
< body_plans
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
        add_option<BODYPLAN>(this);
        add_option<APOP_THRESH>(this);

        
    }
    
    
    virtual void gather_tools() {
        add_tool<analysis::dominant_genetic_graph>(this);
        add_tool<analysis::dominant_causal_graph>(this);
        add_tool<analysis::dominant_reduced_graph>(this);
        
        add_tool<ealib::analysis::markov_movie>(this);
        add_tool<ealib::analysis::recalc_fit>(this);
        
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);



