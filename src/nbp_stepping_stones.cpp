//
//  stepping stones
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
#include <ea/metapopulation.h>
#include <ea/island_model.h>
#include <ea/datafiles/metapopulation_fitness.h>

#include "markov_movie.h"
#include "hourglass.h"
#include "new_body_plans2.h"
#include "stepping_stones.h"
#include "stepping_stone_fitness.h"

using namespace std;
using namespace boost::accumulators;
using namespace ealib;
using namespace mkv;

/*! Steppingstones...
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
        
        
        int offset = get<START_POS>(ea,0);
        for (int q=0; q<get<NUM_START_AGENTS>(ea,1); q++) {
            as.push_back(N); //grid_size, N); // my agents or networks
            int pos = q + offset;
            if (pos >=  grid_size) {
                pos -= grid_size;
            }
            agent_pos[pos] = q;
            as[q].reset(rng.seed());
        }
        
        for (int i=0; i<grid_size; ++i) {
            exec_order[i] = i;
        }
        
        
        std::vector< std::vector<int> > cell_color(grid_size, std::vector<int>(2, 0));
        update_world_stigmergic_communication_N(get<WORLD_UPDATES>(ea,10), 0, agent_pos, exec_order, as, cell_color, ind, rng, ea);
        
        
        accumulator_set<double, stats<tag::max> > fs;
        
        
        // given the final state, evaluate on a variety of different patterns
        double f1 = body_plan_a(grid_size, max_x, max_y, agent_pos, as, ea);
        put<F1>(f1, ind);
        fs(get<W1>(ea) * f1);
        
        double f2 = body_plan_b(grid_size, max_x, max_y, agent_pos, as, ea);
        put<F2>(f2, ind);
        fs(get<W2>(ea) * f2);

        double f3 = body_plan_c(grid_size, max_x, max_y, agent_pos, as, ea);
        put<F3>(f3, ind);
        fs(get<W3>(ea) * f3);
        
        double f4 = body_plan_d(grid_size, max_x, max_y, agent_pos, as, ea);
        put<F4>(f4, ind);
        fs(get<W4>(ea) * f4);
        
        double f5 = body_plan_e(grid_size, max_x, max_y, agent_pos, as, ea);
        put<F5>(f5, ind);
        fs(get<W5>(ea) * f5);
        
        double f6 = body_plan_f(grid_size, max_x, max_y, agent_pos, as, ea);
        put<F6>(f6, ind);
        fs(get<W6>(ea) * f6);
        
        double f7 = body_plan_g(grid_size, max_x, max_y, agent_pos, as, ea);
        put<F7>(f7, ind);
        fs(get<W7>(ea) * f7);
        
        double f8 = body_plan_h(grid_size, max_x, max_y, agent_pos, as, ea);
        put<F8>(f8, ind);
        fs(get<W8>(ea) * f8);
        
        double f9 = body_plan_i(grid_size, max_x, max_y, agent_pos, as, ea);
        put<F9>(f9, ind);
        fs(get<W9>(ea) * f9);
        
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
        add_option<HIDDEN_BIT_ERROR_PROB>(this);
        add_option<OUTPUT_BIT_ERROR_PROB>(this);
        
        add_option<CAPABILITIES_OFF>(this);
        add_option<AGENT_DEATH_PROB>(this);
        add_option<NUM_START_AGENTS>(this);
        add_option<BODYPLAN>(this);
        add_option<APOP_THRESH>(this);
        add_option<START_POS>(this);
        add_option<RAND_ORDER>(this);
        add_option<NO_REP_PERIOD>(this);
        add_option<GATE_EPSILON>(this);
        
        add_option<W1>(this);
        add_option<W2>(this);
        add_option<W3>(this);
        add_option<W4>(this);
        add_option<W5>(this);
        add_option<W6>(this);
        add_option<W7>(this);
        add_option<W8>(this);
        add_option<W9>(this);
        
        
        
    }
    
    
    virtual void gather_tools() {
        /*
        add_tool<analysis::dominant_genetic_graph>(this);
        add_tool<analysis::dominant_causal_graph>(this);
        add_tool<analysis::dominant_reduced_graph>(this);
        
        add_tool<ealib::analysis::markov_movie>(this);
        add_tool<ealib::analysis::markov_movie_ko>(this);
        
        add_tool<ealib::analysis::ko>(this);*/
        add_tool<ealib::analysis::ko_stepping_stones>(this);
   
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
        add_event<datafiles::stepping_stone_fitness_dat>(ea);

    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);



