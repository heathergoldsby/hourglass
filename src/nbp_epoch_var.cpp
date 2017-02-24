//
//  epoch variation
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
#include <ea/datafiles/metapopulation_fitness.h>

#include "markov_movie.h"
#include "hourglass.h"
#include "new_body_plans2.h"


using namespace std;
using namespace boost::accumulators;
using namespace ealib;
using namespace mkv;

/*! Epoch variation -- change environment each epoch. 
 
 Each run is divided up into 10 segments. Each segment, we check for what fitness
 function to run. We then change which environment is used accordingly.
 */


LIBEA_MD_DECL(FIT0, "ea.hourglass.epv.fit0", int); //
LIBEA_MD_DECL(FIT1, "ea.hourglass.epv.fit1", int); //
LIBEA_MD_DECL(FIT2, "ea.hourglass.epv.fit2", int); //
LIBEA_MD_DECL(FIT3, "ea.hourglass.epv.fit3", int); //
LIBEA_MD_DECL(FIT4, "ea.hourglass.epv.fit4", int); //
LIBEA_MD_DECL(FIT5, "ea.hourglass.epv.fit5", int); //
LIBEA_MD_DECL(FIT6, "ea.hourglass.epv.fit6", int); //
LIBEA_MD_DECL(FIT7, "ea.hourglass.epv.fit7", int); //
LIBEA_MD_DECL(FIT8, "ea.hourglass.epv.fit8", int); //
LIBEA_MD_DECL(FIT9, "ea.hourglass.epv.fit9", int); //




struct epoch_variation : fitness_function<unary_fitness<double>, constantS, stochasticS> {
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
        
        // what epoch are we in?
        int run_length = get<RUN_UPDATES>(ea);
        int cur_update = ea.current_update();
        int epoch = floor(cur_update / (run_length / 10));
        int ffToUse;
        
        switch(epoch) {
            case 0:
                ffToUse = get<FIT0>(ea);
                break;
            case 1:
                ffToUse = get<FIT1>(ea);
                break;
            case 2:
                ffToUse = get<FIT2>(ea);
                break;
            case 3:
                ffToUse = get<FIT3>(ea);
                break;
            case 4:
                ffToUse = get<FIT4>(ea);
                break;
            case 5:
                ffToUse = get<FIT5>(ea);
                break;
            case 6:
                ffToUse = get<FIT6>(ea);
                break;
            case 7:
                ffToUse = get<FIT7>(ea);
                break;
            case 8:
                ffToUse = get<FIT8>(ea);
                break;
            case 9:
                ffToUse = get<FIT9>(ea);
                break;
                
        }
        
        switch(ffToUse) {
            case 1:
                f = body_plan_a(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 2:
                f = body_plan_b(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 3:
                f = body_plan_c(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 4:
                f = body_plan_d(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 5:
                f = body_plan_e(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 6:
                f = body_plan_f(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 7:
                f = body_plan_g(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 8:
                f = body_plan_h(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 9:
                f = body_plan_i(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
                
        }
        
    
        return (f);
        
    }
};

// Evolutionary algorithm definition.
typedef markov_network_evolution
< epoch_variation
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

        add_option<FIT0>(this);
        add_option<FIT1>(this);
        add_option<FIT2>(this);
        add_option<FIT3>(this);
        add_option<FIT4>(this);
        add_option<FIT5>(this);
        add_option<FIT6>(this);
        add_option<FIT7>(this);
        add_option<FIT8>(this);
        add_option<FIT9>(this);
        
        
        
    }
    
    
    virtual void gather_tools() {
        /*
         add_tool<analysis::dominant_genetic_graph>(this);
         add_tool<analysis::dominant_causal_graph>(this);
         add_tool<analysis::dominant_reduced_graph>(this);
         
         add_tool<ealib::analysis::markov_movie>(this);
         add_tool<ealib::analysis::markov_movie_ko>(this);
         
         add_tool<ealib::analysis::ko>(this);
        add_tool<ealib::analysis::ko_stepping_stones>(this);*/
        
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
        
    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);



