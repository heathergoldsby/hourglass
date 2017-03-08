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
#include "new_body_plans.h"
#include "body_plans.h"


using namespace std;
using namespace boost::accumulators;
using namespace ealib;
using namespace mkv;

/*! Epoch variation -- change environment each epoch. 
 
 Each run is divided up into 10 segments. Each segment, we check for what fitness
 function to run. We then change which environment is used accordingly.
 */


LIBEA_MD_DECL(FITPEREPOCH, "ea.hourglass.epv.fitperepoch", string); //




struct epoch_variation : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    template <typename Individual, typename RNG, typename EA>
    double operator()(Individual& ind, RNG& rng, EA& ea) {
        
        std::vector<int> vect;
        
        std::stringstream ss(get<FITPEREPOCH>(ea));
        
        int i;
        
        while (ss >> i)
        {
            vect.push_back(i);
            
            if (ss.peek() == ',')
                ss.ignore();
        }
        
        
        
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
        int epoch = floor(cur_update / (run_length / vect.size()));
        int ffToUse = vect[epoch];
        
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
            case 10:
                f = body_plan_a1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 11:
                f = body_plan_b1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 12:
                f = body_plan_c1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 13:
                f = body_plan_d1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 14:
                f = body_plan_e1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 15:
                f = body_plan_f1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 16:
                f = body_plan_g1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 17:
                f = body_plan_h1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 18:
                f = body_plan_i1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 19:
                f = body_plan_j1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 20:
                f = body_plan_k1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 21:
                f = body_plan_l1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 22:
                f = body_plan_m1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 23:
                f = body_plan_n1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 24:
                f = body_plan_o1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 25:
                f = body_plan_p1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 26:
                f = body_plan_q1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 27:
                f = body_plan_r1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 28:
                f = body_plan_s1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 29:
                f = body_plan_t1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 30:
                f = body_plan_u1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 31:
                f = body_plan_v1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 32:
                f = body_plan_w1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 33:
                f = body_plan0(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 34:
                f = body_plan1(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 35:
                f = body_plan2(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 36:
                f = body_plan3(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 37:
                f = body_plan4(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 38:
                f = body_plan5(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 39:
                f = body_plan6(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 40:
                f = body_plan7(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 41:
                f = body_plan8(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 42:
                f = body_plan9(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 43:
                f = body_plan10(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 44:
                f = body_plan11(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 45:
                f = body_plan12(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 46:
                f = body_plan13(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 47:
                f = body_plan14(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 48:
                f = body_plan15(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 49:
                f = body_plan16(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 50:
                f = body_plan17(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 51:
                f = body_plan18(grid_size, max_x, max_y, agent_pos, as, ea);
                break;
            case 52:
                f = body_plan19(grid_size, max_x, max_y, agent_pos, as, ea);
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

        add_option<FITPEREPOCH>(this);

        
        
        
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



