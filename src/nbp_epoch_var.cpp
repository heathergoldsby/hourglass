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
#include <ea/analysis/archive.h>
#include <ea/events.h>


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
LIBEA_MD_DECL(TRANSUPDATES, "ea.hourglass.epv.transupdates", int); //
LIBEA_MD_DECL(MARK, "ea.hourglass.mark", int); //

template <typename EA>
struct mark_event : inheritance_event<EA> {
    //! Constructor.
    mark_event(EA& ea) : inheritance_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~mark_event() {
    }
    
    //! Called for every inheritance event.
    virtual void operator()(typename EA::population_type& parents,
                            typename EA::individual_type& offspring,
                            EA& ea) {
        for(typename EA::population_type::iterator i=parents.begin(); i!=parents.end(); ++i) {
            //offspring.traits().lod_parents().push_back(*i);
            if (get<MARK>(**i, 0) == 0) {
                get<MARK>(offspring,0) = get<MARK>(ea, 0);
                get<MARK>(**i, 0) = get<MARK>(ea, 0);
                int z =get<MARK>(**i, 0);
                int q = 0;
                
            } else {
                get<MARK>(offspring,0) = get<MARK>(**i, 0);
                int z =get<MARK>(offspring, 0);
                int q = 0;
            }

        }
    }
};

///*! An organism rotates to face its parent....
// */
//template <typename EA>
//struct mark_birth_event : birth_event<EA> {
//    
//    //! Constructor.
//    mark_birth_event(EA& ea) : birth_event<EA>(ea) {
//    }
//    
//    //! Destructor.
//    virtual ~mark_birth_event() {
//    }
//    
//    /*! Called for every inheritance event. We are using the orientation of the first parent...
//     */
//    virtual void operator()(typename EA::individual_type& offspring, // individual offspring
//                            typename EA::individual_type& parent, // individual parent
//                            EA& ea) {
//        //ea.env().face_org(parent, offspring);
//        //get<GERM_STATUS>(offspring, true) = get<GERM_STATUS>(ind(parents.begin(),ea), true);
//        get<MARK>(offspring) = get<MARK>(parent,0);
//        
//    }
//};
//



template <typename EA>
double eval_body_med (int ffToUse, int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    double f = 0.0;
    
    switch(ffToUse) {

        case 1: // case 2:
            f = body_plan_b(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 2: //case 6:
            f = body_plan_f(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 3: //case 7:
            f = body_plan_g(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 4: //case 8:
            f = body_plan_h(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 5: //case 13:
            f = body_plan_d1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 6: // case 14:
            f = body_plan_e1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 7: //case 17:
            f = body_plan_h1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 8: //case 23:
            f = body_plan_n1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 9: //case 24:
            f = body_plan_o1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 10: //case 25:
            f = body_plan_p1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 11: // case 31:
            f = body_plan_v1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 12: //case 35:
            f = body_plan2(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 13: //case 40:
            f = body_plan7(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 14: //case 45:
            f = body_plan12(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 15: //case 50:
            f = body_plan17(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
    }
    
    return f;
    
}


template <typename EA>
double eval_body_plan (int ffToUse, int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    double f = 0.0;
    
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

    return f;
   
}

struct epoch_variation_gradual : fitness_function<unary_fitness<double>, constantS, stochasticS> {
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
        int epoch_steps = cur_update - (epoch * (run_length/vect.size()));
        int transupdates = get<TRANSUPDATES>(ea,0);
        
        if ((epoch_steps < transupdates) && (epoch > 0)) {
            double w2 = double(epoch_steps) / double(transupdates);
            double w1 = 1 - (w2);
            
            int prevff = vect[epoch-1];
            double f1 = eval_body_med(prevff, grid_size, max_x, max_y, agent_pos, as, ea);
            double f2 = eval_body_med(ffToUse, grid_size, max_x, max_y, agent_pos, as, ea);
            
            f = w1 * f1 + w2 * f2;
            
        } else {
            f = eval_body_med(ffToUse, grid_size, max_x, max_y, agent_pos, as, ea);
        }
        
        
        return (f);
        
    }
};


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
        
        f = eval_body_med(ffToUse, grid_size, max_x, max_y, agent_pos, as, ea);
        
    
        return (f);
        
    }
};

// Evolutionary algorithm definition.
typedef markov_network_evolution
< epoch_variation_gradual
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
        add_option<TRANSUPDATES>(this);
        add_option<ARCHIVE_OUTPUT>(this);
        add_option<ARCHIVE_INPUT>(this);
        add_option<ANALYSIS_IND_NAME>(this);
        add_option<MARK>(this);
        
        
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
        add_tool<ealib::analysis::archive_dominant>(this);
        add_tool<ealib::analysis::copy_individual>(this);
        add_tool<ealib::analysis::merge_archives>(this);
        add_tool<ealib::analysis::recalculate_fitnesses>(this);
        add_tool<ealib::analysis::rename_individuals>(this);

        
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
        add_event<mark_event>(ea);
        
    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);



