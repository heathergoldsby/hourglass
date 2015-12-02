/* markov_french_flag.cpp
 *
 * This file is part of EALib.
 *
 * Copyright 2014 David B. Knoester.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


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
using namespace ealib;
using namespace mkv;
using namespace std;
using namespace boost::accumulators;



/*! Sample fitness function for Markov networks.
 */
struct french_flag_fitness : fitness_function<unary_fitness<double>, constantS, stochasticS> {
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
        
        
        update_world_N(get<WORLD_UPDATES>(ea,10), agent_pos, exec_order, as, ind, rng, ea);
        
        double f1_10 = 1.0;
        double f1_01 = 1.0;
        double f1_11 = 1.0;
        
        double f2_10 = 1.0;
        double f2_01 = 1.0;
        double f2_11 = 1.0;
        
        double f3_10 = 1.0;
        double f3_01 = 1.0;
        double f3_11 = 1.0;
        
        
        // Compute fitness.
        for (int xy = 0; xy<grid_size; xy++) {
            
            // set the input states...
            int agent_x = floor(xy / max_x);
            int agent_y = xy % max_x;
            
            int p = agent_pos[xy];
            
            if (p == -1) {
                continue;
            }
            
            if (agent_x < (floor(max_x) / 3 )) {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f1_10;
                } else if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f1_01;
                } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f1_11;
                }
                
            } else if ((agent_x >= (floor(max_x) / 3))  && (agent_x < ((floor(max_x) / 3 * 2)))) {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f2_10;
                } else if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f2_01;
                } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f2_11;
                }
            } else if (agent_x >= ((floor(max_x) / 3 * 2))) {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f3_10;
                } else if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f3_01;
                } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f3_11;
                }
            }
            
        }
        
        accumulator_set<double, stats<tag::max> > fs;
        
        double fx = f1_10 * f2_01 * f3_11;
        fs(fx);
        
        fx = f1_10 * f2_11 * f3_01;
        fs(fx);
        
        fx = f1_01 * f2_10 * f3_11;
        fs(fx);
        
        fx = f1_01 * f2_11 * f3_10;
        fs(fx);
        
        fx = f1_11 * f2_01 * f3_10;
        fs(fx);
        
        fx = f1_11 * f2_10 * f3_01;
        fs(fx);
        
        f = boost::accumulators::max(fs);
        return f;
    }
};

// Evolutionary algorithm definition.
typedef markov_network_evolution
< french_flag_fitness
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
    }
    
    
    virtual void gather_tools() {
        add_tool<analysis::dominant_genetic_graph>(this);
        add_tool<analysis::dominant_causal_graph>(this);
        add_tool<analysis::dominant_reduced_graph>(this);
        
        add_tool<ealib::analysis::movie_markov_growth_migration>(this);
        
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);

