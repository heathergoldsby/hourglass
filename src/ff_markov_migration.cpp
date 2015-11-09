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
        
        // start with one agent...
        vector<typename EA::phenotype_type> as; //
        as.push_back(N); //grid_size, N); // my agents or networks
        agent_pos[0] = 0;
        as[0].reset(rng.seed());
        
        
        
        
        for (int i=0; i<grid_size; ++i) {
            exec_order[i] = i;
        }
        
        
        // World update... this is where growth may occur.
        // Run agents for X updates
        
        
        // Inputs:
        // (0) north color, (1) north color,
        // (2) east color, (3) east color,
        // (4) south color, (5) south color,
        // (6) west color, (7) west color
        
        // (8) and (9) x
        // (10) and (11) y
        
        
        // Outputs:
        // (0) color, (1) color
        // (2) direction for repro, (3) direction for repro
        // (4) move
        // (5) reproduce
        
        // for this grid, 0,0 is upper left.
        // for this grid, 0,0 is upper left.
        int world_updates = get<WORLD_UPDATES>(ea,10);
        int brain_updates = get<BRAIN_UPDATES>(ea,10);
        for(int t=0;t<world_updates;t++){
            
            // Must randomize agent execution order...
            std::random_shuffle ( exec_order.begin(), exec_order.end() );
            
            // Brain update... there should be about 4-8 brain updates per agent per world update
            for (int j=0; j<exec_order.size(); ++j) {
                // get x,y coord of agent
                int xy = exec_order[j];
                // where we can find the agent itself
                int p = agent_pos[exec_order[j]];
                
                // no agent exists
                if (p == -1) {
                    continue;
                }
                
                
                // set the input states...
                int agent_y = floor(xy / max_x);
                int agent_x = xy % max_x;
                int north = (agent_y - 1) * max_x + agent_x;
                int east = agent_y * max_x + agent_x + 1;
                int south = (agent_y + 1) * max_x + agent_x;
                int west = agent_y * max_x + agent_x - 1;

                
                
                // migrate
                if (as[p].output(4)) {
                    
                    int d1 = as[p].output(2);
                    int d2 = as[p].output(3);
                    
                    if ((as[p].output(2) == 0) && (as[p].output(3)== 0) && (agent_y > 0)) { // 00 north
                        if (agent_pos[north] == -1){
                            agent_pos[north] = p; // move agent
                            agent_pos[xy] = -1;
                            xy = north;
                        }
                    } else if ((as[p].output(2) == 0) && (as[p].output(3)== 1) && (agent_x < (max_x - 1))) {  // 01 east
                        if (agent_pos[east] == -1) {
                            agent_pos[east] = p; // move agent
                            agent_pos[xy] = -1;
                            xy = east;
                        }
                    } else if ((as[p].output(2) == 1) && (as[p].output(3)== 1) && (agent_y < (max_y - 1))) { // 11 south
                        if (agent_pos[south] == -1){
                            agent_pos[south] = p; // move agent
                            agent_pos[xy] = -1;
                            xy = south;
                        }
                    } else if ((as[p].output(2) == 1) && (as[p].output(3)== 0) && (agent_x > 0)) { // 10 west
                        if (agent_pos[west] == -1){
                            agent_pos[west] = p; // move agent
                            agent_pos[xy] = -1;
                            xy = west;
                        }
                    }
                    
                    agent_y = floor(xy / max_x);
                    agent_x = xy % max_x;
                    north = (agent_y - 1) * max_x + agent_x;
                    east = agent_y * max_x + agent_x + 1;
                    south = (agent_y + 1) * max_x + agent_x;
                    west = agent_y * max_x + agent_x - 1;
                }

                
                
                
                // north neighbor
                if (agent_y > 0) {
                    if (agent_pos[north] != -1) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[north]];
                        (as[p]).input(0) = neighbor.output(0);
                        (as[p]).input(1) = neighbor.output(1);
                    }
                }
                
                // east neighbor
                if (agent_x < (max_x - 2)) {
                    if (agent_pos[east] != -1) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[east]];
                        (as[p]).input(2) = neighbor.output(0);
                        (as[p]).input(3) = neighbor.output(1);
                    }
                }
                
                // south neighbor
                if (agent_y < (max_x - 2)) {
                    if (agent_pos[south] != -1) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[south]];
                        (as[p]).input(4) = neighbor.output(0);
                        (as[p]).input(5) = neighbor.output(1);
                    }
                }
                
                // west neighbor
                if (agent_x > 0) {
                    if (agent_pos[west] != -1) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[west]];
                        (as[p]).input(6) = neighbor.output(0);
                        (as[p]).input(7) = neighbor.output(1);
                    }
                }
                
                
                // origin
                if (agent_x == 0 and agent_y == 0) {
                    as[p].input(8) = 1;
                } else {
                    as[p].input(8) = 0;
                }
                
                
                // edge
                if ((agent_x == 0) or (agent_y == 0) or (agent_x == (max_x -1)) or (agent_y == (max_y -1))) {
                    as[p].input(9) = 1;
                } else {
                    as[p].input(9) = 0;
                }
                
                // Give them their coordinates...
                
                int bsize = 10;
                vector<bool> xcoor(bsize);
                vector<bool> ycoor(bsize);
                
                ealib::algorithm::int2range(agent_x, xcoor.begin());
                ealib::algorithm::int2range(agent_y, ycoor.begin());
                
                int cur_input = 10;
                for (int i = 0; i < bsize; ++i) {
                    as[p].input(cur_input) = xcoor[i];
                    as[p].input(cur_input + bsize) = ycoor[i];
                    ++cur_input;
                }
                
                
                
                // reproduce
                if (as[p].output(5)) {
                    
                    int d1 = as[p].output(2);
                    int d2 = as[p].output(3);
                    
                    if ((as[p].output(2) == 0) && (as[p].output(3)== 0) && (agent_y > 0)) { // 00 north
                        if (agent_pos[north] == -1){
                            as.push_back(N); // Add a new agent.
                            agent_pos[north] = (as.size() -1); // This agent is at the end...
                        }
                    } else if ((as[p].output(2) == 0) && (as[p].output(3)== 1) && (agent_x < (max_x - 1))) {  // 01 east
                        if (agent_pos[east] == -1) {
                            as.push_back(N); // Add a new agent.
                            agent_pos[east] = (as.size() -1); // This agent is at the end...
                        }
                    } else if ((as[p].output(2) == 1) && (as[p].output(3)== 1) && (agent_y < (max_y - 1))) { // 11 south
                        if (agent_pos[south] == -1){
                            as.push_back(N); // Add a new agent.
                            agent_pos[south] = (as.size() -1); // This agent is at the end...
                        }
                    } else if ((as[p].output(2) == 1) && (as[p].output(3)== 0) && (agent_x > 0)) { // 10 west
                        if (agent_pos[west] == -1){
                            as.push_back(N); // Add a new agent.
                            agent_pos[west] = (as.size() -1); // This agent is at the end...
                        }
                    }
                }
                
                
                // apply errors if any.
                float error_prob = get<INPUT_BIT_ERROR_PROB>(ea,0);
                assert(0 <= error_prob <= 1.0);
                if (error_prob > 0) {
                    for (int k = 0; k < get<MKV_INPUT_N>(ea,0); k++) {
                        if (ea.rng().p(error_prob)) {
                            as[p].input(k) = ea.rng().bit();
                        }
                    }
                    
                }
                
                
                // update brain_updates times.
                for (int i = 0; i<brain_updates; ++i) {
                    (as[p]).update();
                }
                
            }
        }
        
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

