/* markov_square.cpp
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


/*! Sample fitness function for Markov networks.
 */
struct square_fitness : fitness_function<unary_fitness<double>, constantS, stochasticS> {
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
                int agent_x = floor(xy / max_x);
                int agent_y = xy % max_x;
            
                
                
                // north neighbor
                int north = (agent_y - 1) * max_x + agent_x;
                if (agent_y > 0) {
                    if (agent_pos[north] != -1) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[north]];
                        (as[p]).input(0) = neighbor.output(0);
                        (as[p]).input(1) = neighbor.output(1);
                    }
                }
                
                // east neighbor
                int east = agent_y * max_x + agent_x + 1;
                if (agent_x < (max_x - 2)) {
                    if (agent_pos[east] != -1) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[east]];
                        (as[p]).input(2) = neighbor.output(0);
                        (as[p]).input(3) = neighbor.output(1);
                    }
                }
                
                // south neighbor
                int south = (agent_y + 1) * max_x + agent_x;
                if (agent_y < (max_x - 2)) {
                    if (agent_pos[south] != -1) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[south]];
                        (as[p]).input(4) = neighbor.output(0);
                        (as[p]).input(5) = neighbor.output(1);
                    }
                }
                
                // west neighbor
                int west = agent_y * max_x + agent_x - 1;
                if (agent_x > 0) {
                    if (agent_pos[west] != -1) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[west]];
                        (as[p]).input(6) = neighbor.output(0);
                        (as[p]).input(7) = neighbor.output(1);
                    }
                }
                
                
                // Give them the solution... check that it works
                
                // For a minute, give them the solution...
                if (agent_x == 0 || agent_x == (max_x-1) || agent_y == 0 || agent_y == (max_y-1)) {
                    as[p].input(8) = 1;
                    as[p].input(9) = 0;
                } else if (agent_x == 1 || agent_x == (max_x-2) || agent_y == 1 || agent_y == (max_y-2)) {
                    as[p].input(8) = 1;
                    as[p].input(9) = 1;
                } else if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)) {
                    as[p].input(8) = 0;
                    as[p].input(9) = 1;
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
                
                
                // update brain_updates times.
                for (int i = 0; i<brain_updates; ++i) {
//                    if ((ea.rng().uniform_integer(0,max_x)) > agent_x) {
//                        (as[p]).input(8) = 1;
//                    } else {
//                        (as[p]).input(8) = 0;
//                    }
//                    
//                    if ((ea.rng().uniform_integer(0,max_y)) > agent_y) {
//                        (as[p]).input(10) = 1;
//                    } else {
//                        (as[p]).input(10) = 0;
//                    }
                    
                    (as[p]).update();
                }
                
                
            }
            
        }
        
        double f1 = 1.0;
        double f2 = 1.0;
        double f3 = 1.0;
        
        // Compute fitness.
        for (int xy = 0; xy<grid_size; xy++) {
            
            // set the input states...
            int agent_x = floor(xy / max_x);
            int agent_y = xy % max_x;
            
            int p = agent_pos[xy];
            
            if (p == -1) {
                continue;
            }
            
            // wrong ff... maybe? strange nested issues...
            if (agent_x == 0 || agent_x == (max_x-1) || agent_y == 0 || agent_y == (max_y-1)) {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f1;
                }
            } else if (agent_x == 1 || agent_x == (max_x-2) || agent_y == 1 || agent_y == (max_y-2)) {
                if  (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)) {
                    ++f2;
                }
            } else if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)) {
                ++f3;
            }

            
        }

        
        // and return some measure of fitness:
        // ponder gamma transform
//        double fit_max = grid_size;
//        double fit_min = 0;
        
//        double rescaled_fit = 100*pow((((f-fit_min) / (fit_max - fit_min))), (get<FIT_GAMMA>(ea))) + 1;
        
        f = f1 * f2 * f3;

        
        
        return f;
    }
};

// Evolutionary algorithm definition.
typedef markov_network_evolution
< square_fitness
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
        add_option<FIT_GAMMA>(this);
        
        
    }
    
    
    virtual void gather_tools() {
        add_tool<analysis::dominant_genetic_graph>(this);
        add_tool<analysis::dominant_causal_graph>(this);
        add_tool<analysis::dominant_reduced_graph>(this);
        
        add_tool<ealib::analysis::movie_markov_growth>(this);
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);
