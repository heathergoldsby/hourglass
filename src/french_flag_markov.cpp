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
#include <ea/mkv/markov_network_evolution.h>
#include <ea/generational_models/moran_process.h>
#include <ea/fitness_function.h>
#include <ea/cmdline_interface.h>
#include <ea/datafiles/fitness.h>
using namespace ealib;
using namespace mkv;
using namespace std;


/*! Sample fitness function for Markov networks.
 */
struct french_flag_fitness : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    template <typename Individual, typename RNG, typename EA>
    double operator()(Individual& ind, RNG& rng, EA& ea) {
        
        int max_x = 6; // change to config params later...
        int max_y = 6;
        int grid_size = max_x * max_y;
        
        // Must create a WHOLE bunch of markov networks here...
        
        // get the "prototype" phenotype (markov network):
        typename EA::phenotype_type &N = ealib::phenotype(ind, ea);
        
        // and make a few copies:
        vector<typename EA::phenotype_type> as(grid_size, N); // my agents or networks
        
        
        
        
        vector<int> agent_pos (grid_size, -1);
        vector<int> exec_order (grid_size);
        double f=0.0;
        
        for (int i=0; i<grid_size; ++i) {
            exec_order[i] = i;
        }
        
        
        /* For right now, there isn't a growth process -- just a FULL grid */
        for(int i=0;i<grid_size;i++) {
            agent_pos[i] = i;
            // probably want to reset the RNG for the markov network:
            as[i].reset(rng.seed());
        }
        
        
        // World update... this is where growth may occur.
        // Run agents for X updates
        
        // setup inputs...
        std::vector<int> inputs (10,0);
        
        // Inputs: (0) x, (1) y, // currently 1,1... they may not need them?
        // (2) north color, (3) north color,
        // (4) east color, (5) east color,
        // (6) south color, (7) south color,
        // (8) west color, (9) west color
        
        // Outputs:
        // (0) color, (1) color
        // (2) direction for repro, (3) direction for repro
        // (4) move
        // (5) reproduce
        
        // for this grid, 0,0 is upper left.
        for(int t=0;t<100;t++){
            
            // Must randomize agent execution order...
            std::random_shuffle ( exec_order.begin(), exec_order.end() );
            
            // Brain update... there should be about 4-8 brain updates per agent per world update
            for (int j=0; j<exec_order.size(); ++j) {
                // get x,y coord of agent
                int xy = exec_order[j];
                // where we can find the agent itself
                int p = agent_pos[exec_order[j]];
                
                // no agent exists
                if (xy == -1) {
                    continue;
                }
                
                // set the input states...
                int agent_x = floor(xy / max_x);
                int agent_y = xy % max_x;
                typename EA::phenotype_type me = as[p];
                
                (as[p]).clear();
                
                (as[p]).input(0) = 1;
                (as[p]).input(1) = 1;
                inputs[0];
                inputs[1];
                typename EA::phenotype_type neighbor;
                
                
                // north neighbor
                int north = (agent_y - 1) * max_x + agent_x;
                
                if (agent_y > 0) {
                    neighbor = as[agent_pos[north]];
                    (as[p]).input(2) = neighbor.output(0);
                    (as[p]).input(3) = neighbor.output(1);
                }
                
                
                // east neighbor
                int east = agent_y * max_x + agent_x + 1;
                if (agent_x < (max_x - 2)) {
                    neighbor = as[agent_pos[east]];
                    (as[p]).input(4) = neighbor.output(0);
                    (as[p]).input(5) = neighbor.output(1);
                }
                
                // south neighbor
                int south = (agent_y + 1) * max_x + agent_x;
                
                if (agent_y < (max_x - 2)) {
                    neighbor = as[agent_pos[south]];
                    (as[p]).input(6) = neighbor.output(0);
                    (as[p]).input(7) = neighbor.output(1);
                }
                
                // west neighbor
                int west = agent_y * max_x + agent_x - 1;
                
                if (agent_x > 0) {
                    neighbor = as[agent_pos[west]];
                    (as[p]).input(8) = neighbor.output(0);
                    (as[p]).input(9) = neighbor.output(1);
                    
                }
                
                // update 4 times.
                for (int i = 0; i<4; ++i) {
                    (as[p]).update();
                }
                
                
            }
            
            
            
            
            
        }
        
        
        // Compute fitness.
        for (int x=0; x < max_x; ++x) {
            for (int y=0; y < max_y; ++y) {
                int xy = y * max_x + x; // the agent position...
                int p = agent_pos[xy];
                
                if (p == -1) {
                    continue;
                }
                
                if ((x < (floor(max_x) / 3 )) && (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0))){
                    // blue 10
                    ++f;
                } else if (((x > (floor(max_x) / 3))  && (x < floor(max_x) / 3 * 2))  &&
                           (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1))){
                    // white 11
                    ++f;
                } else if ((((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)) ) {
                    // red 01
                    ++f;
                }
            }
        }
        
        
        as.clear();
        
        
        // and return some measure of fitness:
        // ponder gamma transform
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
    }
    
    
    virtual void gather_tools() {
        add_tool<analysis::dominant_genetic_graph>(this);
        add_tool<analysis::dominant_causal_graph>(this);
        add_tool<analysis::dominant_reduced_graph>(this);
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);
