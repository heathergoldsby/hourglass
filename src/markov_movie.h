#ifndef _HOURGLASS_MARKOV_MOVIE_H_
#define _HOURGLASS_MARKOV_MOVIE_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include "hourglass.h"
using namespace std;

namespace ealib {
    namespace analysis {

        LIBEA_ANALYSIS_TOOL(movie_markov) {
            double max_fit;
            typename EA::individual_type best;
            
            // recalc all fitness values
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                recalculate_fitness(*i, ea);
                double tmp_fit = static_cast<double>(ealib::fitness(*i,ea));
                if (tmp_fit > max_fit) {
                    best = *i;
                }
            }

            int max_x = get<X_SIZE>(ea,10);
            int max_y = get<Y_SIZE>(ea,10);
            int grid_size = max_x * max_y;
            
            // Must create a WHOLE bunch of markov networks here...
            
            // get the "prototype" phenotype (markov network):
            typename EA::phenotype_type &N = ealib::phenotype(best, ea);
            
            // and make a few copies:
            vector<typename EA::phenotype_type> as(grid_size, N); // my agents or networks
            
            vector<int> agent_pos (grid_size, -1);
            vector<int> exec_order (grid_size);
            double f=0.0;
            
            for (int i=0; i<grid_size; ++i) {
                exec_order[i] = i;
            }
            
            datafile df("movie.dat");
            df.write(max_x);
            df.write(max_y);
            df.endl();
            
            
            /* For right now, there isn't a growth process -- just a FULL grid */
            for(int i=0;i<grid_size;i++) {
                agent_pos[i] = i;
                // probably want to reset the RNG for the markov network:
                as[i].reset(ea.rng().seed());
            }
            
            
            // World update... this is where growth may occur.
            // Run agents for X updates
            
            
            
            // Inputs: (0) x, (1) y,
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
                    if (xy == -1) {
                        continue;
                    }
                    
                    // set the input states...
                    int agent_x = floor(xy / max_x);
                    int agent_y = xy % max_x;
                    
                    
                    // north neighbor
                    int north = (agent_y - 1) * max_x + agent_x;
                    
                    if (agent_y > 0) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[north]];
                        (as[p]).input(2) = neighbor.output(0);
                        (as[p]).input(3) = neighbor.output(1);
                    }
                    
                    
                    // east neighbor
                    int east = agent_y * max_x + agent_x + 1;
                    if (agent_x < (max_x - 2)) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[east]];
                        (as[p]).input(4) = neighbor.output(0);
                        (as[p]).input(5) = neighbor.output(1);
                    }
                    
                    // south neighbor
                    int south = (agent_y + 1) * max_x + agent_x;
                    if (agent_y < (max_x - 2)) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[south]];
                        (as[p]).input(6) = neighbor.output(0);
                        (as[p]).input(7) = neighbor.output(1);
                    }
                    
                    // west neighbor
                    int west = agent_y * max_x + agent_x - 1;
                    if (agent_x > 0) {
                        typename EA::phenotype_type& neighbor = as[agent_pos[west]];
                        (as[p]).input(8) = neighbor.output(0);
                        (as[p]).input(9) = neighbor.output(1);
                        
                    }
                    
                    // edge detection.
                    if ((agent_x == 0) || (agent_y == 0) || (agent_x = (max_x-1)) || (agent_y == (max_y-1))) {
                        as[p].input(10) = 1;
                    }
                    
                    // update brain_updates times.
                    for (int i = 0; i<brain_updates; ++i) {
                        if ((ea.rng().uniform_integer(0,max_x)) > agent_x) {
                            (as[p]).input(0) = 1;
                        }
                        
                        if ((ea.rng().uniform_integer(0,max_y)) > agent_y) {
                            (as[p]).input(1) = 1;
                        }
                        
                        (as[p]).update();
                    }
                    
                    
                }
                
                df.write(t);
                // output time point for movie...
                for (int x=0; x < max_x; ++x) {
                    for (int y=0; y < max_y; ++y) {
                        int xy = y * max_x + x; // the agent position...
                        int p = agent_pos[xy];
                        
                        if (p == -1) {
                            df.write("-1");
                            continue;
                        }
                        
                        int o_zero = as[p].output(0);
                        int o_one = as[p].output(1);
                        
                        if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)) {
                            df.write("0");
                        } else if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)) {
                            df.write("1");
                        } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)) {
                            df.write("2");
                        } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)) {
                            df.write("3");
                        }
                        
                    }
                }
                df.endl();

            
        }
        }
        
        
        LIBEA_ANALYSIS_TOOL(movie_markov_growth) {
            double max_fit = 0;
            typename EA::individual_type best;
            
            // recalc all fitness values
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                recalculate_fitness(*i, ea);
                double tmp_fit = static_cast<int>(ealib::fitness(*i,ea));
                if (tmp_fit > max_fit) {
                    best = *i;
                    max_fit = tmp_fit;
                }
            }
            
            // For the best create a movie...
            int max_x = get<X_SIZE>(ea,10);
            int max_y = get<Y_SIZE>(ea,10);
            int grid_size = max_x * max_y;
            vector<int> agent_pos (grid_size, -1);
            vector<int> exec_order (grid_size);
            
            
            // Must create a WHOLE bunch of markov networks here...
            
            // get the "prototype" phenotype (markov network):
            typename EA::phenotype_type &N = ealib::phenotype(best, ea);
            
            // start with one agent...
            vector<typename EA::phenotype_type> as; //
            as.push_back(N); //grid_size, N); // my agents or networks
            agent_pos[0] = 0;
            as[0].reset(ea.rng().seed());
            
            datafile df("movie.dat");
            df.write(max_x);
            df.write(max_y);
            df.endl();
            
            
            
            for (int i=0; i<grid_size; ++i) {
                exec_order[i] = i;
            }
            
            
            // World update... this is where growth may occur.
            // Run agents for X updates
            
            
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
                    
                    // x and y coord.
                    (as[p]).input(0) = 1;
                    (as[p]).input(1) = 1;
                    
                    
                    
                    // north neighbor
                    int north = (agent_y - 1) * max_x + agent_x;
                    if (agent_y > 0) {
                        if (agent_pos[north] != -1) {
                            typename EA::phenotype_type& neighbor = as[agent_pos[north]];
                            (as[p]).input(2) = neighbor.output(0);
                            (as[p]).input(3) = neighbor.output(1);
                        }
                    }
                    
                    // east neighbor
                    int east = agent_y * max_x + agent_x + 1;
                    if (agent_x < (max_x - 2)) {
                        if (agent_pos[east] != -1) {
                            typename EA::phenotype_type& neighbor = as[agent_pos[east]];
                            (as[p]).input(4) = neighbor.output(0);
                            (as[p]).input(5) = neighbor.output(1);
                        }
                    }
                    
                    // south neighbor
                    int south = (agent_y + 1) * max_x + agent_x;
                    if (agent_y < (max_x - 2)) {
                        if (agent_pos[south] != -1) {
                            typename EA::phenotype_type& neighbor = as[agent_pos[south]];
                            (as[p]).input(6) = neighbor.output(0);
                            (as[p]).input(7) = neighbor.output(1);
                        }
                    }
                    
                    // west neighbor
                    int west = agent_y * max_x + agent_x - 1;
                    if (agent_x > 0) {
                        if (agent_pos[west] != -1) {
                            typename EA::phenotype_type& neighbor = as[agent_pos[west]];
                            (as[p]).input(8) = neighbor.output(0);
                            (as[p]).input(9) = neighbor.output(1);
                        }
                    }
                    
                    // reproduce
                    if (as[p].output(5)) {
                        
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
                        if ((ea.rng().uniform_integer(0,max_x)) > agent_x) {
                            (as[p]).input(0) = 1;
                        }
                        
                        if ((ea.rng().uniform_integer(0,max_y)) > agent_y) {
                            (as[p]).input(1) = 1;
                        }
                        
                        (as[p]).update();
                    }
                    
                } // exec order
                
                
                df.write(t);
                // output time point for movie...
                for (int x=0; x < max_x; ++x) {
                    for (int y=0; y < max_y; ++y) {
                        int xy = y * max_x + x; // the agent position...
                        int p = agent_pos[xy];
                        
                        if (p == -1) {
                            df.write("-1");
                            continue;
                        }
                        
                        int o_zero = as[p].output(0);
                        int o_one = as[p].output(1);
                        
                        if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)) {
                            df.write("0");
                        } else if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)) {
                            df.write("1");
                        } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)) {
                            df.write("2");
                        } else if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)) {
                            df.write("3");
                        }
                    
                    }
                }
                df.endl();

            } // world updates

        }
        
    }
}
#endif
