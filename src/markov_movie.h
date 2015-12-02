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


        
        LIBEA_ANALYSIS_TOOL(movie_markov_growth_migration) {
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
            for (int q=0; q<get<NUM_START_AGENTS>(ea,1); q++) {
                as.push_back(N); //grid_size, N); // my agents or networks
                agent_pos[q] = q;
                as[q].reset(ea.rng().seed());
            }
            
            datafile df("movie.dat");
            df.write(max_x);
            df.write(max_y);
            df.endl();
            
            
            
            for (int i=0; i<grid_size; ++i) {
                exec_order[i] = i;
            }
            
            
            
            int world_updates = get<WORLD_UPDATES>(ea,10);
            for(int t=0;t<world_updates;t++){
                update_world_N(1, agent_pos, exec_order, as, best, ea.rng(), ea);
                
                
                df.write(t);
                // output time point for movie...
                for (int xy = 0; xy<grid_size; xy++) {
                    
                    // set the input states...
                    int agent_y = floor(xy / max_x);
                    int agent_x = xy % max_x;
                    
                    
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
                df.endl();
                
            }
        }
    

    
    }
}
#endif
