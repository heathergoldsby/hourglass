#ifndef _HOURGLASS_MARKOV_MOVIE_H_
#define _HOURGLASS_MARKOV_MOVIE_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>

#include <ea/datafile.h>
#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include "hourglass.h"
using namespace std;
using namespace boost::accumulators;


namespace ealib {
    namespace analysis {


        LIBEA_ANALYSIS_TOOL(recalc_fit) {
            accumulator_set<double, stats<tag::max, tag::mean> > fs;

            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                recalculate_fitness(*i, ea);
                fs(static_cast<int>(ealib::fitness(*i,ea)));
                
            }
            
            datafile df("recalc_fit.dat");
            df.write("pop size");
            df.write("mean");
            df.write("max");
            df.endl();
            df.write(ea.size());
            df.write(boost::accumulators::mean(fs));
            df.write(boost::accumulators::max(fs));
            df.endl();

        }
        
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
                    
                    int p = agent_pos[xy];
                    
                    if (p == -1) {
                        df.write("-1");
                        continue;
                    }

                    
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
    
        
        LIBEA_ANALYSIS_TOOL(markov_movie) {
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
            
            
            std::vector< std::vector<int> > cell_color(grid_size, std::vector<int>(2, 0));
            int world_updates = get<WORLD_UPDATES>(ea,10);
            for(int t=0;t<world_updates;t++){
                update_world_stigmergic_communication_N(1, agent_pos, exec_order, as, cell_color, best, ea.rng(), ea);
                
                
                df.write(t);
                // output time point for movie...
                for (int xy = 0; xy<grid_size; xy++) {
                    
                    int p = agent_pos[xy];
                    
                    if (p == -1) {
                        df.write("-1");
                        continue;
                    }
                    
                    
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
