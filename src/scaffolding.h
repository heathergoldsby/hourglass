//
//  scaffolding.h
//  hourglass
//
//  Created by Heather Goldsby on 1/11/16.
//  Copyright © 2016 Michigan State University. All rights reserved.
//

#ifndef scaffolding_h
#define scaffolding_h

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>

using namespace std;
using namespace boost::accumulators;



using namespace ealib;


template <typename EA>
double perfect_split(int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    double f=0.0;
    return f;

    // rescale by maximum possible fitness.
    

}

template <typename EA>
double three_stripes(int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    double f=0.0;
 
    
    // guts of 3 stripes
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
    
    // rescale by maximum possible fitness.
    
    return f;
}

template <typename EA>
double six_rectangles(int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    double f=0.0;
    return f;
}

template <typename EA>
double nine_rectangles(int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    double f=0.0;
    return f;
}

/*! Sample fitness function for Markov networks.
 */
struct scaffolding : fitness_function<unary_fitness<double>, constantS, stochasticS> {
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
        
        accumulator_set<double, stats<tag::max> > fs;
        

        // given the final state, evaluate on a variety of different patterns
        fs(perfect_split(grid_size, max_x, max_y, agent_pos, as, ea));
        fs(three_stripes(grid_size, max_x, max_y, agent_pos, as, ea));
        fs(six_rectangles(grid_size, max_x, max_y, agent_pos, as, ea));
        fs(nine_rectangles(grid_size, max_x, max_y, agent_pos, as, ea));
        
        f = boost::accumulators::max(fs);
        return (f);
        
        
 
    }
};

#endif /* scaffolding_h */
