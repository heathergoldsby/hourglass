//
//  hourglass.h
//  hourglass
//
//  Created by Heather Goldsby on 9/29/15.
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//

#ifndef _HOURGLASS_H_
#define _HOURGLASS_H_
#include <ea/cmdline_interface.h>
#include <boost/algorithm/string/predicate.hpp>


using namespace ealib;


LIBEA_MD_DECL(X_SIZE, "ea.hourglass.x_size", int); // x dimension
LIBEA_MD_DECL(Y_SIZE, "ea.hourglass.y_size", int); // y dimension
LIBEA_MD_DECL(BRAIN_UPDATES, "ea.hourglass.brain_updates", int); // number of brain updates per world update
LIBEA_MD_DECL(WORLD_UPDATES, "ea.hourglass.world_updates", int); // number of world updates per fitness eval
LIBEA_MD_DECL(INPUT_BIT_ERROR_PROB, "ea.hourglass.input_bit_error_prob", float); // probability of an input bit flipping (between 0 and 1).
LIBEA_MD_DECL(AGENT_DEATH_PROB, "ea.hourglass.agent_death_prob", float); // probability of an agent dying.
LIBEA_MD_DECL(CAPABILITIES_OFF, "ea.hourglass.capabilities_off", std::string); // which capabilities are we turning off? (Setting input bits to 0 or ignoring output.)
LIBEA_MD_DECL(NUM_START_AGENTS, "ea.hourglass.num_start_agents", int); // number of agents to start with


// Run the world...

template <typename Individual, typename RNG, typename EA>
void update_world_N(int n, std::vector<int>& agent_pos, std::vector<int>& exec_order, std::vector<typename EA::phenotype_type>& as, Individual& ind, RNG& rng, EA& ea) {
    
    int max_x = get<X_SIZE>(ea,10);
    int max_y = get<Y_SIZE>(ea,10);
    typename EA::phenotype_type &N = ealib::phenotype(ind, ea);

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
    int brain_updates = get<BRAIN_UPDATES>(ea,10);
    const std::string& capabilities_off = get<CAPABILITIES_OFF>(ea, "");

    bool origin = (!boost::algorithm::icontains(capabilities_off, "origin"));
    bool migrate = (!boost::algorithm::icontains(capabilities_off, "migrate"));
    bool reproduce = (!boost::algorithm::icontains(capabilities_off, "reproduce"));
    bool coordinate = (!boost::algorithm::icontains(capabilities_off, "coordinate"));
    bool edge = (!boost::algorithm::icontains(capabilities_off, "edge"));
    
    
    
    // agent_pos, as, exec_order
    for(int t=0;t<n;t++){
        
        // Must randomize agent execution order...
        std::random_shuffle ( exec_order.begin(), exec_order.end(), ea.rng() );
        
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
            
            // is agent killed?
            float death_prob = get<AGENT_DEATH_PROB>(ea,0.0);
            assert(0 <= death_prob <= 1.0);
            if ((death_prob > 0) && ea.rng().p(death_prob)) {
                agent_pos[xy] = -1;
                continue;
                // we don't actually delete the agent from as b/c this
                // causes ordering trouble.
            }
        
            
            // set the input states...
            int agent_y = floor(xy / max_x);
            int agent_x = xy % max_x;
            int north = (agent_y - 1) * max_x + agent_x;
            int east = agent_y * max_x + agent_x + 1;
            int south = (agent_y + 1) * max_x + agent_x;
            int west = agent_y * max_x + agent_x - 1;
            
            
            
            // migrate
            if (migrate && as[p].output(4)) {
                
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
            if (origin) {
                if (agent_x == 0 and agent_y == 0) {
                    as[p].input(8) = 1;
                } else {
                    as[p].input(8) = 0;
                }
            }
            
            
            // edge
            if (edge) {
                if ((agent_x == 0) or (agent_y == 0) or (agent_x == (max_x -1)) or (agent_y == (max_y -1))) {
                    as[p].input(9) = 1;
                } else {
                    as[p].input(9) = 0;
                }
            }
            
            // Give them their coordinates...
            if (coordinate) {

                int bsize = 10;
                std::vector<bool> xcoor(bsize);
                std::vector<bool> ycoor(bsize);
                
                ealib::algorithm::int2range(agent_x, xcoor.begin());
                ealib::algorithm::int2range(agent_y, ycoor.begin());
                
                int cur_input = 10;
                for (int i = 0; i < bsize; ++i) {
                    as[p].input(cur_input) = xcoor[i];
                    as[p].input(cur_input + bsize) = ycoor[i];
                    ++cur_input;
                }
            }
            
            
            // reproduce
            if (reproduce && as[p].output(5)) {
                
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
}





#endif
