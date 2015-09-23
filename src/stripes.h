//
//  stripes.h
//  ealife
//
//  Created by Heather Goldsby on 11/22/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//


#ifndef _EALIFE_STRIPES_H_
#define _EALIFE_STRIPES_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/instruction_set.h>
#include <ea/digital_evolution/environment.h>
#include <ea/digital_evolution/utils/resource_consumption.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <ea/subpopulation_founder.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/metapopulation.h>
#include <ea/selection/random.h>
#include <ea/selection/proportionate.h>
#include <ea/selection/tournament.h>
#include <ea/mutation.h>
#include <ea/recombination.h>

using namespace ealib;


LIBEA_MD_DECL(PATTERN_FIT, "ea.pattern.fit", int); // count the number of organisms that have the right pattern
LIBEA_MD_DECL(NUM_PROPAGULE_CELL, "ea.stripes.num_propagule_cell", int);

LIBEA_MD_DECL(MULTICELL_REP_TIME, "ea.stripes.mc_rep_time", int);



//! Stripe fitness.
struct permute_stripes : public fitness_function<unary_fitness<double>, nonstationaryS> {
    template <typename EA>
    int eval_permute_stripes(EA& ea) {
        // vert stripes
        double one_fit_not = 0;
        double one_fit_nand = 0;
        double two_fit_not = 0;
        double two_fit_nand = 0;
        // horizontal stripes
        double three_fit_not = 0;
        double three_fit_nand = 0;
        double four_fit_not = 0;
        double four_fit_nand = 0;
        // diagonal stripes
        double five_fit_not = 0;
        double five_fit_nand = 0;
        double six_fit_not = 0;
        double six_fit_nand = 0;
        int num_org = 0;
        
        for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
            for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                typename EA::location_type& l = ea.env().location(x,y);
                if (!l.occupied()) {
                    continue;
                }
                
                std::string lt = get<LAST_TASK>(*l.inhabitant(),"");
                // Vertical stripes!
                if((y % 2) == 0) {
                    if (lt == "nand") { ++one_fit_nand; }
                    if (lt == "not") { ++two_fit_not; }
                } else {
                    if(lt == "not") { ++one_fit_not; }
                    if (lt == "nand") { ++two_fit_nand; }
                }
                
                // Horizontal stripes
                if ((x % 2) == 0) {
                    if (lt == "nand") { ++three_fit_nand; }
                    if (lt == "not") { ++four_fit_not; }
                } else {
                    if(lt == "not") { ++three_fit_not; }
                    if (lt == "nand") { ++four_fit_nand; }
                }
                
                // Diagonal stripes
                if ((x % 2) == (y % 2)) {
                    if(lt == "not") { ++five_fit_not; }
                    if (lt == "nand") { ++six_fit_nand; }
                } else {
                    if(lt == "nand") { ++five_fit_nand; }
                    if (lt == "not") { ++six_fit_not; }
                }
            }
        }
        
        double tmp_one_fit = (one_fit_not + 1)  * (one_fit_nand + 1);
        double tmp_two_fit = (two_fit_not + 1)  * (two_fit_nand + 1);
        double tmp_three_fit = (three_fit_not + 1)  * (three_fit_nand + 1);
        double tmp_four_fit = (four_fit_not + 1)  * (four_fit_nand + 1);
        double tmp_five_fit = (five_fit_not + 1)  * (five_fit_nand + 1);
        double tmp_six_fit = (six_fit_not + 1)  * (six_fit_nand + 1);
        double tmp_fit = std::max(tmp_one_fit, tmp_two_fit);
        tmp_fit = std::max(tmp_fit, tmp_three_fit);
        tmp_fit = std::max(tmp_fit, tmp_four_fit);
        tmp_fit = std::max(tmp_fit, tmp_five_fit);
        tmp_fit = std::max(tmp_fit, tmp_six_fit);
        
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_permute_stripes(sea));
        put<PATTERN_FIT>(f,sea);
        return f;
    }
};





template <typename EA>
std::string get_last_task(int x, int y, EA& ea) {
    typename EA::environment_type::location_type* l = &ea.env().location(x,y);
    std::string lt = "";
    if (l->occupied()) {
        lt = get<LAST_TASK>(*(l->inhabitant()),"");
    }
    return lt;
}



// Can we do this without repeating code?
template <typename EA>
void eval_square(EA& ea) {
    double tmp_fit = 0;
    
    int num_org = 0;
    int max_x = get<SPATIAL_X>(ea) - 1; // minus one because SPATIAL_X is the size
    int max_y = get<SPATIAL_Y>(ea) - 1;
    
    for (int x=0; x < max_x; ++x) {
        for (int y=0; y<max_y; ++y){
            typename EA::location_type& l = ea.env().location(x,y);
            if (!l.occupied()) {
                continue;
            }
            
            std::string lt = get<LAST_TASK>(*l.inhabitant(),"");
            if (x == 0 || x == max_x || y == 0 || y == max_y) {
                if (lt == "not") {
                    ++tmp_fit;
                }
            } else if (x == 1 || x == (max_x - 1) || y == 1 || y == (max_y -1)) {
                if (lt == "nand") {
                    ++tmp_fit;
                }
            } else if (lt == "ornot") {
                ++tmp_fit;
            }
            
            
            
        }
    }
    put<PATTERN_FIT>(tmp_fit, ea);

 
    
}


//! Square fitness - not (outside), nand (next), ornot (inside)
struct square : public fitness_function<unary_fitness<double>, nonstationaryS> {
    template <typename EA>
    int eval_square(EA& ea) {
        double tmp_fit = 0;
      
        int num_org = 0;
        int max_x = get<SPATIAL_X>(ea) - 1; // minus one because SPATIAL_X is the size
        int max_y = get<SPATIAL_Y>(ea) - 1;
        
        for (int x=0; x < max_x; ++x) {
            for (int y=0; y<max_y; ++y){
                typename EA::location_type& l = ea.env().location(x,y);
                if (!l.occupied()) {
                    continue;
                }
                
                std::string lt = get<LAST_TASK>(*l.inhabitant(),"");
                if (x == 0 || x == max_x || y == 0 || y == max_y) {
                    if (lt == "not") {
                        ++tmp_fit;
                    }
                } else if (x == 1 || x == (max_x - 1) || y == 1 || y == (max_y -1)) {
                    if (lt == "nand") {
                        ++tmp_fit;
                    }
                } else if (lt == "ornot") {
                    ++tmp_fit;
                }
                
                
                
            }
        }
        

        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_square(sea));
        put<PATTERN_FIT>(f,sea);
        return f;
    }
};


// Can we do this without repeating code?
template <typename EA>
void eval_french_flag(EA& ea) {
    double tmp_fit = 0;
    
    int num_org = 0;
    int max_x = get<SPATIAL_X>(ea);
    int max_y = get<SPATIAL_Y>(ea);
    
    for (int x=0; x < max_x; ++x) {
        for (int y=0; y<max_y; ++y){
            typename EA::location_type& l = ea.env().location(x,y);
            if (!l.occupied()) {
                continue;
            }
            
            std::string lt = get<LAST_TASK>(*l.inhabitant(),"");
            
            if (x < floor(max_x / 3.0)) {
                if (lt == "not") {
                    ++tmp_fit;
                }
            } else if ((x > (floor(max_x/3.0)) && x > ((2* floor(max_x/3.0))-1))) {
                if (lt == "nand") {
                    ++tmp_fit;
                }
            } else if (lt == "ornot") {
                ++tmp_fit;
            }
            
        }
    }
    put<PATTERN_FIT>(tmp_fit, ea);

}

//! French flag fitness - not, nand, ornot
struct french_flag : public fitness_function<unary_fitness<double>, nonstationaryS> {
    template <typename EA>
    int eval_french_flag(EA& ea) {
        double tmp_fit = 0;
        
        int num_org = 0;
        int max_x = get<SPATIAL_X>(ea);
        int max_y = get<SPATIAL_Y>(ea);
        
        for (int x=0; x < max_x; ++x) {
            for (int y=0; y<max_y; ++y){
                typename EA::location_type& l = ea.env().location(x,y);
                if (!l.occupied()) {
                    continue;
                }
                
                std::string lt = get<LAST_TASK>(*l.inhabitant(),"");
                
                if (x < floor(max_x / 3.0)) {
                    if (lt == "not") {
                        ++tmp_fit;
                    }
                } else if ((x > (floor(max_x/3.0)) && x > ((2* floor(max_x/3.0))-1))) {
                    if (lt == "nand") {
                        ++tmp_fit;
                    }
                } else if (lt == "ornot") {
                    ++tmp_fit;
                }
                
            }
        }
        
        
        return tmp_fit;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        double f = static_cast<double>(eval_french_flag(sea));
        put<PATTERN_FIT>(f,sea);
        return f;
    }
};


        
        
/*! Output the dominant individual
 */
LIBEA_ANALYSIS_TOOL(get_dominant) {
    using namespace ealib;
    
    typename EA::individual_type& ind=*analysis::dominant(ea);
    typename EA::individual_type::individual_ptr_type p = ind.population()[0];
    typename EA::individual_type::individual_type::hardware_type::genome_type g = ind.population()[0]->genome();
    
    
    int count = 0;
    for (typename EA::individual_type::genome_type::iterator f=g.begin(); f != g.end(); f++) {
        //name()
        std::string s = (ind.isa()[*f])->name();
        
        std::cout << "repr[" << count << "] = ea.isa()[\"" << s  << "\"];" << std::endl;
        count ++;
        
    }
    
}

#endif
