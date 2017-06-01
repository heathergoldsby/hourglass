

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>

#include <ea/mkv/markov_network_evolution.h>
#include <ea/generational_models/moran_process.h>
#include <ea/fitness_function.h>
#include <ea/cmdline_interface.h>
#include <ea/datafiles/fitness.h>
#include <ea/line_of_descent.h>

#include "markov_movie.h"
#include "hourglass.h"
#include "new_body_plans2.h"

using namespace std;
using namespace boost::accumulators;
using namespace ealib;
using namespace mkv;





/*! Sample fitness function for Markov networks.
 */
struct big_world : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    template <typename Individual, typename RNG, typename EA>
    double operator()(Individual& ind, RNG& rng, EA& ea) {
        // fixed LARGE grid size
        int max_x = 100;
        int max_y = 100;
        int grid_size = max_x * max_y;
        vector<int> agent_pos (grid_size, -1);
        vector<int> exec_order (grid_size);
        double f=0.0;
        
        
        // get the "prototype" phenotype (markov network):
        typename EA::phenotype_type &N = ealib::phenotype(ind, ea);
        vector<typename EA::phenotype_type> as; //
        
        
        int offset = get<START_POS>(ea,0);
        for (int q=0; q<get<NUM_START_AGENTS>(ea,1); q++) {
            as.push_back(N); //grid_size, N); // my agents or networks
            int pos = q + offset;
            if (pos >=  grid_size) {
                pos -= grid_size;
            }
            agent_pos[pos] = q;
            as[q].reset(rng.seed());
        }
        
        for (int i=0; i<grid_size; ++i) {
            exec_order[i] = i;
        }
        
        
        std::vector< std::vector<int> > cell_color(grid_size, std::vector<int>(2, 0));
        update_big_world_N(get<WORLD_UPDATES>(ea,10), 0, agent_pos, exec_order, as, cell_color, ind, rng, ea);
        
        /* translate agent positions... */
        int upper_left = -1; // upper left individual pops up first
        
        for (int i =0; i<agent_pos.size(); i++) {
            int p = agent_pos[i];
            if (p != -1) {
                upper_left = i;
                break;
            }
        }

        // max digital tissue size. Must be much smaller than the world size (100x100)

        int ul_y = floor(upper_left / max_x);
        int ul_x = upper_left % max_x;
        int dt_max_x = get<X_SIZE>(ea,10);
        int dt_max_y = get<Y_SIZE>(ea,10);
        int dt_grid = dt_max_x * dt_max_y;
        
        vector<int> agent_pos_small (dt_grid, -1);
        

        int count = 0;
        for (int j=0; j<dt_max_y; j++){
            for (int i=0; i<dt_max_x; i++) {
                int xy = upper_left + i + (j*max_y);
                if (xy > grid_size) {break;}
                agent_pos_small[count] = agent_pos[xy];
                int z = agent_pos[xy];
                count++;
            }
        }
        
        

        
        
        
        int fit_func = get<BODYPLAN>(ea,0) ;
        switch(fit_func) {
            case 2:
                f = body_plan_a(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 3:
                f = body_plan_b(dt_grid,dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 4:
                f = body_plan_c(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 5:
                f = body_plan_d(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 6:
                f = body_plan_e(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 7:
                f = body_plan_f(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 8:
                f = body_plan_g(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 9:
                f = body_plan_h(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 10:
                f = body_plan_i(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 100:
                f = body_plan_solid(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
            case 101:
                f = body_plan_drift(dt_grid, dt_max_x, dt_max_y, agent_pos_small, as, ea);
                break;
                
        }
        
        if (f == 0) {
            f = 1;
        }
        return (f);
        
        
        
        
    }
};

// Evolutionary algorithm definition.
typedef markov_network_evolution
< big_world
, recombination::asexual
, generational_models::moran_process< >
, dont_stop
, fill_population
, lifecycles::markov_network_lifecycle
, lod_with_fitness_trait
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
        add_option<HIDDEN_BIT_ERROR_PROB>(this);
        add_option<OUTPUT_BIT_ERROR_PROB>(this);
        
        add_option<CAPABILITIES_OFF>(this);
        add_option<AGENT_DEATH_PROB>(this);
        add_option<NUM_START_AGENTS>(this);
        add_option<BODYPLAN>(this);
        add_option<APOP_THRESH>(this);
        add_option<START_POS>(this);
        add_option<RAND_ORDER>(this);
        add_option<NO_REP_PERIOD>(this);
        add_option<GATE_EPSILON>(this);
        add_option<ANALYSIS_INPUT>(this);
        
        
    }
    
    
    virtual void gather_tools() {
        add_tool<analysis::dominant_genetic_graph>(this);
        add_tool<analysis::dominant_causal_graph>(this);
        add_tool<analysis::dominant_reduced_graph>(this);
        
        add_tool<ealib::analysis::markov_movie_big_world>(this);
        add_tool<ealib::analysis::lod_movie>(this);
        add_tool<ealib::analysis::recalc_fit>(this);
        add_tool<ealib::analysis::ko>(this);
        
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
        add_event<lod_event>(ea);
        add_event<datafiles::mrca_lineage>(ea);
    };
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);



