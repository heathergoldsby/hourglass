#ifndef _PROGRESSION_TRACKING_
#define _PROGRESSION_TRACKING_

#include <numeric>

#include <ea/digital_evolution.h>
#include <ea/events.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <ea/traits.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>

using namespace ealib;

LIBEA_MD_DECL(FITNESS_MULTIPLIER, "ea.prog.fitness_multiplier", double);

LIBEA_MD_DECL(JUVENILE_EVAL_PERIOD, "ea.prog.adult_evalution_period", unsigned int);
LIBEA_MD_DECL(ADULT_EVAL_PERIOD, "ea.prog.juvenile_evaluation_period", unsigned int);

LIBEA_MD_DECL(JUVENILE_FITNESS, "ea.prog.juvenile_fitness", double);
LIBEA_MD_DECL(JUVENILE_FITNESS_MAX, "ea.prog.juvenile_fitness_max", double);
LIBEA_MD_DECL(ADULT_FITNESS, "ea.prog.adult_fitness", double);
LIBEA_MD_DECL(ADULT_FITNESS_MAX, "ea.prog.adult_fitness_max", double);

LIBEA_MD_DECL(OVERALL_FITNESS, "ea.prog.overall_fitness", double);

template 
< typename EA
, typename JuvenileFitnessFunction
, typename AdultFitnessFunction
> struct progression_event : end_of_update_event<EA> {
    progression_event(EA& ea) : end_of_update_event<EA>(ea) { }

    virtual ~progression_event() { }

    virtual void operator()(EA& ea) override {

        const auto age = ea.current_update();

        if (age >= get<JUVENILE_EVAL_PERIOD>(ea)) {
            const auto juvenile_val = JuvenileFitnessFunction()(ea);

            if (juvenile_val > get<JUVENILE_FITNESS_MAX>(ea, 0.0)) {
                put<JUVENILE_FITNESS_MAX>(juvenile_val, ea);
            }
            
            if (age >= get<ADULT_EVAL_PERIOD>(ea)) {
                const auto adult_val = AdultFitnessFunction()(ea);

                if (adult_val > get<ADULT_FITNESS_MAX>(ea, 0.0)) {
                    put<ADULT_FITNESS_MAX>(adult_val, ea);
                }
            }
        }
    }
};

struct progression_fitness : public fitness_function<unary_fitness<double>, nonstationaryS> {
    template <typename EA>
    double eval_progression(EA& ea) {
        const auto original_fitness = std::pow(get<JUVENILE_FITNESS_MAX>(ea, 0.0), 1.5) + std::pow(get<ADULT_FITNESS_MAX>(ea, 0.0), 2);
        const auto weighted_fitness = original_fitness + (get<FITNESS_MULTIPLIER>(ea, 0.0) * get<OVERALL_FITNESS>(ea, 0.0));

        put<OVERALL_FITNESS>(weighted_fitness, ea);
        return weighted_fitness;
    }
    
    template <typename SubpopulationEA, typename MetapopulationEA>
    double operator()(SubpopulationEA& sea, MetapopulationEA& mea) {
        return eval_progression(sea);
    }
};

template <typename EA>
struct progression_dat : record_statistics_event<EA> {
    progression_dat(EA& ea) : record_statistics_event<EA>(ea), _df("prog_fitness.dat") {
        _df.add_field("update")
        .add_field("min_juvenile_fitness")
        .add_field("min_adult_fitness")
        .add_field("min_overall_fitness")

        .add_field("max_juvenile_fitness")
        .add_field("max_adult_fitness")
        .add_field("max_overall_fitness")

        .add_field("avg_juvenile_fitness")
        .add_field("avg_adult_fitness")
        .add_field("avg_overall_fitness");
    }
    
    virtual ~progression_dat() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::min, tag::mean, tag::max> > juvenile_fit;
        accumulator_set<double, stats<tag::min, tag::mean, tag::max> > adult_fit;
        accumulator_set<double, stats<tag::min, tag::mean, tag::max> > overall_fit;
        
        for(typename EA::iterator i = ea.begin(); i != ea.end(); ++i) {
            juvenile_fit(get<JUVENILE_FITNESS_MAX>(*i));
            adult_fit(get<ADULT_FITNESS_MAX>(*i));
            overall_fit(get<OVERALL_FITNESS>(*i));
        }
        
        _df.write(ea.current_update())
        .write(min(juvenile_fit))
        .write(min(adult_fit))
        .write(min(overall_fit))

        .write(max(juvenile_fit))
        .write(max(adult_fit))
        .write(max(overall_fit))

        .write(mean(juvenile_fit))
        .write(mean(adult_fit))
        .write(mean(overall_fit))
        .endl();
    }
    
    datafile _df;
};

LIBEA_ANALYSIS_TOOL(movie_progression) {

    double max_fit = 0;

    typename EA::individual_type best;
    
    for(typename EA::iterator i = ea.begin(); i != ea.end(); ++i) {
        const double curr_fitness = get<OVERALL_FITNESS>(*i);
        
        if (curr_fitness > max_fit) {
            max_fit = curr_fitness;
            best = *i;
        }
    }

    datafile df("movie.dat");
    df.write(get<SPATIAL_X>(ea));
    df.write(get<SPATIAL_Y>(ea));
    df.endl();

    int update_max = get<METAPOP_COMPETITION_PERIOD>(ea);
    typename EA::individual_type best_founder(*best.traits().founder());

    for (int j = 0; j <= update_max; ++j) {
        best_founder.update();
        df.write(j);

        for (int x = 0; x < get<SPATIAL_X>(ea); ++x) {
            for (int y = 0; y < get<SPATIAL_Y>(ea); ++y){
                typename EA::individual_type::environment_type::location_type* l = &best_founder.env().location(x,y);

                if (l->occupied()) {
                    const std::string lt = get<LAST_TASK>(*l->inhabitant(),"");
                    
                    if (lt == "not") {
                        df.write("1");
                    }
                    else if (lt == "nand") {
                        df.write("2");
                    }
                    else if (lt == "") {
                        df.write("0");
                    }
                } else {
                    df.write("-1");
                }
            }
        }
        df.endl();
    }
    df.endl();
}


#endif
