#include <ea/digital_evolution.h>
#include <ea/events.h>

using namespace ealib;

LIBEA_MD_DECL(JUVENILE_EVAL_PERIOD, "ea.prog.adult_evalution_period", unsigned int);
LIBEA_MD_DECL(ADULT_EVAL_PERIOD, "ea.prog.juvenile_evaluation_period", unsigned int);

LIBEA_MD_DECL(JUVENILE_FITNESS, "ea.prog.juvenile_fitness", double);
LIBEA_MD_DECL(JUVENILE_FITNESS_MAX, "ea.prog.juvenile_fitness_max", double);
LIBEA_MD_DECL(ADULT_FITNESS, "ea.prog.adult_fitness", double);
LIBEA_MD_DECL(ADULT_FITNESS_MAX, "ea.prog.adult_fitness_max", double);

template 
< typename EA
, typename JuvenileFitnessFunction
, typename AdultFitnessFunction
> struct progression_event : end_of_update_event<EA> {

    progression_event(EA& ea) : end_of_update_event<EA>(ea) { }

    virtual ~progression_event() { }

    virtual void operator()(EA& ea) override {

        const auto age = ea.current_update();

        if (age > get<JUVENILE_EVAL_PERIOD>(ea)) {
            const auto juvenile_val = JuvenileFitnessFunction()(ea);

            if (juvenile_val > get<JUVENILE_FITNESS_MAX>(ea, 0.0)) {
                put<JUVENILE_FITNESS_MAX>(juvenile_val, ea);
            }
            
            if (age > get<ADULT_EVAL_PERIOD>(ea)) {
                const auto adult_val = AdultFitnessFunction()(ea);

                if (adult_val > get<ADULT_FITNESS_MAX>(ea, 0.0)) {
                    put<ADULT_FITNESS_MAX>(adult_val, ea);
                }
            }
        }
    }
};

// template <typename EA>
// struct test_event : periodic_event<1, EA> {

//     test_event(EA& ea) : periodic_event<1, EA>(ea) { }

//     virtual ~test_event() { }

//     virtual void operator()(EA& ea) override {
//     }
// };
