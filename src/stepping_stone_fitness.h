/* fitness.h
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
#ifndef _EA_DATAFILES_STEPPING_STONE_FITNESS_H_
#define _EA_DATAFILES_STEPPING_STONE_FITNESS_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <ea/datafile.h>
#include <ea/traits.h>

#include "stepping_stones.h"


namespace ealib {
    namespace datafiles {
        
        /*! Datafile for mean generation and min, mean, and max fitness.
         */
        template <typename EA>
        struct stepping_stone_fitness_dat : record_statistics_event<EA> {
            stepping_stone_fitness_dat(EA& ea) : record_statistics_event<EA>(ea), _df("fitness.dat") {
                _df.add_field("update")
                .add_field("mean_fitness_1")
                .add_field("max_fitness_1")
                .add_field("mean_fitness_2")
                .add_field("max_fitness_2")
                .add_field("mean_fitness_3")
                .add_field("max_fitness_3")
                .add_field("mean_fitness_4")
                .add_field("max_fitness_4")
                .add_field("mean_fitness_5")
                .add_field("max_fitness_5")
                .add_field("mean_fitness_6")
                .add_field("max_fitness_6")
                .add_field("mean_fitness_7")
                .add_field("max_fitness_7")
                .add_field("mean_fitness_8")
                .add_field("max_fitness_8")
                .add_field("mean_fitness_9")
                .add_field("max_fitness_9");
            }
            
            virtual ~stepping_stone_fitness_dat() {
            }
            
            virtual void operator()(EA& ea) {
                using namespace boost::accumulators;
                accumulator_set<double, stats<tag::mean, tag::max> > fit1;
                accumulator_set<double, stats<tag::mean, tag::max> > fit2;
                accumulator_set<double, stats<tag::mean, tag::max> > fit3;
                accumulator_set<double, stats<tag::mean, tag::max> > fit4;
                accumulator_set<double, stats<tag::mean, tag::max> > fit5;
                accumulator_set<double, stats<tag::mean, tag::max> > fit6;
                accumulator_set<double, stats<tag::mean, tag::max> > fit7;
                accumulator_set<double, stats<tag::mean, tag::max> > fit8;
                accumulator_set<double, stats<tag::mean, tag::max> > fit9;
                
                for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                    fit1(get<F1>(*i));
                    fit2(get<F2>(*i));
                    fit3(get<F3>(*i));
                    fit4(get<F4>(*i));
                    fit5(get<F5>(*i));
                    fit6(get<F6>(*i));
                    fit7(get<F7>(*i));
                    fit8(get<F8>(*i));
                    fit9(get<F9>(*i));
                }
                
                _df.write(ea.current_update())
                .write(mean(fit1))
                .write(boost::accumulators::max(fit1))
                .write(mean(fit2))
                .write(boost::accumulators::max(fit2))
                .write(mean(fit3))
                .write(boost::accumulators::max(fit3))
                .write(mean(fit4))
                .write(boost::accumulators::max(fit4))
                .write(mean(fit5))
                .write(boost::accumulators::max(fit5))
                .write(mean(fit6))
                .write(boost::accumulators::max(fit6))
                .write(mean(fit7))
                .write(boost::accumulators::max(fit7))
                .write(mean(fit8))
                .write(boost::accumulators::max(fit8))
                .write(mean(fit9))
                .write(boost::accumulators::max(fit9))
                .endl();
            }
            
            datafile _df;
        };
        
    } // datafiles
} // ealib

#endif
