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
#include "stepping_stones.h"

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
            df.write("popsize");
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
    
        
        template <typename Individual, typename EA>
        void generate_one_movie(std::string file_name, Individual& best, EA& ea) {
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
            
            datafile df(file_name);
            df.write(max_x);
            df.write(max_y);
            df.endl();
            
            for (int i=0; i<grid_size; ++i) {
                exec_order[i] = i;
            }
            
            
            std::vector< std::vector<int> > cell_color(grid_size, std::vector<int>(2, 0));
            int world_updates = get<WORLD_UPDATES>(ea,10);
            for(int t=0;t<world_updates;t++){
                update_world_stigmergic_communication_N(1, t, agent_pos, exec_order, as, cell_color, best, ea.rng(), ea);
                
                
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
//                    std::string f = "movie_blah.dat";
//                    generate_one_movie (f, *i, ea);
                    
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
            
            
            // get the "prototype" phenotype (markov network):
            
            int offset = get<START_POS>(ea,0);
            
            
            for (int q=0; q<get<NUM_START_AGENTS>(ea,1); q++) {
                as.push_back(N); //grid_size, N); // my agents or networks
                int pos = q + offset;
                if (pos >=  grid_size) {
                    pos -= grid_size;
                }
                agent_pos[pos] = q;
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
                update_world_stigmergic_communication_N(1, t, agent_pos, exec_order, as, cell_color, best, ea.rng(), ea);
                
                
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
    
        LIBEA_ANALYSIS_TOOL(markov_movie_ko) {
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
            
            std::string f = "movie.dat";
            generate_one_movie (f, best, ea);
            
            get<CAPABILITIES_OFF>(ea, "") = "stigmergic";
            f = "movie_ko_stigmergic.dat";
            generate_one_movie (f, best, ea);

            get<CAPABILITIES_OFF>(ea, "") = "edge";
            f = "movie_ko_edge.dat";
            generate_one_movie (f, best, ea);
            
            get<CAPABILITIES_OFF>(ea, "") = "communication";
            f = "movie_ko_communication.dat";
            generate_one_movie (f, best, ea);
            
            get<CAPABILITIES_OFF>(ea, "") = "neighbor";
            f = "movie_ko_neighbor.dat";
            generate_one_movie (f, best, ea);
            
            get<CAPABILITIES_OFF>(ea, "") = "origin";
            f = "movie_ko_origin.dat";
            generate_one_movie (f, best, ea);


        }
        
        
        LIBEA_ANALYSIS_TOOL(markov_movie_island) {
            double max_fit = 0;
            int count = 0;
            
            // recalc all fitness values
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                int z =  get<ISLAND>(*i);
                int q = 0;
                
                int island =  get<ISLAND>(*i);
                stringstream ss;
                ss << count;
                std::string c = ss.str();
                std::string f = "movie_" + c + ".dat";
                
                for (typename EA::individual_type::iterator j=i->begin(); j!=i->end(); ++j) {
                

                    recalculate_fitness(*j, *i);
                    double tmp_fit = static_cast<int>(ealib::fitness(*j,*i));
                    if (tmp_fit > max_fit) {
                        max_fit = tmp_fit;
                        generate_one_movie (f, *j, *i);
                    }
                }
                


                ++count;
                max_fit = 0;
            }
        }
        
        LIBEA_ANALYSIS_TOOL(markov_movie_ko_island) {
            double max_fit = 0;
            int count = 0;
            typename EA::individual_type best;

            // recalc all fitness values
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                int island =  get<ISLAND>(*i);
                stringstream ss;
                ss << island;
                std::string c = ss.str();
                
                for (typename EA::individual_type::iterator j=i->begin(); j!=i->end(); ++j) {
                    recalculate_fitness(*j, *i);
                    double tmp_fit = static_cast<int>(ealib::fitness(*j,*i));
                    if (tmp_fit > max_fit) {
                        max_fit = tmp_fit;
                        std::string f = "movie_" + c + ".dat";
                        generate_one_movie (f, *j, *i);
                        
                        get<CAPABILITIES_OFF>(*i, "") = "stigmergic";
                        f = "movie_ko_stigmergic_" + c + ".dat";
                        generate_one_movie (f, *j, *i);
                        
                        get<CAPABILITIES_OFF>(*i, "") = "edge";
                        f = "movie_ko_edge_" + c + ".dat" ;
                        generate_one_movie (f, *j, *i);
                        
                        get<CAPABILITIES_OFF>(*i, "") = "communication";
                        f = "movie_ko_communication" + c + ".dat";
                        generate_one_movie (f, *j, *i);
                        
                        get<CAPABILITIES_OFF>(*i, "") = "neighbor";
                        f = "movie_ko_neighbor" + c + ".dat";
                        generate_one_movie (f, *j, *i);
                        
                        get<CAPABILITIES_OFF>(*i, "") = "origin";
                        f = "movie_ko_origin" + c + ".dat";
                        generate_one_movie (f, *j, *i);
                    }
                }
                

                
                ++count;
                max_fit = 0;
            }
        }
        
        
        LIBEA_ANALYSIS_TOOL(markov_ko_island) {
            double max_fit = 0;
            int count = 0;
            typename EA::individual_type best;
            
            
            datafile df("ko_island.dat");
            df.write("treatment")
            .write("island")
            .write("treatment_island")
            .write("mean")
            .write("max");
            
            df.endl();
            
            // recalc all fitness values
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                int island =  get<ISLAND>(*i);
                stringstream ss;
                ss << island;
                std::string c = ss.str();
                
                accumulator_set<double, stats<tag::max, tag::mean> > control;
                accumulator_set<double, stats<tag::max, tag::mean> > stig;
                accumulator_set<double, stats<tag::max, tag::mean> > edge;
                accumulator_set<double, stats<tag::max, tag::mean> > comm;
                accumulator_set<double, stats<tag::max, tag::mean> > neigh;
                accumulator_set<double, stats<tag::max, tag::mean> > ori;

                
                
                // recalc all fitness values
                for (typename EA::individual_type::iterator j=i->begin(); j!=i->end(); ++j) {

                    recalculate_fitness(*j, *i);
                    double tmp_fit = static_cast<int>(ealib::fitness(*j,*i));
                    control(tmp_fit);
                    
                    typename EA::individual_type::individual_type ko_stigmergic = *j;
                    put<CAPABILITIES_OFF>("stigmergic", *i);
                    recalculate_fitness(ko_stigmergic, *i);
                    tmp_fit = static_cast<double>(ealib::fitness(ko_stigmergic,*i));
                    stig(tmp_fit);
                    
                    typename EA::individual_type::individual_type ko_edge = *j;
                    put<CAPABILITIES_OFF>("edge", *i);
                    recalculate_fitness(ko_edge, *i);
                    tmp_fit = static_cast<double>(ealib::fitness(ko_edge,*i));
                    edge(tmp_fit);
                    
                    typename EA::individual_type::individual_type ko_comm = *j;
                    put<CAPABILITIES_OFF>("communication", *i);
                    recalculate_fitness(ko_comm, *i);
                    tmp_fit = static_cast<double>(ealib::fitness(ko_comm,*i));
                    comm(tmp_fit);
                    
                    typename EA::individual_type::individual_type ko_neighbor = *j;
                    put<CAPABILITIES_OFF>("neighbor", *i);
                    recalculate_fitness(ko_neighbor, *i);
                    tmp_fit = static_cast<double>(ealib::fitness(ko_neighbor,*i));
                    neigh(tmp_fit);
                    
                    typename EA::individual_type::individual_type ko_origin = *j;
                    put<CAPABILITIES_OFF>("origin", *i);
                    recalculate_fitness(ko_origin, *i);
                    tmp_fit = static_cast<double>(ealib::fitness(ko_origin,*i));
                    ori(tmp_fit);
                    
                }
                
                //                        f = "movie_ko_neighbor" + c + ".dat";

                df.write("control");
                df.write(c);
                df.write("control_" + c);
                df.write(boost::accumulators::mean(control));
                df.write(boost::accumulators::max(control));
                df.endl();

                df.write("ko_communication");
                df.write(c);
                df.write("ko_communication_" + c);
                df.write(boost::accumulators::mean(comm));
                df.write(boost::accumulators::max(comm));
                df.endl();
                
                df.write("ko_stigmergic");
                df.write(c);
                df.write("ko_stigmergic_" + c);
                df.write(boost::accumulators::mean(stig));
                df.write(boost::accumulators::max(stig));
                df.endl();
                
                df.write("ko_edge");
                df.write(c);
                df.write("ko_edge_" + c);
                df.write(boost::accumulators::mean(edge));
                df.write(boost::accumulators::max(edge));
                df.endl();
                
                df.write("ko_neighbor");
                df.write(c);
                df.write("ko_neighbor_" + c);
                df.write(boost::accumulators::mean(neigh));
                df.write(boost::accumulators::max(neigh));
                df.endl();
                
                df.write("ko_origin");
                df.write(c);
                df.write("ko_origin_" + c);
                df.write(boost::accumulators::mean(ori));
                df.write(boost::accumulators::max(ori));
                df.endl();
                
            }
            
            
                
                ++count;
        }
    
        


        
        /* how the full population of solutions uses these techniques.*/
        LIBEA_ANALYSIS_TOOL(ko) {
            
            accumulator_set<double, stats<tag::max, tag::mean> > control;
            accumulator_set<double, stats<tag::max, tag::mean> > stig;
            accumulator_set<double, stats<tag::max, tag::mean> > edge;
            accumulator_set<double, stats<tag::max, tag::mean> > comm;
            accumulator_set<double, stats<tag::max, tag::mean> > neigh;
            accumulator_set<double, stats<tag::max, tag::mean> > ori;

            
            datafile df("ko.dat");
            df.write("treatment")
            .write("mean")
            .write("max");
            
            df.endl();
            
            // recalc all fitness values
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                recalculate_fitness(*i, ea);
                double tmp_fit = static_cast<int>(ealib::fitness(*i,ea));
                control(tmp_fit);
                
                typename EA::individual_type ko_stigmergic = *i;
                put<CAPABILITIES_OFF>("stigmergic", ea);
                recalculate_fitness(ko_stigmergic, ea);
                tmp_fit = static_cast<double>(ealib::fitness(ko_stigmergic,ea));
                stig(tmp_fit);
                
                typename EA::individual_type ko_edge = *i;
                put<CAPABILITIES_OFF>("edge", ea);
                recalculate_fitness(ko_edge, ea);
                tmp_fit = static_cast<double>(ealib::fitness(ko_edge,ea));
                edge(tmp_fit);
                
                typename EA::individual_type ko_comm = *i;
                put<CAPABILITIES_OFF>("communication", ea);
                recalculate_fitness(ko_comm, ea);
                tmp_fit = static_cast<double>(ealib::fitness(ko_comm,ea));
                comm(tmp_fit);
                
                typename EA::individual_type ko_neighbor = *i;
                put<CAPABILITIES_OFF>("neighbor", ea);
                recalculate_fitness(ko_neighbor, ea);
                tmp_fit = static_cast<double>(ealib::fitness(ko_neighbor,ea));
                neigh(tmp_fit);
                
                typename EA::individual_type ko_origin = *i;
                put<CAPABILITIES_OFF>("origin", ea);
                recalculate_fitness(ko_origin, ea);
                tmp_fit = static_cast<double>(ealib::fitness(ko_origin,ea));
                ori(tmp_fit);
                

            }
            

            df.write("control");
            df.write(boost::accumulators::mean(control));
            df.write(boost::accumulators::max(control));
            df.endl();
            
            df.write("ko_communication");
            df.write(boost::accumulators::mean(comm));
            df.write(boost::accumulators::max(comm));
            df.endl();
            
            df.write("ko_stigmergic");
            df.write(boost::accumulators::mean(stig));
            df.write(boost::accumulators::max(stig));
            df.endl();
            
            df.write("ko_edge");
            df.write(boost::accumulators::mean(edge));
            df.write(boost::accumulators::max(edge));
            df.endl();
            
            df.write("ko_neighbor");
            df.write(boost::accumulators::mean(neigh));
            df.write(boost::accumulators::max(neigh));
            df.endl();
            
            df.write("ko_origin");
            df.write(boost::accumulators::mean(ori));
            df.write(boost::accumulators::max(ori));
            df.endl();
            
            
            
            

        }
        
        
        
       
        
        /* how the full population of solutions uses these techniques.*/
        LIBEA_ANALYSIS_TOOL(ko_stepping_stones) {
            
            accumulator_set<double, stats<tag::max, tag::mean> > control6;
            accumulator_set<double, stats<tag::max, tag::mean> > stig6;
            accumulator_set<double, stats<tag::max, tag::mean> > edge6;
            accumulator_set<double, stats<tag::max, tag::mean> > comm6;
            accumulator_set<double, stats<tag::max, tag::mean> > neigh6;
            accumulator_set<double, stats<tag::max, tag::mean> > ori6;
            accumulator_set<double, stats<tag::max, tag::mean> > control7;
            accumulator_set<double, stats<tag::max, tag::mean> > stig7;
            accumulator_set<double, stats<tag::max, tag::mean> > edge7;
            accumulator_set<double, stats<tag::max, tag::mean> > comm7;
            accumulator_set<double, stats<tag::max, tag::mean> > neigh7;
            accumulator_set<double, stats<tag::max, tag::mean> > ori7;
            accumulator_set<double, stats<tag::max, tag::mean> > control8;
            accumulator_set<double, stats<tag::max, tag::mean> > stig8;
            accumulator_set<double, stats<tag::max, tag::mean> > edge8;
            accumulator_set<double, stats<tag::max, tag::mean> > comm8;
            accumulator_set<double, stats<tag::max, tag::mean> > neigh8;
            accumulator_set<double, stats<tag::max, tag::mean> > ori8;
            accumulator_set<double, stats<tag::max, tag::mean> > control9;
            accumulator_set<double, stats<tag::max, tag::mean> > stig9;
            accumulator_set<double, stats<tag::max, tag::mean> > edge9;
            accumulator_set<double, stats<tag::max, tag::mean> > comm9;
            accumulator_set<double, stats<tag::max, tag::mean> > neigh9;
            accumulator_set<double, stats<tag::max, tag::mean> > ori9;

            
            datafile df("ko.dat");
            df.write("treatment")
            .write("mean")
            .write("max");
            
            df.endl();
            
            // recalc all fitness values
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                recalculate_fitness(*i, ea);
                double tmp_fit = static_cast<int>(ealib::fitness(*i,ea));
                control6(get<F6>(*i));
                control7(get<F7>(*i));
                control8(get<F8>(*i));
                control9(get<F9>(*i));
                
                typename EA::individual_type ko_stigmergic = *i;
                put<CAPABILITIES_OFF>("stigmergic", ea);
                recalculate_fitness(ko_stigmergic, ea);
                stig6(get<F6>(*i));
                stig7(get<F7>(*i));
                stig8(get<F8>(*i));
                stig9(get<F9>(*i));
           
                typename EA::individual_type ko_edge = *i;
                put<CAPABILITIES_OFF>("edge", ea);
                recalculate_fitness(ko_edge, ea);
                tmp_fit = static_cast<double>(ealib::fitness(ko_edge,ea));
                edge6(get<F6>(*i));
                edge7(get<F7>(*i));
                edge8(get<F8>(*i));
                edge9(get<F9>(*i));
                
                typename EA::individual_type ko_comm = *i;
                put<CAPABILITIES_OFF>("communication", ea);
                recalculate_fitness(ko_comm, ea);
                tmp_fit = static_cast<double>(ealib::fitness(ko_comm,ea));
                comm6(get<F6>(*i));
                comm7(get<F7>(*i));
                comm8(get<F8>(*i));
                comm9(get<F9>(*i));
                
                typename EA::individual_type ko_neighbor = *i;
                put<CAPABILITIES_OFF>("neighbor", ea);
                recalculate_fitness(ko_neighbor, ea);
                tmp_fit = static_cast<double>(ealib::fitness(ko_neighbor,ea));
                neigh6(get<F6>(*i));
                neigh7(get<F7>(*i));
                neigh8(get<F8>(*i));
                neigh9(get<F9>(*i));
                
                typename EA::individual_type ko_origin = *i;
                put<CAPABILITIES_OFF>("origin", ea);
                recalculate_fitness(ko_origin, ea);
                tmp_fit = static_cast<double>(ealib::fitness(ko_origin,ea));
                ori6(get<F6>(*i));
                ori7(get<F7>(*i));
                ori8(get<F8>(*i));
                ori9(get<F9>(*i));
                
                
            }
            
            
            df.write("control").write("6").write("control_6");
            df.write(boost::accumulators::mean(control6));
            df.write(boost::accumulators::max(control6));
            df.endl();

            df.write("control").write("7").write("control_7");
            df.write(boost::accumulators::mean(control7));
            df.write(boost::accumulators::max(control7));
            df.endl();
            
            df.write("control").write("8").write("control_8");
            df.write(boost::accumulators::mean(control8));
            df.write(boost::accumulators::max(control8));
            df.endl();
            
            df.write("control").write("9").write("control_9");
            df.write(boost::accumulators::mean(control9));
            df.write(boost::accumulators::max(control9));
            df.endl();
            
            df.write("ko_stigmergic").write("6").write("ko_stigmergic_6");
            df.write(boost::accumulators::mean(stig6));
            df.write(boost::accumulators::max(stig6));
            df.endl();
            
            df.write("ko_stigmergic").write("7").write("ko_stigmergic_7");
            df.write(boost::accumulators::mean(stig7));
            df.write(boost::accumulators::max(stig7));
            df.endl();
            
            df.write("ko_stigmergic").write("8").write("ko_stigmergic_8");
            df.write(boost::accumulators::mean(stig8));
            df.write(boost::accumulators::max(stig8));
            df.endl();
            
            df.write("ko_stigmergic").write("9").write("ko_stigmergic_9");
            df.write(boost::accumulators::mean(stig9));
            df.write(boost::accumulators::max(stig9));
            df.endl();
            
            df.write("ko_edge").write("7").write("ko_edge_6");
            df.write(boost::accumulators::mean(edge6));
            df.write(boost::accumulators::max(edge6)).endl();
            
            df.write("ko_edge").write("7").write("ko_edge_7");
            df.write(boost::accumulators::mean(edge7));
            df.write(boost::accumulators::max(edge7)).endl();
            
            df.write("ko_edge").write("8").write("ko_edge_8");
            df.write(boost::accumulators::mean(edge8));
            df.write(boost::accumulators::max(edge8)).endl();
            
            df.write("ko_edge").write("9").write("ko_edge_9");
            df.write(boost::accumulators::mean(edge9));
            df.write(boost::accumulators::max(edge9)).endl();
            
            df.write("ko_neighbor").write("6").write("ko_neighbor_6");
            df.write(boost::accumulators::mean(neigh6));
            df.write(boost::accumulators::max(neigh6)).endl();
            df.write("ko_neighbor").write("7").write("ko_neighbor_7");
            df.write(boost::accumulators::mean(neigh7));
            df.write(boost::accumulators::max(neigh7)).endl();
            df.write("ko_neighbor").write("8").write("ko_neighbor_8");
            df.write(boost::accumulators::mean(neigh8));
            df.write(boost::accumulators::max(neigh8)).endl();
            df.write("ko_neighbor").write("9").write("ko_neighbor_9");
            df.write(boost::accumulators::mean(neigh9));
            df.write(boost::accumulators::max(neigh9)).endl();
            
            df.write("ko_origin").write("6").write("ko_origin_6");
            df.write(boost::accumulators::mean(ori6));
            df.write(boost::accumulators::max(ori6)).endl();
            df.write("ko_origin").write("7").write("ko_origin_7");
            df.write(boost::accumulators::mean(ori7));
            df.write(boost::accumulators::max(ori7)).endl();
            df.write("ko_origin").write("8").write("ko_origin_8");
            df.write(boost::accumulators::mean(ori8));
            df.write(boost::accumulators::max(ori8)).endl();
            df.write("ko_origin").write("9").write("ko_origin_9");
            df.write(boost::accumulators::mean(ori9));
            df.write(boost::accumulators::max(ori9)).endl();
            
            df.write("ko_communication").write("6").write("ko_communication_6");
            df.write(boost::accumulators::mean(comm6));
            df.write(boost::accumulators::max(comm6)).endl();
            df.write("ko_communication").write("7").write("ko_communication_7");
            df.write(boost::accumulators::mean(comm7));
            df.write(boost::accumulators::max(comm7)).endl();
            df.write("ko_communication").write("8").write("ko_communication_8");
            df.write(boost::accumulators::mean(comm8));
            df.write(boost::accumulators::max(comm8)).endl();
            df.write("ko_communication").write("9").write("ko_communication_9");
            df.write(boost::accumulators::mean(comm9));
            df.write(boost::accumulators::max(comm9)).endl();
            
            
            
        }
        
        



    
    }
}
#endif
