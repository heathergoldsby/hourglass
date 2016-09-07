//
//  body_plans.h
//  hourglass
//
//  Created by Heather Goldsby on 3/7/16.
//  Copyright © 2016 Michigan State University. All rights reserved.
//

#ifndef new_body_plans2_h
#define new_body_plans2_h




/* 0 - all blue */
template <typename EA>
double body_plan_0 (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
            ++f_00;
        }
    }
    
    double e = f_00;
    double f = pow(1.5, e);
    
    return f;
}

/* a - half blue / half green */
template <typename EA>
double body_plan_a (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    double f_01 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        if (agent_y < 3) {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                ++f_00;
            }
        } else {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                ++f_01;
            }
        }
        
    }
    
    double e = f_00 + f_01;
    double f = pow(1.5, e);
    
    return f;
}



/* b - half blue / half green; right yellow */
template <typename EA>
double body_plan_b (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    double f_01 = 0.0;
    double f_10 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        if (agent_x < 3) {
            if (agent_y < 3) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                    ++f_00;
                }
            } else {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f_01;
                }
            }
        } else {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                ++f_10;
            }
        }
        
    }
    
    double e = f_00 + f_01 + f_10;
    double f = pow(1.5, e);
    
    return f;
}


/* c - blue, red, yellow, green */
template <typename EA>
double body_plan_c (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    double f_01 = 0.0;
    double f_10 = 0.0;
    double f_11 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        if (agent_x < 3) {
            if (agent_y < 3) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                    ++f_00;
                }
            } else {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f_01;
                }
            }
        } else {
            if (agent_y < 3){
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f_11;
                }
            } else {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f_10;
                }
            }
        }
        
    }
    
    double e = f_00 + f_01 + f_10 + f_11;
    double f = pow(1.5, e);
    
    return f;
}






/* d */
template <typename EA>
double body_plan_d (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    double f_01 = 0.0;
    double f_10 = 0.0;
    double f_11 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        
        if (agent_y == 1) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                ++f_10;
            }
        } else if (agent_y == 4) {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                ++f_00;
            }
        } else if (agent_x < 3) {
            if (agent_y < 3) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                    ++f_00;
                }
            } else {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f_01;
                }
            }
        } else {
            if (agent_y < 3){
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f_11;
                }
            } else {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f_10;
                }
            }
        }
        
    }
    
    double e = f_00 + f_01 + f_10 + f_11;
    double f = pow(1.5, e);
    
    return f;
}



/* e */
template <typename EA>
double body_plan_e (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    double f_01 = 0.0;
    double f_10 = 0.0;
    double f_11 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        
        if (agent_x == 0) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                ++f_11;
            }
        } else if (agent_x == 5) {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                ++f_01;
            }
        } else if (agent_x < 3) {
            if (agent_y < 3) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                    ++f_00;
                }
            } else {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f_01;
                }
            }
        } else {
            if (agent_y < 3){
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f_11;
                }
            } else {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f_10;
                }
            }
        }
        
    }
    
    double e = f_00 + f_01 + f_10 + f_11;
    double f = pow(1.5, e);
    
    return f;
}




/* f */
template <typename EA>
double body_plan_f (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    double f_01 = 0.0;
    double f_10 = 0.0;
    double f_11 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        
        if (agent_x==0 || agent_x==5 || agent_y==0 || agent_y==5) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                ++f_11;
            }
        } else if (agent_y == 1) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                ++f_10;
            }
        } else if (agent_y == 4) {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                ++f_00;
            }
        } else if (agent_x < 3) {
            if (agent_y < 3) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                    ++f_00;
                }
            } else {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f_01;
                }
            }
        } else {
            if (agent_y < 3){
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f_11;
                }
            } else {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f_10;
                }
            }
        }
        
    }
    
    double e = f_00 + f_01 + f_10 + f_11;
    double f = pow(1.5, e);
    
    return f;
}




/* g */
template <typename EA>
double body_plan_g (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    double f_01 = 0.0;
    double f_10 = 0.0;
    double f_11 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        
        
        if (agent_y == 1) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                ++f_10;
            }
        } else if (agent_y == 2 || agent_y == 3) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                ++f_11;
            }
        } else if (agent_y == 4) {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                ++f_00;
            }
        } else if (agent_x < 3) {
            if (agent_y < 3) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                    ++f_00;
                }
            } else {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f_01;
                }
            }
        } else {
            if (agent_y < 3){
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f_11;
                }
            } else {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f_10;
                }
            }
        }
        
    }
    
    double e = f_00 + f_01 + f_10 + f_11;
    double f = pow(1.5, e);
    
    return f;
}


/* h */
template <typename EA>
double body_plan_h (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    double f_01 = 0.0;
    double f_10 = 0.0;
    double f_11 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        
        if ((agent_x==1 && agent_y > 0 && agent_y < 5) ||
            (agent_x==4 && agent_y > 0 && agent_y < 5) ||
            (agent_y==1 && agent_x > 0 && agent_x < 5) ||
            (agent_y==4 && agent_x > 0 && agent_x < 5) ){
            
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                ++f_10;
            }
        } else  if (agent_x == 0) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                ++f_11;
            }
        } else if (agent_x == 5) {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                ++f_01;
            }
        } else if (agent_x < 3) {
            if (agent_y < 3) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                    ++f_00;
                }
            } else {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f_01;
                }
            }
        } else {
            if (agent_y < 3){
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f_11;
                }
            } else {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f_10;
                }
            }
        }
        
    }
    
    double e = f_00 + f_01 + f_10 + f_11;
    double f = pow(1.5, e);
    
    return f;
}



/* i */
template <typename EA>
double body_plan_i (int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    
    assert(max_x == 6);
    assert(max_y == 6);
    
    double f_00 = 0.0;
    double f_01 = 0.0;
    double f_10 = 0.0;
    double f_11 = 0.0;
    
    // Compute fitness. All 00
    for (int xy = 0; xy<grid_size; xy++) {
        
        // set the input states...
        int agent_x = floor(xy / max_x);
        int agent_y = xy % max_x;
        
        int p = agent_pos[xy];
        
        // no agent
        if (p == -1) {
            continue;
        }
        
        if (agent_y == 0) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                ++f_10;
            }
        } else if (agent_y == 5) {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                ++f_00;
            }
        } else if (agent_x == 0) {
            if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                ++f_11;
            }
        } else if (agent_x == 5) {
            if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                ++f_01;
            }
        } else if (agent_x < 3) {
            if (agent_y < 3) {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 0)){
                    ++f_00;
                }
            } else {
                if (((as[p]).output(0) == 0) &&  ((as[p]).output(1) == 1)){
                    ++f_01;
                }
            }
        } else {
            if (agent_y < 3){
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 1)){
                    ++f_11;
                }
            } else {
                if (((as[p]).output(0) == 1) &&  ((as[p]).output(1) == 0)){
                    ++f_10;
                }
            }
        }
        
    }
    
    double e = f_00 + f_01 + f_10 + f_11;
    double f = pow(1.5, e);
    
    return f;
}







#endif /* new_body_plans_h */



