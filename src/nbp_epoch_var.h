//
//  nbp_epoch_var.h
//  hourglass
//
//  Created by Heather Goldsby on 4/10/17.
//  Copyright © 2017 Michigan State University. All rights reserved.
//

#ifndef nbp_epoch_var_h
#define nbp_epoch_var_h
#include "hourglass.h"
#include "new_body_plans2.h"
#include "new_body_plans.h"
#include "body_plans.h"

template <typename EA>
double eval_body_plan (int ffToUse, int grid_size, int max_x, int max_y, std::vector<int>& agent_pos, std::vector<typename EA::phenotype_type>& as, EA& ea) {
    double f = 0.0;
    
    switch(ffToUse) {
        case 1:
            f = body_plan_a(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 2:
            f = body_plan_b(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 3:
            f = body_plan_c(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 4:
            f = body_plan_d(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 5:
            f = body_plan_e(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 6:
            f = body_plan_f(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 7:
            f = body_plan_g(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 8:
            f = body_plan_h(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 9:
            f = body_plan_i(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 10:
            f = body_plan_a1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 11:
            f = body_plan_b1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 12:
            f = body_plan_c1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 13:
            f = body_plan_d1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 14:
            f = body_plan_e1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 15:
            f = body_plan_f1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 16:
            f = body_plan_g1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 17:
            f = body_plan_h1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 18:
            f = body_plan_i1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 19:
            f = body_plan_j1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 20:
            f = body_plan_k1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 21:
            f = body_plan_l1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 22:
            f = body_plan_m1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 23:
            f = body_plan_n1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 24:
            f = body_plan_o1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 25:
            f = body_plan_p1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 26:
            f = body_plan_q1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 27:
            f = body_plan_r1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 28:
            f = body_plan_s1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 29:
            f = body_plan_t1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 30:
            f = body_plan_u1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 31:
            f = body_plan_v1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 32:
            f = body_plan_w1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 33:
            f = body_plan0(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 34:
            f = body_plan1(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 35:
            f = body_plan2(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 36:
            f = body_plan3(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 37:
            f = body_plan4(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 38:
            f = body_plan5(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 39:
            f = body_plan6(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 40:
            f = body_plan7(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 41:
            f = body_plan8(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 42:
            f = body_plan9(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 43:
            f = body_plan10(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 44:
            f = body_plan11(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 45:
            f = body_plan12(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 46:
            f = body_plan13(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 47:
            f = body_plan14(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 48:
            f = body_plan15(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 49:
            f = body_plan16(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 50:
            f = body_plan17(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 51:
            f = body_plan18(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
        case 52:
            f = body_plan19(grid_size, max_x, max_y, agent_pos, as, ea);
            break;
    }
    
    return f;
    
}


#endif /* nbp_epoch_var_h */
