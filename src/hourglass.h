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

using namespace ealib;


LIBEA_MD_DECL(X_SIZE, "ea.hourglass.x_size", int); // x dimension
LIBEA_MD_DECL(Y_SIZE, "ea.hourglass.y_size", int); // y dimension
LIBEA_MD_DECL(BRAIN_UPDATES, "ea.hourglass.brain_updates", int); // number of brain updates per world update
LIBEA_MD_DECL(WORLD_UPDATES, "ea.hourglass.world_updates", int); // number of world updates per fitness eval
LIBEA_MD_DECL(FIT_GAMMA, "ea.hourglass.fit_gamma", double);




#endif
