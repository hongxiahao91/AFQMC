//
// Created by Hao Shi on 1/19/18.
//

#ifndef AFQMCHUBBARDKANAMORIREAL_HUBBARDKANAMORIREALSDOPERATION_H
#define AFQMCHUBBARDKANAMORIREAL_HUBBARDKANAMORIREALSDOPERATION_H

#include "HubbardKanamoriReal.h"
#include "../../../blocks/walker/SD/include/SD.h"

void fillWalkerRandomly(SD &walker, const HubbardKanamoriReal &model);
void fillWalkerFromModel(SD &walker, HubbardKanamoriReal &model);

#endif //AFQMCHUBBARDKANAMORIREAL_HUBBARDKANAMORIREALSDOPERATION_H