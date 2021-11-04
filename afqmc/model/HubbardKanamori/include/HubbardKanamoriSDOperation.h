//
// Created by Hao Shi on 1/19/18.
//

#ifndef AFQMCHUBBARDKANAMORI_HUBBARDKANAMORISDOPERATION_H
#define AFQMCHUBBARDKANAMORI_HUBBARDKANAMORISDOPERATION_H

#include "HubbardKanamori.h"
#include "../../../blocks/walker/SD/include/SD.h"

void fillWalkerRandomly(SD &walker, const HubbardKanamori &model);
void fillWalkerFromModel(SD &walker, HubbardKanamori &model);

#endif //AFQMCHUBBARDKANAMORI_HUBBARDKANAMORISDOPERATION_H