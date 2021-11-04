//
// Created by Hao Shi on 9/7/17.
//

#ifndef AFQMCLAB_REALMATERIALMOLECULEMDCAS2SOPERATION_H
#define AFQMCLAB_REALMATERIALMOLECULEMDCAS2SOPERATION_H

#include "realMaterialMolecule.h"
#include "../../../blocks/walker/MDCas2s/include/MDCas2s.h"

void fillWalkerRandomly(MDCas2s &walker, const RealMaterialMolecule &model);
void fillWalkerFromModel(MDCas2s &walker, RealMaterialMolecule &model);

#endif //AFQMCLAB_REALMATERIALMOLECULEMDCAS2SOPERATION_H
