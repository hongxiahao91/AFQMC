//
// Created by Hao Shi on 10/4/17.
//

#ifndef AFQMCLAB_HOP2ISSD2SOPERATION_H
#define AFQMCLAB_HOP2ISSD2SOPERATION_H

#include "../../../walker/SD2s/include/SD2s.h"
#include "../../../oneBodyOperator/hop2is/include/hop2is.h"

class Hop2isSD2sOperation
{
 public:
    Hop2isSD2sOperation();
    ~Hop2isSD2sOperation();

    void applyToRight(const Hop2is &oneBody, const SD2s &walker, SD2s &walkerNew) const;
    void applyToLeft(const Hop2is &oneBody, const SD2s &walker, SD2s &walkerNew) const;

 private:
    void checkAndResize(const Hop2is &oneBody, const SD2s &walker, SD2s &walkerNew) const;
};

#endif //AFQMCLAB_HOP2ISSD2SOPERATION_H
