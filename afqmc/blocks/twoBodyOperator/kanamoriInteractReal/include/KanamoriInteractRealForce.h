//
// Created by Hao Shi on 1/12/18.
//

#ifndef AFQMCHUBBARDKANAMORI_KANAMORIINTERACTREALFORCE_H
#define AFQMCHUBBARDKANAMORI_KANAMORIINTERACTREALFORCE_H

#include "../../../../../common/tensorHao/include/tensor_all.h"

class KanamoriInteractRealForce
{
 public:
    tensor_hao::TensorHao<double, 1> UForce;
    tensor_hao::TensorHao<double, 1> U1Force;
    tensor_hao::TensorHao<double, 1> U2JForce;
    tensor_hao::TensorHao<double, 1> JForce;

    KanamoriInteractRealForce();
    KanamoriInteractRealForce(size_t L, size_t numberOfKana);
    KanamoriInteractRealForce(const KanamoriInteractRealForce& x);
    KanamoriInteractRealForce(KanamoriInteractRealForce&& x);
    ~KanamoriInteractRealForce();

    KanamoriInteractRealForce & operator  = (const KanamoriInteractRealForce& x);
    KanamoriInteractRealForce & operator  = (KanamoriInteractRealForce&& x);

    size_t getL() const;
    size_t getNumberOfKana() const;
    double getMemory() const;

 private:
    void copy_deep(const KanamoriInteractRealForce &x);
    void move_deep(KanamoriInteractRealForce &x);
};

#endif //AFQMCHUBBARDKANAMORI_KANAMORIINTERACTREALFORCE_H
