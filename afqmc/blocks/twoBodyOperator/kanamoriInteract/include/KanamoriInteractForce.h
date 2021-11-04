//
// Created by Hao Shi on 1/12/18.
//

#ifndef AFQMCHUBBARDKANAMORI_KANAMORIINTERACTFORCE_H
#define AFQMCHUBBARDKANAMORI_KANAMORIINTERACTFORCE_H

#include "../../../../../common/tensorHao/include/tensor_all.h"

class KanamoriInteractForce
{
 public:
    tensor_hao::TensorHao<double, 1> UForce;
    tensor_hao::TensorHao<double, 1> U1Force;
    tensor_hao::TensorHao<double, 1> U2JForce;
    tensor_hao::TensorHao<std::complex<double>, 1> JForce;

    KanamoriInteractForce();
    KanamoriInteractForce(size_t L, size_t numberOfKana);
    KanamoriInteractForce(const KanamoriInteractForce& x);
    KanamoriInteractForce(KanamoriInteractForce&& x);
    ~KanamoriInteractForce();

    KanamoriInteractForce & operator  = (const KanamoriInteractForce& x);
    KanamoriInteractForce & operator  = (KanamoriInteractForce&& x);

    size_t getL() const;
    size_t getNumberOfKana() const;
    double getMemory() const;

 private:
    void copy_deep(const KanamoriInteractForce &x);
    void move_deep(KanamoriInteractForce &x);
};

#endif //AFQMCHUBBARDKANAMORI_KANAMORIINTERACTFORCE_H
