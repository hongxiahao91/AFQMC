//
// Created by Hao Shi on 1/12/18.
//

#ifndef AFQMCHUBBARDKANAMORI_KANAMORIINTERACTAUX_H
#define AFQMCHUBBARDKANAMORI_KANAMORIINTERACTAUX_H

#include "../../../../../common/tensorHao/include/tensor_all.h"

class KanamoriInteractAux
{
 public:
    tensor_hao::TensorHao<int, 1> UAux;
    tensor_hao::TensorHao<int, 1> U1Aux;
    tensor_hao::TensorHao<int, 1> U2JAux;
    tensor_hao::TensorHao<std::complex<double>, 1> JAux;

    KanamoriInteractAux();
    KanamoriInteractAux(size_t L, size_t numberOfKana);
    KanamoriInteractAux(const KanamoriInteractAux& x);
    KanamoriInteractAux(KanamoriInteractAux&& x);
    ~KanamoriInteractAux();

    KanamoriInteractAux & operator  = (const KanamoriInteractAux& x);
    KanamoriInteractAux & operator  = (KanamoriInteractAux&& x);

    size_t getL() const;
    size_t getNumberOfKana() const;
    double getMemory() const;

private:
    void copy_deep(const KanamoriInteractAux &x);
    void move_deep(KanamoriInteractAux &x);
};


#endif //AFQMCHUBBARDKANAMORI_KANAMORIINTERACTAUX_H
