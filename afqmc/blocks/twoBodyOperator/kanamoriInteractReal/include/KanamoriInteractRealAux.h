//
// Created by Hao Shi on 1/12/18.
//

#ifndef AFQMCHUBBARDKANAMORI_KANAMORIINTERACTREALAUX_H
#define AFQMCHUBBARDKANAMORI_KANAMORIINTERACTREALAUX_H

#include "../../../../../common/tensorHao/include/tensor_all.h"

class KanamoriInteractRealAux
{
 public:
    tensor_hao::TensorHao<int, 1> UAux;
    tensor_hao::TensorHao<int, 1> U1Aux;
    tensor_hao::TensorHao<int, 1> U2JAux;
    tensor_hao::TensorHao<double, 1> JAux;

    KanamoriInteractRealAux();
    KanamoriInteractRealAux(size_t L, size_t numberOfKana);
    KanamoriInteractRealAux(const KanamoriInteractRealAux& x);
    KanamoriInteractRealAux(KanamoriInteractRealAux&& x);
    ~KanamoriInteractRealAux();

    KanamoriInteractRealAux & operator  = (const KanamoriInteractRealAux& x);
    KanamoriInteractRealAux & operator  = (KanamoriInteractRealAux&& x);

    size_t getL() const;
    size_t getNumberOfKana() const;
    double getMemory() const;

private:
    void copy_deep(const KanamoriInteractRealAux &x);
    void move_deep(KanamoriInteractRealAux &x);
};


#endif //AFQMCHUBBARDKANAMORI_KANAMORIINTERACTREALAUX_H
