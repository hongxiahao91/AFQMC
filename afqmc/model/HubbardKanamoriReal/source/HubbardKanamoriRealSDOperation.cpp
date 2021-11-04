//
// Created by Hao Shi on 1/19/18.
//

#include "../include/HubbardKanamoriRealSDOperation.h"

using namespace std;
using namespace tensor_hao;

void fillWalkerRandomly(SD &walker, const HubbardKanamoriReal &model)
{
    size_t L = model.getL();
    size_t N = model.getN();
    walker.resize(2*L, N);

    walker.randomFill();
}

void fillWalkerFromModel(SD &walker, HubbardKanamoriReal &model)
{
    size_t L = model.getL();
    size_t N = model.getN();
    walker.resize(2*L, N);

    model.setKpEigenValueAndVector();
    const TensorHao< complex<double>, 2 > &KpEigenVector = model.getKpEigenVector();
    TensorHao<complex<double>,2> &wf = walker.wfRef();

    copy( KpEigenVector.data(), KpEigenVector.data()+2*L*N, wf.data() );
}