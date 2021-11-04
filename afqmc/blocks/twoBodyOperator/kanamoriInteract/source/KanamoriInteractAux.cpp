//
// Created by Hao Shi on 1/12/18.
//

#include "../include/KanamoriInteractAux.h"

using namespace std;

KanamoriInteractAux::KanamoriInteractAux() { }

KanamoriInteractAux::KanamoriInteractAux(size_t L, size_t numberOfKana)
{
    UAux.resize(L);
    U1Aux.resize(2*numberOfKana);
    U2JAux.resize(2*numberOfKana);
    JAux.resize(numberOfKana);
}

KanamoriInteractAux::KanamoriInteractAux(const KanamoriInteractAux &x) { copy_deep(x); }

KanamoriInteractAux::KanamoriInteractAux(KanamoriInteractAux&&x) { move_deep(x); }

KanamoriInteractAux::~KanamoriInteractAux() { }

KanamoriInteractAux &KanamoriInteractAux::operator=(const KanamoriInteractAux &x)  { copy_deep(x); return *this; }

KanamoriInteractAux &KanamoriInteractAux::operator=(KanamoriInteractAux &&x) { move_deep(x); return *this; }

size_t KanamoriInteractAux::getL() const { return UAux.size(); }

size_t KanamoriInteractAux::getNumberOfKana() const  { return JAux.size(); }

double KanamoriInteractAux::getMemory() const
{
    double mem(0.0);
    mem += UAux.getMemory();
    mem += U1Aux.getMemory();
    mem += U2JAux.getMemory();
    mem += JAux.getMemory();
    return mem;
}

void KanamoriInteractAux::copy_deep(const KanamoriInteractAux &x)
{
    UAux = x.UAux;
    U1Aux = x.U1Aux;
    U2JAux = x.U2JAux;
    JAux = x.JAux;
}

void KanamoriInteractAux::move_deep(KanamoriInteractAux &x)
{
    UAux = move( x.UAux );
    U1Aux = move( x.U1Aux);
    U2JAux = move( x.U2JAux );
    JAux = move( x.JAux );
}