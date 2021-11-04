//
// Created by Hao Shi on 1/12/18.
//

#include "../include/KanamoriInteractRealAux.h"

using namespace std;

KanamoriInteractRealAux::KanamoriInteractRealAux() { }

KanamoriInteractRealAux::KanamoriInteractRealAux(size_t L, size_t numberOfKana)
{
    UAux.resize(L);
    U1Aux.resize(2*numberOfKana);
    U2JAux.resize(2*numberOfKana);
    JAux.resize(numberOfKana);
}

KanamoriInteractRealAux::KanamoriInteractRealAux(const KanamoriInteractRealAux &x) { copy_deep(x); }

KanamoriInteractRealAux::KanamoriInteractRealAux(KanamoriInteractRealAux&&x) { move_deep(x); }

KanamoriInteractRealAux::~KanamoriInteractRealAux() { }

KanamoriInteractRealAux &KanamoriInteractRealAux::operator=(const KanamoriInteractRealAux &x)  { copy_deep(x); return *this; }

KanamoriInteractRealAux &KanamoriInteractRealAux::operator=(KanamoriInteractRealAux &&x) { move_deep(x); return *this; }

size_t KanamoriInteractRealAux::getL() const { return UAux.size(); }

size_t KanamoriInteractRealAux::getNumberOfKana() const  { return JAux.size(); }

double KanamoriInteractRealAux::getMemory() const
{
    double mem(0.0);
    mem += UAux.getMemory();
    mem += U1Aux.getMemory();
    mem += U2JAux.getMemory();
    mem += JAux.getMemory();
    return mem;
}

void KanamoriInteractRealAux::copy_deep(const KanamoriInteractRealAux &x)
{
    UAux = x.UAux;
    U1Aux = x.U1Aux;
    U2JAux = x.U2JAux;
    JAux = x.JAux;
}

void KanamoriInteractRealAux::move_deep(KanamoriInteractRealAux &x)
{
    UAux = move( x.UAux );
    U1Aux = move( x.U1Aux);
    U2JAux = move( x.U2JAux );
    JAux = move( x.JAux );
}