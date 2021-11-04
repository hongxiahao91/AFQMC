//
// Created by Hao Shi on 1/12/18.
//

#include "../include/KanamoriInteractRealForce.h"

using namespace std;

KanamoriInteractRealForce::KanamoriInteractRealForce()  { }

KanamoriInteractRealForce::KanamoriInteractRealForce(size_t L, size_t numberOfKana)
{
    UForce.resize(L);
    U1Force.resize(2*numberOfKana);
    U2JForce.resize(2*numberOfKana);
    JForce.resize(numberOfKana);
}

KanamoriInteractRealForce::KanamoriInteractRealForce(const KanamoriInteractRealForce &x) { copy_deep(x); }

KanamoriInteractRealForce::KanamoriInteractRealForce(KanamoriInteractRealForce &&x) { move_deep(x); }

KanamoriInteractRealForce::~KanamoriInteractRealForce() { }

KanamoriInteractRealForce &KanamoriInteractRealForce::operator=(const KanamoriInteractRealForce &x)  { copy_deep(x); return *this; }

KanamoriInteractRealForce &KanamoriInteractRealForce::operator=(KanamoriInteractRealForce &&x) { move_deep(x); return *this; }

size_t KanamoriInteractRealForce::getL() const { return UForce.size(); }

size_t KanamoriInteractRealForce::getNumberOfKana() const  { return JForce.size(); }

double KanamoriInteractRealForce::getMemory() const
{
    double mem(0.0);
    mem += UForce.getMemory();
    mem += U1Force.getMemory();
    mem += U2JForce.getMemory();
    mem += JForce.getMemory();
    return mem;
}

void KanamoriInteractRealForce::copy_deep(const KanamoriInteractRealForce &x)
{
    UForce = x.UForce;
    U1Force = x.U1Force;
    U2JForce = x.U2JForce;
    JForce = x.JForce;
}

void KanamoriInteractRealForce::move_deep(KanamoriInteractRealForce &x)
{
    UForce = move( x.UForce );
    U1Force = move( x.U1Force );
    U2JForce = move( x.U2JForce );
    JForce = move( x.JForce );
}
