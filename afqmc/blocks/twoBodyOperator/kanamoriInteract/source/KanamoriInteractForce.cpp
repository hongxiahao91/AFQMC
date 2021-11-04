//
// Created by Hao Shi on 1/12/18.
//

#include "../include/KanamoriInteractForce.h"

using namespace std;

KanamoriInteractForce::KanamoriInteractForce()  { }

KanamoriInteractForce::KanamoriInteractForce(size_t L, size_t numberOfKana)
{
    UForce.resize(L);
    U1Force.resize(2*numberOfKana);
    U2JForce.resize(2*numberOfKana);
    JForce.resize(numberOfKana);
}

KanamoriInteractForce::KanamoriInteractForce(const KanamoriInteractForce &x) { copy_deep(x); }

KanamoriInteractForce::KanamoriInteractForce(KanamoriInteractForce &&x) { move_deep(x); }

KanamoriInteractForce::~KanamoriInteractForce() { }

KanamoriInteractForce &KanamoriInteractForce::operator=(const KanamoriInteractForce &x)  { copy_deep(x); return *this; }

KanamoriInteractForce &KanamoriInteractForce::operator=(KanamoriInteractForce &&x) { move_deep(x); return *this; }

size_t KanamoriInteractForce::getL() const { return UForce.size(); }

size_t KanamoriInteractForce::getNumberOfKana() const  { return JForce.size(); }

double KanamoriInteractForce::getMemory() const
{
    double mem(0.0);
    mem += UForce.getMemory();
    mem += U1Force.getMemory();
    mem += U2JForce.getMemory();
    mem += JForce.getMemory();
    return mem;
}

void KanamoriInteractForce::copy_deep(const KanamoriInteractForce &x)
{
    UForce = x.UForce;
    U1Force = x.U1Force;
    U2JForce = x.U2JForce;
    JForce = x.JForce;
}

void KanamoriInteractForce::move_deep(KanamoriInteractForce &x)
{
    UForce = move( x.UForce );
    U1Force = move( x.U1Force );
    U2JForce = move( x.U2JForce );
    JForce = move( x.JForce );
}
