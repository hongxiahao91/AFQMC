//
// Created by Hao Shi on 10/4/17.
//

#include "../include/hop2isSD2sOperation.h"

using namespace std;
using namespace tensor_hao;

Hop2isSD2sOperation::Hop2isSD2sOperation() { }

Hop2isSD2sOperation::~Hop2isSD2sOperation() { }

void Hop2isSD2sOperation::applyToRight(const Hop2is &oneBody, const SD2s &walker, SD2s &walkerNew) const
{
    checkAndResize(oneBody, walker, walkerNew);
    BL_NAME(gmm)( oneBody.matrix, walker.getWfUp(), walkerNew.wfUpRef() );
    BL_NAME(gmm)( oneBody.matrix, walker.getWfDn(), walkerNew.wfDnRef() );
    walkerNew.logwRef() = oneBody.logw + walker.getLogw();
}

void Hop2isSD2sOperation::applyToLeft(const Hop2is &oneBody, const SD2s &walker, SD2s &walkerNew) const
{
    checkAndResize(oneBody, walker, walkerNew);
    BL_NAME(gmm)( oneBody.matrix, walker.getWfUp(), walkerNew.wfUpRef(), 'C' );
    BL_NAME(gmm)( oneBody.matrix, walker.getWfDn(), walkerNew.wfDnRef(), 'C' );
    walkerNew.logwRef() = conj( oneBody.logw ) + walker.getLogw();
}

void Hop2isSD2sOperation::checkAndResize(const Hop2is &oneBody, const SD2s &walker, SD2s &walkerNew) const
{
    size_t L = walker.getL(); size_t Nup = walker.getNup(); size_t Ndn = walker.getNdn();
    if( oneBody.getL() !=  L ) {cout<<"Error!!! Hop2is size is not consistent with walker!"<<endl; exit(1); }
    if( walkerNew.getL() != L  ||  walkerNew.getNup() != Nup || walkerNew.getNdn() != Ndn ) walkerNew.resize( L, Nup, Ndn );
}