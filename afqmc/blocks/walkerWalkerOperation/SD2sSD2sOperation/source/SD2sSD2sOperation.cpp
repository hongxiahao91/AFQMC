//
// Created by boruoshihao on 11/14/18.
//

#include "../include/SD2sSD2sOperation.h"

using namespace std;
using namespace tensor_hao;

SD2sSD2sOperation::SD2sSD2sOperation(): walkerLeft(nullptr), walkerRight(nullptr)
{
    reSet();
}

SD2sSD2sOperation::SD2sSD2sOperation(const SD2s &walkerLeft_, const SD2s &walkerRight_)
{
    set(walkerLeft_, walkerRight_);
}

SD2sSD2sOperation::~SD2sSD2sOperation() { }

SD2sSD2sOperationState SD2sSD2sOperation::getState() const { return state; }

const SD2s *SD2sSD2sOperation::getWalkerLeft() const { return walkerLeft; }

const SD2s *SD2sSD2sOperation::getWalkerRight() const { return walkerRight; }

void SD2sSD2sOperation::set(const SD2s &walkerLeft_, const SD2s &walkerRight_)
{
    walkerLeft  = &walkerLeft_;
    walkerRight = &walkerRight_;

    size_t L = walkerLeft->getL(); size_t Nup = walkerLeft->getNup(); size_t Ndn = walkerLeft->getNdn();
    if( L != walkerRight->getL() || Nup != walkerRight->getNup() || Ndn != walkerRight->getNdn() )
    {
        cout<<"Error!!! Find Inconsistency between walkerLeft and walkerRight!"<<endl;
        exit(1);
    }

    reSet();
}

void SD2sSD2sOperation::reSet()
{
    state = SD2sSD2sOperationState::VOID;
    logOverlapIsCalculated = false;
    greenMatrixUpIsCalculated = false;
    greenMatrixDnIsCalculated = false;
    greenDiagonalUpIsCalculated = false;
    greenDiagonalDnIsCalculated = false;
}

const LUDecomp<complex<double>> &SD2sSD2sOperation::returnLUOverlapUp()
{
    calculateLUOverlap();
    return LUOverlapUp;
}

const LUDecomp<complex<double>> &SD2sSD2sOperation::returnLUOverlapDn()
{
    calculateLUOverlap();
    return LUOverlapDn;
}

const TensorHao<complex<double>, 2> &SD2sSD2sOperation::returnThetaUp_T()
{
    calculateLUOverlap();
    calculateTheta_T();
    return thetaUp_T;
}

const TensorHao<complex<double>, 2> &SD2sSD2sOperation::returnThetaDn_T()
{
    calculateLUOverlap();
    calculateTheta_T();
    return thetaDn_T;
}

complex<double> SD2sSD2sOperation::returnLogOverlap()
{
    if(logOverlapIsCalculated) return logOverlap;
    
    calculateLUOverlap();
    logOverlap =conj(walkerLeft->getLogw())+walkerRight->getLogw()+logDeterminant(LUOverlapUp)+logDeterminant(LUOverlapDn);
    
    logOverlapIsCalculated = true;
    
    return logOverlap;
}

const TensorHao<complex<double>, 2> &SD2sSD2sOperation::returnGreenMatrixUp()
{
    if(greenMatrixUpIsCalculated) return greenMatrixUp;
    
    calculateLUOverlap();
    calculateTheta_T();

    size_t L = walkerLeft->getL();
    greenMatrixUp.resize(L,L);
    BL_NAME(gmm)( conj( walkerLeft->getWfUp() ), thetaUp_T, greenMatrixUp );

    greenMatrixUpIsCalculated = true;
    
    return greenMatrixUp;
}

const TensorHao<complex<double>, 2> &SD2sSD2sOperation::returnGreenMatrixDn()
{
    if(greenMatrixDnIsCalculated) return greenMatrixDn;
    
    calculateLUOverlap();
    calculateTheta_T();

    size_t L = walkerLeft->getL();
    greenMatrixDn.resize(L,L);
    BL_NAME(gmm)( conj( walkerLeft->getWfDn() ), thetaDn_T, greenMatrixDn );

    greenMatrixDnIsCalculated = true;
    
    return greenMatrixDn;
}

const TensorHao<complex<double>, 1> &SD2sSD2sOperation::returnGreenDiagonalUp()
{
    if(greenDiagonalUpIsCalculated) return greenDiagonalUp;
    
    calculateLUOverlap();
    calculateTheta_T();

    size_t L = walkerLeft->getL(); size_t Nup = walkerLeft->getNup();
    const TensorHao<complex<double>, 2> &wfLeftUp = walkerLeft->getWfUp();
    greenDiagonalUp.resize(L); greenDiagonalUp = complex<double>(0,0);
    for(size_t j = 0; j < Nup; ++j)
    {
        for(size_t i = 0; i < L; ++i)
        {
            greenDiagonalUp(i) += conj( wfLeftUp(i, j) ) * thetaUp_T(j,i);
        }
    }
    
    greenDiagonalUpIsCalculated= true;
    
    return greenDiagonalUp;
}

const TensorHao<complex<double>, 1> &SD2sSD2sOperation::returnGreenDiagonalDn()
{
    if(greenDiagonalDnIsCalculated) return greenDiagonalDn;
    
    calculateLUOverlap();
    calculateTheta_T();

    size_t L = walkerLeft->getL(); size_t Ndn = walkerLeft->getNdn();
    const TensorHao<complex<double>, 2> &wfLeftDn = walkerLeft->getWfDn();
    greenDiagonalDn.resize(L); greenDiagonalDn = complex<double>(0,0);
    for(size_t j = 0; j < Ndn; ++j)
    {
        for(size_t i = 0; i < L; ++i)
        {
            greenDiagonalDn(i) += conj( wfLeftDn(i, j) ) * thetaDn_T(j,i);
        }
    }
    
    greenDiagonalDnIsCalculated = true;
    
    return greenDiagonalDn;
}

double SD2sSD2sOperation::getMemory() const
{
    return 8.0*2+LUOverlapUp.A.getMemory()+LUOverlapUp.ipiv.getMemory()
                +LUOverlapDn.A.getMemory()+LUOverlapDn.ipiv.getMemory()
                +thetaUp_T.getMemory()+thetaDn_T.getMemory()
                +16.0+1.0+greenMatrixUp.getMemory()+1.0+greenMatrixDn.getMemory()+1.0
                +greenDiagonalUp.getMemory()+1.0+greenDiagonalDn.getMemory()+1.0;
}

SD2sSD2sOperation::SD2sSD2sOperation(const SD2sSD2sOperation &x) { }

SD2sSD2sOperation &SD2sSD2sOperation::operator=(const SD2sSD2sOperation &x) { return *this; }

void SD2sSD2sOperation::calculateLUOverlap()
{
    if( state >= SD2sSD2sOperationState ::LUOVERLAP ) return;

    size_t Nup = walkerLeft->getNup(); size_t Ndn = walkerLeft->getNdn();

    TensorHao<complex<double>,2> overlapMatrixUp(Nup, Nup), overlapMatrixDn(Ndn, Ndn);
    BL_NAME(gmm)( walkerLeft->getWfUp(), walkerRight->getWfUp(), overlapMatrixUp, 'C' );
    BL_NAME(gmm)( walkerLeft->getWfDn(), walkerRight->getWfDn(), overlapMatrixDn, 'C' );

    LUOverlapUp = BL_NAME(LUconstruct)( move(overlapMatrixUp) );
    LUOverlapDn = BL_NAME(LUconstruct)( move(overlapMatrixDn) );

    state = SD2sSD2sOperationState ::LUOVERLAP;
}

void SD2sSD2sOperation::calculateTheta_T()
{
    if( state >= SD2sSD2sOperationState ::THETA_T ) return;

    thetaUp_T =  BL_NAME(solve_lineq)( LUOverlapUp, trans(walkerRight->getWfUp()), 'T' );
    thetaDn_T =  BL_NAME(solve_lineq)( LUOverlapDn, trans(walkerRight->getWfDn()), 'T' );

    state = SD2sSD2sOperationState ::THETA_T;
}

void setWalkerFromPhiT(vector<SD2s> &walker, vector<bool> &walkerIsAlive, const SD2s &phiT)
{
    int walkerSizePerThread = walker.size();
    int walkerSize = walkerSizePerThread * MPISize();
    if( walkerSize < 1 ) { cout<<"Error!!! Total walkerSize is smaller than 1:  "<<walkerSize<<endl; exit(1); }

    for(int i = 0; i < walkerSizePerThread; ++i)
    {
        walker[i] = phiT;
        walkerIsAlive[i] = true;
    }
}