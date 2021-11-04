//
// Created by Hao Shi on 10/4/17.
//

#include "../include/MDCas2sSD2sOperation.h"

using namespace std;
using namespace tensor_hao;

MDCas2sSD2sOperation::MDCas2sSD2sOperation() { }

MDCas2sSD2sOperation::MDCas2sSD2sOperation(MDCas2s &walkerLeft_, SD2s &walkerRight_, size_t stablizeStep_)
{
    set(walkerLeft_, walkerRight_, stablizeStep_);
}

MDCas2sSD2sOperation::~MDCas2sSD2sOperation() { }

const MDCas2s *MDCas2sSD2sOperation::getWalkerLeft() const { return walkerLeft; }

const SD2s *MDCas2sSD2sOperation::getWalkerRight() const { return walkerRight; }

const TensorHao<complex<double>, 2> &MDCas2sSD2sOperation::getOvlpFullUp() const { return ovlpFullUp; }

const TensorHao<complex<double>, 2> &MDCas2sSD2sOperation::getOvlpFullDn() const { return ovlpFullDn; }

const vector<complex<double>> &MDCas2sSD2sOperation::getDetUp() const { return detUp; }

const vector<complex<double>> &MDCas2sSD2sOperation::getDetDn() const { return detDn; }

const vector<TensorHao<complex<double>, 2>> &MDCas2sSD2sOperation::getOvlpInvUp() const { return ovlpInvUp; }

const vector<TensorHao<complex<double>, 2>> &MDCas2sSD2sOperation::getOvlpInvDn() const { return ovlpInvDn; }

const vector<TensorHao<complex<double>, 2>> &MDCas2sSD2sOperation::getQUp() const { return QUp; }

const vector<TensorHao<complex<double>, 2>> &MDCas2sSD2sOperation::getQDn() const { return QDn; }

void MDCas2sSD2sOperation::set(MDCas2s &walkerLeft_, SD2s &walkerRight_, size_t stablizeStep_)
{
    walkerLeft  = &walkerLeft_;
    walkerRight = &walkerRight_;
    stablizeStep = stablizeStep_;

    size_t L = walkerLeft->getL(); size_t Nup = walkerLeft->getNup(); size_t Ndn = walkerLeft->getNdn();
    size_t NupMax = walkerLeft->getNupMax(); size_t NdnMax = walkerLeft->getNdnMax();

    if( L != walkerRight->getL() || Nup != walkerRight->getNup() || Ndn != walkerRight->getNdn() )
    {
        cout<<"Error!!! Find Inconsistency between walkerLeft and walkerRight!"<<endl;
        exit(1);
    }

    const TensorHao<complex<double>,2> &wfCasUp = walkerLeft->generateWfCasUp();
    const TensorHao<complex<double>,2> &wfCasDn = walkerLeft->generateWfCasDn();
    ovlpFullUp.resize( NupMax, Nup ); BL_NAME(gmm)( wfCasUp, walkerRight->getWfUp(), ovlpFullUp, 'C' );
    ovlpFullDn.resize( NdnMax, Ndn ); BL_NAME(gmm)( wfCasDn, walkerRight->getWfDn(), ovlpFullDn, 'C' );

    setDetOvlpInvQForUpSpin();
    setDetOvlpInvQForDnSpin();

    reSet();
}

void MDCas2sSD2sOperation::reSet()
{
    logOverlapIsCalculated = false;
    thetaIsCalculated = false;
}

complex<double> MDCas2sSD2sOperation::returnLogOverlap()
{
    if(logOverlapIsCalculated) return logOverlap;

    size_t MDSize = walkerLeft->getMDSize();
    const vector< complex<double> > &coe = walkerLeft->getCoe();
    const vector< size_t > &coeLinkUp = walkerLeft->getCoeLinkUp();
    const vector< size_t > &coeLinkDn = walkerLeft->getCoeLinkDn();
    complex<double> overlap(0,0);
    for(size_t i = 0; i < MDSize; ++i)
    {
        overlap += conj( coe[i] ) * detUp[ coeLinkUp[i] ] * detDn[ coeLinkDn[i] ];
    }

    logOverlap = log(overlap) + conj( walkerLeft->getLogw() ) + walkerRight->getLogw();

    logOverlapIsCalculated = true;

    return logOverlap;
}

const vector< TensorHao<complex<double>, 2> > &MDCas2sSD2sOperation::returnThetaUp()
{
    calculateTheta();
    return thetaUp;
}

const vector< TensorHao<complex<double>, 2> > &MDCas2sSD2sOperation::returnThetaDn()
{
    calculateTheta();
    return thetaDn;
}

double MDCas2sSD2sOperation::getMemory() const
{
    double mem(0.0);

    mem += 8.0+8.0+8.0;
    mem += ovlpFullUp.getMemory() + ovlpFullDn.getMemory();
    mem += 16.0*detUp.size() + 16.0*detDn.size();

    size_t MDUpSize =walkerLeft->getMDUpSize();
    size_t MDDnSize =walkerLeft->getMDDnSize();

    for(size_t i = 0; i < MDUpSize; ++i)
    {
        mem += ovlpInvUp[i].getMemory();
        mem += QUp[i].getMemory();
    }

    for(size_t i = 0; i < MDDnSize; ++i)
    {
        mem += ovlpInvDn[i].getMemory();
        mem += QDn[i].getMemory();
    }

    mem += 16.0+1.0;

    for(size_t i = 0; i < MDUpSize; ++i) mem += thetaUp[i].getMemory();
    for(size_t i = 0; i < MDDnSize; ++i) mem += thetaDn[i].getMemory();
    mem += 1.0;

    return mem;
}

MDCas2sSD2sOperation::MDCas2sSD2sOperation(const MDCas2sSD2sOperation &x) { }

MDCas2sSD2sOperation &MDCas2sSD2sOperation::operator=(const MDCas2sSD2sOperation &x) { return *this; }

void MDCas2sSD2sOperation::setDetOvlpInvQForUpSpin()
{
    size_t MDUpSize = walkerLeft->getMDUpSize(); size_t Nup = walkerLeft->getNup();
    const vector< vector<size_t> > &occupancyUp = walkerLeft->getOccupancyUp();
    const vector<int> &isParentUp = walkerLeft->getIsParentUp();
    const vector<size_t> &treeUp = walkerLeft->getTreeUp();
    const vector<vector<size_t>> &diffOrbitalUp = walkerLeft->getDiffOrbitalUp();

    size_t diffNumber, currentParent, currentIndex, parentUpdateNumber;
    TensorHao<complex<double>,2> VUp, OUp, QUUp, AInvUUp;
    LUDecomp<complex<double> > luDecomp;

    detUp.resize(MDUpSize); ovlpInvUp.resize(MDUpSize); QUp.resize(MDUpSize);

    //Det, ovlpInv for 1st step
    fillCasMatrixByRow(ovlpInvUp[0], ovlpFullUp, occupancyUp[0]);
    luDecomp = BL_NAME(LUconstruct)( ovlpInvUp[0] );
    detUp[0] = determinant(luDecomp);
    ovlpInvUp[0] = BL_NAME(inverse)( luDecomp );
    currentParent=0;
    parentUpdateNumber=0;

    for(size_t i = 1; i < MDUpSize; ++i)
    {
        currentIndex = treeUp[i];

        diffNumber = diffOrbitalUp[currentIndex].size();

        getRowDifference(VUp, ovlpFullUp, occupancyUp[currentIndex], occupancyUp[currentParent], diffOrbitalUp[currentIndex]);

        QUp[currentIndex].resize(diffNumber, Nup);
        BL_NAME(gmm)( VUp, ovlpInvUp[currentParent], QUp[currentIndex] );

        fillCasMatrixByCol( OUp, QUp[currentIndex], diffOrbitalUp[currentIndex] );
        for(size_t j = 0; j < diffNumber; ++j) OUp(j,j) += 1.0;
        luDecomp = BL_NAME(LUconstruct)( OUp );
        QUp[currentIndex] =  BL_NAME(solve_lineq)( luDecomp, QUp[currentIndex] );

        fillCasMatrixByCol( QUUp, QUp[currentIndex], diffOrbitalUp[currentIndex] );
        QUUp *= complex<double>(-1.0, 0.0);
        for(size_t j = 0; j < diffNumber ; ++j) { QUUp(j,j) += 1.0; }
        detUp[currentIndex] = detUp[currentParent] / BL_NAME(determinant)( QUUp );

        if( isParentUp[currentIndex] == 1 )
        {
            if( parentUpdateNumber < stablizeStep )
            {
                ovlpInvUp[currentIndex].resize(Nup, Nup);
                fillCasMatrixByCol( AInvUUp, ovlpInvUp[currentParent], diffOrbitalUp[currentIndex] );
                BL_NAME(gmm)(AInvUUp, QUp[currentIndex], ovlpInvUp[currentIndex], 'N', 'N', -1.0, 0.0);
                ovlpInvUp[currentIndex] += ovlpInvUp[currentParent];

                parentUpdateNumber++;
            }
            else
            {
                fillCasMatrixByRow(ovlpInvUp[currentIndex], ovlpFullUp, occupancyUp[currentIndex]);
                luDecomp = BL_NAME(LUconstruct)( ovlpInvUp[currentIndex] );
                detUp[currentIndex] = determinant(luDecomp);
                ovlpInvUp[currentIndex] = BL_NAME(inverse)( luDecomp );

                parentUpdateNumber=0;
            }

            currentParent = currentIndex;
        }
    }
}

void MDCas2sSD2sOperation::setDetOvlpInvQForDnSpin()
{
    size_t MDDnSize = walkerLeft->getMDDnSize(); size_t Ndn = walkerLeft->getNdn();
    const vector< vector<size_t> > &occupancyDn = walkerLeft->getOccupancyDn();
    const vector<int> &isParentDn = walkerLeft->getIsParentDn();
    const vector<size_t> &treeDn = walkerLeft->getTreeDn();
    const vector<vector<size_t>> &diffOrbitalDn = walkerLeft->getDiffOrbitalDn();

    size_t diffNumber, currentParent, currentIndex, parentUpdateNumber;
    TensorHao<complex<double>,2> VDn, ODn, QUDn, AInvUDn;
    LUDecomp<complex<double> > luDecomp;

    detDn.resize(MDDnSize); ovlpInvDn.resize(MDDnSize); QDn.resize(MDDnSize);

    //Det, ovlpInv for 1st step
    fillCasMatrixByRow(ovlpInvDn[0], ovlpFullDn, occupancyDn[0]);
    luDecomp = BL_NAME(LUconstruct)( ovlpInvDn[0] );
    detDn[0] = determinant(luDecomp);
    ovlpInvDn[0] = BL_NAME(inverse)( luDecomp );
    currentParent=0;
    parentUpdateNumber=0;

    for(size_t i = 1; i < MDDnSize; ++i)
    {
        currentIndex = treeDn[i];

        diffNumber = diffOrbitalDn[currentIndex].size();

        getRowDifference(VDn, ovlpFullDn, occupancyDn[currentIndex], occupancyDn[currentParent], diffOrbitalDn[currentIndex]);

        QDn[currentIndex].resize(diffNumber, Ndn);
        BL_NAME(gmm)( VDn, ovlpInvDn[currentParent], QDn[currentIndex] );

        fillCasMatrixByCol( ODn, QDn[currentIndex], diffOrbitalDn[currentIndex] );
        for(size_t j = 0; j < diffNumber; ++j) ODn(j,j) += 1.0;
        luDecomp = BL_NAME(LUconstruct)( ODn );
        QDn[currentIndex] =  BL_NAME(solve_lineq)( luDecomp, QDn[currentIndex] );

        fillCasMatrixByCol( QUDn, QDn[currentIndex], diffOrbitalDn[currentIndex] );
        QUDn *= complex<double>(-1.0, 0.0);
        for(size_t j = 0; j < diffNumber ; ++j) { QUDn(j,j) += 1.0; }
        detDn[currentIndex] = detDn[currentParent] / BL_NAME(determinant)( QUDn );

        if( isParentDn[currentIndex] == 1 )
        {
            if( parentUpdateNumber < stablizeStep )
            {
                ovlpInvDn[currentIndex].resize(Ndn, Ndn);
                fillCasMatrixByCol( AInvUDn, ovlpInvDn[currentParent], diffOrbitalDn[currentIndex] );
                BL_NAME(gmm)(AInvUDn, QDn[currentIndex], ovlpInvDn[currentIndex], 'N', 'N', -1.0, 0.0);
                ovlpInvDn[currentIndex] += ovlpInvDn[currentParent];

                parentUpdateNumber++;
            }
            else
            {
                fillCasMatrixByRow(ovlpInvDn[currentIndex], ovlpFullDn, occupancyDn[currentIndex]);
                luDecomp = BL_NAME(LUconstruct)( ovlpInvDn[currentIndex] );
                detDn[currentIndex] = determinant(luDecomp);
                ovlpInvDn[currentIndex] = BL_NAME(inverse)( luDecomp );

                parentUpdateNumber=0;
            }

            currentParent = currentIndex;
        }
    }
}

void MDCas2sSD2sOperation::calculateTheta()
{
    if(thetaIsCalculated) return;

    calculateThetaUp();
    calculateThetaDn();

    thetaIsCalculated = true;
}

void MDCas2sSD2sOperation::calculateThetaUp()
{
    size_t MDUpSize = walkerLeft->getMDUpSize(); size_t L = walkerLeft->getL(); size_t Nup = walkerLeft->getNup();
    const vector<int> &isParentUp = walkerLeft->getIsParentUp();
    const vector<size_t> &treeUp = walkerLeft->getTreeUp();
    const vector<vector<size_t>> &diffOrbitalUp = walkerLeft->getDiffOrbitalUp();
    TensorHao<complex<double>,2> thetaUUp;

    size_t currentParent, currentIndex, parentUpdateNumber;
    thetaUp.resize(MDUpSize);

    thetaUp[0].resize(L, Nup);
    BL_NAME(gmm)( walkerRight->getWfUp(), ovlpInvUp[0], thetaUp[0] );

    currentParent = 0; parentUpdateNumber = 0;
    for(size_t i = 1; i < MDUpSize; ++i)
    {
        currentIndex = treeUp[i];

        fillCasMatrixByCol( thetaUUp, thetaUp[currentParent], diffOrbitalUp[currentIndex] );
        thetaUp[currentIndex] = thetaUp[currentParent];
        BL_NAME(gmm)( thetaUUp, QUp[currentIndex], thetaUp[currentIndex], 'N', 'N', -1.0, 1.0 );

        if( isParentUp[currentIndex] == 1 )
        {
            if( parentUpdateNumber < stablizeStep )
            {
                parentUpdateNumber++;
            }
            else
            {
                thetaUp[currentIndex].resize(L, Nup);
                BL_NAME(gmm)( walkerRight->getWfUp(), ovlpInvUp[currentIndex], thetaUp[currentIndex] );

                parentUpdateNumber=0;
            }

            currentParent = currentIndex;
        }
    }
}

void MDCas2sSD2sOperation::calculateThetaDn()
{
    size_t MDDnSize = walkerLeft->getMDDnSize(); size_t L = walkerLeft->getL(); size_t Ndn = walkerLeft->getNdn();
    const vector<int> &isParentDn = walkerLeft->getIsParentDn();
    const vector<size_t> &treeDn = walkerLeft->getTreeDn();
    const vector<vector<size_t>> &diffOrbitalDn = walkerLeft->getDiffOrbitalDn();
    TensorHao<complex<double>,2> thetaUDn;

    size_t currentParent, currentIndex, parentUpdateNumber;
    thetaDn.resize(MDDnSize);

    thetaDn[0].resize(L, Ndn);
    BL_NAME(gmm)( walkerRight->getWfDn(), ovlpInvDn[0], thetaDn[0] );

    currentParent = 0; parentUpdateNumber = 0;
    for(size_t i = 1; i < MDDnSize; ++i)
    {
        currentIndex = treeDn[i];

        fillCasMatrixByCol( thetaUDn, thetaDn[currentParent], diffOrbitalDn[currentIndex] );
        thetaDn[currentIndex] = thetaDn[currentParent];
        BL_NAME(gmm)( thetaUDn, QDn[currentIndex], thetaDn[currentIndex], 'N', 'N', -1.0, 1.0 );

        if( isParentDn[currentIndex] == 1 )
        {
            if( parentUpdateNumber < stablizeStep )
            {
                parentUpdateNumber++;
            }
            else
            {
                thetaDn[currentIndex].resize(L, Ndn);
                BL_NAME(gmm)( walkerRight->getWfDn(), ovlpInvDn[currentIndex], thetaDn[currentIndex] );

                parentUpdateNumber=0;
            }

            currentParent = currentIndex;
        }
    }
}

void setWalkerFromPhiT(vector<SD2s> &walker, vector<bool> &walkerIsAlive, MDCas2s &phiT, double noise)
{
    size_t walkerSizePerThread = walker.size();
    size_t walkerSize = walkerSizePerThread * MPISize();
    if( walkerSize < 1 ) { cout<<"Error!!! Total walkerSize is smaller than 1:  "<<walkerSize<<endl; exit(1); }

    size_t MDSize = phiT.getMDSize();
    if( MDSize>walkerSize ) { if( MPIRank()==0 ) cout<<"Warning!!! Total walkerSize is smaller than MDSize!"<<endl; }

    size_t L = phiT.getL(); size_t Nup = phiT.getNup(); size_t Ndn = phiT.getNdn();
    const vector< complex<double> > &coe = phiT.getCoe();
    const vector< size_t > &coeLinkUp = phiT.getCoeLinkUp();
    const vector< size_t > &coeLinkDn = phiT.getCoeLinkDn();
    const vector< vector<size_t> > &occupancyUp = phiT.getOccupancyUp();
    const vector< vector<size_t> > &occupancyDn = phiT.getOccupancyDn();
    const TensorHao<complex<double>,2> &wfUp = phiT.generateWfCasUp();
    const TensorHao<complex<double>,2> &wfDn = phiT.generateWfCasDn();

    size_t linkIndex, mdIndex, col;
    for(size_t i = 0; i < walkerSizePerThread; ++i)
    {
        walker[i].resize(L, Nup, Ndn);

        mdIndex = i + MPIRank() * walkerSizePerThread;

        if( mdIndex<MDSize )
        {
            walker[i].logwRef() = log( coe[ mdIndex ] );

            linkIndex = coeLinkUp[mdIndex];
            TensorHao<complex<double>, 2> &walkerWfUp = walker[i].wfUpRef();
            for(size_t j = 0; j < Nup; ++j)
            {
                col = occupancyUp[linkIndex][j];
                for(size_t k = 0; k < L; ++k)
                {
                    walkerWfUp(k,j) = wfUp(k, col) + noise * complex<double>( uniformHao(), uniformHao() );
                }
            }

            linkIndex = coeLinkDn[mdIndex];
            TensorHao<complex<double>, 2> &walkerWfDn = walker[i].wfDnRef();
            for(size_t j = 0; j < Ndn; ++j)
            {
                col = occupancyDn[linkIndex][j];
                for(size_t k = 0; k < L; ++k)
                {
                    walkerWfDn(k,j) = wfDn(k, col) + noise * complex<double>( uniformHao(), uniformHao() );
                }
            }

            walkerIsAlive[i] = true;
        }
        else
        {
            walkerIsAlive[i] = false;
        }
    }
}