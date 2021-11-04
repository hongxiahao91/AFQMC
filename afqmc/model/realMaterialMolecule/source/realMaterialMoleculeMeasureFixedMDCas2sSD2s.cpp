//
// Created by Hao Shi on 8/31/17.
//

#include "../include/realMaterialMoleculeMeasureFixedMDCas2sSD2s.h"
#include "../../../utilities/manipulateMCData/include/writeThreadSum.h"

using namespace std;
using namespace tensor_hao;

RealMaterialMoleculeMeasureFixedMDCas2sSD2s::RealMaterialMoleculeMeasureFixedMDCas2sSD2s()
{
    initModelWalkerNullptr();
    reSet();
}

RealMaterialMoleculeMeasureFixedMDCas2sSD2s::RealMaterialMoleculeMeasureFixedMDCas2sSD2s(const RealMaterialMolecule &realMaterialMolecule_, MDCas2s &walkerLeft_)
{
    setModelWalker(realMaterialMolecule_, walkerLeft_);
    reSet();
}

RealMaterialMoleculeMeasureFixedMDCas2sSD2s::~RealMaterialMoleculeMeasureFixedMDCas2sSD2s()
{

}

void RealMaterialMoleculeMeasureFixedMDCas2sSD2s::initModelWalkerNullptr()
{
    realMaterialMolecule = nullptr;
    walkerLeft = nullptr;
}

void RealMaterialMoleculeMeasureFixedMDCas2sSD2s::setModelWalker(const RealMaterialMolecule &realMaterialMolecule_, MDCas2s &walkerLeft_)
{
    if( realMaterialMolecule_.getL() != walkerLeft_.getL() ) {cout<<"Model L does not consistent with walker L!"<<endl; exit(1);}
    if( realMaterialMolecule_.getNup() != walkerLeft_.getNup() ) {cout<<"Model Nup does not consistent with walker Nup!"<<endl; exit(1);}
    if( realMaterialMolecule_.getNdn() != walkerLeft_.getNdn() ) {cout<<"Model Ndn does not consistent with walker Ndn!"<<endl; exit(1);}

    realMaterialMolecule = &realMaterialMolecule_;
    walkerLeft = &walkerLeft_;

    wfCasUpConj = conj( walkerLeft->generateWfCasUp() );
    wfCasDnConj = conj( walkerLeft->generateWfCasDn() );

    initWfDaggerT_T();
    initWfDaggerCholeskyVecs_T();
}

void RealMaterialMoleculeMeasureFixedMDCas2sSD2s::reSet()
{
    complex<double> zero(0,0);

    den = zero;
    TNum = zero;
    choleskyBgNum = zero;
    choleskyExNum = zero;
    HNum = zero;
}

complex<double> RealMaterialMoleculeMeasureFixedMDCas2sSD2s::returnEnergy()
{
    complex<double> Htot   = MPISum(HNum);
    complex<double> denTot = MPISum(den);
    complex<double> energy;
    if( MPIRank() == 0 ) energy = Htot/denTot;
    MPIBcast(energy);
    return energy;
}

TensorHao<double,1> RealMaterialMoleculeMeasureFixedMDCas2sSD2s::returnCholeskyBg()
{
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();

    TensorHao<complex<double>, 1> choleskyBgTot( choleskyNumber );
    MPISum( choleskyNumber, choleskyBgNum.data(), choleskyBgTot.data() );

    complex<double> denTot = MPISum(den);

    TensorHao<double,1> choleskyBg( choleskyNumber );
    for(size_t i = 0; i < choleskyNumber; ++i) choleskyBg(i) = ( choleskyBgTot(i)/denTot ).real();
    MPIBcast(choleskyBg);

    return choleskyBg;
}

void RealMaterialMoleculeMeasureFixedMDCas2sSD2s::addMeasurement(MDCas2sSD2sOperation &mdCas2sSD2sOperation, complex<double> denIncrement)
{
    checkWalkerLeft(mdCas2sSD2sOperation);

    den += denIncrement;

    const SD2s *walkerRight = mdCas2sSD2sOperation.getWalkerRight();
    complex<double> logOverlap = mdCas2sSD2sOperation.returnLogOverlap();
    complex<double> numIncrement = denIncrement * exp( conj( walkerLeft->getLogw() ) + walkerRight->getLogw() - logOverlap );

    //add: ( Local Estimator ) * numIncrement* conj( coe[i] )*detUp[i]*detDn[i]
    addEnergy(mdCas2sSD2sOperation, numIncrement);
}

CholeskyRealForce RealMaterialMoleculeMeasureFixedMDCas2sSD2s::getForce(const CholeskyReal &choleskyReal, MDCas2sSD2sOperation &mdCas2sSD2sOperation, double cap)
{
    checkWalkerLeft(mdCas2sSD2sOperation);

    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();

    const SD2s *walkerRight = mdCas2sSD2sOperation.getWalkerRight();
    complex<double> logOverlap = mdCas2sSD2sOperation.returnLogOverlap();
    complex<double> wholeFactor = exp( conj( walkerLeft->getLogw() ) + walkerRight->getLogw() - logOverlap );

    size_t MDSize = walkerLeft->getMDSize();
    size_t MDUpSize = walkerLeft->getMDUpSize();
    size_t MDDnSize = walkerLeft->getMDDnSize();

    const vector< complex<double> > &coe = walkerLeft->getCoe();
    const vector< size_t > &coeLinkUp = walkerLeft->getCoeLinkUp();
    const vector< size_t > &coeLinkDn = walkerLeft->getCoeLinkDn();

    const vector< complex<double> > &detUp = mdCas2sSD2sOperation.getDetUp();
    const vector< complex<double> > &detDn = mdCas2sSD2sOperation.getDetDn();
    const vector< TensorHao<complex<double>, 2> > &thetaUpVec = mdCas2sSD2sOperation.returnThetaUp();
    const vector< TensorHao<complex<double>, 2> > &thetaDnVec = mdCas2sSD2sOperation.returnThetaDn();

    vector<TensorHao<complex<double>, 1>> choleskyBgUp(MDUpSize), choleskyBgDn(MDDnSize);
    for(size_t i = 0; i <MDUpSize; ++i) choleskyBgUp[i] = calculateCholeskyBgUp(i, thetaUpVec[i]);
    for(size_t i = 0; i <MDDnSize; ++i) choleskyBgDn[i] = calculateCholeskyBgDn(i, thetaDnVec[i]);

    complex<double> scaleFactor;
    TensorHao<complex<double>, 1> choleskyBgOne(choleskyNumber);
    TensorHao<complex<double>, 1> choleskyBg(choleskyNumber); choleskyBg=complex<double>(0.0, 0.0);
    size_t indexUp, indexDn;
    for(size_t i = 0; i < MDSize; ++i)
    {
        indexUp = coeLinkUp[i]; indexDn = coeLinkDn[i];
        choleskyBgOne = choleskyBgUp[indexUp] + choleskyBgDn[indexDn];

        scaleFactor = wholeFactor * conj( coe[i] ) * detUp[indexUp] * detDn[indexDn];
        choleskyBg += ( choleskyBgOne * scaleFactor );
    }

    complex<double> sqrtMinusDt = choleskyReal.getSqrtMinusDt();
    const TensorHao<double,1> & currentBg = realMaterialMolecule->getCholeskyBg();

    CholeskyRealForce force(choleskyNumber); complex<double> oneForce;
    for(size_t i = 0; i < choleskyNumber; ++i)
    {
        oneForce = (choleskyBg(i)-currentBg(i)) * sqrtMinusDt;

        if( abs(oneForce) > cap ) force(i) = oneForce*cap/abs(oneForce);
        else force(i) = oneForce;
    }

    return force;
}

void RealMaterialMoleculeMeasureFixedMDCas2sSD2s::write() const
{
    writeThreadSum(den, "den.dat", ios::app);
    writeThreadSum(TNum, "TNum.dat", ios::app);
    writeThreadSum(choleskyBgNum.size(), choleskyBgNum.data(), "choleskyBgNum.dat", ios::app);
    writeThreadSum(choleskyExNum.size(), choleskyExNum.data(), "choleskyExNum.dat", ios::app);
    writeThreadSum(HNum, "HNum.dat", ios::app);
}

double RealMaterialMoleculeMeasureFixedMDCas2sSD2s::getMemory() const
{
    double mem(0.0);
    mem += 8.0*2;
    mem += wfCasUpConj.getMemory() + wfCasDnConj.getMemory();
    mem += 16.0*3; //den, TNum, HNum
    mem += choleskyBgNum.getMemory()+choleskyExNum.getMemory();
    mem += wfUpDaggerT_T.getMemory() + wfDnDaggerT_T.getMemory();
    mem += wfUpDaggerCholeskyVecs_T.getMemory() + wfDnDaggerCholeskyVecs_T.getMemory();
    return mem;
}

RealMaterialMoleculeMeasureFixedMDCas2sSD2s::RealMaterialMoleculeMeasureFixedMDCas2sSD2s(const RealMaterialMoleculeMeasureFixedMDCas2sSD2s &x)
{

}

RealMaterialMoleculeMeasureFixedMDCas2sSD2s & RealMaterialMoleculeMeasureFixedMDCas2sSD2s::operator=(const RealMaterialMoleculeMeasureFixedMDCas2sSD2s &x)
{
    return *this;
}

void RealMaterialMoleculeMeasureFixedMDCas2sSD2s::initWfDaggerT_T()
{
    size_t L = walkerLeft->getL(); size_t NupMax = walkerLeft->getNupMax(); size_t NdnMax = walkerLeft->getNdnMax();

    wfUpDaggerT_T.resize(L, NupMax); wfDnDaggerT_T.resize(L, NdnMax);

    const TensorHao<double,2> & t = realMaterialMolecule->getT();
    TensorHao<complex<double>, 2> tComplex(L,L);
    for(size_t i = 0; i < L; ++i)
    {
        for(size_t j = 0; j < L; ++j) tComplex(j,i) = t(j,i);
    }

    // ( wfCas^{\dagger} . t )^T = t^T . wfCas^c
    BL_NAME(gmm)(tComplex, wfCasUpConj, wfUpDaggerT_T, 'T');
    BL_NAME(gmm)(tComplex, wfCasDnConj, wfDnDaggerT_T, 'T');
}

void RealMaterialMoleculeMeasureFixedMDCas2sSD2s::initWfDaggerCholeskyVecs_T()
{
    size_t L = walkerLeft->getL(); size_t NupMax = walkerLeft->getNupMax(); size_t NdnMax = walkerLeft->getNdnMax();
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();

    wfUpDaggerCholeskyVecs_T.resize(L, NupMax, choleskyNumber); wfDnDaggerCholeskyVecs_T.resize(L, NdnMax, choleskyNumber);

    const TensorHao<double,3> & choleskyVecs = realMaterialMolecule->getCholeskyVecs();
    TensorHao<complex<double>, 2> choleskyVecComplex(L,L);
    TensorHaoRef<complex<double>, 2> wfDaggerCholeskyVec;
    for(size_t k = 0; k < choleskyNumber ; ++k)
    {
        for(size_t i = 0; i < L; ++i)
        {
            for(size_t j = 0; j < L; ++j) choleskyVecComplex(j,i) = choleskyVecs(j,i,k);
        }
        wfDaggerCholeskyVec=wfUpDaggerCholeskyVecs_T[k];
        BL_NAME(gmm)(choleskyVecComplex, wfCasUpConj, wfDaggerCholeskyVec, 'T');
        wfDaggerCholeskyVec=wfDnDaggerCholeskyVecs_T[k];
        BL_NAME(gmm)(choleskyVecComplex, wfCasDnConj, wfDaggerCholeskyVec, 'T');
    }
}

void RealMaterialMoleculeMeasureFixedMDCas2sSD2s::checkWalkerLeft(const MDCas2sSD2sOperation &mdCas2sSD2sOperation)
{
    if( walkerLeft != mdCas2sSD2sOperation.getWalkerLeft() )
    {
        cout<<"Error!!! RealMaterialMoleculeMeasureFixedMDCas2sSD2s only accept MDCas2sSD2sOperation with fixed MDCas2s!"<<endl;
        exit(1);
    }
}

void RealMaterialMoleculeMeasureFixedMDCas2sSD2s::addEnergy(MDCas2sSD2sOperation &mdCas2sSD2sOperation, std::complex<double> numIncrement)
{
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();

    if( choleskyBgNum.rank(0) != choleskyNumber ) { choleskyBgNum.resize(choleskyNumber); choleskyBgNum = complex<double>(0,0); }
    if( choleskyExNum.rank(0) != choleskyNumber ) { choleskyExNum.resize(choleskyNumber); choleskyExNum = complex<double>(0,0); }

    size_t MDSize = walkerLeft->getMDSize();
    size_t MDUpSize = walkerLeft->getMDUpSize();
    size_t MDDnSize = walkerLeft->getMDDnSize();

    const vector< complex<double> > &coe = walkerLeft->getCoe();
    const vector< size_t > &coeLinkUp = walkerLeft->getCoeLinkUp();
    const vector< size_t > &coeLinkDn = walkerLeft->getCoeLinkDn();

    const vector< complex<double> > &detUp = mdCas2sSD2sOperation.getDetUp();
    const vector< complex<double> > &detDn = mdCas2sSD2sOperation.getDetDn();
    const vector< TensorHao<complex<double>, 2> > &thetaUpVec = mdCas2sSD2sOperation.returnThetaUp();
    const vector< TensorHao<complex<double>, 2> > &thetaDnVec = mdCas2sSD2sOperation.returnThetaDn();

    vector<complex<double>> TenergyUp(MDUpSize), TenergyDn(MDDnSize);
    vector<TensorHao<complex<double>, 1>> choleskyBgUp(MDUpSize), choleskyBgDn(MDDnSize);
    vector<TensorHao<complex<double>, 1>> choleskyExUp(MDUpSize), choleskyExDn(MDDnSize);

    for(size_t i = 0; i <MDUpSize; ++i)
    {
        TenergyUp[i]    = calculateTenergyUp(i, thetaUpVec[i]);
        choleskyBgUp[i] = calculateCholeskyBgUp(i, thetaUpVec[i]);
        choleskyExUp[i] = calculateCholeskyExUp(i, thetaUpVec[i]);
    }

    for(size_t i = 0; i <MDDnSize; ++i)
    {
        TenergyDn[i]    = calculateTenergyDn(i, thetaDnVec[i]);
        choleskyBgDn[i] = calculateCholeskyBgDn(i, thetaDnVec[i]);
        choleskyExDn[i] = calculateCholeskyExDn(i, thetaDnVec[i]);
    }

    complex<double> HenergyOne, TenergyOne; TensorHao<complex<double>, 1> choleskyBgOne(choleskyNumber), choleskyExOne(choleskyNumber);
    complex<double> HenergyAdd, TenergyAdd; TensorHao<complex<double>, 1> choleskyBgAdd(choleskyNumber), choleskyExAdd(choleskyNumber);
    complex<double> zero(0.0); HenergyAdd = zero; TenergyAdd=zero; choleskyBgAdd=zero; choleskyExAdd=zero;
    complex<double> scaleFactor;
    size_t indexUp, indexDn;
    for(size_t i = 0; i < MDSize; ++i)
    {
        indexUp = coeLinkUp[i]; indexDn = coeLinkDn[i];
        TenergyOne = TenergyUp[indexUp] + TenergyDn[indexDn];
        choleskyBgOne = choleskyBgUp[indexUp] + choleskyBgDn[indexDn];
        choleskyExOne = choleskyExUp[indexUp] + choleskyExDn[indexDn];

        HenergyOne=0.0;
        for(size_t k = 0; k < choleskyNumber; ++k)
        {
            HenergyOne += ( choleskyBgOne(k)*choleskyBgOne(k) -choleskyExOne(k) );
        }
        HenergyOne *= 0.5;
        HenergyOne += TenergyOne;

        scaleFactor = numIncrement * conj( coe[i] ) * detUp[indexUp] * detDn[indexDn];

        TenergyAdd    += ( TenergyOne    * scaleFactor );
        choleskyBgAdd += ( choleskyBgOne * scaleFactor );
        choleskyExAdd += ( choleskyExOne * scaleFactor );
        HenergyAdd    += ( HenergyOne    * scaleFactor );
    }

    TNum          += TenergyAdd;
    choleskyBgNum += choleskyBgAdd;
    choleskyExNum += choleskyExAdd;
    HNum          += HenergyAdd;
}

complex<double> RealMaterialMoleculeMeasureFixedMDCas2sSD2s::calculateTenergyUp(size_t casIndex, const TensorHao<complex<double>,2> &thetaUp)
{
    size_t L = walkerLeft->getL(); size_t Nup = walkerLeft->getNup();
    const vector< vector<size_t> > &occupancyUp = walkerLeft->getOccupancyUp();
    size_t r;

    complex<double> TenergyUp(0,0);

    for(size_t i = 0; i < Nup; ++i)
    {
        r = occupancyUp[casIndex][i];
        for(size_t j = 0; j < L; ++j)
        {
            TenergyUp += wfUpDaggerT_T(j,r) * thetaUp(j,i);
        }
    }

    return TenergyUp;
}

complex<double> RealMaterialMoleculeMeasureFixedMDCas2sSD2s::calculateTenergyDn(size_t casIndex, const TensorHao<complex<double>, 2> &thetaDn)
{
    size_t L = walkerLeft->getL(); size_t Ndn = walkerLeft->getNdn();
    const vector< vector<size_t> > &occupancyDn = walkerLeft->getOccupancyDn();
    size_t r;

    complex<double> TenergyDn(0,0);

    for(size_t i = 0; i < Ndn; ++i)
    {
        r = occupancyDn[casIndex][i];
        for(size_t j = 0; j < L; ++j)
        {
            TenergyDn += wfDnDaggerT_T(j,r) * thetaDn(j,i);
        }
    }

    return TenergyDn;
}

TensorHao<complex<double>, 1> RealMaterialMoleculeMeasureFixedMDCas2sSD2s::calculateCholeskyBgUp(size_t casIndex, const TensorHao<complex<double>,2> &thetaUp)
{
    size_t L = walkerLeft->getL(); size_t Nup = walkerLeft->getNup();
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();
    const vector< vector<size_t> > &occupancyUp = walkerLeft->getOccupancyUp();

    size_t r;

    TensorHao<complex<double>, 1> choleskyBgUp(choleskyNumber);

    for(size_t k = 0; k < choleskyNumber; ++k)
    {
        choleskyBgUp(k) = 0.0;

        for(size_t i = 0; i < Nup; ++i)
        {
            r = occupancyUp[casIndex][i];
            for(size_t j = 0; j < L; ++j)
            {
                choleskyBgUp(k) += wfUpDaggerCholeskyVecs_T(j,r,k) * thetaUp(j,i);
            }
        }
    }

    return choleskyBgUp;
}

TensorHao<complex<double>, 1> RealMaterialMoleculeMeasureFixedMDCas2sSD2s::calculateCholeskyBgDn(size_t casIndex, const TensorHao<complex<double>, 2> &thetaDn)
{
    size_t L = walkerLeft->getL(); size_t Ndn = walkerLeft->getNdn();
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();
    const vector< vector<size_t> > &occupancyDn = walkerLeft->getOccupancyDn();

    size_t r;

    TensorHao<complex<double>, 1> choleskyBgDn(choleskyNumber);

    for(size_t k = 0; k < choleskyNumber; ++k)
    {
        choleskyBgDn(k) = 0.0;

        for(size_t i = 0; i < Ndn; ++i)
        {
            r = occupancyDn[casIndex][i];
            for(size_t j = 0; j < L; ++j)
            {
                choleskyBgDn(k) += wfDnDaggerCholeskyVecs_T(j,r,k) * thetaDn(j,i);
            }
        }
    }

    return choleskyBgDn;
}

TensorHao<complex<double>, 1> RealMaterialMoleculeMeasureFixedMDCas2sSD2s::calculateCholeskyExUp(size_t casIndex, const TensorHao<complex<double>,2> &thetaUp)
{
    size_t L = walkerLeft->getL(); size_t Nup = walkerLeft->getNup();
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();
    const vector< vector<size_t> > &occupancyUp = walkerLeft->getOccupancyUp();

    TensorHao<complex<double>, 2> densityUp(Nup, Nup), wfDaggerCholeskyVecUp_T(L, Nup);

    size_t r;

    TensorHao<complex<double>, 1> choleskyExUp(choleskyNumber);

    for(size_t k = 0; k < choleskyNumber; ++k)
    {
        for(size_t i = 0; i < Nup; ++i)
        {
            r = occupancyUp[casIndex][i];
            for(size_t j = 0; j < L; ++j)
            {
                wfDaggerCholeskyVecUp_T(j,i) = wfUpDaggerCholeskyVecs_T(j,r,k);
            }
        }
        BL_NAME(gmm)( wfDaggerCholeskyVecUp_T, thetaUp, densityUp, 'T' );

        choleskyExUp(k)=0.0;
        for(size_t i = 0; i < Nup; ++i) { for(size_t j = 0; j < Nup; ++j) choleskyExUp(k) += densityUp(j,i) * densityUp(i,j); }
    }

    return choleskyExUp;
}

TensorHao<complex<double>, 1> RealMaterialMoleculeMeasureFixedMDCas2sSD2s::calculateCholeskyExDn(size_t casIndex, const TensorHao<complex<double>,2> &thetaDn)
{
    size_t L = walkerLeft->getL(); size_t Ndn = walkerLeft->getNdn();
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();
    const vector< vector<size_t> > &occupancyDn = walkerLeft->getOccupancyDn();

    TensorHao<complex<double>, 2> densityDn(Ndn, Ndn), wfDaggerCholeskyVecDn_T(L, Ndn);

    size_t r;

    TensorHao<complex<double>, 1> choleskyExDn(choleskyNumber);

    for(size_t k = 0; k < choleskyNumber; ++k)
    {
        for(size_t i = 0; i < Ndn; ++i)
        {
            r = occupancyDn[casIndex][i];
            for(size_t j = 0; j < L; ++j)
            {
                wfDaggerCholeskyVecDn_T(j,i) = wfDnDaggerCholeskyVecs_T(j,r,k);
            }
        }
        BL_NAME(gmm)( wfDaggerCholeskyVecDn_T, thetaDn, densityDn, 'T' );

        choleskyExDn(k)=0.0;
        for(size_t i = 0; i < Ndn; ++i) { for(size_t j = 0; j < Ndn; ++j) choleskyExDn(k) += densityDn(j,i) * densityDn(i,j); }
    }

    return choleskyExDn;
}
