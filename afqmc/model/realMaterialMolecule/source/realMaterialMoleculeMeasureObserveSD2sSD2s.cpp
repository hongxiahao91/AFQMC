//
// Created by boruoshihao on 11/14/18.
//

#include "../include/realMaterialMoleculeMeasureObserveSD2sSD2s.h"
#include "../../../utilities/manipulateMCData/include/writeThreadSum.h"

using namespace std;
using namespace tensor_hao;

RealMaterialMoleculeMeasureObserveSD2sSD2s::RealMaterialMoleculeMeasureObserveSD2sSD2s()
{
    initModelNullptr();
    reSet();
}

RealMaterialMoleculeMeasureObserveSD2sSD2s::RealMaterialMoleculeMeasureObserveSD2sSD2s(const RealMaterialMolecule &realMaterialMolecule_)
{
    setModel(realMaterialMolecule_);
    reSet();
}

RealMaterialMoleculeMeasureObserveSD2sSD2s::~RealMaterialMoleculeMeasureObserveSD2sSD2s()
{

}

void RealMaterialMoleculeMeasureObserveSD2sSD2s::initModelNullptr()
{
    realMaterialMolecule = nullptr;
}

void RealMaterialMoleculeMeasureObserveSD2sSD2s::setModel(const RealMaterialMolecule &realMaterialMolecule_)
{
    realMaterialMolecule = &realMaterialMolecule_;
}

void RealMaterialMoleculeMeasureObserveSD2sSD2s::reSet()
{
    complex<double> zero(0,0);

    den = zero;
    TNum = zero;
    choleskyBgNum = zero;
    choleskyExNum = zero;
    HNum = zero;
    greenUpNum = zero;
    greenDnNum = zero;
}

complex<double> RealMaterialMoleculeMeasureObserveSD2sSD2s::returnEnergy()
{
    complex<double> Htot   = MPISum(HNum);
    complex<double> denTot = MPISum(den);
    complex<double> energy;
    if( MPIRank() == 0 ) energy = Htot/denTot;
    MPIBcast(energy);
    return energy;
}

TensorHao<double,1> RealMaterialMoleculeMeasureObserveSD2sSD2s::returnCholeskyBg()
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

void RealMaterialMoleculeMeasureObserveSD2sSD2s::addMeasurement(SD2sSD2sOperation &sd2sSD2sOperation, complex<double> denIncrement)
{
    checkWalkerWithModel(sd2sSD2sOperation);
    initWfDaggerT(sd2sSD2sOperation);
    initWfDaggerCholeskyVecs(sd2sSD2sOperation);

    den += denIncrement;

    addEnergy(sd2sSD2sOperation, denIncrement);
    addGreen(sd2sSD2sOperation, denIncrement);

    //Reset matrix, incase we use them for other SD2sSD2sOperation
    //This is a temporary change, better way is to avoid using these matrix all the time.
    wfUpDaggerT.resize(0, 0); wfDnDaggerT.resize(0, 0);
    wfUpDaggerCholeskyVecs.resize(0, 0, 0); wfDnDaggerCholeskyVecs.resize(0, 0, 0);

}

CholeskyRealForce RealMaterialMoleculeMeasureObserveSD2sSD2s::getForce(const CholeskyReal &choleskyReal,
                                                                      SD2sSD2sOperation &sd2sSD2sOperation,
                                                                      double cap)
{
    checkWalkerWithModel(sd2sSD2sOperation);
    initWfDaggerCholeskyVecs(sd2sSD2sOperation);

    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();
    const TensorHao<double,1> & currentBg = realMaterialMolecule->getCholeskyBg();
    complex<double> sqrtMinusDt = choleskyReal.getSqrtMinusDt();

    TensorHao<complex<double>, 1> choleskyBg = calculateCholeskyBg(sd2sSD2sOperation);
    CholeskyRealForce force(choleskyNumber); complex<double> oneForce;
    for(size_t i = 0; i < choleskyNumber; ++i)
    {
        oneForce = (choleskyBg(i)-currentBg(i)) * sqrtMinusDt;

        if( abs(oneForce) > cap ) force(i) = oneForce*cap/abs(oneForce);
        else force(i) = oneForce;
    }

    //Reset matrix, incase we use them for other SD2sSD2sOperation
    //This is a temporary change, better way is to avoid using these matrix all the time.
    wfUpDaggerCholeskyVecs.resize(0, 0, 0); wfDnDaggerCholeskyVecs.resize(0, 0, 0);

    return force;
}

void RealMaterialMoleculeMeasureObserveSD2sSD2s::write() const
{
    writeThreadSum(den, "den.dat", ios::app);
    writeThreadSum(TNum, "TNum.dat", ios::app);
    writeThreadSum(choleskyBgNum.size(), choleskyBgNum.data(), "choleskyBgNum.dat", ios::app);
    writeThreadSum(choleskyExNum.size(), choleskyExNum.data(), "choleskyExNum.dat", ios::app);
    writeThreadSum(HNum, "HNum.dat", ios::app);
    writeThreadSum(greenUpNum.size(), greenUpNum.data(), "greenUpNum.dat", ios::app);
    writeThreadSum(greenDnNum.size(), greenDnNum.data(), "greenDnNum.dat", ios::app);

}

double RealMaterialMoleculeMeasureObserveSD2sSD2s::getMemory() const
{
    double mem(0.0);
    mem += 8.0;
    mem += 16.0*3; //den, TNum, HNum
    mem += choleskyBgNum.getMemory()+choleskyExNum.getMemory();
    mem += greenUpNum.getMemory() + greenDnNum.getMemory();
    mem += wfUpDaggerT.getMemory() + wfDnDaggerT.getMemory();
    mem += wfUpDaggerCholeskyVecs.getMemory() + wfDnDaggerCholeskyVecs.getMemory();
    return mem;
}

RealMaterialMoleculeMeasureObserveSD2sSD2s::RealMaterialMoleculeMeasureObserveSD2sSD2s(const RealMaterialMoleculeMeasureObserveSD2sSD2s &x)
{

}

RealMaterialMoleculeMeasureObserveSD2sSD2s & RealMaterialMoleculeMeasureObserveSD2sSD2s::operator=(const RealMaterialMoleculeMeasureObserveSD2sSD2s &x)
{
    return *this;
}

void RealMaterialMoleculeMeasureObserveSD2sSD2s::checkWalkerWithModel(const SD2sSD2sOperation &sd2sSD2sOperation)
{
    const SD2s *walkerLeft = sd2sSD2sOperation.getWalkerLeft();
    const SD2s *walkerRight = sd2sSD2sOperation.getWalkerRight();

    if( realMaterialMolecule->getL() != walkerLeft->getL() ) {cout<<"Model L does not consistent with walker left L!"<<endl; exit(1);}
    if( realMaterialMolecule->getNup() != walkerLeft->getNup() ) {cout<<"Model Nup does not consistent with walker left Nup!"<<endl; exit(1);}
    if( realMaterialMolecule->getNdn() != walkerLeft->getNdn() ) {cout<<"Model Ndn does not consistent with walker left Ndn!"<<endl; exit(1);}

    if( realMaterialMolecule->getL() != walkerRight->getL() ) {cout<<"Model L does not consistent with walker right L!"<<endl; exit(1);}
    if( realMaterialMolecule->getNup() != walkerRight->getNup() ) {cout<<"Model Nup does not consistent with walker right Nup!"<<endl; exit(1);}
    if( realMaterialMolecule->getNdn() != walkerRight->getNdn() ) {cout<<"Model Ndn does not consistent with walker right Ndn!"<<endl; exit(1);}
}

void RealMaterialMoleculeMeasureObserveSD2sSD2s::initWfDaggerT(SD2sSD2sOperation &sd2sSD2sOperation)
{
    const SD2s  *walkerLeft = sd2sSD2sOperation.getWalkerLeft();

    size_t L = walkerLeft->getL(); size_t Nup = walkerLeft->getNup(); size_t Ndn = walkerLeft->getNdn();

    wfUpDaggerT.resize(Nup, L); wfDnDaggerT.resize(Ndn, L);

    const TensorHao<double,2> & t = realMaterialMolecule->getT();
    TensorHao<complex<double>, 2> tComplex(L,L);
    for(size_t i = 0; i < L; ++i)
    {
        for(size_t j = 0; j < L; ++j) tComplex(j,i) = t(j,i);
    }

    BL_NAME(gmm)(walkerLeft->getWfUp(), tComplex, wfUpDaggerT, 'C');
    BL_NAME(gmm)(walkerLeft->getWfDn(), tComplex, wfDnDaggerT, 'C');
}

void RealMaterialMoleculeMeasureObserveSD2sSD2s::initWfDaggerCholeskyVecs(SD2sSD2sOperation &sd2sSD2sOperation)
{
    const SD2s  *walkerLeft = sd2sSD2sOperation.getWalkerLeft();

    size_t L = walkerLeft->getL(); size_t Nup = walkerLeft->getNup(); size_t Ndn = walkerLeft->getNdn();
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();

    wfUpDaggerCholeskyVecs.resize(Nup, L, choleskyNumber); wfDnDaggerCholeskyVecs.resize(Ndn, L, choleskyNumber);

    const TensorHao<double,3> & choleskyVecs = realMaterialMolecule->getCholeskyVecs();
    TensorHao<complex<double>, 2> choleskyVecComplex(L,L);
    TensorHaoRef<complex<double>, 2> wfDaggerCholeskyVec;
    for(size_t k = 0; k < choleskyNumber ; ++k)
    {
        for(size_t i = 0; i < L; ++i)
        {
            for(size_t j = 0; j < L; ++j) choleskyVecComplex(j,i) = choleskyVecs(j,i,k);
        }
        wfDaggerCholeskyVec=wfUpDaggerCholeskyVecs[k];
        BL_NAME(gmm)(walkerLeft->getWfUp(), choleskyVecComplex, wfDaggerCholeskyVec, 'C');
        wfDaggerCholeskyVec=wfDnDaggerCholeskyVecs[k];
        BL_NAME(gmm)(walkerLeft->getWfDn(), choleskyVecComplex, wfDaggerCholeskyVec, 'C');
    }
}

void RealMaterialMoleculeMeasureObserveSD2sSD2s::addEnergy(SD2sSD2sOperation &sd2sSD2sOperation, complex<double> denIncrement)
{
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();

    if( choleskyBgNum.rank(0) != choleskyNumber ) { choleskyBgNum.resize(choleskyNumber); choleskyBgNum = complex<double>(0,0); }
    if( choleskyExNum.rank(0) != choleskyNumber ) { choleskyExNum.resize(choleskyNumber); choleskyExNum = complex<double>(0,0); }

    complex<double> Tenergy=calculateTenergy(sd2sSD2sOperation);
    TensorHao<complex<double>, 1> choleskyBg = calculateCholeskyBg(sd2sSD2sOperation);
    TensorHao<complex<double>, 1> choleskyEx = calculateCholeskyEx(sd2sSD2sOperation);

    complex<double> Henergy(0.0);
    for(size_t i = 0; i < choleskyNumber; ++i)
    {
        Henergy += ( choleskyBg(i)*choleskyBg(i) -choleskyEx(i) );
    }
    Henergy *= 0.5;
    Henergy += Tenergy;

    TNum += ( Tenergy * denIncrement );
    choleskyBgNum += ( choleskyBg * denIncrement );
    choleskyExNum += ( choleskyEx * denIncrement );
    HNum += ( Henergy * denIncrement );
}

void RealMaterialMoleculeMeasureObserveSD2sSD2s::addGreen(SD2sSD2sOperation &sd2sSD2sOperation, complex<double> denIncrement)
{
    size_t L = realMaterialMolecule->getL();
    if( greenUpNum.rank(0) != L ) { greenUpNum.resize(L,L); greenUpNum = complex<double>(0,0); }
    if( greenDnNum.rank(0) != L ) { greenDnNum.resize(L,L); greenDnNum = complex<double>(0,0); }

    greenUpNum += ( sd2sSD2sOperation.returnGreenMatrixUp() * denIncrement );
    greenDnNum += ( sd2sSD2sOperation.returnGreenMatrixDn() * denIncrement );
}

complex<double> RealMaterialMoleculeMeasureObserveSD2sSD2s::calculateTenergy(SD2sSD2sOperation &sd2sSD2sOperation)
{
    size_t L = realMaterialMolecule->getL(); size_t Nup = realMaterialMolecule->getNup(); size_t Ndn = realMaterialMolecule->getNdn();
    const TensorHao<complex<double>, 2> &thetaUp_T =  sd2sSD2sOperation.returnThetaUp_T();
    const TensorHao<complex<double>, 2> &thetaDn_T =  sd2sSD2sOperation.returnThetaDn_T();

    complex<double> Tenergy(0,0);
    for(size_t i = 0; i < L; ++i)
    {
        for(size_t j = 0; j < Nup; ++j) Tenergy += wfUpDaggerT(j,i) * thetaUp_T(j,i);
        for(size_t j = 0; j < Ndn; ++j) Tenergy += wfDnDaggerT(j,i) * thetaDn_T(j,i);
    }
    return Tenergy;
}

TensorHao<complex<double>, 1> RealMaterialMoleculeMeasureObserveSD2sSD2s::calculateCholeskyBg(SD2sSD2sOperation &sd2sSD2sOperation)
{
    size_t L = realMaterialMolecule->getL(); size_t Nup = realMaterialMolecule->getNup(); size_t Ndn = realMaterialMolecule->getNdn();
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();
    const TensorHao<complex<double>, 2> &thetaUp_T =  sd2sSD2sOperation.returnThetaUp_T();
    const TensorHao<complex<double>, 2> &thetaDn_T =  sd2sSD2sOperation.returnThetaDn_T();

    TensorHaoRef<complex<double>, 2> leftUp(L*Nup, choleskyNumber), leftDn(L*Ndn, choleskyNumber);
    TensorHaoRef<complex<double>, 1> rightUp(L*Nup), rightDn(L*Ndn);
    leftUp.point( wfUpDaggerCholeskyVecs.data() );
    leftDn.point( wfDnDaggerCholeskyVecs.data() );
    rightUp.point( const_cast<complex<double>*>( thetaUp_T.data() ) );
    rightDn.point( const_cast<complex<double>*>( thetaDn_T.data() ) );

    TensorHao<complex<double>, 1> choleskyBg(choleskyNumber);
    BL_NAME(gemv)(leftUp, rightUp, choleskyBg, 'T' );
    BL_NAME(gemv)(leftDn, rightDn, choleskyBg, 'T', 1.0, 1.0 );

    return choleskyBg;
}

TensorHao<complex<double>, 1> RealMaterialMoleculeMeasureObserveSD2sSD2s::calculateCholeskyEx(SD2sSD2sOperation &sd2sSD2sOperation)
{
    size_t Nup = realMaterialMolecule->getNup(); size_t Ndn = realMaterialMolecule->getNdn();
    size_t choleskyNumber = realMaterialMolecule->getCholeskyNumber();
    const TensorHao<complex<double>, 2> &thetaUp_T =  sd2sSD2sOperation.returnThetaUp_T();
    const TensorHao<complex<double>, 2> &thetaDn_T =  sd2sSD2sOperation.returnThetaDn_T();

    TensorHao<complex<double>, 2> densityUp(Nup, Nup), densityDn(Ndn, Ndn);
    TensorHaoRef<complex<double>, 2> wfDaggerCholeskyVec;
    TensorHao<complex<double>, 1> choleskyEx(choleskyNumber);

    TensorHao<complex<double>, 2> thetaUp, thetaDn;
    thetaUp = trans(thetaUp_T); thetaDn = trans(thetaDn_T);
    for(size_t k = 0; k < choleskyNumber; ++k)
    {
        wfDaggerCholeskyVec = wfUpDaggerCholeskyVecs[k];
        BL_NAME(gmm)( wfDaggerCholeskyVec, thetaUp, densityUp);
        wfDaggerCholeskyVec = wfDnDaggerCholeskyVecs[k];
        BL_NAME(gmm)( wfDaggerCholeskyVec, thetaDn, densityDn);

        choleskyEx(k)=0.0;
        for(size_t i = 0; i < Nup; ++i) { for(size_t j = 0; j < Nup; ++j) choleskyEx(k) += densityUp(j,i) * densityUp(i,j); }
        for(size_t i = 0; i < Ndn; ++i) { for(size_t j = 0; j < Ndn; ++j) choleskyEx(k) += densityDn(j,i) * densityDn(i,j); }
    }

    return choleskyEx;
}