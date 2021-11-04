//
// Created by Hao Shi on 3/23/18.
//

#include "../include/HubbardKanamoriMeasureObserveSDSD.h"
#include "../../../utilities/manipulateMCData/include/writeThreadSum.h"

using namespace std;
using namespace tensor_hao;

HubbardKanamoriMeasureObserveSDSD::HubbardKanamoriMeasureObserveSDSD()
{
    initModelNullptr();
    reSet();
}

HubbardKanamoriMeasureObserveSDSD::HubbardKanamoriMeasureObserveSDSD(const HubbardKanamori &HubbardKanamori_)
{
    setModel(HubbardKanamori_);
    reSet();
}

HubbardKanamoriMeasureObserveSDSD::~HubbardKanamoriMeasureObserveSDSD()
{

}

void HubbardKanamoriMeasureObserveSDSD::initModelNullptr()
{
    hubbardKanamori = nullptr;
}

void HubbardKanamoriMeasureObserveSDSD::setModel(const HubbardKanamori &HubbardKanamori_)
{
    hubbardKanamori = &HubbardKanamori_;
}

void HubbardKanamoriMeasureObserveSDSD::reSet()
{
    complex<double> zero(0,0);
    den = zero;
    TNum = zero; UNum = zero; U1Num = zero; U2Num = zero; JNum = zero; HNum = zero;
    kanamoriBgNum = zero;
    NupNum = zero; NdnNum = zero; SplusNum = zero; SminusNum = zero;
    NupTotNum = zero; NdnTotNum = zero; SplusTotNum = zero; SminusTotNum = zero;

    greenMatrixNum = zero; densityDensityNum = zero; splusSminusNum = zero; sminusSplusNum = zero;
}

complex<double> HubbardKanamoriMeasureObserveSDSD::returnEnergy()
{
    complex<double> Htot   = MPISum(HNum);
    complex<double> denTot = MPISum(den);
    complex<double> energy;
    if( MPIRank() == 0 ) energy = Htot/denTot;
    MPIBcast(energy);
    return energy;
}

TensorHao<double, 1> HubbardKanamoriMeasureObserveSDSD::returnKanamoriBg()
{
    size_t numberOfKana = hubbardKanamori->getNumberOfKana();

    TensorHao<complex<double>, 1> kanamoriBgTot( numberOfKana );
    MPISum( numberOfKana, kanamoriBgNum.data(), kanamoriBgTot.data() );

    complex<double> denTot = MPISum(den);

    TensorHao<double,1> kanamoriBg( numberOfKana );
    for(size_t i = 0; i < numberOfKana; ++i) kanamoriBg(i) = ( kanamoriBgTot(i)/denTot ).real();
    MPIBcast(kanamoriBg);

    return kanamoriBg;
}

void HubbardKanamoriMeasureObserveSDSD::addMeasurement(SDSDOperation &sdsdOperation, std::complex<double> denIncrement)
{
    den += denIncrement;

    addEnergy(sdsdOperation, denIncrement);
    addKanamoriBg(sdsdOperation, denIncrement);
    addNupAndNdn(sdsdOperation, denIncrement);
    addSplusAndSminus(sdsdOperation, denIncrement);

    addGreenMatrix(sdsdOperation, denIncrement);
    addDensityDensity(sdsdOperation, denIncrement);
    addSplusSminus(sdsdOperation, denIncrement);
    addSminusSplus(sdsdOperation, denIncrement);
}

KanamoriInteractForce HubbardKanamoriMeasureObserveSDSD::getForce(const KanamoriInteract &kanamoriInteract, SDSDOperation &sdsdOperation, double cap)
{
    size_t L = hubbardKanamori->getL();
    size_t numberOfKana = hubbardKanamori->getNumberOfKana();
    const string &decompTypeU = kanamoriInteract.getDecompTypeU();
    const string &decompTypeU1 = kanamoriInteract.getDecompTypeU1();
    const string &decompTypeU2J = kanamoriInteract.getDecompTypeU2J();
    const string &decompTypeJ = kanamoriInteract.getDecompTypeJ();
    const TensorHao<int, 1> &site_i = hubbardKanamori->getSite_i();
    const TensorHao<int, 1> &site_j = hubbardKanamori->getSite_j();
    size_t ki,kj;

    KanamoriInteractForce force(L,numberOfKana);

    const TensorHao< complex<double>, 2 > &greenMatrix = sdsdOperation.returnGreenMatrix();

    //U force
    TensorHao< complex<double>, 1 > backGroundU(L);
    if( decompTypeU == "densityCharge" )
    {
        for(size_t i = 0; i < L; ++i) backGroundU(i) = greenMatrix(i, i) + greenMatrix(i+L, i+L) -1.0;
    }
    else if( decompTypeU == "densitySpin" )
    {
        for(size_t i = 0; i < L; ++i) backGroundU(i) = greenMatrix(i,i) - greenMatrix(i+L,i+L);
    }
    else
    {
        cout<<"Error! Can not find the matched decompTypeU! "<<decompTypeU<<endl;
        exit(1);
    }

    const TensorHao<complex<double>, 1> &gammaU = kanamoriInteract.getGammaU();
    for (size_t i = 0; i < L; ++i)
    {
        force.UForce(i) = ( gammaU(i) * backGroundU(i) ).real();
        if( force.UForce(i) >  cap ) force.UForce(i) =  cap;
        if( force.UForce(i) < -cap ) force.UForce(i) = -cap;
    }

    //U1 force
    TensorHao< complex<double>, 1 > backGroundU1(2*numberOfKana);
    if( decompTypeU1 == "densityCharge" )
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki=site_i(i);kj=site_j(i);
            backGroundU1(i) = greenMatrix(ki, ki) + greenMatrix(kj+L, kj+L) -1.0;
            backGroundU1(i+numberOfKana) = greenMatrix(ki+L, ki+L) + greenMatrix(kj, kj) -1.0;
        }
    }
    else if( decompTypeU1 == "densitySpin" )
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki=site_i(i);kj=site_j(i);
            backGroundU1(i) = greenMatrix(ki, ki) - greenMatrix(kj+L, kj+L);
            backGroundU1(i+numberOfKana) = greenMatrix(ki+L, ki+L) - greenMatrix(kj, kj);
        }
    }
    else
    {
        cout<<"Error! Can not find the matched decompTypeU1! "<<decompTypeU1<<endl;
        exit(1);
    }

    const TensorHao<complex<double>, 1> &gammaU1 = kanamoriInteract.getGammaU1();
    for (size_t i = 0; i < numberOfKana; ++i)
    {
        force.U1Force(i) = ( gammaU1(i) * backGroundU1(i) ).real();
        if( force.U1Force(i) >  cap ) force.U1Force(i) =  cap;
        if( force.U1Force(i) < -cap ) force.U1Force(i) = -cap;

        force.U1Force(i+numberOfKana) = ( gammaU1(i) * backGroundU1(i+numberOfKana) ).real();
        if( force.U1Force(i+numberOfKana) >  cap ) force.U1Force(i+numberOfKana) =  cap;
        if( force.U1Force(i+numberOfKana) < -cap ) force.U1Force(i+numberOfKana) = -cap;
    }

    //U2J force
    TensorHao< complex<double>, 1 > backGroundU2J(2*numberOfKana);
    if( decompTypeU2J == "densityCharge" )
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki=site_i(i);kj=site_j(i);
            backGroundU2J(i) = greenMatrix(ki,ki) + greenMatrix(kj,kj) -1.0;
            backGroundU2J(i+numberOfKana) = greenMatrix(ki+L,ki+L) + greenMatrix(kj+L, kj+L) -1.0;
        }
    }
    else if( decompTypeU2J == "densitySpin" )
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki=site_i(i);kj=site_j(i);
            backGroundU2J(i) = greenMatrix(ki,ki) - greenMatrix(kj,kj);
            backGroundU2J(i+numberOfKana) = greenMatrix(ki+L,ki+L) - greenMatrix(kj+L, kj+L);
        }
    }
    else
    {
        cout<<"Error! Can not find the matched decompTypeU2J! "<<decompTypeU2J<<endl;
        exit(1);
    }

    const TensorHao<complex<double>, 1> &gammaU2J = kanamoriInteract.getGammaU2J();
    for (size_t i = 0; i < numberOfKana; ++i)
    {
        force.U2JForce(i) = ( gammaU2J(i) * backGroundU2J(i) ).real();
        if( force.U2JForce(i) >  cap ) force.U2JForce(i) =  cap;
        if( force.U2JForce(i) < -cap ) force.U2JForce(i) = -cap;

        force.U2JForce(i+numberOfKana) = ( gammaU2J(i) * backGroundU2J(i+numberOfKana) ).real();
        if( force.U2JForce(i+numberOfKana) >  cap ) force.U2JForce(i+numberOfKana) =  cap;
        if( force.U2JForce(i+numberOfKana) < -cap ) force.U2JForce(i+numberOfKana) = -cap;
    }

    //JForce
    TensorHao< complex<double>, 1 > kanamoriGroundJ(numberOfKana);
    if(decompTypeJ == "charge")
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki=site_i(i);kj=site_j(i);
            kanamoriGroundJ(i) = greenMatrix(ki,kj) + greenMatrix(kj,ki) + greenMatrix(ki+L,kj+L)+ greenMatrix(kj+L,ki+L);
        }
    }
    else if(decompTypeJ =="spin")
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki=site_i(i);kj=site_j(i);
            kanamoriGroundJ(i) = greenMatrix(ki,kj) + greenMatrix(kj,ki) - greenMatrix(ki+L,kj+L) - greenMatrix(kj+L,ki+L);
        }
    }

    const TensorHao<double,1> &currentBg = hubbardKanamori->getKanamoriBg();
    const TensorHao< complex<double>, 1 > &gammaJ = kanamoriInteract.getGammaJ();
    complex<double> oneForce;
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        oneForce = ( kanamoriGroundJ(i) - currentBg(i) ) * gammaJ(i);

        if( abs(oneForce) > cap ) force.JForce(i) = oneForce*cap/abs(oneForce);
        else force.JForce(i) = oneForce;
    }

    return force;
}

void HubbardKanamoriMeasureObserveSDSD::write() const
{
    writeThreadSum(den, "den.dat", ios::app);

    writeThreadSum(TNum, "TNum.dat", ios::app);
    writeThreadSum(UNum, "UNum.dat", ios::app);
    writeThreadSum(U1Num, "U1Num.dat", ios::app);
    writeThreadSum(U2Num, "U2Num.dat", ios::app);
    writeThreadSum(JNum, "JNum.dat", ios::app);
    writeThreadSum(HNum, "HNum.dat", ios::app);

    writeThreadSum(kanamoriBgNum.size(),    kanamoriBgNum.data(),    "kanamoriBgNum.dat",    ios::app);
    writeThreadSum(NupNum.size(),    NupNum.data(),    "NupNum.dat",    ios::app);
    writeThreadSum(NdnNum.size(),    NdnNum.data(),    "NdnNum.dat",    ios::app);
    writeThreadSum(SplusNum.size(),  SplusNum.data(),  "SplusNum.dat",  ios::app);
    writeThreadSum(SminusNum.size(), SminusNum.data(), "SminusNum.dat", ios::app);

    writeThreadSum(NupTotNum,   "NupTotNum.dat",   ios::app);
    writeThreadSum(NdnTotNum,   "NdnTotNum.dat",   ios::app);
    writeThreadSum(SplusTotNum, "SplusTotNum.dat", ios::app);
    writeThreadSum(SminusTotNum, "SminusTotNum.dat", ios::app);

    writeThreadSum(greenMatrixNum.size(), greenMatrixNum.data(), "greenMatrixNum.dat", ios::app);
    writeThreadSum(densityDensityNum.size(), densityDensityNum.data(), "densityDensityNum.dat", ios::app);
    writeThreadSum(splusSminusNum.size(), splusSminusNum.data(), "splusSminusNum.dat", ios::app);
    writeThreadSum(sminusSplusNum.size(), sminusSplusNum.data(), "sminusSplusNum.dat", ios::app);
}

double HubbardKanamoriMeasureObserveSDSD::getMemory() const
{
    double mem(0.0);
    mem += 8.0;
    mem += 8.0;
    mem += 16.0;
    mem += 16.0*6;
    mem += kanamoriBgNum.getMemory();
    mem += NupNum.getMemory()+NdnNum.getMemory()+SplusNum.getMemory()+SminusNum.getMemory();
    mem += 16.0*4;
    mem += greenMatrixNum.getMemory()+densityDensityNum.getMemory()+splusSminusNum.getMemory()+sminusSplusNum.getMemory();
    return mem;
}

HubbardKanamoriMeasureObserveSDSD::HubbardKanamoriMeasureObserveSDSD(const HubbardKanamoriMeasureObserveSDSD &x)
{

}

HubbardKanamoriMeasureObserveSDSD & HubbardKanamoriMeasureObserveSDSD::operator=(const HubbardKanamoriMeasureObserveSDSD &x)
{
    return *this;
}

void HubbardKanamoriMeasureObserveSDSD::addEnergy(SDSDOperation &sdsdOperation, std::complex<double> denIncrement)
{
    size_t L = hubbardKanamori->getL(); size_t L2 = L*2;
    size_t numberOfKana = hubbardKanamori->getNumberOfKana();
    const TensorHao<int, 1> &site_i = hubbardKanamori->getSite_i();
    const TensorHao<int, 1> &site_j = hubbardKanamori->getSite_j();
    size_t ki,kj;

    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();

    complex<double> Tenergy(0,0), Uenergy(0,0), U1energy(0,0), U2energy(0,0), Jenergy(0,0);

    //Add T
    TensorHao< complex<double>, 2 > t = hubbardKanamori->getT();
    for(size_t i = 0; i < L2; ++i)
    {
        for(size_t j = 0; j < L2; ++j) Tenergy +=  t(j,i) * greenMatrix(j,i);
    }

    //Add U
    const TensorHao< double, 1> &U = hubbardKanamori->getU();
    for(size_t i = 0; i < L; ++i)
    {
        Uenergy += U(i) * ( greenMatrix(i,i)*greenMatrix(i+L,i+L) - greenMatrix(i, i+L)*greenMatrix(i+L,i) );
    }

    //Add U1
    const TensorHao< double, 1> &U1 = hubbardKanamori->getU1();
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        ki = site_i(i); kj = site_j(i);
        U1energy += U1(i) * (  greenMatrix(ki,ki)*greenMatrix(kj+L,kj+L)-greenMatrix(ki,kj+L)*greenMatrix(kj+L, ki)
                               +greenMatrix(ki+L,ki+L)*greenMatrix(kj,kj)-greenMatrix(ki+L,kj)*greenMatrix(kj, ki+L)
        );
    }

    //Add U2
    const TensorHao< double, 1> &U2 = hubbardKanamori->getU2();
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        ki = site_i(i); kj = site_j(i);
        U2energy += U2(i) * (  greenMatrix(ki,ki)*greenMatrix(kj,kj)-greenMatrix(ki,kj)*greenMatrix(kj, ki)
                               +greenMatrix(ki+L,ki+L)*greenMatrix(kj+L,kj+L)-greenMatrix(ki+L,kj+L)*greenMatrix(kj+L, ki+L)
        );
    }

    //Add J
    const TensorHao< double, 1> &J = hubbardKanamori->getJ();
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        ki = site_i(i); kj = site_j(i);
        Jenergy += J(i) * (  greenMatrix(ki,kj)*greenMatrix(kj+L,ki+L)-greenMatrix(ki,ki+L)*greenMatrix(kj+L,kj)
                             +greenMatrix(ki,kj)*greenMatrix(ki+L,kj+L)-greenMatrix(ki,kj+L)*greenMatrix(ki+L,kj)
                             +greenMatrix(kj,ki)*greenMatrix(ki+L,kj+L)-greenMatrix(kj,kj+L)*greenMatrix(ki+L,ki)
                             +greenMatrix(kj,ki)*greenMatrix(kj+L,ki+L)-greenMatrix(kj,ki+L)*greenMatrix(kj+L,ki)
        );
    }

    TNum  += ( Tenergy  * denIncrement );
    UNum  += ( Uenergy  * denIncrement );
    U1Num += ( U1energy * denIncrement );
    U2Num += ( U2energy * denIncrement );
    JNum  += ( Jenergy  * denIncrement );
    HNum  += ( ( Tenergy + Uenergy + U1energy + U2energy + Jenergy ) * denIncrement );
}

void HubbardKanamoriMeasureObserveSDSD::addKanamoriBg(SDSDOperation &sdsdOperation, complex<double> denIncrement)
{
    size_t L = hubbardKanamori->getL();
    size_t numberOfKana = hubbardKanamori->getNumberOfKana();
    const TensorHao<int, 1> &site_i = hubbardKanamori->getSite_i();
    const TensorHao<int, 1> &site_j = hubbardKanamori->getSite_j();
    const std::string &decompTypeJ  = hubbardKanamori->getDecompTypeJ();
    size_t ki,kj;

    if( kanamoriBgNum.rank(0) != numberOfKana ) { kanamoriBgNum.resize(numberOfKana); kanamoriBgNum = complex<double>(0,0); }

    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();
    if(decompTypeJ == "charge")
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i(i); kj = site_j(i);
            kanamoriBgNum(i) += (greenMatrix(ki, kj) + greenMatrix(kj, ki) + greenMatrix(ki + L, kj + L) + greenMatrix(kj + L, ki + L)) * denIncrement;
        }
    }
    else if(decompTypeJ == "spin")
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i(i); kj = site_j(i);
            kanamoriBgNum(i) += (greenMatrix(ki, kj) + greenMatrix(kj, ki) - greenMatrix(ki + L, kj + L) - greenMatrix(kj + L, ki + L)) * denIncrement;
        }
    }
    else
    {
        cout<<"Error!!! Do not know the decompTypeJ "<<decompTypeJ<<endl;
        exit(1);
    }
}

void HubbardKanamoriMeasureObserveSDSD::addNupAndNdn(SDSDOperation &sdsdOperation, complex<double> denIncrement)
{
    size_t L = hubbardKanamori->getL();
    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();

    if( NupNum.rank(0) != L ) { NupNum.resize(L); NupNum = complex<double>(0,0); }
    if( NdnNum.rank(0) != L ) { NdnNum.resize(L); NdnNum = complex<double>(0,0); }

    for (size_t i = 0; i < L; ++i)
    {
        NupNum(i) += greenMatrix(i, i) * denIncrement;
        NdnNum(i) += greenMatrix(i+L, i+L) * denIncrement;
    }

    complex<double> diagUpAll(0,0), diagDnAll(0,0);
    for(size_t i = 0; i < L ; ++i)
    {
        diagUpAll += greenMatrix(i, i);
        diagDnAll += greenMatrix(i+L, i+L);
    }

    NupTotNum += diagUpAll * denIncrement;
    NdnTotNum += diagDnAll * denIncrement;
}

void HubbardKanamoriMeasureObserveSDSD::addSplusAndSminus(SDSDOperation &sdsdOperation,
                                                          complex<double> denIncrement)
{
    size_t L = hubbardKanamori->getL();
    const TensorHao<complex<double>, 2> &greenMatrix = sdsdOperation.returnGreenMatrix();

    if( SplusNum.rank(0) != L )  { SplusNum.resize(L);  SplusNum = complex<double>(0,0); }
    if( SminusNum.rank(0) != L ) { SminusNum.resize(L); SminusNum = complex<double>(0,0); }

    for (size_t i = 0; i < L; ++i)
    {
        SplusNum(i)  += greenMatrix(i, i+L) * denIncrement;
        SminusNum(i) += greenMatrix(i+L, i) * denIncrement;
    }

    complex<double> SplusAll(0,0), SminusAll(0,0);
    for(size_t i = 0; i < L; ++i)
    {
        SplusAll +=  greenMatrix(i, i+L);
        SminusAll += greenMatrix(i+L, i);
    }

    SplusTotNum  += SplusAll * denIncrement;
    SminusTotNum += SminusAll * denIncrement;
}

void HubbardKanamoriMeasureObserveSDSD::addGreenMatrix(SDSDOperation &sdsdOperation, complex<double> denIncrement)
{
    const TensorHao< complex<double>, 2 > &greenMatrix = sdsdOperation.returnGreenMatrix();

    size_t L2 = hubbardKanamori->getL() * 2;

    if( greenMatrixNum.rank(0) != L2 ) { greenMatrixNum.resize(L2, L2); greenMatrixNum = complex<double>(0,0); }

    greenMatrixNum += ( greenMatrix * denIncrement );
}

void HubbardKanamoriMeasureObserveSDSD::addDensityDensity(SDSDOperation &sdsdOperation, std::complex<double> denIncrement)
{
    const TensorHao< complex<double>, 2 > &greenMatrix = sdsdOperation.returnGreenMatrix();

    size_t L2 = hubbardKanamori->getL() * 2;

    if( densityDensityNum.rank(0) != L2 ) { densityDensityNum.resize(L2, L2); densityDensityNum = complex<double>(0,0); }

    complex<double> temp;
    for(size_t j = 0; j < L2; ++j)
    {
        for(size_t i = 0; i < L2; ++i)
        {
            if( i==j ) temp = greenMatrix(i,i);
            else temp = greenMatrix(i,i) * greenMatrix(j,j) - greenMatrix(i,j)*greenMatrix(j,i);
            densityDensityNum(i, j) += ( temp * denIncrement );
        }
    }
}

void HubbardKanamoriMeasureObserveSDSD::addSplusSminus(SDSDOperation &sdsdOperation, complex<double> denIncrement)
{
    const TensorHao< complex<double>, 2 > &greenMatrix = sdsdOperation.returnGreenMatrix();

    size_t L = hubbardKanamori->getL();

    if( splusSminusNum.rank(0) != L ) { splusSminusNum.resize(L, L); splusSminusNum = complex<double>(0,0); }

    complex<double> temp;
    for(size_t j = 0; j < L; ++j)
    {
        for(size_t i = 0; i < L; ++i)
        {
            if( i==j )
            {
                temp  = greenMatrix(i,i);
                temp -= ( greenMatrix(i,i)*greenMatrix(i+L,i+L)-greenMatrix(i,i+L)*greenMatrix(i+L,i) );
            }
            else
            {
                temp = -greenMatrix(i,j)*greenMatrix(j+L,i+L)+greenMatrix(i,i+L)*greenMatrix(j+L, j);
            }
            splusSminusNum(i,j) += ( temp * denIncrement );
        }
    }
}

void HubbardKanamoriMeasureObserveSDSD::addSminusSplus(SDSDOperation &sdsdOperation, complex<double> denIncrement)
{
    const TensorHao< complex<double>, 2 > &greenMatrix = sdsdOperation.returnGreenMatrix();

    size_t L = hubbardKanamori->getL();

    if( sminusSplusNum.rank(0) != L ) { sminusSplusNum.resize(L, L); sminusSplusNum = complex<double>(0,0); }

    complex<double> temp;
    for(size_t j = 0; j < L; ++j)
    {
        for(size_t i = 0; i < L; ++i)
        {
            if( i==j )
            {
                temp  = greenMatrix(i+L, i+L);
                temp -= ( greenMatrix(i,i)*greenMatrix(i+L,i+L)-greenMatrix(i,i+L)*greenMatrix(i+L,i) );
            }
            else
            {
                temp = -greenMatrix(i+L,j+L)*greenMatrix(j,i)+greenMatrix(i+L,i)*greenMatrix(j, j+L);
            }
            sminusSplusNum(i,j) += ( temp * denIncrement );
        }
    }
}
