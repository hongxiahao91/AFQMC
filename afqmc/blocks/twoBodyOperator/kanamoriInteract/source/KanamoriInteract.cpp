//
// Created by Hao Shi on 1/12/18.
//

#include "../../../../../common/common.h"
#include "../include/KanamoriInteract.h"

using namespace std;
using namespace tensor_hao;

#define pi 3.14159265358979324

KanamoriInteract::KanamoriInteract() {}

KanamoriInteract::KanamoriInteract(double dt,
                                   const string & decompTypeU,
                                   const string & decompTypeU1,
                                   const string & decompTypeU2J,
                                   const string & decompTypeJ,
                                   const TensorHao<double,1> &U,
                                   size_t numberOfKana,
                                   const TensorHao<int,1> &site_i,
                                   const TensorHao<int,1> &site_j,
                                   const TensorHao<double,1> &U1,
                                   const TensorHao<double,1> &U2J,
                                   const TensorHao<double,1> &J,
                                   const TensorHao<double, 1> &kanamoriBg)
{
    KanamoriInteract::dt =  dt;
    KanamoriInteract::L = U.size();
    KanamoriInteract::numberOfKana = numberOfKana;
    KanamoriInteract::site_i = &site_i;
    KanamoriInteract::site_j = &site_j;

    KanamoriInteract::decompTypeU = decompTypeU;
    KanamoriInteract::decompTypeU1 = decompTypeU1;
    KanamoriInteract::decompTypeU2J = decompTypeU2J;
    KanamoriInteract::decompTypeJ = decompTypeJ;

    KanamoriInteract::dtUSum = dt*U.sum();
    KanamoriInteract::dtU1Sum = dt*U1.sum();
    KanamoriInteract::dtU2JSum =  dt*U2J.sum();

    KanamoriInteract::dtU   =  dt*U;
    KanamoriInteract::dtU1  = dt*U1;
    KanamoriInteract::dtU2J = dt*U2J;
    KanamoriInteract::dtJ   = dt*J;

    setGammaU();
    setGammaU1();
    setGammaU2J();
    setGammaJ();

    KanamoriInteract::kanamoriBg = &kanamoriBg;
}

KanamoriInteract::KanamoriInteract(const KanamoriInteract &x) { copy_deep(x); }

KanamoriInteract::KanamoriInteract(KanamoriInteract&&x) { move_deep(x); }

KanamoriInteract::~KanamoriInteract() { }

KanamoriInteract &KanamoriInteract::operator=(const KanamoriInteract &x)  { copy_deep(x); return *this; }

KanamoriInteract &KanamoriInteract::operator=(KanamoriInteract &&x) { move_deep(x); return *this; }

const string &KanamoriInteract::getDecompTypeU() const { return decompTypeU; }

const string &KanamoriInteract::getDecompTypeU1() const { return decompTypeU1; }

const string &KanamoriInteract::getDecompTypeU2J() const { return decompTypeU2J; }

const string &KanamoriInteract::getDecompTypeJ() const { return decompTypeJ; }

const TensorHao<complex<double>, 1> &KanamoriInteract::getGammaU() const { return gammaU; }

const TensorHao<complex<double>, 1> &KanamoriInteract::getGammaU1() const { return gammaU1; }

const TensorHao<complex<double>, 1> &KanamoriInteract::getGammaU2J() const { return gammaU2J; }

const TensorHao<complex<double>, 1> &KanamoriInteract::getGammaJ() const { return gammaJ; }

std::complex<double> KanamoriInteract::calculateAuxForce(const KanamoriInteractAux &aux, const KanamoriInteractForce &force)
{
    complex<double> auxForce(0,0);
    for(size_t i = 0; i < L; ++i)
    {
        auxForce += aux.UAux(i) * force.UForce(i);
    }

    for(size_t i = 0; i < 2*numberOfKana; ++i)
    {
        auxForce += aux.U1Aux(i) * force.U1Force(i);
        auxForce += aux.U2JAux(i) * force.U2JForce(i);
    }

    for(size_t i = 0; i < numberOfKana; ++i)
    {
        auxForce += aux.JAux(i) * force.JForce(i);
    }

    return auxForce;
}

KanamoriInteractForce KanamoriInteract::readForce(const std::string &filename) const
{
    KanamoriInteractForce force(L, numberOfKana);

    if( !checkFile(filename) )
    {
        force.UForce=0.0;
        force.U1Force=0.0;
        force.U2JForce=0.0;
        force.JForce=complex<double>(0.0,0.0);
    }
    else
    {
        size_t L_in, numberOfKana_in;

        ifstream file;
        file.open(filename, ios::in);
        if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

        readFile( L_in, file );
        readFile( numberOfKana_in, file );
        if( L_in!=L || numberOfKana_in!=numberOfKana )
        {
            cout<<"Error!!! Input force has different L or numberOfKana "<<L<<" "<<L_in<<" "
                <<" "<<numberOfKana<<" "<<numberOfKana_in<<" !"<<endl;
            exit(1);
        }

        readFile( L, force.UForce.data(), file );
        readFile( 2*numberOfKana, force.U1Force.data(), file );
        readFile( 2*numberOfKana, force.U2JForce.data(), file );
        readFile( numberOfKana, force.JForce.data(), file );

        file.close();
    }

    return force;
}

KanamoriInteractAux KanamoriInteract::sampleAuxFromForce(const KanamoriInteractForce &force) const
{
    if( L != force.getL() || numberOfKana !=force.getNumberOfKana() )
    {
        cout<<"Error!!! In KanamoriInteract::sampleAuxFromForce, force size is not consistent!"<<endl;
        exit(1);
    }

    KanamoriInteractAux aux(L, numberOfKana);

    double expPlus, expMinus, prob;

    for(size_t i = 0; i < L; ++i)
    {
        expMinus = exp( -force.UForce(i) );
        expPlus  = exp(  force.UForce(i) );
        prob = expMinus / (expMinus + expPlus);

        if( uniformHao() < prob ) aux.UAux(i) = -1;
        else aux.UAux(i)=1;
    }

    for(size_t i = 0; i < 2*numberOfKana; ++i)
    {
        expMinus = exp( -force.U1Force(i) );
        expPlus  = exp(  force.U1Force(i) );
        prob = expMinus / (expMinus + expPlus);

        if( uniformHao() < prob ) aux.U1Aux(i) = -1;
        else aux.U1Aux(i)=1;
    }

    for(size_t i = 0; i < 2*numberOfKana; ++i)
    {
        expMinus = exp( -force.U2JForce(i) );
        expPlus  = exp(  force.U2JForce(i) );
        prob = expMinus / (expMinus + expPlus);

        if( uniformHao() < prob ) aux.U2JAux(i) = -1;
        else aux.U2JAux(i)=1;
    }

    for(size_t i = 0; i < numberOfKana; ++i)  aux.JAux(i) = gaussianHao() + force.JForce(i);

    return aux;
}

complex<double> KanamoriInteract::logProbOfAuxFromForce(const KanamoriInteractAux &aux, const KanamoriInteractForce &force) const
{
    if( L != aux.getL() || numberOfKana !=aux.getNumberOfKana() )
    {
        cout<<"Error!!! In KanamoriInteract::logProbOfAuxFromForce, aux size is not consistent!"<<endl;
        exit(1);
    }

    if( L != force.getL() || numberOfKana !=force.getNumberOfKana() )
    {
        cout<<"Error!!! In KanamoriInteract::logProbOfAuxFromForce, force size is not consistent!"<<endl;
        exit(1);
    }

    complex<double> logProb(0, 0);

    //U: exp( x*force ) / ( exp(force)+exp(-force) )
    for(size_t i=0; i<L; i++)
    {
        logProb += aux.UAux(i)*force.UForce(i) - log( exp(force.UForce(i)) + exp(-force.UForce(i)) );
    }

    //U1: exp( x*force ) / ( exp(force)+exp(-force) )
    for(size_t i=0; i<2*numberOfKana; i++)
    {
        logProb += aux.U1Aux(i)*force.U1Force(i) - log( exp(force.U1Force(i)) + exp(-force.U1Force(i)) );
    }

    //U2J: exp( x*force ) / ( exp(force)+exp(-force) )
    for(size_t i=0; i<2*numberOfKana; i++)
    {
        logProb += aux.U2JAux(i)*force.U2JForce(i) - log( exp(force.U2JForce(i)) + exp(-force.U2JForce(i)) );
    }

    // Product: 1/sqrt(2.0*Pi) * Exp( -(x-force)^2 / 2.0 )
    complex<double> auxMinusForceSquare(0,0), tmp;
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        tmp = aux.JAux(i)-force.JForce(i);
        auxMinusForceSquare += tmp*tmp;
    }
    logProb += -0.5*log(2.0*pi)*numberOfKana - 0.5*auxMinusForceSquare;

    return logProb;
}

KanamoriInteractSample KanamoriInteract::getTwoBodySampleFromAux(const KanamoriInteractAux &aux) const
{
    if( L != aux.getL() || numberOfKana !=aux.getNumberOfKana() )
    {
        cout<<"Error!!! In KanamoriInteract::getTwoBodySampleFromAux, aux size is not consistent!"<<endl;
        exit(1);
    }

    KanamoriInteractSample twoBodySample(2*L);

    //0.5*exp(dt*U*0.5)*exp(-aux*gamma)
    complex<double> logwU;
    if( decompTypeU == "densityCharge" )
    {
        KahanData<complex<double>> gammaAuxSum;
        for(size_t i = 0; i < L; ++i) gammaAuxSum += aux.UAux(i)*1.0*gammaU(i);
        logwU= log(0.5)*L + dtUSum*0.5 - gammaAuxSum.returnSum();
    }
    else if( decompTypeU == "densitySpin" )
    {
        logwU = log(0.5)*L;
    }

    //0.5*exp(dt*U*0.5)*exp(-aux*gamma)
    complex<double> logwU1;
    if( decompTypeU1 == "densityCharge" )
    {
        KahanData<complex<double>> gammaAuxSum;
        for(size_t i = 0; i < numberOfKana; ++i) gammaAuxSum += ( aux.U1Aux(i)+aux.U1Aux(i+numberOfKana) )*1.0*gammaU1(i);
        logwU1= log(0.5)*2*numberOfKana + dtU1Sum*2*0.5 - gammaAuxSum.returnSum();
    }
    else if( decompTypeU1 == "densitySpin" )
    {
        logwU1 = log(0.5)*2*numberOfKana;
    }

    //0.5*exp(dt*U*0.5)*exp(-aux*gamma)
    complex<double> logwU2J;
    if( decompTypeU2J == "densityCharge" )
    {
        KahanData<complex<double>> gammaAuxSum;
        for(size_t i = 0; i < numberOfKana; ++i) gammaAuxSum += ( aux.U2JAux(i)+aux.U2JAux(i+numberOfKana) ) *1.0*gammaU2J(i);
        logwU2J= log(0.5)*2*numberOfKana + dtU2JSum*2*0.5 - gammaAuxSum.returnSum();
    }
    else if( decompTypeU2J == "densitySpin" )
    {
        logwU2J = log(0.5)*2*numberOfKana;
    }

    complex<double> logwJ;
    // Product: 1/sqrt(2.0*Pi) * Exp( -x^2 / 2.0 ) * Exp( -gammaJ*x*B )
    complex<double> aux2Sum(0,0), auxBGammaJSum(0,0);
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        aux2Sum += aux.JAux(i)*aux.JAux(i);
        auxBGammaJSum += aux.JAux(i) * gammaJ(i) * kanamoriBg->operator()(i);
    }
    logwJ = -0.5*log(2.0*pi)*numberOfKana - 0.5*aux2Sum - auxBGammaJSum;


    twoBodySample.logw = logwU+logwU1+logwU2J+logwJ;

    setTwoBodySampleMatrix(twoBodySample, aux);

    return twoBodySample;
}

KanamoriInteractSample KanamoriInteract::getTwoBodySampleFromAuxForce(const KanamoriInteractAux &aux, const KanamoriInteractForce &force) const
{
    if( L != aux.getL() || numberOfKana !=aux.getNumberOfKana() )
    {
        cout<<"Error!!! In KanamoriInteract::getTwoBodySampleFromAuxForce, aux size is not consistent!"<<endl;
        exit(1);
    }

    if( L != force.getL() || numberOfKana !=force.getNumberOfKana() )
    {
        cout<<"Error!!! In KanamoriInteract::getTwoBodySampleFromAuxForce, force size is not consistent!"<<endl;
        exit(1);
    }

    KanamoriInteractSample twoBodySample = getTwoBodySampleFromAux(aux);

    twoBodySample.logw = twoBodySample.logw - logProbOfAuxFromForce(aux, force);

    return twoBodySample;
}

double KanamoriInteract::getMemory()
{
    double mem(0.0);
    mem += 8.0;
    mem += 8.0 + 8.0;
    mem += 8.0 + 8.0; //pointer
    mem += decompTypeU.size()+decompTypeU1.size()+decompTypeU2J.size()+decompTypeJ.size();
    mem += 8.0*3;
    mem += dtU.getMemory() + dtU1.getMemory() + dtU2J.getMemory() + dtJ.getMemory();
    mem += gammaU.getMemory() + gammaU1.getMemory() + gammaU2J.getMemory() + gammaJ.getMemory();
    mem += 8.0;
    return mem;
}

void KanamoriInteract::copy_deep(const KanamoriInteract &x)
{
    dt = x.dt;
    L = x.L;
    numberOfKana = x.numberOfKana;
    site_i = x.site_i;
    site_j = x.site_j;
    decompTypeU = x.decompTypeU;
    decompTypeU1 = x.decompTypeU1;
    decompTypeU2J = x.decompTypeU2J;
    decompTypeJ = x.decompTypeJ;
    dtUSum = x.dtUSum;
    dtU1Sum = x.dtU1Sum;
    dtU2JSum = x.dtU2JSum;
    dtU = x.dtU;
    dtU1 = x.dtU1;
    dtU2J = x.dtU2J;
    dtJ = x.dtJ;
    gammaU = x.gammaU;
    gammaU1 = x.gammaU1;
    gammaU2J = x.gammaU2J;
    gammaJ = x.gammaJ;
    kanamoriBg = x.kanamoriBg;
}

void KanamoriInteract::move_deep(KanamoriInteract &x)
{
    dt = move( x.dt );
    L = move( x.L );
    numberOfKana = move( x.numberOfKana );
    site_i = move( x.site_i );
    site_j = move( x.site_j );
    decompTypeU = move( x.decompTypeU );
    decompTypeU1 = move( x.decompTypeU1 );
    decompTypeU2J = move( x.decompTypeU2J );
    decompTypeJ = move( x.decompTypeJ );
    dtUSum = move( x.dtUSum );
    dtU1Sum = move( x.dtU1Sum );
    dtU2JSum = move( x.dtU2JSum );
    dtU = move( x.dtU );
    dtU1 = move( x.dtU1 );
    dtU2J = move( x.dtU2J );
    dtJ = move( x.dtJ );
    gammaU = move( x.gammaU );
    gammaU1 = move( x.gammaU1 );
    gammaU2J = move( x.gammaU2J );
    gammaJ = move( x.gammaJ );
    kanamoriBg = move( x.kanamoriBg );
}

void KanamoriInteract::setGammaU()
{
    gammaU.resize(L);

    if( decompTypeU == "densityCharge" )
    {
        for(size_t i=0; i<L; i++) gammaU(i) = solveCoshxEqExpy(-dtU(i) * 0.5 );
    }
    else if( decompTypeU == "densitySpin" )
    {
        for(size_t i=0; i<L; i++) gammaU(i) = solveCoshxEqExpy( dtU(i) * 0.5 );
    }
    else
    {
        cout<<"Error! Can not find the matched decompTypeU! "<<decompTypeU<<endl;
        exit(1);
    }
}

void KanamoriInteract::setGammaU1()
{
    gammaU1.resize(numberOfKana);

    if( decompTypeU1 == "densityCharge" )
    {
        for(size_t i=0; i<numberOfKana; i++) gammaU1(i) = solveCoshxEqExpy(-dtU1(i) * 0.5 );
    }
    else if( decompTypeU1 == "densitySpin" )
    {
        for(size_t i=0; i<numberOfKana; i++) gammaU1(i) = solveCoshxEqExpy( dtU1(i) * 0.5 );
    }
    else
    {
        cout<<"Error! Can not find the matched decompTypeU1! "<<decompTypeU1<<endl;
        exit(1);
    }
}

void KanamoriInteract::setGammaU2J()
{
    gammaU2J.resize(numberOfKana);

    if( decompTypeU2J == "densityCharge" )
    {
        for(size_t i=0; i<numberOfKana; i++) gammaU2J(i) = solveCoshxEqExpy(-dtU2J(i) * 0.5 );
    }
    else if( decompTypeU2J == "densitySpin" )
    {
        for(size_t i=0; i<numberOfKana; i++) gammaU2J(i) = solveCoshxEqExpy( dtU2J(i) * 0.5 );
    }
    else
    {
        cout<<"Error! Can not find the matched decompTypeU2J! "<<decompTypeU2J<<endl;
        exit(1);
    }
}

void KanamoriInteract::setGammaJ()
{
    gammaJ.resize(numberOfKana);

    if(decompTypeJ == "charge")
    {
        for(size_t i = 0; i < numberOfKana; ++i) gammaJ(i) = sqrt( -complex<double>(1.0, 0.0) * dtJ(i) );
    }
    else if(decompTypeJ =="spin")
    {
        for(size_t i = 0; i < numberOfKana; ++i) gammaJ(i) = sqrt(  complex<double>(1.0, 0.0) * dtJ(i) );
    }
}

void KanamoriInteract::setTwoBodySampleMatrix(KanamoriInteractSample &twoBodySample, const KanamoriInteractAux &aux) const
{
    //Set zero
    TensorHao<complex<double>,2> &matrix = twoBodySample.matrix;
    matrix = complex<double>(0.0, 0.0);

    complex<double> tmp;
    size_t ki, kj;

    //auxU
    //charge: (gamma*x-dtU*0.5) *(ni_up+ni_dn)
    //Spin: (gamma*x-dtU*0.5)*ni_up + (-gamma*x-dtU*0.5)*ni_dn
    if( decompTypeU == "densityCharge" )
    {
        for(size_t i = 0; i < L; ++i)
        {
            tmp = 1.0*aux.UAux(i)*gammaU(i) - dtU(i) * 0.5;
            matrix(i,i)      += tmp;
            matrix(i+L, i+L) += tmp;
        }
    }
    else if (decompTypeU == "densitySpin")
    {
        for(size_t i = 0; i < L; ++i)
        {
            matrix(i,i)      +=  1.0*aux.UAux(i)*gammaU(i) - dtU(i) * 0.5;
            matrix(i+L, i+L) += -1.0*aux.UAux(i)*gammaU(i) - dtU(i) * 0.5;
        }
    }
    else
    {
        cout<<"Error! Can not find the matched decompTypeU! "<<decompTypeU<<endl;
        exit(1);
    }

    //auxU1
    //charge: (gamma*x_1-dtU1*0.5) *(ni_up+nj_dn) + (gamma*x_2-dtU1*0.5) *(ni_dn+nj_up)
    //spin: (gamma*x_1-dtU1*0.5)*ni_up + (-gamma*x_1-dtU1*0.5)*nj_dn + (gamma*x_2-dtU1*0.5)*ni_dn + (-gamma*x_2-dtU1*0.5)*nj_up
    if( decompTypeU1 == "densityCharge" )
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i->operator()(i);
            kj = site_j->operator()(i);

            tmp = 1.0*aux.U1Aux(i)*gammaU1(i) - dtU1(i) * 0.5;
            matrix(ki, ki)         += tmp;
            matrix(kj + L, kj + L) += tmp;

            tmp = 1.0*aux.U1Aux(i + numberOfKana)*gammaU1(i) - dtU1(i) * 0.5;
            matrix(ki + L, ki + L) += tmp;
            matrix(kj, kj)         += tmp;
        }
    }
    else if( decompTypeU1 == "densitySpin" )
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i->operator()(i);
            kj = site_j->operator()(i);

            matrix(ki, ki)         +=  1.0*aux.U1Aux(i)*gammaU1(i) - dtU1(i) * 0.5;
            matrix(kj + L, kj + L) += -1.0*aux.U1Aux(i)*gammaU1(i) - dtU1(i) * 0.5;

            matrix(ki + L, ki + L) +=  1.0*aux.U1Aux(i + numberOfKana)*gammaU1(i) - dtU1(i) * 0.5;
            matrix(kj, kj)         += -1.0*aux.U1Aux(i + numberOfKana)*gammaU1(i) - dtU1(i) * 0.5;
        }
    }
    else
    {
        cout<<"Error! Can not find the matched decompTypeU1! "<<decompTypeU1<<endl;
        exit(1);
    }

    //auxU2J
    //charge: (gamma*x_1-dtU2J*0.5) *(ni_up+nj_up) + (gamma*x_2-dtU1*0.5) *(ni_dn+nj_dn)
    //spin: (gamma*x_1-dtU2J*0.5)*ni_up + (-gamma*x_1-dtU2J*0.5)*nj_up + (gamma*x_2-dtU2J*0.5)*ni_dn + (-gamma*x_2-dtU2J*0.5)*nj_dn
    if( decompTypeU2J == "densityCharge" )
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i->operator()(i);
            kj = site_j->operator()(i);

            tmp = 1.0*aux.U2JAux(i)*gammaU2J(i) - dtU2J(i) * 0.5;
            matrix(ki, ki) += tmp;
            matrix(kj, kj) += tmp;

            tmp = 1.0*aux.U2JAux(i + numberOfKana)*gammaU2J(i) - dtU2J(i) * 0.5;
            matrix(ki + L, ki + L) += tmp;
            matrix(kj + L, kj + L) += tmp;
        }
    }
    else if( decompTypeU2J == "densitySpin" )
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i->operator()(i);
            kj = site_j->operator()(i);

            matrix(ki, ki) +=  1.0*aux.U2JAux(i)*gammaU2J(i) - dtU2J(i) * 0.5;
            matrix(kj, kj) += -1.0*aux.U2JAux(i)*gammaU2J(i) - dtU2J(i) * 0.5;

            matrix(ki + L, ki + L) +=  1.0*aux.U2JAux(i + numberOfKana)*gammaU2J(i) - dtU2J(i) * 0.5;
            matrix(kj + L, kj + L) += -1.0*aux.U2JAux(i + numberOfKana)*gammaU2J(i) - dtU2J(i) * 0.5;
        }
    }
    else
    {
        cout<<"Error! Can not find the matched decompTypeU2J! "<<decompTypeU2J<<endl;
        exit(1);
    }

    if(decompTypeJ == "charge")
    {
        //auxJ * sqrt(-dt*J) * sum_sigma ( Ci^d Cj + Cj^d Ci )
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i->operator()(i);
            kj = site_j->operator()(i);

            tmp = aux.JAux(i) * gammaJ(i);
            matrix(ki, kj) += tmp;
            matrix(kj, ki) += tmp;
            matrix(ki + L, kj + L) += tmp;
            matrix(kj + L, ki + L) += tmp;
        }
    }
    else if(decompTypeJ =="spin")
    {
        //auxJ * sqrt(dt*J) * ( Ciu^d Cju + Cju^d Ciu  - Cid^d Cjd - Cjd^d Cid)
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i->operator()(i);
            kj = site_j->operator()(i);

            tmp = aux.JAux(i) * gammaJ(i);
            matrix(ki, kj) += tmp;
            matrix(kj, ki) += tmp;
            matrix(ki + L, kj + L) -= tmp;
            matrix(kj + L, ki + L) -= tmp;
        }
    }
}