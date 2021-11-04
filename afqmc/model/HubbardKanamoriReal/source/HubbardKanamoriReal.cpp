//
// Created by Hao Shi on 11/1/17.
//

#include "../include/HubbardKanamoriReal.h"

using namespace std;
using namespace H5;
using namespace tensor_hao;

HubbardKanamoriReal::HubbardKanamoriReal():L(0),N(0),numberOfKana(0),KpEigenStatus(0) { }

HubbardKanamoriReal::HubbardKanamoriReal(const string &filename) { read(filename); }

HubbardKanamoriReal::~HubbardKanamoriReal() { }

size_t HubbardKanamoriReal::getL() const { return L; }

size_t HubbardKanamoriReal::getN() const { return N; }

const string &HubbardKanamoriReal::getDecompTypeJ() const { return decompTypeJ; }

const TensorHao<complex<double>, 2> &HubbardKanamoriReal::getT() const  { return t; }

const TensorHao<complex<double>, 2> &HubbardKanamoriReal::getK() const { return K; }

const TensorHao<double,1> &HubbardKanamoriReal::getU() const { return U; }

size_t HubbardKanamoriReal::getNumberOfKana() const { return numberOfKana; }

const TensorHao<int, 1> &HubbardKanamoriReal::getSite_i() const { return site_i; }

const TensorHao<int, 1> &HubbardKanamoriReal::getSite_j() const { return site_j; }

const TensorHao<double,1> &HubbardKanamoriReal::getU1() const { return U1; }

const TensorHao<double,1> &HubbardKanamoriReal::getU2() const { return U2; }

const TensorHao<double,1> &HubbardKanamoriReal::getJ() const { return J; }

const TensorHao<double,1> &HubbardKanamoriReal::getKanamoriBg() const { return kanamoriBg; }

const TensorHao<complex<double>,2> &HubbardKanamoriReal::getKpEigenVector() const{ return KpEigenVector; }

void HubbardKanamoriReal::read(const string &filename)
{
    ifstream file;
    file.open(filename, ios::in);
    if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

    readFile( L, file );
    readFile( N, file );
    readFile(numberOfKana, file);
    readFile(decompTypeJ, file);

    t.resize(2*L,2*L); readFile( t.size(),  t.data(),  file );
    U.resize(L);       readFile( U.size(),  U.data(),  file );

    site_i.resize(numberOfKana); readFile(site_i.size(), site_i.data(), file);
    site_j.resize(numberOfKana); readFile(site_j.size(), site_j.data(), file);
    U1.resize(numberOfKana); readFile(U1.size(), U1.data(), file);
    U2.resize(numberOfKana); readFile(U2.size(), U2.data(), file);
    J.resize(numberOfKana); readFile(J.size(), J.data(), file);

    kanamoriBg.resize(numberOfKana);readFile(kanamoriBg.size(),kanamoriBg.data(),file);

    file.close();

    setK();

    KpEigenStatus = 0;
    Kp.resize(0,0);
    KpEigenValue.resize( static_cast<size_t>(0) );
    KpEigenVector.resize( 0, 0 );
}

void HubbardKanamoriReal::write(const string &filename) const
{
    ofstream file;
    file.open(filename, ios::out|ios::trunc);
    if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

    writeFile( L, file );
    writeFile( N, file );
    writeFile( numberOfKana, file);
    writeFile( decompTypeJ, file);

    writeFile( t.size(),  t.data(),  file );
    writeFile( U.size(),  U.data(),  file );

    writeFile( site_i.size(), site_i.data(), file);
    writeFile( site_j.size(), site_j.data(), file);
    writeFile( U1.size(), U1.data(), file );
    writeFile( U2.size(), U2.data(), file );
    writeFile( J.size(), J.data(), file );

    writeFile(kanamoriBg.size(), kanamoriBg.data(), file);

    file.close();
}

#ifdef MPI_HAO
void MPIBcast(HubbardKanamoriReal &buffer, int root, MPI_Comm const &comm)
{
    MPIBcast( buffer.L, root, comm  );
    MPIBcast( buffer.N, root, comm  );
    MPIBcast( buffer.numberOfKana, root, comm);
    MPIBcast( buffer.decompTypeJ, root, comm);

    MPIBcast( buffer.t, root, comm  );
    MPIBcast( buffer.K, root, comm  );
    MPIBcast( buffer.U, root, comm  );

    MPIBcast( buffer.site_i, root, comm  );
    MPIBcast( buffer.site_j, root, comm  );

    MPIBcast( buffer.U1, root, comm  );
    MPIBcast( buffer.U2, root, comm  );
    MPIBcast( buffer.J, root, comm  );

    MPIBcast( buffer.kanamoriBg, root, comm );

    MPIBcast( buffer.KpEigenStatus, root, comm );
    MPIBcast( buffer.Kp, root, comm );
    MPIBcast( buffer.KpEigenValue, root, comm );
    MPIBcast( buffer.KpEigenVector, root, comm );
}
#endif

void HubbardKanamoriReal::writeBackGround(const std::string &filename) const
{
    //We are not using HDF5 file here, can not write background only, call write function instead.
    write(filename);
}

void HubbardKanamoriReal::updateBackGround(const tensor_hao::TensorHao<double, 1> &background)
{
    if( background.size() != numberOfKana ) {cout<<"Error!!! Background size is not numberOfKana!"<<endl; exit(1);}
    KpEigenStatus = 0;
    kanamoriBg = background;
}

void HubbardKanamoriReal::updateBackGround(tensor_hao::TensorHao<double, 1> &&background)
{
    if( background.size() != numberOfKana ) {cout<<"Error!!! Background size is not numberOfKana!"<<endl; exit(1);}
    KpEigenStatus = 0;
    kanamoriBg = move(background);
}

Hop HubbardKanamoriReal::returnExpMinusAlphaK(double alpha)
{
    setKpEigenValueAndVector();

    Hop hop(2*L);
    double Jbg2(0.0);for(size_t i = 0; i < numberOfKana; ++i) Jbg2 += J(i)*kanamoriBg(i)*kanamoriBg(i);

    if(decompTypeJ == "charge") hop.logw = alpha*0.5*Jbg2;
    else if(decompTypeJ == "spin") hop.logw = -alpha*0.5*Jbg2;
    else { cout<<"Error!!! Do not know the decompTypeJ "<<decompTypeJ<<endl; exit(1); }

    BL_NAME(gmm)( KpEigenVector, dMultiMatrix( exp( -alpha*KpEigenValue ), conjtrans(KpEigenVector) ), hop.matrix );

    return hop;
}

LogHop HubbardKanamoriReal::returnLogExpMinusAlphaK(double alpha)
{
    setKp();

    LogHop logHop(2*L);
    double Jbg2(0.0);for(size_t i = 0; i < numberOfKana; ++i) Jbg2 += J(i)*kanamoriBg(i)*kanamoriBg(i);

    if(decompTypeJ == "charge") logHop.logw = alpha * 0.5 * Jbg2;
    else if(decompTypeJ == "spin") logHop.logw = -alpha * 0.5 * Jbg2;
    else { cout<<"Error!!! Do not know the decompTypeJ "<<decompTypeJ<<endl; exit(1); }

    logHop.matrix = complex<double>(-alpha, 0) * Kp;

    return logHop;
}

KanamoriInteractReal HubbardKanamoriReal::returnExpMinusAlphaV(double alpha,
                                                               const std::string &decompTypeU,
                                                               const std::string &decompTypeU1,
                                                               const std::string &decompTypeU2J)
{
    if(decompTypeJ == "charge")
    {
        return KanamoriInteractReal(alpha, decompTypeU, decompTypeU1, decompTypeU2J, decompTypeJ,
                                    U, numberOfKana, site_i, site_j, U1, U2 + J, J, kanamoriBg);
    }
    else if(decompTypeJ == "spin")
    {
        return KanamoriInteractReal(alpha, decompTypeU, decompTypeU1, decompTypeU2J, decompTypeJ,
                                    U, numberOfKana, site_i, site_j, U1, U2 - J, J, kanamoriBg);
    }
    else
    {
        cout<<"Error!!! Do not know the decompTypeJ "<<decompTypeJ<<endl;
        exit(1);
    }
}

void HubbardKanamoriReal::setKp()
{
    if( KpEigenStatus >=1 ) return;

    Kp = K;

    size_t ki,kj; double tmp;
    if(decompTypeJ == "charge")
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki=site_i(i); kj=site_j(i);

            tmp = J(i)*kanamoriBg(i);
            Kp(ki,kj)     += tmp;
            Kp(kj,ki)     += tmp;
            Kp(ki+L,kj+L) += tmp;
            Kp(kj+L,ki+L) += tmp;
        }
    }
    else if(decompTypeJ == "spin")
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki=site_i(i); kj=site_j(i);

            tmp = J(i)*kanamoriBg(i);
            Kp(ki,kj)     -= tmp;
            Kp(kj,ki)     -= tmp;
            Kp(ki+L,kj+L) -= tmp;
            Kp(kj+L,ki+L) -= tmp;
        }
    }
    else
    {
        cout<<"Error!!! Do not know the decompTypeJ "<<decompTypeJ<<endl;
        exit(1);
    }

    KpEigenStatus=1;
}

void HubbardKanamoriReal::setKpEigenValueAndVector()
{
    if( KpEigenStatus >=2 ) return;

    setKp();
    KpEigenVector = Kp;
    KpEigenValue.resize(2*L);
    BL_NAME(eigen)(KpEigenVector, KpEigenValue);

    KpEigenStatus = 2;
}

double HubbardKanamoriReal::getMemory() const
{
    double mem(0.0);
    mem += 8.0+8.0+8.0+decompTypeJ.size(); //L, N, numberOfKana, decompTypeJ
    mem += t.getMemory() + K.getMemory()+ U.getMemory();
    mem += site_i.getMemory()+ site_j.getMemory();
    mem += U1.getMemory() + U2.getMemory() + J.getMemory();
    mem += kanamoriBg.getMemory();
    mem += 8.0+Kp.getMemory()+KpEigenValue.getMemory()+KpEigenVector.getMemory();
    return mem;
}

HubbardKanamoriReal::HubbardKanamoriReal(const HubbardKanamoriReal &x) { }

HubbardKanamoriReal &HubbardKanamoriReal::operator=(const HubbardKanamoriReal &x) { return *this; }

void HubbardKanamoriReal::setK()
{
    size_t ki,kj;

    K = t;

    if(decompTypeJ == "charge")
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i(i);
            K(ki,   ki  ) -= J(i)*0.5;
            K(ki+L, ki+L) -= J(i)*0.5;

            kj = site_j(i);
            K(kj,   kj  ) -= J(i)*0.5;
            K(kj+L, kj+L) -= J(i)*0.5;
        }
    }
    else if(decompTypeJ == "spin")
    {
        for(size_t i = 0; i < numberOfKana; ++i)
        {
            ki = site_i(i);
            K(ki,   ki  ) += J(i)*0.5;
            K(ki+L, ki+L) += J(i)*0.5;

            kj = site_j(i);
            K(kj,   kj  ) += J(i)*0.5;
            K(kj+L, kj+L) += J(i)*0.5;
        }
    }
    else
    {
        cout<<"Error!!! Do not know the decompTypeJ "<<decompTypeJ<<endl;
        exit(1);
    }
}