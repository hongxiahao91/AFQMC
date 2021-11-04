//
// Created by Hao Shi on 11/1/17.
//

#ifndef AFQMCCONSTRAINTPATH_HUBBARDKANAMORI_H
#define AFQMCCONSTRAINTPATH_HUBBARDKANAMORI_H

#include "../../../../common/common.h"
#include "../../../blocks/oneBodyOperator/hop/include/hop.h"
#include "../../../blocks/oneBodyOperator/logHop/include/logHop.h"
#include "../../../blocks/twoBodyOperator/kanamoriInteract/include/KanamoriInteract.h"

//See details in my note.
//t is non-interacting matrix

//For charge
//K is t - J/2
//Kp is K + JB

//For spin
//K is t + J/2
//Kp is K - JB

#ifdef MPI_HAO
class HubbardKanamori;
void MPIBcast(HubbardKanamori &buffer, int root=0,  const MPI_Comm& comm=MPI_COMM_WORLD);
#endif


class HubbardKanamori
{
 private:
    size_t L, N,numberOfKana;
    std::string decompTypeJ; //"charge", "spin"

    tensor_hao::TensorHao< std::complex<double>, 2 > t, K;
    tensor_hao::TensorHao<double,1> U;

    tensor_hao::TensorHao<int,1> site_i, site_j;
    tensor_hao::TensorHao<double,1> U1, U2, J;

    tensor_hao::TensorHao<double,1> kanamoriBg;

    size_t KpEigenStatus; //0: void, 1: calculated Kp, 2: calculated Kp Eigens.
    tensor_hao::TensorHao<std::complex<double>,2> Kp;
    tensor_hao::TensorHao<double,1> KpEigenValue;
    tensor_hao::TensorHao<std::complex<double>,2> KpEigenVector;

 public:
    HubbardKanamori();
    HubbardKanamori(const std::string &filename);
    ~HubbardKanamori();

    size_t getL() const;
    size_t getN() const;
    const std::string &getDecompTypeJ() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getT() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getK() const;
    const tensor_hao::TensorHao<double,1> &getU() const;
    size_t getNumberOfKana() const;
    const tensor_hao::TensorHao<int, 1> &getSite_i() const;
    const tensor_hao::TensorHao<int, 1> &getSite_j() const;
    const tensor_hao::TensorHao<double,1> &getU1() const;
    const tensor_hao::TensorHao<double,1> &getU2() const;
    const tensor_hao::TensorHao<double,1> &getJ() const;
    const tensor_hao::TensorHao<double,1> &getKanamoriBg() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getKp() const;
    const tensor_hao::TensorHao<std::complex<double>,2> &getKpEigenVector() const;

    void read(const std::string &filename);
    void write(const std::string &filename) const;

#ifdef MPI_HAO
    friend void MPIBcast(HubbardKanamori &buffer, int root,  const MPI_Comm& comm);
#endif

    void writeBackGround(const std::string &filename) const;
    void updateBackGround(const tensor_hao::TensorHao<double,1> &background);
    void updateBackGround(tensor_hao::TensorHao<double,1> &&background);

    Hop returnExpMinusAlphaK(double alpha);
    LogHop returnLogExpMinusAlphaK(double alpha);
    KanamoriInteract returnExpMinusAlphaV(double alpha,
                                          const std::string &decompTypeU,
                                          const std::string &decompTypeU1,
                                          const std::string &decompTypeU2J);

    void setKp();
    void setKpEigenValueAndVector();

    double getMemory() const;

 private:
    HubbardKanamori(const HubbardKanamori& x);
    HubbardKanamori & operator  = (const HubbardKanamori& x);

    void setK();
};

#endif //AFQMCCONSTRAINTPATH_HUBBARDKANAMORI_H
