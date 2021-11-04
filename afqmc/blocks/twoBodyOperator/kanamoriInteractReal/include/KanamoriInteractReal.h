//
// Created by Hao Shi on 1/12/18.
//

#ifndef AFQMCHUBBARDKANAMORI_KANAMORIINTERACTREAL_H
#define AFQMCHUBBARDKANAMORI_KANAMORIINTERACTREAL_H

#include "KanamoriInteractRealAux.h"
#include "KanamoriInteractRealForce.h"
#include "KanamoriInteractRealSample.h"

class KanamoriInteractReal
{
 private:
    double dt;
    size_t L, numberOfKana;
    const tensor_hao::TensorHao<int,1> *site_i, *site_j;

    std::string decompTypeU, decompTypeU1, decompTypeU2J;  //densityCharge, densitySpin
    std::string decompTypeJ; //charge, spin
    double dtUSum, dtU1Sum, dtU2JSum;
    tensor_hao::TensorHao<double,1> dtU, dtU1, dtU2J, dtJ;
    tensor_hao::TensorHao<std::complex<double>,1> gammaU, gammaU1, gammaU2J, gammaJ;

    const tensor_hao::TensorHao<double, 1> *kanamoriBg;

 public:
    KanamoriInteractReal();
    KanamoriInteractReal(double dt,
                         const std::string & decompTypeU,
                         const std::string & decompTypeU1,
                         const std::string & decompTypeU2J,
                         const std::string & decompTypeJ,
                         const tensor_hao::TensorHao<double,1> &U,
                         size_t numberOfKana,
                         const tensor_hao::TensorHao<int,1> &site_i,
                         const tensor_hao::TensorHao<int,1> &site_j,
                         const tensor_hao::TensorHao<double,1> &U1,
                         const tensor_hao::TensorHao<double,1> &U2J,
                         const tensor_hao::TensorHao<double,1> &J,
                         const tensor_hao::TensorHao<double, 1> &kanamoriBg);
    KanamoriInteractReal(const KanamoriInteractReal& x);
    KanamoriInteractReal(KanamoriInteractReal&& x);
    ~KanamoriInteractReal();

    KanamoriInteractReal & operator  = (const KanamoriInteractReal& x);
    KanamoriInteractReal & operator  = (KanamoriInteractReal&& x);

    const std::string &getDecompTypeU() const;
    const std::string &getDecompTypeU1() const;
    const std::string &getDecompTypeU2J() const;
    const std::string &getDecompTypeJ() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getGammaU() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getGammaU1() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getGammaU2J() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getGammaJ() const;

    std::complex<double> calculateAuxForce(const KanamoriInteractRealAux &aux, const KanamoriInteractRealForce &force);
    KanamoriInteractRealForce readForce(const std::string &filename) const;
    KanamoriInteractRealAux sampleAuxFromForce(const KanamoriInteractRealForce &force) const;
    std::complex<double> logProbOfAuxFromForce(const KanamoriInteractRealAux &aux, const KanamoriInteractRealForce &force) const;
    KanamoriInteractRealSample getTwoBodySampleFromAux(const KanamoriInteractRealAux &aux) const;
    KanamoriInteractRealSample getTwoBodySampleFromAuxForce(const KanamoriInteractRealAux &aux, const KanamoriInteractRealForce &force) const;

    double getMemory();
 private:
    void copy_deep(const KanamoriInteractReal &x);
    void move_deep(KanamoriInteractReal &x);

    void setGammaU();
    void setGammaU1();
    void setGammaU2J();
    void setGammaJ();

    void setTwoBodySampleMatrix(KanamoriInteractRealSample &twoBodySample, const KanamoriInteractRealAux &aux) const;

};

#endif //AFQMCHUBBARDKANAMORI_KANAMORIINTERACTREAL_H
