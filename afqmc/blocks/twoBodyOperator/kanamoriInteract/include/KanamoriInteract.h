//
// Created by Hao Shi on 1/12/18.
//

#ifndef AFQMCHUBBARDKANAMORI_KANAMORIINTERACT_H
#define AFQMCHUBBARDKANAMORI_KANAMORIINTERACT_H

#include "KanamoriInteractAux.h"
#include "KanamoriInteractForce.h"
#include "KanamoriInteractSample.h"

class KanamoriInteract
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
    KanamoriInteract();
    KanamoriInteract(double dt,
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
    KanamoriInteract(const KanamoriInteract& x);
    KanamoriInteract(KanamoriInteract&& x);
    ~KanamoriInteract();

    KanamoriInteract & operator  = (const KanamoriInteract& x);
    KanamoriInteract & operator  = (KanamoriInteract&& x);

    const std::string &getDecompTypeU() const;
    const std::string &getDecompTypeU1() const;
    const std::string &getDecompTypeU2J() const;
    const std::string &getDecompTypeJ() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getGammaU() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getGammaU1() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getGammaU2J() const;
    const tensor_hao::TensorHao<std::complex<double>, 1> &getGammaJ() const;

    std::complex<double> calculateAuxForce(const KanamoriInteractAux &aux, const KanamoriInteractForce &force);
    KanamoriInteractForce readForce(const std::string &filename) const;
    KanamoriInteractAux sampleAuxFromForce(const KanamoriInteractForce &force) const;
    std::complex<double> logProbOfAuxFromForce(const KanamoriInteractAux &aux, const KanamoriInteractForce &force) const;
    KanamoriInteractSample getTwoBodySampleFromAux(const KanamoriInteractAux &aux) const;
    KanamoriInteractSample getTwoBodySampleFromAuxForce(const KanamoriInteractAux &aux, const KanamoriInteractForce &force) const;

    double getMemory();
 private:
    void copy_deep(const KanamoriInteract &x);
    void move_deep(KanamoriInteract &x);

    void setGammaU();
    void setGammaU1();
    void setGammaU2J();
    void setGammaJ();

    void setTwoBodySampleMatrix(KanamoriInteractSample &twoBodySample, const KanamoriInteractAux &aux) const;

};

#endif //AFQMCHUBBARDKANAMORI_KANAMORIINTERACT_H
