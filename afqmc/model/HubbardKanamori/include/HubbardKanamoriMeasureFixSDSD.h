//
// Created by Hao Shi on 1/17/18.
//

#ifndef AFQMCHUBBARDKANAMORI_HUBBARDKANAMORIMEASUREFIXSDSD_H
#define AFQMCHUBBARDKANAMORI_HUBBARDKANAMORIMEASUREFIXSDSD_H

#include "HubbardKanamori.h"
#include "../../../blocks/walkerWalkerOperation/SDSDOperation/include/SDSDOperation.h"


class HubbardKanamoriMeasureFixSDSD
{
 private:
    const HubbardKanamori * hubbardKanamori;
    const SD *walkerLeft;

    std::complex<double> den;
    std::complex<double> TNum, UNum, U1Num, U2Num, JNum, HNum;
    tensor_hao::TensorHao<std::complex<double>, 1> kanamoriBgNum;
    tensor_hao::TensorHao<std::complex<double>, 1> NupNum, NdnNum, SplusNum, SminusNum;
    std::complex<double> NupTotNum, NdnTotNum, SplusTotNum, SminusTotNum;

    tensor_hao::TensorHao<std::complex<double>,2> wfDaggerT;

 public:
    HubbardKanamoriMeasureFixSDSD();
    HubbardKanamoriMeasureFixSDSD(const HubbardKanamori& HubbardKanamori_, const SD& walkerLeft_);
    ~HubbardKanamoriMeasureFixSDSD();

    void initModelWalkerNullptr();
    void setModelWalker(const HubbardKanamori &HubbardKanamori_, const SD &walkerLeft_);
    void reSet();

    std::complex<double> returnEnergy();
    tensor_hao::TensorHao<double, 1> returnKanamoriBg();
    void addMeasurement(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    KanamoriInteractForce getForce(const KanamoriInteract &kanamoriInteract, SDSDOperation &sdsdOperation, double cap=1.0);

    void write() const;
    double getMemory() const;

 private:
    HubbardKanamoriMeasureFixSDSD(const HubbardKanamoriMeasureFixSDSD& x);
    HubbardKanamoriMeasureFixSDSD & operator = (const HubbardKanamoriMeasureFixSDSD& x);

    void initWfDaggerT();
    void checkWalkerLeft(const SDSDOperation &sdsdOperation);

    void addEnergy(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addKanamoriBg(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addNupNdn(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addSplusSminus(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
};

#endif //AFQMCHUBBARDKANAMORI_HUBBARDKANAMORIMEASUREFIXSDSD_H
