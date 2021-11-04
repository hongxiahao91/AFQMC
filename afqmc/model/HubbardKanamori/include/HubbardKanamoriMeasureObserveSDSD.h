//
// Created by Hao Shi on 3/23/18.
//

#ifndef AFQMCLAB_HUBBARDKANAMORIMEASUREOBSERVESDSD_H
#define AFQMCLAB_HUBBARDKANAMORIMEASUREOBSERVESDSD_H

#include "HubbardKanamori.h"
#include "../../../blocks/walkerWalkerOperation/SDSDOperation/include/SDSDOperation.h"

class HubbardKanamoriMeasureObserveSDSD
{
 private:
    const HubbardKanamori * hubbardKanamori;

    std::complex<double> den;
    std::complex<double> TNum, UNum, U1Num, U2Num, JNum, HNum;
    tensor_hao::TensorHao<std::complex<double>, 1> kanamoriBgNum;
    tensor_hao::TensorHao<std::complex<double>, 1> NupNum, NdnNum, SplusNum, SminusNum;
    std::complex<double> NupTotNum, NdnTotNum, SplusTotNum, SminusTotNum;

    tensor_hao::TensorHao<std::complex<double>, 2> greenMatrixNum;
    tensor_hao::TensorHao<std::complex<double>, 2> densityDensityNum;
    tensor_hao::TensorHao<std::complex<double>, 2> splusSminusNum;
    tensor_hao::TensorHao<std::complex<double>, 2> sminusSplusNum;

 public:
    HubbardKanamoriMeasureObserveSDSD();
    HubbardKanamoriMeasureObserveSDSD(const HubbardKanamori& HubbardKanamori_);
    ~HubbardKanamoriMeasureObserveSDSD();

    void initModelNullptr();
    void setModel(const HubbardKanamori &HubbardKanamori_);
    void reSet();

    std::complex<double> returnEnergy();
    tensor_hao::TensorHao<double, 1> returnKanamoriBg();
    void addMeasurement(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    KanamoriInteractForce getForce(const KanamoriInteract &kanamoriInteract, SDSDOperation &sdsdOperation, double cap=1.0);

    void write() const;
    double getMemory() const;

 private:
    HubbardKanamoriMeasureObserveSDSD(const HubbardKanamoriMeasureObserveSDSD& x);
    HubbardKanamoriMeasureObserveSDSD & operator = (const HubbardKanamoriMeasureObserveSDSD& x);

    void addEnergy(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addKanamoriBg(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addNupAndNdn(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addSplusAndSminus(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);

    void addGreenMatrix(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addDensityDensity(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addSplusSminus(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
    void addSminusSplus(SDSDOperation &sdsdOperation, std::complex<double> denIncrement);
};

#endif //AFQMCLAB_HUBBARDKANAMORIMEASUREOBSERVESDSD_H
