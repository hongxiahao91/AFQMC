//
// Created by boruoshihao on 11/14/18.
//

#ifndef AFQMCLAB_REALMATERIALMOLECULEOBSERVESD2SSD2S_H
#define AFQMCLAB_REALMATERIALMOLECULEOBSERVESD2SSD2S_H

#include "realMaterialMolecule.h"
#include "../../../blocks/walkerWalkerOperation/SD2sSD2sOperation/include/SD2sSD2sOperation.h"

class RealMaterialMoleculeMeasureObserveSD2sSD2s
{
 private:
    const RealMaterialMolecule *realMaterialMolecule;

    std::complex<double> den;
    std::complex<double> TNum;
    tensor_hao::TensorHao<std::complex<double>, 1> choleskyBgNum;
    tensor_hao::TensorHao<std::complex<double>, 1> choleskyExNum;
    std::complex<double> HNum;
    tensor_hao::TensorHao<std::complex<double>, 2> greenUpNum, greenDnNum;

    tensor_hao::TensorHao<std::complex<double>,2> wfUpDaggerT, wfDnDaggerT;
    tensor_hao::TensorHao<std::complex<double>,3> wfUpDaggerCholeskyVecs, wfDnDaggerCholeskyVecs;

 public:
    RealMaterialMoleculeMeasureObserveSD2sSD2s();
    RealMaterialMoleculeMeasureObserveSD2sSD2s(const RealMaterialMolecule& realMaterialMolecule_);
    ~RealMaterialMoleculeMeasureObserveSD2sSD2s();

    void initModelNullptr();
    void setModel(const RealMaterialMolecule& realMaterialMolecule_);
    void reSet();
    std::complex<double> returnEnergy();
    tensor_hao::TensorHao<double, 1> returnCholeskyBg();
    void addMeasurement(SD2sSD2sOperation &sd2sSD2sOperation, std::complex<double> denIncrement);
    CholeskyRealForce getForce(const CholeskyReal &choleskyReal, SD2sSD2sOperation &sd2sSD2sOperation, double cap=1.0);

    void write() const;
    double getMemory() const;
 private:
    RealMaterialMoleculeMeasureObserveSD2sSD2s(const RealMaterialMoleculeMeasureObserveSD2sSD2s& x);
    RealMaterialMoleculeMeasureObserveSD2sSD2s & operator  = (const RealMaterialMoleculeMeasureObserveSD2sSD2s& x);

    void checkWalkerWithModel(const SD2sSD2sOperation &sd2sSD2sOperation);
    void initWfDaggerT(SD2sSD2sOperation &sd2sSD2sOperation);
    void initWfDaggerCholeskyVecs(SD2sSD2sOperation &sd2sSD2sOperation);
    void addEnergy(SD2sSD2sOperation &sd2sSD2sOperation, std::complex<double> denIncrement);
    void addGreen(SD2sSD2sOperation &sd2sSD2sOperation, std::complex<double> denIncrement);
    std::complex<double> calculateTenergy(SD2sSD2sOperation &sd2sSD2sOperation);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyBg(SD2sSD2sOperation &sd2sSD2sOperation);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyEx(SD2sSD2sOperation &sd2sSD2sOperation);
};

#endif //AFQMCLAB_REALMATERIALMOLECULEOBSERVESD2SSD2S_H