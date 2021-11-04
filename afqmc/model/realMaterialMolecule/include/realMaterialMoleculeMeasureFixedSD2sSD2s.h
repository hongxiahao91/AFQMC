//
// Created by boruoshihao on 11/14/18.
//

#ifndef AFQMCLAB_REALMATERIALMOLECULEFIXEDSD2SSD2S_H
#define AFQMCLAB_REALMATERIALMOLECULEFIXEDSD2SSD2S_H

#include "realMaterialMolecule.h"
#include "../../../blocks/walkerWalkerOperation/SD2sSD2sOperation/include/SD2sSD2sOperation.h"

class RealMaterialMoleculeMeasureFixedSD2sSD2s
{
 private:
    const RealMaterialMolecule *realMaterialMolecule;
    const SD2s  *walkerLeft;

    std::complex<double> den;
    std::complex<double> TNum;
    tensor_hao::TensorHao<std::complex<double>, 1> choleskyBgNum;
    tensor_hao::TensorHao<std::complex<double>, 1> choleskyExNum;
    std::complex<double> HNum;
    double HImag; long HImagAdd;
    tensor_hao::TensorHao<std::complex<double>, 2> greenUpNum, greenDnNum;

    tensor_hao::TensorHao<std::complex<double>,2> wfUpDaggerT, wfDnDaggerT;
    tensor_hao::TensorHao<std::complex<double>,3> wfUpDaggerCholeskyVecs, wfDnDaggerCholeskyVecs;

 public:
    RealMaterialMoleculeMeasureFixedSD2sSD2s();
    RealMaterialMoleculeMeasureFixedSD2sSD2s(const RealMaterialMolecule& realMaterialMolecule_, const SD2s &walkerLeft_);
    ~RealMaterialMoleculeMeasureFixedSD2sSD2s();

    void initModelWalkerNullptr();
    void setModelWalker(const RealMaterialMolecule& realMaterialMolecule_, const SD2s &walkerLeft_);
    void reSet();
    std::complex<double> returnEnergy();
    tensor_hao::TensorHao<double, 1> returnCholeskyBg();
    void addMeasurement(SD2sSD2sOperation &sd2sSD2sOperation, std::complex<double> denIncrement);
    CholeskyRealForce getForce(const CholeskyReal &choleskyReal, SD2sSD2sOperation &sd2sSD2sOperation, double cap=1.0);

    void write() const;
    double getMemory() const;
 private:
    RealMaterialMoleculeMeasureFixedSD2sSD2s(const RealMaterialMoleculeMeasureFixedSD2sSD2s& x);
    RealMaterialMoleculeMeasureFixedSD2sSD2s & operator  = (const RealMaterialMoleculeMeasureFixedSD2sSD2s& x);

    void initWfDaggerT();
    void initWfDaggerCholeskyVecs();
    void checkWalkerLeft(const SD2sSD2sOperation &sd2sSD2sOperation);
    void addEnergy(SD2sSD2sOperation &sd2sSD2sOperation, std::complex<double> denIncrement);
    void addGreen(SD2sSD2sOperation &sd2sSD2sOperation, std::complex<double> denIncrement);
    std::complex<double> calculateTenergy(SD2sSD2sOperation &sd2sSD2sOperation);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyBg(SD2sSD2sOperation &sd2sSD2sOperation);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyEx(SD2sSD2sOperation &sd2sSD2sOperation);
};

#endif //AFQMCLAB_REALMATERIALMOLECULEFIXEDSD2SSD2S_H