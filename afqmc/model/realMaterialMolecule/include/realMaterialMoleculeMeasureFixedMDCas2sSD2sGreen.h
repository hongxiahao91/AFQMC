//
// Created by Hao Shi on 8/31/17.
//

#ifndef AFQMCLAB_REALMATERIALMOLECULEMEASUREFIXEDMDCAS2SSD2SGREEN_H
#define AFQMCLAB_REALMATERIALMOLECULEMEASUREFIXEDMDCAS2SSD2SGREEN_H

#include "realMaterialMolecule.h"
#include "../../../blocks/walkerWalkerOperation/MDCas2sSD2sOperation/include/MDCas2sSD2sOperation.h"

class RealMaterialMoleculeMeasureFixedMDCas2sSD2sGreen
{
 private:
    const RealMaterialMolecule *realMaterialMolecule;
    MDCas2s  *walkerLeft;
    tensor_hao::TensorHao<std::complex<double>, 2> wfCasUpConj, wfCasDnConj;

    std::complex<double> den;
    std::complex<double> TNum;
    tensor_hao::TensorHao<std::complex<double>, 1> choleskyBgNum;
    tensor_hao::TensorHao<std::complex<double>, 1> choleskyExNum;
    std::complex<double> HNum;

    tensor_hao::TensorHao<std::complex<double>,2> wfUpDaggerT_T, wfDnDaggerT_T;
    tensor_hao::TensorHao<std::complex<double>,3> wfUpDaggerCholeskyVecs_T, wfDnDaggerCholeskyVecs_T;

 public:
    RealMaterialMoleculeMeasureFixedMDCas2sSD2sGreen();
    RealMaterialMoleculeMeasureFixedMDCas2sSD2sGreen(const RealMaterialMolecule& realMaterialMolecule_, MDCas2s &walkerLeft_);
    ~RealMaterialMoleculeMeasureFixedMDCas2sSD2sGreen();

    void initModelWalkerNullptr();
    void setModelWalker(const RealMaterialMolecule& realMaterialMolecule_, MDCas2s &walkerLeft_);
    void reSet();
    std::complex<double> returnEnergy();
    tensor_hao::TensorHao<double, 1> returnCholeskyBg();
    void addMeasurement(MDCas2sSD2sOperation &mdCas2sSD2sOperation, std::complex<double> denIncrement);
    CholeskyRealForce getForce(const CholeskyReal &choleskyReal, MDCas2sSD2sOperation &mdCas2sSD2sOperation, double cap=1.0);

    void write() const;
    double getMemory() const;
 private:
    RealMaterialMoleculeMeasureFixedMDCas2sSD2sGreen(const RealMaterialMoleculeMeasureFixedMDCas2sSD2sGreen& x);
    RealMaterialMoleculeMeasureFixedMDCas2sSD2sGreen & operator  = (const RealMaterialMoleculeMeasureFixedMDCas2sSD2sGreen& x);

    void initWfDaggerT_T();
    void initWfDaggerCholeskyVecs_T();
    void checkWalkerLeft(const MDCas2sSD2sOperation &mdCas2sSD2sOperation);
    void addEnergy(MDCas2sSD2sOperation &mdCas2sSD2sOperation, std::complex<double> numIncrement);

    std::complex<double> calculateTenergyUp(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaUp);
    std::complex<double> calculateTenergyDn(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaDn);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyBgUp(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaUp);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyBgDn(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaDn);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyExUp(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaUp);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyExDn(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaDn);
};

#endif //AFQMCLAB_REALMATERIALMOLECULEMEASUREFIXEDMDCAS2SSD2SGREEN_H
