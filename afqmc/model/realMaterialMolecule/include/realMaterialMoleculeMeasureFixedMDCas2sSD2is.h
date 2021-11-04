//
// Created by Hao Shi on 8/31/17.
//

#ifndef AFQMCLAB_REALMATERIALMOLECULEMEASUREFIXEDMDCAS2SSD2IS_H
#define AFQMCLAB_REALMATERIALMOLECULEMEASUREFIXEDMDCAS2SSD2IS_H

#include "realMaterialMolecule.h"
#include "../../../blocks/walkerWalkerOperation/MDCas2sSD2isOperation/include/MDCas2sSD2isOperation.h"

class RealMaterialMoleculeMeasureFixedMDCas2sSD2is
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
    RealMaterialMoleculeMeasureFixedMDCas2sSD2is();
    RealMaterialMoleculeMeasureFixedMDCas2sSD2is(const RealMaterialMolecule& realMaterialMolecule_, MDCas2s &walkerLeft_);
    ~RealMaterialMoleculeMeasureFixedMDCas2sSD2is();

    void initModelWalkerNullptr();
    void setModelWalker(const RealMaterialMolecule& realMaterialMolecule_, MDCas2s &walkerLeft_);
    void reSet();
    std::complex<double> returnEnergy();
    tensor_hao::TensorHao<double, 1> returnCholeskyBg();
    void addMeasurement(MDCas2sSD2isOperation &mdCas2sSD2isOperation, std::complex<double> denIncrement);
    CholeskyRealForce getForce(const CholeskyReal &choleskyReal, MDCas2sSD2isOperation &mdCas2sSD2isOperation, double cap=1.0);

    void write() const;
    double getMemory() const;
 private:
    RealMaterialMoleculeMeasureFixedMDCas2sSD2is(const RealMaterialMoleculeMeasureFixedMDCas2sSD2is& x);
    RealMaterialMoleculeMeasureFixedMDCas2sSD2is & operator  = (const RealMaterialMoleculeMeasureFixedMDCas2sSD2is& x);

    void initWfDaggerT_T();
    void initWfDaggerCholeskyVecs_T();
    void checkWalkerLeft(const MDCas2sSD2isOperation &mdCas2sSD2isOperation);
    void addEnergy(MDCas2sSD2isOperation &mdCas2sSD2isOperation, std::complex<double> numIncrement);

    std::complex<double> calculateTenergyUp(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaUp);
    std::complex<double> calculateTenergyDn(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaDn);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyBgUp(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaUp);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyBgDn(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaDn);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyExUp(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaUp);
    tensor_hao::TensorHao<std::complex<double>, 1> calculateCholeskyExDn(size_t casIndex, const tensor_hao::TensorHao<std::complex<double>,2> &thetaDn);
};

#endif //AFQMCLAB_REALMATERIALMOLECULEMEASUREFIXEDMDCAS2SSD2IS_H
