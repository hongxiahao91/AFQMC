//
// Created by Hao Shi on 8/29/17.
//

#ifndef AFQMCLAB_MDCAS2SSD2ISOPERATION_H
#define AFQMCLAB_MDCAS2SSD2ISOPERATION_H

#include "../../../walker/MDCas2s/include/MDCas2s.h"
#include "../../../walker/SD2is/include/SD2is.h"

class MDCas2sSD2isOperation
{
    MDCas2s *walkerLeft;
    SD2is *walkerRight;
    size_t stablizeStep;
    tensor_hao::TensorHaoRef<std::complex<double>, 2> wfRightUp, wfRightDn;
    tensor_hao::TensorHao<std::complex<double>, 2> ovlpFullUp, ovlpFullDn;
    std::vector< std::complex<double> > detUp, detDn;
    std::vector< tensor_hao::TensorHao<std::complex<double>, 2> > ovlpInvUp, ovlpInvDn, QUp, QDn;

    std::complex<double> logOverlap; bool logOverlapIsCalculated;
    std::vector< tensor_hao::TensorHao<std::complex<double>, 2> > thetaUp, thetaDn; bool thetaIsCalculated;
 public:
    MDCas2sSD2isOperation();
    MDCas2sSD2isOperation(MDCas2s &walkerLeft_, SD2is &walkerRight_, size_t stablizeStep_=200);
    ~MDCas2sSD2isOperation();

    const MDCas2s *getWalkerLeft() const;
    const SD2is *getWalkerRight() const;
    const tensor_hao::TensorHaoRef<std::complex<double>, 2> &getWfRightUp() const;
    const tensor_hao::TensorHaoRef<std::complex<double>, 2> &getWfRightDn() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getOvlpFullUp() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getOvlpFullDn() const;
    const std::vector<std::complex<double>> &getDetUp() const;
    const std::vector<std::complex<double>> &getDetDn() const;
    const std::vector<tensor_hao::TensorHao<std::complex<double>, 2>> &getOvlpInvUp() const;
    const std::vector<tensor_hao::TensorHao<std::complex<double>, 2>> &getOvlpInvDn() const;
    const std::vector<tensor_hao::TensorHao<std::complex<double>, 2>> &getQUp() const;
    const std::vector<tensor_hao::TensorHao<std::complex<double>, 2>> &getQDn() const;

    void set(MDCas2s &walkerLeft_, SD2is &walkerRight_, size_t stablizeStep_=200);
    void reSet();

    std::complex<double> returnLogOverlap();
    const std::vector< tensor_hao::TensorHao<std::complex<double>, 2> > &returnThetaUp();
    const std::vector< tensor_hao::TensorHao<std::complex<double>, 2> > &returnThetaDn();
    double getMemory() const;

 private:
    MDCas2sSD2isOperation(const MDCas2sSD2isOperation& x);
    MDCas2sSD2isOperation & operator  = (const MDCas2sSD2isOperation& x);

    void setDetOvlpInvQForUpSpin();
    void setDetOvlpInvQForDnSpin();

    void calculateTheta();
    void calculateThetaUp();
    void calculateThetaDn();
};

void setWalkerFromPhiT(std::vector<SD2is> &walker, std::vector<bool> &walkerIsAlive, const MDCas2s &phiT, double noise=1e-5);

#endif //AFQMCLAB_MDCAS2SSD2ISOPERATION_H
