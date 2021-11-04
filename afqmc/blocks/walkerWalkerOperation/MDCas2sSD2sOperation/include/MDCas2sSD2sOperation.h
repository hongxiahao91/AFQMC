//
// Created by Hao Shi on 10/4/17.
//

#ifndef AFQMCLAB_MDCAS2SSD2SOPERATION_H
#define AFQMCLAB_MDCAS2SSD2SOPERATION_H

#include "../../../walker/MDCas2s/include/MDCas2s.h"
#include "../../../walker/SD2s/include/SD2s.h"

class MDCas2sSD2sOperation
{
    MDCas2s *walkerLeft;
    SD2s *walkerRight;
    size_t stablizeStep;

    tensor_hao::TensorHao<std::complex<double>, 2> ovlpFullUp, ovlpFullDn;
    std::vector< std::complex<double> > detUp, detDn;
    std::vector< tensor_hao::TensorHao<std::complex<double>, 2> > ovlpInvUp, ovlpInvDn, QUp, QDn;

    std::complex<double> logOverlap; bool logOverlapIsCalculated;
    std::vector< tensor_hao::TensorHao<std::complex<double>, 2> > thetaUp, thetaDn; bool thetaIsCalculated;
 public:
    MDCas2sSD2sOperation();
    MDCas2sSD2sOperation(MDCas2s &walkerLeft_, SD2s &walkerRight_, size_t stablizeStep_=200);
    ~MDCas2sSD2sOperation();

    const MDCas2s *getWalkerLeft() const;
    const SD2s *getWalkerRight() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getOvlpFullUp() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getOvlpFullDn() const;
    const std::vector<std::complex<double>> &getDetUp() const;
    const std::vector<std::complex<double>> &getDetDn() const;
    const std::vector<tensor_hao::TensorHao<std::complex<double>, 2>> &getOvlpInvUp() const;
    const std::vector<tensor_hao::TensorHao<std::complex<double>, 2>> &getOvlpInvDn() const;
    const std::vector<tensor_hao::TensorHao<std::complex<double>, 2>> &getQUp() const;
    const std::vector<tensor_hao::TensorHao<std::complex<double>, 2>> &getQDn() const;

    void set(MDCas2s &walkerLeft_, SD2s &walkerRight_, size_t stablizeStep_=200);
    void reSet();

    std::complex<double> returnLogOverlap();
    const std::vector< tensor_hao::TensorHao<std::complex<double>, 2> > &returnThetaUp();
    const std::vector< tensor_hao::TensorHao<std::complex<double>, 2> > &returnThetaDn();
    double getMemory() const;

 private:
    MDCas2sSD2sOperation(const MDCas2sSD2sOperation& x);
    MDCas2sSD2sOperation & operator  = (const MDCas2sSD2sOperation& x);

    void setDetOvlpInvQForUpSpin();
    void setDetOvlpInvQForDnSpin();

    void calculateTheta();
    void calculateThetaUp();
    void calculateThetaDn();
};

void setWalkerFromPhiT(std::vector<SD2s> &walker, std::vector<bool> &walkerIsAlive, MDCas2s &phiT, double noise=1e-5);

#endif //AFQMCLAB_MDCAS2SSD2SOPERATION_H