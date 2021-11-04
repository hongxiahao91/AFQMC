//
// Created by Hao Shi on 8/29/17.
//
#include "../include/MDCas2sSD2isOperation.h"

using namespace std;
using namespace tensor_hao;

class MDCas2sSD2isOperationTest: public ::testing::Test
{
 public:
    size_t L, Nup, Ndn;
    size_t MDSize, MDUpSize, MDDnSize, NupMax, NdnMax;

    string filename="MDCas2s.dat";

    MDCas2sSD2isOperationTest( )
    {
        L=10; Nup=4; Ndn=3;
        MDSize=12; MDUpSize=8; MDDnSize=6; NupMax=9; NdnMax=7;

        complex<double> logw(1.2,2.3);

        vector< complex<double> > coe = { {1.0, 2.0}, {1.0, 3.0}, {2.0, 1.0}, {1.0, 1.0}, {2.0, 2.0}, {-1.0, 1.0},
                                          {3.0, 2.0}, {2.0, 6.0}, {3.0, 5.0}, {4.0, 2.0}, {3.0, 1.0}, { 0.2, 7.0} };
        std::vector< size_t > coeLinkUp = {0, 1, 6, 2, 5, 3, 2, 4, 5, 6, 2, 7};
        std::vector< size_t > coeLinkDn = {0, 1, 2, 4, 5, 3, 2, 0, 1, 2, 5, 2};
        vector< vector<size_t> > occupancyUp = { {0,1,2,3}, {0,2,3,4}, {0,4,7,8}, {0,1,2,5}, {0,1,3,5}, {0,2,7,8}, {0,2,5,6}, {0,5,6,8} };
        vector< vector<size_t> > occupancyDn = { {0,1,2}, {0,3,4}, {1,3,6}, {0,2,5}, {1,3,5}, {0,2,3} };
        TensorHao<complex<double>,2> wfUp, wfDn;
        TensorHao<complex<double>,2> RUp, RDn;

        wfUp.resize(L, NupMax); randomFill(wfUp); MPIBcast(wfUp);
        wfDn.resize(L, NdnMax); randomFill(wfDn); MPIBcast(wfDn);
        RUp.resize(NupMax, NupMax); RUp=complex<double>(0.0); for(size_t i = 0; i < NupMax; ++i) RUp(i, i) = 1.0;
        RDn.resize(NdnMax, NdnMax); RDn=complex<double>(0.0); for(size_t i = 0; i < NdnMax; ++i) RDn(i, i) = 1.0;

        if( MPIRank() == 0 )
        {
            ofstream file;
            file.open(filename, ios::out|ios::trunc);
            if ( ! file.is_open() ) { cout << "Error opening file in File!!! "<<filename<<endl; exit(1); }

            writeFile(L, file);      writeFile(Nup, file);      writeFile(Ndn, file);
            writeFile(MDSize, file); writeFile(MDUpSize, file); writeFile(MDDnSize, file); writeFile(NupMax, file); writeFile(NdnMax, file);
            writeFile(logw, file);
            writeFile(coe.size(), coe.data(), file);
            writeFile(coeLinkUp.size(), coeLinkUp.data(), file);
            writeFile(coeLinkDn.size(), coeLinkDn.data(), file);
            for(size_t i = 0; i < MDUpSize; ++i) writeFile(Nup, occupancyUp[i].data(), file);
            for(size_t i = 0; i < MDDnSize; ++i) writeFile(Ndn, occupancyDn[i].data(), file);
            writeFile(wfUp.size(), wfUp.data(), file);
            writeFile(wfDn.size(), wfDn.data(), file);
            writeFile(RUp.size(), RUp.data(), file);
            writeFile(RDn.size(), RDn.data(), file);

            file.close();
        }
        MPIBarrier();
    }

    ~MDCas2sSD2isOperationTest( )
    {
        removeFile(filename);
    }
};

TEST_F( MDCas2sSD2isOperationTest, logOverlap )
{
    MDCas2s mdCas2s(filename);
    SD2is sd2is(L, Nup, Ndn); randomFill( sd2is.wfRef() ); sd2is.logwRef()=complex<double>(2.0, 3.0);

    MDCas2sSD2isOperation walkerWalkerOperation( mdCas2s, sd2is );

    const TensorHao<complex<double>,2> &wfUpRef = mdCas2s.getWfUp();
    const TensorHao<complex<double>,2> &wfDnRef = mdCas2s.getWfDn();
    const TensorHao<complex<double>,2> &RUpRef = mdCas2s.getRUp();
    const TensorHao<complex<double>,2> &RDnRef = mdCas2s.getRDn();
    const vector< complex<double> > &coeRef = mdCas2s.getCoe();
    const vector< size_t > &coeLinkUpRef = mdCas2s.getCoeLinkUp();
    const vector< size_t > &coeLinkDnRef = mdCas2s.getCoeLinkDn();
    const vector< vector<size_t> > &occupancyUpRef = mdCas2s.getOccupancyUp();
    const vector< vector<size_t> > &occupancyDnRef = mdCas2s.getOccupancyDn();

    TensorHao<complex<double>,2> wfRUp(L, NupMax), wfRDn(L, NdnMax);
    TensorHao<complex<double>,2> walkerUp(L, Nup), walkerDn(L, Ndn);
    TensorHao<complex<double>,2> phiTUp(L, Nup), phiTDn(L, Ndn);
    TensorHao<complex<double>,2> ovlpUp(Nup, Nup), ovlpDn(Ndn, Ndn);

    gmm_cpu(wfUpRef, RUpRef, wfRUp);
    gmm_cpu(wfDnRef, RDnRef, wfRDn);

    copy( sd2is.getWf().data(),sd2is.getWf().data()+L*Nup, walkerUp.data() );
    copy( sd2is.getWf().data(),sd2is.getWf().data()+L*Ndn, walkerDn.data() );

    size_t linkIndex;
    complex<double> ovlpExpect(0,0);

    for(size_t i = 0; i < MDSize; ++i)
    {
        linkIndex = coeLinkUpRef[i];
        for(size_t j = 0; j < Nup; ++j)
        {
            phiTUp[j].copy_deep(  wfRUp[ occupancyUpRef[ linkIndex ][j] ] );
        }
        gmm_cpu( phiTUp, walkerUp, ovlpUp, 'C' );
        complex<double> detUp = determinant( LUconstruct_cpu( ovlpUp ) );

        linkIndex = coeLinkDnRef[i];
        for(size_t j = 0; j < Ndn; ++j)
        {
            phiTDn[j].copy_deep( wfRDn[ occupancyDnRef[ linkIndex ][j] ] );
        }
        gmm_cpu( phiTDn, walkerDn, ovlpDn, 'C' );
        complex<double> detDn = determinant( LUconstruct_cpu( ovlpDn ) );

        ovlpExpect += detUp*detDn*conj( coeRef[i] * exp( mdCas2s.getLogw() ) ) * exp( sd2is.getLogw() );
    }

    EXPECT_COMPLEX_NEAR(ovlpExpect, exp( walkerWalkerOperation.returnLogOverlap() ), abs(1e-14*ovlpExpect) );
}

TEST_F( MDCas2sSD2isOperationTest, theta )
{
    MDCas2s mdCas2s(filename);
    SD2is sd2is(L, Nup, Ndn); randomFill( sd2is.wfRef() ); sd2is.logwRef()=complex<double>(2.0, 3.0);

    MDCas2sSD2isOperation walkerWalkerOperation( mdCas2s, sd2is );

    const TensorHao<complex<double>,2> &wfUpRef = mdCas2s.getWfUp();
    const TensorHao<complex<double>,2> &wfDnRef = mdCas2s.getWfDn();
    const TensorHao<complex<double>,2> &RUpRef = mdCas2s.getRUp();
    const TensorHao<complex<double>,2> &RDnRef = mdCas2s.getRDn();
    const vector< vector<size_t> > &occupancyUpRef = mdCas2s.getOccupancyUp();
    const vector< vector<size_t> > &occupancyDnRef = mdCas2s.getOccupancyDn();
    const vector< TensorHao<complex<double>, 2> > &thetaUpVec = walkerWalkerOperation.returnThetaUp();
    const vector< TensorHao<complex<double>, 2> > &thetaDnVec = walkerWalkerOperation.returnThetaDn();

    TensorHao<complex<double>,2> wfRUp(L, NupMax), wfRDn(L, NdnMax);
    TensorHao<complex<double>,2> walkerUp(L, Nup), walkerDn(L, Ndn);
    TensorHao<complex<double>,2> phiTUp(L, Nup), phiTDn(L, Ndn);
    TensorHao<complex<double>,2> ovlpUp(Nup, Nup), ovlpDn(Ndn, Ndn);
    TensorHao<complex<double>,2> thetaUp(L, Nup), thetaDn(L, Ndn);

    gmm_cpu(wfUpRef, RUpRef, wfRUp);
    gmm_cpu(wfDnRef, RDnRef, wfRDn);

    copy( sd2is.getWf().data(),sd2is.getWf().data()+L*Nup, walkerUp.data() );
    copy( sd2is.getWf().data(),sd2is.getWf().data()+L*Ndn, walkerDn.data() );

    for(size_t i = 0; i < MDUpSize; ++i)
    {
        for(size_t j = 0; j < Nup; ++j)
        {
            phiTUp[j].copy_deep(wfRUp[occupancyUpRef[i][j]]);
        }
        gmm_cpu(phiTUp, walkerUp, ovlpUp, 'C');
        ovlpUp = inverse_cpu(LUconstruct_cpu(ovlpUp));
        gmm_cpu(walkerUp, ovlpUp, thetaUp);
        EXPECT_FALSE(diff(thetaUp, thetaUpVec[i], 1e-10));
    }

    for(size_t i = 0; i < MDDnSize; ++i)
    {
        for(size_t j = 0; j < Ndn; ++j)
        {
            phiTDn[j].copy_deep( wfRDn[ occupancyDnRef[i][j] ] );
        }
        gmm_cpu( phiTDn, walkerDn, ovlpDn, 'C' );
        ovlpDn = inverse_cpu( LUconstruct_cpu( ovlpDn ) );
        gmm_cpu(walkerDn, ovlpDn, thetaDn);
        EXPECT_FALSE( diff( thetaDn, thetaDnVec[i], 1e-10 ) );
    }
}