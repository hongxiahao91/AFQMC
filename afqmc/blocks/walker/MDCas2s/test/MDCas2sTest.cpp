//
// Created by Hao Shi on 8/21/17.
//

#include "../include/MDCas2s.h"
#include "../../../../../common/testHao/gtest_custom.h"
#include "../../../../../common/readWriteHao/include/readWriteHao.h"

using namespace std;
using namespace tensor_hao;

class MDCas2sTest: public ::testing::Test
{
 public:
    size_t L, Nup, Ndn;
    size_t MDSize, MDUpSize, MDDnSize, NupMax, NdnMax;
    complex<double> logw, logwExpect;

    vector< complex<double> > coe;
    std::vector< size_t > coeLinkUp, coeLinkDn;
    vector< vector<size_t> > occupancyUp, occupancyDn;
    vector< complex<double> > coeExpect;
    vector< vector<size_t> > occupancyUpExpect, occupancyDnExpect;

    TensorHao<complex<double>,2> wfUp, wfDn;
    TensorHao<complex<double>,2> RUp, RDn;

    vector<int> isParentUp, isParentDn;
    vector<size_t> treeUp, treeDn;
    vector<vector<size_t>> diffOrbitalUp, diffOrbitalDn;
    vector<size_t> parentIndexUp, parentIndexDn;
    
    string filename="MDCas2s.dat";

    MDCas2sTest( )
    {
        L=10; Nup=4; Ndn=3;
        MDSize=12; MDUpSize=8; MDDnSize=6; NupMax=9; NdnMax=7;
        logw=complex<double>(1.2,2.3);
        coe = { {1.0, 2.0}, {1.0, 3.0}, {2.0, 1.0}, {1.0, 1.0}, {2.0, 2.0}, {-1.0, 1.0},
                {3.0, 2.0}, {2.0, 6.0}, {3.0, 5.0}, {4.0, 2.0}, {3.0, 1.0}, { 0.2, 7.0} };
        coeLinkUp = {0, 1, 6, 2, 5, 3, 2, 4, 5, 6, 2, 7};
        coeLinkDn = {0, 1, 2, 4, 5, 3, 2, 0, 1, 2, 5, 2};
        occupancyUp = { {0,1,2,3}, {0,2,3,4}, {0,4,7,8}, {0,1,2,5}, {0,1,3,5}, {0,2,7,8}, {0,2,5,6}, {0,5,6,8} };
        occupancyDn = { {0,1,2}, {0,3,4}, {1,3,6}, {0,2,5}, {1,3,5}, {0,2,3} };
        wfUp.resize(L, NupMax); randomFill(wfUp); MPIBcast(wfUp);
        wfDn.resize(L, NdnMax); randomFill(wfDn); MPIBcast(wfDn);
        RUp.resize(NupMax, NupMax); RUp=complex<double>(0.0); for(size_t i = 0; i < NupMax; ++i) RUp(i, i) = 1.0;
        RDn.resize(NdnMax, NdnMax); RDn=complex<double>(0.0); for(size_t i = 0; i < NdnMax; ++i) RDn(i, i) = 1.0;

        if( MPIRank() == 0 )
        {
            ofstream file;
            file.open(filename, ios::out|ios::trunc);
            if ( ! file.is_open() ) { cout << "Error opening file in File!!! "<<filename<<endl; exit(1); }

            writeFile(L, file);      writeFile(Nup, file);    writeFile(Ndn, file);
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

        coeExpect = { {1.0, 2.0}, {1.0, 3.0}, {2.0, 1.0}, {-1.0, -1.0}, {2.0, 2.0}, {1.0, -1.0},
                      {-3.0, -2.0}, {-2.0, -6.0}, {-3.0, -5.0}, {4.0, 2.0}, {3.0, 1.0}, { 0.2, 7.0} };
        double norm(0.0); for(size_t i = 0; i < MDSize; ++i) norm += ( coeExpect[i] * conj( coeExpect[i] ) ).real();
        norm =  sqrt( norm ); for(size_t i = 0; i < MDSize; ++i) coeExpect[i] = coeExpect[i]/norm;
        logwExpect = logw + log(norm);

        occupancyUpExpect={ {0,1,2,3}, {0,4,2,3}, {0,4,8,7}, {0,1,2,5}, {0,1,5,3}, {0,2,8,7}, {0,6,2,5}, {0,6,8,5} };
        occupancyDnExpect={ {0,1,2}, {0,3,4}, {1,3,6}, {0,5,2}, {1,3,5}, {0,3,2} };

        isParentUp = {1, 0, 0, 1, 0, 0, 1, 1};
        isParentDn = {1, 1, 0, 0, 0, 1};
        treeUp = {0, 1, 4, 3, 6, 7, 2, 5 };
        treeDn = {0, 3, 5, 1, 2, 4 };
        diffOrbitalUp = { {}, {1}, {1,3}, {3}, {2}, {1,3}, {1}, {2} };
        diffOrbitalDn = { {}, {2}, {0,2}, {1}, {0,2}, {1} };
        parentIndexUp={0, 0, 7, 0, 0, 7, 3, 6};
        parentIndexDn={0, 5, 1, 0, 1, 0};
    }

    ~MDCas2sTest( )
    {
        removeFile(filename);
    }
};

TEST_F(MDCas2sTest, readWriteBcast)
{
    string filenamePrime="MDCas2sPrime.dat";
    MDCas2s mdCas2sBase(filename), mdCas2s;

    if( MPIRank()==0 ) mdCas2sBase.write(filenamePrime);
    MPIBarrier();

    if( MPIRank() == 0 ) mdCas2s.read(filenamePrime);
    MPIBcast( mdCas2s );

    EXPECT_EQ( L, mdCas2s.getL() );
    EXPECT_EQ( Nup, mdCas2s.getNup() );
    EXPECT_EQ( Ndn, mdCas2s.getNdn() );
    EXPECT_EQ( MDSize, mdCas2s.getMDSize() );
    EXPECT_EQ( NupMax, mdCas2s.getNupMax() );
    EXPECT_EQ( NdnMax, mdCas2s.getNdnMax() );
    EXPECT_COMPLEXDOUBLE_EQ( logwExpect, mdCas2s.getLogw() );

    EXPECT_POINTER_COMPLEXDOUBLE_EQ( MDSize, coeExpect.data(), mdCas2s.getCoe().data() );
    EXPECT_EQ( occupancyUpExpect, mdCas2s.getOccupancyUp() );
    EXPECT_EQ( occupancyDnExpect, mdCas2s.getOccupancyDn() );

    EXPECT_FALSE( diff( wfUp, mdCas2s.getWfUp(), 1e-12 ) );
    EXPECT_FALSE( diff( wfDn, mdCas2s.getWfDn(), 1e-12 ) );
    EXPECT_FALSE( diff( RUp, mdCas2s.getRUp(), 1e-12 ) );
    EXPECT_FALSE( diff( RDn, mdCas2s.getRDn(), 1e-12 ) );

    EXPECT_EQ( isParentUp, mdCas2s.getIsParentUp() );
    EXPECT_EQ( isParentDn, mdCas2s.getIsParentDn() );
    EXPECT_EQ( treeUp, mdCas2s.getTreeUp() );
    EXPECT_EQ( treeDn, mdCas2s.getTreeDn() );
    EXPECT_EQ( diffOrbitalUp, mdCas2s.getDiffOrbitalUp() );
    EXPECT_EQ( diffOrbitalDn, mdCas2s.getDiffOrbitalDn() );
    EXPECT_EQ( parentIndexUp, mdCas2s.getParentIndexUp() );
    EXPECT_EQ( parentIndexDn, mdCas2s.getParentIndexDn() );

    TensorHao<complex<double>,2> wfRUp(L, NupMax), wfRDn(L, NdnMax);
    gmm_cpu(wfUp, RUp, wfRUp);
    gmm_cpu(wfDn, RDn, wfRDn);

    EXPECT_FALSE( diff(wfRUp, mdCas2s.generateWfCasUp(), 1e-12) );
    EXPECT_FALSE( diff(wfRDn, mdCas2s.generateWfCasDn(), 1e-12) );

    removeFile(filenamePrime);
}

TEST_F(MDCas2sTest, stabilize)
{
    MDCas2s mdCas2s(filename);
    const vector< complex<double> > &coeRef = mdCas2s.getCoe();
    const TensorHao<complex<double>,2> &wfRUp = mdCas2s.generateWfCasUp();
    const TensorHao<complex<double>,2> &wfRDn = mdCas2s.generateWfCasDn();

    TensorHao<complex<double>,2> walkerUp(L, Nup), walkerDn(L, Ndn);
    TensorHao<complex<double>,2> phiTUp(L, Nup), phiTDn(L, Ndn);
    TensorHao<complex<double>,2> ovlpUp(Nup, Nup), ovlpDn(Ndn, Ndn);
    size_t linkIndex;

    randomFill(walkerUp); randomFill(walkerDn);

    complex<double> ovlpBefore(0,0);
    for(size_t i = 0; i < MDSize; ++i)
    {
        linkIndex = coeLinkUp[i];
        for(size_t j = 0; j < Nup; ++j)
        {
            phiTUp[j].copy_deep(  wfRUp[ occupancyUp[ linkIndex ][j] ] );
        }
        gmm_cpu( phiTUp, walkerUp, ovlpUp, 'C' );
        complex<double> detUp = determinant( LUconstruct_cpu( ovlpUp ) );

        linkIndex = coeLinkDn[i];
        for(size_t j = 0; j < Ndn; ++j)
        {
            phiTDn[j].copy_deep( wfRDn[ occupancyDn[ linkIndex ][j] ] );
        }
        gmm_cpu( phiTDn, walkerDn, ovlpDn, 'C' );
        complex<double> detDn = determinant( LUconstruct_cpu( ovlpDn ) );

        ovlpBefore += detUp*detDn*conj( coeRef[i] * exp( mdCas2s.getLogw() ) );
    }

    mdCas2s.stabilize();

    complex<double> ovlpAfter(0,0);
    for(size_t i = 0; i < MDSize; ++i)
    {
        linkIndex = coeLinkUp[i];
        for(size_t j = 0; j < Nup; ++j)
        {
            phiTUp[j].copy_deep(  wfRUp[ occupancyUp[ linkIndex ] [j] ] );
        }
        gmm_cpu( phiTUp, walkerUp, ovlpUp, 'C' );
        complex<double> detUp = determinant( LUconstruct_cpu( ovlpUp ) );

        linkIndex = coeLinkDn[i];
        for(size_t j = 0; j < Ndn; ++j)
        {
            phiTDn[j].copy_deep( wfRDn[ occupancyDn[ linkIndex ][j] ] );
        }
        gmm_cpu( phiTDn, walkerDn, ovlpDn, 'C' );
        complex<double> detDn = determinant( LUconstruct_cpu( ovlpDn ) );

        ovlpAfter += detUp*detDn*conj( coeRef[i] * exp( mdCas2s.getLogw() ) );
    }
    EXPECT_COMPLEX_NEAR(ovlpBefore, ovlpAfter, 1e-9);
}
