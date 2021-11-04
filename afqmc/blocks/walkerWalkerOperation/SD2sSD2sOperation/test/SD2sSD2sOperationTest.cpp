//
// Created by boruoshihao on 11/14/18.
//
#include "../include/SD2sSD2sOperation.h"
#include "../../../../../common/testHao/gtest_custom.h"

using namespace std;
using namespace tensor_hao;

class SD2sSD2sOperationTest: public ::testing::Test
{
 public:
    size_t L, Nup, Ndn;
    TensorHao<complex<double>,2> wfLeftUp, wfLeftDn, wfRightUp, wfRightDn;
    SD2s walkerLeft;
    SD2s walkerRight;

    SD2sSD2sOperationTest( )
    {
        L=10; Nup=3; Ndn=5;

        wfLeftUp.resize(L,Nup); wfLeftDn.resize(L,Ndn); wfRightUp.resize(L,Nup); wfRightDn.resize(L,Ndn);
        randomFill(wfLeftUp); randomFill(wfLeftDn); randomFill(wfRightUp); randomFill(wfRightDn);

        walkerLeft.resize(L, Nup, Ndn);
        walkerLeft.logwRef() = complex<double>(1.2, 2.0);
        walkerLeft.wfUpRef() = wfLeftUp;
        walkerLeft.wfDnRef() = wfLeftDn;

        walkerRight.resize(L, Nup, Ndn);
        walkerRight.logwRef() = complex<double>(2.2, 3.0);
        walkerRight.wfUpRef() = wfRightUp;
        walkerRight.wfDnRef() = wfRightDn;
    }

    ~SD2sSD2sOperationTest( )  {}
};

TEST_F(SD2sSD2sOperationTest, voidConstruction)
{
    SD2sSD2sOperation sd2sSD2sOperation;
    EXPECT_EQ( sd2sSD2sOperation.getState(), SD2sSD2sOperationState::VOID );
    EXPECT_FALSE( sd2sSD2sOperation.getWalkerLeft() );
    EXPECT_FALSE( sd2sSD2sOperation.getWalkerRight() );
}

TEST_F(SD2sSD2sOperationTest, walkerPointerConstruction)
{
    SD2sSD2sOperation sd2sSD2sOperation(walkerLeft, walkerRight);

    EXPECT_EQ( sd2sSD2sOperation.getState(), SD2sSD2sOperationState::VOID );
    EXPECT_TRUE( sd2sSD2sOperation.getWalkerLeft() );
    EXPECT_TRUE( sd2sSD2sOperation.getWalkerRight() );
}

TEST_F(SD2sSD2sOperationTest, reset)
{
    SD2sSD2sOperation sd2sSD2sOperation(walkerLeft, walkerRight);

    EXPECT_EQ( sd2sSD2sOperation.getState(), SD2sSD2sOperationState::VOID );
    sd2sSD2sOperation.returnLUOverlapUp();
    EXPECT_EQ( sd2sSD2sOperation.getState(), SD2sSD2sOperationState::LUOVERLAP );
    sd2sSD2sOperation.returnLUOverlapDn();
    EXPECT_EQ( sd2sSD2sOperation.getState(), SD2sSD2sOperationState::LUOVERLAP );
    sd2sSD2sOperation.returnThetaUp_T();
    EXPECT_EQ( sd2sSD2sOperation.getState(), SD2sSD2sOperationState::THETA_T );
    sd2sSD2sOperation.returnThetaDn_T();
    EXPECT_EQ( sd2sSD2sOperation.getState(), SD2sSD2sOperationState::THETA_T );
    sd2sSD2sOperation.reSet();
    EXPECT_EQ( sd2sSD2sOperation.getState(), SD2sSD2sOperationState::VOID );
}

TEST_F(SD2sSD2sOperationTest, getLogOverlap)
{
    complex<double> logOverlap = complex<double>(1.2, -2.0) + complex<double>(2.2, 3.0);
    TensorHao< complex<double>,2 > ovlpUp(Nup,Nup); gmm_cpu(wfLeftUp, wfRightUp, ovlpUp, 'C');
    TensorHao< complex<double>,2 > ovlpDn(Ndn,Ndn); gmm_cpu(wfLeftDn, wfRightDn, ovlpDn, 'C');
    logOverlap += logDeterminant( LUconstruct_cpu(ovlpUp) )+logDeterminant( LUconstruct_cpu(ovlpDn) );

    SD2sSD2sOperation sd2sSD2sOperation(walkerLeft, walkerRight);
    EXPECT_COMPLEXDOUBLE_EQ( logOverlap, sd2sSD2sOperation.returnLogOverlap()  );
}

TEST_F(SD2sSD2sOperationTest, returnGreenMatrixAndDiagonalUp)
{
    TensorHao< complex<double>,2 > ovlpUp(Nup,Nup), wfLeftUpDagger(Nup,L), greenMatrixUp(L,L);
    gmm_cpu(wfLeftUp, wfRightUp, ovlpUp, 'C');
    ovlpUp = inverse_cpu(  LUconstruct_cpu( move(ovlpUp) ) );
    gmm_cpu( ovlpUp, wfLeftUp, wfLeftUpDagger, 'N', 'C' );
    gmm_cpu( wfRightUp, wfLeftUpDagger, greenMatrixUp );
    greenMatrixUp = trans( greenMatrixUp );

    TensorHao< complex<double>,1 > greenDiagonalUp(L);
    for (size_t i = 0; i < L; ++i) greenDiagonalUp(i) = greenMatrixUp(i,i);

    SD2sSD2sOperation sd2sSD2sOperation(walkerLeft, walkerRight);
    EXPECT_FALSE( diff(greenMatrixUp, sd2sSD2sOperation.returnGreenMatrixUp(), 1e-12) );
    EXPECT_FALSE( diff(greenDiagonalUp, sd2sSD2sOperation.returnGreenDiagonalUp(), 1e-12) );

}

TEST_F(SD2sSD2sOperationTest, returnGreenMatrixAndDiagonalDn)
{
    TensorHao< complex<double>,2 > ovlpDn(Ndn,Ndn),  wfLeftDnDagger(Ndn,L), greenMatrixDn(L,L);
    gmm_cpu(wfLeftDn, wfRightDn, ovlpDn, 'C');
    ovlpDn = inverse_cpu(  LUconstruct_cpu( move(ovlpDn) ) );
    gmm_cpu( ovlpDn, wfLeftDn, wfLeftDnDagger, 'N', 'C' );
    gmm_cpu( wfRightDn, wfLeftDnDagger, greenMatrixDn );
    greenMatrixDn = trans( greenMatrixDn );

    TensorHao< complex<double>,1 > greenDiagonalDn(L);
    for (size_t i = 0; i < L; ++i) greenDiagonalDn(i) = greenMatrixDn(i,i);

    SD2sSD2sOperation sd2sSD2sOperation(walkerLeft, walkerRight);
    EXPECT_FALSE( diff(greenMatrixDn, sd2sSD2sOperation.returnGreenMatrixDn(), 1e-12) );
    EXPECT_FALSE( diff(greenDiagonalDn, sd2sSD2sOperation.returnGreenDiagonalDn(), 1e-12) );
}