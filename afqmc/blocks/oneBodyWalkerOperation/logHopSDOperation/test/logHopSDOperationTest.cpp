//
// Created by Hao Shi on 1/13/18.
//
#include "../include/logHopSDOperation.h"
#include "../../../../../common/testHao/gtest_custom.h"

using namespace std;
using namespace tensor_hao;

class LogHopSDOperationTest: public ::testing::Test
{
 public:
    size_t L, N;
    TensorHao<complex<double>,2> matrix, wfOld, wfRightNew, wfLeftNew;


    LogHopSDOperationTest( )
    {
        L=3; N=2;
        matrix.resize(L,L);
        wfOld.resize(L, N);
        wfRightNew.resize(L, N);
        wfLeftNew.resize(L, N);

        matrix = {
                {2,  3}, {4,   2}, {1,   1},
                {1, -1}, {0, 0.3}, {1,  -1},
                {3,  2}, {2,   3}, {0, 0.5}
        };
        matrix = complex<double>(-0.01, 0.0) * matrix;

        wfOld = { {1, 0}, {2,0}, {3, 0}, {3,0}, {4,0}, {5, 0} };

        wfRightNew = { {0.8717255743431891,-0.06695756061260984}, {1.9019539751922434,-0.11272541520548454}, {2.9713563453236675, -0.0038676809054978764},
                       {2.753351116138562, -0.14393925193667553}, {3.7836486105153155,-0.21511811236274694}, {4.932649990968739, -0.012861814634312912} };
        wfLeftNew  = { {0.8708226627963166, 0.09656591603485742}, {1.961995737773784, -0.034019035572233965},{2.9316294741182647, 0.09182786294989435},
                       {2.7315235332357157, 0.21256279549696155}, {3.924363675467334, -0.0680085413131974},  {4.833174488748831, 0.19818901143889675} };

    }

    ~LogHopSDOperationTest( )  {}
};

TEST_F(LogHopSDOperationTest, fixedOrderNumber)
{
    SD sd(L,N), sdRightNew, sdLeftNew;
    LogHop logHop; logHop.logw=complex<double>(1.2,1.5); logHop.matrix=matrix;
    sd.logwRef()=1.6; sd.wfRef() = wfOld;

    LogHopSDOperation oneBodyWalkerOperation("fixedOrder", 8);

    oneBodyWalkerOperation.applyToRight(logHop, sd, sdRightNew);
    EXPECT_EQ( static_cast<size_t>(8), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfRightNew, sdRightNew.getWf(), 1e-12) );

    oneBodyWalkerOperation.applyToLeft(logHop, sd, sdLeftNew);
    EXPECT_EQ( static_cast<size_t>(8), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfLeftNew, sdLeftNew.getWf(), 1e-12) );
}

TEST_F(LogHopSDOperationTest, fixedOrderAcuracy)
{
    SD sd(L,N), sdRightNew, sdLeftNew;
    LogHop logHop; logHop.logw=complex<double>(1.2,1.5); logHop.matrix=matrix;
    sd.logwRef()=1.6; sd.wfRef() = wfOld;

    LogHopSDOperation oneBodyWalkerOperation("fixedOrder", 0, 1e-8);

    oneBodyWalkerOperation.applyToRight(logHop, sd, sdRightNew);
    EXPECT_EQ( static_cast<size_t>(6), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfRightNew, sdRightNew.getWf(), 1e-8) );

    oneBodyWalkerOperation.applyToLeft(logHop, sd, sdLeftNew);
    EXPECT_EQ( static_cast<size_t>(6), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfLeftNew, sdLeftNew.getWf(), 1e-8) );
}

TEST_F(LogHopSDOperationTest, dynamicOrder)
{
    SD sd(L,N), sdRightNew, sdLeftNew;
    LogHop logHop; logHop.logw=complex<double>(1.2,1.5); logHop.matrix=matrix;
    sd.logwRef()=1.6; sd.wfRef() = wfOld;

    LogHopSDOperation oneBodyWalkerOperation("dynamicOrder", 0, 1e-8);

    oneBodyWalkerOperation.applyToRight(logHop, sd, sdRightNew);
    EXPECT_EQ( static_cast<size_t>(6), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfRightNew, sdRightNew.getWf(), 1e-8) );

    oneBodyWalkerOperation.applyToLeft(logHop, sd, sdLeftNew);
    EXPECT_EQ( static_cast<size_t>(6), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfLeftNew, sdLeftNew.getWf(), 1e-8) );
}
