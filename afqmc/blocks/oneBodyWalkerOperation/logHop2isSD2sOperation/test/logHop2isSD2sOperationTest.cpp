//
// Created by Hao Shi on 10/4/17.
//

#include "../include/logHop2isSD2sOperation.h"
#include "../../../../../common/testHao/gtest_custom.h"

using namespace std;
using namespace tensor_hao;

class LogHop2isSD2sOperationTest: public ::testing::Test
{
 public:
    size_t L, Nup, Ndn;
    TensorHao<complex<double>,2> matrix, wfUpOld, wfUpRightNew, wfUpLeftNew, wfDnOld, wfDnRightNew, wfDnLeftNew;


    LogHop2isSD2sOperationTest( )
    {
        L=3; Nup=2; Ndn=1;
        matrix.resize(L,L);
        wfUpOld.resize(L, Nup); wfUpRightNew.resize(L, Nup); wfUpLeftNew.resize(L, Nup);
        wfDnOld.resize(L, Ndn); wfDnRightNew.resize(L, Ndn); wfDnLeftNew.resize(L, Ndn);

        matrix = {
                {2,  3}, {4,   2}, {1,   1},
                {1, -1}, {0, 0.3}, {1,  -1},
                {3,  2}, {2,   3}, {0, 0.5}
        };
        matrix = complex<double>(-0.01, 0.0) * matrix;

        wfUpOld = { {1, 0}, {2,0}, {3, 0}, {3,0}, {4,0}, {5, 0} };

        wfUpRightNew = { {0.8717255743431891,-0.06695756061260984}, {1.9019539751922434,-0.11272541520548454}, {2.9713563453236675, -0.0038676809054978764},
                       {2.753351116138562, -0.14393925193667553}, {3.7836486105153155,-0.21511811236274694}, {4.932649990968739, -0.012861814634312912} };
        wfUpLeftNew  = { {0.8708226627963166, 0.09656591603485742}, {1.961995737773784, -0.034019035572233965},{2.9316294741182647, 0.09182786294989435},
                       {2.7315235332357157, 0.21256279549696155}, {3.924363675467334, -0.0680085413131974},  {4.833174488748831, 0.19818901143889675} };

        wfDnOld = { {1.5, 0}, {2.8,0}, {3.9, 0} };
        wfDnRightNew = { {1.3273654189481012,-0.09089527826133821}, {2.66465278547187, -0.15093497686105645}, {3.8588353397101534, -0.004986335990389258} };
        wfDnLeftNew  = { {1.3200948148397493, 0.1351690256672221}, {2.7488038959252905, -0.04561569450693961}, {3.8012155440416793, 0.12907127412412336} };


    }

    ~LogHop2isSD2sOperationTest( )  {}
};

TEST_F(LogHop2isSD2sOperationTest, fixedOrderNumber)
{
    SD2s sd2s(L,Nup, Ndn), sd2sRightNew, sd2sLeftNew;
    LogHop2is logHop2is; logHop2is.logw=complex<double>(1.2,1.5); logHop2is.matrix=matrix;
    sd2s.logwRef()=1.6; sd2s.wfUpRef() = wfUpOld; sd2s.wfDnRef() = wfDnOld;

    LogHop2isSD2sOperation oneBodyWalkerOperation("fixedOrder", 8, 1e-8);

    oneBodyWalkerOperation.applyToRight(logHop2is, sd2s, sd2sRightNew);
    EXPECT_EQ( static_cast<size_t>(8), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfUpRightNew, sd2sRightNew.getWfUp(), 1e-12) );
    EXPECT_FALSE( diff(wfDnRightNew, sd2sRightNew.getWfDn(), 1e-12) );

    oneBodyWalkerOperation.applyToLeft(logHop2is, sd2s, sd2sLeftNew);
    EXPECT_EQ( static_cast<size_t>(8), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfUpLeftNew, sd2sLeftNew.getWfUp(), 1e-12) );
    EXPECT_FALSE( diff(wfDnLeftNew, sd2sLeftNew.getWfDn(), 1e-12) );
}

TEST_F(LogHop2isSD2sOperationTest, fixedOrderAcuracy)
{
    SD2s sd2s(L,Nup, Ndn), sd2sRightNew, sd2sLeftNew;
    LogHop2is logHop2is; logHop2is.logw=complex<double>(1.2,1.5); logHop2is.matrix=matrix;
    sd2s.logwRef()=1.6; sd2s.wfUpRef() = wfUpOld; sd2s.wfDnRef() = wfDnOld;

    LogHop2isSD2sOperation oneBodyWalkerOperation("fixedOrder", 0, 1e-8);

    oneBodyWalkerOperation.applyToRight(logHop2is, sd2s, sd2sRightNew);
    EXPECT_EQ( static_cast<size_t>(6), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfUpRightNew, sd2sRightNew.getWfUp(), 1e-8) );
    EXPECT_FALSE( diff(wfDnRightNew, sd2sRightNew.getWfDn(), 1e-8) );

    oneBodyWalkerOperation.applyToLeft(logHop2is, sd2s, sd2sLeftNew);
    EXPECT_EQ( static_cast<size_t>(6), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfUpLeftNew, sd2sLeftNew.getWfUp(), 1e-8) );
    EXPECT_FALSE( diff(wfDnLeftNew, sd2sLeftNew.getWfDn(), 1e-8) );
}

TEST_F(LogHop2isSD2sOperationTest, dynamicOrder)
{
    SD2s sd2s(L,Nup, Ndn), sd2sRightNew, sd2sLeftNew;
    LogHop2is logHop2is; logHop2is.logw=complex<double>(1.2,1.5); logHop2is.matrix=matrix;
    sd2s.logwRef()=1.6; sd2s.wfUpRef() = wfUpOld; sd2s.wfDnRef() = wfDnOld;

    LogHop2isSD2sOperation oneBodyWalkerOperation("dynamicOrder", 0, 1e-8);

    oneBodyWalkerOperation.applyToRight(logHop2is, sd2s, sd2sRightNew);
    EXPECT_EQ( static_cast<size_t>(6), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfUpRightNew, sd2sRightNew.getWfUp(), 1e-8) );
    EXPECT_FALSE( diff(wfDnRightNew, sd2sRightNew.getWfDn(), 1e-8) );

    oneBodyWalkerOperation.applyToLeft(logHop2is, sd2s, sd2sLeftNew);
    EXPECT_EQ( static_cast<size_t>(6), oneBodyWalkerOperation.getCurrentOrder() );
    EXPECT_FALSE( diff(wfUpLeftNew, sd2sLeftNew.getWfUp(), 1e-8) );
    EXPECT_FALSE( diff(wfDnLeftNew, sd2sLeftNew.getWfDn(), 1e-8) );
}
