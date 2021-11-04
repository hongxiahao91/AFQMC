//
// Created by Hao Shi on 10/4/17.
//

#include "../include/hop2isSD2sOperation.h"
#include "../../../../../common/testHao/gtest_custom.h"

using namespace std;
using namespace tensor_hao;

class Hop2isSD2sOperationTest: public ::testing::Test
{
 public:
    size_t L, Nup, Ndn;
    TensorHao<complex<double>,2> matrix, wfUpOld, wfDnOld, wfRightUpNew, wfRightDnNew, wfLeftUpNew, wfLeftDnNew;

    Hop2isSD2sOperationTest( )
    {
        L=10; Nup=5; Ndn=3;

        matrix.resize(L,L);
        wfUpOld.resize(L,Nup); wfDnOld.resize(L,Ndn);
        wfRightUpNew.resize(L,Nup); wfRightDnNew.resize(L,Ndn);
        wfLeftUpNew.resize(L,Nup); wfLeftDnNew.resize(L,Ndn);

        gmm_cpu(matrix, wfUpOld, wfRightUpNew);
        gmm_cpu(matrix, wfDnOld, wfRightDnNew);
        gmm_cpu(matrix, wfUpOld, wfLeftUpNew, 'C');
        gmm_cpu(matrix, wfDnOld, wfLeftDnNew, 'C');
    }

    ~Hop2isSD2sOperationTest( )  {}
};

TEST_F(Hop2isSD2sOperationTest, applyToRight)
{
    SD2s sd2s(L,Nup, Ndn), sd2sNew;
    Hop2is hop2is; hop2is.logw=complex<double>(1.2,1.5); hop2is.matrix=matrix;
    sd2s.logwRef()=1.6; sd2s.wfUpRef() = wfUpOld; sd2s.wfDnRef() = wfDnOld;

    Hop2isSD2sOperation hop2isSD2sOperation;
    hop2isSD2sOperation.applyToRight(hop2is, sd2s, sd2sNew);

    EXPECT_FALSE( diff(wfRightUpNew, sd2sNew.getWfUp(), 1e-12) );
    EXPECT_FALSE( diff(wfRightDnNew, sd2sNew.getWfDn(), 1e-12) );
    EXPECT_COMPLEXDOUBLE_EQ( complex<double>(2.8, 1.5), sd2sNew.getLogw() );
}

TEST_F(Hop2isSD2sOperationTest, applyToLeft)
{
    SD2s sd2s(L,Nup, Ndn), sd2sNew;
    Hop2is hop2is; hop2is.logw=complex<double>(1.2,1.5); hop2is.matrix=matrix;
    sd2s.logwRef()=1.6; sd2s.wfUpRef() = wfUpOld; sd2s.wfDnRef() = wfDnOld;

    Hop2isSD2sOperation hop2isSD2sOperation;
    hop2isSD2sOperation.applyToLeft(hop2is, sd2s, sd2sNew);

    EXPECT_FALSE( diff(wfLeftUpNew, sd2sNew.getWfUp(), 1e-12) );
    EXPECT_FALSE( diff(wfLeftDnNew, sd2sNew.getWfDn(), 1e-12) );
    EXPECT_COMPLEXDOUBLE_EQ( complex<double>(2.8, -1.5), sd2sNew.getLogw() );
}