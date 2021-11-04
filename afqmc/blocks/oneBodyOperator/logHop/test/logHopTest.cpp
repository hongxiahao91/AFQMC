//
// Created by Hao Shi on 1/13/18.
//

#include "../include/logHop.h"
#include "../../../../../common/testHao/gtest_custom.h"

using namespace std;
using namespace tensor_hao;

TEST(LogHopTest, voidConstruction)
{
    size_t L(0);
    LogHop logHop;
    EXPECT_COMPLEXDOUBLE_EQ( 0.0, logHop.logw );
    EXPECT_FALSE( logHop.matrix.data() );
    EXPECT_EQ( L, logHop.getL() );
}

TEST(LogHopTest, size_tConstruction)
{
    size_t L(10);
    LogHop logHop(L);
    EXPECT_COMPLEXDOUBLE_EQ( 0.0, logHop.logw );
    EXPECT_TRUE( logHop.matrix.data() );
    EXPECT_EQ( L, logHop.getL() );
}

TEST(LogHopTest, copyConstruction)
{
    size_t L(10);
    TensorHao<complex<double>,2> matrix(L,L);
    randomFill(matrix);
    complex<double> logw(2.0, 3.0);

    LogHop logHopBase; logHopBase.logw=logw; logHopBase.matrix = matrix;
    LogHop logHop(logHopBase);

    EXPECT_COMPLEXDOUBLE_EQ( logw, logHop.logw );
    EXPECT_FALSE( diff( matrix, logHop.matrix, 1e-12 ) );
    EXPECT_EQ( L, logHop.getL() );
    EXPECT_FALSE( diff( matrix, logHopBase.matrix, 1e-12 ) );
}

TEST(LogHopTest, moveConstruction)
{
    size_t L(10);
    TensorHao<complex<double>,2> matrix(L,L);
    randomFill(matrix);
    complex<double> logw(2.0, 3.0);

    LogHop logHopBase; logHopBase.logw=logw; logHopBase.matrix = matrix;
    LogHop logHop( move( logHopBase ) );

    EXPECT_COMPLEXDOUBLE_EQ( logw, logHop.logw );
    EXPECT_FALSE( diff( matrix, logHop.matrix, 1e-12 ) );
    EXPECT_EQ( L, logHop.getL() );
    EXPECT_FALSE( logHopBase.matrix.data() );
}

TEST(LogHopTest, copyAssignment)
{
    size_t L(10);
    TensorHao<complex<double>,2> matrix(L,L);
    randomFill(matrix);
    complex<double> logw(2.0, 3.0);

    LogHop logHopBase; logHopBase.logw=logw; logHopBase.matrix = matrix;
    LogHop logHop; logHop = logHopBase;

    EXPECT_COMPLEXDOUBLE_EQ( logw, logHop.logw );
    EXPECT_FALSE( diff( matrix, logHop.matrix, 1e-12 ) );
    EXPECT_EQ( L, logHop.getL() );
    EXPECT_FALSE( diff( matrix, logHopBase.matrix, 1e-12 ) );
}

TEST(LogHopTest, moveAssignment)
{
    size_t L(10);
    TensorHao<complex<double>,2> matrix(L,L);
    randomFill(matrix);
    complex<double> logw(2.0, 3.0);

    LogHop logHopBase; logHopBase.logw=logw; logHopBase.matrix = matrix;
    LogHop logHop; logHop = move( logHopBase );

    EXPECT_COMPLEXDOUBLE_EQ( logw, logHop.logw );
    EXPECT_FALSE( diff( matrix, logHop.matrix, 1e-12 ) );
    EXPECT_EQ( L, logHop.getL() );
    EXPECT_FALSE( logHopBase.matrix.data() );
}