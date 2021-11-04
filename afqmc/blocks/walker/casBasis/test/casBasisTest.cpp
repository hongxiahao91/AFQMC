//
// Created by Hao Shi on 8/21/17.
//
#include "../include/casBasis.h"

using namespace std;
using namespace tensor_hao;

TEST(casBasisTest, checkOccupancyOrder)
{
    vector< vector<size_t> > occupancy;

    EXPECT_EQ( 1, checkOccupancyOrder(occupancy) );

    occupancy.resize(3);
    EXPECT_EQ( 1, checkOccupancyOrder(occupancy) );

    occupancy[0]={1,2,3,4}; occupancy[1]={2,3,4,5}; occupancy[2]={2,7,6,8};
    EXPECT_EQ(0, checkOccupancyOrder(occupancy) );

    occupancy[0]={1,2,3,4}; occupancy[1]={2,3,4,5}; occupancy[2]={2,6,8,10};
    EXPECT_EQ(1, checkOccupancyOrder(occupancy) );
}

TEST(casBasisTest, calculateNumberOfDifferentCasOrbitals)
{
    vector<size_t> a={1,2,3,4,5};
    vector<size_t> b(5);
    size_t expect, actual;

    b={1,2,3,4,5};
    expect = 0;
    actual = calculateNumberOfDifferentCasOrbitals(a,b);
    EXPECT_EQ( expect, actual);

    b={2,3,4,5,6};
    expect = 1;
    actual = calculateNumberOfDifferentCasOrbitals(a,b);
    EXPECT_EQ( expect, actual);

    b={2,5,7,8,9};
    expect = 3;
    actual = calculateNumberOfDifferentCasOrbitals(a,b);
    EXPECT_EQ( expect, actual);

    b={1,3,5,6,7};
    expect = 2;
    actual = calculateNumberOfDifferentCasOrbitals(a,b);
    EXPECT_EQ( expect, actual);

    b={5,6,7,8,9};
    expect = 4;
    actual = calculateNumberOfDifferentCasOrbitals(a,b);
    EXPECT_EQ( expect, actual);

    b={6,7,8,9,10};
    expect = 5;
    actual = calculateNumberOfDifferentCasOrbitals(a,b);
    EXPECT_EQ( expect, actual);
}

TEST(casBasisTest,findMinDiffNumberAndFrequency)
{
    vector<size_t> diff;
    size_t minDiffNumber, frequency;

    diff = {1,2,3,4,5};
    tie(minDiffNumber, frequency) = findMinDiffNumberAndFrequency(diff);
    EXPECT_EQ(static_cast<size_t>(1), minDiffNumber);
    EXPECT_EQ(static_cast<size_t>(1), frequency);

    diff = {5,2,3,2,2};
    tie(minDiffNumber, frequency) = findMinDiffNumberAndFrequency(diff);
    EXPECT_EQ(static_cast<size_t>(2), minDiffNumber);
    EXPECT_EQ(static_cast<size_t>(3), frequency);

    diff = {5,7,3,4,8,3,4,3,9,3,4,5,3};
    tie(minDiffNumber, frequency) = findMinDiffNumberAndFrequency(diff);
    EXPECT_EQ(static_cast<size_t>(3), minDiffNumber);
    EXPECT_EQ(static_cast<size_t>(5), frequency);
}

TEST(casBasisTest, buildParentTree)
{
    vector<int> isParent, isParentExpect;
    vector<size_t> tree, treeExpect;
    vector<vector<size_t> > occupancy;

    occupancy= { {1,2,3}, {1,2,4}, {1,2,5}, {1,2,6}, {1,2,7} };
    isParentExpect={1,0,0,0,0};
    treeExpect={0,1,2,3,4};
    buildParentTree(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(isParent, isParentExpect);
    EXPECT_VECOTR_EQ(tree, treeExpect);

    occupancy= { {1,2,3}, {1,2,4}, {1,2,5}, {1,3,6}, {1,3,7} };
    isParentExpect={1,0,0,0,0};
    treeExpect={0,1,2,3,4};
    buildParentTree(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(isParent, isParentExpect);
    EXPECT_VECOTR_EQ(tree, treeExpect);

    occupancy= { {1,2,3}, {1,2,4}, {1,2,5}, {1,4,6}, {1,4,7} };
    isParentExpect={1,1,0,0,0};
    treeExpect={0,2,1,3,4};
    buildParentTree(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(isParent, isParentExpect);
    EXPECT_VECOTR_EQ(tree, treeExpect);

    occupancy={ {1,2,3,4}, {1,3,4,5}, {1,5,8,9}, {1,2,3,6}, {1,2,4,8}, {1,3,8,9}, {1,2,5,6}, {1,5,7,9} };
    isParentExpect={1,0,1,1,0,0,1,0};
    treeExpect={0,1,4,3,6,7,2,5};
    buildParentTree(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(isParent, isParentExpect);
    EXPECT_VECOTR_EQ(tree, treeExpect);

    occupancy={ {1,2,3,4}, {2,3,4,5}, {2,5,7,9}, {1,2,5,6}, {2,4,7,8}, {1,2,6,9}, {3,5,6,7}, {2,5,8,9} };
    isParentExpect={1,1,0,1,0,0,0,0};
    treeExpect={0,1,2,7,4,6,3,5};
    buildParentTree(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(isParent, isParentExpect);
    EXPECT_VECOTR_EQ(tree, treeExpect);


    occupancy={ {1,2,3,4}, {1,3,5,6}, {1,2,3,5}, {1,3,4,8}, {1,3,4,5}, {1,7,8,9}, {1,2,4,6}, {2,9,10,11} };
    isParentExpect={1,1,1,0,0,1,0,0};
    treeExpect={0,6,3,4,2,1,5,7};
    buildParentTree(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(isParent, isParentExpect);
    EXPECT_VECOTR_EQ(tree, treeExpect);
}

TEST(casBasisTest, permuteOrbitalsInOccupancy)
{
    vector<int> isParent;
    vector<size_t> tree;
    vector<vector<size_t> > occupancy, occupancyExpect;
    vector<int> sign, signExpect;

    occupancy= {{1, 2, 3}, {1, 2, 4}, {1, 2, 5}, {1, 2, 6}, {1, 2, 7} };
    occupancyExpect = {{1, 2, 3}, {1, 2, 4}, {1, 2, 5}, {1, 2, 6}, {1, 2, 7} };
    signExpect={1,1,1,1,1};
    buildParentTree(isParent, tree, occupancy);
    sign = permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(occupancy, occupancyExpect);
    EXPECT_VECOTR_EQ(sign, signExpect);

    occupancy= { {1,2,3}, {1,2,4}, {1,2,5}, {1,3,6}, {1,3,7} };
    occupancyExpect = { {1,2,3}, {1,2,4}, {1,2,5}, {1,6,3}, {1,7,3} };
    signExpect={1,1,1,-1, -1};
    buildParentTree(isParent, tree, occupancy);
    sign = permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(occupancy, occupancyExpect);
    EXPECT_VECOTR_EQ(sign, signExpect);

    occupancy=  { {1,2,3}, {1,2,4}, {1,2,5}, {1,4,6}, {1,4,7} };
    occupancyExpect = { {1,2,3}, {1,2,4}, {1,2,5}, {1,6,4}, {1,7,4} };
    signExpect={1,1,1,-1, -1};
    buildParentTree(isParent, tree, occupancy);
    sign = permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(occupancy, occupancyExpect);
    EXPECT_VECOTR_EQ(sign, signExpect);

    occupancy=  { {1,2,3,4}, {1,3,4,5}, {1,5,8,9}, {1,2,3,6}, {1,2,4,8}, {1,3,8,9}, {1,2,5,6}, {1,5,7,9} };
    occupancyExpect = { {1,2,3,4}, {1,5,3,4}, {1,8,5,9}, {1,2,3,6}, {1,2,8,4}, {1,8,3,9}, {1,2,5,6}, {1,7,5,9} };
    signExpect={1,1,-1,1,-1,-1,1,-1};
    buildParentTree(isParent, tree, occupancy);
    sign = permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(occupancy, occupancyExpect);
    EXPECT_VECOTR_EQ(sign, signExpect);

    occupancy= { {1,2,3,4}, {2,3,4,5}, {2,5,7,9}, {1,2,5,6}, {2,4,7,8}, {1,2,6,9}, {3,5,6,7}, {2,5,8,9} };
    occupancyExpect = { {1,2,3,4}, {5,2,3,4}, {5,2,7,9}, {5,2,1,6}, {8,2,7,4}, {9,2,1,6}, {5,6,3,7}, {5,2,8,9} };
    signExpect={1,-1,-1,-1,1,1,1,-1};
    buildParentTree(isParent, tree, occupancy);
    sign = permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(occupancy, occupancyExpect);
    EXPECT_VECOTR_EQ(sign, signExpect);

    occupancy= { {1,2,3,4}, {1,3,5,6}, {1,2,3,5}, {1,3,4,8}, {1,3,4,5}, {1,7,8,9}, {1,2,4,6}, {2,9,10,11} };
    occupancyExpect = { {1,2,3,4}, {1,6,3,5}, {1,2,3,5}, {1,8,3,4}, {1,5,3,4}, {1,7,8,9}, {1,2,6,4}, {2,11,10,9} };
    signExpect={1,1,1,1,1,1,-1,-1};
    buildParentTree(isParent, tree, occupancy);
    sign = permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(occupancy, occupancyExpect);
    EXPECT_VECOTR_EQ(sign, signExpect);
}

TEST(casBasisTest, calculateDiffOrbitalsInOccupancy)
{
    vector<int> isParent;
    vector<size_t> tree;
    vector<vector<size_t> > occupancy;
    vector<vector<size_t> > diffOrbital, diffOrbitalExpect;

    occupancy= {{1, 2, 3}, {1, 2, 4}, {1, 2, 5}, {1, 2, 6}, {1, 2, 7} };
    diffOrbitalExpect={{},{2},{2},{2},{2}};
    buildParentTree(isParent, tree, occupancy);
    permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    diffOrbital = calculateDiffOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(diffOrbital, diffOrbitalExpect);

    occupancy= { {1,2,3}, {1,2,4}, {1,2,5}, {1,3,6}, {1,3,7} };
    diffOrbitalExpect={{},{2},{2},{1},{1}};
    buildParentTree(isParent, tree, occupancy);
    permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    diffOrbital = calculateDiffOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(diffOrbital, diffOrbitalExpect);

    occupancy= { {1,2,3}, {1,2,4}, {1,2,5}, {1,4,6}, {1,4,7} };
    diffOrbitalExpect={{},{2},{2},{1},{1}};
    buildParentTree(isParent, tree, occupancy);
    permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    diffOrbital = calculateDiffOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(diffOrbital, diffOrbitalExpect);

    occupancy= { {1,2,3,4}, {1,3,4,5}, {1,5,8,9}, {1,2,3,6}, {1,2,4,8}, {1,3,8,9}, {1,2,5,6}, {1,5,7,9} };
    diffOrbitalExpect={{},{1},{1,3},{3},{2},{2},{2},{1,3}};
    buildParentTree(isParent, tree, occupancy);
    permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    diffOrbital = calculateDiffOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(diffOrbital, diffOrbitalExpect);

    occupancy= { {1,2,3,4}, {2,3,4,5}, {2,5,7,9}, {1,2,5,6}, {2,4,7,8}, {1,2,6,9}, {3,5,6,7}, {2,5,8,9} };
    diffOrbitalExpect={{},{0},{2,3},{2,3},{0,2},{0},{1,3},{2,3}};
    buildParentTree(isParent, tree, occupancy);
    permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    diffOrbital = calculateDiffOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(diffOrbital, diffOrbitalExpect);

    occupancy= { {1,2,3,4}, {1,3,5,6}, {1,2,3,5}, {1,3,4,8}, {1,3,4,5}, {1,7,8,9}, {1,2,4,6}, {2,9,10,11} };
    diffOrbitalExpect={{},{1},{3},{1},{1},{1,2,3},{2},{0,1,2}};
    buildParentTree(isParent, tree, occupancy);
    permuteOrbitalsInOccupancy(isParent, tree, occupancy);
    diffOrbital = calculateDiffOrbitalsInOccupancy(isParent, tree, occupancy);
    EXPECT_VECOTR_EQ(diffOrbital, diffOrbitalExpect);
}

TEST(casBasisTest, sortOrbitalsInOccupancy)
{
    vector<vector<size_t> > occupancy;
    vector<vector<size_t> > occupancySort, occupancySortExpect;
    vector<int> sign,signExpect;

    occupancy     = { {1,4,3,2}, {6,5,3,1}, {2,1,5,3}, {4,8,1,3}, {4,5,1,3}, {1,8,7,9}, {6,1,2,4}, {11,10,9,2} };
    occupancySortExpect = { {1,2,3,4}, {1,3,5,6}, {1,2,3,5}, {1,3,4,8}, {1,3,4,5}, {1,7,8,9}, {1,2,4,6}, {2,9,10,11} };
    signExpect = {-1, 1, 1, 1, 1, -1, -1, 1};

    sign = sortOrbitalsInOccupancy(occupancySort, occupancy);

    EXPECT_EQ(occupancySort, occupancySortExpect);
    EXPECT_EQ(sign, signExpect);
}

TEST(casBasisTest, fillCasMatrixByRow)
{
    TensorHao<complex<double>, 2> sm, bm(5,2), smExpect(3,2);
    vector<size_t> row={1,3,4};

    bm = { {1,10},{2,20},{3,30},{4,40},{5,50},{6,60},{7,70},{8,80},{9,90},{10,100} };
    smExpect = {{2,20},{4,40},{5,50},{7,70},{9,90},{10,100}};
    fillCasMatrixByRow(sm, bm, row);

    EXPECT_FALSE(diff(sm, smExpect, 1e-12));
}

TEST(casBasisTest, fillCasMatrixByCol)
{
    TensorHao<complex<double>, 2> sm, bm(2,5), smExpect(2,3);
    vector<size_t> col={1,3,4};

    bm = { {1,10},{2,20},{3,30},{4,40},{5,50},{6,60},{7,70},{8,80},{9,90},{10,100} };
    smExpect = {{3,30},{4,40},{7,70},{8,80},{9,90},{10,100}};
    fillCasMatrixByCol(sm, bm, col);

    EXPECT_FALSE(diff(sm, smExpect, 1e-12));
}

TEST(casBasisTest, getRowDifference)
{
    TensorHao<complex<double>, 2> V, bm(6,3), vExpect(2,3);
    vector<size_t> row1={1,3,4};
    vector<size_t> row2={2,3,5};
    vector<size_t> difference={0,2};

    bm = { {1,10},{3,20},{5,30},{4,40},{5,50},{8,60},{7,70},{10,80},{9,90},{10,100},
           {11,111},{12,120},{13,130},{14,140},{16,150},{16,160},{19,170},{18,180} };

    vExpect = { {-2,-10},{-3,-10},{1,-10},{-1,-9},{-2,-10},{1,-10} };
    getRowDifference(V, bm, row1, row2, difference);

    EXPECT_FALSE( diff(V, vExpect, 1e-12) );
}
