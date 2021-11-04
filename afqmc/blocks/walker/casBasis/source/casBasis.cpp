//
// Created by Hao Shi on 8/21/17.
//

#include "../include/casBasis.h"

using namespace std;
using namespace tensor_hao;

int checkOccupancyOrder(const vector<vector<size_t> > &occupancy)
{
    size_t MDSize = occupancy.size(); if(MDSize == 0) return 1;
    size_t N = ( occupancy[0] ).size(); if( N == 0 ) return 1;

    for(size_t i = 0; i < MDSize; ++i)
    {
        for(size_t j = 0; j < (N-1); ++j)
        {
            if( occupancy[i][j] >= occupancy[i][j+1] ) return 0;
        }
    }

    return 1;
}

size_t calculateNumberOfDifferentCasOrbitals(const vector<size_t> &o1, const vector<size_t> &o2)
{
    size_t N=o1.size();
    size_t numberOfDiff(0), i(0), j(0);
    while( i!=N && j!=N )
    {
        if( o1[i]<o2[j] ) {numberOfDiff++; i++;}
        else if( o1[i]>o2[j] ) j++;
        else {i++; j++;}
    }

    return numberOfDiff+N-i;
}

tuple<size_t, size_t> findMinDiffNumberAndFrequency(const vector<size_t> &diff)
{
    size_t vecSize=diff.size();

    size_t minDiffNumber(INT_MAX);
    for(size_t i = 0; i < vecSize; ++i)
    {
        if( minDiffNumber>diff[i] ) minDiffNumber=diff[i];
    }

    size_t frequency(0);
    for(size_t i = 0; i < vecSize; ++i)
    {
        if( diff[i]==minDiffNumber ) frequency++;
    }

    return make_tuple(minDiffNumber, frequency);
}

void determineAvailChild(vector<size_t> &availChild, const vector<int> &isAvail)
{
    size_t MDSize = isAvail.size();

    availChild.resize(0);
    for(size_t i = 0; i < MDSize; ++i)
    {
        if( isAvail[i] == 1 ) availChild.push_back(i);
    }
}

int permuteOrbitals(const vector<size_t> &parent, vector<size_t> &child)
{
    size_t N=parent.size();

    int sign(1);
    for(size_t i = 0; i < N; ++i)
    {
        for(size_t j = 0; j < N; ++j)
        {
            if( parent[i] == child[j] )
            {
                if( i != j )
                {
                    swap( child[i], child[j] );
                    sign *= -1;
                }
            }
        }
    }
    return sign;
}

vector<size_t> calculateDiffOrbitals(const vector<size_t> &parent, const vector<size_t> &child)
{
    size_t N=parent.size();

    vector<size_t> diffOrbitals;
    for(size_t i = 0; i < N; ++i)
    {
        if( parent[i] != child[i] ) { diffOrbitals.push_back(i);}
    }

    return diffOrbitals;
}

void buildParentTree(vector<int> &isParent, vector<size_t> &tree, const vector<vector<size_t> > &occupancy)
{
    size_t MDSize = occupancy.size();

    isParent.resize(MDSize); for(size_t i = 0; i < MDSize; ++i) isParent[i] = 0;
    tree.resize(MDSize);

    vector<int> isAvail(MDSize, 1);
    size_t availChildSize;
    vector<size_t> availChild; availChild.reserve(MDSize);

    size_t parentBegin(0);
    size_t parentEnd(1);
    tree[0]    = 0;
    isAvail[0] = 0;

    vector<size_t> diff, diffTmp;
    diff.reserve(MDSize); diffTmp.reserve(MDSize);

    while( parentEnd != MDSize )
    {
        determineAvailChild(availChild, isAvail);
        availChildSize = availChild.size();

        size_t parentChosen(INT_MAX);
        size_t minDiffNumber(INT_MAX), frequency(0);
        size_t minDiffNumberTmp, frequencyTmp;
        diffTmp.resize(availChildSize);
        for(size_t i = parentBegin; i < parentEnd; ++i)
        {
            for(size_t j = 0; j < availChildSize; ++j)
            {
                diffTmp[j] = calculateNumberOfDifferentCasOrbitals(occupancy[tree[i]], occupancy[availChild[j]]);
            }

            tie(minDiffNumberTmp, frequencyTmp) = findMinDiffNumberAndFrequency(diffTmp);

            if(minDiffNumber > minDiffNumberTmp)
            {
                parentChosen = i;
                minDiffNumber = minDiffNumberTmp;
                frequency = frequencyTmp;
                diff = diffTmp;
            }
            else if(minDiffNumber == minDiffNumberTmp)
            {
                if(frequency < frequencyTmp)
                {
                    parentChosen = i;
                    frequency = frequencyTmp;
                    diff = diffTmp;
                }
            }
        }

        if( parentChosen != (parentEnd-1) ) swap( tree[parentChosen], tree[parentEnd- 1] );
        isParent[ tree[parentEnd-1] ] = 1;

        parentBegin = parentEnd;
        for(size_t j = 0; j < availChildSize; ++j)
        {
            if(diff[j] == minDiffNumber)
            {
                tree[parentEnd] = availChild[j];
                isAvail[ availChild[j] ] = 0;
                parentEnd++;
            }
        }
    }
}

vector<int> permuteOrbitalsInOccupancy(const vector<int> &isParent, const vector<size_t> &tree, vector<vector<size_t> > &occupancy)
{
    size_t MDSize = occupancy.size();

    vector<int> sign(MDSize, 1);

    size_t currentIndex, currentParent(0);
    for(size_t i = 1; i < MDSize; ++i)
    {
        currentIndex = tree[i];
        sign[ currentIndex ] = permuteOrbitals( occupancy[currentParent], occupancy[ currentIndex ] );
        if( isParent[ currentIndex ] ) currentParent=currentIndex;
    }

    return sign;
}

vector<vector<size_t>> calculateDiffOrbitalsInOccupancy(const vector<int> &isParent,
                                                        const vector<size_t> &tree,
                                                        const vector<vector<size_t> > &occupancy)
{
    size_t MDSize = occupancy.size();

    vector<vector<size_t>> diffOrbitals(MDSize);

    size_t currentIndex, currentParent(0);
    for(size_t i = 1; i < MDSize; ++i)
    {
        currentIndex = tree[i];
        diffOrbitals[ currentIndex ] = calculateDiffOrbitals( occupancy[currentParent], occupancy[ currentIndex ] );
        if( isParent[ currentIndex ] ) currentParent=currentIndex;
    }

    return diffOrbitals;
}

vector<int> sortOrbitalsInOccupancy(vector<vector<size_t>> &occupancySort, const vector<vector<size_t>> &occupancy)
{
    size_t MDSize = occupancy.size();
    if( MDSize==0 ) { occupancySort.resize(0); return vector<int>(); }

    size_t N = occupancy[0].size();

    occupancySort = occupancy;
    vector<int> signVector(MDSize, 1);

    if( N ==0 ) return signVector;

    int sign;
    for(size_t i = 0; i < MDSize; ++i)
    {
        sign=1;
        for(size_t j = 0; j < N; ++j)
        {
            for(size_t k = 0; k < (N-1); ++k)
            {
                if( occupancySort[i][k] > occupancySort[i][k+1] )
                {
                    swap(occupancySort[i][k], occupancySort[i][k+1]);
                    sign *= -1;
                }
            }
        }
        signVector[i] = sign;
    }
    return signVector;
}

void fillCasMatrixByRow(TensorHao<complex<double>, 2> &smallMatrix,
                        const TensorHao<complex<double>, 2> &fullMatrix,
                        const vector<size_t> &row)
{
    size_t L = row.size(); size_t N = fullMatrix.rank(1);
    size_t r;

    smallMatrix.resize(L, N);
    for(size_t i = 0; i < L; ++i)
    {
        r = row[i];
        for(size_t j = 0; j < N; ++j)
        {
            smallMatrix(i,j) = fullMatrix(r,j);
        }
    }
}

void fillCasMatrixByCol(TensorHao<complex<double>, 2> &smallMatrix,
                        const TensorHao<complex<double>, 2> &fullMatrix,
                        const vector<size_t> &col)
{
    size_t L = fullMatrix.rank(0); size_t N = col.size();
    size_t r;

    smallMatrix.resize(L, N);
    for(size_t i = 0; i < N; ++i)
    {
        r = col[i];
        for(size_t j = 0; j < L; ++j)
        {
            smallMatrix(j,i) = fullMatrix(j,r);
        }
    }
}

void getRowDifference(TensorHao<complex<double>, 2> &V,
                      const TensorHao<complex<double>, 2> &fullMatrix,
                      const vector<size_t> &row1,
                      const vector<size_t> &row2,
                      const vector<size_t> &diff)
{
    size_t L = diff.size(); size_t N = fullMatrix.rank(1);
    size_t r, r1, r2;

    V.resize(L, N);
    for(size_t i = 0; i < L; ++i)
    {
        r = diff[i]; r1 = row1[r]; r2 = row2[r];
        for(size_t j = 0; j < N; ++j)
        {
            V(i,j) = fullMatrix(r1,j) - fullMatrix(r2,j);
        }
    }
}
