//
// Created by Hao Shi on 1/13/18.
//

#include <climits>
#include "../include/logHopSDOperation.h"
#include "../../../../../common/common.h"

using namespace std;
using namespace tensor_hao;

LogHopSDOperation::LogHopSDOperation(std::string flag, size_t taylorOrder, double accuracy, size_t baseTaylorOrder)
{
    reset(flag, taylorOrder, accuracy, baseTaylorOrder);
}

LogHopSDOperation::~LogHopSDOperation() { }

void LogHopSDOperation::reset(string flag, size_t taylorOrder, double accuracy, size_t baseTaylorOrder)
{
    LogHopSDOperation::flag = flag;
    LogHopSDOperation::taylorOrder = taylorOrder;
    LogHopSDOperation::accuracy = accuracy;
    LogHopSDOperation::baseTaylorOrder = baseTaylorOrder;

    operationNumber  = 0;
    minTaylorOrder   = INT_MAX;
    maxTaylorOrder   = 0;
    totalTaylorOrder = 0;

    if( flag != "fixedOrder" && flag != "dynamicOrder" )
    {
        cout<<"Error!!! Do not know the input flag for LogHopSDOperation: "<<flag<<endl;
        exit(1);
    }
}

void LogHopSDOperation::applyToRight(const LogHop &oneBody, const SD &walker, SD &walkerNew)
{
    checkAndResize(oneBody, walker, walkerNew);

    char TRANSOneBody='N';
    addOrders(oneBody, walker, walkerNew, TRANSOneBody);
    walkerNew.logwRef() = oneBody.logw + walker.getLogw();
}

void LogHopSDOperation::applyToLeft(const LogHop &oneBody, const SD &walker, SD &walkerNew)
{
    checkAndResize(oneBody, walker, walkerNew);

    char TRANSOneBody='C';
    addOrders(oneBody, walker, walkerNew, TRANSOneBody);
    walkerNew.logwRef() = conj( oneBody.logw ) + walker.getLogw();
}

void LogHopSDOperation::print()
{
    if( flag == "fixedOrder" )
    {
        if( MPIRank() == 0 )
        {
            if( taylorOrder == 0 ) cout<<"Use fixed order, the order is not determined yet."<<endl;
            else cout<<"Use fixedOrder, the order is "<<taylorOrder<<endl;
        }
    }
    else if( flag == "dynamicOrder" )
    {
        size_t globalNumber, globalMin, globalMax, globalTotal;
        globalNumber = MPISum( operationNumber );
#ifdef MPI_HAO
        MPIReduce(minTaylorOrder, globalMin, MPI_MIN);
        MPIReduce(maxTaylorOrder, globalMax, MPI_MAX);
#else
        globalMin = minTaylorOrder;
        globalMax = maxTaylorOrder;
#endif
        globalTotal  = MPISum( totalTaylorOrder );

        if( MPIRank() == 0 )
        {
            cout<<"Use dynamicOrder. "<<endl;
            cout<<"The base order is "<<baseTaylorOrder<<endl;
            cout<<"The min order is "<<globalMin<<endl;
            cout<<"The max order is "<<globalMax<<endl;
            cout<<"The average order is "<< globalTotal*1.0/globalNumber <<endl;
        }
    }
}

size_t LogHopSDOperation::getCurrentOrder() const { return currentOrder; }

void LogHopSDOperation::checkAndResize(const LogHop &oneBody, const SD &walker, SD &walkerNew) const
{
    size_t L = walker.getL(); size_t N = walker.getN();
    if( oneBody.getL() !=  L ) { cout<<"Error!!! LogHop size is not consistent with walker!"<<endl; exit(1); }
    if( walkerNew.getL() != L  ||  walkerNew.getN() != N ) walkerNew.resize( L, N );
}

void LogHopSDOperation::addOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody)
{
    if( flag == "fixedOrder" )
    {
        if( taylorOrder==0 ) determinantAndAddFixedOrders(oneBody, walker, walkerNew, TRANSOneBody);
        else initialAndAddFixedOrders(oneBody, walker, walkerNew, taylorOrder, TRANSOneBody);
        operationNumber++;
    }
    else if( flag == "dynamicOrder" )
    {
        initialAndAddDynamicOrders(oneBody, walker, walkerNew, TRANSOneBody);
        minTaylorOrder = min( minTaylorOrder, currentOrder );
        maxTaylorOrder = max( maxTaylorOrder, currentOrder );
        totalTaylorOrder += currentOrder;
        operationNumber++;
    }
    else
    {
        cout<<"Error!!! Do not know the input flag for LogHopSDOperation: "<<flag<<endl;
        exit(1);
    }
    clearWfTemp();
}

void LogHopSDOperation::determinantAndAddFixedOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody)
{
    initialAndAddFixedOrders(oneBody, walker, walkerNew, baseTaylorOrder, TRANSOneBody);
    addDynamicOrders(oneBody, walker, walkerNew, TRANSOneBody);

#ifdef MPI_HAO
    MPIReduce(currentOrder, taylorOrder, MPI_MAX);
#else
    taylorOrder = currentOrder;
#endif
}

void LogHopSDOperation::initialAndAddFixedOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, size_t maxOrder, char TRANSOneBody)
{
    initialZeroOrder(oneBody, walker, walkerNew, TRANSOneBody);
    addFixedOrders(oneBody, walker, walkerNew, maxOrder, TRANSOneBody);
}

void LogHopSDOperation::initialAndAddDynamicOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody)
{
    initialAndAddFixedOrders(oneBody, walker, walkerNew, baseTaylorOrder, TRANSOneBody);
    addDynamicOrders(oneBody, walker, walkerNew, TRANSOneBody);
}

void LogHopSDOperation::initialZeroOrder(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody)
{
    size_t L = walker.getL(); size_t N = walker.getN();
    if( wfTempNew.rank(0) != L || wfTempNew.rank(1) != N ) wfTempNew.resize( L, N );
    if( wfTempOld.rank(0) != L || wfTempOld.rank(1) != N ) wfTempOld.resize( L, N );

    walkerNew.wfRef() = walker.getWf();
    wfTempOld = walker.getWf();

    currentOrder = 0;
}

void LogHopSDOperation::addFixedOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, size_t maxOrder, char TRANSOneBody)
{
    while( currentOrder < maxOrder )
    {
        BL_NAME(gmm)( oneBody.matrix, wfTempOld, wfTempNew, TRANSOneBody, 'N', 1.0/(currentOrder+1.0) );
        walkerNew.wfRef() += wfTempNew;
        wfTempOld = move( wfTempNew );

        currentOrder++;
    }
}

void LogHopSDOperation::addDynamicOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody)
{
    do
    {
        BL_NAME(gmm)( oneBody.matrix, wfTempOld, wfTempNew, TRANSOneBody, 'N', 1.0/(currentOrder+1.0) );
        walkerNew.wfRef() += wfTempNew;
        wfTempOld = move( wfTempNew );

        currentOrder++;
    } while( ! isConverged() );
}

bool LogHopSDOperation::isConverged()
{
    size_t L = wfTempOld.rank(0); size_t N = wfTempOld.rank(1);
    for(size_t i = 0; i < N; ++i)
    {
        for(size_t j = 0; j < L; ++j)
        {
            if( abs( wfTempOld(j,i) ) > accuracy  ) return false;
        }
    }

    return true;
}

void LogHopSDOperation::clearWfTemp()
{
    wfTempOld.resize(0,0);
    wfTempNew.resize(0,0);
}

LogHopSDOperation::LogHopSDOperation(const LogHopSDOperation &x) { }

LogHopSDOperation &LogHopSDOperation::operator=(const LogHopSDOperation &x) { return *this; }