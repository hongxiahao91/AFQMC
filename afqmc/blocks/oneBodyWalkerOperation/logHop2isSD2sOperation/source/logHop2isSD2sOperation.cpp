//
// Created by Hao Shi on 10/4/17.
//

#include <climits>
#include "../include/logHop2isSD2sOperation.h"
#include "../../../../../common/common.h"

using namespace std;
using namespace tensor_hao;

LogHop2isSD2sOperation::LogHop2isSD2sOperation(std::string flag, size_t taylorOrder, double accuracy, size_t baseTaylorOrder)
{
    reset(flag, taylorOrder, accuracy, baseTaylorOrder);
}

LogHop2isSD2sOperation::~LogHop2isSD2sOperation() { }

void LogHop2isSD2sOperation::reset(string flag, size_t taylorOrder, double accuracy, size_t baseTaylorOrder)
{
    LogHop2isSD2sOperation::flag = flag;
    LogHop2isSD2sOperation::taylorOrder = taylorOrder;
    LogHop2isSD2sOperation::accuracy = accuracy;
    LogHop2isSD2sOperation::baseTaylorOrder = baseTaylorOrder;

    operationNumber  = 0;
    minTaylorOrder   = INT_MAX;
    maxTaylorOrder   = 0;
    totalTaylorOrder = 0;

    if( flag != "fixedOrder" && flag != "dynamicOrder" )
    {
        cout<<"Error!!! Do not know the input flag for LogHop2isSD2sOperation: "<<flag<<endl;
        exit(1);
    }
}

void LogHop2isSD2sOperation::applyToRight(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew)
{
    checkAndResize(oneBody, walker, walkerNew);

    char TRANSOneBody='N';
    addOrders(oneBody, walker, walkerNew, TRANSOneBody);
    walkerNew.logwRef() = oneBody.logw + walker.getLogw();
}

void LogHop2isSD2sOperation::applyToLeft(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew)
{
    checkAndResize(oneBody, walker, walkerNew);

    char TRANSOneBody='C';
    addOrders(oneBody, walker, walkerNew, TRANSOneBody);
    walkerNew.logwRef() = conj( oneBody.logw ) + walker.getLogw();
}

void LogHop2isSD2sOperation::print()
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

size_t LogHop2isSD2sOperation::getCurrentOrder() const { return currentOrder; }

void LogHop2isSD2sOperation::checkAndResize(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew) const
{
    size_t L = walker.getL(); size_t Nup = walker.getNup(); size_t Ndn = walker.getNdn();
    if( oneBody.getL() !=  L ) {cout<<"Error!!! LogHop2is size is not consistent with walker!"<<endl; exit(1); }
    if( walkerNew.getL() != L  ||  walkerNew.getNup() != Nup || walkerNew.getNdn() != Ndn ) walkerNew.resize( L, Nup, Ndn );
}

void LogHop2isSD2sOperation::addOrders(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew, char TRANSOneBody)
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
        cout<<"Error!!! Do not know the input flag for LogHop2isSD2sOperation: "<<flag<<endl;
        exit(1);
    }
    clearWfTemp();
}

void LogHop2isSD2sOperation::determinantAndAddFixedOrders(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew, char TRANSOneBody)
{
    initialAndAddFixedOrders(oneBody, walker, walkerNew, baseTaylorOrder, TRANSOneBody);
    addDynamicOrders(oneBody, walker, walkerNew, TRANSOneBody);

#ifdef MPI_HAO
    MPIReduce(currentOrder, taylorOrder, MPI_MAX);
#else
    taylorOrder = currentOrder;
#endif
}

void LogHop2isSD2sOperation::initialAndAddFixedOrders(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew, size_t maxOrder, char TRANSOneBody)
{
    initialZeroOrder(oneBody, walker, walkerNew, TRANSOneBody);
    addFixedOrders(oneBody, walker, walkerNew, maxOrder, TRANSOneBody);
}

void LogHop2isSD2sOperation::initialAndAddDynamicOrders(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew, char TRANSOneBody)
{
    initialAndAddFixedOrders(oneBody, walker, walkerNew, baseTaylorOrder, TRANSOneBody);
    addDynamicOrders(oneBody, walker, walkerNew, TRANSOneBody);
}

void LogHop2isSD2sOperation::initialZeroOrder(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew, char TRANSOneBody)
{
    size_t L = walker.getL(); size_t Nup = walker.getNup(); size_t Ndn = walker.getNdn();
    if( wfUpTempNew.rank(0) != L || wfUpTempNew.rank(1) != Nup ) wfUpTempNew.resize( L, Nup );
    if( wfDnTempNew.rank(0) != L || wfDnTempNew.rank(1) != Ndn ) wfDnTempNew.resize( L, Ndn );
    if( wfUpTempOld.rank(0) != L || wfUpTempOld.rank(1) != Nup ) wfUpTempOld.resize( L, Nup );
    if( wfDnTempOld.rank(0) != L || wfDnTempOld.rank(1) != Ndn ) wfDnTempOld.resize( L, Ndn );

    walkerNew.wfUpRef() = walker.getWfUp();
    walkerNew.wfDnRef() = walker.getWfDn();
    wfUpTempOld = walker.getWfUp();
    wfDnTempOld = walker.getWfDn();

    currentOrder = 0;
}

void LogHop2isSD2sOperation::addFixedOrders(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew, size_t maxOrder, char TRANSOneBody)
{
    while( currentOrder < maxOrder )
    {
        BL_NAME(gmm)( oneBody.matrix, wfUpTempOld, wfUpTempNew, TRANSOneBody, 'N', 1.0/(currentOrder+1.0) );
        walkerNew.wfUpRef() += wfUpTempNew;
        wfUpTempOld = move( wfUpTempNew );

        BL_NAME(gmm)( oneBody.matrix, wfDnTempOld, wfDnTempNew, TRANSOneBody, 'N', 1.0/(currentOrder+1.0) );
        walkerNew.wfDnRef() += wfDnTempNew;
        wfDnTempOld = move( wfDnTempNew );

        currentOrder++;
    }
}

void LogHop2isSD2sOperation::addDynamicOrders(const LogHop2is &oneBody, const SD2s &walker, SD2s &walkerNew, char TRANSOneBody)
{
    do
    {
        BL_NAME(gmm)( oneBody.matrix, wfUpTempOld, wfUpTempNew, TRANSOneBody, 'N', 1.0/(currentOrder+1.0) );
        walkerNew.wfUpRef() += wfUpTempNew;
        wfUpTempOld = move( wfUpTempNew );

        BL_NAME(gmm)( oneBody.matrix, wfDnTempOld, wfDnTempNew, TRANSOneBody, 'N', 1.0/(currentOrder+1.0) );
        walkerNew.wfDnRef() += wfDnTempNew;
        wfDnTempOld = move( wfDnTempNew );

        currentOrder++;
    } while( ! isConverged() );
}

bool LogHop2isSD2sOperation::isConverged()
{
    size_t L = wfUpTempOld.rank(0); size_t Nup = wfUpTempOld.rank(1); size_t Ndn = wfDnTempOld.rank(1);

    for(size_t i = 0; i < Nup; ++i)
    {
        for(size_t j = 0; j < L; ++j)
        {
            if( abs( wfUpTempOld(j,i) ) > accuracy  ) return false;
        }
    }

    for(size_t i = 0; i < Ndn; ++i)
    {
        for(size_t j = 0; j < L; ++j)
        {
            if( abs( wfDnTempOld(j,i) ) > accuracy  ) return false;
        }
    }

    return true;
}

void LogHop2isSD2sOperation::clearWfTemp()
{
    wfUpTempOld.resize(0,0);
    wfDnTempOld.resize(0,0);

    wfUpTempNew.resize(0,0);
    wfDnTempNew.resize(0,0);
}

LogHop2isSD2sOperation::LogHop2isSD2sOperation(const LogHop2isSD2sOperation &x) { }

LogHop2isSD2sOperation &LogHop2isSD2sOperation::operator=(const LogHop2isSD2sOperation &x) { return *this; }
