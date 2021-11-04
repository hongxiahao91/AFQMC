//
// Created by Hao Shi on 1/13/18.
//

#ifndef AFQMCLAB_LOGHOPSDOPERATION_H
#define AFQMCLAB_LOGHOPSDOPERATION_H

#include "../../../oneBodyOperator/logHop/include/logHop.h"
#include "../../../walker/SD/include/SD.h"

// Three different type of operation.
// flag: fixedOrder, taylorOrder>0. Always use the input taylorOrder order.
// flag: fixedOrder, taylorOrder=0. Determine the taylor order from accuracy and baseTaylorOrder and always use it.
// flag: dynamicOrder. Determine the taylor order from accuracy and baseTaylorOrder in all the operations.

class LogHopSDOperation
{
 private:
    std::string flag; //"fixedOrder", "dynamicOrder"
    size_t taylorOrder;
    double accuracy;
    size_t baseTaylorOrder;

    size_t operationNumber;
    size_t minTaylorOrder, maxTaylorOrder, totalTaylorOrder;

    size_t currentOrder;
    tensor_hao::TensorHao<std::complex<double>,2> wfTempOld, wfTempNew;

 public:
    LogHopSDOperation(std::string flag="dynamicOrder", size_t taylorOrder=0, double accuracy=1e-10, size_t baseTaylorOrder=3);
    ~LogHopSDOperation();

    void reset(std::string flag="dynamicOrder", size_t taylorOrder=0, double accuracy=1e-10, size_t baseTaylorOrder=3);
    void applyToRight(const LogHop &oneBody, const SD &walker, SD &walkerNew);
    void applyToLeft(const LogHop &oneBody, const SD &walker, SD &walkerNew);
    void print();
    size_t getCurrentOrder() const;

 private:
    void checkAndResize(const LogHop &oneBody, const SD &walker, SD &walkerNew) const;

    void addOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody);
    void determinantAndAddFixedOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody);
    void initialAndAddFixedOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, size_t maxOrder, char TRANSOneBody);
    void initialAndAddDynamicOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody);

    void initialZeroOrder(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody);
    void addFixedOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, size_t maxOrder, char TRANSOneBody);
    void addDynamicOrders(const LogHop &oneBody, const SD &walker, SD &walkerNew, char TRANSOneBody);
    bool isConverged();
    void clearWfTemp();

    LogHopSDOperation(const LogHopSDOperation& x);
    LogHopSDOperation & operator  = (const LogHopSDOperation& x);
};


#endif //AFQMCLAB_LOGHOPSDOPERATION_H
