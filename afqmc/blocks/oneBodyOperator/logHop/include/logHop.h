//
// Created by Hao Shi on 1/13/18.
//

#ifndef AFQMCLAB_LOGHOP_H
#define AFQMCLAB_LOGHOP_H

#include "../../../../../common/tensorHao/include/tensor_all.h"

//One body operator: two identical spin species.
//The operator is exp( M ) = 1 + M + 1/(2!) M^2 + 1/(3!) M^3 + ...
//M matrix is stored in LogHop class.

class LogHop
{
 public:
    std::complex<double> logw;
    tensor_hao::TensorHao<std::complex<double>,2> matrix;

    LogHop();
    LogHop(size_t L);
    LogHop(const LogHop& x);
    LogHop(LogHop&& x);
    ~LogHop();

    LogHop & operator  = (const LogHop& x);
    LogHop & operator  = (LogHop&& x);

    size_t getL() const;
    double getMemory() const;

 private:
    void copy_deep(const LogHop &x);
    void move_deep(LogHop &x);
};


#endif //AFQMCLAB_LOGHOP_H
