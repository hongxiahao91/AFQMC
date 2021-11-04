//
// Created by Hao Shi on 1/13/18.
//

#include "../include/logHop.h"

using namespace std;
using namespace tensor_hao;

LogHop::LogHop():logw(0) { }

LogHop::LogHop(size_t L):logw(0) { matrix.resize(L, L); }

LogHop::LogHop(const LogHop &x) { copy_deep(x); }

LogHop::LogHop(LogHop &&x) { move_deep(x); }

LogHop::~LogHop() { }

LogHop &LogHop::operator=(const LogHop &x)  { copy_deep(x); return *this; }

LogHop &LogHop::operator=(LogHop &&x) { move_deep(x); return *this; }

size_t LogHop::getL() const { return matrix.rank(0); }

double LogHop::getMemory() const
{
    return 16.0+matrix.getMemory();
}

void LogHop::copy_deep(const LogHop &x)
{
    logw = x.logw;
    matrix = x.matrix;
}

void LogHop::move_deep(LogHop &x)
{
    logw = x.logw;
    matrix = move( x.matrix );
}