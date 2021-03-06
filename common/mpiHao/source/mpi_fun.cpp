#include <iostream>
#include <climits>
#include "../include/mpi_fun.h"

#ifdef MPI_HAO
using namespace std;

void MPIInit(int& argc,char** & argv)
{
    MPI_Init(&argc,&argv);
}

void MPIInitFunnel(int& argc,char** & argv) 
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if(provided!=MPI_THREAD_FUNNELED)
    {
        if(MPIRank()==0) cout<<"WARNING!!!!The Marchine does not supply MPI_THREAD_FUNNELED, only supply to: "
                                  <<provided<<"-th level!"<<endl;
    }
}

int MPISize(const MPI_Comm& comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    return size;
}

int MPIRank(const MPI_Comm& comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

void MPIBarrier(const MPI_Comm& comm)
{
    MPI_Barrier(comm);
}

void MPIFinalize()
{
    MPI_Finalize();
}

void MPIBcast(bool &buffer, int root, MPI_Comm const &comm)
{
    MPI_Bcast(&buffer,sizeof(buffer), MPI_BYTE, root, comm);
}

void MPIBcast(int& buffer, int root, const MPI_Comm& comm )
{
    MPI_Bcast(&buffer, 1, MPI_INT, root, comm); 
}

void MPIBcast(long& buffer, int root, const MPI_Comm& comm )
{
    MPI_Bcast(&buffer, 1, MPI_LONG, root, comm);
}

void MPIBcast(long long& buffer, int root, const MPI_Comm& comm )
{
    MPI_Bcast(&buffer, 1, MPI_LONG_LONG, root, comm);
}

void MPIBcast(size_t& buffer, int root, const MPI_Comm& comm )
{
    MPI_Bcast(&buffer,sizeof(buffer), MPI_BYTE, root, comm);
}

void MPIBcast(float& buffer, int root, const MPI_Comm& comm )
{
    MPI_Bcast(&buffer, 1, MPI_FLOAT, root, comm);
}

void MPIBcast(double& buffer, int root, const MPI_Comm& comm )
{
    MPI_Bcast(&buffer, 1, MPI_DOUBLE, root, comm);
}

void MPIBcast(complex<float>& buffer, int root, const MPI_Comm& comm )
{
    MPI_Bcast(&buffer, 1, MPI_COMPLEX, root, comm);
}

void MPIBcast(complex<double>& buffer, int root, const MPI_Comm& comm )
{
    MPI_Bcast(&buffer, 1, MPI_DOUBLE_COMPLEX, root, comm);
}

void MPIBcast(string &buffer, int root, MPI_Comm const &comm)
{
    size_t bufferSize=buffer.size();
    MPI_Bcast(&bufferSize, sizeof(bufferSize), MPI_BYTE, root, comm);
    buffer.resize(bufferSize);
    MPI_Bcast(&buffer[0], buffer.size(), MPI_CHAR, root, comm);
}

void MPIBcast(size_t count, bool* buffer, int root, const MPI_Comm& comm )
{
    if( count==0 ) return;

    bool tmp(true);
    size_t chunkSize = INT_MAX;
    size_t chunkNumber = (count-1)/chunkSize;
    bool * bufferPerChunk = buffer;
    size_t currentChunkSize;
    for(size_t i = 0; i <= chunkNumber; ++i)
    {
        currentChunkSize = (i == chunkNumber) ? count - chunkSize * chunkNumber : chunkSize;
        MPI_Bcast(bufferPerChunk, currentChunkSize*sizeof(tmp), MPI_BYTE, root, comm);
        bufferPerChunk += currentChunkSize;
    }
}

void MPIBcast(size_t count, int* buffer, int root, const MPI_Comm& comm )
{
    if( count==0 ) return;

    size_t chunkSize = INT_MAX;
    size_t chunkNumber = (count-1)/chunkSize;
    int * bufferPerChunk = buffer;
    size_t currentChunkSize;
    for(size_t i = 0; i <= chunkNumber; ++i)
    {
        currentChunkSize = (i == chunkNumber) ? count - chunkSize * chunkNumber : chunkSize;
        MPI_Bcast(bufferPerChunk, currentChunkSize, MPI_INT, root, comm);
        bufferPerChunk += currentChunkSize;
    }
}

void MPIBcast(size_t count, size_t * buffer, int root, const MPI_Comm& comm )
{
    if( count==0 ) return;

    size_t chunkSize = INT_MAX;
    size_t chunkNumber = (count-1)/chunkSize;
    size_t * bufferPerChunk = buffer;
    size_t currentChunkSize;
    for(size_t i = 0; i <= chunkNumber; ++i)
    {
        currentChunkSize = (i == chunkNumber) ? count - chunkSize * chunkNumber : chunkSize;
        MPI_Bcast(bufferPerChunk, currentChunkSize*sizeof(size_t), MPI_BYTE, root, comm);
        bufferPerChunk += currentChunkSize;
    }
}

void MPIBcast(size_t count, double *buffer, int root, const MPI_Comm &comm)
{
    if( count==0 ) return;

    size_t chunkSize = INT_MAX;
    size_t chunkNumber = (count-1)/chunkSize;
    double * bufferPerChunk = buffer;
    size_t currentChunkSize;
    for(size_t i = 0; i <= chunkNumber; ++i)
    {
        currentChunkSize = (i == chunkNumber) ? count - chunkSize * chunkNumber : chunkSize;
        MPI_Bcast(bufferPerChunk, currentChunkSize, MPI_DOUBLE, root, comm);
        bufferPerChunk += currentChunkSize;
    }
}

void MPIBcast(size_t count, complex<double>* buffer, int root, const MPI_Comm& comm )
{
    if( count==0 ) return;

    size_t chunkSize = INT_MAX;
    size_t chunkNumber = (count-1)/chunkSize;
    complex<double> * bufferPerChunk = buffer;
    size_t currentChunkSize;
    for(size_t i = 0; i <= chunkNumber; ++i)
    {
        currentChunkSize = (i == chunkNumber) ? count - chunkSize * chunkNumber : chunkSize;
        MPI_Bcast(bufferPerChunk, currentChunkSize, MPI_DOUBLE_COMPLEX, root, comm);
        bufferPerChunk += currentChunkSize;
    }
}

void MPIBcast(vector<int> &buffer, int root, MPI_Comm const &comm)
{
    size_t bufferSize = buffer.size();
    MPIBcast(bufferSize, root, comm);
    buffer.resize(bufferSize);
    MPIBcast(buffer.size(), buffer.data(), root, comm);
}

void MPIBcast(vector<size_t> &buffer, int root, MPI_Comm const &comm)
{
    size_t bufferSize = buffer.size();
    MPIBcast(bufferSize, root, comm);
    buffer.resize(bufferSize);
    MPIBcast(buffer.size(), buffer.data(), root, comm);
}

void MPIBcast(vector<double> &buffer, int root, MPI_Comm const &comm)
{
    size_t bufferSize = buffer.size();
    MPIBcast(bufferSize, root, comm);
    buffer.resize(bufferSize);
    MPIBcast(buffer.size(), buffer.data(), root, comm);
}

void MPIBcast(vector<complex<double>> &buffer, int root, MPI_Comm const &comm)
{
    size_t bufferSize = buffer.size();
    MPIBcast(bufferSize, root, comm);
    buffer.resize(bufferSize);
    MPIBcast(buffer.size(), buffer.data(), root, comm);
}

void MPIReduce(const size_t &sendbuf, size_t &recvbuf, MPI_Op op, int root, MPI_Comm const &comm)
{
    long long sendbufTmp = sendbuf; long long recvbufTmp=0;
    MPI_Reduce(&sendbufTmp, &recvbufTmp, 1, MPI_LONG_LONG, op, root, comm);
    recvbuf = recvbufTmp;
}

void MPIReduce(const double &sendbuf, double &recvbuf, MPI_Op op, int root, MPI_Comm const &comm)
{
    MPI_Reduce(&sendbuf, &recvbuf, 1, MPI_DOUBLE, op, root, comm);
}

void MPIAllreduce(const double &sendbuf, double &recvbuf, MPI_Op op, const MPI_Comm& comm)
{
    MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_DOUBLE, op, comm);
}

void MPIAllreduce(const complex<double> &sendbuf, complex<double> &recvbuf, MPI_Op op, const MPI_Comm& comm)
{
    MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_DOUBLE_COMPLEX, op, comm);
}

void MPIAllreduce(size_t count, const double *sendbuf, double *recvbuf, MPI_Op op, const MPI_Comm& comm)
{
    if( count==0 ) return;

    size_t chunkSize = INT_MAX;
    size_t chunkNumber = (count-1)/chunkSize;
    const double * sendbufPerChunk = sendbuf;
    double * recvbufPerChunk = recvbuf;
    size_t currentChunkSize;
    for(size_t i = 0; i <= chunkNumber; ++i)
    {
        currentChunkSize = (i == chunkNumber) ? count - chunkSize * chunkNumber : chunkSize;
        MPI_Allreduce(sendbufPerChunk, recvbufPerChunk, currentChunkSize, MPI_DOUBLE, op, comm);
        sendbufPerChunk += currentChunkSize; recvbufPerChunk += currentChunkSize;
    }
}

void MPIAllreduce(size_t count, const complex<double> *sendbuf, complex<double> *recvbuf, MPI_Op op, const MPI_Comm& comm)
{
    if( count==0 ) return;

    size_t chunkSize = INT_MAX;
    size_t chunkNumber = (count-1)/chunkSize;
    const complex<double> * sendbufPerChunk = sendbuf;
    complex<double> * recvbufPerChunk = recvbuf;
    size_t currentChunkSize;
    for(size_t i = 0; i <= chunkNumber; ++i)
    {
        currentChunkSize = (i == chunkNumber) ? count - chunkSize * chunkNumber : chunkSize;
        MPI_Allreduce(sendbufPerChunk, recvbufPerChunk, currentChunkSize, MPI_DOUBLE_COMPLEX, op, comm);
        sendbufPerChunk += currentChunkSize; recvbufPerChunk += currentChunkSize;
    }
}

void MPIGather(const double& sendbuf, double* recvbuf, int root, const MPI_Comm& comm)
{
    MPI_Gather(&sendbuf,1,MPI_DOUBLE,recvbuf,1,MPI_DOUBLE,root,comm);
}

void MPIGather(const complex<double>& sendbuf, complex<double>* recvbuf, int root, const MPI_Comm& comm)
{
    MPI_Gather(&sendbuf,1,MPI_DOUBLE_COMPLEX,recvbuf,1,MPI_DOUBLE_COMPLEX,root,comm);
}


int MPISum(const int& sendbuf, int root,  const MPI_Comm& comm)
{
    int recvbuf=0;
    MPI_Reduce(&sendbuf, &recvbuf, 1 , MPI_INT, MPI_SUM, root, comm);
    return recvbuf;
}

long MPISum(const long& sendbuf, int root,  const MPI_Comm& comm)
{
    long recvbuf=0;
    MPI_Reduce(&sendbuf, &recvbuf, 1 , MPI_LONG, MPI_SUM, root, comm);
    return recvbuf;
}

long long MPISum(const long long& sendbuf, int root,  const MPI_Comm& comm)
{
    long long recvbuf=0;
    MPI_Reduce(&sendbuf, &recvbuf, 1 , MPI_LONG_LONG, MPI_SUM, root, comm);
    return recvbuf;
}

size_t MPISum(const size_t& sendbuf, int root,  const MPI_Comm& comm)
{
    //If sendbuf=600, MPIRank=4, recvbuf=2144! (Use Intel compiler, openmpi, on my Macbook Pro)
    //Size_t is not support well in MPI? Use long long instead!
    long long sendbufTmp = sendbuf; long long recvbufTmp=0;
    MPI_Reduce(&sendbufTmp, &recvbufTmp, 1, MPI_LONG_LONG, MPI_SUM, root, comm);
    return recvbufTmp;
}

float MPISum(const float& sendbuf, int root,  const MPI_Comm& comm)
{
    float recvbuf=0;
    MPI_Reduce(&sendbuf, &recvbuf, 1 , MPI_FLOAT, MPI_SUM, root, comm);
    return recvbuf;
}

double MPISum(const double& sendbuf, int root,  const MPI_Comm& comm)
{
    double recvbuf=0;
    MPI_Reduce(&sendbuf, &recvbuf, 1 , MPI_DOUBLE, MPI_SUM, root, comm);
    return recvbuf;
}

complex<float> MPISum(const complex<float>& sendbuf, int root,  const MPI_Comm& comm)
{
    complex<float> recvbuf={0,0};
    MPI_Reduce(&sendbuf, &recvbuf, 1 , MPI_COMPLEX, MPI_SUM, root, comm);
    return recvbuf;
}

complex<double> MPISum(const complex<double>& sendbuf, int root,  const MPI_Comm& comm)
{
    complex<double> recvbuf={0,0};
    MPI_Reduce(&sendbuf, &recvbuf, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, root, comm);
    return recvbuf;
}

void MPISum(size_t count, const double* sendbuf, double* recvbuf, int root, const MPI_Comm& comm)
{
    if( count==0 ) return;

    size_t chunkSize = INT_MAX;
    size_t chunkNumber = (count-1)/chunkSize;
    const double * sendbufPerChunk = sendbuf;
    double * recvbufPerChunk = recvbuf;
    size_t currentChunkSize;
    for(size_t i = 0; i <= chunkNumber; ++i)
    {
        currentChunkSize = (i == chunkNumber) ? count - chunkSize * chunkNumber : chunkSize;
        MPI_Reduce(sendbufPerChunk, recvbufPerChunk, currentChunkSize , MPI_DOUBLE, MPI_SUM, root, comm);
        sendbufPerChunk += currentChunkSize; recvbufPerChunk += currentChunkSize;
    }
}

void MPISum(size_t count, const complex<double>* sendbuf, complex<double>*recvbuf, int root, const MPI_Comm& comm)
{
    if( count==0 ) return;

    size_t chunkSize = INT_MAX;
    size_t chunkNumber = (count-1)/chunkSize;
    const complex<double> * sendbufPerChunk = sendbuf;
    complex<double> * recvbufPerChunk = recvbuf;
    size_t currentChunkSize;
    for(size_t i = 0; i <= chunkNumber; ++i)
    {
        currentChunkSize = (i == chunkNumber) ? count - chunkSize * chunkNumber : chunkSize;
        MPI_Reduce(sendbufPerChunk, recvbufPerChunk, currentChunkSize , MPI_DOUBLE_COMPLEX, MPI_SUM, root, comm);
        sendbufPerChunk += currentChunkSize; recvbufPerChunk += currentChunkSize;
    }
}
#else
void MPIInit(int& argc,char** & argv){return;}
void MPIInitFunnel(int& argc,char** & argv) {return;}
int MPISize() {return 1;}
int MPIRank() {return 0;}
void MPIBarrier() {return;}
void MPIFinalize() {return;}
#endif
