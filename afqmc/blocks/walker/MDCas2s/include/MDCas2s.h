//
// Created by Hao Shi on 8/21/17.
//

#ifndef AFQMCLAB_MDCAS2S_H
#define AFQMCLAB_MDCAS2S_H

#include "../../../../../common/tensorHao/include/tensor_all.h"
#include "../../casBasis/include/casBasis.h"

//Multi Cas type Determinant, two spin species.

#ifdef MPI_HAO
class MDCas2s;
void MPIBcast(MDCas2s &buffer, int root=0,  const MPI_Comm& comm=MPI_COMM_WORLD);
#endif

class MDCas2s
{
 private:
    size_t L, Nup, Ndn;
    size_t MDSize, MDUpSize, MDDnSize, NupMax, NdnMax;
    std::complex<double> logw;
    std::vector< std::complex<double> > coe;
    std::vector< size_t > coeLinkUp, coeLinkDn;
    std::vector< std::vector<size_t> > occupancyUp, occupancyDn;
    tensor_hao::TensorHao<std::complex<double>,2> wfUp, wfDn;
    tensor_hao::TensorHao<std::complex<double>,2> RUp, RDn;

    std::vector<int> isParentUp, isParentDn;
    std::vector<size_t> treeUp, treeDn;
    std::vector<std::vector<size_t>> diffOrbitalUp, diffOrbitalDn;
    std::vector< size_t > parentIndexUp, parentIndexDn;

    bool wfCasIsCalculated;
    tensor_hao::TensorHao<std::complex<double>,2> wfCasUp, wfCasDn;

 public:
    MDCas2s();
    MDCas2s(const std::string &filename);
    MDCas2s(const MDCas2s& x);
    MDCas2s(MDCas2s&& x);
    ~MDCas2s();

    MDCas2s & operator  = (const MDCas2s& x);
    MDCas2s & operator  = (MDCas2s&& x);

    size_t getL() const;
    size_t getNup() const;
    size_t getNdn() const;
    size_t getMDSize() const;
    size_t getMDUpSize() const;
    size_t getMDDnSize() const;
    size_t getNupMax() const;
    size_t getNdnMax() const;
    const std::complex<double> &getLogw() const;
    const std::vector<std::complex<double>> &getCoe() const;
    const std::vector<size_t> &getCoeLinkUp() const;
    const std::vector<size_t> &getCoeLinkDn() const;
    const std::vector<std::vector<size_t>> &getOccupancyUp() const;
    const std::vector<std::vector<size_t>> &getOccupancyDn() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getWfUp() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getWfDn() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getRUp() const;
    const tensor_hao::TensorHao<std::complex<double>, 2> &getRDn() const;
    const std::vector<int> &getIsParentUp() const;
    const std::vector<int> &getIsParentDn() const;
    const std::vector<size_t> &getTreeUp() const;
    const std::vector<size_t> &getTreeDn() const;
    const std::vector<std::vector<size_t>> &getDiffOrbitalUp() const;
    const std::vector<std::vector<size_t>> &getDiffOrbitalDn() const;
    const std::vector<size_t> &getParentIndexUp() const;
    const std::vector<size_t> &getParentIndexDn() const;

    const tensor_hao::TensorHao<std::complex<double>, 2> &generateWfCasUp();
    const tensor_hao::TensorHao<std::complex<double>, 2> &generateWfCasDn();
    void stabilize();
    std::complex<double> normalize();

    void read(const std::string& filename);
    void write(const std::string& filename) const;
    void printParentNumber();
    double getMemory() const;

#ifdef MPI_HAO
    friend void MPIBcast(MDCas2s &buffer, int root,  const MPI_Comm& comm);
#endif

 private:
    void copy_deep(const MDCas2s &x);
    void move_deep(MDCas2s &x);

    void calculateWfCas();
    void normCoe();
    void scaleRMatrix();
    void setUpParentLinkStructure();
    void setParentIndex();
};

#endif //AFQMCLAB_MDCAS2S_H
