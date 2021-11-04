//
// Created by Hao Shi on 8/21/17.
//

#ifndef AFQMCLAB_CASBASIS_H
#define AFQMCLAB_CASBASIS_H

#include "../../../../../common/common.h"

int checkOccupancyOrder(const std::vector< std::vector<size_t> > &occupancy);

size_t calculateNumberOfDifferentCasOrbitals(const std::vector<size_t> &o1, const std::vector<size_t> &o2);

std::tuple<size_t, size_t> findMinDiffNumberAndFrequency(const std::vector<size_t> &diff);

int permuteOrbitals(const std::vector<size_t> &parent, std::vector<size_t> &child);

std::vector<size_t> calculateDiffOrbitals(const std::vector<size_t> &parent, const std::vector<size_t> &child);

void buildParentTree(std::vector<int> &isParent, std::vector<size_t> &tree, const std::vector<std::vector<size_t> > &occupancy);

std::vector<int> permuteOrbitalsInOccupancy(const std::vector<int> &isParent, const std::vector<size_t> &tree, std::vector<std::vector<size_t> > &occupancy);

std::vector<std::vector<size_t>> calculateDiffOrbitalsInOccupancy(const std::vector<int> &isParent,
                                                                  const std::vector<size_t> &tree,
                                                                  const std::vector<std::vector<size_t> > &occupancy);

std::vector<int> sortOrbitalsInOccupancy(std::vector<std::vector<size_t>> &occupancySort, const std::vector<std::vector<size_t>> &occupancy);

void fillCasMatrixByRow(tensor_hao::TensorHao<std::complex<double>, 2> &smallMatrix,
                        const tensor_hao::TensorHao<std::complex<double>, 2> &fullMatrix,
                        const std::vector<size_t> &row);

void fillCasMatrixByCol(tensor_hao::TensorHao<std::complex<double>, 2> &smallMatrix,
                        const tensor_hao::TensorHao<std::complex<double>, 2> &fullMatrix,
                        const std::vector<size_t> &col);

void getRowDifference(tensor_hao::TensorHao<std::complex<double>, 2> &V,
                      const tensor_hao::TensorHao<std::complex<double>, 2> &fullMatrix,
                      const std::vector<size_t> &row1,
                      const std::vector<size_t> &row2,
                      const std::vector<size_t> &diff);
#endif //AFQMCLAB_CASBASIS_H