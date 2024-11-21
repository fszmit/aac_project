#ifndef GRAPH_CYCLES_HPP
#define GRAPH_CYCLES_HPP

#include "common.hpp"

bool isWeaklyConnected(const std::vector<std::vector<int>> * mat);
int MatTrace(std::vector<std::vector<int>> mat);
std::vector<std::vector<int>> MatMul(std::vector<std::vector<int>> mat1, std::vector<std::vector<int>> mat2);
std::vector<std::vector<int>> MatPow(std::vector<std::vector<int>> mat, int exp);
Graph SubGraph(const Graph * const g, std::vector<int> mask, int * neighbours);
int MaxCycle(Graph * g, int * no_cycles);

#endif //GRAPH_CYCLES_HPP