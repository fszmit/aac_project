#ifndef GRAPH_CYCLES_HPP
#define GRAPH_CYCLES_HPP

#include "common.hpp"

bool isWeaklyConnected(const std::vector<std::vector<int>>* mat);
long long MatTrace(std::vector<std::vector<long long>> mat);
std::vector<std::vector<long long>> MatMul(std::vector<std::vector<long long>> mat1, std::vector<std::vector<long long>> mat2);
std::vector<std::vector<long long>> MatPow(std::vector<std::vector<long long>> mat, int exp);
Graph SubGraph(const Graph* const g, std::vector<int> mask, int* neighbours);
int MaxCycle(const Graph* const g, unsigned long long* no_cycles);
Graph PruneGraph(const Graph* const g);
int APXMaxCycle(const Graph* const g, unsigned long long* no_cycles);

#endif //GRAPH_CYCLES_HPP