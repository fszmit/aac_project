#pragma once
#include <vector>

struct WeightedGraph {
    int vertices;
    std::vector<std::vector<int>> weightMatrix;
};

WeightedGraph transform_to_complete_weighted_graph(
    const std::vector<std::vector<int>>& adjacencyMatrix);
void print_weighted_graph(const WeightedGraph& graph);
int hamiltonian_completion_approximation(const std::vector<std::vector<int>>& graph);