#include "approximation.h"

#include <iostream>
#include <vector>

WeightedGraph transform_to_complete_weighted_graph(
    const std::vector<std::vector<int>>& adjacencyMatrix) {
    int n = adjacencyMatrix.size();
    WeightedGraph result;
    result.vertices = n;
    result.weightMatrix =
        std::vector<std::vector<int>>(n, std::vector<int>(n, 1));  // Initialize with weight 1

    // Copy existing edges with weight 0
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (adjacencyMatrix[i][j] == 1) {
                result.weightMatrix[i][j] = 0;
            }
        }
    }

    return result;
}

void print_weighted_graph(const WeightedGraph& graph) {
    std::cout << "Number of vertices: " << graph.vertices << std::endl;
    std::cout << "Weight Matrix:" << std::endl;

    for (const auto& row : graph.weightMatrix) {
        for (int weight : row) {
            std::cout << weight << " ";
        }
        std::cout << std::endl;
    }
}

WeightedGraph hamiltonian_completion_approximation(const std::vector<std::vector<int>>& graph) {
    return transform_to_complete_weighted_graph(graph);
}
