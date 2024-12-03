#include "approximation.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <queue>


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
            if (adjacencyMatrix[i][j] >= 1) { // >= 1 to handle multigraphss
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

// Helper function to calculate cycle length
int calculate_cycle_length(const std::vector<int>& cycle, const std::vector<std::vector<int>>& weights) {
    int length = 0;
    for (size_t i = 0; i < cycle.size(); i++) {
        length += weights[cycle[i]][cycle[(i + 1) % cycle.size()]];
    }
    return length;
}

// Nearest neighbor algorithm for initial cycle
std::vector<int> nearest_neighbor(const std::vector<std::vector<int>>& weights, int start_vertex) {
    int n = weights.size();
    std::vector<bool> visited(n, false);
    std::vector<int> cycle;
    cycle.reserve(n);
    
    int current = start_vertex;
    cycle.push_back(current);
    visited[current] = true;
    
    while (cycle.size() < n) {
        int nearest = -1;
        int min_weight = std::numeric_limits<int>::max();
        
        for (int i = 0; i < n; i++) {
            if (!visited[i] && weights[current][i] < min_weight) {
                min_weight = weights[current][i];
                nearest = i;
            }
        }
        
        current = nearest;
        cycle.push_back(current);
        visited[current] = true;
    }
    
    return cycle;
}

// 3-Opt move without path reversal
bool three_opt_move(std::vector<int>& cycle, const std::vector<std::vector<int>>& weights) {
    int n = cycle.size();
    bool improved = false;
    int best_gain = 0;
    
    // Store the best move
    int best_i = 0, best_j = 0, best_k = 0;
    int best_perm_index = 0;
    
    // Try all possible combinations of three edges
    for (int i = 0; i < n - 2; i++) { 
        for (int j = i + 1; j < n - 1; j++) {
            for (int k = j + 1; k < n; k++) {
                int i_next = (i + 1) % n;
                int j_next = (j + 1) % n;
                int k_next = (k + 1) % n;
                
                // Current edges
                int old_weight = weights[cycle[i]][cycle[i_next]] +
                               weights[cycle[j]][cycle[j_next]] +
                               weights[cycle[k]][cycle[k_next]];
                
                // Try all permutations of the three segments
                std::vector<int> segments = {i_next, j_next, k_next};
                int perm_index = 0;
                do {
                    int new_weight = weights[cycle[i]][cycle[segments[0]]] +
                                   weights[cycle[j]][cycle[segments[1]]] +
                                   weights[cycle[k]][cycle[segments[2]]];
                    
                    int gain = old_weight - new_weight;
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_i = i;
                        best_j = j;
                        best_k = k;
                        best_perm_index = perm_index;
                        improved = true;
                    }
                    perm_index++;
                } while (std::next_permutation(segments.begin(), segments.end()));
            }
        }
    }
    
    // Apply the best move if found
    if (improved) {
        std::vector<int> new_cycle = cycle;
        std::vector<int> segments = {(best_i + 1) % n, (best_j + 1) % n, (best_k + 1) % n};
        std::vector<int> endpoints = {best_j % n, best_k % n, (best_i == 0 ? n - 1 : best_i) % n};
        
        // Get the best permutation
        for (int i = 0; i < best_perm_index; i++) {
            std::next_permutation(segments.begin(), segments.end());
        }
        
        // Reconstruct the cycle with the new permutation
        int pos = 0;
        for (int i = 0; i < 3; i++) {
            int start = segments[i];
            int end = endpoints[i];
            while (start != end) {
                new_cycle[pos++] = cycle[start];
                start = (start + 1) % n;
            }
        }
        
        cycle = new_cycle;
    }
    
    return improved;
}


int hamiltonian_completion_approximation(const std::vector<std::vector<int>>& graph) {
    // Transform to complete weighted graph
    WeightedGraph weighted_graph = transform_to_complete_weighted_graph(graph);
    int n = weighted_graph.vertices;
    
    // Generate random start vertex
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, n - 1);
    int start_vertex = dis(gen);
    
    // Get initial cycle using nearest neighbor
    std::vector<int> cycle = nearest_neighbor(weighted_graph.weightMatrix, start_vertex);
    
    // Apply 3-Opt until no improvement is found
    bool improved;
    int i = 0;
    do {
        improved = three_opt_move(cycle, weighted_graph.weightMatrix);
        i++;
    } while (improved);
    
    // Count the number of edges that need to be added (edges with weight 1 in the cycle)
    int edges_to_add = 0;
    for (size_t i = 0; i < cycle.size(); i++) {
        int from = cycle[i];
        int to = cycle[(i + 1) % cycle.size()];
        if (weighted_graph.weightMatrix[from][to] == 1) {  // If edge doesn't exist in original graph
            edges_to_add++;
        }
    }
    
    return edges_to_add;
}