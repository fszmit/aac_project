#include "approximation.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
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
            if (adjacencyMatrix[i][j] >= 1) {  // >= 1 to handle multigraphss
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
int calculate_cycle_length(const std::vector<int>& cycle,
                           const std::vector<std::vector<int>>& weights) {
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
bool three_opt_move(std::vector<int>& cycle, const std::vector<std::vector<int>>& weights,
                    bool is_symmetric = true) {
    int n = cycle.size();
    bool improved = false;
    int best_gain = 0;

    // Store the best move
    int best_i = 0, best_j = 0, best_k = 0;
    int best_perm_index = 0;

    int i_new = -1;
    int i_next_new = -1;
    int j_new = -1;
    int j_next_new = -1;
    int k_new = -1;
    int k_next_new = -1;

    // Try all possible combinations of three edges
    for (int i = 0; i < n - 2; i++) {
        for (int j = i + 2; j < n - 1; j++) {  // todo revert to j = i + 1
            for (int k = j + 2; k < n; k++) {  // todo revert to k = j + 1
                // std::cout << "combination" << std::endl;

                int i_next = (i + 1) % n;
                int j_next = (j + 1) % n;
                int k_next = (k + 1) % n;

                // given cycle with edges -> A->B, a path "1" B...C, C->D, a path "2" D...E, E->F, a
                // path "3" F...A we will delete edges A->B, C->D, E->F

                // starting from A, we can select the following vertices that are not yet connected
                // B, C, D, E

                // given

                // edges
                // a = [i, i_next]
                // b = [j, j_next]
                // c = [k, k_next]

                int old_weight = weights[cycle[i]][cycle[i_next]] +
                                 weights[cycle[j]][cycle[j_next]] +
                                 weights[cycle[k]][cycle[k_next]];

                std::vector<int> selectable_vertices;

                if (is_symmetric) {
                    selectable_vertices = {i_next, j, j_next, k};
                } else {
                    selectable_vertices = {i_next, j_next};
                }

                for (int f = 0; f < selectable_vertices.size(); f++) {
                    int first_added_edge_start = i;
                    int first_added_edge_end = selectable_vertices[f];

                    int second_added_edge_start = -1;
                    int second_added_edge_end = -1;

                    int third_added_edge_start = -1;
                    int third_added_edge_end = k_next;

                    if (first_added_edge_end == i_next) {
                        second_added_edge_start = j;
                    } else if (first_added_edge_end == j) {
                        second_added_edge_start = i_next;
                    } else if (first_added_edge_end == j_next) {
                        second_added_edge_start = k;
                    } else if (first_added_edge_end == k) {
                        second_added_edge_start = j_next;
                    }

                    std::vector<int> second_selectable_vertices;

                    if (is_symmetric) {
                        if (first_added_edge_end == i_next || first_added_edge_end == j) {
                            second_selectable_vertices = {j_next, k};
                        } else {
                            second_selectable_vertices = {i_next, j};
                        }

                    } else {
                        if (first_added_edge_end == i_next) {
                            second_selectable_vertices = {j_next};
                        } else {
                            second_selectable_vertices = {i_next};
                        }
                    }

                    for (int g = 0; g < second_selectable_vertices.size(); g++) {
                        second_added_edge_end = second_selectable_vertices[g];

                        if (second_added_edge_end == i_next) {
                            third_added_edge_start = j;
                        } else if (second_added_edge_end == j) {
                            third_added_edge_start = i_next;
                        } else if (second_added_edge_end == j_next) {
                            third_added_edge_start = k;
                        } else if (second_added_edge_end == k) {
                            third_added_edge_start = j_next;
                        }
                    }

                    int new_weight =
                        weights[cycle[first_added_edge_start]][cycle[first_added_edge_end]] +
                        weights[cycle[second_added_edge_start]][cycle[second_added_edge_end]] +
                        weights[cycle[third_added_edge_start]][cycle[third_added_edge_end]];

                    // old 50, new 49
                    int gain = old_weight - new_weight;
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_i = i;
                        best_j = j;
                        best_k = k;

                        i_new = first_added_edge_start;
                        i_next_new = first_added_edge_end;
                        j_new = second_added_edge_start;
                        j_next_new = second_added_edge_end;
                        k_new = third_added_edge_start;
                        k_next_new = third_added_edge_end;

                        improved = true;
                        // std::cout << "improvement " << improved << " gain " << gain << " nw "
                        //           << new_weight << " old " << old_weight << std::endl;
                    }
                }

                // Current edges
                // Try all permutations of the three segments
                std::vector<int> segments = {j, k, i};
                int perm_index = 0;
                do {
                    // std::cout << "permutation" << std::endl;

                    std::cout << char(i + int('a')) << "->" << char(segments[0] + int('a')) << " "
                              << char(j + int('a')) << "->" << char(segments[1] + int('a')) << " "
                              << char(k + int('a')) << "->" << char(segments[2] + int('a'))
                              << std::endl;

                    perm_index++;
                } while (std::next_permutation(segments.begin(), segments.end()));
                break;
            }
        }
    }

    // Apply the best move if found
    if (improved) {
        std::vector<int> new_cycle = cycle;

        int best_i_next = (best_i + 1) % n;
        int best_j_next = (best_j + 1) % n;
        int best_k_next = (best_k + 1) % n;

        // for(int i : endpoints){
        //     std::cout << i << " ";
        // }
        // std::cout<<std::endl;
        //  Get the best permutation
        // for (int i = 0; i < best_perm_index; i++) {
        //     // std::cout << "permute to best permutation" << std::endl;
        //     std::next_permutation(segments.begin(), segments.end());
        // }

        new_cycle[best_i] = cycle[i_new];
        new_cycle[best_i_next] = cycle[i_next_new];
        new_cycle[best_j] = cycle[j_new];
        new_cycle[best_j_next] = cycle[j_next_new];
        new_cycle[best_k] = cycle[k_new];
        new_cycle[best_k_next] = cycle[k_next_new];

        // // Reconstruct the cycle with the new permutation
        // int pos = 0;
        // for (int i = 0; i < 3; i++) {
        //     int start = segments[i];
        //     int end = endpoints[i];
        //     while (start != end) {
        //         new_cycle[(pos++) % n] = cycle[start];
        //         start = (start + 1) % n;
        //     }
        // }

        cycle = new_cycle;
    }

    return improved;
}

int get_random_start_vertex(int n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, n - 1);
    return dis(gen);
}

int get_highest_degree_vertex(const std::vector<std::vector<int>>& graph) {
    int n = graph.size();
    int max_degree = -1;
    int max_degree_vertex = 0;

    for (int i = 0; i < n; i++) {
        int degree = 0;
        for (int j = 0; j < n; j++) {
            if (graph[i][j] >= 1) {  // >= 1 to handle multigraphs
                degree++;
            }
        }
        if (degree > max_degree) {
            max_degree = degree;
            max_degree_vertex = i;
        }
    }

    return max_degree_vertex;
}

int get_consistent_start_vertex(int n) {
    return 28;
    /*static int call_count = 0;
    static int last_vertex = -1;
    call_count++;

    std::cout << "Start vertex selection call #" << call_count;

    if (last_vertex != -1) {
        std::cout << " (previous: " << last_vertex << ")";
    }
    std::cout << std::endl;

    last_vertex = call_count % n;
    return last_vertex;*/
}

int edges_to_add_in_cycle(std::vector<int> cycle, WeightedGraph& weighted_graph) {
    int edges_to_add = 0;
    for (size_t i = 0; i < cycle.size(); i++) {
        int from = cycle[i];
        int to = cycle[(i + 1) % cycle.size()];
        if (weighted_graph.weightMatrix[from][to] ==
            1) {  // If edge doesn't exist in original graph
            edges_to_add++;
        }
    }
    return edges_to_add;
}

int hamiltonian_completion_approximation(const std::vector<std::vector<int>>& graph) {
    //   to complete weighted graph
    WeightedGraph weighted_graph = transform_to_complete_weighted_graph(graph);
    int n = weighted_graph.vertices;

    // Get vertex with highest degree as start vertex
    // int start_vertex = get_highest_degree_vertex(graph);

    // int start_vertex = get_random_start_vertex(weighted_graph.vertices);

    // Use consistent vertex selection for debugging
    int start_vertex = get_consistent_start_vertex(n);

    // Get initial cycle using nearest neighbor
    std::vector<int> cycle = nearest_neighbor(weighted_graph.weightMatrix, start_vertex);
    std::cout << edges_to_add_in_cycle(cycle, weighted_graph) << std::endl;

    // Apply 3-Opt until no improvement is found
    bool improved;
    int iters = 0;

    // this makes it iterated three_opt which is slower
    do {
        improved = three_opt_move(cycle, weighted_graph.weightMatrix);
        iters++;
        // std::cout << iters << std::endl;

        std::cout << edges_to_add_in_cycle(cycle, weighted_graph) << std::endl;
    } while (improved);

    // Count the number of edges that need to be added (edges with weight 1 in the cycle)
    return edges_to_add_in_cycle(cycle, weighted_graph);
    // int edges_to_add = 0;
    // for (size_t i = 0; i < cycle.size(); i++) {
    //     int from = cycle[i];
    //     int to = cycle[(i + 1) % cycle.size()];
    //     if (weighted_graph.weightMatrix[from][to] == 1) {  // If edge doesn't exist in original
    //     graph
    //         edges_to_add++;
    //     }
    // }

    // return edges_to_add;
}