#include "approximation.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <set>
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

int edges_to_add_in_cycle(const std::vector<int> cycle,
                          const std::vector<std::vector<int>>& weights) {
    int edges_to_add = 0;
    for (size_t i = 0; i < cycle.size(); i++) {
        int from = cycle[i];
        int to = cycle[(i + 1) % cycle.size()];
        if (weights[from][to] == 1) {  // If edge doesn't exist in original graph
            edges_to_add++;
        }
    }
    return edges_to_add;
}

void print_cycle(const std::vector<int> cycle) {
    for (int i : cycle) {
        std::cout << i << " -> ";
    }
    std::cout << std::endl;
}

bool contains(const std::set<int> s, int x) { return s.find(x) != s.end(); }

// 3-Opt move without path reversal
bool three_opt_move(std::vector<int>& cycle, const std::vector<std::vector<int>>& weights,
                    bool is_symmetric = true) {
    int n = cycle.size();
    bool improved = false;
    int best_gain = 0;

    //// Store the best move
    // int best_i = 0, best_j = 0, best_k = 0;

    // int i_new = -1;
    // int i_next_new = -1;
    // int j_new = -1;
    // int j_next_new = -1;
    // int k_new = -1;
    // int k_next_new = -1;

    std::vector<edge> best_edges;
    std::vector<cycle_segment> best_segments;

    // Try all possible combinations of three edges
    // TODO im not sure it this is correct
    for (int i = 0; i < n - 2; i++) {
        for (int j = i + 1; j < n - 1; j++) {  // todo revert to j = i + 1
            for (int k = j + 1; k < n; k++) {  // todo revert to k = j + 1
                // std::cout << "combination" << std::endl;

                int i_next = (i + 1) % n;
                int j_next = (j + 1) % n;
                int k_next = (k + 1) % n;

                int old_weight = weights[cycle[i]][cycle[i_next]] +
                                 weights[cycle[j]][cycle[j_next]] +
                                 weights[cycle[k]][cycle[k_next]];

                std::vector<int> selectable_vertices;

                std::vector<edge> edges;
                std::vector<cycle_segment> segments;

                auto get_segment = [i_next, j, j_next, k, k_next](int vertex) {
                    cycle_segment s;
                    s.start = vertex;

                    if (vertex == i_next) {
                        s.end = j;
                    } else if (vertex == j) {
                        s.end = i_next;
                        s.reverse = true;
                    } else if (vertex == j_next) {
                        s.end = k;
                    } else if (vertex == k) {
                        s.end = j_next;
                        s.reverse = true;
                    }
                    return s;
                };

                if (is_symmetric) {
                    selectable_vertices = {i_next, j, j_next, k};
                } else {
                    selectable_vertices = {i_next, j_next};
                }

                for (int f = 0; f < selectable_vertices.size(); f++) {
                    // int first_added_edge_start = i;
                    // int first_added_edge_end = selectable_vertices[f];
                    edge first_edge;
                    first_edge.start = i;
                    first_edge.end = selectable_vertices[f];

                    // int second_added_edge_start = -1;
                    // int second_added_edge_end = -1;

                    // int third_added_edge_start = -1;
                    // int third_added_edge_end = k_next;

                    cycle_segment first_segment = get_segment(first_edge.end);

                    edge second_edge;
                    second_edge.start = first_segment.end;

                    std::vector<int> second_selectable_vertices;
                    if (is_symmetric) {
                        if (first_edge.end == i_next || first_edge.end == j) {
                            second_selectable_vertices = {j_next, k};
                        } else {
                            second_selectable_vertices = {i_next, j};
                        }
                    } else {
                        if (first_edge.end == i_next) {
                            second_selectable_vertices = {j_next};
                        } else {
                            second_selectable_vertices = {i_next};
                        }
                    }

                    for (int g = 0; g < second_selectable_vertices.size(); g++) {
                        second_edge.end = second_selectable_vertices[g];

                        cycle_segment second_segment = get_segment(second_edge.end);

                        edge third_edge;
                        third_edge.start = second_segment.end;

                        third_edge.end = k_next;

                        cycle_segment third_segment;
                        third_segment.start = k_next;
                        third_segment.end = i;

                        // if (second_added_edge_end == i_next) {
                        //     third_added_edge_start = j;
                        // } else if (second_added_edge_end == j) {
                        //     third_added_edge_start = i_next;
                        // } else if (second_added_edge_end == j_next) {
                        //     third_added_edge_start = k;
                        // } else if (second_added_edge_end == k) {
                        //     third_added_edge_start = j_next;
                        // }

                        int new_weight = weights[cycle[first_edge.start]][cycle[first_edge.end]] +
                                         weights[cycle[second_edge.start]][cycle[second_edge.end]] +
                                         weights[cycle[third_edge.start]][cycle[third_edge.end]];

                        // std::cout << char(first_added_edge_start + int('a')) << "->"
                        //           << char(first_added_edge_end + int('a')) << " "
                        //           << char(second_added_edge_start + int('a')) << "->"
                        //           << char(second_added_edge_end + int('a')) << " "
                        //           << char(third_added_edge_start + int('a')) << "->"
                        //           << char(third_added_edge_end + int('a')) << std::endl;

                        int gain = old_weight - new_weight;
                        if (gain > best_gain) {
                            best_gain = gain;

                            best_edges = {first_edge, second_edge, third_edge};
                            best_segments = {first_segment, second_segment, third_segment};

                            std::cout << " prev " << cycle[i] << "->" << cycle[i_next] << " "
                                      << cycle[j] << "->" << cycle[j_next] << " " << cycle[k]
                                      << "->" << cycle[k_next] << std::endl;
                            std::cout << "moves " << cycle[first_edge.start] << "->"
                                      << cycle[first_edge.end] << " " << cycle[second_edge.start]
                                      << "->" << cycle[second_edge.end] << " "
                                      << cycle[third_edge.start] << "->" << cycle[third_edge.end]
                                      << std::endl;

                            improved = true;
                            std::cout << "improvement " << improved << " gain " << gain << " new "
                                      << new_weight << " old " << old_weight << std::endl;
                            std::cout << "";
                        }
                    }
                }

                // Current edges
                // Try all permutations of the three segments
                // std::vector<int> segments = {j, k, i};
                // int perm_index = 0;
                // do {
                //    // std::cout << "permutation" << std::endl;

                //              << std::endl;

                //    perm_index++;
                //} while (std::next_permutation(segments.begin(), segments.end()));
                // break;
            }
        }
    }

    // Apply the best move if found
    if (improved) {
        std::cout << "Before improving: " << edges_to_add_in_cycle(cycle, weights) << std::endl;
        print_cycle(cycle);

        edge first_edge = best_edges[0];  // i.....x
        edge second_edge = best_edges[1];
        edge third_edge = best_edges[2];  // z...k_next

        cycle_segment first_segment = best_segments[0];
        cycle_segment second_segment = best_segments[1];
        cycle_segment third_segment = best_segments[2];

        std::vector<int> new_cycle;
        std::set<int> already_added;

        new_cycle.push_back(cycle[first_edge.start]);
        new_cycle.push_back(cycle[first_edge.end]);

        already_added.insert(first_edge.start);
        already_added.insert(first_edge.end);

        if (first_segment.start != first_segment.end) {
            int first_segment_min =
                !first_segment.reverse ? first_segment.start : first_segment.end;
            int first_segment_max =
                !first_segment.reverse ? first_segment.end : first_segment.start;

            std::vector<int> first_segment_vector(cycle.begin() + first_segment_min + 1,
                                                  cycle.begin() + first_segment_max);

            if (first_segment.reverse) {
                std::reverse(first_segment_vector.begin(), first_segment_vector.end());
            }

            new_cycle.insert(new_cycle.end(), first_segment_vector.begin(),
                             first_segment_vector.end());
        }

        if (!contains(already_added, second_edge.start)) {
            new_cycle.push_back(cycle[second_edge.start]);
            already_added.insert(second_edge.start);
        } else {
            std::cout << "already added " << cycle[second_edge.start] << std::endl;
        }

        if (!contains(already_added, second_edge.end)) {
            new_cycle.push_back(cycle[second_edge.end]);
            already_added.insert(second_edge.end);
        } else {
            std::cout << "already added " << cycle[second_edge.end] << std::endl;
        }

        if (second_segment.start != second_segment.end) {
            int second_segment_min =
                !second_segment.reverse ? second_segment.start : second_segment.end;
            int second_segment_max =
                !second_segment.reverse ? second_segment.end : second_segment.start;

            std::vector<int> second_segment_vector(cycle.begin() + second_segment_min + 1,
                                                   cycle.begin() + second_segment_max);

            if (second_segment.reverse) {
                std::reverse(second_segment_vector.begin(), second_segment_vector.end());
            }

            new_cycle.insert(new_cycle.end(), second_segment_vector.begin(),
                             second_segment_vector.end());
        }

        if (!contains(already_added, third_edge.start)) {
            new_cycle.push_back(cycle[third_edge.start]);
            already_added.insert(third_edge.start);
        } else {
            std::cout << "already added " << cycle[third_edge.start] << std::endl;
        }

        if (!contains(already_added, third_edge.end)) {
            new_cycle.push_back(cycle[third_edge.end]);
            already_added.insert(third_edge.end);
        } else {
            std::cout << "already added " << cycle[third_edge.end] << std::endl;
        }

        int third_segment_min = third_segment.start;

        if (third_segment.start > 0) {
            new_cycle.insert(new_cycle.end(), cycle.begin() + third_segment.start + 1, cycle.end());
        }

        if (third_segment.end > 1) {
            new_cycle.insert(new_cycle.end(), cycle.begin() + 1, cycle.begin() + first_edge.start);
        }

        if (third_segment.start != 0 && third_segment.end != 0) {
            new_cycle.insert(new_cycle.end(), cycle.begin(), cycle.begin() + 1);
        }

        std::cout << "After improving: " << edges_to_add_in_cycle(new_cycle, weights) << std::endl;
        print_cycle(new_cycle);
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
    // return 24;
    static int call_count = 0;
    static int last_vertex = -1;
    call_count++;

    std::cout << "Start vertex selection call #" << call_count;

    if (last_vertex != -1) {
        std::cout << " (previous: " << last_vertex << ")";
    }
    std::cout << std::endl;

    last_vertex = call_count % n;
    return last_vertex;
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
    std::cout << edges_to_add_in_cycle(cycle, weighted_graph.weightMatrix) << std::endl;

    // Apply 3-Opt until no improvement is found
    bool improved;
    int iters = 0;

    // this makes it iterated three_opt which is slower
    do {
        improved = three_opt_move(cycle, weighted_graph.weightMatrix);
        iters++;
        //  std::cout << iters << std::endl;

        // std::cout << edges_to_add_in_cycle(cycle, weighted_graph.weightMatrix) << std::endl;
        // std::cout << "After improving: "
        //           << edges_to_add_in_cycle(cycle, weighted_graph.weightMatrix) << std::endl;
        // print_cycle(cycle);
    } while (improved && iters < 10 * n);

    int result = edges_to_add_in_cycle(cycle, weighted_graph.weightMatrix);
    std::cout << "length of the cycle " << cycle.size() << std::endl;
    print_cycle(cycle);

    // Count the number of edges that need to be added (edges with weight 1 in the cycle)
    return result;
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