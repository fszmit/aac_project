#include "GED.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <set>
using namespace std;


int GraphSize(const Graph& graph) {
    int EdgeSum = 0;
    for (int i = 0; i < graph.vertices; i++) {
        for (int j = i; j < graph.vertices; j++) {
            EdgeSum += graph.adjacencyMatrix.at(i).at(j);
        }
    }

    return EdgeSum;
}



// Cost functions
double nodeSubstitutionCost(int u1, int u2) {
    return (u1 == u2) ? 0 : 1; // Example: labels same = 0 cost, different = 1 cost
}

double edgeSubstitutionCost(int e1, int e2) {
    return abs(e1 - e2); // Example: difference in edge multiplicity
}

double insertionCost() {
    return 1.0; // Example: cost to insert a vertex or edge
}

double deletionCost() {
    return 1.0; // Example: cost to delete a vertex or edge
}

// Heuristic function lb(p)
double computeHeuristic(const Graph& g1, const Graph& g2, const unordered_map<int, int>& mapping) {
    // Use a simple heuristic based on the number of unmapped vertices and edges
    set<int> mappedVertices;
    for (const auto& [u1, u2] : mapping) {
        mappedVertices.insert(u2);
    }

    int unmappedG1 = g1.vertices - mapping.size();
    int unmappedG2 = g2.vertices - mappedVertices.size();

    // Assume costs are proportional to the number of remaining vertices/edges
    return unmappedG1 * insertionCost() + unmappedG2 * deletionCost();
}

// A*GED Algorithm
double astarGED(const Graph& g1, const Graph& g2) {
    priority_queue<EditPath, vector<EditPath>, greater<>> openSet;
    double minCost = numeric_limits<double>::infinity();

    // Initialize the search with the first node
    openSet.push({ {}, 0.0, computeHeuristic(g1, g2, {}) });

    while (!openSet.empty()) {
        EditPath current = openSet.top();
        openSet.pop();

        // If all nodes in g1 are mapped, calculate the cost and terminate if it's the minimum
        if (current.mapping.size() == g1.vertices) {
            double finalCost = current.costSoFar + computeHeuristic(g1, g2, current.mapping);
            minCost = min(minCost, finalCost);
            continue;
        }

        // Expand the current path
        int nextNode = current.mapping.size(); // Next unmapped node in g1

        for (int v2 = 0; v2 < g2.vertices; ++v2) {
            if (current.mapping.find(nextNode) != current.mapping.end()) continue;

            // Create a new mapping and calculate costs
            EditPath newPath = current;
            newPath.mapping[nextNode] = v2;

            // Compute edit costs
            double nodeCost = nodeSubstitutionCost(nextNode, v2);
            double edgeCost = 0.0;

            for (const auto& [u1, u2] : current.mapping) {
                if (u2 != -1) { // Only process valid mappings
                    edgeCost += edgeSubstitutionCost(
                        g1.adjacencyMatrix[nextNode][u1],
                        g2.adjacencyMatrix[v2][u2]
                    );
                }
                else {
                    // Handle deletion case: Add the cost of deleting the edge
                    edgeCost += deletionCost();
                }
            }

            newPath.costSoFar += nodeCost + edgeCost;
            newPath.heuristicCost = computeHeuristic(g1, g2, newPath.mapping);

            // Push the new path to the open set
            openSet.push(newPath);
        }

        // Consider deleting the node
        EditPath deletionPath = current;
        deletionPath.mapping[nextNode] = -1; // Use -1 to denote deletion
        deletionPath.costSoFar += deletionCost();
        deletionPath.heuristicCost = computeHeuristic(g1, g2, deletionPath.mapping);
        openSet.push(deletionPath);
    }

    return minCost;
}

int CalculateGraphEditDistance(Graph& A, Graph& B) {
    int distance = 0;
    distance += abs(B.vertices - A.vertices);
    int bigger = B.vertices;
    int smaller = A.vertices;
    if (B.vertices < A.vertices) {
        bigger = A.vertices;
        smaller = B.vertices;
    }
    if (A.additionalInfo == "directed") {
        for (int i = 0; i < bigger; i++) {
            for (int j = 0; j < bigger; j++) {
                if (i >= A.vertices || j >= A.vertices) {
                    distance += B.adjacencyMatrix.at(i).at(j);
                }
                else if (i >= B.vertices || j >= B.vertices) {
                    distance += A.adjacencyMatrix.at(i).at(j);
                }
                else
                    distance += abs(A.adjacencyMatrix.at(i).at(j) - B.adjacencyMatrix.at(i).at(j));
            }
        }
    }
    else {
        for (int i = 0; i < bigger; i++) {
            for (int j = i; j < bigger; j++) {
                if (i >= A.vertices || j >= A.vertices) {
                    distance += B.adjacencyMatrix.at(i).at(j);
                }
                else if (i >= B.vertices || j >= B.vertices) {
                    distance += A.adjacencyMatrix.at(i).at(j);
                }
                else
                    distance += abs(A.adjacencyMatrix.at(i).at(j) - B.adjacencyMatrix.at(i).at(j));
            }
        }
    }
    return distance;
}
