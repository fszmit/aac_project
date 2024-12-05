#ifndef GED_HPP
#define GEDS_HPP

#include "common.hpp"
#include <unordered_map>

struct EditPath {
    unordered_map<int, int> mapping; // Mapping of vertices from g1 to g2
    double costSoFar;               // g(p)
    double heuristicCost;           // lb(p)

    bool operator>(const EditPath& other) const {
        return (costSoFar + heuristicCost) > (other.costSoFar + other.heuristicCost);
    }
};

double nodeSubstitutionCost(int u1, int u2);
int GraphSize(const Graph& graph);
double edgeSubstitutionCost(int e1, int e2);
double insertionCost();
double deletionCost();
double computeHeuristic(const Graph& g1, const Graph& g2, const unordered_map<int, int>& mapping);
double astarGED(const Graph& g1, const Graph& g2);


int CalculateGraphEditDistance(Graph& A, Graph& B);

#endif //GED_HPP