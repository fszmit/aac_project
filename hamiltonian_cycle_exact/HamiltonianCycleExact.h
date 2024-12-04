#pragma once
#include <vector>
#include <string>

using namespace std;

struct Graph {
    int vertices;
    vector<vector<int>> adjacencyMatrix;
    string additionalInfo;
};

void countHamiltonianCyclesUtil(const Graph& graph, int pos, vector<bool>& visited, int count, int& cycles);
int countHamiltonianCycles(const Graph& graph);
int minimalExtension(Graph graph);
void countHamiltonianCyclesUtilUnDirected(const Graph& graph, int pos, vector<bool>& visited, int count, int& cycles);
int countHamiltonianCyclesUnDirected(const Graph& graph);
int minimalExtensionUnDirected(Graph graph);