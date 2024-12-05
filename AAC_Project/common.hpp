#ifndef COMMON_HPP
#define COMMON_HPP

#include <vector>
#include <string>

using namespace std;

typedef unsigned long ul;
typedef unsigned long long ull;

struct Graph {
    int vertices;
    vector<vector<int>> adjacencyMatrix;
    string additionalInfo;
};

vector<Graph> parseGraphs(const string& filename);
void printGraph(const Graph& graph);
Graph InitGraph(int size);
Graph CloneGraph(const Graph* g);
std::vector<std::vector<int>> FlattenMatrix(std::vector<std::vector<int>> mat);
std::vector<std::vector<int>> comb(int n, int k);
unsigned long long binom(int n, int k);
bool checkIfDirected(const Graph& graph);
#endif //COMMON_HPP