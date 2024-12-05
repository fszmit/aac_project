#include "common.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;



vector<Graph> parseGraphs(const string& filename) {
    ifstream inputFile(filename);
    vector<Graph> graphs;

    if (!inputFile.is_open()) {
        cerr << "Error: Could not open file." << endl;
        return graphs;
    }

    string numberOfGraphs;
    getline(inputFile, numberOfGraphs);

    for (int i = 0; i < stoi(numberOfGraphs); ++i) {
        Graph graph;
        string line;

        // number of vertices
        getline(inputFile, line);
        while (line == "")
            getline(inputFile, line);
        graph.vertices = stoi(line);

        // adjacency matrix
        for (int j = 0; j < graph.vertices; ++j) {
            getline(inputFile, line);
            istringstream rowStream(line);   //Input stream class to operate on strings. Objects of this class use a string buffer that contains a sequence of characters. This sequence of characters can be accessed directly as a string object
            vector<int> row;
            int value;

            while (rowStream >> value) {
                row.push_back(value);
            }

            graph.adjacencyMatrix.push_back(row);
        }

        // additional information (if any)
        if (getline(inputFile, line) && !line.empty()) {
            graph.additionalInfo = line;
        }
        graphs.push_back(graph);
    }

    inputFile.close();
    return graphs;
}

void printGraph(const Graph& graph) {
    cout << "Number of vertices: " << graph.vertices << endl;
    cout << "Adjacency Matrix:" << endl;

    for (const auto& row : graph.adjacencyMatrix) {
        for (int value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    if (!graph.additionalInfo.empty()) {
        cout << "Additional Info: " << graph.additionalInfo << endl;
    }
}

Graph InitGraph(int size) {
    Graph g;
    g.vertices = size;
    g.adjacencyMatrix = std::vector<std::vector<int> >(size, std::vector<int>(size));
    return g;
}

Graph CloneGraph(const Graph* g) {
    Graph res = InitGraph(g->vertices);
    for (int i = 0; i < g->vertices; i++)
    {
        for (int j = 0; j < g->vertices; j++)
        {
            res.adjacencyMatrix[i][j] = g->adjacencyMatrix[i][j];
        }
    }
    res.additionalInfo = g->additionalInfo;
    return res;
}

std::vector<std::vector<int>> FlattenMatrix(std::vector<std::vector<int>> mat) {
    std::vector<std::vector<int>> res(mat.size(), std::vector<int>(mat.size(), 0));
    for (ul i = 0; i < mat.size(); i++)
    {
        for (ul j = 0; j < mat.size(); j++)
        {
            if (j == i) {
                res[i][j] = 0;
                continue;
            }
            res[i][j] = min(mat[i][j], 1);
        }
    }
    return res;
}

std::vector<std::vector<int>> comb(int n, int k) {
    std::vector<std::vector<int>> res;

    std::string bitmask(k, 1);
    bitmask.resize(n, 0);

    do
    {
        std::vector<int> sub;
        for (int i = 0; i < n; i++)
        {
            if (bitmask[i])
            {
                sub.push_back(i);
            }
        }
        res.push_back(sub);

    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    return res;
}

unsigned long long binom(int n, int k) {
    if (n < 0 || k < 0 || k > n) {
        return 0;
    }
    int p = 1;
    for (int i = 1; i <= std::min(k, n - k); i++)
    {
        p *= n + 1 - i;
        p /= i;
    }
    return p;
}

bool checkIfDirected(const Graph& graph) {
    for (int i = 0; i < graph.vertices; i++) {
        for (int j = i; j < graph.vertices; j++) {
            if (graph.adjacencyMatrix.at(i).at(j) != graph.adjacencyMatrix.at(j).at(i))
                return true;
        }
    }
    return false;
}