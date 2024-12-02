#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

struct Graph {
    int vertices;
    vector<vector<int>> adjacencyMatrix;
    string additionalInfo;
};

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
        while(line == "")
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

void countHamiltonianCyclesUtil(const Graph& graph, int pos, vector<bool>& visited, int count, int& cycles) {
    int n = graph.vertices;
    if (count == n) {
        if (graph.adjacencyMatrix[pos][0]) {
            cycles++;
        }
        return;
    }
    for (int v = 0; v < n; ++v) {
        if (!visited[v] && graph.adjacencyMatrix[pos][v]) {
            visited[v] = true;
            countHamiltonianCyclesUtil(graph, v, visited, count + 1, cycles);
            visited[v] = false;
        }
    }
}

int countHamiltonianCycles(const Graph& graph) {
    int n = graph.vertices;
    int cycles = 0;
    vector<bool> visited(n, false);
    visited[0] = true;
    countHamiltonianCyclesUtil(graph, 0, visited, 1, cycles);
    return cycles;
}

int minimalExtension(Graph graph) {
    int n = graph.vertices;
    vector<pair<int, int>> missingEdges;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && graph.adjacencyMatrix[i][j] == 0) {
                missingEdges.push_back({i, j});
            }
        }
    }
    int cyclesFirst =  countHamiltonianCycles(graph);
    if (cyclesFirst > 0) {
        return 0;
    }
    int minEdgesToAdd = -1;
    for (int k = 1; k <= missingEdges.size(); ++k) {
        vector<bool> bitmask(k, true);
        bitmask.resize(missingEdges.size(), false);
        do {
            for (size_t i = 0; i < bitmask.size(); ++i) {
                if (bitmask[i]) {
                    auto edge = missingEdges[i];
                    graph.adjacencyMatrix[edge.first][edge.second] = 1;
                }
            }
            int cycles = countHamiltonianCycles(graph);
            if (cycles > 0) {
                minEdgesToAdd = k;
                return minEdgesToAdd;
            }
            for (size_t i = 0; i < bitmask.size(); ++i) {
                if (bitmask[i]) {
                    auto edge = missingEdges[i];
                    graph.adjacencyMatrix[edge.first][edge.second] = 0;
                }
            }
        } while (prev_permutation(bitmask.begin(), bitmask.end()));
    }
    return minEdgesToAdd;
}

void countHamiltonianCyclesUtilUnDirected(const Graph& graph, int pos, vector<bool>& visited, int count, int& cycles) {
    int n = graph.vertices;
    if (count == n) {
        if (graph.adjacencyMatrix[pos][0]) {
            cycles++;
        }
        return;
    }
    for (int v = 0; v < n; v++) {
        if (graph.adjacencyMatrix[pos][v] && !visited[v]) {
            visited[v] = true;
            countHamiltonianCyclesUtil(graph, v, visited, count + 1, cycles);
            visited[v] = false;
        }
    }
}

int countHamiltonianCyclesUnDirected(const Graph& graph) {
    int n = graph.vertices;
    int cycles = 0;
    vector<bool> visited(n, false);
    visited[0] = true;
    countHamiltonianCyclesUtilUnDirected(graph, 0, visited, 1, cycles);
    return cycles / 2;
}

int minimalExtensionUnDirected(Graph graph) {
    int n = graph.vertices;
    vector<pair<int, int>> missingEdges;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (graph.adjacencyMatrix[i][j] == 0) {
                missingEdges.push_back({i, j});
            }
        }
    }

    int cyclesFirst = countHamiltonianCyclesUnDirected(graph);
    if (cyclesFirst > 0) {
        return 0;
    }

    int minEdgesToAdd = -1;
    int m = missingEdges.size();
    for (int k = 1; k <= m; ++k) {
        vector<bool> bitmask(k, true);
        bitmask.resize(m, false);

        do {
            for (int i = 0; i < m; ++i) {
                if (bitmask[i]) {
                    auto edge = missingEdges[i];
                    graph.adjacencyMatrix[edge.first][edge.second] = 1;
                    graph.adjacencyMatrix[edge.second][edge.first] = 1;
                }
            }

            int cycles = countHamiltonianCyclesUnDirected(graph);
            if (cycles > 0) {
                minEdgesToAdd = k;
                return minEdgesToAdd;
            }

            for (int i = 0; i < m; ++i) {
                if (bitmask[i]) {
                    auto edge = missingEdges[i];
                    graph.adjacencyMatrix[edge.first][edge.second] = 0;
                    graph.adjacencyMatrix[edge.second][edge.first] = 0;
                }
            }
        } while (prev_permutation(bitmask.begin(), bitmask.end()));
    }
    return minEdgesToAdd;
}

int main() {
    string filename;
    cout << "Data file path: ";
    cin >> filename;
    
    vector<Graph> graphs = parseGraphs(filename);
    cout << "Parsed " << graphs.size() << " graphs." << endl;

    for (size_t i = 0; i < graphs.size(); ++i) {
        cout << "\nGraph " << i + 1 << ":" << endl;
        printGraph(graphs[i]);
    }

    int minimalExtensionCount = minimalExtensionUnDirected(graphs[0]);
    cout << "Minimal number of edges to add to make the graph Hamiltonian: " << minimalExtensionCount << endl;
    int minimalExtensionCount2 = minimalExtension(graphs[1]);
    cout << "Minimal number of edges to add to make the graph Hamiltonian: " << minimalExtensionCount2 << endl;

    return 0;
}
