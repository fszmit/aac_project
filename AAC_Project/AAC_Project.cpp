#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "hamiltonian_completion/approximation.h"

using namespace std;

struct Graph {
    int vertices;
    vector<vector<int>> adjacencyMatrix;
    string additionalInfo;
};

vector<Graph> parseGraphs(const string &filename) {
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
        while (line == "") getline(inputFile, line);
        graph.vertices = stoi(line);

        // adjacency matrix
        for (int j = 0; j < graph.vertices; ++j) {
            getline(inputFile, line);
            // Input stream class to operate on strings. Objects of this
            // class use a string buffer that contains a sequence of
            // characters. This sequence of characters can be accessed
            // directly as a string object
            istringstream rowStream(line);
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

void printGraph(const Graph &graph) {
    cout << "Number of vertices: " << graph.vertices << endl;
    cout << "Adjacency Matrix:" << endl;

    for (const auto &row : graph.adjacencyMatrix) {
        for (int value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    if (!graph.additionalInfo.empty()) {
        cout << "Additional Info: " << graph.additionalInfo << endl;
    }
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

        print_weighted_graph(transform_to_complete_weighted_graph(graphs[i].adjacencyMatrix));
    }

    return 0;
}
