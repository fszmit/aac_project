#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>

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

        auto total_start_time = std::chrono::high_resolution_clock::now();
        double total_completion = 0;
        vector<int> all_completions;
        vector<long long> all_durations;
        
        for (int j = 0; j < 50; j++) {
            auto iteration_start = std::chrono::high_resolution_clock::now();
            
            int c = hamiltonian_completion_approximation(graphs[i].adjacencyMatrix);
            
            auto iteration_end = std::chrono::high_resolution_clock::now();
            auto iteration_duration = std::chrono::duration_cast<std::chrono::milliseconds>(iteration_end - iteration_start);
            
            total_completion += c;
            all_completions.push_back(c);
            all_durations.push_back(iteration_duration.count());
            
            cout << "Iteration " << j + 1 << ": completion = " << c 
                 << ", time = " << iteration_duration.count() << " ms" << endl;
        }

        auto total_end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end_time - total_start_time);
        
        // Calculate min and max values
        int min_completion = *std::min_element(all_completions.begin(), all_completions.end());
        int max_completion = *std::max_element(all_completions.begin(), all_completions.end());
        long long min_time = *std::min_element(all_durations.begin(), all_durations.end());
        long long max_time = *std::max_element(all_durations.begin(), all_durations.end());
        
        cout << "\nSummary:" << endl;
        cout << "Completion size - Min: " << min_completion << ", Max: " << max_completion 
             << ", Avg: " << (total_completion / 50.0) << endl;
        cout << "Time (ms) - Min: " << min_time << ", Max: " << max_time 
             << ", Avg: " << (total_duration.count() / 50.0) << endl;
        cout << "Total execution time: " << total_duration.count() << " ms" << endl;
    }

    return 0;
}
