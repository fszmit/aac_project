#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <chrono>

using namespace std;

struct Graph
{
    int vertices;
    vector<vector<int>> adjacencyMatrix;
    string additionalInfo;
};

vector<Graph> parseGraphs(const string& filename)
{
    ifstream inputFile(filename);
    vector<Graph> graphs;

    if (!inputFile.is_open())
    {
        cerr << "Error: Could not open file." << endl;
        return graphs;
    }

    string numberOfGraphs;
    getline(inputFile, numberOfGraphs);

    for (int i = 0; i < stoi(numberOfGraphs); ++i)
    {
        Graph graph;
        string line;

        // number of vertices
        getline(inputFile, line);
        while (line == "")
            getline(inputFile, line);
        graph.vertices = stoi(line);

        // adjacency matrix
        for (int j = 0; j < graph.vertices; ++j)
        {
            getline(inputFile, line);
            istringstream rowStream(line);
            //Input stream class to operate on strings. Objects of this class use a string buffer that contains a sequence of characters. This sequence of characters can be accessed directly as a string object
            vector<int> row;
            int value;

            while (rowStream >> value)
            {
                row.push_back(value);
            }

            graph.adjacencyMatrix.push_back(row);
        }

        // additional information (if any)
        if (getline(inputFile, line) && !line.empty())
        {
            graph.additionalInfo = line;
        }
        graphs.push_back(graph);
    }

    inputFile.close();
    return graphs;
}

void printGraph(const Graph& graph)
{
    cout << "Number of vertices: " << graph.vertices << endl;
    cout << "Adjacency Matrix:" << endl;

    for (const auto& row : graph.adjacencyMatrix)
    {
        for (int value : row)
        {
            cout << value << " ";
        }
        cout << endl;
    }

    if (!graph.additionalInfo.empty())
    {
        cout << "Additional Info: " << graph.additionalInfo << endl;
    }
}

bool hasHamiltonianCycleUtil(const Graph& graph, int pos, vector<bool>& visited, int count, bool isDirected)
{
    int n = graph.vertices;
    if (count == n)
    {
        if (isDirected)
        {
            return graph.adjacencyMatrix[pos][0] > 0;
        }
        return (graph.adjacencyMatrix[pos][0] > 0 || graph.adjacencyMatrix[0][pos] > 0);
    }
    for (int v = 0; v < n; ++v)
    {
        if (graph.adjacencyMatrix[pos][v] && !visited[v])
        {
            visited[v] = true;
            if (hasHamiltonianCycleUtil(graph, v, visited, count + 1, isDirected))
            {
                return true;
            }
            visited[v] = false;
        }
    }
    return false;
}

bool hasHamiltonianCycle(const Graph& graph, bool isDirected)
{
    int n = graph.vertices;
    vector<bool> visited(n, false);
    visited[0] = true;
    return hasHamiltonianCycleUtil(graph, 0, visited, 1, isDirected);
}

int minimalExtension(Graph graph, bool isDirected)
{
    int n = graph.vertices;
    vector<pair<int, int>> missingEdges;

    for (int i = 0; i < n; ++i)
    {
        int jStart = isDirected ? 0 : i + 1;
        for (int j = jStart; j < n; ++j)
        {
            if (i != j && graph.adjacencyMatrix[i][j] == 0)
            {
                missingEdges.push_back({i, j});
            }
        }
    }

    if (hasHamiltonianCycle(graph, isDirected))
    {
        return 0;
    }

    int m = missingEdges.size();
    for (int k = 1; k <= m; ++k)
    {
        vector<bool> bitmask(k, true);
        bitmask.resize(m, false);

        do
        {
            for (int i = 0; i < m; ++i)
            {
                if (bitmask[i])
                {
                    auto edge = missingEdges[i];
                    graph.adjacencyMatrix[edge.first][edge.second] = 1;
                    if (!isDirected)
                    {
                        graph.adjacencyMatrix[edge.second][edge.first] = 1;
                    }
                }
            }

            if (hasHamiltonianCycle(graph, isDirected))
            {
                return k;
            }

            for (int i = 0; i < m; ++i)
            {
                if (bitmask[i])
                {
                    auto edge = missingEdges[i];
                    graph.adjacencyMatrix[edge.first][edge.second] = 0;
                    if (!isDirected)
                    {
                        graph.adjacencyMatrix[edge.second][edge.first] = 0;
                    }
                }
            }
        }
        while (prev_permutation(bitmask.begin(), bitmask.end()));
    }
    return -1;
}


bool checkIfDirected(const Graph& graph)
{
    for (int i = 0; i < graph.vertices; i++)
    {
        for (int j = i; j < graph.vertices; j++)
        {
            if (graph.adjacencyMatrix.at(i).at(j) != graph.adjacencyMatrix.at(j).at(i))
                return true;
        }
    }
    return false;
}

void generateRandomGraph(Graph& graph, int vertices, int density)
{
    graph.vertices = vertices;
    graph.adjacencyMatrix = vector<vector<int>>(vertices, vector<int>(vertices, 0));

    srand(time(0));

    for (int i = 0; i < vertices; ++i)
    {
        for (int j = i + 1; j < vertices; ++j)
        {
            if ((rand() % 100) < density)
            {
                graph.adjacencyMatrix[i][j] = 1;
                graph.adjacencyMatrix[j][i] = 1;
            }
        }
    }
}

void generateSparseGraph(Graph& graph, int vertices)
{
    graph.vertices = vertices;
    graph.adjacencyMatrix = vector<vector<int>>(vertices, vector<int>(vertices, 0));

    for (int i = 0; i < vertices - 1; ++i)
    {
        graph.adjacencyMatrix[i][i + 1] = 1;
        graph.adjacencyMatrix[i + 1][i] = 1;
    }
}

void generateDenseGraph(Graph& graph, int vertices)
{
    graph.vertices = vertices;
    graph.adjacencyMatrix = vector<vector<int>>(vertices, vector<int>(vertices, 1));

    for (int i = 0; i < vertices; ++i)
    {
        graph.adjacencyMatrix[i][i] = 0;
    }
}

void generateAlmostCompleteGraph(Graph& graph, int vertices)
{
    graph.vertices = vertices;
    graph.adjacencyMatrix = vector<vector<int>>(vertices, vector<int>(vertices, 0));
    graph.adjacencyMatrix[0][1] = 1;
}

void runTests()
{
    ofstream results("test_results.txt");
    if (!results.is_open())
    {
        cerr << "Error: Could not open results file." << endl;
        return;
    }

    for (int n = 3; n <= 7; ++n)
    {
        results << "Testing for n = " << n << " vertices:\n";
        vector<Graph> testGraphs(4);

        // Generate test cases
        generateSparseGraph(testGraphs[0], n);
        generateDenseGraph(testGraphs[1], n);
        generateRandomGraph(testGraphs[2], n, 50); // 50% density
        generateAlmostCompleteGraph(testGraphs[3], n);

        //string graphTypes[] = {"Sparse", "Dense", "Random"};
        string graphTypes[] = {"Sparse", "Dense", "Random", "Sparse Random"};

        for (int i = 0; i < 4; ++i)
        {
            auto start = chrono::high_resolution_clock::now();
            int result = minimalExtension(testGraphs[i], checkIfDirected(testGraphs[i]));
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = end - start;

            results << graphTypes[i] << " graph: Minimal extensions = " << result
                << ", Time taken = " << elapsed.count() << " seconds.\n";
        }

        results << "\n";
    }

    results.close();
}

int main()
{
    runTests();
    cout << "Tests completed. Results saved to 'test_results.txt'." << endl;
    return 0;
    // string filename;
    // cout << "Data file path: ";
    // cin >> filename;
    //
    // vector<Graph> graphs = parseGraphs(filename);
    // cout << "Parsed " << graphs.size() << " graphs." << endl;
    //
    // for (size_t i = 0; i < graphs.size(); ++i)
    // {
    //     cout << "\nGraph " << i + 1 << ":" << endl;
    //     printGraph(graphs[i]);
    // }
    //
    // for (size_t i = 0; i < graphs.size(); ++i)
    // {
    //     cout << minimalExtension(graphs[i], checkIfDirected(graphs[i])) << endl;
    // }
    // return 0;
}
