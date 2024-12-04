#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cassert>

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

bool hasHamiltonianCycleUtil(const Graph& graph, int pos, vector<bool>& visited, int count)
{
    int n = graph.vertices;
    if (count == n)
    {
        return graph.adjacencyMatrix[pos][0] == 1;
    }
    for (int v = 0; v < n; ++v)
    {
        if (graph.adjacencyMatrix[pos][v] && !visited[v])
        {
            visited[v] = true;
            if (hasHamiltonianCycleUtil(graph, v, visited, count + 1))
            {
                return true;
            }
            visited[v] = false;
        }
    }
    return false;
}

bool hasHamiltonianCycle(const Graph& graph)
{
    int n = graph.vertices;
    vector<bool> visited(n, false);
    visited[0] = true;
    return hasHamiltonianCycleUtil(graph, 0, visited, 1);
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

    if (hasHamiltonianCycle(graph))
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
                    if (isDirected)
                    {
                        graph.adjacencyMatrix[edge.second][edge.first] = 1;
                    }
                }
            }

            if (hasHamiltonianCycle(graph))
            {
                return k;
            }

            for (int i = 0; i < m; ++i)
            {
                if (bitmask[i])
                {
                    auto edge = missingEdges[i];
                    graph.adjacencyMatrix[edge.first][edge.second] = 0;
                    if (isDirected)
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

int main()
{
    string filename;
    cout << "Data file path: ";
    cin >> filename;

    vector<Graph> graphs = parseGraphs(filename);
    cout << "Parsed " << graphs.size() << " graphs." << endl;

    for (size_t i = 0; i < graphs.size(); ++i)
    {
        cout << "\nGraph " << i + 1 << ":" << endl;
        printGraph(graphs[i]);
    }

    for (size_t i = 0; i < graphs.size(); ++i)
    {
        if (hasHamiltonianCycle(graphs[i]))
        {
            cout << 0 << endl;
        }
        else
        {
            cout << minimalExtension(graphs[i], checkIfDirected(graphs[i])) << endl;
        }
    }
    return 0;
}
