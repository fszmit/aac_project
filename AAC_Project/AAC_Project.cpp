#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

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

Graph InitGraph(int size) {
    Graph g;
    g.vertices = size;
    g.adjacencyMatrix = std::vector<std::vector<int> >(size, std::vector<int>(size));
    return g;
}

Graph SubGraph(Graph * g, std::vector<int> mask, int * neighbours = NULL) {

    int sub_size = g->vertices - mask.size();
    Graph s = InitGraph(sub_size);

    std::vector<int> n = std::vector<int>(mask.size(), 0);

    int si = 0;
    for (int i = 0; i < g->vertices; i++) // rows
    {
        int sj = 0;

        if (i == mask[i - si]) { // masked row
            for (int j = 0; j < g->vertices; j++)
            {
                if (j != mask[j - sj]) { // masked row, unmasked column
                    n[i - si] += g->adjacencyMatrix[i][j] >= 1 ? 1 : 0;
                    sj ++;
                }
                //masked row, masked column
                //nop;
            }
            continue;
        }

        //unmasked row
        for (int j = 0; j < g->vertices; j++) // cols
        {
            if (j == mask[j - sj]){ // unmasked row, masked column
                if (g->adjacencyMatrix[i][j] >= 1) {
                    n[j - sj]++;
                }
                continue;
            }
            //unmasked row, unmasked column
            s.adjacencyMatrix[si][sj] = g->adjacencyMatrix[i][j];
            sj ++;
        }
        si ++;
    }
    if (neighbours != NULL) {
        *neighbours = 0;
        for (int i = 0; i < (int)n.size(); i++)
        {
            *neighbours += n[i] > 0 ? 1 : 0;
        }
    }
    return s;
}

std::vector<std::vector<int>> MatMul(std::vector<std::vector<int>> mat1, std::vector<std::vector<int>> mat2) {
    // assumes the sizes are equal
    std::vector<std::vector<int>> res = std::vector<std::vector<int> >(mat1.size(), std::vector<int>(mat1.size(), 0));

    for (int i = 0; i < mat1.size(); i++)
    {
        for (int j = 0; j < mat1.size(); j++)
        {
            for (int k = 0; k < mat1.size(); k++)
            {
                res[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return res;
}

int binom(int n, int k) {
    if (n < 0 || k < 0 || k > n){
        return 0;
    }
    int p = 1;
    for (int i = 1; i <= std::min(k, n-k); i++)
    {
        p *= n + 1 - i; 
        p /= i;
    }
    return p;
}

int MaxCycle(Graph * g, int * no_cycles = NULL) {
    // based on https://arxiv.org/pdf/1612.05531
    // NOTE:: this returnx the number of **directed** cycles. so a singular cycle in an undirected graph will count for two


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
        cout << "\nSubgraph:\n";
        int n = 0;
        Graph sub = SubGraph(&graphs[0], {3}, &n);
        printGraph(sub);
        cout << "neighbours: " << n << '\n';
    }




    return 0;
}
