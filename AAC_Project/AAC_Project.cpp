#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <queue>
#include <algorithm>

using namespace std;

typedef unsigned long ul;
typedef unsigned long long ull;

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

Graph CloneGraph(const Graph * g) {
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

Graph SubGraph(Graph * g, std::vector<int> mask, int * neighbours = NULL) {
    if (mask.size() == 0) {
        return CloneGraph(g);
    }
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

    for (ul i = 0; i < mat1.size(); i++)
    {
        for (ul j = 0; j < mat1.size(); j++)
        {
            for (ul k = 0; k < mat1.size(); k++)
            {
                res[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return res;
}

std::vector<std::vector<int>> MatPow(std::vector<std::vector<int>> mat, int exp){
    std::vector<std::vector<int>> res;
    bool set = false;
    std::vector<std::vector<int>> curr = mat;
    int i = 1;
    while (exp != 0)
    {
        if (exp % 2 == 1){
            //cout << "at " << i;
            if (!set){
                res = curr;
                set = true;
            }
            else {
                res = MatMul(res, curr);
            }
        }
        exp /= 2;
        i*=2;
        curr = MatMul(curr, curr);
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

bool isWeaklyConnected(const std::vector<std::vector<int>> * mat) {
    std::vector<bool> conn (mat->size(), false);
    conn[0] = true;
    std::queue<int> q;
    int connected = 1;
    q.push(0);

    while (!q.empty())
    {
        int c_vert = q.front(); q.pop();
        for (int i = 0; i < (int)mat->size(); i++)
        {
            if (c_vert == i || conn[i] == true) {continue;}
            if ((*mat)[c_vert][i] >= 1 ||(* mat)[i][c_vert] >= 1) {
                q.push(i);
                conn[i] = true;
                connected++;
                if (connected == (int)mat->size())
                {
                    return true;
                }
            }
        }
    }
    return false;

}

std::vector<std::vector<int>> FlattenMatrix(std::vector<std::vector<int>> mat) {
    std::vector<std::vector<int>> res (mat.size(), std::vector<int>(mat.size(), 0));
    for (ul i = 0; i < mat.size(); i++)
    {
        for (ul j = 0; j < mat.size(); j++)
        {
            res[i][j] = max(mat[i][j], 1);
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
        for (int i = 0; i < n; i++)
        {
            std::vector<int> sub;
            if (bitmask[i])
            {
                sub.push_back(i);
            }
            res.push_back(sub);
        }
        
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    return res;
}

int MatTrace(std::vector<std::vector<int>> mat) {
    int res = 0;
    for (ul i = 0; i < mat.size(); i++)
    {
        res += mat[i][i];
    }
    return res;
}

int MaxCycle(Graph * g, int * no_cycles = NULL) {
    // based on https://arxiv.org/pdf/1612.05531
    // NOTE:: this returns the number of **directed** cycles. so a singular cycle in an undirected graph will count for two

    Graph flat;
    flat.vertices = g->vertices;

    flat.adjacencyMatrix = FlattenMatrix(g->adjacencyMatrix);

    flat.additionalInfo = g->additionalInfo;

    for (int i = flat.vertices; i >= 1; i--)
    {
        
        int gamma = 0;
        for (int j = flat.vertices; j >=1 ; j--)
        {
            std::vector<std::vector<int>> combs = comb(flat.vertices, j);
            for (int k = 0; k < (int)combs.size(); k++)
            {
                int sub_gamma = 0;
                int neigh;
                Graph sub = SubGraph(&flat, combs[k], &neigh);
                if (! isWeaklyConnected(&sub.adjacencyMatrix)) {
                    continue;
                }
                sub_gamma = binom(neigh, i - sub.vertices) * ( sub.vertices % 2 == 1 ? -1 : 1 ) * MatTrace( MatPow(sub.adjacencyMatrix, i) );
                gamma += sub_gamma;
            }
            
        }
        gamma /= i;
        gamma *= (i % 2 == 1 ? -1 : 1);
        if (gamma != 0 )
        {
            if (no_cycles != NULL){
                *no_cycles = gamma;
            }
            return i;
        }   
    }   
    *no_cycles = 0;
    return 0;
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
        Graph sub = SubGraph(&graphs[0], {2}, &n);
        printGraph(sub);
        cout << "neighbours: " << n << '\n';

        cout << "power of 5:\n";
        Graph p;
        p.vertices = graphs[i].vertices;
        p.adjacencyMatrix = MatPow(graphs[i].adjacencyMatrix, 5);
        printGraph(p);

        int no_cycles = 0;
        cout << "Max Cycle size: " << MaxCycle(&graphs[i], &no_cycles) << '\n';
        cout << "No Max Cycles: " << no_cycles << '\n';

        cout << "\n\n";
    }


    return 0;
}
