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
    if (mask.size() == g->vertices) {
        if (neighbours != NULL) {
            *neighbours = 0;
        }
        Graph r;
        r.vertices = 0;
        return r;
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
    if (mat == NULL || mat->size() == 0) {return false;}
    if (mat->size() == 1) {return true;}
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
        for (int j = flat.vertices; j >=0 ; j--)
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
int GraphSize(const Graph& graph) {
    int EdgeSum = 0;
    if (graph.additionalInfo.compare("directed") != string::npos) {
        for (int i = 0; i < graph.vertices; i++) {
            for (int j = 0; j < graph.vertices; j++) {
                EdgeSum += graph.adjacencyMatrix.at(i).at(j);
            }
        }
    }
    else {
        for (int i = 0; i < graph.vertices; i++) {
            for (int j = i; j < graph.vertices; j++) {
                EdgeSum += graph.adjacencyMatrix.at(i).at(j);
            }
        }
    }
    return EdgeSum;
}
// assuming both graphs are direct/undirect
int CalculateGraphEditDistance(Graph& A, Graph& B) {
    int distance = 0;
    distance += abs(B.vertices - A.vertices);
    int bigger = B.vertices;
    int smaller = A.vertices;
    if (B.vertices < A.vertices) {
        bigger = A.vertices;
        smaller = B.vertices;
    }
    if (A.additionalInfo == "directed") {
        for (int i = 0; i < bigger; i++) {
            for (int j = 0; j < bigger; j++) {
                if (i >= A.vertices || j >= A.vertices) {
                    distance += B.adjacencyMatrix.at(i).at(j);
                }
                else if(i >= B.vertices || j >= B.vertices) {
                    distance += A.adjacencyMatrix.at(i).at(j);
                }
                else
                    distance += abs(A.adjacencyMatrix.at(i).at(j) - B.adjacencyMatrix.at(i).at(j));
            }
        }
    }
    else {
        for (int i = 0; i < bigger; i++) {
            for (int j = i; j < bigger; j++) {
                if (i >= A.vertices || j >= A.vertices) {
                    distance += B.adjacencyMatrix.at(i).at(j);
                }
                else if (i >= B.vertices || j >= B.vertices) {
                    distance += A.adjacencyMatrix.at(i).at(j);
                }
                else
                    distance += abs(A.adjacencyMatrix.at(i).at(j) - B.adjacencyMatrix.at(i).at(j));
            }
        }
    }
    return distance;
}


void UserInterace(vector<Graph>& graphs) {
    int option = 0;
    int i = 0;
    int size = graphs.size();
    while (true) {
        // Display menu
        cout << "\n=== Graph Operations Menu ===\n";
        cout << "1. Show Size of Graph\n";
        cout << "2. Calculate Distance Between Graphs\n";
        cout << "3. Find Maximal Cycle in a Graph\n";
        cout << "4. Minimal Extension of a Graph for Hamiltonian Cycle\n";
        cout << "5. Show All Graphs\n";
        cout << "0. Exit\n";
        cout << "Enter your choice: ";

        // Get user option
        cin >> option;

        // Handle invalid input
        if (cin.fail() || option < 0 || option > 5) {
            cin.clear(); // Clear the error flag
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Ignore invalid input
            cout << "Invalid choice. Please try again.\n";
            continue;
        }

        // Process based on option
        switch (option) {
        case 1: {
            int graphNumber;
            cout << "Enter the graph number to show size: ";
            cin >> graphNumber;
            if (graphNumber >= size || graphNumber < 0) {
                cout << "No such graph!";
                continue;
            }
            cout << "Size of Graph " << GraphSize(graphs.at(graphNumber)) << ".\n";
            break;
        }
        case 2: {
            int graph1, graph2;
            cout << "Enter the first graph number: ";
            cin >> graph1;
            if (graph1 >= size || graph1 < 0) {
                cout << "No such graph!";
                continue;
            }
            cout << "Enter the second graph number: ";
            cin >> graph2;
            if (graph2 >= size || graph2 < 0) {
                cout << "No such graph!";
                continue;
            }
            cout << "Edit distance between Graph " << graph1 << " and Graph " << graph2 << " is: " << CalculateGraphEditDistance(graphs.at(graph1), graphs.at(graph2)) << endl;
            break;
        }
        case 3: {
            int graphNumber;
            cout << "Enter the graph number to find maximal cycles: ";
            cin >> graphNumber;
            if (graphNumber >= size || graphNumber < 0) {
                cout << "No such graph!";
                continue;
            }
            cout << "You chose to find maximal cycles in Graph " << graphNumber << ".\n";
            break;
        }
        case 4: {
            int graphNumber;
            cout << "Enter the graph number for Hamiltonian cycle operations: ";
            cin >> graphNumber;
            if (graphNumber >= size || graphNumber < 0) {
                cout << "No such graph!";
                continue;
            }
            cout << "You chose to find minimal extension and number of Hamiltonian cycles for Graph "
                << graphNumber << ".\n";
            break;
        }
        case 5:
            
            cout << "You chose to display all graphs.\n";
            for (Graph g : graphs) {
                cout << "[" << i << "]" << endl;
                printGraph(g);
                i++;
            }
            break;
        case 0:
            cout << "Exiting program. Goodbye!\n";
            return;
        default:
            cout << "Invalid option. Please try again.\n";
        }
    }
}

int main() {
    string filename;
    cout << system("pwd");
    cout << "Data file path: ";
    cin >> filename;
    
    vector<Graph> graphs = parseGraphs(filename);
    cout << "Parsed " << graphs.size() << " graphs." << endl;

    UserInterace(graphs);

    /*for (size_t i = 0; i < graphs.size(); ++i) {
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
    */


    return 0;
}
