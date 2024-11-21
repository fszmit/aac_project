#include "common.hpp"
#include "graph_cycles.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <queue>
#include <algorithm>

using namespace std;

int main() {
    string filename;
    cout << system("pwd");
    cout << "Data file path: ";
    cin >> filename;
    
    vector<Graph> graphs = parseGraphs(filename);
    cout << "Parsed " << graphs.size() << " graphs." << endl;

    for (size_t i = 0; i < graphs.size(); ++i) {
        cout << "====\n";
        cout << "\nGraph " << i + 1 << ":" << endl;
        printGraph(graphs[i]);
        cout << "\nSubgraph:\n";


        int no_cycles = 0;
        cout << "Max Cycle size: " << MaxCycle(&(graphs[i]), &no_cycles) << '\n';
        cout << "No Max Cycles: " << no_cycles << '\n';

        cout << "====";
        cout << "\n\n";
    }


    return 0;
}
