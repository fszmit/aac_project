#include "common.hpp"
#include "graph_cycles.hpp"
#include "GED.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <queue>
#include <unordered_map>
#include <set>
#include <limits>
#include <cmath>

using namespace std;

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
            //cout << "Edit distance between Graph " << graph1 << " and Graph " << graph2 << " is: " << CalculateGraphEditDistance(graphs.at(graph1), graphs.at(graph2)) << endl;
            cout << "Edit distance between Graph " << graph1 << " and Graph " << graph2 << " is: " << astarGED(graphs.at(graph1), graphs.at(graph2)) << endl;


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

    return 0;
}