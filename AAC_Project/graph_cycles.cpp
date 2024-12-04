#include "graph_cycles.hpp"

#include <queue>
#include <iostream>
#include <cmath>
#include <algorithm>

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

int MatTrace(std::vector<std::vector<int>> mat) {
    int res = 0;
    for (ul i = 0; i < mat.size(); i++)
    {
        res += mat[i][i];
    }
    return res;
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

Graph SubGraph(const Graph * const g, std::vector<int> mask, int * neighbours = NULL) {
    //std::cout << "mask : " << mask.size() << '\n';
    if (mask.size() == 0) {
        if (neighbours != NULL) {
            *neighbours = 0;
        }
        return CloneGraph(g);
    }
    if (mask.size() == (unsigned long)g->vertices) {
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
        if ( (ul)(i-si) < mask.size() && i == mask[i - si]) { // masked row // !!-->> VALGRIND
            for (int j = 0; j < g->vertices; j++)
            {
                if ( (ul)(j-sj) == mask.size() || j != mask.at(j - sj)) { // masked row, unmasked column // !!-->> VALGRIND
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
            if ( (ul)(j-sj) < mask.size() && j == mask[j - sj]){ // unmasked row, masked column // !!-->> VALGRIND
                if (g->adjacencyMatrix[i][j] >= 1) {
                    n[j - sj]++;
                }
                continue;
            }
            //std::cout << "s: " << si << " : " << sj << " -- " << "b:" << i << " : " << j << '\n';
            //unmasked row, unmasked column
            s.adjacencyMatrix.at(si).at(sj) = g->adjacencyMatrix.at(i).at(j);
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

int MaxCycle(const Graph * const g, int * no_cycles = NULL) {
    // based on https://arxiv.org/pdf/1612.05531
    // NOTE:: this returns the number of **directed** cycles. so a singular cycle in an undirected graph will count for two

    Graph flat;
    flat.vertices = g->vertices;

    flat.adjacencyMatrix = FlattenMatrix(g->adjacencyMatrix);

    flat.additionalInfo = g->additionalInfo;

    for (int i = flat.vertices; i >= 3; i--)
    {   
        int gamma = 0;
        for (int j = flat.vertices; j >=0 ; j--)
        {
            std::vector<std::vector<int>> combs = comb(flat.vertices, j);
            for (int k = 0; k < (int)combs.size(); k++)
            {
                int sub_gamma = 0;
                int neigh = 0;
                std::vector<int> mask = combs[k];
                Graph sub = SubGraph(&flat, mask, &neigh);
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

Graph PruneGraph(const Graph * const g) {
    bool removed = true;
    Graph sub = CloneGraph(g);

    while (removed && sub.vertices > 0)
    {
        removed = false;
        std::vector<size_t> conn(sub.vertices, 0);
        std::vector<size_t> no_neigh(sub.vertices, 0);
        for (int i = 0; i < sub.vertices; i++)
        {
            for (int j = 0; j < sub.vertices; j++)
            {
                if (i != j && sub.adjacencyMatrix.at(i).at(j) >= 1)
                {
                    conn[i] |= 1;
                    conn[j] |= 2;
                }
                if (i != j && (sub.adjacencyMatrix.at(j).at(i) >=1 || (conn[i] & 1) == 1 )) {
                    no_neigh[i]++;
                }
            }
        }
        std::vector<int> mask;

        for (size_t i = 0; i < conn.size(); i++)
        {
            if (conn[i] != 3 || no_neigh[i] <= (size_t)1) {
                mask.push_back(i);
                removed = true;
            }
        }
        if (!removed) {
            break;
        }
        int neigh = 0;
        Graph tmp = SubGraph(&sub, mask, &neigh);
        sub = tmp;
    }
    return sub;
}

void KosarajuVisit(const Graph * const g, const int vert, std::vector<int> * const vert_list, std::vector<bool> * const visited) {
    if (visited->at(vert)) {
        return;
    }
    visited->at(vert) = true;
    for (int i = 0; i < g->vertices; i++)
    {
        if (i != vert && g->adjacencyMatrix[vert][i] >= 1) {
            KosarajuVisit(g, i, vert_list, visited);
        }
    }
    vert_list->push_back(vert);
    return;
}

void KosarajuAssign(const Graph * const g, const int vert, const int root, std::vector<int> * const roots, std::vector<int> * const neg_mask) {
    if ( roots->at(vert) != -1 ) {
        return;
    }
    roots->at(vert) = root;
    neg_mask->push_back(vert);
    for (int i = 0; i < g->vertices; i++)
    {
        if (i != vert && g->adjacencyMatrix[i][vert] >= 1)
        {
            KosarajuAssign(g, i, root, roots, neg_mask);
        }
    }
    return;
}

std::vector<int> InvertMask(const std::vector<int> * const mask, const int no_verts){
    std::vector<int> res_mask;
    int point = 0;
    for (int i = 0; i < no_verts; i++)
    {
        if ( point < (int)mask->size() && i == mask->at(point) ) {
            point++;
            continue;
        }
        res_mask.push_back(i);
    }
    return res_mask;
}

void APXLongestSimpleCycleTraversal(const Graph * const g, int vertex, std::vector<int> * const col, std::vector<int> * const pred) {
    col->at(vertex) = 1;
    for (int i = 0; i < g->vertices; i++)
    {
        if (i == vertex || g->adjacencyMatrix[vertex][i] == 0)
        {
            continue;
        }
        if (col->at(i)==0)
        {
            pred->at(vertex);
            APXLongestSimpleCycleTraversal(g, i, col, pred);
        }
    }
}

int APXLongestSimpleCycle(const Graph * const g, int max_vert) {
    std::vector<int> col(g->vertices, 0);
    std::vector<int> pred(g->vertices, -1);

    for (int i = 0; i < g->vertices; i++)
    {
        if (max_vert == i || g->adjacencyMatrix[max_vert][i] == 0){
            continue;
        }
        if (col[i] == 0){
            APXLongestSimpleCycleTraversal(g, i, &col, &pred);
        }
    }

    int depth = 1;
    int prev = pred[max_vert];
    while (prev != max_vert)
    {
        prev = pred[prev];
        depth ++;
    }
    
    return depth + 1;
}

int MaxOutDegree(const Graph * const g, double * const avg_degree){
    int max_out = -1;
    int v = -1;
    int deg_total = 0;
    for (int i = 0; i < g->vertices; i++)
    {
        int sum = 0;
        for (int j = 0; j < g->vertices; j++)
        {
            sum += g->adjacencyMatrix[i][j] >= 1 ? 1 : 0;
        }
        deg_total += sum;
        if (sum > max_out)
        {
            max_out = sum;
            v = i;
        }
    }
    if (avg_degree != NULL)
    {
        *avg_degree = ((double)deg_total) / ((double)g->vertices);
    }
    return v;
}

int APXMaxCycle(const Graph * const g, int * no_cycles = NULL){

    Graph sub = PruneGraph(g); // O(V^3)

    if (sub.vertices == 0) { // acyclic graph
        if (no_cycles != NULL) {
            *no_cycles = 0;
        }
        return 0;
    }

    // we find strongly connected components, as a cycle will ever cross a bridge 
    // Kosaraju's algorithm -> O(V + E)

    std::vector<int> vert_stack;
    std::vector<bool> vsited(sub.vertices, false);

    for (int i = 0; i < sub.vertices; i++)
    {
        KosarajuVisit(&sub, i, &vert_stack, &vsited);
    }
    std::reverse(vert_stack.begin(), vert_stack.end());

    std::vector<int> roots(sub.vertices, -1);
    std::vector<std::vector<int> > components;

    for (int i = 0; i < (int)vert_stack.size(); i++)
    {
        std::vector<int> neg_mask;
        KosarajuAssign(&sub, i, i, &roots, &neg_mask);
        if (neg_mask.size() != 0){
            components.push_back(neg_mask);
        }
    }
    
    std::sort(components.begin(), components.end(), [](std::vector<int> &a, std::vector<int> &b) -> bool { return a.size() >= b.size() ? true : false; } );

    // std::cout << "components:\n";    
    // for (int i = 0; i < (int)components.size(); i++)
    // {
    //     for (int j = 0; j < (int)components[i].size(); j++)
    //     {
    //         std::cout << (int)components[i][j] << " ";
    //     }
    //     std::cout << "\n";
    // }

    int max_cycle = 0;
    int cycle_count = 0;

    for (int i = 0; i < (int)components.size(); i++)
    {
        std::sort(components[i].begin(), components[i].end());
        std::vector<int> mask = InvertMask(&(components[i]), sub.vertices);
        Graph comp = SubGraph(&sub, mask, NULL);
        double avg_degree = 0;
        int max = MaxOutDegree(&comp, &avg_degree);

        // Kumar, Parveen & Gupta, Nitin. (2014). A Heuristic Algorithm for Longest Simple Cycle Problem. 
        int max_cycle_len = APXLongestSimpleCycle(g, max);
        if (max_cycle > max_cycle_len) { continue; }
        if (max_cycle < max_cycle_len) {
            max_cycle = max_cycle_len;
            cycle_count = 0;
        }
        cycle_count += (int)std::lround( ((double)binom(comp.vertices, max_cycle)) * std::tgamma(avg_degree) );

    }
    
    if (no_cycles != NULL) {
        *no_cycles = cycle_count;
    }

    //estimation w/ the gamma funciton
    return max_cycle;

}