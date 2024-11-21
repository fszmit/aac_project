#include "graph_cycles.hpp"

#include <queue>
#include <iostream>

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