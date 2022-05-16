// Block 0: 1, 2
// 

#include <iostream>
#include "Graph.h"

int main()
{
    //AdjacencyMatrixGraph graph = GenerateRandomGraph(3, 2, false, true);
    //graph.print();

    AdjacencyStructGraph graph(3, false, false);
    graph.AddEdge(2, 0);
    graph.print();
    
}
