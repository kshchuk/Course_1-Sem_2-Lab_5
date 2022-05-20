// Block 0: 1, 2
// Block 1: 7
// Block 2: 13
// Block 3: 14


#include <iostream>
#include "Graph.h"

int main()
{
    AdjacencyMatrixGraph graph = GenerateRandomAdjacencyMatrixGraph(15, 15, true, false, false);
    graph.print();

    std::cout << "\n\n";
    std::vector<int> bypass = graph.getPathByBFS(5, 10);
    for (auto iter : bypass)
        std::cout << "->"  << iter + 1;
    std::cout << "\n\n";


    AdjacencyStructGraph<false> graph1(graph);
    graph1.print();

    std::cout << "\n\n";
    std::vector<int> bypass1 = graph1.getPathByBFS(5, 10);
    for (auto iter : bypass1)
        std::cout << "->" << iter + 1;
    std::cout << "\n\n";

    std::vector<int> path = graph.dijkstrasAlgorithm(5, 10);
    for (auto iter : path)
        std::cout << "->" << iter + 1;

    std::cout << "\n\n";
    std::vector<std::vector<int>> paths = graph.dijkstrasAlgorithm(5);
    for (auto path : paths) {
        for (auto iter : path)
            std::cout << "->" << iter + 1;
        std::cout << "\n";
    }
    std::cout << "\n\n";
    //std::cout << '\n' << graph.hasСycles();
    //AdjacencyStructGraph<false> graph1(graph);
    //graph1.print();
    //std::cout << '\n' << graph1.hasCycles();

    
}
