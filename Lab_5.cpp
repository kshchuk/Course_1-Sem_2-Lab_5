// Block 0: 1, 2
// Block 1: 7
// Block 2: 13
// Block 3: 14
// Block 4: 17
// Block 5: 19, 20 cannot find algo, may in the free time
// Block 6: 21, some bugs


#include <iostream>
#include "Graph.h"


int main()
{
    AdjacencyMatrixGraph MatrixGraph = GenerateRandomAdjacencyMatrixGraph(15, 15, true, false, true);
    MatrixGraph.print();

    std::cout << "\n\n" << MatrixGraph.hasСycles();

    std::cout << "\n\n";
    std::vector<int> bypass = MatrixGraph.getPathByBFS(5, 10);
    for (auto iter : bypass)
        std::cout << "->"  << iter + 1;
    std::cout << "\n\n";


    AdjacencyStructGraph<true> StructGraph(MatrixGraph);
    StructGraph.print();

    std::cout << "\n\n" << StructGraph.hasCycles();

    std::cout << "\n\n";
    std::vector<int> bypass1 = StructGraph.getPathByBFS(5, 10);
    for (auto iter : bypass1)
        std::cout << "->" << iter + 1;
    std::cout << "\n\n";

    std::vector<int> path = MatrixGraph.dijkstrasAlgorithm(5, 10);
    for (auto iter : path)
        std::cout << "->" << iter + 1;

    std::cout << "\n\n";
    std::vector<std::vector<int>> paths = MatrixGraph.dijkstrasAlgorithm(5);
    for (auto path : paths) {
        for (auto iter : path)
            std::cout << "->" << iter + 1;
        std::cout << "\n";
    }
    std::cout << "\n\n";

    std::vector<int> sortedOrder = StructGraph.topologicalSort();
    for (auto sortedOrder : sortedOrder) {
        std::cout << " " << sortedOrder + 1;
    }
    std::cout << "\n\n";

    AdjacencyMatrixGraph spanningTree = MatrixGraph.kruskal();
    spanningTree.print();
    std::cout << "\n\n";

    std::cout << spanningTree.hasСycles() << "\n\n";
}
