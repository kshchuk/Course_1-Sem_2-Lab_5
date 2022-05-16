#pragma once
#include <random>

// Block 0: 1
class AdjacencyMatrixGraph
{
private:
	int** adjacency_matrix;
	int vertices_number = 0;

public:
	AdjacencyMatrixGraph(size_t vertices_number);

	bool AddEdge(size_t vertix_1, size_t vertix_2, int weight = 1, bool isOriented = true);	
	void print();
};

AdjacencyMatrixGraph GenerateRandomAdjacencyMatrixGraph(size_t verticesNum,
										 size_t edgeNum, 
										 bool isOriented = false, 
										 bool hasLoops = false,
										 bool isWeighted = false, 
										 size_t maxWeight = 10);



// Block 0: 2
class AdjacencyStructGraph
{
private:
	struct Edge
	{
		int destination, weight;
	};
	std::vector<std::vector<Edge>>* adjacency_weighted_graph_list;
	std::vector<std::vector<int>>* adjacency_list;
	int vertices_number = 0;
	bool isWeighted;

public:
	AdjacencyStructGraph(size_t vertices_number, bool isWeighted = false);
	bool AddEdge(size_t vertix_1, size_t vertix_2, int weight = 0);
	void print();
};

AdjacencyStructGraph GenerateRandomAdjacencyStructGraphGraph(size_t verticesNum,
	size_t edgeNum,
	bool isOriented = false,
	bool hasLoops = false,
	bool isWeighted = false,
	size_t maxWeight = 10);
