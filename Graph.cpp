#include "Graph.h"
#include <iostream>

AdjacencyMatrixGraph::AdjacencyMatrixGraph(size_t vertices_number) {
	adjacency_matrix = new int* [vertices_number + 1];
	for (size_t i = 0; i < vertices_number; i++) {
		adjacency_matrix[i] = new int[vertices_number + 1];
		for (size_t j = 0; j < vertices_number; j++)
			adjacency_matrix[i][j] = 0;
	}
	this->vertices_number = vertices_number;
}

bool AdjacencyMatrixGraph::AddEdge(size_t vertix_1, size_t vertix_2, int weight, bool isOriented)
{
		adjacency_matrix[vertix_1][vertix_2] += weight;
		if (!isOriented && vertix_1 != vertix_2)
			adjacency_matrix[vertix_2][vertix_1] += weight;
		return true;
}

void AdjacencyMatrixGraph::print()
{
	for (size_t i = 0; i < vertices_number; i++) {
		std::cout << '\n';
		for (size_t j = 0; j < vertices_number; j++)
			std::cout << adjacency_matrix[i][j] << ' ';
	}
}

AdjacencyMatrixGraph GenerateRandomAdjacencyMatrixGraph(size_t verticesNum, size_t edgeNum, bool isOriented,
	bool hasLoops, bool isWeighted, size_t maxWeight)
{
	AdjacencyMatrixGraph graph(verticesNum);

	std::random_device rd;
	std::mt19937 mersenne(rd());

	while (edgeNum) {
		size_t vertice_1 = mersenne() % verticesNum,
			vertice_2 = mersenne() % verticesNum,
			weight = isWeighted ? mersenne() % maxWeight : 1;
		if (!hasLoops && vertice_1 == vertice_2)
			continue;

		graph.AddEdge(vertice_1, vertice_2, weight);
		edgeNum--;
	}
	return graph;
}

AdjacencyStructGraph GenerateRandomGraph(size_t verticesNum, size_t edgeNum, 
	bool isOriented, bool hasLoops, bool isWeighted, size_t maxWeight)
{
	AdjacencyStructGraph(verticesNum, isOriented, isWeighted);
}

AdjacencyStructGraph::AdjacencyStructGraph(size_t vertices_number, bool isWeighted)
{
	if (isWeighted)
		adjacency_weighted_graph_list.resize(vertices_number);
	else
		adjacency_list.resize(vertices_number);

	this->vertices_number = vertices_number;
	this->isOriented = isOriented;
	this->isWeighted = isWeighted;
}

bool AdjacencyStructGraph::AddEdge(size_t vertix_1, size_t vertix_2, int weight)
{
	if (vertix_1 < vertices_number && vertix_2 < vertices_number) {
		if (isWeighted) {
			Edge edge_1;
			edge_1.destination = vertix_2;
			edge_1.weight = weight;
			adjacency_weighted_graph_list[vertix_1].push_back(edge_1);
			if (!isOriented) {
				Edge edge_2;
				edge_2.destination = vertix_1;
				edge_2.weight = weight;
				adjacency_weighted_graph_list[vertix_2].push_back(edge_2);
			}
		}
		else {
			adjacency_list[vertix_1].push_back(vertix_2);
			if (!isOriented)
				adjacency_list[vertix_2].push_back(vertix_1);
		}
		return true;
	}
	else
		return false;
}

void AdjacencyStructGraph::print()
{
	if (isWeighted)
		for (int i = 0; i < vertices_number; i++) {
			std::cout << '\n' << i + 1;
			for (auto j : adjacency_weighted_graph_list[i])
				std::cout << "->(" << j.destination + 1 << ' ' << j.weight << ')';
		}
	else
		for (int i = 0; i < vertices_number; i++) {
			std::cout << '\n' << i + 1;
			for (auto j : adjacency_list[i])
				std::cout << "->" << j + 1;
		}
}
