#pragma once
#include <vector>
#include <queue>
#include <list>
#include <type_traits>


template<bool isWeight>
class AdjacencyStructGraph;

// Block 0: 1
class AdjacencyMatrixGraph
{
private:
	int** adjacency_matrix;
	int vertices_number = 0;
	bool isOriented, isWeighted;

	template<bool>
	friend class AdjacencyStructGraph;

	// Block 1: 7
	bool hasCyclesByDFS(size_t vertix, std::vector<int>& colors, std::vector<int>& from);

	// Returns path existence
	bool bfs(size_t to, std::vector<bool>& used, std::list<int>& buffer,
					std::vector<int>& parents);

	// Block 2: 13
	// Gets path (Recursive function)
	 void get_path(size_t from, size_t to, std::vector<int> parents, std::vector<int>& path);

public:
	AdjacencyMatrixGraph(size_t vertices_number, bool isOriented = false, bool isWeighted = false);
	AdjacencyMatrixGraph(AdjacencyStructGraph<true> graph);
	AdjacencyMatrixGraph(AdjacencyStructGraph<false> graph);

	bool AddEdge(size_t vertix_1, size_t vertix_2, int weight = 1);
	void print();

	// Block 1: 7
	// Uses hasCyclesByDFS()
	bool hasСycles();

	// Block 2: 13
	// Returns path 
	std::vector<int> getPathByBFS(size_t from, size_t to);

	// Block 3: 14
	// Returns path
	std::vector<int> dijkstrasAlgorithm(size_t from, size_t to);
	// Return all paths
	std::vector<std::vector<int>> dijkstrasAlgorithm(size_t from);
};

AdjacencyMatrixGraph GenerateRandomAdjacencyMatrixGraph(size_t verticesNum,
	size_t edgeNum,
	bool isOriented = false,
	bool hasLoops = false,
	bool isWeighted = false,
	size_t maxWeight = 10);



// Block 0: 2
template <bool isWeighted>
class AdjacencyStructGraph
{
private:
	struct Edge
	{
		int destination, weight;
		Edge(int destination, int weight)
			: destination(destination), weight(weight) {}
	};

	typedef std::conditional_t<isWeighted, Edge, int> edge_data_type;
	std::vector<std::list<edge_data_type>> adjacency_list;
	bool isOriented;

	friend class AdjacencyMatrixGraph;

	// Helping structure for graph passing or another way of graph presentation
	struct Node
	{
		edge_data_type* edge;
		Node* parent;
		std::list<Node*> nodes;
		bool is_visited;

		Node(edge_data_type edge) :
			edge(new edge_data_type(edge)),
			is_visited(false) {}

		~Node() {
			delete edge;
		}
	};

	Node* node_by_number(std::list<Node*>& graph, int edge);

	// Block 2: 13
	// Gets path (Recursive function)
	void get_path(Node* from, Node* to, std::vector<int>& path);

	// Returns path existence
	bool bfs(std::list<Node*>& graph, Node* to, std::list<Node*>& buffer);

	// Block 1: 7
	bool hasCyclesByDFS(size_t vertix, std::vector<int>& colors, std::vector<int>& from);

public:
	AdjacencyStructGraph(size_t vertices_number, bool isOriented = false);
	AdjacencyStructGraph(AdjacencyMatrixGraph graph);

	bool AddEdge(size_t vertix_1, size_t vertix_2, int weight = 1);
	bool hasCycles();

	// Block 2: 13
	// Returns path 
	std::vector<int> getPathByBFS(size_t from, size_t to);

	void print();
};

template<bool isWeighted>
AdjacencyStructGraph<isWeighted> GenerateRandomAdjacencyStructGraph(size_t verticesNum,
	size_t edgeNum,
	bool isOriented = true,
	bool hasLoops = false,
	size_t maxWeight = 10);

#include "Graph.ipp"
