#include "Graph.h"
#include <iostream>
#include <random>
#include <map>

#define INFINITY INT_MAX


AdjacencyMatrixGraph::AdjacencyMatrixGraph(size_t vertices_number, bool isOriented, bool isWeighted) 
{
	for (size_t i = 0; i < vertices_number; i++) {
		adjacency_matrix.push_back(std::vector<int>(vertices_number, 0));
	}
	this->vertices_number = vertices_number;
	this->isOriented = isOriented;
	this->isWeighted = isWeighted;
}

inline AdjacencyMatrixGraph::AdjacencyMatrixGraph(std::vector<std::vector<int>> matrix)
{
	this->adjacency_matrix = matrix;
	this->vertices_number = matrix.size();
}

inline AdjacencyMatrixGraph::AdjacencyMatrixGraph(AdjacencyStructGraph<true> graph)
{
	this->isOriented = graph.isOriented;
	this->isWeighted = true;

	vertices_number = graph.adjacency_list.size();
	for (size_t i = 0; i < vertices_number; i++) {
		adjacency_matrix.push_back(std::vector<int>(vertices_number, 0));
	}

	for (size_t i = 0; i < vertices_number; i++)
		for (auto edge : graph.adjacency_list[i])
			adjacency_matrix[i][edge.destination] += edge.weight;

	// Preventing loop recurrence
	for (size_t i = 0; i < vertices_number; i++)
		adjacency_matrix[i][i] /= 2;
}

inline AdjacencyMatrixGraph::AdjacencyMatrixGraph(AdjacencyStructGraph<false> graph)
{
	this->isWeighted = false;
	this->isOriented = graph.isOriented;

	vertices_number = graph.adjacency_list.size();
	for (size_t i = 0; i < vertices_number; i++) {
		adjacency_matrix.push_back(std::vector<int>(vertices_number, 0));
	}

	for (size_t i = 0; i < vertices_number; i++)
		for (auto edge : graph.adjacency_list[i])
			adjacency_matrix[i][edge] += 1;

	// Preventing loop recurrence
	for (size_t i = 0; i < vertices_number; i++)
		adjacency_matrix[i][i] /= 2;
}

bool AdjacencyMatrixGraph::AddEdge(size_t vertix_1, size_t vertix_2, int weight)
{
	if (adjacency_matrix[vertix_1][vertix_2])
		return false;
	adjacency_matrix[vertix_1][vertix_2] = weight;
	if (!isOriented && vertix_1 != vertix_2)
		AddEdge(vertix_2, vertix_1, weight);
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

inline bool AdjacencyMatrixGraph::hasСycles()
{
	std::vector<int> colors, from;
	colors.resize(vertices_number); from.resize(vertices_number);
	for (size_t i = 0; i < vertices_number; i++)
		if (colors[i] == 0 && hasCyclesByDFS(i, colors, from))
			return true;
	return false;
}

inline std::vector<int> AdjacencyMatrixGraph::getPathByBFS(size_t from, size_t to)
{
	std::vector<bool> used;				used.resize(vertices_number);
	std::list<int> buffer;				buffer.push_back(from);
	std::vector<int> parents;			parents.resize(vertices_number);

	if (bfs(to, used, buffer, parents)) {
		std::vector<int> path;
		get_path(from, to, parents, path);
		return path;
	}
	else {
		return std::vector<int>(NULL);
	}
}

inline std::vector<int> AdjacencyMatrixGraph::dijkstrasAlgorithm(size_t from, size_t to)
{
	int** cost = new int*[vertices_number];
	for (size_t i = 0; i < vertices_number; i++) {
		cost[i] = new int[vertices_number];
	}
	int* visited = new int[vertices_number], * pred = new int[vertices_number], * distance = new int[vertices_number];
	int count, mindistance, i, j;
	int nextnode;

	for (i = 0; i < vertices_number; i++)
		for (j = 0; j < vertices_number; j++)
			if (adjacency_matrix[i][j] == 0)
				cost[i][j] = INFINITY;
			else
				cost[i][j] = adjacency_matrix[i][j];

	for (i = 0; i < vertices_number; i++) {
		distance[i] = cost[from][i];
		pred[i] = from;
		visited[i] = 0;
	}
	distance[from] = 0;
	visited[from] = 1;

	count = 1;
	while (count < vertices_number - 2) {
		mindistance = INFINITY;
		for (i = 0; i < vertices_number; i++)
			if (distance[i] < mindistance && !visited[i]) {
				mindistance = distance[i];
				nextnode = i;
			}
		try {
			visited[nextnode] = 1;
		}
		catch (...) {
			visited[nextnode] = 1;
		}
		for (i = 0; i < vertices_number; i++)
			if (!visited[i])
				if (mindistance + cost[nextnode][i] < distance[i]) {
					distance[i] = mindistance + cost[nextnode][i];
					pred[i] = nextnode;
				}
		count++;
	}

	std::vector<int> path;
	path.push_back(to);
	j = to;
	do {
		j = pred[j];
		path.push_back(j);
	} while (j != from);
	std::reverse(path.begin(), path.end());

	for (size_t i = 0; i < vertices_number; i++) {
		delete[] cost[i];
	}

	delete[] cost, visited, pred, distance;

	return path;
}

inline std::vector<std::vector<int>> AdjacencyMatrixGraph::dijkstrasAlgorithm(size_t from)
{
	int** cost = new int* [vertices_number];
	for (size_t i = 0; i < vertices_number; i++) {
		cost[i] = new int[vertices_number];
	}
	int* visited = new int[vertices_number], * pred = new int[vertices_number], * distance = new int[vertices_number];
	int count, mindistance, nextnode, i, j;

	for (i = 0; i < vertices_number; i++)
		for (j = 0; j < vertices_number; j++)
			if (adjacency_matrix[i][j] == 0)
				cost[i][j] = INFINITY;
			else
				cost[i][j] = adjacency_matrix[i][j];

	for (i = 0; i < vertices_number; i++) {
		distance[i] = cost[from][i];
		pred[i] = from;
		visited[i] = 0;
	}
	distance[from] = 0;
	visited[from] = 1;

	count = 1;
	while (count < vertices_number - 1) {
		mindistance = INFINITY;
		for (i = 0; i < vertices_number; i++)
			if (distance[i] < mindistance && !visited[i]) {
				mindistance = distance[i];
				nextnode = i;
			}
		visited[nextnode] = 1;
		for (i = 0; i < vertices_number; i++)
			if (!visited[i])
				if (mindistance + cost[nextnode][i] < distance[i]) {
					distance[i] = mindistance + cost[nextnode][i];
					pred[i] = nextnode;
				}
		count++;
	}

	std::vector < std::vector<int>> paths;
	std::vector<int> path;
	for (i = 0; i < vertices_number; i++)
		if (i != from) {
			j = i;
			path.push_back(j);

			do {
				j = pred[j];
				path.push_back(j);
			} while (j != from);

			std::reverse(path.begin(), path.end());
			paths.push_back(path);
			path.clear();
		}

	delete[] cost, visited, pred, distance;

	return paths;
}

inline AdjacencyMatrixGraph AdjacencyMatrixGraph::kruskal()
{
	std::vector< std::vector<int> > tempMatrix(this->adjacency_matrix);
	std::vector< std::vector<int> > spanningTreeMatr(vertices_number);
	std::vector<int> belongsToTree(vertices_number);

	for (int i = 0; i < vertices_number; i++) {
		spanningTreeMatr[i] = std::vector<int>(vertices_number, INFINITY);
		belongsToTree[i] = i;
	}

	int nodoA;
	int nodoB;
	int arcos = 1;
	while (arcos < vertices_number - 1) {

		int min = INFINITY;
		for (int i = 0; i < vertices_number; i++)
			for (int j = 0; j < vertices_number; j++)
				if (min > tempMatrix[i][j] && belongsToTree[i] != belongsToTree[j] && tempMatrix[i][j] != 0) {
					min = tempMatrix[i][j];
					nodoA = i;
					nodoB = j;
				}

		if (belongsToTree[nodoA] != belongsToTree[nodoB]) {
			spanningTreeMatr[nodoA][nodoB] = min;
			spanningTreeMatr[nodoB][nodoA] = min;

			int temp = belongsToTree[nodoB];
			belongsToTree[nodoB] = belongsToTree[nodoA];
			for (int k = 0; k < vertices_number; k++)
				if (belongsToTree[k] == temp)
					belongsToTree[k] = belongsToTree[nodoA];
			arcos++;
		}
	}
	for (auto& line : spanningTreeMatr) {
		for (auto& edge : line) {
			if (edge == INFINITY)
				edge = 0;
		}
	}

	return AdjacencyMatrixGraph(spanningTreeMatr);
}

inline bool AdjacencyMatrixGraph::bfs(size_t to, std::vector<bool>& used, std::list<int>& buffer, 
										   std::vector<int>& parents)
{
	if (buffer.empty())
		return false;

	size_t from = buffer.front();
	buffer.pop_front();

	if (from == to)
		return true;

	used[from] = true;

	if (isWeighted) {
		// weight - key, vertex - value
		// Automaticaly sorts by a weight
		std::map<int, int> adjacent_vertices;
		for (size_t i = 0; i < vertices_number; ++i) {
			if (adjacency_matrix[from][i])
				adjacent_vertices.insert({ adjacency_matrix[from][i], i });
		}
		for (auto iter : adjacent_vertices)
		{
			if (iter.first == 0)
				continue;
			if (used[iter.second])
				continue;
			if (std::find(buffer.begin(), buffer.end(), iter.second) != buffer.end())
				continue;
			parents[iter.second] = from;
			buffer.push_back(iter.second);
		}
	}
	else {
		for (size_t i = 0; i < vertices_number; i++) {
			if (adjacency_matrix[from][i] == 0)
				continue;
			if (used[i])
				continue;
			if (std::find(buffer.begin(), buffer.end(), i) != buffer.end())
				continue;
			parents[i] = from;
			buffer.push_back(i);
		}
	}
	return bfs(to, used, buffer, parents);
}

inline void AdjacencyMatrixGraph::get_path(size_t from, size_t to, std::vector<int> parents, std::vector<int>& path)
{
	if (to == from) {
		path.push_back(to);
		return;
	}
	int prefrom = parents[to];
	get_path(from, prefrom, parents, path);
	path.push_back(to);
}


inline bool AdjacencyMatrixGraph::hasCyclesByDFS(size_t vertex, std::vector<int>& colors, std::vector<int>& from)
{
	colors[vertex] = 1;
	for (size_t i = 0; i < vertices_number; i++) {
		if (adjacency_matrix[vertex][i]) {
			if (colors[i] == 0) {
				hasCyclesByDFS(i, colors, from);
				from[i] = vertex;
			}
			else if (colors[i] == 1 && (i != from[vertex] || isOriented))
				return true;
		}
	}
	colors[vertex] = 2;
	return false;
}

AdjacencyMatrixGraph GenerateRandomAdjacencyMatrixGraph(size_t verticesNum, size_t edgeNum, bool isOriented,
	bool hasLoops, bool isWeighted, size_t maxWeight)
{
	AdjacencyMatrixGraph graph(verticesNum, isOriented, isWeighted);

	std::random_device rd;
	std::mt19937 mersenne(rd());

	while (edgeNum) {
		size_t vertice_1 = mersenne() % verticesNum,
			vertice_2 = mersenne() % verticesNum,
			weight = isWeighted ? mersenne() % maxWeight + 1 : 1;
		if (!hasLoops && vertice_1 == vertice_2)
			continue;

		if (graph.AddEdge(vertice_1, vertice_2, weight))
			edgeNum--;
	}
	return graph;
}

template<bool isWeighted>
inline bool AdjacencyStructGraph<isWeighted>::bfs(std::list<Node*>& graph, Node* to, std::list<Node*>& buffer)
{
	if (buffer.empty())
		return false;

	Node* from = buffer.front();
	buffer.pop_front();
	from->is_visited = true;
	if (from == to) {
		return true;
	}
	for (typename std::list<Node*>::iterator it = from->nodes.begin(); it != from->nodes.end(); ++it) {
		Node* node = *it;
		if (node->is_visited == true)
			continue;
		if (find(buffer.begin(), buffer.end(), node) != buffer.end())
			continue;
		node->parent = from;
		buffer.push_back(node);
	}
	return bfs(graph, to, buffer);
}

inline bool AdjacencyStructGraph<false>::hasCyclesByDFS(size_t vertex, std::vector<int>& colors, std::vector<int>& from)
{
	colors[vertex] = 1;
	for (auto to : adjacency_list[vertex]) {
		if (colors[to] == 0) {
			hasCyclesByDFS(to, colors, from);
			from[to] = vertex;
		}
		else if (colors[to] == 1 && (to != from[vertex] || isOriented))
			return true;
	}
	colors[vertex] = 2;
	return false;
}

inline bool AdjacencyStructGraph<true>::hasCyclesByDFS(size_t vertex, std::vector<int>& colors, std::vector<int>& from)
{
	colors[vertex] = 1;
	for (auto to : adjacency_list[vertex]) {
		if (colors[to.destination] == 0) {
			hasCyclesByDFS(to.destination, colors, from);
			from[to.destination] = vertex;
		}
		else if (colors[to.destination] == 1 && (to.destination != from[vertex] || isOriented))
			return true;
	}
	colors[vertex] = 2;
	return false;
}

inline void AdjacencyStructGraph<false>::topologicalSortUtil(int vertex, bool visited[], std::stack<int>& Stack)
{
	visited[vertex] = true;

	for (auto edge : adjacency_list[vertex])
		if (!visited[edge])
			topologicalSortUtil(edge, visited, Stack);

	Stack.push(vertex);
}

inline void AdjacencyStructGraph<true>::topologicalSortUtil(int vertex, bool visited[], std::stack<int>& Stack)
{
	visited[vertex] = true;

	for (auto edge : adjacency_list[vertex])
		if (!visited[edge.destination])
			topologicalSortUtil(edge.destination, visited, Stack);

	Stack.push(vertex);
}

template<bool isWeighted>
inline AdjacencyStructGraph<isWeighted>::AdjacencyStructGraph(size_t vertices_number,  bool isOriented)
{
	adjacency_list.resize(vertices_number);
	this->isOriented = isOriented;
}

inline AdjacencyStructGraph<true>::AdjacencyStructGraph(AdjacencyMatrixGraph graph)
{
	adjacency_list.resize(graph.vertices_number);
	for (size_t i = 0; i < graph.vertices_number; i++)
		for (size_t j = 0; j < graph.vertices_number; j++)
			if (graph.adjacency_matrix[i][j] > 0)
				adjacency_list[i].push_back(Edge(j, graph.adjacency_matrix[i][j]));
	this->isOriented = graph.isOriented;
}

template<bool isWeighted>
inline std::vector<int> AdjacencyStructGraph<isWeighted>::topologicalSort()
{
	size_t vertexNum = adjacency_list.size();
	std::stack<int> Stack;

	bool* visited = new bool[vertexNum];
	for (int i = 0; i < vertexNum; i++)
		visited[i] = false;

	for (int i = 0; i < vertexNum; i++)
		if (visited[i] == false)
			topologicalSortUtil(i, visited, Stack);

	std::vector<int> sorted;
	while (Stack.empty() == false)
	{

		sorted.push_back(Stack.top());
		Stack.pop();
	}

	delete[] visited;
	return sorted;
}

inline AdjacencyStructGraph<false>::Node* AdjacencyStructGraph<false>::node_by_number(std::list<Node*>& graph, int edge)
{
	for (typename std::list<Node*>::iterator it = graph.begin(); it != graph.end(); ++it) {
		if (*(*it)->edge == edge)
			return *it;
	}
	return nullptr;
}

inline void AdjacencyStructGraph<false>::get_path(Node* from, Node* to, std::vector<int>& path)
{
	if (to == from) {
		path.push_back(*to->edge);
		return;
	}

	get_path(from, to->parent, path);
	path.push_back(*to->edge);
}

inline void AdjacencyStructGraph<true>::get_path(Node* from, Node* to, std::vector<int>& path)
{
	if (to == from) {
		path.push_back(to->edge->destination);
		return;
	}

	get_path(from, to->parent, path);
	path.push_back(to->edge->destination);
}

inline AdjacencyStructGraph<true>::Node* AdjacencyStructGraph<true>::node_by_number(std::list<Node*>& graph, int edge)
{
	for (typename std::list<Node*>::iterator it = graph.begin(); it != graph.end(); ++it) {
		if ((*it)->edge->destination == edge)
			return *it;
	}
	return nullptr;
}

inline AdjacencyStructGraph<false>::AdjacencyStructGraph(AdjacencyMatrixGraph graph)
{
	adjacency_list.resize(graph.vertices_number);
	for (size_t i = 0; i < graph.vertices_number; i++)
		for (size_t j = 0; j < graph.vertices_number; j++)
			if (graph.adjacency_matrix[i][j] > 0) {
				int decreasing = graph.adjacency_matrix[i][j];
				while (decreasing > 0) {
					adjacency_list[i].push_back(j);
					decreasing--;
				}
			}
	this->isOriented = graph.isOriented;
}

template<bool isWeighted>
inline bool AdjacencyStructGraph<isWeighted>::hasCycles()
{
	std::vector<int> colors, from;
	colors.resize(adjacency_list.size()); 
	from.resize(adjacency_list.size());
	for (size_t i = 0; i < adjacency_list.size(); i++)
		if (colors[i] == 0 && hasCyclesByDFS(i, colors, from))
			return true;
	return false;
}

inline bool AdjacencyStructGraph<true>::AddEdge(size_t vertix_1, size_t vertix_2, int weight)
{
	for (auto edge : adjacency_list[vertix_1])
		if (edge.destination == vertix_2)
			return false;
	adjacency_list[vertix_1].push_back(Edge(vertix_2, weight));
	if (!isOriented)
		adjacency_list[vertix_2].push_back(Edge(vertix_1, weight));
	return true;
}

inline bool AdjacencyStructGraph<false>::AddEdge(size_t vertix_1, size_t vertix_2, int weight)
{
	for (auto edge : adjacency_list[vertix_1])
		if (edge == vertix_2)
			return false;
	adjacency_list[vertix_1].push_back(vertix_2);
	if (!isOriented)
		adjacency_list[vertix_2].push_back(vertix_1);
	return true;
}

inline void AdjacencyStructGraph<true>::print()
{
	for (int i = 0; i < adjacency_list.size(); i++) {
		std::cout << '\n' << i + 1;
		for (int j = 0; j < adjacency_list[i].size(); j++) {
		for (auto edge : adjacency_list[i])
			std::cout << "->" << "(" << edge.destination + 1 <<
				 ", " << edge.weight << ")";
		}
	}
}

inline void AdjacencyStructGraph<false>::print()
{
	for (int i = 0; i < adjacency_list.size(); i++) {
		std::cout << '\n' << i + 1;
		for (auto j : adjacency_list[i])
			std::cout << "->" << j + 1;
	}
}

inline std::vector<int> AdjacencyStructGraph<false>::getPathByBFS(size_t from, size_t to)
{
	// Filling extra vectors
	std::list<Node*> graph, buffer;

	for (size_t i = 0; i < adjacency_list.size(); i++)
		graph.push_back(new Node(i));

	for (size_t i = 0; i < adjacency_list.size(); i++)
		for (auto edge : adjacency_list[i]) {
			int edge_from_number = i,
				edge_to_number = edge;
			Node* edge_from = node_by_number(graph, edge_from_number);
			Node* edge_to = node_by_number(graph, edge_to_number);
			edge_from->nodes.push_back(edge_to);
			edge_to->nodes.push_back(edge_from);
		}

	buffer.push_back(node_by_number(graph, from));

	std::vector<int> path;
	if (bfs(graph, node_by_number(graph, to), buffer)) 
	{
		get_path(node_by_number(graph, from), node_by_number(graph, to), path);
	}

	for (typename std::list<Node*>::iterator it = graph.begin(); it != graph.end(); ++it)
		delete* it;
	
	return path;
}

inline std::vector<int> AdjacencyStructGraph<true>::getPathByBFS(size_t from, size_t to)
{
	// Filling extra vectors

	std::list<Node*> graph, buffer;

	for (size_t i = 0; i < adjacency_list.size(); i++)
		graph.push_back(new Node(Edge(i, 0)));

	for (size_t i = 0; i < adjacency_list.size(); i++)
		for (auto edge : adjacency_list[i]) {
			int edge_from_number = i,
				edge_to_number = edge.destination;
			Node* edge_from = node_by_number(graph, edge_from_number);
			Node* edge_to = node_by_number(graph, edge_to_number);
			edge_from->nodes.push_back(edge_to);
			edge_to->nodes.push_back(edge_from);
		}

	buffer.push_back(node_by_number(graph, from));

	std::vector<int> path;
	if (bfs(graph, node_by_number(graph, to), buffer)) {
		get_path(node_by_number(graph, from), node_by_number(graph, to), path);
	}

	for (typename std::list<Node*>::iterator it = graph.begin(); it != graph.end(); ++it)
		delete* it;

	return path;
}


template<bool isWeighted>
inline AdjacencyStructGraph<isWeighted> GenerateRandomAdjacencyStructGraph(size_t verticesNum, size_t edgeNum,
	bool isOriented, bool hasLoops, size_t maxWeight)
{
	AdjacencyStructGraph<isWeighted> graph(verticesNum, isOriented);

	std::random_device rd;
	std::mt19937 mersenne(rd());

	while (edgeNum) {
		size_t vertice_1 = mersenne() % verticesNum,
			vertice_2 = mersenne() % verticesNum,
			weight = isWeighted ? mersenne() % maxWeight + 1 : 1;
		if (!hasLoops && vertice_1 == vertice_2)
			continue;

		if (graph.AddEdge(vertice_1, vertice_2, weight))
			edgeNum--;
	}
	return graph;
}