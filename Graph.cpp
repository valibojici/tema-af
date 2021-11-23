#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>

class Graph {
private:
	struct Edge {
		int x, y, cost;
		Edge(int x, int y, int cost = 0) : x(x), y(y), cost(cost) {}
	};

	struct Neighbour {
		int node, cost;
		Neighbour(int x, int cost = 0) : node(x), cost(cost) {}
	};

	std::vector<std::vector<Neighbour> >la;
	std::vector<std::vector<Neighbour> >li;
	std::vector<Edge> edges;

	const int N_NODES;
	int N_EDGES;

	void KosarajuTopSortHelper(std::list<int>& nodes, int n, std::vector<bool>& used);
	void KosarajuAddToComponent(int n, std::vector<bool>& used, std::vector<int>& component);
	std::vector<int> getMaxFlowHelperBFS(int start, int end, const std::vector<std::vector<int> >& cost); // parent vector

	static int Find(int x, std::vector<int>& parent);
	static void Union(int x, int y, std::vector<int>& parent, std::vector<int>& height);
public:
	Graph(int N) : N_NODES(N), N_EDGES(0) { la.resize(N + 1); li.resize(N + 1); }
	void getInfo();
	void addEdge(int from, int to, int cost = 0, bool isDirected = 0);

	std::vector<int> BFS(int start); // distance from start 
	int getComponentCount();
	std::vector<int> getTopologicalSort();
	std::vector<std::vector<int> > Kosaraju(); // vector of SCP

	std::pair<int, std::vector<int> > Prim(); // cmin & parent vector
	std::vector<int> Dijkstra(int start);
	std::vector<int> BellmanFord(int start);

	int MaxDist2Nodes();
	std::vector<std::vector<int> > RoyFloyd();
	int getMaxFlow(int start, int end);
	 

	static bool HavelHakimi(const std::vector<int>& deg); // bool if y/n
	static void DisjointSet(int n_nodes, const std::vector<std::pair<int, std::pair<int, int> > >& info, std::ostream&); // "da" & "nu"
};

void Graph::addEdge(int from, int to, int cost, bool isDirected) {
	N_EDGES++;

	la[from].push_back({ to, cost });
	li[to].push_back({ from, cost });
	edges.push_back({ from, to, cost });

	if (!isDirected) {
		la[to].push_back({ from, cost });
		li[from].push_back({ to, cost });
	}
}

std::vector<int> Graph::BFS(int start) {
	std::vector<int> dist(N_NODES + 1, -1);
	dist[start] = 0;

	std::queue<int> q;
	q.push(start);

	while (!q.empty()) {
		int top = q.front();
		q.pop();

		for (const Neighbour& x : la[top]) {
			if (dist[x.node] == -1) {
				dist[x.node] = dist[top] + 1;
				q.push(x.node);
			}
		}
	}
	return dist;
}

int Graph::getComponentCount() {
	std::stack<int> st;
	std::vector<bool> used(N_NODES + 1, false);
	int count = 0;

	for (int i = 1; i <= N_NODES; ++i) {
		if (!used[i]) {
			count++;

			// DFS
			st.push(i);
			while (!st.empty()) {
				int top = st.top();
				st.pop();
				used[top] = 1;

				for (const Neighbour& n : la[top])
					if (!used[n.node]) st.push(n.node);
			}
		}
	}

	return count;
}

void Graph::getInfo() {
	std::cout << "N: " << N_NODES << '\n';
	std::cout << "M: " << N_EDGES << '\n';
	for (int i = 1; i <= N_NODES; ++i) {
		std::cout << i << ": ";
		for (const Neighbour& x : la[i]) {
			std::cout << "(" << x.node << ", " << x.cost << ") ";
		}
		std::cout << '\n';
	}

	std::cout << "Edges: \n";
	for (const Edge& m : edges) {
		std::cout << m.x << " " << m.y << " " << m.cost << '\n';
	}
}

std::vector<int> Graph::getTopologicalSort() {
	std::vector<int> in_deg(N_NODES + 1, 0);
	std::vector<int> result;
	result.reserve(N_NODES);

	for (int i = 1; i <= N_NODES; ++i) {
		for (const Neighbour& x : la[i]) {
			in_deg[x.node]++;
		}
	}

	std::queue<int> q;
	for (int i = 1; i <= N_NODES; ++i) {
		if (in_deg[i] == 0)q.push(i);
	}

	while (!q.empty()) {
		int n = q.front();
		q.pop();

		result.push_back(n);
		for (const Neighbour& x : la[n]) {
			in_deg[x.node]--;
			if (in_deg[x.node] == 0) {
				q.push(x.node);
			}
		}
	}
	return result.size() == N_NODES ? result : std::vector<int>(); // if cycle => empty
}

void Graph::KosarajuTopSortHelper(std::list<int>& nodes, int n, std::vector<bool>& used) {
	used[n] = 1;
	for (const Neighbour& x : la[n]) {
		if (!used[x.node])KosarajuTopSortHelper(nodes, x.node, used);
	}
	nodes.push_front(n);
}

void Graph::KosarajuAddToComponent(int n, std::vector<bool>& used, std::vector <int>& comp) {
	comp.push_back(n);
	used[n] = 1;

	for (const Neighbour& x : li[n]) {
		if (!used[x.node])
			KosarajuAddToComponent(x.node, used, comp);
	}
}

std::vector<std::vector<int> > Graph::Kosaraju() {

	// part 1 (topological sort-ish)
	std::list<int> nodes;
	std::vector<bool> used(N_NODES + 1, false);
	for (int i = 1; i <= N_NODES; ++i) {
		if (!used[i])
			KosarajuTopSortHelper(nodes, i, used);
	}

	// part 2
	std::fill(used.begin(), used.end(), false);
	std::vector<std::vector<int> > result;

	for (int n : nodes) {
		if (!used[n]) {
			std::vector<int> component;
			KosarajuAddToComponent(n, used, component);
			result.push_back(component);
		}
	}

	return result;
}

bool Graph::HavelHakimi(const std::vector<int>& deg) {
	std::list<int> degrees;

	int sum = 0;
	int N = deg.size();
	for (int i : deg) {
		sum += i;
		if (i >= N)return false; // max deg = N - 1
		degrees.push_back(i);
	}

	if (sum & 1)return false; // M = sum deg / 2
	if (sum / 2 > N * (N - 1) / 2)return false; // M <= N * (N-1) / 2

	degrees.sort([](int x, int y) { return x > y; });

	while (!degrees.empty()) {
		int d = degrees.front();
		degrees.pop_front();

		for (int& i : degrees) {
			if (d == 0)break;
			if (i - 1 < 0) return false;
			i--;
			d--;
		}
		degrees.sort([](int x, int y) { return x > y; });
	}

	return 1;
}

std::pair<int, std::vector<int> > Graph::Prim() {
	const int MAX_COST = 1e9;

	std::vector<int> cost(N_NODES + 1, MAX_COST);
	std::vector<bool> used(N_NODES + 1, false);
	std::vector<int> parent(N_NODES + 1, -1);

	struct Compare {
		bool operator()(const Neighbour& a, const Neighbour& b) { return a.cost > b.cost; }
	};

	std::priority_queue<Neighbour, std::vector<Neighbour>, Compare> pq;

	cost[1] = 0;
	int cmin = 0;
	pq.push({ 1, 0 });

	while (!pq.empty()) {
		int n = pq.top().node;
		pq.pop();

		if (used[n])continue;

		used[n] = 1;
		cmin += cost[n];

		for (const Neighbour& x : la[n]) {
			if (!used[x.node] && x.cost < cost[x.node]) {
				cost[x.node] = x.cost;
				parent[x.node] = n;
				pq.push({ x.node, cost[x.node] });
			}
		}
	}
	return { cmin, parent };
}

int Graph::Find(int x, std::vector<int>& parent) {
	int r = x;
	while (parent[r] != r) r = parent[r];
	while (x != r) {
		int temp = parent[x];
		parent[x] = r;
		x = temp;
	}
	return r;
}

void Graph::Union(int x, int y, std::vector<int>& parent, std::vector<int>& height) {
	x = Graph::Find(x, parent);
	y = Graph::Find(y, parent);
	if (x == y)return;
	if (height[x] < height[y])
		parent[x] = y;
	else {
		parent[y] = x;
		if (height[x] == height[y]) height[x]++;
	}
}

void Graph::DisjointSet(int n_nodes, const std::vector<std::pair<int, std::pair<int, int> > >& info, std::ostream& out) {
	std::vector<int> parent(n_nodes + 1);
	std::vector<int> height(n_nodes + 1, 0);
	for (int i = 1; i <= n_nodes; ++i)
		parent[i] = i;

	for (const auto& p : info) {
		int x = p.second.first;
		int y = p.second.second;

		if (p.first == 1) { // Union
			Graph::Union(x, y, parent, height);
		}
		else if (p.first == 2) { // Find
			out << (Graph::Find(x, parent) == Graph::Find(y, parent) ? "DA\n" : "NU\n");
		}
	}
}

std::vector<int> Graph::Dijkstra(int start) {
	const int INF = -1;

	std::vector<int> dist(N_NODES + 1, INF);
	std::vector<bool> used(N_NODES + 1, false);

	struct Compare {
		bool operator()(const Neighbour& a, const Neighbour& b) { return a.cost > b.cost; }
	};

	std::priority_queue<Neighbour, std::vector<Neighbour>, Compare> pq;

	dist[start] = 0;
	pq.push({ start, 0 });

	while (!pq.empty()) {
		int n = pq.top().node;
		pq.pop();

		if (used[n])continue;
		used[n] = 1;

		for (const Neighbour& x : la[n]) {
			if (used[x.node]) continue;

			if (dist[x.node] == INF || (!used[x.node] && dist[n] + x.cost < dist[x.node])) {
				dist[x.node] = dist[n] + x.cost;
				pq.push({ x.node, dist[x.node] });
			}
		}
	}
	for (int& i : dist)
		if (i == INF) i = 0;

	return dist;
}

std::vector<int> Graph::BellmanFord(int start) {
	const int MAX_COST = 2e9;

	std::vector<int> nodes_to_check;
	nodes_to_check.push_back(start);

	std::vector<bool> used(N_NODES + 1, false); // used[i] = 1 ==> i in modified_nodes
	std::vector<int> dist(N_NODES + 1, MAX_COST);
	dist[start] = 0;

	for (int i = 1; i <= N_NODES - 1 && !nodes_to_check.empty(); ++i) {

		std::fill(used.begin(), used.end(), false);

		std::vector<int> modified_nodes;

		for (const int& n : nodes_to_check) {
			for (const Neighbour& x : la[n])
				if (dist[x.node] == MAX_COST || dist[n] + x.cost < dist[x.node]) {
					dist[x.node] = dist[n] + x.cost;

					if (!used[x.node]) // if x.node not in modified_nodes
					{
						modified_nodes.push_back(x.node);
						used[x.node] = 1;
					}
				}
		}
		nodes_to_check.swap(modified_nodes);
	}

	//check neg cycle
	for (int i = 1; i <= N_NODES; ++i)
		for (const Neighbour& x : la[i])
			if (dist[i] + x.cost < dist[x.node])
				throw std::runtime_error("Ciclu negativ!");

	return dist;
}

int Graph::MaxDist2Nodes() {
	std::vector<int> dist = BFS(1);
	int idx_dMax = 1;
	for (int i = 2; i <= N_NODES; ++i)
		if (dist[i] > dist[idx_dMax]) idx_dMax = i;

	dist = BFS(idx_dMax);
	
	int dMax = 0;
	for (int i : dist)
		if (dMax < i) dMax = i;

	return dMax + 1;
}

std::vector<std::vector<int> > Graph::RoyFloyd() {
	const int INF = 1e9;

	std::vector<std::vector<int> > mat;
	mat.resize(N_NODES + 1, std::vector<int>(N_NODES + 1, INF));

	for (int i = 1; i <= N_NODES; ++i)
		for (const Neighbour& x : la[i])
			mat[i][x.node] = x.cost;

	for(int k=1;k<=N_NODES;k++)
		for (int i = 1; i <= N_NODES; ++i) {
			if (i == k)continue;

			for (int j = 1; j <= N_NODES; ++j) {
				if (j == i || j == k) continue;

				if (mat[i][k] + mat[k][j] < mat[i][j]) 
					mat[i][j] = mat[i][k] + mat[k][j];
			}
		}

	for (auto& v : mat) {
		for (int& i : v)
			if (i == INF) i = 0;
	}
	return mat;
}

std::vector<int> Graph::getMaxFlowHelperBFS(int start, int end, const std::vector<std::vector<int> >& cost) {
	std::queue<int> q;
	std::vector<int> parent(N_NODES + 1, -1);
	std::vector<bool> visited(N_NODES + 1, false);
	bool done = false;
	
	q.push(start);
	visited[start] = 1;
	while (!q.empty() && !done) {
		int n = q.front();
		q.pop();
 
		for (int i = 1; i <= N_NODES && !done; ++i)
			if (cost[n][i] > 0 && !visited[i]) {
				parent[i] = n;

				if (i == end)done = true;

				visited[i] = 1;
				q.push(i);
			}
	}

	if (!done) return std::vector<int>();

	return parent;
}

int Graph::getMaxFlow(int start, int end) {
	const int MAX_FLOW = 1e9;

	std::vector<std::vector<int> > cost;
	cost.resize(N_NODES + 1, std::vector<int>(N_NODES + 1, 0));

	int maxFlow = 0;

	for (int i = 1; i <= N_NODES; ++i) {
		for (const Neighbour& v : la[i]) {
			cost[i][v.node] = v.cost;

			// residual edge -> cost[v.node][i] = 0 
		}
	}

	while (1) {
		// get path from start to end
		std::vector<int> parent = getMaxFlowHelperBFS(start, end, cost);

		// if no path -> exit
		if (parent.empty())break;

		int current_node = end;
		int flow = MAX_FLOW;

		// find bottleneck flow in path
		while (parent[current_node] != -1) {
			if (flow > cost[parent[current_node]][current_node])
				flow = cost[parent[current_node]][current_node];

			current_node = parent[current_node];
		}

		// update global max flow
		maxFlow += flow;

		current_node = end;
		// update cost 
		while (parent[current_node] != -1) {
			// out edge
			cost[parent[current_node]][current_node] -= flow;

			// residual edge
			cost[current_node][parent[current_node]] += flow;

			current_node = parent[current_node];
		}
	}
	return maxFlow;
}

int main() {

	std::ifstream f("maxflow.in");
	std::ofstream g("maxflow.out");

	int N, M;
	f >> N >> M;
	Graph a(N);

	for (int i = 0; i < M; ++i) {
		int x, y, c;
		f >> x >> y >> c;
		a.addEdge(x, y, c, 1);
	}
	
	g << a.getMaxFlow(1, N);
}