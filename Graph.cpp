#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <unordered_set>
#include <unordered_map>

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

	bool isEuler();
	
	static std::vector<int> getMaxFlowHelperBFS(int start, int end, const std::vector<std::vector<int> >& capac, const std::vector<std::unordered_set<int> >& la); // parent vector

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

	std::vector<int> eulerCycle();

	static void getMaxFlowBipartite(int N1, int N2, std::vector < std::pair<int, int> >& edges, int& maxFlow);
	static bool HavelHakimi(const std::vector<int>& deg); // bool if y/n
	static void DisjointSet(int n_nodes, const std::vector<std::pair<int, std::pair<int, int> > >& info, std::ostream&); // "da" & "nu"
};

void Graph::addEdge(int from, int to, int cost, bool isDirected) {
	N_EDGES++;

	la[from].push_back({ to, cost });
	//li[to].push_back({ from, cost });
	//edges.push_back({ from, to, cost });

	if (!isDirected) {
		la[to].push_back({ from, cost });
		//li[from].push_back({ to, cost });
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

	for (int k = 1; k <= N_NODES; k++)
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

std::vector<int> Graph::getMaxFlowHelperBFS(int start, int end, const std::vector<std::vector<int> >& capac, const std::vector<std::unordered_set<int> >& la) {
	std::queue<int> q;
	std::vector<int> parent(la.size(), -1);
	std::vector<bool> visited(la.size(), false);

	q.push(start);
	visited[start] = 1;
	while (!q.empty()) {
		int n = q.front();
		q.pop();

		if (n == end)continue;

		for (int i : la[n])
			if (capac[n][i] > 0 && !visited[i]) {
				parent[i] = n;
				visited[i] = 1;
				q.push(i);
			}
	}

	if (!visited[end])return std::vector<int>();

	return parent;
}

int Graph::getMaxFlow(int start, int end) {
	const int INF = 1e9;

	std::vector<std::unordered_set<int> > la_flow;
	la_flow.resize(N_NODES + 1);  // new adj list both ways

	std::vector<std::vector<int> > capac;
	capac.resize(N_NODES + 1, std::vector<int>(N_NODES + 1, 0));

	int maxFlow = 0;

	for (int i = 1; i <= N_NODES; ++i) {
		for (const Neighbour& v : la[i]) {
			capac[i][v.node] = v.cost;

			// residual edge -> cost[v.node][i] = 0 

			la_flow[i].insert(v.node);
			la_flow[v.node].insert(i);
		}
	}


	while (1) {
		// get path from start to end
		std::vector<int> parent = getMaxFlowHelperBFS(start, end, capac, la_flow);

		// if no path -> exit
		if (parent.empty())break;

		for (int i : la_flow[end]) {
			int current_node = end;
			int flow = INF;
			parent[end] = i;

			// find bottleneck flow
			while (parent[current_node] != -1) {
				if (flow > capac[parent[current_node]][current_node])
					flow = capac[parent[current_node]][current_node];

				current_node = parent[current_node];
			}

			// flow == 0 -> saturated path
			// current_node != start when parent[end] not in bfs
			if (flow == 0 || current_node != start) continue;

			maxFlow += flow;

			// update capacity
			current_node = end;
			while (parent[current_node] != -1) {
				// normal edge
				capac[parent[current_node]][current_node] -= flow;

				// residual edge
				capac[current_node][parent[current_node]] += flow;

				current_node = parent[current_node];
			}
		}
	}
	return maxFlow;
}

void Graph::getMaxFlowBipartite(int N1, int N2, std::vector < std::pair<int, int> >& edges, int& maxFlow) {

	// acceasi rezolvare ca la getMaxFlow doar ca acum adaug si muchiile intr-un vector

	const int INF = 1e9;
	maxFlow = 0;

	const int source = N1 + N2 + 1;
	const int dest = N1 + N2 + 2;
	int total_nodes = N1 + N2 + 2;

	std::vector<std::vector<int> > la_flow;
	la_flow.resize(total_nodes + 1);

	// https://stackoverflow.com/questions/32685540/why-cant-i-compile-an-unordered-map-with-a-pair-as-key
	// hash unordered map 
	struct Pair_hash {
		size_t operator()(const std::pair<int,int>& a) const {
			return (1ULL * a.first * 698237239333) + a.second;
		}
	};

	std::unordered_map<std::pair<int, int>, int, Pair_hash> capac;

	// la/capac N1 -> N2
	for (const auto& p : edges) {
		int x = p.first;
		int y = p.second + N1;

		la_flow[x].push_back(y);
		la_flow[y].push_back(x);
		capac[{x, y}] = 1;
		capac[{y, x}] = 0;
	}

	// la/capac source -> N1
	for (int i = 1; i <= N1; ++i) {
		la_flow[source].push_back(i);
		la_flow[i].push_back(source);
		capac[{source, i}] = 1;
		capac[{i, source}] = 0;
	}

	// la/capac N2 -> dest
	for (int i = 1; i <= N2; ++i) {
		la_flow[i + N1].push_back(dest);
		la_flow[dest].push_back(i + N1);
		capac[{i + N1, dest}] = 1;
		capac[{dest, i + N1}] = 0;
	}

	// bfs
	while (1) {
		// get path from start to end
		std::vector<int> parent(total_nodes + 1, -1);
		std::vector<bool> used(total_nodes + 1, false);

		// bfs to find path source - dest
		std::queue<int> q;
		q.push(source);
		used[source] = 1;
		while (!q.empty()) {
			int n = q.front();
			q.pop();
			if (n == dest)continue;

			for (int i : la_flow[n]) {
				if (capac[{n, i}] == 1 && !used[i]) {
					parent[i] = n;
					used[i] = 1;
					q.push(i);
				}
			}
		}

		// if no path -> exit
		if (!used[dest])break;

		for (int i : la_flow[dest]) {
			if (!used[i])continue; // i not in bfs

			int current_node = dest;
			parent[dest] = i;

			// bottleneck flow 0 or 1
			int flow = 1;
			while (parent[current_node] != -1) {
				if (capac[{parent[current_node], current_node}] == 0) {
					flow = 0;
					break;
				}
				current_node = parent[current_node];
			}

			if (flow == 0)continue;

			maxFlow++;

			current_node = dest;
			// update capacity
			while (parent[current_node] != -1) {
				int& par = parent[current_node];

				capac[{par, current_node}] -= 1;
				capac[{current_node, par}] += 1;
				current_node = parent[current_node];
			}

		}
	}

	// return edges
	edges.erase(edges.begin(), edges.end());
	for (int i = 1; i <= N1; ++i)
		for (int j : la_flow[i])
			if (j != source && capac[{i, j}] == 0) {
				edges.push_back({ i, j - N1 });
			}
}

bool Graph::isEuler() {
	for (int i = 1; i <= N_NODES; ++i) {
		if (la[i].size() & 1)
			return false;
	}

	auto dist = BFS(1);
	for (int i = 1; i <= N_NODES; ++i)
		if (dist[i] == -1)
			return false;
 
	return true;
}

std::vector<int> Graph::eulerCycle() {
	if (!isEuler()) {
		throw std::runtime_error("-1");
	}

	using iPair = std::pair<int, int>;

	// https://stackoverflow.com/questions/32685540/why-cant-i-compile-an-unordered-map-with-a-pair-as-key
	// hash unordered map 
	struct Pair_hash {
		size_t operator()(const iPair& a) const {
			//return (std::hash<int>()(a.first) ^ (std::hash<int>()(a.second) >> 1));
			return (1ULL * a.first * 698237239333) + a.second;
		}
	};

	// count edges
	std::unordered_map<iPair, int, Pair_hash> edges;
	for (int i = 1; i <= N_NODES; ++i)
		for (const Neighbour& x : la[i])
			edges[{std::min(i,x.node), std::max(i,x.node)}]++;


	// dfs
	std::stack<int> s;
	std::vector<int> cycle;
	auto la_copy = la; // la copy
	cycle.reserve(N_EDGES);

	/* pseudocod pt ce e mai jos 
		dfs(i) = 
			for neighbour of i
				if edge[i,neighbour] not used 
					then dfs(neighbour)
			add i to cycle
	*/

	s.push(1);
	while (!s.empty()) {
		int top = s.top();
		s.pop();

		bool found_neighbour = false;

		// iau primul vecin disponibil
		while (la_copy[top].size() && !found_neighbour) {
			int i = la_copy[top].back().node;
			la_copy[top].pop_back();

			const int& p1 = std::min(top, i);
			const int& p2 = std::max(top, i);

			if (edges[{p1, p2}] > 0) {
				found_neighbour = true;
				edges[{p1, p2}] -= 2;

				s.push(top); // ca sa continui cu ceilalti vecini
				s.push(i);
			}
		}

		if (!found_neighbour) {
			cycle.push_back(top);
		}
	}
	return cycle;
}

// bfs https://infoarena.ro/job_detail/2789676
// dfs https://infoarena.ro/job_detail/2789680
// biconex = nu am facut
// ctc https://infoarena.ro/job_detail/2790880
// sortaret https://infoarena.ro/job_detail/2790888
// havel hakimi - da, metoda in clasa
// critical-connections-in-a-network = nu am facut
// apm https://infoarena.ro/job_detail/2799654
// disjoint https://infoarena.ro/job_detail/2793259
// dijkstra https://infoarena.ro/job_detail/2799709
// bellmanford https://infoarena.ro/job_detail/2807122
// maxflow https://infoarena.ro/job_detail/2810858
// royfloyd https://infoarena.ro/job_detail/2803900
// darb https://infoarena.ro/job_detail/2803907
// cuplaj 60 pct TLE https://infoarena.ro/job_detail/2810845
// ciclueuler 50 pct TLE https://www.infoarena.ro/job_detail/2818119
// hamilton nu am facut

int main() {
	std::ifstream f("cuplaj.in");
	std::ofstream g("cuplaj.out");
	int N1, N2, M;

	f >> N1 >> N2 >> M;

	std::vector<std::pair<int, int> >edges;
	edges.reserve(M);

	for (int i = 0; i < M; ++i) {
		int x, y;
		f >> x >> y;
		edges.push_back({ x,y });
	}

	int maxflow = 0;
	Graph::getMaxFlowBipartite(N1, N2, edges, maxflow);
	g << maxflow << '\n';
	for (auto& p : edges) {
		g << p.first << ' ' << p.second << '\n';
	}
}